//========================================================================================
// Athena++ astrophysical MHD code, Kokkos version
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particle_random.cpp
//! \brief Problem generator that initializes random particle positions and velocities.

#include <algorithm>
#include <cmath>
#include <sstream>
#include <iostream>

#include "parameter_input.hpp"
#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "particles/particles.hpp"
#include <Kokkos_Random.hpp>

// Athena++ headers
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"



// User-defined history functions
void ParticleHistory(HistoryData *pdata, Mesh *pm);

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem_()
//! \brief Problem Generator for random particle positions/velocities

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  //if (restart) return;


  // enroll user history function
  user_hist_func = ParticleHistory;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->ppart == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Random particles test requires <particles> block in input file"
              << std::endl;
    exit(EXIT_FAILURE);
  }


  Real cfl_part  = pin->GetOrAddReal("particles","cfl_part",0.05);
  Real B0z  = pin->GetOrAddReal("problem","B0z",1.0);

  Real min_mass   = pin->GetOrAddReal("particles","min_mass", 1.0);
  Real mass_log_spacing   = pin->GetOrAddReal("particles","mass_log_spacing", 1.0);


  // capture variables for the kernel
  auto &mbsize = pmy_mesh_->pmb_pack->pmb->mb_size;
  auto &pr = pmy_mesh_->pmb_pack->ppart->prtcl_rdata;
  auto &pi = pmy_mesh_->pmb_pack->ppart->prtcl_idata;
  auto &npart = pmy_mesh_->pmb_pack->ppart->nprtcl_thispack;
  auto &npart_spec = pmy_mesh_->pmb_pack->ppart->nprtcl_perspec_thispack;
  auto gids = pmy_mesh_->pmb_pack->gids;
  auto gide = pmy_mesh_->pmb_pack->gide;
  auto &nspecies_ = pmy_mesh_->pmb_pack->ppart->nspecies;
  Real x1size = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  Real x2size = pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min;
  Real x3size = pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min;

  auto &nmbtp = pmbp->nmb_thispack;

  // initialize particles
  Kokkos::Random_XorShift64_Pool<> rand_pool64(pmbp->gids);
  par_for("part_update",DevExeSpace(),0,(npart-1),
  KOKKOS_LAMBDA(const int p) {
    auto rand_gen = rand_pool64.get_state();  // get random number state this thread
    // choose parent MeshBlock randomly
    int m = static_cast<int>(rand_gen.frand()*(gide - gids + 1.0));
    m = max(m,0);
    m = min(m,nmbtp-1 );
    pi(PGID,p) = gids + m;
    int spec = p /npart_spec;
    spec = max(spec, 0);
    spec = min(spec, nspecies_ -1);
    pi(PSP,p) = spec;
    Real rand = rand_gen.frand();
    pr(IPX,p) = (1. - rand)*mbsize.d_view(m).x1min + rand*mbsize.d_view(m).x1max;
    pr(IPX,p) = fmin(pr(IPX,p),mbsize.d_view(m).x1max);
    pr(IPX,p) = fmax(pr(IPX,p),mbsize.d_view(m).x1min);
    rand = rand_gen.frand();
    pr(IPY,p) = (1. - rand)*mbsize.d_view(m).x2min + rand*mbsize.d_view(m).x2max;
    pr(IPY,p) = fmin(pr(IPY,p),mbsize.d_view(m).x2max);
    pr(IPY,p) = fmax(pr(IPY,p),mbsize.d_view(m).x2min);

    rand = rand_gen.frand();
    pr(IPZ,p) = (1. - rand)*mbsize.d_view(m).x3min + rand*mbsize.d_view(m).x3max;
    pr(IPZ,p) = fmin(pr(IPZ,p),mbsize.d_view(m).x3max);
    pr(IPZ,p) = fmax(pr(IPZ,p),mbsize.d_view(m).x3min);

    Real mu = 0.99; //2.0*(rand_gen.frand()-0.5);
    Real phi = 2.0*M_PI * rand_gen.frand();
    pr(IPVX,p) = std::sqrt(1.0-mu*mu) * std::cos(phi);
    pr(IPVY,p) = std::sqrt(1.0-mu*mu) * std::sin(phi);
    pr(IPVZ,p) = mu;
    pr(IPM,p) = min_mass * pow(mass_log_spacing, spec   );
    pr(IPDX,p) = 0.0;
    pr(IPDY,p) = 0.0;
    pr(IPDZ,p) = 0.0;
    pr(IPDB,p) = 0.0;
    // Fix below so that DF & track at first timestep is wrt to real B field. 
    // Only an issue for the very first output
    pr(IPBX,p) = 0.0;
    pr(IPBY,p) = 0.0;
    pr(IPBZ,p) = 1.0; 
    rand_pool64.free_state(rand_gen);
  });

  // set timestep (which will remain constant for entire run
  // Assumes uniform mesh (no SMR or AMR)
  // Assumes velocities normalized to one, so dt=min(dx)
  Real &dtnew_ = pmy_mesh_->pmb_pack->ppart->dtnew;
  dtnew_ = std::min(mbsize.h_view(0).dx1, mbsize.h_view(0).dx2);
  dtnew_ = std::min(dtnew_, mbsize.h_view(0).dx3);
  dtnew_ *= pin->GetOrAddReal("time", "cfl_number", 0.8);


  if (restart){
    auto &bcc = pmbp->pmhd->bcc0;
    auto &b = pmbp->pmhd->b0;

    // capture variables for kernel
    auto &indcs = pmy_mesh_->mb_indcs;
    int &is = indcs.is; int &ie = indcs.ie;
    int &js = indcs.js; int &je = indcs.je;
    int &ks = indcs.ks; int &ke = indcs.ke;
    int nx1 = indcs.nx1;
    int nx2 = indcs.nx2;
    int nx3 = indcs.nx3;

    const int nmkji = (pmy_mesh_->pmb_pack->nmb_thispack)*nx3*nx2*nx1;
    const int nkji = nx3*nx2*nx1;
    const int nji  = nx2*nx1;
    Real dtOmega = std::numeric_limits<float>::max();

    Kokkos::parallel_reduce("pgen_restart_w_part",Kokkos::RangePolicy<>(DevExeSpace(),
                            0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &min_dt) {
    // compute m,k,j,i indices of thread and call function
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;
      Real Btemp = std::sqrt(  SQR( b.x1f(m,k,j,i) )
                             + SQR( b.x2f(m,k,j,i) )
                             + SQR( b.x3f(m,k,j,i) )  );
      min_dt = fmin(cfl_part*min_mass/Btemp, min_dt);
    }, Kokkos::Min<Real>(dtOmega));

    dtnew_ = std::min(dtnew_, dtOmega);
    return;
  } else{
    // initialize MHD variables ----------------------------------------------------------
    if (pmbp->pmhd != nullptr) {
      EOS_Data &eos = pmbp->pmhd->peos->eos_data;
      Real gm1 = eos.gamma - 1.0;
      Real p0 = 1.0/eos.gamma;
      auto &u0 = pmbp->pmhd->u0;
      auto &b0 = pmbp->pmhd->b0;
      // capture variables for kernel
      auto &indcs = pmy_mesh_->mb_indcs;
      int &is = indcs.is; int &ie = indcs.ie;
      int &js = indcs.js; int &je = indcs.je;
      int &ks = indcs.ks; int &ke = indcs.ke;
      int nx1 = indcs.nx1;
      int nx2 = indcs.nx2;
      int nx3 = indcs.nx3;
      par_for("pgen_mhd", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        // compute cell-centered conserved variables
        u0(m,IDN,k,j,i) = 1.0;
        u0(m,IM1,k,j,i) = 0.0;
        u0(m,IM2,k,j,i) = 0.0;
        u0(m,IM3,k,j,i) = 0.0;

        b0.x1f(m,k,j,i) = 0.0;
        b0.x2f(m,k,j,i) = 0.0;
        b0.x3f(m,k,j,i) = B0z;

        // Include extra face-component at edge of block in each direction
        if (i==ie) {
          b0.x1f(m,k,j,i+1) = 0.0;
        }
        if (j==je) {
          b0.x2f(m,k,j+1,i) = 0.0;
        }
        if (k==ke) {
          b0.x3f(m,k+1,j,i) = B0z;
        }
        if (eos.is_ideal) {
          u0(m,IEN,k,j,i) = p0/gm1 + 0.5*B0z*B0z + // fix contribution from dB
             0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
             SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
        }
      });
    }  // End initialization MHD variables
    dtnew_ = std::min(dtnew_, cfl_part*min_mass / B0z);
    return;
  }
}

// User defined history output for particle data
// Outputs running diffusion coefficients using all particles

void ParticleHistory(HistoryData *pdata, Mesh *pm) {

  particles::Particles *pp = pm->pmb_pack->ppart;
  //int npart = pm->nprtcl_thisrank;
  auto &npart = pm->pmb_pack->ppart->nprtcl_thispack;
  auto &pr = pp->prtcl_rdata;
  auto &pi = pp->prtcl_idata;
  auto nspecies_ = pp->nspecies;
  auto nprtcl_total_ = pm->nprtcl_total;
  int nfields = 4;
  pdata->nhist = nfields*nspecies_;
  int &nhist_ = pdata->nhist;

  for(int i=0; i < nspecies_; ++i){
    pdata->label[nfields*i+0] = "Dx^2";
    pdata->label[nfields*i+1] = "Dy^2";
    pdata->label[nfields*i+2] = "Dz^2";
    pdata->label[nfields*i+3] = "Db^2";
    //pdata->label[nfields*i+4] = "<mu>";
  }

  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("HistSums",Kokkos::RangePolicy<>(DevExeSpace(), 0, npart-1),
  KOKKOS_LAMBDA(const int &p, array_sum::GlobalSum &mb_sum) {
    // MHD conserved variables:
    array_sum::GlobalSum hvars;
    int spec = pi(PSP,p);
    // fill the_array with zeros , if nhist < NHISTORY_VARIABLES
    for (int n=0; n<NHISTORY_VARIABLES; ++n) {
      hvars.the_array[n] = 0.0;
    } 
    // Now fill the fields relevant to that species   
    hvars.the_array[spec*nfields+0] = SQR(pr(IPDX, p))/nprtcl_total_ * nspecies_;
    hvars.the_array[spec*nfields+1] = SQR(pr(IPDY, p))/nprtcl_total_ * nspecies_;
    hvars.the_array[spec*nfields+2] = SQR(pr(IPDZ, p))/nprtcl_total_ * nspecies_;
    hvars.the_array[spec*nfields+3] = SQR(pr(IPDB, p))/nprtcl_total_ * nspecies_;

    // sum into parallel reduce
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  Kokkos::fence();

  // store data into hdata array
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }

  return;
}


