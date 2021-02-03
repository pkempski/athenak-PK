//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_calc_fluxes.cpp
//  \brief Calculate fluxes of the conserved variables, and electro-motive EMFs, for mhd

#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "mhd.hpp"
#include "eos/eos.hpp"
// include inlined reconstruction methods (yuck...)
#include "reconstruct/dc.cpp"
#include "reconstruct/plm.cpp"
#include "reconstruct/ppm.cpp"
// include inlined Riemann solvers (double yuck...)
#include "mhd/rsolvers/advect_mhd.cpp"
#include "mhd/rsolvers/llf_mhd.cpp"
//#include "mhd/rsolvers/hlld.cpp"
//#include "mhd/rsolvers/roe_mhd.cpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn  void MHD::MHDCalcFlux
//  \brief Calculate fluxes of conserved variables, and face-centered area-averaged EMFs
//  for evolution of magnetic field

TaskStatus MHD::MHDCalcFlux(Driver *pdrive, int stage)
{
  int is = pmy_pack->mb_cells.is; int ie = pmy_pack->mb_cells.ie;
  int js = pmy_pack->mb_cells.js; int je = pmy_pack->mb_cells.je;
  int ks = pmy_pack->mb_cells.ks; int ke = pmy_pack->mb_cells.ke;
  int ncells1 = pmy_pack->mb_cells.nx1 + 2*(pmy_pack->mb_cells.ng);

  int nmhd = nmhd;
  int nvars = nmhd + nscalars;
  int nmb = pmy_pack->nmb_thispack;
  auto recon_method = recon_method_;
  auto rsolver_method = rsolver_method_;
  auto &w0_ = w0;
  auto &b0_ = bcc0;
  auto &eos = peos->eos_data;

  //--------------------------------------------------------------------------------------
  // i-direction

  size_t scr_size = (ScrArray2D<Real>::shmem_size(nvars, ncells1) +
                     ScrArray2D<Real>::shmem_size(3, ncells1)) * 2;
  int scr_level = 0;
  auto flx1 = flux1;
  auto emf1 = emf_x1;
  auto bx = b0.x1f;

  par_for_outer("mhd_flux_x1",DevExeSpace(),scr_size,scr_level,0,(nmb-1), ks, ke, js, je,
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j)
    {
      ScrArray2D<Real> wl(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> wr(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> bl(member.team_scratch(scr_level), 3, ncells1);
      ScrArray2D<Real> br(member.team_scratch(scr_level), 3, ncells1);

      // Reconstruct qR[i] and qL[i+1], for both W and Bcc
      switch (recon_method)
      {
        case ReconstructionMethod::dc:
          DonorCellX1(member, m, k, j, is-1, ie+1, w0_, wl, wr);
          DonorCellX1(member, m, k, j, is-1, ie+1, b0_, bl, br);
          break;
        case ReconstructionMethod::plm:
          PiecewiseLinearX1(member, m, k, j, is-1, ie+1, w0_, wl, wr);
          PiecewiseLinearX1(member, m, k, j, is-1, ie+1, b0_, bl, br);
          break;
        case ReconstructionMethod::ppm:
          PiecewiseParabolicX1(member, m, k, j, is-1, ie+1, w0_, wl, wr);
          PiecewiseParabolicX1(member, m, k, j, is-1, ie+1, b0_, bl, br);
          break;
        default:
          break;
      }
      // Sync all threads in the team so that scratch memory is consistent
      member.team_barrier();

      // compute fluxes over [is,ie+1]
      switch (rsolver_method)
      {
        case MHD_RSolver::advect:
          Advect(member, eos, m, k, j, is, ie+1, IVX, wl, wr, bl, br, bx, flx1, emf1);
          break;
        case MHD_RSolver::llf:
          LLF(member, eos, m, k, j, is, ie+1, IVX, wl, wr, bl, br, bx, flx1, emf1);
          break;
//        case MHD_RSolver::hllc:
//          HLLC(member, eos, is, ie+1, IVX, wl, wr, uflux);
//          break;
//        case MHD_RSolver::roe:
//          Roe(member, eos, is, ie+1, IVX, wl, wr, uflux);
//          break;
        default:
          break;
      }

    }
  );
  if (!(pmy_pack->pmesh->nx2gt1)) return TaskStatus::complete;

  //--------------------------------------------------------------------------------------
  // j-direction

  scr_size = (ScrArray2D<Real>::shmem_size(nvars, ncells1) +
              ScrArray2D<Real>::shmem_size(3, ncells1)) * 3;
  auto flx2 = flux2;
  auto emf2 = emf_x2;
  auto by = b0.x2f;

  par_for_outer("mhd_flux2",DevExeSpace(),scr_size,scr_level,0,(nmb-1), ks, ke,
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k)
    {
      ScrArray2D<Real> scr1(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr2(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr3(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr4(member.team_scratch(scr_level), 3, ncells1);
      ScrArray2D<Real> scr5(member.team_scratch(scr_level), 3, ncells1);
      ScrArray2D<Real> scr6(member.team_scratch(scr_level), 3, ncells1);

      for (int j=js-1; j<=je+1; ++j) {
        // Permute scratch arrays.
        auto wl     = scr1;
        auto wl_jp1 = scr2;
        auto wr     = scr3;
        auto bl     = scr4;
        auto bl_jp1 = scr5;
        auto br     = scr6;
        if ((j%2) == 0) {
          wl     = scr2;
          wl_jp1 = scr1;
          bl     = scr5;
          bl_jp1 = scr4;
        }

        // Reconstruct qR[j] and qL[j+1], for both W and Bcc
        switch (recon_method)
        {
          case ReconstructionMethod::dc:
            DonorCellX2(member, m,k,j,is-1,ie+1, w0_, wl_jp1, wr);
            DonorCellX2(member, m,k,j,is-1,ie+1, b0_, bl_jp1, br);
            break;
          case ReconstructionMethod::plm:
            PiecewiseLinearX2(member, m,k,j,is-1,ie+1, w0_, wl_jp1, wr);
            PiecewiseLinearX2(member, m,k,j,is-1,ie+1, b0_, bl_jp1, br);
            break;
          case ReconstructionMethod::ppm:
            PiecewiseParabolicX2(member, m,k,j,is-1,ie+1, w0_, wl_jp1, wr);
            PiecewiseParabolicX2(member, m,k,j,is-1,ie+1, b0_, bl_jp1, br);
            break;
          default:
            break;
        }
        member.team_barrier();

        // compute fluxes over [js,je+1].
        if (j>(js-1)) {
          switch (rsolver_method)
          {
            case MHD_RSolver::advect:
              Advect(member, eos, m,k,j,is-1,ie+1, IVY, wl, wr, bl, br, by, flx2, emf2);
              break;
            case MHD_RSolver::llf:
              LLF(member, eos, m,k,j,is-1,ie+1, IVY, wl, wr, bl, br, by, flx2, emf2);
              break;
//            case MHD_RSolver::hllc:
//              HLLC(member, eos, is, ie, IVY, wl, wr, uf);
//              break;
//            case MHD_RSolver::roe:
//              Roe(member, eos, is, ie, IVY, wl, wr, uf);
//              break;
            default:
              break;
          }
        }

      } // end of loop over j
    }
  );
  if (!(pmy_pack->pmesh->nx3gt1)) return TaskStatus::complete;

  //--------------------------------------------------------------------------------------
  // k-direction. Note order of k,j loops switched

  scr_size = (ScrArray2D<Real>::shmem_size(nvars, ncells1) +
              ScrArray2D<Real>::shmem_size(3, ncells1)) * 3;
  auto flx3 = flux3;
  auto emf3 = emf_x3;
  auto bz = b0.x3f;

  par_for_outer("divflux_x3",DevExeSpace(), scr_size, scr_level, 0, (nmb-1), js, je,
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int j)
    {
      ScrArray2D<Real> scr1(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr2(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr3(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr4(member.team_scratch(scr_level), 3, ncells1);
      ScrArray2D<Real> scr5(member.team_scratch(scr_level), 3, ncells1);
      ScrArray2D<Real> scr6(member.team_scratch(scr_level), 3, ncells1);

      for (int k=ks-1; k<=ke+1; ++k) {
        // Permute scratch arrays.
        auto wl     = scr1;
        auto wl_kp1 = scr2;
        auto wr     = scr3;
        auto bl     = scr4;
        auto bl_kp1 = scr5;
        auto br     = scr6;
        if ((k%2) == 0) {
          wl     = scr2;
          wl_kp1 = scr1;
          bl     = scr5;
          bl_kp1 = scr4;
        }

        // Reconstruct qR[k] and qL[k+1], for both W and Bcc
        switch (recon_method)
        {
          case ReconstructionMethod::dc:
            DonorCellX3(member, m,k,j,is-1,ie+1, w0_, wl_kp1, wr);
            DonorCellX3(member, m,k,j,is-1,ie+1, b0_, bl_kp1, br);
            break;
          case ReconstructionMethod::plm:
            PiecewiseLinearX3(member, m,k,j,is-1,ie+1, w0_, wl_kp1, wr);
            PiecewiseLinearX3(member, m,k,j,is-1,ie+1, b0_, bl_kp1, br);
            break;
          case ReconstructionMethod::ppm:
            PiecewiseParabolicX3(member, m, k, j, is-1, ie+1, w0_, wl_kp1, wr);
            PiecewiseParabolicX3(member, m, k, j, is-1, ie+1, b0_, bl_kp1, br);
            break;
          default:
            break;
        }
        member.team_barrier();

        // compute fluxes over [ks,ke+1].  RS returns flux in input wr array
        if (k>(ks-1)) {
          switch (rsolver_method)
          {
            case MHD_RSolver::advect:
              Advect(member, eos, m,k,j,is-1,ie+1, IVZ, wl, wr, bl, br, bz, flx3, emf3);
              break;
            case MHD_RSolver::llf:
              LLF(member, eos, m,k,j,is-1,ie+1, IVZ, wl, wr, bl, br, bz, flx3, emf3);
              break;
//            case MHD_RSolver::hllc:
//              HLLC(member, eos, is, ie, IVZ, wl, wr, uf);
//              break;
//            case MHD_RSolver::roe:
//              Roe(member, eos, is, ie, IVZ, wl, wr, uf);
//              break;
            default:
              break;
          }
          member.team_barrier();
        }
      } // end loop over k
    }
  );
  return TaskStatus::complete;
}

} // namespace mhd