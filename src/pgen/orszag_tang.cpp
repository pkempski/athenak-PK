//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file orszag-tang.c
//  \brief Problem generator for Orszag-Tang vortex problem.
//
// REFERENCE: For example, see: G. Toth,  "The div(B)=0 constraint in shock capturing
//   MHD codes", JCP, 161, 605 (2000)
//========================================================================================

// C++ headers
#include <math.h> 
#include <iostream>   // endl
#include <sstream>    // stringstream

// Athena++ headers
#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "utils/grid_locations.hpp"
#include "pgen.hpp"

// global varibale shared with vector potential function
Real B0;

//----------------------------------------------------------------------------------------
//! \fn Real A3(const Real x1,const Real x2,const Real x3)
//  \brief A3: 3-component of vector potential

KOKKOS_INLINE_FUNCTION
Real A3(const Real x1, const Real x2) {
  return (B0/(4.0*M_PI))*(std::cos(4.0*M_PI*x1) - 2.0*std::cos(2.0*M_PI*x2));
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Orszag-Tang test.  The initial conditions are
//  constructed assuming the domain extends over [-0.5x0.5, -0.5x0.5], so that exact
//  symmetry can be enforced across x=0 and y=0.

void ProblemGenerator::OrszagTang_(MeshBlockPack *pmbp, ParameterInput *pin)
{
  if (pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Orszag-Tang test can only be run in MHD, but no <mhd> block "
              << "in input file" << std::endl;
    exit(EXIT_FAILURE);
  }

  B0 = 1.0/std::sqrt(4.0*M_PI);  // defined as global variable shared with A3() function
  Real d0 = 25.0/(36.0*M_PI);
  Real v0 = 1.0;
  Real p0 = 5.0/(12.0*M_PI);

  // capture variables for kernel
  int &nx1 = pmbp->mb_cells.nx1;
  int &nx2 = pmbp->mb_cells.nx2;
  int &is = pmbp->mb_cells.is, &ie = pmbp->mb_cells.ie;
  int &js = pmbp->mb_cells.js, &je = pmbp->mb_cells.je;
  int &ks = pmbp->mb_cells.ks, &ke = pmbp->mb_cells.ke;

  EOS_Data &eos = pmbp->pmhd->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &size = pmbp->pmb->mbsize;

  par_for("pgen_ot1", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i)
    {
      Real x1v = CellCenterX(i-is, nx1, size.x1min.d_view(m), size.x1max.d_view(m));
      Real x2v = CellCenterX(j-js, nx2, size.x2min.d_view(m), size.x2max.d_view(m));

      // compute cell-centered conserved variables
      u0(m,IDN,k,j,i) = d0;
      u0(m,IM1,k,j,i) =  d0*v0*std::sin(2.0*M_PI*x2v);
      u0(m,IM2,k,j,i) = -d0*v0*std::sin(2.0*M_PI*x1v);
      u0(m,IM3,k,j,i) = 0.0;

      // Compute face-centered fields from curl(A).
      Real x1f   = LeftEdgeX(i  -is, nx1, size.x1min.d_view(m), size.x1max.d_view(m));
      Real x1fp1 = LeftEdgeX(i+1-is, nx1, size.x1min.d_view(m), size.x1max.d_view(m));
      Real x2f   = LeftEdgeX(j  -js, nx2, size.x2min.d_view(m), size.x2max.d_view(m));
      Real x2fp1 = LeftEdgeX(j+1-js, nx2, size.x2min.d_view(m), size.x2max.d_view(m));
      Real dx1 = size.dx1.d_view(m);
      Real dx2 = size.dx2.d_view(m);

      b0.x1f(m,k,j,i) =  (A3(x1f,  x2fp1) - A3(x1f,x2f))/dx2;
      b0.x2f(m,k,j,i) = -(A3(x1fp1,x2f  ) - A3(x1f,x2f))/dx1;
      b0.x3f(m,k,j,i) = 0.0;

      // Include extra face-component at edge of block in each direction
      if (i==ie) {
        b0.x1f(m,k,j,i+1) =  (A3(x1fp1,x2fp1) - A3(x1fp1,x2f))/dx2;
      }
      if (j==je) {
        b0.x2f(m,k,j+1,i) = -(A3(x1fp1,x2fp1) - A3(x1f,x2fp1))/dx1;
      }
      if (k==ke) {
        b0.x3f(m,k+1,j,i) = 0.0;
      }
    }
  );

  // initialize total energy (requires B to be defined across entire grid first)
  par_for("pgen_ot2", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i)
    {
      u0(m,IEN,k,j,i) = p0/gm1 + (0.5/u0(m,IDN,k,j,i))*
           (SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) + SQR(u0(m,IM3,k,j,i))) +
            0.5*(SQR(0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1))) +
                 SQR(0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i))) +
                 SQR(0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i))));
    }
  );

  return;
}