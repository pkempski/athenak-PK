#ifndef MESH_MESHBLOCK_HPP_
#define MESH_MESHBLOCK_HPP_
//==================================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//==================================================================================================
//! \file meshblock.hpp
//  \brief defines MeshBlock class, and various structs used in them
//  The Mesh is the overall grid structure, and MeshBlocks are local patches of data
//  (potentially on different levels) that tile the entire domain.

#include <cstdint>     // int64_t

//--------------------------------------------------------------------------------------------------
//! \struct GridIndices
//  \brief  structure to store number and indices of grid cells in a MeshBlock

struct GridIndices {
 public:
  GridIndices() {};
  ~GridIndices() {};

  int is,ie,js,je,ks,ke;   // indices of ACTIVE cells
  int nghost;              // number of ghost zones
  int nx1, nx2, nx3;       // number of ACTIVE cells in each dir
  int ncells1, ncells2, ncells3; // total # of cells each dir (i.e. ncells2=nx2+2*NGHOST if nx2>1)
};

//--------------------------------------------------------------------------------------------------
//! \class MeshBlock
//  \brief data/functions associated with a single block

class MeshBlock {

 public:
  MeshBlock(Mesh *pm, ParameterInput *pin, RegionSize input_size);
  ~MeshBlock();

  // data
  Mesh *pmy_mesh;  // ptr to Mesh containing this MeshBlock
  RegionSize block_size;

  // on 1x coarser level MeshBlock (i.e. ncc2=nx2/2 + 2*NGHOST, if nx2>1)
  GridIndices indx, cindx;

  // functions
  int GetNumberOfMeshBlockCells() { return block_size.nx1 * block_size.nx2 * block_size.nx3; }

 private:
  // data

  // functions

};

#endif // MESH_MESHBLOCK_HPP_