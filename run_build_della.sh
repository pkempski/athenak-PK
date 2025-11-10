module purge
module load nvhpc/24.11 openmpi/cuda-12.6/nvhpc-24.11/4.1.6 cudatoolkit/12.6


athenak=/home/pkempski/athenak-PK
build=/home/pkempski/athenak-PK/build
mkdir $build


cd ${athenak}
cmake \
  -D CMAKE_CXX_COMPILER=$athenak/kokkos/bin/nvcc_wrapper \
  -D Kokkos_ENABLE_CUDA=On \
  -D Kokkos_ARCH_AMPERE80=On \
  -D Athena_ENABLE_MPI=On \
  -D PROBLEM=turb \
  -B $build
cd $build
make -j 4
