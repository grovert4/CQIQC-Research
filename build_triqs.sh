# #!/bin/bash

MODULES=" CCEnv StdEnv2023 gcc flexiblas openmpi cmake fftw  hdf5 boost python/3.10.13 llvm/16 eigen mpfr clang"
module purge
module load ${MODULES}

export CC=clang
export CXX=clang++
export CFLAGS="-march=broadwell"
export CXXFLAGS="-stdlib=libc++ -Wno-register -march=broadwell"
export FC=gfortran

export BLA_VENDOR=FlexiBLAS

# set up MKL / OpenMP:
export MKL_INTERFACE_LAYER=GNU,LP64
export MKL_THREADING_LAYER=SEQUENTIAL
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

# set number of threads for compiling and testing
NCORES=10

BUILDDIR=$(pwd)
# set installation directory (default pwd/install)
INSTALLDIR=$(pwd)/install

export TRIQS_ROOT=${INSTALLDIR}
export PATH=${INSTALLDIR}/bin:$PATH
export CPLUS_INCLUDE_PATH=${INSTALLDIR}/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=${INSTALLDIR}/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=${INSTALLDIR}/lib:$LD_LIBRARY_PATH
export PYTHONPATH=${INSTALLDIR}/lib/python3.10/site-packages:$PYTHONPATH
export CMAKE_PREFIX_PATH=${INSTALLDIR}/lib/cmake/triqs:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${INSTALLDIR}/lib/cmake/cpp2py:$CMAKE_PREFIX_PATH
# choose triqs + app to be installed
packages="triqs tprf"


for pkg in ${packages} ; do 
    cd ${BUILDDIR}
    git clone -b unstable --depth 1 https://github.com/TRIQS/$pkg $pkg.src
    # fetch latest changes
    cd $pkg.src && git pull
    mkdir -p build && cd build
    cmake ../ -DCMAKE_INSTALL_PREFIX=${INSTALLDIR} -DMPIEXEC_PREFLAGS='--allow-run-as-root'
    make -j$NCORES
    # some test may use mpi
    ctest -j1 2>&1 >> ${testlog}
    make install
done
