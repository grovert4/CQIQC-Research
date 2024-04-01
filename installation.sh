#!/bin/bash

# installation script for triqs3 stable branch with clang OpenMPI toolchain with new spack modules

# load modules
MODULES="CCEnv StdEnv/2023 gcc/12.3 flexiblas/3.3.1 openmpi/4.1.5 cmake/3.27.7 fftw/3.3.10 nfft hdf5-mpi/1.14.2 boost/1.82.0 python/3.10.13 mpi4py/3.1.4 imkl/2023.2.0 llvm/16.0.6 eigen/3.4.0"
module purge
module load ${MODULES}

export CC=gcc
export CXX=g++
export CFLAGS="-march=broadwell"
export CXXFLAGS="-stdlib=libc++ -Wno-register -march=broadwell"
export FC=gfortran

export BLA_VENDOR=FlexiBLAS

# set up flexiblas:
export MKL_INTERFACE_LAYER=GNU,LP64
export MKL_THREADING_LAYER=SEQUENTIAL
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=12
export NEVANLINNA_NUM_THREADS=4
NCORES=20

BUILDINFO=3.2.x_nix2.2_llvm
BUILDDIR=$(pwd)/triqs${BUILDINFO}_build
INSTALLDIR=$(pwd)/installation
MODULEDIR=$(pwd)/installation/modules
mkdir -p $BUILDDIR
mkdir -p $INSTALLDIR/lib/python3.10/site-packages

export ITENSOR_ROOT=${INSTALLDIR}
export TRIQS_ROOT=${INSTALLDIR}
export PATH=${INSTALLDIR}/bin:$PATH
export CPLUS_INCLUDE_PATH=${INSTALLDIR}/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=${INSTALLDIR}/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=${INSTALLDIR}/lib:$LD_LIBRARY_PATH
export PYTHONPATH=${INSTALLDIR}/lib/python3.10/site-packages:$PYTHONPATH
export CMAKE_PREFIX_PATH=${INSTALLDIR}/lib/cmake/triqs:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${INSTALLDIR}/lib/cmake/cpp2py:$CMAKE_PREFIX_PATH

cd ${BUILDDIR}

module list

# install triqs
cd ${BUILDDIR}
# fetch latest changes
cd triqs.src && git pull && cd ..
rm -rf triqs.build && mkdir -p triqs.build && cd triqs.build

cmake ../triqs.src -DCMAKE_INSTALL_PREFIX=${INSTALLDIR} -DBuild_Deps=Always
# make / test / install
make -j$NCORES
ctest -j$NCORES &>> ${testlog}
make install

cd ${BUILDDIR}
# install TPRF
git clone -b 3.2.x --depth 1 https://github.com/TRIQS/tprf.git tprf.src
# fetch latest changes
cd tprf.src && git pull && cd ..
rm -rf tprf.build && mkdir -p tprf.build && cd tprf.build

cmake ../tprf.src
# make / test / install
make -j$NCORES
ctest -j$NCORES &>> ${testlog}
make install
################



# mkdir -p $MODULEDIR/triqs
# # make the template a proper module
# echo '#%Module' > $MODULEDIR/triqs/$BUILDINFO
# # update module template
# sed "s|REPLACEDIR|${INSTALLDIR}|g;s|MODULES|${MODULES}|g" < src.module >> $MODULEDIR/triqs/$BUILDINFO
