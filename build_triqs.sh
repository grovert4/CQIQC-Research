# #!/bin/bash

MODULES=" CCEnv StdEnv/2020 gcc/11.3.0 flexiblas openmpi cmake fftw/3.3.10 hdf5 python/3.10.2 llvm/16 eigen clang"
module purge
module load ${MODULES}
source ~/triqsenv/bin/activate

export CC=gcc
export CXX=g++
export PYVER=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
export CXXFLAGS="-stdlib=libc++ -Wno-register -march=native"
export FC=gfortran
# compiler flags add stdlib=libc++ for clang

# set blas / lapack Intel10_64_dyn | OpenBLAS | FlexiBLAS
export BLA_VENDOR=Intel10_64_dyn

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
export LD_LIBRARY_PATH=/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/gcccore/11.3.0/lib64:$LD_LIBRARY_PATH
export PYTHONPATH=${INSTALLDIR}/lib/python${PYVER}/site-packages:$PYTHONPATH
export CMAKE_PREFIX_PATH=${INSTALLDIR}/lib/cmake/triqs:${INSTALLDIR}/lib/cmake/cpp2py:$CMAKE_PREFIX_PATH
export Python3_ROOT_DIR = `which python`
packages="triqs cthyb tprf"
printf "hello\nworld\n"
echo $LD_LIBRARY_PATH
printf "hello\nworld\n"
echo $LIBRARY_PATH
for pkg in ${packages} ; do 
    cd ${BUILDDIR}
    git clone -b unstable --depth 1 https://github.com/TRIQS/$pkg $pkg.src
    # fetch latest changes
    cd $pkg.src && git pull
    mkdir -p build && cd build
    cmake ../ -DCMAKE_INSTALL_PREFIX=${INSTALLDIR} -DMPIEXEC_PREFLAGS='--allow-run-as-root' 
    make -j$NCORES
    # some test may use mpi
    ctest --output-on-failure
    make install
done
