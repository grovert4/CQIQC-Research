# #!/bin/bash

# # load modules
# MODULES="CCEnv StdEnv/2023 gcc/12.3 flexiblas/3.3.1 openmpi/4.1.5 cmake/3.27.7 fftw/3.3.10 hdf5-mpi/1.14.2 boost/1.82.0 python/3.10.13 mpi4py/3.1.4 imkl/2023.2.0 llvm/16.0.6 eigen/3.4.0"
# module purge
# module load ${MODULES}
# # choose your compilers here
# BUILDDIR=$(pwd)
# # set installation directory (default pwd/install)
# INSTALLDIR=$(pwd)/install

# export TRIQS_ROOT=${INSTALLDIR}
# export PATH=${INSTALLDIR}/bin:$PATH
# export CPLUS_INCLUDE_PATH=${INSTALLDIR}/include:$CPLUS_INCLUDE_PATH
# export LIBRARY_PATH=${INSTALLDIR}/lib:$LIBRARY_PATH
# export LD_LIBRARY_PATH=${INSTALLDIR}/lib:$LD_LIBRARY_PATH
# export PYTHONPATH=${INSTALLDIR}/lib/python${PYVER}/site-packages:$PYTHONPATH
# export CMAKE_PREFIX_PATH=${INSTALLDIR}/lib/cmake/triqs:${INSTALLDIR}/lib/cmake/cpp2py:$CMAKE_PREFIX_PATH

# # output log files
#!/bin/bash

# choose your compilers here
export CC=clang
export CXX=clang++
export FC=gfortran
#python version detection
export PYVER=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
# compiler flags add stdlib=libc++ for clang
export CXXFLAGS="-stdlib=libc++ -Wno-register -march=native"

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
export PYTHONPATH=${INSTALLDIR}/lib/python${PYVER}/site-packages:$PYTHONPATH
export CMAKE_PREFIX_PATH=${INSTALLDIR}/lib/cmake/triqs:${INSTALLDIR}/lib/cmake/cpp2py:$CMAKE_PREFIX_PATH

# choose triqs + app to be installed
packages="triqs  tprf "


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
