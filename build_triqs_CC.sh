# #!/bin/bash


MODULES="CCEnv StdEnv/2023 gcc/12 openmpi python/3.11 imkl cmake clang mpi4py hdf5 boost fftw"
module purge
module load ${MODULES}
#source ~/triqsenv/bin/activate
virtualenv triqsvenv
source triqsvenv/bin/activate
pip install mako numpy scipy msgpack packaging matplotlib pandas
export CXXFLAGS="-Wno-register -march=native"

mkdir TRIQS
cd TRIQS
mkdir install
BUILDDIR=$(pwd)

INSTALLDIR=$(pwd)/install
export CXXFLAGS="-Wno-register -march=native"
for pkg in ${packages} ; do 
    cd ${BUILDDIR}
    git clone https://github.com/TRIQS/$pkg $pkg.src
    # fetch latest changes
    cd $pkg.src && git pull
    mkdir -p build && cd build
    cmake ../ -DCMAKE_INSTALL_PREFIX=${INSTALLDIR} 
    make install
    LD_LIBRARY_PATH=$(pwd)/lib:$LD_LIBRARY_PATH
    LD_LIBRARY_PATH=$INSTALLDIR/lib:$LD_LIBRARY_PATH
    # some test may use mpi
    ctest --output-on-failure
done

