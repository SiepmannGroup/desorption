#!/bin/bash
path_to_MCCCS_MN="~/MCCCS-MN/"                   # user specific; should have MCCCS-MN/src/
path_to_compilation_directory="~/exe-MCCCS-MN/"  # user specific
 
if [ ! -d "${path_to_compilation_directory}" ]; then
    mkdir -p ${path_to_compilation_directory}
fi

cd ${path_to_compilation_directory}

rm -rf CMake* cmake_install.cmake Makefile defines.h topmon_fortran_interface.h src contrib

FC=mpixlf95_r CC=mpixlc_r CXX=mpixlcxx_r cmake -DDOUBLE_PRECISION=ON -DUSE_MPI=ON -DUSE_OPENMP=ON -DUSE_OWN=ON -DPREFER_STATIC_LIBS=ON -DJOB_FARMING=ON -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_VERBOSE_MAKEFILE=ON  -DBUILD_SHARED_LIBS=OFF -DCMAKE_Fortran_FLAGS="-g -qarch=qp -qtune=qp -qfree=f90 -WF,-qfpp -O3 -qnoipa -qsmp=omp:noauto -qsimd=auto -qhot=level=2 -qprefetch -qunroll=yes -qalias=nostd" -DCMAKE_CXX_FLAGS="-g -qarch=qp -qtune=qp -O3 -qnoipa -qsmp=omp:noauto -qsimd=auto -qhot=level=2 -qprefetch -qunroll=yes" -DCMAKE_C_FLAGS="-g -qarch=qp -qtune=qp -O3 -qnoipa -qsmp=omp:noauto -qsimd=auto -qhot=level=2 -qprefetch -qunroll=yes" -DCMAKE_EXE_LINKER_FLAGS="-g -qnoipa -Wl,-Bstatic -lc" -DMPI_EXTRA_LIBRARY="" ${path_to_MCCCS_MN}

make -j 8

cd src/
rm topmon

mpixlcxx_r -g -qarch=qp -qtune=qp -O3 -qnoipa -qsmp=omp:noauto -qsimd=auto -qhot=level=2 -qprefetch -qunroll=yes  -qsmp=omp  -g -qarch=qp -qtune=qp -O3 -qnoipa -qsmp=omp:noauto -qsimd=auto -qhot=level=2 -qprefetch -qunroll=yes CMakeFiles/topmon.dir/jfmpi.cpp.o  -o topmon -L/soft/compilers/ibmcmp-may2015/xlf/bg/14.1/bglib64 libtopfor.a -lrt -lpthread -lstdc++ -lxlf90_r -lxlomp_ser -lxlfmath -lc
