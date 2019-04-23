list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# the name of the target operating system
SET(CMAKE_SYSTEM_NAME BlueGeneQ-static)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Shared libraries disabled for BlueGene/Q." FORCE)
set(PREFER_STATIC_LIBS ON CACHE BOOL "Prefer static libraries for BlueGene/Q." FORCE)

#__BlueGeneQ_set_static_flags(XL Fortran)
#__BlueGeneQ_set_static_flags(XL C)
#__BlueGeneQ_set_static_flags(XL CXX)

# xl.ndebug is appropriate for production calculations. For debugging,
# use xl to add back error checks and assertions. Using the
# thread-safe compiler version is required, so use (e.g.)
# CMAKE_C_COMPILER=/bgsys/drivers/ppcfloor/comm/xl.ndebug/bin/mpixlc_r
# CMAKE_CXX_COMPILER=/bgsys/drivers/ppcfloor/comm/xl.ndebug/bin/mpixlcxx_r

mark_as_advanced(CMAKE_XL_CreateExportList) # No idea what spams this

# set the compiler
set(CMAKE_Fortran_COMPILER bgxlf90_r)
set(MPI_Fortran_COMPILER mpixlf90_r)
set(CMAKE_C_COMPILER bgxlc_r)
set(MPI_C_COMPILER mpixlc_r)
set(CMAKE_CXX_COMPILER bgxlc++_r)
set(MPI_CXX_COMPILER mpixlcxx_r)

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
# set(CMAKE_FIND_ROOT_PATH
#   /bgsys/drivers/V1R2M2/ppc64/gnu-linux/powerpc64-bgq-linux/
# )

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
# set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
# set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
