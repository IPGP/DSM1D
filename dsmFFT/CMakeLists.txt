enable_language(Fortran)

find_package(HDF5 REQUIRED Fortran)
include_directories(${HDF5_INCLUDE_DIRS})

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")

file(GLOB SRCS "*.f90")

add_executable(dsmFFT.x ${SRCS})
target_link_libraries(dsmFFT.x ${HDF5_LIBRARIES})
