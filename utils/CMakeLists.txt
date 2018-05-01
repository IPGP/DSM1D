enable_language(Fortran)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")

#file(GLOB SRCS "*.f90")
set(SRCS merge.f90)

add_executable(merge.x ${SRCS})
target_link_libraries(merge.x ${HDF5_LIBRARIES})
