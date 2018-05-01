enable_language(Fortran)

find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")

set(SRC
  calmat.f90
  dclib3.f90
  dclib.f90
  formpi.f90
  glu2.f90
  myTraSH.f90
  others.f90
  trial.f90
)

add_executable(dsmSH.x ${SRC})

set(SRC_MPI
  calmat.f90
  dclib3.f90
  dclib.f90
  formpi.f90
  glu2.f90
  myTraSHmpi.f90
  others.f90
  trial.f90
)

add_executable(dsmSHmpi.x ${SRC_MPI})
target_link_libraries(dsmSHmpi.x ${MPI_Fortran_LIBRARIES})

set(SRC_MPI_MS
  calmat.f90
  dclib3.f90
  dclib.f90
  formpi.f90
  glu2.f90
  myTraSHmpi_ms.f90
  others.f90
  trial.f90
)

add_executable(dsmSHmpi_ms.x ${SRC_MPI_MS})
target_link_libraries(dsmSHmpi_ms.x ${MPI_Fortran_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(dsmSHmpi.x PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
  set_target_properties(dsmSHmpi_ms.x PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(dsmSHmpi.x PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
  set_target_properties(dsmSHmpi_ms.x PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
endif()
