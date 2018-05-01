enable_language(Fortran)

find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")

set(SRC
  calmat.f90
  dcsymbdl3.f90
  dcsymbdl.f90
  formpi.f90
  glu2.f90
  myTraPSV.f90
  others.f90
  rk3.f90
  trial.f90
)

add_executable(dsmPSV.x ${SRC})

set(SRC_MPI
  calmat.f90
  dcsymbdl3.f90
  dcsymbdl.f90
  formpi.f90
  glu2.f90
  myTraPSVmpi.f90
  others.f90
  rk3.f90
  trial.f90
)

add_executable(dsmPSVmpi.x ${SRC_MPI})
target_link_libraries(dsmPSVmpi.x ${MPI_Fortran_LIBRARIES})

set(SRC_MPI_MS
  calmat.f90
  dcsymbdl3.f90
  dcsymbdl.f90
  formpi.f90
  glu2.f90
  myTraPSVmpi_ms.f90
  others.f90
  rk3.f90
  trial.f90
)

add_executable(dsmPSVmpi_ms.x ${SRC_MPI_MS})
target_link_libraries(dsmPSVmpi_ms.x ${MPI_Fortran_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(dsmPSVmpi.x PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
  set_target_properties(dsmPSVmpi_ms.x PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(dsmPSVmpi.x PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
  set_target_properties(dsmPSVmpi_ms.x PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
endif()
