set(LEVELSET_STOKES_TARGETS
  PETSc_3Dv1
  PETSc_3Dv2
  PETSc_3Dv3
  PETSc_3Dv4_test
)

add_executable(PETSc_3Dv1
  3Dv1.cpp
)

add_executable(PETSc_3Dv2
  3Dv2.cpp
)

add_executable(PETSc_3Dv3
  3Dv3.cpp
)

add_executable(PETSc_3Dv4_test
  3Dv4_test.cpp
)

foreach(_tgt IN LISTS LEVELSET_STOKES_TARGETS)
  target_link_libraries(${_tgt}
    PUBLIC
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly
    Rodin::Variational
    Rodin::MMG
    Rodin::PETSc
    Rodin::Eikonal
    Rodin::Advection
  )

endforeach()

target_compile_definitions(PETSc_3Dv1
  PRIVATE
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)

target_compile_definitions(PETSc_3Dv2
  PRIVATE
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)

target_compile_definitions(PETSc_3Dv3
  PRIVATE
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)

target_compile_definitions(PETSc_3Dv4_test
  PRIVATE
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)
