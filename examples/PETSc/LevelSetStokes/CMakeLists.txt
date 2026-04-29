set(LEVELSET_STOKES_TARGETS
  PETSc_3Dv8_test
  PETSc_3Dv9
  PETSc_3Dv11_test
  PETSc_3Dv13
  PETSc_3Dv14_test
  PETSc_3Dv15_test
  PETSc_3DPostSmooth
)

add_executable(PETSc_3Dv8_test
  3Dv8_test.cpp
)

add_executable(PETSc_3Dv9
  3Dv9.cpp
)

add_executable(PETSc_3Dv11_test
  3Dv11_test.cpp
)

add_executable(PETSc_3Dv13
  3Dv13.cpp
)

add_executable(PETSc_3Dv14_test
  3Dv14_test.cpp
)

add_executable(PETSc_3Dv15_test
  3Dv15_test.cpp
)

add_executable(PETSc_3DPostSmooth
  3DPostSmooth.cpp
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

target_compile_definitions(PETSc_3Dv8_test
  PRIVATE
  LEVELSET_STOKES_V8_FEATURE_SMOOTH=1
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)

target_compile_definitions(PETSc_3Dv9
  PRIVATE
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)

target_compile_definitions(PETSc_3Dv11_test
  PRIVATE
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)

target_compile_definitions(PETSc_3Dv13
  PRIVATE
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)

target_compile_definitions(PETSc_3Dv14_test
  PRIVATE
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)

target_compile_definitions(PETSc_3Dv15_test
  PRIVATE
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)

target_compile_definitions(PETSc_3DPostSmooth
  PRIVATE
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)
