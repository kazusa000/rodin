set(LEVELSET_STOKES_TARGETS
  PETSc_LevelSetStokes3DObstacle
  PETSc_LevelSetStokes3DObstacle_NS
  PETSc_LevelSetStokes3DObstacle_RV
  PETSc_LevelSetStokes2DObstacle
  PETSc_LevelSetStokes2DObstacle_NS
  PETSc_LevelSetStokes2DObstacle_RV
  PETSc_LevelSetStokes2DObstacle_NSRV
)

add_executable(PETSc_LevelSetStokes3DObstacle
  LevelSetStokes3DObstacle.cpp
)

add_executable(PETSc_LevelSetStokes3DObstacle_NS
  LevelSetStokes3DObstacle_NS.cpp
)

add_executable(PETSc_LevelSetStokes3DObstacle_RV
  LevelSetStokes3DObstacle_RV.cpp
)

add_executable(PETSc_LevelSetStokes2DObstacle
  LevelSetStokes2DObstacle.cpp
)

add_executable(PETSc_LevelSetStokes2DObstacle_NS
  LevelSetStokes2DObstacle_NS.cpp
)

add_executable(PETSc_LevelSetStokes2DObstacle_RV
  LevelSetStokes2DObstacle_RV.cpp
)

add_executable(PETSc_LevelSetStokes2DObstacle_NSRV
  LevelSetStokes2DObstacle_NSRV.cpp
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

target_compile_definitions(PETSc_LevelSetStokes3DObstacle
  PRIVATE
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)

target_compile_definitions(PETSc_LevelSetStokes3DObstacle_NS
  PRIVATE
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)

target_compile_definitions(PETSc_LevelSetStokes3DObstacle_RV
  PRIVATE
  LEVELSET_STOKES_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/3dshapes/sphere_init_aligned.mesh"
)

target_compile_definitions(PETSc_LevelSetStokes2DObstacle
  PRIVATE
  LEVELSET_STOKES_2D_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/2dshapes/circle_init_2d.mesh"
)

target_compile_definitions(PETSc_LevelSetStokes2DObstacle_NS
  PRIVATE
  LEVELSET_STOKES_2D_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/2dshapes/circle_init_2d.mesh"
)

target_compile_definitions(PETSc_LevelSetStokes2DObstacle_RV
  PRIVATE
  LEVELSET_STOKES_2D_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/2dshapes/circle_init_2d.mesh"
)

target_compile_definitions(PETSc_LevelSetStokes2DObstacle_NSRV
  PRIVATE
  LEVELSET_STOKES_2D_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/2dshapes/circle_init_2d.mesh"
)
