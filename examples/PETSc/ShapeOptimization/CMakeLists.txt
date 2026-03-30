add_executable(PETSc_LevelSetCantilever2D LevelSetCantilever2D.cpp)
target_link_libraries(PETSc_LevelSetCantilever2D
  PUBLIC
  Rodin::Math
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::MMG
  Rodin::PETSc
  )

add_executable(PETSc_LevelSetCantilever3D LevelSetCantilever3D.cpp)
target_link_libraries(PETSc_LevelSetCantilever3D
  PUBLIC
  Rodin::Math
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::MMG
  Rodin::PETSc
  )

add_executable(PETSc_LevelSetStokes3DObstacle LevelSetStokes3DObstacle.cpp)
target_link_libraries(PETSc_LevelSetStokes3DObstacle
  PUBLIC
  Rodin::Math
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::MMG
  Rodin::PETSc
  Rodin::IO
  )

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/sphere_init_aligned.mesh
  ${CMAKE_CURRENT_BINARY_DIR}/sphere_init_aligned.mesh
  COPYONLY
)
