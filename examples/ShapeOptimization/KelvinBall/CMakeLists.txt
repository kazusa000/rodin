add_library(KelvinBallSupport STATIC EXCLUDE_FROM_ALL
  Common.cpp
  Configuration.cpp
  LaplacianSmoothing.cpp
  Metrics.cpp
  RotatedNitscheIntegrator.cpp
  SewedOutput.cpp
  Sphere.cpp)
target_link_libraries(KelvinBallSupport PUBLIC
  Rodin::Adaptation
  Rodin::Assembly
  Rodin::Geometry
  Rodin::IO
  Rodin::MMG
  Rodin::Solver
  Rodin::Variational)

add_executable(KelvinBall EXCLUDE_FROM_ALL
  KelvinBall.cpp
  KelvinBallOptimization.cpp)
target_link_libraries(KelvinBall PRIVATE KelvinBallSupport Rodin::Advection)

add_executable(KelvinBallSphere EXCLUDE_FROM_ALL KelvinBallSphere.cpp)
target_link_libraries(KelvinBallSphere PRIVATE KelvinBallSupport)

option(KELVIN_TETRAHEDRAL "Use the 12 proper tetrahedral rotations" OFF)
if(KELVIN_TETRAHEDRAL)
  target_compile_definitions(KelvinBallSupport PUBLIC KELVIN_TETRAHEDRAL)
  add_executable(KelvinBallTetraTest EXCLUDE_FROM_ALL TetrahedralTest.cpp)
  target_link_libraries(KelvinBallTetraTest PRIVATE KelvinBallSupport)
endif()
