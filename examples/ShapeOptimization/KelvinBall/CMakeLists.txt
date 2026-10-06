add_executable(KelvinBall EXCLUDE_FROM_ALL
  KelvinBall.cpp
  KelvinBallOptimization.cpp
  Common.cpp
  Configuration.cpp
  Metrics.cpp
  RotatedNitscheIntegrator.cpp
  SewedOutput.cpp
  Sphere.cpp)
target_link_libraries(KelvinBall
  PRIVATE
  Rodin::Adaptation
  Rodin::Advection
  Rodin::Assembly
  Rodin::Geometry
  Rodin::MMG
  Rodin::Solver
  Rodin::Variational)

add_executable(KelvinBallSphere EXCLUDE_FROM_ALL
  KelvinBallSphere.cpp
  Common.cpp
  Configuration.cpp
  Metrics.cpp
  RotatedNitscheIntegrator.cpp
  SewedOutput.cpp
  Sphere.cpp)
target_link_libraries(KelvinBallSphere
  PRIVATE
  Rodin::Assembly
  Rodin::Geometry
  Rodin::IO
  Rodin::MMG
  Rodin::Solver
  Rodin::Variational)

# Separate executable: the immutable cubic default is not switched at runtime.
get_target_property(kelvin_sources KelvinBall SOURCES)
get_target_property(kelvin_libraries KelvinBall LINK_LIBRARIES)
add_executable(KelvinBallTetra EXCLUDE_FROM_ALL ${kelvin_sources})
target_compile_definitions(KelvinBallTetra PRIVATE KELVIN_TETRAHEDRAL=1)
target_link_libraries(KelvinBallTetra PRIVATE ${kelvin_libraries})

foreach(mode IN ITEMS Cubic Tetra)
  add_executable(KelvinBall${mode}GeometryTest EXCLUDE_FROM_ALL
    SymmetryTest.cpp Common.cpp SewedOutput.cpp Sphere.cpp Configuration.cpp)
  target_link_libraries(KelvinBall${mode}GeometryTest PRIVATE ${kelvin_libraries})
endforeach()
target_compile_definitions(KelvinBallTetraGeometryTest PRIVATE KELVIN_TETRAHEDRAL=1)

foreach(mode IN ITEMS Cubic Tetra)
  add_executable(KelvinBall${mode}ThicknessTest EXCLUDE_FROM_ALL
    SampledThicknessDiagnostic.cpp Common.cpp SewedOutput.cpp)
  target_link_libraries(KelvinBall${mode}ThicknessTest PRIVATE ${kelvin_libraries})
endforeach()
target_compile_definitions(KelvinBallTetraThicknessTest PRIVATE KELVIN_TETRAHEDRAL=1)
