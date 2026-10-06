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
    SymmetryTest.cpp Common.cpp SewedOutput.cpp)
  target_link_libraries(KelvinBall${mode}GeometryTest PRIVATE ${kelvin_libraries})
endforeach()
target_compile_definitions(KelvinBallTetraGeometryTest PRIVATE KELVIN_TETRAHEDRAL=1)

# Opt-in cubic chiral sphere start. Keep KelvinBall and its spherical seed intact.
add_executable(KelvinBallCubicRho EXCLUDE_FROM_ALL ${kelvin_sources})
target_compile_definitions(KelvinBallCubicRho PRIVATE KELVIN_CUBIC_CHIRAL_SEED=1)
target_link_libraries(KelvinBallCubicRho PRIVATE ${kelvin_libraries})
add_executable(KelvinBallCubicRhoGeometryTest EXCLUDE_FROM_ALL
  SymmetryTest.cpp Common.cpp SewedOutput.cpp)
target_compile_definitions(KelvinBallCubicRhoGeometryTest PRIVATE KELVIN_CUBIC_CHIRAL_SEED=1)
target_link_libraries(KelvinBallCubicRhoGeometryTest PRIVATE ${kelvin_libraries})

foreach(mode IN ITEMS Cubic Tetra)
  add_executable(KelvinBall${mode}Thickness EXCLUDE_FROM_ALL
    ThicknessDiagnostic.cpp Common.cpp SewedOutput.cpp)
  target_link_libraries(KelvinBall${mode}Thickness PRIVATE ${kelvin_libraries})
  add_executable(KelvinBall${mode}Boundary EXCLUDE_FROM_ALL
    KelvinBallBoundary.cpp Common.cpp Configuration.cpp Metrics.cpp
    RotatedNitscheIntegrator.cpp SewedOutput.cpp Sphere.cpp)
  target_link_libraries(KelvinBall${mode}Boundary PRIVATE ${kelvin_libraries})
endforeach()
target_compile_definitions(KelvinBallTetraBoundary PRIVATE KELVIN_TETRAHEDRAL=1)
target_compile_definitions(KelvinBallTetraThickness PRIVATE KELVIN_TETRAHEDRAL=1)
add_executable(KelvinBallCubicRhoBoundary EXCLUDE_FROM_ALL
  KelvinBallBoundary.cpp Common.cpp Configuration.cpp Metrics.cpp
  RotatedNitscheIntegrator.cpp SewedOutput.cpp Sphere.cpp)
target_compile_definitions(KelvinBallCubicRhoBoundary PRIVATE KELVIN_CUBIC_CHIRAL_SEED=1)
target_link_libraries(KelvinBallCubicRhoBoundary PRIVATE ${kelvin_libraries})
add_executable(KelvinBallOuterVelocityTest EXCLUDE_FROM_ALL OuterVelocityTest.cpp)
target_link_libraries(KelvinBallOuterVelocityTest PRIVATE ${kelvin_libraries})
add_executable(KelvinBallPressureGaugeTest EXCLUDE_FROM_ALL PressureGaugeTest.cpp Common.cpp)
target_link_libraries(KelvinBallPressureGaugeTest PRIVATE ${kelvin_libraries})
