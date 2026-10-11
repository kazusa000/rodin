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

# Opt-in interfaces and focused numerical contracts. The default KelvinBall
# executable remains the Carlos workflow with zero outer velocity.
foreach(check OuterVelocityTest PressureGaugeTest SeedTest SoftThicknessDiagnostic)
  add_executable(KelvinBall${check} EXCLUDE_FROM_ALL ${check}.cpp)
  target_link_libraries(KelvinBall${check} PRIVATE KelvinBallSupport)
endforeach()
