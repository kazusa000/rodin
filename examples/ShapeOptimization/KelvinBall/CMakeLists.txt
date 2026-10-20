add_library(KelvinBallSupport STATIC EXCLUDE_FROM_ALL
  BEM.cpp
  BEMQuotient.cpp
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
foreach(check OuterVelocityTest PressureGaugeTest SeedTest SoftThicknessDiagnostic QuotientProjectionTest)
  add_executable(KelvinBall${check} EXCLUDE_FROM_ALL ${check}.cpp)
  target_link_libraries(KelvinBall${check} PRIVATE KelvinBallSupport)
endforeach()
target_link_libraries(KelvinBallQuotientProjectionTest PRIVATE Rodin::Advection)

add_executable(KelvinBallBEMValidation EXCLUDE_FROM_ALL BEMValidation.cpp)
target_link_libraries(KelvinBallBEMValidation PRIVATE KelvinBallSupport)

add_executable(KelvinBallBEMMMGValidation EXCLUDE_FROM_ALL BEMMMGValidation.cpp)
target_link_libraries(KelvinBallBEMMMGValidation PRIVATE KelvinBallSupport Rodin::Advection)

add_executable(KelvinBallBEMChamberValidation EXCLUDE_FROM_ALL BEMChamberValidation.cpp)
target_link_libraries(KelvinBallBEMChamberValidation PRIVATE KelvinBallSupport)

add_executable(KelvinBallBEMChamberDiagnostics EXCLUDE_FROM_ALL BEMChamberDiagnostics.cpp)
target_link_libraries(KelvinBallBEMChamberDiagnostics PRIVATE KelvinBallSupport)
