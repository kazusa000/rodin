set(kelvin_support_sources
  Common.cpp
  Configuration.cpp
  LaplacianSmoothing.cpp
  Metrics.cpp
  RotatedNitscheIntegrator.cpp
  SewedOutput.cpp
  Sphere.cpp)
set(kelvin_libraries
  Rodin::Adaptation
  Rodin::Assembly
  Rodin::Geometry
  Rodin::IO
  Rodin::MMG
  Rodin::Solver
  Rodin::Variational)

# Specialize the chamber geometry and load family at compile time.
foreach(example IN ITEMS KelvinBall KelvinBallC3 KelvinBallC4)
  add_library(${example}Support STATIC EXCLUDE_FROM_ALL ${kelvin_support_sources})
  target_link_libraries(${example}Support PUBLIC ${kelvin_libraries})
  if(example STREQUAL "KelvinBallC3")
    target_compile_definitions(${example}Support PUBLIC KELVIN_AXIAL_C3)
  elseif(example STREQUAL "KelvinBallC4")
    target_compile_definitions(${example}Support PUBLIC KELVIN_AXIAL_C4)
  endif()
  add_executable(${example} EXCLUDE_FROM_ALL KelvinBall.cpp KelvinBallOptimization.cpp)
  target_link_libraries(${example} PRIVATE ${example}Support Rodin::Advection)
  add_executable(${example}Sphere EXCLUDE_FROM_ALL KelvinBallSphere.cpp)
  target_link_libraries(${example}Sphere PRIVATE ${example}Support)
endforeach()
