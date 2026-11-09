cmake_minimum_required(VERSION 3.20)
project(KelvinL10Reproduction LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 20)
get_filename_component(RODIN_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../../.." ABSOLUTE)
set(RODIN_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(RODIN_BUILD_UNIT_TESTS OFF CACHE BOOL "" FORCE)
set(RODIN_BUILD_MANUFACTURED_TESTS OFF CACHE BOOL "" FORCE)
set(RODIN_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(RODIN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(RODIN_INSTALL_RESOURCES OFF CACHE BOOL "" FORCE)
set(RODIN_WITH_PY OFF CACHE BOOL "" FORCE)
set(RODIN_USE_PETSC OFF CACHE BOOL "" FORCE)
set(RODIN_USE_MPI OFF CACHE BOOL "" FORCE)
set(RODIN_USE_MUMPS ON CACHE BOOL "" FORCE)
set(RODIN_MULTITHREADED ON CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
add_subdirectory("${RODIN_ROOT}" rodin)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/.." kelvin)
# Packaged OpenMPI headers expose C++ compatibility symbols in this serial build.
find_library(KELVIN_MPI_CXX NAMES mpi_cxx
  HINTS "${MUMPS_ROOT}/lib/x86_64-linux-gnu/openmpi/lib" REQUIRED)
target_link_libraries(KelvinBallSupport PUBLIC "${KELVIN_MPI_CXX}")
set_target_properties(KelvinBall PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(KelvinCheckpointTest CheckpointTest.cpp)
target_link_libraries(KelvinCheckpointTest PRIVATE KelvinBallSupport)
set_target_properties(KelvinCheckpointTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

option(KELVIN_BUILD_MEMORY_TESTS "Build targeted MUMPS/SWIFT memory tests" OFF)
if(KELVIN_BUILD_MEMORY_TESTS)
  add_subdirectory("${RODIN_ROOT}/third-party/googletest" googletest)
  add_executable(KelvinMemoryTest
    "${RODIN_ROOT}/tests/unit/Rodin/Solver/MUMPSTest.cpp"
    "${RODIN_ROOT}/tests/unit/Rodin/Adaptation/SWIFTLinearSolverTest.cpp"
    "${RODIN_ROOT}/tests/unit/Rodin/Adaptation/SWIFTExampleParametersTest.cpp")
  target_link_libraries(KelvinMemoryTest PRIVATE KelvinBallSupport GTest::gtest_main)
  set_target_properties(KelvinMemoryTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endif()
