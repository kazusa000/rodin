window.BENCHMARK_DATA = {
  "lastUpdate": 1772641126684,
  "repoUrl": "https://github.com/kazusa000/rodin",
  "entries": {
    "C++ Rodin Benchmarks": [
      {
        "commit": {
          "author": {
            "email": "carlos.brito524@gmail.com",
            "name": "Carlos Brito-Pacheco",
            "username": "cbritopacheco"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d8c9015de26d8652bc729d670b93b6bed292575d",
          "message": "Develop (#177)\n\n* Add OpenMP specialization for Assembly\n\n* Remove unsafe PETSc assembly casts and harden mesh/FES dispatch (#167)\n\n* Initial plan\n\n* Implement multi-variable assembly specialization\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Refine triplet filtering for Dirichlet DOFs\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Apply Dirichlet elimination to single-variable problems\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix dense assembly fallback without setFromTriplets\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Direct dense assembly paths\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add PETSc single-variable sequential assembly\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix PETSc sequential iteration type deduction\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Guard PETSc sequential for local meshes\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Enforce local-only PETSc sequential and drop Generic fallback\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Remove Generic\n\n* Fix compilation\n\n* Elimination logic\n\n* Fix bug\n\n* Co-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Initialize PETSc solution vectors during assembly\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Replace PETSc assembly void* mesh/FES access with typed visitors\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Remove void* mesh/FES access in PETSc MPI assembly\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Hoist block offsets in PETSc MPI global assembly\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* PETSc multi-field sequential assembly\n\n* Fix P1 MPI IO\n\n* Update\n\n* CI\n\n* Make Stokes3D tests pass\n\n* Update Poisson tests\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\nCo-authored-by: cbritopacheco <carlos.brito524@gmail.com>\n\n* Add getOrder for Min/Max and remove temporary getOrder ops test (#168)\n\n* Initial plan\n\n* Add P1 exact residual manufactured tests\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Clarify mixed residual test comments and naming\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix tests\n\n* Add H1 exact residual manufactured tests\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix Stokes test\n\n* Add getOrder coverage tests for variational ops\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix tests\n\n* Fix tests\n\n* Add getOrder handling and coverage for variational operators\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Remove GetOrderOpsTest per request\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add getOrder for Min/Max and remove temporary getOrder ops test\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\nCo-authored-by: cbritopacheco <carlos.brito524@gmail.com>\n\n* Document and complete P1 quadrature specializations (#169)\n\n* Initial plan\n\n* Add mixed-space support to P1 quadrature specializations\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Refine P1 quadrature documentation and mixed-space handling\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add mixed-space P1 quadrature unit tests\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\nCo-authored-by: Carlos Brito <carlos.brito524@gmail.com>\n\n* Try to fix coverage reporting\n\n* Try to make coverage work pt 2\n\n* Fix SpatialMatrix and Point logic in 3D\n\n* Move MMG into Rodin to support it as first class\n\n* Workon the cantilever 3d example\n\n* Align Copilot setup workflow with required job name (#174)\n\n* Initial plan\n\n* Fix copilot setup steps job name\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Harden `Optional<Attribute>` handling in PETSc assembly backends, add robust `Geometry::MarchingTriangles`, and enforce conforming `MarchingTetrahedra` split consistency (#173)\n\n* Fix PETSc OpenMP/MPI LinearSystem solution-vector sizing, initialization, and RHS sign consistency (#176)\n\n* Initial plan\n\n* fix: correct PETSc solution vector setup in OpenMP and MPI assembly\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* fix: align PETSc OpenMP/MPI block RHS signs with sequential assembly\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Update\n\n---------\n\nCo-authored-by: Copilot <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>",
          "timestamp": "2026-03-03T17:40:43+01:00",
          "tree_id": "c7f723e4b5e7af6daf164ded02446ceb78deb913",
          "url": "https://github.com/kazusa000/rodin/commit/d8c9015de26d8652bc729d670b93b6bed292575d"
        },
        "date": 1772641125075,
        "tool": "googlecpp",
        "benches": [
          {
            "name": "P1Benchmark/UniformTriangular16_Build",
            "value": 0.3147433328842604,
            "unit": "ns/iter",
            "extra": "iterations: 2242038246\ncpu: 0.3147132798732819 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_Build",
            "value": 0.3113649546770902,
            "unit": "ns/iter",
            "extra": "iterations: 2246539421\ncpu: 0.3113395315754845 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular64_Build",
            "value": 0.3112746169902145,
            "unit": "ns/iter",
            "extra": "iterations: 2245677864\ncpu: 0.3112576243482088 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular128_Build",
            "value": 0.31246676725969647,
            "unit": "ns/iter",
            "extra": "iterations: 2250467063\ncpu: 0.3124323908401053 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/2D_Square_GridFunction_Projection_Real_SumOfComponents",
            "value": 648.7442779888377,
            "unit": "ns/iter",
            "extra": "iterations: 1057539\ncpu: 648.3253601049224 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular16_GridFunction_Projection_Real_SumOfComponents",
            "value": 122262.2102123146,
            "unit": "ns/iter",
            "extra": "iterations: 5699\ncpu: 122229.0198280401 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_GridFunction_Projection_Real_SumOfComponents",
            "value": 524099.04915517266,
            "unit": "ns/iter",
            "extra": "iterations: 1302\ncpu: 524055.66743471543 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/2D_Square_GridFunction_Projection_Vector_Components",
            "value": 998.4453749502592,
            "unit": "ns/iter",
            "extra": "iterations: 706187\ncpu: 998.3061625320196 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular16_GridFunction_Projection_Vector_Components",
            "value": 211528.9957377075,
            "unit": "ns/iter",
            "extra": "iterations: 3050\ncpu: 211503.0419672133 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_GridFunction_Projection_Vector_Components",
            "value": 909797.7830065455,
            "unit": "ns/iter",
            "extra": "iterations: 765\ncpu: 909713.4692810453 ns\nthreads: 1"
          },
          {
            "name": "Poisson_UniformGrid_16x16/Assembly_NoCoefficient_ConstantSource",
            "value": 323371.46453089855,
            "unit": "ns/iter",
            "extra": "iterations: 2185\ncpu: 323321.9688787186 ns\nthreads: 1"
          },
          {
            "name": "Poisson_UniformGrid_16x16/Assembly_ConstantCoefficient_ConstantSource",
            "value": 255378.657604709,
            "unit": "ns/iter",
            "extra": "iterations: 2722\ncpu: 255345.00257163856 ns\nthreads: 1"
          },
          {
            "name": "MeshIO/Load_MEDIT_2D_Square",
            "value": 16545.128309765747,
            "unit": "ns/iter",
            "extra": "iterations: 42148\ncpu: 16543.641572553857 ns\nthreads: 1"
          },
          {
            "name": "MeshIO/Load_MEDIT_2D_UniformTriangular64",
            "value": 6469990.694736944,
            "unit": "ns/iter",
            "extra": "iterations: 95\ncpu: 6469267.736842112 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_16x16",
            "value": 121209.83661826377,
            "unit": "ns/iter",
            "extra": "iterations: 5784\ncpu: 121188.36168741345 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_64x64",
            "value": 2236105.788961103,
            "unit": "ns/iter",
            "extra": "iterations: 308\ncpu: 2235907.3246753244 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_128x128",
            "value": 9814653.90277789,
            "unit": "ns/iter",
            "extra": "iterations: 72\ncpu: 9813382.069444425 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_256x256",
            "value": 49732270.24999952,
            "unit": "ns/iter",
            "extra": "iterations: 12\ncpu: 49728645.50000011 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_512x512",
            "value": 325659895.4999958,
            "unit": "ns/iter",
            "extra": "iterations: 2\ncpu: 325633630.5000005 ns\nthreads: 1"
          },
          {
            "name": "Connectivity/Triangular_16x16",
            "value": 48.95318172355231,
            "unit": "ns/iter",
            "extra": "iterations: 14391068\ncpu: 48.94533164599031 ns\nthreads: 1"
          },
          {
            "name": "Connectivity/Triangular_32x32",
            "value": 48.831119159369706,
            "unit": "ns/iter",
            "extra": "iterations: 13915498\ncpu: 48.82791654312332 ns\nthreads: 1"
          },
          {
            "name": "Connectivity/Triangular_64x64",
            "value": 49.231271429578875,
            "unit": "ns/iter",
            "extra": "iterations: 12484789\ncpu: 49.22615360179488 ns\nthreads: 1"
          }
        ]
      }
    ]
  }
}