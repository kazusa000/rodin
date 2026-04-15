window.BENCHMARK_DATA = {
  "lastUpdate": 1776244943582,
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
      },
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
          "id": "d92ff0145dc5ca10428187593ad5c80e383638ef",
          "message": "Develop (#224)\n\n* Implement Jump, Average, and Trace for FunctionBase and ShapeFunctionBase (#214)\n\n* Audit MMG module: fix heap overflow and null-deref bugs, add 72 unit tests (#216)\n\n* Fix bugs in MMG5::copyMesh and copySolution, add MMG unit tests\n\n- Fix xtetra copy guard checking src->tetra instead of src->xtetra\n- Fix off-by-one buffer allocation for namein/nameout strings\n- Add comprehensive unit test suite for the MMG module\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/82d679f6-4d96-4483-a8cb-7e288009c465\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix dangling-else warning, unused variables in MMG tests\n\nAll 32 MMG unit tests pass.\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/82d679f6-4d96-4483-a8cb-7e288009c465\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix tests\n\n* Add 40 new MMG Mesh and GridFunction behavior tests (72 total)\n\nNew Mesh tests: scale, attributes, vertex coordinates, connectivity,\nidempotent setters, copy independence, area, isSurface.\n\nNew GridFunction tests: zero init, lambda/constant projection, data\naccess, operator[], all arithmetic ops (scalar and GF-GF), min/max,\nargmin/argmax, copy/move semantics, point evaluation (scalar and\nvector), linear function exact interpolation, I/O roundtrip (MEDIT\nand MFEM), setData, post-adapt projection workflow.\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/bcccce3e-bca1-4fd3-a909-0034e64478a5\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix copyMesh heap overflow, clear metadata on parent move-assign, explicit delete copy-assign\n\n- Fix heap buffer overflow in MMG5::copyMesh namein/nameout (same pattern\n  as the already-fixed copySolution: strlen+1 for null terminator)\n- MMG::Mesh::operator=(Parent&&) now clears corners/ridges/required sets\n  so stale MMG metadata does not persist after receiving a different mesh\n- Explicitly declare copy assignment as deleted (parent deletes it)\n- Update ParentMoveAssignmentClearsMMGMetadata test to verify metadata\n  is cleared\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/d9417891-0ded-4f34-907d-b1c7290c3f8e\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\nCo-authored-by: Carlos Brito <carlos.brito524@gmail.com>\n\n* Fix quadrature\n\n* Auditing variational module for missing expressions and bugs (#217)\n\n* fix: Variational module bugs (EQ base class, Min traceOf/getOrder, header guards, NEQ implementation, getOrder on Re/Im/Frobenius)\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/d0d3fbcf-c665-4ec2-be66-b89a25d734a7\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* fix: ComplexFunction copy/move constructors and Sum CTAD/alias bugs\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/6b467d6f-8e41-427b-90dc-9557da57a3b7\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* fix: ComplexFunction<FReal,FImag>::getValue() uses -> on value types; fix Potential const-correctness\n\nComplexFunction<FReal,FImag> (callable specialization) called m_re->getValue(p)\nbut m_re/m_imag are stored by value and are plain callables, not pointers.\nFixed to m_re(p)/m_imag(p).\n\nPotential ShapeFunctionBase specialization stored m_kernel as\nreference_wrapper<KernelType> (non-const) but constructor takes const&.\nFixed to reference_wrapper<const KernelType>.\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/2829c624-c980-44cc-ad47-e49f51119b54\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* test: add unit tests for EQ, NEQ, LT, GT, LEQ, GEQ, AND, OR, Cosh, Sinh, Transpose, Zero, IdentityMatrix, Frobenius\n\nAdds 64 new unit tests covering previously untested expression operators:\n- ComparisonTest.cpp: 24 tests for EQ, NEQ, LT, GT, LEQ, GEQ\n- LogicalTest.cpp: 10 tests for AND, OR (with bool constants)\n- HyperbolicTest.cpp: 13 tests for Cosh, Sinh (incl. cosh²-sinh²=1 identity)\n- TransposeTest.cpp: 2 tests for matrix Transpose\n- ZeroTest.cpp: 6 tests for scalar/vector Zero\n- IdentityMatrixTest.cpp: 6 tests for IdentityMatrix\n- FrobeniusTest.cpp: 3 tests for Frobenius norm\n\nAll 1951 unit tests pass (1887 original + 64 new).\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/56c24452-ab0f-46ee-bf88-e1cd19612d16\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* test: add unit tests for UnaryMinus, Conjugate, ComplexFunction, Re, Im\n\nAdds 30 new unit tests across 4 new test files:\n- UnaryMinusTest.cpp: 6 tests (negation, zero, double negation, copy, getOrder)\n- ConjugateTest.cpp: 5 tests (real value, imaginary part, pure imaginary, copy, getOrder)\n- ComplexFunctionTest.cpp: 10 tests (integer/real/complex constants, composite, callable, copy, getOrder)\n- ReImTest.cpp: 9 tests (Re/Im extraction, pure real/imaginary, negatives, copy, Re²+Im²=|f|² identity)\n\nAll 1981 unit tests pass (1951 previous + 30 new).\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/c3bbc2bf-4113-454f-b1ce-c01c1b7a425d\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* test: add tests for F::X/Y coordinate functions, BoundaryNormal, RelativeError audit\n\nAdds 13 new unit tests across 3 new test files:\n- CoordinateFunctionTest.cpp: 8 tests (F::X, F::Y, copy, getOrder, composition)\n- BoundaryNormalTest.cpp: 4 tests (construction, copy, getOrder, unit normal)\n- RelativeErrorTest.cpp: 1 test (documents template deduction bug)\n\nAlso documents RelativeError.h bug: template deduction fails because\nGridFunction<FES> doesn't match the partial specialization GridFunction<FES, Data>.\n\nAll 1994 unit tests pass (1981 previous + 13 new).\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/c3bbc2bf-4113-454f-b1ce-c01c1b7a425d\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* style: rename meshDim to dimension in BoundaryNormalTest for consistency\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/c3bbc2bf-4113-454f-b1ce-c01c1b7a425d\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* fix: P0g.h header guard typo, SparseProblem.h wrong header guard, add EQ/NEQ Number overloads\n\n- P0g.h: header guard `RODIN_VARIATIAL_P0G_H` → `RODIN_VARIATIONAL_P0G_H`\n- SparseProblem.h: header guard `RODIN_VARIATIONAL_DENSEPROBLEM_H` → `RODIN_VARIATIONAL_SPARSEPROBLEM_H`\n- EQ.h: add Number overloads (operator==(Number, FunctionBase) and operator==(FunctionBase, Number))\n- NEQ.h: add Number overloads (operator!=(Number, FunctionBase) and operator!=(FunctionBase, Number))\n- ComparisonTest.cpp: update EQ/NEQ NumberLHS/NumberRHS tests to use new Number overloads\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/e8df599a-f8df-45e2-a2a0-f137d9b362bc\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* test: add Problem and MatrixFunction tests (9 new tests)\n\n- ProblemTest.cpp: 4 tests (construction, assemble, solve Poisson, solution size)\n- MatrixFunctionTest.cpp: 5 tests (constant matrix, rows/columns, copy, getOrder, identity)\n\nAll 2003 unit tests pass.\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/e8df599a-f8df-45e2-a2a0-f137d9b362bc\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* fix: P0/Grad.h swapped GradBase template params and inconsistent header guard\n\n- P0/Grad.h: swap GradBase<Derived, Operand> → GradBase<Operand, Derived> to match\n  the base class declaration GradBase<GridFunction<FES, Data>, Derived>\n- P0/Grad.h: rename header guard RODIN_VARIATIONAL_P0_GRADIENT_H → RODIN_VARIATIONAL_P0_GRAD_H\n  for consistency with P1/Grad.h, H1/Grad.h, P0g/Grad.h\n\nAll 2003 unit tests pass.\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/c14fc487-6ded-4630-8d30-084db40cabb2\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* fix: RelativeError.h template deduction bug and nan ambiguity\n\n- Change GridFunction<FES> to GridFunction<FES, Data> in all 4 methods\n  (l1, l2, lInf, compute) so the compiler can deduce template params\n  from the partial specialization GridFunction<FES, Data>\n- Fix Math::Constants::nan() → Math::nan<Real>() (Constants::nan doesn't exist)\n- Add missing includes (Abs.h, Integral.h, RealFunction.h, Assembly.h, Math/Common.h)\n- Use CTAD for local GridFunction variables instead of explicit GridFunction<FES>\n- Remove invalid setWeights() calls (method doesn't exist on GridFunction)\n- Qualify sqrt → Math::sqrt for proper ADL resolution\n- Hoist getFiniteElementSpace() before switch to reduce duplication\n- Replace placeholder test with 5 real functional tests covering all norms\n\nAll 2007 unit tests pass.\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/0a598b3d-cd7f-41f7-b42b-c73fdff2f432\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* feat: implement H1/Derivative.h specialization for H1 GridFunctions with 7 tests\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/2bc05f89-7eca-4e29-b5b3-e2a63ff55528\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* style: remove unnecessary std::move on primitive size_t in H1/Derivative.h\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/2bc05f89-7eca-4e29-b5b3-e2a63ff55528\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* fix: P0g.h missing Grad/Div/Jacobian includes, P0g/ForwardDecls.h wrong header guard; add P0g and DenseProblem tests\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/f6bb7e02-8104-42a0-9240-2f0e530f5797\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* fix: Min.h header guard, SparseProblem.h CTAD/template bug, missing SparseProblem.h include in Variational.h + tests\n\nBug fixes:\n- Min.h: header guard mixed case RODIN_VARIATIONAL_Min_H → RODIN_VARIATIONAL_MIN_H\n- SparseProblem.h: completely broken CTAD guide and class specialization\n  (wrong template params, wrong inheritance structure) — rewritten to\n  match DenseProblem.h pattern with <LinearSystem, U, V>\n- Variational.h: missing SparseProblem.h include (DenseProblem.h was\n  included but SparseProblem.h was not)\n\nNew tests:\n- SparseProblemTest.cpp: 6 tests (construction, assembly, sparse matrix\n  properties, symmetry, solve Poisson, RHS non-zero)\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/fa48af35-6d04-49f1-9364-4794c258b1cf\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add unit tests for Shard and Sharder classes (partitioning/MPI groundwork) with quadrature and reconcile fixes (#215)\n\n* Initial plan\n\n* Add unit tests for Shard class (partitioning/MPI groundwork)\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/9c0755a8-fc16-411e-ac5b-4a865412f620\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add Sharder unit tests using BalancedCompactPartitioner and register both test targets\n\n- Create SharderTest.cpp (12 tests): basic sharding, cell/vertex ownership\n  uniqueness, vertex coverage, ghost-neighbor validation, owner/halo map\n  consistency, polytope map bidirectionality, quad mesh support, multiple\n  partition counts, isLocal consistency, getContext.\n- Update CMakeLists.txt to register RodinGeometryShardTest and\n  RodinGeometrySharderTest targets.\n- All 1911 unit tests pass (including 24 new Shard + Sharder tests).\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/894e0faa-ed0c-4658-95b0-40892923a53d\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Expand Shard and Sharder tests to cover all geometry types and mixed meshes\n\nShardTest.cpp additions (14 new tests):\n- Parent-mode: Quadrilateral, Segment (1D), Tetrahedron (3D), Hexahedron (3D), Wedge (3D)\n- Parent-mode mixed: 2D tri+quad (full and partial inclusion), 3D tet+wedge\n- Direct-mode: Quadrilateral, Segment (1D), Tetrahedron (3D), Hexahedron (3D), Wedge (3D), mixed 2D tri+quad\n\nSharderTest.cpp additions (18 new tests):\n- 1D Segment: basic sharding, isLocal consistency\n- 3D Tetrahedron: basic sharding, cell ownership uniqueness, ghost neighbors\n- 3D Hexahedron: basic sharding, owner map consistency\n- 3D Wedge: basic sharding, halo map consistency\n- Mixed 2D (tri+quad): basic sharding, vertex ownership, polytope map consistency\n- Mixed 3D (tet+wedge): basic sharding, vertex coverage, isLocal consistency\n- Larger mixed 2D: multiple partition counts\n- All 3D uniform types: parametric multi-partition test\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/36a9b547-1f37-489a-a2c6-59a00b8291e1\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add MPI parallel tests for Sharder with mpiexec (np=1,2,3)\n\nNew file: tests/unit/Rodin/MPI/Geometry/SharderTest.cpp (16 tests)\n- Triangle: cells/vertices present, global cell count, vertex ownership,\n  isLocal invariant, polytope map, owner map, dimension consistency\n- Quadrilateral: global cell count\n- Segment (1D): global cell count\n- Tetrahedron (3D): global cell count, dimension consistency\n- Hexahedron (3D): global cell count\n- Wedge (3D): global cell count\n- Mixed 2D (tri+quad): global cell count, vertex ownership, isLocal\n\nCMakeLists runs tests via mpiexec -n {1,2,3} with Open MPI\n--oversubscribe detection for CI portability.\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/bff95e5a-6b32-4ab5-b5cd-459d9b6c95c0\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add audit checklist gap tests: attribute preservation, connectivity, MPI Mesh API\n\nShardTest.cpp (+3 tests):\n- Builder_ParentMode_AttributePreservation: verify attribute round-trip\n- Builder_ParentMode_ConnectivityPreservation: verify shard-local cell→vertex incidence\n- Builder_DirectMode_Mixed3D: mixed tet+wedge in direct mode\n\nMPI SharderTest.cpp (+10 tests for Mesh<Context::MPI>):\n- GetPolytopeCount_Cells/Vertices/ByType: distributed aggregate queries\n- IndexMapping_RoundTrip/MissingReturnsNullopt: getLocalIndex/getGlobalIndex\n- GetContext: MPIMesh context accessor\n- LocalIterators: cell/vertex iterator counts\n- GetGeometry_Cells: geometry type queries\n- GetVertexCoordinates: coordinate validation\n- Tetrahedron3D_GlobalCellCount: 3D distributed cell count\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/4cf0b652-cc5d-49cb-9b8b-8bb0f754dd70\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix PolytopeQuadratureIndex deserialization and add 12 new MPIMesh API tests\n\nFix pre-existing bug: after MPI deserialization, PolytopeQuadratureIndex\nhad 0 dimension buckets because the cache is not serialized. Added\nm_quadratures.initialize(m_sdim) at end of Mesh::serialize() to\nre-initialize the cache on load.\n\nNew tests: Scale, Flush, IsSubMesh, GetConnectivity, GetSetAttribute,\nSetVertexCoordinates, GetArea_Triangle2D, GetVolume_Tetrahedron3D,\nGetPolytopeTransformation, GetCellByIndex, GetVertexByIndex,\nGetMeasure_MatchesGetArea\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/89060670-4255-4624-9bf5-5e02297e002b\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add 5 new MPIMesh tests: asSubMesh, getQuadrature, getBoundary, isBoundary/isInterface, getPerimeter\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/e9796741-4507-49cb-91de-a37ddcac2c8f\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add UniformGrid, SaveLoad, GetInterface tests and fix np=1 interface skip\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/e9796741-4507-49cb-91de-a37ddcac2c8f\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Address code review: use boost::filesystem for temp paths in SaveLoad test\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/e9796741-4507-49cb-91de-a37ddcac2c8f\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add reconcile() tests (9 tests) and fix pre-existing reconcile communication bug\n\nFix: reconcile() had incomplete participant tracking — entities at ghost\nlayer boundaries (2+ hops from owner) never received distributed IDs.\nTwo fixes applied:\n1. Owner sends IDs to ALL neighbors (not just tracked participants)\n2. Added forwarding round so intermediate ranks relay IDs to their neighbors\n\nTests: 9 reconcile tests covering edges (2D tri/quad), faces (3D tet),\nstate exclusivity, polytope map bidirectionality, owner/halo map\nconsistency, no-duplicate ownership, and single-rank correctness.\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/2e0b0e44-1fa0-4811-aa14-802e75cb6036\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix em dash to double hyphen in reconcile comment\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/2e0b0e44-1fa0-4811-aa14-802e75cb6036\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add 3 reconcile tests: 3D state exclusivity, 3D polytope map bidirectional, idempotency\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/20d3371f-7836-4963-8524-fbaaa2d47d68\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Address code review: replace magic numbers with named enum in idempotency test\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/20d3371f-7836-4963-8524-fbaaa2d47d68\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix compile\n\n* CI\n\n* Fix tests\n\n* Optimize reconcile(): remove per-entity candidateNeighbors, use received flags, eliminate allHolders\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/0e94910f-bba3-403a-bf56-573b93b4438c\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Optimize reconcile(): index-based convergence messages (5-6x smaller), per-neighbor structures replace per-entity sets\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/c18ce694-2747-454f-aef0-0144d3133843\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Address code review: fix halo build-then-erase, remove unnecessary dedup, clarify KeyMsg type\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/c18ce694-2747-454f-aef0-0144d3133843\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add Optional<size_t> maxRounds parameter to reconcile() with unbounded default, plus bounded test\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/9679ffd7-dadf-4fa9-8e35-915afadc3f96\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix bounded reconcile test: use explicit Optional<size_t> bound and maxRounds=1 for codimension-1\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/9679ffd7-dadf-4fa9-8e35-915afadc3f96\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Document and add reconcile option presets (#220)\n\n* Fix compile\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\nCo-authored-by: Carlos Brito <carlos.brito524@gmail.com>\n\n* Audit Math module: fix bugs, add missing functions, complete headers, add unit tests (#222)\n\n* Fix bugs and complete Math module headers: dot() assertions, RK2 semantics, umbrella header, CMakeLists, header guards, missing math functions, documentation\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/97c7fc1b-adba-45d4-8a4c-f0f59e3dc16f\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add comprehensive unit tests for Math/Common functions\n\nAdd 18 new TEST_F cases to CommonTest covering:\n- pow2, pow<N>, pow(base,exp), sqrt\n- isNaN, isInf\n- Trig: cos, sin, tan\n- Hyperbolic: cosh, sinh, tanh\n- Inverse trig: acos, asin, atan, atan2\n- Log: log, log2, log10\n- sgn, binom, factorial, permutation\n- nan<T> factory\n- Form language helpers: sum, minus, mult, division\n- dot product for Real, Complex, and Eigen::Vector3d\n\nCo-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add comprehensive SpatialVector unit tests\n\nAdd 19 new TEST_F tests covering SpatialVector API:\n- Initializer list, Eigen constructor\n- Compound operators (+=, -=, *=, /=), unary negation\n- Resize, named accessors (x/y/z), setZero/setConstant\n- Cross product, dot product, transpose\n- value(), normalize/normalized, norms (squared/stable/blue/lp)\n- conjugate, getData, 2D/1D/0D edge cases\n\nFix SpatialVector::conjugate() for real scalar types by guarding\nstd::conj calls with if constexpr to avoid complex-to-real conversion.\n\nCo-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add NewtonRaphson, RungeKutta test files and update test CMakeLists.txt\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/97c7fc1b-adba-45d4-8a4c-f0f59e3dc16f\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Fix build: add SpatialVector.h include to MatrixTest.cpp, fix NewtonRaphson boundary test\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/dfb90f89-4721-440b-95c6-8ce920779040\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add min/max/clamp functions, SpatialMatrix operator-, documentation, LinearSystem tests, additional Common/Matrix tests\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/dfb90f89-4721-440b-95c6-8ce920779040\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add Deg unit type, document SpatialVector storage design, update umbrella header and CMakeLists\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/dfb90f89-4721-440b-95c6-8ce920779040\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Add pow<N> and integral_constant overload documentation, remove duplicate orphan doc block\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/dfb90f89-4721-440b-95c6-8ce920779040\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Geometry module audit: ~182 new tests, 4 bug fixes, Doxygen, umbrella headers (#223)\n\n* Geometry module audit: tests, docs, umbrella headers, displace() flush fix\n\n- Add Polytope::Key unit tests (construction, iteration, hash, equality, serialization)\n- Add Mesh::UniformGrid tests for Quadrilateral, Tetrahedron, Hexahedron\n- Add Mesh geometric measure tests (getArea, getPerimeter, getVolume, getMeasure)\n- Add Mesh operation tests (scale, copy/move, isEmpty, name, setAttribute, setVertexCoordinates, flush)\n- Add SubMesh tests (Builder, getParent, getAncestors, polytopeMap, skin, trim, keep)\n- Add Polytope::Traits extended tests (isTensorProduct, getCentroid values, half-space validation)\n- Fix Mesh::displace() to call flush() after displacing vertices (stale cache bug)\n- Fix Geometry.h umbrella header: add Point, PointCloud, Connectivity, Region, etc.\n- Update Euclidean.h umbrella header documentation (Magnum-dependent headers noted)\n- Add Doxygen to Polytope::Key, SymmetricEquality, SymmetricHash\n- Add Doxygen to Polytope::Traits::isTensorProduct(), getCentroid()\n- Add Doxygen to Mesh::Box(), trace() overloads, ccl() overloads\n- Add Doxygen to Mesh::flush(), getName/setName, getAttributeIndex, etc.\n- Add Doxygen to Mesh::getDefaultPolytopeTransformation()\n- Add @tparam to Mesh::displace()\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/dd344326-80f0-4012-9205-d8d3fb1e1ca6\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Batch 2: Fix IdentityTransformation, Vertex bounds safety, Project::vertex(), 72 new tests\n\n- Fix IdentityTransformation.h: correct virtual signatures (parameter order matched\n  to base class), remove non-existent getJacobianOrder() override, add missing #endif\n- Fix Vertex::y()/z()/operator() bounds safety: add assert(dimension) guards,\n  move implementations to .cpp to avoid incomplete-type errors\n- Implement Polytope::Project::vertex() (was declared but never implemented)\n- Fix Mesh::Box() documentation: clarify it takes face type, not cell type\n- Add Geometry.h umbrella: include IdentityTransformation.h (now fixed)\n- New tests:\n  - PolytopeProjectTest (20): cell/face/boundary/vertex projections for\n    Segment, Triangle, Quadrilateral, Tetrahedron, Hexahedron\n  - PointCloudTest (16): construction, push_back, access, copy/move,\n    resize/clear, matrix views\n  - IdentityTransformationTest (9): construction, transform, jacobian,\n    copy, move for 2D/3D\n  - MeshOpsTest (15): Box tests (face-type API), CCL, trace, isSurface,\n    getPolytope\n  - ConnectivityExtTest (8): multi-step chains, local incidence,\n    clear/recompute, polytope counts\n  - AttributeIndexTest (5): set/get cell/boundary/vertex attributes\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/0e47f34d-2ee9-4edd-b5a6-1497ab86ec5f\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n* Batch 3: Add PolytopeExtTest (Vertex/Face/Cell extended tests), setPolytopeTransformation docs\n\n- Add PolytopeExtTest.cpp: Vertex coordinate access (2D/3D), Face boundary/interface,\n  Cell iteration/adjacency, geometry type checks, isCell/isFace/isVertex, getVertices\n- Add Doxygen for MeshBase::setPolytopeTransformation() (ownership semantics)\n\nAgent-Logs-Url: https://github.com/cbritopacheco/rodin/sessions/0e47f34d-2ee9-4edd-b5a6-1497ab86ec5f\n\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>\n\n---------\n\nCo-authored-by: Copilot <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: cbritopacheco <6352283+cbritopacheco@users.noreply.github.com>",
          "timestamp": "2026-04-08T11:20:01+02:00",
          "tree_id": "030605bb60b109dca2cf1d4e507c43401c11a97c",
          "url": "https://github.com/kazusa000/rodin/commit/d92ff0145dc5ca10428187593ad5c80e383638ef"
        },
        "date": 1776244940793,
        "tool": "googlecpp",
        "benches": [
          {
            "name": "P1Benchmark/UniformTriangular16_Build",
            "value": 0.35170361979269466,
            "unit": "ns/iter",
            "extra": "iterations: 1947372010\ncpu: 0.35165395542477784 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_Build",
            "value": 0.3518662484536186,
            "unit": "ns/iter",
            "extra": "iterations: 1989656786\ncpu: 0.35182129899261927 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular64_Build",
            "value": 0.3517942448160498,
            "unit": "ns/iter",
            "extra": "iterations: 1989866248\ncpu: 0.35174187797973044 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular128_Build",
            "value": 0.3523593383955874,
            "unit": "ns/iter",
            "extra": "iterations: 1988268664\ncpu: 0.3523064285441055 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/2D_Square_GridFunction_Projection_Real_SumOfComponents",
            "value": 591.9388033599478,
            "unit": "ns/iter",
            "extra": "iterations: 1186405\ncpu: 591.4999245620174 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular16_GridFunction_Projection_Real_SumOfComponents",
            "value": 120850.8873579466,
            "unit": "ns/iter",
            "extra": "iterations: 5806\ncpu: 120838.19600413373 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_GridFunction_Projection_Real_SumOfComponents",
            "value": 521392.17602964886,
            "unit": "ns/iter",
            "extra": "iterations: 1335\ncpu: 521325.0224719101 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/2D_Square_GridFunction_Projection_Vector_Components",
            "value": 1003.0902134519777,
            "unit": "ns/iter",
            "extra": "iterations: 713275\ncpu: 1002.9180554484608 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular16_GridFunction_Projection_Vector_Components",
            "value": 206966.3379697552,
            "unit": "ns/iter",
            "extra": "iterations: 3379\ncpu: 206946.82509618276 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_GridFunction_Projection_Vector_Components",
            "value": 870220.9575000097,
            "unit": "ns/iter",
            "extra": "iterations: 800\ncpu: 870145.3150000017 ns\nthreads: 1"
          },
          {
            "name": "Poisson_UniformGrid_16x16/Assembly_NoCoefficient_ConstantSource",
            "value": 178064.94371526243,
            "unit": "ns/iter",
            "extra": "iterations: 3962\ncpu: 178039.92200908632 ns\nthreads: 1"
          },
          {
            "name": "Poisson_UniformGrid_16x16/Assembly_ConstantCoefficient_ConstantSource",
            "value": 174844.79015416303,
            "unit": "ns/iter",
            "extra": "iterations: 4022\ncpu: 174816.47339632016 ns\nthreads: 1"
          },
          {
            "name": "MeshIO/Load_MEDIT_2D_Square",
            "value": 15811.929881482216,
            "unit": "ns/iter",
            "extra": "iterations: 44382\ncpu: 15810.25474291381 ns\nthreads: 1"
          },
          {
            "name": "MeshIO/Load_MEDIT_2D_UniformTriangular64",
            "value": 7922186.652175025,
            "unit": "ns/iter",
            "extra": "iterations: 92\ncpu: 7920785.206521739 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_16x16",
            "value": 148975.07761974534,
            "unit": "ns/iter",
            "extra": "iterations: 4638\ncpu: 148952.33570504538 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_64x64",
            "value": 2686202.9807699188,
            "unit": "ns/iter",
            "extra": "iterations: 260\ncpu: 2685809.9500000062 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_128x128",
            "value": 15573262.711111763,
            "unit": "ns/iter",
            "extra": "iterations: 45\ncpu: 15571096.00000004 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_256x256",
            "value": 120974953.00001054,
            "unit": "ns/iter",
            "extra": "iterations: 6\ncpu: 120959151.1666665 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_512x512",
            "value": 544873009.0001845,
            "unit": "ns/iter",
            "extra": "iterations: 1\ncpu: 544720193.9999999 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_AllPairs",
            "value": 1299323.9794430288,
            "unit": "ns/iter",
            "extra": "iterations: 535\ncpu: 1299181.876635508 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_AllPairs",
            "value": 17.79249229632145,
            "unit": "ns/iter",
            "extra": "iterations: 38129182\ncpu: 17.785988826091234 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_0_0",
            "value": 1295871.6555675454,
            "unit": "ns/iter",
            "extra": "iterations: 540\ncpu: 1295694.942592606 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_0_1",
            "value": 1014215.9101305049,
            "unit": "ns/iter",
            "extra": "iterations: 701\ncpu: 1013719.2924393622 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_1_0",
            "value": 475210.2108067742,
            "unit": "ns/iter",
            "extra": "iterations: 1480\ncpu: 474968.07702705223 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_1_1",
            "value": 994129.5960236302,
            "unit": "ns/iter",
            "extra": "iterations: 703\ncpu: 994115.716927499 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_0_0",
            "value": 4.287935800249057,
            "unit": "ns/iter",
            "extra": "iterations: 165667930\ncpu: 4.287410279104712 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_0_1",
            "value": 4.229576572715754,
            "unit": "ns/iter",
            "extra": "iterations: 165746263\ncpu: 4.229051607637173 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_1_0",
            "value": 2.4616032736165487,
            "unit": "ns/iter",
            "extra": "iterations: 283725607\ncpu: 2.461210235422985 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_1_1",
            "value": 4.221081233755391,
            "unit": "ns/iter",
            "extra": "iterations: 165694050\ncpu: 4.220690157552417 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_AllPairs",
            "value": 10061606.559999442,
            "unit": "ns/iter",
            "extra": "iterations: 75\ncpu: 10060698.026667107 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_AllPairs",
            "value": 42.12155273699114,
            "unit": "ns/iter",
            "extra": "iterations: 16615718\ncpu: 42.1164609317515 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_2_1",
            "value": 6783033.631168631,
            "unit": "ns/iter",
            "extra": "iterations: 122\ncpu: 6781484.55737687 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_1_2",
            "value": 6838473.255322938,
            "unit": "ns/iter",
            "extra": "iterations: 94\ncpu: 6838151.265957763 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_2_2",
            "value": 2432990.124102197,
            "unit": "ns/iter",
            "extra": "iterations: 282\ncpu: 2432598.599290722 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_1_1",
            "value": 8494171.569641884,
            "unit": "ns/iter",
            "extra": "iterations: 79\ncpu: 8494000.151899166 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_2_1",
            "value": 4.925172713527388,
            "unit": "ns/iter",
            "extra": "iterations: 141968612\ncpu: 4.924454075806547 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_1_2",
            "value": 4.92675030454788,
            "unit": "ns/iter",
            "extra": "iterations: 142153847\ncpu: 4.9261347883184845 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_2_2",
            "value": 4.221516205695852,
            "unit": "ns/iter",
            "extra": "iterations: 165745097\ncpu: 4.220411636067885 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_1_1",
            "value": 5.633260529761141,
            "unit": "ns/iter",
            "extra": "iterations: 124142719\ncpu: 5.631675378400561 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Build_1",
            "value": 4592984.32690211,
            "unit": "ns/iter",
            "extra": "iterations: 156\ncpu: 4591693.467948342 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Build_1",
            "value": 1377960.9196838038,
            "unit": "ns/iter",
            "extra": "iterations: 498\ncpu: 1378106.5401605517 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Transpose_1_2",
            "value": 3083722.6838010694,
            "unit": "ns/iter",
            "extra": "iterations: 253\ncpu: 3082647.2015810595 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Transpose_1_2",
            "value": 567139.8624038269,
            "unit": "ns/iter",
            "extra": "iterations: 1221\ncpu: 567377.0614251258 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Intersection_2_2_via_0",
            "value": 2266778.7231218903,
            "unit": "ns/iter",
            "extra": "iterations: 307\ncpu: 2266061.944625471 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Intersection_2_2_via_0",
            "value": 1106510.708473893,
            "unit": "ns/iter",
            "extra": "iterations: 638\ncpu: 1106559.6755487686 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Intersection_1_1_via_0",
            "value": 4035155.230772384,
            "unit": "ns/iter",
            "extra": "iterations: 182\ncpu: 4034496.13186766 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Intersection_1_1_via_0",
            "value": 1580312.6388231695,
            "unit": "ns/iter",
            "extra": "iterations: 443\ncpu: 1580400.1376974266 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_AllPairs",
            "value": 6034656.983320019,
            "unit": "ns/iter",
            "extra": "iterations: 120\ncpu: 6034525.39999978 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_AllPairs",
            "value": 42.08529091548648,
            "unit": "ns/iter",
            "extra": "iterations: 16661540\ncpu: 42.08117490940212 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_2_1",
            "value": 3558667.658179432,
            "unit": "ns/iter",
            "extra": "iterations: 196\ncpu: 3558795.8316329974 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_1_2",
            "value": 3893859.32044583,
            "unit": "ns/iter",
            "extra": "iterations: 181\ncpu: 3893899.2762430958 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_2_2",
            "value": 1345433.3627617585,
            "unit": "ns/iter",
            "extra": "iterations: 521\ncpu: 1345205.9731285346 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_1_1",
            "value": 4884498.965774994,
            "unit": "ns/iter",
            "extra": "iterations: 146\ncpu: 4884934.684932065 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_2_1",
            "value": 4.949338786419986,
            "unit": "ns/iter",
            "extra": "iterations: 142177585\ncpu: 4.948978954734657 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_1_2",
            "value": 4.955615694762035,
            "unit": "ns/iter",
            "extra": "iterations: 140993623\ncpu: 4.9547732311268415 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_2_2",
            "value": 4.226387408590353,
            "unit": "ns/iter",
            "extra": "iterations: 165192955\ncpu: 4.225653630325797 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_1_1",
            "value": 5.666897498516167,
            "unit": "ns/iter",
            "extra": "iterations: 124114877\ncpu: 5.6660602016308514 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Build_1",
            "value": 2553189.2717326367,
            "unit": "ns/iter",
            "extra": "iterations: 276\ncpu: 2552877.076087558 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Build_1",
            "value": 878164.3692266821,
            "unit": "ns/iter",
            "extra": "iterations: 845\ncpu: 878312.2568046633 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Transpose_1_2",
            "value": 1654134.9255804075,
            "unit": "ns/iter",
            "extra": "iterations: 430\ncpu: 1653680.3488371442 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Transpose_1_2",
            "value": 366545.9215199953,
            "unit": "ns/iter",
            "extra": "iterations: 1924\ncpu: 366763.62318069476 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Intersection_2_2_via_0",
            "value": 1192208.2902659492,
            "unit": "ns/iter",
            "extra": "iterations: 596\ncpu: 1191909.0838922316 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Intersection_2_2_via_0",
            "value": 479852.6880618258,
            "unit": "ns/iter",
            "extra": "iterations: 1449\ncpu: 480059.59006202716 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Intersection_1_1_via_0",
            "value": 2310757.653203867,
            "unit": "ns/iter",
            "extra": "iterations: 297\ncpu: 2310398.952861723 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Intersection_1_1_via_0",
            "value": 854312.0354141258,
            "unit": "ns/iter",
            "extra": "iterations: 819\ncpu: 854484.9474966117 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_AllPairs",
            "value": 86986156.50003679,
            "unit": "ns/iter",
            "extra": "iterations: 6\ncpu: 86963411.99999817 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_AllPairs",
            "value": 78.53066561804641,
            "unit": "ns/iter",
            "extra": "iterations: 8921017\ncpu: 78.51800114269433 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_1",
            "value": 23162074.26667158,
            "unit": "ns/iter",
            "extra": "iterations: 30\ncpu: 23157391.83333439 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_2",
            "value": 27608788.153866746,
            "unit": "ns/iter",
            "extra": "iterations: 26\ncpu: 27603415.230768856 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_2_3",
            "value": 31582614.04548119,
            "unit": "ns/iter",
            "extra": "iterations: 22\ncpu: 31578388.636363152 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_3",
            "value": 12256325.245648202,
            "unit": "ns/iter",
            "extra": "iterations: 57\ncpu: 12255191.578947632 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_1",
            "value": 4.9259364632744225,
            "unit": "ns/iter",
            "extra": "iterations: 142038005\ncpu: 4.9253174951309475 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_2",
            "value": 4.922384063595173,
            "unit": "ns/iter",
            "extra": "iterations: 137806609\ncpu: 4.922065493970654 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_2_3",
            "value": 4.927135024283008,
            "unit": "ns/iter",
            "extra": "iterations: 141793844\ncpu: 4.92653061158283 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_3",
            "value": 4.219935063652577,
            "unit": "ns/iter",
            "extra": "iterations: 165681317\ncpu: 4.219440312633507 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Build_1",
            "value": 12585212.071396004,
            "unit": "ns/iter",
            "extra": "iterations: 56\ncpu: 12584232.428572254 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Build_1",
            "value": 7649386.419372252,
            "unit": "ns/iter",
            "extra": "iterations: 93\ncpu: 7648950.354840136 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Build_2",
            "value": 18337371.674977023,
            "unit": "ns/iter",
            "extra": "iterations: 40\ncpu: 18336638.925001837 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Build_2",
            "value": 8211317.011366439,
            "unit": "ns/iter",
            "extra": "iterations: 88\ncpu: 8210733.011363658 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Transpose_2_3",
            "value": 9259482.974050423,
            "unit": "ns/iter",
            "extra": "iterations: 77\ncpu: 9258105.844156088 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Transpose_2_3",
            "value": 1989050.7690190494,
            "unit": "ns/iter",
            "extra": "iterations: 355\ncpu: 1989115.4281687408 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Intersection_3_3_via_0",
            "value": 11595123.897983953,
            "unit": "ns/iter",
            "extra": "iterations: 49\ncpu: 11593400.551020335 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Intersection_3_3_via_0",
            "value": 8604563.555556165,
            "unit": "ns/iter",
            "extra": "iterations: 81\ncpu: 8603922.9753077 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Intersection_3_3_via_2",
            "value": 9138833.808215356,
            "unit": "ns/iter",
            "extra": "iterations: 73\ncpu: 9137109.753423993 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Intersection_3_3_via_2",
            "value": 2027469.494202896,
            "unit": "ns/iter",
            "extra": "iterations: 344\ncpu: 2027831.0872090585 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_AllPairs",
            "value": 19093614.32431192,
            "unit": "ns/iter",
            "extra": "iterations: 37\ncpu: 19092670.000000488 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_AllPairs",
            "value": 78.5028106996759,
            "unit": "ns/iter",
            "extra": "iterations: 8926425\ncpu: 78.49318019251787 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_1",
            "value": 5779875.111134711,
            "unit": "ns/iter",
            "extra": "iterations: 126\ncpu: 5779787.634921488 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_2",
            "value": 4686378.6308725765,
            "unit": "ns/iter",
            "extra": "iterations: 149\ncpu: 4686485.476509368 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_2_3",
            "value": 5350288.2230707025,
            "unit": "ns/iter",
            "extra": "iterations: 130\ncpu: 5350314.307691607 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_3",
            "value": 1517841.7775377946,
            "unit": "ns/iter",
            "extra": "iterations: 463\ncpu: 1517642.3023756854 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_1",
            "value": 4.922658101520495,
            "unit": "ns/iter",
            "extra": "iterations: 142123315\ncpu: 4.922030280534933 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_2",
            "value": 4.928314909692151,
            "unit": "ns/iter",
            "extra": "iterations: 142108728\ncpu: 4.927919810808556 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_2_3",
            "value": 4.929312724298858,
            "unit": "ns/iter",
            "extra": "iterations: 142107061\ncpu: 4.928974338579866 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_3",
            "value": 4.221642174729155,
            "unit": "ns/iter",
            "extra": "iterations: 165815924\ncpu: 4.221130631579111 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Build_1",
            "value": 4234437.89089409,
            "unit": "ns/iter",
            "extra": "iterations: 165\ncpu: 4233912.9212114755 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Build_1",
            "value": 2204955.500000833,
            "unit": "ns/iter",
            "extra": "iterations: 318\ncpu: 2205286.7138385954 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Build_2",
            "value": 3550553.9394027363,
            "unit": "ns/iter",
            "extra": "iterations: 198\ncpu: 3549965.8535358817 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Build_2",
            "value": 1342302.8354196732,
            "unit": "ns/iter",
            "extra": "iterations: 480\ncpu: 1342481.0041666292 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Transpose_2_3",
            "value": 1855349.6530067045,
            "unit": "ns/iter",
            "extra": "iterations: 366\ncpu: 1854851.6803287333 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Transpose_2_3",
            "value": 484929.602373326,
            "unit": "ns/iter",
            "extra": "iterations: 1431\ncpu: 485211.2844165946 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Intersection_3_3_via_0",
            "value": 1323658.3527134436,
            "unit": "ns/iter",
            "extra": "iterations: 533\ncpu: 1323338.998123458 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Intersection_3_3_via_0",
            "value": 712553.9122813173,
            "unit": "ns/iter",
            "extra": "iterations: 969\ncpu: 712911.6749226653 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Intersection_3_3_via_2",
            "value": 1992484.192420259,
            "unit": "ns/iter",
            "extra": "iterations: 343\ncpu: 1991667.5889208233 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Intersection_3_3_via_2",
            "value": 412391.3955423662,
            "unit": "ns/iter",
            "extra": "iterations: 1704\ncpu: 412603.4906104098 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_AllPairs",
            "value": 32049732.59084909,
            "unit": "ns/iter",
            "extra": "iterations: 22\ncpu: 32047131.63636086 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_AllPairs",
            "value": 78.48957160025029,
            "unit": "ns/iter",
            "extra": "iterations: 8924188\ncpu: 78.47488656671159 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_1",
            "value": 8911429.012170363,
            "unit": "ns/iter",
            "extra": "iterations: 82\ncpu: 8910900.524390038 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_2",
            "value": 8870021.844140628,
            "unit": "ns/iter",
            "extra": "iterations: 77\ncpu: 8869309.935065147 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_2_3",
            "value": 10192196.957164664,
            "unit": "ns/iter",
            "extra": "iterations: 70\ncpu: 10192144.685715983 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_3",
            "value": 3054698.108711739,
            "unit": "ns/iter",
            "extra": "iterations: 230\ncpu: 3054342.3695642897 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_1",
            "value": 4.92421107261912,
            "unit": "ns/iter",
            "extra": "iterations: 142068906\ncpu: 4.923883224665597 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_2",
            "value": 4.927073954411028,
            "unit": "ns/iter",
            "extra": "iterations: 142124942\ncpu: 4.926731051893771 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_2_3",
            "value": 4.926471990723788,
            "unit": "ns/iter",
            "extra": "iterations: 142152849\ncpu: 4.926138715658136 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_3",
            "value": 4.222208077801973,
            "unit": "ns/iter",
            "extra": "iterations: 165671943\ncpu: 4.2219070612336935 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Build_1",
            "value": 6306458.499996855,
            "unit": "ns/iter",
            "extra": "iterations: 114\ncpu: 6305502.140351219 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Build_1",
            "value": 3643299.7083295505,
            "unit": "ns/iter",
            "extra": "iterations: 192\ncpu: 3643153.7499990915 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Build_2",
            "value": 6807211.019795522,
            "unit": "ns/iter",
            "extra": "iterations: 101\ncpu: 6806724.039604028 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Build_2",
            "value": 2390814.046672555,
            "unit": "ns/iter",
            "extra": "iterations: 300\ncpu: 2390991.696667015 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Transpose_2_3",
            "value": 3671024.0277999323,
            "unit": "ns/iter",
            "extra": "iterations: 180\ncpu: 3669774.8444453278 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Transpose_2_3",
            "value": 819729.4028297847,
            "unit": "ns/iter",
            "extra": "iterations: 849\ncpu: 819905.2167252348 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Intersection_3_3_via_0",
            "value": 2909049.34999124,
            "unit": "ns/iter",
            "extra": "iterations: 240\ncpu: 2908462.279167523 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Intersection_3_3_via_0",
            "value": 1632872.9509559793,
            "unit": "ns/iter",
            "extra": "iterations: 428\ncpu: 1633002.0490648802 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Intersection_3_3_via_2",
            "value": 3650643.3743549376,
            "unit": "ns/iter",
            "extra": "iterations: 195\ncpu: 3649313.6820516046 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Intersection_3_3_via_2",
            "value": 723244.3838964659,
            "unit": "ns/iter",
            "extra": "iterations: 969\ncpu: 723463.6130035622 ns\nthreads: 1"
          }
        ]
      }
    ]
  }
}