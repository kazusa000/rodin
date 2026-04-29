window.BENCHMARK_DATA = {
  "lastUpdate": 1777497559671,
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
      },
      {
        "commit": {
          "author": {
            "email": "jiajun.wang621@gmail.com",
            "name": "wang jiajun"
          },
          "committer": {
            "email": "jiajun.wang621@gmail.com",
            "name": "wang jiajun"
          },
          "distinct": true,
          "id": "3895b5739a85f0c4a9b71d67d288a15f06e73054",
          "message": "Add LevelSetStokes PETSc examples and runtime controls",
          "timestamp": "2026-04-15T17:57:13+02:00",
          "tree_id": "739582bf27b948c741df270765fdb16c71c0bd6e",
          "url": "https://github.com/kazusa000/rodin/commit/3895b5739a85f0c4a9b71d67d288a15f06e73054"
        },
        "date": 1776269340761,
        "tool": "googlecpp",
        "benches": [
          {
            "name": "P1Benchmark/UniformTriangular16_Build",
            "value": 0.31224543614075057,
            "unit": "ns/iter",
            "extra": "iterations: 2245046937\ncpu: 0.31224063891364423 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_Build",
            "value": 0.3114310262761153,
            "unit": "ns/iter",
            "extra": "iterations: 2121399139\ncpu: 0.3114015518604393 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular64_Build",
            "value": 0.31143957750738427,
            "unit": "ns/iter",
            "extra": "iterations: 2248080262\ncpu: 0.31141062080104687 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular128_Build",
            "value": 0.3116302967555085,
            "unit": "ns/iter",
            "extra": "iterations: 2249343181\ncpu: 0.3115984154487271 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/2D_Square_GridFunction_Projection_Real_SumOfComponents",
            "value": 597.9567668679115,
            "unit": "ns/iter",
            "extra": "iterations: 1179859\ncpu: 597.8974284215317 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular16_GridFunction_Projection_Real_SumOfComponents",
            "value": 121540.31636553003,
            "unit": "ns/iter",
            "extra": "iterations: 5756\ncpu: 121532.34833217514 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_GridFunction_Projection_Real_SumOfComponents",
            "value": 525684.5101579981,
            "unit": "ns/iter",
            "extra": "iterations: 1329\ncpu: 525623.3461249061 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/2D_Square_GridFunction_Projection_Vector_Components",
            "value": 992.313350514234,
            "unit": "ns/iter",
            "extra": "iterations: 705823\ncpu: 992.2214308686472 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular16_GridFunction_Projection_Vector_Components",
            "value": 204825.14123140846,
            "unit": "ns/iter",
            "extra": "iterations: 3427\ncpu: 204802.08958272546 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_GridFunction_Projection_Vector_Components",
            "value": 897255.9987261434,
            "unit": "ns/iter",
            "extra": "iterations: 785\ncpu: 897170.0445859868 ns\nthreads: 1"
          },
          {
            "name": "Poisson_UniformGrid_16x16/Assembly_NoCoefficient_ConstantSource",
            "value": 191851.90842191485,
            "unit": "ns/iter",
            "extra": "iterations: 3669\ncpu: 191842.71899700226 ns\nthreads: 1"
          },
          {
            "name": "Poisson_UniformGrid_16x16/Assembly_ConstantCoefficient_ConstantSource",
            "value": 191835.524807065,
            "unit": "ns/iter",
            "extra": "iterations: 3628\ncpu: 191814.01019845656 ns\nthreads: 1"
          },
          {
            "name": "MeshIO/Load_MEDIT_2D_Square",
            "value": 13896.966525120519,
            "unit": "ns/iter",
            "extra": "iterations: 50396\ncpu: 13896.265874275718 ns\nthreads: 1"
          },
          {
            "name": "MeshIO/Load_MEDIT_2D_UniformTriangular64",
            "value": 7807362.362636997,
            "unit": "ns/iter",
            "extra": "iterations: 91\ncpu: 7806379.758241747 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_16x16",
            "value": 143200.89922003064,
            "unit": "ns/iter",
            "extra": "iterations: 4872\ncpu: 143178.97988505737 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_64x64",
            "value": 2612291.516728591,
            "unit": "ns/iter",
            "extra": "iterations: 269\ncpu: 2612010.7211895892 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_128x128",
            "value": 13038727.568183037,
            "unit": "ns/iter",
            "extra": "iterations: 44\ncpu: 13036288.295454508 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_256x256",
            "value": 94007267.85713686,
            "unit": "ns/iter",
            "extra": "iterations: 7\ncpu: 93999552.28571396 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_512x512",
            "value": 377456965.5000164,
            "unit": "ns/iter",
            "extra": "iterations: 2\ncpu: 377404229.99999976 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_AllPairs",
            "value": 1258705.2513355012,
            "unit": "ns/iter",
            "extra": "iterations: 561\ncpu: 1258570.0855613921 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_AllPairs",
            "value": 19.501981166066763,
            "unit": "ns/iter",
            "extra": "iterations: 35411721\ncpu: 19.50133618188171 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_0_0",
            "value": 1245869.13214344,
            "unit": "ns/iter",
            "extra": "iterations: 560\ncpu: 1245813.9982142465 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_0_1",
            "value": 969980.7430941082,
            "unit": "ns/iter",
            "extra": "iterations: 724\ncpu: 969963.4116022177 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_1_0",
            "value": 468276.9032679814,
            "unit": "ns/iter",
            "extra": "iterations: 1468\ncpu: 468223.33719346026 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_1_1",
            "value": 973636.6166657085,
            "unit": "ns/iter",
            "extra": "iterations: 720\ncpu: 973580.4888888083 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_0_0",
            "value": 5.004560859999856,
            "unit": "ns/iter",
            "extra": "iterations: 100000000\ncpu: 5.004341739999987 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_0_1",
            "value": 4.981680713828496,
            "unit": "ns/iter",
            "extra": "iterations: 140509405\ncpu: 4.981344821722067 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_1_0",
            "value": 2.4909847927798183,
            "unit": "ns/iter",
            "extra": "iterations: 281056102\ncpu: 2.4908138482615194 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_1_1",
            "value": 4.668185943852157,
            "unit": "ns/iter",
            "extra": "iterations: 149993257\ncpu: 4.667892437324696 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_AllPairs",
            "value": 9055325.610391902,
            "unit": "ns/iter",
            "extra": "iterations: 77\ncpu: 9054890.97402677 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_AllPairs",
            "value": 46.51405275090084,
            "unit": "ns/iter",
            "extra": "iterations: 14963725\ncpu: 46.5101857324964 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_2_1",
            "value": 5551113.650792811,
            "unit": "ns/iter",
            "extra": "iterations: 126\ncpu: 5550895.611110868 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_1_2",
            "value": 6206012.705357113,
            "unit": "ns/iter",
            "extra": "iterations: 112\ncpu: 6205655.6785712 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_2_2",
            "value": 2488886.546100259,
            "unit": "ns/iter",
            "extra": "iterations: 282\ncpu: 2488805.851063639 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_1_1",
            "value": 7794012.788891299,
            "unit": "ns/iter",
            "extra": "iterations: 90\ncpu: 7793801.377777972 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_2_1",
            "value": 5.295270556064656,
            "unit": "ns/iter",
            "extra": "iterations: 132442420\ncpu: 5.2949398840642985 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_1_2",
            "value": 5.294558282754647,
            "unit": "ns/iter",
            "extra": "iterations: 132311248\ncpu: 5.294219641855436 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_2_2",
            "value": 4.673088785285226,
            "unit": "ns/iter",
            "extra": "iterations: 149865082\ncpu: 4.672763692879478 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_1_1",
            "value": 5.919234749607853,
            "unit": "ns/iter",
            "extra": "iterations: 118088930\ncpu: 5.918988740096159 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Build_1",
            "value": 3887512.4861884033,
            "unit": "ns/iter",
            "extra": "iterations: 181\ncpu: 3887079.3922647377 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Build_1",
            "value": 1410597.2500013877,
            "unit": "ns/iter",
            "extra": "iterations: 500\ncpu: 1410504.6220001043 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Transpose_1_2",
            "value": 2348868.3010024764,
            "unit": "ns/iter",
            "extra": "iterations: 299\ncpu: 2348627.6387963803 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Transpose_1_2",
            "value": 564629.3728955863,
            "unit": "ns/iter",
            "extra": "iterations: 1247\ncpu: 564613.3961507706 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Intersection_2_2_via_0",
            "value": 2211771.170490107,
            "unit": "ns/iter",
            "extra": "iterations: 305\ncpu: 2211591.501639544 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Intersection_2_2_via_0",
            "value": 1300482.296098718,
            "unit": "ns/iter",
            "extra": "iterations: 564\ncpu: 1300366.964538964 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Intersection_1_1_via_0",
            "value": 3584447.963918271,
            "unit": "ns/iter",
            "extra": "iterations: 194\ncpu: 3584209.278350722 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Intersection_1_1_via_0",
            "value": 1640786.7796177443,
            "unit": "ns/iter",
            "extra": "iterations: 422\ncpu: 1640684.964454893 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_AllPairs",
            "value": 5622014.47200141,
            "unit": "ns/iter",
            "extra": "iterations: 125\ncpu: 5621673.21600009 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_AllPairs",
            "value": 46.518653136621474,
            "unit": "ns/iter",
            "extra": "iterations: 15050257\ncpu: 46.51664798813729 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_2_1",
            "value": 3274227.481311525,
            "unit": "ns/iter",
            "extra": "iterations: 214\ncpu: 3273917.2149534407 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_1_2",
            "value": 3703234.640211098,
            "unit": "ns/iter",
            "extra": "iterations: 189\ncpu: 3703170.0582012245 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_2_2",
            "value": 1343773.2072946124,
            "unit": "ns/iter",
            "extra": "iterations: 521\ncpu: 1343609.3819575412 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_1_1",
            "value": 4582918.500001857,
            "unit": "ns/iter",
            "extra": "iterations: 154\ncpu: 4582535.24675336 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_2_1",
            "value": 5.295904190357508,
            "unit": "ns/iter",
            "extra": "iterations: 132356615\ncpu: 5.29556126076509 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_1_2",
            "value": 5.299515754606415,
            "unit": "ns/iter",
            "extra": "iterations: 132058045\ncpu: 5.299073698993443 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_2_2",
            "value": 4.682755421371409,
            "unit": "ns/iter",
            "extra": "iterations: 149913405\ncpu: 4.682457756195988 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_1_1",
            "value": 5.923309964509896,
            "unit": "ns/iter",
            "extra": "iterations: 118349704\ncpu: 5.922991324084819 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Build_1",
            "value": 2411675.162067954,
            "unit": "ns/iter",
            "extra": "iterations: 290\ncpu: 2411495.5793103143 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Build_1",
            "value": 907105.2539066925,
            "unit": "ns/iter",
            "extra": "iterations: 768\ncpu: 907094.9114582209 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Transpose_1_2",
            "value": 1481002.4236570855,
            "unit": "ns/iter",
            "extra": "iterations: 465\ncpu: 1480912.4150544098 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Transpose_1_2",
            "value": 365445.9314494664,
            "unit": "ns/iter",
            "extra": "iterations: 1911\ncpu: 365522.37048669596 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Intersection_2_2_via_0",
            "value": 1150404.1976733964,
            "unit": "ns/iter",
            "extra": "iterations: 602\ncpu: 1150327.2990035561 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Intersection_2_2_via_0",
            "value": 535256.1769576366,
            "unit": "ns/iter",
            "extra": "iterations: 1328\ncpu: 535330.0323794095 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Intersection_1_1_via_0",
            "value": 2106257.308383534,
            "unit": "ns/iter",
            "extra": "iterations: 334\ncpu: 2106052.988023568 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Intersection_1_1_via_0",
            "value": 842888.6603324435,
            "unit": "ns/iter",
            "extra": "iterations: 842\ncpu: 842826.9714962391 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_AllPairs",
            "value": 83597236.83334626,
            "unit": "ns/iter",
            "extra": "iterations: 6\ncpu: 83588855.66666885 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_AllPairs",
            "value": 86.12410008979514,
            "unit": "ns/iter",
            "extra": "iterations: 8144007\ncpu: 86.12030662547332 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_1",
            "value": 22256412.45161297,
            "unit": "ns/iter",
            "extra": "iterations: 31\ncpu: 22253625.93548328 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_2",
            "value": 25580485.17856751,
            "unit": "ns/iter",
            "extra": "iterations: 28\ncpu: 25578256.214284778 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_2_3",
            "value": 27697888.19230493,
            "unit": "ns/iter",
            "extra": "iterations: 26\ncpu: 27695586.461541887 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_3",
            "value": 12094705.655165438,
            "unit": "ns/iter",
            "extra": "iterations: 58\ncpu: 12094589.551723849 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_1",
            "value": 5.293612839528062,
            "unit": "ns/iter",
            "extra": "iterations: 132363435\ncpu: 5.293321482628479 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_2",
            "value": 5.292450257188465,
            "unit": "ns/iter",
            "extra": "iterations: 132394457\ncpu: 5.292258096575751 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_2_3",
            "value": 5.297337256294338,
            "unit": "ns/iter",
            "extra": "iterations: 132044815\ncpu: 5.297080638872492 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_3",
            "value": 4.673755016797985,
            "unit": "ns/iter",
            "extra": "iterations: 149703807\ncpu: 4.673416575170999 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Build_1",
            "value": 12011961.08332662,
            "unit": "ns/iter",
            "extra": "iterations: 60\ncpu: 12010988.500001218 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Build_1",
            "value": 6955544.227720215,
            "unit": "ns/iter",
            "extra": "iterations: 101\ncpu: 6954863.693070451 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Build_2",
            "value": 14658496.291673372,
            "unit": "ns/iter",
            "extra": "iterations: 48\ncpu: 14657214.187499434 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Build_2",
            "value": 5732862.975604788,
            "unit": "ns/iter",
            "extra": "iterations: 123\ncpu: 5732544.520324793 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Transpose_2_3",
            "value": 8050913.011231203,
            "unit": "ns/iter",
            "extra": "iterations: 89\ncpu: 8049932.797752258 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Transpose_2_3",
            "value": 1976838.7086847946,
            "unit": "ns/iter",
            "extra": "iterations: 357\ncpu: 1976739.100840621 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Intersection_3_3_via_0",
            "value": 11577390.499995442,
            "unit": "ns/iter",
            "extra": "iterations: 62\ncpu: 11575233.338709084 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Intersection_3_3_via_0",
            "value": 8655263.585357781,
            "unit": "ns/iter",
            "extra": "iterations: 82\ncpu: 8653936.121951332 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Intersection_3_3_via_2",
            "value": 8122927.170458897,
            "unit": "ns/iter",
            "extra": "iterations: 88\ncpu: 8121626.977273458 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Intersection_3_3_via_2",
            "value": 1913982.1810770973,
            "unit": "ns/iter",
            "extra": "iterations: 370\ncpu: 1913913.2459457642 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_AllPairs",
            "value": 18733468.00000736,
            "unit": "ns/iter",
            "extra": "iterations: 37\ncpu: 18731781.756757583 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_AllPairs",
            "value": 86.21499341703742,
            "unit": "ns/iter",
            "extra": "iterations: 8098941\ncpu: 86.20613336978282 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_1",
            "value": 5264321.789479453,
            "unit": "ns/iter",
            "extra": "iterations: 133\ncpu: 5264273.406015761 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_2",
            "value": 4593078.230262162,
            "unit": "ns/iter",
            "extra": "iterations: 152\ncpu: 4592895.407894858 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_2_3",
            "value": 5181226.455880673,
            "unit": "ns/iter",
            "extra": "iterations: 136\ncpu: 5181054.632352196 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_3",
            "value": 1598460.8627023576,
            "unit": "ns/iter",
            "extra": "iterations: 437\ncpu: 1598349.0389015754 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_1",
            "value": 5.291653921074693,
            "unit": "ns/iter",
            "extra": "iterations: 132437952\ncpu: 5.291155166760659 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_2",
            "value": 5.319073653394334,
            "unit": "ns/iter",
            "extra": "iterations: 131924443\ncpu: 5.318641049710612 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_2_3",
            "value": 5.2945654317175155,
            "unit": "ns/iter",
            "extra": "iterations: 132407132\ncpu: 5.294121860444813 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_3",
            "value": 4.668641491597812,
            "unit": "ns/iter",
            "extra": "iterations: 149883497\ncpu: 4.668391824351341 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Build_1",
            "value": 4102885.3333297153,
            "unit": "ns/iter",
            "extra": "iterations: 171\ncpu: 4102468.9941520835 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Build_1",
            "value": 2228088.1305708075,
            "unit": "ns/iter",
            "extra": "iterations: 314\ncpu: 2227907.009554001 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Build_2",
            "value": 3409751.921944695,
            "unit": "ns/iter",
            "extra": "iterations: 205\ncpu: 3409602.5951218577 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Build_2",
            "value": 1391715.1073542074,
            "unit": "ns/iter",
            "extra": "iterations: 503\ncpu: 1391641.5487076095 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Transpose_2_3",
            "value": 1791779.6751895072,
            "unit": "ns/iter",
            "extra": "iterations: 391\ncpu: 1791598.1508946042 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Transpose_2_3",
            "value": 496874.44389618764,
            "unit": "ns/iter",
            "extra": "iterations: 1417\ncpu: 496891.65067031665 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Intersection_3_3_via_0",
            "value": 1373626.729940313,
            "unit": "ns/iter",
            "extra": "iterations: 511\ncpu: 1373453.4559688705 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Intersection_3_3_via_0",
            "value": 802481.3525421476,
            "unit": "ns/iter",
            "extra": "iterations: 885\ncpu: 802456.9220344095 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Intersection_3_3_via_2",
            "value": 1867963.4986386236,
            "unit": "ns/iter",
            "extra": "iterations: 367\ncpu: 1867619.8038156421 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Intersection_3_3_via_2",
            "value": 469459.439810976,
            "unit": "ns/iter",
            "extra": "iterations: 1487\ncpu: 469528.45864129637 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_AllPairs",
            "value": 31177066.363635603,
            "unit": "ns/iter",
            "extra": "iterations: 22\ncpu: 31174193.909089867 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_AllPairs",
            "value": 86.23660726443437,
            "unit": "ns/iter",
            "extra": "iterations: 8115799\ncpu: 86.22722580980712 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_1",
            "value": 8545295.297624633,
            "unit": "ns/iter",
            "extra": "iterations: 84\ncpu: 8544570.011902612 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_2",
            "value": 8890287.935897935,
            "unit": "ns/iter",
            "extra": "iterations: 78\ncpu: 8889405.423076734 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_2_3",
            "value": 9748095.414288726,
            "unit": "ns/iter",
            "extra": "iterations: 70\ncpu: 9747291.785715057 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_3",
            "value": 3035954.644858351,
            "unit": "ns/iter",
            "extra": "iterations: 214\ncpu: 3035512.462616376 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_1",
            "value": 5.295703684614965,
            "unit": "ns/iter",
            "extra": "iterations: 131919063\ncpu: 5.295342652638538 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_2",
            "value": 5.297154374050507,
            "unit": "ns/iter",
            "extra": "iterations: 132227463\ncpu: 5.296367653972313 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_2_3",
            "value": 5.512815439067897,
            "unit": "ns/iter",
            "extra": "iterations: 130763487\ncpu: 5.512115794219935 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_3",
            "value": 4.676589930346163,
            "unit": "ns/iter",
            "extra": "iterations: 149902463\ncpu: 4.676317193000461 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Build_1",
            "value": 5807448.885242534,
            "unit": "ns/iter",
            "extra": "iterations: 122\ncpu: 5806577.024590051 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Build_1",
            "value": 3592081.209185543,
            "unit": "ns/iter",
            "extra": "iterations: 196\ncpu: 3591902.8061213936 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Build_2",
            "value": 6226473.292031187,
            "unit": "ns/iter",
            "extra": "iterations: 113\ncpu: 6226116.371682509 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Build_2",
            "value": 2383602.7585066776,
            "unit": "ns/iter",
            "extra": "iterations: 294\ncpu: 2383357.9931976437 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Transpose_2_3",
            "value": 3076838.4035043614,
            "unit": "ns/iter",
            "extra": "iterations: 228\ncpu: 3076337.982454453 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Transpose_2_3",
            "value": 803370.3606015416,
            "unit": "ns/iter",
            "extra": "iterations: 868\ncpu: 803328.6221198037 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Intersection_3_3_via_0",
            "value": 2809991.715416656,
            "unit": "ns/iter",
            "extra": "iterations: 253\ncpu: 2809486.972331264 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Intersection_3_3_via_0",
            "value": 1824871.7812493118,
            "unit": "ns/iter",
            "extra": "iterations: 384\ncpu: 1824704.6796873168 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Intersection_3_3_via_2",
            "value": 3285646.35922756,
            "unit": "ns/iter",
            "extra": "iterations: 206\ncpu: 3285147.7330107735 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Intersection_3_3_via_2",
            "value": 832627.3540701148,
            "unit": "ns/iter",
            "extra": "iterations: 836\ncpu: 832523.4784687035 ns\nthreads: 1"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "jiajun.wang621@gmail.com",
            "name": "wang jiajun"
          },
          "committer": {
            "email": "jiajun.wang621@gmail.com",
            "name": "wang jiajun"
          },
          "distinct": true,
          "id": "af148b5e4bf2c22695e3f637da56b4ac0bc6fb7f",
          "message": "Add new 3D LevelSetStokes shapes and fix fluid-mesh objectives",
          "timestamp": "2026-04-16T14:29:55+02:00",
          "tree_id": "c6f944935d440f2c903ccef60f7a2e718581c2c4",
          "url": "https://github.com/kazusa000/rodin/commit/af148b5e4bf2c22695e3f637da56b4ac0bc6fb7f"
        },
        "date": 1776343458766,
        "tool": "googlecpp",
        "benches": [
          {
            "name": "P1Benchmark/UniformTriangular16_Build",
            "value": 0.31154026541904334,
            "unit": "ns/iter",
            "extra": "iterations: 2236637457\ncpu: 0.31152074638621235 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_Build",
            "value": 0.3117107117995505,
            "unit": "ns/iter",
            "extra": "iterations: 2245341355\ncpu: 0.31164336613752874 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular64_Build",
            "value": 0.31113442158357113,
            "unit": "ns/iter",
            "extra": "iterations: 2249870980\ncpu: 0.3111205630111288 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular128_Build",
            "value": 0.3116934204203959,
            "unit": "ns/iter",
            "extra": "iterations: 2238980329\ncpu: 0.3116835203780835 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/2D_Square_GridFunction_Projection_Real_SumOfComponents",
            "value": 592.1713952289147,
            "unit": "ns/iter",
            "extra": "iterations: 1185103\ncpu: 592.1357611954404 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular16_GridFunction_Projection_Real_SumOfComponents",
            "value": 120662.47639784485,
            "unit": "ns/iter",
            "extra": "iterations: 5741\ncpu: 120654.31666956974 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_GridFunction_Projection_Real_SumOfComponents",
            "value": 524198.0911144311,
            "unit": "ns/iter",
            "extra": "iterations: 1328\ncpu: 524169.6814759039 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/2D_Square_GridFunction_Projection_Vector_Components",
            "value": 986.3030102626338,
            "unit": "ns/iter",
            "extra": "iterations: 709174\ncpu: 986.2640550838025 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular16_GridFunction_Projection_Vector_Components",
            "value": 204932.4176401964,
            "unit": "ns/iter",
            "extra": "iterations: 3424\ncpu: 204919.22984813084 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_GridFunction_Projection_Vector_Components",
            "value": 893687.1554139868,
            "unit": "ns/iter",
            "extra": "iterations: 785\ncpu: 893601.7248407638 ns\nthreads: 1"
          },
          {
            "name": "Poisson_UniformGrid_16x16/Assembly_NoCoefficient_ConstantSource",
            "value": 187688.84466019078,
            "unit": "ns/iter",
            "extra": "iterations: 3708\ncpu: 187681.01645091688 ns\nthreads: 1"
          },
          {
            "name": "Poisson_UniformGrid_16x16/Assembly_ConstantCoefficient_ConstantSource",
            "value": 190816.63801517675,
            "unit": "ns/iter",
            "extra": "iterations: 3688\ncpu: 190804.45851409974 ns\nthreads: 1"
          },
          {
            "name": "MeshIO/Load_MEDIT_2D_Square",
            "value": 13832.487700704774,
            "unit": "ns/iter",
            "extra": "iterations: 50572\ncpu: 13831.627303646279 ns\nthreads: 1"
          },
          {
            "name": "MeshIO/Load_MEDIT_2D_UniformTriangular64",
            "value": 7044873.149999944,
            "unit": "ns/iter",
            "extra": "iterations: 100\ncpu: 7044341.650000003 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_16x16",
            "value": 144303.48148910564,
            "unit": "ns/iter",
            "extra": "iterations: 4862\ncpu: 144292.93994241033 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_64x64",
            "value": 2609927.3358208793,
            "unit": "ns/iter",
            "extra": "iterations: 268\ncpu: 2609895.9328358183 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_128x128",
            "value": 13029440.31481464,
            "unit": "ns/iter",
            "extra": "iterations: 54\ncpu: 13028775.814814826 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_256x256",
            "value": 96624567.42857525,
            "unit": "ns/iter",
            "extra": "iterations: 7\ncpu: 96615656.28571467 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_512x512",
            "value": 355233533.5000123,
            "unit": "ns/iter",
            "extra": "iterations: 2\ncpu: 355179082.4999994 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_AllPairs",
            "value": 1245128.5275317812,
            "unit": "ns/iter",
            "extra": "iterations: 563\ncpu: 1244705.6873890404 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_AllPairs",
            "value": 19.4982833374091,
            "unit": "ns/iter",
            "extra": "iterations: 35513094\ncpu: 19.49752251944028 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_0_0",
            "value": 1245390.6589709632,
            "unit": "ns/iter",
            "extra": "iterations: 563\ncpu: 1245335.4262877216 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_0_1",
            "value": 965335.7662505672,
            "unit": "ns/iter",
            "extra": "iterations: 723\ncpu: 964810.8810512957 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_1_0",
            "value": 464086.4559311652,
            "unit": "ns/iter",
            "extra": "iterations: 1509\ncpu: 463957.4188205104 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_1_1",
            "value": 966144.9350840476,
            "unit": "ns/iter",
            "extra": "iterations: 724\ncpu: 966015.6477900032 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_0_0",
            "value": 4.989396488496178,
            "unit": "ns/iter",
            "extra": "iterations: 140468938\ncpu: 4.989206688527836 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_0_1",
            "value": 4.982165487444017,
            "unit": "ns/iter",
            "extra": "iterations: 140455647\ncpu: 4.981700030900121 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_1_0",
            "value": 2.491386123012411,
            "unit": "ns/iter",
            "extra": "iterations: 280837131\ncpu: 2.4911342902160625 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_1_1",
            "value": 4.671864268488471,
            "unit": "ns/iter",
            "extra": "iterations: 149911894\ncpu: 4.671468889586581 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_AllPairs",
            "value": 9059618.090908734,
            "unit": "ns/iter",
            "extra": "iterations: 77\ncpu: 9058920.363637105 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_AllPairs",
            "value": 46.60588616940614,
            "unit": "ns/iter",
            "extra": "iterations: 15040104\ncpu: 46.60343193105551 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_2_1",
            "value": 5542009.992064138,
            "unit": "ns/iter",
            "extra": "iterations: 126\ncpu: 5541493.682539762 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_1_2",
            "value": 6202068.672568589,
            "unit": "ns/iter",
            "extra": "iterations: 113\ncpu: 6201587.938053424 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_2_2",
            "value": 2443244.143356858,
            "unit": "ns/iter",
            "extra": "iterations: 286\ncpu: 2443169.1853150823 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_1_1",
            "value": 7766566.855557737,
            "unit": "ns/iter",
            "extra": "iterations: 90\ncpu: 7766083.155555887 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_2_1",
            "value": 5.293042697933195,
            "unit": "ns/iter",
            "extra": "iterations: 132419290\ncpu: 5.292625220993105 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_1_2",
            "value": 5.29620113109041,
            "unit": "ns/iter",
            "extra": "iterations: 132005529\ncpu: 5.29585034275345 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_2_2",
            "value": 4.670612875296656,
            "unit": "ns/iter",
            "extra": "iterations: 150011920\ncpu: 4.670311205936172 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_1_1",
            "value": 5.918685820565149,
            "unit": "ns/iter",
            "extra": "iterations: 118297240\ncpu: 5.9183562439833715 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Build_1",
            "value": 3881018.967032932,
            "unit": "ns/iter",
            "extra": "iterations: 182\ncpu: 3880816.2912088716 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Build_1",
            "value": 1401475.923075954,
            "unit": "ns/iter",
            "extra": "iterations: 481\ncpu: 1401457.2141373488 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Transpose_1_2",
            "value": 2332287.432341706,
            "unit": "ns/iter",
            "extra": "iterations: 303\ncpu: 2332043.8547852533 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Transpose_1_2",
            "value": 572052.2298281122,
            "unit": "ns/iter",
            "extra": "iterations: 1227\ncpu: 572032.5093723383 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Intersection_2_2_via_0",
            "value": 2250923.8424430997,
            "unit": "ns/iter",
            "extra": "iterations: 311\ncpu: 2250761.3697751598 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Intersection_2_2_via_0",
            "value": 1262915.801456411,
            "unit": "ns/iter",
            "extra": "iterations: 549\ncpu: 1262874.3151184903 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Intersection_1_1_via_0",
            "value": 3587034.198956332,
            "unit": "ns/iter",
            "extra": "iterations: 191\ncpu: 3586720.0837693796 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Intersection_1_1_via_0",
            "value": 1652018.830508084,
            "unit": "ns/iter",
            "extra": "iterations: 413\ncpu: 1651947.9128331146 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_AllPairs",
            "value": 5605491.741934451,
            "unit": "ns/iter",
            "extra": "iterations: 124\ncpu: 5605153.4354840135 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_AllPairs",
            "value": 46.49436707140852,
            "unit": "ns/iter",
            "extra": "iterations: 15070136\ncpu: 46.49213749630371 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_2_1",
            "value": 3288808.704226587,
            "unit": "ns/iter",
            "extra": "iterations: 213\ncpu: 3288673.366197411 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_1_2",
            "value": 3730391.4361731503,
            "unit": "ns/iter",
            "extra": "iterations: 188\ncpu: 3730284.005319453 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_2_2",
            "value": 1323203.855787488,
            "unit": "ns/iter",
            "extra": "iterations: 527\ncpu: 1323114.3036053712 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_1_1",
            "value": 4570257.320263415,
            "unit": "ns/iter",
            "extra": "iterations: 153\ncpu: 4570206.287581497 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_2_1",
            "value": 5.289776953659675,
            "unit": "ns/iter",
            "extra": "iterations: 132266954\ncpu: 5.28947532125071 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_1_2",
            "value": 5.295201258206561,
            "unit": "ns/iter",
            "extra": "iterations: 132031734\ncpu: 5.294683859866553 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_2_2",
            "value": 4.678702974986108,
            "unit": "ns/iter",
            "extra": "iterations: 150022550\ncpu: 4.67830724114477 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_1_1",
            "value": 5.915449409859107,
            "unit": "ns/iter",
            "extra": "iterations: 118392574\ncpu: 5.915220941137782 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Build_1",
            "value": 2449283.513983033,
            "unit": "ns/iter",
            "extra": "iterations: 286\ncpu: 2449021.300699257 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Build_1",
            "value": 926241.568210559,
            "unit": "ns/iter",
            "extra": "iterations: 755\ncpu: 926167.7271522157 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Transpose_1_2",
            "value": 1467181.565127337,
            "unit": "ns/iter",
            "extra": "iterations: 476\ncpu: 1467060.3613444464 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Transpose_1_2",
            "value": 365950.67030620866,
            "unit": "ns/iter",
            "extra": "iterations: 1923\ncpu: 366020.6318253442 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Intersection_2_2_via_0",
            "value": 1158417.5801641692,
            "unit": "ns/iter",
            "extra": "iterations: 605\ncpu: 1158313.6479342238 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Intersection_2_2_via_0",
            "value": 533565.4549692078,
            "unit": "ns/iter",
            "extra": "iterations: 1288\ncpu: 533591.4510867314 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Intersection_1_1_via_0",
            "value": 2080366.587535782,
            "unit": "ns/iter",
            "extra": "iterations: 337\ncpu: 2080145.2017802382 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Intersection_1_1_via_0",
            "value": 843275.9134038882,
            "unit": "ns/iter",
            "extra": "iterations: 843\ncpu: 843200.4531431072 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_AllPairs",
            "value": 81721478.11111093,
            "unit": "ns/iter",
            "extra": "iterations: 9\ncpu: 81713282.44444315 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_AllPairs",
            "value": 86.15314218850487,
            "unit": "ns/iter",
            "extra": "iterations: 8135890\ncpu: 86.14430320469023 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_1",
            "value": 21508050.36362761,
            "unit": "ns/iter",
            "extra": "iterations: 33\ncpu: 21505018.121212576 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_2",
            "value": 24345476.068962682,
            "unit": "ns/iter",
            "extra": "iterations: 29\ncpu: 24344461.827587534 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_2_3",
            "value": 26923828.576922935,
            "unit": "ns/iter",
            "extra": "iterations: 26\ncpu: 26920558.38461344 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_3",
            "value": 11078895.555556301,
            "unit": "ns/iter",
            "extra": "iterations: 63\ncpu: 11075921.31746169 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_1",
            "value": 5.292007551448926,
            "unit": "ns/iter",
            "extra": "iterations: 132338581\ncpu: 5.291311367468879 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_2",
            "value": 5.306231168847099,
            "unit": "ns/iter",
            "extra": "iterations: 132468325\ncpu: 5.305723870215787 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_2_3",
            "value": 5.293212154483271,
            "unit": "ns/iter",
            "extra": "iterations: 132451276\ncpu: 5.292841006680754 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_3",
            "value": 4.666933491379418,
            "unit": "ns/iter",
            "extra": "iterations: 149984858\ncpu: 4.66671290911241 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Build_1",
            "value": 11622328.133335222,
            "unit": "ns/iter",
            "extra": "iterations: 60\ncpu: 11621758.133334007 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Build_1",
            "value": 6989857.445544225,
            "unit": "ns/iter",
            "extra": "iterations: 101\ncpu: 6989624.198020083 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Build_2",
            "value": 14224835.285711389,
            "unit": "ns/iter",
            "extra": "iterations: 49\ncpu: 14223675.591837576 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Build_2",
            "value": 5682917.206611308,
            "unit": "ns/iter",
            "extra": "iterations: 121\ncpu: 5682644.15702534 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Transpose_2_3",
            "value": 7413953.829786153,
            "unit": "ns/iter",
            "extra": "iterations: 94\ncpu: 7412556.851064024 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Transpose_2_3",
            "value": 1907641.4347812964,
            "unit": "ns/iter",
            "extra": "iterations: 368\ncpu: 1907542.1385872273 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Intersection_3_3_via_0",
            "value": 10863947.707693674,
            "unit": "ns/iter",
            "extra": "iterations: 65\ncpu: 10862688.353846623 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Intersection_3_3_via_0",
            "value": 7968285.511361599,
            "unit": "ns/iter",
            "extra": "iterations: 88\ncpu: 7967200.943182333 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Intersection_3_3_via_2",
            "value": 7674380.945655047,
            "unit": "ns/iter",
            "extra": "iterations: 92\ncpu: 7673904.4673913615 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Intersection_3_3_via_2",
            "value": 1886383.3844080858,
            "unit": "ns/iter",
            "extra": "iterations: 372\ncpu: 1886316.3306452855 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_AllPairs",
            "value": 19117246.416668206,
            "unit": "ns/iter",
            "extra": "iterations: 36\ncpu: 19116992.249998473 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_AllPairs",
            "value": 86.19749678161331,
            "unit": "ns/iter",
            "extra": "iterations: 8140720\ncpu: 86.19381111252949 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_1",
            "value": 5660711.346773577,
            "unit": "ns/iter",
            "extra": "iterations: 124\ncpu: 5660307.935484171 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_2",
            "value": 5021659.728575401,
            "unit": "ns/iter",
            "extra": "iterations: 140\ncpu: 5021356.885713959 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_2_3",
            "value": 5586816.064516326,
            "unit": "ns/iter",
            "extra": "iterations: 124\ncpu: 5586627.016128199 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_3",
            "value": 1580067.4469533202,
            "unit": "ns/iter",
            "extra": "iterations: 443\ncpu: 1579929.541760911 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_1",
            "value": 5.288163179527996,
            "unit": "ns/iter",
            "extra": "iterations: 131910812\ncpu: 5.287896635796626 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_2",
            "value": 5.2922728631831735,
            "unit": "ns/iter",
            "extra": "iterations: 131797744\ncpu: 5.291843037920291 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_2_3",
            "value": 5.293614000419253,
            "unit": "ns/iter",
            "extra": "iterations: 132300040\ncpu: 5.293374484240563 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_3",
            "value": 4.6716059220495785,
            "unit": "ns/iter",
            "extra": "iterations: 149808006\ncpu: 4.6713983496983165 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Build_1",
            "value": 4717170.899329826,
            "unit": "ns/iter",
            "extra": "iterations: 149\ncpu: 4716982.087248617 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Build_1",
            "value": 2697992.373076265,
            "unit": "ns/iter",
            "extra": "iterations: 260\ncpu: 2697833.476922977 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Build_2",
            "value": 4134198.6272155037,
            "unit": "ns/iter",
            "extra": "iterations: 169\ncpu: 4133945.159762876 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Build_2",
            "value": 1843195.1023609238,
            "unit": "ns/iter",
            "extra": "iterations: 381\ncpu: 1843110.6955383436 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Transpose_2_3",
            "value": 1735356.3876554288,
            "unit": "ns/iter",
            "extra": "iterations: 405\ncpu: 1735118.5432105777 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Transpose_2_3",
            "value": 494002.7191008815,
            "unit": "ns/iter",
            "extra": "iterations: 1424\ncpu: 493976.1903089143 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Intersection_3_3_via_0",
            "value": 1361062.8958746914,
            "unit": "ns/iter",
            "extra": "iterations: 509\ncpu: 1361009.4361485073 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Intersection_3_3_via_0",
            "value": 796471.5671649441,
            "unit": "ns/iter",
            "extra": "iterations: 871\ncpu: 796480.5074631431 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Intersection_3_3_via_2",
            "value": 1821110.1727761175,
            "unit": "ns/iter",
            "extra": "iterations: 382\ncpu: 1820916.445025645 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Intersection_3_3_via_2",
            "value": 469127.46684510366,
            "unit": "ns/iter",
            "extra": "iterations: 1493\ncpu: 469164.4012054501 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_AllPairs",
            "value": 32297638.63636363,
            "unit": "ns/iter",
            "extra": "iterations: 22\ncpu: 32294655.63636066 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_AllPairs",
            "value": 86.11831596671058,
            "unit": "ns/iter",
            "extra": "iterations: 8131692\ncpu: 86.1097316524064 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_1",
            "value": 9007943.076923288,
            "unit": "ns/iter",
            "extra": "iterations: 78\ncpu: 9007345.141026527 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_2",
            "value": 9580748.8356141,
            "unit": "ns/iter",
            "extra": "iterations: 73\ncpu: 9579718.780823803 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_2_3",
            "value": 10523245.348481983,
            "unit": "ns/iter",
            "extra": "iterations: 66\ncpu: 10521916.742425133 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_3",
            "value": 3040417.7339066,
            "unit": "ns/iter",
            "extra": "iterations: 233\ncpu: 3040154.682402587 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_1",
            "value": 5.305519046726196,
            "unit": "ns/iter",
            "extra": "iterations: 100804383\ncpu: 5.305216143230524 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_2",
            "value": 5.289633960166933,
            "unit": "ns/iter",
            "extra": "iterations: 132399443\ncpu: 5.289360590436901 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_2_3",
            "value": 5.292067822992227,
            "unit": "ns/iter",
            "extra": "iterations: 132407547\ncpu: 5.291558690381839 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_3",
            "value": 4.714340118267377,
            "unit": "ns/iter",
            "extra": "iterations: 150012654\ncpu: 4.714278183492399 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Build_1",
            "value": 7027797.424243981,
            "unit": "ns/iter",
            "extra": "iterations: 99\ncpu: 7027282.3434347315 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Build_1",
            "value": 4533401.88311788,
            "unit": "ns/iter",
            "extra": "iterations: 154\ncpu: 4533440.285713707 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Build_2",
            "value": 7301718.697915405,
            "unit": "ns/iter",
            "extra": "iterations: 96\ncpu: 7301150.166666067 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Build_2",
            "value": 3018748.0565231177,
            "unit": "ns/iter",
            "extra": "iterations: 230\ncpu: 3018644.8260886446 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Transpose_2_3",
            "value": 2910583.8298725956,
            "unit": "ns/iter",
            "extra": "iterations: 241\ncpu: 2910441.7427387843 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Transpose_2_3",
            "value": 806524.4994243138,
            "unit": "ns/iter",
            "extra": "iterations: 869\ncpu: 806495.4407361989 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Intersection_3_3_via_0",
            "value": 2712388.5984534114,
            "unit": "ns/iter",
            "extra": "iterations: 259\ncpu: 2712207.3397673685 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Intersection_3_3_via_0",
            "value": 1833842.9107602725,
            "unit": "ns/iter",
            "extra": "iterations: 381\ncpu: 1833747.1522315277 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Intersection_3_3_via_2",
            "value": 3121857.433033678,
            "unit": "ns/iter",
            "extra": "iterations: 224\ncpu: 3121623.526784388 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Intersection_3_3_via_2",
            "value": 830312.1523818884,
            "unit": "ns/iter",
            "extra": "iterations: 840\ncpu: 830282.0285705149 ns\nthreads: 1"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "jiajun.wang621@gmail.com",
            "name": "wang jiajun"
          },
          "committer": {
            "email": "jiajun.wang621@gmail.com",
            "name": "wang jiajun"
          },
          "distinct": true,
          "id": "9e31d4f3fa3cdedd8739f10452e44951c8dc83fe",
          "message": "Unify LevelSetStokes 3D versions",
          "timestamp": "2026-04-21T15:58:40+02:00",
          "tree_id": "3f53ed2c892078e079fb7f23e30087d98a3927e4",
          "url": "https://github.com/kazusa000/rodin/commit/9e31d4f3fa3cdedd8739f10452e44951c8dc83fe"
        },
        "date": 1776780417438,
        "tool": "googlecpp",
        "benches": [
          {
            "name": "P1Benchmark/UniformTriangular16_Build",
            "value": 0.3128824425947241,
            "unit": "ns/iter",
            "extra": "iterations: 2148580954\ncpu: 0.3126962634334139 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_Build",
            "value": 0.3114409984895785,
            "unit": "ns/iter",
            "extra": "iterations: 2249392448\ncpu: 0.3114241063727445 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular64_Build",
            "value": 0.3111697589298059,
            "unit": "ns/iter",
            "extra": "iterations: 2247249252\ncpu: 0.31115732506204435 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular128_Build",
            "value": 0.3114956734736527,
            "unit": "ns/iter",
            "extra": "iterations: 2208990370\ncpu: 0.31145698113659065 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/2D_Square_GridFunction_Projection_Real_SumOfComponents",
            "value": 597.8283677718488,
            "unit": "ns/iter",
            "extra": "iterations: 1179155\ncpu: 597.7740102022204 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular16_GridFunction_Projection_Real_SumOfComponents",
            "value": 120772.35624244306,
            "unit": "ns/iter",
            "extra": "iterations: 5791\ncpu: 120759.96788119494 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_GridFunction_Projection_Real_SumOfComponents",
            "value": 524639.6881559058,
            "unit": "ns/iter",
            "extra": "iterations: 1334\ncpu: 524613.7196401795 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/2D_Square_GridFunction_Projection_Vector_Components",
            "value": 1009.402390893957,
            "unit": "ns/iter",
            "extra": "iterations: 708187\ncpu: 1009.3887913785481 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular16_GridFunction_Projection_Vector_Components",
            "value": 204939.2810668257,
            "unit": "ns/iter",
            "extra": "iterations: 3412\ncpu: 204933.64361078563 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_GridFunction_Projection_Vector_Components",
            "value": 890254.8184143698,
            "unit": "ns/iter",
            "extra": "iterations: 782\ncpu: 890219.8081841425 ns\nthreads: 1"
          },
          {
            "name": "Poisson_UniformGrid_16x16/Assembly_NoCoefficient_ConstantSource",
            "value": 190471.86786786714,
            "unit": "ns/iter",
            "extra": "iterations: 3663\ncpu: 190444.30466830454 ns\nthreads: 1"
          },
          {
            "name": "Poisson_UniformGrid_16x16/Assembly_ConstantCoefficient_ConstantSource",
            "value": 195186.37142069172,
            "unit": "ns/iter",
            "extra": "iterations: 3632\ncpu: 195174.52560572643 ns\nthreads: 1"
          },
          {
            "name": "MeshIO/Load_MEDIT_2D_Square",
            "value": 13879.867924156222,
            "unit": "ns/iter",
            "extra": "iterations: 50683\ncpu: 13876.010910956338 ns\nthreads: 1"
          },
          {
            "name": "MeshIO/Load_MEDIT_2D_UniformTriangular64",
            "value": 7431927.923077042,
            "unit": "ns/iter",
            "extra": "iterations: 91\ncpu: 7431414.901098893 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_16x16",
            "value": 144152.5337309754,
            "unit": "ns/iter",
            "extra": "iterations: 4862\ncpu: 144137.3377211024 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_64x64",
            "value": 2618722.8501873068,
            "unit": "ns/iter",
            "extra": "iterations: 267\ncpu: 2618447.3745318297 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_128x128",
            "value": 13852565.660000665,
            "unit": "ns/iter",
            "extra": "iterations: 50\ncpu: 13851342.680000016 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_256x256",
            "value": 95557293.4285718,
            "unit": "ns/iter",
            "extra": "iterations: 7\ncpu: 95536608.14285742 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_512x512",
            "value": 389489407.9999744,
            "unit": "ns/iter",
            "extra": "iterations: 2\ncpu: 389453929.5000001 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_AllPairs",
            "value": 1248665.3071433245,
            "unit": "ns/iter",
            "extra": "iterations: 560\ncpu: 1248647.4178571391 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_AllPairs",
            "value": 19.581107256207506,
            "unit": "ns/iter",
            "extra": "iterations: 35878399\ncpu: 19.5763623120419 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_0_0",
            "value": 1242851.2477858972,
            "unit": "ns/iter",
            "extra": "iterations: 565\ncpu: 1242480.7327434388 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_0_1",
            "value": 964206.2731034552,
            "unit": "ns/iter",
            "extra": "iterations: 725\ncpu: 964185.3917241028 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_1_0",
            "value": 471132.2394936775,
            "unit": "ns/iter",
            "extra": "iterations: 1499\ncpu: 471109.05870575545 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_1_1",
            "value": 968726.676838422,
            "unit": "ns/iter",
            "extra": "iterations: 721\ncpu: 968568.0790569405 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_0_0",
            "value": 4.988202575724925,
            "unit": "ns/iter",
            "extra": "iterations: 140540847\ncpu: 4.987671384960413 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_0_1",
            "value": 4.98796261242635,
            "unit": "ns/iter",
            "extra": "iterations: 140535892\ncpu: 4.987888247082106 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_1_0",
            "value": 2.49227970244027,
            "unit": "ns/iter",
            "extra": "iterations: 280930234\ncpu: 2.4920298752892527 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_1_1",
            "value": 4.67736137298186,
            "unit": "ns/iter",
            "extra": "iterations: 149819947\ncpu: 4.676863849110829 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_AllPairs",
            "value": 9021503.07791919,
            "unit": "ns/iter",
            "extra": "iterations: 77\ncpu: 9020704.11688287 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_AllPairs",
            "value": 46.52653062933914,
            "unit": "ns/iter",
            "extra": "iterations: 15044008\ncpu: 46.52152451660454 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_2_1",
            "value": 5534771.8888903,
            "unit": "ns/iter",
            "extra": "iterations: 126\ncpu: 5534484.507936133 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_1_2",
            "value": 6257029.362833253,
            "unit": "ns/iter",
            "extra": "iterations: 113\ncpu: 6256571.575221289 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_2_2",
            "value": 2488136.149999929,
            "unit": "ns/iter",
            "extra": "iterations: 280\ncpu: 2487905.5000000047 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_1_1",
            "value": 7944443.494385738,
            "unit": "ns/iter",
            "extra": "iterations: 89\ncpu: 7944063.516854152 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_2_1",
            "value": 5.402978685155224,
            "unit": "ns/iter",
            "extra": "iterations: 132267487\ncpu: 5.402709004367796 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_1_2",
            "value": 5.424801986471022,
            "unit": "ns/iter",
            "extra": "iterations: 128141623\ncpu: 5.424138782759153 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_2_2",
            "value": 4.779781689129037,
            "unit": "ns/iter",
            "extra": "iterations: 146673044\ncpu: 4.779463552962062 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_1_1",
            "value": 6.06864780818565,
            "unit": "ns/iter",
            "extra": "iterations: 113961672\ncpu: 6.067218432878043 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Build_1",
            "value": 3851892.0219785487,
            "unit": "ns/iter",
            "extra": "iterations: 182\ncpu: 3850993.4505494167 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Build_1",
            "value": 1422875.0740001032,
            "unit": "ns/iter",
            "extra": "iterations: 500\ncpu: 1422609.1679998804 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Transpose_1_2",
            "value": 2406823.8807957782,
            "unit": "ns/iter",
            "extra": "iterations: 302\ncpu: 2406528.0099337217 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Transpose_1_2",
            "value": 567293.286740906,
            "unit": "ns/iter",
            "extra": "iterations: 1252\ncpu: 567203.0734823712 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Intersection_2_2_via_0",
            "value": 2230273.7760235597,
            "unit": "ns/iter",
            "extra": "iterations: 317\ncpu: 2230030.9337540837 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Intersection_2_2_via_0",
            "value": 1259940.0996448505,
            "unit": "ns/iter",
            "extra": "iterations: 562\ncpu: 1259858.5943060212 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Intersection_1_1_via_0",
            "value": 3625639.6062150802,
            "unit": "ns/iter",
            "extra": "iterations: 193\ncpu: 3625208.5077718087 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Intersection_1_1_via_0",
            "value": 1646336.244705387,
            "unit": "ns/iter",
            "extra": "iterations: 425\ncpu: 1645778.2588234441 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_AllPairs",
            "value": 5626060.791998043,
            "unit": "ns/iter",
            "extra": "iterations: 125\ncpu: 5624782.431999961 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_AllPairs",
            "value": 46.55614896299156,
            "unit": "ns/iter",
            "extra": "iterations: 15056921\ncpu: 46.55215232915157 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_2_1",
            "value": 3282513.6150280177,
            "unit": "ns/iter",
            "extra": "iterations: 213\ncpu: 3282429.352112616 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_1_2",
            "value": 3747822.569151145,
            "unit": "ns/iter",
            "extra": "iterations: 188\ncpu: 3747366.6595743494 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_2_2",
            "value": 1354568.6473980558,
            "unit": "ns/iter",
            "extra": "iterations: 519\ncpu: 1353661.2061658024 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_1_1",
            "value": 4596612.901314762,
            "unit": "ns/iter",
            "extra": "iterations: 152\ncpu: 4596451.217105724 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_2_1",
            "value": 5.298102109795057,
            "unit": "ns/iter",
            "extra": "iterations: 131358916\ncpu: 5.297267937259713 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_1_2",
            "value": 5.295782397313782,
            "unit": "ns/iter",
            "extra": "iterations: 132371454\ncpu: 5.295261053791814 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_2_2",
            "value": 4.669086629204451,
            "unit": "ns/iter",
            "extra": "iterations: 149928354\ncpu: 4.668555849015737 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_1_1",
            "value": 5.917373720642906,
            "unit": "ns/iter",
            "extra": "iterations: 118323239\ncpu: 5.916864234928463 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Build_1",
            "value": 2402130.0481098387,
            "unit": "ns/iter",
            "extra": "iterations: 291\ncpu: 2401814.178693412 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Build_1",
            "value": 916524.1484889674,
            "unit": "ns/iter",
            "extra": "iterations: 761\ncpu: 916487.6044678973 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Transpose_1_2",
            "value": 1480325.3873683866,
            "unit": "ns/iter",
            "extra": "iterations: 475\ncpu: 1480164.4715790094 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Transpose_1_2",
            "value": 366685.1257926196,
            "unit": "ns/iter",
            "extra": "iterations: 1892\ncpu: 366740.19873162557 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Intersection_2_2_via_0",
            "value": 1165112.3199991768,
            "unit": "ns/iter",
            "extra": "iterations: 600\ncpu: 1164885.5583337317 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Intersection_2_2_via_0",
            "value": 549696.2168218,
            "unit": "ns/iter",
            "extra": "iterations: 1296\ncpu: 549672.881944428 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Intersection_1_1_via_0",
            "value": 2094678.0991017453,
            "unit": "ns/iter",
            "extra": "iterations: 333\ncpu: 2094318.8708704344 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Intersection_1_1_via_0",
            "value": 836798.4793197325,
            "unit": "ns/iter",
            "extra": "iterations: 822\ncpu: 836663.8917275427 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_AllPairs",
            "value": 80587585.00000194,
            "unit": "ns/iter",
            "extra": "iterations: 9\ncpu: 80579782.22222395 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_AllPairs",
            "value": 86.61263600589214,
            "unit": "ns/iter",
            "extra": "iterations: 8133931\ncpu: 86.61016303187247 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_1",
            "value": 21729409.437496018,
            "unit": "ns/iter",
            "extra": "iterations: 32\ncpu: 21727767.78124952 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_2",
            "value": 25237106.214289285,
            "unit": "ns/iter",
            "extra": "iterations: 28\ncpu: 25236029.71428624 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_2_3",
            "value": 27577039.1199967,
            "unit": "ns/iter",
            "extra": "iterations: 25\ncpu: 27573962.67999866 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_3",
            "value": 11466004.566665333,
            "unit": "ns/iter",
            "extra": "iterations: 60\ncpu: 11464821.616666162 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_1",
            "value": 5.290252838806396,
            "unit": "ns/iter",
            "extra": "iterations: 132068176\ncpu: 5.2900767100772805 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_2",
            "value": 5.294416903347364,
            "unit": "ns/iter",
            "extra": "iterations: 132095116\ncpu: 5.293889026146852 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_2_3",
            "value": 5.295279920395772,
            "unit": "ns/iter",
            "extra": "iterations: 132208194\ncpu: 5.2950732161124705 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_3",
            "value": 4.674563930223061,
            "unit": "ns/iter",
            "extra": "iterations: 150000650\ncpu: 4.673902426422794 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Build_1",
            "value": 11416995.087718837,
            "unit": "ns/iter",
            "extra": "iterations: 57\ncpu: 11416755.929824008 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Build_1",
            "value": 6919087.861383492,
            "unit": "ns/iter",
            "extra": "iterations: 101\ncpu: 6917975.623762536 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Build_2",
            "value": 13882201.50979979,
            "unit": "ns/iter",
            "extra": "iterations: 51\ncpu: 13881616.647058407 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Build_2",
            "value": 5669872.679998662,
            "unit": "ns/iter",
            "extra": "iterations: 125\ncpu: 5669316.295999806 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Transpose_2_3",
            "value": 7800154.11111271,
            "unit": "ns/iter",
            "extra": "iterations: 90\ncpu: 7799023.611111503 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Transpose_2_3",
            "value": 1908688.473974495,
            "unit": "ns/iter",
            "extra": "iterations: 365\ncpu: 1908573.994520823 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Intersection_3_3_via_0",
            "value": 11094632.873012178,
            "unit": "ns/iter",
            "extra": "iterations: 63\ncpu: 11091988.714285087 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Intersection_3_3_via_0",
            "value": 8030998.0919550415,
            "unit": "ns/iter",
            "extra": "iterations: 87\ncpu: 8030930.19540264 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Intersection_3_3_via_2",
            "value": 7772407.752807422,
            "unit": "ns/iter",
            "extra": "iterations: 89\ncpu: 7770869.202246899 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Intersection_3_3_via_2",
            "value": 1862296.375662107,
            "unit": "ns/iter",
            "extra": "iterations: 378\ncpu: 1862088.1322751755 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_AllPairs",
            "value": 17955222.000004586,
            "unit": "ns/iter",
            "extra": "iterations: 39\ncpu: 17954814.30769049 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_AllPairs",
            "value": 86.31544610891821,
            "unit": "ns/iter",
            "extra": "iterations: 8127309\ncpu: 86.30693615808204 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_1",
            "value": 5180058.301472122,
            "unit": "ns/iter",
            "extra": "iterations: 136\ncpu: 5179648.139706551 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_2",
            "value": 4491626.974196048,
            "unit": "ns/iter",
            "extra": "iterations: 155\ncpu: 4490800.774194136 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_2_3",
            "value": 5096568.1956544835,
            "unit": "ns/iter",
            "extra": "iterations: 138\ncpu: 5096665.746377358 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_3",
            "value": 1573129.60673983,
            "unit": "ns/iter",
            "extra": "iterations: 445\ncpu: 1573039.9842701328 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_1",
            "value": 5.2894244235204875,
            "unit": "ns/iter",
            "extra": "iterations: 132405671\ncpu: 5.288940025839292 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_2",
            "value": 5.2946570417869685,
            "unit": "ns/iter",
            "extra": "iterations: 131467620\ncpu: 5.294039155801271 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_2_3",
            "value": 5.292531353400674,
            "unit": "ns/iter",
            "extra": "iterations: 132037751\ncpu: 5.292195055639835 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_3",
            "value": 4.6685056425012235,
            "unit": "ns/iter",
            "extra": "iterations: 149824336\ncpu: 4.668221095937286 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Build_1",
            "value": 4041815.8728335723,
            "unit": "ns/iter",
            "extra": "iterations: 173\ncpu: 4041469.62427833 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Build_1",
            "value": 2205130.8301873836,
            "unit": "ns/iter",
            "extra": "iterations: 318\ncpu: 2204989.610062917 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Build_2",
            "value": 3350478.301438727,
            "unit": "ns/iter",
            "extra": "iterations: 209\ncpu: 3350134.5119609395 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Build_2",
            "value": 1362459.1378647361,
            "unit": "ns/iter",
            "extra": "iterations: 515\ncpu: 1362442.8601946246 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Transpose_2_3",
            "value": 1780402.1072326312,
            "unit": "ns/iter",
            "extra": "iterations: 401\ncpu: 1780077.0872820425 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Transpose_2_3",
            "value": 486858.06123863533,
            "unit": "ns/iter",
            "extra": "iterations: 1437\ncpu: 486803.109255621 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Intersection_3_3_via_0",
            "value": 1362776.0136985914,
            "unit": "ns/iter",
            "extra": "iterations: 511\ncpu: 1362585.015655754 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Intersection_3_3_via_0",
            "value": 794057.1936615282,
            "unit": "ns/iter",
            "extra": "iterations: 852\ncpu: 794124.9530514311 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Intersection_3_3_via_2",
            "value": 1815937.5958552002,
            "unit": "ns/iter",
            "extra": "iterations: 386\ncpu: 1815739.83160739 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Intersection_3_3_via_2",
            "value": 472266.5621613174,
            "unit": "ns/iter",
            "extra": "iterations: 1480\ncpu: 472311.1033783119 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_AllPairs",
            "value": 30001819.304353062,
            "unit": "ns/iter",
            "extra": "iterations: 23\ncpu: 30001146.91304297 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_AllPairs",
            "value": 86.27165609279787,
            "unit": "ns/iter",
            "extra": "iterations: 7504978\ncpu: 86.25986378640845 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_1",
            "value": 8037245.022987578,
            "unit": "ns/iter",
            "extra": "iterations: 87\ncpu: 8036640.758618482 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_2",
            "value": 8607003.91463527,
            "unit": "ns/iter",
            "extra": "iterations: 82\ncpu: 8605992.365853969 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_2_3",
            "value": 9384260.640000168,
            "unit": "ns/iter",
            "extra": "iterations: 75\ncpu: 9383803.066668104 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_3",
            "value": 3000391.545063983,
            "unit": "ns/iter",
            "extra": "iterations: 233\ncpu: 2999914.163089422 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_1",
            "value": 5.309572744794473,
            "unit": "ns/iter",
            "extra": "iterations: 132369224\ncpu: 5.309215229667097 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_2",
            "value": 5.2906399471537116,
            "unit": "ns/iter",
            "extra": "iterations: 131723658\ncpu: 5.290366549036971 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_2_3",
            "value": 5.295258780350535,
            "unit": "ns/iter",
            "extra": "iterations: 132406247\ncpu: 5.294764664691376 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_3",
            "value": 4.666643810489936,
            "unit": "ns/iter",
            "extra": "iterations: 145722826\ncpu: 4.666483760066388 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Build_1",
            "value": 5718633.504062193,
            "unit": "ns/iter",
            "extra": "iterations: 123\ncpu: 5718145.780487871 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Build_1",
            "value": 3569195.841836047,
            "unit": "ns/iter",
            "extra": "iterations: 196\ncpu: 3568992.1479601306 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Build_2",
            "value": 6133797.1217401745,
            "unit": "ns/iter",
            "extra": "iterations: 115\ncpu: 6133081.165217924 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Build_2",
            "value": 2348043.1638810374,
            "unit": "ns/iter",
            "extra": "iterations: 299\ncpu: 2347837.428095089 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Transpose_2_3",
            "value": 2994253.7711857576,
            "unit": "ns/iter",
            "extra": "iterations: 236\ncpu: 2993486.5805079467 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Transpose_2_3",
            "value": 806302.1843651027,
            "unit": "ns/iter",
            "extra": "iterations: 857\ncpu: 806227.7491248653 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Intersection_3_3_via_0",
            "value": 2722780.2343752217,
            "unit": "ns/iter",
            "extra": "iterations: 256\ncpu: 2722582.550780839 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Intersection_3_3_via_0",
            "value": 1818408.7722513285,
            "unit": "ns/iter",
            "extra": "iterations: 382\ncpu: 1818317.4554978912 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Intersection_3_3_via_2",
            "value": 3178995.162893104,
            "unit": "ns/iter",
            "extra": "iterations: 221\ncpu: 3178680.416290129 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Intersection_3_3_via_2",
            "value": 831351.5758617405,
            "unit": "ns/iter",
            "extra": "iterations: 870\ncpu: 831282.7126441734 ns\nthreads: 1"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "jiajun.wang621@gmail.com",
            "name": "wang jiajun"
          },
          "committer": {
            "email": "jiajun.wang621@gmail.com",
            "name": "wang jiajun"
          },
          "distinct": true,
          "id": "ba52f2b7530cbcfee1406caabeec8eaf205c92ec",
          "message": "Add 3Dv14 adapt and 3Dv15 thickness variants",
          "timestamp": "2026-04-29T23:08:30+02:00",
          "tree_id": "eeef1438c45693911e28a6a59c52bb1bf84915f9",
          "url": "https://github.com/kazusa000/rodin/commit/ba52f2b7530cbcfee1406caabeec8eaf205c92ec"
        },
        "date": 1777497556446,
        "tool": "googlecpp",
        "benches": [
          {
            "name": "P1Benchmark/UniformTriangular16_Build",
            "value": 0.32337925240025933,
            "unit": "ns/iter",
            "extra": "iterations: 2245079663\ncpu: 0.32334981691916914 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_Build",
            "value": 0.31133339144873134,
            "unit": "ns/iter",
            "extra": "iterations: 2246840208\ncpu: 0.3113241691640583 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular64_Build",
            "value": 0.3113384923901724,
            "unit": "ns/iter",
            "extra": "iterations: 2246650315\ncpu: 0.3113131693571993 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular128_Build",
            "value": 0.31174479402588906,
            "unit": "ns/iter",
            "extra": "iterations: 2250143250\ncpu: 0.31170205941332846 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/2D_Square_GridFunction_Projection_Real_SumOfComponents",
            "value": 593.5257223581716,
            "unit": "ns/iter",
            "extra": "iterations: 1181443\ncpu: 593.4679125442358 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular16_GridFunction_Projection_Real_SumOfComponents",
            "value": 121679.33991341981,
            "unit": "ns/iter",
            "extra": "iterations: 5775\ncpu: 121669.27220779221 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_GridFunction_Projection_Real_SumOfComponents",
            "value": 524780.8456928896,
            "unit": "ns/iter",
            "extra": "iterations: 1335\ncpu: 524730.2846441952 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/2D_Square_GridFunction_Projection_Vector_Components",
            "value": 993.1546018312994,
            "unit": "ns/iter",
            "extra": "iterations: 704513\ncpu: 993.1070044129763 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular16_GridFunction_Projection_Vector_Components",
            "value": 212667.48801170848,
            "unit": "ns/iter",
            "extra": "iterations: 3420\ncpu: 212648.29707602342 ns\nthreads: 1"
          },
          {
            "name": "P1Benchmark/UniformTriangular32_GridFunction_Projection_Vector_Components",
            "value": 899084.8237548118,
            "unit": "ns/iter",
            "extra": "iterations: 783\ncpu: 899053.9157088115 ns\nthreads: 1"
          },
          {
            "name": "Poisson_UniformGrid_16x16/Assembly_NoCoefficient_ConstantSource",
            "value": 190843.18243796032,
            "unit": "ns/iter",
            "extra": "iterations: 3667\ncpu: 190835.48595582205 ns\nthreads: 1"
          },
          {
            "name": "Poisson_UniformGrid_16x16/Assembly_ConstantCoefficient_ConstantSource",
            "value": 190602.35614943618,
            "unit": "ns/iter",
            "extra": "iterations: 3667\ncpu: 190587.90128170128 ns\nthreads: 1"
          },
          {
            "name": "MeshIO/Load_MEDIT_2D_Square",
            "value": 14053.533333332749,
            "unit": "ns/iter",
            "extra": "iterations: 47505\ncpu: 14052.646205662539 ns\nthreads: 1"
          },
          {
            "name": "MeshIO/Load_MEDIT_2D_UniformTriangular64",
            "value": 7019142.200000487,
            "unit": "ns/iter",
            "extra": "iterations: 100\ncpu: 7018550.829999999 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_16x16",
            "value": 143720.3089979542,
            "unit": "ns/iter",
            "extra": "iterations: 4890\ncpu: 143705.20347648268 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_64x64",
            "value": 2617506.9477612376,
            "unit": "ns/iter",
            "extra": "iterations: 268\ncpu: 2617114.138059702 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_128x128",
            "value": 13091959.188679479,
            "unit": "ns/iter",
            "extra": "iterations: 53\ncpu: 13091048.622641478 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_256x256",
            "value": 96825469.28571258,
            "unit": "ns/iter",
            "extra": "iterations: 7\ncpu: 96818434.85714278 ns\nthreads: 1"
          },
          {
            "name": "UniformGrid/Triangular_512x512",
            "value": 345335833.99999654,
            "unit": "ns/iter",
            "extra": "iterations: 2\ncpu: 345325191.0000006 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_AllPairs",
            "value": 1243988.858657416,
            "unit": "ns/iter",
            "extra": "iterations: 566\ncpu: 1243912.5176679017 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_AllPairs",
            "value": 19.656028249076236,
            "unit": "ns/iter",
            "extra": "iterations: 35811153\ncpu: 19.652575916782148 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_0_0",
            "value": 1237577.7367485296,
            "unit": "ns/iter",
            "extra": "iterations: 566\ncpu: 1237174.4876324707 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_0_1",
            "value": 964675.5172410975,
            "unit": "ns/iter",
            "extra": "iterations: 725\ncpu: 964629.9075861891 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_1_0",
            "value": 467062.7943920537,
            "unit": "ns/iter",
            "extra": "iterations: 1498\ncpu: 467007.31108146865 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Cold_Compute_1_1",
            "value": 968428.1533160352,
            "unit": "ns/iter",
            "extra": "iterations: 724\ncpu: 968393.2154697052 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_0_0",
            "value": 4.984491835391008,
            "unit": "ns/iter",
            "extra": "iterations: 140423258\ncpu: 4.9834071575237235 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_0_1",
            "value": 4.9823673302108,
            "unit": "ns/iter",
            "extra": "iterations: 140406758\ncpu: 4.982156428681311 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_1_0",
            "value": 2.4916666590226777,
            "unit": "ns/iter",
            "extra": "iterations: 279087224\ncpu: 2.4915369110554466 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Edge_Warm_Compute_1_1",
            "value": 4.670904851760112,
            "unit": "ns/iter",
            "extra": "iterations: 149636378\ncpu: 4.670318837843019 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_AllPairs",
            "value": 9040737.94871844,
            "unit": "ns/iter",
            "extra": "iterations: 78\ncpu: 9038687.94871826 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_AllPairs",
            "value": 46.566323544906226,
            "unit": "ns/iter",
            "extra": "iterations: 15025667\ncpu: 46.558152393501395 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_2_1",
            "value": 5536101.380951348,
            "unit": "ns/iter",
            "extra": "iterations: 126\ncpu: 5535276.85714257 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_1_2",
            "value": 6246030.874997644,
            "unit": "ns/iter",
            "extra": "iterations: 112\ncpu: 6244973.785713981 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_2_2",
            "value": 2433543.6215275526,
            "unit": "ns/iter",
            "extra": "iterations: 288\ncpu: 2433384.9826393784 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Compute_1_1",
            "value": 7766828.522224185,
            "unit": "ns/iter",
            "extra": "iterations: 90\ncpu: 7766375.68888895 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_2_1",
            "value": 5.290445535318403,
            "unit": "ns/iter",
            "extra": "iterations: 132253825\ncpu: 5.290348040973498 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_1_2",
            "value": 5.2953493832149965,
            "unit": "ns/iter",
            "extra": "iterations: 132399105\ncpu: 5.295075793752541 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_2_2",
            "value": 4.681400171919215,
            "unit": "ns/iter",
            "extra": "iterations: 149894993\ncpu: 4.681173059596424 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Compute_1_1",
            "value": 5.917245023380464,
            "unit": "ns/iter",
            "extra": "iterations: 118272452\ncpu: 5.916958202574519 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Build_1",
            "value": 3853732.5580080366,
            "unit": "ns/iter",
            "extra": "iterations: 181\ncpu: 3853423.464088192 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Build_1",
            "value": 1396262.7219996192,
            "unit": "ns/iter",
            "extra": "iterations: 500\ncpu: 1396184.5040001464 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Transpose_1_2",
            "value": 2307225.4098349093,
            "unit": "ns/iter",
            "extra": "iterations: 305\ncpu: 2307036.1245902823 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Transpose_1_2",
            "value": 557456.1165596466,
            "unit": "ns/iter",
            "extra": "iterations: 1244\ncpu: 557437.1374598484 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Intersection_2_2_via_0",
            "value": 2205088.217664745,
            "unit": "ns/iter",
            "extra": "iterations: 317\ncpu: 2204892.643533252 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Intersection_2_2_via_0",
            "value": 1249037.730357421,
            "unit": "ns/iter",
            "extra": "iterations: 560\ncpu: 1248874.819643004 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Cold_Intersection_1_1_via_0",
            "value": 3541462.5692268577,
            "unit": "ns/iter",
            "extra": "iterations: 195\ncpu: 3541247.2666667677 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Triangle_Warm_Intersection_1_1_via_0",
            "value": 1658244.5704218966,
            "unit": "ns/iter",
            "extra": "iterations: 426\ncpu: 1658127.6478874276 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_AllPairs",
            "value": 5608772.887101515,
            "unit": "ns/iter",
            "extra": "iterations: 124\ncpu: 5608623.524193582 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_AllPairs",
            "value": 46.549680653663614,
            "unit": "ns/iter",
            "extra": "iterations: 15043542\ncpu: 46.547795658761714 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_2_1",
            "value": 3272669.855141324,
            "unit": "ns/iter",
            "extra": "iterations: 214\ncpu: 3272499.785046959 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_1_2",
            "value": 3695948.687829704,
            "unit": "ns/iter",
            "extra": "iterations: 189\ncpu: 3695611.0370370634 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_2_2",
            "value": 1338433.353728162,
            "unit": "ns/iter",
            "extra": "iterations: 523\ncpu: 1338359.1701718546 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Compute_1_1",
            "value": 4546716.954546089,
            "unit": "ns/iter",
            "extra": "iterations: 154\ncpu: 4546355.214284563 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_2_1",
            "value": 5.305567476805629,
            "unit": "ns/iter",
            "extra": "iterations: 132403623\ncpu: 5.3051067718894105 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_1_2",
            "value": 5.300428484803238,
            "unit": "ns/iter",
            "extra": "iterations: 131371987\ncpu: 5.300134327723921 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_2_2",
            "value": 4.671832148867623,
            "unit": "ns/iter",
            "extra": "iterations: 149782725\ncpu: 4.671379900452447 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Compute_1_1",
            "value": 5.970966733168338,
            "unit": "ns/iter",
            "extra": "iterations: 118122429\ncpu: 5.9704227636564955 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Build_1",
            "value": 2488326.8831613865,
            "unit": "ns/iter",
            "extra": "iterations: 291\ncpu: 2488064.3298970945 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Build_1",
            "value": 917150.0839898264,
            "unit": "ns/iter",
            "extra": "iterations: 762\ncpu: 917139.5721785206 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Transpose_1_2",
            "value": 1474406.6890760418,
            "unit": "ns/iter",
            "extra": "iterations: 476\ncpu: 1474251.1008399953 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Transpose_1_2",
            "value": 365615.8224538672,
            "unit": "ns/iter",
            "extra": "iterations: 1915\ncpu: 365703.4778066173 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Intersection_2_2_via_0",
            "value": 1159612.2925622643,
            "unit": "ns/iter",
            "extra": "iterations: 605\ncpu: 1159530.0446278227 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Intersection_2_2_via_0",
            "value": 541071.6821309146,
            "unit": "ns/iter",
            "extra": "iterations: 1315\ncpu: 541112.8798479376 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Cold_Intersection_1_1_via_0",
            "value": 2079282.476191258,
            "unit": "ns/iter",
            "extra": "iterations: 336\ncpu: 2079129.6904767652 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Quadrilateral_Warm_Intersection_1_1_via_0",
            "value": 845557.6516988697,
            "unit": "ns/iter",
            "extra": "iterations: 824\ncpu: 845513.7936892443 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_AllPairs",
            "value": 82599241.3333299,
            "unit": "ns/iter",
            "extra": "iterations: 9\ncpu: 82588875.3333304 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_AllPairs",
            "value": 86.16659081997973,
            "unit": "ns/iter",
            "extra": "iterations: 8128239\ncpu: 86.15746891792945 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_1",
            "value": 21503746.212121397,
            "unit": "ns/iter",
            "extra": "iterations: 33\ncpu: 21500717.39394014 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_2",
            "value": 24275989.413783446,
            "unit": "ns/iter",
            "extra": "iterations: 29\ncpu: 24273547.517241985 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_2_3",
            "value": 26601788.50000359,
            "unit": "ns/iter",
            "extra": "iterations: 26\ncpu: 26599621.999998555 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Compute_3_3",
            "value": 11034501.28571438,
            "unit": "ns/iter",
            "extra": "iterations: 63\ncpu: 11033335.269841937 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_1",
            "value": 5.293116648723811,
            "unit": "ns/iter",
            "extra": "iterations: 131945627\ncpu: 5.292640998249955 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_2",
            "value": 5.290772588854536,
            "unit": "ns/iter",
            "extra": "iterations: 132052358\ncpu: 5.29012627703321 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_2_3",
            "value": 5.29211092121541,
            "unit": "ns/iter",
            "extra": "iterations: 131700396\ncpu: 5.2914668684822255 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Compute_3_3",
            "value": 4.685108026401786,
            "unit": "ns/iter",
            "extra": "iterations: 149899883\ncpu: 4.684387799021842 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Build_1",
            "value": 11439078.803276017,
            "unit": "ns/iter",
            "extra": "iterations: 61\ncpu: 11438305.672129678 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Build_1",
            "value": 6973899.821783822,
            "unit": "ns/iter",
            "extra": "iterations: 101\ncpu: 6973256.643563823 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Build_2",
            "value": 14069344.760000603,
            "unit": "ns/iter",
            "extra": "iterations: 50\ncpu: 14068145.480000284 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Build_2",
            "value": 5688352.081303891,
            "unit": "ns/iter",
            "extra": "iterations: 123\ncpu: 5688154.146342202 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Transpose_2_3",
            "value": 7368590.042559884,
            "unit": "ns/iter",
            "extra": "iterations: 94\ncpu: 7367859.999999806 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Transpose_2_3",
            "value": 1910560.5449573714,
            "unit": "ns/iter",
            "extra": "iterations: 367\ncpu: 1910513.3623980843 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Intersection_3_3_via_0",
            "value": 10567384.38806358,
            "unit": "ns/iter",
            "extra": "iterations: 67\ncpu: 10566344.597014537 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Intersection_3_3_via_0",
            "value": 7798591.306820672,
            "unit": "ns/iter",
            "extra": "iterations: 88\ncpu: 7798465.056816927 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Cold_Intersection_3_3_via_2",
            "value": 7690644.879122313,
            "unit": "ns/iter",
            "extra": "iterations: 91\ncpu: 7689197.329669779 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Tetrahedron_Warm_Intersection_3_3_via_2",
            "value": 1878244.4048227645,
            "unit": "ns/iter",
            "extra": "iterations: 373\ncpu: 1878059.5227881733 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_AllPairs",
            "value": 17950649.461535923,
            "unit": "ns/iter",
            "extra": "iterations: 39\ncpu: 17948953.205127977 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_AllPairs",
            "value": 86.18538851393477,
            "unit": "ns/iter",
            "extra": "iterations: 8127688\ncpu: 86.1753910829256 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_1",
            "value": 5181073.659266784,
            "unit": "ns/iter",
            "extra": "iterations: 135\ncpu: 5180243.281481012 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_2",
            "value": 4538251.806449749,
            "unit": "ns/iter",
            "extra": "iterations: 155\ncpu: 4537528.399999891 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_2_3",
            "value": 5098375.262776204,
            "unit": "ns/iter",
            "extra": "iterations: 137\ncpu: 5097982.9197088815 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Compute_3_3",
            "value": 1576878.6966278192,
            "unit": "ns/iter",
            "extra": "iterations: 445\ncpu: 1576747.928089457 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_1",
            "value": 5.293445982920458,
            "unit": "ns/iter",
            "extra": "iterations: 132288456\ncpu: 5.293038645790855 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_2",
            "value": 5.296960762314508,
            "unit": "ns/iter",
            "extra": "iterations: 132374576\ncpu: 5.296540681648786 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_2_3",
            "value": 5.301492583413729,
            "unit": "ns/iter",
            "extra": "iterations: 132416720\ncpu: 5.301397104534804 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Compute_3_3",
            "value": 4.670186991806217,
            "unit": "ns/iter",
            "extra": "iterations: 149780199\ncpu: 4.669881744515512 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Build_1",
            "value": 4030146.3620698215,
            "unit": "ns/iter",
            "extra": "iterations: 174\ncpu: 4029650.4770110943 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Build_1",
            "value": 2217212.9056618973,
            "unit": "ns/iter",
            "extra": "iterations: 318\ncpu: 2217051.60062869 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Build_2",
            "value": 3399269.2815583446,
            "unit": "ns/iter",
            "extra": "iterations: 206\ncpu: 3398639.820388176 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Build_2",
            "value": 1387447.0714279863,
            "unit": "ns/iter",
            "extra": "iterations: 504\ncpu: 1387295.938492007 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Transpose_2_3",
            "value": 1740839.1374993927,
            "unit": "ns/iter",
            "extra": "iterations: 400\ncpu: 1740707.4825005964 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Transpose_2_3",
            "value": 521838.81799976464,
            "unit": "ns/iter",
            "extra": "iterations: 1000\ncpu: 521913.15899946744 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Intersection_3_3_via_0",
            "value": 1363571.3235896165,
            "unit": "ns/iter",
            "extra": "iterations: 513\ncpu: 1363477.9298250591 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Intersection_3_3_via_0",
            "value": 796863.8361956854,
            "unit": "ns/iter",
            "extra": "iterations: 873\ncpu: 796800.9633439471 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Cold_Intersection_3_3_via_2",
            "value": 1822855.6067733474,
            "unit": "ns/iter",
            "extra": "iterations: 384\ncpu: 1822689.3229168856 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Hexahedron_Warm_Intersection_3_3_via_2",
            "value": 470014.80053563346,
            "unit": "ns/iter",
            "extra": "iterations: 1494\ncpu: 470125.99464525585 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_AllPairs",
            "value": 30025999.608695835,
            "unit": "ns/iter",
            "extra": "iterations: 23\ncpu: 30024075.391306426 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_AllPairs",
            "value": 86.17782103677902,
            "unit": "ns/iter",
            "extra": "iterations: 8137035\ncpu: 86.1751216505773 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_1",
            "value": 8060418.64368241,
            "unit": "ns/iter",
            "extra": "iterations: 87\ncpu: 8059985.137932777 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_2",
            "value": 8439916.51807695,
            "unit": "ns/iter",
            "extra": "iterations: 83\ncpu: 8439610.638555313 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_2_3",
            "value": 9402097.932428893,
            "unit": "ns/iter",
            "extra": "iterations: 74\ncpu: 9401589.891890533 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Compute_3_3",
            "value": 2989549.4957222035,
            "unit": "ns/iter",
            "extra": "iterations: 234\ncpu: 2989410.2777779233 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_1",
            "value": 5.289295598771652,
            "unit": "ns/iter",
            "extra": "iterations: 132063015\ncpu: 5.28882289261681 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_2",
            "value": 5.29045693668857,
            "unit": "ns/iter",
            "extra": "iterations: 132107339\ncpu: 5.290254018363104 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_2_3",
            "value": 5.308911408215555,
            "unit": "ns/iter",
            "extra": "iterations: 132254002\ncpu: 5.308633843836498 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Compute_3_3",
            "value": 4.666961641482893,
            "unit": "ns/iter",
            "extra": "iterations: 149945004\ncpu: 4.666758620380594 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Build_1",
            "value": 5745360.942153992,
            "unit": "ns/iter",
            "extra": "iterations: 121\ncpu: 5745087.008265291 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Build_1",
            "value": 3585233.6769296066,
            "unit": "ns/iter",
            "extra": "iterations: 195\ncpu: 3584994.0461530834 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Build_2",
            "value": 6145444.06140208,
            "unit": "ns/iter",
            "extra": "iterations: 114\ncpu: 6145299.552632341 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Build_2",
            "value": 2383735.1501648524,
            "unit": "ns/iter",
            "extra": "iterations: 293\ncpu: 2383604.5597270094 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Transpose_2_3",
            "value": 2920858.075313088,
            "unit": "ns/iter",
            "extra": "iterations: 239\ncpu: 2920694.213388415 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Transpose_2_3",
            "value": 798614.5273961271,
            "unit": "ns/iter",
            "extra": "iterations: 876\ncpu: 798580.2968039622 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Intersection_3_3_via_0",
            "value": 2697180.7490358287,
            "unit": "ns/iter",
            "extra": "iterations: 259\ncpu: 2697080.664092019 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Intersection_3_3_via_0",
            "value": 1821652.477924958,
            "unit": "ns/iter",
            "extra": "iterations: 385\ncpu: 1821586.9038974156 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Cold_Intersection_3_3_via_2",
            "value": 3138933.1696440238,
            "unit": "ns/iter",
            "extra": "iterations: 224\ncpu: 3138620.6026778095 ns\nthreads: 1"
          },
          {
            "name": "ConnectivityBenchmark/Wedge_Warm_Intersection_3_3_via_2",
            "value": 829749.6757069319,
            "unit": "ns/iter",
            "extra": "iterations: 848\ncpu: 829762.3113212814 ns\nthreads: 1"
          }
        ]
      }
    ]
  }
}