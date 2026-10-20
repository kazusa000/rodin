/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "KelvinBallOptimization.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/Adaptation.h>
#include <Rodin/Alert/Info.h>
#include <Rodin/Location.h>
#include <Rodin/Alert/Notation.h>
#include <Rodin/Alert/Raise.h>
#include <Rodin/Alert/Success.h>
#include <Rodin/Alert/Warning.h>
#include <Rodin/Assembly.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/IO/HDF5.h>
#include <Rodin/MMG.h>
#include <Rodin/Variational.h>

#include "Configuration.h"
#include "Metrics.h"
#include "RotatedCharacteristicContinuation.h"
#include "RotatedNitscheIntegrator.h"
#include "SewedOutput.h"
#include "Thickness.h"
#include "BaselineThickness.h"
#include "Sphere.h"
#include "../../WNGIRExampleParameters.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace KelvinBall
{
  class KelvinBallOptimization::Implementation
  {
    public:
      Implementation(int argc, char** argv)
        : m_argc(argc),
          m_argv(argv)
      {}

      int run();

    private:
      // Nonempty only in the saved-field diagnostic entry; normal runs do not
      // write intermediate meshes or change their reconstruction parameters.
      std::string m_mmgReplayPrefix;
      Real m_mmgReplayOuterRadius = 2;
      bool m_mmgReplayLockCuts = false;

      static constexpr Attribute Obstacle = KelvinBall::Obstacle;
      static constexpr Attribute Fluid = KelvinBall::Fluid;
      static constexpr Attribute Outer = KelvinBall::Outer;
      static constexpr Attribute Gamma = KelvinBall::Gamma;
      static constexpr Attribute SigmaPlus = KelvinBall::SigmaPlus;
      static constexpr Attribute SigmaMinus = KelvinBall::SigmaMinus;
      static constexpr Attribute SigmaXYPlus = KelvinBall::SigmaXYPlus;
      static constexpr Attribute SigmaXYMinus = KelvinBall::SigmaXYMinus;

      static constexpr Real mu = KelvinBall::Mu;
      static constexpr Real remeshGradation = 2.0;
      static constexpr Real rotatedTracePhysicalTolerance = 0.01;
      static constexpr Real rotatedTraceReferenceTolerance = 0.25;
      static constexpr Real linearResidualTolerance = KelvinBall::LinearResidualTolerance;
      static constexpr size_t chamberMultiplicity = KelvinBall::ChamberMultiplicity;
#ifdef RODIN_USE_OPENMP
      static constexpr const char* assemblyBackend = "OpenMP";
#else
      static constexpr const char* assemblyBackend = "Sequential";
#endif

      Alert::Text<Alert::CyanT> stageHeading(const std::string& text)
      {
        Alert::Text<Alert::CyanT> heading(Alert::Cyan, text);
        return heading.setBold();
      }

      /// Formats a vector as (x, y, z), with signed zeros printed as 0.
      std::string vectorText(const Math::SpatialVector<Real>& v)
      {
        std::ostringstream text;
        text << '(';
        for (Eigen::Index i = 0; i < v.size(); ++i)
          text << (i ? ", " : "") << v(i) + Real(0);
        text << ')';
        return text.str();
      }

      Alert::Text<Alert::YellowT> substageHeading(const std::string& text)
      {
        Alert::Text<Alert::YellowT> heading(Alert::Yellow, text);
        return heading.setBold();
      }

      std::string diagnosticLabel(const std::string& text)
      {
        constexpr size_t width = 36;
        return text + std::string(text.size() < width ? width - text.size() : 1, ' ');
      }

      struct GradientDiagnostics
      {
          Real residual;
          Real jump;
      };

      struct MeshDiagnostics
      {
          size_t vertices = 0;
          size_t cells = 0;
          size_t obstacleCells = 0;
          size_t fluidCells = 0;
          size_t interfaceTriangles = 0;
          Real minimumQuality = std::numeric_limits<Real>::infinity();
          Real meanQuality = 0;
          Real maximumQuality = 0;
          Real meanElementSize = 0;
      };

      using ReconstructionDiagnostics = KelvinBall::ReconstructionDiagnostics;

      struct MMGReconstruction
      {
          MMG::Mesh mesh;
          ReconstructionDiagnostics diagnostics;
      };

      template <class Displacement>
      void moveMesh(KelvinBall::Mesh& moved, const KelvinBall::Mesh& mesh,
        const Displacement& displacement)
      {
        const auto& space = displacement.getFiniteElementSpace();
        const auto& coefficients = displacement.getData();
        for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
        {
          auto coordinates = mesh.getVertexCoordinates(vertex);
          const auto& dofs = space.getDOFs(0, vertex);
          for (size_t component = 0; component < mesh.getSpaceDimension(); ++component)
            coordinates(component) += coefficients(dofs[component]);
          moved.setVertexCoordinates(vertex, coordinates);
        }
      }

      template <class LevelSet>
      MMG::Mesh classifyLevelSetForWNGIR(
        const MMG::Mesh& background, const LevelSet& levelSet)
      {
        MMG::Mesh classified(background);
        const auto& space = levelSet.getFiniteElementSpace();
        const auto& coefficients = levelSet.getData();
        for (auto cell = classified.getCell(); cell; ++cell)
        {
          Real meanLevelSet = 0;
          for (const Index vertex : cell->getVertices())
          {
            const auto& dofs = space.getDOFs(0, vertex);
            meanLevelSet += coefficients(dofs[0]);
          }
          meanLevelSet /= static_cast<Real>(cell->getVertices().size());
          classified.setAttribute({classified.getDimension(), cell->getIndex()},
            meanLevelSet < 0 ? Obstacle : Fluid);
        }

        classified.getConnectivity().compute(2, 3);
        for (auto face = classified.getFace(); face; ++face)
        {
          const auto& incident =
            classified.getConnectivity().getIncidence({2, 3}, face->getIndex());
          if (incident.size() != 2)
            continue;
          const auto first = classified.getCell(incident[0])->getAttribute();
          const auto second = classified.getCell(incident[1])->getAttribute();
          if (first && second && *first != *second)
            classified.setAttribute({2, face->getIndex()}, Gamma);
        }
        return classified;
      }

      template <class LevelSet>
      MMGReconstruction fitLevelSetWNGIR(const KelvinBall::Mesh& mesh,
        const LevelSet& levelSet, Real h, Real outerRadius, int argc, char** argv)
      {
        P1<Math::SpatialVector<Real>, KelvinBall::Mesh> gradientSpace(mesh, 3);
        GridFunction projectedGradient(gradientSpace);
        TrialFunction gradientTrial(gradientSpace);
        TestFunction gradientTest(gradientSpace);
        Problem gradientProjection(gradientTrial, gradientTest);
        gradientProjection =
          Integral(gradientTrial, gradientTest) - Integral(Grad(levelSet), gradientTest);
        gradientProjection.assemble();
        Solver::CG(gradientProjection).solve();
        projectedGradient.getData() = gradientTrial.getSolution().getData();

        P1<Math::SpatialVector<Real>, KelvinBall::Mesh> displacementSpace(mesh, 3);
        TrialFunction displacementTrial(displacementSpace);
        TestFunction displacementTest(displacementSpace);
        Rodin::Examples::WNGIRExampleDefaults defaults;
        defaults.kappaBulk = Real(8e-4);
        // A fit from the classified staircase converges in seven to nine steps.
        defaults.maxIterations = 12;
        auto parameters =
          Rodin::Examples::makeWNGIRParameters(argc, argv, h, Gamma, defaults);
        if (!Rodin::Examples::findOption(
              argc, argv, "wngir-primal-barrier-iterations", nullptr))
          parameters.primalBarrierIterations = 30;
        if (!Rodin::Examples::findOption(argc, argv, "wngir-cg-max-iters", nullptr))
          parameters.cgMaxIterations = 10000;
        // The outer sphere is curved, so its nodes are pinned: sliding in a
        // facet plane would walk them off the sphere. The cut planes bound the
        // wedge but are not walls, and the rim of the interface lies entirely
        // on them, so they slide within themselves instead of being frozen.
        parameters.fixedBoundaryAttributes = {Outer};
        parameters.slipBoundaryAttributes = {
          SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus};
        // Pinning the outer sphere removes every rigid mode, so no floor is
        // needed; the option still overrides this default.
        if (!Rodin::Examples::findOption(
              argc, argv, "wngir-rigid-stabilisation", nullptr))
          parameters.rigidStabilisationLevel = 0;
        Adaptation::WNGIR fitting(displacementTrial, displacementTest);
        fitting.setParameters(parameters);

        std::vector<Index> interface;
        for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
          if (face->getAttribute() == Attribute{Gamma})
            interface.push_back(face->getIndex());

        RealFunction target(
          [&](const Geometry::Point& point) { return levelSet.getValue(point); });
        Adaptation::AnalyticVectorFunction targetGradient(
          [&](const Geometry::Point& point) { return projectedGradient.getValue(point); },
          3);

        const auto report = fitting.solve(mesh, interface, target, targetGradient);
        Alert::Info() << substageHeading("WNGIR reconstruction") << Alert::NewLine
                      << diagnosticLabel("Iterations:")
                      << Alert::Notation::Number(report.iterations) << Alert::NewLine

                      << diagnosticLabel("Skeleton normal jump RMS:")
                      << Alert::Notation::Number(report.normalJumpRMS) << Alert::NewLine
                      << diagnosticLabel("Exit reason:") << report.exitReason
                      << Alert::NewLine << diagnosticLabel("Active RMS:")
                      << Alert::Notation::Number(report.activeRMS) << Alert::NewLine
                      << diagnosticLabel("Active RMS / level-set mesh scale:")
                      << Alert::Notation::Number(report.levelSetGradientScale > 0
                             ? report.activeRMS /
                               (parameters.h * report.levelSetGradientScale)
                             : 0)
                      << Alert::NewLine << diagnosticLabel("RMS stopping tolerance:")
                      << Alert::Notation::Number(report.effectiveTauRms) << Alert::NewLine
                      << diagnosticLabel("Scaled RMS stopping tolerance:")
                      << Alert::Notation::Number(report.effectiveTauRmsH)
                      << Alert::NewLine << diagnosticLabel("Active supremum:")
                      << Alert::Notation::Number(report.activeSup) << Alert::NewLine
                      << diagnosticLabel("Supremum stopping tolerance:")
                      << Alert::Notation::Number(report.effectiveTauInf) << Alert::NewLine
                      << diagnosticLabel("Scaled supremum stopping tolerance:")
                      << Alert::Notation::Number(report.effectiveTauInfH)
                      << Alert::NewLine << diagnosticLabel("Minimum Jacobian:")
                      << Alert::Notation::Number(report.minJ) << Alert::NewLine
                      << diagnosticLabel("Maximum relative distortion:")
                      << Alert::Notation::Number(report.maxQRel) << Alert::Raise;

        KelvinBall::Mesh moved(mesh);
        moveMesh(moved, mesh, displacementTrial.getSolution());
        checkFixedGeometry(moved, outerRadius);
        checkMaterials(moved);
        const MeshDiagnostics diagnostics = getMeshDiagnostics(moved);
        return {MMG::Mesh(std::move(moved)),
          {h, h, 0, parameters.fixedBoundaryAttributes.size(), diagnostics.cells,
            diagnostics.cells}};
      }

      void printUsage(const char* executable)
      {
        Alert::Info()
          << "Usage" << Alert::NewLine << "  " << executable << " [options]"
          << Alert::NewLine << Alert::Notation("--n=<points>")
          << "              Background points per edge (default: 13)." << Alert::NewLine
          << Alert::Notation("--h=<size>")
          << "                Requested background size; alternative to --n."
          << Alert::NewLine << Alert::Notation("--outer-radius=<value>")
          << "     Chamber outer radius (default: 2)." << Alert::NewLine
          << Alert::Notation("--iterations=<count>")
          << "       Number of evaluated designs (default: 1)." << Alert::NewLine
          << Alert::Notation("--penalty=<value>")
          << "          Rotational Nitsche penalty (default: 320)." << Alert::NewLine
          << Alert::Notation("--stabilization=<value>")
          << "    P1--P1 pressure-stabilization factor (default: 0.05)." << Alert::NewLine
          << Alert::Notation("--regularization=<value>")
          << "   H1 smoothing length in multiples of h (default: 4)." << Alert::NewLine
          << Alert::Notation("--step=<value>")
          << "             Advection time step in multiples of h (default: 0.1)."
          << Alert::NewLine << Alert::Notation("--level-set-penalty=<value>")
          << " Rotated trace penalty of the level set (default: 1)." << Alert::NewLine
          << Alert::Notation("--thickness-min=<value>")
          << "      Minimum body thickness in h (default: 2; 0 disables it)."
          << Alert::NewLine << Alert::Notation("--thickness-weight=<value>")
          << "   Weight of the thickness penalty (default: 1)." << Alert::NewLine
          << Alert::Notation("--thickness-method=<legacy|baseline>")
          << " Thickness strategy (default: legacy; baseline is Cubic-only)."
          << Alert::NewLine << Alert::Notation("--restart-mesh=<path>")
          << " Read an existing MEDIT volume mesh without initial reconstruction."
          << Alert::NewLine << Alert::Notation("--restart-snapshot=<path>")
          << " Original HDF5 supplying exact double coordinates."
          << Alert::NewLine << Alert::Notation("--restart-target-volume=<value>")
          << " Original chamber target volume; required for restart."
          << Alert::NewLine << Alert::Notation("--initial-global-step=<count>")
          << " Global CSV offset for a restart (default: 0)." << Alert::NewLine
          << Alert::Notation("--motion-every=<count>")
          << "      Write the rigid motion every count iterates (default: 0, off)."
          << Alert::NewLine << Alert::Notation("--motion-force=<fx,fy,fz>")
          << "  Force driving the motion (default: 0,0,1)." << Alert::NewLine
          << Alert::Notation("--motion-frames=<count>")
          << "     Frames over one revolution (default: 24)." << Alert::NewLine
          << Alert::Notation("--advection-quadrature=<order>")
          << " Quadrature order for the transported distance (default: 8)."
          << Alert::NewLine << Alert::Notation("--reconstruction=<method>")
          << "  Interface reconstruction: mmg or wngir (default: mmg)." << Alert::NewLine
          << Alert::Notation("--background-hmin=<value>")
          << "  WNGIR background minimum size, in h (default: 0.1)." << Alert::NewLine
          << Alert::Notation("--background-hmax=<value>")
          << "  WNGIR background maximum size, in h (default: 1)." << Alert::NewLine
          << Alert::Notation("--background-hausdorff=<value>")
          << " WNGIR background Hausdorff tolerance, in h (default: 0.05)."
          << Alert::NewLine << Alert::Notation("--background-gradation=<value>")
          << " WNGIR background gradation (default: 2)." << Alert::NewLine
          << Alert::Notation("--mmg-adapt")
          << "                MMG path: replace the optimization pass after each"
          << Alert::NewLine
          << "                              cut by adaptation to a size map."
          << Alert::NewLine << Alert::Notation("--mmg-adapt-interface-size=<value>")
          << " Size on Gamma, in h (default: 1)." << Alert::NewLine
          << Alert::Notation("--mmg-adapt-far-size=<value>")
          << "   Size away from Gamma, in h (default: 1)." << Alert::NewLine
          << Alert::Notation("--mmg-adapt-width=<value>")
          << "      Distance over which the size grows, in h (default: 3)."
          << Alert::NewLine << Alert::Notation("--mmg-adapt-gradation=<value>")
          << "  Adaptation gradation (default: 1.3)." << Alert::NewLine
          << Alert::Notation("--mmg-snap=<value>")
          << "         MMG path: snap edge crossings closer than this fraction"
          << Alert::NewLine
          << "                              of the edge to a vertex (default: 0, off)."
          << Alert::NewLine << Alert::Notation("--mmg-retries=<count>")
          << "      MMG path: retries of a failed reconstruction, each at half"
          << Alert::NewLine
          << "                              the previous MMG scale (default: 2)."
          << Alert::NewLine << Alert::Notation("--wngir-*=<value>")
          << "        WNGIR fitting parameters (--wngir-steps defaults to 12;"
          << Alert::NewLine
          << "                              also --trace, --j-safe, --j-ls, --j-min)."
          << Alert::NewLine << Alert::Notation("--geometry-only")
          << "             Stop after initial MMG reconstruction." << Alert::NewLine
          << Alert::Notation("--state-only")
          << "                Stop after the first Stokes evaluation." << Alert::NewLine
          << Alert::Notation("--replay-mmg=<directory>")
          << "   Diagnose saved advected fields without optimization." << Alert::NewLine
          << Alert::Notation("--replay-mmg-count=<count>")
          << " Replay one to three saved fields (default: 3)." << Alert::NewLine
          << Alert::Notation("--replay-mmg-lock-cuts")
          << "     Diagnostic control: require cut faces during adaptation." << Alert::NewLine
          << Alert::Notation("--save-mesh")
          << "                 Save KelvinBallInitial.mesh." << Alert::NewLine
          << Alert::Notation("--help") << "                      Show this message."
          << Alert::Raise;
      }

      void announce(const char* stage)
      {
        Alert::Info() << stageHeading(stage) << Alert::Raise;
      }

      using Clock = std::chrono::steady_clock;

      Real elapsedSeconds(const Clock::time_point& start) const
      {
        return std::chrono::duration<Real>(Clock::now() - start).count();
      }

      void reportStageTiming(size_t stage, Real seconds)
      {
        Alert::Info() << substageHeading("Stage " + std::to_string(stage) + " timing")
                      << Alert::NewLine << diagnosticLabel("Elapsed time:")
                      << Alert::Notation::Number(seconds) << " s" << Alert::Raise;
      }

      using LocalMesh = Geometry::Mesh<Context::Local>;
      using VelocitySpace = P1<Math::SpatialVector<Real>, LocalMesh>;
      using PressureSpace = P1<Real, LocalMesh>;

      VelocitySpace makeVelocitySpace(const LocalMesh& mesh)
      {
        return VelocitySpace(mesh, 3);
      }

      Real planeResidual(
        Attribute attribute, const Math::SpatialPoint& x, Real outerRadius)
      {
        if (KelvinBall::Tetrahedral)
        {
          switch (attribute)
          {
            case Outer: return std::min(std::abs(x(0) - outerRadius),
              std::abs(x(1) - outerRadius));
            case SigmaPlus: return std::abs(x(1) - x(2));
            case SigmaMinus: return std::abs(x(0) - x(2));
            case SigmaXYPlus: return std::abs(x(0) + x(2));
            case SigmaXYMinus: return std::abs(x(1) + x(2));
            default: return 0;
          }
        }
        switch (attribute)
        {
          case Outer:
            return std::abs(x(0) - outerRadius);
          case SigmaPlus:
            return std::abs(x(1) - x(2));
          case SigmaMinus:
            return std::abs(x(1) + x(2));
          case SigmaXYPlus:
          case SigmaXYMinus:
            return std::abs(x(0) - x(1));
          default:
            return 0;
        }
      }

      void recordMMGStage(const MMG::Mesh& mesh, const char* stage, Real outerRadius)
      {
        if (m_mmgReplayPrefix.empty())
          return;
        const std::string name = m_mmgReplayPrefix + "-" + stage;
        mesh.save(name + ".mesh", IO::FileFormat::MEDIT);
        IO::XDMF output(name);
        auto grid = output.grid("Chamber");
        grid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
        output.write(0).flush();
        std::ofstream csv(name + "-boundary.csv");
        csv.precision(std::numeric_limits<Real>::max_digits10);
        csv << "attribute,faces,bad_faces,max_plane_residual,worst_vertex,x,y,z\n";
        for (const Attribute attribute : {Outer, SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus})
        {
          size_t count = 0, bad = 0;
          Real maximum = 0;
          Index worst = 0;
          for (auto face = mesh.getPolytope(2); face; ++face)
          {
            if (face->getAttribute() != attribute)
              continue;
            if (count == 0)
              worst = face->getVertices()[0];
            ++count;
            Real faceResidual = 0;
            for (const Index vertex : face->getVertices())
            {
              const Real residual = planeResidual(attribute, mesh.getVertexCoordinates(vertex), outerRadius);
              faceResidual = std::max(faceResidual, residual);
              if (residual > maximum)
              {
                maximum = residual;
                worst = vertex;
              }
            }
            bad += faceResidual > 1e-6;
          }
          const auto& x = mesh.getVertexCoordinates(worst);
          csv << attribute << ',' << count << ',' << bad << ',' << maximum << ','
              << worst << ',' << x(0) << ',' << x(1) << ',' << x(2) << '\n';
        }
      }

      void checkFixedGeometry(const LocalMesh& mesh, Real outerRadius)
      {
        const FlatSet<Attribute> fixed{
          Outer, SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus};
        std::map<Attribute, size_t> counts;
        Real residual = 0;
        for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
        {
          const auto attribute = face->getAttribute();
          if (!attribute || !fixed.contains(*attribute))
            continue;
          ++counts[*attribute];
          for (const Index vertex : face->getVertices())
          {
            const auto& x = mesh.getVertexCoordinates(vertex);
            residual = std::max(residual, planeResidual(*attribute, x, outerRadius));
          }
        }
        for (const Attribute attribute : fixed)
        {
          if (counts[attribute] == 0)
            throw std::runtime_error("A fixed chamber boundary has no triangles.");
        }
        Alert::Info() << substageHeading("Mesh validation") << Alert::NewLine
                      << diagnosticLabel("Vertices:")
                      << Alert::Notation::Number(mesh.getVertexCount()) << Alert::NewLine
                      << diagnosticLabel("Cells:")
                      << Alert::Notation::Number(mesh.getCellCount()) << Alert::NewLine
                      << diagnosticLabel("Fixed-boundary planarity residual:")
                      << Alert::Notation::Number(residual) << Alert::Raise;
        if (residual > 1e-6)
          throw std::runtime_error(
            "Mesh reconstruction changed the fixed boundary geometry.");
      }

      MeshDiagnostics getMeshDiagnostics(
        const LocalMesh& mesh, bool requireCompletePartition = true)
      {
        MeshDiagnostics diagnostics;
        diagnostics.vertices = mesh.getVertexCount();
        diagnostics.cells = mesh.getCellCount();
        Real qualitySum = 0;
        Real elementSizeSum = 0;
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          if (cell->getAttribute() == Obstacle)
            ++diagnostics.obstacleCells;
          else if (cell->getAttribute() == Fluid)
            ++diagnostics.fluidCells;
          else
            throw std::runtime_error("A cell has an unexpected material attribute.");

          const auto& vertices = cell->getVertices();
          if (vertices.size() != 4)
            throw std::runtime_error("Kelvin-ball mesh quality requires tetrahedra.");
          Real squaredEdgeLengthSum = 0;
          Real edgeLengthSum = 0;
          for (size_t i = 0; i < vertices.size(); ++i)
          {
            for (size_t j = i + 1; j < vertices.size(); ++j)
            {
              const Real edgeLength = (mesh.getVertexCoordinates(vertices[i]) -
                mesh.getVertexCoordinates(vertices[j]))
                                        .norm();
              edgeLengthSum += edgeLength;
              squaredEdgeLengthSum += edgeLength * edgeLength;
            }
          }
          const Real quality =
            12 * std::pow(3 * cell->getMeasure(), Real(2) / 3) / squaredEdgeLengthSum;
          diagnostics.minimumQuality = std::min(diagnostics.minimumQuality, quality);
          diagnostics.maximumQuality = std::max(diagnostics.maximumQuality, quality);
          qualitySum += quality;
          elementSizeSum += edgeLengthSum / 6;
        }
        diagnostics.meanQuality = qualitySum / static_cast<Real>(diagnostics.cells);
        diagnostics.meanElementSize =
          elementSizeSum / static_cast<Real>(diagnostics.cells);
        for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
        {
          if (face->getAttribute() == Gamma)
            ++diagnostics.interfaceTriangles;
        }
        if (requireCompletePartition &&
          (diagnostics.obstacleCells == 0 || diagnostics.fluidCells == 0 ||
            diagnostics.interfaceTriangles == 0))
          throw std::runtime_error("The fitted body/fluid partition is incomplete.");
        return diagnostics;
      }

      void checkMaterials(const LocalMesh& mesh)
      {
        const auto diagnostics = getMeshDiagnostics(mesh);
        Alert::Info() << substageHeading("Material partition") << Alert::NewLine
                      << diagnosticLabel("Obstacle cells:")
                      << Alert::Notation::Number(diagnostics.obstacleCells)
                      << Alert::NewLine << diagnosticLabel("Fluid cells:")
                      << Alert::Notation::Number(diagnostics.fluidCells) << Alert::NewLine
                      << diagnosticLabel("Interface triangles:")
                      << Alert::Notation::Number(diagnostics.interfaceTriangles)
                      << Alert::NewLine << diagnosticLabel("Minimum tetrahedron quality:")
                      << Alert::Notation::Number(diagnostics.minimumQuality)
                      << Alert::NewLine << diagnosticLabel("Mean tetrahedron quality:")
                      << Alert::Notation::Number(diagnostics.meanQuality)
                      << Alert::NewLine << diagnosticLabel("Maximum tetrahedron quality:")
                      << Alert::Notation::Number(diagnostics.maximumQuality)
                      << Alert::NewLine << diagnosticLabel("Mean element size:")
                      << Alert::Notation::Number(diagnostics.meanElementSize)
                      << Alert::Raise;
      }

      void reportSphereGeometry(const LocalMesh& mesh)
      {
        std::set<Index> vertices;
        for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
          if (face->getAttribute() == Gamma)
            vertices.insert(face->getVertices().begin(), face->getVertices().end());
        Real squared = 0;
        Real supremum = 0;
        for (const Index vertex : vertices)
        {
          const auto& x = mesh.getVertexCoordinates(vertex);
          const Real error = std::abs(x.norm() - KelvinBall::initialRadius(x));
          squared += error * error;
          supremum = std::max(supremum, error);
        }
        const Real rms = vertices.empty()
          ? std::numeric_limits<Real>::quiet_NaN()
          : std::sqrt(squared / static_cast<Real>(vertices.size()));
        Alert::Info() << substageHeading("Sphere geometry") << Alert::NewLine
                      << diagnosticLabel("Chamber volume:")
                      << Alert::Notation::Number(mesh.getVolume(Obstacle))
                      << Alert::NewLine << diagnosticLabel("Vertex RMS error:")
                      << Alert::Notation::Number(rms) << Alert::NewLine
                      << diagnosticLabel("Vertex maximum error:")
                      << Alert::Notation::Number(supremum) << Alert::Raise;
      }

      Math::SpatialVector<Real> boundingBoxSize(const LocalMesh& mesh)
      {
        Math::SpatialPoint minimum(3);
        Math::SpatialPoint maximum(3);
        minimum.setConstant(std::numeric_limits<Real>::infinity());
        maximum.setConstant(-std::numeric_limits<Real>::infinity());
        for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
        {
          const auto& point = mesh.getVertexCoordinates(vertex);
          for (size_t component = 0; component < 3; ++component)
          {
            minimum(component) = std::min(minimum(component), point(component));
            maximum(component) = std::max(maximum(component), point(component));
          }
        }
        return maximum - minimum;
      }

      Math::SpatialVector<Real> sewnBoundingBoxSize(const LocalMesh& mesh)
      {
        Math::SpatialPoint minimum(3);
        Math::SpatialPoint maximum(3);
        minimum.setConstant(std::numeric_limits<Real>::infinity());
        maximum.setConstant(-std::numeric_limits<Real>::infinity());
        const auto& rotations = KelvinBall::SewedOutput::getCubeRotations();
        for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
        {
          for (const auto& rotation : rotations)
          {
            const Math::SpatialPoint point = rotation * mesh.getVertexCoordinates(vertex);
            for (size_t component = 0; component < 3; ++component)
            {
              minimum(component) = std::min(minimum(component), point(component));
              maximum(component) = std::max(maximum(component), point(component));
            }
          }
        }
        return maximum - minimum;
      }

      template <class LevelSet>
      MMGReconstruction discretizeLevelSetMMG(MMG::Mesh& mesh, const LevelSet& levelSet,
        Real h, const Sphere& sphere, bool adapt, Real snap)
      {
        const size_t previousCells = mesh.getCellCount();
        const Real hmin = 0.1 * h;
        const Real hmax = 10 * h;
        const Real hausdorff = 0.1 * h * h;
        // A retried reconstruction receives the mesh the failed attempt already
        // reduced to a single material, so the input need not be partitioned.
        const MeshDiagnostics inputDiagnostics = getMeshDiagnostics(mesh, false);

        Alert::Info() << substageHeading("MMG input") << Alert::NewLine
                      << diagnosticLabel("Level-set minimum:")
                      << Alert::Notation::Number(levelSet.min()) << Alert::NewLine
                      << diagnosticLabel("Level-set maximum:")
                      << Alert::Notation::Number(levelSet.max()) << Alert::NewLine
                      << diagnosticLabel("Vertices:")
                      << Alert::Notation::Number(inputDiagnostics.vertices)
                      << Alert::NewLine << diagnosticLabel("Cells:")
                      << Alert::Notation::Number(inputDiagnostics.cells) << Alert::NewLine
                      << diagnosticLabel("Interface triangles:")
                      << Alert::Notation::Number(inputDiagnostics.interfaceTriangles)
                      << Alert::NewLine << diagnosticLabel("Minimum tetrahedron quality:")
                      << Alert::Notation::Number(inputDiagnostics.minimumQuality)
                      << Alert::NewLine << diagnosticLabel("Mean tetrahedron quality:")
                      << Alert::Notation::Number(inputDiagnostics.meanQuality)
                      << Alert::NewLine << diagnosticLabel("Mean element size:")
                      << Alert::Notation::Number(inputDiagnostics.meanElementSize)
                      << Alert::Raise;

        // The level set alone defines the new partition. MMG keeps the faces
        // between differently labelled cells as internal boundaries through
        // the cut, where they end up inside one material; the optimization
        // pass (opnbdy) then preserves them, so every past interface would
        // accumulate in the mesh. The cut therefore starts, as the initial
        // one does, from a single material and only the fixed boundary.
        const FlatSet<Attribute> fixed{
          Outer, SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus};
        for (auto cell = mesh.getCell(); cell; ++cell)
          mesh.setAttribute({mesh.getDimension(), cell->getIndex()}, Fluid);
        for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
        {
          const auto attribute = face->getAttribute();
          if (attribute && !fixed.contains(*attribute))
            mesh.setAttribute({mesh.getDimension() - 1, face->getIndex()}, {});
        }
        const size_t requiredTriangles = sphere.protectFixedGeometry(mesh, false);

        MMG::LevelSetDiscretizer discretizer;
        discretizer.split(Fluid, {Obstacle, Fluid})
          .setHMin(hmin)
          .setHMax(hmax)
          .setHausdorff(hausdorff)
          .setGradation(remeshGradation)
          .setBaseReferences(FlatSet<Attribute>{Fluid})
          .setBoundaryReference(Gamma)
          .setAngleDetection(false);
        // A crossed edge (i, j) is cut at t = phi_i / (phi_i - phi_j). A cut
        // with t near 0 or 1 puts the new vertex next to an existing one and
        // leaves a sliver. Snapping the near endpoint to zero instead makes the
        // interface pass through that vertex, so every remaining cut satisfies
        // snap <= t <= 1 - snap. For a distance-like level set the interface
        // moves by at most snap times the edge length.
        std::decay_t<LevelSet> sanitized(levelSet.getFiniteElementSpace());
        sanitized.getData() = levelSet.getData();
        const auto crossingFraction = [&](Index i, Index j) {
          const Real a = sanitized[i], b = sanitized[j];
          return (a < 0) != (b < 0) && a != 0 && b != 0
            ? a / (a - b)
            : std::numeric_limits<Real>::quiet_NaN();
        };
        const auto scanCrossings = [&]() {
          size_t edges = 0, nearVertex = 0;
          Real minimum = 0.5;
          for (auto cell = mesh.getCell(); cell; ++cell)
          {
            const auto& vertices = cell->getVertices();
            for (size_t a = 0; a < vertices.size(); ++a)
              for (size_t b = a + 1; b < vertices.size(); ++b)
              {
                const Real t = crossingFraction(vertices[a], vertices[b]);
                if (std::isnan(t))
                  continue;
                ++edges;
                const Real distanceToVertex = std::min(t, 1 - t);
                minimum = std::min(minimum, distanceToVertex);
                if (distanceToVertex < Real(1e-3))
                  ++nearVertex;
              }
          }
          return std::tuple{edges, nearVertex, minimum};
        };
        const auto [cutEdges, nearVertexCuts, minimumCrossing] = scanCrossings();
        Alert::Info crossingInfo;
        crossingInfo << substageHeading("Level-set crossings") << Alert::NewLine
                     << diagnosticLabel("Crossed cell edges:")
                     << Alert::Notation::Number(cutEdges) << Alert::NewLine
                     << diagnosticLabel("Crossings within 1e-3 of a vertex:")
                     << Alert::Notation::Number(nearVertexCuts) << Alert::NewLine
                     << diagnosticLabel("Minimum crossing fraction:")
                     << Alert::Notation::Number(minimumCrossing);
        size_t snappedCount = 0;
        if (snap > 0)
        {
          std::vector<Real> original(
            sanitized.getData().begin(), sanitized.getData().end());
          std::vector<char> snapped(mesh.getVertexCount(), 0);
          for (auto cell = mesh.getCell(); cell; ++cell)
          {
            const auto& vertices = cell->getVertices();
            for (size_t a = 0; a < vertices.size(); ++a)
              for (size_t b = a + 1; b < vertices.size(); ++b)
              {
                const Index i = vertices[a], j = vertices[b];
                const Real t = original[i] * original[j] < 0
                  ? original[i] / (original[i] - original[j])
                  : Real(0.5);
                if (t < snap)
                  snapped[i] = 1;
                else if (t > 1 - snap)
                  snapped[j] = 1;
              }
          }
          for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
            if (snapped[vertex])
              sanitized[vertex] = 0;
          // MMG requires every tetrahedron to keep a vertex off the level set;
          // the vertex farthest from it keeps its value.
          size_t restored = 0;
          for (auto cell = mesh.getCell(); cell; ++cell)
          {
            const auto& vertices = cell->getVertices();
            Index farthest = vertices[0];
            bool allZero = true;
            for (const Index vertex : vertices)
            {
              allZero = allZero && sanitized[vertex] == 0;
              if (std::abs(original[vertex]) > std::abs(original[farthest]))
                farthest = vertex;
            }
            if (allZero)
            {
              sanitized[farthest] = original[farthest];
              ++restored;
            }
          }
          Real displacement = 0;
          size_t snappedVertices = 0;
          for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
          {
            if (sanitized[vertex] == 0 && original[vertex] != 0)
            {
              ++snappedVertices;
              ++snappedCount;
              displacement = std::max(displacement, std::abs(original[vertex]));
            }
          }
          const auto [snappedEdges, snappedNearVertex, snappedMinimum] = scanCrossings();
          crossingInfo << Alert::NewLine << diagnosticLabel("Snapping fraction:")
                       << Alert::Notation::Number(snap) << Alert::NewLine
                       << diagnosticLabel("Snapped vertices:")
                       << Alert::Notation::Number(snappedVertices) << Alert::NewLine
                       << diagnosticLabel("Restored in all-zero cells:")
                       << Alert::Notation::Number(restored) << Alert::NewLine
                       << diagnosticLabel("Maximum snapped level set:")
                       << Alert::Notation::Number(displacement) << " = "
                       << Alert::Notation::Number(displacement / h) << " h"
                       << Alert::NewLine
                       << diagnosticLabel("Minimum crossing after snapping:")
                       << Alert::Notation::Number(snappedMinimum);
        }
        crossingInfo << Alert::Raise;
        recordMMGStage(mesh, "before-cut", m_mmgReplayOuterRadius);
        if (!m_mmgReplayPrefix.empty())
          sanitized.save(m_mmgReplayPrefix + "-sanitized.sol", IO::FileFormat::MEDIT);
        MMG::Mesh reconstructed = discretizer.discretize(sanitized);

        splitSelfPairedCut(reconstructed);
        recordMMGStage(reconstructed, "after-cut", m_mmgReplayOuterRadius);
        const MeshDiagnostics reconstructionDiagnostics =
          getMeshDiagnostics(reconstructed, false);
        Alert::Info()
          << substageHeading("MMG reconstruction") << Alert::NewLine
          << diagnosticLabel("Minimum size:") << Alert::Notation::Number(hmin)
          << Alert::NewLine << diagnosticLabel("Maximum size:")
          << Alert::Notation::Number(hmax) << Alert::NewLine
          << diagnosticLabel("Hausdorff tolerance:") << Alert::Notation::Number(hausdorff)
          << Alert::NewLine << diagnosticLabel("Required boundary triangles:")
          << Alert::Notation::Number(requiredTriangles) << Alert::NewLine
          << diagnosticLabel("Cell count:") << Alert::Notation::Number(previousCells)
          << " -> " << Alert::Notation::Number(reconstructionDiagnostics.cells)
          << Alert::NewLine << diagnosticLabel("Obstacle cells:")
          << Alert::Notation::Number(reconstructionDiagnostics.obstacleCells)
          << Alert::NewLine << diagnosticLabel("Fluid cells:")
          << Alert::Notation::Number(reconstructionDiagnostics.fluidCells)
          << Alert::NewLine << diagnosticLabel("Interface triangles:")
          << Alert::Notation::Number(inputDiagnostics.interfaceTriangles) << " -> "
          << Alert::Notation::Number(reconstructionDiagnostics.interfaceTriangles)
          << Alert::NewLine << diagnosticLabel("Minimum tetrahedron quality:")
          << Alert::Notation::Number(inputDiagnostics.minimumQuality) << " -> "
          << Alert::Notation::Number(reconstructionDiagnostics.minimumQuality)
          << Alert::NewLine << diagnosticLabel("Mean tetrahedron quality:")
          << Alert::Notation::Number(inputDiagnostics.meanQuality) << " -> "
          << Alert::Notation::Number(reconstructionDiagnostics.meanQuality)
          << Alert::NewLine << diagnosticLabel("Mean element size:")
          << Alert::Notation::Number(inputDiagnostics.meanElementSize) << " -> "
          << Alert::Notation::Number(reconstructionDiagnostics.meanElementSize)
          << Alert::Raise;

        const size_t requiredTrianglesBeforeOptimization =
          sphere.protectFixedGeometry(reconstructed, false);
        recordMMGStage(reconstructed, "before-adapt", m_mmgReplayOuterRadius);
        if (adapt)
        {
          sphere.adapt(reconstructed, h, m_mmgReplayLockCuts);
        }
        else
        {
          MMG::Optimizer()
            .setHMin(hmin)
            .setHMax(hmax)
            .setHausdorff(hausdorff)
            .setGradation(remeshGradation)
            .setAngleDetection(false)
            .optimize(reconstructed);
          splitSelfPairedCut(reconstructed);
        }
        const size_t requiredTrianglesAfterOptimization =
          sphere.protectFixedGeometry(reconstructed, false);
        recordMMGStage(reconstructed, "after-adapt", m_mmgReplayOuterRadius);
        const MeshDiagnostics outputDiagnostics =
          getMeshDiagnostics(reconstructed, false);
        Alert::Info()
          << substageHeading(adapt ? "MMG adaptation" : "MMG optimization")
          << Alert::NewLine << diagnosticLabel("Required boundary triangles:")
          << Alert::Notation::Number(requiredTrianglesBeforeOptimization) << " -> "
          << Alert::Notation::Number(requiredTrianglesAfterOptimization) << Alert::NewLine
          << diagnosticLabel("Cell count:")
          << Alert::Notation::Number(reconstructionDiagnostics.cells) << " -> "
          << Alert::Notation::Number(outputDiagnostics.cells) << Alert::NewLine
          << diagnosticLabel("Minimum tetrahedron quality:")
          << Alert::Notation::Number(reconstructionDiagnostics.minimumQuality) << " -> "
          << Alert::Notation::Number(outputDiagnostics.minimumQuality) << Alert::NewLine
          << diagnosticLabel("Mean tetrahedron quality:")
          << Alert::Notation::Number(reconstructionDiagnostics.meanQuality) << " -> "
          << Alert::Notation::Number(outputDiagnostics.meanQuality) << Alert::NewLine
          << diagnosticLabel("Mean element size:")
          << Alert::Notation::Number(reconstructionDiagnostics.meanElementSize) << " -> "
          << Alert::Notation::Number(outputDiagnostics.meanElementSize) << Alert::Raise;
        ReconstructionDiagnostics diagnostics{hmin, hmax, hausdorff,
          requiredTrianglesAfterOptimization, previousCells, outputDiagnostics.cells};
        diagnostics.minimumCrossing = minimumCrossing;
        diagnostics.snappedVertices = static_cast<Real>(snappedCount);
        return {std::move(reconstructed), diagnostics};
      }

      template <class Space, class Density, class Output>
      GradientDiagnostics identifyGradient(const Space& shapeSpace,
        const Density& density, const KelvinBall::RotatedNitscheIntegrator& coupling,
        Output& gradient, Real regularizationLength, Real nitschePenalty,
        const char* name, const Math::Vector<Real>* additionalLoad = nullptr)
      {
        TrialFunction gradientTrial(shapeSpace);
        TestFunction testField(shapeSpace);
        auto interfaceNormal = FaceNormal(shapeSpace.getMesh());
        interfaceNormal.traceOf(Fluid);
        Problem identification(gradientTrial, testField);
        identification =
          Integral(regularizationLength * regularizationLength * Jacobian(gradientTrial),
            Jacobian(testField)) +
          Integral(gradientTrial, testField) -
          FaceIntegral(density, Dot(interfaceNormal, testField)).over(Gamma) +
          DirichletBC(gradientTrial, VectorFunction{0, 0, 0}).on(Outer);
        identification.assemble();
        auto& system = identification.getLinearSystem();
        coupling.assembleVector(shapeSpace, system,
          regularizationLength * regularizationLength, nitschePenalty,
          FlatSet<Attribute>{Outer});
        if (additionalLoad)
          system.getVector() += *additionalLoad;
        KelvinBall::solveDirect(identification);
        const Real linearResidual =
          (system.getOperator() * system.getSolution() - system.getVector()).norm() /
          std::max(system.getVector().norm(), Real(1));
        if (!std::isfinite(linearResidual) || linearResidual > linearResidualTolerance)
          throw std::runtime_error(
            "The gradient identification problem did not converge.");
        gradient = gradientTrial.getSolution();
        const Real jump = coupling.vectorJump(gradient);
        Alert::Info() << substageHeading(std::string(name) + " gradient")
                      << Alert::NewLine << diagnosticLabel("Linear residual:")
                      << Alert::Notation::Number(linearResidual) << Alert::NewLine
                      << diagnosticLabel("Rotated jump:") << Alert::Notation::Number(jump)
                      << Alert::Raise;
        return {linearResidual, jump};
      }

      int m_argc;
      char** m_argv;
  };
}

KelvinBall::KelvinBallOptimization::KelvinBallOptimization(int argc, char** argv)
  : m_argc(argc),
    m_argv(argv)
{}

int KelvinBall::KelvinBallOptimization::Implementation::run()
{
  const int argc = m_argc;
  char** argv = m_argv;
  Configuration configuration;
  size_t maxIterations = 1;
  bool saveMeshDiagnostic = false;
  bool geometryOnly = false;
  bool stateOnly = false;
  Real regularizationFactor = 4.0;
  Real stepFactor = 0.1;
  Real levelSetPenalty = 1;
  Real thicknessFactor = 2;
  Real thicknessWeight = 1;
  std::string thicknessMethod = "legacy";
  std::string restartMesh;
  std::string restartSnapshot;
  std::string replayMMG;
  size_t replayMMGCount = 3;
  Real restartTargetVolume = 0;
  size_t initialGlobalStep = 0;
  size_t motionEvery = 0;
  size_t motionFrames = 24;
  Math::SpatialVector<Real> motionForce(3);
  motionForce(0) = 0;
  motionForce(1) = 0;
  motionForce(2) = 1;
  size_t advectionQuadratureOrder = 8;
  std::string reconstructionMethod = "mmg";
  for (int argument = 1; argument < argc; ++argument)
  {
    const std::string_view mode(argv[argument]);
    if (configuration.parse(mode))
      continue;
    if (mode.rfind("--iterations=", 0) == 0)
      maxIterations = std::stoul(std::string(mode.substr(13)));
    else if (mode.rfind("--regularization=", 0) == 0)
      regularizationFactor = std::stod(std::string(mode.substr(17)));
    else if (mode.rfind("--step=", 0) == 0)
      stepFactor = std::stod(std::string(mode.substr(7)));
    else if (mode.rfind("--level-set-penalty=", 0) == 0)
      levelSetPenalty = std::stod(std::string(mode.substr(20)));
    else if (mode.rfind("--thickness-min=", 0) == 0)
      thicknessFactor = std::stod(std::string(mode.substr(16)));
    else if (mode.rfind("--thickness-weight=", 0) == 0)
      thicknessWeight = std::stod(std::string(mode.substr(19)));
    else if (mode.rfind("--thickness-method=", 0) == 0)
      thicknessMethod = std::string(mode.substr(19));
    else if (mode.rfind("--restart-mesh=", 0) == 0)
      restartMesh = std::string(mode.substr(15));
    else if (mode.rfind("--restart-snapshot=", 0) == 0)
      restartSnapshot = std::string(mode.substr(19));
    else if (mode.rfind("--replay-mmg=", 0) == 0)
      replayMMG = std::string(mode.substr(std::string_view("--replay-mmg=").size()));
    else if (mode.rfind("--replay-mmg-count=", 0) == 0)
      replayMMGCount = std::stoul(std::string(mode.substr(std::string_view("--replay-mmg-count=").size())));
    else if (mode == "--replay-mmg-lock-cuts")
      m_mmgReplayLockCuts = true;
    else if (mode.rfind("--restart-target-volume=", 0) == 0)
      restartTargetVolume = std::stod(std::string(mode.substr(24)));
    else if (mode.rfind("--initial-global-step=", 0) == 0)
      initialGlobalStep = std::stoul(std::string(mode.substr(22)));
    else if (mode.rfind("--motion-every=", 0) == 0)
      motionEvery = std::stoul(std::string(mode.substr(15)));
    else if (mode.rfind("--motion-frames=", 0) == 0)
      motionFrames = std::stoul(std::string(mode.substr(16)));
    else if (mode.rfind("--motion-force=", 0) == 0)
    {
      std::istringstream components{std::string(mode.substr(15))};
      std::string component;
      for (Eigen::Index i = 0; i < 3; ++i)
      {
        if (!std::getline(components, component, ','))
          throw std::runtime_error("--motion-force needs three components fx,fy,fz.");
        motionForce(i) = std::stod(component);
      }
    }
    else if (mode.rfind("--advection-quadrature=", 0) == 0)
      advectionQuadratureOrder = std::stoul(std::string(mode.substr(23)));
    else if (mode.rfind("--reconstruction=", 0) == 0)
      reconstructionMethod = std::string(mode.substr(17));
    else if (mode.rfind("--wngir-", 0) == 0 || mode.rfind("--quad-order=", 0) == 0 ||
      mode == "--trace" || mode.rfind("--trace=", 0) == 0 ||
      mode.rfind("--j-safe=", 0) == 0 || mode.rfind("--j-ls=", 0) == 0 ||
      mode.rfind("--j-min=", 0) == 0)
      continue;
    else if (mode == "--save-mesh")
      saveMeshDiagnostic = true;
    else if (mode == "--geometry-only")
      geometryOnly = true;
    else if (mode == "--state-only")
      stateOnly = true;
    else if (mode == "--help")
    {
      printUsage(argv[0]);
      return 0;
    }
    else
      throw std::runtime_error("Unknown KelvinBall option: " + std::string(mode));
  }
  configuration.finalize();
  const size_t points = configuration.points;
  const Real outerRadius = configuration.outerRadius;
  const Real nitschePenalty = configuration.nitschePenalty;
  const Real stabilizationFactor = configuration.stabilizationFactor;
  if (maxIterations == 0)
    throw std::runtime_error("The iteration count must be positive.");
  if (!(regularizationFactor > 0))
    throw std::runtime_error("The regularization factor must be positive.");
  if (!(stepFactor > 0))
    throw std::runtime_error("The advection step factor must be positive.");
  if (!(levelSetPenalty > 0))
    throw std::runtime_error("The level-set trace penalty must be positive.");
  if (!std::isfinite(thicknessFactor) || thicknessFactor < 0 ||
      !std::isfinite(thicknessWeight) || !(thicknessWeight > 0))
    throw std::runtime_error("The thickness options must satisfy --thickness-min >= 0 "
                             "and --thickness-weight > 0.");
  if (thicknessMethod != "legacy" && thicknessMethod != "baseline")
    throw std::runtime_error("--thickness-method must be legacy or baseline.");
  if (thicknessMethod == "baseline" && (KelvinBall::Tetrahedral || thicknessWeight != 1))
    throw std::runtime_error("Baseline thickness is Cubic-only and uses an adaptive multiplier, not --thickness-weight.");
  if (!restartMesh.empty() && (restartSnapshot.empty() || reconstructionMethod != "mmg" ||
      !std::isfinite(restartTargetVolume) || !(restartTargetVolume > 0)))
    throw std::runtime_error("Restart requires MMG, original --restart-snapshot and positive --restart-target-volume.");
  if (motionFrames == 0)
    throw std::runtime_error("The motion needs at least one frame.");
  if (!(motionForce.norm() > 0))
    throw std::runtime_error("The motion force must be nonzero.");
  if (advectionQuadratureOrder == 0)
    throw std::runtime_error("The advection quadrature order must be positive.");
  if (reconstructionMethod != "mmg" && reconstructionMethod != "wngir")
    throw std::runtime_error("The reconstruction method must be mmg or wngir.");
  if (configuration.adapt && reconstructionMethod != "mmg")
    throw std::runtime_error("--mmg-adapt applies only to --reconstruction=mmg.");
  if (configuration.mmgSnap > 0 && reconstructionMethod != "mmg")
    throw std::runtime_error("--mmg-snap applies only to --reconstruction=mmg.");
  if (KelvinBall::Tetrahedral && reconstructionMethod != "mmg")
    throw std::runtime_error("The tetrahedral extension supports MMG only.");
  if (geometryOnly && stateOnly)
    throw std::runtime_error("Use either --geometry-only or --state-only, not both.");
  if (!replayMMG.empty() && (restartMesh.empty() || geometryOnly || stateOnly ||
      reconstructionMethod != "mmg" || replayMMGCount == 0 || replayMMGCount > 3))
    throw std::runtime_error("MMG replay requires a saved restart and one to three saved advected fields.");
  if (m_mmgReplayLockCuts && replayMMG.empty())
    throw std::runtime_error("Cut locking is a diagnostic replay control, not an optimization option.");
  const Real h = configuration.getH();
  const Real hilbertLength = regularizationFactor * h;
  const Real dt = stepFactor * h;
  Alert::Info configurationInfo;
  configurationInfo << substageHeading("Configuration") << Alert::NewLine
                    << diagnosticLabel("Symmetry copies:")
                    << Alert::Notation::Number(KelvinBall::ChamberMultiplicity)
                    << Alert::NewLine << diagnosticLabel("Seed perturbation:")
                    << (KelvinBall::Tetrahedral ? "0.02*(H3+H6)/sqrt(2)" : "0.02*H3*H6")
                    << Alert::NewLine
                    << diagnosticLabel("Grid points:") << Alert::Notation::Number(points)
                    << Alert::NewLine << diagnosticLabel("Outer radius:")
                    << Alert::Notation::Number(outerRadius) << Alert::NewLine
                    << diagnosticLabel("Effective h:") << Alert::Notation::Number(h)
                    << Alert::NewLine << diagnosticLabel("Evaluated designs:")
                    << Alert::Notation::Number(maxIterations) << Alert::NewLine
                    << diagnosticLabel("Nitsche penalty:")
                    << Alert::Notation::Number(nitschePenalty) << Alert::NewLine
                    << diagnosticLabel("Stabilization factor:")
                    << Alert::Notation::Number(stabilizationFactor) << Alert::NewLine
                    << diagnosticLabel("Assembly backend:") << assemblyBackend
                    << Alert::NewLine << diagnosticLabel("Direct solver:")
                    << KelvinBall::DirectSolverName << Alert::NewLine
                    << diagnosticLabel("Regularization length:")
                    << Alert::Notation::Number(hilbertLength) << " = "
                    << Alert::Notation::Number(regularizationFactor) << " h"
                    << Alert::NewLine << diagnosticLabel("Advection step:")
                    << Alert::Notation::Number(dt) << " = "
                    << Alert::Notation::Number(stepFactor) << " h" << Alert::NewLine
                    << diagnosticLabel("Level-set trace penalty:")
                    << Alert::Notation::Number(levelSetPenalty) << Alert::NewLine
                    << diagnosticLabel("Advection quadrature order:")
                    << Alert::Notation::Number(advectionQuadratureOrder) << Alert::NewLine
                    << diagnosticLabel("Reconstruction method:") << reconstructionMethod;
  if (reconstructionMethod == "wngir")
  {
    configurationInfo << Alert::NewLine << diagnosticLabel("Background minimum size:")
                      << Alert::Notation::Number(configuration.backgroundHMin) << " h"
                      << Alert::NewLine << diagnosticLabel("Background maximum size:")
                      << Alert::Notation::Number(configuration.backgroundHMax) << " h"
                      << Alert::NewLine << diagnosticLabel("Background Hausdorff:")
                      << Alert::Notation::Number(configuration.backgroundHausdorff)
                      << " h" << Alert::NewLine
                      << diagnosticLabel("Background gradation:")
                      << Alert::Notation::Number(configuration.backgroundGradation);
  }
  else if (configuration.adapt)
  {
    configurationInfo << Alert::NewLine << diagnosticLabel("Adaptation interface size:")
                      << Alert::Notation::Number(configuration.adaptInterfaceSize) << " h"
                      << Alert::NewLine << diagnosticLabel("Adaptation far size:")
                      << Alert::Notation::Number(configuration.adaptFarSize) << " h"
                      << Alert::NewLine << diagnosticLabel("Adaptation width:")
                      << Alert::Notation::Number(configuration.adaptWidth) << " h"
                      << Alert::NewLine << diagnosticLabel("Adaptation gradation:")
                      << Alert::Notation::Number(configuration.adaptGradation);
  }
  if (configuration.mmgSnap > 0)
    configurationInfo << Alert::NewLine << diagnosticLabel("Level-set snapping fraction:")
                      << Alert::Notation::Number(configuration.mmgSnap);
  configurationInfo << Alert::Raise;
  const Real nan = std::numeric_limits<Real>::quiet_NaN();
  const auto stage1Start = Clock::now();
  announce(!restartMesh.empty()
      ? "Stage 1: Loading the saved volume mesh without reconstruction."
      : reconstructionMethod == "wngir"
      ? "Stage 1: Preparing the background mesh and fitting the initial sphere with "
        "WNGIR."
      : "Stage 1: Discretizing the initial sphere with MMG.");
  ReconstructionDiagnostics reconstruction;
  Optional<MMG::Mesh> wngirBackground;
  MMG::Mesh mesh;
  Sphere sphere(configuration);
  if (!restartMesh.empty())
  {
    mesh.load(restartMesh, IO::FileFormat::MEDIT);
    // Boost.Spirit's decimal parser can differ by a final bit. Restore the
    // exact original doubles before constructing any state/distance spaces.
    const auto originalFile = IO::HDF5::File(
      H5Fopen(restartSnapshot.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    if (!originalFile)
      throw std::runtime_error("Cannot read original restart snapshot.");
    const auto [count, dimension] = IO::HDF5::readMatrixShape(
      originalFile.get(), "/Mesh/Geometry/Vertices");
    const auto coordinates = IO::HDF5::readVectorDataset<IO::HDF5::F64>(
      originalFile.get(), "/Mesh/Geometry/Vertices");
    if (count != mesh.getVertexCount() || dimension != 3 || coordinates.size() != 3 * count)
      throw std::runtime_error("Restart snapshot has different vertices.");
    for (Index vertex = 0; vertex < count; ++vertex)
      mesh.setVertexCoordinates(vertex, Math::SpatialPoint{
        coordinates[3 * vertex], coordinates[3 * vertex + 1], coordinates[3 * vertex + 2]});
    // Complete lazy topology before recording all edge/face attribute slots.
    mesh.getConnectivity().compute(3, 1);
    mesh.getConnectivity().compute(3, 2);
    IO::XDMF restartGeometry("RestartGeometry");
    auto original = restartGeometry.grid("Chamber");
    original.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
    restartGeometry.write(0).flush();
  }
  else
  {
    SphereDiscretization initial = reconstructionMethod == "wngir"
      ? sphere.prepareWNGIRBackground()
      : sphere.discretize();
    reconstruction = initial.diagnostics;
    if (reconstructionMethod == "wngir")
    {
      wngirBackground.emplace(std::move(initial.mesh));
      P1 sphereSpace(*wngirBackground);
      GridFunction sphereLevelSet(sphereSpace);
      sphereLevelSet = RealFunction([](const Geometry::Point& point) {
        return point.getPhysicalCoordinates().norm() - Real(1);
      });
      MMG::Mesh classified = classifyLevelSetForWNGIR(*wngirBackground, sphereLevelSet);
      P1 classifiedSphereSpace(classified);
      GridFunction classifiedSphereLevelSet(classifiedSphereSpace);
      classifiedSphereLevelSet.getData() = sphereLevelSet.getData();
      MMGReconstruction fitted =
        fitLevelSetWNGIR(classified, classifiedSphereLevelSet, h, outerRadius, argc, argv);
      mesh = std::move(fitted.mesh);
    }
    else
    {
      mesh = std::move(initial.mesh);
    }
  }
  checkFixedGeometry(mesh, outerRadius);
  checkMaterials(mesh);
  if (saveMeshDiagnostic)
    mesh.save("KelvinBallInitial.mesh", IO::FileFormat::MEDIT);
  if (restartMesh.empty())
    reportSphereGeometry(mesh);
  const Real stage1Seconds = elapsedSeconds(stage1Start);
  reportStageTiming(1, stage1Seconds);
  if (!replayMMG.empty())
  {
    m_mmgReplayOuterRadius = outerRadius;
    // Replay recorded transport fields only: never solve Stokes, identify a
    // direction, advect, or continue optimization. Stop before using a field
    // if reconstructed vertex/cell ordering differs from its saved input.
    for (size_t iteration = 0; iteration < replayMMGCount; ++iteration)
    {
      char suffix[32];
      std::snprintf(suffix, sizeof(suffix), "%06zu.h5", iteration);
      const std::string expectedPath = replayMMG + "/KelvinBall.Chamber.mesh." + suffix;
      const auto expected = IO::HDF5::File(H5Fopen(expectedPath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
      if (!expected)
        throw std::runtime_error("Missing saved replay mesh: " + expectedPath);
      const auto coordinates = IO::HDF5::readVectorDataset<IO::HDF5::F64>(
        expected.get(), "/Mesh/Geometry/Vertices");
      const auto topology = IO::HDF5::readVectorDataset<IO::HDF5::U64>(
        expected.get(), "/Mesh/XDMF/Topology");
      bool identical = coordinates.size() == 3 * mesh.getVertexCount() &&
        topology.size() == 5 * mesh.getCellCount();
      for (Index vertex = 0; identical && vertex < mesh.getVertexCount(); ++vertex)
        for (size_t d = 0; d < 3; ++d)
          identical = identical && mesh.getVertexCoordinates(vertex)(d) == coordinates[3 * vertex + d];
      for (auto cell = mesh.getCell(); identical && cell; ++cell)
      {
        const auto& vertices = cell->getVertices();
        for (size_t k = 0; k < 4; ++k)
          identical = identical && vertices[k] == topology[5 * cell->getIndex() + 1 + k];
      }
      m_mmgReplayPrefix = "replay-" + std::to_string(iteration);
      recordMMGStage(mesh, "input", outerRadius);
      std::ofstream(m_mmgReplayPrefix + "-input-exact.txt") << identical << '\n';
      if (!identical)
        throw std::runtime_error("Replay mesh differs from saved field ordering; refusing to use field.");
      auto& connectivity = mesh.getConnectivity();
      connectivity.discover(3, 2);
      connectivity.discover(3, 1);
      connectivity.restrict(1, 0);
      connectivity.restrict(2, 0);
      connectivity.restrict(2, 3);
      connectivity.discover(0, 0);
      P1 levelSetSpace(mesh);
      MMG::RealGridFunction advected(levelSetSpace);
      const std::string fieldPath = replayMMG + "/KelvinBall.Chamber.Advected." + suffix;
      const auto field = IO::HDF5::File(H5Fopen(fieldPath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
      if (!field)
        throw std::runtime_error("Missing saved advected field: " + fieldPath);
      const auto values = IO::HDF5::readVectorDataset<IO::HDF5::F64>(field.get(), "/GridFunction/Values/Data");
      if (values.size() != mesh.getVertexCount())
        throw std::runtime_error("Saved advected field has different vertices.");
      std::copy(values.begin(), values.end(), advected.getData().data());
      auto result = discretizeLevelSetMMG(mesh, advected, h, sphere,
        configuration.adapt, configuration.mmgSnap);
      mesh = std::move(result.mesh);
      recordMMGStage(mesh, "output", outerRadius);
      checkFixedGeometry(mesh, outerRadius);
    }
    return 0;
  }
  if (geometryOnly)
    return 0;
  const Real targetVolume = restartMesh.empty() ? mesh.getVolume(Obstacle) : restartTargetVolume;

  std::ofstream history("kelvin-ball.csv");
  history.precision(17);
  history << "iteration,outer_radius,h,dt,advection_quadrature_order,"
             "regularization_length,nitsche_penalty,"
             "stabilization_factor,assembly_backend,direct_solver,"
             "vertices,cells,obstacle_cells,fluid_cells,"
             "interface_triangles,mesh_quality_min,mesh_quality_mean,mesh_quality_max,"
             "mean_element_size,"
             "remesh_hmin,remesh_hmax,remesh_hausdorff,"
             "required_boundary_triangles,remesh_cells_before,remesh_cells_after,"
             "k,c,q,rho,coupling_symmetry,nitsche_jump,volume,target_volume,"
             "volume_error,chamber_bbox_x,chamber_bbox_y,chamber_bbox_z,"
             "sewn_bbox_x,sewn_bbox_y,sewn_bbox_z,actual_delta_rho,"
             "predicted_delta_rho,actual_delta_volume,predicted_delta_volume,"
             "null_space_multiplier,xi_rho_inf_norm,theta_inf_norm,d_rho_theta,"
             "d_volume_theta,required_d_volume_theta,rho_gradient_residual,"
             "rho_gradient_jump,volume_gradient_residual,"
             "volume_gradient_jump,stage_1_seconds,stage_2_seconds,"
             "stage_3_seconds,stage_4_seconds,stage_5_seconds,"
             "stage_6_seconds,stage_7_seconds,stage_8_seconds,"
             "stage_9_seconds,"
             "translational_mobility,coupling_mobility,rotational_mobility,"
             "revolution_period,pitch,resistance_radius,"
             "level_set_penalty,thickness_min,thickness_weight,thickness_penalty,"
             "thickness_violating_rays,thickness_deepest_exit,"
             "eikonal_rotated_jump,projected_rotated_jump,distance_correction,"
             "interface_shift_max,advection_increment,advected_rotated_jump,"
             "min_crossing_fraction,snapped_vertices,reconstruction_scale,"
             "z_x,z_y,z_z,omega_x,omega_y,omega_z,global_iteration,"
             "thickness_nominal_penalty,thickness_guard_distance,thickness_minimum_exit,"
             "thickness_maximum_deficit,thickness_active_multiplier,thickness_active_feasible,"
             "thickness_repair_d_volume,d_thickness_theta\n";
  IO::XDMF xdmf("KelvinBall");
  auto chamber = xdmf.grid("Chamber");
  chamber.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
  auto fluidState = xdmf.grid("Fluid");
  IO::XDMF sewedXdmf("KelvinBallSewed");
  auto sewedDesignOutput = sewedXdmf.grid("Design");
  auto sewedFluidOutput = sewedXdmf.grid("Fluid");
  const std::string reconstructionName =
    reconstructionMethod == "wngir" ? "KelvinBallWNGIR" : "KelvinBallMMG";
  IO::XDMF reconstructionXdmf(reconstructionName);
  auto reconstructionOutput = reconstructionXdmf.grid("Reconstructed");
  Optional<Real> previousRho;
  Optional<Real> predictedRhoChange;
  Optional<Real> previousVolume;
  Optional<Real> predictedVolumeChange;
  Optional<MMG::Mesh> nextMesh;

  for (size_t iteration = 0; iteration < maxIterations; ++iteration)
  {
    std::array<Real, 9> stageSeconds;
    stageSeconds.fill(nan);
    if (iteration == 0)
      stageSeconds[0] = stage1Seconds;
    if (nextMesh)
    {
      mesh = std::move(*nextMesh);
      nextMesh.reset();
    }
    Alert::Info() << " ------------------------------------------------------------"
                  << Alert::NewLine << " Iteration "
                  << Alert::Notation::Number(iteration + 1) << " of "
                  << Alert::Notation::Number(maxIterations) << Alert::NewLine
                  << " ------------------------------------------------------------"
                  << Alert::Raise;
    const auto stage2Start = Clock::now();
    announce("Stage 2: Preparing the chamber fluid mesh.");
    auto& connectivity = mesh.getConnectivity();
    connectivity.discover(3, 2);
    connectivity.discover(3, 1);
    connectivity.restrict(1, 0);
    connectivity.restrict(2, 0);
    connectivity.restrict(2, 3);
    connectivity.discover(0, 0);
    SubMesh fluid = mesh.trim(Obstacle);
    splitSelfPairedCut(fluid);
    auto& fluidConnectivity = fluid.getConnectivity();
    fluidConnectivity.discover(2, 1);
    fluidConnectivity.restrict(1, 0);
    fluidConnectivity.restrict(2, 3);

    VelocitySpace Vh = makeVelocitySpace(fluid);
    PressureSpace Qh(fluid);
    KelvinBall::RotatedNitscheIntegrator fluidCoupling(
      fluid, MasterCuts, rotatedTracePhysicalTolerance, rotatedTraceReferenceTolerance);
    GridFunction uT0(Vh), uT1(Vh), uT2(Vh), uR0(Vh), uR1(Vh), uR2(Vh);
    GridFunction pT0(Qh), pT1(Qh), pT2(Qh), pR0(Qh), pR1(Qh), pR2(Qh);
    stageSeconds[1] = elapsedSeconds(stage2Start);
    reportStageTiming(2, stageSeconds[1]);

    const auto stage3Start = Clock::now();
    announce("Stage 3: Solving the translational and rotational Stokes states.");
    const KelvinBall::Parameters metricParameters{h, nitschePenalty, stabilizationFactor};
    const KelvinBall::Metrics resistanceMetrics(metricParameters);
    const KelvinBall::Values metrics = resistanceMetrics.evaluateChamber(
      Vh, Qh, fluidCoupling, uT0, uT1, uT2, uR0, uR1, uR2, pT0, pT1, pT2, pR0, pR1, pR2);
    stageSeconds[2] = elapsedSeconds(stage3Start);
    reportStageTiming(3, stageSeconds[2]);
    const Real k = metrics.k;
    const Real c = metrics.c;
    const Real q = metrics.q;
    const Real rho = metrics.rho;
    const Real couplingSymmetry = metrics.couplingSymmetry;
    const Real volume = mesh.getVolume(Obstacle);
    const Real nitscheJump = metrics.nitscheJump;
    const auto chamberBoundingBox = boundingBoxSize(mesh);
    const auto sewnBoundingBox = sewnBoundingBoxSize(mesh);
    const auto meshDiagnostics = getMeshDiagnostics(mesh);
    const Real actualRhoChange = previousRho ? rho - *previousRho : nan;
    const Real incomingPredictedRhoChange =
      predictedRhoChange ? *predictedRhoChange : nan;
    const Real actualVolumeChange = previousVolume ? volume - *previousVolume : nan;
    const Real incomingPredictedVolumeChange =
      predictedVolumeChange ? *predictedVolumeChange : nan;
    // Diagnostics of the later stages, filled as they run.
    struct
    {
        Real thicknessPenalty = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessViolating = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessDeepest = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessNominalPenalty = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessGuardDistance = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessMinimumExit = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessMaximumDeficit = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessActiveMultiplier = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessActiveFeasible = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessRepairVolume = std::numeric_limits<Real>::quiet_NaN();
        Real dThicknessTheta = std::numeric_limits<Real>::quiet_NaN();
        Real eikonalJump = std::numeric_limits<Real>::quiet_NaN();
        Real projectedJump = std::numeric_limits<Real>::quiet_NaN();
        Real distanceCorrection = std::numeric_limits<Real>::quiet_NaN();
        Real interfaceShift = std::numeric_limits<Real>::quiet_NaN();
        Real advectionIncrement = std::numeric_limits<Real>::quiet_NaN();
        Real advectedJump = std::numeric_limits<Real>::quiet_NaN();
    } stageDiagnostics;
    const auto writeHistory = [&](Real nullSpaceMultiplier, Real xiRhoInfinityNorm,
                                Real thetaInfinityNorm, Real dRhoTheta, Real dVolumeTheta,
                                Real requiredDVolumeTheta,
                                const GradientDiagnostics& rhoGradientDiagnostics,
                                const GradientDiagnostics& volumeGradientDiagnostics) {
      history << iteration << ',' << outerRadius << ',' << h << ',' << dt << ','
              << advectionQuadratureOrder << ',' << hilbertLength << ',' << nitschePenalty
              << ',' << stabilizationFactor << ',' << assemblyBackend << ','
              << KelvinBall::DirectSolverName << ',' << meshDiagnostics.vertices << ','
              << meshDiagnostics.cells << ',' << meshDiagnostics.obstacleCells << ','
              << meshDiagnostics.fluidCells << ',' << meshDiagnostics.interfaceTriangles
              << ',' << meshDiagnostics.minimumQuality << ','
              << meshDiagnostics.meanQuality << ',' << meshDiagnostics.maximumQuality
              << ',' << meshDiagnostics.meanElementSize << ','
              << reconstruction.minimumSize << ',' << reconstruction.maximumSize << ','
              << reconstruction.hausdorffTolerance << ','
              << reconstruction.requiredBoundaryTriangles << ','
              << reconstruction.cellsBefore << ',' << reconstruction.cellsAfter << ','
              << k << ',' << c << ',' << q << ',' << rho << ',' << couplingSymmetry << ','
              << nitscheJump << ',' << volume << ',' << targetVolume << ','
              << volume - targetVolume << ',' << chamberBoundingBox(0) << ','
              << chamberBoundingBox(1) << ',' << chamberBoundingBox(2) << ','
              << sewnBoundingBox(0) << ',' << sewnBoundingBox(1) << ','
              << sewnBoundingBox(2) << ',' << actualRhoChange << ','
              << incomingPredictedRhoChange << ',' << actualVolumeChange << ','
              << incomingPredictedVolumeChange << ',' << nullSpaceMultiplier << ','
              << xiRhoInfinityNorm << ',' << thetaInfinityNorm << ',' << dRhoTheta << ','
              << dVolumeTheta << ',' << requiredDVolumeTheta << ','
              << rhoGradientDiagnostics.residual << ',' << rhoGradientDiagnostics.jump
              << ',' << volumeGradientDiagnostics.residual << ','
              << volumeGradientDiagnostics.jump;
      for (const Real seconds : stageSeconds)
        history << ',' << seconds;
      // Free motion under a unit force without torque, R [Z; w] = [F; 0].
      const Real determinant = k * q - c * c;
      history << ',' << q / determinant << ',' << -c / determinant << ','
              << k / determinant << ','
              << (c != 0 ? 2 * M_PI * determinant / std::abs(c) : nan) << ','
              << (c != 0 ? 2 * M_PI * q / std::abs(c) : nan) << ',' << std::sqrt(q / k)
              << ',' << levelSetPenalty << ',' << thicknessFactor * h << ','
              << thicknessWeight << ',' << stageDiagnostics.thicknessPenalty << ','
              << stageDiagnostics.thicknessViolating << ','
              << stageDiagnostics.thicknessDeepest << ',' << stageDiagnostics.eikonalJump
              << ',' << stageDiagnostics.projectedJump << ','
              << stageDiagnostics.distanceCorrection << ','
              << stageDiagnostics.interfaceShift << ','
              << stageDiagnostics.advectionIncrement << ','
              << stageDiagnostics.advectedJump << ',' << reconstruction.minimumCrossing
              << ',' << reconstruction.snappedVertices << ',' << reconstruction.scale;
      // Z and omega themselves, for the unit force along --motion-force.
      const Math::SpatialVector<Real> unitForce = motionForce / motionForce.norm();
      for (Eigen::Index component = 0; component < 3; ++component)
        history << ',' << q / determinant * unitForce(component);
      for (Eigen::Index component = 0; component < 3; ++component)
        history << ',' << -c / determinant * unitForce(component);
      history << ',' << initialGlobalStep + iteration
              << ',' << stageDiagnostics.thicknessNominalPenalty
              << ',' << stageDiagnostics.thicknessGuardDistance
              << ',' << stageDiagnostics.thicknessMinimumExit
              << ',' << stageDiagnostics.thicknessMaximumDeficit
              << ',' << stageDiagnostics.thicknessActiveMultiplier
              << ',' << stageDiagnostics.thicknessActiveFeasible
              << ',' << stageDiagnostics.thicknessRepairVolume
              << ',' << stageDiagnostics.dThicknessTheta << '\n';
      history.flush();
    };
    const auto stage4Start = Clock::now();
    announce("Stage 4: Evaluating the resistance objective.");
    if (previousRho)
    {
      Alert::Info() << substageHeading("Previous update") << Alert::NewLine
                    << diagnosticLabel("Actual delta rho:")
                    << Alert::Notation::Number(actualRhoChange) << Alert::NewLine
                    << diagnosticLabel("Predicted delta rho:")
                    << Alert::Notation::Number(incomingPredictedRhoChange)
                    << Alert::NewLine << diagnosticLabel("Actual delta volume:")
                    << Alert::Notation::Number(actualVolumeChange) << Alert::NewLine
                    << diagnosticLabel("Predicted delta volume:")
                    << Alert::Notation::Number(incomingPredictedVolumeChange)
                    << Alert::Raise;
    }
    Alert::Info() << substageHeading("Resistance metrics") << Alert::NewLine
                  << diagnosticLabel("Translational resistance k:")
                  << Alert::Notation::Number(k) << Alert::NewLine
                  << diagnosticLabel("Coupling coefficient c:")
                  << Alert::Notation::Number(c) << Alert::NewLine
                  << diagnosticLabel("Rotational resistance q:")
                  << Alert::Notation::Number(q) << Alert::NewLine
                  << diagnosticLabel("Objective rho:") << Alert::Notation::Number(rho)
                  << Alert::NewLine << diagnosticLabel("Volume:")
                  << Alert::Notation::Number(volume) << Alert::NewLine
                  << diagnosticLabel("Chamber bounding box:")
                  << vectorText(chamberBoundingBox) << Alert::NewLine
                  << diagnosticLabel("Sewn mesh bounding box:")
                  << vectorText(sewnBoundingBox) << Alert::NewLine
                  << diagnosticLabel("Nitsche jump:")
                  << Alert::Notation::Number(nitscheJump) << Alert::NewLine
                  << diagnosticLabel("Coupling symmetry residual:")
                  << Alert::Notation::Number(couplingSymmetry) << Alert::Raise;
    {
      // Free motion under a unit force along --motion-force, without torque.
      const Real determinant = k * q - c * c;
      const Math::SpatialVector<Real> force = motionForce / motionForce.norm();
      const Math::SpatialVector<Real> translation = (q / determinant) * force;
      const Math::SpatialVector<Real> angular = (-c / determinant) * force;
      Alert::Info() << substageHeading("Free motion under a unit force") << Alert::NewLine
                    << diagnosticLabel("Force direction:") << vectorText(force)
                    << Alert::NewLine << diagnosticLabel("Translational velocity Z:")
                    << vectorText(translation) << Alert::NewLine
                    << diagnosticLabel("Angular velocity omega:") << vectorText(angular)
                    << Alert::NewLine << diagnosticLabel("Period of one revolution:")
                    << Alert::Notation::Number(
                         c != 0 ? 2 * M_PI * determinant / std::abs(c) : nan)
                    << Alert::NewLine << diagnosticLabel("Pitch:")
                    << Alert::Notation::Number(c != 0 ? 2 * M_PI * q / std::abs(c) : nan)
                    << Alert::NewLine << diagnosticLabel("Resistance length sqrt(q/k):")
                    << Alert::Notation::Number(std::sqrt(q / k)) << Alert::Raise;
    }
    stageSeconds[3] = elapsedSeconds(stage4Start);
    reportStageTiming(4, stageSeconds[3]);
    if (stateOnly)
    {
      const GradientDiagnostics unavailable{nan, nan};
      writeHistory(nan, nan, nan, nan, nan, nan, unavailable, unavailable);
      return 0;
    }

    const auto stage5Start = Clock::now();
    announce("Stage 5: Computing the shape derivative.");
    auto gradUT0 = Jacobian(uT0);
    gradUT0.traceOf(Fluid);
    auto gradUT1 = Jacobian(uT1);
    gradUT1.traceOf(Fluid);
    auto gradUT2 = Jacobian(uT2);
    gradUT2.traceOf(Fluid);
    auto gradUR0 = Jacobian(uR0);
    gradUR0.traceOf(Fluid);
    auto gradUR1 = Jacobian(uR1);
    gradUR1.traceOf(Fluid);
    auto gradUR2 = Jacobian(uR2);
    gradUR2.traceOf(Fluid);
    const auto strainUT0 = 0.5 * (gradUT0 + gradUT0.T());
    const auto strainUT1 = 0.5 * (gradUT1 + gradUT1.T());
    const auto strainUT2 = 0.5 * (gradUT2 + gradUT2.T());
    const auto strainUR0 = 0.5 * (gradUR0 + gradUR0.T());
    const auto strainUR1 = 0.5 * (gradUR1 + gradUR1.T());
    const auto strainUR2 = 0.5 * (gradUR2 + gradUR2.T());
    const Real signC = c < 0 ? -1.0 : 1.0;
    const auto kDensity = chamberMultiplicity * 2.0 * mu / 3.0 *
      (Dot(strainUT0, strainUT0) + Dot(strainUT1, strainUT1) + Dot(strainUT2, strainUT2));
    const auto cDensity = chamberMultiplicity * 2.0 * mu / 3.0 *
      (Dot(strainUT0, strainUR0) + Dot(strainUT1, strainUR1) + Dot(strainUT2, strainUR2));
    const auto qDensity = chamberMultiplicity * 2.0 * mu / 3.0 *
      (Dot(strainUR0, strainUR0) + Dot(strainUR1, strainUR1) + Dot(strainUR2, strainUR2));
    const auto rhoDensity =
      signC / std::sqrt(k * q) * cDensity - 0.5 * rho * (kDensity / k + qDensity / q);

    P1 levelSetSpace(mesh);
    P1 shapeSpace(mesh, 3);
    KelvinBall::RotatedNitscheIntegrator shapeCoupling(mesh,
      FlatSet<Attribute>{SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus},
      rotatedTracePhysicalTolerance, rotatedTraceReferenceTolerance);
    GridFunction rhoGradient(shapeSpace), volumeGradient(shapeSpace);
    stageSeconds[4] = elapsedSeconds(stage5Start);
    reportStageTiming(5, stageSeconds[4]);
    const auto stage6Start = Clock::now();
    announce("Stage 6: Regularizing and constraining the shape direction.");
    // Legacy mixes the load before projection. The baseline method instead
    // identifies and volume-projects its repair independently below.
    Math::Vector<Real> thicknessLoad;
    GridFunction thicknessDescent(shapeSpace);
    if (thicknessFactor > 0)
    {
      thicknessLoad = Math::Vector<Real>::Zero(shapeSpace.getSize());
      mesh.getConnectivity().compute(mesh.getDimension() - 1, mesh.getDimension());
      if (thicknessMethod == "baseline")
      {
        const Real nominal = thicknessFactor * h;
        const Real guard = nominal + Real(2) * dt;
        const KelvinBall::BaselineThicknessPenalty penalty(mesh, guard);
        const auto thickness = penalty.evaluate(mesh, shapeSpace, thicknessLoad, nominal);
        stageDiagnostics.thicknessPenalty = thickness.penalty;
        stageDiagnostics.thicknessViolating = static_cast<Real>(thickness.violating);
        stageDiagnostics.thicknessNominalPenalty = thickness.nominalPenalty;
        stageDiagnostics.thicknessGuardDistance = guard;
        stageDiagnostics.thicknessMinimumExit = thickness.minimumExit;
        stageDiagnostics.thicknessMaximumDeficit = thickness.maximumDeficit;
        stageDiagnostics.thicknessDeepest = thickness.maximumDeficit;
      }
      else
      {
        const KelvinBall::ThicknessPenalty thicknessPenalty(mesh, thicknessFactor * h);
        const Location::AABB<MMG::Mesh> chamberLocator(mesh);
        const auto thickness = thicknessPenalty.evaluate(
          mesh, chamberLocator, shapeSpace, thicknessWeight, thicknessLoad);
        stageDiagnostics.thicknessPenalty = thickness.penalty;
        stageDiagnostics.thicknessViolating = static_cast<Real>(thickness.violating);
        stageDiagnostics.thicknessDeepest = thickness.deepest;
      }
      Alert::Info() << substageHeading("Thickness penalty") << Alert::NewLine
                    << diagnosticLabel("Minimum thickness:")
                    << Alert::Notation::Number(thicknessFactor * h) << " = "
                    << Alert::Notation::Number(thicknessFactor) << " h" << Alert::NewLine
                    << diagnosticLabel("Weight:")
                    << Alert::Notation::Number(thicknessWeight) << Alert::NewLine
                    << diagnosticLabel("Penalty P:")
                    << Alert::Notation::Number(stageDiagnostics.thicknessPenalty) << Alert::NewLine
                    << diagnosticLabel(thicknessMethod == "baseline" ? "Violating ray samples:" : "Rays leaving the body:")
                    << Alert::Notation::Number(stageDiagnostics.thicknessViolating) << Alert::NewLine
                    << diagnosticLabel(thicknessMethod == "baseline" ? "Maximum thickness deficit:" : "Deepest exit:")
                    << Alert::Notation::Number(stageDiagnostics.thicknessDeepest) << Alert::Raise;
    }
    if (thicknessFactor > 0 && thicknessMethod == "baseline")
      identifyGradient(shapeSpace, RealFunction{0}, shapeCoupling, thicknessDescent,
        hilbertLength, nitschePenalty, "Thickness descent", &thicknessLoad);
    const GradientDiagnostics rhoGradientDiagnostics =
      identifyGradient(shapeSpace, -rhoDensity, shapeCoupling, rhoGradient, hilbertLength,
        nitschePenalty, "Rho", thicknessFactor > 0 && thicknessMethod == "legacy" ? &thicknessLoad : nullptr);
    const GradientDiagnostics volumeGradientDiagnostics =
      identifyGradient(shapeSpace, RealFunction{-1}, shapeCoupling, volumeGradient,
        hilbertLength, nitschePenalty, "Volume");
    TestFunction testField(shapeSpace);
    auto interfaceNormal = FaceNormal(mesh);
    interfaceNormal.traceOf(Fluid);
    LinearForm dRho(testField), dVolume(testField);
    dRho = FaceIntegral(-rhoDensity, Dot(interfaceNormal, testField)).over(Gamma);
    dVolume = FaceIntegral(RealFunction{-1}, Dot(interfaceNormal, testField)).over(Gamma);
    dRho.assemble();
    dVolume.assemble();
    const Real volumeMetric = dVolume(volumeGradient);
    if (!std::isfinite(volumeMetric) || std::abs(volumeMetric) <= Real(1e-20))
      throw std::runtime_error("The volume constraint derivative vanished.");
    Real nullSpaceMultiplier = -dVolume(rhoGradient) / volumeMetric;
    GridFunction xiRho(shapeSpace);
    xiRho = rhoGradient + nullSpaceMultiplier * volumeGradient;
    if (thicknessFactor > 0 && thicknessMethod == "baseline")
    {
      const Real thicknessVolumeMultiplier = -dVolume(thicknessDescent) / volumeMetric;
      GridFunction thicknessRepair(shapeSpace);
      thicknessRepair = thicknessDescent + thicknessVolumeMultiplier * volumeGradient;
      stageDiagnostics.thicknessRepairVolume = dVolume(thicknessRepair);
      GridFunction candidate(shapeSpace);
      GridFunction candidateMagnitude(levelSetSpace);
      const auto thicknessRate = [&](Real multiplier) {
        candidate = xiRho + multiplier * thicknessRepair;
        candidateMagnitude = Frobenius(candidate);
        const Real magnitude = std::max(candidateMagnitude.max(), Real(1e-30));
        return -thicknessLoad.dot(candidate.getData()) / magnitude;
      };
      const Real targetRate =
        -stageDiagnostics.thicknessPenalty / (thicknessFactor * h);
      Real multiplier = 0;
      bool feasible = true;
      if (thicknessRate(0) > targetRate)
      {
        candidateMagnitude = Frobenius(thicknessRepair);
        const Real repairMagnitude = std::max(candidateMagnitude.max(), Real(1e-30));
        const Real maximumRepairRate =
          -thicknessLoad.dot(thicknessRepair.getData()) / repairMagnitude;
        if (maximumRepairRate >= targetRate)
        {
          feasible = false;
          multiplier = std::numeric_limits<Real>::infinity();
        }
        else
        {
          Real lower = 0, upper = 1;
          while (thicknessRate(upper) > targetRate && upper < Real(1e6))
            upper *= Real(2);
          if (thicknessRate(upper) > targetRate)
          {
            feasible = false;
            multiplier = std::numeric_limits<Real>::infinity();
          }
          else
          {
            for (size_t search = 0; search < 40; ++search)
            {
              const Real middle = (lower + upper) / Real(2);
              if (thicknessRate(middle) > targetRate)
                lower = middle;
              else
                upper = middle;
            }
            multiplier = upper;
          }
        }
      }
      if (feasible)
      {
        xiRho = rhoGradient + nullSpaceMultiplier * volumeGradient +
          multiplier * thicknessRepair;
        nullSpaceMultiplier += multiplier * thicknessVolumeMultiplier;
      }
      else
      {
        xiRho = thicknessRepair;
        nullSpaceMultiplier = thicknessVolumeMultiplier;
      }
      stageDiagnostics.thicknessActiveMultiplier = multiplier;
      stageDiagnostics.thicknessActiveFeasible = feasible ? Real(1) : Real(0);
      Alert::Info() << substageHeading("Aggregate thickness correction") << Alert::NewLine
                    << diagnosticLabel("Target D guard P [theta]:")
                    << Alert::Notation::Number(targetRate) << Alert::NewLine
                    << diagnosticLabel("Selected multiplier:")
                    << (feasible ? std::to_string(multiplier) : "Repair-only fallback")
                    << Alert::NewLine
                    << diagnosticLabel("Target feasible in tested span:")
                    << (feasible ? "Yes" : "No") << Alert::Raise;
    }
    GridFunction xiRhoNorm(levelSetSpace);
    xiRhoNorm = Frobenius(xiRho);
    const Real xiRhoInfinityNorm = std::max(xiRhoNorm.max(), Real(1e-30));
    GridFunction theta(shapeSpace);
    theta = xiRho;
    theta /= xiRhoInfinityNorm;
    xiRhoNorm = Frobenius(theta);
    const Real thetaInfinityNorm = xiRhoNorm.max();
    const Real dRhoTheta = dRho(theta);
    const Real dVolumeTheta = dVolume(theta);
    if (thicknessFactor > 0 && thicknessMethod == "baseline")
      stageDiagnostics.dThicknessTheta = -thicknessLoad.dot(theta.getData());
    const Real requiredDVolumeTheta = 0;
    Alert::Info() << substageHeading("Constrained direction") << Alert::NewLine
                  << diagnosticLabel("Null-space multiplier:")
                  << Alert::Notation::Number(nullSpaceMultiplier) << Alert::NewLine
                  << diagnosticLabel("Xi rho infinity norm:")
                  << Alert::Notation::Number(xiRhoInfinityNorm) << Alert::NewLine
                  << diagnosticLabel("Theta infinity norm:")
                  << Alert::Notation::Number(thetaInfinityNorm) << Alert::NewLine
                  << diagnosticLabel("D rho [theta]:")
                  << Alert::Notation::Number(dRhoTheta) << Alert::NewLine
                  << diagnosticLabel("D volume [theta]:")
                  << Alert::Notation::Number(dVolumeTheta) << Alert::NewLine
                  << diagnosticLabel("Required D volume [theta]:")
                  << Alert::Notation::Number(requiredDVolumeTheta) << Alert::Raise;
    previousRho = rho;
    predictedRhoChange = dt * dRhoTheta;
    previousVolume = volume;
    predictedVolumeChange = dt * dVolumeTheta;
    GridFunction eikonalDistance(levelSetSpace);
    Distance::Eikonal(eikonalDistance)
      .setInterior(Obstacle)
      .setInterface(Gamma)
      .solve()
      .sign();
    TrialFunction periodicDistance(levelSetSpace);
    TestFunction distanceTest(levelSetSpace);
    Problem distanceProjection(periodicDistance, distanceTest);
    auto eikonalDistanceLoad = Integral(eikonalDistance, distanceTest);
    eikonalDistanceLoad.setOrder(2);
    // The projection only reconciles the traces on the cuts: Gamma is pinned,
    // so it cannot move the interface the transport starts from.
    distanceProjection = Integral(periodicDistance, distanceTest) - eikonalDistanceLoad +
      DirichletBC(periodicDistance, RealFunction(0)).on(Gamma);
    distanceProjection.assemble();
    shapeCoupling.assembleScalarTracePenalty(levelSetSpace,
      distanceProjection.getLinearSystem(), levelSetPenalty, FlatSet<Attribute>{Gamma});
    Solver::CG(distanceProjection).solve();
    GridFunction distance(levelSetSpace);
    distance = periodicDistance.getSolution();
    const auto& distanceSystem = distanceProjection.getLinearSystem();
    const Real distanceProjectionResidual =
      (distanceSystem.getOperator() * distanceSystem.getSolution() -
        distanceSystem.getVector())
        .norm() /
      std::max(distanceSystem.getVector().norm(), Real(1));
    const Real distanceProjectionCorrection =
      (distance.getData() - eikonalDistance.getData()).lpNorm<Eigen::Infinity>();
    // The Eikonal distance vanishes on Gamma, so the projected value at an
    // interface vertex is how far the projection moves the interface there.
    std::vector<char> onInterface(mesh.getVertexCount(), 0);
    for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
    {
      if (face->getAttribute() == Gamma)
        for (const Index vertex : face->getVertices())
          onInterface[vertex] = 1;
    }
    Real interfaceShiftMaximum = 0;
    Real interfaceShiftSquares = 0;
    size_t interfaceVertices = 0;
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
    {
      if (!onInterface[vertex])
        continue;
      const Real shift = std::abs(distance[vertex]);
      interfaceShiftMaximum = std::max(interfaceShiftMaximum, shift);
      interfaceShiftSquares += shift * shift;
      ++interfaceVertices;
    }
    const Real interfaceShiftRMS = interfaceVertices
      ? std::sqrt(interfaceShiftSquares / static_cast<Real>(interfaceVertices))
      : Real(0);
    stageDiagnostics.eikonalJump = shapeCoupling.scalarJump(eikonalDistance);
    stageDiagnostics.projectedJump = shapeCoupling.scalarJump(distance);
    stageDiagnostics.distanceCorrection = distanceProjectionCorrection;
    stageDiagnostics.interfaceShift = interfaceShiftMaximum;
    Alert::Info() << substageHeading("Distance trace projection") << Alert::NewLine
                  << diagnosticLabel("Linear residual:")
                  << Alert::Notation::Number(distanceProjectionResidual) << Alert::NewLine
                  << diagnosticLabel("Eikonal rotated jump:")
                  << Alert::Notation::Number(stageDiagnostics.eikonalJump)
                  << Alert::NewLine << diagnosticLabel("Projected rotated jump:")
                  << Alert::Notation::Number(stageDiagnostics.projectedJump)
                  << Alert::NewLine << diagnosticLabel("Infinity correction:")
                  << Alert::Notation::Number(distanceProjectionCorrection)
                  << Alert::NewLine << diagnosticLabel("Interface shift maximum:")
                  << Alert::Notation::Number(interfaceShiftMaximum) << " = "
                  << Alert::Notation::Number(interfaceShiftMaximum / h) << " h"
                  << Alert::NewLine << diagnosticLabel("Interface shift RMS:")
                  << Alert::Notation::Number(interfaceShiftRMS) << " = "
                  << Alert::Notation::Number(interfaceShiftRMS / h) << " h"
                  << Alert::Raise;

    stageSeconds[5] = elapsedSeconds(stage6Start);
    reportStageTiming(6, stageSeconds[5]);
    const auto stage7Start = Clock::now();
    announce("Stage 7: Writing the chamber and sewn fields.");
    chamber.clear();
    chamber.add("Distance", distance, IO::XDMF::Center::Node);
    chamber.add("Theta", theta, IO::XDMF::Center::Node);
    fluidState.clear();
    fluidState.setMesh(fluid, IO::XDMF::MeshPolicy::Transient);
    fluidState.add("Translation_0", uT0, IO::XDMF::Center::Node);
    fluidState.add("Translation_1", uT1, IO::XDMF::Center::Node);
    fluidState.add("Translation_2", uT2, IO::XDMF::Center::Node);
    fluidState.add("Rotation_0", uR0, IO::XDMF::Center::Node);
    fluidState.add("Rotation_1", uR1, IO::XDMF::Center::Node);
    fluidState.add("Rotation_2", uR2, IO::XDMF::Center::Node);
    fluidState.add("Pressure_Translation_0", pT0, IO::XDMF::Center::Node);
    fluidState.add("Pressure_Translation_1", pT1, IO::XDMF::Center::Node);
    fluidState.add("Pressure_Translation_2", pT2, IO::XDMF::Center::Node);
    fluidState.add("Pressure_Rotation_0", pR0, IO::XDMF::Center::Node);
    fluidState.add("Pressure_Rotation_1", pR1, IO::XDMF::Center::Node);
    fluidState.add("Pressure_Rotation_2", pR2, IO::XDMF::Center::Node);

    KelvinBall::SewedOutput sewedDesign(mesh, FlatSet<Attribute>{Gamma, Outer});
    P1 sewedDesignScalar(sewedDesign.getMesh());
    P1 sewedDesignVector(sewedDesign.getMesh(), 3);
    GridFunction sewedDistance(sewedDesignScalar);
    GridFunction sewedVelocity(sewedDesignVector);
    sewedDesign.setScalar(sewedDistance, distance);
    sewedDesign.setVector(sewedVelocity, theta);
    sewedDesignOutput.clear();
    sewedDesignOutput.setMesh(sewedDesign.getMesh(), IO::XDMF::MeshPolicy::Transient);
    sewedDesignOutput.add("Distance", sewedDistance, IO::XDMF::Center::Node);
    sewedDesignOutput.add("Theta", sewedVelocity, IO::XDMF::Center::Node);

    KelvinBall::SewedOutput sewedFluid(fluid, FlatSet<Attribute>{Gamma, Outer});
    VelocitySpace sewedVelocitySpace = makeVelocitySpace(sewedFluid.getMesh());
    PressureSpace sewedPressureSpace(sewedFluid.getMesh());
    GridFunction sewedUT0(sewedVelocitySpace), sewedUT1(sewedVelocitySpace),
      sewedUT2(sewedVelocitySpace), sewedUR0(sewedVelocitySpace),
      sewedUR1(sewedVelocitySpace), sewedUR2(sewedVelocitySpace);
    GridFunction sewedPT0(sewedPressureSpace), sewedPT1(sewedPressureSpace),
      sewedPT2(sewedPressureSpace), sewedPR0(sewedPressureSpace),
      sewedPR1(sewedPressureSpace), sewedPR2(sewedPressureSpace);
    const auto translations = std::array{&uT0, &uT1, &uT2};
    const auto rotations = std::array{&uR0, &uR1, &uR2};
    const auto translationPressures = std::array{&pT0, &pT1, &pT2};
    const auto rotationPressures = std::array{&pR0, &pR1, &pR2};
    sewedFluid.setVectorLoad(sewedUT0, translations, 0);
    sewedFluid.setVectorLoad(sewedUT1, translations, 1);
    sewedFluid.setVectorLoad(sewedUT2, translations, 2);
    sewedFluid.setVectorLoad(sewedUR0, rotations, 0);
    sewedFluid.setVectorLoad(sewedUR1, rotations, 1);
    sewedFluid.setVectorLoad(sewedUR2, rotations, 2);
    sewedFluid.setScalarLoad(sewedPT0, translationPressures, 0);
    sewedFluid.setScalarLoad(sewedPT1, translationPressures, 1);
    sewedFluid.setScalarLoad(sewedPT2, translationPressures, 2);
    sewedFluid.setScalarLoad(sewedPR0, rotationPressures, 0);
    sewedFluid.setScalarLoad(sewedPR1, rotationPressures, 1);
    sewedFluid.setScalarLoad(sewedPR2, rotationPressures, 2);
    sewedFluidOutput.clear();
    sewedFluidOutput.setMesh(sewedFluid.getMesh(), IO::XDMF::MeshPolicy::Transient);
    sewedFluidOutput.add("Translation_0", sewedUT0, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Translation_1", sewedUT1, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Translation_2", sewedUT2, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Rotation_0", sewedUR0, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Rotation_1", sewedUR1, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Rotation_2", sewedUR2, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Pressure_Translation_0", sewedPT0, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Pressure_Translation_1", sewedPT1, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Pressure_Translation_2", sewedPT2, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Pressure_Rotation_0", sewedPR0, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Pressure_Rotation_1", sewedPR1, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Pressure_Rotation_2", sewedPR2, IO::XDMF::Center::Node);

    if (motionEvery > 0 &&
      ((iteration + 1) % motionEvery == 0 || iteration + 1 == maxIterations))
    {
      // The body is free and pushed by F without torque, so R [Z; w] = [F; 0]
      // with K = kI, C = cI, Q = qI gives a screw motion about F.
      const Real determinant = k * q - c * c;
      const Math::SpatialVector<Real> translation = (q / determinant) * motionForce;
      const Math::SpatialVector<Real> angular = (-c / determinant) * motionForce;
      const Real angularSpeed = angular.norm();
      const Real period = angularSpeed > 0 ? 2 * M_PI / angularSpeed : Real(1);
      Alert::Info() << substageHeading("Rigid motion") << Alert::NewLine
                    << diagnosticLabel("Force:") << vectorText(motionForce)
                    << Alert::NewLine << diagnosticLabel("Translational velocity Z:")
                    << vectorText(translation) << Alert::NewLine
                    << diagnosticLabel("Angular velocity omega:") << vectorText(angular)
                    << Alert::NewLine << diagnosticLabel("Period of one revolution:")
                    << Alert::Notation::Number(period) << Alert::NewLine
                    << diagnosticLabel("Pitch (travel per revolution):")
                    << Alert::Notation::Number(translation.norm() * period)
                    << Alert::Raise;

      // Velocity of the fluid for that motion, by linearity of the states.
      std::vector<Math::SpatialVector<Real>> motionVelocity(
        sewedFluid.getMesh().getVertexCount(), Math::SpatialVector<Real>::Zero(3));
      const auto sewedTranslations = std::array{&sewedUT0, &sewedUT1, &sewedUT2};
      const auto sewedRotations = std::array{&sewedUR0, &sewedUR1, &sewedUR2};
      for (Index vertex = 0; vertex < sewedFluid.getMesh().getVertexCount(); ++vertex)
      {
        const auto dofs = sewedVelocitySpace.getDOFs(0, vertex);
        for (size_t load = 0; load < 3; ++load)
          for (Eigen::Index component = 0; component < 3; ++component)
            motionVelocity[vertex](component) +=
              translation(load) * sewedTranslations[load]->getData()(dofs(component)) +
              angular(load) * sewedRotations[load]->getData()(dofs(component));
      }

      char name[64];
      std::snprintf(name, sizeof(name), "KelvinBallMotion-%06zu", iteration + 1);
      IO::XDMF motionXdmf(name);
      auto bodyOutput = motionXdmf.grid("Body");
      auto fluidOutput = motionXdmf.grid("Fluid");
      const Math::SpatialVector<Real> axis = angularSpeed > 0
        ? Math::SpatialVector<Real>(angular / angularSpeed)
        : Math::SpatialVector<Real>(motionForce / motionForce.norm());
      for (size_t frame = 0; frame < motionFrames; ++frame)
      {
        const Real time =
          period * static_cast<Real>(frame) / static_cast<Real>(motionFrames);
        const Math::SpatialMatrix<Real> rotation = Eigen::AngleAxis<Real>(
          angularSpeed * time, Eigen::Matrix<Real, 3, 1>(axis(0), axis(1), axis(2)))
                                                     .toRotationMatrix();
        const Math::SpatialVector<Real> shift = time * translation;
        LocalMesh body = sewedDesign.getMesh();
        for (Index vertex = 0; vertex < body.getVertexCount(); ++vertex)
          body.setVertexCoordinates(
            vertex, rotation * body.getVertexCoordinates(vertex) + shift);
        LocalMesh moving = sewedFluid.getMesh();
        for (Index vertex = 0; vertex < moving.getVertexCount(); ++vertex)
          moving.setVertexCoordinates(
            vertex, rotation * moving.getVertexCoordinates(vertex) + shift);
        VelocitySpace movingSpace = makeVelocitySpace(moving);
        GridFunction velocity(movingSpace);
        for (Index vertex = 0; vertex < moving.getVertexCount(); ++vertex)
        {
          const auto dofs = movingSpace.getDOFs(0, vertex);
          const Math::SpatialVector<Real> value = rotation * motionVelocity[vertex];
          for (Eigen::Index component = 0; component < 3; ++component)
            velocity.getData()(dofs(component)) = value(component);
        }
        bodyOutput.clear();
        bodyOutput.setMesh(body, IO::XDMF::MeshPolicy::Transient);
        fluidOutput.clear();
        fluidOutput.setMesh(moving, IO::XDMF::MeshPolicy::Transient);
        fluidOutput.add("Velocity", velocity, IO::XDMF::Center::Node);
        // The period is long while the coupling is weak, so the series is
        // indexed by the fraction of a revolution rather than by time.
        motionXdmf.write(time / period).flush();
      }
    }

    if (iteration + 1 == maxIterations)
    {
      xdmf.write(static_cast<Real>(iteration)).flush();
      sewedXdmf.write(static_cast<Real>(iteration)).flush();
      stageSeconds[6] = elapsedSeconds(stage7Start);
      reportStageTiming(7, stageSeconds[6]);
      writeHistory(nullSpaceMultiplier, xiRhoInfinityNorm, thetaInfinityNorm, dRhoTheta,
        dVolumeTheta, requiredDVolumeTheta, rhoGradientDiagnostics,
        volumeGradientDiagnostics);
      continue;
    }

    stageSeconds[6] = elapsedSeconds(stage7Start);
    reportStageTiming(7, stageSeconds[6]);
    const auto stage8Start = Clock::now();
    announce("Stage 8: Advecting the level set.");
    MMG::Mesh& advectionMesh = wngirBackground ? *wngirBackground : mesh;
    auto& advectionConnectivity = advectionMesh.getConnectivity();
    advectionConnectivity.discover(3, 2);
    advectionConnectivity.discover(3, 1);
    advectionConnectivity.restrict(1, 0);
    advectionConnectivity.restrict(2, 0);
    advectionConnectivity.restrict(2, 3);
    advectionConnectivity.discover(0, 0);
    P1 advectionLevelSetSpace(advectionMesh);
    P1 advectionShapeSpace(advectionMesh, 3);
    GridFunction advectionDistance(advectionLevelSetSpace);
    GridFunction advectionDirection(advectionShapeSpace);
    if (wngirBackground)
    {
      // The fitted mesh is the background with its vertices moved: the same
      // numbering at different positions. Copying nodal values would carry
      // each one from x + u(x) back to x and undo the fit near the interface,
      // so both fields are evaluated at the positions of the background.
      const Location::AABB<MMG::Mesh> fittedLocator(mesh);
      std::size_t unlocated = 0;
      advectionDistance = RealFunction([&](const Geometry::Point& point) {
        const auto located = fittedLocator.locate(3, point.getPhysicalCoordinates());
        if (!located)
        {
          ++unlocated;
          return Real(0);
        }
        return distance.getValue(*located);
      });
      advectionDirection =
        VectorFunction(static_cast<size_t>(3), [&](const Geometry::Point& point) {
          Math::SpatialVector<Real> value(3);
          value.setZero();
          const auto located = fittedLocator.locate(3, point.getPhysicalCoordinates());
          if (!located)
          {
            ++unlocated;
            return value;
          }
          const auto sample = theta.getValue(*located);
          for (Eigen::Index component = 0; component < 3; ++component)
            value(component) = sample(component);
          return value;
        });
      if (unlocated > 0)
        throw std::runtime_error(
          "Background points fell outside the fitted mesh during the transfer.");
      Alert::Info()
        << substageHeading("Fitted-to-background transfer") << Alert::NewLine
        << diagnosticLabel("Nodal-copy error (distance):")
        << Alert::Notation::Number(
             (advectionDistance.getData() - distance.getData()).lpNorm<Eigen::Infinity>())
        << Alert::NewLine << diagnosticLabel("Nodal-copy error (direction):")
        << Alert::Notation::Number(
             (advectionDirection.getData() - theta.getData()).lpNorm<Eigen::Infinity>())
        << Alert::Raise;
    }
    else
    {
      advectionDistance.getData() = distance.getData();
      advectionDirection.getData() = theta.getData();
    }
    KelvinBall::RotatedNitscheIntegrator advectionCoupling(advectionMesh,
      FlatSet<Attribute>{SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus},
      rotatedTracePhysicalTolerance, rotatedTraceReferenceTolerance);
    TrialFunction advected(advectionLevelSetSpace);
    TestFunction test(advectionLevelSetSpace);
    KelvinBall::RotatedCharacteristicContinuation rotationalContinuation(
      -dt, advectionMesh, advectionCoupling.getLocator(), RotationPairs);
    Problem transport(advected, test);
    auto transportedDistance =
      Integral(Flow(-dt, advectionDistance, advectionDirection, Math::RungeKutta::RK4{},
                 rotationalContinuation),
        test);
    transportedDistance.setOrder(advectionQuadratureOrder);
    transport = Integral(advected, test) - transportedDistance;
    transport.assemble();
    advectionCoupling.assembleScalarTracePenalty(advectionLevelSetSpace,
      transport.getLinearSystem(), levelSetPenalty, advectionDistance);
    Solver::CG(transport).solve();
    const auto& advectedDistance = advected.getSolution();
    const auto& transportSystem = transport.getLinearSystem();
    const Real transportResidual =
      (transportSystem.getOperator() * transportSystem.getSolution() -
        transportSystem.getVector())
        .norm() /
      std::max(transportSystem.getVector().norm(), Real(1));
    const Real advectionIncrement =
      (advectedDistance.getData() - advectionDistance.getData())
        .lpNorm<Eigen::Infinity>();
    const Real distanceJump = advectionCoupling.scalarJump(advectionDistance);
    const Real advectedJump = advectionCoupling.scalarJump(advectedDistance);
    stageDiagnostics.advectionIncrement = advectionIncrement;
    stageDiagnostics.advectedJump = advectedJump;
    Alert::Info() << substageHeading("Advected distance") << Alert::NewLine
                  << diagnosticLabel("Advection step:") << Alert::Notation::Number(dt)
                  << Alert::NewLine << diagnosticLabel("Linear residual:")
                  << Alert::Notation::Number(transportResidual) << Alert::NewLine
                  << diagnosticLabel("Minimum:")
                  << Alert::Notation::Number(advectedDistance.min()) << Alert::NewLine
                  << diagnosticLabel("Maximum:")
                  << Alert::Notation::Number(advectedDistance.max()) << Alert::NewLine
                  << diagnosticLabel("Increment infinity norm:")
                  << Alert::Notation::Number(advectionIncrement) << Alert::NewLine
                  << diagnosticLabel("Distance rotated jump:")
                  << Alert::Notation::Number(distanceJump) << Alert::NewLine
                  << diagnosticLabel("Advected rotated jump:")
                  << Alert::Notation::Number(advectedJump) << Alert::Raise;
    GridFunction advectedOutput(levelSetSpace);
    advectedOutput.getData() = advectedDistance.getData();
    chamber.add("Advected", advectedOutput, IO::XDMF::Center::Node);
    GridFunction sewedAdvected(sewedDesignScalar);
    sewedDesign.setScalar(sewedAdvected, advectedOutput);
    sewedDesignOutput.add("Advected", sewedAdvected, IO::XDMF::Center::Node);

    xdmf.write(static_cast<Real>(iteration)).flush();
    sewedXdmf.write(static_cast<Real>(iteration)).flush();
    stageSeconds[7] = elapsedSeconds(stage8Start);
    reportStageTiming(8, stageSeconds[7]);

    const auto stage9Start = Clock::now();
    announce(reconstructionMethod == "wngir"
        ? "Stage 9: Fitting the advected interface with WNGIR."
        : "Stage 9: Reconstructing the advected interface with MMG.");
    MMGReconstruction result = [&]() {
      if (reconstructionMethod == "wngir")
      {
        MMG::Mesh classified =
          classifyLevelSetForWNGIR(*wngirBackground, advectedDistance);
        P1 classifiedLevelSetSpace(classified);
        GridFunction classifiedLevelSet(classifiedLevelSetSpace);
        classifiedLevelSet.getData() = advectedDistance.getData();
        return fitLevelSetWNGIR(
          classified, classifiedLevelSet, h, outerRadius, argc, argv);
      }
      // A failed MMG stage is retried on the same advected level set with
      // every MMG size computed from half the scale; the next iteration
      // starts again from h.
      Real scale = h;
      for (size_t attempt = 0;; ++attempt)
      {
        try
        {
          auto reconstructed = discretizeLevelSetMMG(mesh, advectedDistance, scale,
            sphere, configuration.adapt, configuration.mmgSnap);
          reconstructed.diagnostics.scale = scale / h;
          return reconstructed;
        }
        catch (const std::exception& error)
        {
          // Keep what MMG received, so that the failure can be reproduced.
          const std::string failure = "mmg-failure-" + std::to_string(iteration + 1) +
            "-" + std::to_string(attempt);
          mesh.save(failure + ".mesh", IO::FileFormat::MEDIT);
          advectedDistance.save(failure + ".sol", IO::FileFormat::MEDIT);
          Alert::Warning() << "Saved the failed MMG input to " << failure << ".mesh and "
                           << failure << ".sol." << Alert::Raise;
          if (attempt == configuration.mmgRetries)
            throw;
          Alert::Warning() << "MMG reconstruction failed at scale "
                           << Alert::Notation::Number(scale / h) << " h: " << error.what()
                           << Alert::NewLine << "Retrying at scale "
                           << Alert::Notation::Number(scale / (2 * h)) << " h."
                           << Alert::Raise;
          scale /= 2;
        }
      }
    }();
    reconstruction = result.diagnostics;
    reconstructionOutput.clear();
    reconstructionOutput.setMesh(result.mesh, IO::XDMF::MeshPolicy::Transient);
    reconstructionXdmf.write(static_cast<Real>(iteration + 1)).flush();
    checkFixedGeometry(result.mesh, outerRadius);
    checkMaterials(result.mesh);
    nextMesh.emplace(std::move(result.mesh));
    stageSeconds[8] = elapsedSeconds(stage9Start);
    reportStageTiming(9, stageSeconds[8]);
    writeHistory(nullSpaceMultiplier, xiRhoInfinityNorm, thetaInfinityNorm, dRhoTheta,
      dVolumeTheta, requiredDVolumeTheta, rhoGradientDiagnostics,
      volumeGradientDiagnostics);
  }

  Alert::Success() << "Wrote KelvinBall.xdmf, KelvinBallSewed.xdmf, "
                   << reconstructionName
                   << ".xdmf, and "
                      "kelvin-ball.csv"
                   << Alert::Raise;
  return 0;
}

int KelvinBall::KelvinBallOptimization::run()
{
  return Implementation(m_argc, m_argv).run();
}
