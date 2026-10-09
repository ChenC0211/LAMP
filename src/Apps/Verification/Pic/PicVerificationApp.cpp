#include "Verification/Pic/PicVerificationApp.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <vector>

#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Field/CellField.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/IO/OutputLayout.h"
#include "Verification/Pic/DictionaryPicVerification.h"
#include "Solvers/Poisson/EdgeCutLengthHelper.h"
#include "Modules/Plasma/ExplicitPIC/Geometry/LevelSetQuery.h"
#include "Solvers/Poisson/CutCellMetrics.h"
#include "Solvers/Poisson/PoissonFVSolverCellCentered.h"
#include "Solvers/Poisson/PolarizationPoissonFVSolverCellCentered.h"
#include "Solvers/SolversConfig.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Verification/Pic/RegressionHarness.h"
#include "Modules/Plasma/ExplicitPIC/Field/ElectricFieldEvaluator.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/DepositQK.h"
#include "Modules/Plasma/ExplicitPIC/Field/MatchedGather.h"
#include "Modules/Plasma/ExplicitPIC/Diagnostics/PicDiagnostics.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Core/Units/UnitsContext.h"

using namespace std;

namespace
{
using Solvers::Poisson::CutCellMetrics;

static inline std::string lowerCopy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

static bool isLeafCell(const Cell& c)
{
    return c.son[0] == -1;
}

static inline double minCellSize(const Cell& c)
{
    const Vector2D s = c.size();
    return std::min(s.x, s.y);
}

// Coarse-fine composite band predicate (mirrors CouplingWeightsCache.cpp).
static bool isCompositeBandCell(const QuadtreeMesh& mesh, const Cell& c)
{
    // Non-leaf treated as in-band.
    if (!(c.son[0] == -1 && c.son[1] == -1 && c.son[2] == -1 && c.son[3] == -1))
        return true;

    const int lvl = c.level();
    for (int dir = 0; dir < 4; ++dir)
    {
        const int nb0 = c.neiCell[dir][0];
        const int nb1 = c.neiCell[dir][1];
        if (nb1 != -1)
            return true;
        if (nb0 != -1)
        {
            const Cell& nb = mesh.cells(nb0);
            const bool nbLeaf = (nb.son[0] == -1 && nb.son[1] == -1 && nb.son[2] == -1 && nb.son[3] == -1);
            if (!nbLeaf || nb.level() != lvl)
                return true;
        }
    }
    return false;
}

struct BandCellBox
{
    Vector2D min;
    Vector2D max;
};

static std::vector<BandCellBox> collectCompositeBandCells(const QuadtreeMesh& mesh)
{
    std::vector<BandCellBox> out;
    out.reserve(static_cast<size_t>(mesh.cellIndexNum()));
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        if (!isLeafCell(c))
            continue;
        if (!isCompositeBandCell(mesh, c))
            continue;
        const Vector2D cc = c.centerCoord();
        const Vector2D sz = c.size();
        BandCellBox b;
        b.min = Vector2D(cc.x - 0.5 * sz.x, cc.y - 0.5 * sz.y);
        b.max = Vector2D(cc.x + 0.5 * sz.x, cc.y + 0.5 * sz.y);
        out.push_back(b);
    }
    return out;
}

static inline double distancePointToBox(const Vector2D& p, const BandCellBox& b)
{
    const double dx = (p.x < b.min.x) ? (b.min.x - p.x) : (p.x > b.max.x ? p.x - b.max.x : 0.0);
    const double dy = (p.y < b.min.y) ? (b.min.y - p.y) : (p.y > b.max.y ? p.y - b.max.y : 0.0);
    return std::sqrt(dx * dx + dy * dy);
}

static inline double wrapCoord(const double x, const double xmin, const double xmax)
{
    const double L = xmax - xmin;
    if (!(L > 0.0))
        return x;
    double v = x;
    while (v < xmin)
        v += L;
    while (v >= xmax)
        v -= L;
    return v;
}

static inline Vector2D wrapPeriodic(const Vector2D& p, const QuadtreeMesh& mesh, bool px, bool py)
{
    Vector2D out = p;
    if (px)
        out.x = wrapCoord(out.x, mesh.minCoord.x, mesh.maxCoord.x);
    if (py)
        out.y = wrapCoord(out.y, mesh.minCoord.y, mesh.maxCoord.y);
    return out;
}

static void fillCircleLevelSet(QuadtreeMesh& mesh,
                               ImmersedGeometry<QuadtreeMesh>& iGeometry,
                               int levelSetI,
                               const Vector2D& center,
                               double radius)
{
    if (levelSetI < 0 || levelSetI >= (int)iGeometry.levelSets.size())
        return;

    auto& psi = iGeometry.levelSets.at(levelSetI);
    auto& psiDx = iGeometry.LS_dx(levelSetI);
    auto& psiDy = iGeometry.LS_dy(levelSetI);

    if (psiDx.size() != mesh.vertIndexNum())
        psiDx.initialize(mesh, 0.0);
    if (psiDy.size() != mesh.vertIndexNum())
        psiDy.initialize(mesh, 0.0);

    const double eps = 1e-14;
    for (int vi = 0; vi < mesh.vertIndexNum(); ++vi)
    {
        const int vLbl = mesh.vertIndexToLabel[vi];
        const Vert& v = mesh.verts(vLbl);

        const Vector2D xc = v.coord();
        const double dx = xc.x - center.x;
        const double dy = xc.y - center.y;
        const double r = std::sqrt(dx * dx + dy * dy);
        const double invr = 1.0 / std::max(r, eps);

        psi.value[v.index()] = r - radius;
        psiDx.value[v.index()] = dx * invr;
        psiDy.value[v.index()] = dy * invr;
    }
}

static inline bool isInsidePsi(const double psi, const bool insideIsNegative)
{
    return insideIsNegative ? (psi <= 0.0) : (psi >= 0.0);
}

struct EmbeddedContext
{
    std::unique_ptr<DictionaryGeometry> dict;
    std::unique_ptr<ImmersedGeometry<QuadtreeMesh>> geom;
    Solvers::Poisson::CutCellMetrics metrics;
    std::vector<double> alphaByCellIndex;
    int levelSetI = 0;
    bool insideIsNegative = true;
};

static double computeInsideAreaFromAlpha(const QuadtreeMesh& mesh, const CutCellMetrics& metrics)
{
    double area = 0.0;
    const auto& cm = metrics.cellMetricsByIndex();
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        if (!isLeafCell(c))
            continue;
        area += cm[ci].alpha * c.area();
    }
    return area;
}

static double computeInterfacePerimeterFromAgamma(const QuadtreeMesh& mesh, const CutCellMetrics& metrics)
{
    double per = 0.0;
    const auto& cm = metrics.cellMetricsByIndex();
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        if (!isLeafCell(c))
            continue;
        per += cm[ci].A_gamma;
    }
    return per;
}

static void dumpPhiCSV(const std::string& filePath,
                       const QuadtreeMesh& mesh,
                       const CutCellMetrics& metrics,
                       const CellField<double>& phi,
                       const CellField<double>& phiExact)
{
    std::ofstream out(filePath);
    if (!out)
        return;

    out << "cellIndex,cellLabel,x,y,alpha,phi,phiExact,err\n";
    const auto& cm = metrics.cellMetricsByIndex();

    out << std::setprecision(16);
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        if (!isLeafCell(c))
            continue;
        const double alpha = cm[ci].alpha;
        if (alpha <= 0.0)
            continue;

        const Vector2D cc = c.centerCoord();
        const double x = cc.x;
        const double y = cc.y;
        const double p = phi.value[ci];
        const double pe = phiExact.value[ci];
        out << ci << ',' << cLbl << ',' << x << ',' << y << ',' << alpha << ',' << p << ',' << pe << ',' << (p - pe) << "\n";
    }
}

static void dumpPhiRhoCSV(const std::string& filePath,
                          const QuadtreeMesh& mesh,
                          const CutCellMetrics& metrics,
                          const CellField<double>& phi,
                          const CellField<double>& rho)
{
    std::ofstream out(filePath);
    if (!out)
        return;

    out << "cellIndex,cellLabel,x,y,alpha,phi,rho\n";
    const auto& cm = metrics.cellMetricsByIndex();

    out << std::setprecision(16);
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        if (!isLeafCell(c))
            continue;
        const double alpha = cm[ci].alpha;
        if (alpha <= 0.0)
            continue;

        const Vector2D cc = c.centerCoord();
        out << ci << ',' << cLbl << ',' << cc.x << ',' << cc.y << ','
            << alpha << ',' << phi.value[ci] << ',' << rho.value[ci] << "\n";
    }
}

static void dumpEnCSV(const std::string& filePath,
                      const QuadtreeMesh& mesh,
                      const CutCellMetrics& metrics,
                      const std::vector<double>& EnGamma)
{
    std::ofstream out(filePath);
    if (!out)
        return;

    out << "cellIndex,cellLabel,x_gamma,y_gamma,nx,ny,En,A_gamma,d_gamma\n";
    const auto& cm = metrics.cellMetricsByIndex();

    out << std::setprecision(16);
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        if (!metrics.isActiveCellIndex(ci))
            continue;
        const auto& m = cm[ci];
        if (!(m.A_gamma > 0.0))
            continue;

        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        const Vector2D cc = c.centerCoord();
        const Vector2D n = m.n_gamma;
        const double d = m.d_gamma;
        const Vector2D xg = Vector2D(cc.x - n.x * d, cc.y - n.y * d);
        const double En = (ci < static_cast<int>(EnGamma.size())) ? EnGamma[ci] : 0.0;

        out << ci << ',' << cLbl << ',' << xg.x << ',' << xg.y << ','
            << n.x << ',' << n.y << ',' << En << ','
            << m.A_gamma << ',' << d << "\n";
    }
}

// -----------------------------------------------------------------------------
// ElectricFieldEvaluator regression: analytic phi with y-periodic boundaries.
// -----------------------------------------------------------------------------
static inline double phiSinCosPeriodicY(const double x, const double y)
{
    // Unit-box analytic field:
    //   phi = sin(pi x) cos(2 pi y)
    const double pi = std::acos(-1.0);
    return std::sin(pi * x) * std::cos(2.0 * pi * y);
}

static inline Vector2D gradPhiSinCosPeriodicY(const double x, const double y)
{
    const double pi = std::acos(-1.0);
    const double dphidx = pi * std::cos(pi * x) * std::cos(2.0 * pi * y);
    const double dphidy = -2.0 * pi * std::sin(pi * x) * std::sin(2.0 * pi * y);
    return Vector2D(dphidx, dphidy);
}

static inline double exactEnAtFaceCenter(const FaceSegment& seg)
{
    const Vector2D g = gradPhiSinCosPeriodicY(seg.faceCenter.x, seg.faceCenter.y);
    const double dphi_dn = g.x * seg.n.x + g.y * seg.n.y;
    return -dphi_dn;
}

static inline void accumWeightedError(const double err, const double w,
                                      double& sumW, double& sumErr2, double& linf)
{
    sumW += w;
    sumErr2 += w * err * err;
    linf = std::max(linf, std::abs(err));
}

static inline Vector2D periodicShiftForSide(const QuadtreeMesh& mesh, const FaceSide side)
{
    const double Lx = mesh.maxCoord.x - mesh.minCoord.x;
    const double Ly = mesh.maxCoord.y - mesh.minCoord.y;
    switch (side)
    {
    case FaceSide::West:  return Vector2D(+Lx, 0.0);
    case FaceSide::East:  return Vector2D(-Lx, 0.0);
    case FaceSide::South: return Vector2D(0.0, +Ly);
    case FaceSide::North: return Vector2D(0.0, -Ly);
    default:              return Vector2D(0.0, 0.0);
    }
}

static inline int findPeriodicMateCellLabel(const QuadtreeMesh& mesh,
                                            const FaceSegment& seg,
                                            const Vector2D& shift,
                                            const double probeEps)
{
    const Vector2D probe = seg.faceCenter + shift + seg.n * probeEps;
    return mesh.findCellLabel(probe);
}

} // namespace

PicVerificationApp::PicVerificationApp(const RunConfig& cfg)
    : cfg_(cfg)
{
}

int PicVerificationApp::run()
{
    OutputLayout::fields();
    OutputLayout::meta();

    DictionaryPicVerification picDict;
    std::string dictPath = "./system/pic.dict";
    const bool hasPicDict = picDict.readDict(dictPath);
    if (!hasPicDict)
    {
        std::cerr << "[PicVerification] ERROR: missing " << dictPath << std::endl;
        return 2;
    }

    const auto& cfg = picDict.physics();
    const auto& verification = picDict.config();

    RegressionHarness harness(verification.regressionHarness);
    RegressionHarness::Metrics regMetrics;

    const bool doUnitBox = verification.regressionUnitBoxSinSin;
    const bool doUnitBoxPeriodicY = verification.regressionUnitBoxSinCosPeriodicY;
    const bool doEFieldPeriodicY = verification.regressionElectricFieldSinCosPeriodicY;
    const bool doCutCircleMetrics = verification.regressionCutCircle.enabled;
    const bool doCutCirclePoisson = verification.regressionCutCirclePoisson.enabled;
    const bool doCutCirclePolPoisson = verification.regressionCutCirclePolarizationPoisson.enabled;
    const bool doEmbeddedConductor = verification.regressionEmbeddedConductor.enabled;
    const bool doSelfForceTest3 = verification.selfForceTest3.enabled;

    int requiredLevelSetCount = 1;
    if (doCutCircleMetrics)
        requiredLevelSetCount = std::max(requiredLevelSetCount, verification.regressionCutCircle.levelSetIndex + 1);
    if (doCutCirclePoisson)
        requiredLevelSetCount = std::max(requiredLevelSetCount, verification.regressionCutCirclePoisson.levelSetIndex + 1);
    if (doCutCirclePolPoisson)
        requiredLevelSetCount = std::max(requiredLevelSetCount, verification.regressionCutCirclePolarizationPoisson.levelSetIndex + 1);
    if (doEmbeddedConductor)
        requiredLevelSetCount = std::max(requiredLevelSetCount, verification.regressionEmbeddedConductor.levelSetIndex + 1);

    if (!doUnitBox && !doUnitBoxPeriodicY && !doEFieldPeriodicY && !doCutCircleMetrics && !doCutCirclePoisson && !doCutCirclePolPoisson &&
        !doEmbeddedConductor && !doSelfForceTest3)
    {
        std::cout << "[PicVerification] No verification branch enabled in pic.dict." << std::endl;
        std::cout << "  Enable one of: RegressionUnitBoxSinSin / RegressionUnitBoxSinCosPeriodicY / RegressionElectricFieldSinCosPeriodicY / RegressionCutCircle / RegressionCutCirclePoisson / RegressionCutCirclePolarizationPoisson / RegressionEmbeddedConductor / SelfForceTest3" << std::endl;
        return 0;
    }

    // ---------------------------------------------------------------------
    // Common mesh
    // ---------------------------------------------------------------------
    DictionaryMesh dictMesh;
    dictMesh.readDict();
    dictMesh.ensureLevelSetCount(requiredLevelSetCount);
    QuadtreeMesh mesh(dictMesh);

    // Deterministic meta-metrics (stable across platforms, useful for smoke/regress).
    regMetrics["meta.cellIndexNum"] = static_cast<double>(mesh.cellIndexNum());
    regMetrics["meta.vertIndexNum"] = static_cast<double>(mesh.vertIndexNum());

    // Unit-box sin-sin Poisson regression (no EB)
    if (doUnitBox)
    {
        std::cout << "[PicVerification] Regression: unit-box sin-sin (Dirichlet), FV Poisson." << std::endl;
        regMetrics["unitbox.enabled"] = 1.0;
        FaceSegmentGraph graph;
        graph.build(mesh);

        PoissonFVSolverCellCentered solver;
        PoissonFVSolverCellCentered::Options opt;
        if (cfg_.solvers)
            opt = cfg_.solvers->poissonCellCentered;

        CellField<double> phi;
        const auto norms = solver.regressionUnitBoxSinSin(mesh, graph, phi, opt);
        regMetrics["unitbox.L2"] = norms.L2;
        regMetrics["unitbox.Linf"] = norms.Linf;
    }

    // Unit-box regression with y-periodic wrap (exercises Periodic BC assembly)
    if (doUnitBoxPeriodicY)
    {
        std::cout << "[PicVerification] Regression: unit-box sin(pi x) cos(2 pi y), periodic-Y (bottom/top Periodic), FV Poisson."
                  << std::endl;
        regMetrics["unitboxPeriodicY.enabled"] = 1.0;
        FaceSegmentGraph graph;
        graph.build(mesh);

        PoissonFVSolverCellCentered solver;
        PoissonFVSolverCellCentered::Options opt;
        if (cfg_.solvers)
            opt = cfg_.solvers->poissonCellCentered;

        CellField<double> phi;
        const auto norms = solver.regressionUnitBoxSinCosPeriodicY(mesh, graph, phi, opt);
        regMetrics["unitboxPeriodicY.L2"] = norms.L2;
        regMetrics["unitboxPeriodicY.Linf"] = norms.Linf;
    }

    // ElectricFieldEvaluator regression: analytic phi, compare periodic vs non-periodic Y boundary handling.
    if (doEFieldPeriodicY)
    {
        std::cout << "[PicVerification] Regression: ElectricFieldEvaluator, phi=sin(pi x) cos(2 pi y), periodic-Y (bottom/top Periodic)." << std::endl;
        regMetrics["efieldPeriodicY.enabled"] = 1.0;

        FaceSegmentGraph graph;
        graph.build(mesh);

        // Fill analytic phi on cell centers.
        CellField<double> phi;
        phi.initialize(mesh, 0.0);
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (!isLeafCell(c))
                continue;
            const Vector2D cc = c.centerCoord();
            phi.value[ci] = phiSinCosPeriodicY(cc.x, cc.y);
        }

        PicDomainBC bcPeriodic;
        bcPeriodic.left.type = PicBCType::Dirichlet;
        bcPeriodic.right.type = PicBCType::Dirichlet;
        bcPeriodic.left.value = 0.0;
        bcPeriodic.right.value = 0.0;
        bcPeriodic.bottom.type = PicBCType::Periodic;
        bcPeriodic.top.type = PicBCType::Periodic;

        PicDomainBC bcNonPeriodic = bcPeriodic;
        bcNonPeriodic.bottom.type = PicBCType::Dirichlet;
        bcNonPeriodic.top.type = PicBCType::Dirichlet;
        bcNonPeriodic.bottom.value = 0.0;
        bcNonPeriodic.top.value = 0.0;

        ExplicitPIC::Field::ElectricFieldEvaluator::Options eopt;

        const auto EnPeriodic = ExplicitPIC::Field::ElectricFieldEvaluator::computeFaceNormalField(mesh, graph, bcPeriodic, phi, eopt);
        const auto EnNonPeriodic = ExplicitPIC::Field::ElectricFieldEvaluator::computeFaceNormalField(mesh, graph, bcNonPeriodic, phi, eopt);

        const auto& segs = graph.segments();
        double sumW_all = 0.0, sumErr2_all = 0.0, linf_all = 0.0;
        double sumW_pY = 0.0, sumErr2_pY = 0.0, linf_pY = 0.0;
        double sumW_pY_np = 0.0, sumErr2_pY_np = 0.0, linf_pY_np = 0.0;
        double maxPhiJump = 0.0;

        const double probeEps = 1e-12;
        for (size_t si = 0; si < segs.size(); ++si)
        {
            const FaceSegment& seg = segs[si];
            if (seg.K < 0)
                continue;

            // Only use leaf segments with positive area.
            const double w = std::max(seg.Af, 0.0);
            if (w <= 0.0)
                continue;

            const double EnExact = exactEnAtFaceCenter(seg);
            const double errP = EnPeriodic[si] - EnExact;
            accumWeightedError(errP, w, sumW_all, sumErr2_all, linf_all);

            const bool isSouthNorthBoundary = (seg.boundary != Boundary::non) && (seg.side == FaceSide::South || seg.side == FaceSide::North);
            if (isSouthNorthBoundary)
            {
                accumWeightedError(errP, w, sumW_pY, sumErr2_pY, linf_pY);
                const double errNP = EnNonPeriodic[si] - EnExact;
                accumWeightedError(errNP, w, sumW_pY_np, sumErr2_pY_np, linf_pY_np);

                const Vector2D shift = periodicShiftForSide(mesh, seg.side);
                const int lLbl = findPeriodicMateCellLabel(mesh, seg, shift, probeEps);
                if (lLbl >= 0)
                {
                    const int kIndex = mesh.cells(seg.K).index();
                    const int lIndex = mesh.cells(lLbl).index();
                    if (kIndex >= 0 && kIndex < phi.size() && lIndex >= 0 && lIndex < phi.size())
                    {
                        maxPhiJump = std::max(maxPhiJump, std::abs(phi.value[lIndex] - phi.value[kIndex]));
                    }
                }
            }
        }

        const double L2_all = std::sqrt(sumErr2_all / std::max(sumW_all, 1e-30));
        const double L2_pY = std::sqrt(sumErr2_pY / std::max(sumW_pY, 1e-30));
        const double L2_pY_np = std::sqrt(sumErr2_pY_np / std::max(sumW_pY_np, 1e-30));
        const double ratio = (L2_pY > 0.0) ? (L2_pY_np / L2_pY) : 0.0;

        std::cout << "  E_n error (all segments):     L2=" << L2_all << "  Linf=" << linf_all << std::endl;
        std::cout << "  E_n error (Y boundary, Periodic):     L2=" << L2_pY << "  Linf=" << linf_pY << std::endl;
        std::cout << "  E_n error (Y boundary, non-Periodic): L2=" << L2_pY_np << "  Linf=" << linf_pY_np
                  << "  (ratio=" << ratio << ")" << std::endl;
        std::cout << "  max |phi_jump| across periodic-Y wrap: " << maxPhiJump << std::endl;

        regMetrics["efieldPeriodicY.L2_all"] = L2_all;
        regMetrics["efieldPeriodicY.Linf_all"] = linf_all;
        regMetrics["efieldPeriodicY.L2_boundary_periodic"] = L2_pY;
        regMetrics["efieldPeriodicY.Linf_boundary_periodic"] = linf_pY;
        regMetrics["efieldPeriodicY.L2_boundary_nonperiodic"] = L2_pY_np;
        regMetrics["efieldPeriodicY.Linf_boundary_nonperiodic"] = linf_pY_np;
        regMetrics["efieldPeriodicY.boundary_ratio_np_over_p"] = ratio;
        regMetrics["efieldPeriodicY.maxPhiJump"] = maxPhiJump;
    }

    // Cut-circle metrics regression
    if (doCutCircleMetrics)
    {
        std::cout << "[PicVerification] Regression: cut-circle cut-cell metrics." << std::endl;
        const auto& ccfg = verification.regressionCutCircle;

        regMetrics["cutCircle.enabled"] = 1.0;
        regMetrics["cutCircle.centerX"] = ccfg.center.x;
        regMetrics["cutCircle.centerY"] = ccfg.center.y;
        regMetrics["cutCircle.radius"] = ccfg.radius;
        regMetrics["cutCircle.insideIsNegative"] = ccfg.insideIsNegative ? 1.0 : 0.0;
        regMetrics["cutCircle.alphaSamplingN"] = static_cast<double>(ccfg.alphaSamplingN);

        const int levelSetI = ccfg.levelSetIndex;

        if (levelSetI < 0)
        {
            std::cerr << "[RegressionCutCircle] Invalid levelSetIndex=" << levelSetI << std::endl;
            return 2;
        }

        DictionaryGeometry dictGeometry;
        dictGeometry.points.clear();
        dictGeometry.lines.clear();
        dictGeometry.levelSetName.assign(levelSetI + 1, "verification_circle");
        dictGeometry.innerPoints.assign(levelSetI + 1, ccfg.center);
        dictGeometry.moving.assign(levelSetI + 1, false);
        dictGeometry.linesOnInterface.assign(levelSetI + 1, std::vector<SignInt>{});
        dictGeometry.phaseName.clear();
        dictGeometry.levelSetOperation.clear();
        dictGeometry.levelSetIndexA.clear();
        dictGeometry.levelSetIndexB.clear();
        dictGeometry.boundaryLines.clear();

        ImmersedGeometry<QuadtreeMesh> iGeometry(dictGeometry, mesh);

        fillCircleLevelSet(mesh, iGeometry, levelSetI, ccfg.center, ccfg.radius);
        iGeometry.checkCutInfo(mesh, true);
        Solvers::Poisson::EdgeCutLengthHelper::computeLinear(mesh,
                                                                  iGeometry.levelSets.at(levelSetI),
                                                                  levelSetI,
                                                                  ccfg.cutLengthClampEps);

        FaceSegmentGraph graph;
        graph.build(mesh);

        Solvers::Poisson::CutCellMetrics cutMetrics;
        Solvers::Poisson::CutCellMetrics::BuildOptions mopt;
        mopt.levelSetI = levelSetI;
        mopt.insideIsNegative = ccfg.insideIsNegative;
        mopt.alphaSampleN = ccfg.alphaSamplingN;
        mopt.isActive = [](const QuadtreeMesh&, const Cell& c) { return c.son[0] == -1; };
        cutMetrics.build(mesh, iGeometry, graph, mopt);

        const double areaEst = computeInsideAreaFromAlpha(mesh, cutMetrics);
        const double perimEst = computeInterfacePerimeterFromAgamma(mesh, cutMetrics);
        const double pi = 3.141592653589793238462643383279502884;
        const double areaExact = pi * ccfg.radius * ccfg.radius;
        const double perimExact = 2.0 * pi * ccfg.radius;
        const double areaRel = (areaExact > 0.0) ? std::abs(areaEst - areaExact) / areaExact : 0.0;
        const double perimRel = (perimExact > 0.0) ? std::abs(perimEst - perimExact) / perimExact : 0.0;

        std::cout << std::setprecision(16);
        std::cout << "  levelSetIndex     = " << levelSetI << std::endl;
        std::cout << "  center            = (" << ccfg.center.x << ", " << ccfg.center.y << ")" << std::endl;
        std::cout << "  radius            = " << ccfg.radius << std::endl;
        std::cout << "  insideIsNegative  = " << (ccfg.insideIsNegative ? "true" : "false") << std::endl;
        std::cout << "  alphaSamplingN    = " << ccfg.alphaSamplingN << std::endl;
        std::cout << "  [Area]  est = " << areaEst << " , exact = " << areaExact << " , relErr = " << areaRel << std::endl;
        std::cout << "  [Perim] est = " << perimEst << " , exact = " << perimExact << " , relErr = " << perimRel << std::endl;

        regMetrics["cutCircle.areaEst"] = areaEst;
        regMetrics["cutCircle.areaExact"] = areaExact;
        regMetrics["cutCircle.areaRel"] = areaRel;
        regMetrics["cutCircle.perimEst"] = perimEst;
        regMetrics["cutCircle.perimExact"] = perimExact;
        regMetrics["cutCircle.perimRel"] = perimRel;

        if (ccfg.dumpCSV)
        {
            cutMetrics.dumpCellsCSV(mesh, ccfg.cellsCSV, true);
            cutMetrics.dumpSegmentsCSV(mesh, graph, ccfg.segmentsCSV, false);
            std::cout << "  CSV written: " << ccfg.cellsCSV << " and " << ccfg.segmentsCSV << std::endl;
        }
    }

    // Cut-circle Poisson regression
    if (doCutCirclePoisson)
    {
        std::cout << "[PicVerification] Regression: cut-circle FV Poisson (manufactured solution)." << std::endl;
        const auto& pcfg = verification.regressionCutCirclePoisson;
        const bool bigMesh = (mesh.cellIndexNum() > 200000);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePoisson] big-mesh mode: nCell=" << mesh.cellIndexNum() << "\n";

        regMetrics["cutCirclePoisson.enabled"] = 1.0;
        regMetrics["cutCirclePoisson.radius"] = pcfg.radius;
        regMetrics["cutCirclePoisson.centerX"] = pcfg.center.x;
        regMetrics["cutCirclePoisson.centerY"] = pcfg.center.y;
        regMetrics["cutCirclePoisson.insideIsNegative"] = pcfg.insideIsNegative ? 1.0 : 0.0;
        regMetrics["cutCirclePoisson.alphaSamplingN"] = static_cast<double>(pcfg.alphaSamplingN);
        const int levelSetI = pcfg.levelSetIndex;

        if (levelSetI < 0)
        {
            std::cerr << "[RegressionCutCirclePoisson] Invalid levelSetIndex=" << levelSetI << std::endl;
            return 2;
        }

        DictionaryGeometry dictGeometry;
        dictGeometry.points.clear();
        dictGeometry.lines.clear();
        dictGeometry.levelSetName.assign(levelSetI + 1, "verification_circle");
        dictGeometry.innerPoints.assign(levelSetI + 1, pcfg.center);
        dictGeometry.moving.assign(levelSetI + 1, false);
        dictGeometry.linesOnInterface.assign(levelSetI + 1, std::vector<SignInt>{});
        dictGeometry.phaseName.clear();
        dictGeometry.levelSetOperation.clear();
        dictGeometry.levelSetIndexA.clear();
        dictGeometry.levelSetIndexB.clear();
        dictGeometry.boundaryLines.clear();

        ImmersedGeometry<QuadtreeMesh> iGeometry(dictGeometry, mesh);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePoisson] geometry init done.\n";

        fillCircleLevelSet(mesh, iGeometry, levelSetI, pcfg.center, pcfg.radius);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePoisson] fillCircleLevelSet done.\n";
        iGeometry.checkCutInfo(mesh, true);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePoisson] checkCutInfo done.\n";
        Solvers::Poisson::EdgeCutLengthHelper::computeLinear(mesh,
                                                                  iGeometry.levelSets.at(levelSetI),
                                                                  levelSetI,
                                                                  pcfg.cutLengthClampEps);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePoisson] edge cut-length done.\n";

        const auto isActive = [&](const QuadtreeMesh& m, const Cell& c) -> bool {
            if (!isLeafCell(c))
                return false;

            bool anyInside = false;
            for (int k = 0; k < 4; ++k)
            {
                const int vLbl = c.cellToVert[k];
                const int vi = m.verts(vLbl).index();
                const double psiV = iGeometry.levelSets[levelSetI].value[vi];
                if (isInsidePsi(psiV, pcfg.insideIsNegative))
                {
                    anyInside = true;
                    break;
                }
            }
            if (anyInside)
                return true;
            if (c.isCut[levelSetI])
                return true;

            const double psiC = ExplicitPIC::Geometry::LevelSetQuery::psiAt(m, iGeometry, levelSetI, c.centerCoord());
            return isInsidePsi(psiC, pcfg.insideIsNegative);
        };

        FaceSegmentGraph graph;
        graph.build(mesh, isActive);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePoisson] face graph done.\n";

        Solvers::Poisson::CutCellMetrics cutMetrics;
        Solvers::Poisson::CutCellMetrics::BuildOptions mopt;
        mopt.levelSetI = levelSetI;
        mopt.insideIsNegative = pcfg.insideIsNegative;
        mopt.alphaSampleN = pcfg.alphaSamplingN;
        mopt.isActive = isActive;
        cutMetrics.build(mesh, iGeometry, graph, mopt);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePoisson] cut-metrics done.\n";

        if (pcfg.dumpMetricsCSV)
        {
            cutMetrics.dumpCellsCSV(mesh, pcfg.cellsCSV, true);
            cutMetrics.dumpSegmentsCSV(mesh, graph, pcfg.segmentsCSV, false);
            std::cout << "[RegressionCutCirclePoisson] Metrics CSV written: " << pcfg.cellsCSV
                      << " and " << pcfg.segmentsCSV << std::endl;
        }

        const Vector2D minC = mesh.minCoord;
        const Vector2D maxC = mesh.maxCoord;
        const double Lx = maxC.x - minC.x;
        const double Ly = maxC.y - minC.y;
        const double pi = 3.141592653589793238462643383279502884;
        const double kx = pi / Lx;
        const double ky = pi / Ly;
        const double R2 = pcfg.radius * pcfg.radius;
        const double kk2 = kx * kx + ky * ky;

        CellField<double> rhs;
        rhs.initialize(mesh, 0.0);
        CellField<double> phi;
        phi.initialize(mesh, 0.0);
        CellField<double> phiExact;
        phiExact.initialize(mesh, 0.0);

        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (!isActive(mesh, c))
                continue;

            const Vector2D x = c.centerCoord();
            const double dx = x.x - pcfg.center.x;
            const double dy = x.y - pcfg.center.y;
            const double q = dx * dx + dy * dy - R2;

            const double ax = kx * (x.x - minC.x);
            const double ay = ky * (x.y - minC.y);
            const double sx = std::sin(ax);
            const double cx = std::cos(ax);
            const double sy = std::sin(ay);
            const double cy = std::cos(ay);

            const double q2 = q * q;
            const double st = sx * sy;

            const double phiE = pcfg.coeff * q2 * st;

            const double lap = (16.0 * q + 8.0 * R2 - kk2 * q2) * st
                               + 8.0 * dx * q * kx * cx * sy
                               + 8.0 * dy * q * ky * cy * sx;
            const double f = -pcfg.coeff * lap;

            phiExact.value[ci] = phiE;
            rhs.value[ci] = f;
        }

        // Domain BC is configured in ./system/pic.dict under:
        //   BC { domain { left/right/bottom/top <Type> <Value> } }
        // Default is homogeneous Dirichlet (0) if not specified.
        PicDomainBC domainBC = cfg.domainBC;

        // The manufactured solution used in this regression satisfies
        // homogeneous Dirichlet on the outer box. Warn if user config differs.
        const auto warnIfNonHomogDirichlet = [&](const char* side, const PicBoundarySide& bc) {
            if (bc.type != PicBCType::Dirichlet || std::abs(bc.value) > 0.0)
            {
                std::cout << "[RegressionCutCirclePoisson] WARNING: domain BC on " << side
                          << " is not homogeneous Dirichlet; manufactured-solution error norms may be meaningless." << std::endl;
            }
        };
        warnIfNonHomogDirichlet("left", domainBC.left);
        warnIfNonHomogDirichlet("right", domainBC.right);
        warnIfNonHomogDirichlet("bottom", domainBC.bottom);
        warnIfNonHomogDirichlet("top", domainBC.top);

        PoissonFVSolverCellCentered solver;
        PoissonFVSolverCellCentered::Options sopt;
        if (cfg_.solvers)
            sopt = cfg_.solvers->poissonCellCentered;

        const PicEmbeddedBC* ebBC = nullptr;
        if (pcfg.useEmbeddedBC && cfg.embeddedBC.enabled)
            ebBC = &cfg.embeddedBC;

        PoissonFVSolverCellCentered::SolveReport report;
        solver.solve(mesh, graph, domainBC, rhs, phi, &report, sopt, &cutMetrics, ebBC);

        double sumW = 0.0;
        double sumE2 = 0.0;
        double sumW_cf = 0.0;
        double sumE2_cf = 0.0;
        double sumW_eb = 0.0;
        double sumE2_eb = 0.0;
        double linf = 0.0;
        double hMin = std::numeric_limits<double>::infinity();
        const auto& cm = cutMetrics.cellMetricsByIndex();
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (!isActive(mesh, c))
                continue;

            const double alpha = cm[ci].alpha;
            const double w = alpha * c.area();
            if (w <= 0.0)
                continue;

            const double e = phi.value[ci] - phiExact.value[ci];
            sumW += w;
            sumE2 += e * e * w;
            linf = std::max(linf, std::abs(e));
            const Vector2D h = c.size();
            hMin = std::min(hMin, std::min(h.x, h.y));
        }

        const double l2 = (sumW > 0.0) ? std::sqrt(sumE2 / sumW) : 0.0;

        std::cout << "[RegressionCutCirclePoisson] hMin=" << hMin
                  << "  L2=" << l2
                  << "  Linf=" << linf << std::endl;

        regMetrics["cutCirclePoisson.hMin"] = hMin;
        regMetrics["cutCirclePoisson.L2"] = l2;
        regMetrics["cutCirclePoisson.Linf"] = linf;
        regMetrics["cutCirclePoisson.iterations"] = static_cast<double>(report.iterations);
        regMetrics["cutCirclePoisson.finalResidual"] = report.finalResidual;

        if (pcfg.dumpPhiCSV)
        {
            dumpPhiCSV(pcfg.phiCSV, mesh, cutMetrics, phi, phiExact);
            std::cout << "[RegressionCutCirclePoisson] Wrote: " << pcfg.phiCSV << std::endl;
        }

        if (pcfg.computeEError)
        {
            CellField<Vector2D> Ecell;
            Ecell.initialize(mesh, Vector2D(0.0, 0.0));

            ExplicitPIC::Field::ElectricFieldEvaluator::Options eOpt;
            ExplicitPIC::Field::ElectricFieldEvaluator::evaluateCellCenteredField(
                mesh, graph, &cutMetrics,
                domainBC, ebBC, nullptr,
                phi, Ecell,
                /*EnBySegmentOut=*/nullptr,
                /*EnGammaByCellIndexOut=*/nullptr,
                eOpt);

            double sumW_E = 0.0;
            double sumEx2 = 0.0;
            double sumEy2 = 0.0;
            double sumEmag2 = 0.0;
            double linfEx = 0.0;
            double linfEy = 0.0;
            double linfEmag = 0.0;

            for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
            {
                const int cLbl = mesh.cellIndexToLabel[ci];
                const Cell& c = mesh.cells(cLbl);
                if (!isActive(mesh, c))
                    continue;

                const double alpha = cm[ci].alpha;
                const double w = alpha * c.area();
                if (w <= 0.0)
                    continue;

                const Vector2D x = c.centerCoord();
                const double dx = x.x - pcfg.center.x;
                const double dy = x.y - pcfg.center.y;
                const double q = dx * dx + dy * dy - R2;
                const double q2 = q * q;

                const double ax = kx * (x.x - minC.x);
                const double ay = ky * (x.y - minC.y);
                const double sx = std::sin(ax);
                const double cx = std::cos(ax);
                const double sy = std::sin(ay);
                const double cy = std::cos(ay);

                const double st = sx * sy;
                const double stx = kx * cx * sy;
                const double sty = ky * cy * sx;

                const double dphi_dx = pcfg.coeff * (2.0 * q * (2.0 * dx) * st + q2 * stx);
                const double dphi_dy = pcfg.coeff * (2.0 * q * (2.0 * dy) * st + q2 * sty);
                const Vector2D Eexact(-dphi_dx, -dphi_dy);

                const Vector2D E = Ecell.value[ci];
                const double ex = E.x - Eexact.x;
                const double ey = E.y - Eexact.y;
                const double emag = std::sqrt(ex * ex + ey * ey);

                sumW_E += w;
                sumEx2 += ex * ex * w;
                sumEy2 += ey * ey * w;
                sumEmag2 += emag * emag * w;
                linfEx = std::max(linfEx, std::abs(ex));
                linfEy = std::max(linfEy, std::abs(ey));
                linfEmag = std::max(linfEmag, emag);
            }

            const double ExL2 = (sumW_E > 0.0) ? std::sqrt(sumEx2 / sumW_E) : 0.0;
            const double EyL2 = (sumW_E > 0.0) ? std::sqrt(sumEy2 / sumW_E) : 0.0;
            const double EL2 = (sumW_E > 0.0) ? std::sqrt(sumEmag2 / sumW_E) : 0.0;

            std::cout << "[RegressionCutCirclePoisson] E-error: L2=(" << ExL2 << "," << EyL2
                      << "), |E| L2=" << EL2 << " Linf=(" << linfEx << "," << linfEy
                      << "), |E| Linf=" << linfEmag << std::endl;

            regMetrics["cutCirclePoisson.Ex_L2"] = ExL2;
            regMetrics["cutCirclePoisson.Ey_L2"] = EyL2;
            regMetrics["cutCirclePoisson.E_L2"] = EL2;
            regMetrics["cutCirclePoisson.Ex_Linf"] = linfEx;
            regMetrics["cutCirclePoisson.Ey_Linf"] = linfEy;
            regMetrics["cutCirclePoisson.E_Linf"] = linfEmag;
        }
    }

    // Cut-circle polarization Poisson regression (variable-coefficient Poisson).
    if (doCutCirclePolPoisson)
    {
        std::cout << "[PicVerification] Regression: cut-circle FV polarization Poisson (manufactured solution)." << std::endl;
        const auto& pcfg = verification.regressionCutCirclePolarizationPoisson;
        const bool bigMesh = (mesh.cellIndexNum() > 200000);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePolPoisson] big-mesh mode: nCell=" << mesh.cellIndexNum() << "\n";

        std::string kappaMode = lowerCopy(pcfg.kappaMode);
        if (kappaMode.empty())
            kappaMode = "constant";
        if (kappaMode != "constant" && kappaMode != "sin2d")
        {
            std::cout << "[RegressionCutCirclePolPoisson] WARNING: unknown kappaMode='" << pcfg.kappaMode
                      << "', fallback to constant." << std::endl;
            kappaMode = "constant";
        }
        const bool kappaSin2d = (kappaMode == "sin2d");

        regMetrics["cutCirclePolPoisson.enabled"] = 1.0;
        regMetrics["cutCirclePolPoisson.radius"] = pcfg.radius;
        regMetrics["cutCirclePolPoisson.centerX"] = pcfg.center.x;
        regMetrics["cutCirclePolPoisson.centerY"] = pcfg.center.y;
        regMetrics["cutCirclePolPoisson.insideIsNegative"] = pcfg.insideIsNegative ? 1.0 : 0.0;
        regMetrics["cutCirclePolPoisson.alphaSamplingN"] = static_cast<double>(pcfg.alphaSamplingN);
        regMetrics["cutCirclePolPoisson.kappa"] = pcfg.kappa;
        regMetrics["cutCirclePolPoisson.kappaMode"] = kappaSin2d ? 1.0 : 0.0;
        const int levelSetI = pcfg.levelSetIndex;

        if (levelSetI < 0)
        {
            std::cerr << "[RegressionCutCirclePolPoisson] Invalid levelSetIndex=" << levelSetI << std::endl;
            return 2;
        }

        DictionaryGeometry dictGeometry;
        dictGeometry.points.clear();
        dictGeometry.lines.clear();
        dictGeometry.levelSetName.assign(levelSetI + 1, "verification_circle");
        dictGeometry.innerPoints.assign(levelSetI + 1, pcfg.center);
        dictGeometry.moving.assign(levelSetI + 1, false);
        dictGeometry.linesOnInterface.assign(levelSetI + 1, std::vector<SignInt>{});
        dictGeometry.phaseName.clear();
        dictGeometry.levelSetOperation.clear();
        dictGeometry.levelSetIndexA.clear();
        dictGeometry.levelSetIndexB.clear();
        dictGeometry.boundaryLines.clear();

        ImmersedGeometry<QuadtreeMesh> iGeometry(dictGeometry, mesh);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePolPoisson] geometry init done.\n";

        fillCircleLevelSet(mesh, iGeometry, levelSetI, pcfg.center, pcfg.radius);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePolPoisson] fillCircleLevelSet done.\n";
        iGeometry.checkCutInfo(mesh, true);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePolPoisson] checkCutInfo done.\n";
        Solvers::Poisson::EdgeCutLengthHelper::computeLinear(mesh,
                                                                  iGeometry.levelSets.at(levelSetI),
                                                                  levelSetI,
                                                                  pcfg.cutLengthClampEps);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePolPoisson] edge cut-length done.\n";

        const auto isActive = [&](const QuadtreeMesh& m, const Cell& c) -> bool {
            if (!isLeafCell(c))
                return false;

            bool anyInside = false;
            for (int k = 0; k < 4; ++k)
            {
                const int vLbl = c.cellToVert[k];
                const int vi = m.verts(vLbl).index();
                const double psiV = iGeometry.levelSets[levelSetI].value[vi];
                if (isInsidePsi(psiV, pcfg.insideIsNegative))
                {
                    anyInside = true;
                    break;
                }
            }
            if (anyInside)
                return true;
            if (c.isCut[levelSetI])
                return true;

            const double psiC = ExplicitPIC::Geometry::LevelSetQuery::psiAt(m, iGeometry, levelSetI, c.centerCoord());
            return isInsidePsi(psiC, pcfg.insideIsNegative);
        };

        FaceSegmentGraph graph;
        graph.build(mesh, isActive);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePolPoisson] face graph done.\n";

        Solvers::Poisson::CutCellMetrics cutMetrics;
        Solvers::Poisson::CutCellMetrics::BuildOptions mopt;
        mopt.levelSetI = levelSetI;
        mopt.insideIsNegative = pcfg.insideIsNegative;
        mopt.alphaSampleN = pcfg.alphaSamplingN;
        mopt.isActive = isActive;
        cutMetrics.build(mesh, iGeometry, graph, mopt);
        if (bigMesh)
            std::cerr << "[RegressionCutCirclePolPoisson] cut-metrics done.\n";

        if (pcfg.dumpMetricsCSV)
        {
            cutMetrics.dumpCellsCSV(mesh, pcfg.cellsCSV, true);
            cutMetrics.dumpSegmentsCSV(mesh, graph, pcfg.segmentsCSV, false);
            std::cout << "[RegressionCutCirclePolPoisson] Metrics CSV written: " << pcfg.cellsCSV
                      << " and " << pcfg.segmentsCSV << std::endl;
        }

        const Vector2D minC = mesh.minCoord;
        const Vector2D maxC = mesh.maxCoord;
        const double Lx = maxC.x - minC.x;
        const double Ly = maxC.y - minC.y;
        const double pi = 3.141592653589793238462643383279502884;
        const double kx = pi / Lx;
        const double ky = pi / Ly;
        const double R2 = pcfg.radius * pcfg.radius;
        const double kk2 = kx * kx + ky * ky;

        CellField<double> rhs;
        rhs.initialize(mesh, 0.0);
        CellField<double> phi;
        phi.initialize(mesh, 0.0);
        CellField<double> phiExact;
        phiExact.initialize(mesh, 0.0);

        const double kappaConst =
            (!std::isfinite(pcfg.kappa) || pcfg.kappa < 0.0) ? 0.0 : pcfg.kappa;
        const double kappa0 =
            (!std::isfinite(pcfg.kappa0) || pcfg.kappa0 < 0.0) ? 0.0 : pcfg.kappa0;
        const double kappaAmp = std::isfinite(pcfg.kappaAmp) ? pcfg.kappaAmp : 0.0;
        const double kappaKx = std::isfinite(pcfg.kappaKx) ? pcfg.kappaKx : 0.0;
        const double kappaKy = std::isfinite(pcfg.kappaKy) ? pcfg.kappaKy : 0.0;

        auto evalKappa = [&](const Vector2D& x) -> double {
            if (kappaSin2d)
            {
                const double s =
                    std::sin(kappaKx * x.x) * std::sin(kappaKy * x.y);
                const double raw = kappa0 * (1.0 + kappaAmp * s);
                if (!std::isfinite(raw))
                    return 0.0;
                return std::max(0.0, raw);
            }
            return kappaConst;
        };

        auto evalGradKappa = [&](const Vector2D& x) -> Vector2D {
            if (!kappaSin2d)
                return Vector2D(0.0, 0.0);
            const double dkx = kappa0 * kappaAmp * kappaKx *
                               std::cos(kappaKx * x.x) * std::sin(kappaKy * x.y);
            const double dky = kappa0 * kappaAmp * kappaKy *
                               std::sin(kappaKx * x.x) * std::cos(kappaKy * x.y);
            return Vector2D(std::isfinite(dkx) ? dkx : 0.0,
                            std::isfinite(dky) ? dky : 0.0);
        };

        std::vector<double> kappaByCellIndex(static_cast<size_t>(mesh.cellIndexNum()), 1.0);
        double kappaMin = std::numeric_limits<double>::infinity();
        double kappaMax = -std::numeric_limits<double>::infinity();

        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            const Vector2D x = c.centerCoord();
            const double kappaCell = evalKappa(x);
            kappaByCellIndex[static_cast<size_t>(ci)] = kappaCell;
            kappaMin = std::min(kappaMin, kappaCell);
            kappaMax = std::max(kappaMax, kappaCell);
            if (!isActive(mesh, c))
                continue;

            const double dx = x.x - pcfg.center.x;
            const double dy = x.y - pcfg.center.y;
            const double q = dx * dx + dy * dy - R2;

            const double ax = kx * (x.x - minC.x);
            const double ay = ky * (x.y - minC.y);
            const double sx = std::sin(ax);
            const double cx = std::cos(ax);
            const double sy = std::sin(ay);
            const double cy = std::cos(ay);

            const double q2 = q * q;
            const double st = sx * sy;
            const double stx = kx * cx * sy;
            const double sty = ky * cy * sx;
            const double stxx = -(kx * kx) * st;
            const double styy = -(ky * ky) * st;

            const double phiE = pcfg.coeff * q2 * st;
            const double dphi_dx = pcfg.coeff * (4.0 * dx * q * st + q2 * stx);
            const double dphi_dy = pcfg.coeff * (4.0 * dy * q * st + q2 * sty);
            const double d2phi_dxx = pcfg.coeff * ((4.0 * q + 8.0 * dx * dx) * st
                                                   + 8.0 * dx * q * stx
                                                   + q2 * stxx);
            const double d2phi_dyy = pcfg.coeff * ((4.0 * q + 8.0 * dy * dy) * st
                                                   + 8.0 * dy * q * sty
                                                   + q2 * styy);
            const double lapPhi = d2phi_dxx + d2phi_dyy;

            const Vector2D gradKappa = evalGradKappa(x);
            const double rhsMms = -(gradKappa.x * dphi_dx + gradKappa.y * dphi_dy + kappaCell * lapPhi);
            phiExact.value[ci] = phiE;
            rhs.value[ci] = rhsMms;
        }

        if (!std::isfinite(kappaMin))
            kappaMin = 0.0;
        if (!std::isfinite(kappaMax))
            kappaMax = 0.0;
        regMetrics["cutCirclePolPoisson.kappaMin"] = kappaMin;
        regMetrics["cutCirclePolPoisson.kappaMax"] = kappaMax;
        if (kappaMin < 0.0)
        {
            std::cout << "[RegressionCutCirclePolPoisson] WARNING: kappaMin=" << kappaMin
                      << " < 0 (non-elliptic regime)." << std::endl;
        }

        PicDomainBC domainBC = cfg.domainBC;
        const auto warnIfNonHomogDirichlet = [&](const char* side, const PicBoundarySide& bc) {
            if (bc.type != PicBCType::Dirichlet || std::abs(bc.value) > 0.0)
            {
                std::cout << "[RegressionCutCirclePolPoisson] WARNING: domain BC on " << side
                          << " is not homogeneous Dirichlet; manufactured-solution error norms may be meaningless." << std::endl;
            }
        };
        warnIfNonHomogDirichlet("left", domainBC.left);
        warnIfNonHomogDirichlet("right", domainBC.right);
        warnIfNonHomogDirichlet("bottom", domainBC.bottom);
        warnIfNonHomogDirichlet("top", domainBC.top);

        PolarizationPoissonFVSolverCellCentered solver;
        PolarizationPoissonFVSolverCellCentered::Options sopt;
        if (cfg_.solvers)
        {
            const auto& popt = cfg_.solvers->poissonCellCentered;
            sopt.useSuperLU = popt.useSuperLU;
            sopt.useCG = popt.useCG;
            sopt.useIC = popt.useIC;
            sopt.maxIter = popt.maxIter;
            sopt.relTol = popt.relTol;
            sopt.verbose = popt.verbose;
            sopt.embeddedBCMissingPolicy =
                (popt.embeddedBCMissingPolicy == PoissonFVSolverCellCentered::Options::EmbeddedBCMissingPolicy::Skip)
                    ? PolarizationPoissonFVSolverCellCentered::Options::EmbeddedBCMissingPolicy::Skip
                    : PolarizationPoissonFVSolverCellCentered::Options::EmbeddedBCMissingPolicy::Error;
        }

        const PicEmbeddedBC* ebBC = nullptr;
        if (pcfg.useEmbeddedBC && cfg.embeddedBC.enabled)
            ebBC = &cfg.embeddedBC;

        std::vector<double> kappaBySegment(graph.segments().size(), 1.0);
        const auto& segs = graph.segments();
        for (std::size_t sid = 0; sid < segs.size(); ++sid)
        {
            const FaceSegment& s = segs[sid];
            const int kIdx = (s.K >= 0) ? mesh.cells(s.K).index() : -1;
            const int lIdx = (s.L >= 0) ? mesh.cells(s.L).index() : -1;

            double kK = 1.0;
            double kL = 1.0;
            if (kIdx >= 0 && static_cast<std::size_t>(kIdx) < kappaByCellIndex.size())
                kK = kappaByCellIndex[static_cast<std::size_t>(kIdx)];
            if (lIdx >= 0 && static_cast<std::size_t>(lIdx) < kappaByCellIndex.size())
                kL = kappaByCellIndex[static_cast<std::size_t>(lIdx)];

            double kf = kK;
            if (lIdx >= 0)
            {
                const double denom = std::max(1e-300, kK + kL);
                kf = 2.0 * kK * kL / denom;
            }
            if (!(kf >= 0.0) || !std::isfinite(kf))
                kf = 0.0;
            kappaBySegment[sid] = kf;
        }

        PolarizationPoissonFVSolverCellCentered::SolveReport report;
        const bool ok = solver.solve(mesh, graph, domainBC, rhs, phi,
                                     &kappaBySegment, &kappaByCellIndex,
                                     &report, sopt, &cutMetrics, ebBC);
        if (!ok)
        {
            std::cerr << "[RegressionCutCirclePolPoisson] ERROR: solver failed: " << report.message << std::endl;
            return 3;
        }

        double sumW = 0.0;
        double sumE2 = 0.0;
        double sumW_cf = 0.0;
        double sumE2_cf = 0.0;
        double sumW_eb = 0.0;
        double sumE2_eb = 0.0;
        double linf = 0.0;
        double hMin = std::numeric_limits<double>::infinity();
        const auto& cm = cutMetrics.cellMetricsByIndex();
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (!isActive(mesh, c))
                continue;

            const double alpha = cm[ci].alpha;
            const double w = alpha * c.area();
            if (w <= 0.0)
                continue;

            const double e = phi.value[ci] - phiExact.value[ci];
            sumW += w;
            sumE2 += e * e * w;
            const bool inCfBand = isCompositeBandCell(mesh, c);
            const bool inEbBand = (cm[ci].A_gamma > 0.0);
            if (inCfBand)
            {
                sumW_cf += w;
                sumE2_cf += e * e * w;
            }
            if (inEbBand)
            {
                sumW_eb += w;
                sumE2_eb += e * e * w;
            }
            linf = std::max(linf, std::abs(e));
            const Vector2D h = c.size();
            hMin = std::min(hMin, std::min(h.x, h.y));
        }

        const double l2 = (sumW > 0.0) ? std::sqrt(sumE2 / sumW) : 0.0;
        const double l2_cf = (sumW_cf > 0.0) ? std::sqrt(sumE2_cf / sumW_cf) : 0.0;
        const double l2_eb = (sumW_eb > 0.0) ? std::sqrt(sumE2_eb / sumW_eb) : 0.0;

        std::cout << "[RegressionCutCirclePolPoisson] hMin=" << hMin
                  << "  L2=" << l2
                  << "  Linf=" << linf
                  << "  L2_CFband=" << l2_cf
                  << "  L2_EBband=" << l2_eb << std::endl;

        regMetrics["cutCirclePolPoisson.hMin"] = hMin;
        regMetrics["cutCirclePolPoisson.L2"] = l2;
        regMetrics["cutCirclePolPoisson.Linf"] = linf;
        regMetrics["cutCirclePolPoisson.L2_CFband"] = l2_cf;
        regMetrics["cutCirclePolPoisson.L2_EBband"] = l2_eb;
        regMetrics["cutCirclePolPoisson.iterations"] = static_cast<double>(report.iterations);
        regMetrics["cutCirclePolPoisson.finalResidual"] = report.finalResidual;
        regMetrics["cutCirclePolPoisson.N"] = static_cast<double>(report.N);
        regMetrics["cutCirclePolPoisson.nnz"] = static_cast<double>(report.nnz);
        regMetrics["cutCirclePolPoisson.usedSuperLU"] = report.usedSuperLU ? 1.0 : 0.0;

        if (pcfg.dumpPhiCSV)
        {
            dumpPhiCSV(pcfg.phiCSV, mesh, cutMetrics, phi, phiExact);
            std::cout << "[RegressionCutCirclePolPoisson] Wrote: " << pcfg.phiCSV << std::endl;
        }

        if (pcfg.computeEError)
        {
            CellField<Vector2D> Ecell;
            Ecell.initialize(mesh, Vector2D(0.0, 0.0));

            ExplicitPIC::Field::ElectricFieldEvaluator::Options eOpt;
            ExplicitPIC::Field::ElectricFieldEvaluator::evaluateCellCenteredField(
                mesh, graph, &cutMetrics,
                domainBC, ebBC, nullptr,
                phi, Ecell,
                /*EnBySegmentOut=*/nullptr,
                /*EnGammaByCellIndexOut=*/nullptr,
                eOpt);

            double sumW_E = 0.0;
            double sumEx2 = 0.0;
            double sumEy2 = 0.0;
            double sumEmag2 = 0.0;
            double sumW_E_cf = 0.0;
            double sumEmag2_cf = 0.0;
            double sumW_E_eb = 0.0;
            double sumEmag2_eb = 0.0;
            double linfEx = 0.0;
            double linfEy = 0.0;
            double linfEmag = 0.0;

            for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
            {
                const int cLbl = mesh.cellIndexToLabel[ci];
                const Cell& c = mesh.cells(cLbl);
                if (!isActive(mesh, c))
                    continue;

                const double alpha = cm[ci].alpha;
                const double w = alpha * c.area();
                if (w <= 0.0)
                    continue;

                const Vector2D x = c.centerCoord();
                const double dx = x.x - pcfg.center.x;
                const double dy = x.y - pcfg.center.y;
                const double q = dx * dx + dy * dy - R2;
                const double q2 = q * q;

                const double ax = kx * (x.x - minC.x);
                const double ay = ky * (x.y - minC.y);
                const double sx = std::sin(ax);
                const double cx = std::cos(ax);
                const double sy = std::sin(ay);
                const double cy = std::cos(ay);

                const double st = sx * sy;
                const double stx = kx * cx * sy;
                const double sty = ky * cy * sx;

                const double dphi_dx = pcfg.coeff * (2.0 * q * (2.0 * dx) * st + q2 * stx);
                const double dphi_dy = pcfg.coeff * (2.0 * q * (2.0 * dy) * st + q2 * sty);
                const Vector2D Eexact(-dphi_dx, -dphi_dy);

                const Vector2D E = Ecell.value[ci];
                const double ex = E.x - Eexact.x;
                const double ey = E.y - Eexact.y;
                const double emag = std::sqrt(ex * ex + ey * ey);

                sumW_E += w;
                sumEx2 += ex * ex * w;
                sumEy2 += ey * ey * w;
                sumEmag2 += emag * emag * w;
                const bool inCfBand = isCompositeBandCell(mesh, c);
                const bool inEbBand = (cm[ci].A_gamma > 0.0);
                if (inCfBand)
                {
                    sumW_E_cf += w;
                    sumEmag2_cf += emag * emag * w;
                }
                if (inEbBand)
                {
                    sumW_E_eb += w;
                    sumEmag2_eb += emag * emag * w;
                }
                linfEx = std::max(linfEx, std::abs(ex));
                linfEy = std::max(linfEy, std::abs(ey));
                linfEmag = std::max(linfEmag, emag);
            }

            const double ExL2 = (sumW_E > 0.0) ? std::sqrt(sumEx2 / sumW_E) : 0.0;
            const double EyL2 = (sumW_E > 0.0) ? std::sqrt(sumEy2 / sumW_E) : 0.0;
            const double EL2 = (sumW_E > 0.0) ? std::sqrt(sumEmag2 / sumW_E) : 0.0;
            const double EL2_cf = (sumW_E_cf > 0.0) ? std::sqrt(sumEmag2_cf / sumW_E_cf) : 0.0;
            const double EL2_eb = (sumW_E_eb > 0.0) ? std::sqrt(sumEmag2_eb / sumW_E_eb) : 0.0;

            std::cout << "[RegressionCutCirclePolPoisson] E-error: L2=(" << ExL2 << "," << EyL2
                      << "), |E| L2=" << EL2 << " Linf=(" << linfEx << "," << linfEy
                      << "), |E| Linf=" << linfEmag
                      << " |E|L2_CFband=" << EL2_cf
                      << " |E|L2_EBband=" << EL2_eb << std::endl;

            regMetrics["cutCirclePolPoisson.Ex_L2"] = ExL2;
            regMetrics["cutCirclePolPoisson.Ey_L2"] = EyL2;
            regMetrics["cutCirclePolPoisson.E_L2"] = EL2;
            regMetrics["cutCirclePolPoisson.E_L2_CFband"] = EL2_cf;
            regMetrics["cutCirclePolPoisson.E_L2_EBband"] = EL2_eb;
            regMetrics["cutCirclePolPoisson.Ex_Linf"] = linfEx;
            regMetrics["cutCirclePolPoisson.Ey_Linf"] = linfEy;
            regMetrics["cutCirclePolPoisson.E_Linf"] = linfEmag;
        }
    }

    // ---------------------------------------------------------------------
    // paperA Test 5: embedded conductor benchmark (Dirichlet image response)
    // ---------------------------------------------------------------------
    if (doEmbeddedConductor)
    {
        std::cout << "[PicVerification] Test 5: embedded conductor benchmark (Dirichlet response)." << std::endl;
        const auto& tcfg = verification.regressionEmbeddedConductor;

        regMetrics["test5.enabled"] = 1.0;
        regMetrics["test5.centerX"] = tcfg.center.x;
        regMetrics["test5.centerY"] = tcfg.center.y;
        regMetrics["test5.radius"] = tcfg.radius;
        regMetrics["test5.insideIsNegative"] = tcfg.insideIsNegative ? 1.0 : 0.0;
        regMetrics["test5.alphaSamplingN"] = static_cast<double>(tcfg.alphaSamplingN);
        regMetrics["test5.charge"] = tcfg.charge;
        regMetrics["test5.sigma"] = tcfg.sigma;

        const int levelSetI = tcfg.levelSetIndex;
        if (levelSetI < 0)
        {
            std::cerr << "[Test5] Invalid levelSetIndex=" << levelSetI << std::endl;
            return 2;
        }

        DictionaryGeometry dictGeometry;
        dictGeometry.points.clear();
        dictGeometry.lines.clear();
        dictGeometry.levelSetName.assign(levelSetI + 1, "verification_conductor");
        dictGeometry.innerPoints.assign(levelSetI + 1, tcfg.center);
        dictGeometry.moving.assign(levelSetI + 1, false);
        dictGeometry.linesOnInterface.assign(levelSetI + 1, std::vector<SignInt>{});
        dictGeometry.phaseName.clear();
        dictGeometry.levelSetOperation.clear();
        dictGeometry.levelSetIndexA.clear();
        dictGeometry.levelSetIndexB.clear();
        dictGeometry.boundaryLines.clear();

        ImmersedGeometry<QuadtreeMesh> iGeometry(dictGeometry, mesh);
        fillCircleLevelSet(mesh, iGeometry, levelSetI, tcfg.center, tcfg.radius);
        iGeometry.checkCutInfo(mesh, true);
        Solvers::Poisson::EdgeCutLengthHelper::computeLinear(mesh,
                                                                  iGeometry.levelSets.at(levelSetI),
                                                                  levelSetI,
                                                                  tcfg.cutLengthClampEps);

        const bool plasmaInsideIsNegative = !tcfg.insideIsNegative;
        const auto isActive = [&](const QuadtreeMesh& m, const Cell& c) -> bool {
            if (!isLeafCell(c))
                return false;

            bool anyPlasma = false;
            for (int k = 0; k < 4; ++k)
            {
                const int vLbl = c.cellToVert[k];
                const int vi = m.verts(vLbl).index();
                const double psiV = iGeometry.levelSets[levelSetI].value[vi];
                if (isInsidePsi(psiV, plasmaInsideIsNegative))
                {
                    anyPlasma = true;
                    break;
                }
            }
            if (anyPlasma)
                return true;
            if (c.isCut[levelSetI])
                return true;

            const double psiC = ExplicitPIC::Geometry::LevelSetQuery::psiAt(m, iGeometry, levelSetI, c.centerCoord());
            return isInsidePsi(psiC, plasmaInsideIsNegative);
        };

        FaceSegmentGraph graph;
        graph.build(mesh, isActive);

        Solvers::Poisson::CutCellMetrics cutMetrics;
        Solvers::Poisson::CutCellMetrics::BuildOptions mopt;
        mopt.levelSetI = levelSetI;
        mopt.insideIsNegative = plasmaInsideIsNegative;
        mopt.alphaSampleN = tcfg.alphaSamplingN;
        mopt.isActive = isActive;
        cutMetrics.build(mesh, iGeometry, graph, mopt);

        if (tcfg.dumpMetricsCSV)
        {
            cutMetrics.dumpCellsCSV(mesh, tcfg.cellsCSV, true);
            cutMetrics.dumpSegmentsCSV(mesh, graph, tcfg.segmentsCSV, false);
            std::cout << "[Test5] Metrics CSV written: " << tcfg.cellsCSV
                      << " and " << tcfg.segmentsCSV << std::endl;
        }

        CellField<double> rhs;
        rhs.initialize(mesh, 0.0);
        CellField<double> Q;
        Q.initialize(mesh, 0.0);
        CellField<double> rho;
        rho.initialize(mesh, 0.0);
        CellField<double> phi;
        phi.initialize(mesh, 0.0);

        const auto& cm = cutMetrics.cellMetricsByIndex();
        const Vector2D xp = tcfg.particle;
        std::string sampleMode = lowerCopy(tcfg.sampleMode);
        if (sampleMode != "hostcell" && sampleMode != "matched")
        {
            std::cout << "[Test5] WARNING: unknown sampleMode='" << tcfg.sampleMode
                      << "', fallback to matched." << std::endl;
            sampleMode = "matched";
        }
        std::string sourceShape = lowerCopy(tcfg.sourceShape);
        if (sourceShape != "delta" && sourceShape != "gaussian")
        {
            std::cout << "[Test5] WARNING: unknown sourceShape='" << tcfg.sourceShape
                      << "', fallback to gaussian." << std::endl;
            sourceShape = "gaussian";
        }
        std::string chargeUnit = lowerCopy(tcfg.chargeUnit);
        if (chargeUnit != "code" && chargeUnit != "si")
        {
            std::cout << "[Test5] WARNING: unknown chargeUnit='" << tcfg.chargeUnit
                      << "', fallback to code." << std::endl;
            chargeUnit = "code";
        }
        if (chargeUnit == "si")
        {
            std::cout << "[Test5] WARNING: chargeUnit=si is currently interpreted as code charge in Test5." << std::endl;
        }

        const bool useMatchedSample = (sampleMode == "matched");
        bool useGaussianSource = (sourceShape == "gaussian");
        const double chargeCode = tcfg.charge;

        double hFine = std::numeric_limits<double>::infinity();
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            if (!cutMetrics.isActiveCellIndex(ci))
                continue;
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (!isLeafCell(c))
                continue;
            const double alpha = cm[ci].alpha;
            if (!(alpha > 0.0))
                continue;
            hFine = std::min(hFine, minCellSize(c));
        }
        if (!std::isfinite(hFine))
            hFine = 0.0;

        double sigma = tcfg.sigma;
        if (useGaussianSource && tcfg.sigmaAutoFromFine && hFine > 0.0)
            sigma = std::max(1e-12, tcfg.sigmaFineFactor * hFine);
        if (useGaussianSource && !(sigma > 0.0))
        {
            std::cout << "[Test5] WARNING: gaussian source has non-positive sigma; fallback to delta source." << std::endl;
            useGaussianSource = false;
            sigma = 0.0;
        }
        if (!useGaussianSource)
            sigma = 0.0;

        regMetrics["test5.sampleModeMatched"] = useMatchedSample ? 1.0 : 0.0;
        regMetrics["test5.sourceShapeGaussian"] = useGaussianSource ? 1.0 : 0.0;
        regMetrics["test5.sigmaUsed"] = sigma;
        regMetrics["test5.hFine"] = hFine;

        std::vector<double> alphaByCellIndex(static_cast<size_t>(mesh.cellIndexNum()), 1.0);
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            double a = cm[ci].alpha;
            if (!std::isfinite(a))
                a = 1.0;
            if (a < 0.0)
                a = 0.0;
            if (a > 1.0)
                a = 1.0;
            alphaByCellIndex[static_cast<size_t>(ci)] = a;
        }

        CouplingWeightsCache matchedCache;
        bool matchedCacheReady = false;
        auto buildMatchedCache = [&](double hpIn) -> bool {
            std::vector<Vector2D> positions(1, xp);
            std::vector<double> hp(1, hpIn);
            CouplingWeightsCache::BuildOptions wopt;
            wopt.levelSetI = levelSetI;
            wopt.levelSetVert = &iGeometry.levelSets.at(levelSetI);
            wopt.insideIsNegative = tcfg.insideIsNegative;
            wopt.alphaByCellIndex = &alphaByCellIndex;
            wopt.enableCompositeBand = true;
            wopt.enableMomentExactCut = true;
            wopt.parallelBuild = false;
            wopt.deterministicOrdering = true;
            matchedCache.build(mesh, positions, hp, wopt);
            matchedCacheReady = (matchedCache.numParticles() == 1 && matchedCache.nnz() > 0);
            return matchedCacheReady;
        };

        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            Q.value[ci] = 0.0;
            rho.value[ci] = 0.0;
            rhs.value[ci] = 0.0;
        }

        if (useGaussianSource)
        {
            double sumW = 0.0;
            std::vector<double> w(mesh.cellIndexNum(), 0.0);
            const double inv2s2 = 1.0 / (2.0 * sigma * sigma);
            for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
            {
                if (!cutMetrics.isActiveCellIndex(ci))
                    continue;
                const int cLbl = mesh.cellIndexToLabel[ci];
                const Cell& c = mesh.cells(cLbl);
                if (!isLeafCell(c))
                    continue;
                const double alpha = cm[ci].alpha;
                if (!(alpha > 0.0))
                    continue;
                const Vector2D cc = c.centerCoord();
                const double dx = cc.x - xp.x;
                const double dy = cc.y - xp.y;
                const double wi = std::exp(-(dx * dx + dy * dy) * inv2s2);
                w[ci] = wi;
                sumW += wi * alpha * c.area();
            }
            if (sumW > 0.0)
            {
                for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
                {
                    if (!cutMetrics.isActiveCellIndex(ci))
                        continue;
                    const int cLbl = mesh.cellIndexToLabel[ci];
                    const Cell& c = mesh.cells(cLbl);
                    if (!isLeafCell(c))
                        continue;
                    const double alpha = cm[ci].alpha;
                    if (!(alpha > 0.0))
                        continue;
                    const double area = alpha * c.area();
                    const double rho_i = chargeCode * w[ci] / sumW;
                    Q.value[ci] = rho_i * area;
                }
            }
            else
            {
                const int hostLbl = mesh.findCellLabel(xp);
                if (hostLbl >= 0)
                {
                    const int hostCi = mesh.cells(hostLbl).index();
                    if (hostCi >= 0 && hostCi < mesh.cellIndexNum())
                        Q.value[hostCi] = chargeCode;
                }
            }
        }
        else
        {
            if (useMatchedSample && buildMatchedCache(0.0))
            {
                DepositQK depositor;
                std::vector<double> qp(1, chargeCode);
                (void)depositor.deposit(qp, matchedCache, &Q);
            }
            else
            {
                const int cLbl = mesh.findCellLabel(xp);
                if (cLbl >= 0)
                {
                    const Cell& c = mesh.cells(cLbl);
                    const int ci = c.index();
                    if (ci >= 0 && ci < mesh.cellIndexNum())
                        Q.value[ci] = chargeCode;
                }
            }
        }

        double sumCharge = 0.0;
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            if (!cutMetrics.isActiveCellIndex(ci))
                continue;
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (!isLeafCell(c))
                continue;
            const double alpha = cm[ci].alpha;
            if (!(alpha > 0.0))
                continue;
            const double area = std::max(alpha * c.area(), 1e-30);
            rho.value[ci] = Q.value[ci] / area;
            rhs.value[ci] = rho.value[ci];
            sumCharge += Q.value[ci];
        }

        regMetrics["test5.chargeSum"] = sumCharge;

        PicDomainBC domainBC = cfg.domainBC;
        const auto warnIfNotDirichlet0 = [&](const char* side, const PicBoundarySide& bc) {
            if (bc.type != PicBCType::Dirichlet || std::abs(bc.value) > 0.0)
            {
                std::cout << "[Test5] WARNING: domain BC on " << side
                          << " is not homogeneous Dirichlet; image-response interpretation may include box effects." << std::endl;
            }
        };
        warnIfNotDirichlet0("left", domainBC.left);
        warnIfNotDirichlet0("right", domainBC.right);
        warnIfNotDirichlet0("bottom", domainBC.bottom);
        warnIfNotDirichlet0("top", domainBC.top);

        const PicEmbeddedBC* ebBC = nullptr;
        if (cfg.embeddedBC.enabled)
            ebBC = &cfg.embeddedBC;
        else
            std::cout << "[Test5] WARNING: embedded BC is not enabled; conductor response will be missing." << std::endl;

        PoissonFVSolverCellCentered solver;
        PoissonFVSolverCellCentered::Options sopt;
        if (cfg_.solvers)
            sopt = cfg_.solvers->poissonCellCentered;

        PoissonFVSolverCellCentered::SolveReport report;
        solver.solve(mesh, graph, domainBC, rhs, phi, &report, sopt, &cutMetrics, ebBC);

        // Electric field reconstruction.
        CellField<Vector2D> Ecell;
        Ecell.initialize(mesh, Vector2D(0.0, 0.0));
        ExplicitPIC::Field::ElectricFieldEvaluator::Options eOpt;
        std::vector<double> EnGamma;
        ExplicitPIC::Field::ElectricFieldEvaluator::evaluateCellCenteredField(
            mesh, graph, &cutMetrics,
            domainBC, ebBC, nullptr,
            phi, Ecell,
            /*EnBySegmentOut=*/nullptr,
            /*EnGammaByCellIndexOut=*/&EnGamma,
            eOpt);

        // Field at particle (matched gather by default; host-cell center for legacy mode).
        Vector2D Ep(0.0, 0.0);
        if (useMatchedSample)
        {
            if (!matchedCacheReady)
            {
                const double hpForSample = (sigma > 0.0) ? sigma : 0.0;
                (void)buildMatchedCache(hpForSample);
            }
            if (matchedCacheReady)
            {
                std::vector<Vector2D> EpTmp;
                ExplicitPIC::Field::MatchedGather::gatherVector2D(matchedCache, Ecell, EpTmp, false);
                if (!EpTmp.empty())
                    Ep = EpTmp[0];
            }
        }
        if (!useMatchedSample || (!std::isfinite(Ep.x) || !std::isfinite(Ep.y)))
        {
            const int cLbl = mesh.findCellLabel(xp);
            if (cLbl >= 0)
            {
                const int ci = mesh.cells(cLbl).index();
                if (ci >= 0 && ci < mesh.cellIndexNum())
                    Ep = Ecell.value[ci];
            }
        }

        const double Emag = std::sqrt(Ep.x * Ep.x + Ep.y * Ep.y);
        const Vector2D toCenter = Vector2D(tcfg.center.x - xp.x, tcfg.center.y - xp.y);
        const double r = std::sqrt(toCenter.x * toCenter.x + toCenter.y * toCenter.y);
        double cosAlign = 0.0;
        if (Emag > 0.0 && r > 0.0)
            cosAlign = (Ep.x * toCenter.x + Ep.y * toCenter.y) / (Emag * r);

        regMetrics["test5.px"] = xp.x;
        regMetrics["test5.py"] = xp.y;
        regMetrics["test5.Ex"] = Ep.x;
        regMetrics["test5.Ey"] = Ep.y;
        regMetrics["test5.Emag"] = Emag;
        regMetrics["test5.cosAlign"] = cosAlign;
        regMetrics["test5.iterations"] = static_cast<double>(report.iterations);
        regMetrics["test5.finalResidual"] = report.finalResidual;

        double EnMax = -std::numeric_limits<double>::infinity();
        double EnMin = std::numeric_limits<double>::infinity();
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            if (!cutMetrics.isActiveCellIndex(ci))
                continue;
            const auto& m = cm[ci];
            if (!(m.A_gamma > 0.0))
                continue;
            const double En = (ci < static_cast<int>(EnGamma.size())) ? EnGamma[ci] : 0.0;
            EnMax = std::max(EnMax, En);
            EnMin = std::min(EnMin, En);
        }
        if (!std::isfinite(EnMax)) EnMax = 0.0;
        if (!std::isfinite(EnMin)) EnMin = 0.0;
        regMetrics["test5.EnMax"] = EnMax;
        regMetrics["test5.EnMin"] = EnMin;

        if (tcfg.dumpPhiCSV)
        {
            dumpPhiRhoCSV(tcfg.phiCSV, mesh, cutMetrics, phi, rho);
            std::cout << "[Test5] Wrote: " << tcfg.phiCSV << std::endl;
        }
        if (tcfg.dumpEnCSV)
        {
            dumpEnCSV(tcfg.enCSV, mesh, cutMetrics, EnGamma);
            std::cout << "[Test5] Wrote: " << tcfg.enCSV << std::endl;
        }
    }

    // ---------------------------------------------------------------------
    // Test 3: single-particle numerical self-force under periodic neutrality
    // ---------------------------------------------------------------------
    if (doSelfForceTest3)
    {
        std::cout << "[PicVerification] Test 3: single-particle numerical self-force (periodic neutrality)." << std::endl;

        const auto& sf = verification.selfForceTest3;
        const double q = sf.particle.q;

        // Write base mesh snapshot (Tecplot FEPOINT .dat).
        mesh.output("self_force_test3_base");

        auto wantSubtest = [&](const std::string& name) -> bool {
            if (sf.subtests.empty())
                return true;
            const std::string key = lowerCopy(name);
            for (const auto& s : sf.subtests)
            {
                if (lowerCopy(s) == key)
                    return true;
            }
            return false;
        };

        std::vector<std::string> ablations = sf.ablations;
        if (ablations.empty())
            ablations = {"A0", "A1"};

        const auto isPeriodic = [&](const PicBoundarySide& bc) -> bool {
            return bc.type == PicBCType::Periodic;
        };
        const bool periodicX = isPeriodic(cfg.domainBC.left) && isPeriodic(cfg.domainBC.right);
        const bool periodicY = isPeriodic(cfg.domainBC.bottom) && isPeriodic(cfg.domainBC.top);

        if (!(periodicX && periodicY))
        {
            std::cout << "[PicVerification][Test3] WARNING: domain BC is not fully periodic; self-force test expects periodic x/y." << std::endl;
        }

        auto buildEmbeddedContext = [&](QuadtreeMesh& m, FaceSegmentGraph& graph) -> std::unique_ptr<EmbeddedContext> {
            auto ctx = std::make_unique<EmbeddedContext>();
            ctx->dict = std::make_unique<DictionaryGeometry>();
            ctx->dict->readDict(cfg_.units.get());
            const auto embeddedBoundaryIds = ctx->dict->immersedBoundaryIds();
            if (embeddedBoundaryIds.empty())
                return nullptr;
            ctx->geom = std::make_unique<ImmersedGeometry<QuadtreeMesh>>(*ctx->dict, m);

            ctx->geom->initializeLevelSets(m);
            ctx->geom->reinitializeLevelSets(m, 0.0);
            ctx->geom->checkCutInfo(m, true);

            ctx->levelSetI = -1;
            for (const int boundaryId : embeddedBoundaryIds)
            {
                auto it = ctx->dict->boundaryById.find(boundaryId);
                if (it == ctx->dict->boundaryById.end())
                    continue;
                if (it->second.levelSetIndex >= 0)
                {
                    ctx->levelSetI = it->second.levelSetIndex;
                    break;
                }
            }
            if (ctx->levelSetI < 0)
                ctx->levelSetI = 0;
            if (ctx->levelSetI < 0 || ctx->levelSetI >= static_cast<int>(ctx->geom->levelSets.size()))
            {
                std::cerr << "[PicVerification][Test3] WARNING: invalid embedded levelSetIndex=" << ctx->levelSetI << std::endl;
                return nullptr;
            }
            if (!ctx->dict->inferLevelSetInsideIsNegative(ctx->levelSetI, ctx->insideIsNegative))
            {
                std::cerr << "[PicVerification][Test3] WARNING: failed to infer embedded insideIsNegative from geometry.dict. "
                          << "Using runtime default true." << std::endl;
                ctx->insideIsNegative = true;
            }

            Solvers::Poisson::EdgeCutLengthHelper::computeLinear(
                m, ctx->geom->levelSets.at(ctx->levelSetI), ctx->levelSetI, 1e-12);

            const bool plasmaInsideIsNegative = !ctx->insideIsNegative;
            const auto isActive = [&](const QuadtreeMesh& mesh, const Cell& c) -> bool {
                if (!isLeafCell(c))
                    return false;

                bool anyPlasma = false;
                for (int k = 0; k < 4; ++k)
                {
                    const int vLbl = c.cellToVert[k];
                    const int vi = mesh.verts(vLbl).index();
                    const double psiV = ctx->geom->levelSets[ctx->levelSetI].value[vi];
                    if (isInsidePsi(psiV, plasmaInsideIsNegative))
                    {
                        anyPlasma = true;
                        break;
                    }
                }
                if (anyPlasma)
                    return true;
                if (c.isCut[ctx->levelSetI])
                    return true;

                const double psiC = ExplicitPIC::Geometry::LevelSetQuery::psiAt(mesh, *ctx->geom, ctx->levelSetI, c.centerCoord());
                return isInsidePsi(psiC, plasmaInsideIsNegative);
            };

            graph.build(m, isActive);

            Solvers::Poisson::CutCellMetrics::BuildOptions mopt;
            mopt.levelSetI = ctx->levelSetI;
            mopt.insideIsNegative = plasmaInsideIsNegative;
            mopt.alphaSampleN = 6;
            mopt.isActive = isActive;
            ctx->metrics.build(m, *ctx->geom, graph, mopt);

            ctx->alphaByCellIndex.resize(static_cast<size_t>(m.cellIndexNum()), 1.0);
            const auto& cm = ctx->metrics.cellMetricsByIndex();
            for (int ci = 0; ci < m.cellIndexNum(); ++ci)
            {
                double a = cm[ci].alpha;
                if (!ctx->metrics.isActiveCellIndex(ci))
                    a = 0.0;
                ctx->alphaByCellIndex[static_cast<size_t>(ci)] = a;
            }

            return ctx;
        };

        // Helper: build a linear sample path.
        auto buildPath = [](const Vector2D& a, const Vector2D& b, int n) {
            std::vector<Vector2D> pts;
            if (n < 2)
                n = 2;
            pts.reserve(static_cast<size_t>(n));
            for (int i = 0; i < n; ++i)
            {
                const double t = (n == 1) ? 0.0 : static_cast<double>(i) / static_cast<double>(n - 1);
                pts.push_back(Vector2D(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)));
            }
            return pts;
        };

        // Summary accumulators.
        struct SummaryStats
        {
            int n = 0;
            double maxFstar = 0.0;
            double sumFstar2 = 0.0;
            double maxFsi = 0.0;
            double sumFsi2 = 0.0;
            double maxAbsSumQAfter = 0.0;
            double maxAbsPhiMean = 0.0;
            double maxAbsRGhost = 0.0;
            double sumRGhost2 = 0.0;
        };
        std::map<std::string, SummaryStats> uniformStatsByN;
        std::map<std::string, SummaryStats> traverseStatsByAbl;
        std::map<std::string, SummaryStats> traverseCutStatsByAbl;
        std::map<std::string, SummaryStats> cutCrossStatsByAbl;
        std::map<std::string, SummaryStats> crossStatsByAbl;

        auto updateStats = [](SummaryStats& s,
                              double absFstar,
                              double absFsi,
                              double absSumQAfter,
                              double absPhiMean,
                              double absRGhost) {
            s.n += 1;
            s.maxFstar = std::max(s.maxFstar, absFstar);
            s.sumFstar2 += absFstar * absFstar;
            s.maxFsi = std::max(s.maxFsi, absFsi);
            s.sumFsi2 += absFsi * absFsi;
            s.maxAbsSumQAfter = std::max(s.maxAbsSumQAfter, absSumQAfter);
            s.maxAbsPhiMean = std::max(s.maxAbsPhiMean, absPhiMean);
            s.maxAbsRGhost = std::max(s.maxAbsRGhost, absRGhost);
            s.sumRGhost2 += absRGhost * absRGhost;
        };

        auto writeSelfForceHeader = [](std::ostream& out) {
            out << "sampleId,ablTag,testId,meshN,x,y,Ex,Ey,absE,"
                   "Fx_SI,Fy_SI,absF_SI,Fstar,"
                   "sumQ_before,sumQ_after,rhoBar,phiMean,h_local,hp_eff,cellLevel,"
                   "inCompositeBand,distanceToTransitionBand,"
                   "rG_host_abs,rG_L2,rG_Linf\n";
        };

        auto runSamplesToStream = [&](std::ostream& out,
                                      const QuadtreeMesh& m,
                                      const FaceSegmentGraph& graph,
                                      const std::vector<BandCellBox>& bandCells,
                                      const std::string& testId,
                                      const std::vector<Vector2D>& samples,
                                      const std::vector<std::string>& abls,
                                      std::map<std::string, SummaryStats>& statsByAbl,
                                      bool wrapPos,
                                      const EmbeddedContext* embed)
        {
            out << std::setprecision(16);
            const UnitsContext* units = cfg_.units.get();
            const double eps0 = PhysicalConstants::epsilon0_Fpm;

            const int meshN = m.cellIndexNum();
            const auto alphaForCell = [&](int ci) -> double {
                if (!embed || embed->alphaByCellIndex.empty())
                    return 1.0;
                if (ci < 0 || ci >= static_cast<int>(embed->alphaByCellIndex.size()))
                    return 1.0;
                double a = embed->alphaByCellIndex[static_cast<size_t>(ci)];
                if (!std::isfinite(a))
                    a = 1.0;
                if (a < 0.0)
                    a = 0.0;
                if (a > 1.0)
                    a = 1.0;
                return a;
            };

            const double totalVolume = [&]() {
                double v = 0.0;
                for (int ci = 0; ci < m.cellIndexNum(); ++ci)
                {
                    const int cLbl = m.cellIndexToLabel[ci];
                    const Cell& c = m.cells(cLbl);
                    if (!isLeafCell(c))
                        continue;
                    v += c.area() * alphaForCell(ci);
                }
                return v;
            }();

            PoissonFVSolverCellCentered solver;
            PoissonFVSolverCellCentered::Options sopt;
            if (cfg_.solvers)
                sopt = cfg_.solvers->poissonCellCentered;

            CellField<double> Q;
            Q.initialize(m, 0.0);
            CellField<double> rhs;
            rhs.initialize(m, 0.0);
            CellField<double> phi;
            phi.initialize(m, 0.0);
            CellField<Vector2D> Ecell;
            Ecell.initialize(m, Vector2D(0.0, 0.0));

            for (const auto& ablTagRaw : abls)
            {
                const std::string ablTag = ablTagRaw;
                const std::string ablLower = lowerCopy(ablTagRaw);

                const bool useComposite = (ablLower == "a1" || ablLower == "a3");
                const bool useMomentExact = (ablLower == "a2" || ablLower == "a3");

                for (size_t si = 0; si < samples.size(); ++si)
                {
                    Vector2D xp = samples[si];
                    if (wrapPos)
                        xp = wrapPeriodic(xp, m, periodicX, periodicY);

                    // Build coupling weights for a single particle.
                    CouplingWeightsCache cache;
                    std::vector<Vector2D> positions(1, xp);
                    std::vector<double> hp(1, 0.0);

                    CouplingWeightsCache::BuildOptions wopt;
                    wopt.periodicX = periodicX;
                    wopt.periodicY = periodicY;
                    wopt.enableCompositeBand = useComposite;
                    wopt.enableMomentExactCut = useMomentExact;
                    if (embed && embed->geom)
                    {
                        wopt.levelSetI = embed->levelSetI;
                        wopt.levelSetVert = &embed->geom->levelSets.at(embed->levelSetI);
                        wopt.insideIsNegative = embed->insideIsNegative;
                        wopt.alphaByCellIndex = embed->alphaByCellIndex.empty() ? nullptr : &embed->alphaByCellIndex;
                    }
                    wopt.parallelBuild = false;
                    wopt.deterministicOrdering = true;
                    cache.build(m, positions, hp, wopt);

                    // Deposit Q_K (charge amount).
                    for (int ci = 0; ci < m.cellIndexNum(); ++ci)
                        Q.value[ci] = 0.0;
                    DepositQK depositor;
                    std::vector<double> qp(1, q);
                    (void)depositor.deposit(qp, cache, &Q);

                    double sumQ_before = 0.0;
                    for (int ci = 0; ci < m.cellIndexNum(); ++ci)
                    {
                        const int cLbl = m.cellIndexToLabel[ci];
                        const Cell& c = m.cells(cLbl);
                        if (!isLeafCell(c))
                            continue;
                        sumQ_before += Q.value[ci];
                    }

                    double rhoBar = 0.0;
                    if (sf.neutrality.enabled && totalVolume > 0.0)
                    {
                        rhoBar = sumQ_before / totalVolume;
                        for (int ci = 0; ci < m.cellIndexNum(); ++ci)
                        {
                            const int cLbl = m.cellIndexToLabel[ci];
                            const Cell& c = m.cells(cLbl);
                            if (!isLeafCell(c))
                                continue;
                            const double alpha = alphaForCell(ci);
                            if (alpha <= 0.0)
                                continue;
                            Q.value[ci] -= rhoBar * c.area() * alpha;
                        }
                    }

                    double sumQ_after = 0.0;
                    for (int ci = 0; ci < m.cellIndexNum(); ++ci)
                    {
                        const int cLbl = m.cellIndexToLabel[ci];
                        const Cell& c = m.cells(cLbl);
                        if (!isLeafCell(c))
                            continue;
                        sumQ_after += Q.value[ci];
                    }

                    // Build RHS = rho = Q / volume (code units).
                    for (int ci = 0; ci < m.cellIndexNum(); ++ci)
                    {
                        const int cLbl = m.cellIndexToLabel[ci];
                        const Cell& c = m.cells(cLbl);
                        if (!isLeafCell(c))
                        {
                            rhs.value[ci] = 0.0;
                            continue;
                        }
                        const double alpha = alphaForCell(ci);
                        if (alpha <= 0.0)
                        {
                            rhs.value[ci] = 0.0;
                            continue;
                        }
                        const double vol = std::max(1e-300, c.area() * alpha);
                        rhs.value[ci] = Q.value[ci] / vol;
                    }

                    PoissonFVSolverCellCentered::SolveReport report;
                    const Solvers::Poisson::CutCellMetrics* cm = (embed ? &embed->metrics : nullptr);
                    const PicEmbeddedBC* eb = (cfg.embeddedBC.enabled ? &cfg.embeddedBC : nullptr);
                    solver.solve(m, graph, cfg.domainBC, rhs, phi, &report, sopt, cm, eb, nullptr);

                    double phiMean = 0.0;
                    if (sf.gaugeFix.enabled && totalVolume > 0.0)
                    {
                        double sumPhi = 0.0;
                        for (int ci = 0; ci < m.cellIndexNum(); ++ci)
                        {
                            const int cLbl = m.cellIndexToLabel[ci];
                            const Cell& c = m.cells(cLbl);
                            if (!isLeafCell(c))
                                continue;
                            const double alpha = alphaForCell(ci);
                            if (alpha <= 0.0)
                                continue;
                            sumPhi += phi.value[ci] * c.area() * alpha;
                        }
                        phiMean = sumPhi / totalVolume;
                        for (int ci = 0; ci < m.cellIndexNum(); ++ci)
                            phi.value[ci] -= phiMean;
                    }

                    // Evaluate cell-centered E and gather to particle.
                    const Solvers::Poisson::CutCellMetrics* ecm = (embed ? &embed->metrics : nullptr);
                    const PicEmbeddedBC* eeb = (cfg.embeddedBC.enabled ? &cfg.embeddedBC : nullptr);
                    ExplicitPIC::Field::ElectricFieldEvaluator::evaluateCellCenteredField(
                        m, graph, ecm, cfg.domainBC, eeb, nullptr, phi, Ecell);

                    // Same-operator Gauss residual diagnostics.
                    CellField<double> rG;
                    rG.initialize(m, 0.0);
                    const auto diag = ExplicitPIC::Diagnostics::PicDiagnostics::compute(
                        m, graph, cfg.domainBC, rhs, phi, &Q, q, &rG, cm, eb);

                    std::vector<Vector2D> Ep;
                    ExplicitPIC::Field::MatchedGather::gatherVector2D(cache, Ecell, Ep, false);
                    Vector2D E = Ep.empty() ? Vector2D(0.0, 0.0) : Ep[0];

                    const double absE = std::sqrt(E.x * E.x + E.y * E.y);
                    const double Ex_SI = units ? units->ElectricFieldCodeToSI(E.x) : E.x;
                    const double Ey_SI = units ? units->ElectricFieldCodeToSI(E.y) : E.y;
                    const double Fx_SI = q * Ex_SI;
                    const double Fy_SI = q * Ey_SI;
                    const double absF_SI = std::sqrt(Fx_SI * Fx_SI + Fy_SI * Fy_SI);

                    // Local cell info.
                    const int hostLbl = m.findCellLabel(xp);
                    double h_local = 0.0;
                    int level = -1;
                    if (hostLbl >= 0)
                    {
                        const Cell& hc = m.cells(hostLbl);
                        h_local = minCellSize(hc);
                        level = hc.level();
                    }
                    const int hostCi = (hostLbl >= 0) ? m.cells(hostLbl).index() : -1;
                    double absRGhost = 0.0;
                    if (hostCi >= 0 && hostCi < m.cellIndexNum())
                        absRGhost = std::abs(rG.value[hostCi]);
                    double hp_eff = h_local;
                    const auto& hpEffList = cache.effectiveHp();
                    if (!hpEffList.empty() && std::isfinite(hpEffList[0]) && hpEffList[0] > 0.0)
                        hp_eff = hpEffList[0];
                    const double hp_si = units ? units->LengthCodeToSI(hp_eff) : hp_eff;
                    const double denom = (hp_si > 0.0) ? (q * q / (eps0 * hp_si * hp_si)) : 0.0;
                    const double Fstar = (denom > 0.0) ? (absF_SI / denom) : 0.0;

                    bool inBand = false;
                    const auto& inBandList = cache.inCompositeBand();
                    if (!inBandList.empty() && inBandList[0] != 0)
                        inBand = true;

                    double distToBand = -1.0;
                    if (!bandCells.empty())
                    {
                        if (inBand)
                        {
                            distToBand = 0.0;
                        }
                        else
                        {
                            double dmin = std::numeric_limits<double>::infinity();
                            for (const auto& bc : bandCells)
                                dmin = std::min(dmin, distancePointToBox(xp, bc));
                            distToBand = std::isfinite(dmin) ? dmin : -1.0;
                        }
                    }

                    out << si << "," << ablTag << "," << testId << "," << meshN << ","
                        << xp.x << "," << xp.y << ","
                        << E.x << "," << E.y << ","
                        << absE << ","
                        << Fx_SI << "," << Fy_SI << "," << absF_SI << "," << Fstar << ","
                        << sumQ_before << "," << sumQ_after << "," << rhoBar << "," << phiMean << ","
                        << h_local << "," << hp_eff << "," << level << ","
                        << (inBand ? 1 : 0) << "," << distToBand << ","
                        << absRGhost << "," << diag.rG_L2 << "," << diag.rG_Linf
                        << "\n";

                    updateStats(statsByAbl[ablTag],
                                Fstar,
                                absF_SI,
                                std::abs(sumQ_after),
                                std::abs(phiMean),
                                absRGhost);
                }
            }
        };

        auto runSamplesToFile = [&](const std::string& outPath,
                                    const QuadtreeMesh& m,
                                    const FaceSegmentGraph& graph,
                                    const std::vector<BandCellBox>& bandCells,
                                    const std::string& testId,
                                    const std::vector<Vector2D>& samples,
                                    const std::vector<std::string>& abls,
                                    std::map<std::string, SummaryStats>& statsByAbl,
                                    bool wrapPos,
                                    const EmbeddedContext* embed)
        {
            std::ofstream out(outPath);
            if (!out)
            {
                std::cerr << "[PicVerification][Test3] ERROR: failed to open " << outPath << std::endl;
                return;
            }
            writeSelfForceHeader(out);
            runSamplesToStream(out, m, graph, bandCells, testId, samples, abls, statsByAbl, wrapPos, embed);
        };

        // Uniform scaling subtest (use uniform mesh, no refine boxes).
        if (wantSubtest("uniform"))
        {
            std::vector<int> Nlist = sf.uniform.Nlist;
            if (Nlist.empty())
                Nlist.push_back(dictMesh.Nx);

            std::ofstream out(sf.outputUniform);
            if (!out)
            {
                std::cerr << "[PicVerification][Test3] ERROR: failed to open " << sf.outputUniform << std::endl;
            }
            else
            {
                writeSelfForceHeader(out);
            }

            for (int N : Nlist)
            {
                DictionaryMesh dm = dictMesh;
                dm.Nx = N;
                dm.Ny = N;
                dm.refineBoxes.clear();

                QuadtreeMesh umesh(dm);
                umesh.output(std::string("self_force_uniform_N") + std::to_string(N));
                FaceSegmentGraph graph;
                graph.build(umesh);
                const std::vector<BandCellBox> bandCells; // none on uniform mesh

                const std::string testId = std::string("uniform_N") + std::to_string(N);
                const std::vector<std::string> abls = {"A0"};

                if (out)
                {
                    std::map<std::string, SummaryStats> statsTmp;
                    runSamplesToStream(out, umesh, graph, bandCells, testId, sf.uniform.samplePoints, abls, statsTmp, /*wrapPos=*/true, nullptr);

                    SummaryStats s = statsTmp["A0"];
                    uniformStatsByN[testId] = s;
                }
            }
        }

        const bool wantTraverse = wantSubtest("traverse");
        const bool wantTraverseCut = wantSubtest("traversecut");
        const bool wantCutCross = wantSubtest("cutcross");
        if (wantTraverse || wantTraverseCut || wantCutCross)
        {
            std::ofstream out(sf.outputTraverse);
            if (!out)
            {
                std::cerr << "[PicVerification][Test3] ERROR: failed to open " << sf.outputTraverse << std::endl;
            }
            else
            {
                writeSelfForceHeader(out);
                if (wantTraverse)
                {
                    FaceSegmentGraph graph;
                    graph.build(mesh);
                    const auto bandCells = collectCompositeBandCells(mesh);
                    const auto samples = buildPath(sf.traverse.start, sf.traverse.end, sf.traverse.numSamples);
                    runSamplesToStream(out, mesh, graph, bandCells, "traverse", samples, ablations, traverseStatsByAbl, /*wrapPos=*/true, nullptr);
                }

                if (wantTraverseCut || wantCutCross)
                {
                    FaceSegmentGraph graph;
                    auto embed = buildEmbeddedContext(mesh, graph);
                    if (!embed)
                    {
                        std::cerr << "[PicVerification][Test3] WARNING: cut-cell subtest(s) requested but embedded geometry is not enabled/valid." << std::endl;
                    }
                    else
                    {
                        const auto bandCells = collectCompositeBandCells(mesh);
                        if (wantTraverseCut)
                        {
                            const auto samples = buildPath(sf.traverseCut.start, sf.traverseCut.end, sf.traverseCut.numSamples);
                            runSamplesToStream(out, mesh, graph, bandCells, "traverseCut", samples, ablations, traverseCutStatsByAbl, /*wrapPos=*/true, embed.get());
                        }
                        if (wantCutCross)
                        {
                            const auto samples = buildPath(sf.cutCross.start, sf.cutCross.end, sf.cutCross.numSamples);
                            runSamplesToStream(out, mesh, graph, bandCells, "cutCross", samples, ablations, cutCrossStatsByAbl, /*wrapPos=*/true, embed.get());
                        }
                    }
                }
            }
        }

        // Periodic crossing sanity subtest.
        if (wantSubtest("periodiccross"))
        {
            FaceSegmentGraph graph;
            graph.build(mesh);
            const auto bandCells = collectCompositeBandCells(mesh);
            const auto samples = buildPath(sf.periodicCross.start, sf.periodicCross.end, sf.periodicCross.numSamples);
            runSamplesToFile(sf.outputPeriodicCross, mesh, graph, bandCells, "periodicCross", samples, ablations, crossStatsByAbl, /*wrapPos=*/true, nullptr);
        }

        // Write summary JSON (minimal, stable).
        std::ofstream summary(sf.outputSummary);
        if (summary)
        {
            summary << std::setprecision(12);
            summary << "{\n";
            summary << "  \"uniform\": {\n";
            bool firstUniform = true;
            for (const auto& kv : uniformStatsByN)
            {
                const SummaryStats& s = kv.second;
                const double rmsFstar = (s.n > 0) ? std::sqrt(s.sumFstar2 / static_cast<double>(s.n)) : 0.0;
                const double rmsFsi = (s.n > 0) ? std::sqrt(s.sumFsi2 / static_cast<double>(s.n)) : 0.0;
                const double rmsRGhost = (s.n > 0) ? std::sqrt(s.sumRGhost2 / static_cast<double>(s.n)) : 0.0;
                if (!firstUniform)
                    summary << ",\n";
                firstUniform = false;
                summary << "    \"" << kv.first << "\": { \"maxFstar\": " << s.maxFstar
                        << ", \"rmsFstar\": " << rmsFstar
                        << ", \"maxF_SI\": " << s.maxFsi
                        << ", \"rmsF_SI\": " << rmsFsi
                        << ", \"maxAbsSumQAfter\": " << s.maxAbsSumQAfter
                        << ", \"maxAbsPhiMean\": " << s.maxAbsPhiMean
                        << ", \"maxAbsRGhost\": " << s.maxAbsRGhost
                        << ", \"rmsRGhost\": " << rmsRGhost << " }";
            }
            summary << "\n  },\n";
            summary << "  \"traverse\": {\n";
            bool firstTrav = true;
            for (const auto& kv : traverseStatsByAbl)
            {
                const SummaryStats& s = kv.second;
                const double rmsFstar = (s.n > 0) ? std::sqrt(s.sumFstar2 / static_cast<double>(s.n)) : 0.0;
                const double rmsRGhost = (s.n > 0) ? std::sqrt(s.sumRGhost2 / static_cast<double>(s.n)) : 0.0;
                if (!firstTrav)
                    summary << ",\n";
                firstTrav = false;
                summary << "    \"" << kv.first << "\": { \"peakFstar\": " << s.maxFstar
                        << ", \"peakF_SI\": " << s.maxFsi
                        << ", \"rmsFstar\": " << rmsFstar
                        << ", \"maxAbsSumQAfter\": " << s.maxAbsSumQAfter
                        << ", \"maxAbsPhiMean\": " << s.maxAbsPhiMean
                        << ", \"maxAbsRGhost\": " << s.maxAbsRGhost
                        << ", \"rmsRGhost\": " << rmsRGhost << " }";
            }
            summary << "\n  },\n";
            summary << "  \"traverseCut\": {\n";
            bool firstTravCut = true;
            for (const auto& kv : traverseCutStatsByAbl)
            {
                const SummaryStats& s = kv.second;
                const double rmsFstar = (s.n > 0) ? std::sqrt(s.sumFstar2 / static_cast<double>(s.n)) : 0.0;
                const double rmsRGhost = (s.n > 0) ? std::sqrt(s.sumRGhost2 / static_cast<double>(s.n)) : 0.0;
                if (!firstTravCut)
                    summary << ",\n";
                firstTravCut = false;
                summary << "    \"" << kv.first << "\": { \"peakFstar\": " << s.maxFstar
                        << ", \"peakF_SI\": " << s.maxFsi
                        << ", \"rmsFstar\": " << rmsFstar
                        << ", \"maxAbsSumQAfter\": " << s.maxAbsSumQAfter
                        << ", \"maxAbsPhiMean\": " << s.maxAbsPhiMean
                        << ", \"maxAbsRGhost\": " << s.maxAbsRGhost
                        << ", \"rmsRGhost\": " << rmsRGhost << " }";
            }
            summary << "\n  },\n";
            summary << "  \"cutCross\": {\n";
            bool firstCutCross = true;
            for (const auto& kv : cutCrossStatsByAbl)
            {
                const SummaryStats& s = kv.second;
                const double rmsFstar = (s.n > 0) ? std::sqrt(s.sumFstar2 / static_cast<double>(s.n)) : 0.0;
                const double rmsRGhost = (s.n > 0) ? std::sqrt(s.sumRGhost2 / static_cast<double>(s.n)) : 0.0;
                if (!firstCutCross)
                    summary << ",\n";
                firstCutCross = false;
                summary << "    \"" << kv.first << "\": { \"peakFstar\": " << s.maxFstar
                        << ", \"peakF_SI\": " << s.maxFsi
                        << ", \"rmsFstar\": " << rmsFstar
                        << ", \"maxAbsSumQAfter\": " << s.maxAbsSumQAfter
                        << ", \"maxAbsPhiMean\": " << s.maxAbsPhiMean
                        << ", \"maxAbsRGhost\": " << s.maxAbsRGhost
                        << ", \"rmsRGhost\": " << rmsRGhost << " }";
            }
            summary << "\n  },\n";
            summary << "  \"periodicCross\": {\n";
            bool firstCross = true;
            for (const auto& kv : crossStatsByAbl)
            {
                const SummaryStats& s = kv.second;
                const double rmsFstar = (s.n > 0) ? std::sqrt(s.sumFstar2 / static_cast<double>(s.n)) : 0.0;
                const double rmsRGhost = (s.n > 0) ? std::sqrt(s.sumRGhost2 / static_cast<double>(s.n)) : 0.0;
                if (!firstCross)
                    summary << ",\n";
                firstCross = false;
                summary << "    \"" << kv.first << "\": { \"peakFstar\": " << s.maxFstar
                        << ", \"peakF_SI\": " << s.maxFsi
                        << ", \"rmsFstar\": " << rmsFstar
                        << ", \"maxAbsSumQAfter\": " << s.maxAbsSumQAfter
                        << ", \"maxAbsPhiMean\": " << s.maxAbsPhiMean
                        << ", \"maxAbsRGhost\": " << s.maxAbsRGhost
                        << ", \"rmsRGhost\": " << rmsRGhost << " }";
            }
            summary << "\n  }\n";
            summary << "}\n";
        }
    }

    const int rc = harness.process(regMetrics);
    std::cout << "[PicVerification] Done." << std::endl;
    return rc;
}
