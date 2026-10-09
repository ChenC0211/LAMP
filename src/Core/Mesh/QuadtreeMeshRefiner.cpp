#include "Core/Mesh/QuadtreeMeshRefiner.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "Core/Geometry/BSpline.h"

namespace
{
int sourceLeafIndex(const QuadtreeMesh& mesh, Vector2D point)
{
    const int label = mesh.findCellLabel(point);
    if (label < 0 || label >= mesh.cellLabelCapacity())
        throw std::invalid_argument("refinement target point is outside source mesh");
    const int index = mesh.cells(label).index();
    if (index < 0 || index >= mesh.cellIndexNum())
        throw std::invalid_argument("refinement source lookup did not return a leaf");
    return index;
}

std::string toLowerCopy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

static double pointSegmentDist2(const Vector2D& p, const Vector2D& a, const Vector2D& b)
{
    const double vx = b.x - a.x;
    const double vy = b.y - a.y;
    const double wx = p.x - a.x;
    const double wy = p.y - a.y;
    const double c1 = vx * wx + vy * wy;
    if (c1 <= 0.0)
    {
        const double dx = p.x - a.x;
        const double dy = p.y - a.y;
        return dx * dx + dy * dy;
    }
    const double c2 = vx * vx + vy * vy;
    if (c2 <= c1)
    {
        const double dx = p.x - b.x;
        const double dy = p.y - b.y;
        return dx * dx + dy * dy;
    }
    const double t = c1 / c2;
    const double projx = a.x + t * vx;
    const double projy = a.y + t * vy;
    const double dx = p.x - projx;
    const double dy = p.y - projy;
    return dx * dx + dy * dy;
}

static double pointPolylineDist2(const Vector2D& p, const std::vector<Vector2D>& poly)
{
    if (poly.size() < 2)
        return std::numeric_limits<double>::infinity();
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t i = 1; i < poly.size(); ++i)
    {
        const double d2 = pointSegmentDist2(p, poly[i - 1], poly[i]);
        if (d2 < best)
            best = d2;
    }
    return best;
}

static bool aabbOverlaps(const Vector2D& aMin, const Vector2D& aMax,
                         const Vector2D& bMin, const Vector2D& bMax)
{
    if (aMax.x < bMin.x || bMax.x < aMin.x)
        return false;
    if (aMax.y < bMin.y || bMax.y < aMin.y)
        return false;
    return true;
}

static std::vector<Vector2D> discretizeGeometryLine(const Line& L, int numPoints = 32)
{
    std::vector<Vector2D> poly;
    numPoints = std::max(1, numPoints);
    poly.reserve(static_cast<std::size_t>(numPoints) + 1u);

    switch (L.type())
    {
        case LineType::straight:
        {
            poly.push_back(L.points[0]);
            poly.push_back(L.points[1]);
            break;
        }
        case LineType::arc:
        {
            const Vector2D start = L.points[0];
            const double radius = L.radius();
            const double radian = L.radian();
            const Vector2D center = L.arcCenter();
            const double startAngle = std::atan2(start.y - center.y, start.x - center.x);
            for (int i = 0; i <= numPoints; ++i)
            {
                const double t = static_cast<double>(i) / static_cast<double>(numPoints);
                const double angle = startAngle + t * radian * Pi;
                poly.push_back(center + Vector2D(radius * std::cos(angle), radius * std::sin(angle)));
            }
            break;
        }
        case LineType::spline:
        {
            BSplineCurve bspline(L.featurePoints(), 3);
            for (int i = 0; i <= numPoints; ++i)
            {
                const double t = static_cast<double>(i) / static_cast<double>(numPoints);
                poly.push_back(bspline.evaluate(t));
            }
            break;
        }
        default:
            poly.push_back(L.points[0]);
            poly.push_back(L.points[1]);
            break;
    }

    return poly;
}
} // namespace

QuadtreeMeshRefiner::QuadtreeMeshRefiner(const DictionaryMesh& dictMesh)
    : dictMesh_(dictMesh)
{
}

bool QuadtreeMeshRefiner::attachGeometry(const DictionaryGeometry* dictGeometry,
                                         ImmersedGeometry<QuadtreeMesh>* geometry)
{
    dictGeometry_ = dictGeometry;
    geometry_ = geometry;
    geoBoundaryRules_.clear();
    geoRegionRules_.clear();
    geoPhaseRules_.clear();
    regionActive_ = false;
    boundaryActive_ = false;
    phaseActive_ = false;
    geometryActive_ = false;

    if (!dictMesh_.geometryRefinement.enabled)
        return false;

    for (const auto& rule : dictMesh_.geometryRefinement.regionRules)
    {
        GeoRegionRuleRuntime rr;
        rr.name = rule.name;
        rr.shape = rule.shape;
        rr.min = rule.min;
        rr.max = rule.max;
        rr.maxLevel = rule.maxLevel;
        geoRegionRules_.push_back(rr);
    }
    regionActive_ = !geoRegionRules_.empty();
    geometryActive_ = regionActive_;

    const bool needsGeometryContext =
        !dictMesh_.geometryRefinement.boundaryRules.empty() ||
        !dictMesh_.geometryRefinement.phaseRules.empty();
    if (!needsGeometryContext)
        return geometryActive_;

    if (!dictGeometry_ || !geometry_)
    {
        if (!warnedMissingGeometry_)
        {
            std::cerr << "[MeshRefiner] WARNING: GeometryRefinement enabled but geometry is unavailable.\n";
            warnedMissingGeometry_ = true;
        }
        return geometryActive_;
    }

    const int nLS = static_cast<int>(geometry_->levelSets.size());
    if (nLS <= 0)
    {
        if (!warnedMissingGeometry_)
        {
            std::cerr << "[MeshRefiner] WARNING: GeometryRefinement enabled but no level sets are defined.\n";
            warnedMissingGeometry_ = true;
        }
        return geometryActive_;
    }

    for (const auto& rule : dictMesh_.geometryRefinement.boundaryRules)
    {
        auto it = dictGeometry_->boundaryById.find(rule.boundaryId);
        if (it == dictGeometry_->boundaryById.end())
        {
            std::cerr << "[MeshRefiner] WARNING: refineBoundary boundaryId "
                      << rule.boundaryId << " not found in geometry.dict; skipping.\n";
            continue;
        }
        const auto& ent = it->second;
        if (ent.kind == GeometricBoundary::domain)
        {
            std::cerr << "[MeshRefiner] WARNING: refineBoundary boundaryId "
                      << rule.boundaryId << " is a domain boundary; skipping.\n";
            continue;
        }
        if (ent.levelSetIndex < 0 || ent.levelSetIndex >= nLS)
        {
            std::cerr << "[MeshRefiner] WARNING: refineBoundary boundaryId "
                      << rule.boundaryId << " has invalid levelSetIndex; skipping.\n";
            continue;
        }

        GeoBoundaryRuleRuntime rr;
        rr.boundaryId = rule.boundaryId;
        rr.levelSetIndex = ent.levelSetIndex;
        rr.cGamma = (rule.cGamma > 0.0) ? rule.cGamma : 1.0;
        rr.maxLevel = rule.maxLevel;
        rr.lineIds = ent.lineIds;
        rr.linePolylines.reserve(rr.lineIds.size());
        for (int lineId : rr.lineIds)
        {
            if (lineId < 0 || lineId >= static_cast<int>(dictGeometry_->lines.size()))
                continue;
            std::vector<Vector2D> poly = discretizeGeometryLine(
                dictGeometry_->lines[static_cast<std::size_t>(lineId)]);
            if (poly.size() >= 2u)
                rr.linePolylines.push_back(std::move(poly));
        }
        geoBoundaryRules_.push_back(rr);
    }
    boundaryActive_ = !geoBoundaryRules_.empty();

    for (const auto& rule : dictMesh_.geometryRefinement.phaseRules)
    {
        if (rule.phase <= 0)
            continue;
        DictionaryMesh::GeometryRefinementConfig::PhaseRule pr = rule;
        pr.phase = rule.phase - 1; // convert to 0-based
        geoPhaseRules_[pr.phase] = pr;
    }
    phaseActive_ = !geoPhaseRules_.empty();
    geometryActive_ = regionActive_ || boundaryActive_ || phaseActive_;
    return geometryActive_;
}

void QuadtreeMeshRefiner::clearGeometry()
{
    dictGeometry_ = nullptr;
    geometry_ = nullptr;
    geoBoundaryRules_.clear();
    geoRegionRules_.clear();
    geoPhaseRules_.clear();
    regionActive_ = false;
    boundaryActive_ = false;
    phaseActive_ = false;
    geometryActive_ = false;
}

void QuadtreeMeshRefiner::refreshBoundaryPolylines()
{
    if (!boundaryActive_ || !geometry_)
        return;

    for (auto& rule : geoBoundaryRules_)
    {
        rule.linePolylines.clear();
        for (int lineId : rule.lineIds)
        {
            if (lineId < 0 || lineId >= static_cast<int>(geometry_->lines.size()))
                continue;
            std::vector<Vector2D> poly = discretizeGeometryLine(
                geometry_->lines[static_cast<std::size_t>(lineId)]);
            if (poly.size() >= 2u)
                rule.linePolylines.push_back(std::move(poly));
        }
    }
}

void QuadtreeMeshRefiner::clearBoundaryPolylines()
{
    for (auto& rule : geoBoundaryRules_)
        rule.linePolylines.clear();
}

void QuadtreeMeshRefiner::replaceBoundaryPolylines(int boundaryId,
    const std::vector<std::vector<Vector2D>>& polylines)
{
    std::vector<GeoBoundaryRuleRuntime*> selected;
    for (auto& rule : geoBoundaryRules_)
        if (rule.boundaryId == boundaryId) selected.push_back(&rule);
    if (selected.empty())
        throw std::invalid_argument("current boundary polylines require an attached GeometryRefinement boundary rule");
    if (polylines.empty())
        throw std::invalid_argument("current boundary refinement polylines must not be empty");
    for (const auto& line : polylines)
    {
        if (line.size() < 2)
            throw std::invalid_argument("a boundary refinement polyline needs at least two points");
        for (std::size_t i = 0; i < line.size(); ++i)
        {
            if (!std::isfinite(line[i].x) || !std::isfinite(line[i].y))
                throw std::invalid_argument("non-finite current boundary refinement coordinate");
            if (i > 0 && line[i].x == line[i-1].x && line[i].y == line[i-1].y)
                throw std::invalid_argument("degenerate current boundary refinement segment");
        }
    }
    // Stage all allocations, including duplicate rules for the same boundary,
    // before publishing. A rejected update cannot partially move a boundary.
    std::vector<std::vector<std::vector<Vector2D>>> staged(selected.size(), polylines);
    for (std::size_t i = 0; i < selected.size(); ++i)
        selected[i]->linePolylines.swap(staged[i]);
}

bool QuadtreeMeshRefiner::buildMappedTargetLevels(const QuadtreeMesh& target, const CurrentSource& source,
    CellField<double>& targetLevel, CellField<double>* indicatorOut)
{
    const auto& mesh = source.mesh;
    if (target.minCoord.x != mesh.minCoord.x || target.minCoord.y != mesh.minCoord.y
        || target.maxCoord.x != mesh.maxCoord.x || target.maxCoord.y != mesh.maxCoord.y
        || source.protectedCells.size() != static_cast<std::size_t>(mesh.cellIndexNum()))
        throw std::invalid_argument("mapped refinement requires one domain and source-sized protection");
    CurrentFields current;
    current.cellLabels.assign(target.cellIndexToLabel, target.cellIndexToLabel + target.cellIndexNum());
    current.phaseByCell.assign(target.cellIndexNum(), -1);
    for (const auto& field : source.indicators)
    {
        if (field.cellValue.size() != static_cast<std::size_t>(mesh.cellIndexNum()))
            throw std::invalid_argument("mapped refinement indicator does not match source leaves");
        for (double value : field.cellValue)
            if (!std::isfinite(value) || value < 0.0)
                throw std::invalid_argument("mapped refinement source indicator must be finite and nonnegative");
        // Copy only metadata; the source's full field is never duplicated.
        current.indicators.push_back({field.name, std::vector<double>(target.cellIndexNum(), 0.0),
            field.scale, field.reference, field.percentile});
    }
    for (int i = 0; i < target.cellIndexNum(); ++i)
    {
        const auto point = target.cells(target.cellIndexToLabel[i]).centerCoord();
        const int old = sourceLeafIndex(mesh, point);
        if (source.phaseAt) current.phaseByCell[i] = source.phaseAt(point, old);
        for (std::size_t k = 0; k < current.indicators.size(); ++k)
            current.indicators[k].cellValue[i] = source.indicators[k].cellValue[old];
    }
    // A prospective parent center alone can miss a peak on another old leaf.
    if (!current.indicators.empty())
        for (int old = 0; old < mesh.cellIndexNum(); ++old)
        {
            const int i = sourceLeafIndex(target, mesh.cells(mesh.cellIndexToLabel[old]).centerCoord());
            for (std::size_t k = 0; k < current.indicators.size(); ++k)
                current.indicators[k].cellValue[i] = std::max(current.indicators[k].cellValue[i],
                    source.indicators[k].cellValue[old]);
        }
    buildTargetLevels(target, targetLevel, indicatorOut, nullptr, true, &current);
    bool requested = false;
    for (int i = 0; i < target.cellIndexNum(); ++i)
    {
        const auto& cell = target.cells(target.cellIndexToLabel[i]);
        const int old = sourceLeafIndex(mesh, cell.centerCoord());
        if (source.protectedCells[old])
        {
            if (cell.level() != mesh.cells(mesh.cellIndexToLabel[old]).level())
                throw std::invalid_argument("2:1 balancing changed a protected leaf; enlarge configured fine bands");
            targetLevel.value[i] = cell.level();
        }
        requested = requested || targetLevel.value[i] != cell.level();
    }
    return requested && target.mayAdaptToTargets(targetLevel);
}

bool QuadtreeMeshRefiner::buildTargetLevels(const QuadtreeMesh& mesh,
                                            CellField<double>& targetLevel,
                                            CellField<double>* indicatorOut,
                                            const BudgetOptions* budget,
                                            bool coarsenOutsideGeometry,
                                            const CurrentFields* current,
                                            const std::vector<RefinementIndicatorField>* picIndicators)
{
    if (current && picIndicators)
        throw std::invalid_argument("flow and PIC refinement inputs must not be combined");
    if (coarsenOutsideGeometry && (!geometryActive_ || (dictMesh_.physicsRefinement.enabled && !current)))
        throw std::invalid_argument("runtime coarsening requires active geometry and explicit current physics fields when enabled");
    const int nCells = mesh.cellIndexNum();
    if (nCells <= 0)
        return false;

    const auto currentCellLabels = [&](const CellField<double>& field)
    {
        if (field.size() != nCells) return false;
        for (int i = 0; i < nCells; ++i)
            if (field.label(i) != mesh.cellIndexToLabel[i]) return false;
        return true;
    };
    if (!currentCellLabels(targetLevel))
        targetLevel.initialize(mesh, 0.0);
    if (indicatorOut && !currentCellLabels(*indicatorOut))
        indicatorOut->initialize(mesh, 0.0);

    if (current)
    {
        if (!coarsenOutsideGeometry || static_cast<int>(current->cellLabels.size()) != nCells
            || !std::equal(current->cellLabels.begin(), current->cellLabels.end(), mesh.cellIndexToLabel)
            || static_cast<int>(current->phaseByCell.size()) != nCells)
            throw std::invalid_argument("runtime refinement inputs do not identify the current leaf mesh");
        for (int phase : current->phaseByCell)
            if (phase < -1 || (dictGeometry_ && phase >= static_cast<int>(dictGeometry_->phaseName.size())))
                throw std::invalid_argument("invalid current refinement phase ID");
    }
    const bool flowPhysics = current && dictMesh_.physicsRefinement.enabled;
    if (flowPhysics)
    {
        const auto& config = dictMesh_.physicsRefinement;
        if (current->indicators.size() != config.indicators.size() || config.indicators.empty()
            || !std::isfinite(config.hysteresis) || config.hysteresis <= 0 || config.hysteresis >= 1)
            throw std::invalid_argument("incomplete runtime flow refinement inputs");
        for (std::size_t k = 0; k < config.indicators.size(); ++k)
        {
            const auto& field = current->indicators[k];
            const auto key = toLowerCopy(config.indicators[k]);
            if ((key != "velocity" && key != "velocitygradient") || toLowerCopy(field.name) != key
                || !std::isfinite(field.reference) || field.reference <= 0
                || static_cast<int>(field.cellValue.size()) != nCells || !config.maxLevelCaps.count(key))
                throw std::invalid_argument("missing or invalid current flow refinement indicator");
            for (double value : field.cellValue)
                if (!std::isfinite(value) || value < 0 || !std::isfinite(value / field.reference))
                    throw std::invalid_argument("nonfinite or negative flow refinement indicator");
        }
    }
    else if (current && !current->indicators.empty())
        throw std::invalid_argument("inactive runtime refinement indicators supplied");

    const bool wantPhysics = !current && dictMesh_.physicsRefinement.enabled &&
        !dictMesh_.physicsRefinement.indicators.empty();

    RefinementIndicatorField debyeField;
    RefinementIndicatorField eField;
    bool hasDebye = false;
    bool hasEField = false;

    if (wantPhysics)
    {
        for (const auto& name : dictMesh_.physicsRefinement.indicators)
        {
            const RefinementIndicatorField* supplied = nullptr;
            if (picIndicators)
                for (const auto& candidate : *picIndicators)
                    if (toLowerCopy(candidate.name) == toLowerCopy(name))
                    {
                        supplied = &candidate;
                        break;
                    }
            if (!supplied)
            {
                const std::string key = toLowerCopy(name);
                if (warnedMissingIndicators_.insert(key).second)
                {
                    std::cerr << "[MeshRefiner] WARNING: indicator '" << name
                              << "' not provided by any module.\n";
                }
                continue;
            }

            RefinementIndicatorField field = *supplied;
            const std::string key = toLowerCopy(field.name.empty() ? name : field.name);
            if (key == "debyelength")
            {
                debyeField = std::move(field);
                hasDebye = true;
            }
            else if (key == "electricfield")
            {
                eField = std::move(field);
                hasEField = true;
            }
        }
    }

    if (hasDebye && static_cast<int>(debyeField.cellValue.size()) != nCells)
    {
        const std::string key = "size:debyelength";
        if (warnedMissingIndicators_.insert(key).second)
        {
            std::cerr << "[MeshRefiner] WARNING: DebyeLength indicator size mismatch (expected "
                      << nCells << ", got " << debyeField.cellValue.size() << "); ignoring.\n";
        }
        hasDebye = false;
    }
    if (hasEField && static_cast<int>(eField.cellValue.size()) != nCells)
    {
        const std::string key = "size:electricfield";
        if (warnedMissingIndicators_.insert(key).second)
        {
            std::cerr << "[MeshRefiner] WARNING: ElectricField indicator size mismatch (expected "
                      << nCells << ", got " << eField.cellValue.size() << "); ignoring.\n";
        }
        hasEField = false;
    }

    const bool hasPhysics = hasDebye || hasEField;
    const bool hasGeometry = geometryActive_;

    if (!hasPhysics && !hasGeometry)
    {
        for (int ci = 0; ci < nCells; ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            targetLevel.value[ci] = static_cast<double>(c.level());
            if (indicatorOut)
                indicatorOut->value[ci] = 0.0;
        }
        return false;
    }

    const double cd = hasDebye ? std::max(1e-12, debyeField.scale) : 1.0;
    const double hysteresis = std::max(0.0, dictMesh_.physicsRefinement.hysteresis);
    const double coarsenThreshold = std::max(0.0, 1.0 - hysteresis);

    double Eref = 0.0;
    if (hasEField && eField.reference > 0.0)
        Eref = eField.reference;

    std::vector<double> eSamples;
    if (hasEField && Eref <= 0.0 && eField.percentile >= 0.0)
        eSamples.reserve(nCells / 4);

    if (hasEField && Eref <= 0.0)
    {
        double Emax = 0.0;
        const bool wantPercentile = (eField.percentile >= 0.0);
        for (int ci = 0; ci < nCells; ++ci)
        {
            const double val = eField.cellValue[ci];
            if (!std::isfinite(val) || val <= 0.0)
                continue;
            if (val > Emax)
                Emax = val;
            if (wantPercentile)
                eSamples.push_back(val);
        }
        if (wantPercentile && !eSamples.empty())
        {
            double pct = eField.percentile;
            if (pct < 0.0) pct = 0.0;
            if (pct > 100.0) pct = 100.0;
            const size_t n = eSamples.size();
            const size_t k = static_cast<size_t>(
                std::ceil((pct / 100.0) * static_cast<double>(n - 1)));
            if (k < n)
            {
                std::nth_element(eSamples.begin(), eSamples.begin() + k, eSamples.end());
                Eref = eSamples[k];
            }
        }
        if (!(Eref > 0.0))
            Eref = Emax;
    }

    const int maxLevel = std::max(0, mesh.maxLevel());
    int maxAllowedLevel = maxLevel;
    if (hasPhysics)
    {
        auto capFor = [&](const std::string& name) -> int {
            auto it = dictMesh_.physicsRefinement.maxLevelCaps.find(toLowerCopy(name));
            if (it == dictMesh_.physicsRefinement.maxLevelCaps.end())
                return maxAllowedLevel;
            return std::max(0, it->second);
        };
        if (hasDebye)
            maxAllowedLevel = std::min(maxAllowedLevel, capFor("debyelength"));
        if (hasEField)
            maxAllowedLevel = std::min(maxAllowedLevel, capFor("electricfield"));
    }

    std::vector<double> r_use(nCells, hasPhysics ? 0.0 : 1.0);
    if (hasPhysics)
    {
        for (int ci = 0; ci < nCells; ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            const Vector2D sz = c.size();
            const double h_code = std::min(sz.x, sz.y);

        double r_debye = 0.0;
        bool debyeValid = false;
        if (hasDebye)
        {
            const double lambdaD = debyeField.cellValue[ci];
            if (std::isfinite(lambdaD) && lambdaD > 0.0)
            {
                r_debye = h_code / (cd * lambdaD);
                debyeValid = true;
            }
        }

        double eScore = 0.0;
        bool eValid = false;
        if (hasEField && Eref > 0.0)
        {
            const double val = eField.cellValue[ci];
            if (std::isfinite(val) && val > 0.0)
            {
                eScore = val / Eref;
                eValid = true;
            }
        }

        if (!debyeValid && !eValid)
        {
            r_use[ci] = 1.0;
            continue;
        }

        double r = 0.0;
        if (debyeValid)
            r = r_debye;
        if (eValid)
        {
            if (debyeValid)
                r = std::max(r, r_debye * eScore);
            else
                r = std::max(r, eScore);
        }
            r_use[ci] = r;
        }
    }

    double refineCut = -1.0;
    double coarsenCut = -1.0;
    if (hasPhysics && budget && budget->cellBudget > 0)
    {
        std::vector<double> refineScores;
        std::vector<double> coarsenScores;
        refineScores.reserve(nCells / 8);
        coarsenScores.reserve(nCells / 8);

        for (int ci = 0; ci < nCells; ++ci)
        {
            if (budget->excludeMask && ci < static_cast<int>(budget->excludeMask->size()) &&
                (*budget->excludeMask)[ci] != 0)
                continue;

            const double r = r_use[ci];
            if (r > 1.0)
                refineScores.push_back(r);
            if (r < coarsenThreshold)
            {
                const double score = (r > 0.0) ? (1.0 / r) : 1.0e30;
                coarsenScores.push_back(score);
            }
        }

        const int budgetCount = std::max(0, budget->cellBudget);
        if (static_cast<int>(refineScores.size()) > budgetCount)
        {
            auto nth = refineScores.begin() + budgetCount;
            std::nth_element(refineScores.begin(), nth, refineScores.end(),
                             [](double a, double b) { return a > b; });
            refineCut = *nth;
        }
        if (static_cast<int>(coarsenScores.size()) > budgetCount)
        {
            auto nth = coarsenScores.begin() + budgetCount;
            std::nth_element(coarsenScores.begin(), nth, coarsenScores.end(),
                             [](double a, double b) { return a > b; });
            coarsenCut = *nth;
        }
    }

    auto cellPhase = [&](const Cell& c) -> int {
        if (current) return current->phaseByCell[c.index()];
        if (!phaseActive_ || !geometry_)
            return -1;
        int p = -1;
        for (int k = 0; k < 4; ++k)
        {
            const Vert& v = mesh.verts(c.cellToVert[k]);
            const int vp = v.phase();
            if (vp < 0)
                return -1;
            if (p < 0)
                p = vp;
            else if (p != vp)
                return -1;
        }
        return p;
    };

    for (int ci = 0; ci < nCells; ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        const int currentLevel = c.level();
        const double r = r_use[ci];

        // The same ratio/hysteresis criterion as the existing refiner, with
        // independent per-indicator caps. Geometry minima are combined below.
        // One level of coarsening per sweep; zero flow is valid, not missing.
        int flowAim = 0;
        double flowScore = 0;
        if (flowPhysics)
            for (const auto& field : current->indicators)
            {
                const double score = field.cellValue[ci] / field.reference;
                const int cap = dictMesh_.physicsRefinement.maxLevelCaps.at(toLowerCopy(field.name));
                int aim = currentLevel - 1;
                if (score > 1)
                    aim = currentLevel + static_cast<int>(std::min<double>(maxLevel + 1, std::ceil(std::log2(score))));
                else if (score >= coarsenThreshold) aim = currentLevel;
                flowAim = std::max(flowAim, std::min(cap, aim));
                flowScore = std::max(flowScore, score);
            }

        bool wantRefine = (r > 1.0);
        bool wantCoarsen = coarsenOutsideGeometry || (r < coarsenThreshold);
        bool allowRefine = wantRefine;
        bool allowCoarsen = wantCoarsen;

        if (refineCut > 0.0 && wantRefine && r < refineCut)
            allowRefine = false;
        if (coarsenCut > 0.0 && wantCoarsen)
        {
            const double score = (r > 0.0) ? (1.0 / r) : 1.0e30;
            if (score < coarsenCut)
                allowCoarsen = false;
        }

        int geomRegionCap = -1;
        if (regionActive_)
        {
            const Vector2D cc = c.centerCoord();
            const Vector2D sz = c.size();
            const Vector2D cellMin(cc.x - 0.5 * sz.x, cc.y - 0.5 * sz.y);
            const Vector2D cellMax(cc.x + 0.5 * sz.x, cc.y + 0.5 * sz.y);

            // User-defined region refinement enforces a minimum level inside
            // configured geometric boxes without replacing physics/boundary AMR.
            for (const auto& rule : geoRegionRules_)
            {
                if (rule.shape != DictionaryMesh::GeometryRefinementConfig::RegionShape::kBox)
                    continue;
                if (!aabbOverlaps(cellMin, cellMax, rule.min, rule.max))
                    continue;

                const int ruleMax = (rule.maxLevel >= 0) ? rule.maxLevel : maxLevel;
                if (geomRegionCap < 0)
                    geomRegionCap = ruleMax;
                else
                    geomRegionCap = std::max(geomRegionCap, ruleMax);
            }
        }

        if (geomRegionCap >= 0)
        {
            geomRegionCap = std::max(0, std::min(maxLevel, geomRegionCap));
            if (currentLevel < geomRegionCap)
                allowRefine = true;
            if (currentLevel <= geomRegionCap)
                allowCoarsen = false;
        }

        bool inGeomBand = false;
        int geomBandCap = -1;
        if (boundaryActive_ && geometry_)
        {
            for (const auto& rule : geoBoundaryRules_)
            {
                double minDist = std::numeric_limits<double>::infinity();
                bool usedBoundaryLines = false;
                if (!rule.linePolylines.empty())
                {
                    usedBoundaryLines = true;
                    for (const auto& poly : rule.linePolylines)
                    {
                        for (int k = 0; k < 4; ++k)
                        {
                            const Vert& v = mesh.verts(c.cellToVert[k]);
                            const double d2 = pointPolylineDist2(v.coord(), poly);
                            if (d2 < minDist)
                                minDist = d2;
                        }
                    }
                }

                double minAbsPsi = std::numeric_limits<double>::max();
                if (!usedBoundaryLines)
                {
                    if (rule.levelSetIndex < 0 ||
                        rule.levelSetIndex >= static_cast<int>(geometry_->levelSets.size()))
                        continue;
                    for (int k = 0; k < 4; ++k)
                    {
                        const Vert& v = mesh.verts(c.cellToVert[k]);
                        const double psi = geometry_->levelSets[rule.levelSetIndex].value[v.index()];
                        minAbsPsi = std::min(minAbsPsi, std::abs(psi));
                    }
                }
                else
                {
                    if (std::isfinite(minDist))
                        minAbsPsi = std::sqrt(minDist);
                }

                const double h_band = c.size().mod();
                const double band = rule.cGamma * h_band;
                // A parent is twice as wide as its children. In runtime
                // AMR, keep a child until it is outside the
                // PARENT refinement band; otherwise coarsening immediately
                // puts the parent back inside the next sweep's refine band.
                // The legacy/PIC target policy is unchanged when opt-out.
                const int ruleMax = (rule.maxLevel >= 0) ? rule.maxLevel : maxLevel;
                if (coarsenOutsideGeometry && currentLevel <= ruleMax && minAbsPsi <= 2.0 * band)
                    allowCoarsen = false;
                if (minAbsPsi <= band)
                {
                    inGeomBand = true;
                    if (geomBandCap < 0)
                        geomBandCap = ruleMax;
                    else
                        geomBandCap = std::max(geomBandCap, ruleMax);
                }
            }
        }

        if (inGeomBand)
        {
            geomBandCap = std::max(0, std::min(maxLevel, geomBandCap));
            if (currentLevel < geomBandCap)
                allowRefine = true;
            if (currentLevel <= geomBandCap)
                allowCoarsen = false;
            if (currentLevel > geomBandCap)
            {
                allowRefine = false;
                allowCoarsen = true;
            }
        }

        int phaseMin = -1;
        int phaseMax = -1;
        if (phaseActive_)
        {
            const int p = cellPhase(c);
            if (p >= 0)
            {
                auto it = geoPhaseRules_.find(p);
                if (it != geoPhaseRules_.end())
                {
                    phaseMin = it->second.minLevel;
                    phaseMax = it->second.maxLevel;
                }
            }
        }
        if (phaseMin >= 0 && phaseMax >= 0 && phaseMin > phaseMax)
            std::swap(phaseMin, phaseMax);
        if (phaseMin >= 0 && currentLevel < phaseMin)
        {
            allowRefine = true;
            allowCoarsen = false;
        }
        if (phaseMax >= 0 && currentLevel > phaseMax)
        {
            allowRefine = false;
            allowCoarsen = true;
        }

        int aimLevel = currentLevel;
        if (allowRefine)
        {
            double dl = std::ceil(std::log2(std::max(1.0, r)));
            if (!std::isfinite(dl))
                dl = 0.0;
            const int delta = std::max(1, static_cast<int>(dl));
            aimLevel = currentLevel + delta;
        }
        else if (allowCoarsen)
        {
            const double denom = std::max(1e-12, coarsenThreshold);
            double dl = std::ceil(std::log2(std::max(1e-12, r / denom)));
            if (!std::isfinite(dl))
                dl = 0.0;
            const int delta = coarsenOutsideGeometry ? -1 : static_cast<int>(dl);
            aimLevel = currentLevel + delta;
        }

        if (hasPhysics)
            aimLevel = std::min(aimLevel, maxAllowedLevel);

        if (flowPhysics) aimLevel = std::max(aimLevel, flowAim);

        if (inGeomBand)
        {
            if (currentLevel < geomBandCap)
                aimLevel = std::max(aimLevel, currentLevel + 1);
            aimLevel = std::min(aimLevel, geomBandCap);
        }
        if (geomRegionCap >= 0)
            aimLevel = std::max(aimLevel, geomRegionCap);
        if (phaseMin >= 0)
            aimLevel = std::max(aimLevel, phaseMin);
        if (phaseMax >= 0)
            aimLevel = std::min(aimLevel, phaseMax);

        aimLevel = std::max(0, std::min(maxLevel, aimLevel));

        targetLevel.value[ci] = static_cast<double>(aimLevel);
        if (indicatorOut)
            indicatorOut->value[ci] = flowPhysics ? flowScore : r;
    }

    return true;
}
