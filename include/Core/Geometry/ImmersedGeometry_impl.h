#include "Core/IO/FileManagement.h"
#pragma once

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <optional>
#include <stdexcept>
#include "Core/Base/Base.h"
#include "Core/Discrete/FiniteDifference.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Geometry/ExactCurveDistance.h"
#include "Core/Base/Interpolation.h"
#include "Core/Base/Structures.h"
#include "Core/Base/Utility.h"
#include "Core/Base/Vector2D.h"
#include "Core/IO/ConsoleControl.h"
using namespace std;

namespace ImmersedGeometryCache
{
constexpr char kMagic[4] = {'I', 'G', 'C', '2'};
constexpr std::uint32_t kVersion = 2;

inline bool MagicEquals(const char a[4], const char b[4])
{
    return std::memcmp(a, b, 4) == 0;
}
} // namespace ImmersedGeometryCache

namespace ImmersedGeometryDetail
{
inline double clamp01(const double t)
{
    if (t < 0.0) return 0.0;
    if (t > 1.0) return 1.0;
    return t;
}

inline double pointSegmentDist2(const Vector2D& p, const Vector2D& a, const Vector2D& b)
{
    const Vector2D ab = b - a;
    const Vector2D ap = p - a;
    const double denom = ab.x * ab.x + ab.y * ab.y;
    if (!(denom > 0.0))
        return (p - a).mod2();

    const double t = clamp01((ap.x * ab.x + ap.y * ab.y) / denom);
    const Vector2D q = a + ab * t;
    return (p - q).mod2();
}

inline double pointPolylineDist2(const Vector2D& p, const std::vector<Vector2D>& poly)
{
    if (poly.empty())
        return std::numeric_limits<double>::infinity();
    if (poly.size() == 1)
        return (p - poly.front()).mod2();

    double best = std::numeric_limits<double>::infinity();
    for (std::size_t i = 1; i < poly.size(); ++i)
        best = std::min(best, pointSegmentDist2(p, poly[i - 1], poly[i]));
    return best;
}

// The sign path and the distance primitive must describe the same curve.
// Straight lines are exact; splines retain their existing sampled-polyline
// representation. Arcs are never chordal approximations in this path.
struct InitialCurve
{
    std::optional<ExactCurveDistance::CircularArc> arc;
    std::vector<Vector2D> polyline;

    Vector2D start() const { return arc ? arc->start() : polyline.front(); }
    Vector2D end() const { return arc ? arc->end() : polyline.back(); }
    double distance(const Vector2D& point) const
    {
        return arc ? arc->distance(point) : std::sqrt(pointPolylineDist2(point, polyline));
    }
    Vector2D inwardEndpointTangent(bool atEnd) const
    {
        if (arc) return arc->inwardEndpointTangent(atEnd);
        if (polyline.size() < 2)
            throw std::invalid_argument("[ImmersedGeometry] point primitive has no endpoint tangent");
        const Vector2D delta = atEnd ? polyline[polyline.size() - 2] - polyline.back() : polyline[1] - polyline.front();
        const double magnitude = std::hypot(delta.x, delta.y);
        if (!(magnitude > 0.0) || !std::isfinite(magnitude))
            throw std::invalid_argument("[ImmersedGeometry] degenerate endpoint segment has no tangent");
        return Vector2D(delta.x / magnitude, delta.y / magnitude);
    }
    bool rayCrossingParity(const Vector2D& point) const
    {
        if (arc) return arc->rayCrossingParity(point);
        bool parity = false;
        for (std::size_t i = 1; i < polyline.size(); ++i)
            parity ^= ExactCurveDistance::segmentRayCrossingParity(point, polyline[i - 1], polyline[i]);
        return parity;
    }
};

inline bool initialCurveParity(const std::vector<InitialCurve>& curves, Vector2D point,
    Vector2D lower, Vector2D upper, bool explicitClosure)
{
    if (explicitClosure)
    {
        if (point.x == lower.x) point.x = std::nextafter(point.x, upper.x);
        else if (point.x == upper.x) point.x = std::nextafter(point.x, lower.x);
        if (point.y == lower.y) point.y = std::nextafter(point.y, upper.y);
        else if (point.y == upper.y) point.y = std::nextafter(point.y, lower.y);
    }
    bool inside = false;
    for (std::size_t i = 0; i < curves.size(); ++i)
    {
        inside ^= curves[i].rayCrossingParity(point);
        inside ^= ExactCurveDistance::segmentRayCrossingParity(
            point, curves[i].end(), curves[(i + 1) % curves.size()].start());
    }
    return inside;
}

// Smooth-minimum of two non-negative distances.  Returns a C-infinity
// approximation of min(a,b) with maximum deviation eps/2 at a==b.
// Formula: smin_eps(a,b) = (a + b - sqrt((a-b)^2 + eps^2)) / 2
inline double smoothMinDist(const double a, const double b, const double eps2)
{
    const double diff = a - b;
    return 0.5 * (a + b - std::sqrt(diff * diff + eps2));
}

// Diagnostic: dump level-set vertex values to a Tecplot-compatible DAT file.
// Format: TITLE / VARIABLES / ZONE header, then x y psi (space-separated).
// Each (levelSetI, tag) pair is tracked with a call counter.  Only the first
// maxDumps invocations per key are written, to capture base mesh through AMR
// iterations without flooding disk on every time step.
template<typename Mesh, typename VertLS>
inline void dumpLevelSetScatter(const Mesh& mesh, const VertLS& ls,
                                int levelSetI, const std::string& tag)
{
    constexpr int maxDumps = 8; // base + up to 7 AMR/step iterations
    // Track call count per (levelSetI, tag).
    static std::map<std::string, int> callCount;
    std::string key = std::to_string(levelSetI) + "_" + tag;
    int& cnt = callCount[key];
    if (cnt >= maxDumps)
        return;
    int callIdx = cnt++;

    const std::filesystem::path outputDir("./output/meta/ls_diag");
    std::filesystem::create_directories(outputDir);
    const std::filesystem::path fname = outputDir /
        ("ls_diag_LS" + std::to_string(levelSetI) + "_" + tag + "_c" + std::to_string(callIdx) + ".dat");
    std::ofstream ofs(fname.string());
    if (!ofs.is_open())
    {
        std::cerr << "[LS-diag] Cannot open " << fname << "\n";
        return;
    }
    const int nVerts = mesh.vertIndexNum();
    ofs << "TITLE = \"LS" << levelSetI << "_" << tag
        << "_c" << callIdx << "\"\n";
    ofs << "VARIABLES = \"x\", \"y\", \"psi\", \"dpsi_dx\", \"dpsi_dy\", \"grad_mag\"\n";
    ofs << "ZONE T=\"LS" << levelSetI << "_" << tag
        << "_c" << callIdx << "\", I=" << nVerts
        << ", F=POINT\n";
    ofs << std::setprecision(12);
    for (int i = 0; i < nVerts; ++i)
    {
        const int vLbl = mesh.vertIndexToLabel[i];
        const auto& v = mesh.verts(vLbl);
        const double gx = FiniteDifference::dx(ls, mesh, v);
        const double gy = FiniteDifference::dy(ls, mesh, v);
        const double gm = std::sqrt(gx * gx + gy * gy);
        ofs << v.coord().x << " " << v.coord().y << " " << ls.value[i]
            << " " << gx << " " << gy << " " << gm << "\n";
    }
    ofs.close();
}

} // namespace ImmersedGeometryDetail

// Public helper: dump the composite plasma-region level set (min of all LS) with the given tag.
// psi_composite[v] = min_i( levelSets[i].value[v] ) — positive inside plasma, negative inside any solid.
template<typename Mesh>
void ImmersedGeometry<Mesh>::dumpAllLevelSetsScatter(const Mesh& mesh, const std::string& tag)
{
    const int nVerts = mesh.vertIndexNum();
    const int nLS = static_cast<int>(levelSets.size());
    if (nLS == 0 || nVerts == 0)
        return;
    // Build a thin wrapper VertField-like object holding the composite values.
    struct CompositeLS { std::vector<double> value; } composite;
    composite.value.resize(static_cast<std::size_t>(nVerts));
    for (int i = 0; i < nVerts; ++i)
    {
        double v = levelSets[0].value[i];
        for (int j = 1; j < nLS; ++j)
            v = std::min(v, levelSets[j].value[i]);
        composite.value[static_cast<std::size_t>(i)] = v;
    }
    ImmersedGeometryDetail::dumpLevelSetScatter(mesh, composite, 0, tag);
}

namespace ImmersedGeometryDetail
{

template<typename Mesh>
inline bool validVertLabel(const Mesh& mesh, const int label)
{
    if (label < 0 || label >= mesh.vertLabelCapacity())
        return false;
    const int vertIndex = mesh.verts(label).index();
    return vertIndex >= 0 && vertIndex < mesh.vertIndexNum();
}

template<typename Mesh>
inline bool validCellLabel(const Mesh& mesh, const int label)
{
    if (label < 0 || label >= mesh.cellLabelCapacity())
        return false;
    const int cellIndex = mesh.cells(label).index();
    return cellIndex >= 0 && cellIndex < mesh.cellIndexNum();
}

template<typename Mesh>
inline bool validHEdgeLabel(const Mesh& mesh, const int label)
{
    if (label < 0 || label >= mesh.hEdgeLabelCapacity())
        return false;
    const int edgeIndex = mesh.hEdges(label).index();
    return edgeIndex >= 0 && edgeIndex < mesh.hEdgeIndexNum();
}

template<typename Mesh>
inline bool validVEdgeLabel(const Mesh& mesh, const int label)
{
    if (label < 0 || label >= mesh.vEdgeLabelCapacity())
        return false;
    const int edgeIndex = mesh.vEdges(label).index();
    return edgeIndex >= 0 && edgeIndex < mesh.vEdgeIndexNum();
}

template<typename Mesh>
inline bool validCellForInterpolation(const Mesh& mesh, const Cell& cell)
{
    for (int i = 0; i < 4; ++i)
    {
        if (!validVertLabel(mesh, cell.cellToVert[i]))
            return false;
    }

    if (cell.parent != -1)
    {
        if (!validCellLabel(mesh, cell.parent))
            return false;
        const Cell& parentCell = mesh.cells(cell.parent);
        for (int i = 0; i < 4; ++i)
        {
            if (!validCellLabel(mesh, parentCell.son[i]))
                return false;
            const Cell& sonCell = mesh.cells(parentCell.son[i]);
            for (int j = 0; j < 4; ++j)
            {
                if (!validVertLabel(mesh, sonCell.cellToVert[j]))
                    return false;
            }
        }
    }

    return true;
}

template<typename Mesh>
inline int recoverSupportCellLabel(const Mesh& mesh, const Vert& vert)
{
    if (validCellLabel(mesh, vert.supCell))
        return vert.supCell;

    const int preferred[4] = {
        vert.vertToCell[0],
        vert.vertToCell[1],
        vert.vertToCell[2],
        vert.vertToCell[3]
    };
    for (const int label : preferred)
    {
        if (validCellLabel(mesh, label))
            return label;
    }
    return -1;
}

template<typename Mesh>
inline void validateRadialAxisVelocity(const Mesh& mesh, const VertField<double>& vx,
    const VertField<double>& previousVx)
{
    // Once per public entry, never once per vertex/backtrace iteration.
    if (mesh.minCoord.x != 0.0 || !mesh.vertIndexToLabel || !vx.value || !previousVx.value ||
        vx.size() != mesh.vertIndexNum() || previousVx.size() != mesh.vertIndexNum())
        throw std::runtime_error("[ImmersedGeometry] radial-axis velocity requires min.r=0 and matching nodal storage");
    bool foundAxis = false;
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
    {
        const int label = mesh.vertIndexToLabel[k];
        if (!validVertLabel(mesh, label) || mesh.verts(label).index() != k)
            throw std::runtime_error("[ImmersedGeometry] radial-axis velocity has an invalid vertex map");
        if (mesh.verts(label).coord().x != 0.0) continue;
        foundAxis = true;
        if (vx.value[k] != 0.0 || previousVx.value[k] != 0.0)
            throw std::runtime_error("[ImmersedGeometry] actual old/new radial velocity on the symmetry axis must be exactly zero");
    }
    if (!foundAxis) throw std::runtime_error("[ImmersedGeometry] radial-axis mesh has no actual axis vertices");
}
} // namespace ImmersedGeometryDetail

#include "Core/Geometry/ImmersedGeometryBoundary_impl.h"

template<typename Mesh>
ImmersedGeometry<Mesh>::ImmersedGeometry(DictionaryGeometry& dictGeometry, Mesh& mesh) :
points(dictGeometry.points),
lines(dictGeometry.lines)
{
    int levelSetSize = static_cast<int>(dictGeometry.levelSetName.size());
    if (!dictGeometry.signClosureLines.empty()
        && dictGeometry.signClosureLines.size() != static_cast<std::size_t>(levelSetSize))
        throw std::invalid_argument("[ImmersedGeometry] signClosure count must match the level sets");
    signClosureLines_ = dictGeometry.signClosureLines;
    signClosureLines_.resize(levelSetSize);
    levelSetBoundaryConditions_.resize(levelSetSize);
    periodicGeometry_ = !dictGeometry.periodicPairs.empty();
    levelSets.resize(levelSetSize);
    levelSetsPre.resize(levelSetSize);
    levelSetsPre_dx.resize(levelSetSize);
    levelSetsPre_dy.resize(levelSetSize);
    levelSetsPrePre.resize(levelSetSize);
    levelSetsPrePre_dx.resize(levelSetSize);
    levelSetsPrePre_dy.resize(levelSetSize);
    LSs_dx_.resize(levelSetSize);
    LSs_dy_.resize(levelSetSize);
    LSs_dxx_.resize(levelSetSize);
    LSs_dyy_.resize(levelSetSize);
    LSs_dxy_.resize(levelSetSize);
    LSs_dx_backward_.resize(levelSetSize);
    LSs_dy_backward_.resize(levelSetSize);
    LSs_dx_forward_.resize(levelSetSize);
    LSs_dy_forward_.resize(levelSetSize);
    levelSetUpdateDiagnostics_.resize(levelSetSize);
    for (int i = 0; i < levelSetSize; i++)
    {
        levelSets[i].initialize(mesh, 0.0);
        levelSets[i].setInnerPoint(dictGeometry.innerPoints[i]);
        levelSets[i].setLines(dictGeometry.linesOnInterface[i]);
        levelSets[i].setMoving(dictGeometry.moving[i]);

        levelSetsPre[i].initialize(mesh, 0.0);
        levelSetsPre_dx[i].initialize(mesh, 0.0);
        levelSetsPre_dy[i].initialize(mesh, 0.0);
    }

    int phaseSize = static_cast<int>(dictGeometry.phaseName.size());
    phases.resize(phaseSize);
    for (int i = 0; i < phaseSize; i++)
    {
        phases[i].setName(dictGeometry.phaseName[i]);
        phases[i].setOperation(dictGeometry.levelSetOperation[i]);
        phases[i].setLevelSetIndexA(dictGeometry.levelSetIndexA[i]);
        phases[i].setLevelSetIndexB(dictGeometry.levelSetIndexB[i]);
        phases[i].setBoundaryLines(dictGeometry.boundaryLines[i]);
    }
    
    mesh.negFaceVerts.resize(levelSets.size());
    mesh.posFaceVerts.resize(levelSets.size());
}

template<typename Mesh>
ImmersedGeometry<Mesh>::~ImmersedGeometry() {}

template<typename Mesh>
ImmersedGeometryDetail::InitialCurve ImmersedGeometry<Mesh>::makeInitialCurve(SignInt reference, double hMin)
{
    int lineIndex = reference.value;
    // Retain the legacy out-of-range 1-based fallback used by discreteLines.
    if (lineIndex >= static_cast<int>(lines.size()) && lineIndex > 0 &&
        lineIndex - 1 < static_cast<int>(lines.size())) --lineIndex;
    if (lineIndex < 0 || lineIndex >= static_cast<int>(lines.size()))
        throw std::invalid_argument("[ImmersedGeometry] level set references an invalid line index");
    const Line& source = lines[lineIndex];
    ImmersedGeometryDetail::InitialCurve curve;
    if (source.type() == LineType::arc)
        curve.arc.emplace(source.arcCenter(), source.points[0], source.radius(),
                          source.radian(), reference.isNegative);
    else if (source.type() == LineType::straight)
        curve.polyline = reference.isNegative
            ? std::vector<Vector2D>{source.points[1], source.points[0]}
            : std::vector<Vector2D>{source.points[0], source.points[1]};
    else if (source.type() == LineType::spline)
    {
        ExactCurveDistance::requireFinitePoint(source.points[0]);
        ExactCurveDistance::requireFinitePoint(source.points[1]);
        const double samples = std::ceil(ExactCurveDistance::pointDistance(
            source.points[0], source.points[1]) / (0.5 * hMin));
        if (!std::isfinite(samples) || samples >= std::numeric_limits<int>::max())
            throw std::invalid_argument("[ImmersedGeometry] spline sampling count exceeds the supported range");
        curve.polyline = discreteLine(reference, std::max(50, static_cast<int>(samples)));
        if (curve.polyline.size() < 2)
            throw std::invalid_argument("[ImmersedGeometry] spline has insufficient sample points");
    }
    else throw std::invalid_argument("[ImmersedGeometry] unsupported line type during level-set initialization");
    for (const Vector2D& point : curve.polyline) ExactCurveDistance::requireFinitePoint(point);
    return curve;
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::initializeUnsignedDistanceToLines(const Mesh& mesh,
    const std::vector<SignInt>& selectedLines, VertField<double>& distance,
    std::vector<std::array<Vector2D, 2>>* inwardEndpointTangents)
{
    if (selectedLines.empty()) throw std::invalid_argument("[ImmersedGeometry] selected physical line path is empty");
    double hMin = std::numeric_limits<double>::infinity();
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
        hMin = std::min(hMin, mesh.hEdges(mesh.hEdgeIndexToLabel[i]).length());
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
        hMin = std::min(hMin, mesh.vEdges(mesh.vEdgeIndexToLabel[i]).length());
    if (!(hMin > 0.0) || !std::isfinite(hMin))
        throw std::invalid_argument("[ImmersedGeometry] selected distance requires finite positive native spacing");
    std::vector<ImmersedGeometryDetail::InitialCurve> curves;
    std::vector<std::array<Vector2D, 2>> tangents;
    for (const auto reference : selectedLines)
    {
        // This new API has only the documented zero-based SignInt convention.
        if (reference.value < 0 || reference.value >= static_cast<int>(lines.size()))
            throw std::invalid_argument("[ImmersedGeometry] selected line index is outside geometry");
        curves.push_back(makeInitialCurve(reference, hMin));
        if (inwardEndpointTangents)
            tangents.push_back({curves.back().inwardEndpointTangent(false), curves.back().inwardEndpointTangent(true)});
    }
    std::vector<double> staged(mesh.vertIndexNum());
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const int label = mesh.vertIndexToLabel[i];
        if (label < 0 || label >= mesh.vertLabelCapacity() || mesh.verts(label).index() != i)
            throw std::invalid_argument("[ImmersedGeometry] selected distance has invalid vertex mapping");
        const Vector2D point = mesh.verts(label).coord();
        ExactCurveDistance::requireFinitePoint(point);
        double value = std::numeric_limits<double>::infinity();
        for (const auto& curve : curves) value = std::min(value, curve.distance(point));
        if (!std::isfinite(value) || value < 0.0)
            throw std::runtime_error("[ImmersedGeometry] non-finite selected physical distance");
        staged[i] = value;
    }
    distance.initialize(mesh, 0.0);
    for (int i = 0; i < mesh.vertIndexNum(); ++i) distance.value[i] = staged[i];
    if (inwardEndpointTangents) *inwardEndpointTangents = std::move(tangents);
}

template<typename Mesh>
std::function<typename ImmersedGeometry<Mesh>::InitialCurveSample(Vector2D)>
ImmersedGeometry<Mesh>::initialCurveSampler(const Mesh& mesh,
    const std::vector<SignInt>& distanceLines, int signLevelSet)
{
    if (distanceLines.empty() || signLevelSet < 0 || signLevelSet >= static_cast<int>(levelSets.size()))
        throw std::invalid_argument("[ImmersedGeometry] invalid initial sampler line/level-set selection");
    double hMin = std::numeric_limits<double>::infinity();
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
        hMin = std::min(hMin, mesh.hEdges(mesh.hEdgeIndexToLabel[i]).length());
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
        hMin = std::min(hMin, mesh.vEdges(mesh.vEdgeIndexToLabel[i]).length());
    if (!(hMin > 0) || !std::isfinite(hMin))
        throw std::invalid_argument("[ImmersedGeometry] initial sampler requires finite native spacing");
    auto references = levelSets[signLevelSet].lines();
    const auto& closure = signClosureLines_.at(signLevelSet);
    if (references.empty()) throw std::invalid_argument("[ImmersedGeometry] empty initial sign path");
    references.insert(references.end(), closure.begin(), closure.end());
    for (const auto& list : {distanceLines, references})
        for (const auto reference : list)
            if (reference.value < 0 || reference.value >= static_cast<int>(lines.size()))
                throw std::invalid_argument("[ImmersedGeometry] initial sampler requires zero-based line indices");
    if (!closure.empty())
        for (std::size_t i = 0; i < references.size(); ++i)
        {
            const auto a = references[i], b = references[(i + 1) % references.size()];
            const auto end = lines[a.value].points[a.isNegative ? 0 : 1];
            const auto begin = lines[b.value].points[b.isNegative ? 1 : 0];
            if (end.x != begin.x || end.y != begin.y)
                throw std::invalid_argument("[ImmersedGeometry] initial sampler sign closure must join exactly");
        }
    std::vector<ImmersedGeometryDetail::InitialCurve> distanceCurves, signCurves;
    for (auto reference : distanceLines) distanceCurves.push_back(makeInitialCurve(reference, hMin));
    for (auto reference : references) signCurves.push_back(makeInitialCurve(reference, hMin));
    const auto lower = mesh.minCoord, upper = mesh.maxCoord;
    const bool explicitClosure = !closure.empty();
    ExactCurveDistance::requireFinitePoint(levelSets[signLevelSet].innerPoint());
    const bool innerInside = ImmersedGeometryDetail::initialCurveParity(signCurves,
        levelSets[signLevelSet].innerPoint(), lower, upper, explicitClosure);
    return [distanceCurves = std::move(distanceCurves), signCurves = std::move(signCurves),
            lower, upper, explicitClosure, innerInside](Vector2D point) {
        ExactCurveDistance::requireFinitePoint(point);
        double distance = std::numeric_limits<double>::infinity();
        for (const auto& curve : distanceCurves) distance = std::min(distance, curve.distance(point));
        if (!std::isfinite(distance) || distance < 0)
            throw std::runtime_error("[ImmersedGeometry] nonfinite initial curve sample");
        return InitialCurveSample{distance, ImmersedGeometryDetail::initialCurveParity(
            signCurves, point, lower, upper, explicitClosure) == innerInside};
    };
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::initializeLevelSets(Mesh &mesh, double time, bool writeInitialScatter)
{
    std::vector<std::vector<std::optional<double>>> boundaryValues;
    if (hasLevelSetBoundaryConditions())
    {
        validateBoundaryMesh(mesh);
        ImmersedGeometryBoundaryDetail::require(!boundaryHistoryInitialized_,
            "boundary initialization cannot replace existing time history");
        boundaryValues.resize(levelSets.size());
        // Validate every callback before distance/sign initialization writes.
        for (int ls = 0; ls < static_cast<int>(levelSets.size()); ++ls)
        {
            boundaryValues[ls].resize(mesh.vertIndexNum());
            for (int i = 0; i < mesh.vertIndexNum(); ++i)
                boundaryValues[ls][i] = prescribedBoundaryValue(mesh, ls, mesh.verts(mesh.vertIndexToLabel[i]).coord(), time);
        }
    }
    const bool consoleDebug = ConsoleControl::enabled(ConsoleControl::Level::kDebug);
    if (consoleDebug)
        cout << endl << "[ImmersedGeometry] Initializing level set..." << endl;

    // Only spline sampling depends on mesh size. Straight/arc distances and
    // sign classification below use the defining geometric primitives.
    double hMin = std::numeric_limits<double>::max();
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
        hMin = std::min(hMin, mesh.hEdges(mesh.hEdgeIndexToLabel[i]).length());
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
        hMin = std::min(hMin, mesh.vEdges(mesh.vEdgeIndexToLabel[i]).length());
    if (!std::isfinite(hMin) || hMin <= 0.0)
        hMin = 1.0; // fallback; should never happen on a valid mesh

    if (consoleDebug)
        std::cout << "[ImmersedGeometry] hMin=" << hMin << std::endl;

    // Validate once outside the parallel region: exceptions must not escape
    // an OpenMP worker, and an invalid query must not silently become phi=0.
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const int label = mesh.vertIndexToLabel[i];
        if (label < 0 || label >= mesh.vertLabelCapacity())
            throw std::runtime_error("[ImmersedGeometry] invalid vertex label during level-set initialization");
        ExactCurveDistance::requireFinitePoint(mesh.verts(label).coord());
    }

    // Iterate over each level set
    for (int levelSetI = 0; levelSetI < levelSets.size(); levelSetI++)
    {
        const std::size_t physicalCurveCount = levelSets[levelSetI].lines().size();
        auto signReferences = levelSets[levelSetI].lines();
        const auto& closure = signClosureLines_.at(levelSetI);
        if (!closure.empty())
        {
            if (physicalCurveCount == 0)
                throw std::invalid_argument("[ImmersedGeometry] signClosure requires an actual interface path");
            signReferences.insert(signReferences.end(), closure.begin(), closure.end());
            // Check declared endpoints, not trigonometric evaluations of an
            // arc endpoint (which need not reproduce its dictionary bits).
            // The closure must follow the directed true interface, connect its
            // end back to its start, and contain no implicit diagonal gaps.
            for (std::size_t i = 0; i < signReferences.size(); ++i)
            {
                const auto a = signReferences[i], b = signReferences[(i + 1) % signReferences.size()];
                if (a.value < 0 || a.value >= static_cast<int>(lines.size())
                    || b.value < 0 || b.value >= static_cast<int>(lines.size()))
                    throw std::invalid_argument("[ImmersedGeometry] signClosure path has an invalid line index");
                const auto end = lines[a.value].points[a.isNegative ? 0 : 1];
                const auto begin = lines[b.value].points[b.isNegative ? 1 : 0];
                ExactCurveDistance::requireFinitePoint(end);
                ExactCurveDistance::requireFinitePoint(begin);
                if (end.x != begin.x || end.y != begin.y)
                    throw std::invalid_argument("[ImmersedGeometry] directed signClosure and interface endpoints must join exactly");
            }
        }
        std::vector<ImmersedGeometryDetail::InitialCurve> curves;
        curves.reserve(signReferences.size());
        ExactCurveDistance::requireFinitePoint(levelSets[levelSetI].innerPoint());
        for (int localLineI = 0; localLineI < static_cast<int>(signReferences.size()); localLineI++)
        {
            curves.push_back(makeInitialCurve(signReferences[localLineI], hMin));
        }
        if (curves.empty())
            throw std::invalid_argument("[ImmersedGeometry] level-set initialization requires at least one curve");
        if (consoleDebug)
            std::cout << "[ImmersedGeometry] computing distances..." << std::endl;
        // Hard minimum is the distance to the union of the supplied curves.
        // Smoothing that minimum would move corners and perturb the zero set.
#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < mesh.vertIndexNum(); i++)
        {
            const Vector2D coordinate = mesh.verts(mesh.vertIndexToLabel[i]).coord();
            double result = std::numeric_limits<double>::infinity();
            for (std::size_t k = 0; k < physicalCurveCount; ++k)
                result = std::min(result, curves[k].distance(coordinate));
            levelSets[levelSetI].value[i] = result;
        }
        for (int i = 0; i < mesh.vertIndexNum(); ++i)
            if (!std::isfinite(levelSets[levelSetI].value[i]))
                throw std::runtime_error("[ImmersedGeometry] non-finite curve distance during level-set initialization");
        if (consoleDebug)
            std::cout << "[ImmersedGeometry] distances done." << std::endl;

        if (consoleDebug)
            std::cout << "[ImmersedGeometry] assigning sign..." << std::endl;
        // Use the very same arcs for parity; a sampled polygon would leave a
        // thin band with correct distance magnitude but the wrong sign.
        auto pointInsideCurves = [&](const Vector2D& p) -> bool {
            // Same physical primitives, flattened connectors and sign-only
            // domain-wall trace are also used by the owned point sampler.
            return ImmersedGeometryDetail::initialCurveParity(curves, p,
                mesh.minCoord, mesh.maxCoord, !closure.empty());
        };

        const bool innerInside = pointInsideCurves(levelSets[levelSetI].innerPoint());
#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < mesh.vertIndexNum(); i++)
        {
            const Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);
            const Vector2D coord = vert.coord();
            bool inside = pointInsideCurves(coord);
            if (!innerInside)
                inside = !inside;
            if (inside)
                levelSets[levelSetI].value[i] = -fabs(levelSets[levelSetI].value[i]);
            else
                levelSets[levelSetI].value[i] = fabs(levelSets[levelSetI].value[i]);
        }
        if (consoleDebug)
            std::cout << "[ImmersedGeometry] sign done." << std::endl;

        // Legacy zero displacement is retained only for unconstrained nodes.
#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < mesh.vertIndexNum(); i++)
        {
            if (!boundaryValues.empty() && boundaryValues[levelSetI][i])
            {
                levelSets[levelSetI].value[i] = *boundaryValues[levelSetI][i];
                continue;
            }
            if (fabs(levelSets[levelSetI].value[i]) < 1e-10)
            {
                levelSets[levelSetI].value[i] = 1e-10 * (levelSets[levelSetI].value[i] <= 0 ? -1.0 : 1.0);
            }
        }
        if (consoleDebug)
            std::cout << "[ImmersedGeometry] zero check done." << std::endl;

        checkLevelSet(mesh, levelSetI);
        if (consoleDebug)
            std::cout << "[ImmersedGeometry] levelSetI done." << std::endl;
    }

    // --- LS diagnostic: composite plasma region (min of all LS) ---
    if (writeInitialScatter) dumpAllLevelSetsScatter(mesh, "composite_plasma");
}

template<typename Mesh>
vector<Vector2D> ImmersedGeometry<Mesh>::discreteLines(const int levelSetI, const int localLineI, int numSamples)
{
    if (numSamples < 1 || numSamples == std::numeric_limits<int>::max())
        throw std::invalid_argument("[ImmersedGeometry] line sampling count must be positive and below INT_MAX");
    return discreteLine(levelSets[levelSetI].lines()[localLineI], numSamples);
}

template<typename Mesh>
vector<Vector2D> ImmersedGeometry<Mesh>::discreteLine(SignInt line, int numSamples)
{
    vector<Vector2D> discretePoints;
    if (numSamples < 1 || numSamples == std::numeric_limits<int>::max())
        throw std::invalid_argument("[ImmersedGeometry] line sampling count must be positive and below INT_MAX");

    int lineIdx = line.value;
    // Accept either 0-based or 1-based indexing from geometry.dict.
    if (lineIdx < 0 || lineIdx >= static_cast<int>(lines.size()))
    {
        if (lineIdx > 0 && (lineIdx - 1) >= 0 && (lineIdx - 1) < static_cast<int>(lines.size()))
            lineIdx -= 1;
    }
    if (lineIdx < 0 || lineIdx >= static_cast<int>(lines.size()))
    {
        std::cerr << "[ImmersedGeometry] WARNING: line index out of range: " << line.value
                  << " (lines=" << lines.size() << ").\n";
        return discretePoints;
    }
    int numPoints = numSamples;  // Number of points to discretize the line

    // Discretize different types of line elements into discrete points
    switch (lines[lineIdx].type())
    {
        case LineType::straight:
        {
            // Discretize straight line
            for (int i = 0; i <= numPoints; ++i)
            {
                double t = static_cast<double>(i) / numPoints;
                Vector2D p = line.isNegative 
                    ? lines[lineIdx].points[1] + (lines[lineIdx].points[0] - lines[lineIdx].points[1]) * t
                    : lines[lineIdx].points[0] + (lines[lineIdx].points[1] - lines[lineIdx].points[0]) * t ;
                discretePoints.push_back(p);
            }
        }
        break;

        case LineType::arc:
        {
            // Share the exact arc's parameterization. In particular reversing
            // a line reverses t, not an end angle with mismatched pi units.
            const auto& source = lines[lineIdx];
            const ExactCurveDistance::CircularArc arc(source.arcCenter(), source.points[0],
                source.radius(), source.radian(), line.isNegative);
            for (int i = 0; i <= numPoints; ++i)
            {
                const double t = static_cast<double>(i) / numPoints;
                discretePoints.push_back(arc.point(t));
            }
        }
        break;

        case LineType::spline:
        {
            // Discretize spline
            BSplineCurve bspline(lines[lineIdx].featurePoints(), 3);
            for (int i = 0; i <= numPoints; ++i)
            {
                double t = static_cast<double>(i) / numPoints;
                Vector2D p = line.isNegative ? bspline.evaluate(1 - t) : bspline.evaluate(t);
                discretePoints.push_back(p);
            }
        }
        break;

        default:
            cerr << "[ImmersedGeometry] ERROR: unknown Line type." << endl;
            exit(1);
    }

    return discretePoints;
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::markPhases(Mesh& mesh)
{
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        cout << "[ImmersedGeometry] Marking phases..." << endl;

    if (phases.empty()) return;
    // Start from unassigned tags so a changed interface cannot leave stale
    // membership. Stage all tags: malformed expressions publish no partial map.
    std::vector<int> nextPhase(mesh.vertIndexNum(), -1);
    int phaseIndex = 0;  // Phase index
    const int levelSetSize = static_cast<int>(levelSets.size());
    for (Phase &phase : phases)
    {
        const auto operations = phase.operation();
        const auto indicesA = phase.levelSetIndexA(), indicesB = phase.levelSetIndexB();
        if (operations.empty() || indicesA.size() != operations.size() || indicesB.size() != operations.size())
            throw std::invalid_argument("[ImmersedGeometry] phase requires matching nonempty operations and references");
        vector<LevelSet> appLevelSets(operations.size());
        for (std::size_t i = 0; i < operations.size(); i++)
        {
            const int levelSetIdxA = indicesA[i].value;
            const int levelSetIdxB = indicesB[i].value;
            if (levelSetIdxA < 0 || levelSetIdxB < 0
                || levelSetIdxA >= levelSetSize + static_cast<int>(i)
                || levelSetIdxB >= levelSetSize + static_cast<int>(i))
                throw std::invalid_argument("[ImmersedGeometry] phase references an invalid or uncomputed level set");
            if (operations[i] != LevelSetOperation::non && operations[i] != LevelSetOperation::min
                && operations[i] != LevelSetOperation::max)
                throw std::invalid_argument("[ImmersedGeometry] unknown LevelSetOperation");
            if (operations[i] == LevelSetOperation::non
                && (levelSetIdxA != levelSetIdxB || indicesA[i].isNegative != indicesB[i].isNegative))
                throw std::invalid_argument("[ImmersedGeometry] invalid non phase operation");
            for (const int index : {levelSetIdxA, levelSetIdxB})
                if (index < levelSetSize && (!levelSets[index].value
                    || levelSets[index].size() != mesh.vertIndexNum()))
                    throw std::invalid_argument("[ImmersedGeometry] phase level-set storage must match the mesh");
            appLevelSets[i].initialize(mesh, 0.0);
            for (int v = 0; v < mesh.vertIndexNum(); v++)
            {
                double valueA;
                if (levelSetIdxA >= levelSetSize)
                {
                    valueA = (indicesA[i].isNegative ? -1.0 : 1.0) * appLevelSets[levelSetIdxA - levelSetSize].value[v];
                }
                else
                {
                    valueA = (indicesA[i].isNegative ? -1.0 : 1.0) * levelSets[levelSetIdxA].value[v];
                }

                double valueB;
                if (levelSetIdxB >= levelSetSize)
                {
                    valueB = (indicesB[i].isNegative ? -1.0 : 1.0) * appLevelSets[levelSetIdxB - levelSetSize].value[v];
                }
                else
                {
                    valueB = (indicesB[i].isNegative ? -1.0 : 1.0) * levelSets[levelSetIdxB].value[v];
                }
                
                if (!std::isfinite(valueA) || !std::isfinite(valueB))
                    throw std::invalid_argument("[ImmersedGeometry] phase level-set value must be finite");
                if (operations[i] == LevelSetOperation::non)
                {
                    appLevelSets[i].value[v] = valueA;
                }
                else if (operations[i] == LevelSetOperation::min)
                {
                    appLevelSets[i].value[v] = std::min(valueA, valueB);
                }
                else if (operations[i] == LevelSetOperation::max)
                {
                    appLevelSets[i].value[v] = std::max(valueA, valueB);
                }
            }
        }
        // Membership belongs only to the final composite expression. An
        // intermediate constraint must not overwrite another complete phase.
        for (int v = 0; v < mesh.vertIndexNum(); ++v)
            if (appLevelSets.back().value[v] < 0)
                nextPhase[v] = phaseIndex;
        phaseIndex++;
    }
    for (int v = 0; v < mesh.vertIndexNum(); ++v)
        mesh.verts(mesh.vertIndexToLabel[v]).setPhase(nextPhase[v]);
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::diffForReinitLS(const Mesh& mesh, const int levelSetI)
{
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);
        
        LSs_dxx_[levelSetI].value[i] = FiniteDifference::dxx(levelSets[levelSetI], mesh, vert);
        LSs_dyy_[levelSetI].value[i] = FiniteDifference::dyy(levelSets[levelSetI], mesh, vert);
    }

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);
        
        LSs_dx_backward_[levelSetI].value[i] = FiniteDifferenceForLevelSet::dx_backward(levelSets[levelSetI], mesh, vert, LSs_dxx_[levelSetI]);
        LSs_dy_backward_[levelSetI].value[i] = FiniteDifferenceForLevelSet::dy_backward(levelSets[levelSetI], mesh, vert, LSs_dyy_[levelSetI]);
        LSs_dx_forward_[levelSetI].value[i] = FiniteDifferenceForLevelSet::dx_forward(levelSets[levelSetI], mesh, vert, LSs_dxx_[levelSetI]);
        LSs_dy_forward_[levelSetI].value[i] = FiniteDifferenceForLevelSet::dy_forward(levelSets[levelSetI], mesh, vert, LSs_dyy_[levelSetI]);
    }
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::searchCutEdges(const Mesh& mesh, const int levelSetI)
{
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
    {
        Edge& hEdge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        Vert& vert0 = mesh.verts(hEdge.edgeToVert[0]);
        Vert& vert1 = mesh.verts(hEdge.edgeToVert[1]);

        hEdge.setCut(false);
        hEdge.setCutLength(0.0);

        double value0 = levelSets[levelSetI].value[vert0.index()];
        double value1 = levelSets[levelSetI].value[vert1.index()];

        if (value0 * value1 <= 0.0)
        {
            hEdge.setCut(true);
        }
    }

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
    {
        Edge& vEdge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        Vert& vert0 = mesh.verts(vEdge.edgeToVert[0]);
        Vert& vert1 = mesh.verts(vEdge.edgeToVert[1]);

        vEdge.setCut(false);
        vEdge.setCutLength(0.0);

        double value0 = levelSets[levelSetI].value[vert0.index()];
        double value1 = levelSets[levelSetI].value[vert1.index()];

        if (value0 * value1 <= 0.0)
        {
            vEdge.setCut(true);
        }
    }
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::computeCutLength(const Mesh& mesh, const int levelSetI)
{
    searchCutEdges(mesh, levelSetI);

    // Windows Release crashes in this routine have been traced to the OpenMP
    // path around computeCutLength$omp$1. Keep the cut-length rebuild serial
    // until the underlying mesh/edge concurrency issue is fully resolved.
    for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
    {
        Edge& hEdge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        if (hEdge.cut())
        {
            Vert& vert0 = mesh.verts(hEdge.edgeToVert[0]);
            Vert& vert1 = mesh.verts(hEdge.edgeToVert[1]);

            double value0 = levelSets[levelSetI].value[vert0.index()];
            double value1 = levelSets[levelSetI].value[vert1.index()];
            double dxx0 = LSs_dxx_[levelSetI].value[vert0.index()];
            double dxx1 = LSs_dxx_[levelSetI].value[vert1.index()];

            double c2 = 0.5 * minmod(dxx0, dxx1);
            double c1 = (value1 - value0) / hEdge.length();
            double c0 = 0.5 * (value0 + value1) - 0.25 * c2 * hEdge.length() * hEdge.length();

            if (fabs(c2) < MINVAL)
            {
                hEdge.setCutLength(0.5 * hEdge.length() - c0 / c1);
            }
            else
            {
                if (value0 < 0.0)
                {
                    hEdge.setCutLength(0.5 * (hEdge.length() + (-c1 + sqrt(c1 * c1 - 4.0 * c2 *c0)) / c2));
                }
                else
                {
                    hEdge.setCutLength(0.5 * (hEdge.length() + (-c1 - sqrt(c1 * c1 - 4.0 * c2 *c0)) / c2));
                }
            }

            if (hEdge.cutLength() < MINVAL)
            {
                hEdge.setCutLength(MINVAL);
            }
            else if (hEdge.cutLength() > hEdge.length() - MINVAL)
            {
                hEdge.setCutLength(hEdge.length() - MINVAL);
            }
        }
    }

    for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
    {
        Edge& vEdge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        if (vEdge.cut())
        {
            Vert& vert0 = mesh.verts(vEdge.edgeToVert[0]);
            Vert& vert1 = mesh.verts(vEdge.edgeToVert[1]);

            double value0 = levelSets[levelSetI].value[vert0.index()];
            double value1 = levelSets[levelSetI].value[vert1.index()];
            double dyy0 = LSs_dyy_[levelSetI].value[vert0.index()];
            double dyy1 = LSs_dyy_[levelSetI].value[vert1.index()];

            double c2 = 0.5 * minmod(dyy0, dyy1);
            double c1 = (value1 - value0) / vEdge.length();
            double c0 = 0.5 * (value0 + value1) - 0.25 * c2 * vEdge.length() * vEdge.length();

            if (fabs(c2) < MINVAL)
            {
                vEdge.setCutLength(0.5 * vEdge.length() - c0 / c1);
            }
            else
            {
                if (value0 < 0.0)
                {
                    vEdge.setCutLength(0.5 * (vEdge.length() + (-c1 + sqrt(c1 * c1 - 4.0 * c2 *c0)) / c2));
                }
                else
                {
                    vEdge.setCutLength(0.5 * (vEdge.length() + (-c1 - sqrt(c1 * c1 - 4.0 * c2 *c0)) / c2));
                }
            }

            if (vEdge.cutLength() < MINVAL)
            {
                vEdge.setCutLength(MINVAL);
            }
            else if (vEdge.cutLength() > vEdge.length() - MINVAL)
            {
                vEdge.setCutLength(vEdge.length() - MINVAL);
            }
        }
    }
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::reinitializeLevelSets(const Mesh& mesh, const double& t)
{
    pdeReinitializeMovingLevelSets(mesh, t);
    refreshLevelSetDerivedData(mesh, t, /*reinitializeMovingLevelSets=*/false);
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::pdeReinitializeMovingLevelSets(const Mesh& mesh, const double& t)
{
    const std::vector<int> activeLevelSetIndices;
    pdeReinitializeMovingLevelSets(mesh, t, activeLevelSetIndices);
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::pdeReinitializeMovingLevelSets(
    const Mesh& mesh,
    const double& t,
    const std::vector<int>& activeLevelSetIndices)
{
    if (hasLevelSetBoundaryConditions())
        throw std::runtime_error("[ImmersedGeometry] scalar boundary constraints do not support legacy PDE reinitialization");
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[ImmersedGeometry] PDE reinitializing moving level sets...\n";

    for (auto& lsDiag : levelSetUpdateDiagnostics_)
    {
        lsDiag.reinitialized = false;
        lsDiag.reinitIterations = 0;
        lsDiag.reinitFineResidual = 0.0;
        lsDiag.reinitCoarseResidual = 0.0;
    }

    std::vector<unsigned char> shouldReinit;
    if (!activeLevelSetIndices.empty())
    {
        shouldReinit.assign(levelSets.size(), static_cast<unsigned char>(0));
        for (const int levelSetI : activeLevelSetIndices)
        {
            if (levelSetI >= 0 && levelSetI < static_cast<int>(levelSets.size()))
                shouldReinit[static_cast<std::size_t>(levelSetI)] = 1u;
        }
    }

    for (int levelSetI = 0; levelSetI < levelSets.size(); levelSetI++)
    {
        if (!levelSets[levelSetI].moving())
            continue;
        if (!shouldReinit.empty() && shouldReinit[static_cast<std::size_t>(levelSetI)] == 0u)
            continue;

        prepareLevelSetDerivedStorage(mesh, levelSetI);
        computeSecondDerivatives(mesh, levelSetI);
        computeCutLength(mesh, levelSetI);
        reinitLevelSet(mesh, levelSetI, t);
    }
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::resetUpdateDiagnostics()
{
    traceDiagnostics_ = TraceDiagnostics{};
    levelSetUpdateDiagnostics_.assign(levelSets.size(), LevelSetUpdateDiagnostics{});
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::syncHistoricalLevelSetStorageToCurrent(const Mesh& mesh, const int levelSetI)
{
    const bool needsResize =
        (levelSetsPre[levelSetI].size() != mesh.vertIndexNum()) ||
        (levelSetsPre_dx[levelSetI].size() != mesh.vertIndexNum()) ||
        (levelSetsPre_dy[levelSetI].size() != mesh.vertIndexNum()) ||
        (levelSetsPrePre[levelSetI].size() != mesh.vertIndexNum()) ||
        (levelSetsPrePre_dx[levelSetI].size() != mesh.vertIndexNum()) ||
        (levelSetsPrePre_dy[levelSetI].size() != mesh.vertIndexNum());

    if (!needsResize)
        return;

    levelSetsPre[levelSetI].initialize(mesh, 0.0);
    levelSetsPre_dx[levelSetI].initialize(mesh, 0.0);
    levelSetsPre_dy[levelSetI].initialize(mesh, 0.0);
    levelSetsPrePre[levelSetI].initialize(mesh, 0.0);
    levelSetsPrePre_dx[levelSetI].initialize(mesh, 0.0);
    levelSetsPrePre_dy[levelSetI].initialize(mesh, 0.0);

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        levelSetsPre[levelSetI].value[i] = levelSets[levelSetI].value[i];
        levelSetsPre_dx[levelSetI].value[i] = LSs_dx_[levelSetI].value[i];
        levelSetsPre_dy[levelSetI].value[i] = LSs_dy_[levelSetI].value[i];
        levelSetsPrePre[levelSetI].value[i] = levelSets[levelSetI].value[i];
        levelSetsPrePre_dx[levelSetI].value[i] = LSs_dx_[levelSetI].value[i];
        levelSetsPrePre_dy[levelSetI].value[i] = LSs_dy_[levelSetI].value[i];
    }

}

template<typename Mesh>
void ImmersedGeometry<Mesh>::summarizeLevelSetDiagnostics(const Mesh& mesh, const int levelSetI)
{
    LevelSetUpdateDiagnostics lsDiag = levelSetUpdateDiagnostics_[levelSetI];
    lsDiag.minValue = std::numeric_limits<double>::infinity();
    lsDiag.maxValue = -std::numeric_limits<double>::infinity();
    lsDiag.nearZeroCount = 0;
    lsDiag.nonFiniteCount = 0;
    lsDiag.narrowBandCount = 0;

    double minEdgeLength = std::numeric_limits<double>::max();
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);
        if (vert.vertToHEdge[0] != -1)
            minEdgeLength = std::min(minEdgeLength, mesh.hEdges(vert.vertToHEdge[0]).length());
        if (vert.vertToHEdge[1] != -1)
            minEdgeLength = std::min(minEdgeLength, mesh.hEdges(vert.vertToHEdge[1]).length());
        if (vert.vertToVEdge[0] != -1)
            minEdgeLength = std::min(minEdgeLength, mesh.vEdges(vert.vertToVEdge[0]).length());
        if (vert.vertToVEdge[1] != -1)
            minEdgeLength = std::min(minEdgeLength, mesh.vEdges(vert.vertToVEdge[1]).length());
    }
    if (!std::isfinite(minEdgeLength))
        minEdgeLength = 0.0;
    const double narrowBandWidth = levelSets[levelSetI].narrowbandScale() * minEdgeLength;

    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        const double value = levelSets[levelSetI].value[i];
        if (!std::isfinite(value))
        {
            lsDiag.nonFiniteCount += 1;
            continue;
        }
        lsDiag.minValue = std::min(lsDiag.minValue, value);
        lsDiag.maxValue = std::max(lsDiag.maxValue, value);
        if (std::abs(value) <= 10.0 * MINVAL)
            lsDiag.nearZeroCount += 1;
        if (std::abs(value) < narrowBandWidth)
            lsDiag.narrowBandCount += 1;
    }
    if (!std::isfinite(lsDiag.minValue))
        lsDiag.minValue = 0.0;
    if (!std::isfinite(lsDiag.maxValue))
        lsDiag.maxValue = 0.0;
    levelSetUpdateDiagnostics_[levelSetI] = lsDiag;
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::prepareLevelSetDerivedStorage(const Mesh& mesh, const int levelSetI)
{
    const int vertNum = mesh.vertIndexNum();
    LSs_dx_[levelSetI].initialize(mesh, 0.0);
    LSs_dy_[levelSetI].initialize(mesh, 0.0);
    LSs_dxx_[levelSetI].initialize(vertNum, 0.0);
    LSs_dyy_[levelSetI].initialize(vertNum, 0.0);
    LSs_dxy_[levelSetI].initialize(vertNum, 0.0);
    LSs_dx_backward_[levelSetI].initialize(vertNum, 0.0);
    LSs_dy_backward_[levelSetI].initialize(vertNum, 0.0);
    LSs_dx_forward_[levelSetI].initialize(vertNum, 0.0);
    LSs_dy_forward_[levelSetI].initialize(vertNum, 0.0);
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::computeSecondDerivatives(const Mesh& mesh, const int levelSetI)
{
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);
        LSs_dxx_[levelSetI].value[i] = FiniteDifference::dxx(levelSets[levelSetI], mesh, vert);
        LSs_dyy_[levelSetI].value[i] = FiniteDifference::dyy(levelSets[levelSetI], mesh, vert);
    }
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::computeFirstDerivativesAndMixed(const Mesh& mesh, const int levelSetI)
{
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);
        LSs_dx_[levelSetI].value[i] = FiniteDifference::dx(levelSets[levelSetI], mesh, vert);
        LSs_dy_[levelSetI].value[i] = FiniteDifference::dy(levelSets[levelSetI], mesh, vert);
    }

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);
        LSs_dxy_[levelSetI].value[i] = FiniteDifference::dxy(LSs_dx_[levelSetI], mesh, vert);
    }
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::refreshLevelSetDerivedData(
    const Mesh& mesh,
    const double& t,
    const bool reinitializeMovingLevelSets)
{
    if (hasLevelSetBoundaryConditions())
        throw std::runtime_error("[ImmersedGeometry] legacy cut lengths do not support declared scalar boundaries; use refreshLevelSetDerivatives and shared Hermite geometry");
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
    {
        std::cout << "[ImmersedGeometry] Refreshing level-set derived data"
                  << " (reinitializeMoving=" << (reinitializeMovingLevelSets ? "true" : "false") << ")...\n";
    }

    for (int levelSetI = 0; levelSetI < levelSets.size(); levelSetI++)
    {
        prepareLevelSetDerivedStorage(mesh, levelSetI);
        computeSecondDerivatives(mesh, levelSetI);

        // Exact rigid-shift paths can rebuild a fresh signed field directly
        // from moved geometry and only need derivatives/cut-info refreshed.
        if (reinitializeMovingLevelSets && levelSets[levelSetI].moving())
        {
            computeCutLength(mesh, levelSetI);

            reinitLevelSet(mesh, levelSetI, t);

        }
        else
        {
            computeCutLength(mesh, levelSetI);
        }

        computeFirstDerivativesAndMixed(mesh, levelSetI);
        syncHistoricalLevelSetStorageToCurrent(mesh, levelSetI);
        summarizeLevelSetDiagnostics(mesh, levelSetI);
    }

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[ImmersedGeometry] Level-set derived data refresh completed.\n";
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::refreshLevelSetScalarDiagnostics(const Mesh& mesh)
{
    if (levelSetUpdateDiagnostics_.size() != levelSets.size())
        levelSetUpdateDiagnostics_.assign(levelSets.size(), LevelSetUpdateDiagnostics{});
    for (int levelSetI = 0; levelSetI < static_cast<int>(levelSets.size()); ++levelSetI)
        summarizeLevelSetDiagnostics(mesh, levelSetI);
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::reinitLevelSet(const Mesh& mesh, const int levelSetI, const double& t)
{
    VertField<bool> needReinit(mesh, true);
    VertField<double> preL(mesh, 0.0);

    diffForReinitLS(mesh, levelSetI);

    double fineResidual = 0.0;
    double coarseResidual = 0.0;
    double tStep = 0.5;
    double fictitiousTime = 0.0;

    int iter = 0;
    do
    {
        double minEdgeLength = std::numeric_limits<double>::max();
        double minCutLength = minEdgeLength;

#pragma omp parallel for reduction(min:minCutLength, minEdgeLength) num_threads(thread_num)
        for (int i = 0; i < mesh.vertIndexNum(); i++)
        {
            if (needReinit.value[i])
            {
                Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);

                if (vert.vertToHEdge[0] != -1)
                {
                    Edge& hEdge = mesh.hEdges(vert.vertToHEdge[0]);
                    minEdgeLength = std::min(minEdgeLength, hEdge.length());
                    if (hEdge.cut())
                    {
                        minCutLength = std::min(minCutLength, hEdge.length() - hEdge.cutLength());
                    }
                    else
                    {
                        minCutLength = std::min(minCutLength, hEdge.length());
                    }
                }
                if (vert.vertToHEdge[1] != -1)
                {
                    Edge& hEdge = mesh.hEdges(vert.vertToHEdge[1]);
                    minEdgeLength = std::min(minEdgeLength, hEdge.length());
                    if (hEdge.cut())
                    {
                        minCutLength = std::min(minCutLength, hEdge.cutLength());
                    }
                    else
                    {
                        minCutLength = std::min(minCutLength, hEdge.length());
                    }
                }

                if (vert.vertToVEdge[0] != -1)
                {
                    Edge& vEdge = mesh.vEdges(vert.vertToVEdge[0]);
                    minEdgeLength = std::min(minEdgeLength, vEdge.length());
                    if (vEdge.cut())
                    {
                        minCutLength = std::min(minCutLength, vEdge.length() - vEdge.cutLength());
                    }
                    else
                    {
                        minCutLength = std::min(minCutLength, vEdge.length());
                    }
                }
                if (vert.vertToVEdge[1] != -1)
                {
                    Edge& vEdge = mesh.vEdges(vert.vertToVEdge[1]);
                    minEdgeLength = std::min(minEdgeLength, vEdge.length());
                    if (vEdge.cut())
                    {
                        minCutLength = std::min(minCutLength, vEdge.cutLength());
                    }
                    else
                    {
                        minCutLength = std::min(minCutLength, vEdge.length());
                    }
                }
            }
        }
    
        fictitiousTime = tStep * minCutLength;

#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < mesh.vertIndexNum(); i++)
        {
            if (needReinit.value[i])
            {
                preL.value[i] = levelSets[levelSetI].value[i];

                levelSets[levelSetI].value[i] -= fictitiousTime * static_cast<double>(sgn(preL.value[i])) * (operatorHG(LSs_dx_forward_[levelSetI].value[i], LSs_dx_backward_[levelSetI].value[i], LSs_dy_forward_[levelSetI].value[i], LSs_dy_backward_[levelSetI].value[i], preL.value[i]) - 1.0);
            }
        }

        diffForReinitLS(mesh, levelSetI);

#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < mesh.vertIndexNum(); i++)
        {
            if (needReinit.value[i])
            {
                levelSets[levelSetI].value[i] -= fictitiousTime * static_cast<double>(sgn(preL.value[i])) * (operatorHG(LSs_dx_forward_[levelSetI].value[i], LSs_dx_backward_[levelSetI].value[i], LSs_dy_forward_[levelSetI].value[i], LSs_dy_backward_[levelSetI].value[i], preL.value[i]) - 1.0);

                levelSets[levelSetI].value[i] = 0.5 * (preL.value[i] + levelSets[levelSetI].value[i]);
            }
        }

        diffForReinitLS(mesh, levelSetI);

        fineResidual = 0.0;
        coarseResidual = 0.0;
#pragma omp parallel for reduction(max:fineResidual, coarseResidual) num_threads(thread_num)
        for (int i = 0; i < mesh.vertIndexNum(); i++)
        {
            double _residual = fabs(operatorHG(LSs_dx_forward_[levelSetI].value[i], LSs_dx_backward_[levelSetI].value[i], LSs_dy_forward_[levelSetI].value[i], LSs_dy_backward_[levelSetI].value[i], levelSets[levelSetI].value[i]) - 1.0);

            if (fabs(levelSets[levelSetI].value[i]) < levelSets[levelSetI].narrowbandScale() * minEdgeLength)
            {
                if (_residual > fineResidual)
                {
                    fineResidual = _residual;
                }
            }
            else
            {
                if (_residual > coarseResidual)
                {
                    coarseResidual = _residual;
                }
            }

            if (_residual < levelSets[levelSetI].reinitFineEps())
            {
                needReinit.value[i] = false;
            }
            else
            {
                needReinit.value[i] = true;
            }
        }

        if (iter % 100 == 0 && ConsoleControl::enabled(ConsoleControl::Level::kTrace))
        {
            std::cout << "[ImmersedGeometry] Iter=" << iter << ", narrow-band residual=" << fineResidual
                      << ", outside residual=" << coarseResidual << ", fictitious dt=" << fictitiousTime << "\n";
        }

        iter++;

    } while ((fineResidual > levelSets[levelSetI].reinitFineEps() || coarseResidual > levelSets[levelSetI].reinitCoarseEps()) && iter < 500);
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
    {
        std::cout << "[ImmersedGeometry] Iter=" << iter << ", narrow-band residual=" << fineResidual
                  << ", outside residual=" << coarseResidual << "\n";
    }

    checkLevelSet(mesh, levelSetI);
    LevelSetUpdateDiagnostics& lsDiag = levelSetUpdateDiagnostics_[levelSetI];
    lsDiag.reinitialized = true;
    lsDiag.reinitIterations = iter;
    lsDiag.reinitFineResidual = fineResidual;
    lsDiag.reinitCoarseResidual = coarseResidual;
    summarizeLevelSetDiagnostics(mesh, levelSetI);
}

// ---------------------------------------------------------------------------
// exactRedistanceLevelSet — recompute exact signed distance from the current
// zero contour.  The zero contour is extracted as polyline segments from cut
// edges (linear interpolation of the current level-set field).  Then every
// vertex gets its exact unsigned distance to that polyline, with the sign
// taken from the existing advected field.  This eliminates the far-field
// roughness that accumulates from repeated advection + PDE reinit cycles.
// ---------------------------------------------------------------------------
template<typename Mesh>
void ImmersedGeometry<Mesh>::exactRedistanceLevelSet(const Mesh& mesh, const int levelSetI)
{
    // Step 1: Extract current zero-contour segments from cut edges.
    std::vector<Vector2D> zeroContourSegs; // pairs of points (p0,p1,p0,p1,...)

    auto extractFromEdge = [&](const Edge& edge) {
        const int v0 = edge.edgeToVert[0];
        const int v1 = edge.edgeToVert[1];
        if (v0 < 0 || v1 < 0) return;
        const double psi0 = levelSets[levelSetI].value[mesh.verts(v0).index()];
        const double psi1 = levelSets[levelSetI].value[mesh.verts(v1).index()];
        if (psi0 * psi1 > 0.0) return; // no sign change
        const double denom = psi1 - psi0;
        if (std::abs(denom) < 1e-30) return;
        const double t = -psi0 / denom;
        const double tc = (t < 0.0) ? 0.0 : ((t > 1.0) ? 1.0 : t);
        const Vector2D p0 = mesh.verts(v0).coord();
        const Vector2D p1 = mesh.verts(v1).coord();
        zeroContourSegs.push_back(p0 + (p1 - p0) * tc);
    };

    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
        extractFromEdge(mesh.hEdges(mesh.hEdgeIndexToLabel[i]));
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
        extractFromEdge(mesh.vEdges(mesh.vEdgeIndexToLabel[i]));

    if (zeroContourSegs.empty())
    {
        if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
            std::cout << "[ImmersedGeometry] exactRedistance: no zero-contour points found for LS "
                      << levelSetI << "; skipping.\n";
        return;
    }

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[ImmersedGeometry] exactRedistance: LS " << levelSetI
                  << " zero-contour points=" << zeroContourSegs.size() << "\n";

    // Step 2: Build polyline segments from cell-by-cell zero crossings.
    // Each cell with 2 crossing points contributes one segment.
    std::vector<std::pair<Vector2D, Vector2D>> segments;

    // Collect crossing points per cell
    // For each cell, gather zero-crossing points on its 4 edges.
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& cell = mesh.cells(cLbl);
        std::vector<Vector2D> crossings;
        crossings.reserve(4);

        auto checkEdge = [&](int vLabel0, int vLabel1) {
            if (vLabel0 < 0 || vLabel1 < 0) return;
            const int vi0 = mesh.verts(vLabel0).index();
            const int vi1 = mesh.verts(vLabel1).index();
            if (vi0 < 0 || vi1 < 0) return;
            const double psi0 = levelSets[levelSetI].value[vi0];
            const double psi1 = levelSets[levelSetI].value[vi1];
            if (psi0 * psi1 > 0.0) return;
            const double denom = psi1 - psi0;
            if (std::abs(denom) < 1e-30) return;
            double t = -psi0 / denom;
            t = (t < 0.0) ? 0.0 : ((t > 1.0) ? 1.0 : t);
            const Vector2D a = mesh.verts(vLabel0).coord();
            const Vector2D b = mesh.verts(vLabel1).coord();
            crossings.push_back(a + (b - a) * t);
        };

        // Cell vertices are SW(0), SE(1), NE(2), NW(3)
        // Edges: S(0-1), E(1-2), N(2-3), W(3-0)
        checkEdge(cell.cellToVert[0], cell.cellToVert[1]);
        checkEdge(cell.cellToVert[1], cell.cellToVert[2]);
        checkEdge(cell.cellToVert[2], cell.cellToVert[3]);
        checkEdge(cell.cellToVert[3], cell.cellToVert[0]);

        if (crossings.size() == 2)
        {
            segments.push_back({crossings[0], crossings[1]});
        }
        else if (crossings.size() == 4)
        {
            // Ambiguous case — pair by proximity
            segments.push_back({crossings[0], crossings[1]});
            segments.push_back({crossings[2], crossings[3]});
        }
    }

    if (segments.empty())
    {
        if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
            std::cout << "[ImmersedGeometry] exactRedistance: no segments extracted; skipping.\n";
        return;
    }

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[ImmersedGeometry] exactRedistance: " << segments.size()
                  << " zero-contour segments for LS " << levelSetI << "\n";

    // Step 3: Compute exact unsigned distance from each vertex to nearest segment,
    //         then apply the sign from the current (advected) field.
    //         We use hard min (not smooth-min) because cascading smooth-min over
    //         hundreds of cell-sized segments produces a large negative bias that
    //         flips the sign of vertices near the boundary, creating bumps.

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);
        const Vector2D p = vert.coord();

        double bestDist2 = ImmersedGeometryDetail::pointSegmentDist2(p, segments[0].first, segments[0].second);
        for (std::size_t k = 1; k < segments.size(); ++k)
        {
            const double dk2 = ImmersedGeometryDetail::pointSegmentDist2(p, segments[k].first, segments[k].second);
            if (dk2 < bestDist2) bestDist2 = dk2;
        }
        const double dist = std::sqrt(bestDist2);
        const double oldSign = (levelSets[levelSetI].value[i] >= 0.0) ? 1.0 : -1.0;
        levelSets[levelSetI].value[i] = oldSign * dist;

        // Guard against exact zero
        if (std::abs(levelSets[levelSetI].value[i]) < 1e-10)
        {
            const double s = (levelSets[levelSetI].value[i] >= 0.0) ? 1.0 : -1.0;
            levelSets[levelSetI].value[i] = s * 1e-10;
        }
    }

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[ImmersedGeometry] exactRedistance: LS " << levelSetI << " done.\n";
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::exactRedistanceMovingLevelSets(const Mesh& mesh)
{
    if (hasLevelSetBoundaryConditions())
        throw std::runtime_error("[ImmersedGeometry] scalar boundary constraints do not support legacy polygonal redistancing");
    for (int levelSetI = 0; levelSetI < static_cast<int>(levelSets.size()); ++levelSetI)
    {
        if (!levelSets[levelSetI].moving())
            continue;
        exactRedistanceLevelSet(mesh, levelSetI);
    }
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::exactRedistanceMovingLevelSets(
    const Mesh& mesh,
    const std::vector<int>& activeLevelSetIndices)
{
    if (hasLevelSetBoundaryConditions())
        throw std::runtime_error("[ImmersedGeometry] scalar boundary constraints do not support legacy polygonal redistancing");
    if (activeLevelSetIndices.empty())
    {
        exactRedistanceMovingLevelSets(mesh);
        return;
    }
    for (int levelSetI : activeLevelSetIndices)
    {
        if (levelSetI < 0 || levelSetI >= static_cast<int>(levelSets.size()))
            continue;
        if (!levelSets[levelSetI].moving())
            continue;
        exactRedistanceLevelSet(mesh, levelSetI);
    }
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::updateLevelSets(const Mesh &mesh, const VertField<double>& Vx, const VertField<double>& Vy, const VertField<double>& preVx, const VertField<double>& preVy, const double &dt, bool radialAxisSymmetry, double previousTime,
    const LevelSetAdvectionVelocity* advectionVelocity, const LevelSetAdvectionScalar* advectionScalar)
{
    if (advectionVelocity)
        ImmersedGeometryBoundaryDetail::require(hasLevelSetBoundaryConditions()
            && advectionVelocity->current && advectionVelocity->previous,
            "velocity callbacks require explicit scalar boundaries and both old/new interpolants");
    if (advectionScalar)
        ImmersedGeometryBoundaryDetail::require(hasLevelSetBoundaryConditions() && advectionScalar->previous,
            "scalar interpolation requires explicit scalar boundaries and a retained old scalar");
    if (radialAxisSymmetry) ImmersedGeometryDetail::validateRadialAxisVelocity(mesh, Vx, preVx);
    std::vector<std::vector<double>> boundaryValues;
    std::vector<std::vector<TraceQueryDiagnostics>> boundaryDiagnostics;
    std::vector<LevelSetUpdateDiagnostics> freshDiagnostics;
    if (hasLevelSetBoundaryConditions())
    {
        validateBoundaryMesh(mesh);
        ImmersedGeometryBoundaryDetail::require(std::isfinite(dt) && dt >= 0.
            && std::isfinite(previousTime) && std::isfinite(previousTime+dt)
            && (dt == 0. || previousTime+dt > previousTime),
            "boundary advancement requires a resolved finite time interval");
        boundaryValues.resize(levelSets.size()); boundaryDiagnostics.resize(levelSets.size());
        freshDiagnostics.resize(levelSets.size());
        // All new boundary/trajectory callbacks and checks precede every
        // moving scalar, history and diagnostic write, including other LSs.
        for (int i = 0; i < static_cast<int>(levelSets.size()); ++i)
            if (levelSets[i].moving()) stageBoundaryLevelSet(mesh, i, Vx, Vy, preVx, preVy,
                dt, radialAxisSymmetry, previousTime, boundaryValues[i], boundaryDiagnostics[i],
                advectionVelocity, advectionScalar);
    }
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "\n[ImmersedGeometry] Updating level set...\n";

    if (!boundaryValues.empty())
    {
        traceDiagnostics_ = TraceDiagnostics{};
        levelSetUpdateDiagnostics_.swap(freshDiagnostics);
    }
    else resetUpdateDiagnostics();

    for (int levelSetI = 0; levelSetI < levelSets.size(); levelSetI++)
    {
        if (levelSets[levelSetI].moving())
        {
            updateLevelSet(mesh, levelSetI, Vx, Vy, preVx, preVy, dt, radialAxisSymmetry,
                boundaryValues.empty() ? nullptr : &boundaryValues[levelSetI],
                boundaryDiagnostics.empty() ? nullptr : &boundaryDiagnostics[levelSetI]);
        }
    }

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[ImmersedGeometry] Level set updated.\n";
}

template<typename Mesh>
template<typename Field1, typename Field2>
void ImmersedGeometry<Mesh>::updateLevelSet(const Mesh &mesh, const int levelSetI, const Field1 &Vx, const Field2 &Vy, const Field1& preVx, const Field2& preVy, const double &dt, bool radialAxisSymmetry,
    const std::vector<double>* stagedValues, const std::vector<TraceQueryDiagnostics>* stagedDiagnostics)
{
    std::vector<double> axisValues;
    std::vector<TraceQueryDiagnostics> axisDiagnostics;
    if (radialAxisSymmetry && !stagedValues)
    {
        // All new throwing checks occur before history/current-field writes
        // and outside OpenMP. Interpolate the current scalar, which is exactly
        // the history snapshot used by the original update below.
        if (!levelSets[levelSetI].value || levelSets[levelSetI].size() != mesh.vertIndexNum())
            throw std::runtime_error("[ImmersedGeometry] axis update scalar storage does not match the mesh");
        axisValues.resize(mesh.vertIndexNum());
        axisDiagnostics.resize(mesh.vertIndexNum());
        for (int i = 0; i < mesh.vertIndexNum(); ++i)
        {
            const auto point = mesh.verts(mesh.vertIndexToLabel[i]).coord();
            const auto previous = getPreCoordDetailed(mesh, point, Vx, Vy, preVx, preVy, dt,
                &axisDiagnostics[i], false, true);
            const int label = mesh.findCellLabel(previous);
            if (!ImmersedGeometryDetail::validCellLabel(mesh, label) ||
                !ImmersedGeometryDetail::validCellForInterpolation(mesh, mesh.cells(label)))
                throw std::runtime_error("[ImmersedGeometry] axis update lacks a valid scalar interpolation cell");
            axisValues[i] = Interpolation::biQuadInterpolation(mesh, levelSets[levelSetI], mesh.cells(label), previous);
            if (!std::isfinite(axisValues[i]))
                throw std::runtime_error("[ImmersedGeometry] axis update scalar interpolation is nonfinite");
        }
    }
    if (levelSetsPre[levelSetI].size() != mesh.vertIndexNum() ||
        levelSetsPre_dx[levelSetI].size() != mesh.vertIndexNum() ||
        levelSetsPre_dy[levelSetI].size() != mesh.vertIndexNum() ||
        levelSetsPrePre[levelSetI].size() != mesh.vertIndexNum() ||
        levelSetsPrePre_dx[levelSetI].size() != mesh.vertIndexNum() ||
        levelSetsPrePre_dy[levelSetI].size() != mesh.vertIndexNum() ||
        LSs_dx_[levelSetI].size() != mesh.vertIndexNum() ||
        LSs_dy_[levelSetI].size() != mesh.vertIndexNum())
    {
        syncHistoricalLevelSetStorageToCurrent(mesh, levelSetI);
    }

    if (!stagedValues)
    {
        levelSetsPrePre[levelSetI].initialize(mesh, 0.0);
        levelSetsPrePre_dx[levelSetI].initialize(mesh, 0.0);
        levelSetsPrePre_dy[levelSetI].initialize(mesh, 0.0);
    }
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        Vert &vert = mesh.verts(mesh.vertIndexToLabel[i]);

        levelSetsPrePre[levelSetI].value[i] = levelSetsPre[levelSetI].value[i];
        levelSetsPrePre_dx[levelSetI].value[i] = levelSetsPre_dx[levelSetI].value[i];
        levelSetsPrePre_dy[levelSetI].value[i] = levelSetsPre_dy[levelSetI].value[i];
        
        levelSetsPre[levelSetI].value[i] = levelSets[levelSetI].value[i];
        levelSetsPre_dx[levelSetI].value[i] = LSs_dx_[levelSetI].value[i];
        levelSetsPre_dy[levelSetI].value[i] = LSs_dy_[levelSetI].value[i];
    }

    int clampCount = 0;
    int fallbackCellCount = 0;
    int longTraceCount = 0;
    int interpolationFallbackCount = 0;
    int nonFiniteInterpolationCount = 0;
    int invalidCellTopologyCount = 0;
    int inflowBoundaryCount = 0;
    int prescribedBoundaryVertexCount = 0;
    double maxTraceRatio = 0.0;

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        Vert &vert = mesh.verts(mesh.vertIndexToLabel[i]);

        TraceQueryDiagnostics traceDiag;
        double interp;
        if (stagedValues)
        {
            traceDiag = (*stagedDiagnostics)[i];
            interp = (*stagedValues)[i];
        }
        else if (radialAxisSymmetry)
        {
            traceDiag = axisDiagnostics[i];
            interp = axisValues[i];
        }
        else
        {
            Vector2D preCoord = getPreCoordDetailed(mesh, vert.coord(), Vx, Vy, preVx, preVy, dt, &traceDiag);
            int preCellLabel = mesh.findCellLabel(preCoord);
            if (preCellLabel < 0)
            {
                preCellLabel = mesh.findCellLabel(vert.coord());
                if (preCellLabel < 0)
                    preCellLabel = 0;
#pragma omp atomic
                interpolationFallbackCount += 1;
            }
            if (!ImmersedGeometryDetail::validCellLabel(mesh, preCellLabel))
            {
                levelSets[levelSetI].value[i] = levelSetsPre[levelSetI].value[i];
#pragma omp atomic
                invalidCellTopologyCount += 1;
                continue;
            }
            const Cell& preCell = mesh.cells(preCellLabel);
            if (!ImmersedGeometryDetail::validCellForInterpolation(mesh, preCell))
            {
                levelSets[levelSetI].value[i] = levelSetsPre[levelSetI].value[i];
#pragma omp atomic
                invalidCellTopologyCount += 1;
                continue;
            }
            interp = Interpolation::biQuadInterpolation(mesh, levelSetsPre[levelSetI], preCell, preCoord);
        }
        if (std::isfinite(interp))
        {
            levelSets[levelSetI].value[i] = interp;
        }
        else
        {
            levelSets[levelSetI].value[i] = levelSetsPre[levelSetI].value[i];
#pragma omp atomic
            nonFiniteInterpolationCount += 1;
        }

#pragma omp atomic
        clampCount += traceDiag.clampCount;
#pragma omp atomic
        fallbackCellCount += traceDiag.fallbackCellCount;
#pragma omp atomic
        longTraceCount += traceDiag.longTraceCount;
#pragma omp atomic
        inflowBoundaryCount += traceDiag.inflowBoundaryCount;
#pragma omp atomic
        prescribedBoundaryVertexCount += traceDiag.prescribedBoundaryVertexCount;
#pragma omp critical
        {
            maxTraceRatio = std::max(maxTraceRatio, traceDiag.traceRatio);
        }
    }

    checkLevelSet(mesh, levelSetI);

    LevelSetUpdateDiagnostics& lsDiag = levelSetUpdateDiagnostics_[levelSetI];
    lsDiag.reinitialized = false;
    lsDiag.reinitIterations = 0;
    lsDiag.reinitFineResidual = 0.0;
    lsDiag.reinitCoarseResidual = 0.0;
    summarizeLevelSetDiagnostics(mesh, levelSetI);

    traceDiagnostics_.clampCount += clampCount;
    traceDiagnostics_.fallbackCellCount += fallbackCellCount;
    traceDiagnostics_.longTraceCount += longTraceCount;
    traceDiagnostics_.interpolationFallbackCount += interpolationFallbackCount;
    traceDiagnostics_.nonFiniteInterpolationCount += nonFiniteInterpolationCount;
    traceDiagnostics_.maxTraceRatio = std::max(traceDiagnostics_.maxTraceRatio, maxTraceRatio);
    traceDiagnostics_.inflowBoundaryCount += inflowBoundaryCount;
    traceDiagnostics_.prescribedBoundaryVertexCount += prescribedBoundaryVertexCount;

}

template<typename Mesh>
void ImmersedGeometry<Mesh>::checkLevelSet(const Mesh& mesh, const int levelSetI)
{
    const bool boundaryMode = hasLevelSetBoundaryConditions();
    if (boundaryMode)
    {
        validateBoundaryMesh(mesh);
        for (const auto& zero : levelSetBoundaryConditions_[levelSetI].fixedZeroVertices)
            for (int i = 0; i < mesh.vertIndexNum(); ++i)
                if (ImmersedGeometryBoundaryDetail::same(zero, mesh.verts(mesh.vertIndexToLabel[i]).coord()))
                    ImmersedGeometryBoundaryDetail::require(levelSets[levelSetI].value[i] == 0.,
                        "a declared fixed-zero vertex is not exactly zero");
    }
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        if (boundaryMode && prescribedBoundaryVertex(mesh, levelSetI, mesh.verts(mesh.vertIndexToLabel[i]).coord())) continue;
        if (fabs(levelSets[levelSetI].value[i]) < MINVAL)
        {
            if (levelSets[levelSetI].value[i] < 0.0)
            {
                levelSets[levelSetI].value[i] = -MINVAL;
            }
            else
            {
                levelSets[levelSetI].value[i] = MINVAL;
            }
        }
    }
}

template<typename Mesh>
Vector2D ImmersedGeometry<Mesh>::getPreCoord(const Mesh& mesh, const Vector2D& coord, const VertField<double>& Vx, const VertField<double>& Vy, const VertField<double>& preVx, const VertField<double>& preVy, const double& dt, bool radialAxisSymmetry)
{
    if (radialAxisSymmetry) ImmersedGeometryDetail::validateRadialAxisVelocity(mesh, Vx, preVx);
    return getPreCoordDetailed(mesh, coord, Vx, Vy, preVx, preVy, dt, nullptr, false, radialAxisSymmetry);
}

template<typename Mesh>
Vector2D ImmersedGeometry<Mesh>::getNextCoord(const Mesh& mesh, const Vector2D& coord, const VertField<double>& Vx, const VertField<double>& Vy, const VertField<double>& preVx, const VertField<double>& preVy, const double& dt, bool radialAxisSymmetry) const
{
    const auto finitePoint = [](const Vector2D& p) { return std::isfinite(p.x) && std::isfinite(p.y); };
    if (!finitePoint(coord) || !std::isfinite(dt) || dt < 0.0)
        throw std::runtime_error("[ImmersedGeometry] forward trace requires a finite point and nonnegative dt");
    if (radialAxisSymmetry) ImmersedGeometryDetail::validateRadialAxisVelocity(mesh, Vx, preVx);
    // The first strict evaluation validates the initial point, fields and
    // interpolation cells through the shared backtrace before accessing them.
    Vector2D next = coord;
    Vector2D traced = getPreCoordDetailed(mesh, next, Vx, Vy, preVx, preVy, dt, nullptr, true, radialAxisSymmetry);
    const Cell& originCell = mesh.cells(mesh.findCellLabel(coord));
    const double h = std::min(originCell.size().x, originCell.size().y);
    const double coordinateScale = std::max({h, std::abs(mesh.minCoord.x), std::abs(mesh.minCoord.y),
        std::abs(mesh.maxCoord.x), std::abs(mesh.maxCoord.y)});
    const double tolerance = std::max(128.0 * std::numeric_limits<double>::epsilon() * coordinateScale, 1e-12 * h);
    if (!std::isfinite(tolerance) || tolerance <= 0.0 || tolerance > 1e-8 * h)
        throw std::runtime_error("[ImmersedGeometry] forward trace coordinate scale is unresolved");
    for (int iteration = 0; iteration < 32; ++iteration)
    {
        if (iteration != 0)
            traced = getPreCoordDetailed(mesh, next, Vx, Vy, preVx, preVy, dt, nullptr, true, radialAxisSymmetry);
        const Vector2D residual = coord - traced;
        if (std::hypot(residual.x, residual.y) <= tolerance)
            return next;
        // Preserve the actual numerical backward map; this is its local
        // fixed-point inverse, not a separate forward advection formula.
        const Vector2D candidate = next + residual;
        if (!finitePoint(candidate) || std::hypot(candidate.x - coord.x, candidate.y - coord.y) > 0.5 * h)
            throw std::runtime_error("[ImmersedGeometry] forward trace exceeds half a cell displacement");
        if (candidate.x == next.x && candidate.y == next.y)
            throw std::runtime_error("[ImmersedGeometry] forward trace inverse stagnated above tolerance");
        next = candidate;
    }
    throw std::runtime_error("[ImmersedGeometry] forward trace inverse did not converge in 32 evaluations");
}

template<typename Mesh>
Vector2D ImmersedGeometry<Mesh>::getPreCoordDetailed(const Mesh& mesh, const Vector2D& coord, const VertField<double>& Vx, const VertField<double>& Vy, const VertField<double>& preVx, const VertField<double>& preVy, const double& dt, TraceQueryDiagnostics* diag, bool strict, bool radialAxisSymmetry) const
{
    const auto finitePoint = [](const Vector2D& p) { return std::isfinite(p.x) && std::isfinite(p.y); };
    if (strict || radialAxisSymmetry)
    {
        if (!finitePoint(coord) || !finitePoint(mesh.minCoord) || !finitePoint(mesh.maxCoord)
            || mesh.minCoord.x >= mesh.maxCoord.x || mesh.minCoord.y >= mesh.maxCoord.y
            || !std::isfinite(dt) || dt < 0.0)
            throw std::runtime_error("[ImmersedGeometry] strict backtrace has an invalid point/domain/dt");
        if (!Vx.value || !Vy.value || !preVx.value || !preVy.value
            || Vx.size() != mesh.vertIndexNum() || Vy.size() != mesh.vertIndexNum()
            || preVx.size() != mesh.vertIndexNum() || preVy.size() != mesh.vertIndexNum())
            throw std::runtime_error("[ImmersedGeometry] strict backtrace velocity storage does not match the mesh");
        if (radialAxisSymmetry && (mesh.minCoord.x != 0.0 || coord.x < 0.0))
            throw std::runtime_error("[ImmersedGeometry] radial-axis trace requires min.r=0 and a nonnegative radius");
    }
    auto clampToDomain = [&](Vector2D p) {
        const double epsX = 1e-12 * std::max(1.0, std::abs(mesh.maxCoord.x - mesh.minCoord.x));
        const double epsY = 1e-12 * std::max(1.0, std::abs(mesh.maxCoord.y - mesh.minCoord.y));
        const double px = p.x;
        const double py = p.y;
        if (radialAxisSymmetry && (!finitePoint(p) || p.x < 0.0))
            throw std::runtime_error("[ImmersedGeometry] radial-axis trajectory reached a negative/nonfinite radius; no clipping or reflection is allowed");
        if (strict && (!finitePoint(p) || !std::isfinite(epsX) || !std::isfinite(epsY)
            || mesh.minCoord.x + epsX >= mesh.maxCoord.x - epsX
            || mesh.minCoord.y + epsY >= mesh.maxCoord.y - epsY
            || p.x < (radialAxisSymmetry ? 0.0 : mesh.minCoord.x + epsX) || p.x > mesh.maxCoord.x - epsX
            || p.y < mesh.minCoord.y + epsY || p.y > mesh.maxCoord.y - epsY))
            throw std::runtime_error("[ImmersedGeometry] strict backtrace would clip a nonfinite or boundary point");
        p.x = radialAxisSymmetry ? std::min(p.x, mesh.maxCoord.x - epsX)
            : std::clamp(p.x, mesh.minCoord.x + epsX, mesh.maxCoord.x - epsX);
        p.y = std::clamp(p.y, mesh.minCoord.y + epsY, mesh.maxCoord.y - epsY);
        if (diag != nullptr && (std::abs(px - p.x) > 0.0 || std::abs(py - p.y) > 0.0))
            diag->clampCount += 1;
        return p;
    };
    auto safeFindCellLabel = [&](const Vector2D& p, int fallbackLabel) {
        int label = mesh.findCellLabel(clampToDomain(p));
        if ((strict || radialAxisSymmetry) && (!ImmersedGeometryDetail::validCellLabel(mesh, label)
            || !ImmersedGeometryDetail::validCellForInterpolation(mesh, mesh.cells(label))))
            throw std::runtime_error("[ImmersedGeometry] strict backtrace cannot find a valid interpolation cell");
        if (label < 0)
        {
            label = fallbackLabel;
            if (diag != nullptr)
                diag->fallbackCellCount += 1;
        }
        return label;
    };

    const int coordCellLabel = safeFindCellLabel(coord, 0);
    Cell& cell = mesh.cells(coordCellLabel);

    Vector2D vel = Vector2D(Interpolation::biQuadInterpolation(mesh, Vx, cell, coord), Interpolation::biQuadInterpolation(mesh, Vy, cell, coord));

    const auto checkStrictVelocity = [&](const Vector2D& velocity, const Cell& sampleCell) {
        if (!strict && !radialAxisSymmetry) return;
        const double h = std::min(sampleCell.size().x, sampleCell.size().y);
        if (!finitePoint(velocity) || !std::isfinite(h) || h <= 0.0)
            throw std::runtime_error("[ImmersedGeometry] strict backtrace has nonfinite velocity or invalid cell size");
        if (!strict) return;
        const double travel = std::hypot(velocity.x, velocity.y) * dt;
        if (!std::isfinite(travel) || travel > 0.5 * h)
            throw std::runtime_error("[ImmersedGeometry] strict backtrace sampled speed*dt exceeds half a cell");
    };
    checkStrictVelocity(vel, cell);
    if (radialAxisSymmetry && coord.x == 0.0 && vel.x != 0.0)
        throw std::runtime_error("[ImmersedGeometry] axis velocity interpolation is not exactly radial-zero");

    Vector2D halfCoord = clampToDomain(coord - vel * 0.5 * dt);
    const Cell& halfCell = mesh.cells(safeFindCellLabel(halfCoord, coordCellLabel));

    Vector2D preVel = Vector2D(Interpolation::biQuadInterpolation(mesh, Vx, halfCell, halfCoord), Interpolation::biQuadInterpolation(mesh, Vy, halfCell, halfCoord));
    Vector2D prePreVel = Vector2D(Interpolation::biQuadInterpolation(mesh, preVx, halfCell, halfCoord), Interpolation::biQuadInterpolation(mesh, preVy, halfCell, halfCoord));
    checkStrictVelocity(preVel, halfCell);
    checkStrictVelocity(prePreVel, halfCell);
    if (radialAxisSymmetry && halfCoord.x == 0.0 && (preVel.x != 0.0 || prePreVel.x != 0.0))
        throw std::runtime_error("[ImmersedGeometry] old/new midpoint axis velocity is not exactly radial-zero");
    Vector2D halfVel = (preVel + prePreVel) * 0.5;
    checkStrictVelocity(halfVel, halfCell);
    Vector2D traced = clampToDomain(coord - halfVel * dt);

    if (diag != nullptr)
    {
        const double localScale = std::max(cell.size().x, cell.size().y);
        if (localScale > 0.0)
        {
            diag->traceRatio = (coord - traced).mod() / localScale;
            if (diag->traceRatio > 4.0)
                diag->longTraceCount += 1;
        }
    }

    return traced;
}

template<typename Mesh>
Vector2D ImmersedGeometry<Mesh>::getPrePreCoord(const Mesh& mesh, const Vector2D& coord, const VertField<double>& Vx, const VertField<double>& Vy, const VertField<double>& preVx, const VertField<double>& preVy, const double& dt, const double& preDt)
{
    auto clampToDomain = [&](Vector2D p) {
        const double epsX = 1e-12 * std::max(1.0, std::abs(mesh.maxCoord.x - mesh.minCoord.x));
        const double epsY = 1e-12 * std::max(1.0, std::abs(mesh.maxCoord.y - mesh.minCoord.y));
        p.x = std::clamp(p.x, mesh.minCoord.x + epsX, mesh.maxCoord.x - epsX);
        p.y = std::clamp(p.y, mesh.minCoord.y + epsY, mesh.maxCoord.y - epsY);
        return p;
    };
    auto safeFindCellLabel = [&](const Vector2D& p, int fallbackLabel) {
        int label = mesh.findCellLabel(clampToDomain(p));
        if (label < 0)
            label = fallbackLabel;
        return label;
    };

    const int coordCellLabel = safeFindCellLabel(coord, 0);
    Cell& cell = mesh.cells(coordCellLabel);

    Vector2D vel = Vector2D(Interpolation::biQuadInterpolation(mesh, Vx, cell, coord), Interpolation::biQuadInterpolation(mesh, Vy, cell, coord));

    Vector2D preCoord = clampToDomain(coord - vel * dt);

    const Cell& preCell = mesh.cells(safeFindCellLabel(preCoord, coordCellLabel));

    Vector2D preVel = Vector2D(Interpolation::biQuadInterpolation(mesh, preVx, preCell, preCoord), Interpolation::biQuadInterpolation(mesh, preVy, preCell, preCoord));

    return clampToDomain(coord - preVel * (dt + preDt));
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::checkCutInfo(Mesh& mesh, const bool isInit)
{
    if (hasLevelSetBoundaryConditions())
        throw std::runtime_error("[ImmersedGeometry] legacy cut topology does not support declared scalar boundaries; use shared Hermite geometry");
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[ImmersedGeometry] Checking cut-cell/interface information...\n";
    const std::size_t nLS = static_cast<std::size_t>(mesh.levelSetCount());
    for (int i = 0; i < mesh.cellIndexNum(); i++)
    {
        const int cellLabel = mesh.cellIndexToLabel[i];
        if (!ImmersedGeometryDetail::validCellLabel(mesh, cellLabel))
        {
            std::ostringstream oss;
            oss << "[ImmersedGeometry] invalid cellIndexToLabel entry while resetting cut info: "
                << "cellIndex=" << i << " cellLabel=" << cellLabel;
            std::cerr << oss.str() << std::endl;
            throw std::runtime_error(oss.str());
        }
        Cell& cell = mesh.cells(cellLabel);
        cell.ensureLevelSetSize(nLS);
        for (std::size_t j = 0; j < nLS; j++)
            cell.isCut[j] = 0u;
    }
    for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
    {
        const int hEdgeLabel = mesh.hEdgeIndexToLabel[i];
        if (!ImmersedGeometryDetail::validHEdgeLabel(mesh, hEdgeLabel))
        {
            std::ostringstream oss;
            oss << "[ImmersedGeometry] invalid hEdgeIndexToLabel entry while resetting cut info: "
                << "edgeIndex=" << i << " edgeLabel=" << hEdgeLabel;
            std::cerr << oss.str() << std::endl;
            throw std::runtime_error(oss.str());
        }
        Edge& hEdge = mesh.hEdges(hEdgeLabel);
        hEdge.ensureLevelSetSize(nLS);
        for (std::size_t j = 0; j < nLS; j++)
            hEdge.isCut[j] = 0u;
    }
    for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
    {
        const int vEdgeLabel = mesh.vEdgeIndexToLabel[i];
        if (!ImmersedGeometryDetail::validVEdgeLabel(mesh, vEdgeLabel))
        {
            std::ostringstream oss;
            oss << "[ImmersedGeometry] invalid vEdgeIndexToLabel entry while resetting cut info: "
                << "edgeIndex=" << i << " edgeLabel=" << vEdgeLabel;
            std::cerr << oss.str() << std::endl;
            throw std::runtime_error(oss.str());
        }
        Edge& vEdge = mesh.vEdges(vEdgeLabel);
        vEdge.ensureLevelSetSize(nLS);
        for (std::size_t j = 0; j < nLS; j++)
            vEdge.isCut[j] = 0u;
    }
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        const int vertLabel = mesh.vertIndexToLabel[i];
        if (!ImmersedGeometryDetail::validVertLabel(mesh, vertLabel))
        {
            std::ostringstream oss;
            oss << "[ImmersedGeometry] invalid vertIndexToLabel entry while resetting cut info: "
                << "vertIndex=" << i << " vertLabel=" << vertLabel;
            std::cerr << oss.str() << std::endl;
            throw std::runtime_error(oss.str());
        }
        Vert& vert = mesh.verts(vertLabel);
        vert.ensureLevelSetSize(nLS);
        for (std::size_t j = 0; j < nLS; j++)
        {
            vert.onCutEdge[j] = 0u;
            vert.onCutCell[j] = 0u;
            vert.nearestNodeCoord[j] = Vector2D(0.0, 0.0);
            vert.nearestNodeLabel[j] = -1;
        }
    }

    for (int levelSetI = 0; levelSetI < levelSets.size(); levelSetI++)
    {
        if (!levelSets[levelSetI].moving() && !isInit)
        {
            continue;
        }
        mesh.posFaceVerts[levelSetI].clear();
        mesh.negFaceVerts[levelSetI].clear();

        for (int i = 0; i < mesh.cellIndexNum(); i++)
        {
            const int cellLabel = mesh.cellIndexToLabel[i];
            if (!ImmersedGeometryDetail::validCellLabel(mesh, cellLabel))
            {
                std::ostringstream oss;
                oss << "[ImmersedGeometry] invalid cellIndexToLabel entry while checking cut info: "
                    << "levelSet=" << levelSetI
                    << " cellIndex=" << i
                    << " cellLabel=" << cellLabel;
                std::cerr << oss.str() << std::endl;
                throw std::runtime_error(oss.str());
            }
            Cell& cell = mesh.cells(cellLabel);

            if (!ImmersedGeometryDetail::validVertLabel(mesh, cell.cellToVert[0]) ||
                !ImmersedGeometryDetail::validVertLabel(mesh, cell.cellToVert[1]) ||
                !ImmersedGeometryDetail::validVertLabel(mesh, cell.cellToVert[2]) ||
                !ImmersedGeometryDetail::validVertLabel(mesh, cell.cellToVert[3]))
            {
                std::ostringstream oss;
                oss << "[ImmersedGeometry] invalid cellToVert mapping while checking cut info: "
                    << "levelSet=" << levelSetI
                    << " cellLabel=" << cell.label()
                    << " cellIndex=" << cell.index()
                    << " verts=(" << cell.cellToVert[0] << ',' << cell.cellToVert[1]
                    << ',' << cell.cellToVert[2] << ',' << cell.cellToVert[3] << ')';
                std::cerr << oss.str() << std::endl;
                throw std::runtime_error(oss.str());
            }

            Vert& vert0 = mesh.verts(cell.cellToVert[0]);
            Vert& vert1 = mesh.verts(cell.cellToVert[1]);
            Vert& vert2 = mesh.verts(cell.cellToVert[2]);
            Vert& vert3 = mesh.verts(cell.cellToVert[3]);

            if (!((levelSets[levelSetI].value[vert0.index()] > 0.0 
                && levelSets[levelSetI].value[vert1.index()] > 0.0 
                && levelSets[levelSetI].value[vert2.index()] > 0.0 
                && levelSets[levelSetI].value[vert3.index()] > 0.0) 
                ||
                (levelSets[levelSetI].value[vert0.index()] <= 0.0 
                && levelSets[levelSetI].value[vert1.index()] <= 0.0 
                && levelSets[levelSetI].value[vert2.index()] <= 0.0 
                && levelSets[levelSetI].value[vert3.index()] <= 0.0)))
            {
                cell.isCut[levelSetI] = true;
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        const int hEdgeLabel = cell.cellToHEdge[j][k];
                        if (hEdgeLabel != -1)
                        {
                            if (!ImmersedGeometryDetail::validHEdgeLabel(mesh, hEdgeLabel))
                            {
                                std::ostringstream oss;
                                oss << "[ImmersedGeometry] invalid cellToHEdge mapping while checking cut info: "
                                    << "levelSet=" << levelSetI
                                    << " cellLabel=" << cell.label()
                                    << " cellIndex=" << cell.index()
                                    << " slot=(" << j << ',' << k << ')'
                                    << " edgeLabel=" << hEdgeLabel;
                                std::cerr << oss.str() << std::endl;
                                throw std::runtime_error(oss.str());
                            }

                            Edge& hEdge = mesh.hEdges(hEdgeLabel);
                            if (!ImmersedGeometryDetail::validVertLabel(mesh, hEdge.edgeToVert[0]) ||
                                !ImmersedGeometryDetail::validVertLabel(mesh, hEdge.edgeToVert[1]))
                            {
                                std::ostringstream oss;
                                oss << "[ImmersedGeometry] invalid hEdge.edgeToVert mapping from cut cell: "
                                    << "levelSet=" << levelSetI
                                    << " cellLabel=" << cell.label()
                                    << " edgeLabel=" << hEdge.label()
                                    << " edgeIndex=" << hEdge.index()
                                    << " verts=(" << hEdge.edgeToVert[0] << ',' << hEdge.edgeToVert[1] << ')';
                                std::cerr << oss.str() << std::endl;
                                throw std::runtime_error(oss.str());
                            }

                            mesh.verts(hEdge.edgeToVert[0]).onCutCell[levelSetI] = true;
                            mesh.verts(hEdge.edgeToVert[1]).onCutCell[levelSetI] = true;
                        }

                        const int vEdgeLabel = cell.cellToVEdge[j][k];
                        if (vEdgeLabel != -1)
                        {
                            if (!ImmersedGeometryDetail::validVEdgeLabel(mesh, vEdgeLabel))
                            {
                                std::ostringstream oss;
                                oss << "[ImmersedGeometry] invalid cellToVEdge mapping while checking cut info: "
                                    << "levelSet=" << levelSetI
                                    << " cellLabel=" << cell.label()
                                    << " cellIndex=" << cell.index()
                                    << " slot=(" << j << ',' << k << ')'
                                    << " edgeLabel=" << vEdgeLabel;
                                std::cerr << oss.str() << std::endl;
                                throw std::runtime_error(oss.str());
                            }

                            Edge& vEdge = mesh.vEdges(vEdgeLabel);
                            if (!ImmersedGeometryDetail::validVertLabel(mesh, vEdge.edgeToVert[0]) ||
                                !ImmersedGeometryDetail::validVertLabel(mesh, vEdge.edgeToVert[1]))
                            {
                                std::ostringstream oss;
                                oss << "[ImmersedGeometry] invalid vEdge.edgeToVert mapping from cut cell: "
                                    << "levelSet=" << levelSetI
                                    << " cellLabel=" << cell.label()
                                    << " edgeLabel=" << vEdge.label()
                                    << " edgeIndex=" << vEdge.index()
                                    << " verts=(" << vEdge.edgeToVert[0] << ',' << vEdge.edgeToVert[1] << ')';
                                std::cerr << oss.str() << std::endl;
                                throw std::runtime_error(oss.str());
                            }

                            mesh.verts(vEdge.edgeToVert[0]).onCutCell[levelSetI] = true;
                            mesh.verts(vEdge.edgeToVert[1]).onCutCell[levelSetI] = true;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
        {
            const int hEdgeLabel = mesh.hEdgeIndexToLabel[i];
            if (!ImmersedGeometryDetail::validHEdgeLabel(mesh, hEdgeLabel))
            {
                std::ostringstream oss;
                oss << "[ImmersedGeometry] invalid hEdgeIndexToLabel entry while checking cut info: "
                    << "levelSet=" << levelSetI
                    << " edgeIndex=" << i
                    << " edgeLabel=" << hEdgeLabel;
                std::cerr << oss.str() << std::endl;
                throw std::runtime_error(oss.str());
            }
            Edge& hEdge = mesh.hEdges(hEdgeLabel);

            if (!ImmersedGeometryDetail::validVertLabel(mesh, hEdge.edgeToVert[0]) ||
                !ImmersedGeometryDetail::validVertLabel(mesh, hEdge.edgeToVert[1]))
            {
                std::ostringstream oss;
                oss << "[ImmersedGeometry] invalid hEdge.edgeToVert mapping while checking cut info: "
                    << "levelSet=" << levelSetI
                    << " edgeLabel=" << hEdge.label()
                    << " edgeIndex=" << hEdge.index()
                    << " verts=(" << hEdge.edgeToVert[0] << ',' << hEdge.edgeToVert[1] << ')';
                std::cerr << oss.str() << std::endl;
                throw std::runtime_error(oss.str());
            }

            Vert& vert0 = mesh.verts(hEdge.edgeToVert[0]);
            Vert& vert1 = mesh.verts(hEdge.edgeToVert[1]);

            double value0 = levelSets[levelSetI].value[vert0.index()];
            double value1 = levelSets[levelSetI].value[vert1.index()];

            if (value0 * value1 <= 0.0)
            {
                hEdge.isCut[levelSetI] = true;

                vert0.onCutEdge[levelSetI] = true;
                vert1.onCutEdge[levelSetI] = true;
            }
        }

        for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
        {
            const int vEdgeLabel = mesh.vEdgeIndexToLabel[i];
            if (!ImmersedGeometryDetail::validVEdgeLabel(mesh, vEdgeLabel))
            {
                std::ostringstream oss;
                oss << "[ImmersedGeometry] invalid vEdgeIndexToLabel entry while checking cut info: "
                    << "levelSet=" << levelSetI
                    << " edgeIndex=" << i
                    << " edgeLabel=" << vEdgeLabel;
                std::cerr << oss.str() << std::endl;
                throw std::runtime_error(oss.str());
            }
            Edge& vEdge = mesh.vEdges(vEdgeLabel);

            if (!ImmersedGeometryDetail::validVertLabel(mesh, vEdge.edgeToVert[0]) ||
                !ImmersedGeometryDetail::validVertLabel(mesh, vEdge.edgeToVert[1]))
            {
                std::ostringstream oss;
                oss << "[ImmersedGeometry] invalid vEdge.edgeToVert mapping while checking cut info: "
                    << "levelSet=" << levelSetI
                    << " edgeLabel=" << vEdge.label()
                    << " edgeIndex=" << vEdge.index()
                    << " verts=(" << vEdge.edgeToVert[0] << ',' << vEdge.edgeToVert[1] << ')';
                std::cerr << oss.str() << std::endl;
                throw std::runtime_error(oss.str());
            }

            Vert& vert0 = mesh.verts(vEdge.edgeToVert[0]);
            Vert& vert1 = mesh.verts(vEdge.edgeToVert[1]);

            double value0 = levelSets[levelSetI].value[vert0.index()];
            double value1 = levelSets[levelSetI].value[vert1.index()];

            if (value0 * value1 <= 0.0)
            {
                vEdge.isCut[levelSetI] = true;

                vert0.onCutEdge[levelSetI] = true;
                vert1.onCutEdge[levelSetI] = true;
            }
        }

        for (int i = 0; i < mesh.vertIndexNum(); i++)
        {
            const int vertLabel = mesh.vertIndexToLabel[i];
            if (!ImmersedGeometryDetail::validVertLabel(mesh, vertLabel))
            {
                std::ostringstream oss;
                oss << "[ImmersedGeometry] invalid vertIndexToLabel entry while checking cut info: "
                    << "levelSet=" << levelSetI
                    << " vertIndex=" << i
                    << " vertLabel=" << vertLabel;
                std::cerr << oss.str() << std::endl;
                throw std::runtime_error(oss.str());
            }
            Vert& vert = mesh.verts(vertLabel);

            if (vert.onCutEdge[levelSetI])
            {
                if (levelSets[levelSetI].value[vert.index()] > 0.0)
                {
                    mesh.posFaceVerts[levelSetI].push_back(vert.label());
                }
                else
                {
                    mesh.negFaceVerts[levelSetI].push_back(vert.label());
                    double dxL = LSs_dx_[levelSetI].value[vert.index()];
                    double dyL = LSs_dy_[levelSetI].value[vert.index()];
                    double mod_inv = 1.0 / sqrt(dxL * dxL + dyL * dyL);
                    Vector2D normal(dxL * mod_inv, dyL * mod_inv);
                    vert.nearestNodeCoord[levelSetI] = vert.coord() - normal * levelSets[levelSetI].value[vert.index()];
                }
            }
            else if (vert.supCell != -1)
            {
                const int supportCellLabel = ImmersedGeometryDetail::recoverSupportCellLabel(mesh, vert);
                if (!ImmersedGeometryDetail::validCellLabel(mesh, supportCellLabel))
                    continue;
                vert.supCell = supportCellLabel;
                Cell& cell = mesh.cells(supportCellLabel);
                if (cell.isCut[levelSetI])
                {
                    if (levelSets[levelSetI].value[vert.index()] > 0.0)
                    {
                        mesh.posFaceVerts[levelSetI].push_back(vert.label());
                    }
                    else
                    {
                        mesh.negFaceVerts[levelSetI].push_back(vert.label());
                        double dxL = LSs_dx_[levelSetI].value[vert.index()];
                        double dyL = LSs_dy_[levelSetI].value[vert.index()];
                        const double grad2 = dxL * dxL + dyL * dyL;
                        if (std::isfinite(grad2) && grad2 > MINVAL * MINVAL)
                        {
                            double mod_inv = 1.0 / sqrt(grad2);
                            Vector2D normal(dxL * mod_inv, dyL * mod_inv);
                            vert.nearestNodeCoord[levelSetI] =
                                vert.coord() - normal * levelSets[levelSetI].value[vert.index()];
                        }
                        else
                        {
                            // Degenerate local gradient: keep a safe fallback instead of
                            // propagating invalid geometry coordinates into later cut metrics.
                            vert.nearestNodeCoord[levelSetI] = vert.coord();
                        }
                    }
                }
            }
        }
    }

    cutInfoDiagnostics_ = CutInfoDiagnostics{};
    for (int i = 0; i < mesh.cellIndexNum(); i++)
    {
        const int cellLabel = mesh.cellIndexToLabel[i];
        if (!ImmersedGeometryDetail::validCellLabel(mesh, cellLabel))
            continue;
        const Cell& cell = mesh.cells(cellLabel);
        bool anyCut = false;
        for (std::size_t ls = 0; ls < cell.isCut.size(); ++ls)
        {
            if (cell.isCut[ls])
            {
                anyCut = true;
                break;
            }
        }
        if (anyCut)
            cutInfoDiagnostics_.cutCellCount += 1;
    }
    for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
    {
        const int hEdgeLabel = mesh.hEdgeIndexToLabel[i];
        if (!ImmersedGeometryDetail::validHEdgeLabel(mesh, hEdgeLabel))
            continue;
        const Edge& hEdge = mesh.hEdges(hEdgeLabel);
        bool anyCut = false;
        for (std::size_t ls = 0; ls < hEdge.isCut.size(); ++ls)
        {
            if (hEdge.isCut[ls])
            {
                anyCut = true;
                break;
            }
        }
        if (anyCut)
            cutInfoDiagnostics_.cutHEdgeCount += 1;
    }
    for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
    {
        const int vEdgeLabel = mesh.vEdgeIndexToLabel[i];
        if (!ImmersedGeometryDetail::validVEdgeLabel(mesh, vEdgeLabel))
            continue;
        const Edge& vEdge = mesh.vEdges(vEdgeLabel);
        bool anyCut = false;
        for (std::size_t ls = 0; ls < vEdge.isCut.size(); ++ls)
        {
            if (vEdge.isCut[ls])
            {
                anyCut = true;
                break;
            }
        }
        if (anyCut)
            cutInfoDiagnostics_.cutVEdgeCount += 1;
    }

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[ImmersedGeometry] Cut-cell/interface information check completed.\n";
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::wirteCache(const std::string& filePath) const
{
    if (hasLevelSetBoundaryConditions())
        throw std::runtime_error("[ImmersedGeometry] geometry cache v2 cannot store scalar boundary callbacks/constraints");
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        cout << "[ImmersedGeometry] Writing geometry cache..." << endl;

    ostringstream fileName;
	fileName << filePath << "iGeometry_cache.bin";

    auto cache = FileIO::openOutputFile(fileName.str(), ios::binary | ios::out);

    cache.write(reinterpret_cast<const char*>(ImmersedGeometryCache::kMagic), 4);
    const std::uint32_t version = ImmersedGeometryCache::kVersion;
    cache.write(reinterpret_cast<const char*>(&version), sizeof(std::uint32_t));
    const int levelSetCount = static_cast<int>(levelSets.size());
    cache.write(reinterpret_cast<const char*>(&levelSetCount), sizeof(int));

    // points
    int pointNum = points.size();
    cache.write(reinterpret_cast<const char*>(&pointNum), sizeof(int));
    cache.write(reinterpret_cast<const char*>(points.data()), pointNum * sizeof(Vector2D));

    // lines
    int lineNum = lines.size();
    cache.write(reinterpret_cast<const char*>(&lineNum), sizeof(int));
    for (const Line& line : lines)
    {
        cache.write(reinterpret_cast<const char*>(&line.points), 2 * sizeof(Vector2D));
        cache.write(reinterpret_cast<const char*>(&line.levelSetIndex), sizeof(int));
        cache.write(reinterpret_cast<const char*>(&line.phaseIndex), sizeof(int));
        cache.write(reinterpret_cast<const char*>(&line.opposedPhaseIndex), sizeof(int));
        LineType lineType = line.type();
        cache.write(reinterpret_cast<const char*>(&lineType), sizeof(LineType));
        Vector2D arcCenter = line.arcCenter();
        cache.write(reinterpret_cast<const char*>(&arcCenter), sizeof(Vector2D));
        double radius = line.radius();
        cache.write(reinterpret_cast<const char*>(&radius), sizeof(double));
        double radian = line.radian();
        cache.write(reinterpret_cast<const char*>(&radian), sizeof(double));
        int featurePointsNum = line.featurePoints().size();
        cache.write(reinterpret_cast<const char*>(&featurePointsNum), sizeof(int));
        cache.write(reinterpret_cast<const char*>(line.featurePoints().data()), featurePointsNum * sizeof(Vector2D));
        GeometricBoundary geoBoundary = line.geoBoundary();
        cache.write(reinterpret_cast<const char*>(&geoBoundary), sizeof(GeometricBoundary));
        bool ignored = line.ignored();
        cache.write(reinterpret_cast<const char*>(&ignored), sizeof(bool));
        bool immersed = line.immersed();
        cache.write(reinterpret_cast<const char*>(&immersed), sizeof(bool));
    }

    // levelsets
    for (const LevelSet& levelSet : levelSets)
    {
        int lineNum = levelSet.lines().size();
        bool isMoving = levelSet.moving();
        cache.write(reinterpret_cast<const char*>(&lineNum), sizeof(int));
        cache.write(reinterpret_cast<const char*>(levelSet.lines().data()), lineNum * sizeof(SignInt));
        cache.write(reinterpret_cast<const char*>(&isMoving), sizeof(bool));
        int size = levelSet.size();
        cache.write(reinterpret_cast<const char*>(&size), sizeof(int));
        for (int i = 0; i < size; i++)
        {
            cache.write(reinterpret_cast<const char*>(&levelSet.value[i]), sizeof(double));
        }
    }
    // 场
    int fieldNum = levelSetsPre.size();
    cache.write(reinterpret_cast<const char*>(&fieldNum), sizeof(int));
    for (int i = 0; i < levelSetsPre.size(); i++)
    {
        int vertIndexNumPre = levelSetsPre[i].size();
        cache.write(reinterpret_cast<const char*>(&vertIndexNumPre), sizeof(int));
        for (int j = 0; j < levelSetsPre[i].size(); j++)
        {
            int label = levelSetsPre[i].label(j);
            cache.write(reinterpret_cast<const char*>(&label), sizeof(int));
            cache.write(reinterpret_cast<const char*>(&levelSetsPre[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&levelSetsPre_dx[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&levelSetsPre_dy[i].value[j]), sizeof(double));
        }
        int vertIndexNumPrePre = levelSetsPrePre[i].size();
        cache.write(reinterpret_cast<const char*>(&vertIndexNumPrePre), sizeof(int));
        for (int j = 0; j < levelSetsPrePre[i].size(); j++)
        {
            int label = levelSetsPrePre[i].label(j);
            cache.write(reinterpret_cast<const char*>(&label), sizeof(int));
            cache.write(reinterpret_cast<const char*>(&levelSetsPrePre[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&levelSetsPrePre_dx[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&levelSetsPrePre_dy[i].value[j]), sizeof(double));
        }
        int vertIndexNum = LSs_dx_[i].size();
        cache.write(reinterpret_cast<const char*>(&vertIndexNum), sizeof(int));
        for (int j = 0; j < vertIndexNum; j++)
        {
            int label = LSs_dx_[i].label(j);
            cache.write(reinterpret_cast<const char*>(&label), sizeof(int));
            cache.write(reinterpret_cast<const char*>(&LSs_dx_[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&LSs_dy_[i].value[j]), sizeof(double));
        }
    }

    // phases
    int phaseNum = phases.size();
    cache.write(reinterpret_cast<const char*>(&phaseNum), sizeof(int));
    for (const Phase& phase : phases)
    {
        const std::string& name = phase.name();
        int nameLen = name.size();
        cache.write(reinterpret_cast<const char*>(&nameLen), sizeof(int));
        cache.write(name.data(), nameLen);
        
        int operationSize = phase.operation().size();
        cache.write(reinterpret_cast<const char*>(&operationSize), sizeof(int));
        cache.write(reinterpret_cast<const char*>(phase.operation().data()), operationSize * sizeof(LevelSetOperation));

        vector<SignInt> levelSetIndexA = phase.levelSetIndexA();
        int levelSetIndexASize = levelSetIndexA.size();
        cache.write(reinterpret_cast<const char*>(&levelSetIndexASize), sizeof(int));
        cache.write(reinterpret_cast<const char*>(levelSetIndexA.data()), levelSetIndexASize * sizeof(SignInt));

        vector<SignInt> levelSetIndexB = phase.levelSetIndexB();
        int levelSetIndexBSize = levelSetIndexB.size();
        cache.write(reinterpret_cast<const char*>(&levelSetIndexBSize), sizeof(int));
        cache.write(reinterpret_cast<const char*>(levelSetIndexB.data()), levelSetIndexBSize * sizeof(SignInt));

        int boundaryLineNum = phase.boundaryLines().size();
        cache.write(reinterpret_cast<const char*>(&boundaryLineNum), sizeof(int));
        cache.write(reinterpret_cast<const char*>(phase.boundaryLines().data()), boundaryLineNum * sizeof(SignInt));
    }

    for (int i = 0; i < static_cast<int>(LSs_dxx_.size()); i++)
    {
        int diff_size = LSs_dxx_[i].size();
        cache.write(reinterpret_cast<const char*>(&diff_size), sizeof(int));
        for (int j = 0; j < diff_size; j++)
        {
            cache.write(reinterpret_cast<const char*>(&LSs_dxx_[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&LSs_dyy_[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&LSs_dxy_[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&LSs_dx_backward_[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&LSs_dy_backward_[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&LSs_dx_forward_[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&LSs_dy_forward_[i].value[j]), sizeof(double));
        }
    }

    FileIO::finishOutputFile(cache, fileName.str());
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        cout << "[ImmersedGeometry] Geometry cache written." << endl;
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::readCache(Mesh& mesh, const std::string& filePath)
{
    if (hasLevelSetBoundaryConditions())
        throw std::runtime_error("[ImmersedGeometry] geometry cache v2 cannot restore scalar boundary callbacks/constraints");
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        cout << "[ImmersedGeometry] Reading geometry cache..." << endl;

    ostringstream fileName;
    fileName << filePath << "iGeometry_cache.bin";

    ifstream cache(fileName.str(), ios::binary | ios::in);
    if (cache.fail())
    {
        cout << "Reading cache \"" << fileName.str() << "\" failed!" << endl;
        exit(-1);
    }

    char magic[4] = {0, 0, 0, 0};
    cache.read(magic, 4);
    if (!ImmersedGeometryCache::MagicEquals(magic, ImmersedGeometryCache::kMagic))
    {
        cout << "[ImmersedGeometry] ERROR: incompatible cache format in \"" << fileName.str()
             << "\" (delete cache)." << endl;
        exit(-1);
    }
    std::uint32_t version = 0;
    cache.read(reinterpret_cast<char*>(&version), sizeof(std::uint32_t));
    if (version != ImmersedGeometryCache::kVersion)
    {
        cout << "[ImmersedGeometry] ERROR: unsupported cache version " << version
             << " in \"" << fileName.str() << "\"." << endl;
        exit(-1);
    }
    int levelSetCount = 0;
    cache.read(reinterpret_cast<char*>(&levelSetCount), sizeof(int));
    if (levelSetCount < 0)
        levelSetCount = 0;

    // points
    int pointNum = 0;
    cache.read(reinterpret_cast<char*>(&pointNum), sizeof(int));
    points.resize(pointNum);
    cache.read(reinterpret_cast<char*>(points.data()), pointNum * sizeof(Vector2D));

    // lines
    int lineNum = 0;
    cache.read(reinterpret_cast<char*>(&lineNum), sizeof(int));
    lines.resize(lineNum);
    for (int i = 0; i < lineNum; i++)
    {
        Line& line = lines[i];
        cache.read(reinterpret_cast<char*>(&line.points), 2 * sizeof(Vector2D));
        cache.read(reinterpret_cast<char*>(&line.levelSetIndex), sizeof(int));
        cache.read(reinterpret_cast<char*>(&line.phaseIndex), sizeof(int));
        cache.read(reinterpret_cast<char*>(&line.opposedPhaseIndex), sizeof(int));
        LineType lineType;
        cache.read(reinterpret_cast<char*>(&lineType), sizeof(LineType));
        line.setType(lineType);
        Vector2D arcCenter;
        cache.read(reinterpret_cast<char*>(&arcCenter), sizeof(Vector2D));
        line.setArcCenter(arcCenter);
        double radius = 0.0;
        cache.read(reinterpret_cast<char*>(&radius), sizeof(double));
        line.setRadius(radius);
        double radian = 0.0;
        cache.read(reinterpret_cast<char*>(&radian), sizeof(double));
        line.setRadian(radian);
        int featurePointsNum = 0;
        cache.read(reinterpret_cast<char*>(&featurePointsNum), sizeof(int));
        line.featurePoints().resize(featurePointsNum);
        cache.read(reinterpret_cast<char*>(line.featurePoints().data()), featurePointsNum * sizeof(Vector2D));
        GeometricBoundary geoBoundary;
        cache.read(reinterpret_cast<char*>(&geoBoundary), sizeof(GeometricBoundary));
        line.setGeoBoundary(geoBoundary);
        bool ignored = false;
        cache.read(reinterpret_cast<char*>(&ignored), sizeof(bool));
        line.setIgnored(ignored);
        bool immersed = false;
        cache.read(reinterpret_cast<char*>(&immersed), sizeof(bool));
        line.setImmersed(immersed);
    }

    // levelsets
    levelSets.resize(levelSetCount);
    mesh.negFaceVerts.resize(levelSetCount);
    mesh.posFaceVerts.resize(levelSetCount);
    for (int i = 0; i < levelSetCount; i++)
    {
        int lineNum = 0;
        bool isMoving = false;
        cache.read(reinterpret_cast<char*>(&lineNum), sizeof(int));
        vector<SignInt> lines(lineNum);
        cache.read(reinterpret_cast<char*>(lines.data()), lineNum * sizeof(SignInt));
        cache.read(reinterpret_cast<char*>(&isMoving), sizeof(bool));
        levelSets[i].setLines(lines);
        levelSets[i].setMoving(isMoving);
        int size = 0;
        cache.read(reinterpret_cast<char*>(&size), sizeof(int));
        levelSets[i].initialize(mesh, 0.0);
        for (int j = 0; j < size; j++)
        {
            cache.read(reinterpret_cast<char*>(&levelSets[i].value[j]), sizeof(double));
        }
    }

    // 场
    int fieldNum = 0;
    cache.read(reinterpret_cast<char*>(&fieldNum), sizeof(int));
    levelSetsPre.resize(fieldNum);
    levelSetsPre_dx.resize(fieldNum);
    levelSetsPre_dy.resize(fieldNum);
    levelSetsPrePre.resize(fieldNum);
    levelSetsPrePre_dx.resize(fieldNum);
    levelSetsPrePre_dy.resize(fieldNum);
    LSs_dx_.resize(fieldNum);
    LSs_dy_.resize(fieldNum);
    LSs_dxx_.resize(fieldNum);
    LSs_dyy_.resize(fieldNum);
    LSs_dxy_.resize(fieldNum);
    LSs_dx_backward_.resize(fieldNum);
    LSs_dy_backward_.resize(fieldNum);
    LSs_dx_forward_.resize(fieldNum);
    LSs_dy_forward_.resize(fieldNum);
    for (int i = 0; i < fieldNum; i++)
    {
        int vertIndexNumPre = 0;
        cache.read(reinterpret_cast<char*>(&vertIndexNumPre), sizeof(int));
        levelSetsPre[i].initialize(mesh, 0.0);
        levelSetsPre_dx[i].initialize(mesh, 0.0);
        levelSetsPre_dy[i].initialize(mesh, 0.0);
        for (int j = 0; j < vertIndexNumPre; j++)
        {
            int label = 0;
            cache.read(reinterpret_cast<char*>(&label), sizeof(int));
            levelSetsPre[i].setLabel(j, label);
            levelSetsPre_dx[i].setLabel(j, label);
            levelSetsPre_dy[i].setLabel(j, label);
            cache.read(reinterpret_cast<char*>(&levelSetsPre[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&levelSetsPre_dx[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&levelSetsPre_dy[i].value[j]), sizeof(double));
        }
        int vertIndexNumPrePre = 0;
        cache.read(reinterpret_cast<char*>(&vertIndexNumPrePre), sizeof(int));
        levelSetsPrePre[i].initialize(mesh, 0.0);
        levelSetsPrePre_dx[i].initialize(mesh, 0.0);
        levelSetsPrePre_dy[i].initialize(mesh, 0.0);
        for (int j = 0; j < vertIndexNumPrePre; j++)
        {
            int label = 0;
            cache.read(reinterpret_cast<char*>(&label), sizeof(int));
            levelSetsPrePre[i].setLabel(j, label);
            levelSetsPrePre_dx[i].setLabel(j, label);
            levelSetsPrePre_dy[i].setLabel(j, label);
            cache.read(reinterpret_cast<char*>(&levelSetsPrePre[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&levelSetsPrePre_dx[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&levelSetsPrePre_dy[i].value[j]), sizeof(double));
        }
        int vertIndexNum = 0;
        cache.read(reinterpret_cast<char*>(&vertIndexNum), sizeof(int));
        LSs_dx_[i].initialize(mesh, 0.0);
        LSs_dy_[i].initialize(mesh, 0.0);
        for (int j = 0; j < vertIndexNum; j++)
        {
            int label = 0;
            cache.read(reinterpret_cast<char*>(&label), sizeof(int));
            LSs_dx_[i].setLabel(j, label);
            LSs_dy_[i].setLabel(j, label);
            cache.read(reinterpret_cast<char*>(&LSs_dx_[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&LSs_dy_[i].value[j]), sizeof(double));
        }
    }

    // phases
    int phaseNum = 0;
    cache.read(reinterpret_cast<char*>(&phaseNum), sizeof(int));
    phases.resize(phaseNum);
    for (Phase& phase : phases)
    {
        int nameLen = 0;
        cache.read(reinterpret_cast<char*>(&nameLen), sizeof(int));
        std::string name(nameLen, '\0');
        cache.read(&name[0], nameLen);
        phase.setName(name);

        int operationSize = 0;
        cache.read(reinterpret_cast<char*>(&operationSize), sizeof(int));
        vector<LevelSetOperation> operations(operationSize);
        cache.read(reinterpret_cast<char*>(operations.data()), operationSize * sizeof(LevelSetOperation));
        phase.setOperation(operations);

        vector<SignInt> levelSetIndexA;
        int levelSetIndexASize = 0;
        cache.read(reinterpret_cast<char*>(&levelSetIndexASize), sizeof(int));
        levelSetIndexA.resize(levelSetIndexASize);
        cache.read(reinterpret_cast<char*>(levelSetIndexA.data()), levelSetIndexASize * sizeof(SignInt));
        phase.setLevelSetIndexA(levelSetIndexA);

        vector<SignInt> levelSetIndexB;
        int levelSetIndexBSize = 0;
        cache.read(reinterpret_cast<char*>(&levelSetIndexBSize), sizeof(int));
        levelSetIndexB.resize(levelSetIndexBSize);
        cache.read(reinterpret_cast<char*>(levelSetIndexB.data()), levelSetIndexBSize * sizeof(SignInt));
        phase.setLevelSetIndexB(levelSetIndexB);

        vector<SignInt> boundaryLines;
        int boundaryLineNum = 0;
        cache.read(reinterpret_cast<char*>(&boundaryLineNum), sizeof(int));
        boundaryLines.resize(boundaryLineNum);
        cache.read(reinterpret_cast<char*>(boundaryLines.data()), boundaryLineNum * sizeof(SignInt));
        phase.setBoundaryLines(boundaryLines);
    }

    for (int i = 0; i < static_cast<int>(LSs_dxx_.size()); i++)
    {
        int diff_size = 0;
        cache.read(reinterpret_cast<char*>(&diff_size), sizeof(int));
        LSs_dxx_[i].initialize(diff_size, 0.0);
        LSs_dyy_[i].initialize(diff_size, 0.0);
        LSs_dxy_[i].initialize(diff_size, 0.0);
        LSs_dx_backward_[i].initialize(diff_size, 0.0);
        LSs_dy_backward_[i].initialize(diff_size, 0.0);
        LSs_dx_forward_[i].initialize(diff_size, 0.0);
        LSs_dy_forward_[i].initialize(diff_size, 0.0);
        for (int j = 0; j < diff_size; j++)
        {
            cache.read(reinterpret_cast<char*>(&LSs_dxx_[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&LSs_dyy_[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&LSs_dxy_[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&LSs_dx_backward_[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&LSs_dy_backward_[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&LSs_dx_forward_[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&LSs_dy_forward_[i].value[j]), sizeof(double));
        }
    }

    cache.close();
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        cout << "[ImmersedGeometry] Geometry cache read." << endl;
}
