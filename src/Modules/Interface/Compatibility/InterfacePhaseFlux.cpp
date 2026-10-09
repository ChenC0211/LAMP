#include "Modules/Interface/Compatibility/InterfacePhaseFlux.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Geometry/AxisymmetricBoundaryIntegral.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <map>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>

namespace InterfacePhaseFlux
{
namespace
{
void require(bool valid, const char* reason)
{
    if (!valid) throw std::runtime_error(std::string("[InterfacePhaseFlux] ") + reason);
}
}

static FluxBudget phaseFluxBudgetImpl(const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport,
    const LevelSetHermiteGeometry::Snapshot& snapshot, const InterfaceNormalFluxView& velocity,
    const InterfacePhaseGeometry::MappedFacePartitions* mapped)
{
    const auto& domain = transport.geometry();
    InterfacePhaseGeometry::validateLayout(mesh, domain, snapshot);
    require(mapped || !snapshot.usesNativeLeaves(), "native phase flux requires retained mapped face partitions");
    FluxBudget result;
    const long double nan = std::numeric_limits<long double>::quiet_NaN();
    result.positiveFluxByNativeFace.assign(domain.faces(mesh).size(), nan);
    result.outwardFluxByCell.assign(mesh.cellIndexNum(), nan);
    for (int k : domain.cellIndexByRow()) result.outwardFluxByCell[k] = 0.0L;
    for (const auto& face : domain.faces(mesh))
    {
        if (face.negativeCell < 0 && face.positiveCell < 0) continue;
        try
        {
            const auto chart = velocity.faceChart(face.nativeRow);
            const auto& edge = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
            auto a = mesh.verts(edge.edgeToVert[0]).coord(), b = mesh.verts(edge.edgeToVert[1]).coord();
            if (face.horizontal ? a.x > b.x : a.y > b.y) std::swap(a, b);
            const auto same = [](Vector2D p, Vector2D q) { return p.x == q.x && p.y == q.y; };
            require(same(chart.begin, a) && same(chart.end, b) && same(chart.midpoint, face.center),
                "velocity face chart differs from the actual fitted mesh");
            long double flux = 0.0L;
            for (int half = 0; half < 2; ++half)
            {
                const auto begin = half == 0 ? chart.begin : chart.midpoint;
                const auto end = half == 0 ? chart.midpoint : chart.end;
                const bool radialAxis = domain.axisymmetric() && !face.horizontal && chart.begin.x == 0.0
                    && chart.midpoint.x == 0.0 && chart.end.x == 0.0;
                if (transport.boundary(face.nativeRow) == InterfaceTransportDomain::Boundary::StationaryWall || radialAxis)
                {
                    require(velocity.integratedNormalFlux(face.nativeRow, half, 0.0, 1.0) == 0.0L,
                        "wall/axis half-face has nonzero normal flux");
                    continue;
                }
                if (mapped)
                {
                    const auto& query = mapped->halfFace(face.nativeRow,half);
                    flux += velocity.integratedNormalFlux(face.nativeRow,half,query);
                    result.partitionedHalfFaces += query.gridPieceCount() > 1 || query.ranges().size() > 1;
                    continue;
                }
                const int sign = snapshot.certifiedSegmentSign(begin, end);
                if (sign != 0)
                {
                    if (sign < 0) flux += velocity.integratedNormalFlux(face.nativeRow, half, 0.0, 1.0);
                    continue;
                }
                const bool useX = std::abs(end.x-begin.x) >= std::abs(end.y-begin.y);
                const double origin = useX ? begin.x : begin.y;
                const double length = (useX ? end.x : end.y)-origin;
                require(std::isfinite(length) && length > 0.0, "unresolved oriented half-face chart");
                const auto parameter = [&](Vector2D p) {
                    if (same(p, begin)) return 0.0;
                    if (same(p, end)) return 1.0;
                    const double t = ((useX ? p.x : p.y)-origin)/length;
                    require(std::isfinite(t) && t >= 0.0 && t <= 1.0, "phase anchor left its native half-face");
                    return t;
                };
                const auto pieces = snapshot.partitionSegmentParametric(begin, end,
                    LevelSetHermiteGeometry::SegmentPartitionPolicy::CertifyUnresolvedGridPieces);
                ++result.partitionedHalfFaces;
                for (const auto& piece : pieces)
                {
                    if (!piece.negative) continue;
                    const double t0 = parameter(piece.anchorBegin), t1 = parameter(piece.anchorEnd);
                    require(t0 != t1, "distinct phase anchors collapsed in the native half-face chart");
                    // Compose the retained root chart algebraically inside
                    // the integral. Forming two global/half-face root points
                    // first would discard a real sub-ULP wet interval.
                    flux += velocity.integratedNormalFlux(face.nativeRow, half,
                        t0, t1, piece.tBegin, piece.tEnd);
                }
            }
            require(std::isfinite(flux), "nonfinite negative-phase native flux");
            result.positiveFluxByNativeFace[face.nativeRow] = flux;
            if (face.negativeCell >= 0) result.outwardFluxByCell[face.negativeCell] += flux;
            if (face.positiveCell >= 0) result.outwardFluxByCell[face.positiveCell] -= flux;
            result.boundaryOutward += static_cast<long double>(face.outwardSign)*flux;
            result.absoluteFaceFlux += std::abs(flux);
            ++result.activeFaces;
        }
        catch (const std::exception& error)
        {
            std::ostringstream message;
            message << std::setprecision(17) << "[InterfacePhaseFlux] native face=" << face.nativeRow
                << ", center=(" << face.center.x << ',' << face.center.y << "): " << error.what();
            throw std::runtime_error(message.str());
        }
    }
    for (int k : domain.cellIndexByRow())
    {
        require(std::isfinite(result.outwardFluxByCell[k]), "nonfinite phase flux cell budget");
        result.summedCellOutward += result.outwardFluxByCell[k];
    }
    result.assemblyClosure = result.summedCellOutward-result.boundaryOutward;
    require(std::isfinite(result.boundaryOutward) && std::isfinite(result.summedCellOutward)
        && std::isfinite(result.absoluteFaceFlux) && std::isfinite(result.assemblyClosure),
        "nonfinite assembled phase flux budget");
    return result;
}

template<bool Quartic>
static auto phaseFluxSensitivityImpl(const QuadtreeMesh &mesh, const InterfaceTransportDomain& transport,
                                     const InterfacePhaseGeometry::MappedFacePartitions &partitions,
                                     const InterfaceNormalFluxView &velocity,
                                     Quantum::Geometry::CornerSectorVariation sectors)
{
    const auto& domain = transport.geometry();
    partitions.validate(mesh, domain);
    require(partitions.source().usesNativeLeaves(), "phase flux sensitivity requires native leaves");
    require((partitions.source().reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic)==Quartic,
        "phase flux sensitivity port does not match scalar data layout");
    std::conditional_t<Quartic,QuarticFluxSensitivity,FluxSensitivity> result;
    result.positiveFluxByNativeFace.resize(domain.faces(mesh).size());
    result.outwardFluxByCell.resize(mesh.cellIndexNum());
    using Key=std::conditional_t<Quartic,std::tuple<int,int,int>,int>;
    using Values=std::array<long double,Quartic?1:4>;
    using Columns = std::map<Key,Values>;
    std::vector<Columns> cellColumns(mesh.cellIndexNum());
    const auto accumulate = [](Columns &target, const Key& key, const Values &values,
                               long double factor) {
        auto &column = target[key];
        for (std::size_t k = 0; k < values.size(); ++k)
        {
            const long double term = factor * values[k];
            require(std::isfinite(term) && (factor == 0 || values[k] == 0 || term != 0),
                    "nonfinite or underflowed phase-flux sensitivity");
            column[k] += term;
            require(std::isfinite(column[k]), "nonfinite accumulated phase-flux sensitivity");
        }
    };
    const auto finish = [](const Columns &input, auto &output) {
        for (const auto &column : input)
            if (std::any_of(column.second.begin(), column.second.end(),
                            [](long double value) { return value != 0; }))
            {
                if constexpr(Quartic)
                    output.push_back({static_cast<Quantum::Geometry::ActiveQuarticHermiteSensitivity::Kind>(std::get<0>(column.first)),
                        std::get<1>(column.first),std::get<2>(column.first),column.second[0]});
                else output.push_back({column.first, column.second});
            }
    };
    const auto accumulateSensitivity=[&](Columns& target,const auto& column,long double factor) {
        if constexpr(Quartic)
            accumulate(target,{static_cast<int>(column.kind),column.index,column.component},{{column.coefficient}},factor);
        else accumulate(target,column.vertexIndex,column.coefficients,factor);
    };
    for (const auto &face : domain.faces(mesh))
    {
        if (face.negativeCell < 0 && face.positiveCell < 0)
            continue;
        const auto chart = velocity.faceChart(face.nativeRow);
        const bool axis = domain.axisymmetric() && !face.horizontal && chart.begin.x == 0 && chart.midpoint.x == 0 && chart.end.x == 0;
        Columns columns;
        for (int half = 0; half < 2; ++half)
        {
            const auto &query = partitions.halfFace(face.nativeRow, half);
            const auto begin = half == 0 ? chart.begin : chart.midpoint,
                       end = half == 0 ? chart.midpoint : chart.end;
            require(query.begin().x == begin.x && query.begin().y == begin.y && query.end().x == end.x &&
                        query.end().y == end.y,
                    "phase variation and velocity use different actual half-face charts");
            if (axis || transport.boundary(face.nativeRow) == InterfaceTransportDomain::Boundary::StationaryWall)
            {
                require(velocity.integratedNormalFlux(face.nativeRow, half, 0, 1) == 0,
                        "wall/axis phase variation requires identically zero normal flux");
                continue;
            }
            const auto variations=[&] {
                if constexpr(Quartic) return query.nativeQuarticNegativePhaseVariation(sectors);
                else return query.nativeNegativePhaseVariation(sectors);
            }();
            for (const auto &root : variations)
            {
                const auto density =
                    velocity.normalFluxDensity(face.nativeRow, half, root.fromEnd, root.offset);
                for (const auto &column : root.independent)
                    accumulateSensitivity(columns,column,density);
                ++result.roots;
                result.endpointRoots += root.endpoint;
            }
        }
        auto &faceColumns = result.positiveFluxByNativeFace[face.nativeRow];
        finish(columns, faceColumns);
        for (const auto &column : faceColumns)
        {
            if (face.negativeCell >= 0)
                accumulateSensitivity(cellColumns[face.negativeCell],column,1);
            if (face.positiveCell >= 0)
                accumulateSensitivity(cellColumns[face.positiveCell],column,-1);
        }
    }
    for (int cell : domain.cellIndexByRow())
        finish(cellColumns[cell], result.outwardFluxByCell[cell]);
    return result;
}

FluxSensitivity phaseFluxSensitivity(const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions, const InterfaceNormalFluxView& velocity,
    Quantum::Geometry::CornerSectorVariation sectors)
{
    return phaseFluxSensitivityImpl<false>(mesh,transport,partitions,velocity,sectors);
}

QuarticFluxSensitivity quarticPhaseFluxSensitivity(const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions, const InterfaceNormalFluxView& velocity,
    Quantum::Geometry::CornerSectorVariation sectors)
{
    return phaseFluxSensitivityImpl<true>(mesh,transport,partitions,velocity,sectors);
}

FluxBudget phaseFluxBudget(const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport,
    const LevelSetHermiteGeometry::Snapshot& snapshot, const InterfaceNormalFluxView& velocity)
{
    return phaseFluxBudgetImpl(mesh,transport,snapshot,velocity,nullptr);
}
FluxBudget phaseFluxBudget(const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions, const InterfaceNormalFluxView& velocity)
{
    const auto& domain = transport.geometry();
    partitions.validate(mesh,domain);
    return phaseFluxBudgetImpl(mesh,transport,partitions.source(),velocity,&partitions);
}
}
