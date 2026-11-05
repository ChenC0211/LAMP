#include "Modules/Interface/HermiteGalsReflux.h"
#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"

#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <type_traits>

namespace LAMP::Interface
{
namespace
{
void require(bool ok, HermiteGalsIdentity id, const std::string& message)
{
    if (!ok) throw std::invalid_argument("[Interface/HermiteGalsReflux geometry="
        + std::to_string(id.value()) + "] " + message);
}
void add(HermiteGalsSideFlux& side, const HermiteGalsPhaseFaceFlux& face)
{
    side.totalPositiveRate += face.totalPositiveRate;
    side.totalNegativeRate += face.totalNegativeRate;
    side.liquidPositiveRate += face.liquidPositiveRate;
    side.liquidNegativeRate += face.liquidNegativeRate;
    side.totalRate += face.totalRate;
    side.liquidRate += face.liquidRate;
    side.liquidVolume += face.liquidVolume;
    side.integrationErrorEstimate += face.integrationErrorEstimate;
}
}
struct HermiteGalsRefluxResult::Data
{
    explicit Data(const HermiteGalsPhaseFluxResult& input) : source(input) {}
    HermiteGalsPhaseFluxResult source;
    std::vector<HermiteGalsSideFlux> sides;
    HermiteGalsRefluxReport report;
};
HermiteGalsRefluxResult::HermiteGalsRefluxResult(std::shared_ptr<const Data> data)
    : data_(std::move(data)) {}
const HermiteGalsRefluxResult::Data& HermiteGalsRefluxResult::checked() const
{
    require(static_cast<bool>(data_), {}, "result: moved-from handle");
    return *data_;
}
const HermiteGalsPhaseFluxResult& HermiteGalsRefluxResult::source() const { return checked().source; }
const std::vector<HermiteGalsSideFlux>& HermiteGalsRefluxResult::sides() const { return checked().sides; }
const HermiteGalsRefluxReport& HermiteGalsRefluxResult::report() const { return checked().report; }
void HermiteGalsRefluxResult::validateMesh(const QuadtreeMesh& mesh) const { source().validateMesh(mesh); }


namespace
{
struct FaceOwner {int nativeFace,negativeCell,positiveCell;};
FaceOwner owner(const HermiteGalsPhaseFaceFlux& f){return {f.nativeFace,f.negativeCell,f.positiveCell};}
FaceOwner owner(const HermiteGalsIntegratedFaceTransfer& f)
{return {f.face().nativeRow,f.face().negativeCell,f.face().positiveCell};}
HermiteGalsIdentity refluxIdentity(const HermiteGalsPhaseFluxResult& f){return f.geometryIdentity();}
HermiteGalsIdentity refluxIdentity(const HermiteGalsIntegratedInventoryResult& f){return f.interval().previousGeometryIdentity();}
bool finiteSide(const HermiteGalsSideFlux& s)
{
    for(long double value:{s.totalPositiveRate,s.totalNegativeRate,s.liquidPositiveRate,s.liquidNegativeRate,
        s.totalRate,s.liquidRate,s.liquidVolume,s.integrationErrorEstimate})if(!std::isfinite(value))return false;
    return true;
}
bool finiteSide(const HermiteGalsIntegratedSideTransfer& s)
{
    for(long double value:{s.totalPositiveVolume,s.totalNegativeVolume,
        s.liquidPositiveVolume,s.liquidNegativeVolume,s.totalVolume,s.liquidVolume})
        if(!std::isfinite(value))return false;
    return true;
}
void add(HermiteGalsIntegratedSideTransfer& side,const HermiteGalsIntegratedFaceTransfer& face)
{
    const auto& f=face.integratedVolumes();
    side.totalPositiveVolume+=f.integratedTotalPositiveVolume;
    side.totalNegativeVolume+=f.integratedTotalNegativeVolume;
    side.liquidPositiveVolume+=f.integratedLiquidPositiveVolume;
    side.liquidNegativeVolume+=f.integratedLiquidNegativeVolume;
    side.totalVolume+=face.signedTotalVolume();side.liquidVolume+=face.signedLiquidVolume();
}
std::array<long double,3> refluxValues(const HermiteGalsSideFlux& s){return {{s.totalRate,s.liquidRate,s.liquidVolume}};}
std::array<long double,2> refluxValues(const HermiteGalsIntegratedSideTransfer& s){return {{s.totalVolume,s.liquidVolume}};}
void recordResidual(HermiteGalsRefluxReport& r,const std::array<long double,3>& v)
{
    r.maxTotalRateResidual=std::max(r.maxTotalRateResidual,v[0]);
    r.maxLiquidRateResidual=std::max(r.maxLiquidRateResidual,v[1]);
    r.maxLiquidVolumeResidual=std::max(r.maxLiquidVolumeResidual,v[2]);
}
void recordResidual(HermiteGalsIntegratedRefluxReport& r,const std::array<long double,2>& v)
{
    r.maxTotalVolumeResidual=std::max(r.maxTotalVolumeResidual,v[0]);
    r.maxLiquidVolumeResidual=std::max(r.maxLiquidVolumeResidual,v[1]);
}
template<class Flux,class Data>void assembleReflux(const QuadtreeMesh& mesh,const Flux& flux,Data& result)
{
    flux.validateMesh(mesh);
    const auto id = refluxIdentity(flux);
    const auto& faces = flux.faces();
    const auto& volumes = flux.phaseVolumes();
    std::vector<int> cellRow(mesh.cellIndexNum(), -1);
    for (std::size_t k = 0; k < volumes.size(); ++k) cellRow[volumes[k].cell] = static_cast<int>(k);
    FaceSegmentGraph graph;
    graph.build(mesh, [&](const QuadtreeMesh&, const Cell& cell) { return cellRow[cell.index()] >= 0; });
    const int nv = mesh.vEdgeIndexNum();
    std::vector<int> faceRow(nv + mesh.hEdgeIndexNum(), -1);
    for (std::size_t k = 0; k < faces.size(); ++k) faceRow[owner(faces[k]).nativeFace] = static_cast<int>(k);
    std::vector<int> incidence(faces.size(), 0);
    using Side=typename std::decay_t<decltype(result.sides)>::value_type;
    result.sides.reserve(4 * volumes.size());
    for (const auto& volume : volumes)
    {
        const int cell = volume.cell;
        std::array<Side, 4> sides;
        for (int s = 0; s < 4; ++s)
        {
            sides[s].cell = cell;
            sides[s].side = static_cast<FaceSide>(s);
            sides[s].outwardSign = s < 2 ? -1 : 1;
        }
        const auto range = graph.segmentRangeForCellIndex(cell);
        for (auto k = range.first; k < range.second; ++k)
        {
            const auto& segment = graph.segments()[k];
            const auto& edge = segment.isHorizontal ? mesh.hEdges(segment.edgeLabel) : mesh.vEdges(segment.edgeLabel);
            const int row = faceRow[edge.index() + (segment.isHorizontal ? nv : 0)];
            require(row >= 0, id, "assembly: missing final native face");
            const auto& face = faces[row];
            auto& side = sides[static_cast<int>(segment.side)];
            require((side.outwardSign > 0 ? owner(face).negativeCell : owner(face).positiveCell) == cell,
                id, "assembly: native face orientation/owner mismatch");
            require(side.subfaceCount < 2, id, "assembly: side exceeds 2:1 subface count");
            side.subfaces[side.subfaceCount++] = static_cast<std::size_t>(row);
            ++incidence[row];
            add(side, face);
        }
        for (const auto& side : sides)
        {
            require(side.subfaceCount > 0, id, "assembly: incomplete leaf side");
            require(finiteSide(side),id,"assembly: unrepresentable side sum");
            result.sides.push_back(side);
        }
    }
    for (std::size_t k = 0; k < faces.size(); ++k)
        require(incidence[k] == int(owner(faces[k]).negativeCell >= 0) + int(owner(faces[k]).positiveCell >= 0),
            id, "assembly: missing/duplicate native face incidence");

    // Structural gates precede diagnostics. Every selected fine neighbor must
    // point back to the exact SAME actual face, with the opposite orientation.
    for (const auto& side : result.sides) if (side.subfaceCount == 2)
    {
        const auto& coarse = mesh.cells(mesh.cellIndexToLabel[side.cell]);
        for (const auto row : side.subfaces)
        {
            const auto& face = faces[row];
            const bool horizontal = owner(face).nativeFace >= nv;
            const auto& edge = horizontal ? mesh.hEdges(mesh.hEdgeIndexToLabel[owner(face).nativeFace - nv])
                : mesh.vEdges(mesh.vEdgeIndexToLabel[owner(face).nativeFace]);
            const int fineLabel = side.outwardSign > 0 ? edge.masterCell : edge.slaveCell;
            require(fineLabel >= 0 && mesh.cells(fineLabel).level() == coarse.level() + 1,
                id, "assembly: split side requires a 2:1 fine neighbor");
            const int fine = side.outwardSign > 0 ? owner(face).positiveCell : owner(face).negativeCell;
            if (fine < 0) continue; // selected-domain exterior still uses actual subfaces
            const auto& other = result.sides[4 * static_cast<std::size_t>(cellRow[fine])
                + (static_cast<int>(side.side) + 2) % 4];
            require(other.subfaceCount == 1 && other.subfaces[0] == row
                && other.outwardSign == -side.outwardSign, id, "assembly: inconsistent coarse/fine pairing");
        }
    }
    auto& report = result.report;
    report.cellSides = result.sides.size();
    for (const auto& side : result.sides) if (side.subfaceCount == 2)
    {
        ++report.coarseSides;
        decltype(refluxValues(side)) total{};
        std::size_t pairs = 0;
        for (const auto row : side.subfaces)
        {
            const auto& face = faces[row];
            const int fine = side.outwardSign > 0 ? owner(face).positiveCell : owner(face).negativeCell;
            if (fine < 0) continue;
            const auto& other = result.sides[4 * static_cast<std::size_t>(cellRow[fine])
                + (static_cast<int>(side.side) + 2) % 4];
            const auto values=refluxValues(other);
            for(std::size_t k=0;k<total.size();++k)total[k]+=other.outwardSign*values[k];
            ++pairs;
        }
        report.coarseFinePairs += pairs;
        if (pairs != 2) continue;
        ++report.interiorCoarseSides;
        auto residual=refluxValues(side);
        for(std::size_t k=0;k<total.size();++k)residual[k]=std::abs(side.outwardSign*residual[k]+total[k]);
        recordResidual(report,residual);
    }

}
}
HermiteGalsRefluxResult assembleHermiteGalsReflux(const QuadtreeMesh& mesh,const HermiteGalsPhaseFluxResult& flux)
{
    auto result=std::make_shared<HermiteGalsRefluxResult::Data>(flux);
    assembleReflux(mesh,flux,*result);
    return HermiteGalsRefluxResult(std::move(result));
}
struct HermiteGalsIntegratedRefluxResult::Data
{
    explicit Data(const HermiteGalsIntegratedInventoryResult& input):source(input){}
    HermiteGalsIntegratedInventoryResult source;
    std::vector<HermiteGalsIntegratedSideTransfer> sides;
    HermiteGalsIntegratedRefluxReport report;
};
HermiteGalsIntegratedRefluxResult::HermiteGalsIntegratedRefluxResult(std::shared_ptr<const Data> data):data_(std::move(data)){}
const HermiteGalsIntegratedRefluxResult::Data& HermiteGalsIntegratedRefluxResult::checked() const
{require(bool(data_),{},"integrated reflux: moved-from handle");return *data_;}
const HermiteGalsIntegratedInventoryResult& HermiteGalsIntegratedRefluxResult::source() const {return checked().source;}
const std::vector<HermiteGalsIntegratedSideTransfer>& HermiteGalsIntegratedRefluxResult::sides() const {return checked().sides;}
const HermiteGalsIntegratedRefluxReport& HermiteGalsIntegratedRefluxResult::report() const {return checked().report;}
void HermiteGalsIntegratedRefluxResult::validateMesh(const QuadtreeMesh& mesh) const {source().validateMesh(mesh);}
long double HermiteGalsIntegratedRefluxResult::averageTotalFlux(std::size_t row) const
{
    const auto value=sides().at(row).totalVolume/source().interval().dt();
    require(std::isfinite(value),refluxIdentity(source()),"integrated reflux: unrepresentable derived average rate");
    return value;
}
long double HermiteGalsIntegratedRefluxResult::averageLiquidFlux(std::size_t row) const
{
    const auto value=sides().at(row).liquidVolume/source().interval().dt();
    require(std::isfinite(value),refluxIdentity(source()),"integrated reflux: unrepresentable derived average rate");
    return value;
}
HermiteGalsIntegratedRefluxResult assembleHermiteGalsReflux(const QuadtreeMesh& mesh,const HermiteGalsIntegratedInventoryResult& flux)
{
    auto result=std::make_shared<HermiteGalsIntegratedRefluxResult::Data>(flux);
    assembleReflux(mesh,flux,*result);
    return HermiteGalsIntegratedRefluxResult(std::move(result));
}
}
