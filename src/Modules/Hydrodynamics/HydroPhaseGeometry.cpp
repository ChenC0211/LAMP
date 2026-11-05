#include "Modules/Hydrodynamics/HydroPhaseGeometry.h"

#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/PlanarFaceTimeIntegral.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Geometry/Detail/ExactRationalRound.h"
#include "Core/Geometry/UnitPolynomialRoots.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/HydroNormalRatePolynomial.h"
#include "Modules/Interface/HermiteGalsFirstJetStep.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <optional>
#include <unordered_map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace LAMP::Hydrodynamics
{
namespace
{
void require(bool ok,const char* reason)
{
    if(!ok)throw std::invalid_argument(std::string("[HydroPhaseGeometry] ")+reason);
}
void require(bool ok, const std::string& reason)
{
    if (!ok) throw std::invalid_argument("[HydroPhaseGeometry] "+reason);
}
bool sameTime(Interface::HermiteGalsTimeLayer a, Interface::HermiteGalsTimeLayer b)
{
    a.time(); b.time(); // Validate every field before comparison, including NaN.
    return a.step==b.step && a.begin==b.begin && a.end==b.end && a.position==b.position;
}
HydroPhaseMeasure measure(Vector2D a, Vector2D b, const std::vector<LevelSetMappedSegment::Range>& ranges, bool radial)
{
    const long double dx=static_cast<long double>(b.x)-a.x, dy=static_cast<long double>(b.y)-a.y;
    const long double length=std::hypot(dx,dy);
    require(std::isfinite(length) && length>0,"invalid retained half-face length");
    HydroPhaseMeasure result;
    for (const auto& range : ranges) {
        require(std::isfinite(range.offset) && range.offset>=0 &&
            std::isfinite(range.length) && range.length>0,"invalid retained phase range");
        // Integrate r(t) from the independently retained offset and width.
        // Never construct two global endpoints or subtract large primitives.
        const long double origin=range.fromEnd ? b.x : a.x, slope=range.fromEnd ? -dx : dx;
        const long double radius=origin+slope*range.offset;
        const long double meanRadius=radius+slope*(range.length/2);
        require(!radial || (std::isfinite(meanRadius) && meanRadius>=0),"invalid phase radius");
        const long double area=length*range.length*(radial ? CartesianMeasure::twoPi*meanRadius : 1.L);
        require(std::isfinite(area) && area>=0,"nonfinite phase face measure");
        require(area>0 || (radial && a.x==0 && b.x==0),"positive phase measure underflow");
        (range.negative ? result.liquid : result.gas)+=area;
    }
    require(std::isfinite(result.liquid) && std::isfinite(result.gas),"phase face sum overflow");
    return result;
}
HydroPhaseMeasure measure(const LevelSetMappedSegment& segment, bool radial)
{ return measure(segment.begin(),segment.end(),segment.ranges(),radial); }
void add(HydroPhaseMeasure& a,HydroPhaseMeasure b){a.liquid+=b.liquid;a.gas+=b.gas;}
constexpr long double vofClosureTolerance=512.L*std::numeric_limits<double>::epsilon();
struct PositiveVofSum
{
    long double value=0,correction=0;
    void add(long double x)
    {
        require(std::isfinite(x) && x>=0,"invalid positive VOF contribution");
        const auto y=x-correction,t=value+y;
        correction=(t-value)-y;value=t;
        require(std::isfinite(value),"VOF phase sum overflow");
    }
};
long double vofRelativeClosure(HydroPhaseMeasure phases,long double capacity)
{
    require(std::isfinite(capacity) && capacity>0,"invalid VOF rectangle capacity");
    require(std::isfinite(phases.liquid) && std::isfinite(phases.gas) &&
        phases.liquid>=0 && phases.gas>=0,"invalid measured VOF phases");
    const auto closure=std::abs(phases.liquid+phases.gas-capacity)/capacity;
    require(std::isfinite(closure) && closure<=vofClosureTolerance,"VOF phase capacity closure failed");
    return closure;
}
constexpr long double localPlicTolerance=1e-12L;
long double localPlicClosure(HydroPhaseMeasure phases,long double capacity)
{
    require(std::isfinite(capacity)&&capacity>0&&std::isfinite(phases.liquid)&&std::isfinite(phases.gas)
        &&phases.liquid>=0&&phases.gas>=0,"invalid actual local PLIC phases");
    const auto defect=std::abs(phases.liquid+phases.gas-capacity)/capacity;
    require(std::isfinite(defect)&&defect<=localPlicTolerance,"local PLIC phase capacity closure failed");return defect;
}
HydroPhaseMeasure localPlicBox(const Interface::GeometricVofPlicCell& cell,
    const LAMP::Geometry::PlanarCellCut::Rectangle& box,bool rz,bool fromPlane=false)
{
    namespace Cut=LAMP::Geometry::PlanarCellCut;
    const auto& rectangle=cell.rectangle;
    require(box.lower.x>=rectangle.lower.x&&box.lower.y>=rectangle.lower.y
        &&box.upper.x<=rectangle.upper.x&&box.upper.y<=rectangle.upper.y,
        "local PLIC box leaves its original native rectangle");
    const auto capacity=Cut::capacity(box,rz);
    if(!fromPlane&&box.lower.x==rectangle.lower.x&&box.lower.y==rectangle.lower.y
        &&box.upper.x==rectangle.upper.x&&box.upper.y==rectangle.upper.y){
        const HydroPhaseMeasure measured{cell.measured[0],cell.measured[1]};
        localPlicClosure(measured,capacity);return measured;
    }
    // Actual four corners in the original cell-centre frame. Clip by the
    // SAME published plane independently for liquid and gas; no stock mapping.
    const Cut::Point origin{static_cast<long double>(rectangle.lower.x)
        +(static_cast<long double>(rectangle.upper.x)-rectangle.lower.x)/2,
        static_cast<long double>(rectangle.lower.y)+(static_cast<long double>(rectangle.upper.y)-rectangle.lower.y)/2};
    const auto point=[&](double x,double y){return Cut::Point{static_cast<long double>(x)-origin.x,static_cast<long double>(y)-origin.y};};
    const std::vector<Cut::Point> polygon{point(box.lower.x,box.lower.y),point(box.upper.x,box.lower.y),
        point(box.upper.x,box.upper.y),point(box.lower.x,box.upper.y)};
    const auto liquid=Cut::overlapInFrame(rectangle,polygon,origin,cell.normal,cell.offset,rz);
    const auto gas=Cut::overlapInFrame(rectangle,polygon,origin,Vector2D(-cell.normal.x,-cell.normal.y),-cell.offset,rz);
    const HydroPhaseMeasure measured{liquid.liquid,gas.liquid};localPlicClosure(measured,capacity);return measured;
}
HydroPhaseMeasure localPlicTrace(const Interface::GeometricVofPlicCell& cell,
    LAMP::Geometry::PlanarCellCut::Point a,LAMP::Geometry::PlanarCellCut::Point b,bool rz)
{
    namespace Cut=LAMP::Geometry::PlanarCellCut;
    const auto& rectangle=cell.rectangle;
    const auto length=std::hypot(static_cast<long double>(cell.normal.x),static_cast<long double>(cell.normal.y));
    require(std::isfinite(length)&&length>0,"invalid original local PLIC trace normal");
    const auto nx=cell.normal.x/length,ny=cell.normal.y/length;
    const Cut::Point center{static_cast<long double>(rectangle.lower.x)+(static_cast<long double>(rectangle.upper.x)-rectangle.lower.x)/2,
        static_cast<long double>(rectangle.lower.y)+(static_cast<long double>(rectangle.upper.y)-rectangle.lower.y)/2};
    const auto value=[&](Cut::Point p){return nx*(p.x-center.x)+ny*(p.y-center.y)-cell.offset;};
    const auto physical=[&](Cut::Point first,Cut::Point last){
        const auto distance=std::hypot(last.x-first.x,last.y-first.y);
        require(std::isfinite(distance)&&distance>=0&&(!rz||(first.x>=0&&last.x>=0)),"invalid local PLIC physical trace");
        const auto area=distance*(rz?CartesianMeasure::twoPi*(first.x/2+last.x/2):1.L);
        require(std::isfinite(area)&&area>=0,"nonfinite local PLIC trace measure");return area;};
    const auto va=value(a),vb=value(b);require(std::isfinite(va)&&std::isfinite(vb),"nonfinite local PLIC plane trace");
    if(va==0&&vb==0){
        // The retained plane lies on this native boundary. Its actual inside
        // rectangle side supplies the trace; there is no global shared Gamma.
        require(cell.offset!=0,"unresolved local PLIC whole-face side");
        const auto area=physical(a,b);return cell.offset>0?HydroPhaseMeasure{area,0}:HydroPhaseMeasure{0,area};
    }
    const auto capacity=physical(a,b);
    if(va<=0&&vb<=0)return {capacity,0};
    if(va>=0&&vb>=0)return {0,capacity};
    if(capacity==0)return {0,0}; // A true zero-area radial-axis face.
    // Integrate the same affine plane in the segment parameter. A positive
    // tiny interval can be represented even when its WORLD intersection
    // rounds to an endpoint. Keep each phase's own fraction, never 1-small.
    const auto scale=std::max(std::abs(va),std::abs(vb));
    const auto firstWeight=std::abs(va)/scale,lastWeight=std::abs(vb)/scale;
    const auto denominator=firstWeight+lastWeight;
    const auto distance=std::hypot(b.x-a.x,b.y-a.y),dx=b.x-a.x;
    const auto portion=[&](long double fraction,bool fromFirst){
        require(std::isfinite(fraction)&&fraction>0&&fraction<=1,"unrepresentable positive local PLIC trace parameter");
        const auto radius=fromFirst?std::fma(fraction/2,dx,a.x):std::fma(-fraction/2,dx,b.x);
        const auto area=distance*fraction*(rz?CartesianMeasure::twoPi*radius:1.L);
        require(std::isfinite(area)&&area>0&&(!rz||radius>0),"unrepresentable positive local PLIC trace measure");return area;};
    const auto first=portion(firstWeight/denominator,true),last=portion(lastWeight/denominator,false);
    return va<0?HydroPhaseMeasure{first,last}:HydroPhaseMeasure{last,first};
}
bool vofPolygonInsideBox(const std::vector<LAMP::Geometry::PlanarCellCut::AnchoredPoint>& polygon,
    LAMP::Geometry::PlanarCellCut::Point origin,const LAMP::Geometry::PlanarCellCut::Rectangle& box)
{
    return LAMP::Geometry::PlanarCellCut::anchoredPolygonInsideBox(polygon,origin,box);
}
}
struct HydroPhaseGeometry::Data
{
    Data(HydroPhysicalDivergence m, Interface::HermiteGalsPublication p)
        : metrics(std::move(m)), publication(std::move(p)) {}
    Data(HydroPhysicalDivergence m, Interface::GeometricVofPublication p)
        : metrics(std::move(m)), vof(std::move(p)) {}
    HydroPhysicalDivergence metrics;
    std::optional<Interface::HermiteGalsPublication> publication;
    std::optional<Interface::GeometricVofPublication> vof;
    std::vector<HydroPhaseMeasure> cells, faces;
    std::vector<HydroPhaseFaceTrace> traces;
    std::vector<double> alpha;
    HydroPhaseGeometryReport report;
};
HydroPhaseGeometry::HydroPhaseGeometry(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
HydroPhaseGeometry HydroPhaseGeometry::build(const QuadtreeMesh& mesh, HydroPhysicalDivergence metrics,
    Interface::HermiteGalsPublication publication)
{
    metrics.layout().validateMesh(mesh); publication.validateMesh(mesh);
    const auto& primary=publication.state().primary();
    require(primary.transport==Interface::HermiteGalsTransport::FirstJet && primary.phaseVolumes.empty(),
        "requires FirstJet publication without a conservative inventory");
    require(metrics.layout().domain().axisymmetric()==publication.domain().axisymmetric(),
        "different XY/RZ measure");
    const auto view=publication.state().geometry();
    require(publication.volumes().sourceIdentity()==view.identity() &&
        sameTime(publication.volumes().timeLayer(),view.timeLayer()),"phase volume provenance mismatch");
    require(&publication.facePartitions().source()==&view.snapshot(),"different phase partition scalar");
    auto data=std::make_shared<Data>(std::move(metrics),std::move(publication));
    const auto& layout=data->metrics.layout();
    std::vector<const Interface::HermiteGalsPhaseVolumeCell*> cellByIndex(
        data->publication->domain().cellRowByIndex().size(),nullptr);
    for (const auto& cell : data->publication->volumes().cells()) {
        require(cell.cell>=0 && static_cast<std::size_t>(cell.cell)<cellByIndex.size() &&
            !cellByIndex[cell.cell],"duplicate/invalid publication volume cell");
        cellByIndex[cell.cell]=&cell;
    }
    for (int p=0;p<layout.pressureCount();++p) {
        const int native=layout.nativeCell(p);
        require(static_cast<std::size_t>(native)<cellByIndex.size() && cellByIndex[native],
            "missing selected-cell phase geometry");
        const auto& cell=*cellByIndex[native];
        const double capacity=data->metrics.cellVolumes()[p];
        require(cell.capacity==capacity,"different true cell capacity");
        // The publication's certified pure phase uses exact true capacity;
        // rawLiquid is its independently checked Green-ledger diagnostic.
        // It can differ by rounding and is not the pure-phase authority.
        require(cell.certifiedPhase>=-1 && cell.certifiedPhase<=1,"invalid pure-phase certificate");
        const long double liquid=cell.certifiedPhase<0 ? capacity : cell.certifiedPhase>0 ? 0 : cell.rawLiquid;
        require(std::isfinite(liquid) && liquid>=0 && liquid<=capacity,
            "invalid geometric phase volume at cell="+std::to_string(native));
        data->cells.push_back({liquid,capacity-liquid});
        const double alpha=static_cast<double>(liquid/capacity);
        require(std::isfinite(alpha) && alpha>=0 && alpha<=1,"invalid geometric liquid fraction");
        data->alpha.push_back(alpha);
        data->report.cutCells+=cell.certifiedPhase==0;
    }
    for (int u=0;u<layout.velocityCount();++u) {
        const int native=layout.nativeFace(u);
        const auto& partitions=data->publication->facePartitions();
        const auto face=layout.domain().face(mesh,native);
        HydroPhaseMeasure area;HydroPhaseFaceTrace trace;
        trace.negativeCell=face.negativeCell;trace.positiveCell=face.positiveCell;
        for (int half=0;half<2;++half) {
            const auto& segment=partitions.traceHalfFace(native,half);
            require(segment.usesSource(view.snapshot()),"different half-face geometry");
            if(segment.coincidentRanges().empty()) {
                const auto m=measure(segment.ordinary(),layout.domain().axisymmetric());add(area,m);
                if(trace.minusAvailable())add(trace.minus,m);if(trace.plusAvailable())add(trace.plus,m);
            } else {
                ++trace.coincidentHalfFaces;
                const auto ordinary=measure(segment.begin(),segment.end(),segment.ordinaryRanges(),layout.domain().axisymmetric());
                HydroPhaseMeasure minus=ordinary,plus=ordinary;
                for(const auto& zero:segment.coincidentRanges()) {
                    if(trace.minusAvailable()) {
                        require(zero.minus.available && zero.minus.nativeCellIndex==trace.negativeCell,"missing/foreign selected minus trace");
                        add(minus,measure(segment.begin(),segment.end(),zero.minus.ranges,layout.domain().axisymmetric()));
                    }
                    if(trace.plusAvailable()) {
                        require(zero.plus.available && zero.plus.nativeCellIndex==trace.positiveCell,"missing/foreign selected plus trace");
                        add(plus,measure(segment.begin(),segment.end(),zero.plus.ranges,layout.domain().axisymmetric()));
                    }
                }
                if(trace.minusAvailable())add(trace.minus,minus);if(trace.plusAvailable())add(trace.plus,plus);
                if(trace.minusAvailable() && trace.plusAvailable())add(area,{(minus.liquid+plus.liquid)*.5L,(minus.gas+plus.gas)*.5L});
                else add(area,trace.minusAvailable()?minus:plus);
            }
        }
        if(trace.coincidentHalfFaces) {
            ++data->report.coincidentFaces;
            trace.kind=trace.minusAvailable() && trace.plusAvailable()?HydroPhaseFaceTraceKind::InteriorCoincident:HydroPhaseFaceTraceKind::BoundaryCoincident;
            data->report.boundaryCoincidentFaces+=trace.kind==HydroPhaseFaceTraceKind::BoundaryCoincident;
        }
        const double capacity=data->metrics.faceMeasures()[u];
        const long double defect=capacity==0 ? std::abs(area.liquid+area.gas) :
            std::abs(area.liquid+area.gas-capacity)/capacity;
        require(std::isfinite(defect) && defect<=128*std::numeric_limits<double>::epsilon(),
            "native face phase measure closure failed at face="+std::to_string(native));
        require(area.liquid>=0 && area.gas>=0,"negative phase face measure");
        data->report.maxRelativeFaceClosure=std::max(data->report.maxRelativeFaceClosure,defect);
        data->report.splitFaces+=area.liquid>0 && area.gas>0;
        data->faces.push_back(area);
        data->traces.push_back(trace);
    }
    return HydroPhaseGeometry(std::move(data));
}
HydroPhaseGeometry HydroPhaseGeometry::build(const QuadtreeMesh& mesh, HydroPhysicalDivergence metrics,
    Interface::GeometricVofPublication publication)
{
    metrics.layout().validateMesh(mesh);
    publication.validateFor(mesh,publication.identity(),publication.timeLayer());
    require(metrics.layout().domain().axisymmetric()==publication.domain().axisymmetric(),
        "different VOF XY/RZ measure");
    require(metrics.layout().domain().meshIdentity()==publication.domain().meshIdentity(),
        "different retained VOF native mesh binding");
    LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
    auto data=std::make_shared<Data>(std::move(metrics),std::move(publication));
    HydroPhaseGeometry result(data);
    const auto& source=*data->vof;const auto& layout=data->metrics.layout();
    const auto& sourceCells=source.cells();
    require(sourceCells.size()==source.alpha().size(),"different VOF alpha/cell support");
    if(source.geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic){
        const auto& local=source.plicCells();require(local.size()==sourceCells.size(),"different local PLIC native support");
        const bool rz=layout.domain().axisymmetric();
        for(int p=0;p<layout.pressureCount();++p){
            const int native=layout.nativeCell(p);require(native>=0&&static_cast<std::size_t>(native)<local.size()
                &&sourceCells[native].nativeCell==native,"missing selected local PLIC cell");
            const auto capacity=data->metrics.cellVolumes()[p];const auto& cell=local[native];
            require(std::abs(sourceCells[native].capacity-capacity)/capacity<=localPlicTolerance,"different true local PLIC cell capacity");
            const HydroPhaseMeasure phases{cell.measured[0],cell.measured[1]};
            data->report.maxRelativeCellClosure=std::max(data->report.maxRelativeCellClosure,localPlicClosure(phases,capacity));
            for(unsigned phase=0;phase<2;++phase){const auto actual=phase==0?phases.liquid:phases.gas;
                const auto defect=std::abs(actual-cell.inventory[phase])/capacity;
                require(std::isfinite(defect)&&defect<=localPlicTolerance,"local PLIC measured/stock capacity-scaled gate failed");
                data->report.maxRelativeInventoryMeasurement=std::max(data->report.maxRelativeInventoryMeasurement,defect);}
            data->cells.push_back(phases);data->alpha.push_back(source.alpha()[native]);
            data->report.cutCells+=cell.inventory[0]>0&&cell.inventory[1]>0;
        }
        for(int u=0;u<layout.velocityCount();++u){
            const int native=layout.nativeFace(u);const auto face=layout.domain().face(mesh,native);
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            const auto first=mesh.verts(edge.edgeToVert[0]).coord(),last=mesh.verts(edge.edgeToVert[1]).coord();
            const LAMP::Geometry::PlanarCellCut::Point a{first.x,first.y},b{last.x,last.y};
            HydroPhaseFaceTrace trace;trace.negativeCell=face.negativeCell;trace.positiveCell=face.positiveCell;
            if(trace.minusAvailable())trace.minus=localPlicTrace(local[trace.negativeCell],a,b,rz);
            if(trace.plusAvailable())trace.plus=localPlicTrace(local[trace.positiveCell],a,b,rz);
            const auto area=trace.minusAvailable()&&trace.plusAvailable()
                ?HydroPhaseMeasure{(trace.minus.liquid+trace.plus.liquid)/2,(trace.minus.gas+trace.plus.gas)/2}
                :(trace.minusAvailable()?trace.minus:trace.plus);
            const auto capacity=data->metrics.faceMeasures()[u];
            const auto closure=capacity==0?std::abs(area.liquid+area.gas):std::abs(area.liquid+area.gas-capacity)/capacity;
            require(std::isfinite(closure)&&closure<=localPlicTolerance,"native local PLIC face phase closure failed");
            data->report.maxRelativeFaceClosure=std::max(data->report.maxRelativeFaceClosure,closure);
            data->report.splitFaces+=area.liquid>0&&area.gas>0;data->faces.push_back(area);data->traces.push_back(trace);
        }
        return result;
    }
    std::vector<HydroVofBox> boxes;boxes.reserve(layout.pressureCount());
    for(int p=0;p<layout.pressureCount();++p) {
        const int native=layout.nativeCell(p);
        require(native>=0 && static_cast<std::size_t>(native)<sourceCells.size() &&
            sourceCells[native].nativeCell==native,"missing selected VOF cell");
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
        boxes.push_back({native,{mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()}});
    }
    data->cells=result.measureVofBoxes(mesh,boxes);
    for(int p=0;p<layout.pressureCount();++p) {
        const int native=layout.nativeCell(p);const auto& cell=sourceCells[native];
        const auto phases=data->cells[p];const auto capacity=data->metrics.cellVolumes()[p];
        require(std::abs(cell.capacity-capacity)/capacity<=vofClosureTolerance,
            "different true VOF cell capacity");
        data->report.maxRelativeCellClosure=std::max(data->report.maxRelativeCellClosure,
            vofRelativeClosure(phases,capacity));
        const auto alpha=source.alpha()[native];data->alpha.push_back(alpha);
        if(alpha>0 && alpha<1) {
            require(phases.liquid>0 && phases.gas>0,"unresolved mixed VOF phase");
            const auto minor=std::min(cell.inventory,cell.capacity-cell.inventory);
            require(std::isfinite(minor) && minor>0,"unresolved original VOF minor inventory");
            const auto defect=std::abs(phases.liquid-cell.inventory)/minor;
            if(!std::isfinite(defect) || !(defect<=1e-12L)) {
                const auto clock=source.timeLayer();std::ostringstream detail;
                detail.precision(std::numeric_limits<long double>::max_digits10);
                detail<<"actual VOF fan differs from original liquid inventory owner="<<source.identity().value()
                    <<" nativeCell="<<native<<" clockStep="<<clock.step<<" clockBegin="<<clock.begin
                    <<" clockEnd="<<clock.end<<" clockPosition="<<static_cast<int>(clock.position)
                    <<" actualLiquid="<<phases.liquid<<" actualGas="<<phases.gas
                    <<" sourceMeasuredLiquid="<<cell.measured<<" inventory="<<cell.inventory
                    <<" capacity="<<cell.capacity<<" minor="<<minor<<" relative="<<defect;
                require(false,detail.str());
            }
            data->report.maxRelativeInventoryMeasurement=std::max(data->report.maxRelativeInventoryMeasurement,defect);
            ++data->report.cutCells;
        } else {
            require(alpha==0 || alpha==1,"invalid original VOF alpha");
            require(alpha==0 ? phases.liquid==0 : phases.gas==0,"nonzero opposite pure VOF phase");
        }
    }
    namespace Cut=LAMP::Geometry::PlanarCellCut;
    const auto& sourceFaces=source.faces();const auto& vertices=source.anchoredVertices();
    const auto& scalar=source.scalarValues();const auto origin=source.frameOrigin();
    const bool rz=layout.domain().axisymmetric();
    for(int u=0;u<layout.velocityCount();++u) {
        const int native=layout.nativeFace(u);
        require(native>=0 && static_cast<std::size_t>(native)<sourceFaces.size(),"missing selected VOF face");
        const auto& retained=sourceFaces[native];const auto face=layout.domain().face(mesh,native);
        require(retained.nativeRow==native && retained.rangeCount>0 && retained.rangeCount<=retained.ranges.size(),
            "invalid retained VOF face/range identity");
        HydroPhaseFaceTrace trace;trace.negativeCell=face.negativeCell;trace.positiveCell=face.positiveCell;
        require(!trace.minusAvailable() || retained.negativeCell==trace.negativeCell,"foreign VOF minus trace");
        require(!trace.plusAvailable() || retained.positiveCell==trace.positiveCell,"foreign VOF plus trace");
        for(auto vertex:retained.vertices)require(vertex<vertices.size() && vertex<scalar.size(),"invalid VOF face vertex");
        const auto a=vertices[retained.vertices[0]],b=vertices[retained.vertices[1]];
        const auto fullDelta=Cut::anchoredDifference(b,a);
        const auto fullLength=std::hypot(fullDelta.x,fullDelta.y);
        require(std::isfinite(fullLength) && fullLength>0,"invalid original VOF face length");
        require(face.horizontal ? (fullDelta.x>0 && fullDelta.y==0) : (fullDelta.y>0 && fullDelta.x==0),
            "native VOF face is not increasing axis-aligned geometry");
        const auto normalCoordinate=[&](Cut::AnchoredPoint point) {
            return face.horizontal ? Cut::AnchoredPoint{{0,point.anchor.y},{0,point.offset.y}}
                : Cut::AnchoredPoint{{point.anchor.x,0},{point.offset.x,0}};
        };
        const Cut::AnchoredPoint physicalAxis{{-origin.x,0},{0,0}};
        const auto radialCoordinate=[](Cut::AnchoredPoint point) {
            return Cut::AnchoredPoint{{point.anchor.x,0},{point.offset.x,0}};
        };
        const auto physicalRadius=[&](Cut::AnchoredPoint point) {
            const auto radius=Cut::anchoredDifference(radialCoordinate(point),physicalAxis).x;
            require(std::isfinite(radius) && radius>=0,"invalid original VOF trace radius");
            return radius;
        };
        std::array<PositiveVofSum,4> sums;auto cursor=retained.vertices[0];
        for(unsigned j=0;j<retained.rangeCount;++j) {
            const auto& range=retained.ranges[j];
            require(range.vertices[0]==cursor && range.vertices[0]<vertices.size() &&
                range.vertices[1]<vertices.size() && range.vertices[0]!=range.vertices[1],
                "invalid original VOF actual range adjacency");
            const auto first=vertices[range.vertices[0]],last=vertices[range.vertices[1]];
            require(Cut::sameAnchoredPoint(normalCoordinate(first),normalCoordinate(a)) &&
                Cut::sameAnchoredPoint(normalCoordinate(last),normalCoordinate(a)),
                "original VOF actual range leaves its native line");
            const auto delta=Cut::anchoredDifference(last,first);
            const auto fromBegin=Cut::anchoredDifference(first,a),toEnd=Cut::anchoredDifference(b,last);
            require(face.horizontal ? (delta.x>0 && delta.y==0 && fromBegin.x>=0 && toEnd.x>=0)
                : (delta.y>0 && delta.x==0 && fromBegin.y>=0 && toEnd.y>=0),
                "invalid original VOF strict actual range order or support");
            const auto length=std::hypot(delta.x,delta.y);
            const auto radius=rz ? (physicalRadius(first)/2+physicalRadius(last)/2) : 1.L;
            const bool axis=rz && Cut::sameAnchoredPoint(radialCoordinate(first),physicalAxis) &&
                Cut::sameAnchoredPoint(radialCoordinate(last),physicalAxis);
            const auto measure=length*(rz ? CartesianMeasure::twoPi*radius : 1.L);
            require(std::isfinite(length) && length>0 && std::isfinite(measure) && measure>=0 &&
                (measure>0 || (axis && face.area==0)),"unrepresentable positive VOF face trace");
            if(trace.minusAvailable())sums[range.minusLiquid?0:1].add(measure);
            if(trace.plusAvailable())sums[range.plusLiquid?2:3].add(measure);
            cursor=range.vertices[1];
        }
        require(cursor==retained.vertices[1],"incomplete original VOF actual face ranges");
        trace.minus={sums[0].value,sums[1].value};trace.plus={sums[2].value,sums[3].value};
        if(retained.coincident) {
            require(scalar[retained.vertices[0]]==0 && scalar[retained.vertices[1]]==0 &&
                retained.negativeCell>=0 && retained.positiveCell>=0 && face.area>0 &&
                retained.rangeCount==1 && retained.ranges[0].minusLiquid!=retained.ranges[0].plusLiquid,
                "invalid original VOF coincident incidence");
            // VOF publishes the actual whole native trace, not FirstJet half
            // partitions. coincidentHalfFaces stays zero (not applicable).
            trace.kind=trace.minusAvailable() && trace.plusAvailable() ?
                HydroPhaseFaceTraceKind::InteriorCoincident:HydroPhaseFaceTraceKind::BoundaryCoincident;
            ++data->report.coincidentFaces;
            data->report.boundaryCoincidentFaces+=trace.kind==HydroPhaseFaceTraceKind::BoundaryCoincident;
        }
        const HydroPhaseMeasure area=trace.minusAvailable() && trace.plusAvailable() ?
            HydroPhaseMeasure{(trace.minus.liquid+trace.plus.liquid)*.5L,(trace.minus.gas+trace.plus.gas)*.5L} :
            trace.minusAvailable()?trace.minus:trace.plus;
        const auto capacity=data->metrics.faceMeasures()[u];
        const auto closure=capacity==0 ? std::abs(area.liquid+area.gas) :
            std::abs(area.liquid+area.gas-capacity)/capacity;
        require(std::isfinite(closure) && closure<=vofClosureTolerance,"native VOF face phase closure failed");
        data->report.maxRelativeFaceClosure=std::max(data->report.maxRelativeFaceClosure,closure);
        data->report.splitFaces+=area.liquid>0 && area.gas>0;
        data->faces.push_back(area);data->traces.push_back(trace);
    }
    return result;
}
HydroPhaseGeometry::Model HydroPhaseGeometry::model() const noexcept
{return data_->vof ? Model::GeometricVof : Model::FirstJet;}
const HydroPhysicalDivergence& HydroPhaseGeometry::metrics() const noexcept {return data_->metrics;}
const HydroDofLayout& HydroPhaseGeometry::layout() const noexcept {return data_->metrics.layout();}
const Interface::HermiteGalsPublication& HydroPhaseGeometry::publication() const
{require(model()==Model::FirstJet,"FirstJet publication port requires FirstJet");return *data_->publication;}
const Interface::GeometricVofPublication& HydroPhaseGeometry::vofPublication() const
{require(model()==Model::GeometricVof,"GeometricVof publication port requires GeometricVof");return *data_->vof;}
Interface::GeometricVofIdentity HydroPhaseGeometry::vofIdentity() const {return vofPublication().identity();}
Interface::HermiteGalsIdentity HydroPhaseGeometry::geometryIdentity() const {return publication().state().identity();}
Interface::HermiteGalsTimeLayer HydroPhaseGeometry::timeLayer() const
{return data_->vof ? data_->vof->timeLayer() : publication().state().timeLayer();}
const std::vector<HydroPhaseMeasure>& HydroPhaseGeometry::cellVolumes() const noexcept {return data_->cells;}
const std::vector<HydroPhaseMeasure>& HydroPhaseGeometry::faceMeasures() const noexcept {return data_->faces;}
const std::vector<HydroPhaseFaceTrace>& HydroPhaseGeometry::faceTraces() const noexcept {return data_->traces;}
const std::vector<double>& HydroPhaseGeometry::liquidFractions() const noexcept {return data_->alpha;}
const HydroPhaseGeometryReport& HydroPhaseGeometry::report() const noexcept {return data_->report;}
void HydroPhaseGeometry::validateFor(const QuadtreeMesh& mesh, Interface::HermiteGalsIdentity identity,
    Interface::HermiteGalsTimeLayer time) const
{
    require(identity==geometryIdentity(),"different geometry identity: expected="+
        std::to_string(geometryIdentity().value())+" supplied="+std::to_string(identity.value()));
    require(sameTime(time,timeLayer()),"different full time layer for geometry="+
        std::to_string(identity.value()));
    layout().validateMesh(mesh); publication().validateMesh(mesh);
}
void HydroPhaseGeometry::validateVofFor(const QuadtreeMesh& mesh,Interface::GeometricVofIdentity identity,
    Interface::GeometricVofTimeLayer time) const
{
    const auto& source=vofPublication();
    // build() requires the layout and publication to share one immutable
    // mesh binding. The source still validates the LIVE mesh on every call.
    source.validateFor(mesh,identity,time);
}
std::vector<HydroPhaseMeasure> HydroPhaseGeometry::measureVofBoxes(const QuadtreeMesh& mesh,
    const std::vector<HydroVofBox>& boxes) const
{return measureVofBoxesWithMeasure(mesh,boxes,false);}
std::vector<HydroPhaseMeasure> HydroPhaseGeometry::measureVofMeridionalBoxes(const QuadtreeMesh& mesh,
    const std::vector<HydroVofBox>& boxes) const
{return measureVofBoxesWithMeasure(mesh,boxes,true);}
std::vector<HydroPhaseMeasure> HydroPhaseGeometry::measureVofBoxesWithMeasure(const QuadtreeMesh& mesh,
    const std::vector<HydroVofBox>& boxes,bool meridional) const
{
    namespace Cut=LAMP::Geometry::PlanarCellCut;
    const auto& source=vofPublication();validateVofFor(mesh,source.identity(),source.timeLayer());
    LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
    const auto origin=source.frameOrigin();const bool rz=!meridional&&layout().domain().axisymmetric();
    require(!meridional||source.geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic,
        "meridional box area requires actual LocalPlic geometry");
    if(source.geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic){
        const auto& local=source.plicCells();
        const auto prepareBox=[&](const HydroVofBox& query){
            const int native=query.nativeCell;
            require(native>=0 && static_cast<std::size_t>(native)<source.cells().size() &&
                source.cells()[native].nativeCell==native && layout().pressureDof(native)>=0,
                "VOF box requires a selected original native leaf");
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
            const auto lower=mesh.verts(cell.cellToVert[0]).coord(),upper=mesh.verts(cell.cellToVert[2]).coord();
            const auto& box=query.rectangle;
            (void)Cut::capacity(box,rz); // Preserve the original capacity gate before native bounds.
            require(box.lower.x>=lower.x && box.lower.y>=lower.y && box.upper.x<=upper.x && box.upper.y<=upper.y,
                "VOF query box leaves its original native leaf");
            require(static_cast<std::size_t>(native)<local.size(),"missing local PLIC box source");
            return localPlicBox(local[native],box,rz,meridional);
        };
        bool parallel=false;
#ifdef _OPENMP
        parallel=boxes.size()>=32&&!omp_in_parallel()&&omp_get_max_threads()>1;
#endif
        std::vector<HydroPhaseMeasure> result;result.reserve(boxes.size());
        if(!parallel){
            for(const auto& query:boxes)result.push_back(prepareBox(query));
            return result;
        }
        result.resize(boxes.size());std::vector<std::exception_ptr> failures(boxes.size());
        const auto count=static_cast<std::ptrdiff_t>(boxes.size());
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            // Each actual worker authenticates its own rounding/FTZ/DAZ mode.
            // Neither a caller's query scope nor another worker grants permission.
            std::exception_ptr environmentFailure;
            try{LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();}
            catch(...){environmentFailure=std::current_exception();}
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for(std::ptrdiff_t i=0;i<count;++i){
                try{
                    if(environmentFailure)std::rethrow_exception(environmentFailure);
                    result[static_cast<std::size_t>(i)]=prepareBox(boxes[static_cast<std::size_t>(i)]);
                }catch(...){failures[static_cast<std::size_t>(i)]=std::current_exception();}
            }
        }
        // Query results are independent; only the original order selects failure.
        for(const auto& failure:failures)if(failure)std::rethrow_exception(failure);
        return result;
    }
    using Polygons=std::vector<std::vector<Cut::AnchoredPoint>>;
    std::unordered_map<int,std::array<Polygons,2>> cached;
    std::vector<HydroPhaseMeasure> result;result.reserve(boxes.size());
    for(const auto& query:boxes) {
        const int native=query.nativeCell;
        require(native>=0 && static_cast<std::size_t>(native)<source.cells().size() &&
            source.cells()[native].nativeCell==native && layout().pressureDof(native)>=0,
            "VOF box requires a selected original native leaf");
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
        const auto lower=mesh.verts(cell.cellToVert[0]).coord(),upper=mesh.verts(cell.cellToVert[2]).coord();
        const auto& box=query.rectangle;const auto capacity=Cut::capacity(box,rz);
        require(box.lower.x>=lower.x && box.lower.y>=lower.y && box.upper.x<=upper.x && box.upper.y<=upper.y,
            "VOF query box leaves its original native leaf");
        auto entry=cached.try_emplace(native);
        if(entry.second) {
            entry.first->second[0]=source.anchoredPhasePolygons(native,true);
            entry.first->second[1]=source.anchoredPhasePolygons(native,false);
        }
        std::array<PositiveVofSum,2> sum;
        for(unsigned phase=0;phase<2;++phase)for(const auto& polygon:entry.first->second[phase]) {
            // Select by exact original geometric support before integration.
            // Convexity of the query box proves this intersection is the whole
            // original polygon; this path constructs no new cut points.
            if(vofPolygonInsideBox(polygon,origin,box)) {
                sum[phase].add(Cut::anchoredPolygonMeasureInFrame(polygon,origin,rz));
                continue;
            }
            // Static Fubini coordinate integration: the top edge departs
            // vertically through this actual rectangle. Original y endpoints
            // and fan frame translation stay exact inside the Core kernel;
            // no rounded intersection points or physical time state are made.
            const auto integral=LAMP::Geometry::PlanarFaceTime::integrateAnchoredPolygon(polygon,
                {box.lower.x,box.upper.y},{box.upper.x,box.upper.y},
                {0,1},{0,1},box.lower.y,box.upper.y,{0,0},rz,512,{1,1,origin});
            require(integral.halfNegative==0 && integral.fullNegative==0,
                "negative static VOF coordinate integral");
            require(std::isfinite(integral.maximumRelativeUncertainty) &&
                integral.maximumRelativeUncertainty<=32.L*std::numeric_limits<double>::epsilon(),
                "uncertified static VOF coordinate integral");
            sum[phase].add(integral.fullPositive);
        }
        const HydroPhaseMeasure phases{sum[0].value,sum[1].value};
        vofRelativeClosure(phases,capacity);result.push_back(phases);
    }
    return result;
}
struct HydroPhaseRate::Data
{
    Data(HydroPhaseGeometry g,HydroSpatialAdvector a) : geometry(std::move(g)),spatial(std::move(a)) {}
    HydroPhaseGeometry geometry;
    std::vector<HydroPhaseMeasure> faces, cells;
    std::vector<double> total;
    HydroSpatialAdvector spatial;
    std::vector<HydroPhaseMeasure> atoms, dual;
    HydroSpatialPhaseReport report;
};
HydroPhaseRate::HydroPhaseRate(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
const HydroPhaseGeometry& HydroPhaseRate::geometry() const noexcept {return data_->geometry;}
const std::vector<HydroPhaseMeasure>& HydroPhaseRate::faceRates() const noexcept {return data_->faces;}
const std::vector<double>& HydroPhaseRate::totalFaceRates() const noexcept {return data_->total;}
const std::vector<HydroPhaseMeasure>& HydroPhaseRate::cellOutwardRates() const noexcept {return data_->cells;}
bool HydroPhaseRate::isSpatial() const noexcept {return true;}
const HydroSpatialAdvector& HydroPhaseRate::spatialAdvector() const
{ return data_->spatial; }
const HydroMomentumGeometry& HydroPhaseRate::momentumGeometry() const {return spatialAdvector().momentumGeometry();}
HydroAdvectorIdentity HydroPhaseRate::advectorIdentity() const {return spatialAdvector().identity();}
const std::vector<HydroPhaseMeasure>& HydroPhaseRate::atomRates() const {spatialAdvector();return data_->atoms;}
const std::vector<HydroPhaseMeasure>& HydroPhaseRate::dualFaceRates() const {spatialAdvector();return data_->dual;}
const HydroSpatialPhaseReport& HydroPhaseRate::spatialReport() const {spatialAdvector();return data_->report;}

namespace Detail
{
using Term=LAMP::Geometry::UnitPolynomialTerm;
// Metric density p(t)=a(t)*dA/dt from the ORIGINAL advector polynomial.
// Separate factors retain tiny coefficients and translated cancellation in
// Core's existing expansion primitive. No scalar/root data are reconstructed.
std::vector<Term> normalRatePolynomial(const HydroNormalPolynomial& p)
{
    std::vector<Term> result;long double power=p.length;
    for(std::size_t k=0;k<p.beginCoefficients.size();++k) {
        const double coefficient=static_cast<double>(p.beginCoefficients[k]);
        if(coefficient!=0) {
            Term term;term.coefficients[k]=coefficient;
            term.factor=power*(p.radialMeasure?CartesianMeasure::twoPi*p.radiusBegin:1.L);
            result.push_back(term);
            if(p.radialTangent) {
                term.coefficients[k]=0;term.coefficients[k+1]=coefficient;
                term.factor=power*p.length*CartesianMeasure::twoPi;result.push_back(term);
            }
        }
        power*=p.length;
    }
    return result;
}
}
namespace
{
using Term=LAMP::Geometry::UnitPolynomialTerm;
using Trace=LevelSetMappedSegment::NativeQuarticTrace;
constexpr long double spatialTolerance=1e-13L;
struct Sum
{
    long double value=0,correction=0;
    void add(long double x){const auto y=x-correction,t=value+y;correction=(t-value)-y;value=t;}
};
struct AtomIntegral
{
    HydroPhaseMeasure rates,work;
    long double error=0;
};
template<class PhaseTrace>
AtomIntegral integrateTrace(const PhaseTrace& trace,const HydroDofLayout& layout,
    const std::vector<Term>& terms,long double bound)
{
    AtomIntegral out;Sum liquid,gas;
    const auto addRate=[&](bool negative,long double value,long double work) {
        (negative?liquid:gas).add(value);(negative?out.work.liquid:out.work.gas)+=work;
    };
    for(const auto& r:trace.ordinaryRanges()) {
        const auto v=LAMP::Geometry::integrateUnitPolynomialRange(terms,0,1,r.fromEnd,r.offset,r.length);
        require(std::isfinite(v),"nonfinite ordinary spatial integral");addRate(r.negative,v,bound*r.length);
    }
    for(const auto& zero:trace.coincidentRanges())for(int side:{-1,1}) {
        auto signedTerms=terms;if(side>0)for(auto& term:signedTerms)term.weight=-term.weight;
        const auto integrate=[&](const auto& r) {
            const auto v=LAMP::Geometry::integratePositiveUnitPolynomialRange(signedTerms,0,1,
                r.fromEnd,r.offset,r.length,spatialTolerance*bound*r.length);
            require(std::isfinite(v.value) && v.value>=0 && std::isfinite(v.absoluteError),"invalid causal spatial integral");
            out.error+=v.absoluteError;return v;
        };
        const auto& s=side<0?zero.minus:zero.plus;
        const bool selected=s.available && s.nativeCellIndex>=0 && layout.pressureDof(s.nativeCellIndex)>=0;
        if(!selected) {
            const auto incoming=integrate(zero.interval);
            require(incoming.value==0 && incoming.absoluteError==0,
                "zero trace has incoming flow without a formal selected/exterior phase port");
        } else for(const auto& r:s.ranges) {
            const auto v=integrate(r);addRate(r.negative,side<0?v.value:-v.value,bound*r.length);
        }
    }
    out.rates={liquid.value,gas.value};return out;
}
long double phaseChange(const AtomIntegral& a,const AtomIntegral& b)
{
    long double change=0;
    for(int phase=0;phase<2;++phase) {
        const auto x=phase==0?a.rates.liquid:a.rates.gas,y=phase==0?b.rates.liquid:b.rates.gas;
        const auto scale=phase==0?std::max(a.work.liquid,b.work.liquid):std::max(a.work.gas,b.work.gas);
        require(scale>0 || x==y,"unscaled spatial phase refinement change");
        if(scale>0)change=std::max(change,std::abs(x-y)/scale);
    }
    return change;
}
template<class Operation> decltype(auto) visitPhaseTrace(Trace& trace,const Operation& operation)
{return operation(trace);}
template<class Operation> decltype(auto) visitPhaseTrace(
    Interface::HermiteGalsFirstJetTraceSample::Trace& trace,const Operation& operation)
{return std::visit(operation,trace);}
template<class Polynomial,class Modal,class TraceAt,class OwnsTrace>
std::vector<HydroPhaseMeasure> integrateAtomRatesPrepared(
    const HydroAdvectorReconstruction& reconstruction,const HydroDofLayout& layout,HydroSpatialPhaseReport& report,
    const std::vector<int>* requested,const Polynomial& polynomialAt,const Modal& modalAt,
    const TraceAt& traceAt,const OwnsTrace& ownsTrace,bool parallel=false)
{
    report.atoms=requested?requested->size():reconstruction.atoms().size();
    report.tangentUnidentified=reconstruction.report().tangentUnidentified;
    std::vector<HydroPhaseMeasure> rates(report.atoms);
    const auto atomRate=[&](std::size_t q,HydroSpatialPhaseReport& atomReport) {
        const std::size_t i=requested?static_cast<std::size_t>((*requested)[q]):q;
        const auto polynomial=polynomialAt(static_cast<int>(i));
        auto retained=traceAt(q,static_cast<int>(i));
        require(ownsTrace(retained),"spatial trace changed original scalar");++atomReport.sourceQueries;
        return visitPhaseTrace(retained,[&](auto& trace)->HydroPhaseMeasure {
        // A retained strict single-phase trace makes the phase indicator
        // constant on the complete ORIGINAL atom. Preserve that certificate
        // and physical total; no unit-polynomial enclosure/integration is
        // needed when the independently queried modal contribution is zero.
        const auto& ranges=trace.ordinaryRanges();
        if(ranges.size()==1&&ranges[0].offset==0&&ranges[0].length==1
            &&trace.events().empty()&&trace.coincidentRanges().empty()) {
            LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();
            if(modalAt(static_cast<int>(i),trace)==0) {
                const auto total=polynomial.integral(false,0,polynomial.length);
                require(std::isfinite(total),"nonfinite single-phase spatial total");
                ++atomReport.ordinaryRanges;
                const HydroPhaseMeasure pureRate=ranges[0].negative?HydroPhaseMeasure{total,0}:HydroPhaseMeasure{0,total};
                LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();
                return pureRate;
            }
        }
        const auto terms=Detail::normalRatePolynomial(polynomial);const auto bounds=LAMP::Geometry::encloseUnitPolynomial(terms);
        const long double bound=std::max(std::abs(bounds[0]),std::abs(bounds[1]));
        require(std::isfinite(bound),"nonfinite spatial polynomial enclosure");
        auto integral=integrateTrace(trace,layout,terms,bound);bool stable=false;long double change=0;
        for(int sweep=0;sweep<4;++sweep) {
            bool refined=false;
            // Re-read the CURRENT event inventory after every replacement.
            // No ordinal equality between independent side charts is used.
            for(std::size_t k=0;k<trace.events().size();++k) {
                const auto root=trace.events()[k].root;
                if(!root.pointOnly && root.lower<root.upper) {
                    trace=trace.refined(k);++atomReport.refinements;refined=true;
                    require(ownsTrace(retained),"refinement changed original scalar");
                }
            }
            if(!refined){stable=true;break;}
            auto next=integrateTrace(trace,layout,terms,bound);change=phaseChange(integral,next);integral=next;
            if(change<=spatialTolerance){stable=true;break;}
        }
        require(stable,"spatial phase refinement exceeds bounded work/accuracy");
        atomReport.maxRelativeRefinementChange=std::max(atomReport.maxRelativeRefinementChange,change);
        atomReport.positiveIntegrationErrorEstimate+=integral.error;
        atomReport.ordinaryRanges+=trace.ordinaryRanges().size();atomReport.coincidentRanges+=trace.coincidentRanges().size();atomReport.events+=trace.events().size();
        const auto total=polynomial.integral(false,0,polynomial.length);
        const auto error=std::abs(integral.rates.liquid+integral.rates.gas-total);
        require(std::isfinite(error) && (bound>0?error<=1e-12L*bound:error==0),"spatial phase sum differs from actual total");
        if(bound>0)atomReport.maxRelativeAtomClosure=std::max(atomReport.maxRelativeAtomClosure,error/bound);
        const long double modal=modalAt(static_cast<int>(i),trace);
        integral.rates.gas+=modal;
        require(std::isfinite(integral.rates.gas),"modal spatial gas rate overflow");
        const long double combinedError=std::abs(integral.rates.liquid+integral.rates.gas-total-modal),combinedScale=bound+std::abs(modal);
        require(combinedScale>0?combinedError<=1e-12L*combinedScale:combinedError==0,"combined modal phase sum differs from physical total");
        if(combinedScale>0)atomReport.maxRelativeAtomClosure=std::max(atomReport.maxRelativeAtomClosure,combinedError/combinedScale);
        if(modal!=0){++atomReport.tangentialAtoms;atomReport.maxAbsTangentialRate=std::max(atomReport.maxAbsTangentialRate,std::abs(modal));}
        return integral.rates;
        });
    };
    bool runParallel=false;
#ifdef _OPENMP
    runParallel=parallel&&report.atoms>1&&!omp_in_parallel()&&omp_get_max_threads()>1;
#else
    (void)parallel;
#endif
    if(!runParallel) {
        // Retain the original ordered report arithmetic on all old ports.
        for(std::size_t q=0;q<report.atoms;++q)rates[q]=atomRate(q,report);
        return rates;
    }
    std::vector<HydroSpatialPhaseReport> atomReports(report.atoms);
    std::vector<std::exception_ptr> failures(report.atoms);
    const auto prepareAtom=[&](std::size_t q) {
        try {
            LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();
            rates[q]=atomRate(q,atomReports[q]);
        } catch(...) { failures[q]=std::current_exception(); }
    };
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for(std::ptrdiff_t q=0;q<static_cast<std::ptrdiff_t>(report.atoms);++q)
        prepareAtom(static_cast<std::size_t>(q));
    // Immutable requested atoms share no result/report/exception slots.
    // Commit in the original atom order, including the first failure.
    for(std::size_t q=0;q<report.atoms;++q) {
        if(failures[q])std::rethrow_exception(failures[q]);
        const auto& a=atomReports[q];
        require(std::isfinite(a.positiveIntegrationErrorEstimate)&&std::isfinite(a.maxRelativeRefinementChange)
            &&std::isfinite(a.maxRelativeAtomClosure)&&std::isfinite(a.maxAbsTangentialRate),"spatial atom report overflow");
        report.sourceQueries+=a.sourceQueries;report.ordinaryRanges+=a.ordinaryRanges;
        report.coincidentRanges+=a.coincidentRanges;report.events+=a.events;report.refinements+=a.refinements;
        report.tangentialAtoms+=a.tangentialAtoms;
        report.positiveIntegrationErrorEstimate+=a.positiveIntegrationErrorEstimate;
        report.maxRelativeRefinementChange=std::max(report.maxRelativeRefinementChange,a.maxRelativeRefinementChange);
        report.maxRelativeAtomClosure=std::max(report.maxRelativeAtomClosure,a.maxRelativeAtomClosure);
        report.maxAbsTangentialRate=std::max(report.maxAbsTangentialRate,a.maxAbsTangentialRate);
    }
    return rates;
}
std::vector<HydroPhaseMeasure> integrateAtomRates(
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& scalar,
    const HydroSpatialAdvector& advector,const HydroDofLayout& layout,HydroSpatialPhaseReport& report,
    const std::vector<int>* requested=nullptr,bool allowParallel=false)
{
    return integrateAtomRatesPrepared(advector.reconstruction(),layout,report,requested,
        [&](int i){return advector.regularAtomPolynomials()[i];},
        [&](int i,const auto& trace){return advector.tangentialAtomRate(i,trace);},
        [&](std::size_t,int i){const auto& atom=advector.reconstruction().atoms()[i];
            return LevelSetMappedSegment::nativeQuarticTrace(scalar,atom.begin,atom.end,atom.normalAxis);},
        [&](const auto& trace){return trace.usesSource(*scalar);},allowParallel);
}
}
HydroPhaseRate HydroPhaseRate::prepareSpatial(const QuadtreeMesh& mesh,HydroPhaseGeometry geometry,
    const HydroSpatialAdvector& advector)
{
    geometry.validateFor(mesh,advector.geometryIdentity(),advector.timeLayer());
    advector.validateFor(mesh,advector.momentumGeometry(),geometry.geometryIdentity(),geometry.timeLayer());
    require(geometry.layout().sameIdentity(advector.values().layout()),"different spatial native layout identity");
    require(&geometry.metrics().divergence().matrix()==&advector.momentumGeometry().metrics().divergence().matrix(),
        "different original D_I instance");
    auto d=std::make_shared<Data>(std::move(geometry),advector);
    const auto view=std::make_shared<const Interface::HermiteGalsGeometryView>(d->geometry.publication().state().geometry());
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar(view,&view->snapshot());
    const auto& reconstruction=advector.reconstruction();const auto& layout=d->geometry.layout();
    d->atoms=integrateAtomRates(scalar,advector,layout,d->report,nullptr,true);
    for(int u=0;u<layout.velocityCount();++u) {
        Sum liquid,gas;for(int atom:reconstruction.primalAtoms()[u]){liquid.add(d->atoms[atom].liquid);gas.add(d->atoms[atom].gas);}
        require(std::isfinite(liquid.value) && std::isfinite(gas.value),"spatial primal phase rate overflow");
        d->faces.push_back({liquid.value,gas.value});
        const auto total=d->geometry.metrics().faceMeasures()[u]*advector.values().normalVelocity()[u];
        require(std::isfinite(total),"nonfinite original native total rate");d->total.push_back(total);
    }
    for(const auto& ids:reconstruction.dualAtoms()) {
        Sum liquid,gas;for(int atom:ids){liquid.add(d->atoms[atom].liquid);gas.add(d->atoms[atom].gas);}
        require(std::isfinite(liquid.value) && std::isfinite(gas.value),"spatial dual phase rate overflow");
        d->dual.push_back({liquid.value,gas.value});
    }
    d->cells.resize(layout.pressureCount());const auto& di=d->geometry.metrics().divergence().matrix();
    for(Eigen::Index p=0;p<di.outerSize();++p)for(HydroIntegralDivergence::Matrix::InnerIterator it(di,p);it;++it) {
        const auto& r=d->faces[static_cast<std::size_t>(it.col())];const int sign=it.value()>0?1:-1;
        d->cells[static_cast<std::size_t>(p)].liquid+=sign*r.liquid;d->cells[static_cast<std::size_t>(p)].gas+=sign*r.gas;
    }
    for(const auto& r:d->cells)require(std::isfinite(r.liquid) && std::isfinite(r.gas),"spatial cell rate overflow");
    require(std::isfinite(d->report.positiveIntegrationErrorEstimate),"spatial integration estimate overflow");
    return HydroPhaseRate(std::move(d));
}
std::vector<HydroPhaseMeasure> HydroPhaseRate::sampleAtomRates(const QuadtreeMesh& mesh,
    const Interface::HermiteGalsFirstJetSample& sample,const HydroSpatialAdvector& advector)
{
    // Retain the original spatial owner's mesh, full clock and FP checks.
    advector.validateFor(mesh,advector.momentumGeometry(),sample.identity(),sample.timeLayer());
    sample.validateMesh(mesh);
    const auto& layout=advector.momentumGeometry().layout();
    advector.values().validateFor(layout,mesh);
    require(layout.domain().axisymmetric()==(sample.measure()==Interface::HermiteGalsMeasure::AxisymmetricVolume),
        "different XY/RZ measure");
    require(sample.sourceIdentity().value()!=0 && sample.arrivalVelocityIdentity().value()!=0,
        "unbound sample provenance");
    const auto view=std::make_shared<const Interface::HermiteGalsGeometryView>(sample.geometry());
    require(view->identity()==sample.identity() && sameTime(view->timeLayer(),sample.timeLayer()),
        "sample scalar provenance mismatch");
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar(view,&view->snapshot());
    require(scalar->usesNativeLeaves() && (scalar->reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic||scalar->usesC2Macro()),
        "sample requires native Q4 scalar");
    for(int p=0;p<layout.pressureCount();++p)
        require(scalar->nativeCellActive(layout.nativeCell(p)),"missing selected-cell scalar geometry");
    HydroSpatialPhaseReport report;
    auto rates=integrateAtomRates(scalar,advector,layout,report);
    require(std::isfinite(report.positiveIntegrationErrorEstimate),"spatial integration estimate overflow");
    return rates;
}
std::vector<int> HydroPhaseRate::requiredNativeCells(const QuadtreeMesh& mesh,
    const HydroAdvectorReconstruction& reconstruction,const std::vector<int>& atomIds)
{
    const auto& layout=reconstruction.momentumGeometry().layout();layout.validateMesh(mesh);
    require(!layout.domain().axisymmetric(),"local rates require uniform XY mesh");
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    require(grid.nx()>=2&&grid.ny()>=2,"local rates require compact uniform support");
    require(!atomIds.empty(),"empty requested atoms");
    std::vector<unsigned char> seen(reconstruction.atoms().size(),0),cells(mesh.cellIndexNum(),0);
    for(int id:atomIds){require(id>=0&&static_cast<std::size_t>(id)<seen.size(),"invalid atom");
        require(!seen[id],"duplicate atom");seen[id]=1;}
    const double infinity=std::numeric_limits<double>::infinity();
    std::vector<double> lowerX(grid.nx(),infinity),upperX(grid.nx(),-infinity),
        lowerY(grid.ny(),infinity),upperY(grid.ny(),-infinity);
    // Retain the actual corner coordinates, including any row/column ULP
    // differences allowed by the validated uniform lattice. Its resolved
    // spacing makes these column/row envelopes ordered without rounding an
    // atom coordinate to a logical index.
    for(int j=0;j<grid.ny();++j)for(int i=0;i<grid.nx();++i){
        const auto& cell=mesh.cells(grid.cells()[j*grid.nx()+i]);
        const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
        lowerX[i]=std::min(lowerX[i],std::min(a.x,b.x));upperX[i]=std::max(upperX[i],std::max(a.x,b.x));
        lowerY[j]=std::min(lowerY[j],std::min(a.y,b.y));upperY[j]=std::max(upperY[j],std::max(a.y,b.y));
    }
    for(int id:atomIds){const auto& atom=reconstruction.atoms()[id];
        require((atom.normalAxis==0||atom.normalAxis==1)
            &&std::isfinite(atom.begin.x)&&std::isfinite(atom.begin.y)&&std::isfinite(atom.end.x)&&std::isfinite(atom.end.y),
            "invalid actual atom geometry");
        const Vector2D lo{std::min(atom.begin.x,atom.end.x),std::min(atom.begin.y,atom.end.y)},
            hi{std::max(atom.begin.x,atom.end.x),std::max(atom.begin.y,atom.end.y)};
        // Closed bounds include both adjacent cells at a coincident side or
        // endpoint; at the outer domain boundary only the legal side exists.
        const int firstX=static_cast<int>(std::lower_bound(upperX.begin(),upperX.end(),lo.x)-upperX.begin()),
            endX=static_cast<int>(std::upper_bound(lowerX.begin(),lowerX.end(),hi.x)-lowerX.begin()),
            firstY=static_cast<int>(std::lower_bound(upperY.begin(),upperY.end(),lo.y)-upperY.begin()),
            endY=static_cast<int>(std::upper_bound(lowerY.begin(),lowerY.end(),hi.y)-lowerY.begin());
        bool found=false;
        for(int j=firstY;j<endY;++j)for(int i=firstX;i<endX;++i){
            const auto& cell=mesh.cells(grid.cells()[j*grid.nx()+i]);const int c=cell.index();
            const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
            const Vector2D lower{std::min(a.x,b.x),std::min(a.y,b.y)},upper{std::max(a.x,b.x),std::max(a.y,b.y)};
            if(hi.x>=lower.x&&lo.x<=upper.x&&hi.y>=lower.y&&lo.y<=upper.y){cells[c]=1;found=true;}}
        require(found,"actual atom has no native support");}
    std::vector<int> result;for(int c=0;c<mesh.cellIndexNum();++c)if(cells[c])result.push_back(c);
    return result;
}

std::vector<HydroPhaseMeasure> HydroPhaseRate::sampleAtomRates(const QuadtreeMesh& mesh,
    const Interface::HermiteGalsFirstJetLocalSample& sample,const HydroSpatialAdvector& advector,
    const std::vector<int>& atomIds)
{
    advector.validateFor(mesh,advector.momentumGeometry(),sample.identity(),sample.timeLayer());
    sample.validateMesh(mesh);const auto& layout=advector.momentumGeometry().layout();
    advector.values().validateFor(layout,mesh);
    require(!layout.domain().axisymmetric()&&sample.measure()==Interface::HermiteGalsMeasure::CartesianArea,
        "local rates require uniform XY measure");
    require(sample.sourceIdentity().value()!=0&&sample.arrivalVelocityIdentity().value()!=0,"unbound local sample provenance");
    const auto required=requiredNativeCells(mesh,advector.reconstruction(),atomIds);
    const auto view=std::make_shared<const Interface::HermiteGalsGeometryView>(sample.geometry());
    require(view->identity()==sample.identity()&&sameTime(view->timeLayer(),sample.timeLayer()),"local sample scalar provenance mismatch");
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar(view,&view->snapshot());
    require(scalar->usesNativeLeaves()&&!scalar->hasCornerSectors()
        &&scalar->reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,"local sample requires native Q4 scalar");
    // Check before the per-atom kernel, even if its actual normal work is zero.
    for(int c:required)require(scalar->nativeCellActive(c),"missing local native support");
    HydroSpatialPhaseReport report;auto rates=integrateAtomRates(scalar,advector,layout,report,&atomIds);
    require(std::isfinite(report.positiveIntegrationErrorEstimate),"spatial integration estimate overflow");
    return rates;
}

std::vector<HydroPhaseMeasure> HydroPhaseRate::sampleAtomRates(
    const Interface::HermiteGalsFirstJetLocalSample& sample,const HydroAdvectorTimeSlab& slab,
    const std::vector<int>& atomIds)
{
    slab.validateSample(sample.identity(),sample.sourceIdentity(),sample.arrivalVelocityIdentity(),sample.timeLayer());
    return sampleAtomRates(sample,slab,slab.prepareAtomSupport(atomIds));
}
std::vector<HydroPhaseMeasure> HydroPhaseRate::sampleAtomRates(
    const Interface::HermiteGalsFirstJetLocalSample& sample,const HydroAdvectorTimeSlab& slab,
    const HydroAdvectorAtomSupport& support)
{
    const auto& atomIds=support.atomIds();
    slab.validateSample(sample.identity(),sample.sourceIdentity(),sample.arrivalVelocityIdentity(),sample.timeLayer());
    const auto& reconstruction=slab.reconstruction();const auto& layout=reconstruction.momentumGeometry().layout();
    require(!layout.domain().axisymmetric()&&sample.measure()==Interface::HermiteGalsMeasure::CartesianArea,
        "local rates require uniform XY measure");
    slab.validateAtomSupport(support);
    const auto& required=support.nativeCells();
    const auto view=std::make_shared<const Interface::HermiteGalsGeometryView>(sample.geometry());
    require(view->identity()==sample.identity()&&sameTime(view->timeLayer(),sample.timeLayer()),"local sample scalar provenance mismatch");
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar(view,&view->snapshot());
    require(scalar->usesNativeLeaves()&&!scalar->hasCornerSectors()
        &&scalar->reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,"local sample requires native Q4 scalar");
    require(static_cast<std::size_t>(scalar->nativeCellCount())==slab.nativeCellCount(),"different local native mesh");
    for(int c:required){require(scalar->nativeCellActive(c),"missing local native support");
        const auto a=scalar->nativeCellBounds(c),b=slab.nativeCellBounds(c);
        require(a[0].x==b[0].x&&a[0].y==b[0].y&&a[1].x==b[1].x&&a[1].y==b[1].y,"different local native bounds");}
    HydroSpatialPhaseReport report;
    auto rates=integrateAtomRatesPrepared(reconstruction,layout,report,&atomIds,
        [&](int i){return slab.atomPolynomial(i,sample.timeLayer());},
        [](int,const auto&){return 0.L;},
        [&](std::size_t,int i){const auto& atom=reconstruction.atoms()[i];
            return LevelSetMappedSegment::nativeQuarticTrace(scalar,atom.begin,atom.end,atom.normalAxis);},
        [&](const auto& trace){return trace.usesSource(*scalar);},true); // endpoint slab certified zero modal coefficients
    require(std::isfinite(report.positiveIntegrationErrorEstimate),"spatial integration estimate overflow");
    return rates;
}

std::vector<HydroPhaseMeasure> HydroPhaseRate::sampleAtomRates(
    const Interface::HermiteGalsFirstJetTraceSample& sample,const HydroAdvectorTimeSlab& slab,
    const std::vector<int>& atomIds)
{
    slab.validateSample(sample.sourceIdentity(),sample.sourceIdentity(),sample.arrivalVelocityIdentity(),sample.timeLayer());
    return sampleAtomRates(sample,slab,slab.prepareAtomSupport(atomIds));
}
std::vector<HydroPhaseMeasure> HydroPhaseRate::sampleAtomRates(
    const Interface::HermiteGalsFirstJetTraceSample& sample,const HydroAdvectorTimeSlab& slab,
    const HydroAdvectorAtomSupport& support)
{
    const auto& atomIds=support.atomIds();
    LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();
    // The immutable trace sample is derived from this accepted source, with
    // no new full-domain scalar/state identity. Its owner authenticates every
    // actual trace, including a retained full-Q4 fallback and refinements.
    slab.validateSample(sample.sourceIdentity(),sample.sourceIdentity(),
        sample.arrivalVelocityIdentity(),sample.timeLayer());
    const auto& reconstruction=slab.reconstruction();const auto& layout=reconstruction.momentumGeometry().layout();
    require(!layout.domain().axisymmetric()&&sample.measure()==Interface::HermiteGalsMeasure::CartesianArea,
        "local trace rates require uniform XY measure");
    const auto sourceTime=sample.sourceTimeLayer(),endpoint=slab.endpointTimeLayer();
    require(sourceTime.step>=0&&sourceTime.step<endpoint.step&&endpoint.step-sourceTime.step==1
        &&(sourceTime.position==Interface::HermiteGalsTimeLayer::Position::Initial
            ||sourceTime.position==Interface::HermiteGalsTimeLayer::Position::End)
        &&sourceTime.time()==endpoint.begin,"different original source clock");
    slab.validateAtomSupport(support);
    const auto& required=support.nativeCells();
    require(sample.nativeCellCount()==slab.nativeCellCount(),"different local native mesh");
    const auto& available=sample.nativeCells();
    const bool ordered=std::is_sorted(available.begin(),available.end());
    auto cursor=available.begin();
    for(int c:required){
        if(ordered) {
            while(cursor!=available.end()&&*cursor<c)++cursor;
            require(cursor!=available.end()&&*cursor==c,"missing local native support");
        }else require(std::find(available.begin(),available.end(),c)!=available.end(),"missing local native support");
        const auto a=sample.nativeCellBounds(c),b=slab.nativeCellBounds(c);
        require(a[0].x==b[0].x&&a[0].y==b[0].y&&a[1].x==b[1].x&&a[1].y==b[1].y,"different local native bounds");
    }
    HydroSpatialPhaseReport report;
    auto rates=integrateAtomRatesPrepared(reconstruction,layout,report,&atomIds,
        [&](int i){return slab.atomPolynomial(i,sample.timeLayer());},
        [](int,const auto&){return 0.L;},
        [&](std::size_t q,int i){
            auto trace=sample.trace(q);const auto& atom=reconstruction.atoms()[i];
            std::visit([&](const auto& actual){
                const auto a=actual.begin(),b=actual.end();
                require(a.x==atom.begin.x&&a.y==atom.begin.y&&b.x==atom.end.x&&b.y==atom.end.y
                    &&actual.normalAxis()==atom.normalAxis,"different requested physical atom trace");
            },trace);
            return trace;
        },[&](const auto& trace){return sample.ownsTrace(trace);},true);
    require(std::isfinite(report.positiveIntegrationErrorEstimate),"spatial integration estimate overflow");
    LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();
    return rates;
}
}
