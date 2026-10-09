#include "Modules/Hydrodynamics/CchPhaseGeometry.h"

#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/UnitPolynomialRoots.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Modules/Hydrodynamics/CchSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/CchNormalRatePolynomial.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Hydrodynamics
{
namespace
{
void require(bool ok, const std::string& reason)
{
    if (!ok) throw std::invalid_argument("[CchPhaseGeometry] "+reason);
}
bool sameTime(Interface::HermiteGalsTimeLayer a, Interface::HermiteGalsTimeLayer b)
{
    a.time(); b.time(); // Validate every field before comparison, including NaN.
    return a.step==b.step && a.begin==b.begin && a.end==b.end && a.position==b.position;
}
CchPhaseMeasure measure(Vector2D a, Vector2D b, const std::vector<LevelSetMappedSegment::Range>& ranges, bool radial)
{
    const long double dx=static_cast<long double>(b.x)-a.x, dy=static_cast<long double>(b.y)-a.y;
    const long double length=std::hypot(dx,dy);
    require(std::isfinite(length) && length>0,"invalid retained half-face length");
    CchPhaseMeasure result;
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
CchPhaseMeasure measure(const LevelSetMappedSegment& segment, bool radial)
{ return measure(segment.begin(),segment.end(),segment.ranges(),radial); }
void add(CchPhaseMeasure& a,CchPhaseMeasure b){a.liquid+=b.liquid;a.gas+=b.gas;}
}
struct CchPhaseGeometry::Data
{
    Data(CchPhysicalDivergence m, Interface::HermiteGalsPublication p)
        : metrics(std::move(m)), publication(std::move(p)) {}
    CchPhysicalDivergence metrics;
    Interface::HermiteGalsPublication publication;
    std::vector<CchPhaseMeasure> cells, faces;
    std::vector<CchPhaseFaceTrace> traces;
    std::vector<double> alpha;
    CchPhaseGeometryReport report;
};
CchPhaseGeometry::CchPhaseGeometry(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
CchPhaseGeometry CchPhaseGeometry::build(const QuadtreeMesh& mesh, CchPhysicalDivergence metrics,
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
        data->publication.domain().cellRowByIndex().size(),nullptr);
    for (const auto& cell : data->publication.volumes().cells()) {
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
        const auto& partitions=data->publication.facePartitions();
        const auto face=layout.domain().face(mesh,native);
        CchPhaseMeasure area;CchPhaseFaceTrace trace;
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
                CchPhaseMeasure minus=ordinary,plus=ordinary;
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
            trace.kind=trace.minusAvailable() && trace.plusAvailable()?CchPhaseFaceTraceKind::InteriorCoincident:CchPhaseFaceTraceKind::BoundaryCoincident;
            data->report.boundaryCoincidentFaces+=trace.kind==CchPhaseFaceTraceKind::BoundaryCoincident;
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
    return CchPhaseGeometry(std::move(data));
}
const CchPhysicalDivergence& CchPhaseGeometry::metrics() const noexcept {return data_->metrics;}
const CchDofLayout& CchPhaseGeometry::layout() const noexcept {return data_->metrics.layout();}
const Interface::HermiteGalsPublication& CchPhaseGeometry::publication() const noexcept {return data_->publication;}
Interface::HermiteGalsIdentity CchPhaseGeometry::geometryIdentity() const {return publication().state().identity();}
Interface::HermiteGalsTimeLayer CchPhaseGeometry::timeLayer() const {return publication().state().timeLayer();}
const std::vector<CchPhaseMeasure>& CchPhaseGeometry::cellVolumes() const noexcept {return data_->cells;}
const std::vector<CchPhaseMeasure>& CchPhaseGeometry::faceMeasures() const noexcept {return data_->faces;}
const std::vector<CchPhaseFaceTrace>& CchPhaseGeometry::faceTraces() const noexcept {return data_->traces;}
const std::vector<double>& CchPhaseGeometry::liquidFractions() const noexcept {return data_->alpha;}
const CchPhaseGeometryReport& CchPhaseGeometry::report() const noexcept {return data_->report;}
void CchPhaseGeometry::validateFor(const QuadtreeMesh& mesh, Interface::HermiteGalsIdentity identity,
    Interface::HermiteGalsTimeLayer time) const
{
    require(identity==geometryIdentity(),"different geometry identity: expected="+
        std::to_string(geometryIdentity().value())+" supplied="+std::to_string(identity.value()));
    require(sameTime(time,timeLayer()),"different full time layer for geometry="+
        std::to_string(identity.value()));
    layout().validateMesh(mesh); publication().validateMesh(mesh);
}
struct CchPhaseRate::Data
{
    Data(CchPhaseGeometry g,CchSpatialAdvector a) : geometry(std::move(g)),spatial(std::move(a)) {}
    CchPhaseGeometry geometry;
    std::vector<CchPhaseMeasure> faces, cells;
    std::vector<double> total;
    CchSpatialAdvector spatial;
    std::vector<CchPhaseMeasure> atoms, dual;
    CchSpatialPhaseReport report;
};
CchPhaseRate::CchPhaseRate(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
const CchPhaseGeometry& CchPhaseRate::geometry() const noexcept {return data_->geometry;}
const std::vector<CchPhaseMeasure>& CchPhaseRate::faceRates() const noexcept {return data_->faces;}
const std::vector<double>& CchPhaseRate::totalFaceRates() const noexcept {return data_->total;}
const std::vector<CchPhaseMeasure>& CchPhaseRate::cellOutwardRates() const noexcept {return data_->cells;}
bool CchPhaseRate::isSpatial() const noexcept {return true;}
const CchSpatialAdvector& CchPhaseRate::spatialAdvector() const
{ return data_->spatial; }
const CchMomentumGeometry& CchPhaseRate::momentumGeometry() const {return spatialAdvector().momentumGeometry();}
CchAdvectorIdentity CchPhaseRate::advectorIdentity() const {return spatialAdvector().identity();}
const std::vector<CchPhaseMeasure>& CchPhaseRate::atomRates() const {spatialAdvector();return data_->atoms;}
const std::vector<CchPhaseMeasure>& CchPhaseRate::dualFaceRates() const {spatialAdvector();return data_->dual;}
const CchSpatialPhaseReport& CchPhaseRate::spatialReport() const {spatialAdvector();return data_->report;}

namespace Detail
{
using Term=Quantum::Geometry::UnitPolynomialTerm;
// Metric density p(t)=a(t)*dA/dt from the ORIGINAL advector polynomial.
// Separate factors retain tiny coefficients and translated cancellation in
// Core's existing expansion primitive. No scalar/root data are reconstructed.
std::vector<Term> normalRatePolynomial(const CchNormalPolynomial& p)
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
using Term=Quantum::Geometry::UnitPolynomialTerm;
using Trace=LevelSetMappedSegment::NativeQuarticTrace;
constexpr long double spatialTolerance=1e-13L;
struct Sum
{
    long double value=0,correction=0;
    void add(long double x){const auto y=x-correction,t=value+y;correction=(t-value)-y;value=t;}
};
struct AtomIntegral
{
    CchPhaseMeasure rates,work;
    long double error=0;
};
AtomIntegral integrateTrace(const Trace& trace,const CchDofLayout& layout,
    const std::vector<Term>& terms,long double bound)
{
    AtomIntegral out;Sum liquid,gas;
    const auto addRate=[&](bool negative,long double value,long double work) {
        (negative?liquid:gas).add(value);(negative?out.work.liquid:out.work.gas)+=work;
    };
    for(const auto& r:trace.ordinaryRanges()) {
        const auto v=Quantum::Geometry::integrateUnitPolynomialRange(terms,0,1,r.fromEnd,r.offset,r.length);
        require(std::isfinite(v),"nonfinite ordinary spatial integral");addRate(r.negative,v,bound*r.length);
    }
    for(const auto& zero:trace.coincidentRanges())for(int side:{-1,1}) {
        auto signedTerms=terms;if(side>0)for(auto& term:signedTerms)term.weight=-term.weight;
        const auto integrate=[&](const auto& r) {
            const auto v=Quantum::Geometry::integratePositiveUnitPolynomialRange(signedTerms,0,1,
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
}
CchPhaseRate CchPhaseRate::prepareSpatial(const QuadtreeMesh& mesh,CchPhaseGeometry geometry,
    const CchSpatialAdvector& advector)
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
    d->report.atoms=reconstruction.atoms().size();d->report.tangentUnidentified=reconstruction.report().tangentUnidentified;
    for(std::size_t i=0;i<reconstruction.atoms().size();++i) {
        const auto& atom=reconstruction.atoms()[i];const auto& polynomial=advector.regularAtomPolynomials()[i];
        auto trace=LevelSetMappedSegment::nativeQuarticTrace(scalar,atom.begin,atom.end,atom.normalAxis);
        require(trace.usesSource(*scalar),"spatial trace changed original scalar");++d->report.sourceQueries;
        const auto terms=Detail::normalRatePolynomial(polynomial);const auto bounds=Quantum::Geometry::encloseUnitPolynomial(terms);
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
                    trace=trace.refined(k);++d->report.refinements;refined=true;
                    require(trace.usesSource(*scalar),"refinement changed original scalar");
                }
            }
            if(!refined){stable=true;break;}
            auto next=integrateTrace(trace,layout,terms,bound);change=phaseChange(integral,next);integral=next;
            if(change<=spatialTolerance){stable=true;break;}
        }
        require(stable,"spatial phase refinement exceeds bounded work/accuracy");
        d->report.maxRelativeRefinementChange=std::max(d->report.maxRelativeRefinementChange,change);
        d->report.positiveIntegrationErrorEstimate+=integral.error;
        d->report.ordinaryRanges+=trace.ordinaryRanges().size();d->report.coincidentRanges+=trace.coincidentRanges().size();d->report.events+=trace.events().size();
        const auto total=polynomial.integral(false,0,polynomial.length);
        const auto error=std::abs(integral.rates.liquid+integral.rates.gas-total);
        require(std::isfinite(error) && (bound>0?error<=1e-12L*bound:error==0),"spatial phase sum differs from actual total");
        if(bound>0)d->report.maxRelativeAtomClosure=std::max(d->report.maxRelativeAtomClosure,error/bound);
        const long double modal=advector.tangentialAtomRate(static_cast<int>(i),trace);
        integral.rates.gas+=modal;
        require(std::isfinite(integral.rates.gas),"modal spatial gas rate overflow");
        const long double combinedError=std::abs(integral.rates.liquid+integral.rates.gas-total-modal),combinedScale=bound+std::abs(modal);
        require(combinedScale>0?combinedError<=1e-12L*combinedScale:combinedError==0,"combined modal phase sum differs from physical total");
        if(combinedScale>0)d->report.maxRelativeAtomClosure=std::max(d->report.maxRelativeAtomClosure,combinedError/combinedScale);
        if(modal!=0){++d->report.tangentialAtoms;d->report.maxAbsTangentialRate=std::max(d->report.maxAbsTangentialRate,std::abs(modal));}
        d->atoms.push_back(integral.rates);
    }
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
    for(Eigen::Index p=0;p<di.outerSize();++p)for(CchIntegralDivergence::Matrix::InnerIterator it(di,p);it;++it) {
        const auto& r=d->faces[static_cast<std::size_t>(it.col())];const int sign=it.value()>0?1:-1;
        d->cells[static_cast<std::size_t>(p)].liquid+=sign*r.liquid;d->cells[static_cast<std::size_t>(p)].gas+=sign*r.gas;
    }
    for(const auto& r:d->cells)require(std::isfinite(r.liquid) && std::isfinite(r.gas),"spatial cell rate overflow");
    require(std::isfinite(d->report.positiveIntegrationErrorEstimate),"spatial integration estimate overflow");
    return CchPhaseRate(std::move(d));
}
}
