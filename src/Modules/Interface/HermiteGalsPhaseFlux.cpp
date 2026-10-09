#include "Modules/Interface/HermiteGalsPhaseFlux.h"
#include "Modules/Interface/Detail/HermiteGalsInventoryAssembly.h"

#include "Modules/Interface/InterfacePhaseGeometry.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/UnitPolynomialRoots.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Interface
{
namespace
{
using Term=Quantum::Geometry::UnitPolynomialTerm;
using Boundary=HermiteGalsFaceVelocity::Boundary;
constexpr long double roundoffBudget=Detail::hermiteInventoryRoundoffBudget;
void require(bool ok,HermiteGalsIdentity id,const std::string& message)
{if(!ok)throw std::invalid_argument("[Interface/HermiteGalsPhaseFlux geometry="+std::to_string(id.value())+"] "+message);}
struct Sum {long double value=0,correction=0;void add(long double x){const auto a=x-correction,b=value+a;correction=(b-value)-a;value=b;}};
NativeLeafDomain domainFor(const QuadtreeMesh& mesh,const HermiteGalsState& state)
{
    state.validateMesh(mesh);const auto& p=state.primary();const auto id=state.identity();
    require(!p.geometryCells.empty() && p.geometryCells.size()==p.phaseVolumes.size(),id,"domain: inventory requires complete scalar support");
    for(std::size_t i=0;i<p.geometryCells.size();++i)require(p.geometryCells[i]==p.phaseVolumes[i].cell,id,"domain: inventory/scalar cell mismatch");
    const auto& scalar=state.geometry().snapshot();
    require(scalar.usesNativeLeaves() && scalar.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic
        && scalar.radialAxisSymmetry()==(p.measure==HermiteGalsMeasure::AxisymmetricVolume),id,"domain: native Q4/measure mismatch");
    FaceSegmentGraph graph;graph.build(mesh,[&](const QuadtreeMesh&,const Cell& cell){return std::binary_search(p.geometryCells.begin(),p.geometryCells.end(),cell.index());});
    auto result=NativeLeafDomain::build(mesh,graph,p.measure==HermiteGalsMeasure::AxisymmetricVolume);
    InterfacePhaseGeometry::validateLayout(mesh,result,scalar);return result;
}
std::vector<Term> velocityTerms(const std::array<double,11>& input)
{Term term;term.factor=1;std::copy(input.begin(),input.end(),term.coefficients.begin());return {term};}
std::vector<Term> densityTerms(const std::array<double,11>& input,Vector2D a,Vector2D b,bool radial)
{
    const double length=std::hypot(b.x-a.x,b.y-a.y);auto result=velocityTerms(input);
    result[0].factor=length;
    if(radial){result[0].factor*=2*std::acos(-1.);result[0].weight=a.x;
        Term linear;linear.factor=result[0].factor;linear.weight=b.x-a.x;
        for(std::size_t k=0;k<input.size();++k)linear.coefficients[k+1]=input[k];result.push_back(linear);}
    return result;
}
void validateFaces(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const std::vector<HermiteGalsFaceVelocity>& input,HermiteGalsIdentity id)
{
    std::size_t index=0;
    for(const auto face:domain.faces(mesh)){
        if(face.negativeCell<0 && face.positiveCell<0)continue;
        require(index<input.size() && input[index].nativeFace==face.nativeRow,id,"velocity: missing/duplicate/unordered active face="+std::to_string(face.nativeRow));
        const auto& f=input[index++];const bool internal=face.negativeCell>=0 && face.positiveCell>=0;
        require(internal?f.boundary==Boundary::Interior:(f.boundary==Boundary::Impermeable || f.boundary==Boundary::Open),id,"velocity: boundary kind/owner mismatch");
        const bool axis=domain.axisymmetric() && !face.horizontal && face.center.x==0;
        require(!axis || f.boundary==Boundary::Impermeable,id,"velocity: axis must be impermeable");
        for(const auto& half:f.halfPolynomial){for(double value:half)require(std::isfinite(value),id,"velocity: nonfinite coefficient");
            if(f.boundary==Boundary::Impermeable){const auto bounds=Quantum::Geometry::encloseUnitPolynomial(velocityTerms(half));
                require(bounds[0]==0 && bounds[1]==0,id,"velocity: impermeable polynomial must be identically zero");}}
    }
    require(index==input.size(),id,"velocity: foreign/inactive face entry");
}
}
struct HermiteGalsFluxVelocity::Data
{
    Data(const HermiteGalsState& state,const std::vector<HermiteGalsFaceVelocity>& input):source(state),faces(input){}
    HermiteGalsState source;std::vector<HermiteGalsFaceVelocity> faces;
};
HermiteGalsFluxVelocity::HermiteGalsFluxVelocity(const QuadtreeMesh& mesh,const HermiteGalsState& midpoint,const std::vector<HermiteGalsFaceVelocity>& faces)
{
    require(midpoint.timeLayer().position==HermiteGalsTimeLayer::Position::Midpoint,midpoint.identity(),"velocity: requires midpoint layer");
    const auto domain=domainFor(mesh,midpoint);validateFaces(mesh,domain,faces,midpoint.identity());
    data_=std::make_shared<const Data>(midpoint,faces);
}
const HermiteGalsFluxVelocity::Data& HermiteGalsFluxVelocity::checked() const
{require(static_cast<bool>(data_),{},"velocity: moved-from input");return *data_;}
HermiteGalsIdentity HermiteGalsFluxVelocity::geometryIdentity() const {return checked().source.identity();}
HermiteGalsTimeLayer HermiteGalsFluxVelocity::timeLayer() const {return checked().source.timeLayer();}
const std::vector<HermiteGalsFaceVelocity>& HermiteGalsFluxVelocity::faces() const {return checked().faces;}
void HermiteGalsFluxVelocity::validateMesh(const QuadtreeMesh& mesh) const {checked().source.validateMesh(mesh);}
struct HermiteGalsPhaseFluxResult::Data
{
    Data(const HermiteGalsState& old,const HermiteGalsState& mid):previous(old),midpoint(mid){}
    HermiteGalsState previous,midpoint;
    std::vector<HermiteGalsPhaseFaceFlux> faces;
    std::vector<HermiteGalsCellVolume> volumes;
    HermiteGalsPhaseFluxReport report;
};
HermiteGalsPhaseFluxResult::HermiteGalsPhaseFluxResult(std::shared_ptr<const Data> data):data_(std::move(data)){}
const HermiteGalsPhaseFluxResult::Data& HermiteGalsPhaseFluxResult::checked() const
{require(static_cast<bool>(data_),{},"result: moved-from handle");return *data_;}
HermiteGalsIdentity HermiteGalsPhaseFluxResult::previousIdentity() const {return checked().previous.identity();}
void HermiteGalsPhaseFluxResult::validateMesh(const QuadtreeMesh& mesh) const {checked().midpoint.validateMesh(mesh);}
HermiteGalsIdentity HermiteGalsPhaseFluxResult::geometryIdentity() const {return checked().midpoint.identity();}
HermiteGalsTimeLayer HermiteGalsPhaseFluxResult::timeLayer() const {return checked().midpoint.timeLayer();}
HermiteGalsGeometryView HermiteGalsPhaseFluxResult::geometry() const {return checked().midpoint.geometry();}
const std::vector<HermiteGalsPhaseFaceFlux>& HermiteGalsPhaseFluxResult::faces() const {return checked().faces;}
const std::vector<HermiteGalsCellVolume>& HermiteGalsPhaseFluxResult::phaseVolumes() const {return checked().volumes;}
const HermiteGalsPhaseFluxReport& HermiteGalsPhaseFluxResult::report() const {return checked().report;}
HermiteGalsPhaseFluxResult advanceHermiteGalsPhaseVolumes(const QuadtreeMesh& mesh,const HermiteGalsState& previous,
    const HermiteGalsState& midpoint,const HermiteGalsFluxVelocity& velocity,const HermiteGalsPhaseFluxConfig& config)
{
    const auto id=midpoint.identity();
    try {
        previous.validateMesh(mesh);velocity.validateMesh(mesh);
        const auto t=midpoint.timeLayer(),old=previous.timeLayer(),vt=velocity.timeLayer();
        require(t.position==HermiteGalsTimeLayer::Position::Midpoint && t.step>0 && old.step==t.step-1
            && old.position!=HermiteGalsTimeLayer::Position::Midpoint && old.time()==t.begin,id,"advance: incompatible previous/midpoint time layer");
        require(velocity.geometryIdentity()==id && vt.step==t.step && vt.begin==t.begin && vt.end==t.end && vt.position==t.position,id,"advance: foreign velocity geometry/time");
        require(previous.primary().measure==midpoint.primary().measure && previous.primary().geometryCells==midpoint.primary().geometryCells
            && previous.primary().phaseVolumes.size()==midpoint.primary().phaseVolumes.size(),id,"advance: previous/midpoint inventory layout differs");
        for(double tolerance:{config.relativeIntegrationTolerance,config.relativeDivergenceTolerance,config.relativeConservationTolerance})
            require(std::isfinite(tolerance) && tolerance>0 && tolerance<1,id,"advance: invalid relative tolerance");
        const auto domain=domainFor(mesh,midpoint);const auto& cells=domain.cellIndexByRow();
        const long double dt=t.end-t.begin;const auto count=static_cast<std::size_t>(mesh.cellIndexNum());
        const Detail::HermiteInventoryAssembly inventory(count,cells,domain.cellVolumes(),previous.primary().phaseVolumes,id);
        const auto& capacity=inventory.capacity;
        auto owner=std::make_shared<const HermiteGalsGeometryView>(midpoint.geometry());
        const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar(owner,&owner->snapshot());
        InterfacePhaseGeometry::MappedFacePartitions partitions(mesh,domain,scalar);
        auto output=std::make_shared<HermiteGalsPhaseFluxResult::Data>(previous,midpoint);
        std::size_t index=0,rangeCount=0,internalCount=0,counterflowCount=0,subsetReconciliations=0,directionalIdentities=0;Sum errorEstimate;
        long double maxDirectionalAdjustment=0;
        for(const auto face:domain.faces(mesh)){
            if(face.negativeCell<0 && face.positiveCell<0)continue;
            const auto& input=velocity.faces()[index++];HermiteGalsPhaseFaceFlux flux;
            flux.nativeFace=face.nativeRow;flux.negativeCell=face.negativeCell;flux.positiveCell=face.positiveCell;
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
            for(int half=0;half<2;++half){
                auto terms=densityTerms(input.halfPolynomial[half],half==0?a:face.center,half==0?face.center:b,domain.axisymmetric());
                auto negative=terms;for(auto& term:negative)term.factor=-term.factor;
                const auto bounds=Quantum::Geometry::encloseUnitPolynomial(terms);
                const auto scale=std::max(std::abs(bounds[0]),std::abs(bounds[1]));
                const long double tolerance=config.relativeIntegrationTolerance*scale;
                const auto integrate=[&](const std::vector<Term>& polynomial,bool fromEnd,long double offset,long double length,long double* error=nullptr){
                    const auto value=Quantum::Geometry::integratePositiveUnitPolynomialRange(polynomial,0,1,fromEnd,offset,length,tolerance*length);
                    require(std::isfinite(value.value) && value.value>=0,id,"integrate: invalid gross flux face="+std::to_string(face.nativeRow));
                    if(error)*error+=value.absoluteError;
                    flux.integrationErrorEstimate+=value.absoluteError;return value.value;};
                long double fullPositiveError=0,fullNegativeError=0,liquidPositiveError=0,liquidNegativeError=0;
                const auto fullPositive=integrate(terms,false,0,1,&fullPositiveError),fullNegative=integrate(negative,false,0,1,&fullNegativeError);
                flux.totalPositiveRate+=fullPositive;flux.totalNegativeRate+=fullNegative;
                const auto& cut=partitions.halfFace(face.nativeRow,half);
                const bool allLiquid=std::all_of(cut.ranges().begin(),cut.ranges().end(),[](const auto& range){return range.negative;});
                long double liquidPositive=0,liquidNegative=0;
                // A certified full-liquid half may contain isolated tangent
                // roots. Use the identical total primitive, preserving the
                // pure-phase identity without independently rounding sub-sums.
                if(allLiquid){liquidPositive=fullPositive;liquidNegative=fullNegative;rangeCount+=cut.ranges().size();}
                else for(const auto& range:cut.ranges())if(range.negative){++rangeCount;
                    liquidPositive+=integrate(terms,range.fromEnd,range.offset,range.length,&liquidPositiveError);
                    liquidNegative+=integrate(negative,range.fromEnd,range.offset,range.length,&liquidNegativeError);}
                flux.integratedLiquidPositiveRate+=liquidPositive;flux.integratedLiquidNegativeRate+=liquidNegative;
                const auto consistentSubset=[&](long double part,long double partError,long double whole,long double wholeError,const std::vector<Term>& polynomial){
                    if(part<=whole)return part;
                    const auto excess=part-whole,gross=fullPositive+fullNegative;
                    require(std::isfinite(excess) && excess<=partError+wholeError && excess<=roundoffBudget*gross,id,
                        "integrate: subset excess exceeds quadrature/roundoff budget face="+std::to_string(face.nativeRow));
                    // Separate positive-part estimates need not be monotone
                    // in their domains. Measure the complement and reconcile
                    // the pair to the UNCHANGED whole measure. Zero gas gives
                    // the exact transport-direction identity L=total; otherwise
                    // retain the relative nonnegative phase weights. BOTH
                    // quadrature and existing roundoff budgets must cover the
                    // entire adjustment, not merely the initial excess.
                    Sum gas;long double gasError=0;
                    for(const auto& range:cut.ranges())if(!range.negative)
                        gas.add(integrate(polynomial,range.fromEnd,range.offset,range.length,&gasError));
                    const auto partition=part+gas.value;
                    require(std::isfinite(partition) && partition>0,id,"integrate: invalid complementary gross partition");
                    const auto value=gas.value==0?whole:whole*(part/partition);
                    const auto adjustment=part-value;
                    require(std::isfinite(value) && value>=0 && value<=whole && adjustment>=0
                        && adjustment<=partError+wholeError+gasError && adjustment<=roundoffBudget*gross,id,
                        "integrate: complementary reconciliation exceeds quadrature/roundoff budget face="+std::to_string(face.nativeRow));
                    require(partError+adjustment<=tolerance,id,"integrate: consistent subset exceeds integration budget face="+std::to_string(face.nativeRow));
                    flux.integrationErrorEstimate+=adjustment;
                    ++subsetReconciliations;directionalIdentities+=gas.value==0 && gasError==0;
                    maxDirectionalAdjustment=std::max(maxDirectionalAdjustment,adjustment/gross);
                    return value;
                };
                liquidPositive=consistentSubset(liquidPositive,liquidPositiveError,fullPositive,fullPositiveError,terms);
                liquidNegative=consistentSubset(liquidNegative,liquidNegativeError,fullNegative,fullNegativeError,negative);
                // Add the two nonnegative half-face pairs in the SAME order,
                // preserving each proven subset inequality on the full face.
                flux.liquidPositiveRate+=liquidPositive;flux.liquidNegativeRate+=liquidNegative;
            }
            require(std::isfinite(flux.totalPositiveRate) && std::isfinite(flux.totalNegativeRate)
                && flux.liquidPositiveRate<=flux.totalPositiveRate && flux.liquidNegativeRate<=flux.totalNegativeRate,id,
                "integrate: gross phase flux outside total bounds face="+std::to_string(face.nativeRow));
            const auto localCapacity=face.negativeCell<0?capacity[face.positiveCell]:face.positiveCell<0?capacity[face.negativeCell]:
                std::min(capacity[face.negativeCell],capacity[face.positiveCell]);
            const auto subsetAdjustment=std::max(flux.integratedLiquidPositiveRate-flux.liquidPositiveRate,
                flux.integratedLiquidNegativeRate-flux.liquidNegativeRate);
            require(subsetAdjustment>=0 && dt*subsetAdjustment<=roundoffBudget*localCapacity,id,
                "integrate: subset reconciliation exceeds existing inventory roundoff budget face="+std::to_string(face.nativeRow));
            flux.rawLiquidPositiveRate=flux.liquidPositiveRate;flux.rawLiquidNegativeRate=flux.liquidNegativeRate;
            flux.totalRate=flux.totalPositiveRate-flux.totalNegativeRate;
            internalCount+=face.negativeCell>=0 && face.positiveCell>=0;
            counterflowCount+=flux.totalPositiveRate>0 && flux.totalNegativeRate>0;
            errorEstimate.add(flux.integrationErrorEstimate);output->faces.push_back(flux);
        }
        require(std::isfinite(errorEstimate.value),id,"integrate: nonfinite accumulated flux estimate");
        std::vector<Detail::InventoryFace> working;working.reserve(output->faces.size());
        for(const auto& face:output->faces)working.push_back({face.negativeCell,face.positiveCell,
            face.totalPositiveRate,face.totalNegativeRate,face.liquidPositiveRate,face.liquidNegativeRate,
            face.rawLiquidPositiveRate,face.rawLiquidNegativeRate});
        const auto assembled=inventory.finish(working,dt,config.relativeDivergenceTolerance,config.relativeConservationTolerance,output->volumes);
        for(std::size_t k=0;k<working.size();++k){const auto& from=working[k];auto& face=output->faces[k];
            face.liquidPositiveRate=from.liquidPositive;face.liquidNegativeRate=from.liquidNegative;
            face.liquidRate=from.liquidSigned;face.liquidVolume=from.liquidVolume;}
        auto& report=output->report;report.activeFaces=output->faces.size();report.internalFaces=internalCount;report.counterflowFaces=counterflowCount;
        report.phaseRanges=rangeCount;report.roundoffRepairSweeps=assembled.roundoffRepairSweeps;report.repairedFaces=assembled.repairedFaces;
        report.liquidBefore=assembled.liquidBefore;report.liquidAfter=assembled.liquidAfter;
        report.boundaryOutward=assembled.boundaryOutward;report.closure=assembled.closure;
        report.maxOutflowCfl=assembled.maxOutflowCfl;report.maxRelativeDivergence=assembled.maxRelativeDivergence;
        report.maxRelativeFluxCorrection=assembled.maxRelativeFluxCorrection;
        report.subsetReconciliationCount=subsetReconciliations;report.directionalIdentityCount=directionalIdentities;report.maxRelativeSubsetAdjustment=maxDirectionalAdjustment;
        report.integrationErrorEstimate=errorEstimate.value;
        return HermiteGalsPhaseFluxResult(std::move(output));
    }catch(const std::exception& e){throw std::invalid_argument("[Interface/HermiteGalsPhaseFlux geometry="+std::to_string(id.value())+"] "+e.what());}
}
}
