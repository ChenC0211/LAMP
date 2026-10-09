#include "Modules/Interface/Compatibility/InterfaceResolutionPlan.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceParametricResolution.h"
#include "Modules/Interface/Compatibility/InterfaceCurvature.h"
#include "Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Modules/Interface/Compatibility/InterfaceNeckCap.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace InterfaceResolutionPlan
{
static Result evaluateLayer(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceState& state,
    const Options& o,int cellBudget,int minimumInterfaceLevel=0,std::vector<int>* retentionFloor=nullptr,double hysteresis=1)
{
    const auto require=[](bool ok,const char* why){if(!ok) throw std::invalid_argument(why);};
    state.validateMesh(mesh);domain.validateMeshIdentity(mesh);
    const bool parametric=state.hasParametricGeometry();
    require(parametric||o.parametricNeckRadiusCells==0,"parametric terminal policy requires parametric geometry");
    require(parametric || (state.hasTopology() && state.hasConservativePhase()),"resolution planning requires geometry and independent inventory");
    require(parametric || state.pinnedCoordinates().empty(),"resolution planning requires native topology contact identities");
    require(o.maximumLevel>=0 && o.maximumLevel<=24 && std::isfinite(o.maxCurvatureCells)
        && o.maxCurvatureCells>0 && o.maxCurvatureCells<=.5
        && std::isfinite(o.minimumOpposedDistanceCells) && o.minimumOpposedDistanceCells>=2
        && std::isfinite(o.maximumGeometryFractionDifference) && o.maximumGeometryFractionDifference>0
        && o.maximumGeometryFractionDifference<=1 && std::isfinite(o.distanceBandCells) && o.distanceBandCells>0
        && std::isfinite(o.maximumTravel) && o.maximumTravel>=0
        && std::isfinite(o.reserveBandCells) && o.reserveBandCells>=0 && o.reserveBandCells<=64
        && cellBudget>=0,"invalid conservative resolution policy");
    Result result;result.requiredInterfaceLevel=minimumInterfaceLevel;result.localInterfaceLevels=parametric;
    result.targetLevels.resize(mesh.cellIndexNum());std::vector<unsigned char> seed(mesh.cellIndexNum(),0);
    std::vector<int> required(mesh.cellIndexNum(),0);
    for(int c=0;c<mesh.cellIndexNum();++c)result.targetLevels[c]=retentionFloor?0:mesh.cells(mesh.cellIndexToLabel[c]).level();
    const auto request=[&](int c,double ratio) {
        int level=mesh.cells(mesh.cellIndexToLabel[c]).level();
        // A fine retained leaf must not seed an ever-widening fine band.
        // Recover the required physical scale relative to the root lattice,
        // even when the current leaf already exceeds that requirement.
        if(parametric){ratio=std::ldexp(ratio,level);level=0;}
        // One beyond the cap retains evidence of an unmet finite request;
        // do not loop to an unbounded desired level or erase that request.
        while(ratio>1 && level<=o.maximumLevel) {ratio*=.5;++level;}
        result.requiredInterfaceLevel=std::max(result.requiredInterfaceLevel,level);
        required[c]=std::max(required[c],std::min(level,o.maximumLevel));
    };
    if(parametric)parametricIndicators(mesh,domain,state,o,result,seed,request,hysteresis);
    else {
    const auto& phase=state.conservativePhase();
    for(int c:domain.cellIndexByRow()) {
        const int level=mesh.cells(mesh.cellIndexToLabel[c]).level();result.targetLevels[c]=level;
        require(level<=o.maximumLevel,"accepted mesh exceeds resolution policy maximum level");
        const auto f=phase.liquidVolumes()[c]/phase.capacities()[c];
        if(std::min(f,1-f)>Quantum::Geometry::PlanarCellCut::relativeVolumeResolution) seed[c]=1;
    }
    InterfaceCurvature curvature(mesh,domain,state);
    const auto terminal=o.terminalAxialSection?InterfaceNeckCap::terminalResolutionCells(mesh,domain,state,
        *o.terminalAxialSection,o.maximumLevel,o.maxCurvatureCells):std::vector<int>{};
    for(std::size_t s=0;s<state.topology().surfaces().size();++s) {
        const auto& surface=state.topology().surfaces()[s];
        const auto& chain=surface.geometry.nativeChain();
        for(std::size_t g=0;g<chain.graphs().size();++g) {
            const int c=chain.graphs()[g].nativeCellIndex;seed[c]=1;
            const auto size=mesh.cells(mesh.cellIndexToLabel[c]).size();const double h=std::max(size.x,size.y);
            double magnitude=0,azimuthalMagnitude=0,rawMagnitude=0;bool unavailable=false,needsQuality=false,needsGap=false;
            // Use the SAME height-first/local-fit curvature as physical loads,
            // including segment endpoints. Keep raw Hessian spikes observable;
            // they supply a conservative scale estimate only where the physical
            // provider is unavailable, never overriding a valid provider.
            for(double fraction:{0.,.5,1.}) {
                const auto frame=chain.frameAt(static_cast<int>(g),fraction);
                const double k=frame.curvature;
                require(std::isfinite(k),"nonfinite native curvature resolution indicator");
                // The RZ axis has the regular even limit, not 1/r. Away from
                // it, n_r/r resolves a straight narrow cylinder even though
                // meridional curvature vanishes. Do not sum principal parts:
                // their cancellation at a saddle does not resolve a neck.
                const double rawAzimuthal=!domain.axisymmetric()?0:
                    (frame.point.x==0?k:frame.normal.x/frame.point.x);
                require(std::isfinite(rawAzimuthal),"nonfinite native azimuthal resolution indicator");
                rawMagnitude=std::max(rawMagnitude,std::max(std::abs(k),std::abs(rawAzimuthal))*h);
                const auto sample=curvature.evaluate(s,static_cast<int>(g),fraction);const auto& fit=sample.fit;
                if(sample.valid()) {
                    require(std::isfinite(sample.meridional) && std::isfinite(sample.azimuthal),
                        "nonfinite physical curvature resolution indicator");
                    const double azimuthal=domain.axisymmetric()?std::abs(sample.azimuthal)*h:0;
                    azimuthalMagnitude=std::max(azimuthalMagnitude,azimuthal);
                    magnitude=std::max(magnitude,std::max(std::abs(sample.meridional)*h,azimuthal));
                } else {
                    unavailable=true;
                    azimuthalMagnitude=std::max(azimuthalMagnitude,std::abs(rawAzimuthal)*h);
                    magnitude=std::max(magnitude,std::max(std::abs(k),std::abs(rawAzimuthal))*h);
                    const bool supportOnly=sample.height.status==InterfaceHeightCurvature::Status::MixedResolution
                        || fit.refinementCause==InterfaceCurvature::RefinementCause::MixedResolution;
                    needsQuality|=!supportOnly;
                }
                if(fit.refinementCause==InterfaceCurvature::RefinementCause::CurvatureRadius
                    && fit.resolutionScale>0 && std::isfinite(fit.scaledMeridionalMagnitude))
                    magnitude=std::max(magnitude,fit.scaledMeridionalMagnitude*h/fit.resolutionScale);
                if(fit.refinementCause==InterfaceCurvature::RefinementCause::OpposedArc
                    && std::isfinite(fit.opposedSampleDistance) && fit.opposedSampleDistance>=0) {
                    const double ratio=fit.opposedSampleDistance>0?o.minimumOpposedDistanceCells*h/fit.opposedSampleDistance:
                        std::numeric_limits<double>::infinity();
                    if(ratio>1) {needsGap=true;request(c,ratio);}
                }
            }
            result.maxRawCurvatureCells=std::max(result.maxRawCurvatureCells,rawMagnitude);
            if(!unavailable && magnitude<=o.maxCurvatureCells && rawMagnitude>o.maxCurvatureCells)
                ++result.rawOnlyCurvatureCells;
            result.maxObservedCurvatureCells=std::max(result.maxObservedCurvatureCells,magnitude);
            result.maxObservedAzimuthalCells=std::max(result.maxObservedAzimuthalCells,azimuthalMagnitude);
            if(azimuthalMagnitude>o.maxCurvatureCells)++result.azimuthalRequests;
            if(magnitude>o.maxCurvatureCells) {
                ++result.curvatureRequests;
                const int level=mesh.cells(mesh.cellIndexToLabel[c]).level();
                const double atFinest=std::ldexp(magnitude,level-o.maximumLevel);
                if(atFinest>o.maxCurvatureCells && std::binary_search(terminal.begin(),terminal.end(),c)) {
                    ++result.terminalCurvatureRequests;
                    result.maxTerminalCurvatureCells=std::max(result.maxTerminalCurvatureCells,atFinest);
                    result.requiredInterfaceLevel=std::max(result.requiredInterfaceLevel,o.maximumLevel);
                } else request(c,magnitude/o.maxCurvatureCells);
            }
            if(needsGap)++result.gapRequests;
            if(needsQuality) {++result.qualityRequests;request(c,2);}
        }
    }
    if(state.hasConservativeGeometry()) for(int c:domain.cellIndexByRow()) {
        const auto difference=std::abs(state.conservativeGeometry().geometricVolumeByCell[c]-phase.liquidVolumes()[c])/phase.capacities()[c];
        if(difference>o.maximumGeometryFractionDifference) {
            seed[c]=1;++result.inventoryGeometryRequests;request(c,2);
        }
    }
    }
    for(int c:domain.cellIndexByRow()) if(seed[c]) {
        ++result.interfaceCells;result.requiredInterfaceLevel=std::max(result.requiredInterfaceLevel,
            static_cast<int>(mesh.cells(mesh.cellIndexToLabel[c]).level()));
    }
    require(result.interfaceCells>0||(parametric&&state.parametricGeometry().componentCount()==0),"resolution planning found no resolved interface support");
    result.saturated=result.requiredInterfaceLevel>o.maximumLevel;
    result.interfaceLevel=std::min(result.requiredInterfaceLevel,o.maximumLevel);
    const CartesianLeafLocator locator(mesh);
    const double infinity=std::numeric_limits<double>::infinity();
    const auto extendBand=[&](double reserve) {for(int c:domain.cellIndexByRow()) if(seed[c]) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        // A disappearing child cannot define the next review's band. Use
        // its prospective parent footprint, including the retained reserve.
        const auto& footprint=retentionFloor&&cell.level()>0?mesh.cells(cell.parent):cell;
        const auto a=retentionFloor?footprint.centerCoord()-footprint.size()*.5:mesh.verts(cell.cellToVert[0]).coord();
        const auto b=retentionFloor?footprint.centerCoord()+footprint.size()*.5:mesh.verts(cell.cellToVert[2]).coord();
        const int targetLevel=parametric?required[c]:result.interfaceLevel;
        const auto size=cell.size()*std::ldexp(1.,static_cast<int>(cell.level())-targetLevel);
        // Radius-six transport columns, two extra cells for nodal blending,
        // and the scalar plateau/next-review motion envelope share ONE band.
        const double cells=(parametric?o.distanceBandCells:std::max(8.,o.distanceBandCells+2))+reserve;
        const Vector2D pad{cells*size.x+o.maximumTravel,cells*size.y+o.maximumTravel};
        for(int target:locator.intersectingBox(
            {std::nextafter(a.x-pad.x,-infinity),std::nextafter(a.y-pad.y,-infinity)},
            {std::nextafter(b.x+pad.x,infinity),std::nextafter(b.y+pad.y,infinity)}))
            result.targetLevels[target]=std::max(result.targetLevels[target],targetLevel);
    }};
    // Close targets over actual face neighbors before asking Core to adapt.
    // This is deliberately conservative on coarse source leaves; Core still
    // owns allocation, actual refinement and validation of the candidate.
    const auto balance=[&] {for(int sweep=0;sweep<=o.maximumLevel;++sweep) {
        bool changed=false;
        for(int c=0;c<mesh.cellIndexNum();++c)for(int side=0;side<4;++side)for(int part=0;part<2;++part){
            const int neighbor=mesh.cells(mesh.cellIndexToLabel[c]).neiCell[side][part];if(neighbor<0)continue;
            auto& a=result.targetLevels[c];auto& b=result.targetLevels[mesh.cells(neighbor).index()];
            if(a+1<b) {a=b-1;changed=true;}else if(b+1<a) {b=a-1;changed=true;}}
        if(!changed) break;
        require(sweep<o.maximumLevel,"resolution target balance did not converge");
    }};
    const auto leafCount=[&] {
        std::size_t count=0;
        for(int c=0;c<mesh.cellIndexNum();++c) {
            const int delta=result.targetLevels[c]-static_cast<int>(mesh.cells(mesh.cellIndexToLabel[c]).level());
            require(delta>=0 && delta<=24,"invalid refinement level difference");
            const std::size_t leaves=std::size_t(1)<<(2*delta);
            require(leaves<=std::numeric_limits<std::size_t>::max()-count,"refinement leaf count overflow");
            count+=leaves;
        }
        return count;
    };
    if(retentionFloor){
        extendBand(o.reserveBandCells);balance();*retentionFloor=result.targetLevels;
        return result;
    }
    extendBand(0);balance();
    result.requiredLeafCount=result.plannedLeafCount=leafCount();
    // Never chase a moving reserve envelope on every step. The unchanged
    // mandatory band, including physical travel and all stencil supports,
    // is the trigger. Reserve only amortizes an already necessary rebuild.
    if(!result.saturated && o.reserveBandCells>0
        && result.requiredLeafCount>static_cast<std::size_t>(mesh.cellIndexNum())) {
        auto required=result.targetLevels;
        extendBand(o.reserveBandCells);balance();
        result.plannedLeafCount=leafCount();
        if(cellBudget>0 && result.plannedLeafCount>static_cast<std::size_t>(cellBudget)) {
            result.targetLevels=std::move(required);result.plannedLeafCount=result.requiredLeafCount;
            result.reserveBudgetLimited=true;
        } else result.reserveApplied=result.plannedLeafCount>result.requiredLeafCount;
    }
    for(int c=0;c<mesh.cellIndexNum();++c) if(result.targetLevels[c]>static_cast<int>(mesh.cells(mesh.cellIndexToLabel[c]).level())) ++result.refinedCells;
    return result;
}
Result evaluate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceState& state,const Options& o,int budget)
{ return evaluateLayer(mesh,domain,state,o,budget); }

Result evaluate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceState& current,
    const InterfaceState& previous,const Options& o,int budget)
{
    auto mandatory=o;mandatory.reserveBandCells=0;
    if(current.hasParametricGeometry()!=previous.hasParametricGeometry())
        throw std::invalid_argument("resolution histories use different geometry representations");
    auto a=evaluateLayer(mesh,domain,current,mandatory,budget),b=evaluateLayer(mesh,domain,previous,mandatory,budget);
    const int common=std::max(a.interfaceLevel,b.interfaceLevel);
    if(!a.localInterfaceLevels){
        if(a.interfaceLevel<common)a=evaluateLayer(mesh,domain,current,mandatory,budget,common);
        if(b.interfaceLevel<common)b=evaluateLayer(mesh,domain,previous,mandatory,budget,common);
    }
    const auto combine=[&](const Result& x,const Result& y) {
        Result result=x;
        result.interfaceLevel=common;result.requiredInterfaceLevel=std::max(x.requiredInterfaceLevel,y.requiredInterfaceLevel);
        result.saturated=x.saturated || y.saturated;
        result.interfaceCells+=y.interfaceCells;result.curvatureRequests+=y.curvatureRequests;
        result.gapRequests+=y.gapRequests;result.qualityRequests+=y.qualityRequests;
        result.inventoryGeometryRequests+=y.inventoryGeometryRequests;result.azimuthalRequests+=y.azimuthalRequests;
        result.rawOnlyCurvatureCells+=y.rawOnlyCurvatureCells;
        result.terminalCurvatureRequests+=y.terminalCurvatureRequests;
        result.terminalGapRequests+=y.terminalGapRequests;
        result.terminalGeometricVolume=std::max(x.terminalGeometricVolume,y.terminalGeometricVolume);
        result.terminalRanges.insert(result.terminalRanges.end(),y.terminalRanges.begin(),y.terminalRanges.end());
        result.unresolved.insert(result.unresolved.end(),y.unresolved.begin(),y.unresolved.end());
        result.maxTerminalCurvatureCells=std::max(x.maxTerminalCurvatureCells,y.maxTerminalCurvatureCells);
        result.maxObservedCurvatureCells=std::max(x.maxObservedCurvatureCells,y.maxObservedCurvatureCells);
        result.maxObservedAzimuthalCells=std::max(x.maxObservedAzimuthalCells,y.maxObservedAzimuthalCells);
        result.maxRawCurvatureCells=std::max(x.maxRawCurvatureCells,y.maxRawCurvatureCells);
        result.refinedCells=0;result.plannedLeafCount=0;
        for(int c=0;c<mesh.cellIndexNum();++c) {
            result.targetLevels[c]=std::max(x.targetLevels[c],y.targetLevels[c]);
            const int delta=result.targetLevels[c]-static_cast<int>(mesh.cells(mesh.cellIndexToLabel[c]).level());
            const auto leaves=std::size_t(1)<<(2*delta);
            if(leaves>std::numeric_limits<std::size_t>::max()-result.plannedLeafCount)
                throw std::length_error("history refinement leaf count overflow");
            result.plannedLeafCount+=leaves;result.refinedCells+=delta>0;
        }
        result.requiredLeafCount=result.plannedLeafCount;
        return result;
    };
    auto result=combine(a,b);
    if(result.changesMesh() && !result.saturated && o.reserveBandCells>0) {
        const int minimum=a.localInterfaceLevels?0:common;
        auto reserve=combine(evaluateLayer(mesh,domain,current,o,0,minimum),evaluateLayer(mesh,domain,previous,o,0,minimum));
        if(budget>0 && reserve.plannedLeafCount>static_cast<std::size_t>(budget))result.reserveBudgetLimited=true;
        else {
            reserve.requiredLeafCount=result.requiredLeafCount;
            reserve.reserveApplied=reserve.plannedLeafCount>result.requiredLeafCount;
            result=std::move(reserve);
        }
    }
    return result;
}

CoarseningPlan planCoarsening(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceState& state,
    const Options& options,double hysteresis,const std::vector<Vector2D>& protectedPoints)
{
    if(!state.hasParametricGeometry()||!std::isfinite(hysteresis)||hysteresis<=0||hysteresis>=1)
        throw std::invalid_argument("coarsening requires parametric geometry and hysteresis in (0,1)");
    for(const auto& p:protectedPoints)if(!std::isfinite(p.x)||!std::isfinite(p.y))
        throw std::invalid_argument("nonfinite protected coarsening point");
    std::vector<int> floor;evaluateLayer(mesh,domain,state,options,0,0,&floor,hysteresis);
    CoarseningPlan result;result.plannedLeafCount=mesh.cellIndexNum();result.targetLevels.resize(mesh.cellIndexNum());
    for(int c=0;c<mesh.cellIndexNum();++c)result.targetLevels[c]=mesh.cells(mesh.cellIndexToLabel[c]).level();
    for(int c=0;c<mesh.cellIndexNum();++c){
        const auto& child=mesh.cells(mesh.cellIndexToLabel[c]);if(child.level()==0)continue;
        const auto& parent=mesh.cells(child.parent);if(child.label()!=parent.son[0])continue;
        bool eligible=true;
        for(int label:parent.son){
            if(label<0){eligible=false;break;}
            const auto& sibling=mesh.cells(label);const int index=sibling.index();
            if(index<0||index>=mesh.cellIndexNum()||mesh.cellIndexToLabel[index]!=label
                ||sibling.level()!=child.level()||floor[index]>parent.level()){eligible=false;break;}
            // Match Core's pre-change 2:1 eligibility. Families adjacent to
            // finer leaves can become eligible on a later physical review.
            for(int side=0;side<4;side++)for(int part=0;part<2;part++){
                const int neighbor=sibling.neiCell[side][part];
                if(neighbor>=0&&mesh.cells(neighbor).level()>sibling.level())eligible=false;
            }
        }
        if(!eligible)continue;
        const auto a=parent.centerCoord()-parent.size()*.5,b=parent.centerCoord()+parent.size()*.5;
        const double tolerance=64*std::numeric_limits<double>::epsilon()
            *std::max({std::abs(a.x),std::abs(a.y),std::abs(b.x),std::abs(b.y),parent.size().x,parent.size().y});
        for(const auto& point:protectedPoints){
            const bool inside=point.x>=a.x-tolerance&&point.x<=b.x+tolerance&&point.y>=a.y-tolerance&&point.y<=b.y+tolerance;
            const bool corner=(std::abs(point.x-a.x)<=tolerance||std::abs(point.x-b.x)<=tolerance)
                &&(std::abs(point.y-a.y)<=tolerance||std::abs(point.y-b.y)<=tolerance);
            if(inside&&!corner)eligible=false;
        }
        if(!eligible)continue;
        for(int label:parent.son)result.targetLevels[mesh.cells(label).index()]=parent.level();
        ++result.families;result.plannedLeafCount-=3;
    }
    return result;
}

std::unique_ptr<QuadtreeMesh> coarsen(const QuadtreeMesh& mesh,const CoarseningPlan& plan)
{
    if(plan.targetLevels.size()!=static_cast<std::size_t>(mesh.cellIndexNum())
        ||plan.plannedLeafCount+3*plan.families!=static_cast<std::size_t>(mesh.cellIndexNum()))
        throw std::invalid_argument("coarsening plan does not belong to this mesh");
    auto candidate=mesh.clone();auto& target=*candidate;
    CellField<double> levels(target,0.),indicator(target,0.);
    for(int c=0;c<mesh.cellIndexNum();c++){
        const int delta=mesh.cells(mesh.cellIndexToLabel[c]).level()-plan.targetLevels[c];
        if(delta<0||delta>1||plan.targetLevels[c]<0)throw std::invalid_argument("coarsening must remove at most one generation");
        levels.value[c]=plan.targetLevels[c];
    }
    std::vector<VertField<double>> unused;std::vector<VertField<double>*> vertices;
    std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;std::vector<CellField<double>*> cells{&levels};
    if(plan.families)target.adapt(unused,indicator,vertices,horizontal,vertical,cells,&levels);
    if(static_cast<std::size_t>(target.cellIndexNum())!=plan.plannedLeafCount)
        throw std::runtime_error("Core did not realize the eligible sibling coarsening plan");
    const CartesianLeafLocator source(mesh);
    for(int c=0;c<target.cellIndexNum();c++){
        const auto& cell=target.cells(target.cellIndexToLabel[c]);const auto a=cell.centerCoord()-cell.size()*.5,b=cell.centerCoord()+cell.size()*.5;
        bool covered=false;
        for(int from:source.intersectingBox(a,b)){
            const auto& old=mesh.cells(mesh.cellIndexToLabel[from]);const auto lo=old.centerCoord()-old.size()*.5,hi=old.centerCoord()+old.size()*.5;
            if(std::min(b.x,hi.x)<=std::max(a.x,lo.x)||std::min(b.y,hi.y)<=std::max(a.y,lo.y))continue;
            if(cell.level()!=plan.targetLevels[from])throw std::runtime_error("coarsened leaf disagrees with source-family target");covered=true;
        }
        if(!covered)throw std::runtime_error("coarsened leaf has no physical source overlap");
    }
    return candidate;
}

std::unique_ptr<QuadtreeMesh> refine(const QuadtreeMesh& mesh,const Result& plan,int budget)
{
    if(plan.saturated || budget<=0 || plan.targetLevels.size()!=static_cast<std::size_t>(mesh.cellIndexNum()))
        throw std::invalid_argument("refinement requires a resolved native plan and positive cell budget");
    std::size_t requested=0;int maximumLevel=0;
    for(int c=0;c<mesh.cellIndexNum();++c){const int level=plan.targetLevels[c],delta=level-mesh.cells(mesh.cellIndexToLabel[c]).level();
        if(delta<0 || delta>8)throw std::invalid_argument("native plan requires bounded refinement-only targets");
        maximumLevel=std::max(maximumLevel,level);requested+=std::size_t(1)<<(2*delta);
        if(requested>static_cast<std::size_t>(budget))throw std::length_error("requested native leaf count exceeds cellBudget");}
    auto candidate=mesh.clone();auto& target=*candidate;const CartesianLeafLocator oldLeaves(mesh);
    std::vector<VertField<double>> unusedLevelSets;std::vector<VertField<double>*> vertices;
    std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
    bool complete=false;
    for(int sweep=0;sweep<=maximumLevel+1;++sweep){
        CellField<double> levels(target,0.),indicator(target,0.);complete=true;
        for(int c=0;c<target.cellIndexNum();++c){const auto& cell=target.cells(target.cellIndexToLabel[c]);const auto owner=oldLeaves.incident(cell.centerCoord());
            if(owner.count!=1)throw std::runtime_error("adaptive target has ambiguous source-leaf provenance");
            levels.value[c]=plan.targetLevels[owner.indices[0]];
            complete=complete && static_cast<int>(cell.level())==levels.value[c];}
        if(complete)break;
        if(!target.mayAdaptToTargets(levels))throw std::length_error("native adaptation could not fulfill the requested target levels");
        std::vector<CellField<double>*> cells{&levels};target.adapt(unusedLevelSets,indicator,vertices,horizontal,vertical,cells,&levels);
        if(target.cellIndexNum()>budget)throw std::length_error("native refinement exceeded cellBudget");
    }
    if(!complete)throw std::length_error("bounded native adaptation sweeps did not reach requested levels");
    return candidate;
}
}
