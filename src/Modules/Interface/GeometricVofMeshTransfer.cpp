#include "Modules/Interface/GeometricVofMeshTransfer.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <tuple>

namespace LAMP::Interface
{
namespace
{
namespace Cut=Geometry::PlanarCellCut;
namespace FV=Geometry::Detail::FiniteVolume;
using Overlap=GeometricVofTransferOverlap;
constexpr long double arithmeticTolerance=512*std::numeric_limits<long double>::epsilon();
void require(bool ok,const std::string& message)
{ if(!ok)throw std::invalid_argument("[Interface/GeometricVofMeshTransfer] "+message); }
// Literal success checks need no temporary error string. The original string
// overload constructs the identical failure only after the condition fails.
void require(bool ok,const char* message)
{if(!ok)require(false,std::string(message));}
struct Sum
{
    long double value=0,correction=0;
    void add(long double x){const auto y=x-correction,t=value+y;correction=(t-value)-y;value=t;}
};
Cut::Rectangle rectangle(const QuadtreeMesh& mesh,int ci)
{
    const auto& c=mesh.cells(mesh.cellIndexToLabel[ci]);
    return {mesh.verts(c.cellToVert[0]).coord(),mesh.verts(c.cellToVert[2]).coord()};
}
bool contains(Cut::Rectangle a,Cut::Rectangle b)
{ return a.lower.x<=b.lower.x&&a.lower.y<=b.lower.y&&a.upper.x>=b.upper.x&&a.upper.y>=b.upper.y; }
bool equal(Cut::Rectangle a,Cut::Rectangle b){return contains(a,b)&&contains(b,a);}
std::vector<Cut::AnchoredPoint> clipInFrame(Cut::Rectangle box,
    const std::vector<Cut::AnchoredPoint>& polygon,Cut::Point origin)
{
    // Core reads the original actual sites, world rectangle and common frame
    // together; do not round box-origin offsets or fold anchor+offset to Point.
    return Cut::clipAnchoredPolygonInFrame(box,polygon,origin);
}
long double phaseResidual(long double actual,long double expected,long double capacity,bool pure=false)
{
    require(std::isfinite(actual)&&std::isfinite(expected)&&std::isfinite(capacity)&&capacity>0,
        "nonfinite phase ledger");
    // Pure inventories are exact states. For mixed inventories both minority
    // phases constrain the relative transfer error; no volume/phase floor.
    const auto scale=pure?capacity:std::min(expected,capacity-expected);
    return std::abs(actual-expected)/(scale>0?scale:capacity);
}
struct TransferBuilder
{
    const QuadtreeMesh& sourceMesh;
    const GeometricVofPublication& source;
    const QuadtreeMesh& targetMesh;
    const NativeLeafDomain& targetDomain;
    GeometricVofConfig config;
    std::vector<Cut::Rectangle> sourceBoxes,targetBoxes;
    std::vector<long double> targetCapacity;
    std::vector<Overlap> overlaps;
    std::vector<std::vector<std::size_t>> bySource,byTarget;
    GeometricVofTransferReport report;
    bool rz;
    void coverage()
    {
        sourceBoxes.reserve(source.cells().size());bySource.resize(source.cells().size());
        for(const auto& c:source.cells())sourceBoxes.push_back(rectangle(sourceMesh,c.nativeCell));
        byTarget.resize(targetMesh.cellIndexNum());targetBoxes.reserve(byTarget.size());targetCapacity.reserve(byTarget.size());
        for(int ci=0;ci<targetMesh.cellIndexNum();++ci){
            const auto box=rectangle(targetMesh,ci);targetBoxes.push_back(box);
            const auto volume=Cut::capacity(box,rz);targetCapacity.push_back(volume);
            const auto row=targetDomain.cellRowByIndex()[ci];
            require(row>=0,"target domain omits native cell="+std::to_string(ci));
            require(std::abs(volume-targetDomain.cellVolumes()[row])<=config.relativeVolumeTolerance*volume,
                "target Core/native capacity mismatch cell="+std::to_string(ci));
        }
        std::set<std::pair<int,int>> pairs;
        const auto add=[&](int si,int ti,Overlap::Kind kind){
            require(pairs.emplace(si,ti).second,"duplicate physical leaf overlap");
            const auto box=kind==Overlap::Kind::Refined?targetBoxes[ti]:sourceBoxes[si];
            const auto index=overlaps.size();Overlap overlap;overlap.kind=kind;overlap.sourceCell=si;overlap.targetCell=ti;
            overlap.capacity=Cut::capacity(box,rz);overlaps.push_back(overlap);
            bySource[si].push_back(index);byTarget[ti].push_back(index);
        };
        // Nested leaves on the SAME background lattice. Each lookup walks
        // the actual quadtree; no N_source*N_target rectangle search.
        for(int ti=0;ti<targetMesh.cellIndexNum();++ti){
            const auto& cell=targetMesh.cells(targetMesh.cellIndexToLabel[ti]);
            const auto label=sourceMesh.findCellLabel(cell.centerCoord(),0);
            require(label>=0&&label<sourceMesh.cellLabelCapacity(),"target center has no actual source coverage");
            const auto si=sourceMesh.cells(label).index();
            require(si>=0&&static_cast<std::size_t>(si)<bySource.size(),"source lookup is not an actual leaf");
            if(contains(sourceBoxes[si],targetBoxes[ti]))
                add(si,ti,equal(sourceBoxes[si],targetBoxes[ti])?Overlap::Kind::Unchanged:Overlap::Kind::Refined);
            else require(contains(targetBoxes[ti],sourceBoxes[si]),"non-nested source/target leaves");
        }
        for(int si=0;si<sourceMesh.cellIndexNum();++si){
            const auto& cell=sourceMesh.cells(sourceMesh.cellIndexToLabel[si]);
            const auto label=targetMesh.findCellLabel(cell.centerCoord(),0);
            require(label>=0&&label<targetMesh.cellLabelCapacity(),"source center has no actual target coverage");
            const auto ti=targetMesh.cells(label).index();
            require(ti>=0&&static_cast<std::size_t>(ti)<byTarget.size(),"target lookup is not an actual leaf");
            if(contains(targetBoxes[ti],sourceBoxes[si])&&!equal(targetBoxes[ti],sourceBoxes[si]))
                add(si,ti,Overlap::Kind::Coarsened);
            else require(contains(sourceBoxes[si],targetBoxes[ti]),"non-nested target/source leaves");
        }
        const auto check=[&](const auto& groups,const auto& capacity,const char* side){
            for(std::size_t ci=0;ci<groups.size();++ci){
                require(!groups[ci].empty(),std::string(side)+" cell has no physical overlap cell="+std::to_string(ci));
                Sum total;for(auto index:groups[ci])total.add(overlaps[index].capacity);
                const auto error=std::abs(total.value-capacity(ci))/capacity(ci);
                report.maximumRelativeCoverageResidual=std::max(report.maximumRelativeCoverageResidual,error);
                require(error<=arithmeticTolerance,std::string(side)+" incomplete/duplicate physical coverage cell="+std::to_string(ci));
            }
        };
        check(bySource,[&](std::size_t ci){return source.cells()[ci].capacity;},"source");
        check(byTarget,[&](std::size_t ci){return targetCapacity[ci];},"target");
    }
    void donorInventories()
    {
        const auto origin=source.frameOrigin();
        for(std::size_t si=0;si<bySource.size();++si){
            const auto& cell=source.cells()[si];const auto alpha=source.alpha()[si];
            const bool refined=overlaps[bySource[si].front()].kind==Overlap::Kind::Refined;
            if(refined)++report.refinedDonors;
            std::vector<std::vector<Cut::AnchoredPoint>> liquid,gas;
            if(alpha>0&&alpha<1){liquid=source.anchoredPhasePolygons(static_cast<int>(si),true);gas=source.anchoredPhasePolygons(static_cast<int>(si),false);}
            Sum rawInventory,slack;
            for(auto index:bySource[si]){
                auto& overlap=overlaps[index];
                require((overlap.kind==Overlap::Kind::Refined)==refined,"donor mixes incompatible nested transfer kinds");
                if(alpha==0||alpha==1){
                    overlap.liquidMeasured=alpha==1?overlap.capacity:0;
                    overlap.gasMeasured=alpha==0?overlap.capacity:0;
                }else{
                    const auto integrate=[&](const auto& polygons){
                        Sum volume;
                        for(const auto& polygon:polygons){
                            if(refined)volume.add(Cut::anchoredPolygonMeasureInFrame(
                                clipInFrame(targetBoxes[overlap.targetCell],polygon,origin),origin,rz));
                            else volume.add(Cut::anchoredPolygonMeasureInFrame(polygon,origin,rz));
                        }
                        return volume.value;
                    };
                    overlap.liquidMeasured=integrate(liquid);overlap.gasMeasured=integrate(gas);
                }
                const auto l=overlap.liquidMeasured,g=overlap.gasMeasured,total=l+g;
                require(std::isfinite(l)&&std::isfinite(g)&&l>=0&&g>=0&&total>0,
                    "invalid same-fan two-phase measure source="+std::to_string(si));
                require(std::abs(total-overlap.capacity)<=arithmeticTolerance*overlap.capacity,
                    "same-fan phases do not cover actual overlap source="+std::to_string(si));
                if(!refined)overlap.inventoryContribution=cell.inventory;
                else if(l==0)overlap.inventoryContribution=0;
                else if(g==0)overlap.inventoryContribution=overlap.capacity;
                else {
                    overlap.inventoryContribution=overlap.capacity*(l/total);
                    std::ostringstream representation;representation<<"refined mixed overlap is unrepresentable source="<<si
                        <<" target="<<overlap.targetCell<<std::scientific<<" liquid="<<l<<" gas="<<g<<" capacity="<<overlap.capacity;
                    require(overlap.inventoryContribution>0&&overlap.inventoryContribution<overlap.capacity,representation.str());
                    slack.add(std::min(overlap.inventoryContribution,overlap.capacity-overlap.inventoryContribution));
                }
                rawInventory.add(overlap.inventoryContribution);
            }
            if(refined){
                const auto delta=cell.inventory-rawInventory.value;
                report.maximumRelativeRawInventoryResidual=std::max(report.maximumRelativeRawInventoryResidual,
                    phaseResidual(rawInventory.value,cell.inventory,cell.capacity));
                if(slack.value>0){
                    // Explicit LOCAL conservation closure, bounded by every
                    // actual mixed overlap's minority phase. Pure overlap
                    // states stay exact; no update clipping or phase floor.
                    std::ostringstream context;context<<"fan residual exceeds mixed-overlap conservation budget source="<<si
                        <<" owner="<<source.identity().value()<<std::scientific<<" delta="<<delta<<" slack="<<slack.value
                        <<" inventory="<<cell.inventory<<" capacity="<<cell.capacity<<" factor="<<delta/slack.value;
                    const auto factor=delta/slack.value;
                    const auto donorMinor=std::min(cell.inventory,cell.capacity-cell.inventory);
                    require(std::abs(delta)<=config.relativeVolumeTolerance*donorMinor,context.str());
                    report.maximumRelativeDonorClosureCorrection=std::max(report.maximumRelativeDonorClosureCorrection,std::abs(delta)/donorMinor);
                    require(std::abs(factor)<1,context.str());
                    report.maximumRelativeMinorPhaseCorrection=std::max(report.maximumRelativeMinorPhaseCorrection,std::abs(factor));
                    for(auto index:bySource[si]){
                        auto& overlap=overlaps[index];const auto l=overlap.inventoryContribution,g=overlap.capacity-l;
                        if(l==0||g==0)continue;
                        const auto minor=std::min(l,g),correction=factor*minor;
                        overlap.closureCorrection=correction;overlap.inventoryContribution+=correction;
                        report.maximumRelativeClosureCorrection=std::max(report.maximumRelativeClosureCorrection,std::abs(correction)/overlap.capacity);
                        require(overlap.inventoryContribution>0&&overlap.inventoryContribution<overlap.capacity,
                            "conservation closure cannot retain two representable phases source="+std::to_string(si));
                    }
                }else {
                    const auto minor=std::min(cell.inventory,cell.capacity-cell.inventory);
                    require(std::abs(delta)<=arithmeticTolerance*(minor>0?minor:cell.capacity),
                        "no mixed overlap can close non-arithmetic fan inventory residual source="+std::to_string(si));
                }
                Sum finalInventory;for(auto index:bySource[si])finalInventory.add(overlaps[index].inventoryContribution);
                const auto error=phaseResidual(finalInventory.value,cell.inventory,cell.capacity);
                report.maximumRelativeInventoryResidual=std::max(report.maximumRelativeInventoryResidual,error);
                require(error<=config.relativeVolumeTolerance,"refined donor inventory ledger failed source="+std::to_string(si));
            }
        }
    }
    std::vector<double> targetInventories()
    {
        std::vector<double> alpha(byTarget.size());Sum inventory,gasInventory,capacity,contribution,closure;
        for(std::size_t ti=0;ti<byTarget.size();++ti){
            Sum liquid,gas;bool same=false,refined=false;
            for(auto index:byTarget[ti]){
                const auto& o=overlaps[index];liquid.add(o.inventoryContribution);gas.add(o.capacity-o.inventoryContribution);
                contribution.add(o.inventoryContribution);closure.add(o.closureCorrection);
                same|=o.kind==Overlap::Kind::Unchanged;refined|=o.kind==Overlap::Kind::Refined;
            }
            const auto full=targetCapacity[ti];
            require(liquid.value>=0&&gas.value>=0,
                "negative target inventory cell="+std::to_string(ti));
            if(gas.value==0)require(std::abs(liquid.value-full)<=arithmeticTolerance*full,
                "pure-liquid target does not cover actual capacity cell="+std::to_string(ti));
            else require(liquid.value<full,
                "mixed target inventory outside actual capacity cell="+std::to_string(ti));
            if(same){require(byTarget[ti].size()==1,"unchanged leaf has additional coverage");
                alpha[ti]=source.alpha()[overlaps[byTarget[ti][0]].sourceCell];++report.unchangedTargets;}
            else {
                if(refined)++report.refinedTargets;else ++report.coarsenedTargets;
                if(liquid.value==0)alpha[ti]=0;
                else if(gas.value==0)alpha[ti]=1;
                else {
                    alpha[ti]=static_cast<double>(liquid.value/full);
                    require(alpha[ti]>0&&alpha[ti]<1,"mixed target alpha cannot represent both phases cell="+std::to_string(ti));
                }
            }
            const auto actual=static_cast<long double>(alpha[ti])*full;
            const bool pure=liquid.value==0||gas.value==0;
            const auto error=std::max(phaseResidual(actual,liquid.value,full,pure),
                phaseResidual(full-actual,gas.value,full,pure));
            report.maximumRelativeInventoryResidual=std::max(report.maximumRelativeInventoryResidual,error);
            require(error<=config.relativeVolumeTolerance,"target double inventory conversion failed cell="+std::to_string(ti));
            inventory.add(actual);gasInventory.add(full-actual);capacity.add(full);
        }
        Sum oldInventory,oldGas;for(const auto& c:source.cells()){oldInventory.add(c.inventory);oldGas.add(c.capacity-c.inventory);}
        report.sourceInventory=oldInventory.value;report.targetInventory=inventory.value;
        report.sourceGasInventory=oldGas.value;report.targetGasInventory=gasInventory.value;
        report.closureCorrection=closure.value;report.conversionResidual=inventory.value-contribution.value;
        const auto relative=[&](long double actual,long double expected){
            return std::abs(actual-expected)/(expected>0?expected:capacity.value);
        };
        const auto error=std::max(relative(inventory.value,oldInventory.value),relative(gasInventory.value,oldGas.value));
        report.maximumRelativeInventoryResidual=std::max(report.maximumRelativeInventoryResidual,error);
        require(error<=config.relativeVolumeTolerance,"global liquid/gas inventory transfer failed");
        return alpha;
    }
};
struct LocalPlicTransferResult
{
    std::vector<std::array<long double,2>> inventory;
    std::vector<FV::Dyadic> defects;
    std::vector<Overlap> overlaps;
    GeometricVofTransferReport report;
    std::shared_ptr<const Mesh::NativeLeafOverlap> coverage;
};
long double ownPhaseRelative(long double actual,long double expected)
{
    require(std::isfinite(actual)&&std::isfinite(expected)&&actual>=0&&expected>=0,
        "invalid independent phase ledger");
    if(expected==0){require(actual==0,"zero source phase gained inventory");return 0;}
    require(actual>0,"positive source phase lost inventory");
    return std::abs(actual-expected)/expected;
}
LocalPlicTransferResult transferLocalPlic(const QuadtreeMesh& sourceMesh,
    const GeometricVofPublication& source,const QuadtreeMesh& targetMesh,
    const NativeLeafDomain& targetDomain,GeometricVofTimeLayer time,const GeometricVofConfig& config)
{
    const bool rz=targetDomain.axisymmetric();
    auto coverage=std::make_shared<const Mesh::NativeLeafOverlap>(
        Mesh::NativeLeafOverlap::build(sourceMesh,source.domain(),targetMesh,targetDomain));
    const auto& pieces=coverage->pieces();const auto& bySource=coverage->bySource();
    const auto& byTarget=coverage->byTarget();const auto& cells=source.plicCells();
    require(cells.size()==bySource.size(),"LocalPlic source native count mismatch");
    std::vector<Overlap> overlaps(pieces.size());
    std::vector<FV::Dyadic> pieceDefects(pieces.size(),FV::Dyadic(0));
    GeometricVofTransferReport report;report.localPlic=true;
    const auto& sourceDefects=source.capacityDefects();
    require(sourceDefects.size()==cells.size(),"missing source finite-volume closure evidence");
    // This checks a retained source association; it cannot mint a new root.
    const auto exactDefect=[](const std::array<long double,2>& stock,long double capacity){
        FV::Dyadic result(stock[0]);result.add(stock[1]);result.add(capacity,-1);return result;
    };
    const auto certifiedBound=[](long double value,long double capacity,const FV::Dyadic& defect){
        if(!std::isfinite(value)||value<0)return false;
        FV::Dyadic excess(value);excess.add(capacity,-1);
        const FV::Dyadic zero(0);return excess.compare(zero)<=0||excess.compare(defect)<=0;
    };
    report.maximumRelativeCoverageResidual=coverage->maximumRelativeCoverageResidual();
    const auto geometryBudget=config.relativeVolumeTolerance+arithmeticTolerance;
    std::array<Sum,2> sourceTotals,targetTotals;Sum closureTotal;
    for(std::size_t i=0;i<pieces.size();++i){const auto& p=pieces[i];auto& o=overlaps[i];
        o.sourceCell=p.sourceCell;o.targetCell=p.targetCell;o.capacity=p.capacity;o.localPlic=true;
        using NativeKind=Mesh::NativeLeafOverlapPiece::Kind;
        o.kind=p.kind==NativeKind::Unchanged?Overlap::Kind::Unchanged:
            (p.kind==NativeKind::Refined?Overlap::Kind::Refined:Overlap::Kind::Coarsened);
    }
    for(std::size_t si=0;si<bySource.size();++si){
        const auto& sourceCell=cells[si];
        if(source.domain().cellRowByIndex()[si]<0){
            require(bySource[si].empty()&&source.cells()[si].capacity==0&&
                sourceCell.inventory[0]==0&&sourceCell.inventory[1]==0&&
                sourceDefects[si].compare(FV::Dyadic(0))==0,
                "inactive source cannot carry overlap, capacity, stock or closure evidence");
            continue;
        }
        const auto parentCapacity=Cut::capacity(sourceCell.rectangle,rz);
        require(exactDefect(sourceCell.inventory,parentCapacity).compare(sourceDefects[si])==0,
            "source finite-volume closure evidence differs from retained stocks");
        const auto& indices=bySource[si];require(!indices.empty(),"LocalPlic source lacks actual coverage");
        for(std::size_t phase=0;phase<2;++phase)sourceTotals[phase].add(sourceCell.inventory[phase]);
        const bool refined=overlaps[indices.front()].kind==Overlap::Kind::Refined;
        if(!refined){
            require(indices.size()==1,"unchanged/coarsened source has duplicate coverage");
            auto& o=overlaps[indices.front()];
            require(std::abs(o.capacity-parentCapacity)<=arithmeticTolerance*parentCapacity,
                "nonrefined source piece does not cover actual native capacity");
            o.phaseContribution=sourceCell.inventory;o.liquidMeasured=sourceCell.measured[0];
            o.gasMeasured=sourceCell.measured[1];o.geometricMeasureAvailable=sourceCell.resolved;
            pieceDefects[indices.front()]=sourceDefects[si];
            pieceDefects[indices.front()].add(parentCapacity);
            pieceDefects[indices.front()].add(o.capacity,-1);
        }else{
            ++report.refinedDonors;
            const bool pure=sourceCell.inventory[0]==0||sourceCell.inventory[1]==0;
            if(pure||!sourceCell.resolved){
                // Only the genuine pure/P0 tag selects this declared finite
                // partition. Ordinary resolved cuts never fall back to it.
                if(pure)++report.resolvedRefinedDonors;else ++report.unresolvedRefinedDonors;
                auto ordered=indices;
                std::sort(ordered.begin(),ordered.end(),[&](std::size_t a,std::size_t b){
                    const auto& x=pieces[a].rectangle;const auto& y=pieces[b].rectangle;
                    return std::tie(x.lower.x,x.lower.y,x.upper.x,x.upper.y)<
                        std::tie(y.lower.x,y.lower.y,y.upper.x,y.upper.y);
                });
                std::vector<FV::Dyadic> weights;weights.reserve(ordered.size());FV::Dyadic coverageExact(0);
                for(auto index:ordered){auto weight=FV::rectangleCapacityExact(pieces[index].rectangle,rz);
                    coverageExact.add(weight);weights.push_back(std::move(weight));}
                require(coverageExact.compare(FV::rectangleCapacityExact(sourceCell.rectangle,rz))==0,
                    "finite pure/P0 actual partition does not cover source rectangle");
                const auto liquid=FV::partition(sourceCell.inventory[0],weights);
                const auto gas=FV::partition(sourceCell.inventory[1],weights);
                const FV::Dyadic zero(0);
                require(liquid.values.size()==ordered.size()&&gas.values.size()==ordered.size()
                    &&liquid.exactAllocations.size()==ordered.size()&&gas.exactAllocations.size()==ordered.size()
                    &&liquid.sumError.compare(zero)==0&&gas.sumError.compare(zero)==0,
                    "finite pure/P0 phase partition is not exactly conservative");
                for(std::size_t k=0;k<ordered.size();++k){const auto index=ordered[k];auto& o=overlaps[index];
                    o.phaseContribution={liquid.values[k],gas.values[k]};
                    require(liquid.exactAllocations[k].compare(FV::Dyadic(liquid.values[k]))==0
                        &&gas.exactAllocations[k].compare(FV::Dyadic(gas.values[k]))==0,
                        "finite partition receipt differs from native allocation");
                    // The typed integer allocations derive from the actual
                    // source stock and certified geometric partition above.
                    pieceDefects[index]=liquid.exactAllocations[k];pieceDefects[index].add(gas.exactAllocations[k]);
                    pieceDefects[index].add(o.capacity,-1);
                    if(pure){o.geometricMeasureAvailable=true;
                        o.liquidMeasured=sourceCell.inventory[1]==0?o.capacity:0;
                        o.gasMeasured=sourceCell.inventory[0]==0?o.capacity:0;
                    }else{o.usedP0=true;++report.p0RefinedPieces;}
                }
            }else{
                ++report.resolvedRefinedDonors;std::array<Sum,2> raw;Sum slack;
                const auto cx=(static_cast<long double>(sourceCell.rectangle.lower.x)+sourceCell.rectangle.upper.x)/2;
                const auto cy=(static_cast<long double>(sourceCell.rectangle.lower.y)+sourceCell.rectangle.upper.y)/2;
                const auto nx=static_cast<long double>(sourceCell.normal.x),ny=static_cast<long double>(sourceCell.normal.y);
                const auto length=std::hypot(nx,ny);require(std::isfinite(length)&&length>0,"invalid parent PLIC normal");
                for(auto index:indices){auto& o=overlaps[index];const auto& box=pieces[index].rectangle;
                    o.geometricMeasureAvailable=true;
                    const auto tx=(static_cast<long double>(box.lower.x)+box.upper.x)/2;
                    const auto ty=(static_cast<long double>(box.lower.y)+box.upper.y)/2;
                    const auto offset=sourceCell.offset+(nx*(cx-tx)+ny*(cy-ty))/length;
                    o.liquidMeasured=Cut::volume(box,sourceCell.normal,offset,rz);
                    o.gasMeasured=Cut::volume(box,Vector2D(-sourceCell.normal.x,-sourceCell.normal.y),-offset,rz);
                    report.localPlaneCuts+=2;
                    o.phaseContribution={o.liquidMeasured,o.gasMeasured};
                    pieceDefects[index]=exactDefect(o.phaseContribution,o.capacity);
                    for(std::size_t phase=0;phase<2;++phase){const auto value=o.phaseContribution[phase];
                        require(std::isfinite(value)&&value>=0&&value<=o.capacity,"invalid raw independent PLIC piece");raw[phase].add(value);}
                    const auto joint=std::abs((o.liquidMeasured+o.gasMeasured)-o.capacity)/o.capacity;
                    require(joint<=arithmeticTolerance,"raw PLIC piece does not cover actual capacity");
                    report.maximumRelativeLocalGeometryResidual=std::max(report.maximumRelativeLocalGeometryResidual,joint);
                    if(o.liquidMeasured>0&&o.gasMeasured>0)slack.add(std::min(o.liquidMeasured,o.gasMeasured));
                }
                std::array<long double,2> delta{};bool corrected=false;
                for(std::size_t phase=0;phase<2;++phase){
                    const auto measuredDefect=std::abs(raw[phase].value-sourceCell.measured[phase])/parentCapacity;
                    const auto stockDefect=std::abs(raw[phase].value-sourceCell.inventory[phase])/parentCapacity;
                    report.maximumRelativeLocalGeometryResidual=std::max(report.maximumRelativeLocalGeometryResidual,measuredDefect);
                    report.maximumRelativeRawInventoryResidual=std::max(report.maximumRelativeRawInventoryResidual,stockDefect);
                    require(measuredDefect<=geometryBudget&&stockDefect<=geometryBudget,
                        "resolved parent raw geometry/stock budget failed source="+std::to_string(si));
                    delta[phase]=sourceCell.inventory[phase]-raw[phase].value;
                    if(delta[phase]!=0&&slack.value==0)
                        require(std::abs(delta[phase])<=arithmeticTolerance*sourceCell.inventory[phase],
                            "resolved donor has no mixed overlap slack source="+std::to_string(si));
                    if(slack.value>0)require(std::abs(delta[phase])<slack.value,"resolved correction exceeds donor mixed slack");
                }
                if(slack.value>0)for(auto index:indices){auto& o=overlaps[index];
                    if(o.liquidMeasured==0||o.gasMeasured==0)continue;
                    const auto weight=std::min(o.liquidMeasured,o.gasMeasured)/slack.value;
                    for(std::size_t phase=0;phase<2;++phase){o.phaseCorrection[phase]=weight*delta[phase];
                        const auto before=o.phaseContribution[phase];
                        o.phaseContribution[phase]+=o.phaseCorrection[phase];
                        pieceDefects[index].add(o.phaseContribution[phase]);pieceDefects[index].add(before,-1);
                        require(o.phaseContribution[phase]>0&&o.phaseContribution[phase]<o.capacity,
                            "donor-local prolongation cannot retain both positive phases");
                        const auto relative=std::abs(o.phaseCorrection[phase])/sourceCell.inventory[phase];
                        report.maximumRelativePhaseCorrection[phase]=std::max(report.maximumRelativePhaseCorrection[phase],relative);
                        corrected|=o.phaseCorrection[phase]!=0;}
                }
                if(corrected)++report.correctedRefinedDonors;
            }
            if(sourceCell.resolved&&!pure){
            // Ordinary resolved cuts retain their original own-phase
            // arithmetic-tail law and strict overlap bounds. Finite pure/P0
            // partitions above already telescope exactly to each source phase.
            for(std::size_t phase=0;phase<2;++phase){Sum total;
                for(auto index:indices)total.add(overlaps[index].phaseContribution[phase]);
                const auto stock=sourceCell.inventory[phase],tail=stock-total.value;
                if(stock==0){require(total.value==0,"zero refined phase gained inventory");continue;}
                require(std::abs(tail)<=arithmeticTolerance*stock,
                    "prolongation left a non-arithmetic own-phase tail source="+std::to_string(si));
                if(tail==0)continue;
                std::size_t chosen=pieces.size();long double largest=-1;
                for(auto index:indices){const auto& o=overlaps[index];const auto before=o.phaseContribution[phase];
                    if(before<=0)continue;const auto after=before+tail;
                    const auto other=o.phaseContribution[1-phase];
                    if(!std::isfinite(after)||after<=0||after>o.capacity)continue;
                    if(other>0&&after>=o.capacity)continue;
                    if(std::abs((after+other)-o.capacity)>config.relativeVolumeTolerance*o.capacity)continue;
                    if(before>largest||(before==largest&&(chosen==pieces.size()||o.targetCell<overlaps[chosen].targetCell))){
                        chosen=index;largest=before;}
                }
                if(chosen!=pieces.size()){auto& o=overlaps[chosen];const auto before=o.phaseContribution[phase];
                    o.phaseContribution[phase]=before+tail;o.phaseRoundoffClosure[phase]=o.phaseContribution[phase]-before;
                    pieceDefects[chosen].add(o.phaseContribution[phase]);pieceDefects[chosen].add(before,-1);
                    if(o.phaseRoundoffClosure[phase]!=0){++report.roundoffClosedPhases;
                        report.maximumRelativePhaseRoundoffClosure[phase]=std::max(
                            report.maximumRelativePhaseRoundoffClosure[phase],std::abs(o.phaseRoundoffClosure[phase])/stock);}
                }
            }
            }
        }
        for(auto index:indices)require(exactDefect(overlaps[index].phaseContribution,overlaps[index].capacity).compare(pieceDefects[index])==0,
            "overlap closure receipt differs from declared transfer operation");
        for(std::size_t phase=0;phase<2;++phase){Sum total;
            for(auto index:indices){const auto& o=overlaps[index];const auto value=o.phaseContribution[phase];
                require(certifiedBound(value,o.capacity,pieceDefects[index]),"independent transferred phase outside certified overlap capacity");
                total.add(value);}
            const auto residual=ownPhaseRelative(total.value,sourceCell.inventory[phase]);
            report.maximumRelativeDonorPhaseResidual[phase]=std::max(report.maximumRelativeDonorPhaseResidual[phase],residual);
            require(residual<=config.relativeVolumeTolerance,"donor own-phase conservation failed source="+std::to_string(si));
        }
    }
    std::vector<std::array<long double,2>> targetInventory(byTarget.size());
    std::vector<FV::Dyadic> targetDefects(byTarget.size(),FV::Dyadic(0));
    for(std::size_t ti=0;ti<byTarget.size();++ti){
        if(targetDomain.cellRowByIndex()[ti]<0){
            require(byTarget[ti].empty(),"inactive target cannot carry physical overlap");continue;
        }
        std::array<Sum,2> totals;bool same=false,refined=false;
        std::array<FV::Dyadic,2> exactTotals{FV::Dyadic(0),FV::Dyadic(0)};
        const auto full=Cut::capacity(coverage->targetRectangles()[ti],rz);
        auto& predicted=targetDefects[ti];predicted.add(full,-1);
        for(auto index:byTarget[ti]){auto& o=overlaps[index];
            const auto closure=std::abs((o.phaseContribution[0]+o.phaseContribution[1])-o.capacity)/o.capacity;
            require(closure<=config.relativeVolumeTolerance,"transferred piece joint phase closure failed");
            for(std::size_t phase=0;phase<2;++phase)exactTotals[phase].add(o.phaseContribution[phase]);
            predicted.add(pieceDefects[index]);predicted.add(o.capacity);
            o.inventoryContribution=o.phaseContribution[0];
            o.closureCorrection=o.phaseCorrection[0]+o.phaseRoundoffClosure[0];closureTotal.add(o.closureCorrection);
            same|=o.kind==Overlap::Kind::Unchanged;refined|=o.kind==Overlap::Kind::Refined;
        }
        // Coarsening is one final RNE of each independent exact donor sum.
        // Its errors extend predecessor receipts, not a freshly self-signed q-C.
        for(std::size_t phase=0;phase<2;++phase){const auto rounded=FV::round(exactTotals[phase]);
            totals[phase].value=rounded.value;predicted.add(rounded.error);}
        for(std::size_t phase=0;phase<2;++phase){const auto value=totals[phase].value;
            if(!certifiedBound(value,full,predicted)){
                // Retain the original detailed native-geometry failure
                // context when the derived finite stock bound is rejected.
                std::ostringstream context;context<<"target independent phase outside actual capacity"
                    <<" rz="<<rz<<" source_owner="<<source.identity().value()
                    <<" source_leaves="<<sourceMesh.cellIndexNum()<<" target_leaves="<<targetMesh.cellIndexNum()
                    <<" target="<<ti<<" target_label="<<targetMesh.cellIndexToLabel[ti]
                    <<" target_level="<<targetMesh.cells(targetMesh.cellIndexToLabel[ti]).level()
                    <<" phase="<<phase<<" clock_step="<<time.step<<" clock_position="<<static_cast<int>(time.position);
                const auto scalar=[&](const std::string& name,long double x){
                    context<<' '<<name<<'='<<std::scientific<<std::setprecision(std::numeric_limits<long double>::max_digits10)<<x
                        <<' '<<name<<"_hex="<<std::hexfloat<<x<<std::defaultfloat;
                };
                const auto box=[&](const std::string& name,const Cut::Rectangle& r){
                    scalar(name+"_lo_x",r.lower.x);scalar(name+"_lo_y",r.lower.y);
                    scalar(name+"_hi_x",r.upper.x);scalar(name+"_hi_y",r.upper.y);
                };
                scalar("clock_begin",time.begin);scalar("clock_end",time.end);
                const auto sourceClock=source.timeLayer();context<<" source_step="<<sourceClock.step<<" source_position="<<static_cast<int>(sourceClock.position);
                scalar("source_begin",sourceClock.begin);scalar("source_end",sourceClock.end);
                scalar("q",value);scalar("capacity",full);scalar("excess",value-full);
                scalar("q_liquid",totals[0].value);scalar("q_gas",totals[1].value);
                scalar("sum_correction_liquid",totals[0].correction);scalar("sum_correction_gas",totals[1].correction);
                const auto ulp=std::nextafter(full,std::numeric_limits<long double>::infinity())-full;
                scalar("capacity_ulp_up",ulp);scalar("excess_in_up_ulps",(value-full)/ulp);
                const int row=targetDomain.cellRowByIndex()[ti];context<<" target_domain_row="<<row;
                if(row>=0&&static_cast<std::size_t>(row)<targetDomain.cellVolumes().size())scalar("domain_capacity",targetDomain.cellVolumes()[row]);
                scalar("corner_capacity",Cut::capacity(rectangle(targetMesh,static_cast<int>(ti)),rz));
                box("target_box",coverage->targetRectangles()[ti]);box("target_corner_box",rectangle(targetMesh,static_cast<int>(ti)));
                Sum covered;for(const auto index:byTarget[ti])covered.add(overlaps[index].capacity);
                scalar("piece_capacity_sum",covered.value);scalar("coverage_excess",covered.value-full);
                context<<" donors="<<byTarget[ti].size()<<" donor_details_limit=64";
                std::size_t detail=0;
                for(const auto index:byTarget[ti]){
                    if(detail==64){context<<" donor_details_truncated=1";break;}
                    const auto& o=overlaps[index];const auto& donor=cells[o.sourceCell];
                    const std::string n="d"+std::to_string(detail++);
                    context<<' '<<n<<"_index="<<index<<' '<<n<<"_kind="<<static_cast<int>(o.kind)
                        <<' '<<n<<"_source="<<o.sourceCell<<' '<<n<<"_source_label="<<sourceMesh.cellIndexToLabel[o.sourceCell]
                        <<' '<<n<<"_source_level="<<sourceMesh.cells(sourceMesh.cellIndexToLabel[o.sourceCell]).level()
                        <<' '<<n<<"_resolved="<<donor.resolved<<' '<<n<<"_p0="<<o.usedP0
                        <<' '<<n<<"_geometry_available="<<o.geometricMeasureAvailable;
                    scalar(n+"_source_capacity",Cut::capacity(donor.rectangle,rz));scalar(n+"_piece_capacity",o.capacity);
                    scalar(n+"_stock_l",donor.inventory[0]);scalar(n+"_stock_g",donor.inventory[1]);
                    scalar(n+"_source_measured_l",donor.measured[0]);scalar(n+"_source_measured_g",donor.measured[1]);
                    scalar(n+"_piece_measured_l",o.liquidMeasured);scalar(n+"_piece_measured_g",o.gasMeasured);
                    scalar(n+"_contribution_l",o.phaseContribution[0]);scalar(n+"_contribution_g",o.phaseContribution[1]);
                    scalar(n+"_correction_l",o.phaseCorrection[0]);scalar(n+"_correction_g",o.phaseCorrection[1]);
                    scalar(n+"_tail_l",o.phaseRoundoffClosure[0]);scalar(n+"_tail_g",o.phaseRoundoffClosure[1]);
                    scalar(n+"_normal_x",donor.normal.x);scalar(n+"_normal_y",donor.normal.y);scalar(n+"_offset",donor.offset);
                    box(n+"_source_box",donor.rectangle);box(n+"_piece_box",pieces[index].rectangle);
                }
                require(false,context.str());
            }
            targetInventory[ti][phase]=value;targetTotals[phase].add(value);}
        require(exactDefect(targetInventory[ti],full).compare(predicted)==0,
            "target closure receipt differs from actual independent summation");
        const auto closure=std::abs((totals[0].value+totals[1].value)-full)/full;
        report.maximumRelativeTargetPhaseClosure=std::max(report.maximumRelativeTargetPhaseClosure,closure);
        require(closure<=config.relativeVolumeTolerance,"target independent phase closure failed");
        if(same){require(byTarget[ti].size()==1,"unchanged target has multiple overlaps");
            require(targetInventory[ti]==cells[overlaps[byTarget[ti].front()].sourceCell].inventory,
                "unchanged target did not retain both stocks exactly");++report.unchangedTargets;}
        else if(refined)++report.refinedTargets;else ++report.coarsenedTargets;
    }
    for(std::size_t phase=0;phase<2;++phase){
        report.sourcePhaseInventory[phase]=sourceTotals[phase].value;report.targetPhaseInventory[phase]=targetTotals[phase].value;
        report.relativeGlobalPhaseResidual[phase]=ownPhaseRelative(targetTotals[phase].value,sourceTotals[phase].value);
        require(report.relativeGlobalPhaseResidual[phase]<=config.relativeVolumeTolerance,"global own-phase conservation failed");
        report.maximumRelativeInventoryResidual=std::max({report.maximumRelativeInventoryResidual,
            report.maximumRelativeDonorPhaseResidual[phase],report.relativeGlobalPhaseResidual[phase]});
    }
    report.sourceInventory=report.sourcePhaseInventory[0];report.targetInventory=report.targetPhaseInventory[0];
    report.sourceGasInventory=report.sourcePhaseInventory[1];report.targetGasInventory=report.targetPhaseInventory[1];
    report.sourceMeasured=source.report().measured;report.closureCorrection=closureTotal.value;
    return {std::move(targetInventory),std::move(targetDefects),std::move(overlaps),report,std::move(coverage)};
}

}
GeometricVofMeshTransfer::GeometricVofMeshTransfer(GeometricVofIdentity sourceIdentity,
    GeometricVofPublication candidate,std::vector<GeometricVofTransferOverlap> overlaps,GeometricVofTransferReport report,
    std::shared_ptr<const Mesh::NativeLeafOverlap> nativeOverlaps)
    :sourceIdentity_(sourceIdentity),candidate_(std::move(candidate)),overlaps_(std::move(overlaps)),report_(report),
     nativeOverlaps_(std::move(nativeOverlaps)){}
const Mesh::NativeLeafOverlap& GeometricVofMeshTransfer::nativeOverlaps() const
{
    require(static_cast<bool>(nativeOverlaps_),"native overlap port is available only for LocalPlic transfer");
    return *nativeOverlaps_;
}
GeometricVofMeshTransfer GeometricVofMeshTransfer::prepare(const QuadtreeMesh& sourceMesh,
    const GeometricVofPublication& source,const QuadtreeMesh& targetMesh,NativeLeafDomain targetDomain,
    GeometricVofTimeLayer time,const GeometricVofConfig& config)
{
    source.validateFor(sourceMesh,source.identity(),time);
    require(time.position==GeometricVofTimeLayer::Position::End,"regrid requires the same accepted End clock");
    require(std::isfinite(config.relativeVolumeTolerance)&&config.relativeVolumeTolerance>0&&config.relativeVolumeTolerance<=1e-12L,
        "transfer tolerance must be finite, positive and no looser than 1e-12");
    targetDomain.validateMeshIdentity(targetMesh);
    require(source.domain().axisymmetric()==targetDomain.axisymmetric(),"source/target physical measure mismatch");
    require(sourceMesh.Nx==targetMesh.Nx&&sourceMesh.Ny==targetMesh.Ny&&
        sourceMesh.minCoord.x==targetMesh.minCoord.x&&sourceMesh.minCoord.y==targetMesh.minCoord.y&&
        sourceMesh.maxCoord.x==targetMesh.maxCoord.x&&sourceMesh.maxCoord.y==targetMesh.maxCoord.y&&
        sourceMesh.cellSize.x==targetMesh.cellSize.x&&sourceMesh.cellSize.y==targetMesh.cellSize.y,
        "regrid must retain the actual background lattice/domain");
    require(targetDomain.cellRowByIndex().size()==static_cast<std::size_t>(targetMesh.cellIndexNum()),
        "target selected map must cover all native indices");
    if(source.geometryKind()==GeometricVofGeometryKind::LocalPlic){
        auto local=transferLocalPlic(sourceMesh,source,targetMesh,targetDomain,time,config);
        // Only the typed transfer member can publish the derived association.
        const bool complete=source.domain().cellIndexByRow().size()==source.domain().cellRowByIndex().size()&&
            targetDomain.cellIndexByRow().size()==targetDomain.cellRowByIndex().size();
        auto candidate=[&]{
            if(complete)return GeometricVofPublication::prepareCertifiedLocalPlicInventory(targetMesh,std::move(targetDomain),
                time,source,GeometricVofPublication::CertifiedOperation::MeshTransfer,local.inventory,local.defects,config);
            const GeometricVofPublication::MeshTransferReceipt receipt(sourceMesh,source,targetMesh,*local.coverage);
            return GeometricVofPublication::prepareCertifiedLocalPlicInventory(targetMesh,std::move(targetDomain),
                time,source,GeometricVofPublication::CertifiedOperation::MeshTransfer,local.inventory,local.defects,config,nullptr,&receipt);
        }();
        candidate.validateFor(targetMesh,candidate.identity(),time);
        require(candidate.identity()!=source.identity(),"regrid candidate must have a fresh owner");
        local.report.targetMeasured=candidate.report().measured;
        local.report.maximumRelativeGeometryResidual=candidate.report().maximumRelativeVolumeResidual;
        return GeometricVofMeshTransfer(source.identity(),std::move(candidate),std::move(local.overlaps),
            local.report,std::move(local.coverage));
    }
    require(targetDomain.cellIndexByRow().size()==static_cast<std::size_t>(targetMesh.cellIndexNum()),
        "compatible fan target domain must cover all actual leaves");
    TransferBuilder builder{sourceMesh,source,targetMesh,targetDomain,config,{},{},{},{},{},{},{},targetDomain.axisymmetric()};
    builder.coverage();builder.donorInventories();const auto alpha=builder.targetInventories();
    auto candidate=GeometricVofPublication::prepare(targetMesh,std::move(targetDomain),time,alpha,config);
    candidate.validateFor(targetMesh,candidate.identity(),time);
    auto& report=builder.report;
    report.sourceMeasured=source.report().measured;report.targetMeasured=candidate.report().measured;
    report.maximumRelativeGeometryResidual=candidate.report().maximumRelativeVolumeResidual;
    report.sourceComponents=source.topology().components;report.targetComponents=candidate.topology().components;
    report.sourceJunctions=source.topology().junctionVertices.size();report.targetJunctions=candidate.topology().junctionVertices.size();
    require(report.maximumRelativeGeometryResidual<=config.relativeVolumeTolerance,"target geometry/inventory mismatch");
    require(candidate.identity()!=source.identity(),"regrid candidate must have a fresh owner");
    return GeometricVofMeshTransfer(source.identity(),std::move(candidate),std::move(builder.overlaps),report);
}
}
