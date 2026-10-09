#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Core/Mesh/CartesianLeafMesh.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>

namespace
{
void require(bool valid,const char* message)
{ if (!valid) throw std::invalid_argument(std::string("[InterfaceConservativePhase] ")+message); }
constexpr long double eps=128*std::numeric_limits<long double>::epsilon();
void validateStep(const InterfaceStepProgress& previous,const InterfaceStepTime& time)
{
    require(previous.completedSteps<std::numeric_limits<int>::max()
        && time.index==previous.completedSteps+1 && time.begin==previous.endTime
        && std::isfinite(time.begin) && std::isfinite(time.end) && time.end>time.begin
        && std::isfinite(time.dt) && time.dt>0
        && std::abs(time.end-time.begin-time.dt)<=16*std::numeric_limits<double>::epsilon()
            *std::max({time.dt,std::abs(time.begin),std::abs(time.end)}),"invalid external step clock");
}
}
InterfaceConservativePhase InterfaceConservativePhase::initialize(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const std::vector<long double>& volumes)
{
    require(domain.uniformLeaves(),"direct mixed-leaf initialization requires a conservative mesh transfer");
    return initializeNative(mesh,domain,volumes);
}
InterfaceConservativePhase InterfaceConservativePhase::initializeNative(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const std::vector<long double>& volumes)
{
    domain.validateMeshIdentity(mesh);
    require(volumes.size()==static_cast<std::size_t>(mesh.cellIndexNum()),"volume array has wrong native layout");
    InterfaceConservativePhase state;
    state.revision_=std::make_shared<const Revision>();
    state.identity_=domain.meshIdentity(); state.axisymmetric_=domain.axisymmetric();
    state.cells_=domain.cellIndexByRow(); state.liquid_=volumes;
    state.capacity_.assign(volumes.size(),std::numeric_limits<long double>::quiet_NaN());
    for (std::size_t row=0;row<state.cells_.size();++row) {
        const int cell=state.cells_[row]; const long double capacity=domain.cellVolumes()[row];
        require(std::isfinite(capacity) && capacity>0,"invalid cell capacity");
        require(std::isfinite(volumes[cell]) && volumes[cell]>=0 && volumes[cell]<=capacity,
            "initial liquid volume is outside cell capacity");
        state.capacity_[cell]=capacity; state.report_.liquidBefore+=volumes[cell];
    }
    for (std::size_t cell=0;cell<volumes.size();++cell)
        if (std::isnan(state.capacity_[cell])) require(std::isnan(volumes[cell]),"inactive cell must be NaN");
    state.report_.liquidAfter=state.report_.liquidBefore;
    return state;
}
InterfaceConservativePhase InterfaceConservativePhase::transferMesh(const QuadtreeMesh& oldMesh,
    const NativeLeafDomain& oldDomain,const InterfaceConservativePhase& old,
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const std::vector<Vector2D>& normals,
    const std::vector<std::optional<Quantum::Geometry::QuadraticCellCut::Shape>>& curves,
    const std::vector<long double>& geometricVolumes)
{
    namespace P=Quantum::Geometry::PlanarCellCut;
    namespace Q=Quantum::Geometry::QuadraticCellCut;
    old.validate(oldMesh,oldDomain);domain.validateMeshIdentity(mesh);
    validateCartesianLeafMesh(oldMesh);validateCartesianLeafMesh(mesh);
    require(oldDomain.axisymmetric()==domain.axisymmetric() && oldMesh.Nx==mesh.Nx && oldMesh.Ny==mesh.Ny
        && oldMesh.minCoord.x==mesh.minCoord.x && oldMesh.minCoord.y==mesh.minCoord.y
        && oldMesh.maxCoord.x==mesh.maxCoord.x && oldMesh.maxCoord.y==mesh.maxCoord.y,
        "mesh transfer requires unchanged background lattice and measure");
    require(normals.size()==static_cast<std::size_t>(oldMesh.cellIndexNum())
        && (curves.empty() || curves.size()==normals.size()),"mesh transfer donor layout mismatch");
    const auto box=[](const QuadtreeMesh& m,int c) {
        const auto& cell=m.cells(m.cellIndexToLabel[c]);
        return P::Rectangle{m.verts(cell.cellToVert[0]).coord(),m.verts(cell.cellToVert[2]).coord()};
    };
    const bool rz=domain.axisymmetric();const CartesianLeafLocator locator(mesh);
    const auto count=static_cast<std::size_t>(mesh.cellIndexNum());
    require(geometricVolumes.empty() || geometricVolumes.size()==count,"mesh transfer geometric reference has wrong target layout");
    const auto nan=std::numeric_limits<long double>::quiet_NaN();
    std::vector<long double> liquid(count,nan),coverage(count,0);
    std::vector<int> contributors(count,0);
    std::vector<bool> allLiquid(count,true);
    for(int c:domain.cellIndexByRow()) liquid[c]=0;
    MeshTransferReport report;report.liquidBefore=old.report_.liquidAfter;
    for(int c=0;c<oldMesh.cellIndexNum();++c) {
        const auto parent=box(oldMesh,c);
        struct Child {int cell;P::Rectangle box;long double capacity,liquid=0;};
        std::vector<Child> children;long double covered=0;
        for(int target:locator.intersectingBox(parent.lower,parent.upper)) {
            const auto next=box(mesh,target);
            const P::Rectangle overlap{{std::max(parent.lower.x,next.lower.x),std::max(parent.lower.y,next.lower.y)},
                {std::min(parent.upper.x,next.upper.x),std::min(parent.upper.y,next.upper.y)}};
            if(!(overlap.lower.x<overlap.upper.x && overlap.lower.y<overlap.upper.y)) continue;
            const bool active=oldDomain.cellRowByIndex()[c]>=0;
            require(active==(domain.cellRowByIndex()[target]>=0),"mesh transfer changes active/solid coverage");
            if(!active) continue;
            const bool contains=parent.lower.x<=next.lower.x && parent.lower.y<=next.lower.y
                && parent.upper.x>=next.upper.x && parent.upper.y>=next.upper.y;
            const bool contained=next.lower.x<=parent.lower.x && next.lower.y<=parent.lower.y
                && next.upper.x>=parent.upper.x && next.upper.y>=parent.upper.y;
            require(contains || contained,"mesh transfer overlaps are not native parent/child cells");
            const auto capacity=P::capacity(overlap,rz);covered+=capacity;coverage[target]+=capacity;
            children.push_back({target,overlap,contains?domain.cellVolumes()[domain.cellRowByIndex()[target]]:old.capacity_[c]});
        }
        if(oldDomain.cellRowByIndex()[c]<0) continue;
        const auto capacity=old.capacity_[c],amount=old.liquid_[c];
        const auto tolerance=(P::relativeVolumeResolution+32*std::numeric_limits<long double>::epsilon()*children.size())*capacity;
        require(!children.empty() && std::abs(covered-capacity)<=tolerance,"mesh transfer source partition is incomplete");
        if(children.size()==1) children[0].liquid=amount;
        else {
            ++report.refinedSourceCells;
            const bool tail=std::min(amount,capacity-amount)<=P::relativeVolumeResolution*capacity;
            const auto geometricPartition=[&]() {
                if(geometricVolumes.empty() || amount==0 || amount==capacity || tail) return false;
                std::vector<long double> weights;weights.reserve(children.size());
                long double minimum=0,maximum=0,referenceSum=0;
                for(auto& child:children) {
                    const auto value=geometricVolumes[child.cell];
                    require(std::isfinite(value) && value>=0 && value<=child.capacity,
                        "mesh subdivision reference is outside native capacity");
                    child.liquid=value;referenceSum+=value;
                    const auto fraction=value/child.capacity;
                    weights.push_back(value*(1-fraction)); // Positive only on a geometric cut.
                    if(value==child.capacity) minimum+=child.capacity;
                    if(value>0) maximum+=child.capacity;
                }
                if(amount<minimum-tolerance || amount>maximum+tolerance) {++report.geometricFallbackCells;return false;}
                // Minimize sum (V_child-V_geometry)^2/(2*w_child), with
                // 0<=V_child<=capacity and sum V_child=the OLD parent C.
                // Fix only active bounds; free children keep their original
                // geometric reference until the final common multiplier.
                for(std::size_t iteration=0;iteration<=children.size();++iteration) {
                    long double sum=0,weight=0;
                    for(std::size_t i=0;i<children.size();++i) {sum+=children[i].liquid;weight+=weights[i];}
                    const auto defect=amount-sum;
                    if(std::abs(defect)<=tolerance) break;
                    if(!(weight>0)) {++report.geometricFallbackCells;return false;}
                    bool saturated=false;
                    for(std::size_t i=0;i<children.size();++i) if(weights[i]>0) {
                        auto& child=children[i];const auto proposed=child.liquid+defect*weights[i]/weight;
                        if(proposed<0 || proposed>child.capacity) {
                            child.liquid=proposed<0?0:child.capacity;weights[i]=0;saturated=true;
                        }
                    }
                    if(saturated) continue;
                    for(std::size_t i=0;i<children.size();++i) if(weights[i]>0)
                        children[i].liquid+=defect*weights[i]/weight;
                    break;
                }
                ++report.geometricSourceCells;report.geometricCorrectionL1+=std::abs(amount-referenceSum);return true;
            };
            const bool geometric=geometricPartition();
            std::optional<P::Reconstruction> plane;std::optional<Q::Reconstruction> quadratic;
            if(amount>0 && amount<capacity && !tail && !geometric) {
                if(!curves.empty() && curves[c]) quadratic=Q::reconstruct(parent,*curves[c],amount/capacity,rz);
                else plane=P::reconstruct(parent,normals[c],amount/capacity,rz);
            } else if(amount>0 && amount<capacity && tail) ++report.uniformTailSplits;
            long double sum=0;
            for(auto& child:children) {
                if(!geometric) {
                    if(amount==0) child.liquid=0;
                    else if(amount==capacity) child.liquid=child.capacity;
                    else if(tail) child.liquid=amount/capacity*child.capacity;
                    else {
                        const P::Point origin{parent.lower.x,parent.lower.y};
                        const long double x0=static_cast<long double>(child.box.lower.x)-origin.x,
                            y0=static_cast<long double>(child.box.lower.y)-origin.y,
                            x1=static_cast<long double>(child.box.upper.x)-origin.x,
                            y1=static_cast<long double>(child.box.upper.y)-origin.y;
                        const std::vector<P::Point> polygon{{x0,y0},{x1,y0},{x1,y1},{x0,y1}};
                        const auto cut=quadratic?Q::overlapInFrame(parent,polygon,origin,*quadratic,rz):
                            P::overlapInFrame(parent,polygon,origin,plane->normal,plane->offset,rz);
                        child.liquid=cut.liquid==cut.total?child.capacity:cut.liquid;
                    }
                }
                require(std::isfinite(child.liquid) && child.liquid>=0 && child.liquid<=child.capacity,
                    "mesh transfer child volume is outside its capacity");
                sum+=child.liquid;
            }
            const auto defect=amount-sum;
            require(std::abs(defect)<=tolerance,"mesh transfer donor partition does not close");
            // Correct only a measured roundoff defect, within existing bounds.
            // Pure-capacity representation can have no admissible recipient;
            // retain/report that roundoff closure instead of inventing mass.
            auto best=children.end();long double room=0;
            for(auto it=children.begin();it!=children.end();++it) {
                // A pure parent remains pure in every child. Native RZ child
                // capacities may sum a few ulps away from the native parent
                // capacity; retain that already bounded partition residual
                // instead of making one full child artificially mixed.
                if(amount==0 || amount==capacity)continue;
                // The geometric partition fixes pure children, including in
                // this closure step. Depositing an ulp there seeds a spurious
                // PLIC fragment that subsequent refinement can resolve.
                if(geometric && (geometricVolumes[it->cell]==0
                    || geometricVolumes[it->cell]==it->capacity)) continue;
                const auto available=defect>=0?it->capacity-it->liquid:it->liquid;
                if(available>room){room=available;best=it;}
            }
            if(defect!=0 && best!=children.end() && room>=std::abs(defect)) {
                const auto previous=best->liquid;best->liquid+=defect;
                sum+=best->liquid-previous;
                report.maxPartitionCorrection=std::max(report.maxPartitionCorrection,std::abs(best->liquid-previous));
            }
            report.maxSourceClosure=std::max(report.maxSourceClosure,std::abs(sum-amount));
        }
        for(const auto& child:children) {
            liquid[child.cell]+=child.liquid;++contributors[child.cell];
            allLiquid[child.cell]=allLiquid[child.cell] && amount==capacity;
        }
    }
    for(int c:domain.cellIndexByRow()) {
        const auto capacity=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        require(contributors[c]>0 && std::abs(coverage[c]-capacity)<=eps*(contributors[c]+1)*capacity,
            "mesh transfer target partition is incomplete");
        if(contributors[c]>1) {
            ++report.coarsenedTargetCells;
            if(allLiquid[c]) {
                require(std::abs(liquid[c]-capacity)<=eps*(contributors[c]+1)*capacity,"pure target capacity mismatch");
                report.pureCapacityAdjustment+=capacity-liquid[c];
                liquid[c]=capacity; // Canonical pure measure; global roundoff remains in the ledger.
            }
        }
    }
    auto next=initializeNative(mesh,domain,liquid);
    report.liquidAfter=next.report_.liquidAfter;report.closure=report.liquidAfter-report.liquidBefore;
    long double scale=0;for(int c:old.cells_) scale+=old.capacity_[c];
    require(std::abs(report.closure)<=eps*(old.cells_.size()+domain.cellIndexByRow().size()+1)*scale,
        "mesh transfer inventory closure exceeds summation precision");
    next.progress_=old.progress_;next.meshSourceRevision_=old.revision_;next.meshTransferReport_=report;
    return next;
}
bool InterfaceConservativePhase::transferredFrom(const InterfaceConservativePhase& source) const noexcept
{
    return meshSourceRevision_ && meshSourceRevision_==source.revision_ && meshTransferReport_
        && progress_.completedSteps==source.progress_.completedSteps && progress_.endTime==source.progress_.endTime;
}
void InterfaceConservativePhase::validate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) const
{
    domain.validateMeshIdentity(mesh);
    require(identity_ && identity_==domain.meshIdentity() && axisymmetric_==domain.axisymmetric()
        && cells_==domain.cellIndexByRow(),"inventory domain identity changed");
}
bool InterfaceConservativePhase::follows(const InterfaceConservativePhase& previous,const InterfaceStepTime& time) const noexcept
{
    return revision_ && parentRevision_ && parentRevision_==previous.revision_ && identity_==previous.identity_
        && time.index==stepTime_.index && time.begin==stepTime_.begin && time.end==stepTime_.end && time.dt==stepTime_.dt
        && previous.progress_.completedSteps<std::numeric_limits<int>::max()
        && time.index==previous.progress_.completedSteps+1 && time.begin==previous.progress_.endTime;
}
InterfaceConservativePhase InterfaceConservativePhase::advance(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceConservativePhase& previous,
    const InterfaceStepTime& time,const std::vector<FaceFlux>& fluxes)
{
    previous.validate(mesh,domain);
    validateStep(previous.progress_,time);
    auto next=advanceInterval(mesh,domain,previous,fluxes);
    next.progress_={time.index,time.end};
    next.stepTime_=time;
    next.parentRevision_=previous.revision_;next.revision_=std::make_shared<const Revision>();
    return next;
}
InterfaceConservativePhase InterfaceConservativePhase::advanceInterval(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceConservativePhase& previous,
    const std::vector<FaceFlux>& fluxes)
{
    const auto faceCount=static_cast<std::size_t>(domain.nativeFaceCount());
    require(fluxes.empty() || fluxes.size()==faceCount,"flux array has wrong native layout");
    auto next=previous;
    next.meshSourceRevision_.reset();next.meshTransferReport_.reset();
    next.report_={};
    const auto nan=std::numeric_limits<long double>::quiet_NaN();
    next.totalFlux_.assign(faceCount,nan); next.liquidFlux_.assign(faceCount,nan); next.limiter_.assign(faceCount,nan);
    next.directionalFlux_.assign(faceCount,{});
    next.report_.intervals=1;
    std::vector<long double> low(faceCount,0),anti(faceCount,0);
    std::vector<DirectionalFlux> lowDirections(faceCount),highDirections(faceCount);
    std::vector<bool> bidirectional(faceCount,false);
    const auto blend=[](long double a,long double b,long double theta) {
        if (theta==1) return b;
        if (theta==0 || a==b) return a;
        return b>=a?a+theta*(b-a):b+(1-theta)*(a-b);
    };
    const auto limitedDirections=[&](int row,long double theta) {
        auto result=highDirections[row];
        result.liquidPositive=blend(lowDirections[row].liquidPositive,result.liquidPositive,theta);
        result.liquidNegative=blend(lowDirections[row].liquidNegative,result.liquidNegative,theta);
        return result;
    };
    const auto limitedNet=[&](int row,long double theta) {
        if (bidirectional[row]) {
            const auto value=limitedDirections(row,theta);return value.liquidPositive-value.liquidNegative;
        }
        return theta==1?(fluxes.empty()?0.L:fluxes[row].highOrderLiquidVolume):low[row]+theta*anti[row];
    };
    std::vector<long double> div(previous.liquid_.size(),0),out=div,lowExchange=div,lowVolume=div,highExchange=div;
    std::vector<long double> positive=div,negative=div,rplus=div,rminus=div;
    for (const auto face:domain.faces(mesh)) {
        if (face.negativeCell<0 && face.positiveCell<0) continue;
        const auto f=fluxes.empty()?FaceFlux{}:fluxes[face.nativeRow];
        require(std::isfinite(f.totalVolume) && std::isfinite(f.highOrderLiquidVolume),"nonfinite face flux");
        const auto total=f.totalVolume;
        const int row=face.nativeRow;
        const auto high=f.directional.value_or(DirectionalFlux{std::max(0.L,total),std::max(0.L,-total),
            std::max(0.L,f.highOrderLiquidVolume),std::max(0.L,-f.highOrderLiquidVolume)});
        require(std::isfinite(high.totalPositive) && std::isfinite(high.totalNegative)
            && std::isfinite(high.liquidPositive) && std::isfinite(high.liquidNegative)
            && high.totalPositive>=0 && high.totalNegative>=0
            && high.liquidPositive>=0 && high.liquidPositive<=high.totalPositive
            && high.liquidNegative>=0 && high.liquidNegative<=high.totalNegative,
            "directional high order flux is not a phase subvolume");
        require(total==high.totalPositive-high.totalNegative
            && f.highOrderLiquidVolume==high.liquidPositive-high.liquidNegative,
            "net and directional face volumes disagree");
        const auto upwind=[&](long double magnitude,int donor) {
            if (magnitude==0) return 0.L;
            const long double c=donor>=0?previous.liquid_[donor]/previous.capacity_[donor]:f.boundaryLiquidFraction;
            require(std::isfinite(c) && c>=0 && c<=1,"missing/invalid directional inflow fraction");
            return magnitude*c;
        };
        auto baseline=high;
        baseline.liquidPositive=upwind(high.totalPositive,face.negativeCell);
        baseline.liquidNegative=upwind(high.totalNegative,face.positiveCell);
        highDirections[row]=high;lowDirections[row]=baseline;bidirectional[row]=f.directional.has_value();
        low[row]=baseline.liquidPositive-baseline.liquidNegative;
        anti[row]=f.highOrderLiquidVolume-low[row]; next.totalFlux_[row]=total;
        const auto accumulate=[&](int cell,long double sign) {
            if (cell<0) return;
            const auto outward=sign*total;
            div[cell]+=outward; out[cell]+=sign>0?high.totalPositive:high.totalNegative;
            lowExchange[cell]-=sign*low[row];
            highExchange[cell]-=sign*f.highOrderLiquidVolume;
            const auto increment=-sign*anti[row];
            positive[cell]+=std::max(0.L,increment); negative[cell]+=std::max(0.L,-increment);
        };
        accumulate(face.negativeCell,1); accumulate(face.positiveCell,-1);
    }
    for (int cell:previous.cells_) {
        const auto capacity=previous.capacity_[cell];
        require(std::isfinite(div[cell]) && std::isfinite(out[cell]) && std::isfinite(lowExchange[cell]),"overflowed flux sum");
        require(std::abs(div[cell])<=eps*capacity,"total volume flux is not discretely divergence free");
        require(out[cell]<=capacity,"outflow CFL exceeds one; split the interval");
    }
    // A roundoff-size native divergence can move a pure cell one ulp outside
    // its bounds. Repair shared INTERIOR phase fluxes, never a cell inventory
    // or prescribed boundary flux. Reduce only the direction that contributes
    // to the violation, with both owners consuming the same revised value.
    // The whole correction stays within the existing cell roundoff scale;
    // real CFL/divergence/feasibility defects and exhausted propagation reject.
    const auto originalLow=lowDirections;
    for (int sweep=0;sweep<64;++sweep) {
        bool bounded=true;
        for (int cell:previous.cells_) {
            const auto value=previous.liquid_[cell]+lowExchange[cell],cap=previous.capacity_[cell];
            bounded=bounded && value>=0 && value<=cap;
            require(value>=-eps*cap && value<=cap*(1+eps),"low order violation exceeds roundoff repair budget");
        }
        if (bounded) break;
        ++next.report_.baselineRepairSweeps;
        std::vector<long double> incoming(previous.liquid_.size(),0),outgoing=incoming;
        for (const auto face:domain.faces(mesh)) if (face.negativeCell>=0 && face.positiveCell>=0) {
            const auto& d=lowDirections[face.nativeRow];
            outgoing[face.negativeCell]+=d.liquidPositive;incoming[face.positiveCell]+=d.liquidPositive;
            outgoing[face.positiveCell]+=d.liquidNegative;incoming[face.negativeCell]+=d.liquidNegative;
        }
        const auto factor=[](long double violation,long double available) {
            if (!(violation>0)) return 1.L;
            require(available>violation,"no interior flux capacity for bounded roundoff repair");
            return std::max(0.L,std::nextafter(1-violation/available,0.L));
        };
        for (int cell:previous.cells_) {
            const auto value=previous.liquid_[cell]+lowExchange[cell];
            incoming[cell]=factor(value-previous.capacity_[cell],incoming[cell]);
            outgoing[cell]=factor(-value,outgoing[cell]);
        }
        std::fill(lowExchange.begin(),lowExchange.end(),0);
        for (const auto face:domain.faces(mesh)) {
            const int a=face.negativeCell,b=face.positiveCell,row=face.nativeRow;
            if (a<0 && b<0) continue;
            auto& d=lowDirections[row];
            if (a>=0 && b>=0) {
                d.liquidPositive*=std::min(outgoing[a],incoming[b]);
                d.liquidNegative*=std::min(outgoing[b],incoming[a]);
                const auto& original=originalLow[row];
                const auto correction=std::max(original.liquidPositive-d.liquidPositive,
                    original.liquidNegative-d.liquidNegative);
                require(correction<=eps*std::min(previous.capacity_[a],previous.capacity_[b]),
                    "shared low order flux repair exceeds roundoff budget");
                next.report_.maxBaselineFluxCorrection=std::max(next.report_.maxBaselineFluxCorrection,correction);
            }
            low[row]=d.liquidPositive-d.liquidNegative;
            if (a>=0) lowExchange[a]-=low[row];
            if (b>=0) lowExchange[b]+=low[row];
        }
    }
    if (next.report_.baselineRepairSweeps) {
        std::fill(positive.begin(),positive.end(),0);std::fill(negative.begin(),negative.end(),0);
        for (const auto face:domain.faces(mesh)) {
            const int a=face.negativeCell,b=face.positiveCell,row=face.nativeRow;
            if (a<0 && b<0) continue;
            anti[row]=(fluxes.empty()?0.L:fluxes[row].highOrderLiquidVolume)-low[row];
            if (a>=0) {positive[a]+=std::max(0.L,-anti[row]);negative[a]+=std::max(0.L,anti[row]);}
            if (b>=0) {positive[b]+=std::max(0.L,anti[row]);negative[b]+=std::max(0.L,-anti[row]);}
            if (originalLow[row].liquidPositive!=lowDirections[row].liquidPositive
                || originalLow[row].liquidNegative!=lowDirections[row].liquidNegative) ++next.report_.baselineRepairedFaces;
        }
    }
    bool highAdmissible=true;
    for (int cell:previous.cells_) {
        const auto capacity=previous.capacity_[cell], value=previous.liquid_[cell];
        require(std::isfinite(div[cell]) && std::isfinite(out[cell]) && std::isfinite(lowExchange[cell]),"overflowed flux sum");
        require(std::abs(div[cell])<=eps*capacity,"total volume flux is not discretely divergence free");
        require(out[cell]<=capacity,"outflow CFL exceeds one; split the interval");
        // Use the SAME shared-face arithmetic for both baseline and final
        // inventory. Equal inflow/outflow cancels before adding pure-cell mass;
        // value*(1-out/capacity)+inflow can round ABOVE capacity even at C=1.
        lowVolume[cell]=value+lowExchange[cell];
        if (!(lowVolume[cell]>=0 && lowVolume[cell]<=capacity)) {
            std::ostringstream detail;detail<<std::setprecision(21)
                <<"[InterfaceConservativePhase] low order update lost boundedness; cell="<<cell
                <<", capacity="<<capacity<<", before="<<value<<", exchange="<<lowExchange[cell]
                <<", after="<<lowVolume[cell]<<", divergence="<<div[cell]<<", outflow="<<out[cell];
            for (const auto face:domain.faces(mesh)) if (face.negativeCell==cell || face.positiveCell==cell) {
                const auto row=face.nativeRow;const auto& d=lowDirections[row];
                detail<<"; face="<<row<<", owners="<<face.negativeCell<<'/'<<face.positiveCell
                    <<", total="<<next.totalFlux_[row]<<", low="<<low[row]
                    <<", gross="<<d.totalPositive<<'/'<<d.totalNegative
                    <<", phase="<<d.liquidPositive<<'/'<<d.liquidNegative;
            }
            throw std::invalid_argument(detail.str());
        }
        require(std::isfinite(positive[cell]) && std::isfinite(negative[cell]),"overflowed correction budget");
        rplus[cell]=positive[cell]>0 ? std::min(1.L,(capacity-lowVolume[cell])/positive[cell]) : 1;
        rminus[cell]=negative[cell]>0 ? std::min(1.L,lowVolume[cell]/negative[cell]) : 1;
        const auto high=value+highExchange[cell];
        highAdmissible &= std::isfinite(high) && high>=0 && high<=capacity;
        next.report_.liquidBefore+=value;
        next.report_.maxDivergence=std::max(next.report_.maxDivergence,std::abs(div[cell])/capacity);
        next.report_.maxOutflowCfl=std::max(next.report_.maxOutflowCfl,out[cell]/capacity);
    }
    // Retain cancellation between incoming/outgoing antidiffusive corrections.
    // Start at the high-order candidate and reduce only corrections that push
    // an actually inadmissible cell farther outside its bounds. Both owners
    // share each reduction. A small floating margin belongs to the face factor,
    // never to a clipped cell mass or a relaxed accepted bound.
    std::vector<long double> chosen(faceCount,1),trial=highExchange;
    // Correction propagates through shared faces. Sixteen sweeps prematurely
    // exhausted an otherwise convergent native-AMR case (39 sweeps), making
    // a tiny high-order bound defect trigger broad upwind/FCT diffusion.
    // This is a work cap, not a looser phase/roundoff tolerance: every accepted
    // cell still satisfies the exact bounds below, or the old fallback runs.
    constexpr int maximumCancellationSweeps=128;
    for (int sweep=0;!highAdmissible && sweep<maximumCancellationSweeps;++sweep) {
        ++next.report_.limiterSweeps;
        std::vector<long double> plus(previous.liquid_.size(),0),minus=plus;
        for (const auto face:domain.faces(mesh)) {
            if (face.negativeCell<0 && face.positiveCell<0) continue;
            const auto correction=chosen[face.nativeRow]*anti[face.nativeRow];
            if (face.negativeCell>=0) {
                plus[face.negativeCell]+=std::max(0.L,-correction);
                minus[face.negativeCell]+=std::max(0.L,correction);
            }
            if (face.positiveCell>=0) {
                plus[face.positiveCell]+=std::max(0.L,correction);
                minus[face.positiveCell]+=std::max(0.L,-correction);
            }
        }
        const auto reduction=[](long double violation,long double budget) {
            return budget>0 ? std::max(0.L,1-violation/budget)*(1-eps) : 0.L;
        };
        for (int cell:previous.cells_) {
            const auto value=previous.liquid_[cell]+trial[cell],cap=previous.capacity_[cell];
            plus[cell]=value>cap?reduction(value-cap,plus[cell]):1;
            minus[cell]=value<0?reduction(-value,minus[cell]):1;
        }
        std::fill(trial.begin(),trial.end(),0);
        for (const auto face:domain.faces(mesh)) {
            if (face.negativeCell<0 && face.positiveCell<0) continue;
            const int row=face.nativeRow;
            long double factor=1;
            if (face.negativeCell>=0) factor=std::min(factor,anti[row]>=0?minus[face.negativeCell]:plus[face.negativeCell]);
            if (face.positiveCell>=0) factor=std::min(factor,anti[row]>=0?plus[face.positiveCell]:minus[face.positiveCell]);
            chosen[row]*=factor;
            const auto value=limitedNet(row,chosen[row]);
            if (face.negativeCell>=0) trial[face.negativeCell]-=value;
            if (face.positiveCell>=0) trial[face.positiveCell]+=value;
        }
        highAdmissible=true;
        for (int cell:previous.cells_) {
            const auto value=previous.liquid_[cell]+trial[cell];
            highAdmissible &= std::isfinite(value) && value>=0 && value<=previous.capacity_[cell];
        }
    }
    // The bounded upwind/FCT construction remains a guaranteed fallback if
    // local correction propagation does not finish within the work limit.
    // Pairwise face correction: both neighbours always use the SAME theta.
    // Accumulate the net exchange before adding the inventory. In particular,
    // equal inflow/outflow in a pure cell must cancel before rounding against
    // the much larger stored mass; no post-update clipping is permitted.
    next.report_.limiterFallbacks=highAdmissible?0:1;
    std::vector<long double> exchange(previous.liquid_.size(),0);
    for (const auto face:domain.faces(mesh)) {
        if (face.negativeCell<0 && face.positiveCell<0) continue;
        const int row=face.nativeRow;
        long double theta=chosen[row];
        if (!highAdmissible) {
            if (face.negativeCell>=0) theta=std::min(theta,anti[row]>=0?rminus[face.negativeCell]:rplus[face.negativeCell]);
            if (face.positiveCell>=0) theta=std::min(theta,anti[row]>=0?rplus[face.positiveCell]:rminus[face.positiveCell]);
            if (theta<1) theta*=1-eps;
        }
        if (theta<1) ++next.report_.limitedFaces;
        const auto liquid=limitedNet(row,theta);
        next.limiter_[row]=theta; next.liquidFlux_[row]=liquid;
        if (bidirectional[row]) next.directionalFlux_[row]=limitedDirections(row,theta);
        else next.directionalFlux_[row]={highDirections[row].totalPositive,highDirections[row].totalNegative,
            std::max(0.L,liquid),std::max(0.L,-liquid)};
        const auto& directions=next.directionalFlux_[row];
        require(directions.liquidPositive>=0 && directions.liquidPositive<=directions.totalPositive
            && directions.liquidNegative>=0 && directions.liquidNegative<=directions.totalNegative,
            "limited directional phase flux lost boundedness");
        if (face.negativeCell>=0) exchange[face.negativeCell]-=liquid;
        if (face.positiveCell>=0) exchange[face.positiveCell]+=liquid;
        next.report_.boundaryOutward+=face.outwardSign*liquid;
    }
    for (int cell:previous.cells_) {
        next.liquid_[cell]=previous.liquid_[cell]+exchange[cell];
        const auto value=next.liquid_[cell], capacity=previous.capacity_[cell];
        require(std::isfinite(value) && value>=0 && value<=capacity,"limited update lost boundedness");
        next.report_.liquidAfter+=value;
    }
    next.report_.closure=next.report_.liquidAfter-next.report_.liquidBefore+next.report_.boundaryOutward;
    require(std::isfinite(next.report_.closure),"nonfinite global balance");
    return next;
}
InterfaceConservativePhase InterfaceConservativePhase::advanceSubsteps(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceConservativePhase& previous,
    const InterfaceStepTime& time,const std::vector<double>& ends,const FluxProvider& provider)
{
    previous.validate(mesh,domain);validateStep(previous.progress_,time);
    require(!ends.empty() && ends.back()==1 && static_cast<bool>(provider),"missing interval partition/provider");
    std::vector<Interval> intervals;intervals.reserve(ends.size());
    double fraction=0,begin=time.begin;
    for (double endFraction:ends) {
        require(std::isfinite(endFraction) && endFraction>fraction && endFraction<=1,"invalid interval partition");
        const auto end=endFraction==1?time.end:std::fma(endFraction,time.dt,time.begin);
        const auto dt=(endFraction-fraction)*time.dt;
        require(std::isfinite(end) && end>begin && std::isfinite(dt) && dt>0,"unrepresentable substep time");
        intervals.push_back({begin,end,dt});begin=end;fraction=endFraction;
    }
    auto next=previous;
    Report aggregate;
    std::vector<long double> total(domain.nativeFaceCount(),0),liquid=total,theta(total.size(),1);
    std::vector<DirectionalFlux> directional(total.size());
    for (std::size_t i=0;i<intervals.size();++i) {
        const auto flux=provider(i,intervals[i],next);
        // The caller can query intermediate C but does not own it. Identity
        // is rechecked after the callback in case an aliased mesh was changed.
        next.validate(mesh,domain);
        next=advanceInterval(mesh,domain,next,flux);
        if (i==0) aggregate.liquidBefore=next.report_.liquidBefore;
        aggregate.liquidAfter=next.report_.liquidAfter;
        aggregate.boundaryOutward+=next.report_.boundaryOutward;
        aggregate.maxDivergence=std::max(aggregate.maxDivergence,next.report_.maxDivergence);
        aggregate.maxOutflowCfl=std::max(aggregate.maxOutflowCfl,next.report_.maxOutflowCfl);
        aggregate.limitedFaces+=next.report_.limitedFaces;++aggregate.intervals;
        aggregate.limiterSweeps+=next.report_.limiterSweeps;
        aggregate.limiterFallbacks+=next.report_.limiterFallbacks;
        aggregate.baselineRepairSweeps+=next.report_.baselineRepairSweeps;
        aggregate.baselineRepairedFaces+=next.report_.baselineRepairedFaces;
        aggregate.maxBaselineFluxCorrection=std::max(aggregate.maxBaselineFluxCorrection,next.report_.maxBaselineFluxCorrection);
        for (std::size_t f=0;f<total.size();++f) {
            total[f]+=next.totalFlux_[f];liquid[f]+=next.liquidFlux_[f];
            theta[f]=std::isnan(next.limiter_[f])?next.limiter_[f]:std::min(theta[f],next.limiter_[f]);
            auto& d=directional[f];const auto& v=next.directionalFlux_[f];
            d.totalPositive+=v.totalPositive;d.totalNegative+=v.totalNegative;
            d.liquidPositive+=v.liquidPositive;d.liquidNegative+=v.liquidNegative;
        }
    }
    aggregate.closure=aggregate.liquidAfter-aggregate.liquidBefore+aggregate.boundaryOutward;
    require(std::isfinite(aggregate.closure),"nonfinite physical step balance");
    for (const auto face:domain.faces(mesh)) if (face.negativeCell>=0 || face.positiveCell>=0) {
        const auto f=face.nativeRow;const auto& d=directional[f];
        require(std::isfinite(total[f]) && std::isfinite(liquid[f]) && std::isfinite(d.totalPositive)
            && std::isfinite(d.totalNegative) && std::isfinite(d.liquidPositive) && std::isfinite(d.liquidNegative),
            "physical step flux accumulation overflowed");
    }
    next.totalFlux_=std::move(total);next.liquidFlux_=std::move(liquid);next.limiter_=std::move(theta);
    next.directionalFlux_=std::move(directional);next.report_=aggregate;next.progress_={time.index,time.end};
    next.stepTime_=time;
    next.parentRevision_=previous.revision_;next.revision_=std::make_shared<const Revision>();
    return next;
}
