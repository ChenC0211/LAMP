#include "Modules/Interface/Compatibility/InterfaceParametricProjection.h"
#include <numeric>
#include <sstream>
#include <set>
#include <optional>

namespace InterfaceParametricProjection {
using namespace Quantum::Geometry::Parametric;
using V=Point;
namespace {
void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
void validateGrid(const Grid& grid){
    require(grid.nx>0&&grid.ny>0&&grid.nx<=65536&&grid.ny<=65536&&grid.size()>0&&grid.width>0&&grid.height>0
        &&std::isfinite(grid.left)&&std::isfinite(grid.bottom)&&std::isfinite(grid.width)&&std::isfinite(grid.height),"invalid parametric grid");
    require(!grid.axisymmetric||grid.left>=0,"negative RZ projection radius");
}
double remainderError(const RetiredComponent& c,const Grid& grid){
    require(c.id&&std::isfinite(c.time)&&c.inventory&&c.inventory->size()==grid.size(),"invalid retired component ledger");double error=0;
    for(std::size_t k=0;k<grid.size();k++){
        require(std::isfinite((*c.inventory)[k]),"nonfinite retired inventory");error=std::max(error,std::abs((*c.inventory)[k])/grid.capacity(k));
    }return error;
}
void validateBudget(const Candidate& candidate,const Options& options){
    const double budget=options.inventoryFractionUncertainty+options.measureFractionReserve;
    require(std::isfinite(budget)&&options.inventoryFractionUncertainty>0&&options.measureFractionReserve>0
        &&std::isfinite(options.inventoryFractionTolerance)&&options.inventoryFractionTolerance>budget,"invalid parametric inventory budget");
    require(candidate.report().maximumInventoryFractionError<=budget,"combined active/retired inventory exceeds original fraction budget");
    require(std::abs(candidate.report().totalGeometryInventoryDifference)<=options.measureFractionReserve*candidate.grid().totalCapacity(),
        "combined active/retired inventory exceeds original total-volume reserve");
}
}
int Grid::locate(Point p) const {
    if(!leaves){const int i=int(std::floor((p.x-left)/(width/nx))),j=int(std::floor((p.y-bottom)/(height/ny)));return i>=0&&i<nx&&j>=0&&j<ny?i*ny+j:-1;}
    const auto hit=leaves->locator.incident({p.x,p.y});
    for(int k=0;k<hit.count;k++){const int slot=leaves->slotByCell.at(hit.indices[k]);if(slot<0)continue;const auto& c=leaves->cells[slot];
        if(p.x>=c.left&&p.x<c.right&&p.y>=c.bottom&&p.y<c.top)return slot;}
    return -1;
}
std::vector<int> Grid::overlapping(const Rectangle& box) const {
    std::vector<int> candidates,result;
    if(leaves){for(int cell:leaves->locator.intersectingBox({box.left,box.bottom},{box.right,box.top}))
        candidates.push_back(leaves->slotByCell.at(cell));}
    else {
        const int i0=std::max(0,int(std::floor((box.left-left)*nx/width))),i1=std::min(nx,int(std::ceil((box.right-left)*nx/width)));
        const int j0=std::max(0,int(std::floor((box.bottom-bottom)*ny/height))),j1=std::min(ny,int(std::ceil((box.top-bottom)*ny/height)));
        for(int i=i0;i<i1;i++)for(int j=j0;j<j1;j++)candidates.push_back(i*ny+j);
    }
    for(int slot:candidates)if(slot>=0){const auto c=rectangle(slot);
        if(std::max(c.left,box.left)<std::min(c.right,box.right)&&std::max(c.bottom,box.bottom)<std::min(c.top,box.top))result.push_back(slot);}
    std::sort(result.begin(),result.end());return result;
}
std::vector<double> Grid::lines(int axis) const {
    std::vector<double> result;
    if(leaves){for(const auto& p:leaves->planes.at(axis))result.push_back(p.coordinate);}
    else {const int n=axis?ny:nx;const double begin=axis?bottom:left,h=(axis?height:width)/n;for(int k=0;k<=n;k++)result.push_back(begin+k*h);}
    return result;
}
ArcPrediction predictArc(const Source& source,double a,double b) {
    const double length=b-a;const auto partition=source.partition(a,b);std::array<V,3> moments{};
    std::vector<PolynomialArc<9>> panels;
    for(std::size_t i=1;i<partition.size();++i) {
        const double lo=partition[i-1],hi=partition[i],width=(hi-lo)/length,origin=(lo-a)/length;
        const auto panel=restrictedArc(source.arcs[source.index(lo)],lo,hi);panels.push_back(panel);
        const std::array<std::array<double,3>,3> test{{{{1,0,0}},{{2*origin-1,2*width,0}},{{6*origin*origin-6*origin+1,(12*origin-6)*width,6*width*width}}}};
        // Integral R_tt L_j = [R_t L_j] - integral R_t L_j'.
        // Evaluate the weak identity coefficientwise on the first-jet
        // polynomial, cancelling the constant mode algebraically. No
        // transported second derivative or velocity Hessian is sampled.
        for(int k=2;k<9;++k)for(int j=0;j<3;++j) {
            double weight=0;for(int m=0;m<=j;++m)weight+=test[j][m]/(k-1+m);
            const V firstJetCoefficient=panel.p[k]*(k/width);
            moments[j]+=firstJetCoefficient*((k-1)*(2*j+1)*weight);
        }
    }
    // Keep the first coefficient in the retained local chart. Converting it
    // to a global derivative and multiplying back loses low bits that an
    // extremely short arc's area constraint would amplify into curvature.
    Arc free{a,b,{{panels.front().p[0],panels.front().p[1]*(length/(panels.front().b-panels.front().a)),(moments[0]-moments[1]+moments[2])*.5,
        moments[1]/3-moments[2],moments[2]*.5}}};
    Arc base=free;base.p[2]=base.p[2]-free.p[4];base.p[3]+=free.p[4]*2;base.p[4]={};
    const V chord=base.increment(a,b),areaTangent=base.p[1]/30+base.p[2]/30+base.p[3]/35,areaNormal{areaTangent.y,-areaTangent.x};
    require(norm(chord)>0 && norm(areaTangent)>1e-12*norm(chord),"first-jet arc has degenerate area direction");
    V accumulated;double areaDefect=0;const auto& rule=gaussEight();
    for(const auto& panel:panels) {
        const auto freePanel=restrictedArc(free,panel.a,panel.b);std::array<V,9> defect{};defect[0]=accumulated;
        for(int k=1;k<9;++k){defect[k]=panel.p[k]-(k<5?freePanel.p[k]:V{});accumulated+=defect[k];}
        // Integrate first-jet differences from the shared left point. Large
        // absolute coordinates are never subtracted to form a small defect.
        for(int q=0;q<8;++q){const double f=.5*(1+rule.x[q]);areaDefect+=.25*rule.w[q]*cross(localValue(defect,f),localValue(panel.p,f,1)+localValue(freePanel.p,f,1));}

    }
    return {free,areaNormal,areaDefect};
}
double Candidate::liquidArea(double left,double right,double bottom,double top) const {
    require(bool(boundary_),"empty parametric candidate");
    double out=0;for(const auto& c:boundary_->components()){PhaseMeasure phase(*c.curve);out+=phase.simpleRectangle(left,right,bottom,top);}return out;
}
Quantum::Geometry::Parametric::DistanceSample Candidate::query(Point point) const {
    require(bool(boundary_),"empty parametric candidate");
    if(boundary_->components().empty()){DistanceSample out;out.value=std::numeric_limits<double>::infinity();return out;}
    auto out=Quantum::Geometry::Parametric::signedDistance(*boundary_->components().front().curve,point);
    for(std::size_t i=1;i<componentCount();i++){auto q=Quantum::Geometry::Parametric::signedDistance(*boundary_->components()[i].curve,point);
        if(q.value<out.value)out=q;else if(q.value==out.value)out.differentialAvailable=false;}return out;
}
double Candidate::liquidVolume(double left,double right,double bottom,double top) const {
    require(bool(boundary_),"empty parametric candidate");double out=0;
    for(const auto& c:boundary_->components()){PhaseMeasure phase(*c.curve);out+=grid_.axisymmetric?phase.rectanglePhaseMoments(left,right,bottom,top,true).values[0]:phase.simpleRectangle(left,right,bottom,top);}return out;
}
std::vector<double> Candidate::geometricVolumes() const {
    std::vector<double> out(grid_.size());
    for(const auto& component:boundary_->components()){PhaseMeasure phase(*component.curve);
    for(std::size_t k=0;k<grid_.size();k++){const auto c=grid_.rectangle(k);
        out[k]+=grid_.axisymmetric?phase.rectanglePhaseMoments(c.left,c.right,c.bottom,c.top,true).values[0]:phase.simpleRectangle(c.left,c.right,c.bottom,c.top);}}
    return out;
}
Candidate Candidate::component(std::size_t i) const {
    if(componentCount()==1&&retired_.empty()){require(i==0,"invalid material component");return *this;}
    Candidate out;out.grid_=grid_;out.boundary_=std::make_shared<const Boundary>(std::vector<Boundary::Component>{boundary_->components().at(i)});
    out.inventory_=componentInventories_.at(i);out.report_=componentReports_.at(i);return out;
}
Candidate combine(std::vector<Candidate> components,const std::vector<std::uint64_t>& ids){
    require(!components.empty()&&components.size()==ids.size(),"invalid material component assembly");
    const auto grid=components.front().grid();return combine(grid,std::move(components),ids);
}
Candidate combine(const Grid& grid,std::vector<Candidate> components,const std::vector<std::uint64_t>& ids,std::vector<RetiredComponent> retired){
    validateGrid(grid);require(components.size()==ids.size(),"invalid material component assembly");
    Candidate out;out.grid_=grid;out.inventory_.assign(grid.size(),0);std::vector<Boundary::Component> curves;std::set<std::uint64_t> identities;
    for(std::size_t i=0;i<components.size();i++){
        const auto& c=components[i];const auto& a=c.grid();const auto& b=out.grid_;
        require(c.componentCount()==1&&c.retired_.empty()&&a.sameLayout(b),"material components have different grids or nested ledgers");
        require(ids[i]&&identities.insert(ids[i]).second,"duplicate material identity");
        curves.push_back({ids[i],c.geometryHandle()});out.componentInventories_.push_back(c.inventory());out.componentReports_.push_back(c.report());
        for(std::size_t k=0;k<out.inventory_.size();k++)out.inventory_[k]+=c.inventory()[k];
        const auto& r=c.report();auto& s=out.report_;s.arcs+=r.arcs;s.refinements+=r.refinements;s.activeCellConstraints+=r.activeCellConstraints;s.endpointIncidences+=r.endpointIncidences;
        s.basisCoefficients+=r.basisCoefficients;s.denseBasisCoefficients+=r.denseBasisCoefficients;
        s.relativeH2=std::max(s.relativeH2,r.relativeH2);s.maximumInventoryFractionError+=r.maximumInventoryFractionError;s.totalGeometryInventoryDifference+=r.totalGeometryInventoryDifference;
        s.maximumGeometryDisplacement=std::max(s.maximumGeometryDisplacement,r.maximumGeometryDisplacement);s.minimumLongitudinalDerivative=std::min(s.minimumLongitudinalDerivative,r.minimumLongitudinalDerivative);
        s.positionJoin=std::max(s.positionJoin,r.positionJoin);s.tangentJoin=std::max(s.tangentJoin,r.tangentJoin);s.stationarityResidual=std::max(s.stationarityResidual,r.stationarityResidual);s.maximumEndpointIncidenceWidth=std::max(s.maximumEndpointIncidenceWidth,r.maximumEndpointIncidenceWidth);
    }
    for(const auto& c:retired){
        require(identities.insert(c.id).second,"retired identity is still active or duplicated");
        out.report_.maximumInventoryFractionError+=remainderError(c,grid);
        for(std::size_t k=0;k<grid.size();k++){out.inventory_[k]+=(*c.inventory)[k];out.report_.totalGeometryInventoryDifference-=(*c.inventory)[k];}
    }
    out.retired_=std::move(retired);out.boundary_=std::make_shared<const Boundary>(std::move(curves));
    if(components.size()==1&&out.retired_.empty()){out.componentInventories_.clear();out.componentReports_.clear();}
    return out;
}
Options activeOptions(const Candidate& candidate,const Options& options){
    double fixed=0,fixedGlobal=0;
    for(const auto& c:candidate.retiredComponents()){
        fixed+=remainderError(c,candidate.grid());
        fixedGlobal+=double(std::abs(std::accumulate(c.inventory->begin(),c.inventory->end(),0.L))/candidate.grid().totalCapacity());
    }
    const double budget=options.inventoryFractionUncertainty+options.measureFractionReserve;
    require(std::isfinite(budget)&&options.inventoryFractionUncertainty>0&&options.measureFractionReserve>0
        &&std::isfinite(options.inventoryFractionTolerance)&&options.inventoryFractionTolerance>budget&&fixed<=budget
        &&fixedGlobal<=options.measureFractionReserve,"retired inventory exceeds original local/global volume budget");
    if(!candidate.componentCount())return options;
    require(fixed<budget&&fixedGlobal<options.measureFractionReserve,"retired inventory exhausted the active geometry budget");
    auto result=options;const double count=double(candidate.componentCount()),remaining=budget-fixed;
    if(fixed==0&&fixedGlobal==0){
        // Preserve the existing arithmetic when there is no frozen error to
        // reserve, including the common closed-domain single-component path.
        result.inventoryFractionUncertainty/=count;result.measureFractionReserve/=count;result.inventoryFractionTolerance/=count;return result;
    }
    const double reserve=std::min(options.measureFractionReserve-fixedGlobal,remaining*(options.measureFractionReserve/budget));
    result.inventoryFractionUncertainty=(remaining-reserve)/count;result.measureFractionReserve=reserve/count;
    result.inventoryFractionTolerance=(options.inventoryFractionTolerance-fixed)/count;return result;
}
Candidate retireExterior(const Candidate& source,double time,const Options& options){
    require(std::isfinite(time),"invalid component retirement time");
    const auto& grid=source.grid();const DomainRectangle box{grid.left,grid.left+grid.width,grid.bottom,grid.bottom+grid.height};
    std::vector<Candidate> active;std::vector<std::uint64_t> ids;auto retired=source.retiredComponents();
    for(const auto& c:retired)require(time>=c.time,"component retirement moves backward in physical time");
    bool changed=false;
    for(std::size_t i=0;i<source.componentCount();i++){
        const auto& c=source.boundaryHandle()->components()[i];const auto& curve=*c.curve;bool absent=false;
        if(!curve.boundary.fixed(0)&&!curve.boundary.fixed(1)){
            const auto restriction=restrictSurface(curve,box);bool hasSurface=false;
            for(const auto& intervals:restriction.intervals)hasSurface=hasSurface||!intervals.empty();
            if(!hasSurface){PhaseMeasure phase(curve);absent=phase.rectanglePhaseMoments(box.left,box.right,box.bottom,box.top,grid.axisymmetric).phase==RectanglePhase::Empty;}
        }
        if(absent){retired.push_back({c.id,time,std::make_shared<const std::vector<double>>(source.componentInventory(i))});changed=true;}
        else {active.push_back(source.component(i));ids.push_back(c.id);}
    }
    if(!changed){validateBudget(source,options);return source;}
    auto result=combine(grid,std::move(active),ids,std::move(retired));validateBudget(result,options);
    // Regrouping active/retired ledgers can change addition order, but cannot
    // create a physical flux or consume the FV remainder.
    for(std::size_t k=0;k<grid.size();k++)require(std::abs(result.inventory()[k]-source.inventory()[k])<=64*std::numeric_limits<double>::epsilon()*grid.capacity(k),"component retirement changed independent inventory");
    return result;
}
Candidate initialize(std::vector<Curve> curves,const Grid& grid){
    std::vector<Candidate> components;std::vector<std::uint64_t> ids;
    for(auto& curve:curves){components.push_back(initialize(std::move(curve),grid));ids.push_back(ids.size()+1);}return combine(grid,std::move(components),ids);
}

Candidate initialize(Curve curve,const Grid& grid) {
    validateGrid(grid);
    curve.validateBoundary(grid.left,grid.left+grid.width,grid.bottom,grid.bottom+grid.height,grid.axisymmetric);
    require(!curve.arcs.empty() && curve.arcs.front().a==0 && curve.arcs.back().b==1,"initial parametric boundary needs one material interval");
    for(std::size_t i=0;i<curve.arcs.size();++i) {
        const auto& a=curve.arcs[i];const auto& b=curve.arcs[(i+1)%curve.arcs.size()];
        require(a.b>a.a && (i+1==curve.arcs.size() || a.b==b.a),"initial material partition does not close");
        for(const auto& p:a.p)require(std::isfinite(p.x) && std::isfinite(p.y),"nonfinite initial geometry");
        if(!curve.boundary.open()||i+1<curve.arcs.size())require(norm(a.sampleFraction(1)-b.p[0])<2e-11 && norm(a.sampleFraction(1,1)-b.sampleFraction(0,1))<2e-8,"initial first-jet joins are inconsistent");
    }
    if(curve.boundary.axis(0))sampleSurfaceRZ(curve.arcs.front(),0,true);if(curve.boundary.axis(1))sampleSurfaceRZ(curve.arcs.back(),1,true);
    if(curve.boundary.axisPair())require(curve.sample(1).y>curve.sample(0).y,"axis meridian must have liquid on its left");
    Candidate out;out.grid_=grid;out.boundary_=std::make_shared<const Boundary>(std::vector<Boundary::Component>{{1,std::make_shared<const Curve>(std::move(curve))}});out.inventory_.resize(grid.size());
    PhaseMeasure phase(out.geometry());
    for(std::size_t k=0;k<grid.size();k++) {
        const auto c=grid.rectangle(k);const double capacity=grid.capacity(k);
        const double volume=grid.axisymmetric?phase.rectanglePhaseMoments(c.left,c.right,c.bottom,c.top,true).values[0]:phase.simpleRectangle(c.left,c.right,c.bottom,c.top);
        if(!(std::isfinite(volume) && volume>=-2e-12*capacity && volume<=capacity*(1+2e-12))){
            std::ostringstream reason;reason<<std::setprecision(17)<<"initial effective-boundary measure is not physical: cell="<<k
                <<" rectangle="<<c.left<<','<<c.right<<','<<c.bottom<<','<<c.top<<" volume="<<volume<<" capacity="<<capacity;
            throw std::runtime_error(reason.str());
        }
        out.inventory_[k]=volume;
    }
    out.report_.arcs=out.geometry().arcs.size();return out;
}

Candidate transferMesh(const Candidate& source,const Grid& grid,const Options& options,MeshTransferReport& report){
    if(source.componentCount()>1||!source.retiredComponents().empty()){
        const auto local=activeOptions(source,options);
        std::vector<Candidate> components;std::vector<std::uint64_t> ids;MeshTransferReport total;
        const auto append=[&](const MeshTransferReport& r){
            if(total.overlaps.empty()){total.overlaps=r.overlaps;total.refinedSources=r.refinedSources;total.coarsenedTargets=r.coarsenedTargets;total.retainedCells=r.retainedCells;}
            total.inventoryBefore+=r.inventoryBefore;total.inventoryAfter+=r.inventoryAfter;total.inventoryBalance+=r.inventoryBalance;total.maximumFractionError+=r.maximumFractionError;total.maximumInheritedFraction+=r.maximumInheritedFraction;
        };
        for(std::size_t i=0;i<source.componentCount();i++){MeshTransferReport r;components.push_back(transferMesh(source.component(i),grid,local,r));ids.push_back(source.boundary_->components()[i].id);
            append(r);}
        auto retired=source.retiredComponents();
        for(auto& c:retired){
            // The same nested overlap remap, with identically zero geometry,
            // transfers this numerical remainder in the volume norm.
            auto ledger=combine(source.grid(),{},{});ledger.inventory_=*c.inventory;MeshTransferReport r;
            c.inventory=std::make_shared<const std::vector<double>>(transferMesh(ledger,grid,options,r).inventory());append(r);
        }
        auto result=combine(grid,std::move(components),ids,std::move(retired));validateBudget(result,options);
        result.boundary_=source.boundary_;report=std::move(total);return result;
    }
    const auto& old=source.grid();
    require(source.boundary_&&grid.size()>0&&grid.nx>0&&grid.ny>0&&old.axisymmetric==grid.axisymmetric
        &&old.left==grid.left&&old.bottom==grid.bottom&&old.width==grid.width&&old.height==grid.height,
        "parametric remap requires identical physical domains and measure");
    const double budget=options.inventoryFractionUncertainty+options.measureFractionReserve;
    require(std::isfinite(budget)&&options.inventoryFractionUncertainty>0&&options.measureFractionReserve>0
        &&options.inventoryFractionTolerance>budget,"invalid parametric remap budget");
    MeshTransferReport measured;
    std::optional<PhaseMeasure> phase;if(source.componentCount())phase.emplace(source.geometry());
    const auto geometricVolume=[&](const Grid::Rectangle& c){return !phase?0.:grid.axisymmetric?
        phase->rectanglePhaseMoments(c.left,c.right,c.bottom,c.top,true).values[0]:phase->simpleRectangle(c.left,c.right,c.bottom,c.top);};
    std::vector<long double> target(grid.size(),0),coverage(grid.size(),0);
    std::vector<int> donorCount(grid.size(),0);
    const auto contains=[](const Grid::Rectangle& a,const Grid::Rectangle& b){
        return a.left<=b.left&&a.right>=b.right&&a.bottom<=b.bottom&&a.top>=b.top;};
    for(std::size_t from=0;from<old.size();from++){
        const auto box=old.rectangle(from);const double capacity=old.capacity(from),inventory=source.inventory().at(from);
        require(std::isfinite(inventory)&&inventory>=-budget*capacity&&inventory<=(1+budget)*capacity,"remap source inventory exceeds original fraction budget");
        const auto children=grid.overlapping(box);require(!children.empty(),"parametric remap source has no target coverage");
        long double covered=0;std::vector<double> volumes(children.size()),geometry(children.size());
        for(std::size_t k=0;k<children.size();k++){
            const int to=children[k];const auto child=grid.rectangle(to);
            require(contains(box,child)||contains(child,box),"parametric remap requires nested Cartesian leaves");
            const Grid::Rectangle overlap{std::max(box.left,child.left),std::min(box.right,child.right),std::max(box.bottom,child.bottom),std::min(box.top,child.top)};
            const double volume=(overlap.right-overlap.left)*(overlap.top-overlap.bottom)*(grid.axisymmetric?std::acos(-1.)*(overlap.left+overlap.right):1.);
            volumes[k]=volume;covered+=volume;coverage[to]+=volume;donorCount[to]++;
            measured.overlaps.push_back({int(from),to,overlap,volume});
            if(children.size()>1)geometry[k]=geometricVolume(overlap);
        }
        require(std::abs(covered-capacity)<=64*std::numeric_limits<double>::epsilon()*capacity,"parametric remap donor coverage is incomplete");
        if(children.size()==1){target[children[0]]+=inventory;if(grid.capacity(children[0])==capacity)measured.retainedCells++;}
        else {
            measured.refinedSources++;const long double sum=std::accumulate(geometry.begin(),geometry.end(),0.L),defect=static_cast<long double>(inventory)-sum;
            measured.maximumInheritedFraction=std::max(measured.maximumInheritedFraction,double(std::abs(defect)/capacity));
            // Minimum sum(delta_child^2 / V_child) subject to the independent
            // parent inventory. This transfers its numerical discrepancy in
            // the same fraction norm; it is not a physical phase reseeding.
            long double assigned=0;const auto largest=std::size_t(std::max_element(volumes.begin(),volumes.end())-volumes.begin());
            for(std::size_t k=0;k<children.size();k++)if(k!=largest){const long double liquid=geometry[k]+defect*volumes[k]/covered;target[children[k]]+=liquid;assigned+=liquid;}
            target[children[largest]]+=static_cast<long double>(inventory)-assigned;
        }
    }
    Candidate result;result.grid_=grid;result.boundary_=source.boundary_;result.inventory_.resize(grid.size());
    long double before=0,after=0,geometryTotal=0;
    for(double v:source.inventory())before+=v;
    for(std::size_t k=0;k<grid.size();k++){
        const double capacity=grid.capacity(k);require(std::abs(coverage[k]-capacity)<=64*std::numeric_limits<double>::epsilon()*capacity,"parametric remap target coverage is incomplete");
        if(donorCount[k]>1)measured.coarsenedTargets++;
        const double inventory=double(target[k]);require(std::isfinite(inventory)&&inventory>=-budget*capacity&&inventory<=(1+budget)*capacity,"remapped independent inventory exceeds original fraction budget");
        result.inventory_[k]=inventory;after+=inventory;const auto c=grid.rectangle(k);
        const double geometric=geometricVolume(c);
        geometryTotal+=geometric;measured.maximumFractionError=std::max(measured.maximumFractionError,std::abs(inventory-geometric)/capacity);
    }
    require(measured.maximumFractionError<=budget,"parametric remap geometry/inventory fraction exceeds original budget");
    require(std::abs(after-before)<=64*std::numeric_limits<double>::epsilon()*grid.totalCapacity(),"parametric remap lost independent liquid");
    measured.inventoryBefore=double(before);measured.inventoryAfter=double(after);measured.inventoryBalance=double(after-before);
    result.report_=source.report_;result.report_.maximumInventoryFractionError=measured.maximumFractionError;
    result.report_.totalGeometryInventoryDifference=double(geometryTotal-after);
    report=std::move(measured);return result;
}

} // namespace InterfaceParametricProjection
