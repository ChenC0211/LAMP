#include "Modules/Interface/Compatibility/Detail/InterfaceParametricResolution.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace InterfaceResolutionPlan {
namespace {
namespace P=Quantum::Geometry::Parametric;
struct Piece {P::Arc arc;P::CurveBoundary boundary;std::uint64_t component=0;int terminal=-1;};
struct CellCurvature {double principal=0,azimuthal=0;P::Point point;bool intersects=false;};
CellCurvature cellCurvature(const QuadtreeMesh& mesh,int index,const Piece& part,bool rz,const std::vector<double>& critical){
    const auto& cell=mesh.cells(mesh.cellIndexToLabel[index]);
    const auto lower=mesh.verts(cell.cellToVert[0]).coord(),upper=mesh.verts(cell.cellToVert[2]).coord();
    const auto& arc=part.arc;auto points=critical;
    // The extrema of a curve restricted to a leaf occur at a stationary
    // point or at a leaf-face intersection (including a one-sided C1 join).
    for(int axis=0;axis<2;axis++)for(double plane:{axis?lower.y:lower.x,axis?upper.y:upper.x}){
        auto equation=arc.coordinate(axis);equation[0]-=plane;const auto roots=P::roots(equation,0,1);
        points.insert(points.end(),roots.begin(),roots.end());
    }
    std::sort(points.begin(),points.end());points.erase(std::unique(points.begin(),points.end()),points.end());
    const double scale=std::max({std::abs(lower.x),std::abs(lower.y),std::abs(upper.x),std::abs(upper.y),cell.size().x,cell.size().y});
    const double roundoff=128*std::numeric_limits<double>::epsilon()*scale;
    CellCurvature result;
    for(double t:points){auto q=P::sampleSurface(arc,t);
        const auto surface=rz?P::sampleSurfaceRZ(arc,t,part.boundary):q;
        q.point=surface.point;
        if(q.point.x<lower.x-roundoff||q.point.x>upper.x+roundoff||q.point.y<lower.y-roundoff||q.point.y>upper.y+roundoff)continue;
        const double hoop=rz?surface.curvature-q.curvature:0;
        if(!std::isfinite(hoop))throw std::invalid_argument("nonfinite local parametric hoop curvature");
        const double magnitude=std::max(std::abs(q.curvature),std::abs(hoop));
        if(!result.intersects||magnitude>result.principal){result.principal=magnitude;result.point=q.point;}
        result.azimuthal=std::max(result.azimuthal,std::abs(hoop));result.intersects=true;
    }
    return result;
}
bool axialGraph(const P::Arc& arc,const P::CurveBoundary& boundary) {
    std::array<double,4> dz{};for(int k=1;k<5;k++)dz[k-1]=k*arc.p[k].y;
    // The ready curve already owns the regular pole jet. Remove its exact
    // endpoint zero before testing strict positivity in the open interval.
    if(boundary.axis(0)&&arc.a==0){std::array<double,3> q{{dz[1],dz[2],dz[3]}};return P::range(q)[0]>0;}
    if(boundary.axis(1)&&arc.b==1){std::array<double,3> q{{dz[0],dz[0]+dz[1],dz[0]+dz[1]+dz[2]}};return P::range(q)[0]>0;}
    return P::range(dz)[0]>0;
}
void terminalRanges(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceState& state,
    const Options& o,Result& result) {
    if(o.parametricNeckRadiusCells==0)return;
    if(!domain.axisymmetric()||mesh.minCoord.x!=0||!std::isfinite(o.parametricNeckRadiusCells)
        ||o.parametricNeckRadiusCells<=0||o.parametricNeckRadiusCells>.5)
        throw std::invalid_argument("parametric terminal resolution requires the declared axisymmetric neck model");
    result.finestNeckScale=std::ldexp(std::min(mesh.cellSize.x,mesh.cellSize.y),-o.maximumLevel);
    result.terminalRadius=std::ldexp(std::max(mesh.cellSize.x,mesh.cellSize.y),-o.maximumLevel)/o.maxCurvatureCells;
    for(const auto& component:state.parametricGeometry().boundaryHandle()->components()){
        const auto& curve=*component.curve;
        if(!curve.boundary.axis(1)||(!curve.boundary.axis(0)&&!curve.boundary.fixed(0)))continue;
        std::vector<double> cuts{0,1};
        for(const auto& arc:curve.arcs){auto r=arc.coordinate(0);r[0]-=result.terminalRadius;
            for(double t:P::roots(r,0,1))cuts.push_back(arc.a+(arc.b-arc.a)*t);
            if(P::polynomial(r,0)==0)cuts.push_back(arc.a);if(P::polynomial(r,1)==0)cuts.push_back(arc.b);}
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        for(std::size_t k=1;k<cuts.size();k++){
            const double begin=cuts[k-1],end=cuts[k];
            if(!(begin<end)||curve.sample(.5*(begin+end)).x>=result.terminalRadius)continue;
            // An entirely subscale component has no resolved outer matching
            // boundary. It needs a droplet model, outside this neck closure.
            if(begin==0&&end==1)continue;
            const bool pole=begin==0||end==1;std::uint64_t birth=0;
            if(pole){const int side=begin==0?0:1;if(!curve.boundary.axis(side))continue;
                const auto event=state.parametricModelEvent(component.id,double(side));if(!event)continue;birth=*event;}
            bool admissible=true;double volume=0;const auto& rule=P::gaussEight();
            for(const auto& arc:curve.arcs){const double lo=std::max(begin,arc.a),hi=std::min(end,arc.b);if(!(lo<hi))continue;
                const auto part=P::restrictedArc(arc,lo,hi);admissible=admissible&&axialGraph(part,curve.boundary);
                for(int j=0;j<8;j++){const double t=.5*(1+rule.x[j]);const auto q=part.sampleFraction(t),d=part.sampleFraction(t,1);
                    volume+=.5*(part.b-part.a)*rule.w[j]*std::acos(-1.)*q.x*q.x*d.y;}}
            if(!admissible)continue;
            if(!(volume>0&&std::isfinite(volume)))throw std::runtime_error("invalid terminal-region geometric measure");
            result.terminalRanges.push_back({component.id,birth,begin,end,curve.sample(begin).y,curve.sample(end).y,state.progress().endTime,pole});
            result.terminalGeometricVolume+=volume;
        }
    }
}
std::vector<Piece> pieces(const InterfaceState& state,const Result& result) {
    std::vector<Piece> out;
    for(const auto& component:state.parametricGeometry().boundaryHandle()->components())for(const auto& arc:component.curve->arcs){
        std::vector<double> cuts{arc.a,arc.b};for(const auto& range:result.terminalRanges)if(range.component==component.id){
            if(range.begin>arc.a&&range.begin<arc.b)cuts.push_back(range.begin);
            if(range.end>arc.a&&range.end<arc.b)cuts.push_back(range.end);}
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        for(std::size_t k=1;k<cuts.size();k++){const double middle=.5*(cuts[k-1]+cuts[k]);int terminal=-1;
            for(std::size_t j=0;j<result.terminalRanges.size();j++){const auto& range=result.terminalRanges[j];
                if(range.component==component.id&&middle>=range.begin&&middle<=range.end){terminal=int(j);break;}}
            out.push_back({P::restrictedArc(arc,cuts[k-1],cuts[k]),component.curve->boundary,component.id,terminal});}
    }
    return out;
}
}
void parametricIndicators(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceState& state,
    const Options& o,Result& result,std::vector<unsigned char>& seed,const std::function<void(int,double)>& request,double hysteresis) {
    InterfaceParametricMotion::validateDomain(mesh,domain,state);
    if(o.terminalAxialSection)throw std::invalid_argument("parametric resolution does not use the legacy terminal neck model");
    terminalRanges(mesh,domain,state,o,result);const auto arcs=pieces(state,result);
    const auto& geometry=state.parametricGeometry();const auto& indices=InterfaceParametricMotion::cellIndices(state);
    const CartesianLeafLocator locator(mesh);const auto volumes=geometry.geometricVolumes();
    const double finestH=std::ldexp(std::max(mesh.cellSize.x,mesh.cellSize.y),-o.maximumLevel);
    const double physicalTime=state.progress().endTime;
    const auto finest=[&](int c){seed[c]=1;const int level=mesh.cells(mesh.cellIndexToLabel[c]).level();request(c,std::ldexp(1.,o.maximumLevel-level));};
    for(std::size_t slot=0;slot<indices.size();++slot){const int c=indices[slot];
        if(mesh.cells(mesh.cellIndexToLabel[c]).level()>o.maximumLevel)throw std::invalid_argument("accepted mesh exceeds resolution policy maximum level");
        const double capacity=geometry.grid().capacity(slot),fraction=volumes[slot]/capacity;
        if(std::min(fraction,1-fraction)>Quantum::Geometry::PlanarCellCut::relativeVolumeResolution)seed[c]=1;
        // Independent inventory defects are never terminal-model exceptions.
        const double difference=std::abs(volumes[slot]-geometry.inventory()[slot])/capacity;
        if(difference>o.maximumGeometryFractionDifference){seed[c]=1;++result.inventoryGeometryRequests;request(c,2);
            if(mesh.cells(mesh.cellIndexToLabel[c]).level()==o.maximumLevel)
                result.unresolved.push_back({Result::UnresolvedFeature::Kind::Inventory,c,0,0,0,0,difference,o.maximumGeometryFractionDifference,physicalTime,mesh.cells(mesh.cellIndexToLabel[c]).centerCoord(),{}});}
    }
    // Resolve each terminal branch and its enclosed liquid at the declared
    // finest scale. This is a physical sublevel region, not a list of cells
    // granted immunity after the allocator happens to exhaust its budget.
    for(const auto& range:result.terminalRanges)
        for(int c:locator.intersectingBox({0,range.lowerZ},{result.terminalRadius,range.upperZ}))finest(c);
    // A prospective split plane partitions the independent parent inventory.
    // Occupied coarse cells intersecting its finest-plane envelope must also
    // be refined; unrelated gas and other components need no artificial cut.
    if(o.parametricNeckRadiusCells>0){const double hz=std::ldexp(mesh.cellSize.y,-o.maximumLevel);
        for(const auto& site:InterfaceParametricProjection::inspectNecks(geometry,o.parametricNeckRadiusCells)){
            bool admitted=false;for(const auto& range:result.terminalRanges)if(!range.pole&&range.component==site.options.parent&&site.position>range.lowerZ&&site.position<range.upperZ)admitted=true;
            if(!admitted)continue;
            std::size_t part=0;while(geometry.boundaryHandle()->components()[part].id!=site.options.parent)++part;
            const auto& inventory=geometry.componentInventory(part);
            for(std::size_t slot=0;slot<indices.size();slot++){const auto box=geometry.grid().rectangle(slot);
                if(inventory[slot]!=0&&box.bottom<site.position+hz&&box.top>site.position-hz)finest(indices[slot]);}
        }
    }
    std::vector<P::PhaseMeasure::Bounds> bounds;for(const auto& part:arcs)bounds.push_back({P::range(part.arc.coordinate(0)),P::range(part.arc.coordinate(1))});
    for(std::size_t index=0;index<arcs.size();index++){
        const auto& part=arcs[index];const auto& arc=part.arc;const auto& box=bounds[index];
        const auto cells=locator.intersectingBox({box.x[0],box.y[0]},{box.x[1],box.y[1]});
        double gap=std::numeric_limits<double>::infinity(),modelGap=gap;
        P::Point gapPoint,opposedPoint;std::uint64_t opposedComponent=0;
        const auto critical=P::surfacePartition(arc,part.boundary,domain.axisymmetric());
        auto samples=critical;for(double t:P::gaussEight().x)samples.push_back(.5*(1+t));
        std::sort(samples.begin(),samples.end());samples.erase(std::unique(samples.begin(),samples.end()),samples.end());
        for(double t:samples){auto q=P::sampleSurface(arc,t);
            const bool pole=(part.boundary.axis(0)&&arc.a==0&&t==0)||(part.boundary.axis(1)&&arc.b==1&&t==1);
            if(pole)q.point.x=0;
            const P::Point tangent{-q.normal.y,q.normal.x};
            for(std::size_t other=0;other<arcs.size();other++){
                bool paired=false;const auto& partner=arcs[other];
                if(part.terminal>=0&&partner.terminal>=0&&part.component!=partner.component){
                    const auto& a=result.terminalRanges[part.terminal];const auto& b=result.terminalRanges[partner.terminal];
                    paired=a.pole&&b.pole&&a.birthParent!=0&&a.birthParent==b.birthParent;
                }
                const auto& b=bounds[other];const double dx=std::max({b.x[0]-q.point.x,0.,q.point.x-b.x[1]}),dy=std::max({b.y[0]-q.point.y,0.,q.point.y-b.y[1]});
                const double bound=paired?modelGap:gap;if(dx*dx+dy*dy>=bound*bound)continue;
                std::array<double,5> line{};for(int k=0;k<5;k++)line[k]=P::dot(partner.arc.p[k],tangent);line[0]-=P::dot(q.point,tangent);
                auto roots=P::roots(line,0,1);if(P::polynomial(line,0)==0)roots.push_back(0);if(P::polynomial(line,1)==0)roots.push_back(1);
                for(double s:roots){const auto r=P::sampleSurface(partner.arc,s);if(P::dot(q.normal,r.normal)>=0)continue;
                    const double distance=P::norm(r.point-q.point);if(distance==0)throw std::runtime_error("parametric refinement encountered an unresolved opposed contact");
                    if(paired)modelGap=std::min(modelGap,distance);else if(distance<gap){gap=distance;gapPoint=q.point;opposedPoint=r.point;opposedComponent=partner.component;}}
            }
        }
        for(int c:cells){
            const auto observed=cellCurvature(mesh,c,part,domain.axisymmetric(),critical);if(!observed.intersects)continue;
            const double curvature=observed.principal,azimuthal=observed.azimuthal;const auto curvaturePoint=observed.point;
            seed[c]=1;const auto size=mesh.cells(mesh.cellIndexToLabel[c]).size();const double h=std::max(size.x,size.y);
            if(part.terminal<0&&curvature*finestH>o.maxCurvatureCells)
                result.unresolved.push_back({Result::UnresolvedFeature::Kind::Curvature,c,part.component,0,arc.a,arc.b,curvature,o.maxCurvatureCells/finestH,physicalTime,{curvaturePoint.x,curvaturePoint.y},{}});
            if(gap<o.minimumOpposedDistanceCells*finestH)
                result.unresolved.push_back({Result::UnresolvedFeature::Kind::OpposedGap,c,part.component,opposedComponent,arc.a,arc.b,gap,o.minimumOpposedDistanceCells*finestH,physicalTime,{gapPoint.x,gapPoint.y},{opposedPoint.x,opposedPoint.y}});
            result.maxObservedCurvatureCells=std::max(result.maxObservedCurvatureCells,curvature*h);
            result.maxObservedAzimuthalCells=std::max(result.maxObservedAzimuthalCells,azimuthal*h);result.maxRawCurvatureCells=result.maxObservedCurvatureCells;
            if(curvature*h>o.maxCurvatureCells){++result.curvatureRequests;if(part.terminal>=0)++result.terminalCurvatureRequests;}
            if(azimuthal*h>o.maxCurvatureCells)++result.azimuthalRequests;
            if(part.terminal>=0){finest(c);result.maxTerminalCurvatureCells=std::max(result.maxTerminalCurvatureCells,curvature*h);}
            else request(c,curvature*h/(hysteresis*o.maxCurvatureCells));
            if(gap<o.minimumOpposedDistanceCells*h)++result.gapRequests;
            if(std::isfinite(gap))request(c,o.minimumOpposedDistanceCells*h/(hysteresis*gap));
            if(modelGap<o.minimumOpposedDistanceCells*h)++result.terminalGapRequests;
        }
    }
}
bool admitsNeck(const QuadtreeMesh& mesh,const InterfaceState& state,const Result& plan,int maximumLevel,
    const InterfaceParametricProjection::NeckSite& site) {
    state.validateMesh(mesh);bool admitted=false;
    for(const auto& range:plan.terminalRanges)if(!range.pole&&range.time==state.progress().endTime&&range.component==site.options.parent
        &&site.position>range.lowerZ&&site.position<range.upperZ)admitted=true;
    if(!admitted)return false;
    const auto& grid=state.parametricGeometry().grid();const auto& indices=InterfaceParametricMotion::cellIndices(state);int count=0;
    for(std::size_t slot=0;slot<grid.size();slot++){const auto box=grid.rectangle(slot);
        if(box.left==0&&(box.top==site.options.section||box.bottom==site.options.section)){
            if(mesh.cells(mesh.cellIndexToLabel[indices[slot]]).level()!=maximumLevel)
                throw std::logic_error("reviewed terminal neck lost its finest native support");++count;
        }}
    if(count!=2)throw std::logic_error("reviewed terminal neck has no unique pair of axis cells");
    return true;
}
}
