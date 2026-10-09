#include "Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"
#include "Modules/Interface/Compatibility/InterfaceNeckCap.h"
#include "Modules/Interface/Compatibility/InterfaceSurfaceTransfer.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceStateInternal.h"
#include "Modules/Interface/Compatibility/InterfacePlicTransport.h"
#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Geometry/CubicGraphCellMeasure.h"
#include "Core/Geometry/ActiveLevelSetTransport.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianLeafLevelIndex.h"
#include "Core/Geometry/InterfaceGeometryResolution.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceHeightMoments.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceSymmetryStencil.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <sstream>
#include <iomanip>
#include <exception>

namespace
{
namespace Cut=Quantum::Geometry::PlanarCellCut;
void require(bool value,const char* message)
{ if (!value) throw std::invalid_argument(std::string("[InterfaceConservativeGeometry] ")+message); }
std::vector<InterfaceComponentInventory> transportComponents(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceState& previous,const InterfaceConservativePhase& next,
    const InterfaceTopology& topology,const InterfaceStepTime* time,InterfaceConservativeGeometryReport& report)
{
    const auto& components=previous.componentInventories();
    if(components.empty())return {};
    require(time && next.follows(previous.conservativePhase(),*time),
        "child inventory requires its actual physical phase step; zero-time remap is not implemented");
    require(components.size()==2 && topology.liquids().size()==2,
        "child inventory currently supports one retained two-child split");
    const auto& old=previous.conservativePhase().liquidVolumes();const auto& volume=next.liquidVolumes();
    const int count=mesh.cellIndexNum();
    for(const auto& component:components)require(component.liquidId!=0
        && component.liquidVolumes.size()==old.size()
        && std::any_of(topology.liquids().begin(),topology.liquids().end(),
            [&](const auto& liquid){return liquid.id==component.liquidId;}),"child inventory lost its liquid identity/layout");
    require(components[0].liquidId!=components[1].liquidId,"duplicate child inventory identity");
    struct Exchange {int donor,receiver;long double volume;};
    std::vector<Exchange> exchanges;
    std::vector<long double> denominator=old,weight(count),rhs(count),updated(count),outgoing(count);
    long double scale=0;
    const long double eps=std::numeric_limits<long double>::epsilon();
    for(int c:domain.cellIndexByRow()) {
        const auto a=components[0].liquidVolumes[c],b=components[1].liquidVolumes[c];
        require(std::isfinite(a+b) && a>=0 && b>=0 && a<=old[c] && b<=old[c]
            && std::abs(a+b-old[c])<=4*eps*old[c],"child inventories do not partition the old liquid");
        weight[c]=old[c]>0?a/old[c]:0;scale+=old[c];
    }
    const auto append=[&](int donor,int receiver,long double amount) {
        require(std::isfinite(amount) && amount>=0,"invalid accepted directional child flux");
        if(amount==0)return;
        require(donor>=0 && receiver>=0,"child inventory boundary material flux is not yet implemented");
        exchanges.push_back({donor,receiver,amount});denominator[receiver]+=amount;
        outgoing[donor]+=amount;scale+=amount;
    };
    for(const auto face:domain.faces(mesh)) {
        const auto& flux=next.directionalFluxes().at(face.nativeRow);
        append(face.negativeCell,face.positiveCell,flux.liquidPositive);
        append(face.positiveCell,face.negativeCell,flux.liquidNegative);
    }
    std::vector<bool> implicit(count,false);
    std::vector<long double> retained=components[0].liquidVolumes,storage=denominator;
    for(int c:domain.cellIndexByRow()) {
        implicit[c]=outgoing[c]>old[c];
        if(implicit[c])++report.componentImplicitCells;
        else {
            const auto remaining=old[c]-outgoing[c];
            // Multiply a bounded concentration by the remaining amount.
            // A_old*(remaining/C_old) can exceed remaining by one ulp even
            // for A_old==C_old; that order loses pure-child boundedness.
            retained[c]=remaining*weight[c];
            storage[c]=remaining;
        }
    }
    for(const auto& exchange:exchanges)if(!implicit[exchange.receiver])storage[exchange.receiver]+=exchange.volume;
    // Donors use old concentration whenever their old liquid can supply all
    // accepted outflow. This prevents unnecessary instantaneous label mixing
    // through arbitrarily long chains. Only donors requiring this step's
    // inflow use local implicit mixing on the ACCEPTED gross liquid fluxes:
    // (C_old + sum F_in) w_new = A_old + sum F_in w_donor_new.
    // Positive sums avoid subtractive denominators; every iterate stays in
    // [0,1] by monotonic arithmetic. No clipping or geometric relabeling.
    // This first-order material mixing does not claim sharp/high-order label
    // transport. Total C and its existing high-order phase flux are unchanged.
    bool converged=false;
    for(int iteration=1;iteration<=512;++iteration) {
        rhs=components[0].liquidVolumes;
        for(const auto& exchange:exchanges)rhs[exchange.receiver]+=exchange.volume*weight[exchange.donor];
        long double change=0;
        for(int c:domain.cellIndexByRow()) {
            updated[c]=implicit[c]?rhs[c]/denominator[c]:weight[c];
            require(std::isfinite(updated[c]) && updated[c]>=0 && updated[c]<=1,
                "implicit child mixing lost boundedness");
            change=std::max(change,std::abs(updated[c]-weight[c]));
        }
        weight.swap(updated);report.componentTransportIterations=iteration;
        if(change<=16*eps){converged=true;break;}
    }
    require(converged,"implicit child inventory transport did not converge");
    auto result=components;std::array<long double,2> massChange{};
    rhs=retained;
    for(const auto& exchange:exchanges)rhs[exchange.receiver]+=exchange.volume*weight[exchange.donor];
    for(int c:domain.cellIndexByRow()) {
        const auto concentration=storage[c]>0?rhs[c]/storage[c]:0;
        require(std::isfinite(concentration) && concentration>=0 && concentration<=1,
            "transported child concentration lost boundedness");
        result[0].liquidVolumes[c]=volume[c]*concentration;
        result[1].liquidVolumes[c]=volume[c]-result[0].liquidVolumes[c];
        const auto residual=result[0].liquidVolumes[c]+outgoing[c]*weight[c]
            -components[0].liquidVolumes[c]-(rhs[c]-retained[c]);
        report.componentRowResidual=std::max(report.componentRowResidual,std::abs(residual));
        for(int k=0;k<2;++k)massChange[k]+=result[k].liquidVolumes[c]-components[k].liquidVolumes[c];
    }
    report.componentMassResidual=std::max(std::abs(massChange[0]),std::abs(massChange[1]));
    // Covers the finite solve residual plus the existing phase-update and
    // native reduction roundoff, without correcting either C or child mass.
    require(report.componentMassResidual<=1024*eps*scale && report.componentRowResidual<=1024*eps*scale,
        "child material transport exceeds its arithmetic conservation budget");
    return result;
}
struct Piece
{
    int cell=-1;
    long double cx=0,cy=0,offset=0;
    Vector2D normal;
    std::array<Cut::Point,2> ends;
    long double planeValue(Vector2D p) const
    {return normal.x*(static_cast<long double>(p.x)-cx)+normal.y*(static_cast<long double>(p.y)-cy)-offset;}
    long double distance(Vector2D p) const
    {
        const long double x=static_cast<long double>(p.x)-cx,y=static_cast<long double>(p.y)-cy;
        const auto a=ends[0],b=ends[1];const auto dx=b.x-a.x,dy=b.y-a.y;
        const auto t=std::clamp(((x-a.x)*dx+(y-a.y)*dy)/(dx*dx+dy*dy),0.L,1.L);
        return std::hypot(x-a.x-t*dx,y-a.y-t*dy);
    }
};
struct ContactColumn {Vector2D point;int cell;};
// Reconstruct a local graph from actual column volumes. This shares the
// moment-fit service with height curvature, with separate stencil and branch
// checks: pure ends, one phase transition, connected mixed cells and coherent
// predictor normals. It does not borrow stale accepted geometry as new C.
class HeightReconstruction
{
    const QuadtreeMesh& mesh;
    const NativeLeafDomain& domain;
    const InterfaceConservativePhase& inventory;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> predictor;
    const std::vector<Piece>& pieces;
    std::vector<ContactColumn> contacts;
    std::vector<LevelSetInterfacePanels::SymmetryPlane> sections;
    CartesianLeafLevelIndex grid;
    std::vector<std::array<int,2>> logical;
    std::vector<int> pieceByCell;
    struct Chart
    {
        InterfaceHeightMoments::Fit fit;
        int axis=0,sign=1;
        long double tlo=0,thi=0;
        double radialWidth=0;
    };
    struct GraphChart
    {
        Vector2D origin,normal,tangent;
        Eigen::Vector4d coefficients;
        double h=0,lo=0,hi=0,rms=0,maxResidual=0;bool cubic=false;
    };
    // Each support owns one fixed fit. Re-selecting columns at the query
    // point would differentiate different functions at neighboring vertices.
    mutable std::vector<std::optional<Chart>> charts;
    mutable std::vector<bool> chartsPrepared;
    mutable std::vector<std::optional<GraphChart>> graphCharts;
    mutable std::vector<bool> graphsPrepared;
    mutable std::vector<int> axisLiquid; // 0 unknown, +/-1 certified/not certified.
    static double component(Vector2D p,int axis) {return axis?p.y:p.x;}
    static bool chartJet(const Chart& chart,Vector2D point,
        Quantum::Geometry::ActiveHermitePatch::NodalDof& node)
    {
        const auto height=chart.fit.evaluate(component(point,1-chart.axis));
        if (!height.valid) return false;
        const auto q=1/std::sqrt(1+height.first*height.first),dq=-height.first*height.second*q*q*q;
        if(chart.radialWidth>0 && point.x==0) {
            // Keep the ordinary value/derivative at r=h, but impose evenness
            // at r=0 without moving a subcell radial zero. With T=R/h, the
            // resulting radial Hermite polynomial is
            // h*(t-T)*(-t*t+2*t+2*T)/(1+2*T), t=r/h.
            // Its second factor is positive on [0,1] for R>0. It has exactly
            // the desired root R, is monotone for r>0, and matches r-R in
            // value/first derivative at h. Only the axis nodal jet changes.
            const auto r=height.height,h=static_cast<long double>(chart.radialWidth),denominator=h+2*r;
            const auto value=-2*r*r/denominator;
            const auto derivative=-4*r*(h+r)/(denominator*denominator);
            node={static_cast<double>(chart.sign*value*q),0,
                static_cast<double>(chart.sign*(derivative*height.first*q+value*dq)),0};
            return std::isfinite(node.value) && std::isfinite(node.dy);
        }
        const auto delta=component(point,chart.axis)-height.height;
        const auto value=chart.sign*delta*q,dt=chart.sign*(-height.first*q+delta*dq),dn=chart.sign*q,mixed=chart.sign*dq;
        if (!std::isfinite(value) || !std::isfinite(dt) || !std::isfinite(dn) || !std::isfinite(mixed)) return false;
        node={static_cast<double>(value),static_cast<double>(chart.axis?dt:dn),
            static_cast<double>(chart.axis?dn:dt),static_cast<double>(mixed)};
        return true;
    }
    static bool graphJet(const GraphChart& chart,Vector2D point,
        Quantum::Geometry::ActiveHermitePatch::NodalDof& node)
    {
        const auto& c=chart.coefficients;const auto normal=chart.normal,tangent=chart.tangent;
        const auto offset=point-chart.origin;const double h=chart.h;
        const double t=(offset.x*tangent.x+offset.y*tangent.y)/h,n=(offset.x*normal.x+offset.y*normal.y)/h;
        if (t<chart.lo || t>chart.hi) return false;
        const auto height=c[0]+c[1]*t+c[2]*t*t+c[3]*t*t*t,
            first=c[1]+2*c[2]*t+3*c[3]*t*t,second=2*c[2]+6*c[3]*t,third=6*c[3];
        const auto delta=n-height,q=1/std::sqrt(1+first*first);
        if (q<.9 || std::abs(delta*second)>.5*(1+first*first)) return false;
        const auto dq=-first*second*q*q*q,ddq=-(second*second+first*third)*q*q*q+3*first*first*second*second*q*q*q*q*q;
        const auto dt=-first*q+delta*dq,dn=q,dtt=(-second*q-2*first*dq+delta*ddq)/h,dtn=dq/h;
        const Quantum::Geometry::ActiveHermitePatch::NodalDof proposed{h*delta*q,tangent.x*dt+normal.x*dn,tangent.y*dt+normal.y*dn,
            tangent.x*tangent.y*dtt+(tangent.x*normal.y+normal.x*tangent.y)*dtn};
        if (!std::isfinite(proposed.value) || !std::isfinite(proposed.dx) || !std::isfinite(proposed.dy) || !std::isfinite(proposed.dxy)) return false;
        node=proposed;return true;
    }
    bool attempt(Vector2D point,const Piece& anchor,int axis,int radius,
        Quantum::Geometry::ActiveHermitePatch::NodalDof& node,Vector2D* surfaceNormal,int stencilColumns=4,
        std::optional<Quantum::Geometry::QuadraticCellCut::Shape>* curve=nullptr,Chart* retained=nullptr,
        bool projectDonor=false) const
    {
        const auto centre=logical[anchor.cell];const int tangent=1-axis;
        const auto level=grid.address(anchor.cell).level;const auto extent=grid.extent(level);
        const int nt=extent[tangent],nn=extent[axis];
        const auto h=grid.spacing(level);const auto hn=component(h,axis),ht=component(h,tangent);
        const auto cell=[&](int t,int n,int normalAxis){return normalAxis?grid.cell(level,t,n):grid.cell(level,n,t);};
        const auto direction=component(anchor.normal,axis);if (direction==0) return false;
        const int sign=direction>0?1:-1;
        const auto& host=mesh.cells(mesh.cellIndexToLabel[anchor.cell]);
        const bool redundant=stencilColumns==7;
        int first=centre[tangent]-(component(point,tangent)<component(host.centerCoord(),tangent)?2:1),count=4;
        if (stencilColumns==3) {first=centre[tangent]-1;count=3;}
        if (stencilColumns==5) {first=centre[tangent]-2;count=5;}
        if (redundant) {first=centre[tangent]-3;count=7;}
        const auto reflection=InterfaceSymmetryStencil::select(mesh.minCoord,mesh.maxCoord,tangent,nt,
            first,count,domain.axisymmetric(),sections);
        const bool even=reflection.even;
        // Endpoint/volume projection does not constrain the derivative at
        // the symmetry axis. Retain the parity-constrained reconstruction.
        if(projectDonor && even) return false;
        if (even && redundant) return false;
        const ContactColumn* contact=nullptr;
        for (const auto& candidate:contacts) {
            if(grid.address(candidate.cell).level!=level) continue;
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[candidate.cell]);
            const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
            const auto wall=component(candidate.point,tangent);
            const int index=logical[candidate.cell][tangent];
            int shifted=first;
            if (wall==component(lo,tangent) && first<index) shifted=index;
            if (wall==component(hi,tangent) && first+count>index+1) shifted=index+1-count;
            if (shifted!=first) {
                if (contact || even) return false;
                contact=&candidate;first=shifted;
            }
        }
        if (first<0 || first+count>nt) return false;
        const int lower=std::max(0,centre[axis]-radius),upper=std::min(nn-1,centre[axis]+radius);
        for(int t=first;t<first+count;++t) for(int n=lower;n<=upper;++n)
            if(cell(t,n,axis)<0) return false;
        std::vector<InterfaceHeightMoments::Column> columns;std::vector<int> selected;
        for (int i=first;i<first+count;++i) {
            int lo,hi;
            if (!InterfaceHeightMoments::activeSpan(centre[axis],lower,upper,[&](int row){
                return domain.cellRowByIndex()[cell(i,row,axis)]>=0;
            },lo,hi)) return false;
            const auto phaseAt=[&](int row) {const int c=cell(i,row,axis);
                return InterfaceHeightMoments::columnPhase(inventory.liquidVolumes()[c],inventory.capacities()[c]);};
            const int lowerPhase=phaseAt(lo);
            if(lowerPhase!=-sign) {
                // An attached RZ neck can fill only part of the first radial
                // cell. The physical axis then replaces a pure-liquid end
                // cell, but only when the ENTIRE predictor axis trace is
                // strictly liquid. An axis pole/gas gap does not qualify.
                if(!domain.axisymmetric() || axis!=0 || sign!=1 || lo!=0
                    || mesh.minCoord.x!=0 || lowerPhase!=0) return false;
                const int c=cell(i,lo,axis);auto& certified=axisLiquid[c];
                if(certified==0) {
                    const auto& axisCell=mesh.cells(mesh.cellIndexToLabel[c]);
                    const auto a=mesh.verts(axisCell.cellToVert[0]).coord();
                    const auto b=mesh.verts(axisCell.cellToVert[3]).coord();
                    certified=-1;
                    if(predictor->negative(a) && predictor->negative(b)) {
                        const LevelSetMappedSegment trace(predictor,a,b);
                        if(trace.roots().empty() && !trace.ranges().empty()
                            && std::all_of(trace.ranges().begin(),trace.ranges().end(),[](const auto& range){return range.negative;}))
                            certified=1;
                    }
                }
                if(certified!=1) return false;
            }
            if(phaseAt(hi)!=sign) return false;
            int phase=0;long double volume=0;
            for (int j=lo;j<=hi;++j) {
                const int c=cell(i,j,axis);const int next=1+sign*phaseAt(j);
                if (next<phase) return false;phase=next;volume+=inventory.liquidVolumes()[c];
                if (pieceByCell[c]>=0) {
                    const auto& p=pieces[pieceByCell[c]];
                    if (sign*component(p.normal,axis)<=0 || p.normal.x*anchor.normal.x+p.normal.y*anchor.normal.y<=.5) return false;
                    selected.push_back(c);
                }
            }
            const auto& lowerCell=mesh.cells(mesh.cellIndexToLabel[cell(i,lo,axis)]);
            const auto& upperCell=mesh.cells(mesh.cellIndexToLabel[cell(i,hi,axis)]);
            const auto lower=mesh.verts(lowerCell.cellToVert[0]).coord(),upper=mesh.verts(upperCell.cellToVert[2]).coord();
            columns.push_back({component(lower,tangent),component(upper,tangent),component(lower,axis),component(upper,axis),volume});
        }
        if (selected.empty() || std::find(selected.begin(),selected.end(),anchor.cell)==selected.end()) return false;
        if (contact && std::find(selected.begin(),selected.end(),contact->cell)==selected.end()) return false;
        std::vector<int> visited{anchor.cell};
        for (std::size_t k=0;k<visited.size();++k) for (int c:selected)
            if (std::find(visited.begin(),visited.end(),c)==visited.end()
                && std::abs(logical[c][0]-logical[visited[k]][0])<=1 && std::abs(logical[c][1]-logical[visited[k]][1])<=1)
                visited.push_back(c);
        if (visited.size()!=selected.size()) return false;
        const auto t=component(point,tangent),n=component(point,axis);
        if (t<columns.front().tlo || t>columns.back().thi) return false;
        const auto fit=InterfaceHeightMoments::fit(columns,even?reflection.origin:contact?component(contact->point,tangent):t,
            contact?component(contact->point,axis):n,ht,hn,axis,sign,domain.axisymmetric(),even,redundant?5:0,contact!=nullptr);
        auto height=fit.evaluate(t);if (!height.valid) return false;
        if(projectDonor) {
            const auto& column=columns.at(centre[tangent]-first);
            const auto a=fit.evaluate(column.tlo),b=fit.evaluate(column.thi);
            if(!a.valid || !b.valid) return false;
            const Cut::Rectangle support=axis?Cut::Rectangle{{static_cast<double>(column.tlo),static_cast<double>(column.nlo)},
                {static_cast<double>(column.thi),static_cast<double>(column.nhi)}}:
                Cut::Rectangle{{static_cast<double>(column.nlo),static_cast<double>(column.tlo)},
                {static_cast<double>(column.nhi),static_cast<double>(column.thi)}};
            const auto projected=Quantum::Geometry::QuadraticCellCut::projectEndpointVolume(support,axis,sign,
                a.height,b.height,column.volume,domain.axisymmetric());
            if(!projected) return false;
            const auto delta=t-(column.tlo+column.thi)/2;
            height={(column.nlo+column.nhi)/2+projected->offset+projected->shape.slope*delta+projected->shape.second*delta*delta/2,
                projected->shape.slope+projected->shape.second*delta,projected->shape.second,true};
        }
        double scaledResidual=fit.scaledResidual*static_cast<double>(fit.scale)/hn;
        if (fit.squaredRadius) scaledResidual/=2*std::max(static_cast<double>(n),hn);
        if (scaledResidual>.025) return false;
        const long double anchorT=tangent?anchor.cy:anchor.cx,anchorN=axis?anchor.cy:anchor.cx;
        const auto planeHeight=anchorN+(anchor.offset-component(anchor.normal,tangent)*(t-anchorT))/direction;
        // A volume column estimates the front's direction even in an almost
        // pure donor left behind by that front. Its local PLIC offset must
        // still represent this cell's own C, not the column height. Requiring
        // those two positions to coincide suppresses a valid normal extension
        // and feeds off-interface predictor distortion back into tiny facets.
        // A fixed, shared column chart represents that same resolved front
        // throughout its support. A near-pure trailing cell can have its own
        // PLIC facet far below/above this front while retaining the correct
        // column moment. Rejecting that chart would instead blend a distorted
        // local facet fit into the front. Keep the column's phase-transition,
        // connectivity, orientation and moment checks; its jets still pass
        // the shared-node phase and Bernstein constraints before publication.
        // Only an isolated point fit needs agreement with its anchor plane.
        // Donor offsets continue to come from each unchanged cell volume.
        if (!surfaceNormal && !retained && std::abs(height.height-planeHeight)>.5*hn) return false;
        const auto q=1/std::sqrt(1+height.first*height.first);
        const double ntangent=static_cast<double>(-sign*height.first*q),nnormal=static_cast<double>(sign*q);
        if (ntangent*component(anchor.normal,tangent)+nnormal*direction<.8) return false;
        const Chart chart{fit,axis,sign,columns.front().tlo,columns.back().thi,
            domain.axisymmetric() && axis==0 && mesh.minCoord.x==0?hn:0};
        if (!chartJet(chart,point,node)) return false;
        if (retained) *retained=chart;
        if (surfaceNormal) *surfaceNormal=axis?Vector2D{ntangent,nnormal}:Vector2D{nnormal,ntangent};
        // A resolved single graph with bounded slope variation across one
        // donor cell. No quadratic is inferred across an unresolved branch.
        if (curve && std::abs(height.first)<=2 && std::abs(height.second)*ht<=.5)
            *curve=Quantum::Geometry::QuadraticCellCut::Shape{axis,sign,height.first,height.second};
        return true;
    }
    const Chart* fixedChart(const Piece& anchor) const
    {
        const int index=pieceByCell[anchor.cell];
        if (!chartsPrepared[index]) {
            chartsPrepared[index]=true;
            const Vector2D centre{static_cast<double>(anchor.cx+(anchor.ends[0].x+anchor.ends[1].x)/2),
                static_cast<double>(anchor.cy+(anchor.ends[0].y+anchor.ends[1].y)/2)};
            const int preferred=std::abs(anchor.normal.y)>=std::abs(anchor.normal.x)?1:0;
            Quantum::Geometry::ActiveHermitePatch::NodalDof node;Chart chart;
            // Reproduce quartic column heights before falling back to the
            // narrower cubic chart. Both interpolate their volume moments;
            // the shared nodal scalar still uses the same compact blending.
            for (int width:{5,4}) for (int axis:{preferred,1-preferred}) for (int radius:{2,3,4})
                if (attempt(centre,anchor,axis,radius,node,nullptr,width,nullptr,&chart)) {
                    charts[index]=std::move(chart);return &*charts[index];
                }
        }
        return charts[index]?&*charts[index]:nullptr;
    }
public:
    struct BlendQuality {int contributors=0,graphs=0,cubicGraphs=0;double spread=0,graphRms=0,graphResidual=0;};
    HeightReconstruction(const QuadtreeMesh& m,const NativeLeafDomain& d,const InterfaceConservativePhase& c,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> prediction,const std::vector<Piece>& p,const std::vector<ContactColumn>& ends={},
        const std::vector<LevelSetInterfacePanels::SymmetryPlane>& planes={}):mesh(m),domain(d),inventory(c),predictor(prediction),pieces(p),contacts(ends),sections(planes),grid(m),
        logical(m.cellIndexNum()),pieceByCell(m.cellIndexNum(),-1),charts(p.size()),chartsPrepared(p.size(),false),
        graphCharts(p.size()),graphsPrepared(p.size(),false),axisLiquid(m.cellIndexNum(),0)
    {
        for(int c=0;c<mesh.cellIndexNum();++c) {
            const auto address=grid.address(c);logical[c]={address.i,address.j};
        }
        for (std::size_t i=0;i<pieces.size();++i) pieceByCell[pieces[i].cell]=static_cast<int>(i);
    }
    void requireDonorSupport() const
    {
        if(domain.uniformLeaves()) return;
        // Union of both permitted five-column, radius-six donor rectangles.
        // Keep every selectable stencil unchanged from the uniform method;
        // a coarse neighbor is not grounds for silently choosing a plane.
        for(const auto& piece:pieces) {
            const auto a=grid.address(piece.cell);const auto extent=grid.extent(a.level);
            for(int di=-6;di<=6;++di) for(int dj=-6;dj<=6;++dj) {
                if(std::abs(di)>2 && std::abs(dj)>2) continue;
                const int i=a.i+di,j=a.j+dj;
                if(i<0 || j<0 || i>=extent[0] || j>=extent[1]) continue; // Physical boundary, not an AMR gap.
                if(grid.cell(a.level,i,j)<0) throw Quantum::Geometry::InterfaceGeometryResolution(
                    Quantum::Geometry::InterfaceGeometryResolution::Reason::MixedResolution,piece.cell,
                    "transport donor needs same-level five-column/radius-six support; enlarge the interface band");
            }
        }
    }
    bool evaluate(Vector2D point,const Piece& anchor,Quantum::Geometry::ActiveHermitePatch::NodalDof& node,
        Vector2D* surfaceNormal=nullptr,std::optional<Quantum::Geometry::QuadraticCellCut::Shape>* curve=nullptr,
        std::optional<InterfaceConservativeGeometry::DonorModel> donorModel={}) const
    {
        const int preferred=std::abs(anchor.normal.y)>=std::abs(anchor.normal.x)?1:0;
        // A quadratic transport graph must retain its column volume moment.
        // Truncating a quartic's derivatives to degree two introduces an
        // amplifying band even for flat, constant tangential translation.
        // Curvature/publication still use their independent higher-order
        // moment stencils; this choice changes only internal flux donors.
        if(surfaceNormal && donorModel) {
            if(*donorModel==InterfaceConservativeGeometry::DonorModel::ProjectedMoments)
                for(int axis:{preferred,1-preferred}) for(int radius:{2,3,4,5,6})
                    if(attempt(point,anchor,axis,radius,node,surfaceNormal,5,curve,nullptr,true)) return true;
            for(int axis:{preferred,1-preferred}) for(int radius:{2,3,4,5,6})
                if(attempt(point,anchor,axis,radius,node,surfaceNormal,3,curve)) return true;
            return false;
        }
        // Five exact column moments give a first-derivative symbol with the
        // correct orientation at every resolved flat-interface wavenumber.
        // The seven-column redundant quartic reverses the high-frequency
        // slope. Keep scalar publication separate: return a unit normal,
        // then invert its offset against the unchanged cell inventory.
        if (surfaceNormal) for (int radius:{2,3,4,5,6})
            if (attempt(point,anchor,preferred,radius,node,surfaceNormal,5,curve)) return true;
        for (int axis:{preferred,1-preferred}) for (int radius:{2,3,4})
            if (attempt(point,anchor,axis,radius,node,surfaceNormal,4,curve)) return true;
        return false;
    }
    // Partition of unity over fixed height/local charts of connected PLIC
    // support. Compact ellipses stay inside the selected tangent intervals.
    // Include all weight derivatives; topology/phase selection remains local.
    // Cell adjacency is not a certificate for unresolved nearby branches.
    bool blended(Vector2D point,const Piece& anchor,Quantum::Geometry::ActiveHermitePatch::NodalDof& node,
        BlendQuality& quality) const
    {
        const auto level=grid.address(anchor.cell).level;const auto spacing=grid.spacing(level);
        const long double h=std::min(spacing.x,spacing.y),radius2=4*h*h;
        std::vector<const Piece*> support;
        for (const auto& p:pieces) {
            const auto x=p.cx+(p.ends[0].x+p.ends[1].x)/2-point.x;
            const auto y=p.cy+(p.ends[0].y+p.ends[1].y)/2-point.y;
            if (x*x+y*y<radius2 && p.normal.x*anchor.normal.x+p.normal.y*anchor.normal.y>.8) {
                if(grid.address(p.cell).level!=level) return false;
                support.push_back(&p);
            }
        }
        if (support.empty()) return false;
        std::vector<int> connected{anchor.cell};
        for (std::size_t k=0;k<connected.size();++k) for (const auto* p:support)
            if (std::find(connected.begin(),connected.end(),p->cell)==connected.end()
                && std::abs(logical[p->cell][0]-logical[connected[k]][0])<=1
                && std::abs(logical[p->cell][1]-logical[connected[k]][1])<=1) connected.push_back(p->cell);
        long double w=0,wx=0,wy=0,wxy=0,n=0,nx=0,ny=0,nxy=0;
        double lo=std::numeric_limits<double>::infinity(),hi=-lo;int count=0;BlendQuality accepted;
        for (const auto* p:support) {
            if (std::find(connected.begin(),connected.end(),p->cell)==connected.end()) continue;
            const auto cx=p->cx+(p->ends[0].x+p->ends[1].x)/2,cy=p->cy+(p->ends[0].y+p->ends[1].y)/2;
            const auto* fixed=fixedChart(*p);const GraphChart* graph=nullptr;
            Vector2D tangent{1,0},normal{0,1};long double rt=0;
            if (fixed) {
                const auto ct=fixed->axis?cx:cy;
                rt=2*h;
                // A reflection section is not an interior support edge.
                if (!(fixed->fit.even && fixed->fit.origin==fixed->thi)) rt=std::min(rt,fixed->thi-ct);
                if (!(fixed->fit.even && fixed->fit.origin==fixed->tlo)) rt=std::min(rt,ct-fixed->tlo);
                if (!fixed->axis) {tangent={0,1};normal={1,0};}
            } else {
                const int index=pieceByCell[p->cell];
                if (!graphsPrepared[index]) {
                    graphsPrepared[index]=true;GraphChart candidate;
                    Quantum::Geometry::ActiveHermitePatch::NodalDof unused;
                    if (localGraph({static_cast<double>(cx),static_cast<double>(cy)},*p,unused,&candidate))
                        graphCharts[index]=std::move(candidate);
                }
                if (!graphCharts[index]) continue;
                graph=&*graphCharts[index];tangent=graph->tangent;normal=graph->normal;
                rt=std::min<long double>({2*h,-graph->lo*graph->h,graph->hi*graph->h});
            }
            if (!(rt>0)) continue;
            const auto x=point.x-cx,y=point.y-cy,t=x*tangent.x+y*tangent.y,d=x*normal.x+y*normal.y;
            const auto rt2=rt*rt,z=1-t*t/rt2-d*d/radius2;
            if (!(z>0)) continue;
            Quantum::Geometry::ActiveHermitePatch::NodalDof chart;
            if (!(fixed?chartJet(*fixed,point,chart):graphJet(*graph,point,chart))) continue;
            const auto zx=-2*(t*tangent.x/rt2+d*normal.x/radius2),zy=-2*(t*tangent.y/rt2+d*normal.y/radius2),
                zxy=-2*(tangent.x*tangent.y/rt2+normal.x*normal.y/radius2);
            const auto a=z*z*z,ax=3*z*z*zx,ay=3*z*z*zy,axy=6*z*zx*zy+3*z*z*zxy;
            w+=a;wx+=ax;wy+=ay;wxy+=axy;
            n+=a*chart.value;nx+=ax*chart.value+a*chart.dx;ny+=ay*chart.value+a*chart.dy;
            nxy+=axy*chart.value+ax*chart.dy+ay*chart.dx+a*chart.dxy;
            // Complete the compact weight's even extension at a declared
            // section. The chart itself already uses even physical column
            // moments; reflect only its weight centre, never the inventory.
            if(fixed && fixed->fit.even && std::any_of(sections.begin(),sections.end(),[&](const auto& section) {
                return section.normalAxis==1-fixed->axis && section.coordinate==fixed->fit.origin;
            })) {
                const auto mirrorX=fixed->axis?2*fixed->fit.origin-cx:cx;
                const auto mirrorY=fixed->axis?cy:2*fixed->fit.origin-cy;
                const auto mx=point.x-mirrorX,my=point.y-mirrorY;
                const auto mt=mx*tangent.x+my*tangent.y,md=mx*normal.x+my*normal.y;
                const auto mz=1-mt*mt/rt2-md*md/radius2;
                if(mz>0) {
                    const auto mzx=-2*(mt*tangent.x/rt2+md*normal.x/radius2);
                    const auto mzy=-2*(mt*tangent.y/rt2+md*normal.y/radius2);
                    const auto ma=mz*mz*mz,max=3*mz*mz*mzx,may=3*mz*mz*mzy;
                    const auto maxy=6*mz*mzx*mzy+3*mz*mz*zxy;
                    w+=ma;wx+=max;wy+=may;wxy+=maxy;
                    n+=ma*chart.value;nx+=max*chart.value+ma*chart.dx;ny+=may*chart.value+ma*chart.dy;
                    nxy+=maxy*chart.value+max*chart.dy+may*chart.dx+ma*chart.dxy;
                }
            }
            lo=std::min(lo,chart.value);hi=std::max(hi,chart.value);++count;
            if (graph) {++accepted.graphs;accepted.cubicGraphs+=graph->cubic;accepted.graphRms=std::max(accepted.graphRms,graph->rms);
                accepted.graphResidual=std::max(accepted.graphResidual,graph->maxResidual);}
        }
        if (count==0 || !(w>0) || hi-lo>.25*h) return false;
        const auto value=n/w,dx=(nx-value*wx)/w,dy=(ny-value*wy)/w;
        const auto dxy=(nxy-dx*wy-dy*wx-value*wxy)/w;
        if (!std::isfinite(value) || !std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dxy)) return false;
        node={static_cast<double>(value),static_cast<double>(dx),static_cast<double>(dy),static_cast<double>(dxy)};
        accepted.spread=(hi-lo)/static_cast<double>(h);accepted.contributors=count;quality=accepted;return true;
    }
    // Pre-publication geometric fallback. PLIC segment moments seed one
    // connected, coherently oriented quadratic graph, with cubic promotion
    // when four effective volume constraints support it. Away from reflected
    // axis support, refine against actual cell volumes, never a new inventory.
    bool localGraph(Vector2D point,const Piece& anchor,
        Quantum::Geometry::ActiveHermitePatch::NodalDof& node,GraphChart* retained=nullptr) const
    {
        const auto level=grid.address(anchor.cell).level;const auto spacing=grid.spacing(level);
        const auto h=std::max(spacing.x,spacing.y);
        const Vector2D normal=anchor.normal,tangent{-normal.y,normal.x};
        const long double ox=anchor.cx+(anchor.ends[0].x+anchor.ends[1].x)/2;
        const long double oy=anchor.cy+(anchor.ends[0].y+anchor.ends[1].y)/2;
        std::vector<const Piece*> support;
        for (const auto& p:pieces) {
            const auto x=p.cx+(p.ends[0].x+p.ends[1].x)/2-ox;
            const auto y=p.cy+(p.ends[0].y+p.ends[1].y)/2-oy;
            if (std::hypot(x,y)<=2.5*h && p.normal.x*normal.x+p.normal.y*normal.y>.5) {
                if(grid.address(p.cell).level!=level) return false;
                support.push_back(&p);
            }
        }
        std::vector<int> connected{anchor.cell};
        for (std::size_t k=0;k<connected.size();++k) for (const auto* p:support)
            if (std::find(connected.begin(),connected.end(),p->cell)==connected.end()
                && std::abs(logical[p->cell][0]-logical[connected[k]][0])<=1
                && std::abs(logical[p->cell][1]-logical[connected[k]][1])<=1) connected.push_back(p->cell);
        if (connected.size()!=support.size()) return false;
        // Axis parity adds reflected geometric samples only when this same
        // patch actually meets the symmetry axis. No second liquid is joined.
        bool meetsAxis=false;
        for (const auto* p:support) for (const auto e:p->ends) if (p->cx+e.x==0) meetsAxis=true;
        const bool reflect=domain.axisymmetric() && mesh.minCoord.x==0 && point.x<2*h
            && std::abs(normal.y)>.8 && meetsAxis;
        struct Sample {double t,n,w,meanT,meanT2,meanN,meanT3;};std::vector<Sample> samples;
        const auto add=[&](long double ax,long double ay,long double bx,long double by) {
            const auto x=(ax+bx)/2,y=(ay+by)/2;
            const auto t=((x-ox)*tangent.x+(y-oy)*tangent.y)/h;
            const auto n=((x-ox)*normal.x+(y-oy)*normal.y)/h;
            if (std::abs(t)>2.5) return;
            const auto dt=((bx-ax)*tangent.x+(by-ay)*tangent.y)/h;
            const auto dn=((bx-ax)*normal.x+(by-ay)*normal.y)/h;
            // A volume-matched PLIC segment constrains the average graph,
            // not its midpoint. Integrate quadratic moments on the segment;
            // in RZ, the first volume variation is weighted by its radius.
            // This is a linearized cell-volume condition: moving endpoints
            // and quadratic radial displacement remain higher-order errors.
            if (domain.axisymmetric() && x==0) return;
            const long double radial=domain.axisymmetric()?(bx-ax)/x:0;
            const auto meanT=t+radial*dt/12,meanT2=t*t+dt*dt/12+radial*t*dt/6;
            const auto meanN=n+radial*dn/12;
            samples.push_back({static_cast<double>(t),static_cast<double>(n),static_cast<double>(std::hypot(bx-ax,by-ay)/h),
                static_cast<double>(meanT),static_cast<double>(meanT2),static_cast<double>(meanN),
                static_cast<double>(t*t*t+t*dt*dt/4+radial*(t*t*dt/4+dt*dt*dt/80))});
        };
        for (const auto* p:support) {
            const auto ax=p->cx+p->ends[0].x,ay=p->cy+p->ends[0].y;
            const auto bx=p->cx+p->ends[1].x,by=p->cy+p->ends[1].y;
            add(ax,ay,bx,by);if (reflect && ax+bx>0) add(-ax,ay,-bx,by);
        }
        if (samples.size()<5) return false;
        const double t=((point.x-ox)*tangent.x+(point.y-oy)*tangent.y)/h;
        double lo=0,hi=0;for (const auto& s:samples) {lo=std::min(lo,s.t);hi=std::max(hi,s.t);}
        if (lo>-.75 || hi<.75 || t<lo || t>hi) return false;
        Eigen::MatrixXd matrix(samples.size(),3);Eigen::VectorXd rhs(samples.size());
        for (std::size_t i=0;i<samples.size();++i) {
            const auto& s=samples[i];const auto weight=std::sqrt(s.w);
            matrix.row(i)<<weight,weight*s.meanT,weight*s.meanT2;rhs[i]=weight*s.meanN;
        }
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(matrix);qr.setThreshold(1e-10);
        if (qr.rank()!=3) return false;
        const auto diagonal=qr.matrixR().topLeftCorner(3,3).diagonal().cwiseAbs().eval();
        if (diagonal.minCoeff()<1e-5*diagonal.maxCoeff()) return false;
        Eigen::Vector3d c=qr.solve(rhs);if (!c.allFinite()) return false;
        // The segment moments are a stable initial estimate. Fit the actual
        // liquid volumes next: curved edge intersections and the quadratic
        // radial moment must not be discarded as linearized PLIC residuals.
        // Axis-reflected support retains its existing parity-constrained fit.
        if (!reflect) {
            const auto seed=c;
            const auto volumes=[&](const Eigen::Vector3d& candidate,Eigen::VectorXd& residual,Eigen::MatrixXd& jacobian) {
                residual.resize(support.size());jacobian.resize(support.size(),3);
                for(std::size_t i=0;i<support.size();++i) {
                    const auto& cell=mesh.cells(mesh.cellIndexToLabel[support[i]->cell]);
                    const Cut::Rectangle box{mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()};
                    const auto measure=Quantum::Geometry::QuadraticCellCut::graphMeasure(box,{ox,oy},normal,h,
                        {{candidate[0],candidate[1],candidate[2]}},domain.axisymmetric());
                    const auto capacity=inventory.capacities()[support[i]->cell];
                    residual[i]=static_cast<double>((measure.volume-inventory.liquidVolumes()[support[i]->cell])/capacity);
                    for(int k=0;k<3;++k) jacobian(i,k)=static_cast<double>(measure.derivative[k]/capacity);
                }
                return residual.allFinite() && jacobian.allFinite();
            };
            Eigen::VectorXd residual;Eigen::MatrixXd jacobian;
            if (volumes(c,residual,jacobian)) for(int iteration=0;iteration<6;++iteration) {
                Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solve(jacobian);solve.setThreshold(1e-10);
                if(solve.rank()!=3) break;
                const auto diag=solve.matrixR().topLeftCorner(3,3).diagonal().cwiseAbs().eval();
                if(diag.minCoeff()<1e-5*diag.maxCoeff()) break;
                const Eigen::Vector3d delta=solve.solve(-residual);if(!delta.allFinite()) break;
                bool accepted=false;
                for(double fraction=1;fraction>=1./128;fraction/=2) {
                    const Eigen::Vector3d candidate=c+fraction*delta,change=candidate-seed;
                    // Keep the same local branch; bound the full polynomial
                    // displacement on its retained tangent interval.
                    const double bound=std::abs(change[0])+std::abs(change[1])*std::max(-lo,hi)
                        +std::abs(change[2])*std::max(lo*lo,hi*hi);
                    if(bound>.25) continue;
                    Eigen::VectorXd nextResidual;Eigen::MatrixXd nextJacobian;
                    if(!volumes(candidate,nextResidual,nextJacobian) || nextResidual.squaredNorm()>=residual.squaredNorm()) continue;
                    c=candidate;residual=std::move(nextResidual);jacobian=std::move(nextJacobian);accepted=true;break;
                }
                if(!accepted || delta.lpNorm<Eigen::Infinity>()<1e-12) break;
            }
        }
        Eigen::Vector4d polynomial(c[0],c[1],c[2],0);bool cubic=false;
        double supportWeight=0,supportWeight2=0;
        for(const auto& s:samples) {supportWeight+=s.w;supportWeight2+=s.w*s.w;}
        // Promote only with four effective, well-conditioned volume constraints.
        // The quadratic solution is the nested starting point and fallback.
        // Exact cell measures, not PLIC point residuals, select the improvement.
        if(!reflect && support.size()>=5 && supportWeight*supportWeight>=4*supportWeight2) {
            const auto seed=polynomial;
            const auto volumes=[&](const Eigen::Vector4d& candidate,Eigen::VectorXd& residual,Eigen::MatrixXd& jacobian) {
                residual.resize(support.size());jacobian.resize(support.size(),4);
                try {
                    for(std::size_t i=0;i<support.size();++i) {
                        const auto& cell=mesh.cells(mesh.cellIndexToLabel[support[i]->cell]);
                        const Cut::Rectangle box{mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()};
                        const auto measure=Quantum::Geometry::CubicGraphCellMeasure::evaluate(box,{ox,oy},normal,h,
                            {{candidate[0],candidate[1],candidate[2],candidate[3]}},domain.axisymmetric());
                        const auto capacity=inventory.capacities()[support[i]->cell];
                        residual[i]=static_cast<double>((measure.volume-inventory.liquidVolumes()[support[i]->cell])/capacity);
                        for(int j=0;j<4;++j) jacobian(i,j)=static_cast<double>(measure.derivative[j]/capacity);
                    }
                } catch(const std::invalid_argument&) {return false;}
                  catch(const std::runtime_error&) {return false;}
                return residual.allFinite() && jacobian.allFinite();
            };
            Eigen::VectorXd residual;Eigen::MatrixXd jacobian;auto candidate=seed;
            if(volumes(candidate,residual,jacobian)) for(int iteration=0;iteration<6;++iteration) {
                Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solve(jacobian);solve.setThreshold(1e-10);
                if(solve.rank()!=4) break;
                const auto diag=solve.matrixR().topLeftCorner(4,4).diagonal().cwiseAbs().eval();
                if(diag.minCoeff()<1e-5*diag.maxCoeff()) break;
                const Eigen::Vector4d delta=solve.solve(-residual);if(!delta.allFinite()) break;
                bool accepted=false;
                for(double fraction=1;fraction>=1./128;fraction/=2) {
                    const Eigen::Vector4d next=candidate+fraction*delta,change=next-seed;
                    const double extent=std::max(-lo,hi);
                    const double bound=std::abs(change[0])+extent*(std::abs(change[1])
                        +extent*(std::abs(change[2])+extent*std::abs(change[3])));
                    if(bound>.25) continue;
                    Eigen::VectorXd nextResidual;Eigen::MatrixXd nextJacobian;
                    if(!volumes(next,nextResidual,nextJacobian) || nextResidual.squaredNorm()>=residual.squaredNorm()) continue;
                    double error=0;
                    for(const auto& s:samples) error+=s.w*std::pow(next[0]+next[1]*s.meanT+next[2]*s.meanT2+next[3]*s.meanT3-s.meanN,2);
                    const GraphChart proposed{{static_cast<double>(ox),static_cast<double>(oy)},normal,tangent,next,h,lo,hi,0,0,true};
                    Quantum::Geometry::ActiveHermitePatch::NodalDof unused;
                    if(error>.075*.075*supportWeight || !graphJet(proposed,point,unused)) continue;
                    candidate=next;residual=std::move(nextResidual);jacobian=std::move(nextJacobian);accepted=true;cubic=true;break;
                }
                if(!accepted || delta.lpNorm<Eigen::Infinity>()<1e-12) break;
            }
            if(cubic) polynomial=candidate;
        }
        double weight=0,weight2=0,error=0,maxResidual=0,left=0,right=0;
        for (const auto& s:samples) {
            const auto residual=polynomial[0]+polynomial[1]*s.meanT+polynomial[2]*s.meanT2+polynomial[3]*s.meanT3-s.meanN;
            weight+=s.w;weight2+=s.w*s.w;error+=s.w*residual*residual;
            // Keep the original midpoint gap visible as a separate geometric
            // diagnostic; the accepted RMS now measures the moment equations.
            maxResidual=std::max(maxResidual,std::abs(polynomial[0]+polynomial[1]*s.t+polynomial[2]*s.t*s.t+polynomial[3]*s.t*s.t*s.t-s.n));
            if (s.t<0) left+=s.w;else if (s.t>0) right+=s.w;
        }
        // Judge the same length-weighted objective that was fitted, while
        // requiring actual support on both sides and enough effective samples.
        // A tiny facet cannot veto a supported graph through its unweighted
        // residual; its full C remains in the geometric discrepancy ledger.
        const auto rms=std::sqrt(error/weight);
        if (!std::isfinite(rms) || rms>.075 || weight*weight<3*weight2 || left<.5 || right<.5) return false;
        const GraphChart graph{{static_cast<double>(ox),static_cast<double>(oy)},normal,tangent,polynomial,h,lo,hi,rms,maxResidual,cubic};
        if (!graphJet(graph,point,node)) return false;
        if (retained) *retained=graph;
        return true;
    }
};
struct PlicReconstruction
{
    std::vector<Piece> pieces;
    std::vector<Vector2D> normals;
    InterfacePlicTransport::CurvedDonors curved;
    std::vector<int> unresolved;
    long double unresolvedVolumeDefect=0;
    std::size_t heightNormals=0;
    double maxNormalChange=0;
};
PlicReconstruction reconstructPieces(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceConservativePhase& inventory,std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> predictor,
    const std::vector<ContactColumn>& contacts={},std::optional<InterfaceConservativeGeometry::DonorModel> donorModel={},
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& sections={})
{
    inventory.validate(mesh,domain);
    const double nan=std::numeric_limits<double>::quiet_NaN();
    PlicReconstruction result;result.normals.assign(mesh.cellIndexNum(),{nan,nan});result.curved.resize(mesh.cellIndexNum());
    for (int c:domain.cellIndexByRow()) result.normals[c]={1,0};
    std::vector<Piece> initial;
    const auto construct=[&](int c,Vector2D normal,bool recordUnresolved) -> std::optional<Piece> {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        const Cut::Rectangle box{mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()};
        const auto volume=inventory.liquidVolumes()[c],capacity=inventory.capacities()[c];
        const auto defect=std::min(volume,capacity-volume);
        const auto volumeResolution=InterfaceHeightMoments::volumeResolution(capacity);
        // Below the inverse solver's own absolute volume uncertainty a
        // minority phase does not determine a trustworthy facet position.
        // Preserve that inventory (including all later fluxes), report its
        // full minority volume, and never amplify a roundoff-sized corner
        // fragment into a high-order interface with grid-sized gradients.
        if (defect<=volumeResolution) {
            if (recordUnresolved) {
                result.unresolved.push_back(c);result.unresolvedVolumeDefect+=defect;
            }
            return std::nullopt;
        }
        const auto plane=Cut::reconstruct(box,normal,volume/capacity,domain.axisymmetric());
        // C remains authoritative. A mixed inventory can be smaller than the
        // representable displacement of a plane from a cell corner. The Core
        // inverse has already bounded its actual geometric volume residual.
        // Do not invent an endpoint, snap C, or publish a zero-length segment.
        const auto ends=plane.volume==0 || plane.volume==Cut::capacity(box,domain.axisymmetric())
            ?std::nullopt:Cut::interfaceSegment(box,plane.normal,plane.offset);
        if (!ends) {
            std::ostringstream message;message<<std::setprecision(17)
                <<"unresolved PLIC segment exceeds inverse precision; cell="<<c<<", fraction="<<volume/capacity;
            throw std::invalid_argument(message.str());
        }
        Piece p;p.cell=c;p.normal=plane.normal;p.offset=plane.offset;
        p.cx=static_cast<long double>(box.lower.x)+(static_cast<long double>(box.upper.x)-box.lower.x)/2;
        p.cy=static_cast<long double>(box.lower.y)+(static_cast<long double>(box.upper.y)-box.lower.y)/2;
        p.ends=*ends;
        return p;
    };
    for (int c:domain.cellIndexByRow()) {
        const auto f=inventory.liquidVolumes()[c]/inventory.capacities()[c];
        if (f==0 || f==1) continue;
        const auto jet=predictor->scalarJetInNativeCell(c,mesh.cells(mesh.cellIndexToLabel[c]).centerCoord());
        result.normals[c]={jet.dx,jet.dy};
        if (const auto p=construct(c,result.normals[c],true)) {
            initial.push_back(*p);result.normals[c]=p->normal;
        }
    }
    // All fits see the SAME immutable initial orientations and volume columns.
    // Updating this support during the pass would make directions depend on
    // native cell traversal order. One fit per cell, no iterative fitting loop.
    const HeightReconstruction heights(mesh,domain,inventory,predictor,initial,contacts,sections);
    if(donorModel || !domain.uniformLeaves()) heights.requireDonorSupport();
    result.pieces=initial;
    for (std::size_t i=0;i<initial.size();++i) {
        const auto& piece=initial[i];Quantum::Geometry::ActiveHermitePatch::NodalDof unused;Vector2D normal;
        const auto point=mesh.cells(mesh.cellIndexToLabel[piece.cell]).centerCoord();
        if (!heights.evaluate(point,piece,unused,&normal,&result.curved[piece.cell],donorModel)) continue;
        // The graph normal is evaluated at its interface height, not the
        // gradient of the off-interface normalized graph scalar at this centre.
        const auto revised=construct(piece.cell,normal,false);
        if (!revised) continue; // Keep the already resolved predictor plane.
        result.pieces[i]=*revised;result.normals[piece.cell]=revised->normal;++result.heightNormals;
        result.maxNormalChange=std::max(result.maxNormalChange,std::hypot(revised->normal.x-piece.normal.x,revised->normal.y-piece.normal.y));
    }
    return result;
}
InterfaceConservativeGeometry::TransportDonors transportDonors(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& source,
    const std::vector<ContactColumn>& contacts={},
    InterfaceConservativeGeometry::DonorModel model=InterfaceConservativeGeometry::defaultDonorModel)
{
    using Model=InterfaceConservativeGeometry::DonorModel;
    using Donors=InterfaceConservativeGeometry::TransportDonors;
    using Status=Donors::NormalFitStatus;
    namespace Q=Quantum::Geometry::QuadraticCellCut;
    auto reconstruction=reconstructPieces(mesh,domain,source.conservativePhase(),source.geometry(),contacts,
        model==Model::NormalFrameVolumeFit?Model::ProjectedMoments:model,source.topology().symmetryPlanes());
    Donors result{std::move(reconstruction.normals),std::move(reconstruction.curved),reconstruction.heightNormals};
    if(model!=Model::NormalFrameVolumeFit) return result;
    const auto& inventory=source.conservativePhase();
    const auto& surfaces=source.topology().surfaces();
    struct Branch {int surface=-1,graph=-1,count=0;};
    std::vector<Branch> branches(mesh.cellIndexNum());
    for(std::size_t s=0;s<surfaces.size();++s) {
        const auto& lookup=surfaces[s].geometry.nativeChain().cellToGraphs();
        for(std::size_t c=0;c<lookup.size();++c) for(int g:lookup[c]) {
            auto& b=branches[c];b.surface=static_cast<int>(s);b.graph=g;++b.count;
        }
    }
    const CartesianLeafLevelIndex grid(mesh);
    InterfaceCurvature curvature(mesh,domain,source);
    const auto rectangle=[&](int c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        return Cut::Rectangle{mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()};};
    result.normalFits.resize(reconstruction.pieces.size());
    const auto fitPiece=[&](std::size_t index) {
        const auto& piece=reconstruction.pieces[index];
        const int c=piece.cell;Donors::NormalFit report;report.cell=c;
        const auto fit=[&]() {
            const auto branch=branches[c];if(branch.count!=1) return Status::MissingBranch;
            const auto address=grid.address(c);const auto extent=grid.extent(address.level);
            const auto& chain=surfaces[branch.surface].geometry.nativeChain();
            // The 3x3 objective may only consume one local, coherently
            // oriented branch. A second strand of the same closed surface
            // is not admitted merely because its component id agrees.
            std::vector<int> support;
            for(int di=-1;di<=1;++di) for(int dj=-1;dj<=1;++dj) {
                const int i=address.i+di,j=address.j+dj;
                if(i<0 || j<0 || i>=extent[0] || j>=extent[1]) return Status::BoundarySupport;
                const int neighbor=grid.cell(address.level,i,j);
                if(neighbor<0 || domain.cellRowByIndex()[neighbor]<0) return Status::BoundarySupport;
                const auto b=branches[neighbor];
                const auto v=inventory.liquidVolumes()[neighbor],cap=inventory.capacities()[neighbor];
                const bool resolved=std::min(v,cap-v)>InterfaceHeightMoments::volumeResolution(cap);
                if(b.count>1 || (resolved && b.count!=1)) return Status::AmbiguousBranch;
                if(b.count==1) {
                    if(b.surface!=branch.surface) return Status::AmbiguousBranch;
                    int distance=std::abs(b.graph-branch.graph);
                    if(chain.closed()) distance=std::min(distance,static_cast<int>(chain.graphs().size())-distance);
                    if(distance>4) return Status::AmbiguousBranch;
                    const auto a=chain.frameAt(branch.graph,.5).normal,n=chain.frameAt(b.graph,.5).normal;
                    if(a.x*n.x+a.y*n.y<=.5) return Status::AmbiguousBranch;
                }
                support.push_back(neighbor);
            }
            const auto sample=curvature.evaluate(branch.surface,branch.graph);
            if(!sample.valid() || !std::isfinite(sample.meridional)) return Status::CurvatureUnavailable;
            report.meridional=sample.meridional;
            report.regularizedCurvature=sample.regularized;
            const auto box=rectangle(c);const auto centre=mesh.cells(mesh.cellIndexToLabel[c]).centerCoord();
            const long double h=std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y);
            const auto capacity=inventory.capacities()[c],fraction=inventory.liquidVolumes()[c]/capacity;
            const double base=std::atan2(sample.normal.y,sample.normal.x);
            struct Trial {long double error;Q::Reconstruction curve;Vector2D normal;double angle;};
            const auto evaluate=[&](double angle) {
                const Vector2D n{std::cos(base+angle),std::sin(base+angle)};
                const auto curve=Q::reconstruct(box,{1,1,0,-sample.meridional,n},fraction,domain.axisymmetric());
                long double error=0;
                for(int neighbor:support) {
                    const auto measure=Q::graphMeasure(rectangle(neighbor),{centre.x,centre.y},n,h,
                        {{curve.offset/h,0,-sample.meridional*h/2}},domain.axisymmetric());
                    const auto residual=(measure.volume-inventory.liquidVolumes()[neighbor])/inventory.capacities()[neighbor];
                    error+=residual*residual;
                }
                return Trial{error,curve,n,angle};
            };
            // Bracket a local direction minimum inside +/-30 degrees, then
            // minimize the unchanged native-volume residual. A boundary
            // minimum reports insufficient orientation support; no angle is
            // clipped and no source C is ever replaced by fitted volumes.
            constexpr int intervals=8;const double bound=std::acos(-1.)/6,step=2*bound/intervals;
            std::vector<Trial> trials;for(int k=0;k<=intervals;++k) trials.push_back(evaluate(-bound+k*step));
            int best=0;for(int k=1;k<=intervals;++k) if(trials[k].error<trials[best].error)best=k;
            report.initialRms=std::sqrt(trials[intervals/2].error/support.size());
            if(best==0 || best==intervals) return Status::OrientationBoundary;
            auto accepted=trials[best];double lo=trials[best-1].angle,hi=trials[best+1].angle;
            constexpr double ratio=.6180339887498948482;
            auto left=evaluate(hi-ratio*(hi-lo)),right=evaluate(lo+ratio*(hi-lo));
            for(int k=0;k<48;++k) {
                if(left.error<accepted.error) accepted=left;if(right.error<accepted.error) accepted=right;
                if(left.error<right.error){hi=right.angle;right=left;left=evaluate(hi-ratio*(hi-lo));}
                else {lo=left.angle;left=right;right=evaluate(lo+ratio*(hi-lo));}
            }
            report.finalRms=std::sqrt(accepted.error/support.size());report.angleChange=accepted.angle;
            report.inverseResidual=accepted.curve.residual;
            result.curved[c]=accepted.curve.shape;result.normals[c]=accepted.normal;
            return Status::Fitted;
        };
        report.status=fit();result.normalFits[index]=report;
    };
    // All fits consume the same immutable source, branch lookup and eager
    // curvature knots. A task writes only its own cell and ordered report;
    // never parallelize reconstructPieces' separate lazy height-chart pass.
    const auto count=reconstruction.pieces.size();
    if(count) fitPiece(0); // Establish the strict geometry/FP contract on the caller.
    std::vector<std::exception_ptr> failures(count);
#pragma omp parallel for schedule(dynamic, 1) if(count>16)
    for(std::ptrdiff_t index=1;index<static_cast<std::ptrdiff_t>(count);++index) {
        try {fitPiece(static_cast<std::size_t>(index));}
        catch(...) {failures[index]=std::current_exception();}
    }
    for(const auto& failure:failures) if(failure) std::rethrow_exception(failure);
    return result;
}
void recordTransportDonors(InterfaceConservativeGeometryReport& report,
    const InterfaceConservativeGeometry::TransportDonors& donors)
{
    report.curvedTransportCells=std::count_if(donors.curved.begin(),donors.curved.end(),[](const auto& c){return c.has_value();});
    for(const auto& fit:donors.normalFits) {
        ++report.normalFrameDonorStatus[static_cast<std::size_t>(fit.status)];
        if(fit.status==InterfaceConservativeGeometry::TransportDonors::NormalFitStatus::Fitted) {
            report.maxNormalFrameDonorRms=std::max(report.maxNormalFrameDonorRms,fit.finalRms);
            report.regularizedDonorCells+=fit.regularizedCurvature;
        }
    }
}
void scaleSectorJets(Quantum::Geometry::ActiveHermitePatch& patch,const std::vector<double>& scale,
    const std::vector<std::uint8_t>& fixed)
{
    for (auto& sector:patch.cornerSectors) {
        const auto factor=scale[sector.vertexIndex];if (factor==1) continue;
        const auto mask=fixed.empty()?0:fixed[sector.vertexIndex];
        if (!(mask&1)) sector.dof.dx*=factor;if (!(mask&2)) sector.dof.dy*=factor;sector.dof.dxy*=factor;
    }
}
void limitPureCellJets(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceConservativePhase& inventory,
    Quantum::Geometry::ActiveHermitePatch& patch,InterfaceConservativeGeometryReport& report,
    const std::vector<std::uint8_t>& fixed={})
{
    // Each cubic Bernstein control adjacent to a corner is phi + theta*D.
    // A single shared factor for that node's derivative jet makes ALL incident
    // pure-cell controls have the inventory sign. Values and C never change;
    // mixed neighbours consume the same limited jet, preserving shared traces.
    std::vector<double> scale(patch.nodes.size(),1);
    for (int c:domain.cellIndexByRow()) {
        const auto volume=inventory.liquidVolumes()[c],capacity=inventory.capacities()[c];
        const int sign=InterfaceHeightMoments::columnPhase(volume,capacity);
        if (sign==0) continue;
        if (volume!=0 && volume!=capacity) ++report.unresolvedPhaseConstrainedCells;
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        const auto h=mesh.verts(cell.cellToVert[2]).coord()-mesh.verts(cell.cellToVert[0]).coord();
        for (int k=0;k<4;++k) {
            const int v=mesh.verts(cell.cellToVert[k]).index();
            const auto* dof=&patch.nodes[v];
            for (const auto& sector:patch.cornerSectors)
                if (sector.cellIndex==c && sector.vertexIndex==v) dof=&sector.dof;
            const auto& node=*dof;
            const int sx=(k==0 || k==3)?1:-1,sy=k<2?1:-1;
            const long double value=sign*node.value;
            require(value>=0,"pure-cell inventory contradicts a shared nodal phase");
            const long double dx=sign*sx*static_cast<long double>(h.x)*node.dx;
            const long double dy=sign*sy*static_cast<long double>(h.y)*node.dy;
            const long double dxy=sign*sx*sy*static_cast<long double>(h.x)*h.y*node.dxy;
            const auto mask=fixed.empty()?0:fixed[v];
            const long double fixedX=(mask&1)?dx:0,fixedY=(mask&2)?dy:0;
            const long double freeX=(mask&1)?0:dx,freeY=(mask&2)?0:dy;
            const auto constrain=[&](long double base,long double derivative,long double magnitude) {
                // Integer-cleared controls avoid forming 1/3 or 1/9. Reserve
                // a floating arithmetic margin before rounding the factor
                // toward zero; this is not a volume/fraction cutoff.
                require(base>=0,"fixed wall trace conflicts with pure-cell phase constraints");
                const auto lower=derivative-32*std::numeric_limits<double>::epsilon()*magnitude;
                if (lower < -base) scale[v]=std::min(scale[v],std::nextafter(static_cast<double>(base/(-lower)),0.));
            };
            constrain(3*value+fixedX,freeX,std::abs(freeX));constrain(3*value+fixedY,freeY,std::abs(freeY));
            constrain(9*value+3*fixedX+3*fixedY,3*freeX+3*freeY+dxy,
                3*std::abs(freeX)+3*std::abs(freeY)+std::abs(dxy));
        }
    }
    for (std::size_t v=0;v<scale.size();++v) if (scale[v]<1) {
        auto& node=patch.nodes[v];const auto mask=fixed.empty()?0:fixed[v];
        if (!(mask&1)) node.dx*=scale[v];if (!(mask&2)) node.dy*=scale[v];node.dxy*=scale[v];
        ++report.pureLimitedVertices;report.minPureJetScale=std::min(report.minPureJetScale,scale[v]);
    }
    scaleSectorJets(patch,scale,fixed);
}
bool recoverPurePhaseFit(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceConservativePhase& inventory,int vertex,int sign,double fallback,
    Quantum::Geometry::ActiveHermitePatch::NodalDof& proposed)
{
    // For fixed derivatives, each pure-cell Bernstein control is s*phi + D.
    // Their intersection is a half-line. Project only a phase-conflicting fit
    // onto its nearest endpoint, instead of replacing an accurate chart by
    // distance to a finite PLIC segment endpoint. No C or flux is modified.
    // Keep the existing fallback when this correction would put the node
    // farther into the pure phase; no adjustable geometric cutoff is added.
    auto node=proposed;const auto& point=mesh.verts(mesh.vertIndexToLabel[vertex]);
    if(domain.axisymmetric() && point.coord().x==0) node.dx=node.dxy=0;
    long double bound=0;bool incident=false;
    for(int label:point.vertToCell) if(label>=0) {
        const auto& cell=mesh.cells(label);const int c=cell.index();
        if(c<0 || domain.cellRowByIndex()[c]<0) continue;
        const int phase=InterfaceHeightMoments::columnPhase(inventory.liquidVolumes()[c],inventory.capacities()[c]);
        if(!phase) continue;if(phase!=sign) return false;
        const auto h=mesh.verts(cell.cellToVert[2]).coord()-mesh.verts(cell.cellToVert[0]).coord();
        for(int k=0;k<4;++k) if(mesh.verts(cell.cellToVert[k]).index()==vertex) {
            incident=true;const int sx=(k==0 || k==3)?1:-1,sy=k<2?1:-1;
            const long double dx=sign*sx*static_cast<long double>(h.x)*node.dx,
                dy=sign*sy*static_cast<long double>(h.y)*node.dy,
                dxy=sign*sx*sy*static_cast<long double>(h.x)*h.y*node.dxy;
            const auto eps=64*std::numeric_limits<double>::epsilon();
            bound=std::max({bound,(-dx+eps*std::abs(dx))/3,(-dy+eps*std::abs(dy))/3,
                (-(3*dx+3*dy+dxy)+eps*(3*std::abs(dx)+3*std::abs(dy)+std::abs(dxy)))/9});
        }
    }
    if(!incident || !std::isfinite(bound)) return false;
    // Preserve an exact zero endpoint. Replacing it by the smallest subnormal
    // would introduce an artificial arithmetic scale into exact predicates.
    const double value=bound>0?std::nextafter(static_cast<double>(bound),std::numeric_limits<double>::infinity()):0;
    if(value>sign*fallback) return false;
    node.value=sign*value;proposed=node;return true;
}
void limitMixedCellJets(const QuadtreeMesh& mesh,const std::vector<Piece>& pieces,
    Quantum::Geometry::ActiveHermitePatch& patch,InterfaceConservativeGeometryReport& report,
    const std::vector<std::uint8_t>& fixed={})
{
    std::vector<double> scale(patch.nodes.size(),1);
    for (const auto& piece:pieces) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[piece.cell]);
        const auto h=mesh.verts(cell.cellToVert[2]).coord()-mesh.verts(cell.cellToVert[0]).coord();
        const int axis=std::abs(piece.normal.y)>=std::abs(piece.normal.x)?1:0;
        const int sign=(axis?piece.normal.y:piece.normal.x)>0?1:-1;
        struct Control {int vertex;long double value,jet,magnitude;};
        const auto control=[&](int i,int j) {
            const bool right=i>=2,top=j>=2;
            const int k=top?(right?2:3):(right?1:0),v=mesh.verts(cell.cellToVert[k]).index();
            const auto& n=patch.nodes[v];const int sx=right?-1:1,sy=top?-1:1;
            const auto dx=(i==1 || i==2)?3*sx*static_cast<long double>(h.x)*n.dx:0;
            const auto dy=(j==1 || j==2)?3*sy*static_cast<long double>(h.y)*n.dy:0;
            const auto dxy=(i==1 || i==2) && (j==1 || j==2)
                ?sx*sy*static_cast<long double>(h.x)*h.y*n.dxy:0;
            const auto mask=fixed.empty()?0:fixed[v];
            const auto freeX=(mask&1)?0:dx,freeY=(mask&2)?0:dy;
            return Control{v,9*static_cast<long double>(n.value)+((mask&1)?dx:0)+((mask&2)?dy:0),
                freeX+freeY+dxy,std::abs(freeX)+std::abs(freeY)+std::abs(dxy)};
        };
        const auto a=control(0,0),b=control(axis?0:3,axis?3:0),
            c=control(axis?3:0,axis?0:3),d=control(3,3);
        if (sign*(b.value-a.value)<0 || sign*(d.value-c.value)<0) {
            ++report.nonmonotoneNodalCells;continue;
        }
        bool feasible=true;
        if (!fixed.empty()) for (int row=0;row<4;++row) for (int col=0;col<3;++col)
            if (sign*(control(axis?row:col+1,axis?col+1:row).value-
                control(axis?row:col,axis?col:row).value)<0) feasible=false;
        if (!feasible) {++report.nonmonotoneNodalCells;continue;}
        report.monotoneCells.push_back({piece.cell,axis,sign});
        // Every Bernstein derivative is a difference of adjacent controls.
        // Bound negative contributions separately: later shared-vertex or
        // pure-cell reductions must not invalidate this cell's inequality.
        for (int row=0;row<4;++row) for (int col=0;col<3;++col) {
            const auto lo=control(axis?row:col,axis?col:row),
                hi=control(axis?row:col+1,axis?col+1:row);
            const auto base=sign*(hi.value-lo.value);
            if (lo.vertex==hi.vertex) {
                const auto derivative=sign*(hi.jet-lo.jet);
                if (derivative < -base) scale[lo.vertex]=std::min(scale[lo.vertex],
                    base>0?std::nextafter(static_cast<double>(base/(-derivative)),0.):0.);
                continue;
            }
            const auto eps=32*std::numeric_limits<double>::epsilon();
            const auto negativeLo=std::max(0.L,sign*lo.jet+eps*lo.magnitude),
                negativeHi=std::max(0.L,-sign*hi.jet+eps*hi.magnitude),negative=negativeLo+negativeHi;
            if (negative>base) {
                const auto factor=base>0?std::nextafter(static_cast<double>(base/negative),0.):0.;
                if (negativeLo>0) scale[lo.vertex]=std::min(scale[lo.vertex],factor);
                if (negativeHi>0) scale[hi.vertex]=std::min(scale[hi.vertex],factor);
            }
        }
    }
    for (std::size_t v=0;v<scale.size();++v) if (scale[v]<1) {
        auto& node=patch.nodes[v];const auto mask=fixed.empty()?0:fixed[v];
        if (!(mask&1)) node.dx*=scale[v];if (!(mask&2)) node.dy*=scale[v];node.dxy*=scale[v];
        ++report.monotoneLimitedVertices;report.minMonotoneJetScale=std::min(report.minMonotoneJetScale,scale[v]);
    }
    scaleSectorJets(patch,scale,fixed);
}
struct ConservativeBoundary
{
    using Kind=InterfaceTransportDomain::Boundary;
    struct Face {
        Kind kind;Boundary side;int row,cell,phase;bool horizontal;Vector2D a,b;
        int slipPhase=0; // Pure outer slip layer; zero on axis/ports/no-slip walls.
        bool contains(Vector2D p) const {
            return horizontal?p.y==a.y && p.x>=a.x && p.x<=b.x:p.x==a.x && p.y>=a.y && p.y<=b.y;
        }
    };
    Vector2D lower,upper;bool radial=false,hasPorts=false,hasInternalWalls=false;
    std::vector<Face> faces;
    std::vector<int> phases;
    std::vector<std::uint8_t> fixed;
    std::vector<Vector2D> contacts;
    std::vector<ContactColumn> heightContacts;
    std::vector<LevelSetInterfacePanels::SymmetryPlane> sections;
    std::vector<ActiveLevelSetTransport::FixedWallCornerExtension> wetExtensions;
    const Quantum::Geometry::ActiveHermitePatch* scalar=nullptr;
    ConservativeBoundary(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceState& previous,
        const InterfaceConservativeGeometry::BoundaryInput& input)
        :lower(mesh.minCoord),upper(mesh.maxCoord),radial(domain.axisymmetric()),phases(domain.nativeFaceCount(),0),
        fixed(mesh.vertIndexNum(),0)
    {
        require(input.domain.geometry().meshIdentity()==domain.meshIdentity()
            && input.domain.geometry().axisymmetric()==domain.axisymmetric()
            && input.domain.geometry().cellIndexByRow()==domain.cellIndexByRow(),"boundary input belongs to another domain");
        scalar=std::get_if<Quantum::Geometry::ActiveHermitePatch>(&previous.scalarPatch());
        require(scalar,"boundary synchronization requires a native cubic scalar");
        sections=previous.topology().symmetryPlanes();
        for(const auto& section:sections)
            require(!radial || section.normalAxis==1,
                "RZ physical section transport supports axial planes; a finite-radius meridional reflection is not a physical planar symmetry");
        // Pure ports, the radial axis and separated outer slip layers use
        // actual native face/cell identities only. They need no uniform row
        // layout; the donor band and scalar-transition checks below still
        // enforce their existing same-level reconstruction requirements.
        // Wet corner sectors retain an explicit local chart. Fixed wall
        // endpoints use native faces; mixed transitions must remain in
        // the certified pure scalar plateau checked by prepareImpl.
        // Native scalar construction already validates the locally same-level
        // reentrant corner, its C0 seam and exact pure extension. The remote
        // coarse/fine boundary remains subject to the scalar plateau gate.
        for (const auto face:domain.faces(mesh)) {
            const auto kind=input.domain.boundary(face.nativeRow);
            const int phase=input.portPhase?input.portPhase(face.nativeRow):0;
            if ((face.negativeCell>=0)==(face.positiveCell>=0)) {
                require(kind==Kind::Interior && phase==0,"boundary input mislabels an interior or inactive face");continue;
            }
            const bool exterior=face.horizontal?(face.center.y==lower.y || face.center.y==upper.y):
                (face.center.x==lower.x || face.center.x==upper.x);
            hasInternalWalls|=!exterior;
            const bool axis=radial && !face.horizontal && face.center.x==0;
            const bool section=std::any_of(sections.begin(),sections.end(),[&](const auto& p){
                return p.normalAxis==(face.horizontal?1:0) && p.coordinate==(face.horizontal?face.center.y:face.center.x);
            });
            if(section) require(exterior && kind==Kind::Symmetry && phase==0,
                "declared section requires an impermeable symmetry law on every actual boundary face");
            if (axis) require(kind==Kind::Symmetry && phase==0,"radial axis requires symmetry without a port phase");
            else if (kind==Kind::SinglePhasePort) {
                require(exterior && (phase==-1 || phase==1) && face.area>0,
                    "single-phase port needs an explicit pure phase at an outer boundary");hasPorts=true;
            } else if (kind==Kind::Symmetry)
                require(exterior && phase==0,"slip boundary requires an outer impermeable face");
            else require(kind==Kind::StationaryWall && phase==0,"unsupported conservative boundary law");
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            if (face.horizontal?a.x>b.x:a.y>b.y) std::swap(a,b);
            const auto side=face.horizontal?(face.negativeCell<0?Boundary::south:Boundary::north):
                (face.negativeCell<0?Boundary::west:Boundary::east);
            faces.push_back({kind,side,face.nativeRow,face.negativeCell<0?face.positiveCell:face.negativeCell,phase,face.horizontal,a,b});
            if(kind==Kind::Symmetry && !axis && !section) {
                auto& slip=faces.back();slip.slipPhase=previous.topology().cellPhase()[slip.cell];
                require((slip.slipPhase==-1 || slip.slipPhase==1)
                    && InterfaceHeightMoments::columnPhase(previous.conservativePhase().liquidVolumes()[slip.cell],
                        previous.conservativePhase().capacities()[slip.cell])==slip.slipPhase,
                    "outer slip wall requires a pure separated interface layer; moving contact is unsupported");
            }
            phases[face.nativeRow]=phase;
            if (kind==Kind::StationaryWall)
                for (int label:edge.edgeToVert) fixed[mesh.verts(label).index()]|=face.horizontal?1:2;
        }
        validateGeometry(previous);
        for (const auto& face:faces) if (face.phase)
            require(previous.conservativePhase().liquidVolumes()[face.cell]==
                (face.phase<0?previous.conservativePhase().capacities()[face.cell]:0),"old port inventory is not the declared pure phase");
        for (const auto& surface:previous.topology().surfaces()) {
            const auto& chain=surface.geometry.nativeChain();if (chain.closed()) continue;
            for (bool end:{false,true}) {
                const auto kind=end?chain.endEndpointKind():chain.beginEndpointKind();
                if (kind==LevelSetInterfacePanels::EndpointKind::DeclaredBoundary) {
                    const auto p=chain.frameAt(end?static_cast<int>(chain.graphs().size())-1:0,end?1:0).point;
                    require(std::any_of(faces.begin(),faces.end(),[&](const auto& f){return f.kind==Kind::StationaryWall && f.contains(p);}),
                        "declared contact does not belong to a stationary wall");contacts.push_back(p);
                    heightContacts.push_back({p,chain.graphs()[end?chain.graphs().size()-1:0].nativeCellIndex});
                } else if(kind==LevelSetInterfacePanels::EndpointKind::SymmetryPlane) {
                    const auto p=chain.frameAt(end?static_cast<int>(chain.graphs().size())-1:0,end?1:0).point;
                    require(chain.endpointSymmetryPlane(end) && std::any_of(faces.begin(),faces.end(),[&](const auto& f){
                        return f.kind==Kind::Symmetry && f.contains(p);}),"section endpoint has no physical symmetry face");
                } else require(kind==LevelSetInterfacePanels::EndpointKind::Axis,"open endpoint has no conservative boundary identity");
            }
        }
        for (const auto& sector:scalar->cornerSectors) {
            const auto point=mesh.verts(mesh.vertIndexToLabel[sector.vertexIndex]).coord();
            const Face* wet=nullptr;
            for (const auto& face:faces) if (face.kind==Kind::StationaryWall && face.cell==sector.cellIndex && face.contains(point)) {
                require(!wet,"wet corner has multiple candidate stationary walls");wet=&face;
            }
            require(wet,"wet sector has no actual stationary wall");
            const auto& base=scalar->nodes[sector.vertexIndex];
            require(sector.dof.value==0 && base.value==0 && sector.dof.dxy==base.dxy
                && (wet->horizontal?(sector.dof.dx==0 && sector.dof.dy==base.dy):
                    (sector.dof.dy==0 && sector.dof.dx==base.dx)),"wet sector does not match the flat wall-tangent extension");
            bool bound=false;
            for (const auto& surface:previous.topology().surfaces()) {
                const auto& chain=surface.geometry.nativeChain();if (chain.closed()) continue;
                for (bool end:{false,true}) {
                    if ((end?chain.endEndpointKind():chain.beginEndpointKind())!=LevelSetInterfacePanels::EndpointKind::DeclaredBoundary) continue;
                    const int g=end?static_cast<int>(chain.graphs().size())-1:0;
                    const auto p=chain.frameAt(g,end?1:0).point;
                    if (p.x!=point.x || p.y!=point.y) continue;
                    require(!bound,"wet corner belongs to multiple interface ends");bound=true;
                }
            }
            require(bound,"wet sector is not a declared interface endpoint");
            wetExtensions.push_back({sector.cellIndex,sector.vertexIndex,sector.phaseSign,wet->horizontal});
            // The wet wall owns the substituted derivative. Applying that
            // constraint to the base cap jet would freeze the contact angle.
            fixed[sector.vertexIndex]&=~(wet->horizontal?1:2);
            require(fixed[sector.vertexIndex]!=0,"wet corner lacks its complementary stationary wall trace");
            require(previous.conservativePhase().liquidVolumes()[sector.cellIndex]==
                (sector.phaseSign<0?previous.conservativePhase().capacities()[sector.cellIndex]:0),
                "old wet-corner inventory is not its declared pure phase");
        }
    }
    void validateGeometry(const InterfaceState& state) const {
        require(state.topology().symmetryPlanes()==sections,"conservative synchronization changed symmetry declarations");
        for (const auto& face:faces) {
            const int phase=face.phase?face.phase:face.slipPhase;
            if(phase) require(state.topology().cellPhase()[face.cell]==phase,"interface reached the required pure port/slip layer");
        }
    }
    void validateInventory(const InterfaceConservativePhase& inventory,bool physicalStep=true) const {
        for (const auto& wet:wetExtensions)
            require(inventory.liquidVolumes()[wet.cellIndex]==(wet.phaseSign<0?inventory.capacities()[wet.cellIndex]:0),
                "wet-corner inventory lost its pure extension; a contact event is required");
        for (const auto& face:faces) {
            if(face.slipPhase) require(InterfaceHeightMoments::columnPhase(inventory.liquidVolumes()[face.cell],
                inventory.capacities()[face.cell])==face.slipPhase,"resolved interface inventory reached the separated slip layer");
            if (face.phase) {
                require(inventory.liquidVolumes()[face.cell]==(face.phase<0?inventory.capacities()[face.cell]:0),
                    "transported port layer is no longer pure");
            }
            if(!physicalStep)continue; // A remap has no physical-step flux ledger.
            const auto& f=inventory.directionalFluxes().at(face.row);
            if (face.phase) {
                require(f.liquidPositive==(face.phase<0?f.totalPositive:0) && f.liquidNegative==(face.phase<0?f.totalNegative:0),
                    "port phase flux disagrees with its declared history");
            } else require(f.totalPositive==0 && f.totalNegative==0 && f.liquidPositive==0 && f.liquidNegative==0
                && inventory.totalFluxes().at(face.row)==0 && inventory.liquidFluxes().at(face.row)==0,
                "impermeable boundary has nonzero gross volume or phase flux");
        }
    }
    void transportOptions(const QuadtreeMesh& mesh,ActiveLevelSetTransport::Options& options) const {
        for (const auto& wet:wetExtensions) {
            options.fixedZeroVertices.push_back(wet.vertexIndex);options.fixedWallCornerExtensions.push_back(wet);
            const auto point=mesh.verts(mesh.vertIndexToLabel[wet.vertexIndex]).coord();
            const auto center=mesh.cells(mesh.cellIndexToLabel[wet.cellIndex]).centerCoord();
            const auto coordinate=[&](Vector2D p){return wet.wallHorizontal?p.x:p.y;};
            const auto origin=coordinate(point);const bool below=coordinate(center)<origin;
            std::vector<std::uint8_t> cap(mesh.vertIndexNum(),0),extension(mesh.vertIndexNum(),0);
            for (int c=0;c<mesh.cellIndexNum();++c) if (scalar->activeCells[c])
                for (int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert) {
                    const auto& vertex=mesh.verts(label);const double x=coordinate(vertex.coord());
                    const bool seam=x==origin,wetSide=(x<origin)==below;
                    cap[vertex.index()]=seam || !wetSide;extension[vertex.index()]=seam || wetSide;
                }
            if (options.oneSidedSamples.empty()) {
                options.oneSidedSamples=std::move(cap);options.extensionSideSamples=std::move(extension);
            } else require(options.oneSidedSamples==cap && options.extensionSideSamples==extension,
                "multiple wet corners require incompatible scalar charts");
        }
    }
    void constrainSectionJets(const QuadtreeMesh& mesh,Quantum::Geometry::ActiveHermitePatch& patch) const {
        // Homogeneous Neumann data for the even scalar. Tangential value and
        // derivative remain transported; unlike a wall contact, no phi value
        // or endpoint radius is pinned. The published complete trace is then
        // certified by Core, including its mixed derivative along the edge.
        for(int v=0;v<mesh.vertIndexNum();++v) {
            auto& node=patch.nodes[v];if(!std::isfinite(node.value))continue;
            const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();
            for(const auto& section:sections) if((section.normalAxis?p.y:p.x)==section.coordinate) {
                if(section.normalAxis)node.dy=0;else node.dx=0;
                node.dxy=0;
            }
        }
    }
    void appendSectors(Quantum::Geometry::ActiveHermitePatch& patch) const {
        for (const auto& wet:wetExtensions) {
            auto dof=patch.nodes[wet.vertexIndex];
            if (wet.wallHorizontal) dof.dx=0;else dof.dy=0;
            patch.cornerSectors.push_back({wet.cellIndex,wet.vertexIndex,dof,wet.phaseSign});
        }
    }
    Vector2D velocity(Vector2D p,Vector2D u) const {
        if (!hasInternalWalls && p.x!=lower.x && p.x!=upper.x && p.y!=lower.y && p.y!=upper.y) return u;
        for (const auto& face:faces) if (face.contains(p)) {
            if (face.kind==Kind::StationaryWall) require(u.x==0 && u.y==0,"stationary wall velocity is not zero");
            if (face.kind==Kind::Symmetry) require((face.horizontal?u.y:u.x)==0,"symmetry/slip normal velocity is not zero");
        }
        return u;
    }
    void validateScalar(Boundary side,Vector2D p,double value) const {
        int phase=0;
        for (const auto& face:faces) if (face.side==side && face.contains(p)) {
            require(face.phase!=0 && (phase==0 || phase==face.phase),"scalar inflow has no unique single-phase port history");phase=face.phase;
        }
        require(phase!=0 && std::isfinite(value) && phase*value>0,"scalar inflow disagrees with the conservative port phase");
    }
};
void resolvedComponents(const InterfaceState& state,int cells,bool fixedWalls=false)
{
    require(state.topology().symmetryPlanes().empty() || fixedWalls,
        "symmetry-section transport requires explicit physical boundary laws");
    std::vector<int> ownership(cells,0);
    for (const auto& liquid:state.topology().liquids()) {
        require(fixedWalls || !liquid.touchesBoundary,"non-axis contact requires explicit fixed-wall synchronization");
        for (int c:liquid.cells) require(++ownership.at(c)==1,
            "several liquids share a cell; refine before single-plane reconstruction");
    }
}
}
void InterfaceConservativeGeometry::prepare(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& previous,InterfaceConservativePhase inventory,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> predictor,
    const InterfaceStepTime& time,InterfaceState& destination,double distanceBandCells,
    const BoundaryInput* boundaries,CurvatureDiagnostics curvatureDiagnostics,const NeckEventPolicy* neckEvents)
{
    require(previous.hasConservativePhase() && inventory.follows(previous.conservativePhase(),time),
        "inventory does not descend from this physical source step");
    prepareImpl(mesh,domain,previous,std::move(inventory),std::move(predictor),&time,destination,distanceBandCells,boundaries,curvatureDiagnostics,neckEvents);
}
void InterfaceConservativeGeometry::resynchronize(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& previous,InterfaceState& destination,double distanceBandCells,const BoundaryInput* boundaries,
    CurvatureDiagnostics curvatureDiagnostics)
{
    require(&previous!=&destination && destination.stage()==InterfaceState::Stage::Empty,"fresh zero-time destination required");
    InterfacePreparation preparation(InterfaceStateAccess::get(destination));
    require(previous.stage()==InterfaceState::Stage::Ready && previous.hasConservativePhase()
        && previous.hasTopology() && (boundaries || (previous.pinnedCoordinates().empty() && previous.cornerSectors().empty())),
        "zero-time synchronization requires a ready state and explicit contact boundary laws");
    require(previous.progress().completedSteps==previous.conservativePhase().progress().completedSteps
        && previous.progress().endTime==previous.conservativePhase().progress().endTime,
        "zero-time source geometry/inventory clocks differ");
    prepareImpl(mesh,domain,previous,previous.conservativePhase(),previous.geometry(),nullptr,destination,distanceBandCells,boundaries,curvatureDiagnostics,nullptr);
    preparation.complete(InterfaceState::Stage::Ready);
}
void InterfaceConservativeGeometry::prepareImpl(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& previous,InterfaceConservativePhase inventory,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> predictor,
    const InterfaceStepTime* time,InterfaceState& destination,double distanceBandCells,
    const BoundaryInput* boundaries,CurvatureDiagnostics curvatureDiagnostics,const NeckEventPolicy* neckEvents)
{
    require(destination.stage()==InterfaceState::Stage::Empty,"fresh destination required");
    previous.validateMesh(mesh);inventory.validate(mesh,domain);
    require(previous.hasConservativePhase() && previous.hasTopology(),"previous complete geometry and inventory required");
    require(boundaries || domain.cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "selected-domain synchronization requires explicit boundary laws");
    require(predictor && predictor->usesNativeLeaves(),"native scalar predictor required");
    InterfacePhaseGeometry::validateLayout(mesh,domain,*predictor);
    require(std::isfinite(distanceBandCells) && distanceBandCells>=2,"distance band must cover at least two cells");
    resolvedComponents(previous,mesh.cellIndexNum(),boundaries!=nullptr);
    std::vector<Vector2D> contacts;
    std::vector<std::uint8_t> fixed;
    const Quantum::Geometry::ActiveHermitePatch* wallSource=nullptr;
    std::optional<ConservativeBoundary> boundary;
    if (boundaries) {
        boundary.emplace(mesh,domain,previous,*boundaries);boundary->validateInventory(inventory,time!=nullptr);
        wallSource=boundary->scalar;fixed=boundary->fixed;contacts=boundary->contacts;
    }
    auto& target=InterfaceStateAccess::get(destination);
    InterfacePreparation preparation(target);
    InterfaceConservativeGeometryReport report;report.predictor=predictor;
    const std::vector<ContactColumn> heightContacts=boundary?boundary->heightContacts:std::vector<ContactColumn>{};
    auto reconstruction=reconstructPieces(mesh,domain,inventory,predictor,heightContacts,{},previous.topology().symmetryPlanes());
    report.heightPlicCells=reconstruction.heightNormals;report.maxPlicNormalChange=reconstruction.maxNormalChange;
    report.unresolvedPlicCells=std::move(reconstruction.unresolved);
    report.unresolvedPlicVolumeDefect=reconstruction.unresolvedVolumeDefect;
    const auto& pieces=reconstruction.pieces;
    require(!pieces.empty(),"no mixed-cell PLIC pieces; pure-face interface reconstruction is not implemented");
    const auto& interfaceCell=mesh.cells(mesh.cellIndexToLabel[pieces.front().cell]);
    const auto h=domain.uniformLeaves()?domain.spacing():interfaceCell.size();
    const double band=distanceBandCells*std::max(h.x,h.y);
    require(std::isfinite(band) && band>0,"unrepresentable distance band");report.distanceBandWidth=band;
    if(!domain.uniformLeaves()) for(const auto& piece:pieces)
        if(mesh.cells(mesh.cellIndexToLabel[piece.cell]).level()!=interfaceCell.level())
            throw Quantum::Geometry::InterfaceGeometryResolution(
                Quantum::Geometry::InterfaceGeometryResolution::Reason::MixedResolution,piece.cell,
                "synchronization requires one common interface-band level");
    report.plicSegments=pieces.size();
    for(const auto& piece:pieces) report.plicDirections.push_back({piece.cell,piece.normal});
    const HeightReconstruction heights(mesh,domain,inventory,predictor,pieces,heightContacts,previous.topology().symmetryPlanes());
    // Phase evidence belongs to cells actually incident on this vertex. A
    // nearest PLIC segment in another cell has no authority to extend its
    // liquid half-space over this point. At conflicting local PLIC traces,
    // retain the transported sign; never vote by the number of neighbours.
    std::vector<const Piece*> pieceByCell(mesh.cellIndexNum(),nullptr);
    for (const auto& piece:pieces) pieceByCell[piece.cell]=&piece;
    std::vector<unsigned char> localSides(mesh.vertIndexNum(),0),inventorySides(mesh.vertIndexNum(),0);
    for (int c:domain.cellIndexByRow()) {
        const auto volume=inventory.liquidVolumes()[c],capacity=inventory.capacities()[c];
        // PLIC and column moments already classify sub-inverse-resolution
        // minorities as geometrically pure. Shared node phase and Bernstein
        // constraints must use that same decision, or a retained tiny C can
        // create a much larger polynomial fragment in its unresolved cell.
        const int geometricPhase=InterfaceHeightMoments::columnPhase(volume,capacity);
        const unsigned char phase=geometricPhase>0?1:geometricPhase<0?2:4;
        for (int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert) {
            const auto& vertex=mesh.verts(label);const int v=vertex.index();inventorySides[v]|=phase;
            if (phase!=4) localSides[v]|=phase;
            else if (pieceByCell[c]) {
                const auto side=pieceByCell[c]->planeValue(vertex.coord());
                localSides[v]|=side>0?1:side<0?2:4;
            } else localSides[v]|=8; // Explicitly unrepresentable local plane.
        }
    }
    Quantum::Geometry::ActiveHermitePatch patch;
    patch.activeCells.assign(mesh.cellIndexNum(),0);
    for (int c:domain.cellIndexByRow()) patch.activeCells[c]=1;
    const double nan=std::numeric_limits<double>::quiet_NaN();
    patch.nodes.assign(mesh.vertIndexNum(),{nan,nan,nan,nan});
    report.closestPlicCellByVertex.assign(mesh.vertIndexNum(),-1);
    for (int v=0;v<mesh.vertIndexNum();++v) {
        if (inventorySides[v]==0) continue;
        patch.nodes[v]={};
        const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();
        long double distance=std::numeric_limits<long double>::infinity();const Piece* closest=nullptr;
        for (const auto& piece:pieces) {
            const auto d=piece.distance(p);
            if (d<distance) { distance=d;closest=&piece; }
        }
        require(closest && std::isfinite(distance),"unresolved PLIC distance");
        const auto predicted=predictor->value(p);
        require(std::isfinite(predicted),"nonfinite scalar prediction");
        const auto sides=localSides[v];const int pureSides=inventorySides[v]&3;
        const bool conflict=(sides&3)==3;
        const int sign=pureSides==1?1:pureSides==2?-1:!(sides&8) && (sides&3)==1?1:
            !(sides&8) && (sides&3)==2?-1:predicted<0?-1:1;
        auto& node=patch.nodes[v];node.value=sign*static_cast<double>(std::min(distance,static_cast<long double>(band)));
        if (pureSides==3) node.value=0; // Opposite pure cells require a shared zero.
        if (distance>=band) {
            // A local PLIC supporting plane cannot classify the far field of
            // a closed surface. Its half-space can flip across donor changes
            // even in uniform gas, creating artificial zero contours. Outside
            // the reconstructed band preserve the transported scalar phase.
            require(predicted!=0,"predictor zero outside reconstructed PLIC distance band");
            node.value=predicted<0?-band:band;
            continue;
        }
        report.closestPlicCellByVertex[v]=closest->cell;
        InterfaceConservativeGeometryReport::NodalJet trace;trace.vertex=v;
        trace.rejectedValues={nan,nan,nan};int attempt=0;
        if (conflict) ++report.localPhaseDisagreementNodes;
        const auto phaseAdmissible=[&](Quantum::Geometry::ActiveHermitePatch::NodalDof& fit) {
            const double value=fit.value;
            const bool valid=pureSides==1?value>=0:pureSides==2?value<=0:pureSides==3?value==0:true;
            if (!valid) {++report.purePhaseRejectedFits;trace.rejectedValues[attempt]=value;}
            if(!valid && (pureSides==1 || pureSides==2) && (fixed.empty() || !fixed[v])
                && recoverPurePhaseFit(mesh,domain,inventory,v,pureSides==1?1:-1,node.value,fit)) {
                ++report.phaseCorrectedFits;trace.scalarCorrection=fit.value-value;
                report.maxPhaseScalarCorrection=std::max(report.maxPhaseScalarCorrection,std::abs(trace.scalarCorrection));
                return true;
            }
            return valid;
        };
        // Each shared node takes ONE local branch. Opposed nearby normals
        // signal a distance ridge, where no averaged derivative is legitimate.
        bool ridge=false;
        for (const auto& piece:pieces) if (piece.normal.x*closest->normal.x+piece.normal.y*closest->normal.y<0
            && std::abs(piece.distance(p)-distance)<.125*std::min(h.x,h.y)) { ridge=true;break; }
        auto proposed=node;
        HeightReconstruction::BlendQuality blend;
        if (ridge) ++report.ambiguousDistanceNodes;
        else if (heights.blended(p,*closest,proposed,blend) && phaseAdmissible(proposed)) {
            node=proposed;
            if (blend.graphs) {
                trace.source=2;
                ++report.localGraphNodes;++report.blendedLocalGraphNodes;
                if(blend.cubicGraphs) ++report.cubicLocalGraphNodes;
                report.maxLocalGraphRms=std::max(report.maxLocalGraphRms,blend.graphRms);
                report.maxLocalGraphResidual=std::max(report.maxLocalGraphResidual,blend.graphResidual);
                report.maxLocalGraphChartSpread=std::max(report.maxLocalGraphChartSpread,blend.spread);
            } else {
                trace.source=1;
                ++report.heightReconstructedNodes;report.blendedHeightNodes+=blend.contributors>1;
                report.maxHeightChartSpread=std::max(report.maxHeightChartSpread,blend.spread);
            }
        }
        else if ((attempt=1,heights.evaluate(p,*closest,proposed)) && phaseAdmissible(proposed)) {
            node=proposed;trace.source=3;++report.heightReconstructedNodes;
        }
        else if ((attempt=2,heights.localGraph(p,*closest,proposed)) && phaseAdmissible(proposed)) {
            node=proposed;trace.source=4;++report.localGraphNodes;
        }
        else { node.dx=closest->normal.x;node.dy=closest->normal.y;trace.source=5;++report.distanceFallbackNodes; }
        // Scalar axis parity; no division by r and no artificial negative-r
        // sample. Height and local graph jets retain their mixed derivatives.
        if (domain.axisymmetric() && p.x==0) node.dx=node.dxy=0;
        report.maxNodalScalarChange=std::max(report.maxNodalScalarChange,std::abs(node.value-predicted));
        report.nodalJets.push_back(trace);
    }
    if (wallSource) for (std::size_t v=0;v<fixed.size();++v) if (fixed[v]) {
        auto& node=patch.nodes[v];const auto& old=wallSource->nodes[v];
        if(!domain.uniformLeaves() && std::abs(node.value)==band && node.dx==0 && node.dy==0 && node.dxy==0) {
            // Away from every interface, reinitialization may retain the
            // signed plateau on a stationary wall. Do not carry a varying
            // old amplitude into a hanging transition. Every active incident
            // cell must independently have this same pure geometric phase.
            const int sign=node.value>0?1:-1;bool pure=true;
            for(int label:mesh.verts(mesh.vertIndexToLabel[static_cast<int>(v)]).vertToCell)if(label>=0){
                const int c=mesh.cells(label).index();
                if(domain.cellRowByIndex()[c]>=0)pure=pure && previous.topology().cellPhase()[c]==sign;
            }
            if(pure)continue;
        }
        for (const auto& wet:boundary->wetExtensions) if (wet.vertexIndex==static_cast<int>(v)) {
            // A zero-surface gradient has an arbitrary scalar amplitude.
            // Match the retained complementary wall tangent by scaling the
            // whole proposed jet, rather than rotating its normal by fixing
            // just one component of a differently normalized graph scalar.
            const double before=wet.wallHorizontal?old.dy:old.dx;
            const double proposed=wet.wallHorizontal?node.dy:node.dx;
            const double scale=before/proposed;
            require(std::isfinite(scale) && scale>0,"wet corner cannot reconcile its cap normal with the fixed wall trace");
            node.dx*=scale;node.dy*=scale;node.dxy*=scale;
        }
        node.value=old.value;
        if (fixed[v]&1) node.dx=old.dx;if (fixed[v]&2) node.dy=old.dy;
        ++report.fixedWallTraceVertices;
    }
    if (boundary) {boundary->appendSectors(patch);boundary->constrainSectionJets(mesh,patch);}
    if(!domain.uniformLeaves()) {
        // Existing jet limiters act on independent same-level nodes. The
        // coarse/fine transition must lie in a constant signed plateau,
        // outside the reconstructed distance band. Certify the entire local
        // polynomial data, not only the sign at its centre, before reusing
        // those limiters and the native hanging constraints.
        for(int c:domain.cellIndexByRow()) {
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
            bool transition=cell.level()!=interfaceCell.level();
            for(int label:cell.cellToVert) transition=transition || mesh.verts(label).supCell>=0;
            if(!transition) continue;
            const int phase=InterfaceHeightMoments::columnPhase(inventory.liquidVolumes()[c],inventory.capacities()[c]);
            bool plateau=phase!=0;
            for(int label:cell.cellToVert) {
                const auto& node=patch.nodes[mesh.verts(label).index()];
                plateau=plateau && node.value==phase*band && node.dx==0 && node.dy==0 && node.dxy==0;
            }
            if(!plateau) throw Quantum::Geometry::InterfaceGeometryResolution(
                Quantum::Geometry::InterfaceGeometryResolution::Reason::MixedResolution,c,
                "coarse/fine transition intersects the reconstructed scalar band; enlarge the common-level band");
        }
    }
    if(domain.axisymmetric() && mesh.minCoord.x==0) for(int v=0;v<mesh.vertIndexNum();++v)
        if(mesh.verts(mesh.vertIndexToLabel[v]).coord().x==0 && std::isfinite(patch.nodes[v].value))
            report.axisJets.push_back({v,patch.nodes[v].value,patch.nodes[v].dy,patch.nodes[v].dy,patch.nodes[v].dy});
    for(auto& jet:report.nodalJets) {
        const auto& node=patch.nodes[jet.vertex];jet.value=node.value;jet.before={node.dx,node.dy,node.dxy};
    }
    limitMixedCellJets(mesh,pieces,patch,report,fixed);
    for(auto& jet:report.nodalJets) {
        const auto& node=patch.nodes[jet.vertex];jet.mixed={node.dx,node.dy,node.dxy};
    }
    for(auto& jet:report.axisJets) jet.dyMixed=patch.nodes[jet.vertex].dy;
    limitPureCellJets(mesh,domain,inventory,patch,report,fixed);
    for(auto& jet:report.nodalJets) {
        const auto& node=patch.nodes[jet.vertex];jet.pure={node.dx,node.dy,node.dxy};
    }
    for(auto& jet:report.axisJets) jet.dyPure=patch.nodes[jet.vertex].dy;
    const bool detachedTerminal=std::any_of(previous.topology().events().begin(),previous.topology().events().end(),
        [&](const auto& event) {
            return event.kind==InterfaceTopologyEvent::Kind::SectionDetachment
                && std::any_of(event.detachedSections.begin(),event.detachedSections.end(),
                    [&](const auto& section){return section.normalAxis==1 && section.coordinate==mesh.maxCoord.y;});
        });
    if(!previous.eventCurvatureCells().empty() && detachedTerminal) {
        require(time,"terminal-event geometry requires its physical predictor; zero-time remeshing is not implemented");
        // The predictor-only cap is necessary while C cannot resolve the
        // newly opened gap. It must not become a permanent alternative to
        // conservative synchronization. A full gas strip at the old section
        // proves separation using C itself; pure-cell Bernstein constraints
        // above also impose this separation on the reconstructed scalar.
        double axialWidth=0;
        for(int c:previous.eventCurvatureCells())
            axialWidth=std::max(axialWidth,mesh.cells(mesh.cellIndexToLabel[c]).size().y);
        bool separated=axialWidth>0;
        const double gasBegin=mesh.maxCoord.y-2*axialWidth;
        for(int c:domain.cellIndexByRow()) {
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
            if(mesh.verts(cell.cellToVert[2]).coord().y<=gasBegin)continue;
            ++report.terminalGasStripCells;
            separated=separated && InterfaceHeightMoments::columnPhase(
                inventory.liquidVolumes()[c],inventory.capacities()[c])==1;
        }
        report.terminalInventoryHandoff=separated && report.terminalGasStripCells>0;
        if(!report.terminalInventoryHandoff) {
            // Below a cell width C cannot distinguish the attached neck from
            // an equal-volume detached cap. Keep the transported scalar in
            // this explicitly modeled region until C certifies the gas gap.
            // Shared boundary nodes also touch adjacent cells, all included
            // in the reported model scope. C is never overwritten here.
            std::vector<bool> modeledVertex(mesh.vertIndexNum(),false);
            for(int c:previous.eventCurvatureCells())for(int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert) {
                const int v=mesh.verts(label).index();
                require(mesh.verts(label).supCell<0,"terminal-event transport requires independent local nodes");
                const auto jet=predictor->scalarJetInNativeCell(c,mesh.verts(label).coord());
                require(std::isfinite(jet.value+jet.dx+jet.dy+jet.dxy),"nonfinite terminal-event predictor");
                patch.nodes[v]={jet.value,jet.dx,jet.dy,jet.dxy};modeledVertex[v]=true;
            }
            for(int c:domain.cellIndexByRow())for(int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert)
                if(modeledVertex[mesh.verts(label).index()]){target.eventCurvatureCells.push_back(c);break;}
            // Limiter provenance belongs to pre-model reconstruction.
            report.nodalJets.clear();report.axisJets.clear();
        }
    }
    const bool modeledInterior=std::any_of(previous.topology().events().begin(),previous.topology().events().end(),
        [](const auto& event){return !event.interiorNecks.empty();});
    if(modeledInterior && !previous.eventCurvatureCells().empty()) {
        require(time,"interior-event zero-time remeshing requires its own model transfer");
        // A sub-cell opened gap cannot be inferred from C alone. Release its
        // transported model only when an entire native axial row BETWEEN the
        // current child tips is independently pure gas. This uses moving
        // endpoints, not the original fixed split coordinate.
        std::vector<std::pair<double,double>> intervals;
        for(const auto& surface:previous.topology().surfaces()) {
            const auto& chain=surface.geometry.nativeChain();
            require(!chain.closed() && chain.axisymmetric(),"interior-gap transport requires RZ axis components");
            std::array<double,2> z{};
            for(bool end:{false,true}) {
                const auto* section=chain.endpointSymmetryPlane(end);
                const auto kind=end?chain.endEndpointKind():chain.beginEndpointKind();
                require(kind==LevelSetInterfacePanels::EndpointKind::Axis || (section && section->normalAxis==1),
                    "interior-gap component has an unsupported endpoint");
                z[end]=chain.frameAt(end?static_cast<int>(chain.graphs().size())-1:0,end?1:0).point.y;
            }
            intervals.push_back(std::minmax(z[0],z[1]));
        }
        std::sort(intervals.begin(),intervals.end());
        for(std::size_t i=1;i<intervals.size();++i) {
            require(intervals[i-1].second<intervals[i].first,"interior-gap axis intervals overlap");
            bool gasRow=false;
            for(int axis:domain.cellIndexByRow()) {
                const auto box=previous.geometry()->nativeCellBounds(axis);
                if(box[0].x!=0 || box[0].y<intervals[i-1].second || box[1].y>intervals[i].first)continue;
                bool pure=true;
                for(int c:domain.cellIndexByRow()) {
                    const auto other=previous.geometry()->nativeCellBounds(c);
                    if(other[1].y<=box[0].y || other[0].y>=box[1].y)continue;
                    pure=pure && InterfaceHeightMoments::columnPhase(inventory.liquidVolumes()[c],inventory.capacities()[c])==1;
                }
                gasRow=gasRow || pure;
            }
            report.interiorGapProtected=report.interiorGapProtected || !gasRow;
        }
        if(report.interiorGapProtected) {
            std::vector<bool> modeledVertex(mesh.vertIndexNum(),false);
            for(int c:previous.eventCurvatureCells())for(int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert) {
                const auto& vertex=mesh.verts(label);const int v=vertex.index();
                require(vertex.supCell<0,"interior-event transport requires independent model nodes");
                const auto jet=predictor->scalarJetInNativeCell(c,vertex.coord());
                require(std::isfinite(jet.value+jet.dx+jet.dy+jet.dxy),"nonfinite interior-event predictor");
                patch.nodes[v]={jet.value,jet.dx,jet.dy,jet.dxy};modeledVertex[v]=true;
            }
            for(int c:domain.cellIndexByRow())for(int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert)
                if(modeledVertex[mesh.verts(label).index()]){target.eventCurvatureCells.push_back(c);break;}
            std::sort(target.eventCurvatureCells.begin(),target.eventCurvatureCells.end());
            target.eventCurvatureCells.erase(std::unique(target.eventCurvatureCells.begin(),target.eventCurvatureCells.end()),target.eventCurvatureCells.end());
            report.nodalJets.clear();report.axisJets.clear();
        }
    }
    if(modeledInterior && !previous.eventCurvatureCells().empty())
        correctEventInventory(mesh,domain,previous,inventory,patch,report);
    if(!domain.uniformLeaves()) for(int v=0;v<mesh.vertIndexNum();++v)
        if(mesh.verts(mesh.vertIndexToLabel[v]).supCell>=0) patch.nodes[v]={nan,nan,nan,nan};
    const bool interiorSplit=neckEvents && time && InterfaceNeckCap::prepareInterior(mesh,domain,previous,
        predictor,inventory,*time,patch,neckEvents->radiusCells,neckEvents->maximumLevel,destination,
        neckEvents->retainMaterial?&report.interiorMaterialTransfer:nullptr);
    report.interiorNeckSplit=interiorSplit;
    try {
        if(!interiorSplit)
        InterfaceGeometryReconstruction::prepare(mesh,std::move(patch),domain.axisymmetric(),contacts,time?&previous:nullptr,time,destination,
            previous.topology().symmetryPlanes());
    } catch(const Quantum::Geometry::InterfaceGeometryResolution& failure) {
        if(!detachedTerminal)throw;
        throw Quantum::Geometry::InterfaceGeometryResolution(failure.reason(),failure.cellIndex(),
            "post-section synchronization, inventory_handoff="+std::to_string(report.terminalInventoryHandoff)
            +", gas_strip_cells="+std::to_string(report.terminalGasStripCells)+": "+failure.what());
    }
    if(!time) {
        std::vector<std::vector<int>> identity(mesh.cellIndexNum());
        for(int c:domain.cellIndexByRow())identity[c].push_back(c);
        InterfaceMeshTransferTopology::retain(previous.topology(),*target.topology,identity);
        target.progress=previous.progress();
        const auto checkShift=[&](const InterfaceTopology& from,const InterfaceTopology& to) {
            for(const auto& surface:from.surfaces()) {
                const auto other=std::find_if(to.surfaces().begin(),to.surfaces().end(),
                    [&](const auto& candidate){return candidate.id==surface.id;});
                require(other!=to.surfaces().end(),"zero-time projection lost its surface identity");
                const auto& chain=surface.geometry.nativeChain();
                for(std::size_t g=0;g<chain.graphs().size();++g) {
                    const auto size=mesh.cells(mesh.cellIndexToLabel[chain.graphs()[g].nativeCellIndex]).size();
                    for(double fraction:{0.,.5,1.}) {
                        const auto frame=chain.frameAt(static_cast<int>(g),fraction);
                        const auto projection=other->geometry.projectedFrame(frame.point);
                        require(std::isfinite(projection.distance) && projection.distance<=.25*std::min(size.x,size.y)
                            && frame.normal.x*projection.frame.normal.x+frame.normal.y*projection.frame.normal.y>0,
                            "zero-time scalar synchronization exceeds local interface displacement/orientation bound");
                        report.maxZeroTimeInterfaceShift=std::max(report.maxZeroTimeInterfaceShift,projection.distance);
                    }
                }
            }
        };
        checkShift(previous.topology(),destination.topology());checkShift(destination.topology(),previous.topology());
    }
    require(destination.topology().events().size()==previous.topology().events().size()+(interiorSplit?1:0),
        "synchronization changed topology; conservative event/material reconstruction is not yet implemented");
    resolvedComponents(destination,mesh.cellIndexNum(),boundaries!=nullptr);
    if (boundary) boundary->validateGeometry(destination);
    InterfacePhaseGeometry::MappedFacePartitions partitions(mesh,domain,destination.geometry());
    auto measure=InterfaceCellMeasure::evaluateVolumes(mesh,domain,partitions,destination);
    report.geometricVolume=measure.summedNegativeVolume;
    report.geometricVolumeByCell=std::move(measure.negativeVolumeByCell);
    for (int c:domain.cellIndexByRow()) {
        const auto volume=inventory.liquidVolumes()[c],difference=report.geometricVolumeByCell[c]-volume;
        require(std::isfinite(difference),"nonfinite geometry/inventory mismatch");
        report.inventoryVolume+=volume;report.signedDifference+=difference;report.l1Difference+=std::abs(difference);
        report.maxFractionDifference=std::max(report.maxFractionDifference,std::abs(difference)/inventory.capacities()[c]);
    }
    std::vector<std::uint8_t> assigned(mesh.cellIndexNum(),0);
    for (const auto& liquid:destination.topology().liquids()) {
        InterfaceConservativeGeometryReport::ComponentDifference difference;difference.liquidId=liquid.id;
        for (int c:liquid.cells) {
            assigned[c]=1;
            const auto volume=inventory.liquidVolumes()[c],geometric=report.geometricVolumeByCell[c];
            difference.inventoryVolume+=volume;difference.geometricVolume+=geometric;
            difference.signedDifference+=geometric-volume;difference.l1Difference+=std::abs(geometric-volume);
        }
        report.components.push_back(difference);
    }
    for (int c:domain.cellIndexByRow()) if (!assigned[c]) report.unassignedInventoryVolume+=inventory.liquidVolumes()[c];
    // No old-style instantaneous phase budget is invented for this state.
    // New consumers obtain accepted step fluxes from conservativePhase() and
    // independent geometric measures from conservativeGeometry().
    if(!interiorSplit)target.componentInventories=transportComponents(mesh,domain,previous,inventory,
        destination.topology(),time,report);
    target.conservativePhase.emplace(std::move(inventory));
    if(curvatureDiagnostics==CurvatureDiagnostics::Record) {
    InterfaceCurvature curvature(mesh,domain,destination);
    const auto evaluate=[&](std::size_t surface,int graph,double fraction) {
        auto sample=curvature.evaluate(surface,graph,fraction);
        if (sample.height.valid()) ++report.validHeightSamples;else ++report.invalidHeightSamples;
        if (sample.method==InterfaceCurvature::Method::LocalFit) ++report.validLocalFitSamples;
        if (!sample.valid()) ++report.unavailableCurvatureSamples;
        return sample;
    };
    for (std::size_t s=0;s<destination.topology().surfaces().size();++s) {
        const auto& surface=destination.topology().surfaces()[s];
        const auto& chain=surface.geometry.nativeChain();
        InterfaceConservativeGeometryReport::SurfaceCurvature samples;samples.surfaceId=surface.id;
        for (std::size_t g=0;g<chain.graphs().size();++g) {
            samples.graphMidpoints.push_back(evaluate(s,static_cast<int>(g),.5));
        }
        if (!chain.closed()) for (bool end:{false,true}) {
            const int graph=end?static_cast<int>(chain.graphs().size())-1:0;
            const auto kind=end?chain.endEndpointKind():chain.beginEndpointKind();
            if (kind==LevelSetInterfacePanels::EndpointKind::Axis) samples.axisEnds.push_back(evaluate(s,graph,end?1:0));
            if (kind==LevelSetInterfacePanels::EndpointKind::SymmetryPlane) samples.symmetryEnds.push_back(evaluate(s,graph,end?1:0));
        }
        report.curvature.push_back(std::move(samples));
    }
    }
    target.conservativeGeometry.emplace(std::move(report));
    preparation.complete(InterfaceState::Stage::Ready);
}
InterfaceConservativeGeometry::TransportDonors InterfaceConservativeGeometry::reconstructTransportDonors(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceState& source,const BoundaryInput* boundaries,DonorModel model)
{
    source.validateMesh(mesh);domain.validateMeshIdentity(mesh);
    require(source.hasConservativePhase() && source.hasTopology(),"complete conservative donor source required");
    require(model==DonorModel::QuadraticMoments || model==DonorModel::ProjectedMoments
        || model==DonorModel::NormalFrameVolumeFit,"unknown transport donor model");
    resolvedComponents(source,mesh.cellIndexNum(),boundaries!=nullptr);
    if(boundaries) {
        const ConservativeBoundary boundary(mesh,domain,source,*boundaries);
        return transportDonors(mesh,domain,source,boundary.heightContacts,model);
    }
    return transportDonors(mesh,domain,source,{},model);
}
void InterfaceConservativeGeometry::translate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& previous,const InterfaceStepTime& time,Vector2D velocity,
    InterfaceState& destination,double distanceBandCells,bool retainSurfaceTransfer)
{
    require(destination.stage()==InterfaceState::Stage::Empty,"fresh destination required");
    previous.validateMesh(mesh);domain.validateMeshIdentity(mesh);
    require(previous.hasConservativePhase() && previous.hasTopology(),"closed conservative source required");
    require(!retainSurfaceTransfer || previous.eventCurvatureCells().empty(),
        "terminal neck model has no surface-material/charge continuation yet");
    require(std::isfinite(velocity.x) && std::isfinite(velocity.y),"nonfinite velocity");
    for (const auto face:domain.faces(mesh)) {
        if (face.negativeCell>=0 && face.positiveCell>=0) continue;
        if (domain.axisymmetric() && !face.horizontal && face.center.x==0) continue;
        const auto cell=face.negativeCell>=0?face.negativeCell:face.positiveCell;
        require(cell>=0 && previous.topology().cellPhase()[cell]==1,
            "closed translation requires a certified pure-gas exterior layer");
    }
    const Vector2D displacement=velocity*time.dt;
    const auto source=previous.geometry();
    const auto donors=transportDonors(mesh,domain,previous);
    InterfacePreparation preparation(InterfaceStateAccess::get(destination));
    const auto flux=InterfacePlicTransport::translation(mesh,domain,previous.conservativePhase(),donors.normals,
        displacement,[](const Cut::Rectangle&){return InterfacePlicTransport::CellData{0,{1,0}};},donors.curved);
    auto inventory=InterfaceConservativePhase::advance(mesh,domain,previous.conservativePhase(),time,flux);
    Quantum::Geometry::ActiveHermitePatch patch;
    patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());
    for (int v=0;v<mesh.vertIndexNum();++v) {
        const auto point=mesh.verts(mesh.vertIndexToLabel[v]).coord()-displacement;
        const Vector2D inside{std::clamp(point.x,mesh.minCoord.x,mesh.maxCoord.x),
            std::clamp(point.y,mesh.minCoord.y,mesh.maxCoord.y)};
        const auto jet=source->scalarJetNative(inside);
        auto& node=patch.nodes[v];
        if (inside.x==point.x && inside.y==point.y) node={jet.value,jet.dx,jet.dy,jet.dxy};
        else {
            // Outside departure points are in the explicitly prescribed gas
            // extension. They cannot provide liquid or interior material data.
            require(jet.value>0,"exterior scalar continuation is not gas");
            node={jet.value,0,0,0};
        }
        if (domain.axisymmetric() && mesh.verts(mesh.vertIndexToLabel[v]).coord().x==0) {
            node.dx=0;node.dxy=0;
        }
    }
    auto predictor=std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,domain.axisymmetric()));
    prepare(mesh,domain,previous,std::move(inventory),std::move(predictor),time,destination,distanceBandCells);
    recordTransportDonors(*InterfaceStateAccess::get(destination).conservativeGeometry,donors);
    if (retainSurfaceTransfer) InterfaceStateAccess::get(destination).conservativeGeometry->surfaceTransfer=
        InterfaceSurfaceTransfer::create(previous,destination,time,[displacement](Vector2D p){return p+displacement;});
    preparation.complete(InterfaceState::Stage::Ready);
}
void InterfaceConservativeGeometry::advance(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& previous,const InterfaceStepTime& time,const InterfaceMotionVelocityLayers& velocities,
    const InterfaceNormalFluxView& previousFlux,const InterfaceNormalFluxView& currentFlux,
    InterfaceState& destination,double distanceBandCells,const std::vector<std::vector<double>>& partitions,
    AdvanceReport* report,bool retainSurfaceTransfer,const BoundaryInput* boundaries,bool parallelFlux,
    CurvatureDiagnostics curvatureDiagnostics,const NeckEventPolicy* neckEvents)
{
    require(destination.stage()==InterfaceState::Stage::Empty,"fresh destination required");
    previous.validateMesh(mesh);domain.validateMeshIdentity(mesh);
    require(previous.hasConservativePhase() && previous.hasTopology(),"closed conservative source required");
    require(!neckEvents || !retainSurfaceTransfer || neckEvents->retainMaterial,
        "interior-neck material motion requires its retained local closure");
    // A terminal curvature tag does not alter material topology. The ordinary
    // transfer still checks unchanged surface parents/endpoints and ordered
    // bounded projections. Actual detachment follows in a separate zero-time
    // closure, owned and consumed by the joint physical caller.
    resolvedComponents(previous,mesh.cellIndexNum(),boundaries!=nullptr);
    auto physical=velocities;
    if (domain.cellIndexByRow().size()!=static_cast<std::size_t>(mesh.cellIndexNum()))
        require(boundaries && partitions.empty() && previousFlux.canInterpolateRoots(currentFlux),
            "selected-domain motion requires explicit boundary laws and space-time flux");
    std::shared_ptr<const ConservativeBoundary> boundary;
    InterfacePlicTransport::SinglePhasePorts ports;
    if (boundaries) {
        boundary=std::make_shared<const ConservativeBoundary>(mesh,domain,previous,*boundaries);
        require(velocities.current && velocities.previous,"fixed wall motion requires both velocity layers");
        if (boundary->hasPorts) {
            require(boundaries->inflowScalar && partitions.empty() && previousFlux.canInterpolateRoots(currentFlux),
                "single-phase ports require scalar history and the space-time flux path");
            ports=[boundary](int row){return boundary->phases.at(row);};
        }
        // Every actual wall query must obey the declared no-slip model. Do
        // not clamp a violating velocity, move the pin, or erase a wall flux.
        const auto wrap=[boundary](auto velocity) {
            return [boundary,velocity](const Vector2D& p){return boundary->velocity(p,velocity(p));};
        };
        physical.current=wrap(velocities.current);physical.previous=wrap(velocities.previous);
        for (const auto& surface:previous.topology().surfaces()) {
            const auto& chain=surface.geometry.nativeChain();if (chain.closed()) continue;
            for (bool end:{false,true}) {
                const auto p=chain.frameAt(end?static_cast<int>(chain.graphs().size())-1:0,end?1:0).point;
                physical.current(p);physical.previous(p);
            }
        }
    }
    InterfacePreparation preparation(InterfaceStateAccess::get(destination));
    AdvanceReport quality;
    auto stageBegin=std::chrono::steady_clock::time_point{};
    if(report)stageBegin=std::chrono::steady_clock::now();
    const auto timeStage=[&](int index) {
        if(!report)return;
        const auto end=std::chrono::steady_clock::now();
        quality.stageSeconds[index]=std::chrono::duration<double>(end-stageBegin).count();stageBegin=end;
    };
    const auto donors=transportDonors(mesh,domain,previous,
        boundary?boundary->heightContacts:std::vector<ContactColumn>{});
    quality.heightPlicCells=donors.heightCells;
    timeStage(0);
    const auto flux=[&]() {
        if (partitions.empty() && previousFlux.canInterpolateRoots(currentFlux)) {
            quality.usedSpaceTime=true;
            InterfacePlicTransport::SpaceTimeOptions fluxOptions;fluxOptions.parallelFaces=parallelFlux;
            return InterfacePlicTransport::spaceTime(mesh,domain,previous.conservativePhase(),donors.normals,time,physical,
                previousFlux,currentFlux,fluxOptions,&quality.spaceTime,donors.curved,ports);
        }
        const auto pieces=InterfacePlicTransport::traceSweeps(mesh,domain,time,physical,previousFlux,currentFlux,partitions,&quality.traces);
        return InterfacePlicTransport::sweptSubfaces(mesh,domain,previous.conservativePhase(),donors.normals,pieces,{},&quality.sweeps,donors.curved);
    }();
    timeStage(1);
    auto inventory=InterfaceConservativePhase::advance(mesh,domain,previous.conservativePhase(),time,flux);
    timeStage(2);
    ActiveLevelSetTransport::Options options;
    std::function<Vector2D(Vector2D)> forward;
    options.physicalVelocity=true;options.radialAxisSymmetry=domain.axisymmetric() && mesh.minCoord.x==0;
    if (boundary) boundary->transportOptions(mesh,options);
    if (boundary && boundary->hasPorts)
        options.inflowScalar=[boundary,history=boundaries->inflowScalar,end=time.end](Boundary side,Vector2D p,double elapsed) {
            const double value=history(side,p,end-elapsed);boundary->validateScalar(side,p,value);return value;
        };
    auto predictor=std::visit([&](const auto& patch) -> std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> {
        auto candidate=[&]() {
            if constexpr (std::is_same_v<std::decay_t<decltype(patch)>,Quantum::Geometry::ActiveHermitePatch>)
                return ActiveLevelSetTransport::advanceNativeJet(mesh,patch,time.dt,physical.current,physical.previous,options);
            else return ActiveLevelSetTransport::advanceQuarticNativeJet(mesh,patch,time.dt,physical.current,physical.previous,options);
        }();
        if constexpr (std::is_same_v<std::decay_t<decltype(patch)>,Quantum::Geometry::ActiveHermitePatch>) {
            if(boundary && !boundary->sections.empty()) {
                boundary->constrainSectionJets(mesh,candidate.patch);
                candidate.geometry=std::make_unique<LevelSetHermiteGeometry::Snapshot>(
                    LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,candidate.patch,options.radialAxisSymmetry));
            }
        }
        quality.predictedNodes=candidate.activeNodes;quality.maxTraceRatio=candidate.maxTraceRatio;
        if (retainSurfaceTransfer) forward=std::move(candidate.forwardMaterialPoint);
        return std::move(candidate.geometry);
    },previous.scalarPatch());
    timeStage(3);
    prepare(mesh,domain,previous,std::move(inventory),std::move(predictor),time,destination,distanceBandCells,boundaries,curvatureDiagnostics,neckEvents);
    recordTransportDonors(*InterfaceStateAccess::get(destination).conservativeGeometry,donors);
    timeStage(4);
    if (retainSurfaceTransfer) InterfaceStateAccess::get(destination).conservativeGeometry->surfaceTransfer=
        InterfaceSurfaceTransfer::create(previous,destination.conservativeGeometry().interiorMaterialTransfer
            ? *destination.conservativeGeometry().interiorMaterialTransfer->sourceState:destination,time,forward);
    timeStage(5);
    if (report) *report=std::move(quality);
    preparation.complete(InterfaceState::Stage::Ready);
}
