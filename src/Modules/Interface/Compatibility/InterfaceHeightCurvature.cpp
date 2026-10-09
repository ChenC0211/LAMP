#include "Modules/Interface/Compatibility/InterfaceHeightCurvature.h"
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianLeafLevelIndex.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceHeightMoments.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceSymmetryStencil.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <set>

namespace
{
void require(bool ok,const char* why) { if (!ok) throw std::invalid_argument(why); }
double coordinate(Vector2D p,int axis) { return axis==0?p.x:p.y; }

}
struct InterfaceHeightCurvature::Impl
{
    const QuadtreeMesh& mesh;
    const NativeLeafDomain& domain;
    const InterfaceState& state;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source;
    const InterfaceConservativePhase* inventory;
    CartesianLeafLevelIndex grid;
    Impl(const QuadtreeMesh& m,const NativeLeafDomain& d,const InterfaceState& s)
        :mesh(m),domain(d),state(s),source(s.geometry()),inventory(&s.conservativePhase()),
         grid(m)
    {
        s.validateMesh(m);inventory->validate(m,d);
        require(s.hasTopology(),"height curvature requires conservative topology");
    }
    Sample attempt(std::size_t surface,int graph,const LevelSetHermiteGeometry::Frame& query,int axis,int radius) const
    {
        Sample result;result.point=query.point;result.normalAxis=axis;result.halfHeight=radius;
        const auto& chain=state.topology().surfaces()[surface].geometry.nativeChain();
        const auto address=grid.address(chain.graphs()[graph].nativeCellIndex);
        const std::array<int,2> centre{{address.i,address.j}};
        const auto extent=grid.extent(address.level);
        const int tangent=1-axis,nt=extent[tangent],nn=extent[axis];
        const auto h=grid.spacing(address.level);const double hn=coordinate(h,axis),ht=coordinate(h,tangent);
        const auto cell=[&](int t,int n,int normalAxis) {
            return normalAxis==1?grid.cell(address.level,t,n):grid.cell(address.level,n,t);
        };
        const auto point=query.point;
        const double normalComponent=coordinate(query.normal,axis);
        if (normalComponent==0) { result.status=Status::PoorOrientation;return result; }
        const int sign=normalComponent>0?1:-1;
        int first=centre[tangent]-2,count=5;
        const auto reflection=InterfaceSymmetryStencil::select(mesh.minCoord,mesh.maxCoord,tangent,nt,
            first,count,domain.axisymmetric(),chain.symmetryPlanes());
        const bool even=reflection.even;
        int boundaryGraph=-1;Vector2D boundaryAnchor;
        if (!chain.closed()) for (bool end:{false,true}) {
                const auto kind=end?chain.endEndpointKind():chain.beginEndpointKind();
                if(kind==LevelSetInterfacePanels::EndpointKind::SymmetryPlane && even && !reflection.axis) {
                    const auto* plane=chain.endpointSymmetryPlane(end);
                    if(plane && plane->normalAxis==tangent && plane->coordinate==reflection.origin) {
                        if(boundaryGraph>=0) {result.status=Status::MixedBranches;return result;}
                        boundaryGraph=end?static_cast<int>(chain.graphs().size())-1:0;
                        boundaryAnchor=chain.frameAt(boundaryGraph,end?1:0).point;
                    }
                    continue;
                }
                if (kind!=LevelSetInterfacePanels::EndpointKind::DeclaredBoundary) continue;
                const int g=end?static_cast<int>(chain.graphs().size())-1:0;
                const auto endpoint=chain.frameAt(g,end?1:0).point;
                const int owner=chain.graphs()[g].nativeCellIndex;
                const auto& cell=mesh.cells(mesh.cellIndexToLabel[owner]);
                const auto lower=mesh.verts(cell.cellToVert[0]).coord(),upper=mesh.verts(cell.cellToVert[2]).coord();
                const auto endpointAddress=grid.address(owner);
                // A remote endpoint at another level cannot anchor this
                // local chart. Any cross-level cells in its actual rectangle
                // are rejected below before looking at phase values.
                if(endpointAddress.level!=address.level) continue;
                const int index=tangent==0?endpointAddress.i:endpointAddress.j;int shifted=first;
                if (coordinate(endpoint,tangent)==coordinate(lower,tangent) && first<index) shifted=index;
                if (coordinate(endpoint,tangent)==coordinate(upper,tangent) && first+count>index+1) shifted=index+1-count;
                if (shifted==first) continue;
                if (even) return result;
                if (boundaryGraph>=0) { result.status=Status::MixedBranches;return result; }
                boundaryGraph=g;boundaryAnchor=endpoint;first=shifted;
            result.oneSidedBoundary=true;result.boundaryAnchor=boundaryAnchor;
        }
        if (first<0 || first+count>nt) return result;
        if(even && !reflection.axis && boundaryGraph<0) {result.status=Status::MixedBranches;return result;}
        result.columns=count;result.evenAxisFit=reflection.axis;result.evenSectionFit=even && !reflection.axis;
        const int lo=std::max(0,centre[axis]-radius),hi=std::min(nn-1,centre[axis]+radius);
        std::set<int> supportCells;
        for(int t=first;t<first+count;++t) for(int n=lo;n<=hi;++n) {
            const int c=cell(t,n,axis);
            if(c<0) {result.status=Status::MixedResolution;return result;}
            supportCells.insert(c);
        }
        // All graphs intersecting the rectangle must belong to ONE local arc,
        // not merely to one surface/liquid. No nearby opposite branch enters
        // either the volume sum or derivative stencil unnoticed.
        std::vector<bool> selected(chain.graphs().size(),false);
        for (std::size_t s=0;s<state.topology().surfaces().size();++s) {
            const auto& surfaceGeometry=state.topology().surfaces()[s].geometry;
            const auto& other=surfaceGeometry.nativeChain();
            for (std::size_t g=0;g<other.graphs().size();++g) {
                if(!supportCells.count(other.graphs()[g].nativeCellIndex)) continue;
                if (s!=surface) { result.status=Status::MixedBranches;return result; }
                selected[g]=true;result.branchGraphs.push_back(static_cast<int>(g));
                const auto& panel=surfaceGeometry.panels()[g];
                for (auto normal:{panel.beginFrame.normal,panel.frame.normal,panel.endFrame.normal})
                    // A merely positive component admits nearly vertical
                    // heights and unbounded derivatives at remote columns.
                    // Bound the sampled chart slope by sqrt(3) over the
                    // whole supporting arc, not only at the query point.
                    // The regular even axis stencil retains its original
                    // sign guard: it fits parity moments about r=0 and is
                    // also the sole source of the azimuthal axis limit.
                    if (sign*coordinate(normal,axis)<=0 || (!reflection.axis && sign*coordinate(normal,axis)<.5)) {
                    result.status=Status::PoorOrientation;return result;
                }
            }
        }
        int blocks=0;
        for (std::size_t g=0;g<selected.size();++g) if (selected[g]
            && (g==0?(!chain.closed() || !selected.back()):!selected[g-1])) ++blocks;
        if (!selected[graph] || blocks!=1 || (boundaryGraph>=0 && !selected[boundaryGraph])) {
            result.status=Status::MixedBranches;return result;
        }
        std::vector<InterfaceHeightMoments::Column> columns;
        const long double origin=even?reflection.origin:coordinate(boundaryGraph>=0?boundaryAnchor:point,tangent);
        const long double queryNormal=coordinate(point,axis);
        const long double reference=boundaryGraph>=0?coordinate(boundaryAnchor,axis):queryNormal;
        for (int column=0;column<count;++column) {
            const int ci=first+column;
            int columnLo,columnHi;
            if (!InterfaceHeightMoments::activeSpan(centre[axis],lo,hi,[&](int row){
                return domain.cellRowByIndex()[cell(ci,row,axis)]>=0;
            },columnLo,columnHi)) return result;
            const auto phaseAt=[&](int row) {
                const int c=cell(ci,row,axis);
                return InterfaceHeightMoments::columnPhase(inventory->liquidVolumes()[c],inventory->capacities()[c]);
            };
            const auto low=phaseAt(columnLo),high=phaseAt(columnHi);
            const int lowCell=cell(ci,columnLo,axis),highCell=cell(ci,columnHi,axis);
            if (low!=-sign || high!=sign
                || state.topology().cellPhase()[lowCell]!=(sign>0?-1:1)
                || state.topology().cellPhase()[highCell]!=(sign>0?1:-1)) {
                result.status=Status::MissingPureEnds;return result;
            }
            int phase=0;long double volume=0;
            for (int row=columnLo;row<=columnHi;++row) {
                const int next=1+sign*phaseAt(row);
                if (next<phase) { result.status=Status::MultipleTransitions;return result; }
                const int c=cell(ci,row,axis);const auto liquid=inventory->liquidVolumes()[c],capacity=inventory->capacities()[c];
                if (liquid>0 && liquid<capacity && next!=1) {
                    ++result.unresolvedColumnCells;result.unresolvedColumnVolume+=std::min(liquid,capacity-liquid);
                }
                phase=next;volume+=liquid;
            }
            const auto& a=mesh.cells(mesh.cellIndexToLabel[lowCell]);
            const auto& b=mesh.cells(mesh.cellIndexToLabel[highCell]);
            const auto lower=mesh.verts(a.cellToVert[0]).coord(),upper=mesh.verts(b.cellToVert[2]).coord();
            const long double tlo=coordinate(lower,tangent),thi=coordinate(upper,tangent);
            const long double nlo=coordinate(lower,axis),nhi=coordinate(upper,axis);
            columns.push_back({tlo,thi,nlo,nhi,volume});
            if (column==0) {result.supportLower=lower;result.supportUpper=upper;}
            else {
                result.supportLower={std::min(result.supportLower.x,lower.x),std::min(result.supportLower.y,lower.y)};
                result.supportUpper={std::max(result.supportUpper.x,upper.x),std::max(result.supportUpper.y,upper.y)};
            }
        }
        auto fit=InterfaceHeightMoments::fit(columns,origin,reference,ht,hn,axis,sign,domain.axisymmetric(),even,
            0,boundaryGraph>=0 && !even);
        if(reflection.axis) fit=InterfaceHeightMoments::axisCurvatureFit(columns,fit,sign);
        result.implicitAxisFit=fit.implicitEven;
        result.scaledRankRatio=fit.rankRatio;
        if (!fit.valid) {result.status=Status::SingularFit;return result;}
        // Convert squared-radius residual to a height scale as well. This
        // bounds a geometric column mismatch, not an inventory correction.
        result.scaledMomentResidual=fit.scaledResidual*static_cast<double>(fit.scale)/hn;
        if (fit.squaredRadius) result.scaledMomentResidual/=2*std::max(static_cast<double>(queryNormal),hn);
        if (result.scaledMomentResidual>.025) {result.status=Status::PoorMomentFit;return result;}
        const auto height=fit.evaluate(coordinate(point,tangent));
        if (!height.valid) {result.status=Status::GeometryMismatch;return result;}
        const auto f=height.height,df=height.first,ddf=height.second;
        result.height=static_cast<double>(f);result.first=static_cast<double>(df);result.second=static_cast<double>(ddf);
        result.heightMismatch=static_cast<double>(std::abs(f-queryNormal));
        if (result.heightMismatch>hn*.5 || !std::isfinite(df) || !std::isfinite(ddf)) {
            result.status=Status::GeometryMismatch;return result;
        }
        const auto length=std::sqrt(1+df*df);
        const double ntangent=static_cast<double>(-sign*df/length),nnormal=static_cast<double>(sign/length);
        result.normal=axis==1?Vector2D{ntangent,nnormal}:Vector2D{nnormal,ntangent};
        result.meridional=static_cast<double>(-sign*ddf/(length*length*length));
        result.azimuthal=0;
        if (domain.axisymmetric()) {
            if (reflection.axis) {
                const auto x=(coordinate(point,tangent)-origin)/ht;
                const long double radialSlopeRatio=fit.implicitEven
                    ?fit.scale/(ht*ht)*2*fit.coefficients[1]/std::sqrt(1-4*fit.coefficients[2]*(fit.coefficients[0]+fit.coefficients[1]*x*x))
                    :fit.scale/(ht*ht)*(2*fit.coefficients[1]+4*fit.coefficients[2]*x*x);
                result.azimuthal=static_cast<double>(-sign*radialSlopeRatio/length);
                // At the exact axis the even graph gives H'(0)=0 and
                // lim H'(r)/r=H''(0). Evaluate that shared principal limit
                // once: algebraically equivalent products/divisions above
                // can otherwise differ by one ulp and spuriously reject an
                // entire axis-to-axis regularization branch.
                if(point.x==0)result.azimuthal=result.meridional;
            } else {
                if (point.x<=0) { result.status=Status::GeometryMismatch;return result; }
                result.azimuthal=result.normal.x/point.x;
            }
        }
        result.curvature=result.meridional+result.azimuthal;
        if (!std::isfinite(result.curvature)) { result.status=Status::SingularFit;return result; }
        result.status=Status::Valid;return result;
    }
};
InterfaceHeightCurvature::InterfaceHeightCurvature(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& state):impl_(std::make_unique<Impl>(mesh,domain,state)) {}
InterfaceHeightCurvature::~InterfaceHeightCurvature()=default;
InterfaceHeightCurvature::Sample InterfaceHeightCurvature::evaluate(std::size_t surface,int graph,double fraction) const
{
    const auto& p=*impl_;
    require(p.state.geometry()==p.source && p.state.hasConservativePhase() && &p.state.conservativePhase()==p.inventory,
        "height curvature source state was replaced");
    const auto& chain=p.state.topology().surfaces().at(surface).geometry.nativeChain();
    const auto frame=chain.frameAt(graph,fraction);
    const int preferred=std::abs(frame.normal.y)>=std::abs(frame.normal.x)?1:0;
    Sample failure;
    // Five exact column moments give the centred flat-interface H'' stencil
    // [-1/8,3/2,-11/4,3/2,-1/8]/h^2, with Fourier symbol
    // -4*sin(theta/2)^2*(1+sin(theta/2)^2/2)/h^2 <= 0.
    // Seven-column quartic least squares instead has +34/33 at Nyquist:
    // fitting smooth data well does not preserve capillary restoring signs.
    // Keep the compact moment operator and the dedicated even-axis limit.
    for (int axis:{preferred,1-preferred}) for (int radius:{2,3,4,5}) {
        auto sample=p.attempt(surface,graph,frame,axis,radius);
        if (sample.valid()) return sample;
        failure=std::move(sample);
    }
    return failure;
}
