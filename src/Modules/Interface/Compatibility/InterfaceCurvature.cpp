#include "Modules/Interface/Compatibility/InterfaceCurvature.h"
#include <iomanip>
#include <sstream>
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Geometry/ArcScalarRegularization.h"
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
double dot(Vector2D a,Vector2D b) { return a.x*b.x+a.y*b.y; }
double length(Vector2D a) { return std::hypot(a.x,a.y); }
}
struct InterfaceCurvature::Impl
{
    const InterfaceState& state;
    const QuadtreeMesh& mesh;
    const NativeLeafDomain& domain;
    InterfaceHeightCurvature heights;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source;
    const InterfaceConservativePhase* inventory;
    bool rz;
    // Eager immutable knots: force faces repeatedly query the same accepted
    // geometry/inventory. No lazy mutation, cross-state cache or time history.
    std::vector<std::vector<Source>> midpoints;
    std::vector<std::array<Source,2>> endpoints;
    // Ordinary-side endpoint knots at an explicit event-model boundary.
    // Prepared eagerly so parallel force queries never refit or mutate them.
    std::vector<std::vector<std::array<Source,2>>> modelBoundaries;
    std::vector<std::vector<std::array<double,2>>> halfLengths;
    std::vector<RegularizationStatus> regularization;
    std::vector<std::vector<std::uint8_t>> radialCharts;
    bool eventCell(int cell) const {
        const auto& cells=state.eventCurvatureCells();
        return std::find(cells.begin(),cells.end(),cell)!=cells.end();
    }
    Impl(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceState& s)
        :state(s),mesh(mesh),domain(domain),heights(mesh,domain,s),source(s.geometry()),inventory(&s.conservativePhase()),rz(domain.axisymmetric()) {}
    void validate() const
    {
        if (state.geometry()!=source || !state.hasConservativePhase() || &state.conservativePhase()!=inventory)
            throw std::invalid_argument("curvature source state was replaced");
    }
};
InterfaceCurvature::InterfaceCurvature(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceState& state,
    Regularization regularization)
    :impl_(std::make_unique<Impl>(mesh,domain,state))
{
    auto& p=*impl_;const auto& surfaces=state.topology().surfaces();
    if(regularization!=Regularization::None && regularization!=Regularization::ArcThirdDerivative)
        throw std::invalid_argument("unknown curvature regularization");
    p.midpoints.resize(surfaces.size());p.endpoints.resize(surfaces.size());p.halfLengths.resize(surfaces.size());
    p.modelBoundaries.resize(surfaces.size());
    p.regularization.assign(surfaces.size(),RegularizationStatus::NotRequested);
    p.radialCharts.resize(surfaces.size());
    std::vector<int> radialCells(mesh.cellIndexNum(),-1);
    for (std::size_t s=0;s<surfaces.size();++s) {
        const auto& chain=surfaces[s].geometry.nativeChain();const auto count=chain.graphs().size();
        p.radialCharts[s].resize(count);
        for(std::size_t g=0;g<count;++g) {
            const auto& graph=chain.graphs()[g];bool radial=!graph.coordinateIsX;
            if(!radial && p.rz && !state.eventCurvatureCells().empty()) {
                auto& certified=radialCells[graph.nativeCellIndex];
                if(certified<0)certified=p.source->nativeZeroSetPartialResolved(graph.nativeCellIndex,0)?1:0;
                radial=certified==1;
                // A cap cell may also contain the axis pole, where phi_r=0.
                // Certify the actual retained subchart before declaring a
                // remote tail's connected radial arc unresolved. This uses
                // outward scalar bounds, never sampled endpoint normals.
                if(!radial && (graph.localLower.x!=0 || graph.localLower.y!=0
                    || graph.localUpper.x!=1 || graph.localUpper.y!=1))
                    radial=p.source->nativeZeroSetPartialResolved(graph.nativeCellIndex,0,
                        graph.localLower,graph.localUpper);
            }
            p.radialCharts[s][g]=radial;
        }
        p.midpoints[s].reserve(count);p.halfLengths[s].resize(count);p.modelBoundaries[s].resize(count);
        for (std::size_t g=0;g<count;++g) {
            p.midpoints[s].push_back(evaluatePointwise(s,static_cast<int>(g),.5));
            if(!p.eventCell(chain.graphs()[g].nativeCellIndex))for(int side:{0,1}) {
                const int neighbor=static_cast<int>(g)+(side?1:-1);
                if((chain.closed() || (neighbor>=0 && neighbor<static_cast<int>(count)))
                    && p.eventCell(chain.graphs()[(neighbor+static_cast<int>(count))%count].nativeCellIndex))
                    p.modelBoundaries[s][g][side]=evaluatePointwise(s,static_cast<int>(g),side);
            }
            for (int half:{0,1}) for (const auto& q:chain.quadrature(static_cast<int>(g),half*.5,(half+1)*.5))
                p.halfLengths[s][g][half]+=q.weight;
        }
        if (!chain.closed()) {
            p.endpoints[s][0]=evaluatePointwise(s,0,0);
            p.endpoints[s][1]=evaluatePointwise(s,static_cast<int>(count)-1,1);
        }
        if(regularization==Regularization::None)continue;
        auto& status=p.regularization[s];status=RegularizationStatus::UnsupportedGeometry;
        const bool event=std::any_of(chain.graphs().begin(),chain.graphs().end(),
            [&](const auto& graph){return p.eventCell(graph.nativeCellIndex);});
        if(event) {
            // Split at event cells: the terminal cap must not smooth its
            // model curvature into the resolved part of the surface. Each
            // clean run retains its end knots and uses only its own ds.
            status=RegularizationStatus::EventPartition;
            for(std::size_t first=0;first<count;) {
                if(p.midpoints[s][first].method==Method::EventModel){++first;continue;}
                std::size_t last=first;
                while(last<count && p.midpoints[s][last].method!=Method::EventModel)++last;
                bool valid=last-first>=4;double h=0;
                std::vector<double> edges,meridional,azimuthal;
                for(std::size_t g=first;g<last;++g) {
                    const auto& raw=p.midpoints[s][g];
                    const auto spacing=mesh.cells(mesh.cellIndexToLabel[chain.graphs()[g].nativeCellIndex]).size();
                    const double scale=std::max(spacing.x,spacing.y);
                    if(g==first)h=scale;
                    valid&=raw.valid() && h==scale;
                    if(g>first){valid&=dot(p.midpoints[s][g-1].normal,raw.normal)>0;
                        edges.push_back(p.halfLengths[s][g-1][1]+p.halfLengths[s][g][0]);}
                    meridional.push_back(raw.meridional);azimuthal.push_back(raw.azimuthal);
                }
                if(valid) {
                    namespace Arc=Quantum::Geometry::ArcScalarRegularization;
                    const auto smooth=Arc::openThirdDerivative(edges,meridional,h,Arc::Endpoint::FixedValue,Arc::Endpoint::FixedValue);
                    const auto ring=p.rz?Arc::openThirdDerivative(edges,azimuthal,h,Arc::Endpoint::FixedValue,Arc::Endpoint::FixedValue):Arc::Result{};
                    for(std::size_t g=first;g<last;++g){auto& raw=p.midpoints[s][g];
                        raw.meridional=smooth.values[g-first];if(p.rz)raw.azimuthal=ring.values[g-first];
                        raw.curvature=raw.meridional+raw.azimuthal;raw.regularized=true;}
                }
                first=last;
            }
            continue;
        }
        using Kind=LevelSetInterfacePanels::EndpointKind;
        const bool open=!chain.closed();
        const auto admittedEnd=[&](Kind kind){return kind==Kind::DeclaredBoundary || kind==Kind::SymmetryPlane || (p.rz && kind==Kind::Axis);};
        if((open && (!admittedEnd(chain.beginEndpointKind()) || !admittedEnd(chain.endEndpointKind()))) || count<8)continue;
        std::vector<double> edges(chain.closed()?count:count+1),values,azimuthal;
        double h=0;status=RegularizationStatus::Applied;
        for(std::size_t g=0;g<count;++g) {
            const auto& raw=p.midpoints[s][g];const auto& next=p.midpoints[s][g+1<count?g+1:(chain.closed()?0:g)];
            if(!raw.valid() || !next.valid()){status=RegularizationStatus::IncompleteCurvature;break;}
            if(dot(raw.normal,next.normal)<=0){status=RegularizationStatus::UnresolvedTurn;break;}
            const auto spacing=mesh.cells(mesh.cellIndexToLabel[chain.graphs()[g].nativeCellIndex]).size();
            const double scale=std::max(spacing.x,spacing.y);if(g==0)h=scale;
            if(scale!=h){status=RegularizationStatus::MixedScale;break;}
            // Actual native half-arc quadrature, not midpoint chord length or
            // half of a whole graph's measure. No missing knot is skipped.
            if(chain.closed() || g+1<count)edges[g+(open?1:0)]=p.halfLengths[s][g][1]+p.halfLengths[s][(g+1)%count][0];
        }
        if(status!=RegularizationStatus::Applied)continue;
        if(open) {
            for(int end:{0,1}) {const auto& point=p.endpoints[s][end];const auto& neighbor=p.midpoints[s][end?count-1:0];
                if(!point.valid()){status=RegularizationStatus::IncompleteCurvature;break;}
                const bool axis=(end?chain.endEndpointKind():chain.beginEndpointKind())==Kind::Axis;
                if((axis && (point.point.x!=0 || point.normal.x!=0 || point.meridional!=point.azimuthal))
                    || dot(point.normal,neighbor.normal)<=0){status=RegularizationStatus::UnresolvedTurn;break;}}
            if(status!=RegularizationStatus::Applied)continue;
            edges.front()=p.halfLengths[s].front()[0];edges.back()=p.halfLengths[s].back()[1];
            values.push_back(p.endpoints[s][0].meridional);azimuthal.push_back(p.endpoints[s][0].azimuthal);
        }
        for(const auto& raw:p.midpoints[s]){values.push_back(raw.meridional);azimuthal.push_back(raw.azimuthal);}
        if(open){values.push_back(p.endpoints[s][1].meridional);azimuthal.push_back(p.endpoints[s][1].azimuthal);}
        namespace Arc=Quantum::Geometry::ArcScalarRegularization;
        const auto apply=[&](const std::vector<double>& field){return open
            ?Arc::openThirdDerivative(edges,field,h,
                chain.beginEndpointKind()!=Kind::DeclaredBoundary?Arc::Endpoint::FixedEven:Arc::Endpoint::FixedValue,
                chain.endEndpointKind()!=Kind::DeclaredBoundary?Arc::Endpoint::FixedEven:Arc::Endpoint::FixedValue)
            :Arc::closedThirdDerivative(edges,field,h);};
        const auto smoothed=apply(values);const auto ring=p.rz?apply(azimuthal):Arc::Result{};
        for(std::size_t g=0;g<count;++g){auto& raw=p.midpoints[s][g];const auto i=g+(open?1:0);
            raw.meridional=smoothed.values[i];if(p.rz)raw.azimuthal=ring.values[i];
            raw.curvature=raw.meridional+raw.azimuthal;raw.regularized=true;}
        if(open)for(auto& end:p.endpoints[s])end.regularized=true; // Original endpoint values are hard constraints.
    }
}
InterfaceCurvature::~InterfaceCurvature()=default;
InterfaceCurvature::RegularizationStatus InterfaceCurvature::regularizationStatus(std::size_t surface) const
{impl_->validate();return impl_->regularization.at(surface);}

InterfaceCurvature::Fit InterfaceCurvature::localFit(std::size_t surface,int graph,double fraction,LocalModel model,
    double maximumArcCells) const
{
    const auto& p=*impl_;p.validate();
    if(maximumArcCells==0) {
        auto fit=localFit(surface,graph,fraction,model,3.5);
        if(model==defaultLocalModel && (fit.status==FitStatus::GeometryMismatch
            || (fit.status==FitStatus::PoorOrientation && fit.orientationCause==OrientationCause::SupportTurn))) {
            // A finite wide stencil may reach a turn that the re-integrated
            // GL nodes resolve more clearly than the initial support search.
            // Retry on the same branch with the existing two-cell stencil;
            // its support, normal cone, rank and accuracy gates all remain.
            // Never retry a fitted fold/normal mismatch, proximity, resolution,
            // endpoint or rank failure by selecting a more convenient fit.
            const auto shorter=localFit(surface,graph,fraction,model,2.);
            if(shorter.valid())fit=shorter;
            else {
                std::ostringstream detail;detail<<std::setprecision(17)
                    <<"graph="<<graph<<", fraction="<<fraction
                    <<", status="<<static_cast<int>(shorter.status)
                    <<", orientation="<<static_cast<int>(shorter.orientationCause)
                    <<", arc_cells="<<shorter.maximumArcCells
                    <<", support="<<shorter.negativeSupport<<'/'<<shorter.positiveSupport
                    <<", residual="<<shorter.scaledResidual<<", offset="<<shorter.scaledPositionError
                    <<", rank="<<shorter.scaledRankRatio<<", samples="<<shorter.samples
                    <<", kappa_h="<<shorter.scaledMeridionalMagnitude
                    <<", refinement="<<static_cast<int>(shorter.refinementCause);
                fit.rejectedShorterFit=detail.str();
            }
        }
        return fit;
    }
    if(!std::isfinite(maximumArcCells) || maximumArcCells<2 || maximumArcCells>3.5)
        throw std::invalid_argument("local curvature arc bound must be in [2,3.5] cells");
    if(model!=LocalModel::IsotropicQuadratic && model!=LocalModel::TangentQuartic && model!=LocalModel::ArcMeasureQuartic
        && model!=LocalModel::SymmetricArcQuartic && model!=LocalModel::ClippedArcQuartic)
        throw std::invalid_argument("unknown local curvature model");
    const auto& geometry=p.state.topology().surfaces().at(surface).geometry;
    const auto& chain=geometry.nativeChain();const auto query=chain.frameAt(graph,fraction);
    Fit result;result.point=query.point;result.maximumArcCells=maximumArcCells;
    const auto& queryCell=p.mesh.cells(p.mesh.cellIndexToLabel[chain.graphs()[graph].nativeCellIndex]);
    const auto spacing=queryCell.size();const double h=std::max(spacing.x,spacing.y);
    result.resolutionScale=h;
    // RZ axis limits come from the even height stencil. Do not divide a noisy
    // fitted radial normal by a vanishing radius or invent endpoint data.
    if (p.rz && query.point.x<h) { result.status=FitStatus::UnsupportedEndpoint;return result; }
    const Vector2D normal=query.normal,tangent{normal.y,-normal.x};
    // A sampled resolution guard, not a certified minimum branch distance.
    // Opposed arcs within two cells need refinement even when a geometric fit
    // on just one of them would look smooth. No whole-liquid identity shortcut.
    for (const auto& other:p.state.topology().surfaces()) for (const auto& panel:other.geometry.panels())
        for (const auto& frame:{panel.beginFrame,panel.frame,panel.endFrame})
            if (dot(frame.normal,normal)<=0 && length(frame.point-query.point)<2*h) {
                result.refinementCause=RefinementCause::OpposedArc;
                result.opposedSampleDistance=length(frame.point-query.point);
                result.status=FitStatus::RefinementRequired;return result;
            }
    struct Datum { double x,y,weight; };
    std::vector<Datum> data;
    const int count=static_cast<int>(chain.graphs().size());
    std::vector<bool> selected(count,false);
    // Follow BOTH actual ordered arcs, using positive ds quadrature weights.
    // The finite stencil ends at the requested sampled arc length, a normal turn of
    // 60 degrees, an endpoint or a tangent fold; no spatial-neighbor search.
    for (int direction:{-1,1}) {
        int current=graph;double first=fraction,arc=0,previousX=0;int tangentSign=0;
        Vector2D previous=query.point;bool finished=false;
        for (int visited=0;visited<count && !finished;++visited) {
            const double a=direction>0?first:0,b=direction>0?1:first;
            auto quadrature=chain.quadrature(current,a,b);
            if (direction<0) std::reverse(quadrature.begin(),quadrature.end());
            for (const auto& sample:quadrature) {
                const auto delta=sample.frame.point-query.point;
                arc+=length(sample.frame.point-previous);previous=sample.frame.point;
                if (arc>maximumArcCells*h || dot(sample.frame.normal,normal)<.5) { finished=true;break; }
                if(p.mesh.cells(p.mesh.cellIndexToLabel[chain.graphs()[current].nativeCellIndex]).level()!=queryCell.level()) {
                    result.status=FitStatus::RefinementRequired;result.refinementCause=RefinementCause::MixedResolution;return result;
                }
                const double x=dot(delta,tangent)/h;
                if (!tangentSign && x!=0) tangentSign=x>0?1:-1;
                if (!tangentSign || tangentSign*x<=tangentSign*previousX) { finished=true;break; }
                previousX=x;
                data.push_back({x,dot(delta,normal)/h,sample.weight/h});selected[current]=true;
                result.negativeSupport=std::max(result.negativeSupport,-x*h);
                result.positiveSupport=std::max(result.positiveSupport,x*h);
            }
            current+=direction;
            if (current<0 || current>=count) {
                if (!chain.closed()) break;
                current=(current+count)%count;
            }
            first=direction>0?0:1;
        }
    }
    result.samples=data.size();
    for (int g=0;g<count;++g) if (selected[g]) result.branchGraphs.push_back(g);
    if (data.size()<10 || result.negativeSupport<h || result.positiveSupport<h) {
        result.status=FitStatus::InsufficientArc;return result;
    }
    // Reintegrate the actual native graphs at the admitted support endpoints;
    // do not retain whole quadrature weights after cutting through a graph.
    // Production retains both available sides. Symmetric clipping reduces
    // truncation bias on smooth data, but discards useful support on transported
    // geometry. Keep it as an explicit control, not a forced production rule.
    // Branch/turn/resolution gates remain identical for all clipped variants.
    if(model!=LocalModel::ArcMeasureQuartic) {
    const double extent=std::min(result.negativeSupport,result.positiveSupport)/h;
    const double lowerBound=model==LocalModel::ClippedArcQuartic?-result.negativeSupport/h:-extent;
    const double upperBound=model==LocalModel::ClippedArcQuartic?result.positiveSupport/h:extent;
    const bool arcMeasure=model==LocalModel::SymmetricArcQuartic || model==LocalModel::ClippedArcQuartic;
    data.clear();result.branchGraphs.clear();
    for(int g=0;g<count;++g) if(selected[g]) {
        const auto coordinate=[&](double f){return dot(chain.frameAt(g,f).point-query.point,tangent)/h;};
        const double x0=coordinate(0),x1=coordinate(1);
        const double lower=std::max(lowerBound,std::min(x0,x1)),upper=std::min(upperBound,std::max(x0,x1));
        if(!(lower<upper)) continue;
        const auto inverse=[&](double x) {
            if(x==x0) return 0.;if(x==x1) return 1.;
            double lo=0,hi=1;
            for(int k=0;k<52;++k) {const double mid=lo+(hi-lo)*.5;
                if(mid==lo || mid==hi) break;const double value=coordinate(mid);
                if(value==x) return mid;
                if((value<x)==(x1>x0)) lo=mid;else hi=mid;
            }
            return lo+(hi-lo)*.5;
        };
        const double first=inverse(x1>x0?lower:upper),last=inverse(x1>x0?upper:lower);
        if(!(first<last)) {result.status=FitStatus::InsufficientArc;return result;}
        for(const auto& sample:chain.quadrature(g,first,last)) {
            const double cosine=dot(sample.frame.normal,normal);
            if(cosine<.5) {result.status=FitStatus::PoorOrientation;
                result.orientationCause=OrientationCause::SupportTurn;return result;}
            const auto delta=sample.frame.point-query.point;
            data.push_back({dot(delta,tangent)/h,dot(delta,normal)/h,
                sample.weight*(arcMeasure?1:cosine)/h});
        }
        result.branchGraphs.push_back(g);
    }
    result.samples=data.size();result.negativeSupport=-lowerBound*h;result.positiveSupport=upperBound*h;
    if(data.size()<10) {result.status=FitStatus::InsufficientArc;return result;}
    }
    // Local implicit graph in the query's normal/tangent frame, scaled by h:
    // y = c0+c1*x+c2*(x*x+y*y)+c3*x^3+c4*x^4.
    // The isotropic quadratic reproduces circles without truncating their
    // square-root height. Cubic/quartic corrections retain variable curvature;
    // this is a local five-parameter fit, not a component-wide circle prior.
    // Keep the ordered branch, positive weights and the stable QR solve.
    const int columns=model==LocalModel::IsotropicQuadratic?3:5;
    Eigen::MatrixXd matrix(data.size(),columns);Eigen::VectorXd rhs(data.size());double weights=0;
    for (std::size_t i=0;i<data.size();++i) {
        const auto& v=data[i];const double w=std::sqrt(v.weight);weights+=v.weight;
        matrix(i,0)=w;matrix(i,1)=w*v.x;matrix(i,2)=w*(v.x*v.x+v.y*v.y);
        if(columns==5) {matrix(i,3)=w*v.x*v.x*v.x;matrix(i,4)=w*v.x*v.x*v.x*v.x;}rhs[i]=w*v.y;
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(matrix);qr.setThreshold(1e-10);
    if (qr.rank()!=columns) { result.status=FitStatus::SingularFit;return result; }
    const auto diagonal=qr.matrixR().topLeftCorner(columns,columns).diagonal().cwiseAbs().eval();
    result.scaledRankRatio=diagonal.minCoeff()/diagonal.maxCoeff();
    Eigen::VectorXd coefficients=Eigen::VectorXd::Zero(5);coefficients.head(columns)=qr.solve(rhs);
    if (result.scaledRankRatio<1e-8 || !coefficients.allFinite()) {
        result.status=FitStatus::SingularFit;return result;
    }
    // Select the graph root continuous with c0 as c2 -> 0, without subtracting
    // nearly equal numbers. A vanishing normal derivative is not a graph.
    const double discriminant=1-4*coefficients[2]*coefficients[0];
    if (!(discriminant>0) || !std::isfinite(discriminant)) { result.status=FitStatus::GeometryMismatch;return result; }
    const double normalDerivative=std::sqrt(discriminant);
    const double offset=2*coefficients[0]/(1+normalDerivative);
    result.scaledPositionError=std::abs(offset);
    // Compare a normal-distance residual, rather than an arbitrary implicit
    // polynomial scaling. Reject a fitted fold anywhere in the sampled arc.
    double squaredResidual=0;
    for (const auto& v:data) {
        const double gx=-coefficients[1]-2*coefficients[2]*v.x
            -3*coefficients[3]*v.x*v.x-4*coefficients[4]*v.x*v.x*v.x;
        const double gy=1-2*coefficients[2]*v.y,gradient=std::hypot(gx,gy);
        if (!(gy>0) || !std::isfinite(gradient)) {result.status=FitStatus::PoorOrientation;
            result.orientationCause=OrientationCause::FittedFold;return result;}
        const double residual=(v.y-coefficients[0]-coefficients[1]*v.x
            -coefficients[2]*(v.x*v.x+v.y*v.y)-coefficients[3]*v.x*v.x*v.x
            -coefficients[4]*v.x*v.x*v.x*v.x)/gradient;
        squaredResidual+=v.weight*residual*residual;
    }
    result.scaledResidual=std::sqrt(squaredResidual/weights);
    // These are stencil admissibility gates, NOT a curvature accuracy claim.
    if (result.scaledResidual>.025 || result.scaledPositionError>.05) {
        result.status=FitStatus::GeometryMismatch;return result;
    }
    const double slope=coefficients[1]/normalDerivative,norm=std::sqrt(1+slope*slope);
    const auto fittedNormal=(normal-tangent*slope)*(1/norm);
    if (dot(fittedNormal,normal)<.98) { result.status=FitStatus::PoorOrientation;
        result.orientationCause=OrientationCause::NormalMismatch;return result; }
    // Implicit differentiation at (0,offset): y''=2*c2*(1+y'^2)/F_y.
    const double meridional=-2*coefficients[2]/(h*normalDerivative*norm);
    result.scaledMeridionalMagnitude=std::abs(meridional)*h;
    if (result.scaledMeridionalMagnitude>.5) {
        result.refinementCause=RefinementCause::CurvatureRadius;
        result.status=FitStatus::RefinementRequired;return result;
    }
    const double azimuthal=p.rz?fittedNormal.x/query.point.x:0;
    if (!std::isfinite(meridional+azimuthal)) { result.status=FitStatus::SingularFit;return result; }
    result.normal=fittedNormal;result.meridional=meridional;result.azimuthal=azimuthal;
    result.curvature=meridional+azimuthal;result.status=FitStatus::Valid;return result;
}
InterfaceCurvature::Source InterfaceCurvature::evaluatePointwise(std::size_t surface,int graph,double fraction) const
{
    impl_->validate();Source result;
    const auto& p=*impl_;const auto& chain=p.state.topology().surfaces().at(surface).geometry.nativeChain();
    if(p.eventCell(chain.graphs().at(graph).nativeCellIndex)) {
        const auto frame=chain.frameAt(graph,fraction);
        result.point=frame.point;result.normal=frame.normal;result.meridional=frame.curvature;
        result.azimuthal=p.rz?(frame.point.x==0?frame.curvature:frame.normal.x/frame.point.x):0;
        result.curvature=result.meridional+result.azimuthal;
        if(std::isfinite(result.curvature))result.method=Method::EventModel;
        return result;
    }
    result.height=impl_->heights.evaluate(surface,graph,fraction);
    result.point=result.height.point;
    if (result.height.valid()) {
        result.method=Method::HeightFunction;result.normal=result.height.normal;
        result.meridional=result.height.meridional;result.azimuthal=result.height.azimuthal;
        result.curvature=result.height.curvature;return result;
    }
    result.fit=localFit(surface,graph,fraction);
    if (result.fit.valid()) {
        result.method=Method::LocalFit;result.normal=result.fit.normal;
        result.meridional=result.fit.meridional;result.azimuthal=result.fit.azimuthal;
        result.curvature=result.fit.curvature;
    }
    return result;
}
InterfaceCurvature::Sample InterfaceCurvature::evaluate(std::size_t surface,int graph,double fraction) const
{
    const auto& p=*impl_;p.validate();
    const auto& chain=p.state.topology().surfaces().at(surface).geometry.nativeChain();
    const auto query=chain.frameAt(graph,fraction);
    const int count=static_cast<int>(chain.graphs().size());
    if(p.eventCell(chain.graphs().at(graph).nativeCellIndex)) {
        Sample result;static_cast<Source&>(result)=evaluatePointwise(surface,graph,fraction);
        result.primaryMethod=result.method;result.primaryPoint=result.point;return result;
    }
    const auto arc=[&](int g,double a,double b) {
        if (a==0 && b==.5) return p.halfLengths[surface][g][0];
        if (a==.5 && b==1) return p.halfLengths[surface][g][1];
        // A remapped query can differ from a retained knot by a few ulps.
        // Such a span need not have a representable quadrature interior.
        // Its endpoint chord has the same arc-length limit (cubic error in
        // span length), and is used only for this interpolation weight.
        if(a<b && b-a<=64*std::numeric_limits<double>::epsilon())
            return length(chain.frameAt(g,b).point-chain.frameAt(g,a).point);
        double sum=0;
        if (a<b) for (const auto& q:chain.quadrature(g,a,b)) sum+=q.weight;
        return sum;
    };
    const auto knot=[&](int g,double f) -> const Source& {
        return f==.5?p.midpoints[surface][g]:p.endpoints[surface][f==0?0:1];
    };
    Sample result;
    if (fraction==.5 || (!chain.closed() && ((graph==0 && fraction==0)
        || (graph==count-1 && fraction==1)))) {
        static_cast<Source&>(result)=knot(graph,fraction);
        result.primaryMethod=result.method;result.primaryPoint=result.point;return result;
    }
    int left=graph,right=graph;double leftFraction=.5,rightFraction=.5,travel=0,total=0;
    if (fraction<.5) {
        if (graph==0 && !chain.closed()) leftFraction=0;
        else {left=(graph+count-1)%count;travel=arc(left,.5,1);}
        total=travel+arc(graph,0,.5);travel+=arc(graph,0,fraction);
    } else {
        if (graph==count-1 && !chain.closed()) rightFraction=1;
        else {right=(graph+1)%count;total=arc(right,0,.5);}
        total+=arc(graph,.5,1);travel=arc(graph,.5,fraction);
    }
    const Source* first=&knot(left,leftFraction);const Source* second=&knot(right,rightFraction);
    // Keep the ordinary region's fixed-knot representation all the way to
    // the model boundary. Its own endpoint replaces the neighboring model
    // midpoint; no model curvature is blended into the ordinary graph and no
    // new pointwise fit is attempted at each projected force query.
    if(first->method==Method::EventModel) {
        first=&p.modelBoundaries[surface][graph][0];
        total=arc(graph,0,.5);travel=arc(graph,0,fraction);
    }
    if(second->method==Method::EventModel) {
        second=&p.modelBoundaries[surface][graph][1];
        total=arc(graph,.5,1);travel=arc(graph,.5,fraction);
    }
    const auto& a=*first;const auto& b=*second;
    static_cast<Source&>(result)=a;result.primaryMethod=a.method;result.primaryPoint=a.point;
    result.point=query.point;result.secondary=b;result.interpolated=true;
    // No crossing a component, skipping an unavailable knot, or extrapolating
    // over an unresolved turn. Both estimates retain their own branch gates.
    if (!a.valid() || !b.valid() || !(total>0) || !std::isfinite(total)
        || travel<0 || travel>total || dot(a.normal,b.normal)<=0) {
        result.method=Method::Unavailable;
        result.curvature=result.meridional=result.azimuthal=std::numeric_limits<double>::quiet_NaN();
        return result;
    }
    const double w=travel/total;result.secondaryWeight=w;
    double cw=w;
    if(p.rz && a.regularized && b.regularized && !chain.closed()) {
        if(leftFraction==0 && chain.beginEndpointKind()==LevelSetInterfacePanels::EndpointKind::Axis)cw=w*w;
        if(rightFraction==1 && chain.endEndpointKind()==LevelSetInterfacePanels::EndpointKind::Axis)cw=1-(1-w)*(1-w);
    }
    result.secondaryCurvatureWeight=cw;
    const auto normal=a.normal*(1-w)+b.normal*w;
    result.normal=normal*(1/length(normal));
    result.meridional=(1-cw)*a.meridional+cw*b.meridional;
    result.regularized=a.regularized && b.regularized;
    result.azimuthal=(1-cw)*a.azimuthal+cw*b.azimuthal;
    result.curvature=result.meridional+result.azimuthal;
    result.method=a.method==b.method?a.method:Method::Mixed;
    return result;
}
InterfaceCurvature::Extension InterfaceCurvature::extend(Vector2D query,double halfWidth) const
{ return extendImpl(query,halfWidth,false); }
InterfaceCurvature::Extension InterfaceCurvature::extendConservativeFace(int nativeFace,double halfWidth) const
{
    const auto face=impl_->domain.face(impl_->mesh,nativeFace);
    std::vector<int> owners;
    for(int c:{face.negativeCell,face.positiveCell})if(c>=0)owners.push_back(c);
    return extendConservativeFace(nativeFace,halfWidth,owners);
}
InterfaceCurvature::Extension InterfaceCurvature::extendConservativeFace(int nativeFace,double halfWidth,
    const std::vector<int>& pressureSupport,std::uint64_t liquidId) const
{
    const auto& p=*impl_;p.validate();
    const auto face=p.domain.face(p.mesh,nativeFace);
    bool tail=face.negativeCell>=0 && face.positiveCell>=0;int samePhase=0;double minorityMax=0;
    const std::vector<long double>* volumes=&p.inventory->liquidVolumes();
    const std::vector<int>* componentCells=nullptr;
    if(liquidId) {
        const auto& components=p.state.componentInventories();
        const auto material=std::find_if(components.begin(),components.end(),
            [&](const auto& component){return component.liquidId==liquidId;});
        const auto& liquids=p.state.topology().liquids();
        const auto liquid=std::find_if(liquids.begin(),liquids.end(),
            [&](const auto& component){return component.id==liquidId;});
        if(material==components.end() || liquid==liquids.end())
            throw std::invalid_argument("scoped curvature requires transported child inventory and geometry");
        volumes=&material->liquidVolumes;componentCells=&liquid->cells;
    }
    for(int owner:{face.negativeCell,face.positiveCell})if(owner>=0
        && std::find(pressureSupport.begin(),pressureSupport.end(),owner)==pressureSupport.end())
        throw std::invalid_argument("capillary pressure support omitted a native face owner");
    for (int cell:pressureSupport) {
        if(cell<0 || cell>=p.mesh.cellIndexNum() || p.domain.cellRowByIndex()[cell]<0)
            throw std::invalid_argument("capillary pressure support contains an inactive cell");
        const auto volume=volumes->at(cell),capacity=p.inventory->capacities()[cell];
        const auto minority=std::min(volume,capacity-volume);
        const int phase=volume<=capacity/2?1:-1;
        // A cell occupied by another child is geometrically gas for this
        // component. Certify the SAME C_k whose pressure gradient is loaded;
        // total C would incorrectly turn a sub-resolution material tail into
        // a resolved interface. The original tail budget/branch gates remain.
        const int geometricPhase=componentCells
            && std::find(componentCells->begin(),componentCells->end(),cell)==componentCells->end()
            ?1:p.state.topology().cellPhase()[cell];
        tail=tail && minority<=Quantum::Geometry::PlanarCellCut::relativeVolumeResolution*capacity
            && geometricPhase==phase && (!samePhase || samePhase==phase);
        samePhase=phase;minorityMax=std::max(minorityMax,static_cast<double>(minority/capacity));
    }
    auto result=extendImpl(face.center,halfWidth,tail && minorityMax>0,liquidId);
    if (result.unresolvedInventoryTail) result.maxMinorityFraction=minorityMax;
    return result;
}
InterfaceCurvature::Extension InterfaceCurvature::extendImpl(Vector2D query,double halfWidth,bool permitTail,std::uint64_t liquidId) const
{
    const auto& p=*impl_;p.validate();
    if (!std::isfinite(query.x) || !std::isfinite(query.y) || !std::isfinite(halfWidth) || halfWidth<=0)
        throw std::invalid_argument("curvature extension requires a finite point and positive physical width");
    Extension result;result.query=query;
    struct Candidate { std::size_t surface;LevelSetInterfacePanels::LocalProjection projection; };
    std::vector<Candidate> candidates;
    const auto& surfaces=p.state.topology().surfaces();
    for (std::size_t s=0;s<surfaces.size();++s)if(!liquidId || surfaces[s].liquidId==liquidId)
        for (auto candidate:surfaces[s].geometry.nativeChain().loadProjectionCandidates(query))
            candidates.push_back({s,std::move(candidate)});
    if (candidates.empty()) return result;
    const auto best=std::min_element(candidates.begin(),candidates.end(),[](const Candidate& a,const Candidate& b) {
        return a.projection.distance<b.projection.distance;
    });
    result.surface=best->surface;result.surfaceId=surfaces[result.surface].id;
    result.projection=best->projection;result.distance=best->projection.distance;
    result.signedDistance=p.source->value(query)<0?-result.distance:result.distance;
    // Inventory resolution is certified by C and the complete pressure
    // support, not by the arbitrary projection search width. A retained
    // pure-phase tail keeps the same chart policy on either side of the band.
    result.unresolvedInventoryTail=permitTail;
    if (result.distance>halfWidth) {
        if (!permitTail) { result.status=ExtensionStatus::OutsideBand;return result; }
    }
    const auto spacing=surfaces[result.surface].geometry.nativeChain().minimumActiveSpacing();
    const double h=std::min(spacing.x,spacing.y);
    if (result.unresolvedInventoryTail)
        result.sample=evaluate(result.surface,result.projection.panelIndex,result.projection.nativeGraphFraction);
    const auto sameTailChart=[&](const Candidate& other) {
        if (!result.unresolvedInventoryTail || other.surface!=result.surface || !result.sample.valid()
            || result.sample.method!=Method::HeightFunction) return false;
        const auto containsBoth=[&](const InterfaceHeightCurvature::Sample& height) {
            return height.valid()
                && std::find(height.branchGraphs.begin(),height.branchGraphs.end(),result.projection.panelIndex)!=height.branchGraphs.end()
                && std::find(height.branchGraphs.begin(),height.branchGraphs.end(),other.projection.panelIndex)!=height.branchGraphs.end();
        };
        // A common admitted column stencil is a local connected height arc,
        // not merely a whole-surface ID or a Euclidean distance threshold.
        // Both knots of an interpolation must give the same chart witness.
        return containsBoth(result.sample.height) && (!result.sample.interpolated
            || (result.sample.secondary.height.normalAxis==result.sample.height.normalAxis
                && containsBoth(result.sample.secondary.height)));
    };
    const auto sameEventChart=[&](const Candidate& other) {
        if(other.surface!=result.surface)return false;
        const auto& chain=surfaces[result.surface].geometry.nativeChain();
        const int first=result.projection.panelIndex,last=other.projection.panelIndex;
        const auto radialArc=[&](int g,double from,double to) {
            if(from==to)return true;
            const auto a=chain.frameAt(g,from),b=chain.frameAt(g,to);
            if(!(b.point.y>a.point.y) || !(a.normal.x>0) || !(b.normal.x>0))return false;
            return p.radialCharts[result.surface][g]
                || chain.zeroSetPartialResolved(g,0,from,to);
        };
        if(first==last && p.eventCell(chain.graphs()[first].nativeCellIndex)) {
            if(p.radialCharts[result.surface][first])return true;
            return radialArc(first,std::min(result.projection.nativeGraphFraction,other.projection.nativeGraphFraction),
                std::max(result.projection.nativeGraphFraction,other.projection.nativeGraphFraction));
        }
        if(p.state.eventCurvatureCells().empty())return false;
        if(!p.rz || chain.closed())return false;
        if(!result.unresolvedInventoryTail) {
            const auto liquid=surfaces[result.surface].liquidId;
            const auto& events=p.state.topology().events();
            if(!std::any_of(events.begin(),events.end(),[&](const auto& event) {
                return (event.kind==InterfaceTopologyEvent::Kind::SectionDetachment
                    || (event.kind==InterfaceTopologyEvent::Kind::Split && !event.interiorNecks.empty()))
                    && std::find(event.children.begin(),event.children.end(),liquid)!=event.children.end();
            }))return false;
        }
        // Inside the force band, this policy is limited to one wholly modeled
        // connecting arc, as for the same-graph event case above. Artificial
        // graph boundaries must not change that local model's branch policy.
        // A remote unresolved tail may also cross the MODEL'S boundary
        // without extending model curvature into the resolved arc.
        // Geometry, not provider tags, certifies one connected radial chart.
        // The actual candidates must each have an available provider below;
        // no curvature is interpolated through the intervening segments.
        // A surface ID alone cannot merge two sides of a folded curve.
        const bool forward=first<last || (first==last
            && result.projection.nativeGraphFraction<=other.projection.nativeGraphFraction);
        const auto& begin=forward?result.projection:other.projection;
        const auto& end=forward?other.projection:result.projection;
        for(int g=begin.panelIndex;g<=end.panelIndex;++g) {
            if(!result.unresolvedInventoryTail && !p.eventCell(chain.graphs()[g].nativeCellIndex))return false;
            // Prove only the connecting arc. The unused remainder of an
            // endpoint graph may contain the pole or a turn; including it
            // makes acceptance depend on the integration partition.
            if(!radialArc(g,g==begin.panelIndex?begin.nativeGraphFraction:0,
                g==end.panelIndex?end.nativeGraphFraction:1))return false;
        }
        return true;
    };
    // Ignore duplicate roots/nearby stationary coordinates on the same local
    // arc; a distinct point with almost the same distance is a force-band ridge.
    // The 0.01h tie layer is a numerical resolution policy, not an exact medial axis.
    for (const auto& candidate:candidates) {
        if(candidate.projection.distance<=result.distance+.01*h && sameEventChart(candidate)) {
            const auto a=evaluate(result.surface,result.projection.panelIndex,result.projection.nativeGraphFraction);
            const auto b=evaluate(candidate.surface,candidate.projection.panelIndex,candidate.projection.nativeGraphFraction);
            if(a.valid() && b.valid() && dot(a.normal,b.normal)>0) {
                if(length(candidate.projection.frame.point-best->projection.frame.point)>.25*h) {
                    ++result.eventChartCandidates;
                    result.eventCurvatureSpread=std::max(result.eventCurvatureSpread,std::abs(a.curvature-b.curvature));
                }
                continue;
            }
        }
        if (candidate.projection.distance<=result.distance+.01*h && sameTailChart(candidate)) {
            if (length(candidate.projection.frame.point-best->projection.frame.point)>.25*h)
                ++result.sameChartTailCandidates;
            continue;
        }
        if (length(candidate.projection.frame.point-best->projection.frame.point)>.25*h
            && candidate.projection.distance<result.competingDistance) {
            result.competingDistance=candidate.projection.distance;
            result.competingSurface=candidate.surface;result.competingProjection=candidate.projection;
        }
    }
    if (result.competingDistance<=result.distance+.01*h) {
        if(result.unresolvedInventoryTail && result.sample.valid() && std::isfinite(result.sample.curvature)) {
            double lower=result.sample.curvature,upper=lower;
            for(const auto& candidate:candidates) {
                if(candidate.projection.distance>result.distance+.01*h
                    || length(candidate.projection.frame.point-best->projection.frame.point)<=.25*h)continue;
                const auto sample=evaluate(candidate.surface,candidate.projection.panelIndex,candidate.projection.nativeGraphFraction);
                if(!sample.valid() || !std::isfinite(sample.curvature)) {
                    result.tailAlternatives.clear();break;
                }
                result.tailAlternatives.push_back({candidate.surface,candidate.projection,sample.curvature});
                lower=std::min(lower,sample.curvature);upper=std::max(upper,sample.curvature);
            }
            if(!result.tailAlternatives.empty())result.tailCurvatureSpread=upper-lower;
        }
        result.status=ExtensionStatus::AmbiguousBranch;return result;
    }
    if (!result.unresolvedInventoryTail)
        result.sample=evaluate(result.surface,result.projection.panelIndex,result.projection.nativeGraphFraction);
    result.status=result.sample.valid()?ExtensionStatus::Valid:ExtensionStatus::CurvatureUnavailable;
    return result;
}
