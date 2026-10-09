#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroFittedVelocity.h"
#include "Modules/Hydrodynamics/HydroSplineVelocity.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpline.h"
#include "Modules/Hydrodynamics/HydroBoundaryVelocity.h"
#include "Modules/Hydrodynamics/Detail/HydroNativeVelocityTrace.h"
#include "Modules/Interface/InterfacePhaseGeometry.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/LevelSetGaussTrace.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Mesh/CartesianLeafLocator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace HydroInterfaceTransport
{
namespace
{
void require(bool valid, const char* reason)
{ if (!valid) throw std::runtime_error(std::string("[HydroInterfaceTransport] ")+reason); }
const HydroBoundaryPatch* patchAt(const NativeLeafFace& face, const std::vector<HydroBoundaryPatch>& patches)
{
    const auto side=face.horizontal ? (face.outwardSign<0?HydroBoundarySide::South:HydroBoundarySide::North)
        : (face.outwardSign<0?HydroBoundarySide::West:HydroBoundarySide::East);
    const double position=face.horizontal?face.center.x:face.center.y;
    for(const auto& patch:patches)
        if(patch.side==side && position>patch.start && position<patch.end) return &patch;
    return nullptr;
}
void completeDomain(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const std::vector<HydroBoundaryPatch>& patches)
{
    domain.validateMeshIdentity(mesh);
    require(domain.cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "complete velocity bridge requires every native leaf active");
    if(!domain.uniformLeaves())require(patches.empty(),"mixed complete velocity currently requires closed boundaries");
    else {const UniformLeafGrid grid(mesh);validateHydroBoundaryPatches(patches,mesh.minCoord,mesh.maxCoord,grid.spacing(),domain.axisymmetric());}
}
}

HydroFittedVelocitySampler completeVelocity(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const VEdgeField<double>& ux,const HEdgeField<double>& uy,const std::vector<HydroBoundaryPatch>& patches)
{
    completeDomain(mesh,domain,patches);
    if(!domain.uniformLeaves()) {
        require(ux.value && uy.value && ux.size()==mesh.vEdgeIndexNum() && uy.size()==mesh.hEdgeIndexNum(),"mixed native velocity storage mismatch");
        std::vector<HydroFittedVelocitySampler::VertexTrace> traces(mesh.vertIndexNum());
        for(int v=0;v<mesh.vertIndexNum();++v) {
            const auto& vertex=mesh.verts(mesh.vertIndexToLabel[v]);std::vector<int> seed;
            for(int label:vertex.vertToCell)if(label>=0)seed.push_back(mesh.cells(label).index());
            const auto support=HydroFittedDetail::nativeTraceSupport(mesh,domain,seed);
            for(int component=0;component<2;++component) {
                const auto fit=HydroFittedDetail::nativeVelocityTrace(mesh,domain,support,component,vertex.coord());
                const auto evaluate=[&](const auto& trace) {
                    double value=0;for(const auto& term:trace)value+=term.second*(component?uy.value[term.first-mesh.vEdgeIndexNum()]:ux.value[term.first]);
                    require(std::isfinite(value),"nonfinite mixed velocity trace");return value;
                };
                if(component){traces[v].velocity.y=evaluate(fit.value);traces[v].axialDerivative=evaluate(fit.gradient[1]);}
                else traces[v].velocity.x=evaluate(fit.value);
            }
        }
        for(const auto face:domain.faces(mesh))if(face.outwardSign!=0) {
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            const double normal=face.horizontal?uy.value[edge.index()]:ux.value[edge.index()];
            require(normal==0,"mixed closed velocity has nonzero boundary flux");
            for(int label:edge.edgeToVert) {
                auto& trace=traces[mesh.verts(label).index()];
                if(face.horizontal)trace.velocity.y=0;else trace.velocity.x=0;
            }
        }
        return HydroFittedVelocitySampler(mesh,domain,ux,uy,traces);
    }
    const UniformLeafGrid grid(mesh);
    const auto native=hydroBoundaryAdvectionVelocity(mesh,ux,uy,patches,domain.axisymmetric() && mesh.minCoord.x==0);
    std::vector<HydroFittedVelocitySampler::VertexTrace> traces(mesh.vertIndexNum());
    for(int k=0;k<mesh.vertIndexNum();++k) traces[k].velocity=native(mesh.verts(mesh.vertIndexToLabel[k]).coord());
    const auto& vertices=grid.vertices();
    for(int j=0;j<=grid.ny();++j) for(int i=0;i<=grid.nx();++i)
    {
        const int low=j==0?0:j==grid.ny()?j-2:j-1;
        const auto& a=mesh.verts(vertices[i+(grid.nx()+1)*low]);
        const auto& b=mesh.verts(vertices[i+(grid.nx()+1)*(low+1)]);
        const auto& c=mesh.verts(vertices[i+(grid.nx()+1)*(low+2)]);
        const long double h0=static_cast<long double>(b.coord().y)-a.coord().y,
            h1=static_cast<long double>(c.coord().y)-b.coord().y;
        const long double s0=(static_cast<long double>(traces[b.index()].velocity.y)-traces[a.index()].velocity.y)/h0,
            s1=(static_cast<long double>(traces[c.index()].velocity.y)-traces[b.index()].velocity.y)/h1;
        const long double derivative=j==0?s0-h0*(s1-s0)/(h0+h1):j==grid.ny()?s1+h1*(s1-s0)/(h0+h1)
            :(h1*s0+h0*s1)/(h0+h1);
        traces[mesh.verts(vertices[i+(grid.nx()+1)*j]).index()].axialDerivative=static_cast<double>(derivative);
    }
    for(const auto face:domain.faces(mesh)) if(face.outwardSign!=0)
    {
        const auto* patch=patchAt(face,patches);
        if(patch && patch->type==HydroBoundaryType::PressureOutlet) continue;
        const double normal=patch?face.outwardSign*patch->value:0;
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        for(int label:edge.edgeToVert)
        {
            const auto v=traces[mesh.verts(label).index()].velocity;
            require((face.horizontal?v.y:v.x)==normal,
                "shared vertex cannot represent discontinuous prescribed normal boundary traces");
        }
    }
    return HydroFittedVelocitySampler(mesh,domain,ux,uy,traces);
}

InterfaceTransportDomain domainView(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const std::vector<HydroBoundaryPatch>& patches)
{
    completeDomain(mesh,domain,patches);
    return InterfaceTransportDomain(domain,[&mesh,&domain,patches](int row) {
        using Boundary=InterfaceTransportDomain::Boundary;
        const auto face=domain.face(mesh,row);
        if(face.outwardSign==0) return Boundary::Interior;
        const auto* patch=patchAt(face,patches);
        return patch && (patch->type==HydroBoundaryType::PressureOutlet || patch->value!=0)
            ?Boundary::SinglePhasePort:Boundary::Symmetry;
    });
}

InterfaceTransportDomain domainView(const HydroFittedDomain& domain)
{
    return InterfaceTransportDomain(domain.geometry(), [&domain](int row) {
        using Hydro = HydroFittedDomain::BoundaryType;
        using Interface = InterfaceTransportDomain::Boundary;
        switch (domain.faces().at(row).boundary.type)
        {
        case Hydro::Interior: return Interface::Interior;
        case Hydro::NoSlipWall: return Interface::StationaryWall;
        case Hydro::Symmetry: return Interface::Symmetry;
        case Hydro::PressureReservoir: case Hydro::PressureOutlet: return Interface::SinglePhasePort;
        }
        throw std::invalid_argument("Hydro interface transport received an unsupported boundary");
    });
}

InterfaceNormalFluxView normalFlux(const HydroFittedVelocitySampler& velocity)
{
    InterfaceNormalFluxView::Queries queries;
    queries.chart = [&velocity](int row) {
        const auto chart = velocity.faceChart(row);
        return InterfaceNormalFluxView::FaceChart{chart.begin,chart.midpoint,chart.end};
    };
    queries.integral = [&velocity](int row,int half,double a,double b,double c,double d) {
        return velocity.integratedNormalFlux(row,half,a,b,c,d);
    };
    queries.phaseIntegral = [&velocity](int row,int half,const LevelSetMappedSegment& segment) {
        return velocity.integratedNormalFlux(row,half,segment);
    };
    queries.density = [&velocity](int row,int half,bool fromEnd,long double offset) {
        return velocity.normalFluxDensity(row,half,fromEnd,offset);
    };
    queries.roots = [&velocity](int row,int half) {return velocity.normalFluxRoots(row,half);};
    queries.polynomial = [&velocity](int row,int half) {return velocity.normalFluxPolynomial(row,half);};
    return InterfaceNormalFluxView(std::move(queries));
}

InterfaceParametricTransport::Flow frozenParametricFlow(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const HydroFittedVelocitySampler& velocity,int blocks)
{
    completeDomain(mesh,domain,{});
    require(!domain.axisymmetric() && domain.uniformLeaves(),
        "parametric flow currently requires a complete uniform XY domain");
    require(blocks>0,"parametric flow requires positive Gauss block count");
    // Smooth pieces are defined by the actual native half-face coordinates,
    // exactly as in the velocity reconstruction; there is no resampled lattice.
    using Patch=LevelSetGaussTrace::SmoothPatch;
    std::vector<Patch> patches;
    std::vector<std::vector<std::size_t>> byCell(mesh.cellIndexNum());
    for(int index:domain.cellIndexByRow()) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[index]);
        const auto lower=mesh.verts(cell.cellToVert[0]).coord(),upper=mesh.verts(cell.cellToVert[2]).coord();
        std::vector<double> xs{lower.x,upper.x},ys{lower.y,upper.y};
        for(int side=0;side<4;++side) {
            const bool horizontal=side>=2;
            const auto& labels=horizontal?cell.cellToHEdge[side%2]:cell.cellToVEdge[side%2];
            for(int part=0;part<2;++part)if(const int label=labels[part];label>=0) {
                const auto& edge=horizontal?mesh.hEdges(label):mesh.vEdges(label);
                const auto chart=velocity.faceChart(edge.index()+(horizontal?mesh.vEdgeIndexNum():0));
                const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
                const auto m=edge.midpointCoord();
                require(chart.begin.x==a.x && chart.begin.y==a.y && chart.end.x==b.x && chart.end.y==b.y &&
                    chart.midpoint.x==m.x && chart.midpoint.y==m.y,"parametric velocity chart differs from native mesh");
                for(auto p:{chart.begin,chart.midpoint,chart.end})(horizontal?xs:ys).push_back(horizontal?p.x:p.y);
            }
        }
        const auto unique=[](auto& coordinates){std::sort(coordinates.begin(),coordinates.end());
            coordinates.erase(std::unique(coordinates.begin(),coordinates.end()),coordinates.end());};
        unique(xs);unique(ys);
        for(std::size_t i=1;i<xs.size();++i)for(std::size_t j=1;j<ys.size();++j) {
            byCell[index].push_back(patches.size());patches.push_back({{xs[i-1],ys[j-1]},{xs[i],ys[j]}});
        }
    }
    const LevelSetGaussTrace::PatchQuery select=[patches=std::move(patches),byCell=std::move(byCell),locator=CartesianLeafLocator(mesh)]
        (Vector2D p,Vector2D direction) {
        std::size_t selected=patches.size();const auto cells=locator.incident(p);
        for(int i=0;i<cells.count;++i)for(auto id:byCell[cells.indices[i]]) {
            const auto& patch=patches[id];
            if(id>=selected || !LevelSetGaussTrace::contains(patch,p))continue;
            if((p.x==patch.lower.x && direction.x<0) || (p.x==patch.upper.x && direction.x>0) ||
               (p.y==patch.lower.y && direction.y<0) || (p.y==patch.upper.y && direction.y>0))continue;
            selected=id;
        }
        require(selected<patches.size(),"no native smooth patch in the flow direction");
        return patches[selected];
    };
    return [&velocity,select,blocks](Vector2D p,double duration) {
        const auto map=LevelSetGaussTrace::advance([&velocity](Vector2D x) {
            const auto jet=velocity.volumeFirstJet(x);
            return LevelSetGaussTrace::VelocityJet{jet.value,jet.dx,jet.dy};
        },select,p,duration,blocks);
        return InterfaceParametricTransport::MapJet{map.point,
            {map.derivative(0,0),map.derivative(1,0)},{map.derivative(0,1),map.derivative(1,1)},
            map.timeDerivative,map.displacement,map.coordinateResolution};
    };
}

InterfaceNormalFluxView normalFlux(const HydroSplineVelocity& velocity)
{
    InterfaceNormalFluxView::Queries q;
    q.chart=[&velocity](int row){const auto& c=velocity.faceFlux(row).chart();return InterfaceNormalFluxView::FaceChart{c.begin,c.midpoint,c.end};};
    q.integral=[&velocity](int row,int half,double a,double b,double c,double d){return velocity.faceFlux(row).integrate(half,a,b,c,d);};
    q.phaseIntegral=[&velocity](int row,int half,const LevelSetMappedSegment& segment){return velocity.faceFlux(row).integrate(half,segment);};
    q.density=[&velocity](int row,int half,bool end,long double t){return velocity.faceFlux(row).density(half,end,t);};
    q.roots=[&velocity](int row,int half){return velocity.faceFlux(row).densityRoots(half);};
    q.polynomial=[&velocity](int row,int half){return velocity.faceFlux(row).densityPolynomial(half);};
    return InterfaceNormalFluxView(std::move(q));
}

template<class Spline> static auto transportJet(const Spline& velocity,Vector2D p){return velocity.firstJet(p);}
static auto transportJet(const Quantum::Hydrodynamics::HydroCompatibleMomentSpline::NativeVelocity& velocity,Vector2D p){return velocity.transportJet(p);}
template<class Spline> static LevelSetGaussTrace::SmoothPatch transportDomain(const Spline& velocity){return {velocity.lower(),velocity.upper()};}
static LevelSetGaussTrace::SmoothPatch transportDomain(const Quantum::Hydrodynamics::HydroCompatibleMomentSpline::NativeVelocity& velocity){return {velocity.transportLower(),velocity.transportUpper()};}
template<class Spline>
static InterfaceParametricTransport::Flow frozenSplineFlow(const Spline& velocity,int blocks)
{
    require(blocks>0 && blocks<=4096,"spline flow requires a valid Gauss block count");
    const auto chart=transportDomain(velocity);
    return [&velocity,chart,blocks](Vector2D p,double duration) {
        const auto map=LevelSetGaussTrace::advance([&velocity](Vector2D x){const auto j=transportJet(velocity,x);
            return LevelSetGaussTrace::VelocityJet{j.value,j.dx,j.dy};},[chart](Vector2D,Vector2D){return chart;},p,duration,blocks);
        return InterfaceParametricTransport::MapJet{map.point,{map.derivative(0,0),map.derivative(1,0)},
            {map.derivative(0,1),map.derivative(1,1)},map.timeDerivative,map.displacement,map.coordinateResolution};
    };
}

template<class Spline>
static InterfaceParametricTransport::Flow timeLinearSplineFlow(const Spline& previous,
    const Spline& current,double stepDuration,int blocks)
{
    require(std::isfinite(stepDuration) && stepDuration>0 && blocks>0 && blocks<=4096,
        "time-linear spline flow requires positive physical duration and valid Gauss blocks");
    require(previous.sameLayout(current),"time-linear spline snapshots have different native layouts");
    const auto chart=transportDomain(previous);
    return [&previous,&current,stepDuration,chart,blocks](Vector2D p,double elapsed) {
        require(std::isfinite(elapsed) && elapsed>=0 && elapsed<=stepDuration,"spline flow query is outside its physical time interval");
        const auto map=LevelSetGaussTrace::advanceTimeDependent([&](Vector2D x,double t){
            const auto a=transportJet(previous,x),b=transportJet(current,x);const double fraction=t/stepDuration;
            return LevelSetGaussTrace::TimeVelocityJet{a.value+(b.value-a.value)*fraction,
                a.dx+(b.dx-a.dx)*fraction,a.dy+(b.dy-a.dy)*fraction,(b.value-a.value)*(1/stepDuration)};
        },chart,p,elapsed,blocks);
        return InterfaceParametricTransport::MapJet{map.point,{map.derivative(0,0),map.derivative(1,0)},
            {map.derivative(0,1),map.derivative(1,1)},map.timeDerivative,map.displacement,map.coordinateResolution};
    };
}

InterfaceParametricTransport::Flow frozenParametricFlow(const HydroSplineVelocity& velocity,int blocks){return frozenSplineFlow(velocity,blocks);}
InterfaceParametricTransport::Flow timeLinearParametricFlow(const HydroSplineVelocity& previous,const HydroSplineVelocity& current,double duration,int blocks){return timeLinearSplineFlow(previous,current,duration,blocks);}
InterfaceNormalFluxView normalFlux(const Quantum::Hydrodynamics::HydroCompatibleMomentSpline& velocity){return velocity.normalFlux();}
InterfaceParametricTransport::Flow frozenParametricFlow(const Quantum::Hydrodynamics::HydroCompatibleMomentSpline& velocity,int blocks){return frozenSplineFlow(velocity,blocks);}
InterfaceParametricTransport::Flow frozenParametricFlow(const Quantum::Hydrodynamics::HydroCompatibleMomentSpline::NativeVelocity& velocity,int blocks){return frozenSplineFlow(velocity,blocks);}
InterfaceParametricTransport::Flow timeLinearParametricFlow(const Quantum::Hydrodynamics::HydroCompatibleMomentSpline::NativeVelocity& previous,
    const Quantum::Hydrodynamics::HydroCompatibleMomentSpline::NativeVelocity& current,double duration,int blocks){
    require(std::abs(current.time()-previous.time()-duration)<=32*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(current.time()),std::abs(previous.time()),std::abs(duration)}),
        "complete native velocity snapshots do not span the physical step");
    return timeLinearSplineFlow(previous,current,duration,blocks);
}
InterfaceParametricTransport::Flow timeLinearParametricFlow(const Quantum::Hydrodynamics::HydroCompatibleMomentSpline& previous,
    const Quantum::Hydrodynamics::HydroCompatibleMomentSpline& current,double duration,int blocks){
    require(std::abs(current.time()-previous.time()-duration)<=32*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(current.time()),std::abs(previous.time()),std::abs(duration)}),
        "moment spline snapshots do not span the physical step");
    return timeLinearSplineFlow(previous,current,duration,blocks);
}

long double outwardPhaseFlux(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& snapshot,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy)
{
    InterfacePhaseGeometry::validateLayout(mesh, domain.geometry(), snapshot);
    require(ux.value && uy.value && ux.size() == mesh.vEdgeIndexNum() && uy.size() == mesh.hEdgeIndexNum(),
        "native velocity storage mismatch");
    long double flux = 0.0L;
    for (const auto& face : domain.faces())
    {
        if (face.boundary.type != HydroFittedDomain::BoundaryType::PressureReservoir
            && face.boundary.type != HydroFittedDomain::BoundaryType::PressureOutlet) continue;
        const int owner = face.negativeCell >= 0 ? face.negativeCell : face.positiveCell;
        require(owner >= 0 && owner < mesh.cellIndexNum()
            && (face.negativeCell < 0) != (face.positiveCell < 0), "pressure face lacks one unique active owner");
        const int sign = snapshot.closedCellSign(mesh.cells(mesh.cellIndexToLabel[owner]).centerCoord());
        require(sign != 0, "pressure owner cell is not certified single phase");
        if (sign > 0) continue;
        const int k = face.horizontal ? face.nativeRow - mesh.vEdgeIndexNum() : face.nativeRow;
        require(k >= 0 && k < (face.horizontal ? uy.size() : ux.size())
            && (face.horizontal ? uy.label(k) : ux.label(k)) == face.edgeLabel,
            "pressure native velocity index/label mismatch");
        const double u = face.horizontal ? uy.value[k] : ux.value[k];
        require(std::isfinite(u), "nonfinite negative-phase pressure velocity");
        flux += static_cast<long double>(face.outwardSign) * face.area * u;
    }
    require(std::isfinite(flux), "nonfinite negative-phase boundary volume rate");
    return flux;
}

}
