#include "Modules/Hydrodynamics/HydroUniformLoad.h"
#include "Modules/Hydrodynamics/HydroInterfaceLoadGeometry.h"
#include "Modules/Hydrodynamics/HydroConservativePressure.h"
#include "Modules/Hydrodynamics/HydroSplineVelocity.h"
#include "Modules/Hydrodynamics/Detail/ConservativeSurfaceMoment.h"
#include "Modules/Interface/Compatibility/InterfaceCurvature.h"
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>
#include <exception>
#include <optional>
#include <limits>
#include <Eigen/Cholesky>

namespace Quantum::Hydrodynamics
{
namespace {
void require(bool ok, const char* reason)
{ if (!ok) throw std::invalid_argument(std::string("[HydroUniformLoad] ") + reason); }
double checked(long double value, const char* reason)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result), reason);
    return result;
}
void cellAverages(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    bool axisymmetric, HydroUniformLoad& result)
{
    for (const auto& face : graph.segments())
    {
        const int owner = mesh.cells(face.K).index();
        if (face.isHorizontal)
            result.forceY.value[owner] += 0.5 * result.faceForceY.value[mesh.hEdges(face.edgeLabel).index()];
        else if (axisymmetric)
            result.forceX.value[owner] += 0.5 * (face.faceCenter.x / mesh.cells(face.K).centerCoord().x)
                * result.faceForceX.value[mesh.vEdges(face.edgeLabel).index()];
        else result.forceX.value[owner] += 0.5 * result.faceForceX.value[mesh.vEdges(face.edgeLabel).index()];
    }
}
void integrate(const QuadtreeMesh& mesh,const FaceSegmentGraph& graph,bool axisymmetric,HydroUniformLoad& result)
{
    cellAverages(mesh, graph, axisymmetric, result);
    long double totalX = 0.0L, totalY = 0.0L;
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const double fx = result.forceX.value[i], fy = result.forceY.value[i];
        require(std::isfinite(fx) && std::isfinite(fy), "non-finite cell force");
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        const double area = CartesianMeasure::cell(cell.area(), cell.centerCoord().x, axisymmetric);
        totalX += static_cast<long double>(fx) * area;
        totalY += static_cast<long double>(fy) * area;
        result.maxForce = std::max(result.maxForce, std::hypot(fx, fy));
    }
    result.totalForceX = checked(totalX, "integrated interface force x overflow");
    result.totalForceY = checked(totalY, "integrated interface force y overflow");
    require(std::isfinite(result.maxForce), "non-finite maximum force");
}
void integrateNative(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroConservativePressure& pressure,HydroUniformLoad& result)
{
    long double totalX=0,totalY=0;
    for(const auto face:domain.faces(mesh)) {
        const auto& volume=pressure.duals().volumes()[face.nativeRow];
        const double force=face.horizontal?result.faceForceY.value[face.nativeRow-mesh.vEdgeIndexNum()]
            :result.faceForceX.value[face.nativeRow];
        (face.horizontal?totalY:totalX)+=static_cast<long double>(force)*volume.measure;
        for(int owner:{face.negativeCell,face.positiveCell})if(owner>=0) {
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[owner]);
            const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
            const Vector2D a{std::max(lo.x,volume.lower.x),std::max(lo.y,volume.lower.y)},
                b{std::min(hi.x,volume.upper.x),std::min(hi.y,volume.upper.y)};
            const double measure=CartesianMeasure::cell((b.x-a.x)*(b.y-a.y),(a.x+b.x)/2,domain.axisymmetric());
            (face.horizontal?result.forceY:result.forceX).value[owner]+=force*measure/domain.cellVolumes()[domain.cellRowByIndex()[owner]];
        }
    }
    for(int cell:domain.cellIndexByRow()) {
        const double magnitude=std::hypot(result.forceX.value[cell],result.forceY.value[cell]);
        require(std::isfinite(magnitude),"nonfinite native capillary cell force");result.maxForce=std::max(result.maxForce,magnitude);
    }
    result.totalForceX=checked(totalX,"native capillary total x overflow");
    result.totalForceY=checked(totalY,"native capillary total y overflow");
}
template<std::size_t N,class Function>
std::array<long double,N> surfaceIntegral(const Function& function,double lo,double hi,
    const std::array<double,N>& scale,ParametricCapillaryReport& report,int depth=0)
{
    namespace P=Quantum::Geometry::Parametric;
    const auto rule=[&](const auto& x,const auto& w){std::array<long double,N> result{};
        for(std::size_t k=0;k<x.size();++k){const auto value=function(.5*(lo+hi+(hi-lo)*x[k]));++report.surfaceSamples;
            for(std::size_t i=0;i<N;++i)result[i]+=.5L*(hi-lo)*w[k]*value[i];}return result;};
    const auto coarse=rule(P::gx,P::gw),fine=rule(P::gaussEight().x,P::gaussEight().w);
    double error=0;
    for(std::size_t i=0;i<N;++i){const auto denominator=std::max(std::abs(fine[i]),static_cast<long double>((hi-lo)*scale[i]));
        require(std::isfinite(fine[i]) && denominator>0,"invalid parameter surface integral scale");
        error=std::max(error,static_cast<double>(std::abs(fine[i]-coarse[i])/denominator));}
    if(error<=512*std::numeric_limits<double>::epsilon()){report.quadratureError=std::max(report.quadratureError,error);return fine;}
    const double middle=.5*(lo+hi);require(depth<24 && middle>lo && middle<hi,"parameter surface integral is unresolved");
    auto result=surfaceIntegral<N>(function,lo,middle,scale,report,depth+1);const auto right=surfaceIntegral<N>(function,middle,hi,scale,report,depth+1);
    for(std::size_t i=0;i<N;++i)result[i]+=right[i];return result;
}
}

HydroUniformLoad::HydroUniformLoad(const QuadtreeMesh& mesh)
    : forceX(mesh, 0.0), forceY(mesh, 0.0), faceForceX(mesh, 0.0), faceForceY(mesh, 0.0) {}

HydroUniformLoad buildCompatibleScalarGradientLoad(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const HydroConservativePressure& pressure,
    const std::vector<double>& scalar,const std::vector<double>& coefficient)
{
    pressure.validate(mesh,domain);
    require(pressure.divergenceAdjoint(),"compatible normal load requires the adjoint pressure operator");
    require(scalar.size()==static_cast<std::size_t>(mesh.cellIndexNum())
        && coefficient.size()==static_cast<std::size_t>(domain.nativeFaceCount()),"compatible load storage mismatch");
    std::vector<double> boundary(domain.nativeFaceCount(),0);
    for(const auto face:domain.faces(mesh)) {
        require(std::isfinite(coefficient[face.nativeRow]),"nonfinite compatible normal coefficient");
        if(face.outwardSign!=0) {
            require(coefficient[face.nativeRow]==0,"closed normal boundary cannot receive a compatible interface load");
            boundary[face.nativeRow]=scalar[std::max(face.negativeCell,face.positiveCell)];
        }
    }
    const auto gradient=pressure.apply(mesh,domain,scalar,boundary);HydroUniformLoad result(mesh);
    for(const auto face:domain.faces(mesh)) {
        const double load=checked(static_cast<long double>(coefficient[face.nativeRow])*gradient[face.nativeRow]
            /pressure.duals().volumes()[face.nativeRow].measure,"compatible normal force density overflow");
        (face.horizontal?result.faceForceY.value[face.nativeRow-mesh.vEdgeIndexNum()]:result.faceForceX.value[face.nativeRow])=load;
    }
    integrateNative(mesh,domain,pressure,result);return result;
}

HydroUniformLoad buildParametricCapillaryLoad(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroConservativePressure& pressure,const InterfaceState& state,double sigma,ParametricCapillaryReport* output)
{
    namespace P=Quantum::Geometry::Parametric;
    InterfaceParametricMotion::validateDomain(mesh,domain,state);pressure.validate(mesh,domain);
    require(pressure.divergenceAdjoint() && !domain.axisymmetric(),"parameter capillarity requires compatible XY pressure");
    require(std::isfinite(sigma) && sigma>=0,"invalid parameter surface tension");
    const auto& candidate=state.parametricGeometry();const auto& curve=candidate.geometry();const auto& grid=candidate.grid();
    const double h=std::min(grid.width/grid.nx,grid.height/grid.ny);
    ParametricCapillaryReport report;std::array<long double,3> totals{};
    const auto tangentBound=[](const P::Arc& arc){double bound=0;for(int k=1;k<5;++k)bound+=k*P::norm(arc.p[k]);return bound;};
    for(const auto& arc:curve.arcs){for(int axis=0;axis<2;++axis){const auto extent=P::range(arc.coordinate(axis));
            require(extent[0]>(axis?grid.bottom:grid.left) && extent[1]<(axis?grid.bottom+grid.height:grid.left+grid.width),"parameter capillarity does not supply a wall contact law");}
        const double bound=tangentBound(arc);
        const auto values=surfaceIntegral<3>([&](double t){const auto q=P::sampleSurface(arc,t);
            return std::array<double,3>{q.speed,q.curvature*q.speed,q.curvature*q.curvature*q.speed};},0,1,{bound,bound/h,bound/(h*h)},report);
        for(int i=0;i<3;++i)totals[i]+=values[i];}
    require(totals[0]>0,"empty parameter surface");report.length=static_cast<double>(totals[0]);
    report.meanCurvature=static_cast<double>(totals[1]/totals[0]);report.curvatureSquaredIntegral=static_cast<double>(totals[2]);
    std::vector<double> alpha(mesh.cellIndexNum()),coefficient(domain.nativeFaceCount());
    const auto& cells=InterfaceParametricMotion::cellIndices(state);
    for(std::size_t slot=0;slot<cells.size();++slot){const int cell=cells[slot];alpha[cell]=candidate.inventory()[slot]/domain.cellVolumes()[domain.cellRowByIndex()[cell]];}
    for(const auto face:domain.faces(mesh))coefficient[face.nativeRow]=face.outwardSign?0:sigma*report.meanCurvature;
    auto result=buildCompatibleScalarGradientLoad(mesh,domain,pressure,alpha,coefficient);
    HydroSplineVelocityAdjoint adjoint(mesh,domain);
    const std::array<std::vector<double>,2> planes{{grid.lines(0),grid.lines(1)}};
    for(const auto& arc:curve.arcs){std::vector<double> cuts{0,1};
        for(int axis=0;axis<2;++axis)for(std::size_t line=1;line+1<planes[axis].size();++line){auto equation=arc.coordinate(axis);equation[0]-=planes[axis][line];
            for(double t:P::roots(equation,0,1))cuts.push_back(t);}
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        const double bound=tangentBound(arc);
        for(std::size_t part=1;part<cuts.size();++part){const double lo=cuts[part-1],hi=cuts[part];const auto middle=arc.sampleFraction(.5*(lo+hi));
            const int slot=grid.locate(middle);require(slot>=0,"parameter surface piece left the native domain");
            const int cell=cells[slot];const auto box=grid.rectangle(slot);
            const double hx=box.right-box.left,hy=box.top-box.bottom;
            std::array<double,32> scale;scale.fill(bound/(h*h));
            const auto values=surfaceIntegral<32>([&](double t){const auto q=P::sampleSurface(arc,t);const double x=(q.point.x-box.left)/hx,y=(q.point.y-box.bottom)/hy;
                const double enclosure=128*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(box.left)/hx,std::abs(box.right)/hx,std::abs(box.bottom)/hy,std::abs(box.top)/hy});
                require(x>=-enclosure && x<=1+enclosure && y>=-enclosure && y<=1+enclosure,"parameter surface quadrature crossed its cell partition");
                // Only terminal basis evaluation uses the floating coordinate
                // enclosure. The retained curve and partition are unchanged.
                const double variation=q.curvature-report.meanCurvature;
                return adjoint.moment(cell,{std::clamp(x,0.,1.),std::clamp(y,0.,1.)},
                    {-variation*q.derivative.y,variation*q.derivative.x});},lo,hi,scale,report);
            ++report.surfacePieces;adjoint.accumulate(cell,values);
        }
    }
    const auto variation=adjoint.integratedForce();
    for(const auto face:domain.faces(mesh)){const int row=face.nativeRow;
        auto& force=face.horizontal?result.faceForceY.value[row-mesh.vEdgeIndexNum()]:result.faceForceX.value[row];
        force+=checked(sigma*variation[row]/pressure.duals().volumes()[row].measure,"nonfinite parameter capillary variation");
    }
    result.forceX.initialize(mesh,0);result.forceY.initialize(mesh,0);result.maxForce=0;
    integrateNative(mesh,domain,pressure,result);if(output)*output=report;return result;
}

HydroUniformLoad buildUniformScalarGradientLoad(const QuadtreeMesh& mesh,
    const NativeMeshIdentity& source, const std::vector<double>& scalar,
    const std::vector<double>& coefficient)
{
    source.validate(mesh);
    const UniformLeafGrid grid(mesh, CartesianGridScope::RootOnly);
    require(scalar.size() == static_cast<std::size_t>(mesh.cellIndexNum())
        && coefficient.size() == static_cast<std::size_t>(mesh.vEdgeIndexNum() + mesh.hEdgeIndexNum()),
        "scalar/coefficient storage does not match the source mesh");
    for (double value : scalar) require(std::isfinite(value), "non-finite cell scalar");
    for (double value : coefficient) require(std::isfinite(value), "non-finite native coefficient");
    FaceSegmentGraph graph; graph.build(mesh);
    HydroUniformLoad result(mesh);
    for (const auto& face : graph.segments())
    {
        const int edge = face.isHorizontal ? mesh.hEdges(face.edgeLabel).index() : mesh.vEdges(face.edgeLabel).index();
        const int row = edge + (face.isHorizontal ? mesh.vEdgeIndexNum() : 0);
        if (face.boundary != Boundary::non) {
            require(coefficient[row] == 0.0, "prescribed boundary velocity cannot receive a scalar-gradient load");
            continue;
        }
        if (face.n.x + face.n.y < 0.0) continue;
        const int low = mesh.cells(face.K).index(), high = mesh.cells(face.L).index();
        const double force = coefficient[row] * ((scalar[high] - scalar[low]) / face.df);
        require(std::isfinite(force), "non-finite scalar-gradient force");
        if (face.isHorizontal) result.faceForceY.value[edge] = force;
        else result.faceForceX.value[edge] = force;
    }
    cellAverages(mesh, graph, false, result);
    // Preserve the original diffuse double-precision cell reduction order.
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const double fx = result.forceX.value[i], fy = result.forceY.value[i];
        require(std::isfinite(fx) && std::isfinite(fy), "non-finite cell force");
        const double area = mesh.cells(mesh.cellIndexToLabel[i]).area();
        result.maxForce = std::max(result.maxForce, std::hypot(fx, fy));
        result.totalForceX += fx * area;
        result.totalForceY += fy * area;
    }
    require(std::isfinite(result.maxForce) && std::isfinite(result.totalForceX)
        && std::isfinite(result.totalForceY), "integrated scalar-gradient load overflow");
    return result;
}

HydroUniformLoad buildUniformInterfaceTractionLoad(const QuadtreeMesh& mesh,
    const VertField<double>& scalar, const HydroInterfaceLoadGeometry& geometry,
    const std::shared_ptr<const InterfaceTraceGeometry>& queries,
    const std::vector<Vector2D>& traction, bool includeTangential)
{
    geometry.validateSource(mesh, scalar);
    const UniformLeafGrid grid(mesh, CartesianGridScope::RootOnly);
    require(queries && queries == geometry.traceGeometry(), "traction query source differs from the receiving load geometry");
    require(traction.size() == geometry.faces().size(), "traction count does not match the complete receiving support");
    require(!includeTangential || geometry.options().includeTangentialLoads, "tangential traction requires transverse support");
    HydroUniformLoad result(mesh);
    for (std::size_t i = 0; i < traction.size(); ++i)
    {
        const auto& face = geometry.faces()[i];
        const double tau = traction[i].x, tangent = traction[i].y;
        require(std::isfinite(tau) && std::isfinite(tangent), "non-finite interface traction");
        require(includeTangential || tangent == 0.0, "normal-only loading cannot discard tangential traction");
        double force = 0.0;
        if (includeTangential)
        {
            const double normalLoad = checked(static_cast<long double>(tau) * face.gradient, "non-finite normal face load");
            const double tangentLoad = checked(static_cast<long double>(tangent) * face.transverseGradient, "non-finite transverse face load");
            require(tau == 0.0 || face.gradient == 0.0 || normalLoad != 0.0, "nonzero normal face load underflow");
            require(tangent == 0.0 || face.transverseGradient == 0.0 || tangentLoad != 0.0, "nonzero transverse face load underflow");
            force = checked(static_cast<long double>(normalLoad) + (face.horizontal ? tangentLoad : -tangentLoad),
                "non-finite full interface face load");
        }
        else
        {
            force = checked(static_cast<long double>(tau) * face.gradient, "non-finite normal face load");
            require(tau == 0.0 || face.gradient == 0.0 || force != 0.0, "nonzero normal face load underflow");
        }
        const double integral = checked(static_cast<long double>(force) * face.dualMeasure, "dual-face integral overflow");
        require(force == 0.0 || integral != 0.0, "nonzero dual-face integral underflow");
        if (face.horizontal) result.faceForceY.value[face.edgeIndex] = force;
        else result.faceForceX.value[face.edgeIndex] = force;
    }
    FaceSegmentGraph graph; graph.build(mesh);
    const bool axisymmetric = geometry.options().axisymmetric;
    integrate(mesh,graph,axisymmetric,result);
    return result;
}
static HydroUniformLoad conservativeCapillaryLoad(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceState& state,double surfaceTension,
    double projectionHalfWidth,ConservativeCapillaryReport* report,const HydroConservativePressure* pressure,
    std::shared_ptr<const ConservativeTractionGeometry>* tractionGeometry=nullptr,
    const std::vector<long double>* componentFractions=nullptr,std::uint64_t liquidId=0,
    const InterfaceCurvature* sharedCurvature=nullptr)
{
    require(std::isfinite(surfaceTension) && surfaceTension>=0,"invalid surface tension");
    require(std::isfinite(projectionHalfWidth) && projectionHalfWidth>0,"invalid capillary projection width");
    domain.validateMeshIdentity(mesh);
    require(domain.cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "uniform conservative capillary load requires the complete fluid domain");
    if(pressure)pressure->validate(mesh,domain);
    else { const UniformLeafGrid grid(mesh,CartesianGridScope::RootOnly); }
    const auto& phase=state.conservativePhase();
    const bool splitComponents=state.topology().liquids().size()>1
        && std::any_of(state.topology().events().begin(),state.topology().events().end(),
            [](const auto& event){return !event.interiorNecks.empty();});
    if(splitComponents && !componentFractions) {
        // C is unchanged. Child inventories follow accepted directional phase
        // fluxes, so gas-gap tails retain their material provenance.
        // Linearity preserves G(C)=sum G(C_k), including a pressure face where
        // two different caps contribute and the total C gradient cancels.
        const auto& components=state.componentInventories();
        require(components.size()==state.topology().liquids().size(),"split force lacks transported child inventories");
        for(const auto& component:components)require(component.liquidVolumes.size()==phase.liquidVolumes().size()
            && std::any_of(state.topology().liquids().begin(),state.topology().liquids().end(),
                [&](const auto& liquid){return liquid.id==component.liquidId;}),"split force lost a child inventory identity/layout");
        for(int c:domain.cellIndexByRow()) {
            long double sum=0;
            for(const auto& component:components) {
                require(std::isfinite(component.liquidVolumes[c]) && component.liquidVolumes[c]>=0,
                    "split force has an invalid material inventory");sum+=component.liquidVolumes[c];
            }
            require(std::abs(sum-phase.liquidVolumes()[c])<=4*std::numeric_limits<long double>::epsilon()*phase.liquidVolumes()[c],
                "split force material inventories do not partition C");
        }
        HydroUniformLoad combined(mesh);ConservativeCapillaryReport quality;
        auto combinedTraction=tractionGeometry?std::make_shared<ConservativeTractionGeometry>():nullptr;
        // Every child uses the same immutable geometry/C and curvature knots.
        // Prepare them once for this load only; scoped projection still uses
        // each child's actual material inventory and liquid identity.
        std::optional<InterfaceCurvature> curvature;
        for(const auto& component:components) {
            std::vector<long double> fractions(mesh.cellIndexNum(),0);bool present=false;
            for(int c:domain.cellIndexByRow()) {
                fractions[c]=component.liquidVolumes[c]/phase.capacities()[c];present=present || fractions[c]!=0;
            }
            if(!present)continue;
            if(!curvature)curvature.emplace(mesh,domain,state);
            ConservativeCapillaryReport partial;
            std::shared_ptr<const ConservativeTractionGeometry> childTraction;
            const auto load=conservativeCapillaryLoad(mesh,domain,state,surfaceTension,projectionHalfWidth,
                &partial,pressure,tractionGeometry?&childTraction:nullptr,&fractions,component.liquidId,&*curvature);
            if(combinedTraction)combinedTraction->components.push_back(std::move(childTraction));
            for(int i=0;i<combined.faceForceX.size();++i)combined.faceForceX.value[i]+=load.faceForceX.value[i];
            for(int i=0;i<combined.faceForceY.size();++i)combined.faceForceY.value[i]+=load.faceForceY.value[i];
            quality.heightFaces+=partial.heightFaces;quality.localFitFaces+=partial.localFitFaces;
            quality.mixedFaces+=partial.mixedFaces;quality.eventModelFaces+=partial.eventModelFaces;
            quality.unresolvedTailFaces+=partial.unresolvedTailFaces;quality.boundedTailRidgeFaces+=partial.boundedTailRidgeFaces;
            quality.maxEventModelForce=std::max(quality.maxEventModelForce,partial.maxEventModelForce);
            quality.maxEventProjectionForceSpread=std::max(quality.maxEventProjectionForceSpread,partial.maxEventProjectionForceSpread);
            quality.maxUnresolvedTailForce=std::max(quality.maxUnresolvedTailForce,partial.maxUnresolvedTailForce);
            quality.maxProjectionDistance=std::max(quality.maxProjectionDistance,partial.maxProjectionDistance);
            quality.maxScaledTailForceUncertainty=std::max(quality.maxScaledTailForceUncertainty,partial.maxScaledTailForceUncertainty);
            quality.maxTailForceUncertainty=std::max(quality.maxTailForceUncertainty,partial.maxTailForceUncertainty);
            quality.surfaceIntegratedFaces+=partial.surfaceIntegratedFaces;
            quality.maxNormalMomentDefect=std::max(quality.maxNormalMomentDefect,partial.maxNormalMomentDefect);
            if(report)for(auto& face:partial.faces)quality.faces.push_back(std::move(face));
        }
        if(pressure)integrateNative(mesh,domain,*pressure,combined);
        else {FaceSegmentGraph graph;graph.build(mesh);integrate(mesh,graph,domain.axisymmetric(),combined);}
        if(combinedTraction) {
            combinedTraction->meshIdentity=domain.meshIdentity();combinedTraction->nativePressure=pressure!=nullptr;
            combinedTraction->liquidVolumes=phase.liquidVolumes();*tractionGeometry=std::move(combinedTraction);
        }
        if(report)*report=std::move(quality);return combined;
    }
    std::optional<InterfaceCurvature> localCurvature;
    if(!sharedCurvature)localCurvature.emplace(mesh,domain,state);
    const auto& curvature=sharedCurvature?*sharedCurvature:*localCurvature;
    const auto fraction=[&](int cell) {return componentFractions?componentFractions->at(cell):phase.liquidVolumes()[cell]/phase.capacities()[cell];};
    std::vector<double> nativeGradient;
    if(pressure) {
        std::vector<double> values(mesh.cellIndexNum()),boundary(domain.nativeFaceCount());
        for(int cell:domain.cellIndexByRow())values[cell]=static_cast<double>(fraction(cell));
        for(const auto face:domain.faces(mesh))if(face.negativeCell<0 || face.positiveCell<0)
            boundary[face.nativeRow]=values[std::max(face.negativeCell,face.positiveCell)];
        nativeGradient=pressure->apply(mesh,domain,values,boundary);
    }
    std::vector<double> transverse(tractionGeometry?domain.nativeFaceCount():0);
    std::vector<std::vector<int>> transverseRows(tractionGeometry?domain.nativeFaceCount():0);
    if(tractionGeometry) {
        // Fixed incidence stencil of the two actual owners, including their
        // subfaces. Affine interpolation uses face gradients per dual volume.
        // Prescribed axis/reflection boundary velocities receive zero loads.
        std::vector<double> gradient(domain.nativeFaceCount());
        for(const auto f:domain.faces(mesh)) if(f.negativeCell>=0 && f.positiveCell>=0)
            gradient[f.nativeRow]=pressure?nativeGradient[f.nativeRow]/pressure->duals().volumes()[f.nativeRow].measure
                :static_cast<double>((fraction(f.positiveCell)-fraction(f.negativeCell))/f.normalSpacing);
        for(const auto f:domain.faces(mesh)) {
            if(f.negativeCell<0 || f.positiveCell<0)continue;
            auto& rows=transverseRows[f.nativeRow];Vector2D h{0,0};
            for(int owner:{f.negativeCell,f.positiveCell}) {
                const auto& cell=mesh.cells(mesh.cellIndexToLabel[owner]);
                h.x=std::max(h.x,cell.size().x);h.y=std::max(h.y,cell.size().y);
                for(int side=0;side<2;++side)for(int part=0;part<2;++part) {
                    const int label=f.horizontal?cell.cellToVEdge[side][part]:cell.cellToHEdge[side][part];
                    if(label>=0)rows.push_back(f.horizontal?mesh.vEdges(label).index():mesh.vEdgeIndexNum()+mesh.hEdges(label).index());
                }
            }
            std::sort(rows.begin(),rows.end());rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
            Eigen::Matrix3d gram=Eigen::Matrix3d::Zero();std::vector<Eigen::Vector3d> basis;
            for(int row:rows) {
                const auto p=domain.face(mesh,row).center;
                basis.emplace_back(1,(p.x-f.center.x)/h.x,(p.y-f.center.y)/h.y);
                gram.noalias()+=basis.back()*basis.back().transpose();
            }
            Eigen::LDLT<Eigen::Matrix3d> fit(gram);
            require(fit.info()==Eigen::Success && fit.isPositive() && fit.vectorD().minCoeff()>1e-12,
                "conservative transverse stencil is rank deficient");
            const Eigen::Vector3d coefficient=fit.solve(Eigen::Vector3d(1,0,0));
            Eigen::Vector3d reproduction=Eigen::Vector3d::Zero();long double value=0;
            for(std::size_t i=0;i<rows.size();++i) {
                const double weight=basis[i].dot(coefficient);reproduction+=weight*basis[i];
                value+=static_cast<long double>(weight)*gradient[rows[i]];
            }
            require((reproduction-Eigen::Vector3d(1,0,0)).cwiseAbs().maxCoeff()<1e-10,
                "conservative transverse stencil lost affine reproduction");
            transverse[f.nativeRow]=checked(value,"conservative transverse gradient overflow");
        }
    }
    // Local event-model support is determined before any projection succeeds
    // or fails. Include every arc of this liquid in the two owner cells, not
    // just the nearest graph. Non-event owner neighbours are included too.
    std::vector<std::vector<std::pair<std::size_t,int>>> graphsByCell(mesh.cellIndexNum());
    std::vector<unsigned char> eventCell(mesh.cellIndexNum(),0),integratedRow(domain.nativeFaceCount(),0);
    if(!pressure && !state.eventCurvatureCells().empty()) {
        for(int c:state.eventCurvatureCells())eventCell.at(c)=1;
        const auto& surfaces=state.topology().surfaces();
        for(std::size_t s=0;s<surfaces.size();++s)if(!liquidId || surfaces[s].liquidId==liquidId) {
            const auto& graphs=surfaces[s].geometry.nativeChain().graphs();
            for(std::size_t g=0;g<graphs.size();++g)graphsByCell.at(graphs[g].nativeCellIndex).push_back({s,static_cast<int>(g)});
        }
        for(const auto f:domain.faces(mesh))if(f.negativeCell>=0 && f.positiveCell>=0)
            integratedRow[f.nativeRow]=(eventCell[f.negativeCell] && !graphsByCell[f.negativeCell].empty())
                || (eventCell[f.positiveCell] && !graphsByCell[f.positiveCell].empty());
    }
    struct SurfaceNode {const LevelSetInterfacePanels::QuadratureNode* geometry;InterfaceCurvature::Sample curvature;};
    std::vector<SurfaceNode> surfaceNodes;
    std::vector<std::vector<std::size_t>> nodesByCell(mesh.cellIndexNum());
    std::vector<unsigned char> preparedCell(mesh.cellIndexNum(),0);
    const auto prepareSurfaceCell=[&](int cell) {
        if(preparedCell[cell])return;preparedCell[cell]=1;
        for(const auto [s,g]:graphsByCell[cell]) {
            const auto& panels=state.topology().surfaces()[s].geometry.panels();
            require(panels.at(g).cellIndex==cell,"surface moment lost its native panel identity");
            for(const auto& q:panels[g].quadratureNodes) {
                auto k=curvature.evaluate(s,g,q.nativeGraphFraction);
                require(k.valid(),"surface moment has unavailable branch curvature");
                require(q.frame.hasCellChartPoint && q.frame.nativeCellIndex==cell && q.physicalWeight>0,
                    "surface moment needs retained positive native quadrature");
                nodesByCell[cell].push_back(surfaceNodes.size());surfaceNodes.push_back({&q,std::move(k)});
            }
        }
    };
    // Projection gets expensive when nonzero unresolved inventory tails
    // spread to more pressure rows. Queries share only immutable geometry,
    // curvature knots and C. Precompute each needed row independently, then
    // retain the original serial force/report/diagnostic assembly order.
    struct PreparedProjection {
        int nativeRow;
        std::vector<int> pressureSupport;
        std::optional<InterfaceCurvature::Extension> extension;
        std::exception_ptr failure;
        bool boundedTailRidge=false;
        double resolutionScale=0;
        long double scaledTailForceUncertainty=0;
    };
    std::vector<PreparedProjection> projections;
    std::vector<int> projectionByRow(domain.nativeFaceCount(),-1);
    if(surfaceTension!=0 || tractionGeometry) for(const auto face:domain.faces(mesh)) {
        if(face.negativeCell<0 || face.positiveCell<0) continue;
        if(integratedRow[face.nativeRow])continue;
        const long double difference=pressure?static_cast<long double>(nativeGradient[face.nativeRow])
            :fraction(face.positiveCell)-fraction(face.negativeCell);
        if(difference==0 && (!tractionGeometry || transverse[face.nativeRow]==0)) continue;
        projectionByRow[face.nativeRow]=static_cast<int>(projections.size());
        projections.push_back({face.nativeRow});
    }
    const auto prepareProjection=[&](std::size_t index) {
        auto& prepared=projections[index];
        try {
            const auto face=domain.face(mesh,prepared.nativeRow);
            auto& support=prepared.pressureSupport;
            support={face.negativeCell,face.positiveCell};
            if(pressure || tractionGeometry) {
                if(pressure)for(HydroConservativePressure::Matrix::InnerIterator term(pressure->gradientLoad(),face.nativeRow);term;++term)
                    support.push_back(static_cast<int>(term.col()));
                if(tractionGeometry)for(int row:transverseRows[face.nativeRow]) {
                    if(pressure)for(HydroConservativePressure::Matrix::InnerIterator term(pressure->gradientLoad(),row);term;++term)
                        support.push_back(static_cast<int>(term.col()));
                    else {
                        const auto f=domain.face(mesh,row);
                        if(f.negativeCell>=0)support.push_back(f.negativeCell);
                        if(f.positiveCell>=0)support.push_back(f.positiveCell);
                    }
                }
                std::sort(support.begin(),support.end());
                support.erase(std::unique(support.begin(),support.end()),support.end());
            }
            prepared.extension=curvature.extendConservativeFace(face.nativeRow,projectionHalfWidth,support,liquidId);
            const auto& extension=*prepared.extension;
            if(extension.status==InterfaceCurvature::ExtensionStatus::AmbiguousBranch
                && extension.unresolvedInventoryTail && !extension.tailAlternatives.empty()) {
                double h=std::numeric_limits<double>::infinity();
                for(int cell:support) {
                    const auto width=mesh.cells(mesh.cellIndexToLabel[cell]).size();
                    h=std::min({h,width.x,width.y});
                }
                const long double difference=pressure?static_cast<long double>(nativeGradient[face.nativeRow]):
                    fraction(face.positiveCell)-fraction(face.negativeCell);
                const double spacing=pressure?pressure->duals().volumes()[face.nativeRow].measure:face.normalSpacing;
                prepared.resolutionScale=h;
                prepared.scaledTailForceUncertainty=std::abs(difference/spacing)*extension.tailCurvatureSpread*h*h;
                prepared.boundedTailRidge=std::isfinite(prepared.scaledTailForceUncertainty)
                    && prepared.scaledTailForceUncertainty<=Quantum::Geometry::PlanarCellCut::relativeVolumeResolution;
            }
        } catch(...) {prepared.failure=std::current_exception();}
    };
    // First actual load query stays on the caller for the strict geometry/FP
    // contract. Even its exception is consumed in the original face order,
    // after any earlier boundary checks. Never throw through an OpenMP team.
    if(!projections.empty()) prepareProjection(0);
#pragma omp parallel for schedule(dynamic, 4) if(projections.size()>16)
    for(std::ptrdiff_t index=1;index<static_cast<std::ptrdiff_t>(projections.size());++index)
        prepareProjection(static_cast<std::size_t>(index));
    ConservativeCapillaryReport quality;HydroUniformLoad result(mesh);
    auto traction=tractionGeometry?std::make_shared<ConservativeTractionGeometry>():nullptr;
    std::vector<InterfaceTraceGeometry::Point> queries;
    std::vector<InterfaceTraceGeometry::Point> alternateQueries;
    for (const auto face:domain.faces(mesh)) {
        if (face.negativeCell<0 || face.positiveCell<0) {
            const int cell=std::max(face.negativeCell,face.positiveCell);
            const auto& sections=state.topology().symmetryPlanes();
            const bool symmetry=std::any_of(sections.begin(),sections.end(),[&](const auto& section) {
                return section.normalAxis==(face.horizontal?1:0)
                    && section.coordinate==(face.horizontal?face.center.y:face.center.x)
                    && (!domain.axisymmetric() || section.normalAxis==1);
            });
            require((domain.axisymmetric() && face.center.x==0)
                || symmetry
                || (fraction(cell)<=Quantum::Geometry::PlanarCellCut::relativeVolumeResolution
                    && state.topology().cellPhase()[cell]==1),
                "conservative capillary load requires separation from boundaries without axis/section symmetry");
            continue;
        }
        const long double difference=pressure?static_cast<long double>(nativeGradient[face.nativeRow])
            :fraction(face.positiveCell)-fraction(face.negativeCell);
        const double spacing=pressure?pressure->duals().volumes()[face.nativeRow].measure:face.normalSpacing;
        if(integratedRow[face.nativeRow] && (surfaceTension!=0 || tractionGeometry)) {
            std::vector<std::size_t> indices;std::vector<Detail::SurfaceMomentNode> nodes;
            for(int owner:{face.negativeCell,face.positiveCell}) {
                prepareSurfaceCell(owner);
                const auto bounds=state.geometry()->nativeCellBounds(owner);
                for(std::size_t id:nodesByCell[owner]) {
                    const auto& q=*surfaceNodes[id].geometry;
                    const auto basis=Detail::surfaceFaceBasis(q.frame.cellChartPoint,bounds,phase.capacities()[owner],
                        spacing,face.horizontal,owner==face.negativeCell,domain.axisymmetric());
                    indices.push_back(id);nodes.push_back({basis*q.physicalWeight,
                        face.horizontal?q.frame.normal.y:q.frame.normal.x,
                        face.horizontal?q.frame.normal.x:-q.frame.normal.y});
                }
            }
            const long double gradient=difference/spacing;
            const auto balanced=Detail::balanceSurfaceMoments(nodes,-gradient,
                tractionGeometry?(face.horizontal?-1.L:1.L)*transverse[face.nativeRow]:0);
            const auto& first=surfaceNodes.at(indices.front());
            const double reference=first.curvature.curvature;
            long double force=reference*gradient,mean=0,meridional=0,azimuthal=0,area=0;
            bool event=false;
            for(std::size_t j=0;j<indices.size();++j) {
                const auto& k=surfaceNodes[indices[j]].curvature;
                force-=balanced.normal[j]*(k.curvature-reference);
                area+=nodes[j].weight;mean+=nodes[j].weight*k.curvature;
                meridional+=nodes[j].weight*k.meridional;azimuthal+=nodes[j].weight*k.azimuthal;
                event=event || k.method==InterfaceCurvature::Method::EventModel;
            }
            const double value=checked(surfaceTension*force,"integrated capillary force overflow");
            if(traction) {
                ConservativeTractionGeometry::IntegratedFace receiving{{face.nativeRow,face.horizontal,
                    static_cast<double>(gradient),transverse[face.nativeRow]},{}};
                for(std::size_t j=0;j<indices.size();++j)
                    receiving.moments.push_back({indices[j],balanced.normal[j],balanced.tangent[j]});
                traction->integratedFaces.push_back(std::move(receiving));
            }
            ++quality.surfaceIntegratedFaces;
            quality.maxNormalMomentDefect=std::max(quality.maxNormalMomentDefect,
                checked(std::abs(balanced.normalDefect),"surface moment defect overflow"));
            if(event) {++quality.eventModelFaces;quality.maxEventModelForce=std::max(quality.maxEventModelForce,std::abs(value));}
            else ++quality.mixedFaces;
            if(report) {
                ConservativeCapillaryReport::Face entry;
                entry.nativeRow=face.nativeRow;entry.center=face.center;entry.interfacePoint=first.geometry->frame.point;
                entry.gradient=static_cast<double>(gradient);entry.curvature=checked(mean/area,"mean curvature overflow");
                entry.meridional=checked(meridional/area,"mean meridional curvature overflow");
                entry.azimuthal=checked(azimuthal/area,"mean azimuthal curvature overflow");entry.force=value;
                entry.primaryMethod=static_cast<int>(first.curvature.method);entry.liquidId=liquidId;
                entry.surfaceIntegrated=true;entry.surfaceSamples=indices.size();
                entry.normalMomentDefect=checked(balanced.normalDefect,"surface moment defect overflow");
                quality.faces.push_back(entry);
            }
            if(face.horizontal)result.faceForceY.value[mesh.hEdges(face.edgeLabel).index()]=value;
            else result.faceForceX.value[mesh.vEdges(face.edgeLabel).index()]=value;
            continue;
        }
        if(tractionGeometry) {if(difference==0 && transverse[face.nativeRow]==0)continue;}
        else if (difference==0 || surfaceTension==0) continue;
        const auto& prepared=projections.at(projectionByRow[face.nativeRow]);
        if(prepared.failure) std::rethrow_exception(prepared.failure);
        const auto& pressureSupport=prepared.pressureSupport;
        const auto& extension=*prepared.extension;
        if (!extension.valid() && !prepared.boundedTailRidge) {
            std::ostringstream detail;detail<<std::setprecision(21)
                <<"[HydroUniformLoad] conservative capillary face "<<face.nativeRow
                <<" has unavailable curvature: extension="<<static_cast<int>(extension.status)
                <<", height="<<static_cast<int>(extension.sample.height.status)
                <<", fit="<<static_cast<int>(extension.sample.fit.status)
                <<", secondary_height="<<static_cast<int>(extension.sample.secondary.height.status)
                <<", secondary_fit="<<static_cast<int>(extension.sample.secondary.fit.status)
                <<", point=("<<face.center.x<<','<<face.center.y<<"), owners="<<face.negativeCell<<'/'<<face.positiveCell
                <<", C="<<fraction(face.negativeCell)<<'/'<<fraction(face.positiveCell)
                <<", distance="<<extension.distance<<", half_width="<<projectionHalfWidth
                <<", unresolved_tail="<<extension.unresolvedInventoryTail
                <<", tail_alternatives="<<extension.tailAlternatives.size()
                <<", scaled_tail_force_uncertainty="<<prepared.scaledTailForceUncertainty;
            if(extension.status==InterfaceCurvature::ExtensionStatus::CurvatureUnavailable) {
                for(const auto* sample:{static_cast<const InterfaceCurvature::Source*>(&extension.sample),&extension.sample.secondary}) {
                    detail<<"; fit_at=("<<sample->fit.point.x<<','<<sample->fit.point.y<<')'
                        <<", arc_cells="<<sample->fit.maximumArcCells
                        <<", support="<<sample->fit.negativeSupport<<'/'<<sample->fit.positiveSupport
                        <<", residual="<<sample->fit.scaledResidual<<", offset="<<sample->fit.scaledPositionError
                        <<", rank="<<sample->fit.scaledRankRatio<<", refinement="<<static_cast<int>(sample->fit.refinementCause);
                    detail<<", orientation="<<static_cast<int>(sample->fit.orientationCause);
                    if(!sample->fit.rejectedShorterFit.empty())
                        detail<<"; rejected_shorter_fit={"<<sample->fit.rejectedShorterFit<<'}';
                }
            }
            if(pressure) {
                detail<<", integrated_C_gradient="<<difference<<", support=";
                for(int cell:pressureSupport)detail<<cell<<':'<<fraction(cell)<<';';
            }
            if (extension.status==InterfaceCurvature::ExtensionStatus::AmbiguousBranch) {
                for (bool other:{false,true}) {
                    const auto& candidate=other?extension.competingProjection:extension.projection;
                    const auto provider=curvature.evaluate(other?extension.competingSurface:extension.surface,
                        candidate.panelIndex,candidate.nativeGraphFraction);
                    detail<<(other?"; competitor":"; nearest")<<" surface="<<(other?extension.competingSurface:extension.surface)
                        <<", graph="<<candidate.panelIndex<<", fraction="<<candidate.nativeGraphFraction
                        <<", point=("<<candidate.frame.point.x<<','<<candidate.frame.point.y<<')'
                        <<", normal=("<<candidate.frame.normal.x<<','<<candidate.frame.normal.y<<')'
                        <<", raw_curvature="<<candidate.frame.curvature<<", distance="<<candidate.distance
                        <<", stationary_residual="<<candidate.kktResidual
                        <<", provider="<<static_cast<int>(provider.method)<<", provider_curvature="<<provider.curvature;
                }
                if(extension.competingSurface==extension.surface) {
                    const auto& chain=state.topology().surfaces()[extension.surface].geometry.nativeChain();
                    for(int g=std::min(extension.projection.panelIndex,extension.competingProjection.panelIndex);
                        g<=std::max(extension.projection.panelIndex,extension.competingProjection.panelIndex);++g) {
                        const auto& graph=chain.graphs()[g];const auto a=chain.frameAt(g,0),b=chain.frameAt(g,1);
                        detail<<"; connecting_graph="<<g<<", cell="<<graph.nativeCellIndex
                            <<", coordinate_x="<<graph.coordinateIsX
                            <<", chart=("<<graph.localLower.x<<','<<graph.localLower.y<<")-("
                            <<graph.localUpper.x<<','<<graph.localUpper.y<<")"
                            <<", endpoint_r="<<a.point.x<<'/'<<b.point.x
                            <<", endpoint_z="<<a.point.y<<'/'<<b.point.y
                            <<", endpoint_nr="<<a.normal.x<<'/'<<b.normal.x;
                    }
                }
            }
            throw std::invalid_argument(detail.str());
        }
        if(tractionGeometry) {
            const auto& frame=extension.projection.frame;
            queries.push_back({frame,state.geometry()->curvatureTangentialDerivativeAtNativeCellChart(
                frame.nativeCellIndex,frame.cellChartPoint),extension.projection.kktResidual});
            traction->faces.push_back({face.nativeRow,face.horizontal,static_cast<double>(difference/spacing),transverse[face.nativeRow]});
            if(prepared.boundedTailRidge)for(const auto& alternative:extension.tailAlternatives) {
                const auto& frame=alternative.projection.frame;
                traction->tailAlternatives.push_back({traction->faces.size()-1,alternateQueries.size(),prepared.resolutionScale});
                alternateQueries.push_back({frame,state.geometry()->curvatureTangentialDerivativeAtNativeCellChart(
                    frame.nativeCellIndex,frame.cellChartPoint),alternative.projection.kktResidual});
            }
            if(difference==0 || surfaceTension==0)continue;
        }
        if (extension.sample.method==InterfaceCurvature::Method::HeightFunction) ++quality.heightFaces;
        else if (extension.sample.method==InterfaceCurvature::Method::LocalFit) ++quality.localFitFaces;
        else if (extension.sample.method==InterfaceCurvature::Method::EventModel) ++quality.eventModelFaces;
        else ++quality.mixedFaces;
        quality.maxProjectionDistance=std::max(quality.maxProjectionDistance,extension.distance);
        const double force=checked(surfaceTension*static_cast<long double>(extension.sample.curvature)
            *difference/spacing,"conservative capillary force overflow");
        require(extension.sample.curvature==0 || force!=0,"conservative capillary force underflow");
        if(extension.sample.method==InterfaceCurvature::Method::EventModel)
            quality.maxEventModelForce=std::max(quality.maxEventModelForce,std::abs(force));
        quality.maxEventProjectionForceSpread=std::max(quality.maxEventProjectionForceSpread,
            checked(surfaceTension*std::abs(difference/spacing)*extension.eventCurvatureSpread,"event force sensitivity overflow"));
        if (extension.unresolvedInventoryTail) {
            ++quality.unresolvedTailFaces;quality.maxUnresolvedTailForce=std::max(quality.maxUnresolvedTailForce,std::abs(force));
        }
        if(prepared.boundedTailRidge) {
            ++quality.boundedTailRidgeFaces;
            quality.maxScaledTailForceUncertainty=std::max(quality.maxScaledTailForceUncertainty,
                static_cast<double>(prepared.scaledTailForceUncertainty));
            quality.maxTailForceUncertainty=std::max(quality.maxTailForceUncertainty,
                checked(surfaceTension*std::abs(difference/spacing)*extension.tailCurvatureSpread,"tail force sensitivity overflow"));
        }
        if(report) quality.faces.push_back({face.nativeRow,face.center,extension.sample.point,
            static_cast<double>(difference/spacing),extension.sample.curvature,
            extension.sample.meridional,extension.sample.azimuthal,force,
            extension.sample.primaryMethod==InterfaceCurvature::Method::HeightFunction,
            extension.sample.height.valid() && extension.sample.height.evenAxisFit,
            extension.sample.height.valid() && extension.sample.height.implicitAxisFit,
            extension.sample.interpolated,static_cast<int>(extension.sample.primaryMethod),
            static_cast<int>(extension.sample.secondary.method),extension.sample.secondaryWeight,
            extension.sample.secondary.height.valid() && extension.sample.secondary.height.evenAxisFit,
            extension.sample.secondary.height.valid() && extension.sample.secondary.height.implicitAxisFit,
            extension.unresolvedInventoryTail,extension.distance,extension.maxMinorityFraction,extension.sameChartTailCandidates,
            extension.sample.regularized,extension.sample.secondaryCurvatureWeight,
            static_cast<int>(curvature.regularizationStatus(extension.surface)),
            extension.eventChartCandidates,extension.eventCurvatureSpread,prepared.boundedTailRidge,
            extension.tailCurvatureSpread,static_cast<double>(prepared.scaledTailForceUncertainty),extension.tailAlternatives.size(),liquidId});
        if (face.horizontal) result.faceForceY.value[mesh.hEdges(face.edgeLabel).index()]=force;
        else result.faceForceX.value[mesh.vEdges(face.edgeLabel).index()]=force;
    }
    if(pressure)integrateNative(mesh,domain,*pressure,result);
    else { FaceSegmentGraph graph;graph.build(mesh);integrate(mesh,graph,domain.axisymmetric(),result); }
    if (report) *report=std::move(quality);
    if(tractionGeometry) {
        traction->meshIdentity=domain.meshIdentity();
        traction->nativePressure=pressure!=nullptr;
        traction->liquidVolumes=phase.liquidVolumes();
        traction->liquidId=liquidId;
        if(liquidId) {
            const auto found=std::find_if(state.componentInventories().begin(),state.componentInventories().end(),
                [&](const auto& component){return component.liquidId==liquidId;});
            require(found!=state.componentInventories().end(),"child traction geometry lost its material inventory");
            traction->liquidVolumes=found->liquidVolumes;
        }
        for(auto& alternative:traction->tailAlternatives)alternative.queryIndex+=queries.size();
        queries.insert(queries.end(),alternateQueries.begin(),alternateQueries.end());
        for(auto& face:traction->integratedFaces)for(auto& moment:face.moments)moment.queryIndex+=queries.size();
        traction->integrationQueryCount=surfaceNodes.size();
        for(const auto& node:surfaceNodes) {
            const auto& frame=node.geometry->frame;
            queries.push_back({frame,state.geometry()->curvatureTangentialDerivativeAtNativeCellChart(
                frame.nativeCellIndex,frame.cellChartPoint),0});
        }
        traction->trace=std::make_shared<const InterfaceTraceGeometry>(state.geometry(),std::move(queries));
        *tractionGeometry=std::move(traction);
    }
    return result;
}
HydroUniformLoad buildUniformConservativeCapillaryLoad(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceState& state,double surfaceTension,
    double projectionHalfWidth,ConservativeCapillaryReport* report,
    std::shared_ptr<const ConservativeTractionGeometry>* tractionGeometry)
{
    return conservativeCapillaryLoad(mesh,domain,state,surfaceTension,projectionHalfWidth,report,nullptr,tractionGeometry);
}
HydroUniformLoad buildNativeConservativeCapillaryLoad(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceState& state,const HydroConservativePressure& pressure,
    double surfaceTension,double projectionHalfWidth,ConservativeCapillaryReport* report,
    std::shared_ptr<const ConservativeTractionGeometry>* tractionGeometry)
{
    return conservativeCapillaryLoad(mesh,domain,state,surfaceTension,projectionHalfWidth,report,&pressure,tractionGeometry);
}
HydroUniformLoad buildConservativeInterfaceTractionLoad(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceState& state,const HydroConservativePressure* pressure,
    const ConservativeTractionGeometry& geometry,const std::shared_ptr<const InterfaceTraceGeometry>& queries,
    const std::vector<Vector2D>& traction,ConservativeTractionReport* report)
{
    domain.validateMeshIdentity(mesh);state.validateMesh(mesh);
    if(pressure)pressure->validate(mesh,domain);else {const UniformLeafGrid grid(mesh,CartesianGridScope::RootOnly);}
    require(geometry.nativePressure==(pressure!=nullptr),"conservative traction uses a different pressure discretization");
    const std::vector<long double>* receivingInventory=&state.conservativePhase().liquidVolumes();
    if(geometry.liquidId) {
        const auto found=std::find_if(state.componentInventories().begin(),state.componentInventories().end(),
            [&](const auto& component){return component.liquidId==geometry.liquidId;});
        require(found!=state.componentInventories().end(),"traction query has no matching child inventory");
        receivingInventory=&found->liquidVolumes;
    }
    require(geometry.components.empty() && geometry.meshIdentity==domain.meshIdentity() && geometry.trace && queries==geometry.trace
        && queries->source()==state.geometry() && geometry.liquidVolumes==*receivingInventory,
        "conservative traction source/queries/inventory mismatch");
    require(traction.size()==geometry.faces.size()+geometry.tailAlternatives.size()+geometry.integrationQueryCount
        && traction.size()==queries->points().size(),"incomplete conservative tractions");
    ConservativeTractionReport quality;
    const auto forceValue=[](const ConservativeTractionGeometry::Face& f,Vector2D t) {
        require(std::isfinite(t.x) && std::isfinite(t.y),"nonfinite conservative traction");
        return -static_cast<long double>(t.x)*f.gradient+(f.horizontal?-1.L:1.L)*t.y*f.transverse;
    };
    for(std::size_t i=0;i<geometry.tailAlternatives.size();++i) {
        const auto& alternative=geometry.tailAlternatives[i];
        require(alternative.faceIndex<geometry.faces.size()
            && alternative.queryIndex==geometry.faces.size()+i
            && std::isfinite(alternative.resolutionScale) && alternative.resolutionScale>0,
            "incomplete or reordered tail traction witnesses");
        const auto& f=geometry.faces[alternative.faceIndex];
        const auto a=traction[alternative.faceIndex],b=traction[alternative.queryIndex];
        const auto spread=std::abs(forceValue(f,a)-forceValue(f,b));
        const long double scale=std::max({std::abs(a.x),std::abs(a.y),std::abs(b.x),std::abs(b.y)});
        const long double scaled=scale>0?spread*alternative.resolutionScale/scale:0;
        require(std::isfinite(scaled) && scaled<=Quantum::Geometry::PlanarCellCut::relativeVolumeResolution,
            "unresolved inventory tail exceeds the Maxwell traction branch-error budget");
        ++quality.tailAlternatives;
        quality.maxScaledTailForceUncertainty=std::max(quality.maxScaledTailForceUncertainty,static_cast<double>(scaled));
        quality.maxTailForceUncertainty=std::max(quality.maxTailForceUncertainty,checked(spread,"tail Maxwell sensitivity overflow"));
    }
    HydroUniformLoad result(mesh);
    std::vector<unsigned char> usedRow(domain.nativeFaceCount(),0);
    int previousRow=-1;
    for(std::size_t i=0;i<geometry.faces.size();++i) {
        const auto& f=geometry.faces[i];const auto t=traction[i];
        require(f.row>previousRow && f.row<domain.nativeFaceCount(),"unordered/invalid conservative receiving row");
        const auto native=domain.face(mesh,f.row);previousRow=f.row;
        usedRow[f.row]=1;
        require(native.negativeCell>=0 && native.positiveCell>=0 && native.horizontal==f.horizontal
            && std::isfinite(f.gradient) && std::isfinite(f.transverse),"invalid conservative receiving face/gradient");
        require(std::isfinite(t.x) && std::isfinite(t.y),"nonfinite conservative traction");
        // n points from liquid to gas; grad(C)=-n delta, t=(-ny,nx).
        const double value=checked(-static_cast<long double>(t.x)*f.gradient
            +(f.horizontal?-1.L:1.L)*t.y*f.transverse,"conservative traction overflow");
        if(f.horizontal)result.faceForceY.value[f.row-mesh.vEdgeIndexNum()]=value;
        else result.faceForceX.value[f.row]=value;
    }
    previousRow=-1;
    const auto integrationBegin=geometry.faces.size()+geometry.tailAlternatives.size();
    std::vector<unsigned char> usedQuery(geometry.integrationQueryCount,0);
    for(const auto& integrated:geometry.integratedFaces) {
        const auto& f=integrated.face;
        require(f.row>previousRow && f.row<domain.nativeFaceCount() && !usedRow[f.row]
            && !integrated.moments.empty(),"invalid/duplicate integrated receiving row");
        previousRow=f.row;usedRow[f.row]=1;
        const auto native=domain.face(mesh,f.row);
        require(!pressure && native.negativeCell>=0 && native.positiveCell>=0 && native.horizontal==f.horizontal,
            "integrated receiving face is not a uniform interior face");
        const auto reference=traction.at(integrated.moments.front().queryIndex);
        long double force=forceValue(f,reference);
        for(const auto& m:integrated.moments) {
            require(m.queryIndex>=integrationBegin && m.queryIndex<traction.size()
                && std::isfinite(m.normal) && std::isfinite(m.tangent),"invalid surface moment query/coefficient");
            usedQuery[m.queryIndex-integrationBegin]=1;
            const auto t=traction[m.queryIndex];
            require(std::isfinite(t.x) && std::isfinite(t.y),"nonfinite integrated traction");
            force+=m.normal*(t.x-reference.x)+m.tangent*(t.y-reference.y);
        }
        const double value=checked(force,"integrated interface traction overflow");
        if(f.horizontal)result.faceForceY.value[f.row-mesh.vEdgeIndexNum()]=value;
        else result.faceForceX.value[f.row]=value;
    }
    require(std::all_of(usedQuery.begin(),usedQuery.end(),[](auto used){return used!=0;}),"unused integrated traction query");
    if(pressure)integrateNative(mesh,domain,*pressure,result);
    else {FaceSegmentGraph graph;graph.build(mesh);integrate(mesh,graph,domain.axisymmetric(),result);}
    if(report)*report=quality;
    return result;
}
}
