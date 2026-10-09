#include "Modules/Hydrodynamics/HydroConservativeStress.h"
#include "Modules/Hydrodynamics/HydroConservativePressure.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

namespace Quantum::Hydrodynamics
{
namespace {
using Trace=HydroFittedDetail::MomentumTrace;
using Type=HydroFittedDomain::BoundaryType;
void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[HydroConservativeStress] ")+why);}
double at(Vector2D p,int axis) {return axis?p.y:p.x;}
void add(Trace& target,const Trace& source,double scale)
{
    for(const auto& term:source) {
        const auto it=std::find_if(target.begin(),target.end(),[&](const auto& old){return old.first==term.first;});
        if(it==target.end())target.emplace_back(term.first,scale*term.second);else it->second+=scale*term.second;
    }
}
}
Eigen::SparseMatrix<double> HydroConservativeStress::assembleParametricEnergy(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceState& state,double liquidMu,double gasMu,EnergyReport* output)
{
    InterfaceParametricMotion::validateDomain(mesh,domain,state);
    require(!domain.axisymmetric() && domain.uniformLeaves(),"energy viscosity requires uniform XY parameter geometry");
    require(std::isfinite(liquidMu) && liquidMu>=0 && std::isfinite(gasMu) && gasMu>=0,"invalid phase viscosity");
    Quantum::Geometry::Parametric::PhaseMeasure phase(state.parametricGeometry().geometry());EnergyReport report;
    const auto material=[&](double x0,double x1,double y0,double y1){const double volume=(x1-x0)*(y1-y0);
        const double liquid=phase.simpleRectangle(x0,x1,y0,y1);
        require(std::isfinite(liquid) && liquid>=-2e-12*volume && liquid<=volume*(1+2e-12),"phase viscosity measure exceeds geometry budget");
        const double bounded=std::clamp(liquid,0.,volume);
        report.maximumMeasureFractionAdjustment=std::max(report.maximumMeasureFractionAdjustment,std::abs(bounded-liquid)/volume);
        // Terminal positive quadrature measure only; never writes FV inventory.
        return liquidMu*bounded+gasMu*(volume-bounded);};
    auto matrix=assembleEnergy(mesh,domain,material);
    if(output)*output=report;return matrix;
}
Eigen::SparseMatrix<double> HydroConservativeStress::assembleCellEnergy(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const std::vector<double>& viscosity)
{
    domain.validateMeshIdentity(mesh);
    require(domain.uniformLeaves() && viscosity.size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "cell energy viscosity requires uniform native cell coefficients");
    for(double value:viscosity)require(std::isfinite(value) && value>=0,"invalid cell energy viscosity");
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    const auto h=grid.spacing(),origin=grid.minCoord();
    const auto material=[&](double x0,double x1,double y0,double y1){
        const int i0=std::max(0,static_cast<int>(std::floor((x0-origin.x)/h.x))),
            i1=std::min(grid.nx()-1,static_cast<int>(std::floor((x1-origin.x)/h.x)));
        const int j0=std::max(0,static_cast<int>(std::floor((y0-origin.y)/h.y))),
            j1=std::min(grid.ny()-1,static_cast<int>(std::floor((y1-origin.y)/h.y)));
        long double sum=0;
        for(int j=j0;j<=j1;++j)for(int i=i0;i<=i1;++i){
            const double a=std::max(x0,origin.x+i*h.x),b=std::min(x1,origin.x+(i+1)*h.x),
                c=std::max(y0,origin.y+j*h.y),d=std::min(y1,origin.y+(j+1)*h.y);
            if(b>a && d>c)sum+=viscosity[mesh.cells(grid.cells()[j*grid.nx()+i]).index()]
                *static_cast<long double>(CartesianMeasure::cell((b-a)*(d-c),.5*(a+b),domain.axisymmetric()));
        }
        return static_cast<double>(sum);};
    return assembleEnergy(mesh,domain,material);
}
Eigen::SparseMatrix<double> HydroConservativeStress::assembleEnergy(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const ViscosityMeasure& material)
{
    domain.validateMeshIdentity(mesh);
    require(domain.uniformLeaves(),"energy viscosity requires uniform leaves");
    require(static_cast<bool>(material),"missing integrated viscosity measure");
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);const int nx=grid.nx(),ny=grid.ny();
    const auto h=grid.spacing(),origin=grid.minCoord();const int radial=mesh.vEdgeIndexNum();
    const auto xrow=[&](int i,int j){return mesh.vEdges(grid.vertical()[j*(nx+1)+i]).index();};
    const auto yrow=[&](int i,int j){return radial+mesh.hEdges(grid.horizontal()[j*nx+i]).index();};
    std::vector<Eigen::Triplet<double>> entries;
    const auto addEnergy=[&](const Trace& strain,double weight){require(weight>=0 && std::isfinite(weight),"negative/nonfinite strain energy weight");
        if(weight==0)return;for(const auto& a:strain)for(const auto& b:strain)entries.emplace_back(a.first,b.first,weight*a.second*b.second);};
    for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){const double x=origin.x+i*h.x,y=origin.y+j*h.y;
        const double weight=2*material(x,x+h.x,y,y+h.y);
        addEnergy({{xrow(i,j),-1/h.x},{xrow(i+1,j),1/h.x}},weight);
        addEnergy({{yrow(i,j),-1/h.y},{yrow(i,j+1),1/h.y}},weight);
        if(domain.axisymmetric()) {
            // Cell-centred u_r/r. The first cell uses u_r(0)=0 and a
            // positive centre radius: its regular affine limit is exact.
            // This strain participates in the same energy and transpose,
            // rather than being appended as a separate radial force.
            const double inverseDiameter=1/(2*x+h.x);
            addEnergy({{xrow(i,j),inverseDiameter},{xrow(i+1,j),inverseDiameter}},weight);
        }}
    for(int j=1;j<ny;++j)for(int i=1;i<nx;++i){const double x=origin.x+i*h.x,y=origin.y+j*h.y;
        const double weight=material(x-.5*h.x,x+.5*h.x,y-.5*h.y,y+.5*h.y);
        addEnergy({{xrow(i,j-1),-1/h.y},{xrow(i,j),1/h.y},{yrow(i-1,j),-1/h.x},{yrow(i,j),1/h.x}},weight);}
    Eigen::SparseMatrix<double> matrix(domain.nativeFaceCount(),domain.nativeFaceCount());matrix.setFromTriplets(entries.begin(),entries.end());
    return matrix;
}
HydroConservativeStress HydroConservativeStress::build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroConservativePressure& pressure,const std::vector<HydroFittedDomain::Boundary>& boundary)
{
    pressure.validate(mesh,domain);
    require(boundary.size()==static_cast<std::size_t>(domain.nativeFaceCount()),"boundary map size mismatch");
    for(const auto face:domain.faces(mesh)) {
        const auto type=boundary[face.nativeRow].type;
        if(face.outwardSign==0)require(type==Type::Interior,"boundary type on an interior face");
        else {
            require(type==Type::NoSlipWall || type==Type::Symmetry || type==Type::PressureOutlet || type==Type::PressureReservoir,
                "missing physical stress boundary");
            if(domain.axisymmetric() && !face.horizontal && face.center.x==0)
                require(type==Type::Symmetry,"true radial axis requires its symmetry law");
        }
    }
    HydroConservativeStress result;result.identity_=domain.meshIdentity();result.axisymmetric_=domain.axisymmetric();
    result.faces_=domain.nativeFaceCount();result.cells_=mesh.cellIndexNum();result.vertices_=mesh.vertIndexNum();
    const auto& duals=pressure.duals();
    std::map<int,HydroFittedDetail::NativeVelocityTrace> fits;
    const auto fit=[&](int row)->const HydroFittedDetail::NativeVelocityTrace& {
        auto found=fits.find(row);if(found!=fits.end())return found->second;
        const auto face=domain.face(mesh,row);
        const auto cells=HydroFittedDetail::nativeTraceSupport(mesh,domain,{face.negativeCell,face.positiveCell});
        return fits.emplace(row,HydroFittedDetail::nativeVelocityTrace(mesh,domain,cells,face.horizontal?1:0,face.center)).first->second;
    };
    for(const auto& segment:duals.segments()) {
        if(segment.measure==0)continue; // True axis; no finite surface traction.
        Connection connection;connection.negative=segment.negativeNative;connection.positive=segment.positiveNative;
        connection.area=segment.measure;
        const int row=segment.negativeNative>=0?segment.negativeNative:segment.positiveNative;
        const auto a=domain.face(mesh,row);const int component=a.horizontal?1:0,normal=segment.normalAxis;
        Vector2D center{(segment.begin.x+segment.end.x)/2,(segment.begin.y+segment.end.y)/2};
        if(result.axisymmetric_ && normal==1)center.x=(2./3.)*(segment.begin.x*segment.begin.x
            +segment.begin.x*segment.end.x+segment.end.x*segment.end.x)/(segment.begin.x+segment.end.x);
        if(segment.negativeNative<0 || segment.positiveNative<0) {
            const auto external=domain.face(mesh,segment.primalNative);const auto type=boundary[segment.primalNative].type;
            // Prescribed complete outlet traction belongs to the pressure/load
            // boundary, and zero shear has no additional viscous traction.
            if(type!=Type::NoSlipWall || component==normal)continue;
            const double distance=external.outwardSign*(at(external.center,normal)-at(a.center,normal));
            require(distance>0,"wall shear lacks an inward native sample");
            connection.strain={{row,-external.outwardSign/distance}};
        } else {
            const auto b=domain.face(mesh,segment.positiveNative);
            const double distance=at(b.center,normal)-at(a.center,normal),offset=at(b.center,1-normal)-at(a.center,1-normal);
            require(distance>0,"nonpositive physical stress centre distance");
            connection.strain={{b.nativeRow,1/distance},{a.nativeRow,-1/distance}};
            if(offset!=0) {
                add(connection.strain,fit(a.nativeRow).gradient[1-normal],-.5*offset/distance);
                add(connection.strain,fit(b.nativeRow).gradient[1-normal],-.5*offset/distance);
            }
            if(component==normal)for(auto& term:connection.strain)term.second*=2;
            else {
                const auto cells=HydroFittedDetail::nativeTraceSupport(mesh,domain,
                    {a.negativeCell,a.positiveCell,b.negativeCell,b.positiveCell});
                const auto cross=HydroFittedDetail::nativeVelocityTrace(mesh,domain,cells,normal,center);
                add(connection.strain,cross.gradient[component],1);
            }
        }
        if(segment.ownerCell>=0)connection.cellMaterial={{segment.ownerCell,1}};
        else {
            const auto face=domain.face(mesh,segment.primalNative);
            if(face.outwardSign!=0)connection.cellMaterial={{std::max(face.negativeCell,face.positiveCell),1}};
            else {
                const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
                const auto& start=mesh.verts(edge.edgeToVert[0]);const auto& end=mesh.verts(edge.edgeToVert[1]);
                const int tangent=face.horizontal?0:1;
                const double fraction=(at(center,tangent)-at(start.coord(),tangent))/(at(end.coord(),tangent)-at(start.coord(),tangent));
                require(fraction>=0 && fraction<=1,"stress material trace leaves its real primal face");
                connection.vertexMaterial={{start.index(),1-fraction},{end.index(),fraction}};
            }
        }
        result.connections_.push_back(std::move(connection));
    }
    if(result.axisymmetric_)for(const auto face:domain.faces(mesh))if(!face.horizontal) {
        const auto& volume=duals.volumes()[face.nativeRow];const auto r=face.center.x;
        for(int owner:{face.negativeCell,face.positiveCell})if(owner>=0) {
            const auto c=mesh.cells(mesh.cellIndexToLabel[owner]).centerCoord().x;
            const long double lo=std::min(c,r),hi=std::max(c,r);
            require(hi>lo,"empty radial viscous half-piece");
            Hoop hoop;hoop.row=face.nativeRow;hoop.owner=owner;hoop.derivative=fit(face.nativeRow).gradient[0];
            const long double factor=2*CartesianMeasure::twoPi*(volume.upper.y-volume.lower.y);
            if(lo==0) {
                // u_r(0)=0 is constrained. Its smooth odd extension gives a
                // finite integral of u_r/r; do not manufacture log(0).
                require(r==0 && boundary[face.nativeRow].type==Type::Symmetry,"singular unconstrained radial stress volume");
                hoop.slope=static_cast<double>(factor*(hi-lo));
            } else {
                const auto logarithm=std::log1p((hi-lo)/lo);
                hoop.constant=static_cast<double>(factor*logarithm);
                hoop.slope=static_cast<double>(factor*((hi-lo)-r*logarithm));
            }
            result.hoop_.push_back(std::move(hoop));
        }
    }
    return result;
}
Eigen::SparseMatrix<double> HydroConservativeStress::assemble(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const std::vector<double>& mu,const std::vector<double>& vertexMu) const
{
    domain.validateMeshIdentity(mesh);
    require(identity_ && identity_==domain.meshIdentity() && axisymmetric_==domain.axisymmetric(),"stress source/domain differs");
    require(mu.size()==static_cast<std::size_t>(cells_) && vertexMu.size()==static_cast<std::size_t>(vertices_),"viscosity storage mismatch");
    for(double value:mu)require(std::isfinite(value) && value>=0,"invalid cell viscosity");
    for(double value:vertexMu)require(std::isfinite(value) && value>=0,"invalid vertex viscosity");
    std::vector<Eigen::Triplet<double>> entries;
    const auto append=[&](int row,int column,double value) {
        require(std::isfinite(value),"viscous coefficient overflow");if(value!=0)entries.emplace_back(row,column,value);
    };
    for(const auto& connection:connections_) {
        double viscosity=0;
        for(const auto& term:connection.cellMaterial)viscosity+=term.second*mu[term.first];
        for(const auto& term:connection.vertexMaterial)viscosity+=term.second*vertexMu[term.first];
        const double scale=viscosity*connection.area;
        for(const auto& term:connection.strain) {
            if(connection.negative>=0)append(connection.negative,term.first,scale*term.second);
            if(connection.positive>=0)append(connection.positive,term.first,-scale*term.second);
        }
    }
    for(const auto& hoop:hoop_) {
        append(hoop.row,hoop.row,-mu[hoop.owner]*hoop.constant);
        for(const auto& term:hoop.derivative)append(hoop.row,term.first,-mu[hoop.owner]*hoop.slope*term.second);
    }
    Eigen::SparseMatrix<double> matrix(faces_,faces_);matrix.setFromTriplets(entries.begin(),entries.end());return matrix;
}
}
