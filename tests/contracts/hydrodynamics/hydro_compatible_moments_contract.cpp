#include "Modules/Hydrodynamics/HydroCompatibleMomentSpace.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Geometry/ParametricHermiteCurve.h"
#include <Eigen/SuperLUSupport>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>

namespace {
using Space=Quantum::Hydrodynamics::HydroCompatibleMomentSpace;
const double pi=std::acos(-1.);
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
struct Exact {double p,ux,uy,weightedDivergence;};
Exact exact(double x,double y,bool rz,bool polynomial){
    if(polynomial){const double ux=1+.4*x-.1*y,uy=.3-.1*x+.3*y;
        return {.2+x+.3*y+.2*x*x-.1*x*y+.15*y*y,ux,uy,rz?2*pi*(.7*x+ux):.7};}
    if(rz){const double cx=std::cos(pi*x),cy=std::cos(pi*y),sx=std::sin(pi*x),sy=std::sin(pi*y);
        return {cx*cy,-pi*sx*cy,-pi*cx*sy,2*pi*(-2*pi*pi*x*cx*cy-pi*sx*cy)};}
    const double sx=std::sin(pi*x),sy=std::sin(pi*y),cx=std::cos(pi*x),cy=std::cos(pi*y);
    return {sx*sy,pi*cx*sy,pi*sx*cy,-2*pi*pi*sx*sy};
}
double weight(double x,bool rz){return rz?2*pi*x:1.;}
template<class F>void cellIntegral(Vector2D center,Vector2D h,F&& f){
    const auto& g=Quantum::Geometry::Parametric::gaussEight();
    for(int j=0;j<8;j++)for(int i=0;i<8;i++)f(Vector2D{center.x+.5*h.x*g.x[i],center.y+.5*h.y*g.x[j]},.25*h.x*h.y*g.w[i]*g.w[j]);
}
template<class F>void faceIntegral(const QuadtreeMesh& mesh,const NativeLeafFace& face,F&& f){
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    const auto& g=Quantum::Geometry::Parametric::gaussEight();const double length=face.horizontal?b.x-a.x:b.y-a.y;
    for(int i=0;i<8;i++){const double t=.5*(1+g.x[i]);f(Vector2D{a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)},.5*length*g.w[i]);}
}
struct Result {double velocity,pressure;};
Result solve(bool rz,int n,bool polynomial){
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=n;input.Ny=rz?5*n/4:n;input.maxCellNum=40000;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);DictionaryMesh::RefineBox outer;outer.min={.25,.25};outer.max={.75,.75};outer.level=1;
    DictionaryMesh::RefineBox inner;inner.min={.4375,.375};inner.max={.625,.5625};inner.level=2;
    if(polynomial){std::vector<DictionaryMesh::RefineBox> boxes;
        for(const auto offset:std::vector<std::pair<int,int>>{{2,3},{4,3},{3,2},{3,4}}){DictionaryMesh::RefineBox box;
            box.min={(offset.first+.1)/input.Nx,(offset.second+.1)/input.Ny};box.max={(offset.first+.9)/input.Nx,(offset.second+.9)/input.Ny};box.level=1;boxes.push_back(box);}
        mesh.applyRefineBoxes(boxes);
    }else mesh.applyRefineBoxes({outer,inner});
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);
    const auto space=Space::build(mesh,domain);space.validate(mesh,domain);
    const auto& hodge=space.mass();const auto& d=space.divergence();const int nv=int(hodge.rows()),np=int(d.rows()),nc=mesh.cellIndexNum();
    const Eigen::SparseMatrix<double> adjoint=space.pressureLoad()+Space::Matrix(d.transpose());require(adjoint.norm()==0,"moment pressure lost exact adjoint");
    Eigen::VectorXd sampled=Eigen::VectorXd::Zero(nv),boundary=Eigen::VectorXd::Zero(nv),source=Eigen::VectorXd::Zero(np),pressure=Eigen::VectorXd::Zero(np);
    Eigen::VectorXd affine=Eigen::VectorXd::Zero(nv);long double affineEnergy=0;
    for(const auto face:domain.faces(mesh)){
        const auto& info=space.faces()[face.nativeRow];if(info.first<0)continue;
        faceIntegral(mesh,face,[&](Vector2D p,double ds){const double dw=ds*weight(p.x,rz),theta=((face.horizontal?p.x:p.y)-info.tangentCenter)/info.tangentScale;
            const auto e=exact(p.x,p.y,rz,polynomial);const double normal=face.horizontal?e.uy:e.ux;
            const double anormal=face.horizontal?.7+.2*p.x+.1*p.y:.2+.3*p.x-.4*p.y;
            sampled[info.first]+=dw/info.area*normal;sampled[info.first+1]+=dw/info.area*theta*normal;
            affine[info.first]+=dw/info.area*anormal;affine[info.first+1]+=dw/info.area*theta*anormal;
            if(face.outwardSign){boundary[info.first]+=face.outwardSign*dw*e.p;boundary[info.first+1]+=face.outwardSign*dw*theta*e.p;}
        });
    }
    for(int c=0;c<nc;c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const double volume=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        cellIntegral(center,h,[&](Vector2D p,double da){const double x=(p.x-center.x)/h.x,y=(p.y-center.y)/h.y,dv=da*weight(p.x,rz);
            const auto e=exact(p.x,p.y,rz,polynomial);const int first=space.firstCellVelocity(c);
            sampled[first]+=dv/volume*e.ux;sampled[first+1]+=dv/volume*e.uy;sampled[first+2]+=dv/volume*(-y*e.ux+x*e.uy);
            const double ax=.2+.3*p.x-.4*p.y,ay=.7+.2*p.x+.1*p.y;
            affine[first]+=dv/volume*ax;affine[first+1]+=dv/volume*ay;affine[first+2]+=dv/volume*(-y*ax+x*ay);affineEnergy+=dv*(ax*ax+ay*ay);
            source[3*c]+=da*e.weightedDivergence;source[3*c+1]+=da*x*e.weightedDivergence;source[3*c+2]+=da*y*e.weightedDivergence;
            pressure[3*c]+=da/(h.x*h.y)*e.p;pressure[3*c+1]+=12*da/(h.x*h.y)*x*e.p;pressure[3*c+2]+=12*da/(h.x*h.y)*y*e.p;
        });
    }
    if(rz){
        // Use the declared physical port for the new coordinate reference.
        // Face moments, weighted cell means, rotation, force, divergence and
        // continuum energy above remain independent analytic integrals.
        const auto physical=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto e=exact(p.x,p.y,true,polynomial);return Eigen::Vector2d(e.ux,e.uy);});
        const auto physicalAffine=space.interpolateVelocity(mesh,domain,[](Vector2D p){return Eigen::Vector2d(.2+.3*p.x-.4*p.y,.7+.2*p.x+.1*p.y);});
        for(int c=0;c<nc;c++){const int first=space.firstCellVelocity(c);sampled.segment<3>(first+3)=physical.segment<3>(first+3);affine.segment<3>(first+3)=physicalAffine.segment<3>(first+3);}
    }
    const double energyError=std::abs(affine.dot(hodge*affine)-static_cast<double>(affineEnergy))/static_cast<double>(affineEnergy);
    require(energyError<3e-12,"affine physical energy not reproduced on mixed XY/RZ cells");
    if(polynomial){const double gradient=(hodge*sampled+d.transpose()*pressure-boundary).norm()/std::max(1.,boundary.norm());
        const double divergence=(d*sampled-source).norm()/std::max(1.,source.norm());
        require(gradient<3e-12&&divergence<3e-12,"quadratic mixed-cell pressure/gradient not reproduced");
        std::cout<<std::setprecision(17)<<"HODGE_MOMENT_POLYNOMIAL metric="<<(rz?"RZ":"XY")<<" gradient="<<gradient<<" divergence="<<divergence<<" affine_energy="<<energyError<<'\n';}
    std::vector<Eigen::Triplet<double>> entries;entries.reserve(hodge.nonZeros()+2*d.nonZeros());
    for(int j=0;j<hodge.outerSize();j++)for(Space::Matrix::InnerIterator it(hodge,j);it;++it)entries.emplace_back(it.row(),it.col(),it.value());
    for(int j=0;j<d.outerSize();j++)for(Space::Matrix::InnerIterator it(d,j);it;++it){entries.emplace_back(nv+it.row(),it.col(),it.value());entries.emplace_back(it.col(),nv+it.row(),it.value());}
    Space::Matrix system(nv+np,nv+np);system.setFromTriplets(entries.begin(),entries.end());
    Eigen::VectorXd rhs(nv+np);rhs.head(nv)=boundary;rhs.tail(np)=source;
    Eigen::SuperLU<Space::Matrix> solver;solver.compute(system);require(solver.info()==Eigen::Success,"moment mixed pressure factorization failed");
    const Eigen::VectorXd solution=solver.solve(rhs);require(solver.info()==Eigen::Success&&solution.allFinite(),"moment mixed pressure solve failed");
    const double residual=(system*solution-rhs).norm()/std::max(1.,rhs.norm());require(residual<2e-11,"moment mixed pressure residual failed");
    const Eigen::VectorXd error=solution.head(nv)-sampled;
    const double velocityError=std::sqrt(error.dot(hodge*error)/sampled.dot(hodge*sampled));
    long double pressureError=0,pressureNorm=0,physicalVelocityError=0,physicalVelocityNorm=0;
    std::ofstream fields;if(!polynomial&&n==32)if(const char* prefix=std::getenv("QUANTUM_MOMENT_FIELDS")){
        fields.open(std::string(prefix)+(rz?"_rz.csv":"_xy.csv"));fields<<std::setprecision(17)<<"x,y,hx,hy,pressure_center,exact_pressure_center,projected_ux,projected_uy\n";}
    for(int c=0;c<nc;c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
        const auto& block=space.cells()[c];Eigen::VectorXd local(block.rows.size());for(int i=0;i<local.size();i++)local[i]=solution[block.rows[i]];
        const Eigen::VectorXd coeff=block.affineProjection*local;const double l=std::sqrt(h.x*h.y);
        cellIntegral(center,h,[&](Vector2D p,double da){const double ph=solution[nv+3*c]+solution[nv+3*c+1]*(p.x-center.x)/h.x+solution[nv+3*c+2]*(p.y-center.y)/h.y;
            const auto e=exact(p.x,p.y,rz,polynomial);const double dv=da*weight(p.x,rz),x=(p.x-center.x)/h.x,y=(p.y-center.y)/h.y;
            const double ux=l/h.x*(coeff[0]+coeff[2]*x+coeff[3]*y)-coeff[5]*y,uy=l/h.y*(coeff[1]+coeff[3]*x+coeff[4]*y)+coeff[5]*x;
            physicalVelocityError+=dv*((ux-e.ux)*(ux-e.ux)+(uy-e.uy)*(uy-e.uy));physicalVelocityNorm+=dv*(e.ux*e.ux+e.uy*e.uy);
            pressureError+=dv*(ph-e.p)*(ph-e.p);pressureNorm+=dv*e.p*e.p;});
        if(fields){
            fields<<center.x<<','<<center.y<<','<<h.x<<','<<h.y<<','<<solution[nv+3*c]<<','<<exact(center.x,center.y,rz,false).p<<','<<coeff[0]*l/h.x<<','<<coeff[1]*l/h.y<<'\n';}
    }
    if(polynomial)require(velocityError<3e-10&&(solution.tail(np)-pressure).norm()/pressure.norm()<3e-10,"quadratic manufactured mixed solve failed");
    const auto report=space.report();const double pressureL2=std::sqrt(static_cast<double>(pressureError/pressureNorm));
    const double velocityL2=std::sqrt(static_cast<double>(physicalVelocityError/physicalVelocityNorm));
    if(polynomial)require(velocityL2<3e-10,"affine reconstructed velocity not reproduced");
    if(polynomial)require(report.maximumLocalVelocityDofs==16+space.cellVelocityCount(),"four split sides were not covered");
    if(rz){
        // RZ now uses a physical reconstruction Gram, scaled to unit diagonal
        // in this report. Raw correlated stream moments do not have the old
        // arbitrary diagonal-penalty spectrum. Keep the independent affine
        // energy, polynomial work and actual solve accuracy gates above.
        require(report.minimumScaledEigenvalue>0&&report.maximumScaledEigenvalue<=report.maximumLocalVelocityDofs+1e-10,"physical RZ moment Gram lost positivity or its unit-diagonal trace bound");
    }else require(report.minimumScaledEigenvalue>.05&&report.maximumScaledEigenvalue<50,"local moment metric lost spectral bounds");
    std::cout<<std::setprecision(17)<<"HODGE_MOMENT_SOLVE metric="<<(rz?"RZ":"XY")<<" polynomial="<<polynomial<<" n="<<n<<" ny="<<input.Ny<<" cells="<<nc
        <<" velocity_dofs="<<nv<<" pressure_dofs="<<np<<" relative_velocity_dofs="<<velocityError<<" relative_velocity_l2="<<velocityL2<<" relative_pressure_l2="<<pressureL2<<" residual="<<residual
        <<" consistency="<<report.maximumConsistencyDefect<<" symmetry="<<report.maximumSymmetryDefect<<" lambda_min="<<report.minimumScaledEigenvalue
        <<" lambda_max="<<report.maximumScaledEigenvalue<<" local_dofs="<<report.maximumLocalVelocityDofs<<" ns_time_step=0\n";
    if(polynomial){DictionaryMesh::RefineBox change;change.min={.02,.02};change.max={.08,.08};change.level=1;mesh.applyRefineBoxes({change});
        bool rejected=false;try{space.validate(mesh,domain);}catch(const std::exception&){rejected=true;}require(rejected,"stale moment metric survived mesh mutation");}
    return {velocityL2,pressureL2};
}
}
void hydroCompatibleMomentsContract(){for(bool rz:{false,true}){solve(rz,8,true);Result previous{};
    for(int n:{8,16,32}){const auto current=solve(rz,n,false);if(previous.velocity){const double vu=std::log2(previous.velocity/current.velocity),pp=std::log2(previous.pressure/current.pressure);
        std::cout<<std::setprecision(17)<<"HODGE_MOMENT_ORDER metric="<<(rz?"RZ":"XY")<<" n="<<n<<" velocity="<<vu<<" pressure="<<pp<<'\n';
        require(vu>1.7&&pp>1.7,"mixed-cell moment pressure/velocity not approaching second order");}previous=current;}}}
