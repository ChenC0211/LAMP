#include "Modules/Hydrodynamics/HydroCompatibleMomentTransport.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentMaterial.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentBoundary.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Quantum::Hydrodynamics {
namespace {
using Matrix=HydroCompatibleMomentSpace::Matrix;
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
Eigen::VectorXd localValues(const HydroCompatibleMomentSpace::CellBlock& cell,const Eigen::VectorXd& global){
    Eigen::VectorXd result(cell.rows.size());for(int i=0;i<result.size();i++)result[i]=global[cell.rows[i]];return result;
}
void scatter(std::vector<Eigen::Triplet<double>>& entries,const std::vector<int>& rows,const std::vector<int>& cols,const Eigen::MatrixXd& block){
    for(int i=0;i<block.rows();i++)for(int j=0;j<block.cols();j++)if(block(i,j)!=0)entries.emplace_back(rows[i],cols[j],block(i,j));
}
constexpr double gx[3]={-.774596669241483377,0,.774596669241483377};
constexpr double gw[3]={5./9,8./9,5./9};
}
Eigen::Matrix<double,2,12> HydroCompatibleMomentTransport::rtBasis(Vector2D p){
    const double a[3]={1,std::sqrt(12.)*p.x,std::sqrt(180.)*(p.x*p.x-1./12)};
    const double b[3]={1,std::sqrt(12.)*p.y,std::sqrt(180.)*(p.y*p.y-1./12)};
    Eigen::Matrix<double,2,12> result=Eigen::Matrix<double,2,12>::Zero();
    for(int j=0;j<2;j++)for(int i=0;i<3;i++){result(0,3*j+i)=a[i]*b[j];result(1,6+3*j+i)=a[j]*b[i];}return result;
}
HydroCompatibleMomentTransport HydroCompatibleMomentTransport::build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space){
    space.validate(mesh,domain);if(domain.axisymmetric())return buildAxisymmetric(mesh,domain,space);
    HydroCompatibleMomentTransport result;result.identity_=domain.meshIdentity();result.velocityCount_=int(space.mass().rows());result.cells_.resize(mesh.cellIndexNum());
    std::vector<std::vector<int>> incident(mesh.cellIndexNum());
    for(const auto f:domain.faces(mesh)){if(f.negativeCell>=0)incident[f.negativeCell].push_back(f.nativeRow);if(f.positiveCell>=0)incident[f.positiveCell].push_back(f.nativeRow);}
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const auto& bulk=space.cells()[c];
        const int side=incident[c].size()>4?2:1,pieces=side*side,nlocal=int(bulk.rows.size());const Vector2D hs{h.x/side,h.y/side};
        const int constraints=8*side+4*side*(side-1)+4*pieces+3; // boundary, inner faces, divergence, moments
        Eigen::MatrixXd a=Eigen::MatrixXd::Zero(constraints,12*pieces),b=Eigen::MatrixXd::Zero(constraints,nlocal);
        auto& saved=result.cells_[c];saved.pieces.resize(pieces);
        for(int j=0;j<side;j++)for(int i=0;i<side;i++)saved.pieces[j*side+i]={{center.x+h.x*((i+.5)/side-.5),center.y+h.y*((j+.5)/side-.5)},hs,{}};
        int row=0;
        const auto trace=[&](int r,int piece,bool horizontal,int sign,double coefficient){const int first=12*piece+(horizontal?6:0);
            for(int m=0;m<2;m++){a(r+m,first+3*m)+=coefficient;a(r+m,first+3*m+1)+=coefficient*sign*std::sqrt(3.);a(r+m,first+3*m+2)+=coefficient*std::sqrt(5.);}};
        for(int direction=0;direction<2;direction++)for(int sign:{-1,1})for(int s=0;s<side;s++){
            const bool horizontal=direction==1;const int i=horizontal?s:(sign<0?0:side-1),j=horizontal?(sign<0?0:side-1):s,piece=j*side+i;
            const double tcenter=horizontal?saved.pieces[piece].center.x:saved.pieces[piece].center.y;
            const double length=horizontal?hs.x:hs.y;int found=-1;
            for(int k=0;k<int(incident[c].size());k++){
                const auto f=domain.face(mesh,incident[c][k]);if(f.horizontal!=horizontal||(f.negativeCell==c?1:-1)!=sign)continue;
                const double fc=horizontal?f.center.x:f.center.y;
                if(std::abs(tcenter-fc)<.5*f.area+1e-13*length){
                    require(std::abs(f.area-length)<1e-11*length||std::abs(f.area-2*length)<1e-11*length,"moment transport needs at most 2:1 split sides");
                    require(std::abs(tcenter-fc)+.5*length<=.5*f.area+1e-11*length,"RT boundary piece crosses a native subface");found=k;break;
                }
            }
            require(found>=0,"RT reconstruction could not match a native boundary face");const auto& info=space.faces()[incident[c][found]];
            trace(row,piece,horizontal,sign,1);b(row,2*found)=1;b(row,2*found+1)=(tcenter-info.tangentCenter)/info.tangentScale;
            b(row+1,2*found+1)=length/(std::sqrt(12.)*info.tangentScale);row+=2;
        }
        for(int j=0;j<side;j++)for(int i=0;i+1<side;i++){trace(row,j*side+i,false,1,1);trace(row,j*side+i+1,false,-1,-1);row+=2;}
        for(int j=0;j+1<side;j++)for(int i=0;i<side;i++){trace(row,j*side+i,true,1,1);trace(row,(j+1)*side+i,true,-1,-1);row+=2;}
        Eigen::MatrixXd divergence(3,nlocal);
        for(int k=0;k<3;k++)for(int i=0;i<nlocal;i++)divergence(k,i)=space.divergence().coeff(3*c+k,bulk.rows[i])*(k?12:1)/(h.x*h.y);
        const double scale=std::sqrt(hs.x*hs.y),dx=scale/hs.x,dy=scale/hs.y;
        for(int p=0;p<pieces;p++){
            const int first=12*p;const auto pc=saved.pieces[p].center;
            a(row,first+1)=std::sqrt(12.)*dx;a(row,first+7)=std::sqrt(12.)*dy;
            a(row+1,first+2)=std::sqrt(60.)*dx;a(row+1,first+10)=std::sqrt(12.)*dy;
            a(row+2,first+4)=std::sqrt(12.)*dx;a(row+2,first+8)=std::sqrt(60.)*dy;
            a(row+3,first+5)=std::sqrt(60.)*dx;a(row+3,first+11)=std::sqrt(60.)*dy;
            b.row(row)=scale*(divergence.row(0)+(pc.x-center.x)/h.x*divergence.row(1)+(pc.y-center.y)/h.y*divergence.row(2));
            b.row(row+1)=scale/(side*std::sqrt(12.))*divergence.row(1);b.row(row+2)=scale/(side*std::sqrt(12.))*divergence.row(2);row+=4;
        }
        for(int p=0;p<pieces;p++){
            const int first=12*p;const auto pc=saved.pieces[p].center;const double fraction=1./pieces;
            a(row,first)=a(row+1,first+6)=fraction;
            a(row+2,first)=-fraction*(pc.y-center.y)/h.y;a(row+2,first+6)=fraction*(pc.x-center.x)/h.x;
            a(row+2,first+3)=-fraction/(side*std::sqrt(12.));a(row+2,first+9)=fraction/(side*std::sqrt(12.));
        }
        for(int k=0;k<3;k++)b(row+k,nlocal-3+k)=1;row+=3;require(row==constraints,"RT constraint count differs");
        // Equal-volume pieces and orthonormal RT basis: Euclidean coefficient
        // minimum is precisely minimum physical L2 norm. Keep all constraints;
        // the three global divergence/moment dependencies are not deleted.
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(a,Eigen::ComputeThinU|Eigen::ComputeThinV);
        const double tolerance=4096*std::numeric_limits<double>::epsilon();int rank=0;
        for(int k=0;k<svd.singularValues().size();k++)if(svd.singularValues()[k]>tolerance*svd.singularValues()[0])rank++;
        require(rank==(side==1?12:40),"RT moment constraints have unexpected rank");
        const Eigen::MatrixXd coefficients=svd.matrixV().leftCols(rank)*svd.singularValues().head(rank).cwiseInverse().asDiagonal()*svd.matrixU().leftCols(rank).transpose()*b;
        saved.constraintDefect=(a*coefficients-b).norm()/b.norm();require(saved.constraintDefect<8*tolerance,"incompatible RT moment constraints");
        for(int p=0;p<pieces;p++)saved.pieces[p].coefficients=coefficients.middleRows(12*p,12);
    }
    return result;
}
void HydroCompatibleMomentTransport::validate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space) const{
    space.validate(mesh,domain);require(identity_&&identity_==domain.meshIdentity()&&axisymmetric_==domain.axisymmetric()&&velocityCount_==space.mass().rows(),"moment transport/domain source differs");
}
HydroCompatibleMomentTransport::Assembly HydroCompatibleMomentTransport::assemble(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,const Eigen::VectorXd& wind,double rho,Flux flux,const HydroCompatibleMomentBoundary* boundary) const{
    validate(mesh,domain,space);require(wind.size()==velocityCount_&&wind.allFinite()&&std::isfinite(rho)&&rho>0,"invalid moment advecting state");
    if(boundary)boundary->validate(mesh,domain,space);
    require(flux==Flux::Central||flux==Flux::TraceUpwind,"invalid moment convection flux");
    const Eigen::VectorXd divergence=space.divergence()*wind,scale=space.divergence().cwiseAbs()*wind.cwiseAbs();
    require(divergence.stableNorm()<=1e-10*scale.stableNorm(),"moment convection requires incompressible advecting moments");
    Assembly result;const int nc=mesh.cellIndexNum();result.cellMomentum.resize(nc,Eigen::Vector2d::Zero());result.cellBoundaryMomentum=result.cellMomentum;
    result.boundaryLoad=Eigen::VectorXd::Zero(velocityCount_);
    std::vector<Eigen::Triplet<double>> central,dissipation;std::vector<Eigen::VectorXd> loads(nc),values(nc);
    for(int c=0;c<nc;c++){
        const auto& bulk=space.cells()[c];const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const double l=std::sqrt(h.x*h.y);
        values[c]=localValues(bulk,wind);Eigen::Matrix<double,2,6> x=Eigen::Matrix<double,2,6>::Zero(),y=x;
        x(0,2)=l/(h.x*h.x);x(1,3)=l/(h.x*h.y);x(1,5)=1/h.x;
        y(0,3)=l/(h.x*h.y);y(0,5)=-1/h.y;y(1,4)=l/(h.y*h.y);
        const Eigen::MatrixXd gradX=x*bulk.affineProjection,gradY=y*bulk.affineProjection;
        Eigen::MatrixXd local=Eigen::MatrixXd::Zero(bulk.rows.size(),bulk.rows.size());
        for(int j=0;j<3;j++)for(int i=0;i<3;i++){
            const Eigen::MatrixXd basis=HydroCompatibleMomentSpace::affineVelocityBasis({.5*gx[i],.5*gx[j]},h)*bulk.affineProjection;
            const Eigen::Vector2d velocity=basis*values[c];
            // The RT reconstruction has exactly this P1 L2 projection, so
            // these polynomial work integrals require no microcell quadrature.
            const double weight=domain.axisymmetric()?CartesianMeasure::twoPi*(center.x+.5*h.x*gx[i]):1.;
            local-=.25*rho*h.x*h.y*gw[i]*gw[j]*weight*(velocity[0]*gradX+velocity[1]*gradY).transpose()*basis;
        }
        loads[c]=local*values[c];scatter(central,bulk.rows,bulk.rows,local);
    }
    for(const auto f:domain.faces(mesh)){
        const auto& info=space.faces()[f.nativeRow];if(info.first<0)continue;
        const double mean=wind[info.first],linear=wind[info.first+1];
        if(f.outwardSign){
            using Kind=HydroCompatibleMomentBoundary::Kind;
            if(!boundary||HydroCompatibleMomentBoundary::impermeable(boundary->faces()[f.nativeRow].kind)){require(mean==0&&linear==0,"moment convection boundary is not closed");continue;}
            const auto& port=boundary->faces()[f.nativeRow];const int c=f.negativeCell>=0?f.negativeCell:f.positiveCell;
            const auto& bulk=space.cells()[c];const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
            const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            const auto va=mesh.verts(edge.edgeToVert[0]).coord(),vb=mesh.verts(edge.edgeToVert[1]).coord();
            const double a=f.horizontal?va.x:va.y,b=f.horizontal?vb.x:vb.y;
            const double q0=f.outwardSign*(mean+linear*(a-info.tangentCenter)/info.tangentScale);
            const double q1=f.outwardSign*(mean+linear*(b-info.tangentCenter)/info.tangentScale);
            if(port.kind==Kind::PressureOutlet)require(port.incomingVelocity||std::min(q0,q1)>=0,"moment pressure outlet backflow requires explicit incoming momentum data");
            Eigen::MatrixXd faceC=Eigen::MatrixXd::Zero(bulk.rows.size(),bulk.rows.size());
            Eigen::VectorXd supplied=Eigen::VectorXd::Zero(bulk.rows.size());
            for(int g=0;g<3;g++){
                const double s=.5*(a+b)+.5*(b-a)*gx[g],ds=.5*(b-a)*gw[g]*(domain.axisymmetric()?CartesianMeasure::twoPi*(f.horizontal?s:f.center.x):1.);
                const Vector2D point=f.horizontal?Vector2D{s,f.center.y}:Vector2D{f.center.x,s};
                const double q=f.outwardSign*(mean+linear*(s-info.tangentCenter)/info.tangentScale),mass=rho*q;
                const Eigen::MatrixXd projection=HydroCompatibleMomentSpace::affineVelocityBasis({(point.x-center.x)/h.x,(point.y-center.y)/h.y},h)*bulk.affineProjection;
                const Eigen::Vector2d inside=projection*values[c];Eigen::Vector2d upwind=inside;
                if(q<0){
                    require(bool(port.incomingVelocity),"moment inflow has no exterior velocity");upwind<<port.incomingVelocity->x,port.incomingVelocity->y;
                    supplied-=ds*mass*projection.transpose()*upwind;
                    result.boundaryInflowDissipation-=.5*ds*mass*(inside-upwind).squaredNorm();
                }else faceC+=ds*mass*projection.transpose()*projection;
                result.boundaryKineticFlux+=.5*ds*mass*upwind.squaredNorm();result.boundaryVolumeRate+=ds*q;
                result.cellBoundaryMomentum[c]+=ds*mass*upwind;
            }
            loads[c]+=faceC*values[c]-supplied;scatter(central,bulk.rows,bulk.rows,faceC);
            for(int i=0;i<supplied.size();i++)result.boundaryLoad[bulk.rows[i]]+=supplied[i];
            continue;
        }
        const int minus=f.negativeCell,plus=f.positiveCell;const auto& bm=space.cells()[minus];const auto& bp=space.cells()[plus];
        const auto& cm=mesh.cells(mesh.cellIndexToLabel[minus]);const auto& cp=mesh.cells(mesh.cellIndexToLabel[plus]);
        const auto mc=cm.centerCoord(),mh=cm.size(),pc=cp.centerCoord(),ph=cp.size();const int nm=int(bm.rows.size()),np=int(bp.rows.size()),normal=f.horizontal?1:0,tangent=1-normal;
        std::vector<int> rows=bm.rows;rows.insert(rows.end(),bp.rows.begin(),bp.rows.end());
        Eigen::MatrixXd faceC=Eigen::MatrixXd::Zero(nm+np,nm+np),faceS=faceC;
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);const auto va=mesh.verts(edge.edgeToVert[0]).coord(),vb=mesh.verts(edge.edgeToVert[1]).coord();
        const double a=f.horizontal?va.x:va.y,b=f.horizontal?vb.x:vb.y;
        std::vector<double> breaks{a,b};if(linear!=0){const double root=info.tangentCenter-mean*info.tangentScale/linear;if(root>a&&root<b)breaks.insert(breaks.begin()+1,root);}
        for(int part=0;part+1<int(breaks.size());part++)for(int g=0;g<3;g++){
            const double s=.5*(breaks[part]+breaks[part+1])+.5*(breaks[part+1]-breaks[part])*gx[g];
            const double ds=.5*(breaks[part+1]-breaks[part])*gw[g]*(domain.axisymmetric()?CartesianMeasure::twoPi*(f.horizontal?s:f.center.x):1.);
            const Vector2D p=f.horizontal?Vector2D{s,f.center.y}:Vector2D{f.center.x,s};const double theta=(s-info.tangentCenter)/info.tangentScale,mass=rho*(mean+linear*theta);
            const Eigen::MatrixXd pm=HydroCompatibleMomentSpace::affineVelocityBasis({(p.x-mc.x)/mh.x,(p.y-mc.y)/mh.y},mh)*bm.affineProjection;
            const Eigen::MatrixXd pp=HydroCompatibleMomentSpace::affineVelocityBasis({(p.x-pc.x)/ph.x,(p.y-pc.y)/ph.y},ph)*bp.affineProjection;
            Eigen::MatrixXd jump(2,nm+np),average(2,nm+np);jump<<pm,-pp;average<<.5*pm,.5*pp;
            faceC+=ds*mass*jump.transpose()*average;
            Eigen::Vector2d momentum=mass*.5*(pm*values[minus]+pp*values[plus]);
            loads[minus]+=ds*pm.transpose()*momentum;loads[plus]-=ds*pp.transpose()*momentum;
            if(flux==Flux::TraceUpwind){
                const double absolute=std::abs(mass);
                // One face-trace residual: the common normal component is the
                // owned face polynomial; tangential component is the centered
                // physical projection. Sum of its two squared side residuals.
                Eigen::VectorXd rm=-pm.row(normal).transpose(),rp=-pp.row(normal).transpose();
                for(int i=0;i<nm;i++){if(bm.rows[i]==info.first)rm[i]+=1;if(bm.rows[i]==info.first+1)rm[i]+=theta;}
                for(int i=0;i<np;i++){if(bp.rows[i]==info.first)rp[i]+=1;if(bp.rows[i]==info.first+1)rp[i]+=theta;}
                faceS.topLeftCorner(nm,nm)+=ds*absolute*rm*rm.transpose();faceS.bottomRightCorner(np,np)+=ds*absolute*rp*rp.transpose();
                faceS+=.5*ds*absolute*jump.row(tangent).transpose()*jump.row(tangent);
                loads[minus]+=ds*absolute*rm*rm.dot(values[minus]);loads[plus]+=ds*absolute*rp*rp.dot(values[plus]);
                Eigen::Vector2d upwind=Eigen::Vector2d::Zero();upwind[tangent]=.5*absolute*((pm*values[minus])[tangent]-(pp*values[plus])[tangent]);
                loads[minus]+=ds*pm.transpose()*upwind;loads[plus]-=ds*pp.transpose()*upwind;momentum+=upwind;
            }
            result.cellBoundaryMomentum[minus]+=ds*momentum;result.cellBoundaryMomentum[plus]-=ds*momentum;
        }
        scatter(central,rows,rows,faceC);scatter(dissipation,rows,rows,faceS);
    }
    for(int c=0;c<nc;c++){
        const auto& h=mesh.cells(mesh.cellIndexToLabel[c]).size();const double l=std::sqrt(h.x*h.y);const auto& samples=space.cells()[c].affineSamples;
        result.cellMomentum[c]<<h.x/l*samples.col(0).dot(loads[c]),h.y/l*samples.col(1).dot(loads[c]);
        result.maximumLocalMomentumDefect=std::max(result.maximumLocalMomentumDefect,(result.cellMomentum[c]-result.cellBoundaryMomentum[c]).norm());
    }
    result.central.resize(velocityCount_,velocityCount_);result.central.setFromTriplets(central.begin(),central.end());
    result.dissipation.resize(velocityCount_,velocityCount_);result.dissipation.setFromTriplets(dissipation.begin(),dissipation.end());
    return result;
}
HydroCompatibleMomentTransport::Assembly HydroCompatibleMomentTransport::assemblePhase(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,
    const Eigen::VectorXd& wind,const InterfaceState& previous,const InterfaceState& next,const InterfaceStepTime& time,
    const HydroCompatibleMomentMaterial& before,const HydroCompatibleMomentMaterial& midpoint,const HydroCompatibleMomentMaterial& after,
    const HydroCompatibleMomentBoundary* boundary) const {
    validate(mesh,domain,space);
    if(boundary)boundary->validate(mesh,domain,space);
    for(const auto* material:{&before,&midpoint,&after})material->validate(mesh,domain,space);
    require(before.boundaryHandle()->matches(*previous.parametricGeometry().boundaryHandle())&&after.boundaryHandle()->matches(*next.parametricGeometry().boundaryHandle()),
        "moment transport material is not the actual endpoint interface");
    const auto equalTime=[](double a,double b){return std::abs(a-b)<=16*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(a),std::abs(b)});};
    require(time.dt>0&&equalTime(before.geometryTime(),time.begin)&&equalTime(after.geometryTime(),time.end)
        &&equalTime(midpoint.geometryTime(),time.begin+.5*time.dt),"moment material time layers differ from phase transport");
    const auto properties=before.properties();
    for(const auto* material:{&midpoint,&after}){const auto p=material->properties();
        require(p.liquidDensity==properties.liquidDensity&&p.gasDensity==properties.gasDensity
            &&p.liquidViscosity==properties.liquidViscosity&&p.gasViscosity==properties.gasViscosity,"phase properties changed during momentum transport");}
    require(wind.size()==velocityCount_&&wind.allFinite(),"invalid phase midpoint velocity");
    require((space.divergence()*wind).stableNorm()<=1e-10*(space.divergence().cwiseAbs()*wind.cwiseAbs()).stableNorm(),
        "phase momentum requires an incompressible midpoint");
    const auto flux=phaseMassFlux(mesh,domain,space,previous,next,time,properties.liquidDensity,properties.gasDensity);
    Assembly result;result.phaseMassClosure=flux.maximumRelativeClosure;result.boundaryLoad=Eigen::VectorXd::Zero(velocityCount_);
    const int nc=mesh.cellIndexNum();std::vector<Eigen::Matrix3d> faceGram(nc,Eigen::Matrix3d::Zero());
    std::vector<Eigen::VectorXd> values(nc),loads(nc);std::vector<Eigen::Triplet<double>> central,boundaryMetric;
    result.cellMomentum.resize(nc,Eigen::Vector2d::Zero());result.cellBoundaryMomentum=result.cellMomentum;
    for(int c=0;c<nc;c++){values[c]=localValues(space.cells()[c],wind);loads[c]=Eigen::VectorXd::Zero(values[c].size());}
    for(const auto face:domain.faces(mesh)){
        const auto& info=space.faces()[face.nativeRow];const auto& f=flux.face[face.nativeRow];
        if(info.first<0){require(f.isZero(0),"inactive axis has nonzero mass flux");continue;}
        Eigen::MatrixXd pm[2],pp[2];
        const auto trace=[&](int c,Eigen::MatrixXd* b,double sign){
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto h=cell.size(),x=cell.centerCoord();
            const Vector2D point=face.horizontal?Vector2D{info.tangentCenter,face.center.y}:Vector2D{face.center.x,info.tangentCenter};
            const Vector2D a{(point.x-x.x)/h.x,(point.y-x.y)/h.y};
            const Vector2D t=face.horizontal?Vector2D{info.tangentScale/h.x,0}:Vector2D{0,info.tangentScale/h.y};
            const Eigen::Vector3d p0(1,a.x,a.y),p1(0,t.x,t.y);
            faceGram[c]+=sign*(f[0]*p0*p0.transpose()+f[1]*(p0*p1.transpose()+p1*p0.transpose())+f[2]*p1*p1.transpose());
            b[0]=HydroCompatibleMomentSpace::affineVelocityBasis(a,h)*space.cells()[c].affineProjection;
            b[1]=(HydroCompatibleMomentSpace::affineVelocityBasis(a+t,h)*space.cells()[c].affineProjection-b[0]).eval();
        };
        if(face.outwardSign){
            using Kind=HydroCompatibleMomentBoundary::Kind;
            if(!boundary||HydroCompatibleMomentBoundary::impermeable(boundary->faces()[face.nativeRow].kind)){require(f.isZero(0),"phase momentum wall carries mass");continue;}
            require(next.parametricFluxes().hasExteriorDirectionMoments,"phase ports require actual signed/directional face moments");
            const auto& port=boundary->faces()[face.nativeRow];const int c=face.negativeCell>=0?face.negativeCell:face.positiveCell;
            require(port.incomingPhase==1||port.incomingPhase==-1,"two-phase port has no declared exterior phase");
            const auto& volume=InterfaceParametricMotion::faceVolumes(next,face.nativeRow);
            const double totalIn=face.outwardSign>0?volume.totalNegative:volume.totalPositive;
            const double liquidIn=face.outwardSign>0?volume.liquidNegative:volume.liquidPositive;
            const double capacity=domain.cellVolumes()[domain.cellRowByIndex()[c]];
            require(std::abs(liquidIn-(port.incomingPhase<0?totalIn:0.))<=1e-12*capacity,"incoming phase disagrees with the owned boundary geometry");
            require(totalIn==0||port.incomingVelocity,"actual C1 backflow requires explicit exterior momentum");
            Eigen::Vector2d exterior=Eigen::Vector2d::Zero();if(port.incomingVelocity)exterior<<port.incomingVelocity->x,port.incomingVelocity->y;
            const Eigen::Vector3d out=face.outwardSign>0?flux.exteriorPositive[face.nativeRow]:Eigen::Vector3d(flux.exteriorPositive[face.nativeRow]-f);
            const Eigen::Vector3d in=out-face.outwardSign*f;trace(c,pm,face.outwardSign);
            const auto& bulk=space.cells()[c];const int size=int(bulk.rows.size());
            Eigen::MatrixXd block=Eigen::MatrixXd::Zero(size,size),metric=block;Eigen::VectorXd supplied=Eigen::VectorXd::Zero(size);
            Eigen::Vector2d v[2]={pm[0]*values[c],pm[1]*values[c]},difference[2]={v[0]-exterior,v[1]},momentum=-in[0]*exterior/time.dt;
            for(int i=0;i<2;i++){
                supplied+=in[i]*pm[i].transpose()*exterior/time.dt;momentum+=out[i]*v[i]/time.dt;
                for(int j=0;j<2;j++){
                    block+=out[i+j]*pm[i].transpose()*pm[j]/time.dt;metric+=(out[i+j]+in[i+j])*pm[i].transpose()*pm[j];
                    result.boundaryKineticFlux+=.5*out[i+j]*v[i].dot(v[j])/time.dt;
                    result.boundaryInflowDissipation+=.5*in[i+j]*difference[i].dot(difference[j])/time.dt;
                }
            }
            result.boundaryKineticFlux-=.5*in[0]*exterior.squaredNorm()/time.dt;
            result.boundaryVolumeRate+=face.outwardSign*volume.total/time.dt;
            loads[c]+=block*values[c]-supplied;result.cellBoundaryMomentum[c]+=momentum;
            for(int i=0;i<size;i++)result.boundaryLoad[bulk.rows[i]]+=supplied[i];
            scatter(central,bulk.rows,bulk.rows,block);scatter(boundaryMetric,bulk.rows,bulk.rows,metric);continue;
        }
        const int cm=face.negativeCell,cp=face.positiveCell;const auto& bm=space.cells()[cm];const auto& bp=space.cells()[cp];
        const int nm=int(bm.rows.size()),np=int(bp.rows.size());
        trace(cm,pm,1);trace(cp,pp,-1);
        Eigen::MatrixXd jump[2],average[2];for(int k=0;k<2;k++){jump[k].resize(2,nm+np);average[k].resize(2,nm+np);jump[k]<<pm[k],-pp[k];average[k]<<.5*pm[k],.5*pp[k];}
        Eigen::MatrixXd block=Eigen::MatrixXd::Zero(nm+np,nm+np);
        Eigen::VectorXd both(nm+np);both<<values[cm],values[cp];
        Eigen::Vector2d boundary=Eigen::Vector2d::Zero();
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)block+=f[i+j]*jump[i].transpose()*average[j]/time.dt;
        for(int j=0;j<2;j++)boundary+=f[j]*average[j]*both/time.dt;
        const Eigen::VectorXd load=block*both;loads[cm]+=load.head(nm);loads[cp]+=load.tail(np);
        result.cellBoundaryMomentum[cm]+=boundary;result.cellBoundaryMomentum[cp]-=boundary;
        std::vector<int> rows=bm.rows;rows.insert(rows.end(),bp.rows.begin(),bp.rows.end());scatter(central,rows,rows,block);
    }
    for(int c=0;c<nc;c++){
        const auto& bulk=space.cells()[c];const auto h=mesh.cells(mesh.cellIndexToLabel[c]).size();
        const auto& g0=before.cells()[c].density;const auto& gm=midpoint.cells()[c].density;const auto& g1=after.cells()[c].density;
        const Eigen::Matrix3d delta=g1-g0,continuity=delta+faceGram[c];
        const Eigen::MatrixXd b0=HydroCompatibleMomentSpace::affineVelocityBasis({0,0},h)*bulk.affineProjection;
        const Eigen::MatrixXd bx=HydroCompatibleMomentSpace::affineVelocityBasis({1,0},h)*bulk.affineProjection-b0;
        const Eigen::MatrixXd by=HydroCompatibleMomentSpace::affineVelocityBasis({0,1},h)*bulk.affineProjection-b0;
        const Eigen::Vector2d v0=b0*values[c],vx=bx*values[c],vy=by*values[c];
        const Eigen::Vector3d ux(v0[0],vx[0],vy[0]),uy(v0[1],vx[1],vy[1]);
        // Weak continuity for X,Y,X^2,XY,Y^2 determines five bulk moments.
        // In RZ, g and f already contain the SAME 2*pi*r measure. Integration
        // by parts uses div(r*rho*u), so no extra 1/r source or separately
        // reconstructed density flux belongs in these moment identities.
        // Their one rotational complement comes from the SAME midpoint
        // material and velocity projection, not a second face flux.
        const double rotation=time.dt*(-gm.row(2).dot(ux)+gm.row(1).dot(uy));
        const double qxy=(continuity(1,2)-rotation/h.y)/(1/h.x+1/h.y),qyx=qxy+rotation;
        const Eigen::Vector3d qx(h.x*continuity(0,1),.5*h.x*continuity(1,1),qxy);
        const Eigen::Vector3d qy(h.y*continuity(0,2),qyx,.5*h.y*continuity(2,2));
        Eigen::MatrixXd local=-(bx/h.x).transpose()*(qx[0]*b0+qx[1]*bx+qx[2]*by)
            -(by/h.y).transpose()*(qy[0]*b0+qy[1]*bx+qy[2]*by);
        const Eigen::MatrixXd complement=Eigen::MatrixXd::Identity(bulk.rows.size(),bulk.rows.size())-bulk.affineSamples*bulk.affineProjection;
        // Discrete metric transport in the non-polynomial Hodge complement.
        // It annihilates affine velocity and is fixed by the mass product rule.
        local-=.5*delta(0,0)/before.cells()[c].volume(0,0)*complement.transpose()*bulk.mass*complement;
        local/=time.dt;loads[c]+=local*values[c];scatter(central,bulk.rows,bulk.rows,local);
        result.geometricMassClosure=std::max(result.geometricMassClosure,std::abs(continuity(0,0))/(g0(0,0)+g1(0,0)));
        const double l=std::sqrt(h.x*h.y);
        result.cellMomentum[c]<<h.x/l*bulk.affineSamples.col(0).dot(loads[c]),h.y/l*bulk.affineSamples.col(1).dot(loads[c]);
        result.maximumLocalMomentumDefect=std::max(result.maximumLocalMomentumDefect,(result.cellMomentum[c]-result.cellBoundaryMomentum[c]).norm());
    }
    result.central.resize(velocityCount_,velocityCount_);result.central.setFromTriplets(central.begin(),central.end());
    result.dissipation.resize(velocityCount_,velocityCount_);result.dissipation.setZero();
    Matrix outer(velocityCount_,velocityCount_);outer.setFromTriplets(boundaryMetric.begin(),boundaryMetric.end());
    const Matrix defect=time.dt*(result.central+Matrix(result.central.transpose()))+after.kineticMass()-before.kineticMass()-outer;
    result.massMetricDefect=defect.norm()/(before.kineticMass().norm()+after.kineticMass().norm());
    require(std::isfinite(result.massMetricDefect),"nonfinite moving-mass transport");
    return result;
}
HydroCompatibleMomentTransport::PhaseMassFlux HydroCompatibleMomentTransport::phaseMassFlux(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,
    const InterfaceState& previous,const InterfaceState& next,const InterfaceStepTime& time,double liquidDensity,double gasDensity) {
    space.validate(mesh,domain);InterfaceParametricMotion::validateDomain(mesh,domain,previous);InterfaceParametricMotion::validateDomain(mesh,domain,next);
    require(next.parametricFollows(previous,time),"moment phase flux source geometry/time mismatch");
    require(liquidDensity>0 && gasDensity>0 && std::isfinite(liquidDensity) && std::isfinite(gasDensity),"invalid moment phase densities");
    const auto& retained=next.parametricFluxes();require(retained.momentDegree==2,"moment momentum requires the actual first/second phase flux moments");
    const auto& cells=InterfaceParametricMotion::cellIndices(previous);
    std::vector<double> first(mesh.cellIndexNum()),last(mesh.cellIndexNum());
    for(int slot=0;slot<int(cells.size());slot++){
        first[cells[slot]]=previous.parametricGeometry().inventory()[slot];last[cells[slot]]=next.parametricGeometry().inventory()[slot];
    }
    return massFlux(mesh,domain,space,[&](int row)->const auto& {return InterfaceParametricMotion::faceVolumes(next,row);},
        first,last,liquidDensity,gasDensity,retained.hasVelocityProducts,retained.hasExteriorDirectionMoments);
}
HydroCompatibleMomentTransport::PhaseMassFlux HydroCompatibleMomentTransport::backgroundMassFlux(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,
    const std::vector<InterfaceParametricTransport::FaceVolumes::Face>& volumes,double density){
    require(volumes.size()==std::size_t(domain.nativeFaceCount()),"background mass flux changed native faces");
    for(const auto& face:volumes)require(face.liquid==0&&face.liquidPositive==0&&face.liquidNegative==0,
        "background mass flux contains a material interface");
    const std::vector<double> absent(mesh.cellIndexNum(),0.);
    return massFlux(mesh,domain,space,[&](int row)->const auto& {return volumes.at(row);},absent,absent,density,density,true,true);
}
HydroCompatibleMomentTransport::PhaseMassFlux HydroCompatibleMomentTransport::massFlux(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,
    const std::function<const InterfaceParametricTransport::FaceVolumes::Face&(int)>& faceVolumes,
    const std::vector<double>& previousLiquid,const std::vector<double>& nextLiquid,
    double liquidDensity,double gasDensity,bool velocityProducts,bool exteriorDirectionMoments){
    space.validate(mesh,domain);
    require(liquidDensity>0&&gasDensity>0&&std::isfinite(liquidDensity)&&std::isfinite(gasDensity),"invalid mass-flux density");
    require(previousLiquid.size()==std::size_t(mesh.cellIndexNum())&&nextLiquid.size()==previousLiquid.size(),"mass-flux inventory layout differs");
    PhaseMassFlux result;result.face.resize(domain.nativeFaceCount());result.exteriorPositive.assign(domain.nativeFaceCount(),Eigen::Vector3d::Zero());result.massBefore.resize(mesh.cellIndexNum());result.massAfter.resize(mesh.cellIndexNum());
    if(velocityProducts){result.products.resize(domain.nativeFaceCount());result.exteriorPositiveProducts.resize(domain.nativeFaceCount());}
    std::vector<long double> residual(mesh.cellIndexNum()),scale(mesh.cellIndexNum());
    for(int c=0;c<mesh.cellIndexNum();c++){
        const double capacity=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        const long double first=previousLiquid[c],last=nextLiquid[c];
        const auto before=liquidDensity*first+gasDensity*(capacity-first),after=liquidDensity*last+gasDensity*(capacity-last);
        result.massBefore[c]=static_cast<double>(before);result.massAfter[c]=static_cast<double>(after);
        residual[c]=after-before;scale[c]=before+after;
    }
    for(const auto f:domain.faces(mesh)){
        const auto& volumes=faceVolumes(f.nativeRow);
        if(velocityProducts){
            constexpr int degree=InterfaceParametricTransport::FaceVolumes::productDegree;auto& full=result.products[f.nativeRow];auto& positive=result.exteriorPositiveProducts[f.nativeRow];
            full.resize(degree+1);positive=Eigen::VectorXd::Zero(degree+1);
            for(int k=0;k<=degree;k++){
                full[k]=double(static_cast<long double>(liquidDensity)*volumes.liquidProducts[k]+static_cast<long double>(gasDensity)*(volumes.totalProducts[k]-volumes.liquidProducts[k]));
                if(f.outwardSign)positive[k]=double(static_cast<long double>(liquidDensity)*volumes.liquidPositiveProducts[k]+static_cast<long double>(gasDensity)*(volumes.totalPositiveProducts[k]-volumes.liquidPositiveProducts[k]));
            }
            require(full.allFinite()&&positive.allFinite(),"nonfinite complete mass-flux functional");
            const double lower=(f.horizontal?f.center.x:f.center.y)-.5*f.tangentSpacing;
            const auto constant=InterfaceParametricTransport::FaceVolumes::constantProduct(f.horizontal,domain.axisymmetric(),lower,lower+f.tangentSpacing);
            long double represented=0;for(int k=0;k<=degree;k++)represented+=constant[k]*static_cast<long double>(full[k]);
            const long double actual=static_cast<long double>(liquidDensity)*volumes.liquid+static_cast<long double>(gasDensity)*(volumes.total-volumes.liquid);
            double capacity=std::numeric_limits<double>::infinity();for(int cell:{f.negativeCell,f.positiveCell})if(cell>=0)capacity=std::min(capacity,domain.cellVolumes()[domain.cellRowByIndex()[cell]]);
            result.productConstantDefect=std::max(result.productConstantDefect,double(std::abs(represented-actual)/(std::max(liquidDensity,gasDensity)*capacity)));
        }
        // The very same Q_l and Q are combined at every order. No face
        // density, face velocity interpolation or reconstructed alpha enters.
        long double m[3];for(int k=0;k<3;k++){
            const long double total=k?volumes.totalHigher[k-1]:volumes.total,liquid=k?volumes.liquidHigher[k-1]:volumes.liquid;
            m[k]=liquidDensity*liquid+gasDensity*(total-liquid);
        }
        const auto& info=space.faces()[f.nativeRow];
        if(info.first<0){require(domain.axisymmetric()&&!f.horizontal&&f.center.x==0&&m[0]==0&&m[1]==0&&m[2]==0,"inactive RZ axis carries a phase mass flux");
            result.face[f.nativeRow].setZero();continue;}
        require(info.tangentScale>0,"invalid moment phase-flux face scale");
        const double lower=(f.horizontal?f.center.x:f.center.y)-.5*f.tangentSpacing;
        const long double a=(lower-info.tangentCenter)/info.tangentScale,b=f.tangentSpacing/info.tangentScale;
        result.face[f.nativeRow]<<double(m[0]),double(a*m[0]+b*m[1]),double(a*a*m[0]+2*a*b*m[1]+b*b*m[2]);
        if(f.outwardSign&&exteriorDirectionMoments){
            long double positive[3];for(int k=0;k<3;k++){
                const long double total=k?volumes.totalPositiveHigher[k-1]:volumes.totalPositive,liquid=k?volumes.liquidPositiveHigher[k-1]:volumes.liquidPositive;
                positive[k]=liquidDensity*liquid+gasDensity*(total-liquid);
            }
            result.exteriorPositive[f.nativeRow]<<double(positive[0]),double(a*positive[0]+b*positive[1]),double(a*a*positive[0]+2*a*b*positive[1]+b*b*positive[2]);
            require(result.exteriorPositive[f.nativeRow].allFinite(),"nonfinite exterior directional mass moment");
        }
        require(result.face[f.nativeRow].allFinite(),"nonfinite phase mass-flux moment");
        if(f.negativeCell>=0){residual[f.negativeCell]+=m[0];scale[f.negativeCell]+=std::abs(m[0]);}
        if(f.positiveCell>=0){residual[f.positiveCell]-=m[0];scale[f.positiveCell]+=std::abs(m[0]);}
    }
    for(int c=0;c<mesh.cellIndexNum();c++){require(scale[c]>0,"invalid moment inventory mass scale");const auto error=std::abs(residual[c])/scale[c];
        result.maximumRelativeClosure=std::max(result.maximumRelativeClosure,double(error));
        require(error<=4096*std::numeric_limits<double>::epsilon(),"shared phase flux does not close the actual moment cell mass");}
    return result;
}
}
