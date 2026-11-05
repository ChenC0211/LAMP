#include "Modules/Electric/ElectricVofField.h"
#include "Modules/Electric/ElectricVofMaxwellStress.h"
#include "Modules/Electric/ElectricVofChargeTransport.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <tuple>
#include <cfenv>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include <cerrno>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Modules/Electric/ElectricPotential.h"
#include "Modules/Electric/ElectricMaxwellTrace.h"
#include "Modules/Electric/ElectricPhaseGeometry.h"
#include "Modules/Electric/SurfaceChargeField.h"
#include <optional>
#include <utility>
#include <cstring>
#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Modules/Interface/HermiteGalsFirstJetStep.h"
#include <set>
#include "Solvers/Poisson/ElectricLocalOperator.h"
#include "Solvers/Poisson/PolarizationPoissonFVSolverCellCentered.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <algorithm>
#include <string>
#include <sstream>
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Field/CellField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/VertField.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <type_traits>
#include <typeinfo>
#include <limits>
#include <cstdlib>

namespace ElectricVofOhmicQuantumContract { void verify(const QuadtreeMesh&,const LAMP::Electric::ElectricVofField&,bool,bool,bool,bool); }
namespace ElectricVofPostprocessContract { void executionReceipt(const LAMP::Electric::ElectricVofField&,unsigned); }
namespace {
using namespace LAMP::Electric;
void check(bool ok, const char* reason) { if (!ok) throw std::runtime_error(reason); }
template<class F> void rejects(F action) {
    try { action(); } catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("invalid electric source accepted");
}

namespace ElectricContract {
namespace H = LAMP::ElectricDiscretization;
constexpr double pi=3.1415926535897932384626433832795;
constexpr int powers[10][2]={{0,0},{1,0},{0,1},{2,0},{1,1},{0,2},{3,0},{2,1},{1,2},{0,3}};
constexpr double nodes[8]={-.9602898564975363,-.7966664774136267,-.5255324099163290,-.1834346424956498,
    .1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363};
constexpr double weights[8]={.1012285362903763,.2223810344533745,.3137066458778873,.3626837833783620,
    .3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763};
using Coefficients=std::array<double,10>;
struct Jet { double v=0,x=0,y=0,xx=0,xy=0,yy=0; };
double monomial(double x,int order) { return order==0?1.:std::pow(x,order); }
Jet polynomial(Vector2D point,Vector2D center,double scale,const Coefficients& coefficients) {
    const double x=(point.x-center.x)/scale,y=(point.y-center.y)/scale;
    Jet result;
    for(int k=0;k<10;++k) {
        const int a=powers[k][0],b=powers[k][1];const double c=coefficients[k];
        result.v+=c*monomial(x,a)*monomial(y,b);
        if(a)result.x+=c*a*monomial(x,a-1)*monomial(y,b)/scale;
        if(b)result.y+=c*b*monomial(x,a)*monomial(y,b-1)/scale;
        if(a>1)result.xx+=c*a*(a-1)*monomial(x,a-2)*monomial(y,b)/(scale*scale);
        if(a&&b)result.xy+=c*a*b*monomial(x,a-1)*monomial(y,b-1)/(scale*scale);
        if(b>1)result.yy+=c*b*(b-1)*monomial(x,a)*monomial(y,b-2)/(scale*scale);
    }
    return result;
}
template<class F> void rejected(F action,const char* reason) {
    bool caught=false;try {action();}catch(const std::exception&){caught=true;}
    check(caught,reason);
}
void close(double actual,double expected,double tolerance,const char* reason) {
    check(std::abs(actual-expected)<=tolerance*(1+std::abs(expected)),reason);
}
template<class F> double areaIntegral(Vector2D lo,Vector2D hi,F function) {
    const auto h=hi-lo,center=lo+h*.5;long double sum=0;
    for(int y=0;y<8;++y)for(int x=0;x<8;++x)
        sum+=weights[x]*weights[y]*function(center+Vector2D{.5*h.x*nodes[x],.5*h.y*nodes[y]});
    return static_cast<double>(sum*.25*h.x*h.y);
}
template<class F> double sideIntegral(Vector2D lo,Vector2D hi,int side,F function) {
    const auto h=hi-lo;const bool vertical=side==0||side==2;
    const Vector2D center=vertical?Vector2D{side==0?lo.x:hi.x,.5*(lo.y+hi.y)}
        :Vector2D{.5*(lo.x+hi.x),side==1?lo.y:hi.y};
    const double length=vertical?h.y:h.x;long double sum=0;
    for(int k=0;k<8;++k)sum+=weights[k]*function(center+(vertical?Vector2D{0,.5*length*nodes[k]}:
        Vector2D{.5*length*nodes[k],0}),nodes[k]);
    return static_cast<double>(sum*.5*length);
}
Vector2D normal(int side) {
    const Vector2D directions[4]={{-1,0},{0,-1},{1,0},{0,1}};return directions[side];
}
H::LocalGeometry rectangle(Vector2D lo,Vector2D hi,bool axisymmetric) {
    H::LocalGeometry result;const auto h=hi-lo;result.center=lo+h*.5;
    result.scale=std::max(h.x,h.y);result.axisymmetric=axisymmetric;
    for(int y=0;y<8;++y)for(int x=0;x<8;++x)result.volume.push_back({
        result.center+Vector2D{.5*h.x*nodes[x],.5*h.y*nodes[y]},.25*h.x*h.y*weights[x]*weights[y]});
    for(int side=0;side<4;++side) {
        if(axisymmetric&&side==0&&lo.x==0)continue;
        H::FaceQuadrature face;face.length=side==0||side==2?h.y:h.x;
        sideIntegral(lo,hi,side,[&](Vector2D p,double zeta) {
            // This callback is visited in the same eight-point rule; its
            // return does not participate in construction of the input data.
            int k=0;while(k<8&&nodes[k]!=zeta)++k;
            check(k<8,"independent rectangle face node not found");
            face.nodes.push_back({p,normal(side),.5*face.length*weights[k],zeta});return 0.;
        });
        result.faces.push_back(std::move(face));
    }
    return result;
}
double moment(int exponent,double halfExtent) {
    return exponent%2?0.:2*std::pow(halfExtent,exponent+1)/(exponent+1);
}
H::CellVector ordinaryProjection(Vector2D lo,Vector2D hi,double scale,const Coefficients& c) {
    Eigen::Matrix<double,6,6> gram;H::CellVector rhs;const auto h=hi-lo;
    const auto exact=[&](int a,int b){return scale*scale*moment(a,h.x/(2*scale))*moment(b,h.y/(2*scale));};
    for(int i=0;i<6;++i) {
        rhs[i]=0;
        for(int j=0;j<6;++j)gram(i,j)=exact(powers[i][0]+powers[j][0],powers[i][1]+powers[j][1]);
        for(int j=0;j<10;++j)rhs[i]+=c[j]*exact(powers[i][0]+powers[j][0],powers[i][1]+powers[j][1]);
    }
    return gram.ldlt().solve(rhs);
}
double relativeSymmetry(const Eigen::MatrixXd& matrix) {
    return (matrix-matrix.transpose()).norm()/std::max(matrix.norm(),std::numeric_limits<double>::min());
}
void positiveSemidefinite(const Eigen::MatrixXd& matrix,const char* reason) {
    check(relativeSymmetry(matrix)<=1e-12,"electric energy block is not symmetric");
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(.5*(matrix+matrix.transpose()));
    check(eig.info()==Eigen::Success,reason);
    const double scale=std::max(matrix.norm(),std::numeric_limits<double>::min());
    check(eig.eigenvalues().minCoeff()>=-2e-12*scale,reason);
}
void local() {
    double maximumJet=0,maximumCellProjection=0,maximumPde=0,maximumRecovery=0;
    for(int configuration=0;configuration<3;++configuration) {
        const bool axis=configuration!=0;const Vector2D lo=configuration==2?Vector2D{0,-.31}:Vector2D{.27,-.31};
        const Vector2D hi=configuration==2?Vector2D{.86,.83}:Vector2D{1.13,.83};
        const auto geometry=rectangle(lo,hi,axis);constexpr double epsilon=2.3;
        H::LocalOperator op(geometry,epsilon);
        check(op.faceCount()==std::size_t(configuration==2?3:4),"axis face was not omitted from local unknowns");
        const auto& K=op.stiffnessMatrix();const auto& S=op.condensedMatrix();
        check(K.rows()==6+3*int(op.faceCount())&&S.rows()==3*int(op.faceCount()),"local cell-only condensation layout");
        positiveSemidefinite(K,"local K is not PSD");positiveSemidefinite(S,"local skeleton S is not PSD");
        positiveSemidefinite(op.stabilizationMatrix(),"local stabilization is not PSD");
        for(double contrast:{1e-6,1e6}) {
            H::LocalOperator scaled(geometry,contrast);const double ratio=contrast/epsilon;
            check((scaled.reconstructionMatrix()-op.reconstructionMatrix()).norm()
                <=1e-11*(1+op.reconstructionMatrix().norm()),"constant epsilon changed polynomial reconstruction");
            for(const auto pair:{std::make_pair(&scaled.stiffnessMatrix(),&K),
                    std::make_pair(&scaled.condensedMatrix(),&S)}) {
                check((*pair.first-ratio*(*pair.second)).norm()<=1e-11*ratio*pair.second->norm(),
                    "local energy or Schur block does not scale with epsilon");
            }
            check((scaled.gradientGramMatrix()-ratio*op.gradientGramMatrix()).norm()
                <=1e-11*ratio*op.gradientGramMatrix().norm(),"gradient Gram does not scale with epsilon");
        }
        Eigen::VectorXd arbitrary(K.rows());
        for(int i=0;i<arbitrary.size();++i)arbitrary[i]=std::sin(.73*i)+.1*i;
        const Eigen::VectorXd fullLoad=K*arbitrary,arbitraryFaces=arbitrary.tail(S.rows());
        const H::CellVector originalCell=arbitrary.head<6>(),arbitraryLoad=fullLoad.head<6>();
        const auto originalRecovered=op.recoverCell(arbitraryFaces,arbitraryLoad);
        check((originalRecovered-originalCell).norm()<=1e-11*(1+originalCell.norm()),
            "cell static recovery failed the full arbitrary local equation");
        const auto arbitraryCondensed=op.condensedRhs(arbitraryLoad,Eigen::VectorXd::Zero(S.rows()));
        check((S*arbitraryFaces-arbitraryCondensed-fullLoad.tail(S.rows())).norm()
            <=1e-11*(1+fullLoad.norm()),"arbitrary full local residual differs from Schur");
        const auto arbitraryReconstruction=op.reconstruct(originalCell,arbitraryFaces);
        Coefficients r{},t{};
        for(int k=0;k<10;++k)r[k]=arbitraryReconstruction[k];
        for(int k=0;k<6;++k)t[k]=originalCell[k];
        const double meanDefect=areaIntegral(lo,hi,[&](Vector2D point) {
            return (axis?2*pi*point.x:1)*(polynomial(point,geometry.center,geometry.scale,r).v
                -polynomial(point,geometry.center,geometry.scale,t).v);
        });
        close(meanDefect,0,1e-11,"reconstruction lost its physical weighted mean constraint");
        const auto constant=op.reduction([](Vector2D){return 1.;});
        check((K*constant).norm()<=1e-11*K.norm()*constant.norm(),"local energy lost the constant kernel");
        check((S*constant.tail(S.rows())).norm()<=1e-11*S.norm()*constant.tail(S.rows()).norm(),"Schur lost the constant kernel");
        for(int column=0;column<10;++column) {
            Coefficients c{};c[column]=1;
            const auto analytic=[&](Vector2D p){return polynomial(p,geometry.center,geometry.scale,c);};
            const auto reduced=op.reduction([&](Vector2D p){return analytic(p).v;});
            const auto exactCell=ordinaryProjection(lo,hi,geometry.scale,c);
            const double projection=(reduced.head<6>()-exactCell).norm()/(1+exactCell.norm());
            maximumCellProjection=std::max(maximumCellProjection,projection);
            check(projection<=1e-11,"cell projector used the RZ measure instead of ordinary dA");
            const H::PotentialVector reconstructed=op.reconstruct(reduced.head<6>(),reduced.tail(reduced.size()-6));
            Coefficients represented{};for(int k=0;k<10;++k)represented[k]=reconstructed[k];
            for(double x:{.07,.39,.91})for(double y:{.13,.61,.87}) {
                const auto p=lo+Vector2D{x*(hi.x-lo.x),y*(hi.y-lo.y)};
                const auto exact=analytic(p),actual=polynomial(p,geometry.center,geometry.scale,represented);
                const double h=geometry.scale,error=std::max({std::abs(actual.v-exact.v),
                    h*std::abs(actual.x-exact.x),h*std::abs(actual.y-exact.y),
                    h*h*std::abs(actual.xx-exact.xx),h*h*std::abs(actual.xy-exact.xy),h*h*std::abs(actual.yy-exact.yy)});
                maximumJet=std::max(maximumJet,error);check(error<=1e-11,"weighted P3 polynomial reconstruction failed");
            }
        }
        const Coefficients c{{1.7,.23,-.31,.19,.07,-.11,.13,-.17,.09,.21}};
        const auto analytic=[&](Vector2D p){return polynomial(p,geometry.center,geometry.scale,c);};
        const auto reduced=op.reduction([&](Vector2D p){return analytic(p).v;});
        Eigen::VectorXd rhs=Eigen::VectorXd::Zero(reduced.size());
        for(int mode=0;mode<6;++mode)rhs[mode]=areaIntegral(lo,hi,[&](Vector2D p) {
            const auto j=analytic(p);const double x=(p.x-geometry.center.x)/geometry.scale,y=(p.y-geometry.center.y)/geometry.scale;
            const double weightedSource=-epsilon*(axis?2*pi*(p.x*(j.xx+j.yy)+j.x):j.xx+j.yy);
            return weightedSource*monomial(x,powers[mode][0])*monomial(y,powers[mode][1]);
        });
        int ordinal=0;
        for(int side=0;side<4;++side) {
            if(axis&&side==0&&lo.x==0)continue;
            for(int mode=0;mode<3;++mode)rhs[6+3*ordinal+mode]=sideIntegral(lo,hi,side,[&](Vector2D p,double zeta) {
                const auto j=analytic(p);const auto n=normal(side);
                const double basis=mode==0?1.:mode==1?zeta:.5*(3*zeta*zeta-1);
                return epsilon*(j.x*n.x+j.y*n.y)*(axis?2*pi*p.x:1)*basis;
            });
            ++ordinal;
        }
        const double pde=(K*reduced-rhs).norm()/(1+rhs.norm());maximumPde=std::max(maximumPde,pde);
        check(pde<=1e-11,"local energy differs from independent cubic PDE/face moments");
        const Eigen::VectorXd faces=reduced.tail(reduced.size()-6);const H::CellVector load=rhs.head<6>();
        const auto recovered=op.recoverCell(faces,load);
        const double recovery=(recovered-reduced.head<6>()).norm()/(1+reduced.head<6>().norm());
        maximumRecovery=std::max(maximumRecovery,recovery);check(recovery<=1e-11,"static recovery differs from original cell equations");
        const auto condensed=op.condensedRhs(load,Eigen::VectorXd::Zero(faces.size()));
        check((S*faces-condensed-rhs.tail(faces.size())).norm()<=1e-11*(1+rhs.norm()),"Schur residual differs from full face equations");
        const auto residual=op.faceResidual(recovered,faces,Eigen::VectorXd::Zero(faces.size()));
        check((residual-rhs.tail(faces.size())).norm()<=1e-11*(1+rhs.norm()),"face residual differs from independent conormal moments");
        double flux=0;for(int face=0;face<int(op.faceCount());++face)flux+=residual[3*face];
        close(flux+rhs[0],0,1e-11,"local Gauss balance failed");
        auto mutatedGeometry=geometry;mutatedGeometry.volume.front().weight*=2;
        close(op.geometry().volume.front().weight,geometry.volume.front().weight,0,"local operator borrowed mutable quadrature");
        rejected([&]{H::LocalOperator bad(geometry,0.);},"nonpositive epsilon accepted");
        rejected([&]{H::LocalOperator bad(geometry,epsilon,-1.);},"negative stabilization accepted");
    }
    std::cout<<"ELECTRIC_HHO_P0_LOCAL,polynomial_columns=30,jet="<<maximumJet<<",ordinary_projection="<<maximumCellProjection
        <<",independent_pde="<<maximumPde<<",recover="<<maximumRecovery<<",PASS\n";
}
void localFrame() {
    const auto geometry=rectangle({-.125,.25},{.875,1.5},false);
    constexpr double epsilon=2.3;
    H::LocalOperator original(geometry,epsilon);
    auto shifted=geometry;shifted.coordinateOrigin={std::ldexp(1.,40),-std::ldexp(1.,39)};
    H::LocalOperator framed(shifted,epsilon);
    const auto sameBits=[](const auto& a,const auto& b,const char* reason) {
        check(a.rows()==b.rows()&&a.cols()==b.cols(),reason);
        check(std::memcmp(a.data(),b.data(),sizeof(double)*static_cast<std::size_t>(a.size()))==0,reason);
    };
    sameBits(framed.reconstructionMatrix(),original.reconstructionMatrix(),"XY frame altered reconstruction bits");
    sameBits(framed.stiffnessMatrix(),original.stiffnessMatrix(),"XY frame altered stiffness bits");
    sameBits(framed.condensedMatrix(),original.condensedMatrix(),"XY frame altered Schur bits");
    sameBits(framed.gradientGramMatrix(),original.gradientGramMatrix(),"XY frame altered gradient Gram bits");
    sameBits(framed.stabilizationMatrix(),original.stabilizationMatrix(),"XY frame altered stabilization bits");
    const auto origin=shifted.coordinateOrigin;
    const auto worldValue=[&](Vector2D point) {
        return 1.25+.5*(point.x-origin.x)-.75*(point.y-origin.y);
    };
    const auto callbacks=[&](const std::vector<Vector2D>& points,auto actual,auto reference) {
        std::size_t visited=0;
        const auto checked=[&](Vector2D point) {
            check(visited<points.size(),"framed callback visited an extra node");
            const auto expected=origin+points[visited++];
            check(point.x==expected.x&&point.y==expected.y,"framed callback did not receive actual world point");
            return worldValue(point);
        };
        const auto actualResult=actual(checked);
        const auto expectedResult=reference([&](Vector2D point){return worldValue(origin+point);});
        check(visited==points.size(),"framed callback missed actual quadrature nodes");
        sameBits(actualResult,expectedResult,"world callback changed local reduction/load arithmetic");
    };
    std::vector<Vector2D> volumePoints,facePoints;
    for(const auto& node:geometry.volume)volumePoints.push_back(node.point);
    for(const auto& face:geometry.faces)for(const auto& node:face.nodes)facePoints.push_back(node.point);
    callbacks(volumePoints,[&](const H::ScalarFunction& f){return framed.reductionCell(f);},
        [&](const H::ScalarFunction& f){return original.reductionCell(f);});
    callbacks(facePoints,[&](const H::ScalarFunction& f){return framed.reductionFaces(f);},
        [&](const H::ScalarFunction& f){return original.reductionFaces(f);});
    callbacks(volumePoints,[&](const H::ScalarFunction& f){return framed.loadCell(f);},
        [&](const H::ScalarFunction& f){return original.loadCell(f);});
    for(std::size_t face=0;face<geometry.faces.size();++face){
        std::vector<Vector2D> points;for(const auto& node:geometry.faces[face].nodes)points.push_back(node.point);
        callbacks(points,[&](const H::ScalarFunction& f){return framed.loadFace(face,f);},
            [&](const H::ScalarFunction& f){return original.loadFace(face,f);});
    }
    const auto rzGeometry=rectangle({0,-.25},{.875,1.},true);
    H::LocalOperator rz(rzGeometry,epsilon);
    auto zeroOrigin=rzGeometry;zeroOrigin.coordinateOrigin={0,0};
    H::LocalOperator rzZero(zeroOrigin,epsilon);
    sameBits(rz.stiffnessMatrix(),rzZero.stiffnessMatrix(),"RZ zero origin changed stiffness");
    std::size_t rzVisited=0;
    rzZero.loadCell([&](Vector2D point){
        check(rzVisited<rzGeometry.volume.size(),"RZ zero-origin callback visited an extra node");
        const auto expected=rzGeometry.volume[rzVisited++].point;
        check(std::memcmp(&point.x,&expected.x,sizeof(double))==0
            &&std::memcmp(&point.y,&expected.y,sizeof(double))==0,"RZ zero-origin callback changed point bits");
        return 1.;
    });
    check(rzVisited==rzGeometry.volume.size(),"RZ zero-origin callback missed a node");
    for(const auto originValue:{Vector2D{.5,0},Vector2D{0,.5}}){
        auto invalid=rzGeometry;invalid.coordinateOrigin=originValue;
        rejected([&]{H::LocalOperator bad(invalid,epsilon);},"RZ nonzero coordinate origin accepted");
    }
    auto invalid=geometry;invalid.coordinateOrigin.x=std::numeric_limits<double>::infinity();
    rejected([&]{H::LocalOperator bad(invalid,epsilon);},"nonfinite coordinate origin accepted");
    check(geometry.coordinateOrigin.x==0&&geometry.coordinateOrigin.y==0,"framed operator changed caller geometry");
    std::cout<<"ELECTRIC_HHO_FRAME,XY_matrix_bits=PASS,world_callbacks=PASS,RZ_zero_origin=PASS\n";
}
void globalCase(int kind) {
    const bool axis=kind>=3;constexpr double epsilon=2.3;
    DictionaryMesh input;input.Nx=4;input.Ny=3;input.maxInitialIterations=0;input.maxCellNum=128;
    QuadtreeMesh mesh(input);UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    Coefficients c{};
    if(kind==0)c[0]=2.5;
    else if(kind==1||kind==3){c[0]=1.2;c[2]=.7;if(!axis)c[1]=.4;}
    else if(axis){c[0]=1;c[2]=1;c[9]=.15;c[7]=.2;}
    else c={{1.3,.4,-.2,.3,.15,.2,.07,-.08,.09,-.05}};
    const auto analytic=[c](Vector2D p){return polynomial(p,{0,0},1,c);};
    ElectricCartesianProblem problem;problem.axisymmetric=axis;problem.epsilonByCell.assign(mesh.cellIndexNum(),epsilon);
    problem.source=[analytic,axis,epsilon](Vector2D p){const auto j=analytic(p);
        // Axis MMS is even in r; its analytic limit is finite, including r=0.
        return -epsilon*(j.xx+j.yy+(axis?(p.x==0?j.xx:j.x/p.x):0));};
    const bool mixed=kind==1||kind==3;
    for(int side=0;side<4;++side) {
        const bool neumann=(axis&&side==0)||(mixed&&side!=1);
        problem.boundaries[side].type=neumann?ElectricCartesianBoundary::Type::Neumann:ElectricCartesianBoundary::Type::Dirichlet;
        if(neumann)problem.boundaries[side].value=[analytic,side,epsilon](Vector2D p){const auto j=analytic(p);const auto n=normal(side);
            return epsilon*(j.x*n.x+j.y*n.y);};
        else problem.boundaries[side].value=[analytic](Vector2D p){return analytic(p).v;};
    }
    const auto result=solveElectricPotential(mesh,problem);result.validateMesh(mesh);
    const auto& report=result.report();check(result.nativeCellCount()==12&&report.cells==12,"global electric did not use actual native cells");
    check(report.axisFaces==std::size_t(axis?3:0),"global axis face inventory is incorrect");
    check(report.relativeLinearResidual<=1e-11&&report.relativeSymmetryDefect<=1e-12
        &&report.minimumSkeletonPivot>0,"global skeleton failed symmetry/SPD/linear residual");
    const Eigen::MatrixXd skeleton(result.skeletonMatrix());
    check(relativeSymmetry(skeleton)<=1e-12,"owned skeleton is not symmetric");
    Eigen::Map<const Eigen::VectorXd> solution(result.skeletonSolution().data(),result.skeletonSolution().size());
    Eigen::Map<const Eigen::VectorXd> rhs(result.skeletonRhs().data(),result.skeletonRhs().size());
    check((skeleton*solution-rhs).norm()<=1e-11*(1+rhs.norm()),"owned skeleton field does not solve its equations");
    double maximumPotential=0,maximumGradient=0,maximumFlux=0,maximumInterior=0;double boundaryFlux=0,totalSource=0;
    for(int y=0;y<grid.ny();++y)for(int x=0;x<grid.nx();++x) {
        const int native=mesh.cells(grid.cells()[y*grid.nx()+x]).index();
        const Vector2D lo=grid.minCoord()+Vector2D{x*grid.spacing().x,y*grid.spacing().y},hi=lo+grid.spacing();
        for(double sx:{.17,.53,.89})for(double sy:{.19,.57,.83}) {
            const auto p=lo+Vector2D{sx*grid.spacing().x,sy*grid.spacing().y};const auto expected=analytic(p);
            const auto field=result.gradientInCell(native,p);
            maximumPotential=std::max(maximumPotential,std::abs(result.potentialInCell(native,p)-expected.v)/(1+std::abs(expected.v)));
            maximumGradient=std::max({maximumGradient,std::abs(field.x-expected.x)/(1+std::abs(expected.x)),
                std::abs(field.y-expected.y)/(1+std::abs(expected.y))});
        }
        const double source=areaIntegral(lo,hi,[&](Vector2D p){const auto j=analytic(p);
            return -epsilon*(axis?2*pi*(p.x*(j.xx+j.yy)+j.x):j.xx+j.yy);});
        totalSource+=source;double outward=0;
        for(int side=0;side<4;++side) {
            if(axis&&side==0&&x==0) {
                rejected([&]{result.outwardFaceResidual(native,side);},"axis acquired a face residual DOF");continue;
            }
            const bool boundary=side==0?x==0:side==1?y==0:side==2?x==grid.nx()-1:y==grid.ny()-1;
            const auto residual=result.outwardFaceResidual(native,side);
            const double conormal=sideIntegral(lo,hi,side,[&](Vector2D p,double){const auto j=analytic(p);const auto n=normal(side);
                return epsilon*(j.x*n.x+j.y*n.y)*(axis?2*pi*p.x:1);});
            const double prescribed=boundary&&problem.boundaries[side].type==ElectricCartesianBoundary::Type::Neumann?conormal:0;
            const double raw=residual[0]+prescribed;outward+=raw;
            maximumFlux=std::max(maximumFlux,std::abs(raw-conormal)/(1+std::abs(conormal)));
            if(boundary)boundaryFlux+=raw;
            else if(side==2||side==3) {
                const int other=mesh.cells(grid.cells()[(y+(side==3))*grid.nx()+x+(side==2)]).index();
                const auto opposite=result.outwardFaceResidual(other,side==2?0:1);
                for(int mode=0;mode<3;++mode)maximumInterior=std::max(maximumInterior,
                    std::abs(residual[mode]+opposite[mode])/(1+std::abs(residual[mode])+std::abs(opposite[mode])));
            }
        }
        close(outward+source,0,1e-10,"global cell Gauss equation failed independent volume/face integral");
    }
    const double tolerance=(kind==2||kind==4)?1e-9:1e-10;
    check(maximumPotential<=tolerance&&maximumGradient<=tolerance&&maximumFlux<=tolerance
        &&maximumInterior<=1e-10,"actual shared electric polynomial MMS/face balance failed");
    close(boundaryFlux+totalSource,0,1e-10,"global independent boundary Gauss balance failed");
    for(int family=0;family<2;++family) {
        const bool horizontal=family==1;const auto& labels=horizontal?grid.horizontal():grid.vertical();
        for(int label:labels) {
            const int edge=horizontal?mesh.hEdges(label).index():mesh.vEdges(label).index();
            if(axis&&!horizontal&&mesh.vEdges(label).midpointCoord().x==0) {
                rejected([&]{result.nativeFaceTrace(false,edge);},"axis acquired a shared trace DOF");continue;
            }
            const auto box=result.nativeFaceBounds(horizontal,edge);const auto trace=result.nativeFaceTrace(horizontal,edge);
            const auto center=(box[0]+box[1])*.5,h=box[1]-box[0];
            for(int mode=0;mode<3;++mode) {
                long double defect=0;
                for(int k=0;k<8;++k) {
                    const auto p=center+h*(.5*nodes[k]);const double legendre=.5*(3*nodes[k]*nodes[k]-1);
                    const double represented=trace[0]+trace[1]*nodes[k]+trace[2]*legendre;
                    const double basis=mode==0?1.:mode==1?nodes[k]:legendre;
                    defect+=weights[k]*(axis?2*pi*p.x:1)*(represented-analytic(p).v)*basis;
                }
                check(std::abs(defect)<=2e-9,"shared physical face trace is not the weighted P2 projection");
            }
        }
    }
    const auto retained=result;const int first=mesh.cells(grid.cells()[0]).index();
    const Vector2D inFirst{.1,.1};const double saved=retained.potentialInCell(first,inFirst);
    problem.epsilonByCell[0]=0;problem.source=[](Vector2D){return 100.;};
    close(retained.potentialInCell(first,inFirst),saved,0,"solved result borrowed mutable input data");
    auto invalid=problem;invalid.epsilonByCell.assign(mesh.cellIndexNum(),epsilon);
    for(auto& boundary:invalid.boundaries){boundary.type=ElectricCartesianBoundary::Type::Neumann;boundary.value=[](Vector2D){return 0.;};}
    invalid.source={};
    rejected([&]{solveElectricPotential(mesh,invalid);},"unanchored electric skeleton accepted");
    invalid=problem;rejected([&]{solveElectricPotential(mesh,invalid);},"nonpositive global epsilon accepted");
    if(axis) {
        invalid=problem;invalid.epsilonByCell.assign(mesh.cellIndexNum(),epsilon);
        invalid.boundaries[0].type=ElectricCartesianBoundary::Type::Dirichlet;
        rejected([&]{solveElectricPotential(mesh,invalid);},"axis Dirichlet trace accepted");
        invalid.boundaries[0].type=ElectricCartesianBoundary::Type::Neumann;
        invalid.boundaries[0].value=[](Vector2D){return 1.;};
        rejected([&]{solveElectricPotential(mesh,invalid);},"nonzero axis flux accepted");
    }
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]);vertex.setCoord(vertex.coord()+Vector2D{.01,0});
    rejected([&]{retained.validateMesh(mesh);},"stale mesh identity accepted by electric owner");
    close(retained.potentialInCell(first,inFirst),saved,0,"mesh mutation changed immutable electric polynomial");
    std::cout<<"ELECTRIC_HHO_P0_GLOBAL,"<<kind<<",cells=12,potential="<<maximumPotential<<",gradient="<<maximumGradient
        <<",flux="<<maximumFlux<<",interior="<<maximumInterior<<",gauss="<<boundaryFlux+totalSource<<",PASS\n";
}
void run() {localFrame();local();for(int kind=0;kind<5;++kind)globalCase(kind);}
} // namespace ElectricContract


namespace ElectricPhaseContract {
namespace I=LAMP::Interface;
namespace H=LAMP::ElectricDiscretization;
using ElectricContract::pi;
using ElectricContract::rejected;
constexpr double radius=.213,epsilonNegative=5.,epsilonPositive=1.;
Vector2D origin(bool axis) { return {axis?0.:.5,.5}; }
double scalar(Vector2D p,bool axis) {
    const auto q=p-origin(axis);return q.x*q.x+q.y*q.y-radius*radius;
}
I::HermiteGalsPublication publication(const QuadtreeMesh& mesh,bool axis,I::HermiteGalsTimeLayer time={}) {
    I::HermiteGalsPrimaryState primary;primary.transport=I::HermiteGalsTransport::FirstJet;
    primary.measure=axis?I::HermiteGalsMeasure::AxisymmetricVolume:I::HermiteGalsMeasure::CartesianArea;
    for(int v=0;v<mesh.vertIndexNum();++v) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord(),q=p-origin(axis);
        primary.vertices.push_back({v,scalar(p,axis),2*q.x,2*q.y});
    }
    std::set<std::array<int,2>> keys;
    std::vector<double> centers;
    for(int c=0;c<mesh.cellIndexNum();++c) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);primary.geometryCells.push_back(c);
        centers.push_back(scalar(cell.centerCoord(),axis));
        for(int side=0;side<4;++side) {
            std::array<int,2> key{{mesh.verts(cell.cellToVert[side]).index(),mesh.verts(cell.cellToVert[(side+1)%4]).index()}};
            const auto a=mesh.verts(mesh.vertIndexToLabel[key[0]]).coord(),b=mesh.verts(mesh.vertIndexToLabel[key[1]]).coord();
            if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);
            keys.insert(key);
        }
    }
    for(const auto& key:keys) {
        const auto a=mesh.verts(mesh.vertIndexToLabel[key[0]]).coord(),b=mesh.verts(mesh.vertIndexToLabel[key[1]]).coord();
        const auto p=(a+b)*.5,q=p-origin(axis);
        primary.edges.push_back({key,scalar(p,axis),2*(a.y==b.y?q.y:q.x)});
    }
    I::HermiteGalsCandidate candidate(mesh,time,primary);
    // AxisymmetricVolume already supplies the r=0 symmetry identity.
    // Additional declared planes must not register that axis a second time.
    const std::vector<LevelSetInterfacePanels::SymmetryPlane> symmetry;
    const auto mixed=I::closeHermiteGalsMixedAmr(mesh,candidate,symmetry);
    I::HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=symmetry;
    const I::HermiteGalsCenterPrediction prediction{candidate.identity(),centers};
    return I::publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,prediction,config);
}
struct Reference {
    double value=0;
    Vector2D gradient;
};
Reference exact(Vector2D p,bool axis,bool liquid,int kind) {
    if(kind==0)return {2.5,{0,0}};
    const auto q=p-origin(axis);const double squared=q.x*q.x+q.y*q.y;
    const double epsilon=liquid?epsilonNegative:epsilonPositive;
    if(kind==1) {
        const double phi=squared-radius*radius;
        if(axis)return {phi*q.y/epsilon,{2*q.x*q.y/epsilon,(phi+2*q.y*q.y)/epsilon}};
        return {phi*q.x/epsilon,{(phi+2*q.x*q.x)/epsilon,2*q.x*q.y/epsilon}};
    }
    check(!axis,"only the declared XY analytic dipole is part of the P1 short contract");
    const double C=radius*radius*(epsilonPositive-epsilonNegative)/(epsilonNegative+epsilonPositive);
    const double A=2*epsilonPositive/(epsilonNegative+epsilonPositive);
    if(liquid)return {A*q.x,{A,0}};
    return {q.x+C*q.x/squared,
        {1+C*(squared-2*q.x*q.x)/(squared*squared),-2*C*q.x*q.y/(squared*squared)}};
}
ElectricPhaseProblem problem(bool axis,int kind) {
    ElectricPhaseProblem result;result.epsilonByPhase={{epsilonNegative,epsilonPositive}};
    if(kind==1)for(auto& source:result.sourceByPhase)source=[axis](Vector2D p) {
        const auto q=p-origin(axis);return axis?-10*q.y:-8*q.x;
    };
    for(int side=0;side<4;++side) {
        result.boundaries[side].type=axis&&side==0?ElectricCartesianBoundary::Type::Neumann:ElectricCartesianBoundary::Type::Dirichlet;
        if(axis&&side==0)result.boundaries[side].value=[](Vector2D){return 0.;};
        else result.boundaries[side].value=[axis,kind](Vector2D p){return exact(p,axis,false,kind).value;};
    }
    return result;
}
void geometryChecks(const ElectricPhaseGeometry& geometry) {
    const bool axis=geometry.axisymmetric();const auto& report=geometry.report();
    check(report.interfaceFaces>0&&report.elements==geometry.elements().size(),"P1 adapter omitted physical phase elements/interface faces");
    check(report.minimumRelativePhaseArea>=1e-6,"declared healthy phase fixture contains an unsupported small cut");
    check(report.acceptedRefinement>0,"P1 cut matrices have no actual quadrature refinement evidence");
    check(std::max({report.maxRelativeReconstructionDifference,report.maxRelativeGradientGramDifference,
        report.maxRelativeStiffnessDifference,report.maxRelativeCondensedDifference})<=1e-9,
        "P1 actual R/G/K/S quadrature did not converge");
    long double ordinary=0,negativeArea=0,negativeVolume=0,arcLength=0,surface=0;
    int cutElements=0;bool testedKernel=false;
    for(const auto& element:geometry.elements()) {
        bool cut=false;
        for(int face:element.faces)cut=cut||geometry.faces()[face].kind==ElectricPhaseSkeletonFace::Kind::Interface;
        if(cut)++cutElements;
        for(const auto& node:element.geometry.volume) {
            check(node.weight>0&&std::isfinite(node.weight),"P1 volume integration is not a positive physical rule");
            ordinary+=node.weight;
            if(element.liquid){negativeArea+=node.weight;negativeVolume+=node.weight*(axis?2*pi*node.point.x:1);}
            check(element.liquid?scalar(node.point,axis)<1e-12:scalar(node.point,axis)>-1e-12,
                "phase quadrature point belongs to a different continuum phase");
        }
        if(cut&&!testedKernel) {
            H::LocalOperator op(element.geometry,element.liquid?epsilonNegative:epsilonPositive);
            const auto one=op.reduction([](Vector2D){return 1.;});
            const auto& K=op.stiffnessMatrix();const auto& S=op.condensedMatrix();
            check(ElectricContract::relativeSymmetry(K)<=1e-12&&ElectricContract::relativeSymmetry(S)<=1e-12,
                "true curved phase energy/Schur block is not symmetric");
            check((K*one).norm()<=1e-11*K.norm()*one.norm(),"true curved element lost the constant energy kernel");
            const auto faces=one.tail(S.rows());
            check((S*faces).norm()<=1e-11*S.norm()*faces.norm(),"true curved Schur lost the shared constant trace kernel");
            testedKernel=true;
        }
    }
    check(cutElements>0&&testedKernel,"P1 fixture contains no actual curved phase element");
    for(const auto& face:geometry.faces()) {
        check(face.length>0&&!face.nodes.empty(),"P1 retained an empty/zero-area skeleton face");
        for(const auto& node:face.nodes) {
            check(node.weight>0&&std::isfinite(node.coordinate),"P1 face rule is not positive and resolved");
            if(face.kind==ElectricPhaseSkeletonFace::Kind::Interface) {
                const auto radial=(node.point-origin(axis))*(1/radius);
                ElectricContract::close(scalar(node.point,axis),0,1e-10,"P1 interface node is not on the analytic Q2 circle");
                ElectricContract::close(node.outwardNormal.x,radial.x,1e-10,"P1 interface normal differs from same-circle negative-to-positive normal");
                ElectricContract::close(node.outwardNormal.y,radial.y,1e-10,"P1 interface normal differs from same-circle negative-to-positive normal");
                arcLength+=node.weight;surface+=node.weight*(axis?2*pi*node.point.x:1);
            } else if(axis)check(node.point.x>0,"RZ zero-area axis acquired a skeleton quadrature/DOF");
        }
    }
    ElectricContract::close(static_cast<double>(ordinary),1,1e-8,"P1 phase areas do not partition the physical square");
    ElectricContract::close(static_cast<double>(negativeArea),(axis?.5:1)*pi*radius*radius,1e-8,
        "P1 negative phase area differs from independent circle/hemicircle area");
    ElectricContract::close(static_cast<double>(negativeVolume),axis?4*pi*radius*radius*radius/3:pi*radius*radius,1e-8,
        "P1 negative phase metric differs from independent circle/sphere volume");
    ElectricContract::close(static_cast<double>(arcLength),(axis?1:2)*pi*radius,1e-8,
        "P1 actual interface ds differs from independent circle/semicircle length");
    ElectricContract::close(static_cast<double>(surface),axis?4*pi*radius*radius:2*pi*radius,1e-8,
        "P1 physical interface measure differs from independent circle/sphere area");
    std::cout<<"ELECTRIC_HHO_P1_GEOMETRY,"<<axis<<",elements="<<geometry.elements().size()<<",interface_faces="<<report.interfaceFaces
        <<",refinement="<<report.acceptedRefinement<<",minimum_area="<<report.minimumRelativePhaseArea
        <<",R="<<report.maxRelativeReconstructionDifference<<",G="<<report.maxRelativeGradientGramDifference
        <<",K="<<report.maxRelativeStiffnessDifference<<",S="<<report.maxRelativeCondensedDifference<<",PASS\n";
}
struct Error {double value=0,gradient=0;};
Error solveAndCheck(const QuadtreeMesh& mesh,const ElectricPhaseGeometry& geometry,int kind,
    std::optional<ElectricPhaseResult>* retained=nullptr) {
    const bool axis=geometry.axisymmetric();const auto input=problem(axis,kind);
    const auto result=solveElectricPhaseElectrostatics(geometry,input);result.validateSource(mesh,geometry.publication());
    check(&result.geometry().publication().state().geometry().snapshot()==&geometry.publication().state().geometry().snapshot(),
        "P1 solver introduced a second scalar geometry owner");
    const auto& report=result.report();
    check(report.interfaceFaces==geometry.report().interfaceFaces
        &&result.phaseElementCount()==geometry.elements().size(),"P1 result lost actual phase/shared interface layout");
    check(report.skeleton.relativeLinearResidual<=1e-11&&report.skeleton.relativeSymmetryDefect<=1e-12
        &&report.skeleton.minimumSkeletonPivot>0,"P1 skeleton did not pass SPD/symmetry/linear residual");
    check(report.maxRelativeInterfaceBalance<=1e-10,"P1 reported interface transmission moments are unbalanced");
    const Eigen::MatrixXd matrix(result.skeletonMatrix());
    check(ElectricContract::relativeSymmetry(matrix)<=1e-12,"P1 published shared skeleton matrix is not symmetric");
    Eigen::Map<const Eigen::VectorXd> solution(result.skeletonSolution().data(),result.skeletonSolution().size());
    Eigen::Map<const Eigen::VectorXd> rhs(result.skeletonRhs().data(),result.skeletonRhs().size());
    check((matrix*solution-rhs).norm()<=1e-11*(1+rhs.norm()),"P1 returned skeleton field does not solve its original equations");
    std::vector<std::array<long double,3>> sums(geometry.faces().size()),scales(geometry.faces().size());
    std::vector<int> incidences(geometry.faces().size(),0);
    long double valueError=0,gradientError=0,valueNorm=0,gradientNorm=0,sourceIntegral=0,sourceScale=0;
    long double boundaryFlux=0,boundaryScale=0,interfaceFluxNorm=0;
    double maximumValue=0,maximumGradient=0,maximumBalance=0;
    for(std::size_t e=0;e<geometry.elements().size();++e) {
        const auto& element=geometry.elements()[e];
        check(result.nativeCellIndex(e)==element.nativeCell&&result.liquid(e)==element.liquid,
            "phase result replaced native element identity");
        long double cellSource=0,cellSourceScale=0;
        for(const auto& node:element.geometry.volume) {
            const auto reference=exact(node.point,axis,element.liquid,kind);
            const auto gradient=result.gradientInElement(e,node.point);
            const double value=result.potentialInElement(e,node.point);
            const double dv=value-reference.value,dx=gradient.x-reference.gradient.x,dy=gradient.y-reference.gradient.y;
            const double weight=node.weight*(axis?2*pi*node.point.x:1);
            valueError+=weight*dv*dv;gradientError+=weight*(dx*dx+dy*dy);
            valueNorm+=weight*reference.value*reference.value;
            gradientNorm+=weight*(reference.gradient.x*reference.gradient.x+reference.gradient.y*reference.gradient.y);
            maximumValue=std::max(maximumValue,std::abs(dv)/(1+std::abs(reference.value)));
            maximumGradient=std::max(maximumGradient,std::hypot(dx,dy)/(1+std::hypot(reference.gradient.x,reference.gradient.y)));
            if(kind==1) {
                const auto q=node.point-origin(axis);const double f=axis?-10*q.y:-8*q.x;
                cellSource+=weight*f;cellSourceScale+=weight*std::abs(f);
            }
        }
        sourceIntegral+=cellSource;sourceScale+=cellSourceScale;
        long double cellFlux=0,cellFluxScale=0;
        for(std::size_t local=0;local<element.faces.size();++local) {
            const int id=element.faces[local];const auto& face=geometry.faces()[id];
            const auto residual=result.outwardFaceResidual(e,local);++incidences[id];
            for(int mode=0;mode<3;++mode){sums[id][mode]+=residual[mode];scales[id][mode]+=std::abs(residual[mode]);}
            cellFlux+=residual[0];cellFluxScale+=std::abs(residual[0]);
            if(face.boundarySide>=0){boundaryFlux+=residual[0];boundaryScale+=std::abs(residual[0]);}
            if(face.kind==ElectricPhaseSkeletonFace::Kind::Interface)interfaceFluxNorm+=std::abs(residual[0]);
        }
        const double relative=static_cast<double>(std::abs(cellFlux+cellSource)/(1+cellFluxScale+cellSourceScale));
        maximumBalance=std::max(maximumBalance,relative);check(relative<=1e-10,"P1 local phase Gauss residual/source balance failed");
    }
    for(std::size_t face=0;face<geometry.faces().size();++face) {
        const auto& physical=geometry.faces()[face];
        if(physical.boundarySide>=0)check(incidences[face]==1,"physical boundary phase subface has multiple owners");
        else {
            check(incidences[face]==2,"shared phase/interface face has no two actual incidences");
            for(int mode=0;mode<3;++mode) {
                const double relative=static_cast<double>(std::abs(sums[face][mode])/(1+scales[face][mode]));
                maximumBalance=std::max(maximumBalance,relative);check(relative<=1e-10,"P1 shared face/qs0 transmission moment balance failed");
            }
        }
    }
    check(std::abs(boundaryFlux+sourceIntegral)<=1e-10*(1+boundaryScale+sourceScale),
        "P1 global Gauss balance differs from physical source and numerical boundary flux");
    // An independent continuum boundary flux uses the analytic gradient, not
    // the solved residual. The exact whole-square source integral is zero in
    // these symmetric fixtures; both terms are kept and their scales recorded.
    long double referenceBoundary=0,referenceBoundaryScale=0;
    for(const auto& face:geometry.faces())if(face.boundarySide>=0)for(const auto& node:face.nodes) {
        const auto reference=exact(node.point,axis,face.liquid,kind);
        const auto normal=ElectricContract::normal(face.boundarySide);
        const double epsilon=face.liquid?epsilonNegative:epsilonPositive;
        const double flux=node.weight*(axis?2*pi*node.point.x:1)*epsilon*
            (reference.gradient.x*normal.x+reference.gradient.y*normal.y);
        referenceBoundary+=flux;referenceBoundaryScale+=std::abs(flux);
    }
    check(std::abs(referenceBoundary)<=1e-10*(1+referenceBoundaryScale),"analytic fixture violates independent continuous Gauss balance");
    check(std::abs(boundaryFlux-referenceBoundary)<=1e-10*(1+boundaryScale+referenceBoundaryScale),
        "P1 global boundary displacement disagrees with independent analytic Gauss reference");
    Error error{std::sqrt(static_cast<double>(valueError/(kind==0?1:valueNorm))),
        std::sqrt(static_cast<double>(gradientError/(kind==0?1:gradientNorm)))};
    if(kind!=2) {
        const double tolerance=kind==0?1e-10:2e-7;
        check(maximumValue<=tolerance&&maximumGradient<=tolerance,"P1 constant/special transmission MMS potential or E failed");
        if(kind==1)check(interfaceFluxNorm>1e-4,"P1 transmission MMS became a trivial zero-field interface test");
    }
    std::cout<<"ELECTRIC_HHO_P1_SOLVE,"<<axis<<','<<kind<<",elements="<<result.phaseElementCount()<<",V_L2="<<error.value
        <<",E_L2="<<error.gradient<<",V_max="<<maximumValue<<",E_max="<<maximumGradient
        <<",balance="<<maximumBalance<<",interface_flux="<<static_cast<double>(interfaceFluxNorm)<<",PASS\n";
    if(retained)*retained=result;
    return error;
}

Reference chargedExact(Vector2D p,bool axis,bool liquid,bool harmonic) {
    const auto q=p-origin(axis);const double squared=q.x*q.x+q.y*q.y;
    if(!harmonic) {
        const double factor=liquid?1/epsilonNegative:.5;
        return {factor*(squared-radius*radius),q*(2*factor)};
    }
    auto result=exact(p,false,liquid,2);
    if(!liquid) {
        const double coefficient=radius*radius/epsilonPositive;
        result.value-=coefficient*std::log(std::sqrt(squared)/radius);
        result.gradient=result.gradient-q*(coefficient/squared);
    }
    return result;
}
ElectricPhaseProblem chargedProblem(bool axis,bool harmonic) {
    ElectricPhaseProblem result;result.epsilonByPhase={{epsilonNegative,epsilonPositive}};
    if(!harmonic) {
        result.sourceByPhase[0]=[axis](Vector2D){return axis?-6.:-4.;};
        result.sourceByPhase[1]=[axis](Vector2D){return axis?-3.:-2.;};
    }
    for(int side=0;side<4;++side) {
        const bool pole=axis&&side==0;
        result.boundaries[side].type=pole?ElectricCartesianBoundary::Type::Neumann:ElectricCartesianBoundary::Type::Dirichlet;
        if(pole)result.boundaries[side].value=[](Vector2D){return 0.;};
        else result.boundaries[side].value=[axis,harmonic](Vector2D p){return chargedExact(p,axis,false,harmonic).value;};
    }
    return result;
}
bool identicalCharge(const std::vector<double>& a,const std::vector<double>& b) {
    return a.size()==b.size()&&(a.empty()||std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0);
}
// Independently evaluate the documented Field polynomial, not P_Q nor the
// stencil's evaluator. The integration below uses the analytic circular arc.
double densityPolynomial(const SurfaceChargeDensityStencil& stencil,
    const std::vector<double>& charge,double arclength) {
    std::array<double,3> coefficient{};
    for(int row=0;row<3;++row)for(int column=0;column<3;++column)
        coefficient[row]+=stencil.polynomialWeights[3*row+column]*charge.at(stencil.neighbor[column]);
    const double u=(arclength-stencil.origin)/stencil.scale;
    const double coordinate=stencil.evenAtAxis?u*u:u;
    return coefficient[0]+coordinate*(coefficient[1]+coordinate*coefficient[2]);
}
std::array<double,3> circleChargeMoments(const SurfaceChargeField& field,std::size_t panel,
    const std::vector<double>& charge,bool axis) {
    const auto& material=field.geometry().panels().at(panel);
    const auto begin=material.begin-origin(axis),end=material.end-origin(axis);
    const double first=std::atan2(begin.y,begin.x);
    double angle=std::atan2(end.y,end.x)-first;
    if(angle<=0)angle+=2*pi;
    check(angle>0&&angle<pi,"charge reference lost the native circle orientation");
    const double length=radius*angle;
    ElectricContract::close(material.measure,length,1e-10,"material ds differs from independent analytic circle arc");
    const auto stencil=field.densityStencil(panel);
    std::array<double,3> result{};
    constexpr int intervals=32;
    for(int interval=0;interval<intervals;++interval)for(int node=0;node<8;++node) {
        const double s=length*(interval+.5*(1+ElectricContract::nodes[node]))/intervals;
        const double zeta=2*s/length-1;
        const double metric=axis?2*pi*radius*std::cos(first+s/radius):1.;
        const double weight=.5*length/intervals*ElectricContract::weights[node]*metric;
        const double density=densityPolynomial(stencil,charge,s);
        const std::array<double,3> basis{{1,zeta,.5*(3*zeta*zeta-1)}};
        for(int mode=0;mode<3;++mode)result[mode]+=weight*density*basis[mode];
    }
    return result;
}
void chargeMomentChecks(const ElectricPhaseGeometry& geometry,const SurfaceChargeField& field,
    const ElectricSurfaceChargeMomentOperator& op,bool fullReference) {
    const bool axis=geometry.axisymmetric();const auto& stored=field.transport().panelCharge();
    const auto saved=stored;const auto& panels=field.geometry().panels();
    check(op.panelCount()==panels.size()&&op.report().panels==panels.size(),"P_Q omitted a material panel");
    check(op.projectionMatrix().rows()==3*static_cast<int>(panels.size())
        &&op.projectionMatrix().cols()==static_cast<int>(panels.size()),"P_Q shape is not three moments per unique Q owner");
    check(op.interfaceInjectionMatrix().rows()==3*static_cast<int>(geometry.faces().size())
        &&op.interfaceInjectionMatrix().cols()==3*static_cast<int>(panels.size()),"B_Gamma does not cover the original shared skeleton");
    check(op.report().maxRelativeMomentDifference<=1e-9
        &&op.report().maxRelativePanelLengthDifference<=1e-9,"P_Q actual high-moment arclength refinement did not converge");
    const auto load=op.apply(stored);
    check(identicalCharge(load.panelCharge(),stored),"candidate load changed input integrated Q bits");
    double total=0,nonzeroHighMoment=0;
    std::vector<double> first(stored.size()),second(stored.size()),sum(stored.size());
    bool negative=false,positive=false;
    std::set<int> interfaceIds;
    for(std::size_t panel=0;panel<panels.size();++panel) {
        const int face=op.interfaceFace(panel);
        check(interfaceIds.insert(face).second&&geometry.faces().at(face).kind==ElectricPhaseSkeletonFace::Kind::Interface,
            "material panel is not mapped bijectively to an original shared Gamma face");
        check(op.component(panel)==0&&op.graph(panel)==static_cast<int>(panel)
            &&geometry.faces()[face].component==op.component(panel)&&geometry.faces()[face].graph==op.graph(panel),
            "P_Q lost actual component/graph order");
        check(std::memcmp(&load.moments()[panel][0],&stored[panel],sizeof(double))==0,
            "eta0 is not a direct bitwise copy of the panel Q");
        ElectricContract::close(stored[panel],radius*field.transport().panelMeasure()[panel],0,
            "constant initial charge was not integrated over the actual physical material measure");
        total+=stored[panel];nonzeroHighMoment=std::max(nonzeroHighMoment,std::abs(load.moments()[panel][1]));
        first[panel]=radius*field.transport().panelMeasure()[panel]*(.15+std::sin(.91*(panel+1)));
        second[panel]=radius*field.transport().panelMeasure()[panel]*std::cos(.67*(panel+1));
        sum[panel]=first[panel]+second[panel];negative=negative||first[panel]<0;positive=positive||first[panel]>0;
    }
    check(negative&&positive,"P_Q signed trial vector became a constant positive special case");
    check(interfaceIds.size()==geometry.report().interfaceFaces,"P_Q omitted or duplicated a native interface graph");
    ElectricContract::close(total,axis?4*pi*radius*radius*radius:2*pi*radius*radius,1e-10,
        "integrated total Q differs from independent circle/sphere charge");
    if(axis)check(nonzeroHighMoment>1e-6,"RZ weighted high moments were replaced with zero constant-density moments");
    const auto a=op.apply(first),b=op.apply(second),c=op.apply(sum);
    const Eigen::Map<const Eigen::VectorXd> q(first.data(),first.size());
    const Eigen::VectorXd projected=op.projectionMatrix()*q;
    Eigen::VectorXd flat(3*static_cast<int>(panels.size()));
    for(std::size_t panel=0;panel<panels.size();++panel)for(int mode=0;mode<3;++mode) {
        ElectricContract::close(c.moments()[panel][mode],a.moments()[panel][mode]+b.moments()[panel][mode],1e-12,
            "P_Q is not linear for signed nonuniform trial charges");
        ElectricContract::close(projected[3*panel+mode],a.moments()[panel][mode],1e-12,
            "immutable candidate moments differ from the published linear P_Q");
        flat[3*panel+mode]=a.moments()[panel][mode];
    }
    const Eigen::VectorXd injected=op.interfaceInjectionMatrix()*flat;
    for(std::size_t face=0;face<geometry.faces().size();++face)for(int mode=0;mode<3;++mode)
        ElectricContract::close(injected[3*face+mode],a.faceMoments()[face][mode],1e-12,
            "B_Gamma injection changed charge moments on the shared skeleton");
    double maximumReference=0,maximumSource=0;
    if(fullReference) {
        for(std::size_t panel=0;panel<panels.size();++panel) {
            const auto reference=circleChargeMoments(field,panel,first,axis);
            for(int mode=0;mode<3;++mode) {
                const double error=std::abs(a.moments()[panel][mode]-reference[mode])/(1+std::abs(reference[mode]));
                maximumReference=std::max(maximumReference,error);
                check(error<=1e-9,"P_Q high moment differs from independent true-arclength density integral");
            }
        }
        const auto originalSource=field.source(first);
        // These are exact saved interpolation knots of the old Field's own
        // Frame -> ds map, not an assumed equivalence of graph fraction and ds.
        const std::array<std::size_t,3> samples{{0,panels.size()/2,panels.size()-1}};
        for(auto panel:samples) {
            const auto& nodes=panels[panel].quadratureNodes;
            check(!nodes.empty(),"material panel lost its original source quadrature frames");
            const auto& node=nodes[nodes.size()/2];
            const double value=originalSource.evaluate(node.frame).density;
            const double reference=densityPolynomial(field.densityStencil(panel),first,node.localArclength);
            const double error=std::abs(value-reference)/(1+std::abs(reference));
            maximumSource=std::max(maximumSource,error);
            check(error<=1e-10,"density stencil is not the actual old Field.source trial reconstruction");
        }
    }
    check(identicalCharge(field.transport().panelCharge(),saved),"P_Q preparation/apply mutated the sole Field/Transport Q owner");
    auto wrongClock=geometry.timeLayer();++wrongClock.step;
    rejected([&]{ElectricSurfaceChargeMomentOperator::prepare(geometry,{&field},wrongClock);},
        "P_Q accepted a foreign explicit coeval clock declaration");
    rejected([&]{ElectricSurfaceChargeMomentOperator::prepare(geometry,{&field,&field},geometry.timeLayer());},
        "P_Q accepted a foreign component inventory");
    rejected([&]{op.apply(std::vector<double>(stored.size()+1));},"P_Q accepted an invalid charge vector size");
    auto invalid=stored;invalid[0]=std::numeric_limits<double>::quiet_NaN();
    rejected([&]{op.apply(invalid);},"P_Q accepted nonfinite trial charge");
    std::cout<<"ELECTRIC_HHO_P2_MOMENTS,"<<axis<<",panels="<<panels.size()<<",total_Q="<<total
        <<",eta_reference="<<maximumReference<<",legacy_source="<<maximumSource
        <<",quadrature="<<op.report().maxRelativeMomentDifference<<",PASS\n";
}
void zeroChargeRegression(const ElectricPhaseGeometry& geometry,const ElectricSurfaceChargeMomentOperator& op,
    const ElectricPhaseResult& baseline) {
    const auto zero=op.apply(std::vector<double>(op.panelCount(),0.));
    const auto result=solveElectricPhaseElectrostatics(geometry,problem(geometry.axisymmetric(),1),zero);
    check(result.surfaceChargeLoad()!=nullptr&&baseline.surfaceChargeLoad()==nullptr,
        "typed charged and original P1 results lost their explicit load distinction");
    const Eigen::MatrixXd before(baseline.skeletonMatrix()),after(result.skeletonMatrix());
    check((after-before).norm()==0,"Q0 regression changed the shared P1 skeleton matrix");
    check(result.skeletonRhs()==baseline.skeletonRhs(),"Q0 injection changed the original P1 right-hand side");
    for(std::size_t e=0;e<geometry.elements().size();++e) {
        const auto old=baseline.reconstructedPotentialP3(e),current=result.reconstructedPotentialP3(e);
        for(int j=0;j<10;++j)ElectricContract::close(current[j],old[j],1e-12,"Q0 changed the retained P1 reconstructed field");
    }
}

struct CircularPanelReference {
    double first=0,angle=0,physicalMeasure=0,current=0;
};
CircularPanelReference circularPanelCurrent(const LevelSetInterfacePanels::Panel& panel,bool axis,
    bool harmonic,const std::array<double,2>& sigma) {
    const auto begin=panel.begin-origin(axis),end=panel.end-origin(axis);
    const double first=std::atan2(begin.y,begin.x);
    double angle=std::atan2(end.y,end.x)-first;
    if(angle<=0)angle+=2*pi;
    check(angle>0&&angle<pi,"P3 analytic current lost the original oriented circular arc");
    const double last=first+angle;
    const double measure=axis?2*pi*radius*radius*(std::sin(last)-std::sin(first)):radius*angle;
    check(measure>0,"P3 analytic physical panel measure is not positive");
    const double integral=harmonic?
        radius*(5*sigma[1]-sigma[0])/3*(std::sin(last)-std::sin(first))-sigma[1]*radius*radius*angle:
        (sigma[1]-.4*sigma[0])*radius*measure;
    return {first,angle,measure,integral};
}
struct CurrentBalance {
    double total=0,exterior=0,bulk=0,maximum=0;
};
CurrentBalance currentConservationChecks(const ElectricPhaseResult& result,
    const ElectricPanelCurrentOperator& current,const std::array<double,2>& sigma,bool harmonic) {
    const auto& geometry=result.geometry();const bool axis=geometry.axisymmetric();
    const std::array<double,2> rho{{sigma[0]/epsilonNegative,sigma[1]/epsilonPositive}};
    std::vector<std::array<double,3>> outward(geometry.faces().size()),scales(geometry.faces().size());
    std::vector<int> incidence(geometry.faces().size()),panelAtFace(geometry.faces().size(),-1);
    const auto& op=current.momentOperator();
    for(std::size_t panel=0;panel<op.panelCount();++panel)
        panelAtFace[op.interfaceFace(panel)]=static_cast<int>(panel);
    std::vector<double> rawPanelCurrent(op.panelCount(),0.);
    CurrentBalance balance;double exteriorScale=0,bulkScale=0,panelScale=0;
    for(std::size_t element=0;element<geometry.elements().size();++element) {
        const auto& retained=geometry.elements()[element];const double ratio=rho[retained.liquid?0:1];
        double source=0,sourceScale=0,flux=0,fluxScale=0;
        if(!harmonic)for(const auto& node:retained.geometry.volume) {
            const double f=axis?(retained.liquid?-6.:-3.):(retained.liquid?-4.:-2.);
            const double contribution=node.weight*(axis?2*pi*node.point.x:1)*ratio*f;
            source+=contribution;sourceScale+=std::abs(contribution);
        }
        for(std::size_t local=0;local<retained.faces.size();++local) {
            const int id=retained.faces[local];
            const auto displacement=result.physicalOutwardDisplacementMoments(element,local);++incidence[id];
            if(panelAtFace[id]>=0) {
                const auto raw=result.outwardFaceResidual(element,local);
                rawPanelCurrent[static_cast<std::size_t>(panelAtFace[id])]-=ratio*raw[0];
            }
            for(int mode=0;mode<3;++mode) {
                const double value=ratio*displacement[mode];
                outward[id][mode]+=value;scales[id][mode]+=std::abs(value);
            }
            flux+=ratio*displacement[0];fluxScale+=std::abs(ratio*displacement[0]);
            if(geometry.faces()[id].boundarySide>=0) {
                balance.exterior+=ratio*displacement[0];exteriorScale+=std::abs(ratio*displacement[0]);
            }
        }
        const double error=std::abs(flux-source)/(1+fluxScale+sourceScale);
        balance.maximum=std::max(balance.maximum,error);
        check(error<=1e-10,"P3 phase-local physical current Gauss balance failed");
        balance.bulk+=source;bulkScale+=sourceScale;
    }
    for(std::size_t panel=0;panel<op.panelCount();++panel) {
        ElectricContract::close(current.panelCurrent()[panel],rawPanelCurrent[panel],1e-10,
            "P3 panel current is not -rho_minus*raw_r_minus0-rho_plus*raw_r_plus0");
        const auto negative=current.normalCurrentMoments(panel,true),positive=current.normalCurrentMoments(panel,false);
        const auto dnNegative=result.normalDisplacementMoments(op.interfaceFace(panel),true);
        const auto dnPositive=result.normalDisplacementMoments(op.interfaceFace(panel),false);
        const auto pNegative=current.normalCurrentPolynomialP2(panel,true),pPositive=current.normalCurrentPolynomialP2(panel,false);
        const auto dNegative=result.normalDisplacementPolynomialP2(op.interfaceFace(panel),true);
        const auto dPositive=result.normalDisplacementPolynomialP2(op.interfaceFace(panel),false);
        for(int mode=0;mode<3;++mode) {
            ElectricContract::close(negative[mode],rho[0]*dnNegative[mode],1e-12,
                "P3 negative current moment is not based on the actual one-sided displacement");
            ElectricContract::close(positive[mode],rho[1]*dnPositive[mode],1e-12,
                "P3 positive current moment is not based on the actual one-sided displacement");
            ElectricContract::close(pNegative[mode],rho[0]*dNegative[mode],1e-12,"P3 negative current P2 used another electric field");
            ElectricContract::close(pPositive[mode],rho[1]*dPositive[mode],1e-12,"P3 positive current P2 used another electric field");
        }
        ElectricContract::close(current.panelCurrent()[panel],negative[0]-positive[0],1e-10,
            "P3 integrated panel current differs from the original two-sided current moments");
        balance.total+=current.panelCurrent()[panel];panelScale+=std::abs(current.panelCurrent()[panel]);
    }
    for(std::size_t face=0;face<geometry.faces().size();++face) {
        if(geometry.faces()[face].boundarySide>=0){check(incidence[face]==1,"P3 boundary current has multiple owners");continue;}
        check(incidence[face]==2,"P3 shared current face lacks two actual incidences");
        std::array<double,3> expected{};
        if(panelAtFace[face]>=0) {
            const auto panel=static_cast<std::size_t>(panelAtFace[face]);
            const auto negative=current.normalCurrentMoments(panel,true),positive=current.normalCurrentMoments(panel,false);
            for(int mode=0;mode<3;++mode)expected[mode]=negative[mode]-positive[mode];
        }
        for(int mode=0;mode<3;++mode) {
            const double error=std::abs(outward[face][mode]-expected[mode])/(1+scales[face][mode]+std::abs(expected[mode]));
            balance.maximum=std::max(balance.maximum,error);
            check(error<=1e-10,"P3 shared Cartesian/Gamma three-moment current balance failed");
        }
    }
    const double residual=balance.exterior+balance.total-balance.bulk;
    check(std::abs(residual)<=1e-10*(1+exteriorScale+panelScale+bulkScale),
        "P3 global exterior-plus-interface current balance failed against physical volume sources");
    const auto& report=current.report();
    check(report.panels==op.panelCount()&&report.relativeGlobalCurrentBalance<=1e-10
        &&report.maxRelativeAffineResidual<=1e-10,"P3 reported physical current/affine gates failed");
    ElectricContract::close(report.totalPanelCurrent,balance.total,1e-10,"P3 report total current differs from actual original panels");
    ElectricContract::close(report.exteriorCurrent,balance.exterior,1e-10,"P3 report exterior current omitted its physical boundary load");
    ElectricContract::close(report.bulkCurrent,balance.bulk,1e-10,"P3 report bulk current differs from independent analytic source integration");
    ElectricContract::close(report.globalCurrentBalance,residual,1e-10,"P3 report substituted another current conservation equation");
    const auto affine=current.apply(result.skeletonSolution());
    check(current.matrix().rows()==static_cast<int>(op.panelCount())
        &&current.matrix().cols()==static_cast<int>(result.skeletonSolution().size())
        &&current.offset().size()==op.panelCount(),"P3 affine current has a foreign free-skeleton layout");
    const Eigen::Map<const Eigen::VectorXd> u(result.skeletonSolution().data(),result.skeletonSolution().size());
    const Eigen::Map<const Eigen::VectorXd> offset(current.offset().data(),current.offset().size());
    const Eigen::VectorXd explicitAffine=current.matrix()*u+offset;
    for(std::size_t panel=0;panel<op.panelCount();++panel) {
        ElectricContract::close(affine[panel],current.panelCurrent()[panel],1e-10,"P3 C*u+affine_offset does not reproduce the original raw panel current");
        ElectricContract::close(explicitAffine[panel],current.panelCurrent()[panel],1e-10,"P3 exposed C/affine_offset differ from the original raw current");
    }
    return balance;
}
void currentChecks(const ElectricPhaseResult& result,const ElectricSurfaceChargeMomentOperator& op,
    const SurfaceChargeField& field,bool harmonic,bool withTrial,const ElectricSurfaceChargeMomentOperator* foreign) {
    const auto& geometry=result.geometry();const bool axis=geometry.axisymmetric();
    constexpr std::array<double,2> sigma{{.7,.03}};
    const auto chargeBefore=field.transport().panelCharge();
    const auto solvedBefore=result.skeletonSolution();
    const auto current=result.preparePanelCurrent(op,sigma);
    check(current.epsilonByPhase()==std::array<double,2>{{epsilonNegative,epsilonPositive}}
        &&current.conductivityByPhase()==sigma,"P3 conductivity ratio is not bound to the actual solved 5:1 permittivity");
    check(current.geometry().geometryIdentity()==geometry.geometryIdentity()
        &&&current.geometry().publication().state().geometry().snapshot()==&geometry.publication().state().geometry().snapshot(),
        "P3 current introduced another scalar geometry owner");
    const auto layer=current.timeLayer(),expectedLayer=geometry.timeLayer();
    check(layer.step==expectedLayer.step&&layer.begin==expectedLayer.begin&&layer.end==expectedLayer.end
        &&layer.position==expectedLayer.position,"P3 current lost a field of the complete physical clock");
    const auto balance=currentConservationChecks(result,current,sigma,harmonic);
    const auto& panels=field.geometry().panels();double error=0,norm=0,maximum=0;
    bool positive=false,negative=false;
    for(std::size_t panel=0;panel<panels.size();++panel) {
        const auto reference=circularPanelCurrent(panels[panel],axis,harmonic,sigma);
        const double difference=current.panelCurrent()[panel]-reference.current;
        const double normalized=std::abs(difference)/(1+std::abs(reference.current));
        maximum=std::max(maximum,normalized);
        error+=difference*difference/reference.physicalMeasure;
        norm+=reference.current*reference.current/reference.physicalMeasure;
        positive=positive||reference.current>0;negative=negative||reference.current<0;
        if(!harmonic)check(normalized<=2e-7,"P3 panel current differs from independent charged-MMS circular surface integral");
    }
    const double relative=std::sqrt(error/norm);
    if(harmonic) {
        check(relative<=.15,"P3 harmonic true-arc panel-average current failed independent weighted L2 accuracy");
        check(positive&&negative,"P3 harmonic current did not retain both signs");
    }
    const double volume=axis?4*pi*radius*radius*radius/3:pi*radius*radius;
    const double charge=axis?4*pi*radius*radius*radius:2*pi*radius*radius;
    const double totalReference=harmonic?-sigma[1]*charge:(sigma[1]-.4*sigma[0])*(axis?4:2)*pi*radius*(axis?radius*radius:radius);
    const double exteriorReference=harmonic?sigma[1]*charge:(axis?-3*pi*sigma[1]:-2*sigma[1]);
    const double bulkReference=harmonic?0:axis?-3*pi*sigma[1]+(3*sigma[1]-1.2*sigma[0])*volume:
        -2*sigma[1]+(2*sigma[1]-.8*sigma[0])*volume;
    ElectricContract::close(balance.total,totalReference,1e-10,"P3 total panel current differs from independent circle/sphere flux");
    ElectricContract::close(balance.exterior,exteriorReference,1e-10,"P3 exterior current differs from independent continuum Gauss flux");
    ElectricContract::close(balance.bulk,bulkReference,1e-10,"P3 total Ohmic volume source differs from independent material volumes");
    const auto gasZero=result.preparePanelCurrent(op,{{sigma[0],0.}});
    currentConservationChecks(result,gasZero,{{sigma[0],0.}},harmonic);
    for(std::size_t panel=0;panel<panels.size();++panel) {
        const auto plus=gasZero.normalCurrentMoments(panel,false);
        for(double value:plus)check(value==0,"zero gas conductivity produced a nonzero positive current moment");
        if(!harmonic)ElectricContract::close(gasZero.panelCurrent()[panel],
            circularPanelCurrent(panels[panel],axis,false,{{sigma[0],0.}}).current,2e-7,
            "zero gas conductivity lost the independent liquid current reference");
    }
    const auto bothZero=result.preparePanelCurrent(op,{{0.,0.}});
    check(bothZero.matrix().norm()==0&&bothZero.report().totalPanelCurrent==0
        &&bothZero.report().exteriorCurrent==0&&bothZero.report().bulkCurrent==0,"zero conductivity retained a current response");
    for(double value:bothZero.offset())check(value==0,"zero conductivity retained a source current offset");
    for(double value:bothZero.panelCurrent())check(value==0,"zero conductivity retained a panel current");
    const auto zeroApplied=bothZero.apply(result.skeletonSolution());
    for(double value:zeroApplied)check(value==0,"zero conductivity affine apply produced a current");
    rejected([&]{result.preparePanelCurrent(op,{{-1.,sigma[1]}});},"P3 accepted negative liquid conductivity");
    rejected([&]{result.preparePanelCurrent(op,{{sigma[0],-1.}});},"P3 accepted negative gas conductivity");
    rejected([&]{result.preparePanelCurrent(op,{{std::numeric_limits<double>::infinity(),sigma[1]}});},
        "P3 accepted nonfinite conductivity");
    rejected([&]{result.preparePanelCurrent(op,{{std::numeric_limits<double>::denorm_min(),sigma[1]}});},
        "P3 accepted an unrepresentable positive conductivity/permittivity ratio");
    rejected([&]{current.apply(std::vector<double>(solvedBefore.size()+1));},"P3 current accepted a foreign skeleton size");
    auto nonfinite=solvedBefore;nonfinite[0]=std::numeric_limits<double>::quiet_NaN();
    rejected([&]{current.apply(nonfinite);},"P3 current accepted a nonfinite free trace");
    if(foreign)rejected([&]{result.preparePanelCurrent(*foreign,sigma);},"P3 accepted a foreign prepared moment/geometry owner");
    if(withTrial) {
        check(!axis&&!harmonic,"P3 additional signed-Q solve escaped its single XY4 fixture");
        auto trialCharge=chargeBefore;bool trialPositive=false,trialNegative=false;
        for(std::size_t panel=0;panel<trialCharge.size();++panel) {
            trialCharge[panel]*=.15+std::sin(.91*(panel+1));
            trialPositive=trialPositive||trialCharge[panel]>0;trialNegative=trialNegative||trialCharge[panel]<0;
        }
        check(trialPositive&&trialNegative,"P3 affine trial charge has no signed nonuniform perturbation");
        const auto trialLoad=op.apply(trialCharge);
        const auto trial=solveElectricPhaseElectrostatics(geometry,chargedProblem(false,false),trialLoad);
        check(trial.report().skeleton.relativeLinearResidual<=1e-11&&trial.report().skeleton.relativeSymmetryDefect<=1e-12
            &&trial.report().skeleton.minimumSkeletonPivot>0&&trial.report().maxRelativeInterfaceBalance<=1e-10
            &&trial.report().maxRelativeLocalGaussResidual<=1e-10&&trial.report().relativeGlobalGaussResidual<=1e-10,
            "P3 signed-Q affine probe does not solve the original charged equations");
        const auto changed=trial.preparePanelCurrent(op,sigma);
        currentConservationChecks(trial,changed,sigma,false);
        const Eigen::MatrixXd before(current.matrix()),after(changed.matrix());
        check((after-before).norm()==0&&changed.offset()==current.offset(),"P3 C/affine_offset changed with trial Q instead of only its solved trace");
        const auto prediction=current.apply(trial.skeletonSolution());
        std::vector<std::vector<std::size_t>> cellPanels(
            geometry.publication().state().geometry().snapshot().nativeCellCount());
        for(std::size_t panel=0;panel<panels.size();++panel)cellPanels.at(panels[panel].cellIndex).push_back(panel);
        double maximumDensityVariation=0;bool multipleGraphs=false;
        for(const auto& group:cellPanels)for(std::size_t first=0;first<group.size();++first)
            for(std::size_t second=first+1;second<group.size();++second) {
                multipleGraphs=true;
                maximumDensityVariation=std::max(maximumDensityVariation,std::abs(
                    changed.panelCurrent()[group[first]]/panels[group[first]].physicalMeasure
                    -changed.panelCurrent()[group[second]]/panels[group[second]].physicalMeasure));
            }
        check(multipleGraphs&&maximumDensityVariation>1e-6,
            "P3 coarse signed-Q probe assigned one area-uniform current to multiple native graphs in a cell");
        double response=0,traceResponse=0;
        for(std::size_t panel=0;panel<prediction.size();++panel) {
            ElectricContract::close(prediction[panel],changed.panelCurrent()[panel],1e-10,
                "P3 original C/affine_offset cannot predict a second independent signed-Q solve's raw currents");
            response+=std::abs(changed.panelCurrent()[panel]-current.panelCurrent()[panel]);
        }
        for(std::size_t mode=0;mode<solvedBefore.size();++mode)traceResponse+=std::abs(trial.skeletonSolution()[mode]-solvedBefore[mode]);
        check(response>1e-5&&traceResponse>1e-5&&current.matrix().norm()>0,
            "P3 affine response probe became a zero-matrix or unchanged-state single-point check");
        std::cout<<"ELECTRIC_HHO_P3_AFFINE_PROBE,current_response="<<response<<",trace_response="<<traceResponse
            <<",multi_graph_density_variation="<<maximumDensityVariation<<",PASS\n";
    }
    check(identicalCharge(field.transport().panelCharge(),chargeBefore)&&result.skeletonSolution()==solvedBefore,
        "P3 current preparation/trial solve mutated the accepted charge or retained potential owner");
    std::cout<<"ELECTRIC_HHO_P3_CURRENT,"<<axis<<",harmonic="<<harmonic<<",I_max="<<maximum<<",I_L2="<<relative
        <<",total="<<balance.total<<",exterior="<<balance.exterior<<",bulk="<<balance.bulk
        <<",balance="<<balance.maximum<<",PASS\n";
}

Error chargedSolveAndCheck(const QuadtreeMesh& mesh,const ElectricPhaseGeometry& geometry,
    const ElectricSurfaceChargeMomentOperator& op,const std::vector<double>& charge,bool harmonic,
    std::optional<ElectricPhaseResult>* retained=nullptr) {
    const bool axis=geometry.axisymmetric();const auto input=chargedProblem(axis,harmonic);
    const auto load=op.apply(charge);const auto result=solveElectricPhaseElectrostatics(geometry,input,load);
    if(!harmonic && !axis) {
        ElectricCartesianControls invalid;
        invalid.relativeConservationTolerance=1e-8;
        rejected([&]{solveElectricPhaseElectrostatics(geometry,input,load,invalid);},
            "charged solve accepted conservation tolerance above the independent E-R5 bound");
    }
    result.validateSource(mesh,geometry.publication());
    check(result.surfaceChargeLoad()!=nullptr&&identicalCharge(result.surfaceChargeLoad()->panelCharge(),charge),
        "charged result lost immutable Q provenance");
    const auto& report=result.report();
    check(report.skeleton.relativeLinearResidual<=1e-11&&report.skeleton.relativeSymmetryDefect<=1e-12
        &&report.skeleton.minimumSkeletonPivot>0,"charged skeleton failed SPD/symmetry/linear gates");
    check(report.maxRelativeInterfaceBalance<=1e-10&&report.maxRelativeLocalGaussResidual<=1e-10
        &&report.relativeGlobalGaussResidual<=1e-10,"charged report failed jump/local/global Gauss gates");
    const Eigen::MatrixXd matrix(result.skeletonMatrix());
    check(ElectricContract::relativeSymmetry(matrix)<=1e-12,"charged shared skeleton matrix is not symmetric");
    const Eigen::Map<const Eigen::VectorXd> x(result.skeletonSolution().data(),result.skeletonSolution().size());
    const Eigen::Map<const Eigen::VectorXd> rhs(result.skeletonRhs().data(),result.skeletonRhs().size());
    check((matrix*x-rhs).norm()<=1e-11*(1+rhs.norm()),"charged result does not solve its actual skeleton equations");
    std::vector<std::array<double,3>> outward(geometry.faces().size()),scales(geometry.faces().size());
    std::vector<int> incidence(geometry.faces().size());
    double valueError=0,gradientError=0,valueNorm=0,gradientNorm=0,maximumValue=0,maximumGradient=0;
    double boundary=0,boundaryScale=0,totalSource=0,sourceScale=0,maximumBalance=0;
    std::array<double,2> fieldNorm{};
    for(std::size_t e=0;e<geometry.elements().size();++e) {
        const auto& element=geometry.elements()[e];double localSource=0,localSourceScale=0;
        for(const auto& node:element.geometry.volume) {
            const auto reference=chargedExact(node.point,axis,element.liquid,harmonic);
            const auto actual=result.gradientInElement(e,node.point);
            const double dv=result.potentialInElement(e,node.point)-reference.value;
            const double dx=actual.x-reference.gradient.x,dy=actual.y-reference.gradient.y;
            const double weight=node.weight*(axis?2*pi*node.point.x:1);
            valueError+=weight*dv*dv;gradientError+=weight*(dx*dx+dy*dy);
            valueNorm+=weight*reference.value*reference.value;
            gradientNorm+=weight*(reference.gradient.x*reference.gradient.x+reference.gradient.y*reference.gradient.y);
            fieldNorm[element.liquid?0:1]+=weight*(actual.x*actual.x+actual.y*actual.y);
            maximumValue=std::max(maximumValue,std::abs(dv)/(1+std::abs(reference.value)));
            maximumGradient=std::max(maximumGradient,std::hypot(dx,dy)/(1+std::hypot(reference.gradient.x,reference.gradient.y)));
            if(!harmonic) {
                const double source=axis?(element.liquid?-6.:-3.):(element.liquid?-4.:-2.);
                localSource+=weight*source;localSourceScale+=weight*std::abs(source);
            }
        }
        totalSource+=localSource;sourceScale+=localSourceScale;
        double cellOutward=0,cellScale=0;
        for(std::size_t local=0;local<element.faces.size();++local) {
            const int id=element.faces[local];const auto& face=geometry.faces()[id];
            const auto physical=result.physicalOutwardDisplacementMoments(e,local);
            const auto primal=result.outwardFaceResidual(e,local);++incidence[id];
            for(int mode=0;mode<3;++mode) {
                // The only Neumann datum in these fixtures is the absent axis
                // face. Physical D is therefore the negative primal residual.
                ElectricContract::close(physical[mode],-primal[mode],1e-12,"physical flux getter changed the actual local residual");
                outward[id][mode]+=physical[mode];scales[id][mode]+=std::abs(physical[mode]);
            }
            cellOutward+=physical[0];cellScale+=std::abs(physical[0]);
            if(face.boundarySide>=0){boundary+=physical[0];boundaryScale+=std::abs(physical[0]);}
        }
        const double defect=std::abs(cellOutward-localSource)/(1+cellScale+localSourceScale);
        maximumBalance=std::max(maximumBalance,defect);
        check(defect<=1e-10,"charged local physical Gauss balance failed against analytic bulk source");
    }
    double totalCharge=0;for(double value:charge)totalCharge+=value;
    for(std::size_t id=0;id<geometry.faces().size();++id) {
        const auto& face=geometry.faces()[id];if(face.boundarySide>=0){check(incidence[id]==1,"boundary lost unique phase incidence");continue;}
        check(incidence[id]==2,"charged shared face lost its two physical incidences");
        const bool interface=face.kind==ElectricPhaseSkeletonFace::Kind::Interface;
        const auto eta=interface?result.interfaceChargeMoments(id):std::array<double,3>{};
        for(int mode=0;mode<3;++mode) {
            const double defect=std::abs(outward[id][mode]+eta[mode])/(1+scales[id][mode]+std::abs(eta[mode]));
            maximumBalance=std::max(maximumBalance,defect);
            check(defect<=1e-10,"charged shared Cartesian/Gamma three-moment balance failed");
        }
        if(interface) {
            const auto negative=result.normalDisplacementMoments(id,true),positive=result.normalDisplacementMoments(id,false);
            const auto minusPolynomial=result.normalDisplacementPolynomialP2(id,true);
            const auto plusPolynomial=result.normalDisplacementPolynomialP2(id,false);
            Eigen::Matrix3d mass=Eigen::Matrix3d::Zero();std::array<double,3> analyticEta{};
            for(const auto& node:face.nodes) {
                const double zeta=node.coordinate;
                const Eigen::Vector3d basis(1,zeta,.5*(3*zeta*zeta-1));
                const double weight=node.weight*(axis?2*pi*node.point.x:1);
                mass+=weight*basis*basis.transpose();
                for(int mode=0;mode<3;++mode)analyticEta[mode]+=weight*radius*basis[mode];
            }
            for(int mode=0;mode<3;++mode) {
                ElectricContract::close(positive[mode]-negative[mode],eta[mode],1e-10,
                    "numerical one-sided D traces do not satisfy the charged three-moment jump");
                ElectricContract::close(eta[mode],load.faceMoments()[id][mode],0,"result substituted a second interface charge load");
                ElectricContract::close(eta[mode],analyticEta[mode],1e-9,
                    "constant physical qs charge moments differ from independent surface integration");
                double negativeFromPolynomial=0,positiveFromPolynomial=0;
                for(int column=0;column<3;++column) {
                    negativeFromPolynomial+=mass(mode,column)*minusPolynomial[column];
                    positiveFromPolynomial+=mass(mode,column)*plusPolynomial[column];
                }
                ElectricContract::close(negativeFromPolynomial,negative[mode],1e-10,"Dminus P2 is not the physical face-mass solution");
                ElectricContract::close(positiveFromPolynomial,positive[mode],1e-10,"Dplus P2 is not the physical face-mass solution");
            }
        }
    }
    check(fieldNorm[0]>1e-6&&fieldNorm[1]>1e-6,"charged benchmark lost a nonzero field in one material");
    const double negativeVolume=axis?4*pi*radius*radius*radius/3:pi*radius*radius;
    const double exactSource=harmonic?0:axis?-3*pi-3*negativeVolume:-2-2*negativeVolume;
    const double exactCharge=axis?4*pi*radius*radius*radius:2*pi*radius*radius;
    const double exactBoundary=harmonic?exactCharge:(axis?-3*pi:-2.);
    ElectricContract::close(totalSource,exactSource,1e-10,"charged bulk integral differs from independent circle/sphere volume");
    ElectricContract::close(totalCharge,exactCharge,1e-10,"charged solve used an incorrect independent total Q");
    check(std::abs(boundary-totalSource-totalCharge)<=1e-10*(1+boundaryScale+sourceScale+std::abs(totalCharge)),
        "charged whole-domain physical Gauss balance failed");
    ElectricContract::close(boundary,exactBoundary,1e-10,"charged boundary D integral differs from independent analytic PDE Gauss reference");
    const Error error{std::sqrt(valueError/valueNorm),std::sqrt(gradientError/gradientNorm)};
    if(harmonic)check(error.value<=.05&&error.gradient<=.15,"charged dipole-plus-radial harmonic field failed frozen physical L2 accuracy");
    else check(maximumValue<=2e-7&&maximumGradient<=2e-7,"charged P2 transmission MMS failed independent potential/electric-field accuracy");
    std::cout<<"ELECTRIC_HHO_P2_SOLVE,"<<axis<<",harmonic="<<harmonic<<",V_L2="<<error.value
        <<",E_L2="<<error.gradient<<",V_max="<<maximumValue<<",E_max="<<maximumGradient
        <<",balance="<<maximumBalance<<",boundary_D="<<boundary<<",source="<<totalSource
        <<",total_Q="<<totalCharge<<",PASS\n";
    if(retained)*retained=result;
    return error;
}

bool sameChargeClock(const I::HermiteGalsTimeLayer& a,const I::HermiteGalsTimeLayer& b) {
    return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;
}

double maxwellDot(Vector2D a,Vector2D b) {return a.x*b.x+a.y*b.y;}
Vector2D maxwellTensor(Vector2D electric,double epsilon,Vector2D normal) {
    return electric*(epsilon*maxwellDot(electric,normal))
        -normal*(.5*epsilon*maxwellDot(electric,electric));
}
Vector2D maxwellVelocity(Vector2D point) {
    // Prescribed nonuniform mechanical test velocity, independent of the field.
    return {.4+.3*point.x-.2*point.y+.11*point.x*point.y,
        -.2+.17*point.x+.4*point.y+.07*point.x*point.x};
}
Reference sphericalChargedExact(Vector2D point,bool liquid) {
    const auto q=point-origin(true);const double squared=maxwellDot(q,q),rho=std::sqrt(squared);
    constexpr double A=3./7.;
    if(liquid)return {A*q.y,{0,A}};
    const double C=-4*radius*radius*radius/7,monopole=radius*radius*radius;
    return {q.y+C*q.y/(squared*rho)+monopole/rho-radius*radius,
        {-3*C*q.y*q.x/(squared*squared*rho)-monopole*q.x/(squared*rho),
        1+C/(squared*rho)-3*C*q.y*q.y/(squared*squared*rho)-monopole*q.y/(squared*rho)}};
}
ElectricPhaseProblem sphericalChargedProblem() {
    ElectricPhaseProblem input;input.epsilonByPhase={{epsilonNegative,epsilonPositive}};
    for(int side=0;side<4;++side) {
        const bool pole=side==0;
        input.boundaries[side].type=pole?ElectricCartesianBoundary::Type::Neumann:ElectricCartesianBoundary::Type::Dirichlet;
        if(pole)input.boundaries[side].value=[](Vector2D){return 0.;};
        else input.boundaries[side].value=[](Vector2D p){return sphericalChargedExact(p,false).value;};
    }
    return input;
}
struct MaxwellReference {
    std::array<double,2> displacement{};
    double tangential=0,normalTraction=0,tangentialTraction=0;
    Vector2D traction;
};
MaxwellReference maxwellReference(Vector2D point,bool axis,bool harmonic) {
    const auto normal=(point-origin(axis))*(1/radius);const Vector2D tangent{-normal.y,normal.x};
    const auto minus=(axis&&harmonic?sphericalChargedExact(point,true):chargedExact(point,axis,true,harmonic)).gradient*(-1);
    const auto plus=(axis&&harmonic?sphericalChargedExact(point,false):chargedExact(point,axis,false,harmonic)).gradient*(-1);
    MaxwellReference reference;
    reference.displacement={{epsilonNegative*maxwellDot(minus,normal),epsilonPositive*maxwellDot(plus,normal)}};
    reference.tangential=maxwellDot(minus,tangent);
    // Full continuous tensors are the independent physical oracle.
    reference.traction=maxwellTensor(plus,epsilonPositive,normal)-maxwellTensor(minus,epsilonNegative,normal);
    reference.normalTraction=maxwellDot(reference.traction,normal);
    reference.tangentialTraction=maxwellDot(reference.traction,tangent);
    return reference;
}
struct MaxwellError {
    // Dminus, Dplus, Et, normal traction, tangent traction, full vector.
    std::array<double,6> relative{};
    double force=0,power=0,cubic=0,kkt=0;
    bool physicalBudgetPassed=true;
};
double maxwellRelative(double defect,double scale,const char* reason) {
    check(std::isfinite(defect)&&std::isfinite(scale)&&scale>=0,reason);
    if(scale==0){check(defect==0,reason);return 0;}
    return defect/scale;
}
Eigen::Vector4d maxwellBasis(double z) {
    return {1,z,.5*(3*z*z-1),.5*(5*z*z*z-3*z)};
}
double maxwellKkt(const ElectricPhaseResult& field,const ElectricMaxwellTrace& maxwell,std::size_t faceId) {
    const auto& geometry=field.geometry();const auto& face=geometry.faces().at(faceId);
    std::array<std::size_t,2> element{{geometry.elements().size(),geometry.elements().size()}};
    for(std::size_t e=0;e<geometry.elements().size();++e)
        if(std::find(geometry.elements()[e].faces.begin(),geometry.elements()[e].faces.end(),
            static_cast<int>(faceId))!=geometry.elements()[e].faces.end())
            element[geometry.elements()[e].liquid?0:1]=e;
    check(element[0]<geometry.elements().size()&&element[1]<geometry.elements().size(),
        "P5 independent KKT lacks its actual two-sided phase incidence");
    const auto minus=field.reconstructedPotentialP3(element[0]),plus=field.reconstructedPotentialP3(element[1]);
    Eigen::Matrix4d mass=Eigen::Matrix4d::Zero();Eigen::Vector4d target=Eigen::Vector4d::Zero();
    double potentialFamily=0;for(double c:minus)potentialFamily+=std::abs(c);for(double c:plus)potentialFamily+=std::abs(c);
    for(const auto& node:face.nodes) {
        const auto basis=maxwellBasis(node.coordinate);const double weight=node.weight*(geometry.axisymmetric()?2*pi*node.point.x:1);
        const double pMinus=ElectricContract::polynomial(node.point,field.elementCenter(element[0]),field.elementScale(element[0]),minus).v;
        const double pPlus=ElectricContract::polynomial(node.point,field.elementCenter(element[1]),field.elementScale(element[1]),plus).v;
        mass+=weight*basis*basis.transpose();target+=weight*.5*(pMinus+pPlus)*basis;
    }
    const auto trace=field.skeletalFaceTrace(faceId);
    Eigen::Vector3d shared(trace[0],trace[1],trace[2]);
    const Eigen::Matrix<double,3,4> constraint=mass.topRows<3>();
    Eigen::Matrix<double,7,7> system=Eigen::Matrix<double,7,7>::Zero();
    system.topLeftCorner<4,4>()=mass;system.topRightCorner<4,3>()=constraint.transpose();
    system.bottomLeftCorner<3,4>()=constraint;
    Eigen::Matrix<double,7,1> rhs;rhs.head<4>()=target;rhs.tail<3>()=mass.topLeftCorner<3,3>()*shared;
    // Independent complete constrained solve, not the production one-mode formula.
    Eigen::FullPivLU<Eigen::Matrix<double,7,7>> factor(system);
    check(factor.isInvertible(),"P5 weighted cubic KKT has no unique constrained minimizer");
    const Eigen::Matrix<double,7,1> reference=factor.solve(rhs);
    check((system*reference-rhs).norm()<=1e-11*(system.norm()*reference.norm()+rhs.norm()),
        "P5 independent cubic KKT did not solve its original weighted equations");
    const auto coefficient=maxwell.potentialPolynomialP3(faceId);
    Eigen::Vector4d actual(coefficient[0],coefficient[1],coefficient[2],coefficient[3]);
    const double coefficientScale=potentialFamily+shared.norm()+actual.norm()+reference.head<4>().norm();
    const double error=maxwellRelative((actual-reference.head<4>()).norm(),coefficientScale,
        "P5 nonfinite common cubic coefficient comparison");
    const double equationScale=mass.norm()*coefficientScale+target.norm()+rhs.tail<3>().norm();
    const double projection=maxwellRelative((constraint*actual-rhs.tail<3>()).norm(),equationScale,
        "P5 nonfinite weighted common-trace projection");
    const double optimality=maxwellRelative((mass*actual+constraint.transpose()*reference.tail<3>()-target).norm(),
        equationScale,"P5 nonfinite constrained common-trace optimality");
    check(std::max({error,projection,optimality})<=1e-10,
        "P5 common cubic differs from independent KKT/projection/optimality");
    return std::max({error,projection,optimality});
}
std::pair<Vector2D,double> maxwellMechanicalReference(bool axis) {
    Vector2D force;double power=0;
    constexpr int pieces=32;const double begin=axis?-.5*pi:0,end=axis?.5*pi:2*pi;
    for(int piece=0;piece<pieces;++piece)for(int k=0;k<8;++k) {
        const double a=begin+(end-begin)*piece/pieces,b=begin+(end-begin)*(piece+1)/pieces;
        const double theta=.5*(a+b)+.5*(b-a)*ElectricContract::nodes[k];
        const auto point=origin(axis)+Vector2D{radius*std::cos(theta),radius*std::sin(theta)};
        const double weight=.5*(b-a)*ElectricContract::weights[k]*radius*(axis?2*pi*point.x:1);
        const auto reference=maxwellReference(point,axis,true);
        force+=reference.traction*weight;power+=weight*maxwellDot(reference.traction,maxwellVelocity(point));
    }
    if(axis) {
        // Only the axial component is a genuine Cartesian net force in RZ.
        const double exactAxial=-4*pi*radius*radius*radius;
        check(std::abs(force.y-exactAxial)<=1e-11*std::abs(exactAxial),
            "P5 independent spherical tensor force contradicts its analytic axial net force");
    } else {
        const Vector2D exactForce{-2*pi*radius*radius,0};
        check(std::hypot(force.x-exactForce.x,force.y-exactForce.y)<=1e-11*std::abs(exactForce.x),
            "P5 independent continuous full-tensor force contradicts its analytic net force");
    }
    check(std::abs(power)>1e-6,"P5 independent mechanical-power reference became zero by symmetry");
    return {force,power};
}
MaxwellError maxwellChecks(const QuadtreeMesh& mesh,const ElectricPhaseResult& field,bool harmonic,
    bool fine=false,bool checkInvalid=false) {
    const auto oldSkeleton=field.skeletonSolution();
    const auto oldCharge=field.surfaceChargeLoad()?field.surfaceChargeLoad()->panelCharge():std::vector<double>{};
    const auto maxwell=ElectricMaxwellTrace::prepare(field);const auto& geometry=field.geometry();const bool axis=geometry.axisymmetric();
    check(&maxwell.geometry().publication().state().geometry().snapshot()
        ==&geometry.publication().state().geometry().snapshot(),"P5 Maxwell introduced a different geometry owner");
    check(sameChargeClock(maxwell.timeLayer(),field.timeLayer())
        &&sameChargeClock(maxwell.geometrySourceTimeLayer(),geometry.timeLayer()),
        "P5 Maxwell lost the actual field clock or immutable geometry source clock");
    maxwell.validateSource(mesh,geometry.publication(),field.timeLayer());
    check(maxwell.interfaceFaces().size()==geometry.report().interfaceFaces,"P5 Maxwell omitted an actual Gamma face");
    std::array<long double,6> error{},norm{};Vector2D force;long double power=0;
    double maximumD=0,maximumEt=0,maximumTangentTraction=0,maximumCharge=0,maximumKkt=0,maximumCubic=0,maximumAlgebra=0;
    for(std::size_t faceId:maxwell.interfaceFaces()) {
        const auto& face=geometry.faces().at(faceId);const auto& samples=maxwell.quadrature(faceId);
        check(face.kind==ElectricPhaseSkeletonFace::Kind::Interface&&samples.size()==face.nodes.size(),
            "P5 Maxwell replaced its original Gamma quadrature/layout");
        maximumKkt=std::max(maximumKkt,maxwellKkt(field,maxwell,faceId));
        const auto coefficient=maxwell.potentialPolynomialP3(faceId);
        maximumCubic=std::max(maximumCubic,std::abs(coefficient[3]));
        const auto tangential=maxwell.tangentialFieldPolynomialP2(faceId);
        Eigen::Vector3d qMoment=Eigen::Vector3d::Zero(),negative=Eigen::Vector3d::Zero(),positive=Eigen::Vector3d::Zero();
        Vector2D integrated;long double facePower=0;
        for(std::size_t j=0;j<samples.size();++j) {
            const auto& sample=samples[j];const auto& node=face.nodes[j];
            const double expectedWeight=node.weight*(axis?2*pi*node.point.x:1);
            check(sample.point.x==node.point.x&&sample.point.y==node.point.y&&sample.coordinate==node.coordinate
                &&sample.physicalWeight==expectedWeight,"P5 Maxwell sampled a second geometry or physical measure");
            ElectricContract::close(sample.normal.x,node.outwardNormal.x,1e-13,"P5 Maxwell normal lost source orientation");
            ElectricContract::close(sample.normal.y,node.outwardNormal.y,1e-13,"P5 Maxwell normal lost source orientation");
            ElectricContract::close(sample.tangent.x,-node.outwardNormal.y,1e-13,"P5 Maxwell tangent has the wrong native arclength orientation");
            ElectricContract::close(sample.tangent.y,node.outwardNormal.x,1e-13,"P5 Maxwell tangent has the wrong native arclength orientation");
            const auto basis=maxwellBasis(node.coordinate);
            const double derivative=coefficient[1]+3*node.coordinate*coefficient[2]
                +.5*(15*node.coordinate*node.coordinate-3)*coefficient[3];
            ElectricContract::close(sample.tangentialField,-2*derivative/face.length,1e-11,
                "P5 Et is not the arclength derivative of the common cubic potential");
            ElectricContract::close(sample.tangentialField,tangential[0]+tangential[1]*basis[1]+tangential[2]*basis[2],1e-11,
                "P5 exported Et polynomial differs from the same common-trace derivative");
            const auto minus=sample.normal*(sample.normalDisplacement[0]/epsilonNegative)+sample.tangent*sample.tangentialField;
            const auto plus=sample.normal*(sample.normalDisplacement[1]/epsilonPositive)+sample.tangent*sample.tangentialField;
            const auto tensor=maxwellTensor(plus,epsilonPositive,sample.normal)-maxwellTensor(minus,epsilonNegative,sample.normal);
            const double algebra=std::hypot(sample.tractionJump.x-tensor.x,sample.tractionJump.y-tensor.y);
            const double tensorScale=std::hypot(tensor.x,tensor.y)+epsilonNegative*maxwellDot(minus,minus)
                +epsilonPositive*maxwellDot(plus,plus);
            maximumAlgebra=std::max(maximumAlgebra,maxwellRelative(algebra,tensorScale,"P5 nonfinite full Maxwell tensor consistency"));
            check(algebra<=1e-10*tensorScale,"P5 traction is not the full two-sided Maxwell tensor jump");
            const auto velocity=maxwellVelocity(node.point);
            const double work=maxwellDot(sample.tractionJump,velocity);
            const double projected=sample.normalTraction*maxwellDot(sample.normal,velocity)
                +sample.tangentialTraction*maxwellDot(sample.tangent,velocity);
            check(std::abs(work-projected)<=1e-10*(std::abs(work)+std::abs(projected)+tensorScale*std::hypot(velocity.x,velocity.y)),
                "P5 Maxwell vector traction and mechanical power disagree");
            const auto reference=maxwellReference(node.point,axis,harmonic);
            const double actual[5]={sample.normalDisplacement[0],sample.normalDisplacement[1],sample.tangentialField,
                sample.normalTraction,sample.tangentialTraction};
            const double exact[5]={reference.displacement[0],reference.displacement[1],reference.tangential,
                reference.normalTraction,reference.tangentialTraction};
            for(int quantity=0;quantity<5;++quantity) {
                const double difference=actual[quantity]-exact[quantity];
                error[quantity]+=expectedWeight*difference*difference;norm[quantity]+=expectedWeight*exact[quantity]*exact[quantity];
            }
            const auto difference=sample.tractionJump-reference.traction;
            error[5]+=expectedWeight*maxwellDot(difference,difference);norm[5]+=expectedWeight*maxwellDot(reference.traction,reference.traction);
            maximumD=std::max({maximumD,std::abs(actual[0]-exact[0]),std::abs(actual[1]-exact[1])});
            maximumEt=std::max(maximumEt,std::abs(sample.tangentialField));
            maximumTangentTraction=std::max(maximumTangentTraction,std::abs(sample.tangentialTraction));
            maximumCharge=std::max(maximumCharge,std::abs(sample.chargeDensity-radius));
            qMoment+=expectedWeight*sample.chargeDensity*basis.head<3>();
            negative+=expectedWeight*sample.normalDisplacement[0]*basis.head<3>();
            positive+=expectedWeight*sample.normalDisplacement[1]*basis.head<3>();
            integrated+=sample.tractionJump*expectedWeight;facePower+=expectedWeight*work;
        }
        const auto q=field.interfaceChargeMoments(faceId),dnMinus=field.normalDisplacementMoments(faceId,true),
            dnPlus=field.normalDisplacementMoments(faceId,false);
        const Eigen::Vector3d qReference(q[0],q[1],q[2]),minusReference(dnMinus[0],dnMinus[1],dnMinus[2]),
            plusReference(dnPlus[0],dnPlus[1],dnPlus[2]);
        const double momentScale=qReference.norm()+minusReference.norm()+plusReference.norm();
        check((qMoment-qReference).norm()<=1e-10*momentScale
            &&(negative-minusReference).norm()<=1e-10*momentScale
            &&(positive-plusReference).norm()<=1e-10*momentScale
            &&(positive-negative-qMoment).norm()<=1e-10*momentScale,
            "P5 Maxwell q/Dn moments no longer come from the same original charge/flux");
        const auto receipt=maxwell.integratedTractionJump(faceId);
        check(std::hypot(receipt.x-integrated.x,receipt.y-integrated.y)<=1e-10*
            (std::hypot(receipt.x,receipt.y)+std::hypot(integrated.x,integrated.y)),
            "P5 integrated traction receipt differs from the actual physical quadrature");
        force+=integrated;power+=facePower;
    }
    MaxwellError result;result.kkt=maximumKkt;result.cubic=maximumCubic;
    for(int quantity=0;quantity<6;++quantity) {
        // The radial polynomial MMS has exactly zero Et and tangent traction.
        // Their meaningful error is an absolute maximum, not roundoff/roundoff.
        if(!harmonic&&(quantity==2||quantity==4))result.relative[quantity]=0;
        else result.relative[quantity]=norm[quantity]>0?std::sqrt(static_cast<double>(error[quantity]/norm[quantity])):0;
    }
    check(maximumCharge<=2e-7,"P5 derived physical q_h differs from the constant independent surface charge");
    if(harmonic) {
        check(maximumEt>1e-3&&maximumTangentTraction>1e-4&&maximumCubic>1e-10,
            "P5 harmonic Maxwell fixture lost nonzero Et/charge tangential traction/cubic mode");
        const auto mechanical=maxwellMechanicalReference(axis);
        result.force=axis?std::abs(force.y-mechanical.first.y)/std::abs(mechanical.first.y):
            std::hypot(force.x-mechanical.first.x,force.y-mechanical.first.y)/
                std::hypot(mechanical.first.x,mechanical.first.y);
        result.power=std::abs(static_cast<double>(power)-mechanical.second)/std::abs(mechanical.second);
        if(axis) {
            result.physicalBudgetPassed=result.relative[0]<=.15&&result.relative[1]<=.15&&result.relative[2]<=.15
                &&result.relative[4]<=.15&&result.relative[3]<=.30&&result.relative[5]<=.30&&result.force<=.30;
        } else if(fine) {
            result.physicalBudgetPassed=result.relative[0]<=.05&&result.relative[1]<=.05&&result.relative[2]<=.05
                &&result.relative[4]<=.05&&result.relative[3]<=.15&&result.relative[5]<=.15&&result.force<=.15&&result.power<=.15;
        }
    } else {
        check(maximumD<=2e-7&&maximumEt<=2e-7&&maximumTangentTraction<=2e-7&&result.relative[3]<=2e-6,
            "P5 nonzero-charge XY/RZ polynomial Maxwell transmission MMS failed");
    }
    if(checkInvalid) {
        ElectricMaxwellControls invalid;invalid.relativeProjectionTolerance=0;
        rejected([&]{ElectricMaxwellTrace::prepare(field,invalid);},"P5 accepted an invalid numerical projection tolerance");
        invalid={};invalid.relativeChargeJumpTolerance=1e-8;
        rejected([&]{ElectricMaxwellTrace::prepare(field,invalid);},"P5 relaxed the independent physical charge-jump target");
        auto wrongClock=field.timeLayer();wrongClock.end+=.01;
        rejected([&]{maxwell.validateSource(mesh,geometry.publication(),wrongClock);},"P5 accepted the wrong electric physical clock");
        const auto foreign=publication(mesh,axis,geometry.timeLayer());
        rejected([&]{maxwell.validateSource(mesh,foreign,field.timeLayer());},"P5 accepted a foreign equal-sampled geometry owner");
        const auto cartesian=std::find_if(geometry.faces().begin(),geometry.faces().end(),
            [](const auto& face){return face.kind!=ElectricPhaseSkeletonFace::Kind::Interface;});
        check(cartesian!=geometry.faces().end(),"P5 invalid face probe has no Cartesian source face");
        const auto faceId=static_cast<std::size_t>(std::distance(geometry.faces().begin(),cartesian));
        rejected([&]{maxwell.quadrature(faceId);},"P5 accepted a noninterface Maxwell query");
    }
    check(field.skeletonSolution()==oldSkeleton&&(!field.surfaceChargeLoad()
        ||identicalCharge(field.surfaceChargeLoad()->panelCharge(),oldCharge)),
        "P5 Maxwell preparation/queries changed the immutable electric result or Q receipt");
    std::cout<<"ELECTRIC_HHO_P5_MAXWELL,"<<axis<<",harmonic="<<harmonic<<",fine="<<fine
        <<",Dn_minus="<<result.relative[0]<<",Dn_plus="<<result.relative[1]<<",Et="<<result.relative[2]
        <<",tau_n="<<result.relative[3]<<",tau_t="<<result.relative[4]<<",vector="<<result.relative[5]
        <<",Et_abs_max="<<maximumEt<<",tau_t_abs_max="<<maximumTangentTraction
        <<",force="<<result.force<<",power="<<result.power<<",cubic="<<result.cubic<<",KKT="<<result.kkt
        <<",tensor="<<maximumAlgebra<<",integrated=("<<force.x<<','<<force.y<<"),work="<<static_cast<double>(power)
        <<(result.physicalBudgetPassed?",PASS\n":",FAIL_PHYSICAL_BUDGET\n");
    return result;
}
void maxwellChargeCandidate(const QuadtreeMesh& mesh,const ElectricFixedGeometryChargeStep& step,const SurfaceChargeField& field) {
    const auto accepted=field.transport().panelCharge(),candidate=step.integratedCharge();
    const auto original=step.potential().skeletonSolution();
    const auto maxwell=ElectricMaxwellTrace::prepare(step.potential());
    check(sameChargeClock(maxwell.timeLayer(),step.timeLayer())
        &&sameChargeClock(maxwell.geometrySourceTimeLayer(),step.geometrySourceTimeLayer())
        &&!sameChargeClock(maxwell.timeLayer(),maxwell.geometrySourceTimeLayer()),
        "P5 fixed-geometry candidate lost its true separate physical and source clocks");
    const auto& geometry=step.potential().geometry();
    maxwell.validateSource(mesh,geometry.publication(),step.timeLayer());
    check(identicalCharge(field.transport().panelCharge(),accepted)&&identicalCharge(step.integratedCharge(),candidate)
        &&step.potential().skeletonSolution()==original,"P5 candidate Maxwell changed accepted/candidate Q or the electric field");
}

I::HermiteGalsTimeLayer nextChargeClock(const ElectricSurfaceChargeTimeState& state,double dt) {
    const auto old=state.timeLayer();const double begin=old.time();
    return {old.step+1,begin,begin+dt,I::HermiteGalsTimeLayer::Position::End};
}
void chargeTimeEquation(const ElectricFixedGeometryChargeStep& step,const std::vector<double>& old,
    const std::vector<double>* previous) {
    const auto& report=step.report();const auto& charge=step.integratedCharge();
    const auto& current=step.current().panelCurrent();
    check(charge.size()==old.size()&&current.size()==charge.size(),"P4 final charge/current inventory changed");
    check(report.relativeSchurResidual<=1e-11&&report.maxRelativePanelEquationResidual<=1e-10
        &&report.relativeTotalEquationResidual<=1e-10,"P4 final reported charge equation gate failed");
    long double total=0,scale=0;double maximum=0;
    for(std::size_t i=0;i<charge.size();++i) {
        const double history=previous?report.a2*(*previous)[i]:0;
        const double residual=report.a0*charge[i]+report.a1*old[i]+history-report.dt*current[i];
        const double localScale=std::abs(report.a0*charge[i])+std::abs(report.a1*old[i])
            +std::abs(history)+std::abs(report.dt*current[i]);
        check(localScale>0||residual==0,"P4 nonzero charge equation has no physical scale");
        maximum=std::max(maximum,localScale>0?std::abs(residual)/localScale:0);
        total+=residual;scale+=localScale;
    }
    check(maximum<=1e-10&&std::abs(total)<=1e-10*scale,"P4 actual final raw current does not close the time equation");
    const auto& electric=step.potential().report();
    check(electric.skeleton.relativeLinearResidual<=1e-11&&electric.skeleton.relativeSymmetryDefect<=1e-12
        &&electric.skeleton.minimumSkeletonPivot>0&&electric.maxRelativeInterfaceBalance<=1e-10
        &&electric.maxRelativeLocalGaussResidual<=1e-10&&electric.relativeGlobalGaussResidual<=1e-10,
        "P4 final charge solve omitted original electric/Gauss gates");
    check(step.potential().surfaceChargeLoad()
        &&identicalCharge(step.potential().surfaceChargeLoad()->panelCharge(),charge),
        "P4 final electric field belongs to another charge iterate");
    check(sameChargeClock(step.timeLayer(),step.potential().timeLayer())
        &&sameChargeClock(step.timeLayer(),step.current().timeLayer()),
        "P4 final Q/V/raw current physical clocks differ");
    check(sameChargeClock(step.geometrySourceTimeLayer(),step.potential().geometry().timeLayer()),
        "P4 time advance retagged the immutable geometry source clock");
}
Eigen::MatrixXd freeChargeInjection(const ElectricPhaseResult& result,
    const ElectricSurfaceChargeMomentOperator& op,const ElectricPhaseProblem& input) {
    const Eigen::MatrixXd physical=op.interfaceInjectionMatrix()*op.projectionMatrix();
    const int n=static_cast<int>(result.skeletonSolution().size());
    Eigen::MatrixXd injection=Eigen::MatrixXd::Zero(n,static_cast<int>(op.panelCount()));
    int offset=0;
    for(std::size_t face=0;face<result.geometry().faces().size();++face) {
        const int side=result.geometry().faces()[face].boundarySide;
        const bool prescribed=side>=0&&input.boundaries[static_cast<std::size_t>(side)].type
            ==ElectricCartesianBoundary::Type::Dirichlet;
        if(!prescribed) {
            check(offset+3<=n,"P4 reference charge injection guessed a foreign free layout");
            injection.middleRows(offset,3)=physical.middleRows(3*static_cast<Eigen::Index>(face),3);offset+=3;
        }
    }
    check(offset==n,"P4 reference injection did not consume actual free skeleton inventory");
    return injection;
}
void chargeTimeContract(const QuadtreeMesh& mesh,const ElectricPhaseResult& base,const ElectricSurfaceChargeMomentOperator& op,
    SurfaceChargeField& field,bool harmonic,bool withMaxwell=false) {
    static_assert(!std::is_copy_constructible_v<ElectricSurfaceChargeTimeState>);
    const auto saved=field.transport().panelCharge();const auto sourceClock=base.geometry().timeLayer();
    const auto sourceIdentity=base.geometry().geometryIdentity();
    auto signedCharge=saved;bool positive=false,negative=false;
    for(std::size_t i=0;i<signedCharge.size();++i) {
        signedCharge[i]*=.2+std::sin(.73*(i+1));
        positive=positive||signedCharge[i]>0;negative=negative||signedCharge[i]<0;
    }
    check(positive&&negative,"P4 retained-panel charge fixture has no signed nonuniform data");
    constexpr double rate=.14;
    const auto scalarSolver=ElectricFixedGeometryChargeSolver::prepare(base,op,{{5*rate,rate}});
    double maximumStiffError=0;
    for(double z:{.1,1.,10.,100.}) {
        field.assignIntegratedCharge(signedCharge);
        auto state=ElectricSurfaceChargeTimeState::initialize(scalarSolver,{&field},sourceClock);
        const auto old=state.integratedCharge();const auto target=nextChargeClock(state,z/rate);
        const auto step=scalarSolver.advance(state,target);
        check(step.report().temporalOrder==1,"P4 first physical step did not use BE startup");
        chargeTimeEquation(step,old,nullptr);
        if(withMaxwell)maxwellChargeCandidate(mesh,step,field);
        for(std::size_t i=0;i<old.size();++i) {
            const double reference=old[i]/(1+z);
            const double error=std::abs(step.integratedCharge()[i]-reference)/(1+std::abs(reference));
            maximumStiffError=std::max(maximumStiffError,error);
            check(error<=1e-10,"P4 large-step BE differs from independent dielectric-relaxation reference");
        }
        check(identicalCharge(field.transport().panelCharge(),old)&&sameChargeClock(state.timeLayer(),sourceClock),
            "P4 solving a candidate changed the accepted Field/history");
        state.commit(step);
        check(identicalCharge(field.transport().panelCharge(),step.integratedCharge())
            &&sameChargeClock(state.timeLayer(),target)&&state.previousTimeLayer()
            &&sameChargeClock(*state.previousTimeLayer(),sourceClock)&&state.previousIntegratedCharge()
            &&identicalCharge(*state.previousIntegratedCharge(),old),
            "P4 commit did not atomically retain actual accepted Q history and full clock");
        rejected([&]{state.commit(step);},"P4 accepted a stale/double commit");
        const auto oldSecond=state.integratedCharge(),older=*state.previousIntegratedCharge();
        const auto secondTarget=nextChargeClock(state,.5*z/rate);
        const auto second=scalarSolver.advance(state,secondTarget);
        chargeTimeEquation(second,oldSecond,&older);
        if(withMaxwell)maxwellChargeCandidate(mesh,second,field);
        const double ratio=.5,a0=(1+2*ratio)/(1+ratio),a1=-(1+ratio),a2=ratio*ratio/(1+ratio);
        check(second.report().temporalOrder==2,"P4 accepted two physical histories but did not use variable BDF2");
        ElectricContract::close(second.report().a0,a0,1e-14,"P4 BDF2 a0 ignores the real previous dt");
        ElectricContract::close(second.report().a1,a1,1e-14,"P4 BDF2 a1 ignores the real previous dt");
        ElectricContract::close(second.report().a2,a2,1e-14,"P4 BDF2 a2 ignores the real previous dt");
        for(std::size_t i=0;i<old.size();++i) {
            const double reference=(-a1*oldSecond[i]-a2*older[i])/(a0+.5*z);
            ElectricContract::close(second.integratedCharge()[i],reference,1e-10,
                "P4 stiff variable BDF2 differs from independent two-history formula");
        }
        state.commit(second);
        auto badRatio=nextChargeClock(state,3*second.report().dt);
        rejected([&]{scalarSolver.advance(state,badRatio);},"P4 accepted an unsupported BDF2 history step ratio");
    }
    field.assignIntegratedCharge(signedCharge);
    const auto unequal=ElectricFixedGeometryChargeSolver::prepare(base,op,{{.7,.03}});
    auto state=ElectricSurfaceChargeTimeState::initialize(unequal,{&field},sourceClock);
    rejected([&]{ElectricSurfaceChargeTimeState::initialize(unequal,{&field},sourceClock);},
        "P4 bound two live accepted time owners to one mutable Field");
    const auto input=chargedProblem(base.axisymmetric(),harmonic);
    const auto zero=solveElectricPhaseElectrostatics(base.geometry(),input,op.apply(std::vector<double>(op.panelCount(),0)));
    const auto rawCurrent=zero.preparePanelCurrent(op,{{.7,.03}});
    const auto injection=freeChargeInjection(zero,op,input);
    const Eigen::MatrixXd skeleton=zero.skeletonMatrix(),currentMap=rawCurrent.matrix();
    const int nu=static_cast<int>(zero.skeletonSolution().size()),nq=static_cast<int>(op.panelCount());
    const Eigen::Map<const Eigen::VectorXd> staticRhs(zero.skeletonRhs().data(),nu);
    const Eigen::Map<const Eigen::VectorXd> currentOffset(rawCurrent.offset().data(),nq);
    double maximumBlockError=0;
    for(double dt:{10/rate,5/rate}) {
        const auto old=state.integratedCharge();
        const auto older=state.previousIntegratedCharge()?*state.previousIntegratedCharge():std::vector<double>{};
        const auto target=nextChargeClock(state,dt);
        const auto step=unequal.advance(state,target);
        chargeTimeEquation(step,old,older.empty()?nullptr:&older);
        if(withMaxwell)maxwellChargeCandidate(mesh,step,field);
        const auto& report=step.report();
        Eigen::MatrixXd block=Eigen::MatrixXd::Zero(nu+nq,nu+nq);
        block.topLeftCorner(nu,nu)=skeleton;block.topRightCorner(nu,nq)=-injection;
        block.bottomLeftCorner(nq,nu)=-dt*currentMap;
        block.bottomRightCorner(nq,nq)=report.a0*Eigen::MatrixXd::Identity(nq,nq);
        Eigen::VectorXd rhs(nu+nq);rhs.head(nu)=staticRhs;
        for(int i=0;i<nq;++i)rhs[nu+i]=-report.a1*old[static_cast<std::size_t>(i)]
            -(older.empty()?0:report.a2*older[static_cast<std::size_t>(i)])+dt*currentOffset[i];
        Eigen::FullPivLU<Eigen::MatrixXd> direct(block);
        check(direct.isInvertible(),"P4 independent coupled reference has no unique solution");
        const Eigen::VectorXd reference=direct.solve(rhs);
        check((block*reference-rhs).norm()<=1e-11*(rhs.norm()+block.norm()*reference.norm()),
            "P4 independent coupled block reference failed its actual equations");
        for(int i=0;i<nq;++i) {
            const double error=std::abs(step.integratedCharge()[static_cast<std::size_t>(i)]-reference[nu+i])
                /(1+std::abs(reference[nu+i]));
            maximumBlockError=std::max(maximumBlockError,error);
            check(error<=1e-9,"P4 matrix-free Schur differs from the independent coupled electric/charge solve");
        }
        for(int i=0;i<nu;++i)ElectricContract::close(step.potential().skeletonSolution()[static_cast<std::size_t>(i)],
            reference[i],1e-9,"P4 final coeval potential differs from independent coupled block solution");
        state.commit(step);
    }
    const auto accepted=state.integratedCharge();const auto acceptedClock=state.timeLayer();
    auto badClock=nextChargeClock(state,state.timeLayer().end-state.timeLayer().begin);
    badClock.begin+=.01;
    rejected([&]{unequal.advance(state,badClock);},"P4 accepted a disconnected physical clock interval");
    auto badPosition=nextChargeClock(state,state.timeLayer().end-state.timeLayer().begin);
    badPosition.position=I::HermiteGalsTimeLayer::Position::Midpoint;
    rejected([&]{unequal.advance(state,badPosition);},"P4 accepted a midpoint as a final charge step");
    ElectricChargeTimeControls badControls;badControls.relativeChargeBalanceTolerance=1e-8;
    rejected([&]{unequal.advance(state,nextChargeClock(state,1),badControls);},"P4 relaxed the physical charge gate");
    badControls={};badControls.maximumIterations=0;
    rejected([&]{unequal.advance(state,nextChargeClock(state,1),badControls);},"P4 accepted invalid solve controls");
    check(identicalCharge(field.transport().panelCharge(),accepted)&&sameChargeClock(state.timeLayer(),acceptedClock),
        "P4 rejected solve changed accepted Q or history clock");
    auto changed=accepted;changed[0]+=1e-4;field.assignIntegratedCharge(changed);
    rejected([&]{unequal.advance(state,nextChargeClock(state,1));},"P4 accepted a Field mutated outside its history transaction");
    check(identicalCharge(field.transport().panelCharge(),changed),"P4 external-mutation rejection changed Field data");
    field.assignIntegratedCharge(saved);
    rejected([&]{ElectricSurfaceChargeTimeState::initialize(unequal,{&field,&field},sourceClock);},
        "P4 bound duplicate/foreign component charge owners");
    check(base.geometry().geometryIdentity()==sourceIdentity
        &&sameChargeClock(base.geometry().timeLayer(),sourceClock),"P4 changed its fixed geometry source identity/clock");
    std::cout<<"ELECTRIC_HHO_P4_TIME,"<<base.axisymmetric()<<",stiff_z_max=100,BE_error="<<maximumStiffError
        <<",coupled_block_error="<<maximumBlockError<<",PASS\n";
}

void chargeContract(const QuadtreeMesh& mesh,const ElectricPhaseGeometry& geometry,
    const ElectricPhaseResult* baseline,bool harmonic,bool withCurrent=false,
    std::optional<ElectricSurfaceChargeMomentOperator>* previousMoment=nullptr,bool withTime=false,
    bool withMaxwell=false,MaxwellError* harmonicMaxwell=nullptr) {
    check(geometry.publication().topology().curves.size()==1,"declared P2 benchmark is not one actual interface component");
    SurfaceChargeField field(mesh,geometry.publication().topology().curves[0],radius);
    const auto saved=field.transport().panelCharge();
    const auto identity=geometry.geometryIdentity();
    const auto* owner=&geometry.publication().state().geometry().snapshot();
    const auto op=ElectricSurfaceChargeMomentOperator::prepare(geometry,{&field},geometry.timeLayer());
    chargeMomentChecks(geometry,field,op,!harmonic);
    if(baseline)zeroChargeRegression(geometry,op,*baseline);
    std::optional<ElectricPhaseResult> retained;
    chargedSolveAndCheck(mesh,geometry,op,saved,harmonic,(withCurrent||withMaxwell)?&retained:nullptr);
    if(withCurrent) {
        const auto* foreign=previousMoment&&previousMoment->has_value()?&**previousMoment:nullptr;
        currentChecks(*retained,op,field,harmonic,!harmonic&&!geometry.axisymmetric(),foreign);
        if(previousMoment&&!previousMoment->has_value())*previousMoment=op;
    }
    if(withMaxwell) {
        const auto error=maxwellChecks(mesh,*retained,harmonic,harmonic,!harmonic&&!geometry.axisymmetric());
        if(harmonic&&harmonicMaxwell)*harmonicMaxwell=error;
        if(!harmonic) {
            const auto extra=geometry.axisymmetric()?
                solveElectricPhaseElectrostatics(geometry,sphericalChargedProblem(),op.apply(saved)):
                solveElectricPhaseElectrostatics(geometry,chargedProblem(false,true),op.apply(saved));
            const auto harmonicError=maxwellChecks(mesh,extra,true);
            if(harmonicMaxwell)*harmonicMaxwell=harmonicError;
        }
    }
    check(identicalCharge(field.transport().panelCharge(),saved),"charged solve published trial Q to its original material owner");
    check(geometry.geometryIdentity()==identity&&&geometry.publication().state().geometry().snapshot()==owner,
        "charged solve created or changed the scalar geometry owner");
    if(withTime)chargeTimeContract(mesh,*retained,op,field,harmonic,withMaxwell);
}

// E-M1: true FirstJet/RK2 geometry publications and the production material
// endpoint/physical-overlap factory. The spherical reference below observes
// the transfer; it never supplies a weight, current or candidate charge.
I::HermiteGalsPublication chargeMovingPublication(const QuadtreeMesh& mesh,
    const I::HermiteGalsPublication& source,const I::HermiteGalsTimeLayer& target,
    double displacement,double totalDisplacement) {
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    const LevelSetBoundaryTrace::Domain domain{grid.minCoord(),grid.maxCoord(),grid.spacing(),true};
    I::HermiteGalsCandidate authority(mesh,target,source.state().primary());
    const double dt=target.time()-source.state().timeLayer().time();
    check(dt>0,"E-M1 geometry fixture has no actual positive interval");
    const auto velocity=[speed=displacement/dt](const Vector2D&) {
        return LevelSetBoundaryTrace::VelocityJet{{0,speed},{0,0,0,0}};
    };
    const I::HermiteGalsVelocityView previous(source.state().identity(),source.state().timeLayer(),velocity);
    const I::HermiteGalsVelocityView arrival(authority.identity(),target,velocity);
    const I::HermiteGalsBoundaryFirstJetQuery boundary=[target,totalDisplacement,&mesh](
        Vector2D p,I::HermiteGalsTimeLayer layer,Boundary) {
        check(sameChargeClock(layer,target),"E-M1 prescribed boundary lost its complete target clock");
        check(p.x==mesh.minCoord.x||p.x==mesh.maxCoord.x||p.y==mesh.minCoord.y||p.y==mesh.maxCoord.y,
            "E-M1 prescribed first jet was used for an interior reset");
        const Vector2D q{p.x,p.y-.5-totalDisplacement};
        return I::HermiteGalsBoundaryFirstJet{q.x*q.x+q.y*q.y-radius*radius,{2*q.x,2*q.y}};
    };
    const auto step=I::prepareHermiteGalsFirstJetStep(mesh,source,domain,previous,arrival,boundary);
    check(step.sourceIdentity()==source.state().identity()
        &&sameChargeClock(step.publication().state().timeLayer(),target),
        "E-M1 target is not the real FirstJet publication of the supplied source/motion");
    return step.publication();
}
double chargeIntegralScale(const std::vector<double>& values) {
    long double scale=0;for(double q:values)scale+=std::abs(q);return static_cast<double>(scale);
}
void chargeIntegralConserved(const std::vector<double>& old,const std::vector<double>& mapped) {
    long double before=0,after=0;for(double q:old)before+=q;for(double q:mapped)after+=q;
    check(std::abs(before-after)<=512*std::numeric_limits<double>::epsilon()*chargeIntegralScale(old),
        "E-M1 material remap changed signed integrated total charge");
}
std::vector<double> chargeSphericalOverlapReference(const SurfaceChargeMaterialTransfer& transfer,
    const std::vector<double>& charge,double displacement,double& maximumRelativeBoundaryError) {
    const auto& source=transfer.sourcePanels().panels();
    const auto& target=transfer.targetPanels().panels();
    const auto& boundaries=transfer.materialBoundaries();
    check(source.size()==charge.size()&&boundaries.size()==source.size()+1,
        "E-M1 actual material transfer does not cover every original boundary");
    const double sphereArea=4*pi*radius*radius;
    for(std::size_t i=0;i<boundaries.size();++i) {
        const double moved=(i==source.size()?source.back().end.y:source[i].begin.y)+displacement;
        const auto location=boundaries[i];double expected=0;
        if(location.panelIndex==target.size()) {
            check(location.localMeasure==0,"E-M1 endpoint sentinel has nonzero physical measure");
            ElectricContract::close(moved,target.back().end.y,1e-10,"E-M1 material map lost the final axis pole");
        } else {
            check(location.panelIndex<target.size(),"E-M1 completed material boundary has a foreign panel");
            const auto& panel=target[location.panelIndex];
            // For a sphere meridian, 2*pi*r*ds = 2*pi*R*dz. A pure axial
            // translation preserves material z-center; this reference uses
            // physical endpoint coordinates, independently of local map code.
            expected=2*pi*radius*(moved-panel.begin.y);
            const double error=std::abs(location.localMeasure-expected)/sphereArea;
            maximumRelativeBoundaryError=std::max(maximumRelativeBoundaryError,error);
            check(error<=2e-9,"E-M1 real endpoint map disagrees with independent RZ physical measure");
        }
    }
    std::vector<double> result(target.size(),0);
    for(std::size_t i=0;i<source.size();++i) {
        const double a=source[i].begin.y,b=source[i].end.y;
        check(b>a,"E-M1 sphere chain does not retain its positive bottom-to-top orientation");
        check(std::abs(source[i].physicalMeasure-2*pi*radius*(b-a))<=2e-9*sphereArea,
            "E-M1 source physical area is not the independent spherical RZ area");
        long double totalShare=0;
        for(std::size_t j=0;j<target.size();++j) {
            const double c=target[j].begin.y-displacement,d=target[j].end.y-displacement;
            const double overlap=std::min(b,d)-std::max(a,c);
            if(overlap>0) {const double share=overlap/(b-a);result[j]+=charge[i]*share;totalShare+=share;}
        }
        check(std::abs(totalShare-1)<=2e-9,"E-M1 independent sphere intervals do not cover a material donor");
    }
    return result;
}
double chargeMappedReferenceError(const std::vector<double>& actual,const std::vector<double>& reference,
    const std::vector<double>& source) {
    check(actual.size()==reference.size(),"E-M1 mapped panel count differs from real target geometry");
    long double error=0;for(std::size_t i=0;i<actual.size();++i)error+=std::abs(actual[i]-reference[i]);
    const double scale=chargeIntegralScale(source);
    check(scale>0,"E-M1 signed physical test charge has no nonzero integral scale");
    const double relative=static_cast<double>(error)/scale;
    check(relative<=2e-8,"E-M1 mapped Q differs from independent spherical material overlap");
    return relative;
}
void chargeMovingCoeval(const QuadtreeMesh& mesh,const ElectricMovingGeometryChargeStep& moving,
    const I::HermiteGalsPublication& target,const I::HermiteGalsTimeLayer& layer) {
    const auto& step=moving.step();const auto& potential=step.potential();
    const auto& op=step.current().momentOperator();
    check(moving.componentCount()==1&&moving.field(0).geometry().usesSource(target.state().geometry().snapshot())
        &&op.geometry().geometryIdentity()==target.state().identity()
        &&sameChargeClock(op.timeLayer(),layer)
        &&sameChargeClock(step.timeLayer(),layer)&&sameChargeClock(step.geometrySourceTimeLayer(),layer)
        &&sameChargeClock(moving.geometrySourceTimeLayer(),layer),
        "E-M1 candidate Q/P_Q/field/current is not represented on the same true target geometry/full clock");
    check(potential.surfaceChargeLoad()&&potential.surfaceChargeLoad()->momentOperator().geometry().geometryIdentity()
        ==op.geometry().geometryIdentity(),"E-M1 final potential uses a foreign target P_Q");
    const auto maxwell=ElectricMaxwellTrace::prepare(potential);
    maxwell.validateSource(mesh,target,layer);
    check(sameChargeClock(maxwell.timeLayer(),layer)&&sameChargeClock(maxwell.geometrySourceTimeLayer(),layer),
        "E-M1 Maxwell did not retain the actual coeval moving Q/V/current clock");
}

// E-C1 exercises real native publications driven by prescribed velocity,
// independently of Hydro load deposition or a coupled AC transaction.
struct ChargeStageReference {
    double chargeError=0,potentialError=0,currentError=0,blockResidual=0;
    double panelEquation=0,totalEquation=0,overlapError=0,boundaryError=0;
    double dt=0,previousDt=0,a0=0,a1=0,a2=0;
    std::size_t materialQueries=0;
    int rejections=0;
};
ChargeStageReference chargeStageBlockReference(const ElectricMovingGeometryChargeStep& moving,
    const ElectricPhaseProblem& electric,const std::array<double,2>& sigma) {
    const auto& step=moving.step();const auto& old=moving.mappedIntegratedCharge();
    const auto* previous=moving.mappedPreviousIntegratedCharge();
    const double h=step.timeLayer().time()-moving.sourceTimeLayer().time();
    const double d=moving.previousSourceTimeLayer()
        ?moving.sourceTimeLayer().time()-moving.previousSourceTimeLayer()->time():0;
    check(h>0&&(!previous||d>0),"E-C1 reference has invalid actual physical time nodes");
    ChargeStageReference result;result.dt=h;result.previousDt=d;
    // Derivative of the Lagrange interpolant through physical nodes
    // t_stage, t_n and t_(n-1), multiplied by h. No reported coefficient is used.
    result.a0=previous?(2*h+d)/(h+d):1;
    result.a1=previous?-(h+d)/d:-1;
    result.a2=previous?h*h/(d*(h+d)):0;
    const auto& report=step.report();
    ElectricContract::close(report.dt,h,1e-14,"E-C1 used the full AC interval instead of the actual stage increment");
    ElectricContract::close(report.previousDt,d,1e-14,"E-C1 stage history lost its original accepted physical increment");
    ElectricContract::close(report.a0,result.a0,1e-14,"E-C1 a0 differs from independent physical-node derivative");
    ElectricContract::close(report.a1,result.a1,1e-14,"E-C1 a1 differs from independent physical-node derivative");
    ElectricContract::close(report.a2,result.a2,1e-14,"E-C1 a2 differs from independent physical-node derivative");
    check(report.temporalOrder==(previous?2:1),"E-C1 stage discarded a real history or invented a startup history");
    const auto& op=step.current().momentOperator();
    const auto zero=solveElectricPhaseElectrostatics(step.potential().geometry(),electric,
        op.apply(std::vector<double>(op.panelCount(),0)));
    const auto current=zero.preparePanelCurrent(op,sigma);
    const Eigen::MatrixXd injection=freeChargeInjection(zero,op,electric);
    const int nu=static_cast<int>(zero.skeletonSolution().size()),nq=static_cast<int>(op.panelCount());
    Eigen::MatrixXd block=Eigen::MatrixXd::Zero(nu+nq,nu+nq);
    block.topLeftCorner(nu,nu)=Eigen::MatrixXd(zero.skeletonMatrix());
    block.topRightCorner(nu,nq)=-injection;
    block.bottomLeftCorner(nq,nu)=-h*Eigen::MatrixXd(current.matrix());
    block.bottomRightCorner(nq,nq)=result.a0*Eigen::MatrixXd::Identity(nq,nq);
    Eigen::VectorXd rhs(nu+nq);
    rhs.head(nu)=Eigen::Map<const Eigen::VectorXd>(zero.skeletonRhs().data(),nu);
    for(int i=0;i<nq;++i)rhs[nu+i]=-result.a1*old[static_cast<std::size_t>(i)]
        -(previous?result.a2*(*previous)[static_cast<std::size_t>(i)]:0)+h*current.offset()[static_cast<std::size_t>(i)];
    const Eigen::FullPivLU<Eigen::MatrixXd> direct(block);
    check(direct.isInvertible(),"E-C1 original coupled block reference has no unique solution");
    const Eigen::VectorXd reference=direct.solve(rhs);
    const double blockScale=rhs.norm()+block.norm()*reference.norm();
    result.blockResidual=(block*reference-rhs).norm()/blockScale;
    check(result.blockResidual<=1e-11,"E-C1 independent coupled reference does not solve its original equations");
    const Eigen::VectorXd referenceCurrent=Eigen::MatrixXd(current.matrix())*reference.head(nu)
        +Eigen::Map<const Eigen::VectorXd>(current.offset().data(),nq);
    long double equationTotal=0,equationScale=0,currentError=0,currentScale=0;
    for(int i=0;i<nq;++i) {
        const auto k=static_cast<std::size_t>(i);
        result.chargeError=std::max(result.chargeError,
            std::abs(step.integratedCharge()[k]-reference[nu+i])/(1+std::abs(reference[nu+i])));
        currentError+=std::abs(step.current().panelCurrent()[k]-referenceCurrent[i]);
        currentScale+=std::abs(step.current().panelCurrent()[k])+std::abs(referenceCurrent[i]);
        const double history=previous?result.a2*(*previous)[k]:0;
        const double residual=result.a0*step.integratedCharge()[k]+result.a1*old[k]+history-h*step.current().panelCurrent()[k];
        const double scale=std::abs(result.a0*step.integratedCharge()[k])+std::abs(result.a1*old[k])
            +std::abs(history)+std::abs(h*step.current().panelCurrent()[k]);
        result.panelEquation=std::max(result.panelEquation,scale?std::abs(residual)/scale:0);
        equationTotal+=residual;equationScale+=scale;
    }
    for(int i=0;i<nu;++i)result.potentialError=std::max(result.potentialError,
        std::abs(step.potential().skeletonSolution()[static_cast<std::size_t>(i)]-reference[i])/(1+std::abs(reference[i])));
    check(currentScale>0&&equationScale>0,"E-C1 signed nonuniform fixture lost its physical current/equation scale");
    result.currentError=static_cast<double>(currentError/currentScale);
    result.totalEquation=static_cast<double>(std::abs(equationTotal)/equationScale);
    check(result.chargeError<=1e-9&&result.potentialError<=1e-9&&result.currentError<=1e-9
        &&result.panelEquation<=1e-10&&result.totalEquation<=1e-10,
        "E-C1 Schur Q/V/current or charge-current balance differs from the independent original coupled block");
    check(step.current().conductivityByPhase()==sigma&&step.current().epsilonByPhase()==electric.epsilonByPhase
        &&step.current().report().relativeGlobalCurrentBalance<=1e-10,
        "E-C1 candidate reused another material/current preparation");
    chargeTimeEquation(step,old,previous);
    return result;
}
ElectricFixedGeometryChargeStep chargeACMidpointContract(const QuadtreeMesh& mesh,
    ElectricSurfaceChargeTimeState& state,const I::HermiteGalsPublication& source,
    I::HermiteGalsTimeLayer target,double displacement,double totalDisplacement,
    const ElectricPhaseProblem& electric,const std::array<double,2>& sigma,ChargeStageReference& metric) {
    target.position=I::HermiteGalsTimeLayer::Position::Midpoint;
    const auto old=state.integratedCharge();
    const auto previous=state.previousIntegratedCharge()?*state.previousIntegratedCharge():std::vector<double>{};
    const auto acceptedClock=state.timeLayer(),representationClock=state.geometrySourceTimeLayer();
    const auto* acceptedField=&state.field(0);const auto* acceptedSolver=&state.solver();
    const std::optional<I::HermiteGalsTimeLayer> previousClock=state.previousTimeLayer()
        ?std::optional<I::HermiteGalsTimeLayer>(*state.previousTimeLayer()):std::nullopt;
    const auto unchanged=[&] {
        check(&state.field(0)==acceptedField&&&state.solver()==acceptedSolver
            &&identicalCharge(state.integratedCharge(),old)&&identicalCharge(acceptedField->transport().panelCharge(),old)
            &&sameChargeClock(state.timeLayer(),acceptedClock)&&sameChargeClock(state.geometrySourceTimeLayer(),representationClock)
            &&static_cast<bool>(state.previousTimeLayer())==previousClock.has_value()
            &&(!previousClock||(sameChargeClock(*state.previousTimeLayer(),*previousClock)
                &&state.previousIntegratedCharge()&&identicalCharge(*state.previousIntegratedCharge(),previous))),
            "E-C1 stage preparation/rejection mutated accepted owner/Q/history/full clocks");
    };
    const SurfaceChargeField* retainedField=nullptr;
    // Destroy the moving wrapper and all local publication/transfer handles;
    // the returned numerical step must retain the real target Field lifetime.
    const auto retained=[&] {
        const auto pub=chargeMovingPublication(mesh,source,target,displacement,totalDisplacement);
        const auto geometry=ElectricPhaseGeometry::build(mesh,pub);
        const LevelSetInterfacePanels panels(pub.topology().curves.front());
        const auto transfer=acceptedField->prepareMaterialTransfer(mesh,panels,[&](Vector2D p) {
            ++metric.materialQueries;return p+Vector2D{0,displacement};
        },static_cast<int>(target.step),acceptedClock.time(),target.time());
        const auto queries=metric.materialQueries;
        check(queries==acceptedField->geometry().panels().size()+1,
            "E-C1 actual material endpoints were not traced exactly once");
        const auto stage=state.prepareMovingStage(mesh,geometry,{transfer},target,electric,sigma);
        const auto block=chargeStageBlockReference(stage,electric,sigma);
        metric.chargeError=block.chargeError;metric.potentialError=block.potentialError;
        metric.currentError=block.currentError;metric.blockResidual=block.blockResidual;
        metric.panelEquation=block.panelEquation;metric.totalEquation=block.totalEquation;
        metric.dt=block.dt;metric.previousDt=block.previousDt;
        metric.a0=block.a0;metric.a1=block.a1;metric.a2=block.a2;
        unchanged();chargeMovingCoeval(mesh,stage,pub,target);
        check(sameChargeClock(stage.sourceTimeLayer(),acceptedClock)
            &&sameChargeClock(stage.sourceGeometryTimeLayer(),representationClock)
            &&static_cast<bool>(stage.previousSourceTimeLayer())==previousClock.has_value()
            &&(!previousClock||sameChargeClock(*stage.previousSourceTimeLayer(),*previousClock))
            &&identicalCharge(stage.field(0).transport().panelCharge(),stage.step().integratedCharge()),
            "E-C1 stage altered historical physical clocks or its coeval target Field Q");
        chargeIntegralConserved(old,stage.mappedIntegratedCharge());
        const auto mappedReference=chargeSphericalOverlapReference(transfer,old,displacement,metric.boundaryError);
        metric.overlapError=chargeMappedReferenceError(stage.mappedIntegratedCharge(),mappedReference,old);
        if(previousClock) {
            check(stage.mappedPreviousIntegratedCharge(),"E-C1 midpoint omitted accepted previous charge");
            chargeIntegralConserved(previous,*stage.mappedPreviousIntegratedCharge());
            const auto previousReference=chargeSphericalOverlapReference(transfer,previous,displacement,metric.boundaryError);
            metric.overlapError=std::max(metric.overlapError,
                chargeMappedReferenceError(*stage.mappedPreviousIntegratedCharge(),previousReference,previous));
        }
        const auto reject=[&](auto action,const char* reason) {
            rejected(action,reason);++metric.rejections;unchanged();
        };
        reject([&]{state.commit(stage);},"E-C1 publicly committed a Midpoint moving stage");
        reject([&]{state.commit(stage.step());},"E-C1 publicly committed a Midpoint numerical step");
        reject([&]{state.solver().advance(state,target);},"E-C1 public fixed advance accepted a Midpoint");
        reject([&]{state.prepareMoving(mesh,geometry,{transfer},target,electric,sigma);},
            "E-C1 End preparation accepted a Midpoint publication");
        auto wrongClock=target;wrongClock.begin+=.001;wrongClock.end-=.001;
        reject([&]{state.prepareMovingStage(mesh,geometry,{transfer},wrongClock,electric,sigma);},
            "E-C1 accepted the same midpoint physical time with a foreign full clock");
        auto endClock=target;endClock.position=I::HermiteGalsTimeLayer::Position::End;
        reject([&]{state.prepareMovingStage(mesh,geometry,{transfer},endClock,electric,sigma);},
            "E-C1 stage preparation accepted an End clock");
        const auto wrongInterval=acceptedField->prepareMaterialTransfer(mesh,panels,
            [=](Vector2D p){return p+Vector2D{0,displacement};},static_cast<int>(target.step),
            acceptedClock.time()+.001,target.time());
        reject([&]{state.prepareMovingStage(mesh,geometry,{wrongInterval},target,electric,sigma);},
            "E-C1 stage accepted a disconnected transfer physical interval");
        ElectricChargeTimeControls invalid;invalid.maximumIterations=0;
        reject([&]{state.prepareMovingStage(mesh,geometry,{transfer},target,electric,sigma,{}, {},invalid);},
            "E-C1 stage accepted invalid Schur solve controls");
        check(metric.materialQueries==queries,"E-C1 stage/remapped histories retraced prepared material endpoints");
        retainedField=&stage.field(0);
        return stage.step();
    }();
    check(retainedField&&identicalCharge(retainedField->transport().panelCharge(),retained.integratedCharge())
        &&retainedField->geometry().usesSource(retained.potential().geometry().publication().state().geometry().snapshot()),
        "E-C1 numerical step did not retain actual target Field after local stage destruction");
    const auto& op=retained.current().momentOperator();const auto load=op.apply(retained.integratedCharge());
    check(identicalCharge(load.panelCharge(),retained.integratedCharge()),"E-C1 retained stage P_Q became invalid");
    const auto maxwell=ElectricMaxwellTrace::prepare(retained.potential());
    maxwell.validateSource(mesh,retained.potential().geometry().publication(),target);
    unchanged();
    const auto precision=std::cout.precision();
    std::cout<<std::setprecision(17)<<"ELECTRIC_HHO_EC1_STAGE,"<<(previousClock?"history":"startup")
        <<",h="<<metric.dt<<",full_interval="<<target.end-target.begin<<",previous_h="<<metric.previousDt
        <<",a0="<<metric.a0<<",a1="<<metric.a1<<",a2="<<metric.a2
        <<",block_Q_error="<<metric.chargeError<<",block_V_error="<<metric.potentialError
        <<",block_current_error="<<metric.currentError<<",block_residual="<<metric.blockResidual
        <<",panel_equation="<<metric.panelEquation<<",total_equation="<<metric.totalEquation
        <<",overlap="<<metric.overlapError<<",boundary_area="<<metric.boundaryError
        <<",material_queries="<<metric.materialQueries<<",atomic_checks="<<metric.rejections
        <<",stage_step="<<target.step<<",stage_begin="<<target.begin<<",stage_end="<<target.end
        <<",stage_position="<<static_cast<int>(target.position)
        <<",Q_V_I_Maxwell_fullclock=matched,accepted_history=unchanged,retained_Field=valid,PASS\n";
    std::cout.precision(precision);
    return retained;
}

void chargeMovingContract() {
    static_assert(!std::is_default_constructible_v<SurfaceChargeMaterialTransfer>);
    DictionaryMesh input;input.Nx=input.Ny=8;input.maxInitialIterations=0;input.maxCellNum=512;
    QuadtreeMesh mesh(input);const auto initialPub=publication(mesh,true);
    const auto initialGeometry=ElectricPhaseGeometry::build(mesh,initialPub);
    SurfaceChargeField initialField(mesh,initialPub.topology().curves.front(),0.);
    auto q0=initialField.transport().panelCharge();bool positive=false,negative=false;
    for(std::size_t i=0;i<q0.size();++i) {
        q0[i]=initialField.transport().panelMeasure()[i]*(std::sin(1.7*(i+1))+.2+.01*i);
        positive=positive||q0[i]>0;negative=negative||q0[i]<0;
    }
    check(positive&&negative,"E-M1 real RZ panels do not contain signed nonuniform Q");
    initialField.assignIntegratedCharge(q0);
    const auto initialOp=ElectricSurfaceChargeMomentOperator::prepare(initialGeometry,{&initialField},initialGeometry.timeLayer());
    ElectricPhaseProblem electric;electric.epsilonByPhase={{epsilonNegative,epsilonPositive}};
    for(int side=0;side<4;++side) {
        electric.boundaries[side].type=side==0?ElectricCartesianBoundary::Type::Neumann:ElectricCartesianBoundary::Type::Dirichlet;
        electric.boundaries[side].value=side==0?std::function<double(Vector2D)>([](Vector2D){return 0.;})
            :std::function<double(Vector2D)>([](Vector2D p){return .05*p.y;});
    }
    const auto initialPotential=solveElectricPhaseElectrostatics(initialGeometry,electric,initialOp.apply(q0));
    const auto initialSolver=ElectricFixedGeometryChargeSolver::prepare(initialPotential,initialOp,{{.7,.03}});
    auto state=ElectricSurfaceChargeTimeState::initialize(initialSolver,{&initialField},initialGeometry.timeLayer());
    int rejections=0;const auto reject=[&](auto action,const char* reason) {rejected(action,reason);++rejections;};
    constexpr double dt1=.06,dt2=.09,displacement=.008,rate=.14;
    const auto clock0=state.timeLayer(),clock1=nextChargeClock(state,dt1);
    const auto pub1=chargeMovingPublication(mesh,initialPub,clock1,0,0);
    const auto geometry1=ElectricPhaseGeometry::build(mesh,pub1);
    const LevelSetInterfacePanels panels1(pub1.topology().curves.front());
    const auto identity=initialField.prepareMaterialTransfer(mesh,panels1,[](Vector2D p){return p;},
        static_cast<int>(clock1.step),clock0.time(),clock1.time());
    check(identity.geometryReport().maxMaterialDisplacement==0&&identity.sourcePanels().axisymmetric()
        &&identity.sourcePanels().beginEndpointKind()==LevelSetInterfacePanels::EndpointKind::Axis
        &&identity.sourcePanels().endEndpointKind()==LevelSetInterfacePanels::EndpointKind::Axis,
        "E-M1 startup is not an actual zero-velocity RZ axis-to-axis material transfer");
    ChargeStageReference startupStageMetric;
    const auto startupStage=chargeACMidpointContract(mesh,state,initialPub,clock1,0,0,
        electric,{{.7,.03}},startupStageMetric);
    ElectricContract::close(startupStageMetric.dt,dt1*.5,1e-14,"E-C1 startup midpoint did not use half of the real AC interval");
    const auto fixedReference=initialSolver.advance(state,clock1);
    const auto first=state.prepareMoving(mesh,geometry1,{identity},clock1,electric,{{.7,.03}});
    chargeTimeEquation(first.step(),first.mappedIntegratedCharge(),nullptr);chargeMovingCoeval(mesh,first,pub1,clock1);
    check(first.step().report().temporalOrder==1&&!first.mappedPreviousIntegratedCharge()
        &&!first.previousSourceTimeLayer()&&sameChargeClock(first.sourceTimeLayer(),clock0),
        "E-M1 actual missing history did not preserve BE startup/source clock");
    check(identicalCharge(first.field(0).transport().panelCharge(),first.step().integratedCharge()),
        "E-M1 staged target Field does not contain final candidate Q");
    double boundaryError=0;
    const auto identityReference=chargeSphericalOverlapReference(identity,q0,0,boundaryError);
    const double identityError=chargeMappedReferenceError(first.mappedIntegratedCharge(),identityReference,q0);
    chargeIntegralConserved(q0,first.mappedIntegratedCharge());
    for(std::size_t i=0;i<q0.size();++i)ElectricContract::close(first.step().integratedCharge()[i],
        fixedReference.integratedCharge()[i],1e-9,"E-M1 zero-motion publication changed the existing fixed-geometry BE result");
    check(identicalCharge(initialField.transport().panelCharge(),q0)&&identicalCharge(state.integratedCharge(),q0)
        &&sameChargeClock(state.timeLayer(),clock0)&&!state.previousTimeLayer(),
        "E-M1 identity preparation changed accepted Field/history before commit");
    state.commit(first);
    const auto q1=state.integratedCharge(),oldOn1=*state.previousIntegratedCharge();
    check(!identicalCharge(q1,q0)&&sameChargeClock(state.timeLayer(),clock1)
        &&sameChargeClock(state.geometrySourceTimeLayer(),clock1)&&state.previousTimeLayer()
        &&sameChargeClock(*state.previousTimeLayer(),clock0)
        &&identicalCharge(oldOn1,first.mappedIntegratedCharge())
        &&identicalCharge(state.field(0).transport().panelCharge(),q1)
        &&identicalCharge(initialField.transport().panelCharge(),q0),
        "E-M1 identity acceptance lost the true history clock/mapped Q or changed the old Field");
    positive=negative=false;for(double q:q1){positive=positive||q>0;negative=negative||q<0;}
    check(positive&&negative,"E-M1 real BE evolution lost the signed history fixture");
    reject([&]{state.commit(first);},"E-M1 accepted a previously committed moving candidate");
    reject([&]{state.commit(fixedReference);},"E-M1 accepted the old solver's stale fixed-geometry candidate");
    // The old borrowed Field is writable through a new legitimate owner again;
    // retained first/potential/current handles must not retain its live lease.
    auto foreignState=ElectricSurfaceChargeTimeState::initialize(initialSolver,{&initialField},clock0);
    const auto foreignQ=foreignState.integratedCharge();
    const auto clock2=nextChargeClock(state,dt2);
    const auto pub2=chargeMovingPublication(mesh,pub1,clock2,displacement,displacement);
    const auto geometry2=ElectricPhaseGeometry::build(mesh,pub2);
    const LevelSetInterfacePanels panels2(pub2.topology().curves.front());
    const auto* sourceField=&state.field(0);const auto* sourceSolver=&state.solver();
    std::size_t materialQueries=0;
    const auto transfer=sourceField->prepareMaterialTransfer(mesh,panels2,[&](Vector2D p){
        ++materialQueries;return p+Vector2D{0,displacement};
    },static_cast<int>(clock2.step),clock1.time(),clock2.time());
    const auto queriesAfterFactory=materialQueries;
    check(materialQueries==sourceField->geometry().panels().size()+1
        &&transfer.sourcePanels().usesSource(pub1.state().geometry().snapshot())
        &&transfer.targetPanels().usesSource(pub2.state().geometry().snapshot())
        &&transfer.geometryReport().maxMaterialDisplacement>0,
        "E-M1 nonzero move did not trace the actual native material endpoints exactly once");
    ElectricContract::close(transfer.geometryReport().maxMaterialDisplacement,displacement,1e-14,
        "E-M1 factory material displacement differs from the prescribed physical translation");
    const auto unchanged=[&] {
        check(&state.field(0)==sourceField&&&state.solver()==sourceSolver
            &&identicalCharge(state.integratedCharge(),q1)
            &&state.previousIntegratedCharge()&&identicalCharge(*state.previousIntegratedCharge(),oldOn1)
            &&identicalCharge(sourceField->transport().panelCharge(),q1)
            &&sameChargeClock(state.timeLayer(),clock1)&&sameChargeClock(state.geometrySourceTimeLayer(),clock1)
            &&state.previousTimeLayer()&&sameChargeClock(*state.previousTimeLayer(),clock0),
            "E-M1 rejected/prepared candidate changed accepted Field/binding/history/full clocks");
    };
    ChargeStageReference historyStageMetric;
    const auto historyStage=chargeACMidpointContract(mesh,state,pub1,clock2,displacement*.5,displacement*.5,
        electric,{{.7,.03}},historyStageMetric);
    ElectricContract::close(historyStageMetric.dt,dt2*.5,1e-14,"E-C1 history midpoint did not use half of the real AC interval");
    unchanged();
    const auto nonproportionalEnd=state.prepareMoving(mesh,geometry2,{transfer},clock2,electric,{{.7,.03}});
    const auto nonproportionalEndMetric=chargeStageBlockReference(nonproportionalEnd,electric,{{.7,.03}});
    chargeMovingCoeval(mesh,nonproportionalEnd,pub2,clock2);unchanged();
    check(materialQueries==queriesAfterFactory,"E-C1 readonly End retraced its shared prepared material map");
    const std::array<double,2> sigma{{epsilonNegative*rate,epsilonPositive*rate}};
    auto wrongClock=clock2;++wrongClock.step;
    reject([&]{state.prepareMoving(mesh,geometry2,{transfer},wrongClock,electric,sigma);},
        "E-M1 accepted equal physical time with a foreign step/full target clock");unchanged();
    wrongClock=clock2;wrongClock.position=I::HermiteGalsTimeLayer::Position::Midpoint;
    wrongClock.end=2*clock2.end-clock2.begin; // Same physical time, different full clock.
    reject([&]{state.prepareMoving(mesh,geometry2,{transfer},wrongClock,electric,sigma);},
        "E-M1 accepted equal physical time at a midpoint instead of the actual target End");unchanged();
    const auto wrongInterval=sourceField->prepareMaterialTransfer(mesh,panels2,
        [=](Vector2D p){return p+Vector2D{0,displacement};},static_cast<int>(clock2.step),
        clock1.time()+.001,clock2.time());
    reject([&]{state.prepareMoving(mesh,geometry2,{wrongInterval},clock2,electric,sigma);},
        "E-M1 accepted a disconnected material source physical interval");unchanged();
    const auto foreignSourcePub=chargeMovingPublication(mesh,initialPub,clock1,0,0);
    SurfaceChargeField foreignSource(mesh,foreignSourcePub.topology().curves.front(),0.);
    foreignSource.assignIntegratedCharge(q1);
    const auto foreignSourceTransfer=foreignSource.prepareMaterialTransfer(mesh,panels2,
        [=](Vector2D p){return p+Vector2D{0,displacement};},static_cast<int>(clock2.step),clock1.time(),clock2.time());
    reject([&]{state.prepareMoving(mesh,geometry2,{foreignSourceTransfer},clock2,electric,sigma);},
        "E-M1 accepted equal scalar/panels from a foreign source publication owner");unchanged();
    const auto foreignTargetPub=chargeMovingPublication(mesh,pub1,clock2,displacement,displacement);
    const LevelSetInterfacePanels foreignTargetPanels(foreignTargetPub.topology().curves.front());
    const auto foreignTargetTransfer=sourceField->prepareMaterialTransfer(mesh,foreignTargetPanels,
        [=](Vector2D p){return p+Vector2D{0,displacement};},static_cast<int>(clock2.step),clock1.time(),clock2.time());
    reject([&]{state.prepareMoving(mesh,geometry2,{foreignTargetTransfer},clock2,electric,sigma);},
        "E-M1 accepted equal target panels owned by another target publication");unchanged();
    ElectricChargeTimeControls failedSolve;failedSolve.maximumIterations=0;
    reject([&]{state.prepareMoving(mesh,geometry2,{transfer},clock2,electric,sigma,{}, {},failedSolve);},
        "E-M1 moving candidate accepted invalid Schur controls");unchanged();
    reject([&]{ElectricSurfaceChargeTimeState::initialize(state.solver(),
        {const_cast<SurfaceChargeField*>(sourceField)},clock1);},
        "E-M1 rejected candidates released the live accepted Field lease");unchanged();
    const auto staleFixed=state.solver().advance(state,clock2);
    const auto second=state.prepareMoving(mesh,geometry2,{transfer},clock2,electric,sigma);
    check(materialQueries==queriesAfterFactory,"E-M1 charge/history remap retraced a completed material trajectory");
    unchanged();chargeMovingCoeval(mesh,second,pub2,clock2);
    check(second.mappedPreviousIntegratedCharge()&&second.previousSourceTimeLayer()
        &&sameChargeClock(second.sourceTimeLayer(),clock1)
        &&sameChargeClock(*second.previousSourceTimeLayer(),clock0)
        &&sameChargeClock(second.sourceGeometryTimeLayer(),clock1),
        "E-M1 remap retagged the two historical physical clocks or their old representation clock");
    const auto& mappedOld=second.mappedIntegratedCharge();const auto& mappedOlder=*second.mappedPreviousIntegratedCharge();
    chargeIntegralConserved(q1,mappedOld);chargeIntegralConserved(oldOn1,mappedOlder);
    const auto referenceOld=chargeSphericalOverlapReference(transfer,q1,displacement,boundaryError);
    const auto referenceOlder=chargeSphericalOverlapReference(transfer,oldOn1,displacement,boundaryError);
    const double oldError=chargeMappedReferenceError(mappedOld,referenceOld,q1);
    const double olderError=chargeMappedReferenceError(mappedOlder,referenceOlder,oldOn1);
    check(second.remapReports().size()==1&&second.remapReports()[0].transport.success
        &&second.historyRemapReports().size()==1&&second.historyRemapReports()[0].size()==1
        &&second.historyRemapReports()[0][0].success,
        "E-M1 current/previous histories did not use one successful real component transfer");
    chargeTimeEquation(second.step(),mappedOld,&mappedOlder);
    const auto& report=second.step().report();const double ratio=dt2/dt1;
    const double a0=(1+2*ratio)/(1+ratio),a1=-(1+ratio),a2=ratio*ratio/(1+ratio);
    check(report.temporalOrder==2&&report.a2!=0,"E-M1 moving rebind silently restarted BE and discarded history");
    ElectricContract::close(report.previousDt,dt1,1e-14,"E-M1 BDF2 previous increment lost its original physical clock");
    ElectricContract::close(report.dt,dt2,1e-14,"E-M1 moving BDF2 used another target increment");
    ElectricContract::close(report.a0,a0,1e-14,"E-M1 moving BDF2 a0 ignores unequal actual increments");
    ElectricContract::close(report.a1,a1,1e-14,"E-M1 moving BDF2 a1 ignores unequal actual increments");
    ElectricContract::close(report.a2,a2,1e-14,"E-M1 moving BDF2 a2 ignores unequal actual increments");
    const auto& q2=second.step().integratedCharge();const auto& current=second.step().current().panelCurrent();
    long double referenceError=0,currentError=0,currentScale=0;double independentPanelEquation=0;
    for(std::size_t i=0;i<q2.size();++i) {
        const double reference=(-a1*mappedOld[i]-a2*mappedOlder[i])/(a0+dt2*rate);
        referenceError+=std::abs(q2[i]-reference);currentError+=std::abs(current[i]+rate*q2[i]);
        currentScale+=std::abs(current[i])+rate*std::abs(q2[i]);
        const double residual=report.a0*q2[i]+report.a1*mappedOld[i]+report.a2*mappedOlder[i]-report.dt*current[i];
        const double scale=std::abs(report.a0*q2[i])+std::abs(report.a1*mappedOld[i])
            +std::abs(report.a2*mappedOlder[i])+std::abs(report.dt*current[i]);
        independentPanelEquation=std::max(independentPanelEquation,scale?std::abs(residual)/scale:0);
    }
    const double relativeTimeReference=static_cast<double>(referenceError)/chargeIntegralScale(q2);
    const double relativeScalarCurrent=static_cast<double>(currentError/currentScale);
    check(relativeTimeReference<=1e-10&&relativeScalarCurrent<=1e-10,
        "E-M1 moved BDF2 differs from the independent two-history dielectric relaxation/current reference");
    check(second.step().current().conductivityByPhase()==sigma
        &&second.step().current().epsilonByPhase()==electric.epsilonByPhase
        &&second.step().current().report().relativeGlobalCurrentBalance<=1e-10,
        "E-M1 target preparation reused old conductivity/material or lost the raw current balance");
    check(identicalCharge(second.field(0).transport().panelCharge(),q2),
        "E-M1 target candidate Field and its final coeval Q/V/I disagree");
    reject([&]{foreignState.commit(second);},"E-M1 foreign accepted owner committed another state's moving candidate");
    check(identicalCharge(foreignState.integratedCharge(),foreignQ)
        &&identicalCharge(initialField.transport().panelCharge(),q0)&&sameChargeClock(foreignState.timeLayer(),clock0)
        &&!foreignState.previousTimeLayer(),"E-M1 foreign commit rejection mutated its Field/history");unchanged();
    const auto duplicate=second;
    state.commit(second);
    check(&state.field(0)!=sourceField&&state.solver().momentOperator().geometry().geometryIdentity()==pub2.state().identity()
        &&sameChargeClock(state.timeLayer(),clock2)&&sameChargeClock(state.geometrySourceTimeLayer(),clock2)
        &&sameChargeClock(state.solver().geometrySourceTimeLayer(),clock2)
        &&state.previousTimeLayer()&&sameChargeClock(*state.previousTimeLayer(),clock1)
        &&state.previousIntegratedCharge()&&identicalCharge(*state.previousIntegratedCharge(),mappedOld)
        &&identicalCharge(state.integratedCharge(),q2)&&identicalCharge(state.field(0).transport().panelCharge(),q2)
        &&identicalCharge(sourceField->transport().panelCharge(),q1),
        "E-M1 commit did not atomically replace target Field/solver and retain mapped true physical history");
    reject([&]{state.commit(duplicate);},"E-M1 successful commit left its stale moving candidate reusable");
    reject([&]{state.commit(staleFixed);},"E-M1 successful move accepted an old source-geometry fixed candidate");
    check(identicalCharge(state.integratedCharge(),q2)&&sameChargeClock(state.timeLayer(),clock2)
        &&identicalCharge(*state.previousIntegratedCharge(),mappedOld),
        "E-M1 stale candidate rejection changed newly accepted state");
    rejected([&]{state.commit(startupStage);},"E-C1 committed a retained startup Midpoint after End acceptance");
    rejected([&]{state.commit(historyStage);},"E-C1 committed a retained history Midpoint after End acceptance");
    rejected([&]{state.commit(nonproportionalEnd);},"E-C1 committed a stale readonly nonproportional End candidate");
    check(sameChargeClock(*state.previousTimeLayer(),clock1)
        &&identicalCharge(*state.previousIntegratedCharge(),mappedOld)&&identicalCharge(state.integratedCharge(),q2),
        "E-C1 actual End acceptance stored stage Q/clock in the physical history");
    const auto retainedMaxwell=ElectricMaxwellTrace::prepare(historyStage.potential());
    retainedMaxwell.validateSource(mesh,historyStage.potential().geometry().publication(),historyStage.timeLayer());
    check(sameChargeClock(historyStage.timeLayer(),historyStage.current().timeLayer())
        &&historyStage.current().momentOperator().apply(historyStage.integratedCharge()).panelCharge()==historyStage.integratedCharge(),
        "E-C1 retained midpoint numerical handles lost their actual target after source End acceptance");
    const auto stagePrecision=std::cout.precision();
    std::cout<<std::setprecision(17)<<"ELECTRIC_HHO_EC1_AC_STAGE,RZ,cells="<<mesh.cellIndexNum()
        <<",startup_h="<<startupStageMetric.dt<<",history_h="<<historyStageMetric.dt
        <<",history_previous_dt="<<historyStageMetric.previousDt
        <<",history_a0="<<historyStageMetric.a0<<",history_a1="<<historyStageMetric.a1<<",history_a2="<<historyStageMetric.a2
        <<",startup_charge_block="<<startupStageMetric.chargeError<<",startup_potential_block="<<startupStageMetric.potentialError
        <<",startup_current_block="<<startupStageMetric.currentError<<",history_charge_block="<<historyStageMetric.chargeError
        <<",history_potential_block="<<historyStageMetric.potentialError<<",history_current_block="<<historyStageMetric.currentError
        <<",End_charge_block="<<nonproportionalEndMetric.chargeError<<",End_potential_block="<<nonproportionalEndMetric.potentialError
        <<",End_current_block="<<nonproportionalEndMetric.currentError
        <<",block_residual="<<std::max({startupStageMetric.blockResidual,historyStageMetric.blockResidual,nonproportionalEndMetric.blockResidual})
        <<",panel_equation="<<std::max({startupStageMetric.panelEquation,historyStageMetric.panelEquation,nonproportionalEndMetric.panelEquation})
        <<",total_equation="<<std::max({startupStageMetric.totalEquation,historyStageMetric.totalEquation,nonproportionalEndMetric.totalEquation})
        <<",material_overlap="<<std::max(startupStageMetric.overlapError,historyStageMetric.overlapError)
        <<",boundary_area="<<std::max(startupStageMetric.boundaryError,historyStageMetric.boundaryError)
        <<",midpoint_material_queries="<<startupStageMetric.materialQueries+historyStageMetric.materialQueries
        <<",atomic_stage_rejections="<<startupStageMetric.rejections+historyStageMetric.rejections
        <<",Q_V_I_Maxwell_fullclock=matched,physical_history=End_only,retained_Field=valid"
        <<",nonproportional_End=readonly,actual_End_sigma=0.7:0.14,PASS\n";
    std::cout.precision(stagePrecision);
    const auto relativeTotal=[&](const std::vector<double>& a,const std::vector<double>& b) {
        long double x=0,y=0;for(double q:a)x+=q;for(double q:b)y+=q;
        return static_cast<double>(std::abs(x-y))/chargeIntegralScale(a);
    };
    const auto oldPrecision=std::cout.precision();const auto& gauss=second.step().potential().report();
    std::cout<<std::setprecision(17)<<"ELECTRIC_HHO_EM1_MOVING,RZ,cells="<<mesh.cellIndexNum()
        <<",source_panels="<<q1.size()<<",target_panels="<<q2.size()<<",material_displacement="<<displacement
        <<",dt="<<report.dt<<",previous_dt="<<report.previousDt<<",a0="<<report.a0<<",a1="<<report.a1<<",a2="<<report.a2
        <<",temporal_order="<<report.temporalOrder<<",identity_overlap="<<identityError
        <<",current_history_overlap="<<oldError<<",previous_history_overlap="<<olderError
        <<",boundary_area="<<boundaryError<<",current_history_total="<<relativeTotal(q1,mappedOld)
        <<",previous_history_total="<<relativeTotal(oldOn1,mappedOlder)<<",time_reference="<<relativeTimeReference
        <<",scalar_current="<<relativeScalarCurrent<<",panel_equation="<<independentPanelEquation
        <<",Schur="<<report.relativeSchurResidual<<",global_charge_equation="<<report.relativeTotalEquationResidual
        <<",local_Gauss="<<gauss.maxRelativeLocalGaussResidual<<",global_Gauss="<<gauss.relativeGlobalGaussResidual
        <<",current_balance="<<second.step().current().report().relativeGlobalCurrentBalance
        <<",material_queries="<<materialQueries<<",atomic_owner_rejections="<<rejections
        <<",source_step="<<clock1.step<<",source_begin="<<clock1.begin<<",source_end="<<clock1.end
        <<",previous_step="<<clock0.step<<",previous_position="<<static_cast<int>(clock0.position)
        <<",target_step="<<clock2.step<<",target_begin="<<clock2.begin<<",target_end="<<clock2.end
        <<",target_position="<<static_cast<int>(clock2.position)<<",Q_V_I_Maxwell_fullclock=matched,PASS\n";
    std::cout.precision(oldPrecision);
}

void run(bool withCharge=false,bool withCurrent=false,bool withTime=false,bool withMaxwell=false) {
    if(withTime)chargeMovingContract();
    check(!withMaxwell||withCharge,"P5 Maxwell contract requires its original charged solve");
    check(!withTime||withCurrent,"P4 time contract requires the original current operator");
    check(!withCurrent||withCharge,"P3 current contract requires its retained P2 charged result");
    std::optional<ElectricSurfaceChargeMomentOperator> previousMoment;
    std::array<Error,2> dipole;
    std::array<MaxwellError,2> maxwellError;
    MaxwellError sphericalMaxwell;
    for(bool axis:{false,true}) {
        DictionaryMesh input;input.Nx=input.Ny=4;input.maxInitialIterations=0;input.maxCellNum=512;
        QuadtreeMesh mesh(input);const auto pub=publication(mesh,axis);
        const auto geometry=ElectricPhaseGeometry::build(mesh,pub);geometryChecks(geometry);
        const auto identity=pub.state().identity();const auto* scalarOwner=&pub.state().geometry().snapshot();
        std::optional<ElectricPhaseResult> baseline;
        std::optional<ElectricPhaseResult> constant;
        solveAndCheck(mesh,geometry,0,withMaxwell?&constant:nullptr);
        if(withMaxwell) {
            const auto original=constant->skeletonSolution();
            const auto zero=ElectricMaxwellTrace::prepare(*constant);double maximumEt=0,maximumTraction=0;
            check(!constant->surfaceChargeLoad()&&sameChargeClock(zero.timeLayer(),constant->timeLayer())
                &&sameChargeClock(zero.geometrySourceTimeLayer(),geometry.timeLayer()),
                "P5 uncharged constant field lost its legal source/physical clock provenance");
            zero.validateSource(mesh,pub,constant->timeLayer());
            for(std::size_t face:zero.interfaceFaces())for(const auto& sample:zero.quadrature(face)) {
                maximumEt=std::max(maximumEt,std::abs(sample.tangentialField));
                maximumTraction=std::max(maximumTraction,std::hypot(sample.tractionJump.x,sample.tractionJump.y));
            }
            check(maximumEt<=2e-7&&maximumTraction<=2e-7&&constant->skeletonSolution()==original,
                "P5 Maxwell failed the existing legal zero-field absolute consistency check");
            std::cout<<"ELECTRIC_HHO_P5_CONSTANT,"<<axis<<",Et_max="<<maximumEt
                <<",traction_max="<<maximumTraction<<",PASS\n";
        }
        solveAndCheck(mesh,geometry,1,withCharge?&baseline:nullptr);
        if(!axis)dipole[0]=solveAndCheck(mesh,geometry,2);
        if(withCharge)chargeContract(mesh,geometry,&*baseline,false,withCurrent,withCurrent?&previousMoment:nullptr,
            withTime,withMaxwell,axis?&sphericalMaxwell:&maxwellError[0]);
        const auto actualClock=geometry.timeLayer(),publicationClock=pub.state().timeLayer();
        check(actualClock.step==publicationClock.step&&actualClock.begin==publicationClock.begin
            &&actualClock.end==publicationClock.end&&actualClock.position==publicationClock.position,
            "P1 prepared geometry lost a field of its complete publication clock");
        if(!axis) {
            const auto other=publication(mesh,axis);
            rejected([&]{geometry.validateSource(mesh,other);},"equal sampled geometry substituted for P1 source identity");
            if(withCharge) {
                SurfaceChargeField foreign(mesh,other.topology().curves[0],radius);
                rejected([&]{ElectricSurfaceChargeMomentOperator::prepare(geometry,{&foreign},geometry.timeLayer());},
                    "P_Q accepted a numerically equal curve from a foreign immutable source");
            }
            const I::HermiteGalsTimeLayer later{1,0,.01,I::HermiteGalsTimeLayer::Position::End};
            const auto anotherClock=publication(mesh,axis,later);
            rejected([&]{geometry.validateSource(mesh,anotherClock);},"P1 geometry accepted a foreign full source clock");
        }
        auto invalid=problem(axis,0);invalid.epsilonByPhase[0]=0;
        rejected([&]{solveElectricPhaseElectrostatics(geometry,invalid);},"P1 accepted nonpositive phase permittivity");
        check(pub.state().identity()==identity&&&pub.state().geometry().snapshot()==scalarOwner,
            "P1 solve or failed preflight changed its original publication owner");
        auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]);vertex.setCoord(vertex.coord()+Vector2D{.01,0});
        rejected([&]{geometry.validateMesh(mesh);},"P1 accepted stale native mesh identity");
    }
    DictionaryMesh fineInput;fineInput.Nx=fineInput.Ny=8;fineInput.maxInitialIterations=0;fineInput.maxCellNum=512;
    QuadtreeMesh fine(fineInput);const auto finePub=publication(fine,false);
    const auto fineGeometry=ElectricPhaseGeometry::build(fine,finePub);geometryChecks(fineGeometry);
    dipole[1]=solveAndCheck(fine,fineGeometry,2);
    if(withCharge)chargeContract(fine,fineGeometry,nullptr,true,withCurrent,withCurrent?&previousMoment:nullptr,
        false,withMaxwell,&maxwellError[1]);
    if(withMaxwell) {
        // Collect the two independent physical cases before the same hard gates.
        // Failed cases were printed FAIL_PHYSICAL_BUDGET, never PASS.
        if(!sphericalMaxwell.physicalBudgetPassed||!maxwellError[1].physicalBudgetPassed) {
            std::ostringstream reason;reason.precision(17);
            reason<<"P5 frozen physical Maxwell budget failed";
            for(const auto& named:std::array<std::pair<const char*,const MaxwellError*>,2>{{
                {"RZ_sphere",&sphericalMaxwell},{"XY_fine",&maxwellError[1]}}}) {
                reason<<" "<<named.first<<" [";
                for(int quantity=0;quantity<6;++quantity)reason<<" relative_"<<quantity<<"="<<named.second->relative[quantity];
                reason<<" force="<<named.second->force<<" power="<<named.second->power
                    <<" passed="<<named.second->physicalBudgetPassed<<"]";
            }
            throw std::runtime_error(reason.str());
        }
        for(int quantity=0;quantity<6;++quantity)check(maxwellError[0].relative[quantity]>0
            &&maxwellError[1].relative[quantity]<maxwellError[0].relative[quantity],
            "P5 independent XY Maxwell quantity failed its frozen coarse/fine refinement gate");
        check(maxwellError[1].force<maxwellError[0].force,
            "P5 independent XY full-tensor force did not improve on refinement");
        std::cout<<"ELECTRIC_HHO_P5_REFINEMENT,4_to_8";
        for(int quantity=0;quantity<6;++quantity)std::cout<<",quantity_"<<quantity<<"_ratio="
            <<maxwellError[1].relative[quantity]/maxwellError[0].relative[quantity];
        std::cout<<",force_ratio="<<maxwellError[1].force/maxwellError[0].force<<",PASS\n";
    }
    check(dipole[0].value>0&&dipole[0].gradient>0,"analytic dipole refinement has no nonzero coarse error");
    check(dipole[1].value<.7*dipole[0].value&&dipole[1].gradient<.8*dipole[0].gradient
        &&dipole[1].value<=.05&&dipole[1].gradient<=.15,"P1 genuine nonpolynomial circle dipole did not pass frozen refinement/accuracy gates");
    std::cout<<"ELECTRIC_HHO_P1_DIPOLE,4_to_8,V_ratio="<<dipole[1].value/dipole[0].value
        <<",E_ratio="<<dipole[1].gradient/dipole[0].gradient<<",PASS\n";
}
} // namespace ElectricPhaseContract



namespace ElectricVofContract {
namespace I=LAMP::Interface;
using Clock=LAMP::TimeLayer;
using Publication=I::GeometricVofPublication;
using Boundary=ElectricVofBoundary;
using Type=Boundary::Type;
using Stock=std::array<long double,2>;
using Backend=ElectricVofControls::Backend;
constexpr long double gate=1e-10L;
bool sameClock(Clock a,Clock b) {
    return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;
}
void bounded(long double error,long double scale,long double tolerance,const char* message) {
    check(std::isfinite(error)&&std::isfinite(scale)&&scale>=0,message);
    check(scale>0?std::abs(error)<=tolerance*scale:error==0,message);
}
void backendReceipt(const ElectricVofField& field,Backend backend) {
    check(backend==Backend::Gmres||backend==Backend::SuperLU,"E_BACKEND invalid expected backend");
    check(field.report().usedGMRES==(backend==Backend::Gmres)&&
        field.report().usedSuperLU==(backend==Backend::SuperLU),"E_BACKEND actual solver receipt is not exclusive");
}
template<class F>void rejected(F f,const char* message) {
    bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,message);
}
std::unique_ptr<QuadtreeMesh> meshFor(int n=8) {
    DictionaryMesh cfg;cfg.Nx=cfg.Ny=n;cfg.minCoord={0,0};cfg.maxCoord={1,1};
    cfg.maxInitialIterations=0;cfg.maxCellNum=4096;
    auto mesh=std::make_unique<QuadtreeMesh>(cfg);
    DictionaryMesh::RefineBox box;box.min={.5,0};box.max={1,1};box.level=1;mesh->applyRefineBoxes({box});
    return mesh;
}
NativeLeafDomain domainFor(const QuadtreeMesh& mesh,bool rz) {
    FaceSegmentGraph graph;graph.build(mesh,[](const QuadtreeMesh&,const Cell&){return true;});
    return NativeLeafDomain::build(mesh,graph,rz);
}
std::vector<Stock> phaseStocks(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) {
    std::vector<Stock> q(mesh.cellIndexNum());
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto capacity=domain.cellVolumes()[domain.cellRowByIndex()[ci]];
        const auto center=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord();
        q[ci]=center.x<.5?Stock{capacity,0}:Stock{0,capacity};}
    return q;
}
std::vector<std::array<long double,8>> snapshot(const Publication& p) {
    std::vector<std::array<long double,8>> q;for(const auto& c:p.plicCells())q.push_back({c.inventory[0],c.inventory[1],
        c.measured[0],c.measured[1],c.normal.x,c.normal.y,c.offset,c.resolved?1.L:0.L});return q;
}
double affine(Vector2D p,bool rz) {return .17+(rz?p.y:p.x+.37*p.y);}
struct FaceReference {
    NativeLeafFace face;
    double area=0,epsilon=0,sigma=0,a=0,affineH=0;
};
std::vector<FaceReference> references(const QuadtreeMesh& mesh,const Publication& publication,
    const ElectricVofMaterials& m,double dt,std::vector<Boundary>& boundaries) {
    const auto& domain=publication.domain();std::vector<double> epsilon(mesh.cellIndexNum()),sigma(mesh.cellIndexNum());
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto& p=publication.plicCells()[ci];
        const long double capacity=domain.cellVolumes()[domain.cellRowByIndex()[ci]];
        epsilon[ci]=static_cast<double>((m.epsilonLiquid*p.inventory[0]+m.epsilonGas*p.inventory[1])/capacity);
        sigma[ci]=static_cast<double>((m.conductivityLiquid*p.inventory[0]+m.conductivityGas*p.inventory[1])/capacity);}
    std::vector<FaceReference> result;boundaries.resize(domain.nativeFaceCount());
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto lo=mesh.verts(edge.edgeToVert[0]).coord(),hi=mesh.verts(edge.edgeToVert[1]).coord();
        // Physical measure from the actual edge endpoints, not from the new operator.
        const double length=f.horizontal?std::abs(hi.x-lo.x):std::abs(hi.y-lo.y);
        const double radius=f.horizontal?(lo.x+hi.x)/2:lo.x;
        const double area=length*(domain.axisymmetric()?CartesianMeasure::twoPi*radius:1.);
        const auto distance=[&](int ci){const auto c=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord();
            return f.horizontal?std::abs(c.y-f.center.y):std::abs(c.x-f.center.x);};
        const int a=f.negativeCell,b=f.positiveCell;double ep=0,sg=0;
        if(a>=0&&b>=0){const double da=distance(a),db=distance(b);check(da>0&&db>0,"E1 actual half distances");
            ep=(da+db)/(da/epsilon[a]+db/epsilon[b]);
            sg=sigma[a]==0||sigma[b]==0?0:(da+db)/(da/sigma[a]+db/sigma[b]);
            boundaries[row]={Type::Interior,0};
        }else{const int ci=a>=0?a:b;ep=epsilon[ci];sg=sigma[ci];
            boundaries[row]=area==0?Boundary{Type::HomogeneousNeumann,0}:Boundary{Type::Dirichlet,affine(f.center,domain.axisymmetric())};}
        const double gradient=f.horizontal?(domain.axisymmetric()?1.:.37):(domain.axisymmetric()?0.:1.);
        result.push_back({f,area,ep,sg,ep+dt*sg,-area*gradient});
    }
    return result;
}
void scatter(std::vector<long double>& destination,const NativeLeafFace& f,long double flux) {
    if(f.negativeCell>=0)destination[f.negativeCell]+=flux;
    if(f.positiveCell>=0)destination[f.positiveCell]-=flux;
}
struct Metrics {long double gauss=0,charge=0,global=0,affineField=0,potential=0;};
Metrics verify(const QuadtreeMesh& mesh,const Publication& geometry,const ElectricVofCharge& source,
    const ElectricVofField& field,Clock candidate,const ElectricVofMaterials& material,
    const std::vector<FaceReference>& reference,bool analytic,bool zeroSigma,Backend backend=Backend::Gmres) {
    const auto& domain=geometry.domain();const auto& old=source.integratedCharge();const auto& next=field.chargeCandidate().integratedCharge();
    const auto& potential=field.potential();const double dt=candidate.time()-source.timeLayer().time();
    check(sameClock(field.geometryTime(),geometry.timeLayer())&&sameClock(field.chargeSourceTime(),source.timeLayer())&&
        sameClock(field.candidateTime(),candidate)&&sameClock(field.chargeCandidate().timeLayer(),candidate),"E1 separate full clocks lost");
    check(field.effectiveDt()==dt&&field.geometry().identity()==geometry.identity()&&field.sourceCharge().sameInstance(source),"E1 real dt/source owner lost");
    field.validateFor(mesh,geometry,source,candidate);check(field.chargeCandidate().identity()!=source.identity(),"E1 candidate charge reused source identity");
    check(potential.size()==old.size()&&next.size()==old.size()&&field.faces().size()==reference.size(),"E1 actual native field sizes");
    std::vector<long double> D(old.size()),J(old.size()),absD(old.size()),absJ(old.size());
    long double boundaryJ=0,boundaryAbsJ=0,sourceTotal=0,targetTotal=0,sourceAbs=0,targetAbs=0,potentialError=0;
    Metrics errors;
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto& c=geometry.plicCells()[ci];const auto capacity=domain.cellVolumes()[domain.cellRowByIndex()[ci]];
        const double ep=static_cast<double>((material.epsilonLiquid*c.inventory[0]+material.epsilonGas*c.inventory[1])/capacity);
        const double sg=static_cast<double>((material.conductivityLiquid*c.inventory[0]+material.conductivityGas*c.inventory[1])/capacity);
        bounded(field.epsilonByCell()[ci]-ep,ep,1e-12L,"E1 independent raw-stock epsilon");
        bounded(field.conductivityByCell()[ci]-sg,sg,1e-12L,"E1 independent raw-stock conductivity");
        sourceTotal+=old[ci];targetTotal+=next[ci];sourceAbs+=std::abs(old[ci]);targetAbs+=std::abs(next[ci]);
        if(analytic)potentialError=std::max(potentialError,std::abs(static_cast<long double>(potential[ci])-affine(mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord(),domain.axisymmetric())));
    }
    for(std::size_t row=0;row<reference.size();++row){const auto& f=field.faces()[row];const auto& ref=reference[row];
        check(f.nativeRow==static_cast<int>(row)&&f.negativeCell==ref.face.negativeCell&&f.positiveCell==ref.face.positiveCell,"E1 canonical native orientation changed");
        bounded(f.physicalMeasure-ref.area,ref.area,1e-12L,"E1 actual XY/RZ native-face measure");
        bounded(f.epsilon-ref.epsilon,ref.epsilon,1e-12L,"E1 epsilon half-distance harmonic");
        bounded(f.conductivity-ref.sigma,ref.sigma,1e-12L,"E1 conductivity half-distance harmonic");
        bounded(f.implicitCoefficient-ref.a,ref.a,1e-12L,"E1 a differs from epsilonFace+dt*sigmaFace");
        long double h=f.gradientConstant;for(const auto& coefficient:f.gradientCoefficients){
            check(coefficient.nativeCell>=0&&coefficient.nativeCell<mesh.cellIndexNum(),"E1 public H coefficient references foreign cell");
            h+=static_cast<long double>(coefficient.value)*potential[coefficient.nativeCell];}
        // This checks the same-H public receipt; the affine oracle below is independent.
        bounded(h-f.integratedField,std::abs(h)+std::abs(f.integratedField)+ref.area,1e-12L,"E1 H/stencil receipt mismatch");
        bounded(f.displacement-static_cast<long double>(ref.epsilon)*f.integratedField,
            std::abs(f.displacement)+std::abs(static_cast<long double>(ref.epsilon)*f.integratedField),1e-12L,"E1 D does not consume same H");
        bounded(f.current-static_cast<long double>(ref.sigma)*f.integratedField,
            std::abs(f.current)+std::abs(static_cast<long double>(ref.sigma)*f.integratedField),1e-12L,"E1 J does not consume same H");
        if(ref.area==0)check(f.integratedField==0&&f.displacement==0&&f.current==0,"E1 true axis flux is not exact zero");
        if(analytic&&ref.area>0){const long double scale=ref.area*(domain.axisymmetric()?1.:std::hypot(1.,.37));
            errors.affineField=std::max(errors.affineField,std::abs(static_cast<long double>(f.integratedField)-ref.affineH)/scale);}
        scatter(D,ref.face,f.displacement);scatter(J,ref.face,f.current);
        if(ref.face.negativeCell>=0){absD[ref.face.negativeCell]+=std::abs(f.displacement);absJ[ref.face.negativeCell]+=std::abs(f.current);}
        if(ref.face.positiveCell>=0){absD[ref.face.positiveCell]+=std::abs(f.displacement);absJ[ref.face.positiveCell]+=std::abs(f.current);}
        if(ref.face.negativeCell<0||ref.face.positiveCell<0){const int sign=ref.face.negativeCell>=0?1:-1;
            boundaryJ+=sign*static_cast<long double>(f.current);boundaryAbsJ+=std::abs(f.current);}
        if(zeroSigma)check(f.current==0,"E1 zero conductivity produced current");
    }
    for(std::size_t ci=0;ci<old.size();++ci){const long double dg=D[ci]-next[ci],dq=static_cast<long double>(next[ci])-old[ci]+dt*J[ci];
        const long double gs=std::abs(next[ci])+absD[ci],qs=std::abs(old[ci])+std::abs(next[ci])+dt*absJ[ci];
        bounded(dg,gs,gate,"E1 independent actual-face local Gauss");bounded(dq,qs,gate,"E1 independent integrated-Q conduction equation");
        if(gs>0)errors.gauss=std::max(errors.gauss,std::abs(dg)/gs);if(qs>0)errors.charge=std::max(errors.charge,std::abs(dq)/qs);
    }
    const long double global=targetTotal-sourceTotal+dt*boundaryJ,scale=sourceAbs+std::abs(targetTotal)+dt*boundaryAbsJ;
    bounded(global,scale,gate,"E1 independent signed total charge/boundary current");if(scale>0)errors.global=std::abs(global)/scale;
    backendReceipt(field,backend);
    check(field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=gate,"E1 retained linear/Gauss gates after independent conservation checks");
    check(field.report().nonorthogonalFaces>0,"E1 real split-face nonorthogonal path not consumed");
    if(domain.axisymmetric())check(field.report().axisFaces>0,"E1 RZ witness has no actual axis faces");
    if(analytic){errors.potential=potentialError;check(potentialError<=1e-10&&errors.affineField<=1e-10L,"E1 independent affine potential/physical H patch");}
    if(zeroSigma)check(std::memcmp(old.data(),next.data(),old.size()*sizeof(double))==0,"E1 zero sigma changed an integrated charge bit");
    return errors;
}
void run() {
    std::size_t solved=0,guards=0;
    for(bool rz:{false,true}){
        auto mesh=meshFor();const auto domain=domainFor(*mesh,rz);check(!domain.uniformLeaves(),"E1 fixture lacks true2:1 native AMR");
        const auto stocks=phaseStocks(*mesh,domain);const Clock acceptedClock{0,.2,.2,Clock::Position::Initial};
        const Clock geometryClock{1,.2,.22,Clock::Position::Midpoint};
        const Clock candidate{1,.2,.22,rz?Clock::Position::Midpoint:Clock::Position::End};
        const auto accepted=Publication::prepareLocalPlicInventory(*mesh,domain,acceptedClock,stocks,{});
        const auto geometry=Publication::prepareLocalPlicInventory(*mesh,domain,geometryClock,stocks,{});const auto saved=snapshot(geometry);
        std::vector<double> arbitrary(mesh->cellIndexNum());bool positive=false,negative=false;
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto center=mesh->cells(mesh->cellIndexToLabel[ci]).centerCoord();
            arbitrary[ci]=domain.cellVolumes()[domain.cellRowByIndex()[ci]]*(center.x-.5+.19*(center.y-.5));positive|=arbitrary[ci]>0;negative|=arbitrary[ci]<0;}
        check(positive&&negative,"E1 independent Qi fixture must contain both signs");
        for(int kind=0;kind<3;++kind){ElectricVofMaterials material{2.,5.,.7,.03};if(kind==2)material.conductivityLiquid=material.conductivityGas=0;
            std::vector<Boundary> boundaries;const double dt=candidate.time()-acceptedClock.time();
            const auto reference=references(*mesh,geometry,material,dt,boundaries);std::vector<double> q=arbitrary;
            if(kind==0){std::vector<long double> analytic(q.size());for(const auto& face:reference)scatter(analytic,face.face,static_cast<long double>(face.a)*face.affineH);
                for(std::size_t ci=0;ci<q.size();++ci)q[ci]=static_cast<double>(analytic[ci]);}
            const auto charge=ElectricVofCharge::fromIntegrated(*mesh,accepted,q);check(charge.identity()!=0&&charge.integratedCharge()==q,"E1 signed source charge changed");
            const auto field=ElectricVofField::prepareOhmic(*mesh,geometry,charge,candidate,material,boundaries);
            const auto metric=verify(*mesh,geometry,charge,field,candidate,material,reference,kind==0,kind==2);++solved;
            check(charge.integratedCharge()==q&&snapshot(geometry)==saved,"E1 successful candidate mutated original Q/VOF stocks/plane");
            if(kind==1){
                const auto other=ElectricVofCharge::fromIntegrated(*mesh,accepted,q);
                check(!charge.sameInstance(other),"E1 distinct source charge factory reused instance");
                rejected([&]{field.validateFor(*mesh,geometry,other,candidate);},"E1 foreign charge owner accepted");++guards;
                const auto otherGeometry=Publication::prepareLocalPlicInventory(*mesh,domain,geometryClock,stocks,{});
                rejected([&]{field.validateFor(*mesh,otherGeometry,charge,candidate);},"E1 foreign frozen VOF owner accepted");++guards;
                auto wrong=candidate;++wrong.step;check(wrong.time()==candidate.time(),"E1 clock negative must retain same physical time");
                rejected([&]{field.validateFor(*mesh,geometry,charge,wrong);},"E1 same-time foreign full clock accepted");++guards;
                rejected([&]{ElectricVofField::prepareOhmic(*mesh,geometry,charge,wrong,material,boundaries);},"E1 wrong interval/geometry fullclock accepted");++guards;
                auto wrongSource=acceptedClock;++wrongSource.step;
                rejected([&]{charge.validateFor(*mesh,wrongSource);},"E1 source accepted fullclock ignored");++guards;
                auto foreignMesh=meshFor(4);
                rejected([&]{field.validateFor(*foreignMesh,geometry,charge,candidate);},"E1 different native geometry accepted");++guards;
                rejected([&]{ElectricVofCharge::fromIntegrated(*mesh,geometry,q);},"E1 Midpoint geometry accepted as accepted source charge");++guards;
                auto badMaterial=material;badMaterial.epsilonLiquid=0;
                rejected([&]{ElectricVofField::prepareOhmic(*mesh,geometry,charge,candidate,badMaterial,boundaries);},"E1 nonpositive epsilon accepted");++guards;
                badMaterial=material;badMaterial.conductivityGas=-1;
                rejected([&]{ElectricVofField::prepareOhmic(*mesh,geometry,charge,candidate,badMaterial,boundaries);},"E1 negative conductivity accepted");++guards;
                ElectricVofControls relaxed;relaxed.relativeConservationTolerance=1e-8;
                rejected([&]{ElectricVofField::prepareOhmic(*mesh,geometry,charge,candidate,material,boundaries,relaxed);},"E1 relaxed original conservation gate accepted");++guards;
                auto badBoundary=boundaries;for(std::size_t row=0;row<reference.size();++row)if(reference[row].face.negativeCell>=0&&reference[row].face.positiveCell>=0){badBoundary[row]={Type::Dirichlet,0};break;}
                rejected([&]{ElectricVofField::prepareOhmic(*mesh,geometry,charge,candidate,material,badBoundary);},"E1 interior face boundary accepted");++guards;
                if(rz){badBoundary=boundaries;for(std::size_t row=0;row<reference.size();++row)if(reference[row].area==0){badBoundary[row]={Type::Dirichlet,0};break;}
                    rejected([&]{ElectricVofField::prepareOhmic(*mesh,geometry,charge,candidate,material,badBoundary);},"E1 zero-axis Dirichlet accepted");++guards;}
                check(charge.integratedCharge()==q&&snapshot(geometry)==saved,"E1 failed candidates mutated accepted charge/VOF owner");
                field.validateFor(*mesh,geometry,charge,candidate);
            }
            std::cout<<"ELECTRIC_VOF_OHMIC,"<<(rz?"RZ":"XY")<<",kind,"<<kind<<",cells,"<<mesh->cellIndexNum()<<",native_faces,"<<reference.size()
                <<",nonorthogonal,"<<field.report().nonorthogonalFaces<<",axis,"<<field.report().axisFaces<<",geometry_time,"<<field.geometryTime().time()
                <<",source_time,"<<field.chargeSourceTime().time()<<",candidate_time,"<<field.candidateTime().time()<<",conduction_dt,"<<field.effectiveDt()
                <<",independent_Gauss,"<<metric.gauss<<",charge_equation,"<<metric.charge<<",boundary_total,"<<metric.global
                <<",affine_H,"<<metric.affineField<<",affine_potential,"<<metric.potential<<",Maxwell,NOT_CONNECTED\n";
        }
    }
    check(solved==6&&guards==23,"E1 minimum signed-Q/clock/geometry evidence incomplete");
    std::cout<<"ELECTRIC_VOF_OHMIC_PASS,solves,"<<solved<<",guards,"<<guards<<",Gauss_charge_gate,1e-10,no_charge_advection,no_AMR_charge_transfer,no_Maxwell\n";
}
} // namespace ElectricVofContract

namespace ElectricVofBinary64Contract {
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using D=FV::Dyadic;
unsigned nearestEven(unsigned n,unsigned d) {
    check(d>0,"E2P independent integer denominator");const unsigned q=n/d,r=n%d;
    return q+(2*r>d||(2*r==d&&(q&1))?1:0);
}
void receipt(const FV::PartitionDouble& p,double source) {
    check(p.values.size()==p.exactAllocations.size(),"E2P binary64 allocation receipt sizes");D sum;
    for(std::size_t k=0;k<p.values.size();++k){
        check(std::isfinite(p.values[k])&&p.values[k]>=0,"E2P finite nonnegative binary64 piece");
        check(D(p.values[k]).compare(p.exactAllocations[k])==0,"E2P stored double differs from its exact allocation");sum.add(p.exactAllocations[k]);}
    check(sum.compare(D(source))==0&&p.sumError.compare(D{})==0,"E2P exact own-donor telescoping lost charge magnitude");
}
void tinyOracle(unsigned units,const std::vector<unsigned>& integerWeights,const std::vector<D>& weights) {
    const double source=std::scalbn(static_cast<double>(units),-1074);const auto partition=FV::partitionDouble(source,weights);
    receipt(partition,source);unsigned total=0;for(auto w:integerWeights)total+=w;
    unsigned prefix=0,prior=0;check(integerWeights.size()==partition.values.size(),"E2P oracle weight count");
    for(std::size_t k=0;k<integerWeights.size();++k){prefix+=integerWeights[k];const unsigned current=nearestEven(units*prefix,total);
        const double expected=std::scalbn(static_cast<double>(current-prior),-1074);
        check(partition.values[k]==expected,"E2P independent binary64 integer-CDF piece oracle");prior=current;}
    check(prior==units,"E2P independent integer endpoint");
    for(int sign:{-1,1}){D signedSum;for(double q:partition.values)signedSum.add(static_cast<long double>(sign)*q);
        const auto result=FV::roundDouble(signedSum);check(result.value==sign*source&&result.error.compare(D{})==0,
            "E2P signed donor allocation/coarsen did not preserve source");}
}
template<class F>void reject(F f,const char* reason) {bool threw=false;try{f();}catch(const std::exception&){threw=true;}check(threw,reason);}
void run() {
    const double tiny=std::numeric_limits<double>::denorm_min();
    std::cout<<std::setprecision(17)<<"E2P_BINARY64_MODEL,double_digits,"<<std::numeric_limits<double>::digits
        <<",long_double_digits,"<<std::numeric_limits<long double>::digits<<",double_min_exponent,"<<std::numeric_limits<double>::min_exponent
        <<",double_size,"<<sizeof(double)<<",double_dmin,"<<std::hexfloat<<tiny<<std::defaultfloat<<'\n';
    for(unsigned units:{1u,3u,7u})tinyOracle(units,{1,1,1,1},{D(1),D(1),D(1),D(1)});
    // Same actual D16 RZ parent, split into four real rectangles. Independent
    // radial-square differences are45:47 and the two z widths are equal.
    const double rs[3]={11./64,23./128,3./16},zs[3]={5./16,41./128,21./64};std::vector<D> rz;
    for(int x=0;x<2;++x)for(int z=0;z<2;++z)rz.push_back(FV::rectangleCapacityExact({{rs[x],zs[z]},{rs[x+1],zs[z+1]}},true));
    for(unsigned units:{1u,3u,7u})tinyOracle(units,{45,45,47,47},rz);
    const auto repeated=FV::partitionDouble(1.,{D(0),D(1),D(0),D(3),D(0)});receipt(repeated,1.);
    check(repeated.values==std::vector<double>({0,.25,0,.75,0}),"E2P zero intervals/actual endpoint rather than array-last allocation");
    const auto zero=FV::partitionDouble(0.,{D(1),D(0)});receipt(zero,0.);check(zero.values==std::vector<double>({0,0}),"E2P zero-source pieces");
    const double maximum=std::numeric_limits<double>::max();receipt(FV::partitionDouble(maximum,{D(1),D(0)}),maximum);
    for(double factor:{.5,1.5,2.5,-.5,-1.5}){
        const D exact=D::product({tiny,factor});const auto result=FV::roundDouble(exact);
        const double expected=factor==.5||factor==-.5?std::copysign(0.,factor):std::copysign(2*tiny,factor);
        check(result.value==expected&&std::signbit(result.value)==std::signbit(expected),"E2P independent subnormal ties-to-even/sign oracle");
        D identity=result.error;identity.add(exact);identity.add(result.value,-1);check(identity.compare(D{})==0,"E2P signed rounding error value-minus-exact");}
    D cancellation(std::scalbn(3.,-1074));cancellation.add(std::scalbn(2.,-1074),-1);
    check(FV::roundDouble(cancellation).value==tiny,"E2P exact signed receiver cancellation lost one double quantum");
    D half=D::product({1.L,std::scalbn(1.L,-53)}),tie(1);tie.add(half);
    check(FV::roundDouble(tie).value==1.,"E2P lower-even normal midpoint");
    const double odd=std::nextafter(1.,2.);D upperTie(odd);upperTie.add(half);
    check(FV::roundDouble(upperTie).value==std::nextafter(odd,2.),"E2P upper-even normal midpoint");
    // An exact bit beyond native-LD precision separates direct binary64 RNE
    // from forbidden native-LD-round-then-double-cast on wider platforms.
    D above=tie;above.add(D::product({std::scalbn(1.L,-std::numeric_limits<long double>::digits),std::scalbn(1.L,-64)}));
    check(FV::roundDouble(above).value==odd&&above.roundedDouble()==odd,"E2P direct binary64 RNE double-rounding counterexample");
    if(std::numeric_limits<long double>::digits>53)check(static_cast<double>(above.rounded())==1.,"E2P wider-native double-rounding witness not exposed");
    const double minimum=std::numeric_limits<double>::min();D transition(minimum);transition.add(D::product({tiny,.5L}),-1);
    check(FV::roundDouble(transition).value==minimum,"E2P normal/subnormal even-boundary oracle");
    D belowOverflow(maximum);belowOverflow.add(D(std::scalbn(1.L,969)));
    check(FV::roundDouble(belowOverflow).value==maximum,"E2P finite near-overflow rounding");
    D overflow(maximum);overflow.add(D(std::scalbn(1.L,970)));
    reject([&]{FV::roundDouble(overflow);},"E2P final binary64 overflow accepted");
    reject([&]{FV::partitionDouble(-tiny,{D(1)});},"E2P negative magnitude accepted");
    reject([&]{FV::partitionDouble(1.,{D(0)});},"E2P positive source with zero geometry accepted");
    reject([&]{FV::partitionDouble(0.,{D(-1)});},"E2P invalid geometry ignored for zero source");
    reject([&]{FV::partitionDouble(std::numeric_limits<double>::infinity(),{D(1)});},"E2P infinite source accepted");
    // Original native-LD public entries remain real consumers, not aliases.
    const auto old=FV::partition(3.L,{D(1),D(1)});check(old.values==std::vector<long double>({1.5L,1.5L}),"E2P original native-LD partition changed");
    const auto oldRound=FV::round(D(1.5L));check(oldRound.value==1.5L&&oldRound.error.compare(D{})==0,"E2P original native-LD round changed");
    std::cout<<"E2P_BINARY64_MATH_PASS,integer_CDF_XY_RZ,min_odd_subnormal,signed_cancellation,ties_normal_transition,overflow_rejected,old_LD_preserved\n";
}
} // namespace ElectricVofBinary64Contract

namespace ElectricVofMovingContract {
using namespace ElectricVofContract;
using ChargeTransport=LAMP::Electric::ElectricVofChargeTransport;
using ChargeRemap=LAMP::Electric::ElectricVofChargeMeshTransfer;
using VofTransport=I::GeometricVofTransport;
using Slab=I::GeometricVofAdvectorSlab;
using VofTransfer=I::GeometricVofMeshTransfer;
std::vector<double> normalFlow(const QuadtreeMesh& mesh,const Publication& p,double dt) {
    const auto psi=[](Vector2D v){const long double x=v.x,y=v.y;return x*x*(1-x)*(1-x)*y*y*(1-y)*(1-y);};
    std::vector<double> result(p.domain().nativeFaceCount());
    for(int row=0;row<p.domain().nativeFaceCount();++row){const auto f=p.domain().face(mesh,row);
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const auto volume=(f.horizontal?-1:1)*(psi(b)-psi(a));
        if(f.outwardSign!=0||f.area==0){check(volume==0,"E2 actual closed streamfunction leaks at boundary/axis");result[row]=0;}
        else result[row]=static_cast<double>(volume/f.area);}
    const auto rates=VofTransport::localPlicCourantRates(mesh,p,result);
    const auto peak=std::max({rates.maximumLocalColumn,rates.maximumIncoming,rates.maximumOutflow});check(peak>0,"E2 moving fixture has zero actual flow");
    for(auto& u:result)u=static_cast<double>(u*(.25L/(dt*peak)));return result;
}
void verifyAdvection(const QuadtreeMesh& mesh,const LAMP::Electric::ElectricVofAdvectedCharge& a,bool half) {
    const auto& q=a.sourceCharge().integratedCharge();const auto& next=a.integratedCharge();const auto& domain=a.transport().source().domain();
    a.validateFor(mesh,a.sourceCharge(),a.transport(),a.timeLayer());check(a.geometry().identity()==(half?a.transport().midpoint():a.transport().endpoint()).identity(),"E2 receipt geometry is not actual transport stage");
    std::vector<long double> net(q.size());long double oldSum=0,newSum=0,l1=0,faceL1=0;std::size_t moving=0,opposite=0;
    for(std::size_t row=0;row<a.faces().size();++row){const auto& f=a.faces()[row];const auto& vf=a.transport().faces()[row];const auto t=half?vf.half.total:vf.full.total;
        check(f.nativeRow==vf.nativeRow&&f.negativeCell==vf.negativeCell&&f.positiveCell==vf.positiveCell&&f.totalVolume==t,"E2 did not consume actual shared total T");
        if(f.negativeCell<0||f.positiveCell<0)check(t==0&&f.integratedCharge==0,"E2 closed exterior charge flux");
        else if(t!=0){++moving;const int donor=t>0?f.negativeCell:f.positiveCell;const double magnitude=std::abs(q[donor]);
            const long double expected=t*q[donor]/domain.cellVolumes()[domain.cellRowByIndex()[donor]];
            int e=0;std::frexp(magnitude,&e);const double quantum=magnitude==0?0:std::scalbn(1.,std::max(e-53,-1074));
            check(std::abs(static_cast<long double>(f.integratedCharge)-expected)<=2.L*quantum,"E2 actual donor P0/CDF error exceeds one parent-quantum scale");
            if(q[donor]<0&&f.integratedCharge!=0){check((f.integratedCharge>0)!=(t>0),"E2 donor was selected by charge sign");++opposite;}}
        if(f.negativeCell>=0)net[f.negativeCell]+=f.integratedCharge;if(f.positiveCell>=0)net[f.positiveCell]-=f.integratedCharge;faceL1+=std::abs(f.integratedCharge);}
    for(std::size_t i=0;i<q.size();++i){const auto scale=std::abs(q[i])+std::abs(next[i])+std::abs(net[i]);
        bounded(static_cast<long double>(next[i])-q[i]+net[i],scale,gate,"E2 independent cell signed charge advection");oldSum+=q[i];newSum+=next[i];l1+=std::abs(q[i])+std::abs(next[i]);}
    bounded(newSum-oldSum,l1+faceL1,gate,"E2 independent closed signed global charge");
    bounded(a.report().chargeResidual-(newSum-oldSum),l1+faceL1,1e-12L,"E2 actual transport report global charge");
    check(a.report().relativeChargeResidual<=gate&&a.report().maximumOutflowCfl<=.5,"E2 original physical charge/Courant gates");
    check(moving>0,"E2 advection did not consume moving native faces");
    if(a.report().nonzeroFaceCharges>8)check(opposite>0,"E2 signed-Q donor witness absent");
}
void verifyConduction(const QuadtreeMesh& mesh,const LAMP::Electric::ElectricVofAdvectedCharge& adv,const ElectricVofField& field,bool zeroSigma) {
    const auto& old=adv.sourceCharge().integratedCharge();const auto& input=adv.integratedCharge();const auto& next=field.chargeCandidate().integratedCharge();
    check(field.advection()&&field.advection()->sameInstance(adv)&&field.sourceCharge().sameInstance(adv.sourceCharge()),"E2 field lost actual accepted/advection owners");
    check(field.implicitInputCharge()==input&&sameClock(field.chargeSourceTime(),adv.sourceCharge().timeLayer())&&sameClock(field.candidateTime(),adv.timeLayer()),"E2 advected End relabeled as Begin");
    const double dt=adv.timeLayer().time()-adv.sourceCharge().timeLayer().time();check(field.effectiveDt()==dt,"E2 conduction duration did not use original accepted source");
    std::vector<long double> D(old.size()),J(old.size()),aD(old.size()),aJ(old.size());long double boundary=0,boundaryAbs=0,source=0,rhs=0,target=0,l1=0;
    for(const auto& f:field.faces()){const auto actual=adv.geometry().domain().face(mesh,f.nativeRow);scatter(D,actual,f.displacement);scatter(J,actual,f.current);
        if(actual.negativeCell>=0){aD[actual.negativeCell]+=std::abs(f.displacement);aJ[actual.negativeCell]+=std::abs(f.current);}
        if(actual.positiveCell>=0){aD[actual.positiveCell]+=std::abs(f.displacement);aJ[actual.positiveCell]+=std::abs(f.current);}
        if(actual.negativeCell<0||actual.positiveCell<0){const int sign=actual.negativeCell>=0?1:-1;boundary+=sign*static_cast<long double>(f.current);boundaryAbs+=std::abs(f.current);}}
    for(std::size_t i=0;i<old.size();++i){bounded(D[i]-next[i],std::abs(next[i])+aD[i],gate,"E2 actual native Gauss after advection");
        bounded(static_cast<long double>(next[i])-input[i]+dt*J[i],std::abs(next[i])+std::abs(input[i])+dt*aJ[i],gate,"E2 Ohmic update does not use Qadv RHS");
        source+=old[i];rhs+=input[i];target+=next[i];l1+=std::abs(old[i])+std::abs(input[i])+std::abs(next[i]);}
    bounded(target-rhs+dt*boundary,l1+dt*boundaryAbs,gate,"E2 independent Qadv/Ohmic boundary ledger");
    bounded(target-source+dt*boundary,l1+dt*boundaryAbs,gate,"E2 accepted/advection/Ohmic total ledger");
    check(field.report().usedAdvection,"E2 actual advection request missing after independent Ohmic/combined checks");
    if(zeroSigma)check(std::memcmp(input.data(),next.data(),input.size()*sizeof(double))==0,"E2 sigma0 did not retain Qadv bits");
}
unsigned nearest(unsigned n,unsigned d){const unsigned q=n/d,r=n%d;return q+(2*r>d||(2*r==d&&(q&1))?1:0);}
unsigned rectangleWeight(const LAMP::Geometry::PlanarCellCut::Rectangle& b,bool rz) {
    const auto integer=[](long double x){const long double scaled=x*16;const auto n=static_cast<unsigned>(std::llround(scaled));check(scaled==n,"E2 independent rectangle is not actual1/16 lattice");return n;};
    const auto x0=integer(b.lower.x),x1=integer(b.upper.x),y0=integer(b.lower.y),y1=integer(b.upper.y);
    return (rz?x1*x1-x0*x0:x1-x0)*(y1-y0);
}
void tinyRemapOracle(const ChargeRemap& result) {
    const auto& actual=result.interfaceTransfer().nativeOverlaps();const auto& source=result.sourceCharge().integratedCharge();const bool rz=result.sourceGeometry().domain().axisymmetric();
    std::vector<long long> expected(actual.byTarget().size());const double unit=std::numeric_limits<double>::denorm_min();
    for(std::size_t si=0;si<source.size();++si){auto group=actual.bySource()[si];std::sort(group.begin(),group.end(),[&](auto a,auto b){const auto& x=actual.pieces()[a].rectangle;const auto& y=actual.pieces()[b].rectangle;
            return std::tie(x.lower.x,x.lower.y,x.upper.x,x.upper.y)<std::tie(y.lower.x,y.lower.y,y.upper.x,y.upper.y);});
        const unsigned units=static_cast<unsigned>(std::llround(std::abs(source[si])/unit));check(units<=7&&units*unit==std::abs(source[si]),"E2 tiny oracle source not explicit integer double quantum");
        unsigned total=0;for(auto pi:group)total+=rectangleWeight(actual.pieces()[pi].rectangle,rz);unsigned prefix=0,prior=0;
        for(auto pi:group){prefix+=rectangleWeight(actual.pieces()[pi].rectangle,rz);const unsigned current=nearest(units*prefix,total);const int sign=std::signbit(source[si])?-1:1;
            expected[actual.pieces()[pi].targetCell]+=sign*static_cast<long long>(current-prior);prior=current;}check(prior==units,"E2 independent integer source endpoint");}
    const auto& q=result.chargeCandidate().integratedCharge();for(std::size_t ti=0;ti<q.size();++ti)
        check(q[ti]==expected[ti]*unit,"E2 actual AMR differs from independent XY/RZ integer CDF/coarsen oracle");
    check(result.report().chargeResidual==0&&result.report().roundingResidual==0,"E2 tiny AMR lost a binary64 quantum");
}
void run() {
    std::size_t solves=0,guards=0,remaps=0;
    for(bool rz:{false,true}){
        auto mesh=meshFor();auto coarse=mesh->clone();const auto domain=domainFor(*mesh,rz);const auto stocks=phaseStocks(*mesh,domain);
        const Clock start{0,.2,.2,Clock::Position::Initial},end{1,.2,.202,Clock::Position::End};
        const auto geometry=Publication::prepareLocalPlicInventory(*mesh,domain,start,stocks,{});const auto saved=snapshot(geometry);
        std::vector<double> q(mesh->cellIndexNum());for(int i=0;i<mesh->cellIndexNum();++i){const auto c=mesh->cells(mesh->cellIndexToLabel[i]).centerCoord();q[i]=domain.cellVolumes()[domain.cellRowByIndex()[i]]*(c.x-.5+.19*(c.y-.5));}
        const auto charge=ElectricVofCharge::fromIntegrated(*mesh,geometry,q);const auto normals=normalFlow(*mesh,geometry,end.end-end.begin);
        const auto slab=Slab::prepare(*mesh,geometry,end,normals,normals);const auto flow=VofTransport::prepareLocalPlicSweep(*mesh,geometry,slab,1000,1.226);
        const auto moved=ChargeTransport::prepareClosed(*mesh,charge,flow);verifyAdvection(*mesh,moved.midpoint(),true);verifyAdvection(*mesh,moved.endpoint(),false);
        std::vector<long double> tinyOut(q.size());
        for(const auto& f:flow.faces())if(f.full.total!=0){const int donor=f.full.total>0?f.negativeCell:f.positiveCell;
            check(donor>=0,"E2 tiny donor outside closed domain");tinyOut[donor]+=std::abs(f.full.total)/domain.cellVolumes()[domain.cellRowByIndex()[donor]];}
        const auto tinyDonor=static_cast<std::size_t>(std::max_element(tinyOut.begin(),tinyOut.end())-tinyOut.begin());
        std::vector<double> advTiny(q.size());advTiny[tinyDonor]=-7*std::numeric_limits<double>::denorm_min();
        const auto advTinySource=ElectricVofCharge::fromIntegrated(*mesh,geometry,advTiny);const auto advTinyMoved=ChargeTransport::prepareClosed(*mesh,advTinySource,flow);
        verifyAdvection(*mesh,advTinyMoved.midpoint(),true);verifyAdvection(*mesh,advTinyMoved.endpoint(),false);
        long long tinyUnits=0;for(const auto value:advTinyMoved.endpoint().integratedCharge())tinyUnits+=std::llround(value/std::numeric_limits<double>::denorm_min());
        check(tinyUnits==-7&&advTinyMoved.endpoint().report().nonzeroFaceCharges>0,"E2 actual signed tiny donor was deleted or did not reach any real moving face");
        const ElectricVofMaterials material{2,5,.7,.03},zeroSigma{2,5,0,0};std::vector<Boundary> boundary;
        references(*mesh,flow.midpoint(),material,moved.midpoint().timeLayer().time()-start.time(),boundary);
        const auto half=ElectricVofField::prepareOhmic(*mesh,flow.midpoint(),moved.midpoint(),material,boundary);verifyConduction(*mesh,moved.midpoint(),half,false);++solves;
        references(*mesh,flow.midpoint(),zeroSigma,end.time()-start.time(),boundary);
        const auto full=ElectricVofField::prepareOhmic(*mesh,flow.midpoint(),moved.endpoint(),zeroSigma,boundary);verifyConduction(*mesh,moved.endpoint(),full,true);++solves;
        const auto zero=std::vector<double>(normals.size());const auto stillSlab=Slab::prepare(*mesh,geometry,end,zero,zero);
        const auto stillFlow=VofTransport::prepareLocalPlicSweep(*mesh,geometry,stillSlab,1000,1.226);
        auto bitQ=q;bitQ[0]=-0.;const auto bitCharge=ElectricVofCharge::fromIntegrated(*mesh,geometry,bitQ);const auto still=ChargeTransport::prepareClosed(*mesh,bitCharge,stillFlow);
        for(const auto* a:{&still.midpoint(),&still.endpoint()})check(std::memcmp(bitQ.data(),a->integratedCharge().data(),q.size()*sizeof(double))==0,"E2 zero-flow changed charge bits");
        auto wrong=moved.endpoint().timeLayer();++wrong.step;rejected([&]{moved.endpoint().validateFor(*mesh,charge,flow,wrong);},"E2 same-time wrong fullclock accepted");++guards;
        const auto foreignCharge=ElectricVofCharge::fromIntegrated(*mesh,geometry,q);rejected([&]{moved.endpoint().validateFor(*mesh,foreignCharge,flow,end);},"E2 foreign accepted charge owner accepted");++guards;
        const auto otherFlow=VofTransport::prepareLocalPlicSweep(*mesh,geometry,slab,1000,1.226);
        rejected([&]{moved.endpoint().validateFor(*mesh,charge,otherFlow,end);},"E2 foreign transport owner accepted");++guards;
        const auto fresh=Publication::prepareLocalPlicInventory(*mesh,domain,flow.midpoint().timeLayer(),stocks,{});
        rejected([&]{ElectricVofField::prepareOhmic(*mesh,fresh,moved.endpoint(),zeroSigma,boundary);},"E2 fresh same-clock geometry self-signed advection");++guards;
        const auto endQ=ElectricVofCharge::fromIntegrated(*mesh,flow.endpoint(),moved.endpoint().integratedCharge());
        rejected([&]{ElectricVofField::prepareOhmic(*mesh,flow.midpoint(),endQ,end,zeroSigma,boundary);},"E2 End Qadv was reused as this interval Begin");++guards;
        auto fine=mesh->clone();DictionaryMesh::RefineBox box;box.min={0,0};box.max={1,1};box.level=1;fine->applyRefineBoxes({box});check(fine->cellIndexNum()>mesh->cellIndexNum(),"E2 actual refine changed no native leaves");
        const auto refine=VofTransfer::prepare(*mesh,flow.endpoint(),*fine,domainFor(*fine,rz),end);
        const auto remap=ChargeRemap::prepare(*mesh,flow.endpoint(),full.chargeCandidate(),*fine,refine);++remaps;
        remap.validateFor(*mesh,flow.endpoint(),full.chargeCandidate(),*fine,refine);check(remap.report().refinedDonors>0&&remap.report().unchangedTargets>0,"E2 refine/unchanged paths not real");
        std::vector<double> tiny(q.size());int a=-1,b=-1;
        for(int i=0;i<mesh->cellIndexNum();++i)if(refine.nativeOverlaps().bySource()[i].size()>1){if(a<0)a=i;else if(b<0)b=i;}
        check(a>=0&&b>=0,"E2 true refine lacks two source donors");const auto unit=std::numeric_limits<double>::denorm_min();tiny[a]=7*unit;tiny[b]=-unit;
        const auto tinyCharge=ElectricVofCharge::fromIntegrated(*mesh,flow.endpoint(),tiny);const auto tinyFine=ChargeRemap::prepare(*mesh,flow.endpoint(),tinyCharge,*fine,refine);tinyRemapOracle(tinyFine);++remaps;
        const Clock next{2,end.end,end.end+.002,Clock::Position::End};const auto normal2=normalFlow(*fine,refine.candidate(),next.end-next.begin);
        const auto slab2=Slab::prepare(*fine,refine.candidate(),next,normal2,normal2);const auto flow2=VofTransport::prepareLocalPlicSweep(*fine,refine.candidate(),slab2,1000,1.226);
        const auto moved2=ChargeTransport::prepareClosed(*fine,remap.chargeCandidate(),flow2);verifyAdvection(*fine,moved2.endpoint(),false);
        references(*fine,flow2.midpoint(),zeroSigma,next.time()-end.time(),boundary);
        const auto field2=ElectricVofField::prepareOhmic(*fine,flow2.midpoint(),moved2.endpoint(),zeroSigma,boundary);verifyConduction(*fine,moved2.endpoint(),field2,true);++solves;
        const auto coarsen=VofTransfer::prepare(*fine,flow2.endpoint(),*coarse,domainFor(*coarse,rz),next);
        const auto back=ChargeRemap::prepare(*fine,flow2.endpoint(),field2.chargeCandidate(),*coarse,coarsen);++remaps;
        check(back.report().coarsenedTargets>0&&back.report().unchangedTargets>0,"E2 actual coarsen/unchanged paths absent");
        const auto tinyCoarsen=VofTransfer::prepare(*fine,refine.candidate(),*coarse,domainFor(*coarse,rz),end);
        const auto tinyBack=ChargeRemap::prepare(*fine,refine.candidate(),tinyFine.chargeCandidate(),*coarse,tinyCoarsen);tinyRemapOracle(tinyBack);++remaps;
        check(tinyBack.chargeCandidate().integratedCharge()==tiny,"E2 signed minimum/odd charge refine/coarsen cycle changed source");
        const Clock third{3,next.end,next.end+.002,Clock::Position::End};const auto normal3=normalFlow(*coarse,coarsen.candidate(),third.end-third.begin);
        const auto slab3=Slab::prepare(*coarse,coarsen.candidate(),third,normal3,normal3);const auto flow3=VofTransport::prepareLocalPlicSweep(*coarse,coarsen.candidate(),slab3,1000,1.226);
        const auto moved3=ChargeTransport::prepareClosed(*coarse,back.chargeCandidate(),flow3);references(*coarse,flow3.midpoint(),zeroSigma,third.time()-next.time(),boundary);
        const auto field3=ElectricVofField::prepareOhmic(*coarse,flow3.midpoint(),moved3.endpoint(),zeroSigma,boundary);verifyConduction(*coarse,moved3.endpoint(),field3,true);++solves;
        rejected([&]{back.chargeCandidate().validateFor(*fine,next);},"E2 source-size/foreign native mesh accepted after remap");++guards;
        rejected([&]{remap.validateFor(*mesh,flow.endpoint(),tinyCharge,*fine,refine);},"E2 charge remap foreign source owner accepted");++guards;
        check(charge.integratedCharge()==q&&snapshot(geometry)==saved,"E2 success/failure mutated accepted Q/VOF");
        std::cout<<"ELECTRIC_VOF_MOVING,"<<(rz?"RZ":"XY")<<",source_cells,"<<mesh->cellIndexNum()<<",fine_cells,"<<fine->cellIndexNum()<<",coarse_cells,"<<coarse->cellIndexNum()
            <<",half_dt,"<<half.effectiveDt()<<",full_dt,"<<full.effectiveDt()<<",refine_donors,"<<remap.report().refinedDonors<<",coarsen_targets,"<<back.report().coarsenedTargets
            <<",advection_relative,"<<moved.endpoint().report().relativeChargeResidual<<",Ohmic_equation,independently_verified,Maxwell,NOT_CONNECTED\n";
    }
    check(solves==8&&remaps==8&&guards==14,"E2 minimum actual consumer paths incomplete");
    std::cout<<"ELECTRIC_VOF_MOVING_PASS,solves,"<<solves<<",AMR_receipts,"<<remaps<<",guards,"<<guards<<",physical_gate,1e-10,closed_only,no_Maxwell,no_App_commit\n";
}
} // namespace ElectricVofMovingContract


namespace ElectricVofGradientContract {
using namespace ElectricVofContract;
void run() {
    std::size_t solves=0,gradientCells=0,constantCells=0;
    long double peakGradient=0,peakH=0;
    for(bool rz:{false,true}){
        auto mesh=meshFor();const auto domain=domainFor(*mesh,rz);
        check(!domain.uniformLeaves(),"E3a gradient witness lacks actual 2:1 AMR");
        const auto stocks=phaseStocks(*mesh,domain);
        const Clock acceptedClock{0,.2,.2,Clock::Position::Initial};
        const Clock geometryClock{1,.2,.22,Clock::Position::Midpoint};
        const Clock candidate{1,.2,.22,rz?Clock::Position::Midpoint:Clock::Position::End};
        const auto accepted=Publication::prepareLocalPlicInventory(*mesh,domain,acceptedClock,stocks,{});
        const auto geometry=Publication::prepareLocalPlicInventory(*mesh,domain,geometryClock,stocks,{});
        const auto saved=snapshot(geometry);const ElectricVofMaterials material{2.,5.,.7,.03};
        for(bool constant:{false,true}){
            std::vector<Boundary> boundaries;
            auto reference=references(*mesh,geometry,material,candidate.time()-acceptedClock.time(),boundaries);
            std::vector<double> q(mesh->cellIndexNum());
            if(constant){
                // A genuine zero-potential patch has exact zero physical
                // charge/current scales: no activity floor or relaxed gate.
                for(std::size_t row=0;row<reference.size();++row){reference[row].affineH=0;
                    if(boundaries[row].type==Type::Dirichlet)boundaries[row].value=0;}
            }else{
                std::vector<long double> integrated(q.size());
                for(const auto& face:reference)scatter(integrated,face.face,static_cast<long double>(face.a)*face.affineH);
                for(std::size_t i=0;i<q.size();++i)q[i]=static_cast<double>(integrated[i]);
            }
            const auto charge=ElectricVofCharge::fromIntegrated(*mesh,accepted,q);
            const auto field=ElectricVofField::prepareOhmic(*mesh,geometry,charge,candidate,material,boundaries);
            // Reuse every original independent physical Gauss, current,
            // integrated-Q, boundary and fullclock assertion without changes.
            const auto metric=verify(*mesh,geometry,charge,field,candidate,material,reference,!constant,false);
            const auto& gradients=field.cellPotentialGradients();
            check(gradients.size()==q.size(),"E3a actual cell LS gradient native size");
            const Vector2D exact=constant?Vector2D{0,0}:(rz?Vector2D{0,1}:Vector2D{1,.37});
            long double localPeak=0;
            for(std::size_t i=0;i<gradients.size();++i){const auto g=gradients[i];
                check(std::isfinite(g.x)&&std::isfinite(g.y),"E3a nonfinite full cell LS gradient");
                const long double error=std::hypot(static_cast<long double>(g.x)-exact.x,static_cast<long double>(g.y)-exact.y);
                const long double scale=constant?0:std::hypot(static_cast<long double>(exact.x),static_cast<long double>(exact.y));
                bounded(error,scale,gate,"E3a independent affine/constant complete cell LS gradient");
                localPeak=std::max(localPeak,error);++gradientCells;
                if(constant){check(field.potential()[i]==0&&g.x==0&&g.y==0,"E3a zero patch acquired a gradient");++constantCells;}
            }
            if(constant)for(const auto& face:field.faces())
                check(face.integratedField==0&&face.displacement==0&&face.current==0,"E3a exact zero patch acquired H/D/J");
            const auto retained=field;
            check(&retained.cellPotentialGradients()==&field.cellPotentialGradients(),"E3a copied field lost original gradient owner");
            check(charge.integratedCharge()==q&&snapshot(geometry)==saved,"E3a passive gradient changed source charge/stocks/plane");
            field.validateFor(*mesh,geometry,charge,candidate);
            peakGradient=std::max(peakGradient,localPeak);peakH=std::max(peakH,metric.affineField);++solves;
            std::cout<<"ELECTRIC_VOF_CELL_GRADIENT,"<<(rz?"RZ":"XY")<<",constant_phi,"<<constant
                <<",cells,"<<gradients.size()<<",nonorthogonal,"<<field.report().nonorthogonalFaces<<",axis,"<<field.report().axisFaces
                <<",gradient_error,"<<localPeak<<",independent_H,"<<metric.affineField<<",Gauss,"<<metric.gauss
                <<",charge_equation,"<<metric.charge<<",boundary_total,"<<metric.global<<",Maxwell,NOT_CONNECTED\n";
        }
    }
    check(solves==4&&constantCells>0&&gradientCells==2*constantCells,"E3a minimum actual gradient consumers incomplete");
    std::cout<<"ELECTRIC_VOF_CELL_GRADIENT_PASS,solves,"<<solves<<",gradient_cells,"<<gradientCells
        <<",peak_gradient_error,"<<peakGradient<<",peak_H_error,"<<peakH<<",Gauss_charge_gate,1e-10,no_Maxwell_no_Hydro_load\n";
}
} // namespace ElectricVofGradientContract


namespace ElectricVofMaxwellContract {
using namespace ElectricVofContract;
using Stress=ElectricVofMaxwellStress;
std::unique_ptr<QuadtreeMesh> jumpMesh() {
    DictionaryMesh cfg;cfg.Nx=cfg.Ny=8;cfg.minCoord={0,0};cfg.maxCoord={1,1};cfg.maxInitialIterations=0;cfg.maxCellNum=4096;
    auto mesh=std::make_unique<QuadtreeMesh>(cfg);DictionaryMesh::RefineBox box;
    box.min={.75,0};box.max={1,1};box.level=1;mesh->applyRefineBoxes({box});return mesh;
}
double jumpPotential(Vector2D p) {
    const double displacement=1/(.5/2.+.5/5.);
    return 1-displacement*(p.x<=.5?p.x/2.:.5/2.+(p.x-.5)/5.);
}
void run() {
    std::size_t solves=0,samples=0,guards=0,jumpFaces=0,nonzeroRadial=0;
    long double peakField=0,peakTensor=0,peakSide=0,peakGauss=0,peakSquaredRadius=0;
    for(int kind=0;kind<4;++kind){const bool rz=kind==1||kind==3,jump=kind==2,charged=kind==3;
        auto mesh=jump?jumpMesh():meshFor();const auto domain=domainFor(*mesh,rz);check(!domain.uniformLeaves(),"E3B2a fixture lacks actual 2:1 AMR");
        const auto stocks=phaseStocks(*mesh,domain);const Clock oldClock{0,.2,.2,Clock::Position::Initial},clock{1,.2,.22,Clock::Position::Midpoint};
        const auto accepted=Publication::prepareLocalPlicInventory(*mesh,domain,oldClock,stocks,{});
        const auto geometry=Publication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,{});const auto savedGeometry=snapshot(geometry);
        const ElectricVofMaterials material=jump?ElectricVofMaterials{2,5,0,0}:ElectricVofMaterials{3,3,0,0};
        std::vector<Boundary> boundaries;const auto reference=references(*mesh,geometry,material,clock.time()-oldClock.time(),boundaries);
        std::vector<double> q(mesh->cellIndexNum());
        if(jump){for(std::size_t row=0;row<reference.size();++row)if(boundaries[row].type==Type::Dirichlet)
                boundaries[row].value=jumpPotential(reference[row].face.center);
        }else if(charged){for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto p=mesh->cells(mesh->cellIndexToLabel[ci]).centerCoord();
                q[ci]=domain.cellVolumes()[domain.cellRowByIndex()[ci]]*(p.x-.5+.19*(p.y-.5));}
        }else{std::vector<long double> integrated(q.size());for(const auto& f:reference)scatter(integrated,f.face,static_cast<long double>(f.a)*f.affineH);
            for(std::size_t ci=0;ci<q.size();++ci)q[ci]=static_cast<double>(integrated[ci]);}
        const auto charge=ElectricVofCharge::fromIntegrated(*mesh,accepted,q);
        const auto field=ElectricVofField::prepareOhmic(*mesh,geometry,charge,clock,material,boundaries);++solves;
        const auto metric=verify(*mesh,geometry,charge,field,clock,material,reference,!jump&&!charged,true);
        const auto savedPhi=field.potential();std::vector<std::array<double,3>> savedFaces;
        for(const auto& f:field.faces())savedFaces.push_back({f.integratedField,f.displacement,f.current});
        const auto stress=Stress::prepare(*mesh,field);stress.validateFor(*mesh,field,clock);
        check(sameClock(stress.timeLayer(),clock)&&stress.axisymmetric()==rz&&stress.cells().size()==q.size()&&
            &stress.field().potential()==&field.potential()&&stress.report().wholeCellSidesOnly,"E3B2a actual field/fullclock/RT0 scope lost");
        if(rz)check(stress.report().axisCells>0,"E3B2a RZ has no true-axis cells");
        long double localField=0,localTensor=0,localSide=0,localGauss=0,localRadius=0;
        for(const auto& c:stress.cells()){
            const auto center=(c.lower+c.upper)*.5;
            const std::array<Vector2D,5> points{center,c.lower,Vector2D{c.lower.x,c.upper.y},Vector2D{c.upper.x,c.lower.y},c.upper};
            for(const auto point:points){const auto s=stress.sample(c.nativeCell,point);++samples;
                check(s.nativeCell==c.nativeCell&&s.point.x==point.x&&s.point.y==point.y,"E3B2a actual sample cell/point changed");
                bounded(s.electricField.x-static_cast<long double>(s.displacement.x)/c.epsilon,
                    std::abs(static_cast<long double>(s.displacement.x)/c.epsilon),1e-12L,"E3B2a radial/x E not actual D/epsilon_i");
                bounded(s.electricField.y-static_cast<long double>(s.displacement.y)/c.epsilon,
                    std::abs(static_cast<long double>(s.displacement.y)/c.epsilon),1e-12L,"E3B2a axial/y E not actual D/epsilon_i");
                const long double ex=s.electricField.x,ey=s.electricField.y,energy=.5L*c.epsilon*(ex*ex+ey*ey);
                const std::array<long double,4> expected{c.epsilon*ex*ex-energy,c.epsilon*ey*ey-energy,c.epsilon*ex*ey,-energy};
                const std::array<double,4> actual{s.tensor.xx,s.tensor.yy,s.tensor.xy,s.tensor.hoop};
                for(int k=0;k<4;++k){const auto error=std::abs(actual[k]-expected[k]);bounded(error,energy,1e-12L,"E3B2a represented Maxwell tensor contraction/hoop");
                    if(energy>0)localTensor=std::max(localTensor,error/energy);}
                if(!charged){const long double dx=jump?1/(.5L/2+.5L/5):0;
                    const long double exactX=jump?dx/c.epsilon:(rz?0:-1),exactY=jump?0:(rz?-1:-.37L);
                    const long double norm=std::hypot(exactX,exactY),error=std::hypot(ex-exactX,ey-exactY);
                    bounded(error,norm,gate,"E3B2a independent affine/serial-resistance electric field");localField=std::max(localField,error/norm);
                    const long double e=.5L*c.epsilon*(exactX*exactX+exactY*exactY);
                    const std::array<long double,4> physical{c.epsilon*exactX*exactX-e,c.epsilon*exactY*exactY-e,c.epsilon*exactX*exactY,-e};
                    for(int k=0;k<4;++k)bounded(actual[k]-physical[k],e,gate,"E3B2a independent physical Maxwell tensor");
                }
                if(rz){if(point.x==0)check(s.displacement.x==0&&s.electricField.x==0&&s.tensor.xy==0,"E3B2a axis analytic limit/fake radius");
                    if(charged&&s.displacement.x!=0)++nonzeroRadial;
                    const auto left=stress.sample(c.nativeCell,{c.lower.x,point.y}),right=stress.sample(c.nativeCell,{c.upper.x,point.y});
                    const long double r=point.x,r0=c.lower.x,r1=c.upper.x,t=(r*r-r0*r0)/(r1*r1-r0*r0);
                    const long double expectedFlux=(1-t)*r0*left.displacement.x+t*r1*right.displacement.x;
                    const long double scale=std::abs((1-t)*r0*left.displacement.x)+std::abs(t*r1*right.displacement.x);
                    const long double error=r*s.displacement.x-expectedFlux;bounded(error,scale,1e-12L,"E3B2a actual RZ rDr not affine in squared radius");
                    if(scale>0)localRadius=std::max(localRadius,std::abs(error)/scale);}
            }
            const std::array<Vector2D,4> sidePoints{Vector2D{c.lower.x,center.y},Vector2D{c.upper.x,center.y},Vector2D{center.x,c.lower.y},Vector2D{center.x,c.upper.y}};
            std::array<long double,4> density{};long double sideScale=0;
            for(int side=0;side<4;++side){const auto s=stress.sample(c.nativeCell,sidePoints[side]);density[side]=side<2?s.displacement.x:s.displacement.y;
                const auto represented=density[side]*c.sidePhysicalMeasure[side],error=represented-c.sideDisplacement[side];
                const auto scale=std::abs(represented)+std::abs(c.sideDisplacement[side]);bounded(error,scale,1e-12L,"E3B2a reconstructed whole-side integral differs from actual D");
                if(scale>0)localSide=std::max(localSide,std::abs(error)/scale);sideScale+=std::abs(c.sideDisplacement[side]);}
            const long double dx=static_cast<long double>(c.upper.x)-c.lower.x,dz=static_cast<long double>(c.upper.y)-c.lower.y;
            const long double radial=rz?2*(c.upper.x*density[1]-c.lower.x*density[0])/(dx*(static_cast<long double>(c.upper.x)+c.lower.x)):(density[1]-density[0])/dx;
            const long double divergence=radial+(density[3]-density[2])/dz;
            const long double volume=dx*dz*(rz?CartesianMeasure::twoPi*(static_cast<long double>(c.upper.x)+c.lower.x)/2:1.L);
            const long double error=volume*divergence-field.chargeCandidate().integratedCharge()[c.nativeCell];
            const long double scale=sideScale+std::abs(field.chargeCandidate().integratedCharge()[c.nativeCell]);
            bounded(error,scale,gate,"E3B2a independent RT0 physical divergence differs from actual Gauss charge");
            if(scale>0)localGauss=std::max(localGauss,std::abs(error)/scale);
        }
        if(jump){std::size_t samePhaseHanging=0;const long double d=1/(.5L/2+.5L/5),expectedJump=.5L*d*d*(1/2.L-1/5.L);
            for(const auto& f:field.faces()){if(f.negativeCell<0||f.positiveCell<0)continue;const auto native=domain.face(*mesh,f.nativeRow);
                const auto a=mesh->cells(mesh->cellIndexToLabel[f.negativeCell]).size(),b=mesh->cells(mesh->cellIndexToLabel[f.positiveCell]).size();
                if(a.x!=b.x||a.y!=b.y){check(field.epsilonByCell()[f.negativeCell]==5&&field.epsilonByCell()[f.positiveCell]==5&&native.center.x>.5,
                        "E3B2a hanging transition must be same-phase away from epsilon jump");++samePhaseHanging;}
                if(!native.horizontal&&native.center.x==.5){const auto l=stress.sample(f.negativeCell,native.center),r=stress.sample(f.positiveCell,native.center);
                    bounded((l.tensor.xx-r.tensor.xx)-expectedJump,std::abs(expectedJump),gate,"E3B2a independent one-dimensional Maxwell normal pressure jump");
                    for(const auto& s:{l,r})bounded(s.electricField.x-static_cast<long double>(f.displacement)/(f.physicalMeasure*field.epsilonByCell()[s.nativeCell]),
                        std::abs(d/field.epsilonByCell()[s.nativeCell]),gate,"E3B2a planar jump one-sided actual D/epsilon");++jumpFaces;}}
            check(samePhaseHanging>0&&jumpFaces>0,"E3B2a missing genuine remote AMR transition or epsilon jump");}
        if(kind==0){const auto copied=field;stress.validateFor(*mesh,copied,clock);
            auto wrong=clock;++wrong.step;rejected([&]{stress.validateFor(*mesh,field,wrong);},"E3B2a same-time wrong fullclock accepted");++guards;
            auto otherMesh=meshFor(4);rejected([&]{stress.validateFor(*otherMesh,field,clock);},"E3B2a different mesh accepted");++guards;
            const auto fresh=ElectricVofField::prepareOhmic(*mesh,geometry,charge,clock,material,boundaries);++solves;
            rejected([&]{stress.validateFor(*mesh,fresh,clock);},"E3B2a fresh same-value field owner accepted");++guards;
            auto end=clock;end.position=Clock::Position::End;
            const auto noncoeval=ElectricVofField::prepareOhmic(*mesh,geometry,charge,end,material,boundaries);++solves;
            rejected([&]{Stress::prepare(*mesh,noncoeval);},"E3B2a frozen Midpoint to charge End relabeled as coeval stress");++guards;
            rejected([&]{stress.sample(-1,{0,0});},"E3B2a negative native cell accepted");++guards;
            rejected([&]{stress.sample(mesh->cellIndexNum(),{0,0});},"E3B2a foreign native cell accepted");++guards;
            const auto& c=stress.cells().front();
            rejected([&]{stress.sample(c.nativeCell,{std::nextafter(c.lower.x,-std::numeric_limits<double>::infinity()),c.lower.y});},"E3B2a sample outside actual cell accepted");++guards;
            rejected([&]{stress.sample(c.nativeCell,{std::numeric_limits<double>::quiet_NaN(),c.lower.y});},"E3B2a nonfinite sample accepted");++guards;
        }
        check(savedPhi==field.potential()&&savedGeometry==snapshot(geometry)&&charge.integratedCharge()==q,"E3B2a readonly stress changed phi/Q/VOF");
        for(std::size_t row=0;row<field.faces().size();++row)check(savedFaces[row]==std::array<double,3>{field.faces()[row].integratedField,field.faces()[row].displacement,field.faces()[row].current},"E3B2a readonly stress changed H/D/J");
        peakField=std::max(peakField,localField);peakTensor=std::max(peakTensor,localTensor);peakSide=std::max(peakSide,localSide);
        peakGauss=std::max(peakGauss,localGauss);peakSquaredRadius=std::max(peakSquaredRadius,localRadius);
        std::cout<<"ELECTRIC_VOF_MAXWELL_RT0,"<<(rz?"RZ":"XY")<<",kind,"<<kind<<",cells,"<<mesh->cellIndexNum()<<",axis_cells,"<<stress.report().axisCells
            <<",independent_E,"<<localField<<",tensor,"<<localTensor<<",whole_side,"<<localSide<<",RT0_Gauss,"<<localGauss<<",rD_r2,"<<localRadius
            <<",original_Gauss,"<<metric.gauss<<",whole_sides_only,1,no_Hydro_load\n";
    }
    check(solves==6&&guards==8&&jumpFaces>0&&nonzeroRadial>0,"E3B2a minimum actual RT0/axis/jump/owner witnesses incomplete");
    std::cout<<"ELECTRIC_VOF_MAXWELL_RT0_PASS,solves,"<<solves<<",samples,"<<samples<<",guards,"<<guards<<",jump_faces,"<<jumpFaces
        <<",nonzero_radial,"<<nonzeroRadial<<",E,"<<peakField<<",tensor,"<<peakTensor<<",whole_side,"<<peakSide<<",Gauss,"<<peakGauss
        <<",rD_r2,"<<peakSquaredRadius<<",whole_sides_only,no_Hydro_force_no_curlfree_or_sharp_precision_claim\n";
}
} // namespace ElectricVofMaxwellContract

namespace ElectricVofSnapshotContract {
using namespace ElectricVofContract;
using namespace ElectricVofMovingContract;
void verifySnapshot(const QuadtreeMesh& mesh,const Publication& p,const ElectricVofCharge& charge,
    const ElectricVofField& field,const ElectricVofMaterials& material,Backend backend=Backend::Gmres) {
    const auto clock=p.timeLayer();const auto& q=charge.integratedCharge();
    check(field.isElectrostaticSnapshot()&&field.report().electrostaticSnapshot&&!field.report().usedAdvection&&!field.advection(),
        "E3D0 snapshot is not an explicit passive accepted-state mode");
    check(field.effectiveDt()==0&&sameClock(clock,field.geometryTime())&&sameClock(clock,field.chargeSourceTime())&&
        sameClock(clock,field.candidateTime())&&sameClock(clock,field.chargeCandidate().timeLayer()),"E3D0 snapshot invented a clock/interval");
    check(field.geometry().identity()==p.identity()&&field.sourceCharge().sameInstance(charge)&&field.chargeCandidate().sameInstance(charge),
        "E3D0 snapshot replaced an accepted charge/geometry owner");
    check(field.implicitInputCharge().size()==q.size()&&field.potential().size()==q.size(),"E3D0 actual snapshot native sizes");
    check(std::memcmp(q.data(),field.implicitInputCharge().data(),q.size()*sizeof(double))==0&&
        std::memcmp(q.data(),field.chargeCandidate().integratedCharge().data(),q.size()*sizeof(double))==0,
        "E3D0 snapshot altered signed integrated charge bits");
    field.validateFor(mesh,p,charge,clock);
    ElectricVofOhmicQuantumContract::verify(mesh,field,false,true,false,false);
    ElectricVofPostprocessContract::executionReceipt(field,0);
    std::vector<Boundary> bc;const auto ref=references(mesh,p,material,0,bc);
    check(field.faces().size()==ref.size(),"E3D0 actual snapshot native face count");
    std::vector<long double> d(q.size()),absD(q.size()),j(q.size());long double boundaryJ=0,boundaryAbs=0,sum=0,l1=0;
    for(std::size_t row=0;row<ref.size();++row){const auto& f=field.faces()[row];const auto& r=ref[row];
        check(f.nativeRow==static_cast<int>(row)&&f.negativeCell==r.face.negativeCell&&f.positiveCell==r.face.positiveCell,
            "E3D0 changed actual shared native incidence");
        bounded(f.physicalMeasure-r.area,r.area,1e-12L,"E3D0 actual RZ/XY physical native measure");
        bounded(f.epsilon-r.epsilon,r.epsilon,1e-12L,"E3D0 independent native epsilon resistance");
        bounded(f.conductivity-r.sigma,r.sigma,1e-12L,"E3D0 independent native conductivity resistance");
        check(f.implicitCoefficient==f.epsilon,"E3D0 snapshot matrix is not exactly epsilon");
        long double h=f.gradientConstant,hScale=std::abs(f.gradientConstant);
        for(const auto& term:f.gradientCoefficients){check(term.nativeCell>=0&&static_cast<std::size_t>(term.nativeCell)<q.size(),"E3D0 foreign actual H stencil cell");
            const long double v=static_cast<long double>(term.value)*field.potential()[term.nativeCell];h+=v;hScale+=std::abs(v);}
        bounded(f.integratedField-h,hScale,1e-12L,"E3D0 snapshot did not retain actual geometric H");
        bounded(f.displacement-static_cast<long double>(f.epsilon)*f.integratedField,std::abs(f.displacement),1e-12L,"E3D0 snapshot D differs from epsilon H");
        bounded(f.current-static_cast<long double>(f.conductivity)*f.integratedField,std::abs(f.current),1e-12L,"E3D0 snapshot diagnostic J differs from sigma H");
        scatter(d,r.face,f.displacement);scatter(j,r.face,f.current);
        if(f.negativeCell>=0)absD[f.negativeCell]+=std::abs(f.displacement);
        if(f.positiveCell>=0)absD[f.positiveCell]+=std::abs(f.displacement);
        if(f.negativeCell<0||f.positiveCell<0){const int sign=f.negativeCell>=0?1:-1;boundaryJ+=sign*static_cast<long double>(f.current);boundaryAbs+=std::abs(f.current);}
        if(r.area==0)check(p.domain().axisymmetric()&&r.face.center.x==0&&f.integratedField==0&&f.displacement==0&&f.current==0,
            "E3D0 true axis invented a finite radial flux");
    }
    long double peakGauss=0,netCurrentL1=0;
    for(std::size_t i=0;i<q.size();++i){const auto scale=absD[i]+std::abs(q[i]);bounded(d[i]-q[i],scale,gate,"E3D0 independent accepted signed-charge Gauss");
        if(scale>0)peakGauss=std::max(peakGauss,std::abs(d[i]-q[i])/scale);sum+=q[i];l1+=std::abs(q[i]);netCurrentL1+=std::abs(j[i]);}
    backendReceipt(field,backend);
    check(field.report().relativeLinearResidual<=1e-11&&
        field.report().maximumRelativeGaussResidual<=gate,"E3D0 original linear/Gauss gates or passive-Q semantics failed");
    // Positive conductivity and nonzero divergence of instantaneous J are
    // allowed: this is electrostatics, never a hidden conduction advance.
    check(netCurrentL1>0,"E3D0 positive-sigma snapshot lacks diagnostic current witness");
    std::cout<<"ELECTRIC_VOF_SNAPSHOT,"<<(p.domain().axisymmetric()?"RZ":"XY")<<",step,"<<clock.step<<",time,"<<clock.time()
        <<",cells,"<<q.size()<<",Gauss,"<<peakGauss<<",diagnostic_netJ_L1,"<<netCurrentL1<<",dt,0,same_charge_owner,1\n";
}
void run() {
    std::size_t snapshots=0,remaps=0,guards=0,oldOhmic=0;
    for(bool rz:{false,true}){
        auto mesh=meshFor();const auto domain=domainFor(*mesh,rz);check(!domain.uniformLeaves(),"E3D0 snapshot requires actual 2:1 AMR");
        const auto stocks=phaseStocks(*mesh,domain);const Clock start{0,.2,.2,Clock::Position::Initial},end{1,.2,.202,Clock::Position::End};
        const auto p=Publication::prepareLocalPlicInventory(*mesh,domain,start,stocks,{});const auto saved=snapshot(p);
        const ElectricVofMaterials material{2.,5.,.7,.03},zeroSigma{2.,5.,0.,0.};
        std::vector<double> q(mesh->cellIndexNum());for(int ci=0;ci<mesh->cellIndexNum();++ci)
            q[ci]=static_cast<double>(domain.cellVolumes()[domain.cellRowByIndex()[ci]]*.013L*(static_cast<int>(ci%7)-3));
        q[0]=-0.;const auto charge=ElectricVofCharge::fromIntegrated(*mesh,p,q);std::vector<Boundary> bc;references(*mesh,p,material,0,bc);
        const auto initial=ElectricVofField::prepareElectrostaticSnapshot(*mesh,p,charge,material,bc);verifySnapshot(*mesh,p,charge,initial,material);++snapshots;
        const auto initialStress=LAMP::Electric::ElectricVofMaxwellStress::prepare(*mesh,initial);initialStress.validateFor(*mesh,initial,start);
        const auto fresh=Publication::prepareLocalPlicInventory(*mesh,domain,start,stocks,{});
        rejected([&]{ElectricVofField::prepareElectrostaticSnapshot(*mesh,fresh,charge,material,bc);},"E3D0 foreign same-clock publication self-signed snapshot");++guards;
        const auto foreignCharge=ElectricVofCharge::fromIntegrated(*mesh,p,q);
        rejected([&]{initial.validateFor(*mesh,p,foreignCharge,start);},"E3D0 fresh same-value source replaced actual snapshot owner");++guards;
        const Clock otherClock{7,.19,.2,Clock::Position::End};
        const auto other=Publication::prepareLocalPlicInventory(*mesh,domain,otherClock,stocks,{});
        rejected([&]{ElectricVofField::prepareElectrostaticSnapshot(*mesh,other,charge,material,bc);},"E3D0 same-time different complete accepted clock allowed");++guards;
        rejected([&]{initial.validateFor(*mesh,p,charge,otherClock);},"E3D0 snapshot validation ignored complete clock");++guards;
        rejected([&]{ElectricVofField::prepareOhmic(*mesh,p,charge,start,material,bc);},"E3D0 old Ohmic accepted a fake zero-duration snapshot");++guards;
        const std::vector<double> zero(domain.nativeFaceCount());const auto slab=Slab::prepare(*mesh,p,end,zero,zero);
        const auto flow=VofTransport::prepareLocalPlicSweep(*mesh,p,slab,1000,1.226);
        const auto moved=ChargeTransport::prepareClosed(*mesh,charge,flow);
        check(std::memcmp(q.data(),moved.endpoint().integratedCharge().data(),q.size()*sizeof(double))==0,"E3D0 real zero-advection lost signed Q");
        const auto reference=references(*mesh,flow.endpoint(),zeroSigma,end.time()-start.time(),bc);
        const auto endpoint=ElectricVofField::prepareOhmic(*mesh,flow.endpoint(),moved.endpoint(),zeroSigma,bc);
        verifyConduction(*mesh,moved.endpoint(),endpoint,true);
        verify(*mesh,flow.endpoint(),charge,endpoint,end,zeroSigma,reference,false,true);++oldOhmic;
        check(!endpoint.isElectrostaticSnapshot()&&!endpoint.report().electrostaticSnapshot,"E3D0 old positive-dt Ohmic mode changed");
        if(!rz){references(*mesh,flow.midpoint(),material,flow.midpoint().timeLayer().time()-start.time(),bc);
            const auto mid=ElectricVofField::prepareOhmic(*mesh,flow.midpoint(),moved.midpoint(),material,bc);verifyConduction(*mesh,moved.midpoint(),mid,false);++oldOhmic;
            rejected([&]{ElectricVofField::prepareElectrostaticSnapshot(*mesh,flow.midpoint(),mid.chargeCandidate(),material,bc);},"E3D0 Midpoint was relabeled as an accepted snapshot");++guards;}
        auto fine=mesh->clone();DictionaryMesh::RefineBox box;box.min={0,0};box.max={1,1};box.level=1;fine->applyRefineBoxes({box});
        check(fine->cellIndexNum()>mesh->cellIndexNum(),"E3D0 same-End snapshot lacks actual changed tree");
        const auto transfer=VofTransfer::prepare(*mesh,flow.endpoint(),*fine,domainFor(*fine,rz),end);
        const auto remap=ChargeRemap::prepare(*mesh,flow.endpoint(),endpoint.chargeCandidate(),*fine,transfer);++remaps;
        check(remap.report().refinedDonors>0&&remap.report().unchangedTargets>0&&sameClock(remap.chargeCandidate().timeLayer(),end),"E3D0 actual E2 same-End remap receipt absent");
        long double sourceTotal=0,targetTotal=0,l1=0;for(double v:endpoint.chargeCandidate().integratedCharge()){sourceTotal+=v;l1+=std::abs(v);}
        for(double v:remap.chargeCandidate().integratedCharge()){targetTotal+=v;l1+=std::abs(v);}
        bounded(targetTotal-sourceTotal,l1,gate,"E3D0 independent signed global actual AMR charge ledger");
        references(*fine,transfer.candidate(),material,0,bc);
        const auto target=ElectricVofField::prepareElectrostaticSnapshot(*fine,transfer.candidate(),remap.chargeCandidate(),material,bc);
        verifySnapshot(*fine,transfer.candidate(),remap.chargeCandidate(),target,material);++snapshots;
        const auto stress=LAMP::Electric::ElectricVofMaxwellStress::prepare(*fine,target);stress.validateFor(*fine,target,end);
        rejected([&]{initial.validateFor(*fine,transfer.candidate(),remap.chargeCandidate(),end);},"E3D0 old-layout snapshot survived actual remap");++guards;
        rejected([&]{ElectricVofField::prepareOhmic(*fine,transfer.candidate(),remap.chargeCandidate(),end,material,bc);},"E3D0 same-End refresh weakened old Ohmic next-step guard");++guards;
        check(snapshot(p)==saved&&std::memcmp(q.data(),charge.integratedCharge().data(),q.size()*sizeof(double))==0,
            "E3D0 snapshot success/rejection/AMR mutated accepted VOF/charge source");
    }
    check(snapshots==4&&remaps==2&&guards==15&&oldOhmic==3,"E3D0 minimum actual snapshot/AMR/old-mode witnesses incomplete");
    std::cout<<"ELECTRIC_VOF_SNAPSHOT_PASS,snapshots,"<<snapshots<<",actual_E2_remaps,"<<remaps<<",old_positive_dt_Ohmic,"<<oldOhmic
        <<",guards,"<<guards<<",Gauss_gate,1e-10,no_Q_advance,no_App_commit\n";
}
} // namespace ElectricVofSnapshotContract

namespace ElectricVofLayerRelaxationContract {
using namespace ElectricVofContract;
using namespace ElectricVofMovingContract;
struct Layers {
    long double liquidField=0,gasField=0,charge=0,leftCurrent=0,rightCurrent=0;
};
Layers observe(const QuadtreeMesh& mesh,const ElectricVofField& field) {
    const auto& domain=field.geometry().domain();Layers r;long double liquidH=0,gasH=0,liquidArea=0,gasArea=0;
    // Independent actual plane integrals; no gradient fitted from this test's
    // Q update, matrix coefficients, or continuum answer.
    for(const auto& f:field.faces()){const auto n=domain.face(mesh,f.nativeRow);
        if(n.horizontal){if(n.negativeCell<0||n.positiveCell<0)check(f.integratedField==0&&f.current==0,"E3D layer insulating top/bottom leak");continue;}
        if(n.center.x==.25){liquidH+=f.integratedField;liquidArea+=f.physicalMeasure;}
        if(n.center.x==.875){gasH+=f.integratedField;gasArea+=f.physicalMeasure;}
        if(n.center.x==0)r.leftCurrent+=f.current;
        if(n.center.x==1)r.rightCurrent+=f.current;
    }
    check(liquidArea==1&&gasArea==1,"E3D layer actual section does not cover unit physical area");
    r.liquidField=liquidH/liquidArea;r.gasField=gasH/gasArea;
    for(double q:field.chargeCandidate().integratedCharge())r.charge+=q;
    r.charge/=liquidArea;return r;
}
void run() {
    const long double e1=2,e2=5,k1=1,k2=.25L,l1=.5L,l2=.5L,V=1;
    const long double A=e1*l2+e2*l1,B=k1*l2+k2*l1,tau=A/B;
    const long double el0=e2*V/A,eg0=e1*V/A,elInf=k2*V/B,egInf=k1*V/B,qInf=V*(e2*k1-e1*k2)/B;
    check(e1*k2!=e2*k1&&tau>0&&qInf>0,"E3D layer must exercise nonproportional epsilon/conductivity");
    const ElectricVofMaterials material{static_cast<double>(e1),static_cast<double>(e2),static_cast<double>(k1),static_cast<double>(k2)};
    DictionaryMesh cfg;cfg.Nx=64;cfg.Ny=4;cfg.minCoord={0,0};cfg.maxCoord={1,1};cfg.maxInitialIterations=0;cfg.maxCellNum=4096;
    auto mesh=std::make_unique<QuadtreeMesh>(cfg);DictionaryMesh::RefineBox box;box.min={.75,0};box.max={1,1};box.level=1;mesh->applyRefineBoxes({box});
    const auto domain=domainFor(*mesh,false);check(mesh->cellIndexNum()==448&&!domain.uniformLeaves(),"E3D layer true remote 2:1 fixture changed");
    const auto stocks=phaseStocks(*mesh,domain);const Clock start{0,0,0,Clock::Position::Initial};
    auto geometry=Publication::prepareLocalPlicInventory(*mesh,domain,start,stocks,{});
    std::vector<double> q(mesh->cellIndexNum());auto charge=ElectricVofCharge::fromIntegrated(*mesh,geometry,q);
    std::vector<Boundary> boundary(domain.nativeFaceCount());std::size_t hanging=0;
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(*mesh,row);
        if(f.negativeCell>=0&&f.positiveCell>=0){boundary[row]={Type::Interior,0};
            const auto& a=mesh->cells(mesh->cellIndexToLabel[f.negativeCell]);const auto& b=mesh->cells(mesh->cellIndexToLabel[f.positiveCell]);
            if(a.level()!=b.level()){check(f.center.x==.75&&stocks[f.negativeCell][0]==0&&stocks[f.positiveCell][0]==0,
                    "E3D layer epsilon jump contaminated by hanging transition");++hanging;}
        }else if(f.horizontal)boundary[row]={Type::HomogeneousNeumann,0};
        else{check(f.center.x==0||f.center.x==1,"E3D layer foreign electric wall");boundary[row]={Type::Dirichlet,f.center.x==0?1.:0.};}
    }
    check(hanging>0,"E3D layer no actual same-phase hanging faces");
    const auto initial=ElectricVofField::prepareElectrostaticSnapshot(*mesh,geometry,charge,material,boundary);
    const auto initialObserved=observe(*mesh,initial);
    check(initial.isElectrostaticSnapshot()&&initial.chargeCandidate().sameInstance(charge)&&initial.effectiveDt()==0,
        "E3D layer initial state was not genuine passive electrostatics");
    bounded(initialObserved.liquidField-el0,std::abs(el0),1e-10L,"E3D layer independent initial liquid series resistance");
    bounded(initialObserved.gasField-eg0,std::abs(eg0),1e-10L,"E3D layer independent initial gas series resistance");
    check(initialObserved.charge==0,"E3D layer initial charge must be genuinely zero");
    std::cout<<std::setprecision(17)<<"ELECTRIC_VOF_LAYER_RELAXATION_BEGIN,cells,"<<mesh->cellIndexNum()<<",hanging,"<<hanging
        <<",epsilon_L,"<<e1<<",epsilon_G,"<<e2<<",conductivity_L,"<<k1<<",conductivity_G,"<<k2<<",tau,"<<tau
        <<",E_L_0,"<<el0<<",E_G_0,"<<eg0<<",E_L_inf,"<<elInf<<",E_G_inf,"<<egInf<<",q_inf,"<<qInf
        <<",space_budget,"<<2.L/64<<",method,first_order_BE_and_P0,no_steady_or_high_order_claim\n";
    const double seed=static_cast<double>(tau/16);const std::vector<double> zero(domain.nativeFaceCount());
    long double thetaBE=1,peakGauss=0,peakPhysical=0;std::size_t checkpoints=0;
    for(int n=1;n<=32;++n){
        const auto oldGeometry=snapshot(geometry);const auto oldQ=charge.integratedCharge();const auto source=charge;
        const double begin=charge.timeLayer().time(),end=static_cast<double>(n)*seed;const Clock target{n,begin,end,Clock::Position::End};
        const auto slab=Slab::prepare(*mesh,geometry,target,zero,zero);
        const auto transport=VofTransport::prepareLocalPlicSweep(*mesh,geometry,slab,1000,1.226);
        for(const auto& f:transport.faces())check(f.full.total==0&&f.half.total==0,"E3D layer actual stationary VOF has nonzero T");
        const auto moved=ChargeTransport::prepareClosed(*mesh,charge,transport);
        check(std::memcmp(oldQ.data(),moved.endpoint().integratedCharge().data(),oldQ.size()*sizeof(double))==0,"E3D layer stationary Qadv changed bits");
        const auto field=ElectricVofField::prepareOhmic(*mesh,transport.endpoint(),moved.endpoint(),material,boundary);
        verifyConduction(*mesh,moved.endpoint(),field,false);
        check(!field.isElectrostaticSnapshot()&&field.sourceCharge().sameInstance(source)&&
            field.geometry().identity()==transport.endpoint().identity()&&sameClock(field.geometryTime(),target)&&sameClock(field.candidateTime(),target),
            "E3D layer Ohmic stage lost actual endpoint/source/fullclock");
        check(field.report().usedGMRES&&!field.report().usedSuperLU&&field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=gate&&field.report().nonorthogonalFaces>0,
            "E3D layer changed retained electric solver/Gauss gates");
        const auto actual=observe(*mesh,field);const long double dt=field.effectiveDt();check(dt==target.time()-source.timeLayer().time(),"E3D layer fake nominal conduction dt");
        thetaBE/=1+dt/tau;long double oldSum=0,l1Q=0;for(double v:oldQ){oldSum+=v;l1Q+=std::abs(v);}
        for(double v:field.chargeCandidate().integratedCharge())l1Q+=std::abs(v);
        bounded(actual.charge-oldSum+dt*(actual.rightCurrent-actual.leftCurrent),
            l1Q+dt*(std::abs(actual.leftCurrent)+std::abs(actual.rightCurrent)),gate,"E3D layer independent global Q plus actual boundary current");
        check(snapshot(geometry)==oldGeometry&&std::memcmp(oldQ.data(),source.integratedCharge().data(),oldQ.size()*sizeof(double))==0,
            "E3D layer successful candidate mutated accepted independent phase/Q source");
        check(snapshot(transport.endpoint())==oldGeometry,"E3D layer zero flow changed independent L/G inventory/PLIC plane");
        peakGauss=std::max(peakGauss,static_cast<long double>(field.report().maximumRelativeGaussResidual));
        if(n==8||n==16||n==32){
            const long double t=target.time(),theta=std::exp(-t/tau),expectedQ=qInf*(1-theta);
            const long double expectedL=elInf+(el0-elInf)*theta,expectedG=egInf+(eg0-egInf)*theta;
            const long double errorQ=std::abs(actual.charge-expectedQ)/std::abs(qInf);
            const long double errorL=std::abs(actual.liquidField-expectedL)/std::abs(el0-elInf);
            const long double errorG=std::abs(actual.gasField-expectedG)/std::abs(eg0-egInf);
            const long double timeBudget=std::abs(thetaBE-theta),budget=2.L/64+timeBudget;
            std::cout<<"ELECTRIC_VOF_LAYER_RELAXATION_SAMPLE,step,"<<n<<",time,"<<t<<",time_over_tau,"<<t/tau
                <<",q,"<<actual.charge<<",q_exact,"<<expectedQ<<",E_L,"<<actual.liquidField<<",E_L_exact,"<<expectedL
                <<",E_G,"<<actual.gasField<<",E_G_exact,"<<expectedG<<",error_q,"<<errorQ<<",error_E_L,"<<errorL
                <<",error_E_G,"<<errorG<<",space_budget,"<<2.L/64<<",BE_time_budget,"<<timeBudget<<",physical_budget,"<<budget<<",tau_fit,";
            if(actual.charge>0&&actual.charge<qInf)std::cout<<-t/std::log(1-actual.charge/qInf);else std::cout<<"NA";
            std::cout<<",left_J,"<<actual.leftCurrent<<",right_J,"<<actual.rightCurrent<<std::endl;
            // This is the declared first-order physical validation budget,
            // distinct from the unchanged original 1e-10 algebra gates.
            check(errorQ<=budget&&errorL<=budget&&errorG<=budget,"E3D layer continuous RC relaxation outside declared spatial/BE physical budget");
            peakPhysical=std::max({peakPhysical,errorQ,errorL,errorG});++checkpoints;
        }
        geometry=transport.endpoint();charge=field.chargeCandidate();
    }
    check(checkpoints==3&&sameClock(geometry.timeLayer(),charge.timeLayer())&&charge.timeLayer().step==32,
        "E3D layer did not complete the declared actual time window");
    std::cout<<"ELECTRIC_VOF_LAYER_RELAXATION_PASS,steps,32,checkpoints,"<<checkpoints<<",accepted_time,"<<charge.timeLayer().time()
        <<",T_over_tau,"<<charge.timeLayer().time()/tau<<",Gauss,"<<peakGauss
        <<",physical_exp_error,"<<peakPhysical<<",independent_charge_gate,1e-10,declared_first_order_budget,no_steady_or_high_order_claim\n";
}
} // namespace ElectricVofLayerRelaxationContract

namespace ElectricVofBoundaryChargeContract {
using namespace ElectricVofContract;
using namespace ElectricVofMovingContract;
using ChargeBoundary=LAMP::Electric::ElectricVofChargeBoundary;
using FluidBoundary=Slab::Boundary;
using Mode=FluidBoundary::Mode;
void checkReceipt(const QuadtreeMesh& mesh,const LAMP::Electric::ElectricVofAdvectedCharge& a,
    const ChargeBoundary& boundary,double inletDensity,bool half) {
    a.validateFor(mesh,a.sourceCharge(),a.transport(),a.timeLayer());
    check(a.boundary()&&a.boundary()->sameInstance(boundary)&&a.report().boundaryEnabled,"F4 receipt lost actual inlet boundary owner");
    const auto& domain=a.transport().source().domain();const auto& old=a.sourceCharge().integratedCharge();const auto& next=a.integratedCharge();
    std::vector<long double> net(old.size()),absNet(old.size());long double source=0,target=0,l1=0,boundaryNet=0,incoming=0,outgoing=0,fluxL1=0;
    std::size_t inRows=0,outRows=0,signedOpposite=0;
    for(const auto& f:a.faces()){const auto actual=domain.face(mesh,f.nativeRow);const auto& original=a.transport().faces()[f.nativeRow];
        const long double T=half?original.half.total:original.full.total;
        check(f.totalVolume==T&&f.negativeCell==actual.negativeCell&&f.positiveCell==actual.positiveCell,"F4 actual F1 shared T/native incidence changed");
        const bool exterior=actual.negativeCell<0||actual.positiveCell<0;const int sign=actual.negativeCell>=0?1:-1;
        const long double outward=sign>0?T:-T;long double expected=0;
        if(exterior&&outward<0){expected=static_cast<long double>(inletDensity)*T;++inRows;incoming-=sign*static_cast<long double>(f.integratedCharge);}
        else if(T!=0){const int donor=exterior?(actual.negativeCell>=0?actual.negativeCell:actual.positiveCell):(T>0?actual.negativeCell:actual.positiveCell);
            const auto capacity=domain.cellVolumes()[domain.cellRowByIndex()[donor]];expected=static_cast<long double>(old[donor])*T/capacity;
            if(exterior){++outRows;outgoing+=sign*static_cast<long double>(f.integratedCharge);}}
        bounded(static_cast<long double>(f.integratedCharge)-expected,std::abs(expected)+std::abs(f.integratedCharge),1e-12L,
            "F4 independent reservoir rho*T or actual inside P0 donor charge");
        if(T!=0&&f.integratedCharge!=0&&(T>0)!=(f.integratedCharge>0))++signedOpposite;
        if(actual.negativeCell>=0){net[actual.negativeCell]+=f.integratedCharge;absNet[actual.negativeCell]+=std::abs(f.integratedCharge);}
        if(actual.positiveCell>=0){net[actual.positiveCell]-=f.integratedCharge;absNet[actual.positiveCell]+=std::abs(f.integratedCharge);}
        if(exterior)boundaryNet+=sign*static_cast<long double>(f.integratedCharge);fluxL1+=std::abs(f.integratedCharge);
    }
    for(std::size_t i=0;i<old.size();++i){bounded(static_cast<long double>(next[i])-old[i]+net[i],std::abs(next[i])+std::abs(old[i])+absNet[i],gate,
            "F4 independent cell signed advection ledger");source+=old[i];target+=next[i];l1+=std::abs(old[i])+std::abs(next[i]);}
    const long double residual=target-source+boundaryNet;
    bounded(residual,l1+fluxL1,gate,"F4 independent nonzero boundary charge advection");
    bounded(a.report().chargeResidual-residual,l1+fluxL1,1e-12L,"F4 report is not actual signed boundary-inclusive charge residual");
    bounded(a.report().boundaryIncomingCharge-incoming,fluxL1,1e-12L,"F4 signed incoming diagnostic");
    bounded(a.report().boundaryOutgoingCharge-outgoing,fluxL1,1e-12L,"F4 signed outgoing diagnostic");
    bounded(a.report().boundaryNetCharge-boundaryNet,fluxL1,1e-12L,"F4 signed boundary diagnostic");
    check(inRows>0&&outRows>0&&signedOpposite>0&&boundaryNet!=0&&target!=source,
        "F4 must include signed rho/Q and genuinely nonzero net injection");
    check(a.report().relativeChargeResidual<=gate&&a.report().relativeBoundaryModelChargeResidual<=gate&&a.report().maximumOutflowCfl<=.5,
        "F4 original physical charge/donor gates");
}
void checkField(const QuadtreeMesh& mesh,const LAMP::Electric::ElectricVofAdvectedCharge& adv,const ElectricVofField& field,bool zeroSigma) {
    const auto& old=adv.sourceCharge().integratedCharge();const auto& input=adv.integratedCharge();const auto& next=field.chargeCandidate().integratedCharge();
    check(field.advection()&&field.advection()->sameInstance(adv)&&field.sourceCharge().sameInstance(adv.sourceCharge())&&
        field.geometry().identity()==adv.geometry().identity()&&sameClock(field.geometryTime(),adv.timeLayer())&&sameClock(field.candidateTime(),adv.timeLayer()),
        "F4 Ohmic consumer lost actual open receipt/coeval geometry/source interval");
    field.validateFor(mesh,adv.geometry(),adv.sourceCharge(),adv.timeLayer());
    const long double dt=field.effectiveDt();check(dt==adv.timeLayer().time()-adv.sourceCharge().timeLayer().time(),"F4 fake conduction interval");
    std::vector<long double> d(old.size()),j(old.size()),absD(old.size()),absJ(old.size());long double boundaryJ=0,absBoundaryJ=0,boundaryQt=0,absBoundaryQt=0;
    for(const auto& f:field.faces()){const auto actual=adv.geometry().domain().face(mesh,f.nativeRow);
        bounded(f.displacement-static_cast<long double>(f.epsilon)*f.integratedField,std::abs(f.displacement),1e-12L,"F4 original same-H displacement");
        bounded(f.current-static_cast<long double>(f.conductivity)*f.integratedField,std::abs(f.current),1e-12L,"F4 original same-H current");
        scatter(d,actual,f.displacement);scatter(j,actual,f.current);
        if(actual.negativeCell>=0){absD[actual.negativeCell]+=std::abs(f.displacement);absJ[actual.negativeCell]+=std::abs(f.current);}
        if(actual.positiveCell>=0){absD[actual.positiveCell]+=std::abs(f.displacement);absJ[actual.positiveCell]+=std::abs(f.current);}
        if(actual.negativeCell<0||actual.positiveCell<0){const int sign=actual.negativeCell>=0?1:-1;boundaryJ+=sign*static_cast<long double>(f.current);absBoundaryJ+=std::abs(f.current);}}
    for(const auto& f:adv.faces())if(f.negativeCell<0||f.positiveCell<0){const int sign=f.negativeCell>=0?1:-1;
        boundaryQt+=sign*static_cast<long double>(f.integratedCharge);absBoundaryQt+=std::abs(f.integratedCharge);}
    long double before=0,rhs=0,after=0,l1=0;
    for(std::size_t i=0;i<old.size();++i){bounded(d[i]-next[i],absD[i]+std::abs(next[i]),gate,"F4 independent actual Gauss");
        bounded(static_cast<long double>(next[i])-input[i]+dt*j[i],std::abs(next[i])+std::abs(input[i])+dt*absJ[i],gate,"F4 independent unchanged Qcandidate-Qadv+dtJ");
        before+=old[i];rhs+=input[i];after+=next[i];l1+=std::abs(old[i])+std::abs(input[i])+std::abs(next[i]);}
    bounded(after-rhs+dt*boundaryJ,l1+dt*absBoundaryJ,gate,"F4 original global Ohmic input charge ledger");
    bounded(after-before+boundaryQt+dt*boundaryJ,l1+absBoundaryQt+dt*absBoundaryJ,gate,"F4 independent accepted Q plus actual Qt and J");
    check(boundaryQt!=0&&
        field.report().maximumRelativeGaussResidual<=gate&&field.report().relativeLinearResidual<=1e-11,
        "F4 retained linear/Gauss gates and boundary-charge witness");
    if(zeroSigma)check(std::memcmp(input.data(),next.data(),input.size()*sizeof(double))==0,"F4 sigma0 changed actual Qadv bits");
    std::cout<<"ELECTRIC_VOF_OPEN_CHARGE_FIELD,"<<(adv.geometry().domain().axisymmetric()?"RZ":"XY")<<",stage,"<<(adv.timeLayer().position==Clock::Position::Midpoint?"Mid":"End")
        <<",source_Q,"<<before<<",Qadv,"<<rhs<<",candidate_Q,"<<after<<",boundary_Qt,"<<boundaryQt<<",conduction_dt,"<<dt<<",boundary_J,"<<boundaryJ
        <<",combined_equation,independently_verified,Gauss,"<<field.report().maximumRelativeGaussResidual<<'\n';
}
void sameClosedPayload(const LAMP::Electric::ElectricVofAdvectedCharge& a,const LAMP::Electric::ElectricVofAdvectedCharge& b) {
    const auto& x=a.integratedCharge();const auto& y=b.integratedCharge();check(!a.boundary()&&!b.boundary()&&!a.report().boundaryEnabled&&!b.report().boundaryEnabled,
        "F4 allClosed entry did not delegate the original owner path");
    check(x.size()==y.size()&&std::memcmp(x.data(),y.data(),x.size()*sizeof(double))==0&&a.faces().size()==b.faces().size(),"F4 allClosed changed charge bits");
    for(std::size_t r=0;r<a.faces().size();++r){const auto& p=a.faces()[r];const auto& q=b.faces()[r];
        check(p.nativeRow==q.nativeRow&&p.negativeCell==q.negativeCell&&p.positiveCell==q.positiveCell&&p.totalVolume==q.totalVolume&&p.integratedCharge==q.integratedCharge,
            "F4 allClosed changed actual native face payload");}
    check(a.report().chargeResidual==b.report().chargeResidual&&a.report().roundingResidual==b.report().roundingResidual&&
        a.report().relativeChargeResidual==b.report().relativeChargeResidual&&a.report().maximumOutflowCfl==b.report().maximumOutflowCfl,"F4 allClosed numerical reports changed");
}
void run() {
    std::size_t solves=0,guards=0,delegates=0,varyingScenes=0;
    for(bool rz:{false,true}){
        DictionaryMesh cfg;cfg.Nx=cfg.Ny=8;cfg.minCoord={rz?1.:0.,0};cfg.maxCoord={rz?2.:1.,1};cfg.maxInitialIterations=0;cfg.maxCellNum=4096;
        auto mesh=std::make_unique<QuadtreeMesh>(cfg);DictionaryMesh::RefineBox refine;refine.level=1;refine.min={rz?1.5:.5,0};refine.max=cfg.maxCoord;mesh->applyRefineBoxes({refine});
        const auto domain=domainFor(*mesh,rz);check(!domain.uniformLeaves(),"F4 actual open scene lacks real 2:1 mesh");
        const auto stocks=phaseStocks(*mesh,domain);const Clock start{0,.2,.2,Clock::Position::Initial},end{1,.2,.2625,Clock::Position::End};
        const auto source=Publication::prepareLocalPlicInventory(*mesh,domain,start,stocks,{});const auto saved=snapshot(source);
        std::vector<double> q(mesh->cellIndexNum());const double oldDensity=rz?-.125:.125,inletDensity=rz?.25:-.25;
        for(std::size_t i=0;i<q.size();++i)q[i]=oldDensity*domain.cellVolumes()[domain.cellRowByIndex()[i]];
        const auto charge=ElectricVofCharge::fromIntegrated(*mesh,source,q);std::vector<Mode> rules(domain.nativeFaceCount(),Mode::Interior);
        std::vector<double> normal(rules.size()),rho(rules.size());int firstIn=-1,firstOut=-1,firstInterior=-1;
        const double left=cfg.minCoord.x,right=cfg.maxCoord.x;
        for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(*mesh,row);
            normal[row]=f.horizontal?0:.125/(rz?f.center.x:1.);
            if(f.negativeCell>=0&&f.positiveCell>=0){firstInterior=row;continue;}
            if(f.horizontal)rules[row]=Mode::Closed;
            else if(f.center.x==left){rules[row]=Mode::PureLiquidReservoir;rho[row]=inletDensity;firstIn=row;}
            else{check(f.center.x==right,"F4 foreign external port");rules[row]=Mode::PureGasReservoir;firstOut=row;}
        }
        check(firstIn>=0&&firstOut>=0&&firstInterior>=0,"F4 scene lost actual port/interior rows");
        const auto fluidBoundary=FluidBoundary::prepare(*mesh,source,rules);
        const auto slab=Slab::prepareBoundary(*mesh,source,fluidBoundary,end,normal,normal);
        const auto volume=VofTransport::prepareLocalPlicBoundarySweep(*mesh,source,slab,1000,1.226);
        const auto boundary=ChargeBoundary::prepare(*mesh,charge,volume,rho);boundary.validateFor(*mesh,charge,volume);
        check(!boundary.allClosed()&&boundary.transport().sameInstance(volume)&&boundary.sourceCharge().sameInstance(charge),"F4 detached inlet config lost actual owner");
        const auto moved=ChargeTransport::prepareBoundary(*mesh,charge,volume,boundary);
        checkReceipt(*mesh,moved.midpoint(),boundary,inletDensity,true);checkReceipt(*mesh,moved.endpoint(),boundary,inletDensity,false);
        for(std::size_t row=0;row<moved.endpoint().faces().size();++row){const auto h=moved.midpoint().faces()[row].integratedCharge,f=moved.endpoint().faces()[row].integratedCharge;
            bounded(static_cast<long double>(f)-2.L*h,std::abs(f)+2.L*std::abs(h),1e-12L,"F4 steady signed half/full charge flux");}
        const ElectricVofMaterials material{2.,5.,rz?0.:1.,rz?0.:.25};
        for(const auto* a:{&moved.midpoint(),&moved.endpoint()}){std::vector<Boundary> bc;references(*mesh,a->geometry(),material,a->timeLayer().time()-start.time(),bc);
            const auto field=ElectricVofField::prepareOhmic(*mesh,a->geometry(),*a,material,bc);checkField(*mesh,*a,field,rz);++solves;}
        const std::vector<double> zero(rules.size());const auto closedSlab=Slab::prepare(*mesh,source,end,zero,zero);
        const auto closedVolume=VofTransport::prepareLocalPlicSweep(*mesh,source,closedSlab,1000,1.226);
        const auto closedBoundary=ChargeBoundary::prepare(*mesh,charge,closedVolume,zero);
        const auto closed=ChargeTransport::prepareClosed(*mesh,charge,closedVolume),delegated=ChargeTransport::prepareBoundary(*mesh,charge,closedVolume,closedBoundary);
        sameClosedPayload(closed.midpoint(),delegated.midpoint());sameClosedPayload(closed.endpoint(),delegated.endpoint());++delegates;
        const auto foreign=ElectricVofCharge::fromIntegrated(*mesh,source,q);
        rejected([&]{boundary.validateFor(*mesh,foreign,volume);},"F4 same-value foreign accepted charge allowed");++guards;
        const auto freshVolume=VofTransport::prepareLocalPlicBoundarySweep(*mesh,source,slab,1000,1.226);
        rejected([&]{boundary.validateFor(*mesh,charge,freshVolume);},"F4 fresh same-value transport owner allowed");++guards;
        auto badClock=moved.endpoint().timeLayer();++badClock.step;
        rejected([&]{moved.endpoint().validateFor(*mesh,charge,volume,badClock);},"F4 wrong complete receipt clock allowed");++guards;
        auto bad=rho;bad[firstInterior]=1;rejected([&]{ChargeBoundary::prepare(*mesh,charge,volume,bad);},"F4 nonzero unused interior rho allowed");++guards;
        bad=rho;bad[firstOut]=1;rejected([&]{ChargeBoundary::prepare(*mesh,charge,volume,bad);},"F4 outflow reservoir rho replaced actual donor");++guards;
        bad=rho;bad[firstIn]=std::numeric_limits<double>::infinity();rejected([&]{ChargeBoundary::prepare(*mesh,charge,volume,bad);},"F4 infinite signed rho allowed");++guards;
        bad=zero;bad[firstIn]=1;rejected([&]{ChargeBoundary::prepare(*mesh,charge,closedVolume,bad);},"F4 closed external rho allowed");++guards;
        auto varying=normal;for(auto& u:varying)u*=2;
        const auto varyingSlab=Slab::prepareBoundary(*mesh,source,fluidBoundary,end,normal,varying);
        const auto varyingVolume=VofTransport::prepareLocalPlicBoundarySweep(*mesh,source,varyingSlab,1000,1.226);
        // This is the real F1 linear-normal history, not a guessed half of a full flux.
        const auto varyingBoundary=ChargeBoundary::prepare(*mesh,charge,varyingVolume,rho);
        const auto varyingMoved=ChargeTransport::prepareBoundary(*mesh,charge,varyingVolume,varyingBoundary);
        varyingBoundary.validateFor(*mesh,charge,varyingVolume);
        checkReceipt(*mesh,varyingMoved.midpoint(),varyingBoundary,inletDensity,true);
        checkReceipt(*mesh,varyingMoved.endpoint(),varyingBoundary,inletDensity,false);
        const long double actualDt=static_cast<long double>(end.end)-end.begin;
        std::size_t nonSteadyRows=0;
        for(const auto& f:varyingVolume.faces()){
            const auto actual=domain.face(*mesh,f.nativeRow);
            const long double u0=normal[f.nativeRow],u1=varying[f.nativeRow],area=actual.area;
            const long double expectedHalf=area*actualDt*(3*u0+u1)/8;
            const long double expectedFull=area*actualDt*(u0+u1)/2;
            bounded(f.half.total-expectedHalf,std::abs(f.half.total)+std::abs(expectedHalf),1e-12L,
                "F4b independent half integral of actual linear native normal");
            bounded(f.full.total-expectedFull,std::abs(f.full.total)+std::abs(expectedFull),1e-12L,
                "F4b independent full integral of actual linear native normal");
            if(u0!=0){check(f.full.total!=2*f.half.total,"F4b varying actual T was falsely made steady half/full");++nonSteadyRows;}
        }
        check(nonSteadyRows>0&&varyingMoved.midpoint().faces()[firstIn].integratedCharge*2!=varyingMoved.endpoint().faces()[firstIn].integratedCharge,
            "F4b actual reservoir charge half/full lacks a nonsteady temporal witness");
        for(const auto* a:{&varyingMoved.midpoint(),&varyingMoved.endpoint()}){
            std::vector<Boundary> bc;references(*mesh,a->geometry(),material,a->timeLayer().time()-start.time(),bc);
            const auto field=ElectricVofField::prepareOhmic(*mesh,a->geometry(),*a,material,bc);
            checkField(*mesh,*a,field,rz);++solves;
        }
        ++varyingScenes;
        std::cout<<"ELECTRIC_VOF_LINEAR_NATIVE_CHARGE,"<<(rz?"RZ":"XY")<<",rows,"<<nonSteadyRows
            <<",inlet_half_T,"<<varyingVolume.faces()[firstIn].half.total<<",inlet_full_T,"<<varyingVolume.faces()[firstIn].full.total
            <<",inlet_half_Qt,"<<varyingMoved.midpoint().faces()[firstIn].integratedCharge
            <<",inlet_full_Qt,"<<varyingMoved.endpoint().faces()[firstIn].integratedCharge<<",actual_F1_linear_no_reversal\n";
        auto reverse=normal;reverse[firstIn]=-normal[firstIn];bool reversed=false;
        try{const auto badSlab=Slab::prepareBoundary(*mesh,source,fluidBoundary,end,normal,reverse);
            VofTransport::prepareLocalPlicBoundarySweep(*mesh,source,badSlab,1000,1.226);}
        catch(const std::invalid_argument& e){reversed=std::string(e.what()).find("local sweep requires discrete incompressible endpoint snapshots")!=std::string::npos;}
        check(reversed,"F4b single-face reversal must retain the actual endpoint-divergence gate");++guards;
        check(snapshot(source)==saved&&std::memcmp(q.data(),charge.integratedCharge().data(),q.size()*sizeof(double))==0,"F4 attempts changed original independent source phase/Q");
    }
    check(solves==8&&delegates==2&&guards==16&&varyingScenes==2,"F4b steady/linear native consumers and owner/reversal negatives incomplete");
    std::cout<<"ELECTRIC_VOF_OPEN_CHARGE_PASS,Ohmic_solves,"<<solves<<",closed_delegates,"<<delegates<<",guards,"<<guards
        <<",linear_native_scenes,"<<varyingScenes<<",actual_F1_XY_RZ,signed_inlet_and_donor,nonzero_net_injection,original_gate,1e-10,fixed_mesh_F1_linear_no_reversal,no_App_or_open_AMR_claim\n";
}
} // namespace ElectricVofBoundaryChargeContract

namespace ElectricVofDielectricInclusionContract {
using namespace ElectricVofContract;
constexpr long double radius=.23L,centerZ=.513L,contrast=2.5L,farField=1.L;
struct Rect {long double x0,x1,z0,z1;};
struct Vec {long double x=0,z=0;};
struct Oracle {long double phi=0,epsilon=0;Vec e;std::array<long double,4> tensor{};};
Rect rectangle(const QuadtreeMesh& mesh,int ci,bool rz) {
    const auto& c=mesh.cells(mesh.cellIndexToLabel[ci]);const auto lo=mesh.verts(c.cellToVert[0]).coord(),hi=mesh.verts(c.cellToVert[2]).coord();
    const long double xc=rz?0.L:.47L;return {lo.x-xc,hi.x-xc,lo.y-centerZ,hi.y-centerZ};
}
int phase(const Rect& b) {
    const auto distance=[](long double a,long double z){return a>0?a:(z<0?-z:0.L);};
    const auto x=distance(b.x0,b.x1),z=distance(b.z0,b.z1);
    const auto xFar=std::max(std::abs(b.x0),std::abs(b.x1)),zFar=std::max(std::abs(b.z0),std::abs(b.z1));
    if(xFar*xFar+zFar*zFar<radius*radius)return 0;
    if(x*x+z*z>radius*radius)return 1;
    return -1;
}
std::vector<long double> cuts(const Rect& b) {
    const long double low=std::max(b.z0,-radius),high=std::min(b.z1,radius);if(!(high>low))return {};
    std::vector<long double> c{low,high};
    for(long double x:{b.x0,b.x1})if(std::abs(x)<radius){const auto z=std::sqrt(radius*radius-x*x);
        if(-z>low&&-z<high)c.push_back(-z);if(z>low&&z<high)c.push_back(z);}
    std::sort(c.begin(),c.end());c.erase(std::unique(c.begin(),c.end()),c.end());return c;
}
// Independent closed-form box integrals of the analytic circle/sphere.
// No production initialization sampler, PLIC polygon, epsilon or E is an oracle.
long double liquid(const Rect& b,bool rz) {
    const auto c=cuts(b);long double volume=0;
    const auto primitive=[](long double z){const auto q=std::max(0.L,radius*radius-z*z);
        return .5L*(z*std::sqrt(q)+radius*radius*std::asin(std::max(-1.L,std::min(1.L,z/radius))));};
    for(std::size_t k=1;k<c.size();++k){const auto a=c[k-1],d=c[k],mid=(a+d)/2,extent=std::sqrt(std::max(0.L,radius*radius-mid*mid));
        if(rz){if(extent<=b.x0)continue;
            const auto polynomial=extent>=b.x1?(b.x1*b.x1-b.x0*b.x0)*(d-a):
                (radius*radius-b.x0*b.x0)*(d-a)-(d*d*d-a*a*a)/3;
            volume+=static_cast<long double>(CartesianMeasure::twoPi)*polynomial/2;
        }else{const auto high=std::min(b.x1,extent),low=std::max(b.x0,-extent);if(!(high>low))continue;
            const bool upperArc=extent<b.x1,lowerArc=-extent>b.x0;
            const auto constant=(upperArc?0.L:b.x1)-(lowerArc?0.L:b.x0);
            volume+=constant*(d-a)+(int(upperArc)+int(lowerArc))*(primitive(d)-primitive(a));}}
    check(std::isfinite(volume)&&volume>=0,"dielectric analytic box integral nonfinite/negative");return volume;
}
Oracle analytic(Vector2D point,bool rz) {
    const long double x=static_cast<long double>(point.x)-(rz?0.L:.47L),z=static_cast<long double>(point.y)-centerZ;
    const auto s2=x*x+z*z;const int m=rz?3:2;const auto A=m/(contrast+m-1),B=(contrast-1)/(contrast+m-1);
    Oracle o;
    if(s2<radius*radius){o.e={0,A*farField};o.phi=-A*farField*z;o.epsilon=contrast;}
    else{const auto s=std::sqrt(s2),factor=B*farField*std::pow(radius/s,m);
        o.phi=-farField*z+factor*z;o.e={factor*m*z*x/s2,farField+factor*(m*z*z/s2-1)};o.epsilon=1;}
    const auto xx=.5L*o.epsilon*(o.e.x*o.e.x-o.e.z*o.e.z);
    o.tensor={xx,-xx,o.epsilon*o.e.x*o.e.z,-.5L*o.epsilon*(o.e.x*o.e.x+o.e.z*o.e.z)};return o;
}
struct Error {
    std::size_t count=0;long double weight=0,phi=0,ls=0,e=0,et=0,tensor=0;
    long double maxE=0,maxEt=0,maxTensor=0;int worstE=-1;Vector2D worstPoint{0,0};
    void add(int ci,Vector2D point,long double capacity,const Oracle& ref,Vector2D gradient,
        const ElectricVofMaxwellSample& sample,long double potential,bool rz) {
        const long double x=static_cast<long double>(point.x)-(rz?0.L:.47L),z=static_cast<long double>(point.y)-centerZ,s=std::hypot(x,z);
        check(s>0&&capacity>0,"dielectric error sample has zero radius/capacity");const Vec tangent{z/s,-x/s};
        const auto ex=static_cast<long double>(sample.electricField.x)-ref.e.x,ez=static_cast<long double>(sample.electricField.y)-ref.e.z;
        const auto errE=std::hypot(ex,ez)/farField,errEt=std::abs(ex*tangent.x+ez*tangent.z)/farField;
        const auto errLs=std::hypot(-static_cast<long double>(gradient.x)-ref.e.x,-static_cast<long double>(gradient.y)-ref.e.z)/farField;
        const auto errPhi=(potential-ref.phi)/(farField*radius);
        const std::array<long double,4> actual{sample.tensor.xx,sample.tensor.yy,sample.tensor.xy,sample.tensor.hoop};long double t2=0;
        for(int k=0;k<4;++k){const auto d=(actual[k]-ref.tensor[k])/(farField*farField);t2+=(k==2?2:1)*d*d;}
        const auto errT=std::sqrt(t2);++count;weight+=capacity;phi+=capacity*errPhi*errPhi;ls+=capacity*errLs*errLs;
        e+=capacity*errE*errE;et+=capacity*errEt*errEt;tensor+=capacity*t2;
        if(errE>maxE){maxE=errE;worstE=ci;worstPoint=point;}maxEt=std::max(maxEt,errEt);maxTensor=std::max(maxTensor,errT);
    }
    long double rms(long double squared)const {check(count>0&&weight>0,"dielectric error region has no independent pure cells");return std::sqrt(squared/weight);}
    void print(const char* geometry,int n,const char* region,int ph)const {
        std::cout<<"ELECTRIC_DIELECTRIC_ERROR,"<<geometry<<",N,"<<n<<",region,"<<region<<",phase,"<<ph<<",samples,"<<count
            <<",physical_weight,"<<weight<<",phi_RMS,"<<rms(phi)<<",LS_E_RMS,"<<rms(ls)<<",RT0_E_RMS,"<<rms(e)<<",Et_RMS,"<<rms(et)
            <<",T_RMS,"<<rms(tensor)<<",RT0_E_max,"<<maxE<<",Et_max,"<<maxEt<<",T_max,"<<maxTensor
            <<",worst_E_cell,"<<worstE<<",worst_E_xr,"<<worstPoint.x<<",worst_E_yz,"<<worstPoint.y<<"\n";
    }
};
struct Result {Error ring,bulk;std::array<Error,2> phaseRing,phaseBulk;};
Result scene(int n,bool rz) {
    DictionaryMesh cfg;cfg.Nx=cfg.Ny=n;cfg.minCoord={0,0};cfg.maxCoord={1,1};cfg.maxInitialIterations=0;
    cfg.maxCellNum=n*n*2;auto mesh=std::make_unique<QuadtreeMesh>(cfg);const auto domain=domainFor(*mesh,rz);
    check(domain.uniformLeaves()&&mesh->cellIndexNum()==n*n,"dielectric fixture must be one fixed uniform native mesh");
    std::vector<Stock> stocks(mesh->cellIndexNum());std::vector<int> classification(stocks.size());long double sumL=0;
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto b=rectangle(*mesh,ci,rz);const auto C=static_cast<long double>(domain.cellVolumes()[domain.cellRowByIndex()[ci]]);
        const auto p=phase(b);classification[ci]=p;const auto L=p==0?C:(p==1?0.L:liquid(b,rz));
        check(L>=0&&L<=C,"independent analytic liquid integral outside actual capacity; no clipping allowed");
        stocks[ci]={L,C-L};sumL+=L;}
    const auto measurePi=static_cast<long double>(CartesianMeasure::twoPi)/2;
    const auto exactTotal=rz?4*measurePi*radius*radius*radius/3:measurePi*radius*radius;
    bounded(sumL-exactTotal,exactTotal,1e-12L,"independent analytic circle/sphere integral total");
    const Clock clock{0,0,0,Clock::Position::Initial};
    const auto pub=Publication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,{});const auto saved=snapshot(pub);
    const auto sourceQ=ElectricVofCharge::fromIntegrated(*mesh,pub,std::vector<double>(stocks.size()));
    const auto savedQ=sourceQ.integratedCharge();std::vector<Boundary> boundaries(domain.nativeFaceCount());std::size_t mixed=0;
    for(std::size_t ci=0;ci<classification.size();++ci)if(classification[ci]<0)++mixed;
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(*mesh,row);
        boundaries[row]=f.negativeCell>=0&&f.positiveCell>=0?Boundary{Type::Interior,0}:
            (f.area==0?Boundary{Type::HomogeneousNeumann,0}:Boundary{Type::Dirichlet,static_cast<double>(analytic(f.center,rz).phi)});}
    const ElectricVofMaterials material{double(contrast),1,0,0};
    const auto field=ElectricVofField::prepareElectrostaticSnapshot(*mesh,pub,sourceQ,material,boundaries);
    field.validateFor(*mesh,pub,sourceQ,clock);
    check(field.isElectrostaticSnapshot()&&field.effectiveDt()==0&&field.chargeCandidate().sameInstance(sourceQ)
        &&field.sourceCharge().sameInstance(sourceQ)&&field.geometry().identity()==pub.identity()
        &&sameClock(field.geometryTime(),clock)&&sameClock(field.candidateTime(),clock)&&sameClock(field.chargeSourceTime(),clock),"dielectric snapshot owner/fullclock/dt0 mismatch");
    check(field.report().usedGMRES&&!field.report().usedSuperLU&&field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=1e-10,"dielectric original linear/Gauss gates");
    const auto stress=ElectricVofMaxwellStress::prepare(*mesh,field);stress.validateFor(*mesh,field,clock);
    std::vector<long double> divergence(stocks.size()),gross(stocks.size());long double peakGauss=0,faceH2=0,faceD2=0,faceWeight=0;
    std::size_t pureFaces=0,axis=0;const auto h=1.L/n;
    for(const auto& f:field.faces()){const auto native=domain.face(*mesh,f.nativeRow);scatter(divergence,native,f.displacement);
        if(native.negativeCell>=0)gross[native.negativeCell]+=std::abs(f.displacement);
        if(native.positiveCell>=0)gross[native.positiveCell]+=std::abs(f.displacement);
        check(f.current==0,"pure dielectric has nonzero Ohmic current");
        if(f.physicalMeasure==0){check(rz&&f.integratedField==0&&f.displacement==0,"dielectric true axis D/H not exact zero");++axis;continue;}
        int p=-2;bool eligible=true;for(int ci:{native.negativeCell,native.positiveCell})if(ci>=0){
            const int q=classification[ci];if(q<0||(p>=0&&p!=q))eligible=false;p=q;}
        if(!eligible||p<0)continue;
        const auto& edge=native.horizontal?mesh->hEdges(native.edgeLabel):mesh->vEdges(native.edgeLabel);
        const auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();
        const long double length=std::hypot(static_cast<long double>(b.x)-a.x,static_cast<long double>(b.y)-a.y);
        long double H=0,D=0;
        for(int k=0;k<8;++k){const auto u=(1.L+ElectricContract::nodes[k])/2;
            const Vector2D point{static_cast<double>(a.x+u*(b.x-a.x)),static_cast<double>(a.y+u*(b.y-a.y))};const auto ref=analytic(point,rz);
            const auto Af=length*ElectricContract::weights[k]/2*(rz?CartesianMeasure::twoPi*static_cast<long double>(point.x):1.L);
            const auto normalE=native.horizontal?ref.e.z:ref.e.x;H+=Af*normalE;D+=Af*ref.epsilon*normalE;}
        const auto measure=static_cast<long double>(f.physicalMeasure);
        const auto eh=(static_cast<long double>(f.integratedField)-H)/(measure*farField),ed=(static_cast<long double>(f.displacement)-D)/(measure*farField);
        faceH2+=measure*eh*eh;faceD2+=measure*ed*ed;faceWeight+=measure;++pureFaces;
    }
    for(std::size_t ci=0;ci<divergence.size();++ci){bounded(divergence[ci],gross[ci],1e-10L,"dielectric independent actual-face Gauss");
        if(gross[ci]>0)peakGauss=std::max(peakGauss,std::abs(divergence[ci])/gross[ci]);}
    check(pureFaces>0&&faceWeight>0&&(!rz||axis==static_cast<std::size_t>(n)),"dielectric pure native face/axis samples absent");
    Result result;std::size_t exactPure=0;
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const int p=classification[ci];if(p<0)continue;++exactPure;
        const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);const auto point=cell.centerCoord();const auto ref=analytic(point,rz);
        bounded(field.epsilonByCell()[ci]-(p==0?contrast:1.L),p==0?contrast:1.L,1e-12L,"analytic pure cell has different represented epsilon");
        const long double dx=static_cast<long double>(point.x)-(rz?0.L:.47L),dz=static_cast<long double>(point.y)-centerZ;
        const auto distance=std::abs(std::hypot(dx,dz)-radius);const auto C=static_cast<long double>(domain.cellVolumes()[domain.cellRowByIndex()[ci]]);
        const auto sample=stress.sample(ci,point);const auto gradient=field.cellPotentialGradients()[ci];const auto phi=field.potential()[ci];
        if(distance<=3*h){result.ring.add(ci,point,C,ref,gradient,sample,phi,rz);result.phaseRing[p].add(ci,point,C,ref,gradient,sample,phi,rz);}
        if(distance>2*h){result.bulk.add(ci,point,C,ref,gradient,sample,phi,rz);result.phaseBulk[p].add(ci,point,C,ref,gradient,sample,phi,rz);}
    }
    const auto foreign=ElectricVofCharge::fromIntegrated(*mesh,pub,savedQ);
    rejected([&]{field.validateFor(*mesh,pub,foreign,clock);},"dielectric same-value foreign Q owner accepted");
    const Clock wrong{1,0,1e-6,Clock::Position::End};
    rejected([&]{field.validateFor(*mesh,pub,sourceQ,wrong);},"dielectric wrong complete accepted clock allowed");
    const auto copy=Publication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,{});
    rejected([&]{field.validateFor(*mesh,copy,sourceQ,clock);},"dielectric fresh same-value VOF owner allowed");
    check(snapshot(pub)==saved&&std::memcmp(sourceQ.integratedCharge().data(),savedQ.data(),savedQ.size()*sizeof(double))==0
        &&field.chargeCandidate().sameInstance(sourceQ),"dielectric solve/sample/negative changed source inventories/normal/offset/Q bits");
    const char* geometry=rz?"RZ_sphere":"XY_cylinder";
    std::cout<<std::setprecision(17)<<"ELECTRIC_DIELECTRIC_SCENE,"<<geometry<<",N,"<<n<<",cells,"<<stocks.size()<<",h,"<<h<<",a,"<<radius
        <<",epsilon_ratio,"<<contrast<<",K,0,Q,0,analytic_liquid,"<<sumL<<",analytic_total_error,"<<sumL-exactTotal
        <<",actual_pure_cells,"<<exactPure<<",mixed_cells,"<<mixed<<",pure_faces,"<<pureFaces<<",axis_faces,"<<axis
        <<",H_RMS,"<<std::sqrt(faceH2/faceWeight)<<",D_RMS,"<<std::sqrt(faceD2/faceWeight)<<",Gauss,"<<peakGauss
        <<",linear,"<<field.report().relativeLinearResidual<<",solve_seconds,"<<field.report().totalSeconds
        <<",sharp_Et_jump,NA,body_flow,NOT_RUN\n";
    result.ring.print(geometry,n,"pure_ring_abs_distance_le3h",-1);result.bulk.print(geometry,n,"pure_bulk_abs_distance_gt2h",-1);
    for(int p=0;p<2;++p){result.phaseRing[p].print(geometry,n,"pure_ring_abs_distance_le3h",p);result.phaseBulk[p].print(geometry,n,"pure_bulk_abs_distance_gt2h",p);}
    return result;
}
void run() {
    std::array<std::array<Result,2>,2> results;
    // Complete all four diagnostics before evaluating the frozen physical gate.
    // No fitted solver tolerance, geometry rule or error threshold is adjusted.
    for(int axis=0;axis<2;++axis)for(int level=0;level<2;++level)results[axis][level]=scene(level?64:32,axis!=0);
    bool passed=true;
    for(int axis=0;axis<2;++axis){const auto& coarse=results[axis][0].ring;const auto& fine=results[axis][1].ring;
        const std::array<long double,3> c{{coarse.rms(coarse.e),coarse.rms(coarse.et),coarse.rms(coarse.tensor)}};
        const std::array<long double,3> f{{fine.rms(fine.e),fine.rms(fine.et),fine.rms(fine.tensor)}};
        const std::array<long double,3> gates{{.05L,.05L,.10L}};bool absolute=true,refined=true;
        for(int k=0;k<3;++k){absolute=absolute&&f[k]<=gates[k];refined=refined&&f[k]<=c[k];}
        std::cout<<"ELECTRIC_DIELECTRIC_FROZEN_GATE,"<<(axis?"RZ_sphere":"XY_cylinder")<<",fine_ring_E,"<<f[0]<<",limit_E,0.05,fine_ring_Et,"<<f[1]
            <<",limit_Et,0.05,fine_ring_T,"<<f[2]<<",limit_T,0.10,all_fine_le_coarse,"<<refined<<",absolute_pass,"<<absolute
            <<",sharp_Et_jump,NA,engineering_gate_not_paper_threshold\n";passed=passed&&absolute&&refined;}
    std::cout<<"ELECTRIC_DIELECTRIC_DIAGNOSTIC_COMPLETED,actual_snapshots,4,N,32:64,XY_cylinder_RZ_sphere,epsilon_ratio,2.5,K,0,free_Q,0,"
        <<"analytic_box_stocks,no_production_initialization_sampling,three_owner_guards_per_scene,source_bits_preserved,"
        <<"no_Hydro_no_drop_deformation_no_sharp_Et_jump_qualification\n";
    check(passed,"dielectric pure-ring physical E/Et/T frozen engineering gate FAILED; do not relax it");
    std::cout<<"ELECTRIC_DIELECTRIC_INCLUSION_PASS,physical_ring_engineering_gates_only,no_leaky_drop_or_Taylor_flow_claim\n";
}
} // namespace ElectricVofDielectricInclusionContract

}
namespace ElectricVofCounterflowContract {
using namespace ElectricVofContract;
using namespace ElectricVofMovingContract;
using ChargeBoundary=LAMP::Electric::ElectricVofChargeBoundary;
using Mode=Slab::Boundary::Mode;
struct Gross {long double positive=0,negative=0;};
// Independent integral of a linear normal on [0,s]. All chosen roots are
// exact binary quarters; this oracle does not call FaceTime or the new kernel.
Gross gross(long double first,long double last,long double s,long double dt,long double area) {
    std::vector<long double> cuts{0,s};const auto slope=last-first;
    if(slope!=0){const auto root=-first/slope;if(root>0&&root<s)cuts.push_back(root);}
    std::sort(cuts.begin(),cuts.end());Gross q;
    for(std::size_t k=1;k<cuts.size();++k){const auto l=cuts[k-1],r=cuts[k];
        const auto value=area*dt*(first*(r-l)+slope*(r*r-l*l)/2);
        if(first+slope*(l+r)/2>0)q.positive+=value;else q.negative-=value;
    }return q;
}
void verifyCharge(const QuadtreeMesh& mesh,const LAMP::Electric::ElectricVofAdvectedCharge& a,
    long double u0,long double u1,double reservoirDensity,bool half,std::size_t& zeroNetWitness) {
    const auto& t=a.transport();const auto& d=t.source().domain();const auto& old=a.sourceCharge().integratedCharge();
    const auto& next=a.integratedCharge();std::vector<long double> net(old.size()),activity(old.size());
    const long double dt=t.endTime().end-t.endTime().begin;long double before=0,after=0,in=0,out=0,canonicalBoundary=0,fluxActivity=0;
    for(const auto& f:a.faces()){const auto actual=d.face(mesh,f.nativeRow);const auto& accepted=half?t.faces()[f.nativeRow].half:t.faces()[f.nativeRow].full;
        const bool exterior=actual.negativeCell<0||actual.positiveCell<0;const int sign=actual.negativeCell>=0?1:-1;
        const auto expectedGross=gross(actual.horizontal?u0:0,actual.horizontal?u1:0,half?.5L:1.L,dt,actual.area);
        if(exterior){
            bounded(accepted.positiveTotal-expectedGross.positive,std::abs(accepted.positiveTotal)+std::abs(expectedGross.positive),1e-12L,"F15 independent positive temporal gross");
            bounded(accepted.negativeTotal-expectedGross.negative,std::abs(accepted.negativeTotal)+std::abs(expectedGross.negative),1e-12L,"F15 independent negative temporal gross");
        }else{
            const auto expectedNet=expectedGross.positive-expectedGross.negative;
            bounded(accepted.total-expectedNet,std::abs(accepted.total)+std::abs(expectedNet),1e-12L,"F15 unchanged interior frozen-mean net law");
            check(accepted.positiveTotal==std::max(accepted.total,0.L)&&accepted.negativeTotal==std::max(-accepted.total,0.L),
                "F15 interior was changed to chronological counterflow");
        }
        check(f.totalVolume==accepted.total,"F15 charge replaced actual canonical total volume");
        long double expected=0;
        if(exterior){const int cell=actual.negativeCell>=0?actual.negativeCell:actual.positiveCell;
            const auto c=d.cellVolumes()[d.cellRowByIndex()[cell]];
            const auto outward=sign>0?accepted.positiveTotal:accepted.negativeTotal,inward=sign>0?accepted.negativeTotal:accepted.positiveTotal;
            const auto oldOut=static_cast<long double>(old[cell])*outward/c;
            const auto reservoirIn=static_cast<long double>(reservoirDensity)*inward;
            expected=sign*(oldOut-reservoirIn);out+=oldOut;in+=reservoirIn;
            canonicalBoundary+=sign*static_cast<long double>(f.integratedCharge);
            if(accepted.total==0&&outward>0&&inward>0){check(f.integratedCharge!=0,"F15 net-zero volume erased real signed charge exchange");++zeroNetWitness;}
        }else if(accepted.total!=0){const int donor=accepted.total>0?actual.negativeCell:actual.positiveCell;
            expected=static_cast<long double>(old[donor])*accepted.total/d.cellVolumes()[d.cellRowByIndex()[donor]];}
        bounded(static_cast<long double>(f.integratedCharge)-expected,std::abs(expected)+std::abs(f.integratedCharge),1e-12L,
            "F15 independent old-source outgoing P0 plus signed reservoir incoming charge");
        if(actual.negativeCell>=0){net[actual.negativeCell]+=f.integratedCharge;activity[actual.negativeCell]+=std::abs(f.integratedCharge);}
        if(actual.positiveCell>=0){net[actual.positiveCell]-=f.integratedCharge;activity[actual.positiveCell]+=std::abs(f.integratedCharge);}
        fluxActivity+=std::abs(f.integratedCharge);
    }
    long double sourceL1=0,targetL1=0;
    for(std::size_t i=0;i<old.size();++i){bounded(static_cast<long double>(next[i])-old[i]+net[i],std::abs(next[i])+std::abs(old[i])+activity[i],gate,"F15 independent canonical cell charge ledger");
        before+=old[i];after+=next[i];sourceL1+=std::abs(old[i]);targetL1+=std::abs(next[i]);}
    const auto& report=a.report();const long double scale=sourceL1+targetL1+fluxActivity+std::abs(in)+std::abs(out);
    bounded(after-before+canonicalBoundary,scale,gate,"F15 independent canonical global charge ledger");
    bounded(after-before+out-in,scale,gate,"F15 independent physical gross boundary charge ledger");
    bounded(report.boundaryIncomingCharge-in,scale,1e-12L,"F15 signed incoming gross charge diagnostic");
    bounded(report.boundaryOutgoingCharge-out,scale,1e-12L,"F15 signed outgoing gross charge diagnostic");
    bounded(report.boundaryNetCharge-canonicalBoundary,scale,1e-12L,"F15 canonical boundary Qt diagnostic");
    bounded(canonicalBoundary-out+in-report.boundaryProductRoundingResidual-report.boundaryNetRoundingResidual,scale,1e-12L,
        "F15 gross-product/net RNE signs disagree with physical boundary ledger");
    check(report.relativeChargeResidual<=gate&&report.relativeBoundaryModelChargeResidual<=gate&&report.maximumOutflowCfl<=.5,
        "F15 changed the original physical charge/CFL gate");
    check(after!=before&&canonicalBoundary!=0,"F15 actual gross charge exchange needs nonzero physical injection/extraction");
}
void run() {
    std::size_t scenes=0,records=0,zeroNetWitness=0,solves=0,delegates=0,guards=0;
    const std::array<std::array<double,2>,4> patterns{{{{.125,-.125}},{{-.125,.125}},{{.125,-.375}},{{-.375,.125}}}};
    for(bool rz:{false,true}){
        auto mesh=meshFor();const auto domain=domainFor(*mesh,rz);check(!domain.uniformLeaves(),"F15 electric fixture lacks actual 2:1");
        std::vector<Stock> stocks(mesh->cellIndexNum());
        for(std::size_t i=0;i<stocks.size();++i){const auto c=domain.cellVolumes()[domain.cellRowByIndex()[i]];stocks[i]=rz?Stock{0,c}:Stock{c,0};}
        const Clock initial{0,.2,.2,Clock::Position::Initial},end{1,.2,.23125,Clock::Position::End};
        const auto source=Publication::prepareLocalPlicInventory(*mesh,domain,initial,stocks,{});const auto saved=snapshot(source);
        const double chargeSign=rz?-1.:1.,rhoValue=-chargeSign*.25;
        std::vector<double> q(stocks.size());for(std::size_t i=0;i<q.size();++i)q[i]=chargeSign*(1+static_cast<int>(i%3))/128.;
        const auto charge=ElectricVofCharge::fromIntegrated(*mesh,source,q);std::vector<Mode> rules(domain.nativeFaceCount(),Mode::Interior);
        for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(*mesh,row);
            if(f.negativeCell<0||f.positiveCell<0)rules[row]=f.horizontal?(rz?Mode::PureGasReservoir:Mode::PureLiquidReservoir):Mode::Closed;}
        const auto boundary=Slab::Boundary::prepare(*mesh,source,rules);
        for(std::size_t scene=0;scene<patterns.size();++scene){const auto values=patterns[scene];std::vector<double> first(rules.size()),last(rules.size());
            for(int row=0;row<domain.nativeFaceCount();++row)if(domain.face(*mesh,row).horizontal){first[row]=values[0];last[row]=values[1];}
            const auto slab=Slab::prepareBoundary(*mesh,source,boundary,end,first,last);
            const auto volume=VofTransport::prepareLocalPlicBoundarySweep(*mesh,source,slab,1000,1.226);
            if(scene==0||scene==2){bool actualNetZeroCounterflow=false;
                for(const auto& face:volume.faces())if(face.negativeCell<0||face.positiveCell<0){
                    const auto& selected=scene==0?face.full:face.half;
                    actualNetZeroCounterflow=actualNetZeroCounterflow||(selected.total==0&&selected.positiveTotal>0&&selected.negativeTotal>0);}
                check(actualNetZeroCounterflow,"F15 closed-entry negative lacks actual producer net-zero/nonzero-gross extent");
                bool closedRejected=false;try{ChargeTransport::prepareClosed(*mesh,charge,volume);}
                catch(const std::invalid_argument& error){closedRejected=std::string(error.what()).find("closed charge entry rejects nonzero exterior total flow")!=std::string::npos;}
                check(closedRejected,"F15 closed entry accepted a true open temporal-counterflow owner");++guards;
            }
            std::vector<double> rho(rules.size());for(const auto& f:volume.faces())if(f.negativeCell<0||f.positiveCell<0){const int sign=f.negativeCell>=0?1:-1;
                if((sign>0?f.full.negativeTotal:f.full.positiveTotal)>0)rho[f.nativeRow]=rhoValue;}
            const auto chargeBoundary=ChargeBoundary::prepare(*mesh,charge,volume,rho);
            const auto moved=ChargeTransport::prepareBoundary(*mesh,charge,volume,chargeBoundary);
            for(bool half:{true,false}){const auto& adv=half?moved.midpoint():moved.endpoint();
                verifyCharge(*mesh,adv,values[0],values[1],rhoValue,half,zeroNetWitness);++records;
                if(scene==0){const ElectricVofMaterials material{2,5,0,0};std::vector<LAMP::Electric::ElectricVofBoundary> bc;
                    references(*mesh,adv.geometry(),material,adv.timeLayer().time()-initial.time(),bc);
                    const auto field=ElectricVofField::prepareOhmic(*mesh,adv.geometry(),adv,material,bc);
                    ElectricVofBoundaryChargeContract::checkField(*mesh,adv,field,true);++solves;}
            }
            check(snapshot(source)==saved&&q==charge.integratedCharge(),"F15 successful counterflow changed accepted source phase/Q");
            ++scenes;
        }
        const std::vector<double> zero(rules.size());const auto closedSlab=Slab::prepare(*mesh,source,end,zero,zero);
        const auto closedVolume=VofTransport::prepareLocalPlicSweep(*mesh,source,closedSlab,1000,1.226);
        const auto closedBoundary=ChargeBoundary::prepare(*mesh,charge,closedVolume,zero);
        const auto closed=ChargeTransport::prepareClosed(*mesh,charge,closedVolume),delegated=ChargeTransport::prepareBoundary(*mesh,charge,closedVolume,closedBoundary);
        ElectricVofBoundaryChargeContract::sameClosedPayload(closed.midpoint(),delegated.midpoint());
        ElectricVofBoundaryChargeContract::sameClosedPayload(closed.endpoint(),delegated.endpoint());++delegates;
        auto bad=zero;for(int row=0;row<domain.nativeFaceCount();++row)if(rules[row]==Mode::Closed){bad[row]=1;break;}
        rejected([&]{ChargeBoundary::prepare(*mesh,charge,closedVolume,bad);},"F15 foreign closed/axis density allowed");++guards;
        check(snapshot(source)==saved&&q==charge.integratedCharge(),"F15 counterflow/negative attempts changed accepted source");
    }
    check(scenes==8&&records==16&&zeroNetWitness>0&&solves==4&&delegates==2&&guards==6,"F15 minimum true counterflow consumers incomplete");
    std::cout<<"ELECTRIC_VOF_COUNTERFLOW_PASS,scenes,"<<scenes<<",half_full_records,"<<records<<",zero_net_volume_nonzero_Qt,"<<zeroNetWitness
        <<",Ohmic,"<<solves<<",closed_delegates,"<<delegates<<",guards,"<<guards<<",XY_RZ_true_2to1,positive_negative_Q,roots_quarter_half_threequarter,original_gate,1e-10\n";
}
} // namespace ElectricVofCounterflowContract

namespace ElectricVofValidationBatchContract {
using namespace ElectricVofContract;
using ChargeBoundary=LAMP::Electric::ElectricVofChargeBoundary;
using ChargeTransport=LAMP::Electric::ElectricVofChargeTransport;
using Slab=I::GeometricVofAdvectorSlab;
using Transport=I::GeometricVofTransport;
using Mode=Slab::Boundary::Mode;
void run() {
    // Existing F4/F15 native helpers; only the missing independent public
    // mutable-mesh/clone boundary is new. Old selectors own physical oracles.
    auto mesh=meshFor(4);const auto domain=domainFor(*mesh,true);
    const Clock initial{0,0,0,Clock::Position::Initial},end{1,0,.03125,Clock::Position::End};
    std::vector<Stock> stocks(mesh->cellIndexNum());std::vector<double> q(stocks.size(),1./1024.);
    for(std::size_t i=0;i<stocks.size();++i)stocks[i]={0,domain.cellVolumes()[domain.cellRowByIndex()[i]]};
    const auto source=Publication::prepareLocalPlicInventory(*mesh,domain,initial,stocks,{});const auto saved=snapshot(source);
    const auto charge=ElectricVofCharge::fromIntegrated(*mesh,source,q);
    std::vector<Mode> rules(domain.nativeFaceCount(),Mode::Interior);std::vector<double> normal(rules.size()),rho(rules.size());
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(*mesh,row);
        if(f.negativeCell<0||f.positiveCell<0)rules[row]=f.horizontal?Mode::PureGasReservoir:Mode::Closed;
        if(f.horizontal)normal[row]=.0625;}
    const auto flow=Slab::Boundary::prepare(*mesh,source,rules);
    const auto slab=Slab::prepareBoundary(*mesh,source,flow,end,normal,normal);
    const auto volume=Transport::prepareLocalPlicBoundarySweep(*mesh,source,slab,1000,1.226);
    const auto boundary=ChargeBoundary::prepare(*mesh,charge,volume,rho);
    const auto moved=ChargeTransport::prepareBoundary(*mesh,charge,volume,boundary);const auto& adv=moved.endpoint();
    const auto halfOnly=ChargeTransport::prepareMidpointBoundary(*mesh,charge,volume,boundary);
    const ElectricVofMaterials material{2,5,0,0};std::vector<LAMP::Electric::ElectricVofBoundary> electric;references(*mesh,adv.geometry(),material,end.time(),electric);
    const auto field=ElectricVofField::prepareOhmic(*mesh,adv.geometry(),adv,material,electric);
    check(field.report().maximumRelativeGaussResidual<=gate,
        "Task32c original public fixture gates");
    const auto valid=[&](const QuadtreeMesh& actual){charge.validateFor(actual,initial);boundary.validateFor(actual,charge,volume);
        adv.validateFor(actual,charge,volume,end);field.validateFor(actual,adv.geometry(),charge,end);
        halfOnly.validateFor(actual,charge,volume,halfOnly.timeLayer());};
    valid(*mesh);auto clone=mesh->clone();valid(*clone); // exact contents, different C++ mesh address
    std::fenv_t before{},after{};const int error=errno;check(std::fegetenv(&before)==0,"Task32c native FP capture");
#if defined(_M_X64) || defined(__SSE2__)
    const unsigned csr=_mm_getcsr();
#endif
    valid(*mesh);check(std::fegetenv(&after)==0&&std::memcmp(&before,&after,sizeof before)==0&&errno==error,"Task32c validation changed FP/errno");
#if defined(_M_X64) || defined(__SSE2__)
    check(csr==_mm_getcsr(),"Task32c validation changed MXCSR");
#endif
    std::size_t guards=0,requestedCompletionGuards=0;
    const auto invalid=[&]{const auto reject=[&](auto action){bool caught=false;try{action();}catch(const std::invalid_argument&){caught=true;}
            check(caught,"Task32c independent public entry reused an obsolete live-mesh proof");++guards;};
        reject([&]{charge.validateFor(*mesh,initial);});reject([&]{boundary.validateFor(*mesh,charge,volume);});
        reject([&]{adv.validateFor(*mesh,charge,volume,end);});reject([&]{field.validateFor(*mesh,adv.geometry(),charge,end);});
        bool caught=false;try{(void)ChargeTransport::completeEndpoint(*mesh,halfOnly);}catch(const std::invalid_argument&){caught=true;}
        check(caught,"completeEndpoint reused stale public live-mesh proof");++requestedCompletionGuards;};
    auto& leaf=mesh->cells(mesh->cellIndexToLabel[0]);const auto center=leaf.centerCoord();leaf.setCenterCoord(center+Vector2D{.0009765625,0});
    invalid();leaf.setCenterCoord(center);valid(*mesh);
    int parent=-1;for(int i=0;i<mesh->cellIndexNum()&&parent<0;++i){const int p=mesh->cells(mesh->cellIndexToLabel[i]).parent;
        if(p>=0&&std::find(mesh->cellIndexToLabel,mesh->cellIndexToLabel+mesh->cellIndexNum(),p)==mesh->cellIndexToLabel+mesh->cellIndexNum())parent=p;}
    check(parent>=0,"Task32c fixture lacks a nonindexed ancestor");auto& ancestor=mesh->cells(parent);
    check(ancestor.son[0]>=0&&ancestor.son[1]>=0&&ancestor.son[0]!=ancestor.son[1],"Task32c actual ancestor children");
    std::swap(ancestor.son[0],ancestor.son[1]);invalid();std::swap(ancestor.son[0],ancestor.son[1]);valid(*mesh);
    auto& edge=mesh->vEdges(mesh->vEdgeIndexToLabel[0]);std::swap(edge.edgeToVert[0],edge.edgeToVert[1]);invalid();std::swap(edge.edgeToVert[0],edge.edgeToVert[1]);valid(*mesh);
    auto& vertex=mesh->verts(mesh->vertIndexToLabel[0]);const auto coordinate=vertex.coord();vertex.setCoord(coordinate+Vector2D{0,.0009765625});
    invalid();vertex.setCoord(coordinate);valid(*mesh);
    check(guards==16&&requestedCompletionGuards==4&&snapshot(source)==saved&&q==charge.integratedCharge(),"Task32c/requested completion live-mutation checks/source preservation incomplete");
    std::cout<<"ELECTRIC_VOF_VALIDATION_BATCH_PASS,RZaxis_true2to1,public_mutation_guards,"<<guards<<",content_clone,1,FP_unchanged,1,no_retained_proof"
        <<",requested_completion_live_mutation_guards,"<<requestedCompletionGuards<<"\n";
}
} // namespace ElectricVofValidationBatchContract

namespace ElectricVofPostprocessContract {
using namespace ElectricVofContract;
bool sameBits(double a,double b) {return std::memcmp(&a,&b,sizeof a)==0;}
template<class T>void bits(const std::vector<T>& a,const std::vector<T>& b,const char* message) {
    check(a.size()==b.size()&&(a.empty()||std::memcmp(a.data(),b.data(),a.size()*sizeof(T))==0),message);
}
void payload(const ElectricVofField& a,const ElectricVofField& b,bool withCellGradients=true) {
    check(a.geometry().identity()==b.geometry().identity()&&a.sourceCharge().sameInstance(b.sourceCharge())
        &&sameClock(a.geometryTime(),b.geometryTime())&&sameClock(a.chargeSourceTime(),b.chargeSourceTime())
        &&sameClock(a.candidateTime(),b.candidateTime())&&sameBits(a.effectiveDt(),b.effectiveDt()),"requested Field lost original owner/fullclock/duration");
    if(a.advection()||b.advection())check(a.advection()&&b.advection()&&a.advection()->sameInstance(*b.advection()),"requested Field lost actual advection owner");
    const auto& at=a.potentialTerms();const auto& bt=b.potentialTerms();
    check(at.size()==bt.size()&&std::memcmp(at.data(),bt.data(),at.size()*sizeof(std::array<double,3>))==0,
        "Task34 actual solved potential expansion changed between literal/parallel consumers");
    bits(a.potential(),b.potential(),"Task34 potential bits");bits(a.implicitInputCharge(),b.implicitInputCharge(),"Task34 implicit Q bits");
    bits(a.chargeCandidate().integratedCharge(),b.chargeCandidate().integratedCharge(),"Task34 candidate Q bits");
    bits(a.epsilonByCell(),b.epsilonByCell(),"Task34 epsilon bits");bits(a.conductivityByCell(),b.conductivityByCell(),"Task34 sigma bits");
    if(withCellGradients){check(a.cellPotentialGradients().size()==b.cellPotentialGradients().size(),"Task34 gradient shape");
        for(std::size_t i=0;i<a.cellPotentialGradients().size();++i)
            check(sameBits(a.cellPotentialGradients()[i].x,b.cellPotentialGradients()[i].x)
                &&sameBits(a.cellPotentialGradients()[i].y,b.cellPotentialGradients()[i].y),"Task34 LS gradient bits");}
    check(a.faces().size()==b.faces().size(),"Task34 face shape");
    for(std::size_t i=0;i<a.faces().size();++i){const auto& x=a.faces()[i];const auto& y=b.faces()[i];
        check(x.nativeRow==y.nativeRow&&x.negativeCell==y.negativeCell&&x.positiveCell==y.positiveCell,"Task34 face incidence");
        for(auto member:{&ElectricVofFace::physicalMeasure,&ElectricVofFace::epsilon,&ElectricVofFace::conductivity,
            &ElectricVofFace::implicitCoefficient,&ElectricVofFace::gradientConstant,&ElectricVofFace::integratedField,
            &ElectricVofFace::displacement,&ElectricVofFace::current})check(sameBits(x.*member,y.*member),"Task34 canonical face bits");
        check(x.gradientCoefficients.size()==y.gradientCoefficients.size(),"Task34 face coefficient size");
        for(std::size_t k=0;k<x.gradientCoefficients.size();++k)check(x.gradientCoefficients[k].nativeCell==y.gradientCoefficients[k].nativeCell
            &&sameBits(x.gradientCoefficients[k].value,y.gradientCoefficients[k].value),"Task34 actual gradient coefficient bits");}
    const auto& x=a.report();const auto& y=b.report();
    for(auto member:{&ElectricVofReport::relativeLinearResidual,&ElectricVofReport::maximumRelativeGaussResidual,
        &ElectricVofReport::maximumAffineGradientDefect,&ElectricVofReport::maximumPhaseCapacityDefect})
        check(sameBits(x.*member,y.*member),"Task34 retained physical report bits");
    check(x.cells==y.cells&&x.nativeFaces==y.nativeFaces&&x.nonorthogonalFaces==y.nonorthogonalFaces&&x.axisFaces==y.axisFaces
        &&x.usedSuperLU==y.usedSuperLU&&x.usedGMRES==y.usedGMRES&&x.usedAdvection==y.usedAdvection&&x.electrostaticSnapshot==y.electrostaticSnapshot,"Task34 retained report flags");
}
void executionReceipt(const ElectricVofField& field,unsigned requestedWorkers) {
    const auto& r=field.report();
    check(r.postprocessTimingAvailable&&r.integratedExpandedResidualAvailable&&r.integratedExpandedResidualTimingAvailable,
        "V56 actual postprocess/full-A timing receipt absent");
    for(double t:{r.postprocessFaceSeconds,r.postprocessCellSeconds,r.postprocessFinalizeSeconds,r.postprocessSeconds,
        r.integratedExpandedResidualRouteSeconds,r.integratedExpandedResidualSeconds})
        check(std::isfinite(t)&&t>=0,"V56 passive actual duration is not finite/nonnegative");
    check(r.postprocessSeconds>=std::max({r.postprocessFaceSeconds,r.postprocessCellSeconds,r.postprocessFinalizeSeconds}),
        "V56 nested phase timer exceeds its actual parent");
    check(r.integratedExpandedResidualRows==field.geometry().domain().cellIndexByRow().size(),
        "V56 full-A receipt did not cover the actual selected rows");
    check(r.integratedExpandedResidualScratchBytes<=64u*1024u*1024u,"V56 residual logical scratch exceeds optional cap");
    if(requestedWorkers==1)check(!r.integratedExpandedResidualParallelUsed&&!r.integratedExpandedResidualParallelFallback&&
        r.integratedExpandedResidualWorkers==0,"V56 literal residual claims a parallel join/fallback");
    else if(requestedWorkers>1)check(r.integratedExpandedResidualParallelUsed&&!r.integratedExpandedResidualParallelFallback&&
        r.integratedExpandedResidualWorkers>=2&&r.integratedExpandedResidualWorkers<=requestedWorkers&&
        r.integratedExpandedResidualScratchBytes>0,"V56 actual full-A parallel receipt missing");
    if(r.integratedExpandedResidualParallelUsed)check(r.integratedExpandedResidualWorkers>=2,
        "V56 unjoined team advertised as used");
}
void missingCellGradients(const ElectricVofField& field) {
    bool caught=false;try{(void)field.cellPotentialGradients();}
    catch(const std::invalid_argument& e){caught=std::string(e.what()).find("cell potential gradients were not requested")!=std::string::npos;}
    check(caught,"Omit cell-gradient getter fabricated empty/zero data or lost explicit missing-request rejection");
}
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
struct Runtime {
    std::fenv_t native{};unsigned csr=0;int error=errno,threads=omp_get_max_threads(),dynamic=omp_get_dynamic();
    Runtime(){check(std::fegetenv(&native)==0,"Task34 caller FP capture");csr=_mm_getcsr();}
    ~Runtime(){omp_set_dynamic(dynamic);omp_set_num_threads(threads);std::fesetenv(&native);_mm_setcsr(csr);errno=error;}
};
struct Receipt {
    std::fenv_t native{};unsigned csr=0;int error=0;
    void capture(){error=errno;csr=_mm_getcsr();check(std::fegetenv(&native)==0,"Task34 final FP capture");}
    bool same(const Receipt& b)const{return error==b.error&&csr==b.csr&&std::memcmp(&native,&b.native,sizeof native)==0;}
};
void run() {
    Runtime restore;omp_set_dynamic(0);
    std::size_t solves=0,guards=0,parallel=0,explicitRecover=0,omitted=0,omitParallel=0;
    for(bool rz:{false,true}){
        auto mesh=meshFor();const auto domain=domainFor(*mesh,rz);const auto stocks=phaseStocks(*mesh,domain);
        const Clock initial{0,0,0,Clock::Position::Initial},mid{1,0,.015625,Clock::Position::Midpoint};
        const auto accepted=Publication::prepareLocalPlicInventory(*mesh,domain,initial,stocks,{});
        const auto geometry=Publication::prepareLocalPlicInventory(*mesh,domain,mid,stocks,{});
        std::vector<double> q(mesh->cellIndexNum());
        for(std::size_t i=0;i<q.size();++i){const auto p=mesh->cells(mesh->cellIndexToLabel[static_cast<int>(i)]).centerCoord();
            q[i]=domain.cellVolumes()[domain.cellRowByIndex()[i]]*(p.x-.5+.19*(p.y-.5));}
        q.front()=-0.; // sigma=0 must preserve this accepted signed-zero payload.
        const auto charge=ElectricVofCharge::fromIntegrated(*mesh,accepted,q);const auto saved=snapshot(geometry);
        for(int kind=0;kind<3;++kind){
            ElectricVofMaterials material{2.,5.,kind==1?0.:.7,kind==1?0.:.03};
            std::vector<LAMP::Electric::ElectricVofBoundary> boundary;const auto& pub=kind==2?accepted:geometry;
            const auto refs=references(*mesh,pub,material,kind==2?0.:mid.time(),boundary);
            std::fenv_t entry{};check(std::fegetenv(&entry)==0,"Task34 common FP capture");const auto entryCsr=_mm_getcsr();
            const auto reset=[&]{check(std::fesetenv(&entry)==0,"Task34 common FP reset");_mm_setcsr(entryCsr);errno=173;};
            const auto prepare=[&]{return kind==2?ElectricVofField::prepareElectrostaticSnapshot(*mesh,pub,charge,material,boundary):
                ElectricVofField::prepareOhmic(*mesh,pub,charge,mid,material,boundary);};
            omp_set_num_threads(1);reset();const auto serial=prepare();Receipt sf;sf.capture();
            omp_set_num_threads(4);reset();const auto actual=prepare();Receipt pf;pf.capture();++solves;
            check(!serial.report().postprocessParallelUsed&&!serial.report().postprocessParallelFallback,"Task34 one-thread must be literal");
            check(actual.report().postprocessParallelUsed&&!actual.report().postprocessParallelFallback&&actual.report().postprocessWorkers>=2,
                "Task34 actual parallel execution missing");++parallel;
            check(actual.report().postprocessScratchBytes>0&&actual.report().postprocessScratchBytes<=64u*1024u*1024u,"Task34 scratch bound");
            executionReceipt(serial,1);executionReceipt(actual,4);
            ElectricVofOhmicQuantumContract::verify(*mesh,serial,kind==1,kind==2,false,false);
            ElectricVofOhmicQuantumContract::verify(*mesh,actual,kind==1,kind==2,false,false);
            payload(serial,actual);check(sf.same(pf),"Task34 complete FP/errno differs from literal");
            if(kind==2)check(actual.chargeCandidate().sameInstance(charge)&&serial.chargeCandidate().sameInstance(charge),"Task34 snapshot owner changed");
            else{check(actual.chargeCandidate().identity()==serial.chargeCandidate().identity()+1,"Task34 speculative candidate consumed extra Q identity");
                verify(*mesh,pub,charge,actual,mid,material,refs,false,kind==1);}
            if(kind==0){
                // The immutable geometry/structure offer contains no old
                // material, boundary values, Q, duration, matrix or factor.
                const Clock changedClock{1,0,.03125,Clock::Position::Midpoint};
                const auto changedGeometry=Publication::prepareLocalPlicInventory(*mesh,domain,changedClock,stocks,{});
                auto changedMaterial=material;changedMaterial.epsilonLiquid=3.125;changedMaterial.epsilonGas=4.75;
                changedMaterial.conductivityLiquid=.375;changedMaterial.conductivityGas=.0625;
                auto changedBoundary=boundary;bool changed=false;
                for(auto& bc:changedBoundary)if(bc.type==LAMP::Electric::ElectricVofBoundary::Type::Dirichlet){bc.value+=.0625;changed=true;}
                check(changed,"v80 Field requires actual changed Dirichlet values");
                const auto prepareChanged=[&](const ElectricVofField* offered){return ElectricVofField::prepareOhmic(*mesh,changedGeometry,
                    charge,changedClock,changedMaterial,changedBoundary,ElectricVofControls{},ElectricVofField::CellGradientRecovery::Recover,offered);};
                omp_set_num_threads(1);reset();const auto fresh=prepareChanged(nullptr);
                reset();const auto warmSerial=prepareChanged(&actual);Receipt ws;ws.capture();
                omp_set_num_threads(8);reset();const auto warmParallel=prepareChanged(&actual);Receipt wp;wp.capture();
                payload(fresh,warmSerial);payload(warmSerial,warmParallel);
                check(ws.same(wp),"v80 offered Field numeric1/8 native/rawCSR/errno changed");
                check(fresh.effectiveDt()!=actual.effectiveDt()&&fresh.epsilonByCell()!=actual.epsilonByCell()
                    &&fresh.potential()!=actual.potential(),"v80 Field dynamic material/BC/duration change retained old numeric payload");
                ElectricVofOhmicQuantumContract::verify(*mesh,warmParallel,false,false,false,false);
                warmParallel.validateFor(*mesh,changedGeometry,charge,changedClock);
            }
            actual.validateFor(*mesh,pub,charge,kind==2?initial:mid);
            auto wrong=kind==2?initial:mid;++wrong.step;
            rejected([&]{actual.validateFor(*mesh,pub,charge,wrong);},"Task34 foreign full clock");++guards;
            const auto other=ElectricVofCharge::fromIntegrated(*mesh,accepted,q);
            rejected([&]{actual.validateFor(*mesh,pub,other,kind==2?initial:mid);},"Task34 foreign source owner");++guards;
            using Request=ElectricVofField::CellGradientRecovery;
            const auto requested=[&](Request request){return kind==2?
                ElectricVofField::prepareElectrostaticSnapshot(*mesh,pub,charge,material,boundary,ElectricVofControls{},request):
                ElectricVofField::prepareOhmic(*mesh,pub,charge,mid,material,boundary,ElectricVofControls{},request);};
            // Preserve the original default Recover FP contract above.  Each
            // explicit request is compared only with its own serial FP history.
            omp_set_num_threads(1);reset();const auto recovered=requested(Request::Recover);Receipt rf;rf.capture();
            payload(serial,recovered);check(sf.same(rf),"default Recover differs from explicit original request FP/errno");++explicitRecover;
            reset();const auto omitSerial=requested(Request::Omit);Receipt osf;osf.capture();
            payload(serial,omitSerial,false);missingCellGradients(omitSerial);++omitted;
            check(!omitSerial.report().postprocessParallelUsed&&!omitSerial.report().postprocessParallelFallback,"Omit single thread was not literal");
            for(int workers:{2,8}){omp_set_num_threads(workers);reset();const auto omitParallelField=requested(Request::Omit);Receipt opf;opf.capture();
                executionReceipt(omitSerial,1);executionReceipt(omitParallelField,static_cast<unsigned>(workers));
                ElectricVofOhmicQuantumContract::verify(*mesh,omitParallelField,kind==1,kind==2,false,false);
                payload(omitSerial,omitParallelField,false);check(osf.same(opf),"Omit own-request 1/2/8 FP/errno mismatch");
                missingCellGradients(omitParallelField);++omitted;
                const auto& report=omitParallelField.report();check(report.postprocessParallelUsed&&!report.postprocessParallelFallback
                    &&report.postprocessWorkers==static_cast<std::size_t>(workers)&&report.postprocessScratchBytes>0
                    &&report.postprocessScratchBytes<=64u*1024u*1024u,"Omit actual 2/8-worker receipt/scratch limit");++omitParallel;
                omitParallelField.validateFor(*mesh,pub,charge,kind==2?initial:mid);
                if(kind==2)check(omitParallelField.chargeCandidate().sameInstance(charge),"Omit snapshot lost original accepted charge owner");
                else verify(*mesh,pub,charge,omitParallelField,mid,material,refs,false,kind==1);
            }
            check(charge.integratedCharge()==q&&snapshot(geometry)==saved,"Task34 source mutated");
        }
    }
    check(solves==6&&parallel==6&&guards==12&&explicitRecover==6&&omitted==18&&omitParallel==12,"Task34/requested-stage minimum evidence count");
    std::cout<<"ELECTRIC_VOF_POSTPROCESS_PARALLEL_PASS,scenes,2,paired_solves,"<<solves
        <<",parallel,"<<parallel<<",guards,"<<guards<<",payload_bits,exact,FP_errno,exact,scratch_limit,67108864,PDE_once_per_call"
        <<",explicit_default_Recover,"<<explicitRecover<<",Omit_missing_getters,"<<omitted<<",Omit_own_1_2_8_FP_parallel_pairs,"<<omitParallel<<"\n";
}
#else
void run(){throw std::runtime_error("Task34 parallel contract requires qualified MSVC x64 OpenMP build");}
#endif
} // namespace ElectricVofPostprocessContract

namespace ElectricVofChargeParallelContract {
using namespace ElectricVofContract;
using namespace ElectricVofMovingContract;
using ChargeBoundary=LAMP::Electric::ElectricVofChargeBoundary;
using Advected=LAMP::Electric::ElectricVofAdvectedCharge;
using Report=LAMP::Electric::ElectricVofChargeAdvectionReport;
using Mode=Slab::Boundary::Mode;
using ElectricVofPostprocessContract::sameBits;
using ElectricVofPostprocessContract::bits;
void payload(const Advected& a,const Advected& b) {
    check(a.sourceCharge().sameInstance(b.sourceCharge())&&a.transport().sameInstance(b.transport())
        &&a.geometry().identity()==b.geometry().identity()&&sameClock(a.timeLayer(),b.timeLayer()),"Task36 receipt owner/fullclock");
    if(a.boundary()||b.boundary())check(a.boundary()&&b.boundary()&&a.boundary()->sameInstance(*b.boundary()),"Task36 actual boundary owner");
    bits(a.integratedCharge(),b.integratedCharge(),"Task36 target charge bits");
    check(a.faces().size()==b.faces().size(),"Task36 native face count");
    for(std::size_t i=0;i<a.faces().size();++i){const auto& x=a.faces()[i];const auto& y=b.faces()[i];
        check(x.nativeRow==y.nativeRow&&x.negativeCell==y.negativeCell&&x.positiveCell==y.positiveCell
            &&std::memcmp(&x.totalVolume,&y.totalVolume,sizeof x.totalVolume)==0
            &&sameBits(x.integratedCharge,y.integratedCharge)&&sameBits(x.lowOrderIntegratedCharge,y.lowOrderIntegratedCharge)
            &&sameBits(x.highOrderIntegratedCharge,y.highOrderIntegratedCharge)&&sameBits(x.limiter,y.limiter),"Task36 actual shared face Qt/L/H/limiter bits");}
    const auto& x=a.report();const auto& y=b.report();
    check(x.boundaryEnabled==y.boundaryEnabled&&x.boundaryReservoirFaces==y.boundaryReservoirFaces
        &&x.cells==y.cells&&x.nativeFaces==y.nativeFaces&&x.donors==y.donors
        &&x.nonzeroFaceCharges==y.nonzeroFaceCharges&&x.quantizedZeroPieces==y.quantizedZeroPieces
        &&x.limitedLinearFct==y.limitedLinearFct&&x.reconstructedCells==y.reconstructedCells
        &&x.highOrderFaces==y.highOrderFaces&&x.limitedFaces==y.limitedFaces,"Task36 original report counts/scheme");
    for(auto member:{&Report::boundaryIncomingCharge,&Report::boundaryOutgoingCharge,&Report::boundaryNetCharge,
        &Report::boundaryProductRoundingResidual,&Report::maximumBoundaryProductRoundingError,
        &Report::boundaryNetRoundingResidual,&Report::maximumBoundaryNetRoundingError,&Report::boundaryModelChargeResidual,
        &Report::relativeBoundaryModelChargeResidual,&Report::maximumOutflowCfl,&Report::sourceTotalCharge,
        &Report::candidateTotalCharge,&Report::chargeResidual,&Report::relativeChargeResidual,&Report::roundingResidual,
        &Report::maximumCellRoundingError,&Report::minimumFaceLimiter,&Report::maximumDensityBoundDefect})
        check(sameBits(x.*member,y.*member),"Task36 original physical report bits");
}
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
using ElectricVofPostprocessContract::Runtime;
using ElectricVofPostprocessContract::Receipt;
struct Entry {
    std::fenv_t native{};unsigned csr=0;
    Entry(){check(std::feraiseexcept(FE_INEXACT)==0,"Task36 prime actual public inexact");
        csr=_mm_getcsr();check(std::fegetenv(&native)==0,"Task36 common FP capture");}
    void restore()const{check(std::fesetenv(&native)==0,"Task36 common FP restore");_mm_setcsr(csr);errno=173;}
};
template<class Prepare,class Verify>
void paired(const Prepare& prepare,const Verify& verify,std::size_t& pairs,std::size_t& parallel,bool requireDonor=true,bool requireParallel=true) {
    Entry entry;omp_set_num_threads(1);entry.restore();const auto serial=prepare();Receipt sf;sf.capture();
    omp_set_num_threads(4);entry.restore();const auto actual=prepare();Receipt pf;pf.capture();++pairs;
    check(sf.same(pf),"Task36 complete native/rawCSR/errno differs from literal");
    for(bool half:{true,false}){const auto& x=half?serial.midpoint():serial.endpoint();const auto& y=half?actual.midpoint():actual.endpoint();
        payload(x,y);check(!x.report().donorParallelUsed&&!x.report().targetParallelUsed&&!x.report().parallelFallback
            &&x.report().parallelWorkers==0&&x.report().parallelScratchBytes==0,"Task36 single thread was not literal");
        const auto& r=y.report();
        if(r.donorParallelUsed||r.targetParallelUsed)check(r.parallelWorkers>=2&&r.parallelScratchBytes>0
            &&r.parallelScratchBytes<=64u*1024u*1024u,"Task36 successful optional storage/worker report");
        if(requireParallel){check(!r.parallelFallback&&r.targetParallelUsed&&(!requireDonor||r.donorParallelUsed)
                &&r.parallelWorkers>=2,"Task36 qualified parallel loops not actually used");
            check(r.parallelScratchBytes>0&&r.parallelScratchBytes<=64u*1024u*1024u,"Task36 actual scratch cap");++parallel;}
        else if(r.parallelFallback){check(!r.donorParallelUsed&&!r.targetParallelUsed&&r.parallelWorkers==0&&r.parallelScratchBytes==0,
                "Task36 replay leaked optional statistics");}
        verify(y,half);
    }
}
void run() {
    Runtime restore;omp_set_dynamic(0);std::size_t pairs=0,parallel=0,guards=0,zeroNetWitness=0,signedZeroWitness=0;
    for(bool rz:{false,true}){
        auto mesh=meshFor();const auto domain=domainFor(*mesh,rz);check(!domain.uniformLeaves(),"Task36 fixture lacks true 2:1");
        std::vector<Stock> stocks(mesh->cellIndexNum());
        for(std::size_t i=0;i<stocks.size();++i){const auto C=domain.cellVolumes()[domain.cellRowByIndex()[i]];stocks[i]=rz?Stock{0,C}:Stock{C,0};}
        const Clock initial{0,.2,.2,Clock::Position::Initial},end{1,.2,.23125,Clock::Position::End};
        const auto source=Publication::prepareLocalPlicInventory(*mesh,domain,initial,stocks,{});const auto saved=snapshot(source);
        std::vector<double> mixed(stocks.size());for(std::size_t i=0;i<mixed.size();++i)mixed[i]=(i%2?-1.:1.)*(1+static_cast<int>(i%3))/128.;
        const auto closedQ=ElectricVofCharge::fromIntegrated(*mesh,source,mixed);
        const auto normal=normalFlow(*mesh,source,end.end-end.begin);
        const auto closedSlab=Slab::prepare(*mesh,source,end,normal,normal);
        const auto closed=VofTransport::prepareLocalPlicSweep(*mesh,source,closedSlab,1000,1.226);
        paired([&]{return ChargeTransport::prepareClosed(*mesh,closedQ,closed);},
            [&](const Advected& a,bool half){verifyAdvection(*mesh,a,half);},pairs,parallel);
        // Original public owner/fullclock guards remain independently callable.
        const auto actual=ChargeTransport::prepareClosed(*mesh,closedQ,closed);
        const auto foreign=ElectricVofCharge::fromIntegrated(*mesh,source,mixed);
        rejected([&]{actual.endpoint().validateFor(*mesh,foreign,closed,end);},"Task36 same-value foreign source");++guards;
        auto wrong=end;++wrong.step;
        rejected([&]{actual.endpoint().validateFor(*mesh,closedQ,closed,wrong);},"Task36 foreign complete clock");++guards;
        std::vector<Mode> rules(domain.nativeFaceCount(),Mode::Interior);
        for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(*mesh,row);
            if(f.negativeCell<0||f.positiveCell<0)rules[row]=f.horizontal?(rz?Mode::PureGasReservoir:Mode::PureLiquidReservoir):Mode::Closed;}
        const auto fluid=Slab::Boundary::prepare(*mesh,source,rules);
        const double sign=rz?-1.:1.,rhoValue=-sign*.25;
        std::vector<double> q(stocks.size());for(std::size_t i=0;i<q.size();++i)q[i]=sign*(1+static_cast<int>(i%3))/128.;
        const auto charge=ElectricVofCharge::fromIntegrated(*mesh,source,q);
        for(int pattern=0;pattern<2;++pattern){const double first=.125,last=pattern?-.125:.125;
            std::vector<double> u0(rules.size()),u1(rules.size());
            for(int row=0;row<domain.nativeFaceCount();++row)if(domain.face(*mesh,row).horizontal){u0[row]=first;u1[row]=last;}
            const auto slab=Slab::prepareBoundary(*mesh,source,fluid,end,u0,u1);
            const auto volume=VofTransport::prepareLocalPlicBoundarySweep(*mesh,source,slab,1000,1.226);
            std::vector<double> rho(rules.size());
            for(const auto& f:volume.faces())if(f.negativeCell<0||f.positiveCell<0){const int outward=f.negativeCell>=0?1:-1;
                if((outward>0?f.full.negativeTotal:f.full.positiveTotal)>0)rho[f.nativeRow]=rhoValue;}
            const auto boundary=ChargeBoundary::prepare(*mesh,charge,volume,rho);
            paired([&]{return ChargeTransport::prepareBoundary(*mesh,charge,volume,boundary);},
                [&](const Advected& a,bool half){if(pattern)ElectricVofCounterflowContract::verifyCharge(*mesh,a,first,last,rhoValue,half,zeroNetWitness);
                    else ElectricVofBoundaryChargeContract::checkReceipt(*mesh,a,boundary,rhoValue,half);},pairs,parallel);
        }
        // The actual source -0 survives the original net==0 bypass. Tiny
        // signed stocks additionally exercise finite native-double quanta.
        const double unit=std::numeric_limits<double>::denorm_min();
        auto tiny=mixed;for(std::size_t i=0;i<tiny.size();++i)tiny[i]=(i%2?-3.:5.)*unit;tiny.front()=-0.;
        const auto tinyQ=ElectricVofCharge::fromIntegrated(*mesh,source,tiny);
        const std::vector<double> zero(rules.size());const auto zeroSlab=Slab::prepare(*mesh,source,end,zero,zero);
        const auto still=VofTransport::prepareLocalPlicSweep(*mesh,source,zeroSlab,1000,1.226);
        paired([&]{return ChargeTransport::prepareClosed(*mesh,tinyQ,still);},
            [&](const Advected& a,bool){bits(tiny,a.integratedCharge(),"Task36 zero-net tiny/signed-zero bits");
                check(std::signbit(a.integratedCharge().front())&&a.integratedCharge().front()==0,"Task36 signed zero bypass");
                check(a.report().donors==0&&a.report().relativeChargeResidual==0,"Task36 zero-net physical ledger");++signedZeroWitness;},
            pairs,parallel,false,false);
        // FTZ selects the original literal path, whose existing binary64
        // environment gate REJECTS it. Do not mistake fallback for success.
        const auto ftz=_mm_getcsr();_mm_setcsr(ftz|0x8000u);
        Entry ftzEntry;std::string serialError,actualError;Receipt sf,pf;
        for(int threads:{1,4}){omp_set_num_threads(threads);ftzEntry.restore();
            try{const auto unexpected=ChargeTransport::prepareClosed(*mesh,closedQ,closed);(void)unexpected;}
            catch(const std::invalid_argument& error){(threads==1?serialError:actualError)=error.what();}
            if(threads==1)sf.capture();else pf.capture();
        }
        ++pairs;
        check(serialError=="[FiniteVolume] binary64 gradual underflow required"&&actualError==serialError,
            "Task36 unsupported FTZ did not preserve the original literal first exception");
        check(sf.same(pf),"Task36 rejected FTZ literal FP/errno differs between thread settings");
        _mm_setcsr(ftz);
        bits(mixed,closedQ.integratedCharge(),"Task36 closed accepted source mutated");
        bits(q,charge.integratedCharge(),"Task36 open accepted source mutated");
        check(snapshot(source)==saved,"Task36 accepted VOF geometry/stocks mutated");
    }
    check(pairs==10&&parallel==12&&guards==4&&zeroNetWitness>0&&signedZeroWitness==4,"Task36 required finite scenes missing");
    std::cout<<"ELECTRIC_VOF_CHARGE_PARALLEL_PASS,paired_factories,"<<pairs<<",actual_parallel_half_full,"<<parallel
        <<",foreign_guards,"<<guards<<",zero_net_gross_witness,"<<zeroNetWitness<<",signed_zero_records,"<<signedZeroWitness
        <<",exact_Qt_Q_report_bits,1,FP_errno_exact,1,scratch_limit,67108864,no_PDE\n";
}
#else
void run(){throw std::runtime_error("Task36 parallel contract requires MSVC x64 OpenMP");}
#endif
} // namespace ElectricVofChargeParallelContract

namespace ElectricVofSecondOrderChargeContract {
namespace FctExact=LAMP::Geometry::Detail::FiniteVolume;
using namespace ElectricVofContract;
using namespace ElectricVofMovingContract;
using Scheme=LAMP::Electric::ElectricVofChargeAdvectionScheme;
using Advected=LAMP::Electric::ElectricVofAdvectedCharge;
using Mode=Slab::Boundary::Mode;
using ChargeBoundary=LAMP::Electric::ElectricVofChargeBoundary;
std::size_t stages=0,highFaces=0,limitedFaces=0,parallelPairs=0;
ElectricVofControls controls(Scheme scheme){ElectricVofControls c;c.chargeAdvection=scheme;return c;}
bool sameBits(double a,double b){return std::memcmp(&a,&b,sizeof(double))==0;}
// Independent exact polynomial primitives: compact B(s)=(1-s*s)^4 and s*B(s).
long double primitive(long double s,bool odd){s=std::max(-1.L,std::min(1.L,s));const auto t=s*s;
    return odd?t*(.5L+t*(-1.L+t*(1.L+t*(-.5L+t*.1L)))):s*(1.L+t*(-4.L/3+t*(6.L/5+t*(-4.L/7+t/9))));}
long double integral(long double lo,long double hi,long double center,bool odd){constexpr long double r=.2L;
    return r*(primitive((hi-center)/r,odd)-primitive((lo-center)/r,odd));}
long double exactCell(const QuadtreeMesh& mesh,int ci,bool rz,long double time){
    const auto& c=mesh.cells(mesh.cellIndexToLabel[ci]);const auto x=c.centerCoord(),h=c.size();
    const long double x0=x.x-h.x/2,x1=x.x+h.x/2,z0=x.y-h.y/2,z1=x.y+h.y/2;
    const auto z=integral(z0,z1,.43L+(rz?.125L:.09375L)*time,true);
    return rz?CartesianMeasure::twoPi*(x1*x1-x0*x0)*z/2:integral(x0,x1,.43L+.125L*time,false)*z;
}
std::unique_ptr<QuadtreeMesh> uniformMesh(int n){DictionaryMesh cfg;cfg.Nx=cfg.Ny=n;cfg.minCoord={0,0};cfg.maxCoord={1,1};
    cfg.maxInitialIterations=0;cfg.maxCellNum=32768;return std::make_unique<QuadtreeMesh>(cfg);}
std::vector<Stock> gasStocks(const NativeLeafDomain& d){std::vector<Stock> s(d.cellRowByIndex().size());
    for(std::size_t i=0;i<s.size();++i)s[i]={0,d.cellVolumes()[d.cellRowByIndex()[i]]};return s;}
void receipt(const QuadtreeMesh& mesh,const Advected& a,bool half,bool bounds){
    const auto& old=a.sourceCharge().integratedCharge();const auto& next=a.integratedCharge();const auto& d=a.transport().source().domain();
    a.validateFor(mesh,a.sourceCharge(),a.transport(),a.timeLayer());
    check(a.geometry().identity()==(half?a.transport().midpoint():a.transport().endpoint()).identity(),"FCT actual stage geometry");
    check(a.report().limitedLinearFct,"FCT silently used the old P0 method");
    std::vector<long double> net(old.size()),activity(old.size()),lo(old.size()),hi(old.size()),lowRho(old.size());
    std::vector<FctExact::Dyadic> lowNet(old.size());
    for(const auto& f:a.faces()){
        if(f.negativeCell>=0)lowNet[f.negativeCell].add(f.lowOrderIntegratedCharge);
        if(f.positiveCell>=0)lowNet[f.positiveCell].add(f.lowOrderIntegratedCharge,-1);}
    long double globalLo=std::numeric_limits<long double>::infinity(),globalHi=-globalLo;
    for(std::size_t i=0;i<old.size();++i){const auto volume=d.cellVolumes()[d.cellRowByIndex()[i]];
        FctExact::Dyadic low(old[i]);low.add(lowNet[i],-1);
        lowRho[i]=static_cast<long double>(FctExact::roundDouble(low).value/volume);
        const auto sourceRho=static_cast<long double>(old[i]/volume);
        lo[i]=std::min(sourceRho,lowRho[i]);hi[i]=std::max(sourceRho,lowRho[i]);
        globalLo=std::min(globalLo,sourceRho);globalHi=std::max(globalHi,sourceRho);}
    long double before=0,after=0,boundary=0,l1=0;
    for(const auto& f:a.faces()){const auto& t=a.transport().faces()[f.nativeRow];const auto T=half?t.half.total:t.full.total;
        check(f.totalVolume==T&&f.negativeCell==t.negativeCell&&f.positiveCell==t.positiveCell,"FCT replaced actual shared VOF T/incidence");
        check(std::isfinite(f.limiter)&&f.limiter>=0&&f.limiter<=1,"FCT invalid limiter");
        const long double L=f.lowOrderIntegratedCharge,H=f.highOrderIntegratedCharge,F=f.integratedCharge;
        check(std::isfinite(L)&&std::isfinite(H),"FCT nonfinite diagnostic flux");
        bounded(F-(L+f.limiter*(H-L)),std::abs(L)+std::abs(H)+std::abs(F),32*std::numeric_limits<double>::epsilon(),"FCT reported low+theta*(high-low)");
        check(F>=std::min(L,H)&&F<=std::max(L,H),"FCT flux left low/high segment");
        FctExact::Dyadic difference(f.highOrderIntegratedCharge);difference.add(f.lowOrderIntegratedCharge,-1);
        const auto allowed=difference.multiplied(f.limiter);FctExact::Dyadic predicted(f.lowOrderIntegratedCharge);predicted.add(allowed);
        const double raw=FctExact::roundDouble(predicted).value;
        check(f.integratedCharge==raw||f.integratedCharge==std::nextafter(raw,f.lowOrderIntegratedCharge),"FCT face is not native RNE or its single inward neighbor");
        FctExact::Dyadic actualCorrection(f.integratedCharge);actualCorrection.add(f.lowOrderIntegratedCharge,-1);
        check(actualCorrection.absolute().compare(allowed.absolute())<=0
            &&(actualCorrection.compare(FctExact::Dyadic{})==0||actualCorrection.compare(FctExact::Dyadic{})==difference.compare(FctExact::Dyadic{})),
            "represented FCT shared correction overspends or reverses its exact budget");
        if(f.negativeCell>=0){net[f.negativeCell]+=F;activity[f.negativeCell]+=std::abs(F);}
        if(f.positiveCell>=0){net[f.positiveCell]-=F;activity[f.positiveCell]+=std::abs(F);}
        if(f.negativeCell>=0&&f.positiveCell>=0){const int i=f.negativeCell,j=f.positiveCell;
            const auto ri=static_cast<long double>(old[i]/d.cellVolumes()[d.cellRowByIndex()[i]]),rj=static_cast<long double>(old[j]/d.cellVolumes()[d.cellRowByIndex()[j]]);
            lo[i]=std::min({lo[i],rj,lowRho[j]});hi[i]=std::max({hi[i],rj,lowRho[j]});
            lo[j]=std::min({lo[j],ri,lowRho[i]});hi[j]=std::max({hi[j],ri,lowRho[i]});
        }else{const int cell=f.negativeCell>=0?f.negativeCell:f.positiveCell;boundary+=(f.negativeCell>=0?1:-1)*F;
            check(f.lowOrderIntegratedCharge==f.integratedCharge&&f.highOrderIntegratedCharge==f.integratedCharge,"FCT changed declared P0 gross boundary model");
            const auto& accepted=half?t.half:t.full;
            if(a.boundary()&&(f.negativeCell>=0?accepted.negativeTotal:accepted.positiveTotal)>0){
                const auto rho=static_cast<long double>(a.boundary()->inflowChargeDensities()[f.nativeRow]);
                lo[cell]=std::min(lo[cell],rho);hi[cell]=std::max(hi[cell],rho);globalLo=std::min(globalLo,rho);globalHi=std::max(globalHi,rho);}}
    }
    for(std::size_t i=0;i<old.size();++i){bounded(static_cast<long double>(next[i])-old[i]+net[i],std::abs(old[i])+std::abs(next[i])+activity[i],gate,"FCT signed cell ledger");
        if(bounds){const auto rho=static_cast<long double>(next[i]/d.cellVolumes()[d.cellRowByIndex()[i]]);
            const auto globalSlack=512*std::numeric_limits<double>::epsilon()*std::max(std::abs(globalLo),std::abs(globalHi));
            check(lowRho[i]>=globalLo-globalSlack&&lowRho[i]<=globalHi+globalSlack,"incompressible P0 baseline created a global signed density extremum");
            check(rho>=globalLo-globalSlack&&rho<=globalHi+globalSlack,"FCT created a global signed density extremum");
            const auto slack=512*std::numeric_limits<double>::epsilon()*std::max(std::abs(lo[i]),std::abs(hi[i]));
            if(!(rho>=lo[i]-slack&&rho<=hi[i]+slack)){
                const auto c=mesh.cells(mesh.cellIndexToLabel[i]).centerCoord();
                std::cerr<<std::setprecision(17)<<"FCT_BOUND_FAILURE,stage,"<<stages<<",half,"<<half<<",cells,"<<old.size()<<",cell,"<<i<<",x,"<<c.x<<",y,"<<c.y
                    <<",source_rho,"<<static_cast<long double>(old[i])/d.cellVolumes()[d.cellRowByIndex()[i]]<<",next_rho,"<<rho<<",lo,"<<lo[i]<<",hi,"<<hi[i]<<",slack,"<<slack<<'\n';
                throw std::runtime_error("FCT created a local signed density extremum");}}
        before+=old[i];after+=next[i];l1+=std::abs(old[i])+std::abs(next[i])+activity[i];}
    bounded(after-before+boundary,l1,gate,"FCT boundary-inclusive global charge ledger");
    bounded(a.report().chargeResidual-(after-before+boundary),l1,1e-12L,"FCT actual report charge residual");
    check(a.report().relativeChargeResidual<=gate&&a.report().maximumOutflowCfl<=.5,"FCT changed physical charge/CFL gates");
    check(std::isfinite(a.report().minimumFaceLimiter)&&a.report().minimumFaceLimiter>=0&&a.report().minimumFaceLimiter<=1
        &&std::isfinite(a.report().maximumDensityBoundDefect)&&a.report().maximumDensityBoundDefect>=0,"FCT invalid limiter/bounds report");
    highFaces+=a.report().highOrderFaces;limitedFaces+=a.report().limitedFaces;++stages;
}
struct Error{long double fct=0,p0=0;};
Error smooth(bool rz,int n){auto mesh=uniformMesh(n);const auto d=domainFor(*mesh,rz);const auto stocks=gasStocks(d);
    const Clock initial{0,0,0,Clock::Position::Initial};auto geometry=Publication::prepareLocalPlicInventory(*mesh,d,initial,stocks,{});
    std::vector<double> q(mesh->cellIndexNum());for(int i=0;i<mesh->cellIndexNum();++i)q[i]=static_cast<double>(exactCell(*mesh,i,rz,0));
    auto fctQ=ElectricVofCharge::fromIntegrated(*mesh,geometry,q),p0Q=fctQ;const int steps=n/2;constexpr double finalTime=.5;
    for(int step=1;step<=steps;++step){const Clock end{step,finalTime*(step-1)/steps,finalTime*step/steps,Clock::Position::End};
        std::vector<Mode> modes(d.nativeFaceCount(),Mode::Interior);std::vector<double> normal(d.nativeFaceCount()),rho(d.nativeFaceCount());
        for(int row=0;row<d.nativeFaceCount();++row){const auto f=d.face(*mesh,row);normal[row]=f.horizontal?(rz?.125:.09375):(rz?0.:.125);
            if(f.negativeCell<0||f.positiveCell<0)modes[row]=f.area==0?Mode::Closed:Mode::PureGasReservoir;}
        const auto fluidBoundary=Slab::Boundary::prepare(*mesh,geometry,modes);const auto slab=Slab::prepareBoundary(*mesh,geometry,fluidBoundary,end,normal,normal);
        const auto volume=VofTransport::prepareLocalPlicBoundarySweep(*mesh,geometry,slab,1000,1.226);
        const auto boundary=ChargeBoundary::prepare(*mesh,fctQ,volume,rho),p0Boundary=ChargeBoundary::prepare(*mesh,p0Q,volume,rho);
        const auto fct=ChargeTransport::prepareBoundary(*mesh,fctQ,volume,boundary,controls(Scheme::LimitedLinearFct));
        const auto p0=ChargeTransport::prepareBoundary(*mesh,p0Q,volume,p0Boundary,controls(Scheme::DonorCell));
        for(bool half:{true,false})receipt(*mesh,half?fct.midpoint():fct.endpoint(),half,true);
        geometry=volume.endpoint();fctQ=ElectricVofCharge::fromIntegrated(*mesh,geometry,fct.endpoint().integratedCharge());
        p0Q=ElectricVofCharge::fromIntegrated(*mesh,geometry,p0.endpoint().integratedCharge());}
    Error e;long double scale=0;for(int i=0;i<mesh->cellIndexNum();++i){const auto exact=exactCell(*mesh,i,rz,finalTime);
        e.fct+=std::abs(static_cast<long double>(fctQ.integratedCharge()[i])-exact);e.p0+=std::abs(static_cast<long double>(p0Q.integratedCharge()[i])-exact);scale+=std::abs(exact);}
    check(scale>0,"FCT smooth analytic fixture is empty");e.fct/=scale;e.p0/=scale;
    std::cout<<"ELECTRIC_VOF_FCT_SMOOTH,"<<(rz?"RZ_z":"XY_diagonal")<<",n,"<<n<<",steps,"<<steps<<",L1_FCT,"<<e.fct<<",L1_P0,"<<e.p0<<'\n';return e;
}
void signedAmr(bool rz){auto mesh=meshFor();const auto d=domainFor(*mesh,rz);check(!d.uniformLeaves(),"FCT signed fixture lacks actual2:1");
    const auto stocks=gasStocks(d);const Clock initial{0,.2,.2,Clock::Position::Initial},end{1,.2,.23125,Clock::Position::End};
    const auto geometry=Publication::prepareLocalPlicInventory(*mesh,d,initial,stocks,{});std::vector<double> q(mesh->cellIndexNum());
    for(int i=0;i<mesh->cellIndexNum();++i){const auto c=mesh->cells(mesh->cellIndexToLabel[i]).centerCoord();q[i]=d.cellVolumes()[d.cellRowByIndex()[i]]*(c.x+.37*c.y<.63?-2.5:.75);}
    const auto charge=ElectricVofCharge::fromIntegrated(*mesh,geometry,q);const auto normal=normalFlow(*mesh,geometry,end.end-end.begin);
    const auto slab=Slab::prepare(*mesh,geometry,end,normal,normal);const auto volume=VofTransport::prepareLocalPlicSweep(*mesh,geometry,slab,1000,1.226);
#if defined(_OPENMP)
    const auto oldThreads=omp_get_max_threads();omp_set_dynamic(0);omp_set_num_threads(1);
#endif
    const auto serial=ChargeTransport::prepareClosed(*mesh,charge,volume,controls(Scheme::LimitedLinearFct));
#if defined(_OPENMP)
    omp_set_num_threads(4);
    int actualWorkers=0;
#pragma omp parallel
    {
#pragma omp single
        actualWorkers=omp_get_num_threads();
    }
    check(actualWorkers==4,"FCT parallel comparison did not form its requested actual four-thread team");
#endif
    const auto parallel=ChargeTransport::prepareClosed(*mesh,charge,volume,controls(Scheme::LimitedLinearFct));
#if defined(_OPENMP)
    omp_set_num_threads(oldThreads);
#endif
    for(bool half:{true,false}){const auto& a=half?serial.midpoint():serial.endpoint();const auto& b=half?parallel.midpoint():parallel.endpoint();receipt(*mesh,a,half,true);receipt(*mesh,b,half,true);
        check(a.integratedCharge().size()==b.integratedCharge().size()&&std::memcmp(a.integratedCharge().data(),b.integratedCharge().data(),a.integratedCharge().size()*sizeof(double))==0,"FCT serial/parallel Q differs");
        for(std::size_t row=0;row<a.faces().size();++row){const auto& x=a.faces()[row];const auto& y=b.faces()[row];check(sameBits(x.integratedCharge,y.integratedCharge)&&sameBits(x.lowOrderIntegratedCharge,y.lowOrderIntegratedCharge)
            &&sameBits(x.highOrderIntegratedCharge,y.highOrderIntegratedCharge)&&sameBits(x.limiter,y.limiter),"FCT serial/parallel shared-face values differ");}++parallelPairs;}
    check(std::memcmp(q.data(),charge.integratedCharge().data(),q.size()*sizeof(double))==0,"FCT changed accepted source Q");
    auto wrong=end;++wrong.step;rejected([&]{parallel.endpoint().validateFor(*mesh,charge,volume,wrong);},"FCT wrong fullclock accepted");
    const auto foreign=ElectricVofCharge::fromIntegrated(*mesh,geometry,q);rejected([&]{parallel.endpoint().validateFor(*mesh,foreign,volume,end);},"FCT foreign source owner accepted");
}
// Merge inside ElectricVofSecondOrderChargeContract before run(); call boundaryCounterflow().
// This checks FCT's physical boundary ledger independently of its retained extra Ledger.
// Only the already-qualified finite-arithmetic CDF primitive is reused for the
// declared low-order outgoing branch; no potential solve or second PDE is added.
struct FctBoundaryCounters {
    std::size_t scenes=0,records=0,zeroNetNonzeroQt=0,productRounding=0,netRounding=0,nonzeroExtra=0,distinctModel=0;
};
std::vector<double> boundaryOutgoingOracle(const QuadtreeMesh& mesh,const ElectricVofCharge& q,
    const VofTransport& t,bool half){
    const auto& d=t.source().domain();std::vector<double> result(t.faces().size());
    std::vector<std::vector<std::pair<int,long double>>> outgoing(q.integratedCharge().size());
    for(const auto& f:t.faces()){const auto& v=half?f.half:f.full;
        if(f.negativeCell<0||f.positiveCell<0){const int cell=f.negativeCell>=0?f.negativeCell:f.positiveCell;
            const auto volume=f.negativeCell>=0?v.positiveTotal:v.negativeTotal;if(volume>0)outgoing[cell].push_back({f.nativeRow,volume});}
        else if(v.total!=0)outgoing[v.total>0?f.negativeCell:f.positiveCell].push_back({f.nativeRow,std::abs(v.total)});}
    for(std::size_t ci=0;ci<outgoing.size();++ci){auto& rows=outgoing[ci];if(rows.empty())continue;std::sort(rows.begin(),rows.end());
        std::vector<FctExact::Dyadic> weights;FctExact::Dyadic sum;
        for(const auto& r:rows){weights.emplace_back(r.second);sum.add(weights.back());}
        FctExact::Dyadic retained(d.cellVolumes()[d.cellRowByIndex()[ci]]);retained.add(sum,-1);weights.push_back(retained);
        const auto pieces=FctExact::partitionDouble(std::abs(q.integratedCharge()[ci]),weights);
        for(std::size_t k=0;k<rows.size();++k)result[rows[k].first]=std::copysign(pieces.values[k],q.integratedCharge()[ci]);
    }
    (void)mesh;return result;
}
void boundaryStageOracle(const QuadtreeMesh& mesh,const Advected& a,const Advected& p0,
    long double u0,long double u1,bool half,FctBoundaryCounters& counts){
    const auto& t=a.transport();const auto& d=t.source().domain();const auto& old=a.sourceCharge().integratedCharge();
    const auto& next=a.integratedCharge();const auto& low=p0.integratedCharge();const auto& report=a.report();const auto& baseline=p0.report();
    check(a.boundary()&&p0.boundary()&&a.boundary()->sameInstance(*p0.boundary()),"FCT/P0 do not share one actual signed charge boundary owner");
    check(report.boundaryEnabled==baseline.boundaryEnabled&&report.boundaryReservoirFaces==baseline.boundaryReservoirFaces,"FCT changed actual boundary role/count");
    using Report=LAMP::Electric::ElectricVofChargeAdvectionReport;
    for(auto member:{&Report::boundaryIncomingCharge,&Report::boundaryOutgoingCharge,&Report::boundaryNetCharge,
        &Report::boundaryProductRoundingResidual,&Report::maximumBoundaryProductRoundingError,
        &Report::boundaryNetRoundingResidual,&Report::maximumBoundaryNetRoundingError,&Report::maximumOutflowCfl,&Report::sourceTotalCharge})
        check(sameBits(report.*member,baseline.*member),"FCT changed an authenticated pure-boundary/P0 diagnostic bit");
    const auto out=boundaryOutgoingOracle(mesh,a.sourceCharge(),t,half);const auto& rho=a.boundary()->inflowChargeDensities();
    FctExact::Dyadic before,after,canonical,physicalOut,rawIn,representedIn,rawInMagnitude,counterOutMagnitude,productError,netError;
    FctExact::Dyadic scale;double maxProductError=0,maxNetError=0;
    std::vector<long double> lo(old.size()),hi(old.size());
    for(std::size_t i=0;i<old.size();++i){const auto c=d.cellVolumes()[d.cellRowByIndex()[i]];
        const auto r=static_cast<long double>(old[i])/c,l=static_cast<long double>(low[i])/c;lo[i]=std::min(r,l);hi[i]=std::max(r,l);
        before.add(old[i]);after.add(next[i]);scale.add(std::abs(old[i]));scale.add(std::abs(next[i]));}
    for(const auto& f:a.faces()){const auto& accepted=half?t.faces()[f.nativeRow].half:t.faces()[f.nativeRow].full;
        scale.add(std::abs(f.integratedCharge));
        if(f.negativeCell>=0&&f.positiveCell>=0){const int i=f.negativeCell,j=f.positiveCell;
            const auto ci=d.cellVolumes()[d.cellRowByIndex()[i]],cj=d.cellVolumes()[d.cellRowByIndex()[j]];
            for(const auto r:{static_cast<long double>(old[j])/cj,static_cast<long double>(low[j])/cj}){lo[i]=std::min(lo[i],r);hi[i]=std::max(hi[i],r);}
            for(const auto r:{static_cast<long double>(old[i])/ci,static_cast<long double>(low[i])/ci}){lo[j]=std::min(lo[j],r);hi[j]=std::max(hi[j],r);}continue;}
        const auto actual=d.face(mesh,f.nativeRow);const int sign=f.negativeCell>=0?1:-1,cell=f.negativeCell>=0?f.negativeCell:f.positiveCell;
        const auto incoming=sign>0?accepted.negativeTotal:accepted.positiveTotal,outgoing=sign>0?accepted.positiveTotal:accepted.negativeTotal;
        const auto analytic=ElectricVofCounterflowContract::gross(actual.horizontal?u0:0,actual.horizontal?u1:0,half?.5L:1.L,
            t.endTime().end-t.endTime().begin,actual.area);
        bounded(accepted.positiveTotal-analytic.positive,std::abs(accepted.positiveTotal)+std::abs(analytic.positive),1e-12L,"FCT boundary independent positive gross integral");
        bounded(accepted.negativeTotal-analytic.negative,std::abs(accepted.negativeTotal)+std::abs(analytic.negative),1e-12L,"FCT boundary independent negative gross integral");
        const auto exactIncoming=FctExact::Dyadic::product({rho[f.nativeRow],-sign*incoming});
        const auto inRne=FctExact::roundDouble(exactIncoming);FctExact::Dyadic flux(sign*out[f.nativeRow]);flux.add(inRne.value);
        const auto netRne=FctExact::roundDouble(flux);
        check(sameBits(f.integratedCharge,p0.faces()[f.nativeRow].integratedCharge)&&sameBits(f.integratedCharge,netRne.value),"FCT changed actual boundary Qt or disagrees with independent gross model");
        check(sameBits(f.lowOrderIntegratedCharge,f.integratedCharge)&&sameBits(f.highOrderIntegratedCharge,f.integratedCharge)&&f.limiter==1,"FCT applied an interior high-order correction at physical boundary");
        productError.add(inRne.error,sign);netError.add(netRne.error,sign);
        maxProductError=std::max(maxProductError,FctExact::roundDouble(inRne.error.absolute()).value);
        maxNetError=std::max(maxNetError,FctExact::roundDouble(netRne.error.absolute()).value);
        physicalOut.add(out[f.nativeRow]);rawIn.add(exactIncoming,-sign);representedIn.add(inRne.value,-sign);rawInMagnitude.add(exactIncoming.absolute());
        canonical.add(f.integratedCharge,sign);
        if(incoming>0){lo[cell]=std::min(lo[cell],static_cast<long double>(rho[f.nativeRow]));hi[cell]=std::max(hi[cell],static_cast<long double>(rho[f.nativeRow]));}
        if(incoming>0&&outgoing>0)counterOutMagnitude.add(std::abs(out[f.nativeRow]));
        if(accepted.total==0&&incoming>0&&outgoing>0){check(f.integratedCharge!=0,"FCT net-zero volume erased actual nonzero signed boundary exchange");++counts.zeroNetNonzeroQt;}
    }
    // Global shared Qt and raw physical reservoir products are separate exact ledgers.
    auto global=after;global.add(before,-1);global.add(canonical);
    auto model=after;model.add(before,-1);model.add(physicalOut);model.add(rawIn,-1);
    auto fromErrors=global;fromErrors.add(productError,-1);fromErrors.add(netError,-1);
    check(model.compare(fromErrors)==0,"FCT physical gross ledger and boundary product/net errors disagree exactly");
    auto modelScale=scale;modelScale.add(rawInMagnitude);modelScale.add(counterOutMagnitude);
    check(model.absolute().compare(modelScale.multiplied(1e-10))<=0,"FCT physical signed boundary model exceeds original exact budget");
    check(sameBits(report.chargeResidual,FctExact::roundDouble(global).value)&&sameBits(report.roundingResidual,FctExact::roundDouble(global).value),"FCT global/cell-RNE exact ledger report mismatch");
    check(sameBits(report.boundaryModelChargeResidual,FctExact::roundDouble(model).value),"FCT failed retained boundaryModelExtra residual reconstruction");
    const auto expectedRelative=FctExact::roundDouble(model.absolute()).value/FctExact::roundDouble(modelScale).value;
    check(sameBits(report.relativeBoundaryModelChargeResidual,expectedRelative),"FCT failed retained boundaryModelExtra activity/scale reconstruction");
    check(sameBits(report.boundaryIncomingCharge,FctExact::roundDouble(representedIn).value)&&sameBits(report.boundaryOutgoingCharge,FctExact::roundDouble(physicalOut).value)
        &&sameBits(report.boundaryNetCharge,FctExact::roundDouble(canonical).value),"FCT signed gross boundary diagnostics are not actual branch sums");
    check(sameBits(report.boundaryProductRoundingResidual,FctExact::roundDouble(productError).value)&&sameBits(report.maximumBoundaryProductRoundingError,maxProductError)
        &&sameBits(report.boundaryNetRoundingResidual,FctExact::roundDouble(netError).value)&&sameBits(report.maximumBoundaryNetRoundingError,maxNetError),"FCT exact product/net rounding diagnostics mismatch");
    for(std::size_t i=0;i<old.size();++i){const auto r=static_cast<long double>(next[i])/d.cellVolumes()[d.cellRowByIndex()[i]];
        const auto slack=512*std::numeric_limits<double>::epsilon()*std::max(std::abs(lo[i]),std::abs(hi[i]));
        check(r>=lo[i]-slack&&r<=hi[i]+slack,"FCT signed counterflow exceeds old+low local bounds");}
    if(productError.compare(FctExact::Dyadic{})!=0)++counts.productRounding;
    if(netError.compare(FctExact::Dyadic{})!=0)++counts.netRounding;
    auto extra=model;extra.add(global,-1);if(extra.compare(FctExact::Dyadic{})!=0)++counts.nonzeroExtra;
    if(!sameBits(report.boundaryModelChargeResidual,report.chargeResidual))++counts.distinctModel;
    receipt(mesh,a,half,false); // old+low bounds are checked independently just above.
    ++counts.records;
}
void boundaryCounterflow(){
    FctBoundaryCounters counts;const std::array<std::array<double,2>,4> patterns{{{{.125,-.125}},{{-.125,.125}},{{.125,-.375}},{{-.375,.125}}}};
    for(bool rz:{false,true}){auto mesh=meshFor();const auto d=domainFor(*mesh,rz);check(!d.uniformLeaves()&&mesh->cellIndexNum()==160,"FCT boundary fixture lacks genuine160-leaf2:1 mesh");
        const Clock initial{0,.2,.2,Clock::Position::Initial},end{1,.2,.23125,Clock::Position::End};
        const auto source=Publication::prepareLocalPlicInventory(*mesh,d,initial,gasStocks(d),{});const auto phaseBefore=snapshot(source);
        std::vector<double> q(mesh->cellIndexNum());bool positive=false,negative=false;
        for(std::size_t i=0;i<q.size();++i){const auto c=mesh->cells(mesh->cellIndexToLabel[i]).centerCoord();q[i]=(c.y<.5?-1.:1.)*(1+static_cast<int>(i%3))/128.;positive=positive||q[i]>0;negative=negative||q[i]<0;}
        check(positive&&negative,"FCT boundary source is not genuinely signed");const auto charge=ElectricVofCharge::fromIntegrated(*mesh,source,q);
        std::vector<Mode> rules(d.nativeFaceCount(),Mode::Interior);for(int row=0;row<d.nativeFaceCount();++row){const auto f=d.face(*mesh,row);
            if(f.negativeCell<0||f.positiveCell<0)rules[row]=f.horizontal?Mode::PureGasReservoir:Mode::Closed;}
        const auto fluid=Slab::Boundary::prepare(*mesh,source,rules);
        for(const double reservoir:{-.1,.1})for(const auto& pattern:patterns){std::vector<double> first(rules.size()),last(rules.size());
            for(int row=0;row<d.nativeFaceCount();++row)if(d.face(*mesh,row).horizontal){first[row]=pattern[0];last[row]=pattern[1];}
            const auto slab=Slab::prepareBoundary(*mesh,source,fluid,end,first,last);const auto volume=VofTransport::prepareLocalPlicBoundarySweep(*mesh,source,slab,1000,1.226);
            std::vector<double> rho(rules.size());for(const auto& f:volume.faces())if(f.negativeCell<0||f.positiveCell<0){const int sign=f.negativeCell>=0?1:-1;
                if((sign>0?f.full.negativeTotal:f.full.positiveTotal)>0)rho[f.nativeRow]=reservoir;}
            const auto boundary=ChargeBoundary::prepare(*mesh,charge,volume,rho);
            const auto fct=ChargeTransport::prepareBoundary(*mesh,charge,volume,boundary,controls(Scheme::LimitedLinearFct));
            const auto p0=ChargeTransport::prepareBoundary(*mesh,charge,volume,boundary,controls(Scheme::DonorCell));
            for(bool half:{true,false})boundaryStageOracle(*mesh,half?fct.midpoint():fct.endpoint(),half?p0.midpoint():p0.endpoint(),pattern[0],pattern[1],half,counts);
            check(q==charge.integratedCharge()&&snapshot(source)==phaseBefore,"FCT counterflow mutated accepted signed Q or pure-gas geometry");++counts.scenes;}
    }
    check(counts.scenes==16&&counts.records==32&&counts.zeroNetNonzeroQt>0&&counts.productRounding>0&&counts.netRounding>0&&counts.nonzeroExtra>0&&counts.distinctModel>0,"FCT nonzero signed boundary/model-extra evidence incomplete");
    std::cout<<"ELECTRIC_VOF_FCT_SIGNED_BOUNDARY_PASS,scenes,"<<counts.scenes<<",half_full_records,"<<counts.records<<",net_zero_nonzero_Qt,"<<counts.zeroNetNonzeroQt
        <<",product_rounding_records,"<<counts.productRounding<<",net_rounding_records,"<<counts.netRounding<<",nonzero_extra_records,"<<counts.nonzeroExtra<<",distinct_model_records,"<<counts.distinctModel
        <<",XY_RZ_160leaf_2to1,pure_gas,positive_negative_reservoir,signed_source,exact_physical_model_ledger,boundary_Qt_and_pure_diagnostics_bitwise_unchanged,no_PDE\n";
}
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
void parallelStages(){
    using namespace ElectricVofChargeParallelContract;
    Runtime restore;omp_set_dynamic(0);std::size_t pairs=0,actualStages=0,replayedErrors=0;
    for(bool rz:{false,true}){auto mesh=meshFor();const auto domain=domainFor(*mesh,rz);
        check(mesh->cellIndexNum()==160&&!domain.uniformLeaves(),"FCT staged pair needs actual160-cell2:1 native mesh");
        const Clock initial{0,.2,.2,Clock::Position::Initial},end{1,.2,.23125,Clock::Position::End};
        const auto source=Publication::prepareLocalPlicInventory(*mesh,domain,initial,gasStocks(domain),{});const auto saved=snapshot(source);
        std::vector<double> q(mesh->cellIndexNum());for(std::size_t i=0;i<q.size();++i)q[i]=(i%2?-1.:1.)*(1+static_cast<int>(i%3))/128.;
        const auto charge=ElectricVofCharge::fromIntegrated(*mesh,source,q);const auto normal=normalFlow(*mesh,source,end.end-end.begin);
        const auto slab=Slab::prepare(*mesh,source,end,normal,normal);const auto volume=VofTransport::prepareLocalPlicSweep(*mesh,source,slab,1000,1.226);
        const auto pair=[&](const auto& prepare,bool requireWorkers,bool physicalBounds){
            Entry entry;omp_set_num_threads(1);entry.restore();const auto serial=prepare();Receipt sf;sf.capture();
            omp_set_num_threads(8);entry.restore();const auto candidate=prepare();Receipt pf;pf.capture();
            check(sf.same(pf),"FCT staged complete native/rawCSR/errno differs from literal1thread");
            for(bool half:{true,false}){const auto& a=half?serial.midpoint():serial.endpoint();const auto& b=half?candidate.midpoint():candidate.endpoint();
                payload(a,b);receipt(*mesh,a,half,physicalBounds);receipt(*mesh,b,half,physicalBounds);
                check(!a.report().targetParallelUsed&&!a.report().parallelFallback&&a.report().parallelWorkers==0,"FCT staged serial path was not literal");
                if(requireWorkers){check(!b.report().parallelFallback&&b.report().targetParallelUsed&&b.report().parallelWorkers==8
                        &&b.report().parallelScratchBytes>0&&b.report().parallelScratchBytes<=64u*1024u*1024u,"FCT staged actual8-worker path missing");++actualStages;}
            }
            ++pairs;
        };
        pair([&]{return ChargeTransport::prepareClosed(*mesh,charge,volume,controls(Scheme::LimitedLinearFct));},true,true);
        std::vector<Mode> rules(domain.nativeFaceCount(),Mode::Interior);
        for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(*mesh,row);
            if(f.negativeCell<0||f.positiveCell<0)rules[row]=f.horizontal?Mode::PureGasReservoir:Mode::Closed;}
        const auto fluid=Slab::Boundary::prepare(*mesh,source,rules);std::vector<double> first(rules.size()),last(rules.size());
        for(int row=0;row<domain.nativeFaceCount();++row)if(domain.face(*mesh,row).horizontal){first[row]=.125;last[row]=-.125;}
        const auto openSlab=Slab::prepareBoundary(*mesh,source,fluid,end,first,last);
        const auto openVolume=VofTransport::prepareLocalPlicBoundarySweep(*mesh,source,openSlab,1000,1.226);
        std::vector<double> rho(rules.size());for(const auto& f:openVolume.faces())if(f.negativeCell<0||f.positiveCell<0){const int sign=f.negativeCell>=0?1:-1;
            if((sign>0?f.full.negativeTotal:f.full.positiveTotal)>0)rho[f.nativeRow]=-.1;}
        const auto boundary=ChargeBoundary::prepare(*mesh,charge,openVolume,rho);
        pair([&]{return ChargeTransport::prepareBoundary(*mesh,charge,openVolume,boundary,controls(Scheme::LimitedLinearFct));},true,false);
        std::vector<double> zero(rules.size());const auto stillSlab=Slab::prepare(*mesh,source,end,zero,zero);
        const auto still=VofTransport::prepareLocalPlicSweep(*mesh,source,stillSlab,1000,1.226);
        for(bool tiny:{false,true}){std::vector<double> values(q.size(),0.);
            if(tiny)for(std::size_t i=0;i<values.size();++i)values[i]=(i%2?-3.:5.)*std::numeric_limits<double>::denorm_min();values.front()=-0.;
            const auto stock=ElectricVofCharge::fromIntegrated(*mesh,source,values);
            pair([&]{return ChargeTransport::prepareClosed(*mesh,stock,still,controls(Scheme::LimitedLinearFct));},!tiny,false);
            check(std::signbit(stock.integratedCharge().front()),"FCT staged tiny/zero source lost its negative zero");
        }
        auto tinyTolerance=controls(Scheme::LimitedLinearFct);
        tinyTolerance.relativeConservationTolerance=std::numeric_limits<double>::denorm_min();
        pair([&]{return ChargeTransport::prepareClosed(*mesh,charge,still,tinyTolerance);},true,false);
        // Every worker range sees finite signed stocks whose Q/V overflows.
        // The failed private attempt must replay the literal first error/flags.
        std::vector<double> huge(q.size());for(std::size_t i=0;i<huge.size();++i)
            huge[i]=(i%2?-1.:1.)*std::numeric_limits<double>::max();
        const auto unrepresentable=ElectricVofCharge::fromIntegrated(*mesh,source,huge);
        Entry failedEntry;std::string serialFailure,parallelFailure;Receipt serialFp,parallelFp;
        for(int threads:{1,8}){omp_set_num_threads(threads);failedEntry.restore();
            try{(void)ChargeTransport::prepareClosed(*mesh,unrepresentable,still,controls(Scheme::LimitedLinearFct));}
            catch(const std::invalid_argument& e){(threads==1?serialFailure:parallelFailure)=e.what();}
            if(threads==1)serialFp.capture();else parallelFp.capture();}
        check(serialFailure=="[ElectricVofChargeTransport] unrepresentable signed Q/V for linear charge advection"
            &&parallelFailure==serialFailure&&serialFp.same(parallelFp),"FCT multichunk failure changed literal first error or FP/errno prefix");
        ++replayedErrors;failedEntry.restore();
        const auto csr=_mm_getcsr();_mm_setcsr(csr|0x8000u);Entry rejectedEntry;std::string literalError,parallelError;Receipt sf,pf;
        for(int threads:{1,8}){omp_set_num_threads(threads);rejectedEntry.restore();
            try{(void)ChargeTransport::prepareClosed(*mesh,charge,volume,controls(Scheme::LimitedLinearFct));}
            catch(const std::invalid_argument& e){(threads==1?literalError:parallelError)=e.what();}
            if(threads==1)sf.capture();else pf.capture();}
        check(literalError=="[FiniteVolume] binary64 gradual underflow required"&&parallelError==literalError&&sf.same(pf),
            "FCT staged rejected FTZ changed literal first error or FP/errno prefix");_mm_setcsr(csr);
        check(q==charge.integratedCharge()&&snapshot(source)==saved,"FCT staged attempts changed source Q/geometry");
    }
    check(pairs==10&&actualStages==16&&replayedErrors==2,"FCT staged signed/boundary/zero/tiny pair coverage incomplete");
    std::cout<<"ELECTRIC_VOF_FCT_STAGED_PARALLEL_PASS,pairs,"<<pairs<<",actual8worker_stages,"<<actualStages
        <<",XY_RZ_true2to1,shared_Qt_L_H_limiter_report_bits,FP_errno_exact,typed_boundary_gross,zero_subnormal_FTZ_first_error,denorm_min_tolerance,multichunk_overflow_replay\n";
}
#endif
void run(){check(ElectricVofControls{}.chargeAdvection==Scheme::DonorCell,"FCT changed direct module P0 default");
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    parallelStages();
#endif
    boundaryCounterflow();
    for(bool rz:{false,true}){const auto coarse=smooth(rz,16),middle=smooth(rz,32),fine=smooth(rz,64);const auto order=std::log(middle.fct/fine.fct)/std::log(2.L);
        check(std::isfinite(order)&&order>=1.7L,"FCT smooth fine-grid L1 order below1.7");
        check(fine.fct<fine.p0&&middle.fct<coarse.fct&&fine.fct<middle.fct,"FCT smooth refinement or gain over P0 missing");
        std::cout<<"ELECTRIC_VOF_FCT_ORDER,"<<(rz?"RZ_z":"XY_diagonal")<<",fine_L1_order,"<<order<<'\n';signedAmr(rz);}
    check(stages>0&&highFaces>0&&parallelPairs==4,"FCT high-order/AMR/half-full coverage missing");
    std::cout<<"ELECTRIC_VOF_CHARGE_SECOND_ORDER_PASS,stages,"<<stages<<",high_order_faces,"<<highFaces<<",limited_faces,"<<limitedFaces<<",serial_parallel_stage_pairs,"<<parallelPairs
        <<",actual_VOF_T,XY_diagonal_RZ_z,signed_old_low_local_and_global_bounds,shared_Qt,NO_OHMIC_TEMPORAL_ORDER_CLAIM\n";
}
} // namespace ElectricVofSecondOrderChargeContract
namespace ElectricVofSelectedDomainContract {
using namespace ElectricVofContract;
NativeLeafDomain selectedDomain(const QuadtreeMesh& mesh,bool rz,bool corner) {
    FaceSegmentGraph graph;
    graph.build(mesh,[corner](const QuadtreeMesh&,const Cell& cell){const auto p=cell.centerCoord();
        return p.y>.5||(corner&&p.x<.5);});
    return NativeLeafDomain::build(mesh,graph,rz);
}
std::vector<Stock> selectedStocks(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) {
    std::vector<Stock> stocks(mesh.cellIndexNum());
    for(int ci:domain.cellIndexByRow()){const auto volume=domain.cellVolumes()[domain.cellRowByIndex()[ci]];
        stocks[ci]=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord().y>.75?Stock{volume,0}:Stock{0,volume};}
    return stocks;
}
void run(Backend backend=Backend::Gmres) {
    ElectricVofControls controls;controls.backend=backend;
    std::size_t scenes=0,guards=0,internalFaces=0,excludedFaces=0,axisFaces=0;
    long double maximumPotentialError=0,maximumFieldError=0,maximumCharge=0;
    for(bool rz:{false,true})for(bool corner:{false,true}){
        auto mesh=meshFor();const auto domain=selectedDomain(*mesh,rz,corner);
        check(!domain.uniformLeaves()&&domain.cellIndexByRow().size()<static_cast<std::size_t>(mesh->cellIndexNum()),
            "E_SELECTED fixture must contain both excluded leaves and actual 2:1 split faces");
        const auto stocks=selectedStocks(*mesh,domain);
        const Clock start{0,0,0,Clock::Position::Initial},end{1,0,.002,Clock::Position::End};
        const auto publication=Publication::prepareLocalPlicInventory(*mesh,domain,start,stocks,{});
        const auto saved=snapshot(publication);std::vector<double> q(mesh->cellIndexNum());
        const auto charge=ElectricVofCharge::fromIntegrated(*mesh,publication,q);
        const ElectricVofMaterials material{2,2,.3,.3},insulator{2,2,0,0};
        std::vector<ElectricVofBoundary> bc(domain.nativeFaceCount());std::vector<int> internal,excluded;
        for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(*mesh,row);
            if(f.negativeCell<0&&f.positiveCell<0){excluded.push_back(row);++excludedFaces;continue;}
            if(f.negativeCell>=0&&f.positiveCell>=0)continue;
            if(f.area==0){bc[row]={Type::HomogeneousNeumann,0};++axisFaces;continue;}
            const bool actualInternal=f.center.x>0&&f.center.x<1&&f.center.y>0&&f.center.y<1;
            if(actualInternal){internal.push_back(row);++internalFaces;bc[row]={Type::Dirichlet,6};}
            else bc[row]={Type::Dirichlet,corner?6.5:6+f.center.y-.5};
        }
        check(!internal.empty()&&!excluded.empty(),"E_SELECTED missing actual internal conductor/excluded native faces");
        const auto field=ElectricVofField::prepareElectrostaticSnapshot(*mesh,publication,charge,material,bc,controls);
        check(field.report().cells==domain.cellIndexByRow().size()&&field.potential().size()==q.size()&&
            field.faces().size()==bc.size(),"E_SELECTED compact unknown count or dense public indexing lost");
        backendReceipt(field,backend);
        check(field.chargeCandidate().sameInstance(charge)&&field.report().relativeLinearResidual<=1e-11&&
            field.report().maximumRelativeGaussResidual<=gate,"E_SELECTED original snapshot solver/Gauss/owner gates failed");
        field.validateFor(*mesh,publication,charge,start);
        std::vector<long double> divergence(q.size()),gross(q.size());
        for(std::size_t row=0;row<field.faces().size();++row){const auto& f=field.faces()[row];const auto actual=domain.face(*mesh,static_cast<int>(row));
            check(f.nativeRow==static_cast<int>(row)&&f.negativeCell==actual.negativeCell&&f.positiveCell==actual.positiveCell,
                "E_SELECTED native face incidence changed");
            if(actual.negativeCell<0&&actual.positiveCell<0){check(f.physicalMeasure==0&&f.epsilon==0&&f.conductivity==0&&
                f.implicitCoefficient==0&&f.gradientCoefficients.empty()&&f.gradientConstant==0&&f.integratedField==0&&f.displacement==0&&f.current==0,
                "E_SELECTED both-excluded face fabricated a boundary/material/flux");continue;}
            if(actual.area==0)check(f.integratedField==0&&f.displacement==0&&f.current==0,"E_SELECTED true RZ axis has nonzero flux");
            scatter(divergence,actual,f.displacement);
            if(f.negativeCell>=0)gross[f.negativeCell]+=std::abs(f.displacement);
            if(f.positiveCell>=0)gross[f.positiveCell]+=std::abs(f.displacement);
            for(const auto& c:f.gradientCoefficients)check(domain.cellRowByIndex()[c.nativeCell]>=0,
                "E_SELECTED field gradient reached an excluded solid leaf");
            if(!corner){const long double expected=-actual.area*(actual.horizontal?1.L:0.L);
                maximumFieldError=std::max(maximumFieldError,std::abs(f.integratedField-expected));
                bounded(f.integratedField-expected,actual.area,2e-10L,"E_SELECTED independent affine integrated H or RZ measure failed");}
        }
        for(std::size_t ci=0;ci<q.size();++ci){const auto gradient=field.cellPotentialGradients()[ci];
            if(domain.cellRowByIndex()[ci]<0){check(field.potential()[ci]==0&&field.epsilonByCell()[ci]==0&&field.conductivityByCell()[ci]==0&&
                gradient.x==0&&gradient.y==0&&field.chargeCandidate().integratedCharge()[ci]==0,
                "E_SELECTED excluded cell is not an exact zero sentinel");continue;}
            bounded(divergence[ci],gross[ci],gate,"E_SELECTED independently scattered actual selected Gauss failed");
            if(!corner){const long double expected=6+mesh->cells(mesh->cellIndexToLabel[ci]).centerCoord().y-.5;
                maximumPotentialError=std::max(maximumPotentialError,std::abs(field.potential()[ci]-expected));
                check(std::abs(field.potential()[ci]-expected)<=2e-10L&&std::abs(gradient.x)<=2e-10&&std::abs(gradient.y-1)<=2e-10,
                    "E_SELECTED analytic affine potential/gradient failed");
            }else check(field.potential()[ci]>=6-2e-10&&field.potential()[ci]<=6.5+2e-10,
                "E_SELECTED actual L-shaped conductor solution violates independent maximum principle");
        }
        for(int row:internal){const auto& face=field.faces()[row];const auto actual=domain.face(*mesh,row);
            const int ci=actual.negativeCell>=0?actual.negativeCell:actual.positiveCell;
            const auto center=mesh->cells(mesh->cellIndexToLabel[ci]).centerCoord();
            const double distance=actual.horizontal?std::abs(actual.center.y-center.y):std::abs(actual.center.x-center.x);
            const double sign=actual.negativeCell>=0?1:-1;
            bounded(face.gradientConstant+sign*actual.area*6/distance,std::abs(face.gradientConstant),1e-12L,
                "E_SELECTED conductor distance/value uses inactive center rather than actual boundary face");}
        const auto endpoint=Publication::prepareLocalPlicInventory(*mesh,domain,end,stocks,{});
        const auto ohmic=ElectricVofField::prepareOhmic(*mesh,endpoint,charge,end,corner?insulator:material,bc,controls);
        backendReceipt(ohmic,backend);
        check(!ohmic.isElectrostaticSnapshot()&&sameClock(ohmic.chargeCandidate().timeLayer(),end)&&
            ohmic.report().maximumRelativeGaussResidual<=gate,"E_SELECTED positive-duration field/Gauss consumer failed");
        for(std::size_t ci=0;ci<q.size();++ci){const double next=ohmic.chargeCandidate().integratedCharge()[ci];
            if(domain.cellRowByIndex()[ci]<0||corner)check(next==0,"E_SELECTED inactive/zero-conductivity Q changed");
            else{maximumCharge=std::max(maximumCharge,static_cast<long double>(std::abs(next)));
                check(std::abs(next)<=1e-11*domain.cellVolumes()[domain.cellRowByIndex()[ci]],
                    "E_SELECTED charge-free uniform-conductivity affine solution created resolved charge");}}
        auto badQ=q;const auto inactive=std::find(domain.cellRowByIndex().begin(),domain.cellRowByIndex().end(),-1)-domain.cellRowByIndex().begin();
        badQ[inactive]=.125;
        rejected([&]{ElectricVofCharge::fromIntegrated(*mesh,publication,badQ);},"E_SELECTED nonzero solid charge accepted");++guards;
        auto bad=bc;bad[excluded.front()]={Type::Dirichlet,6};
        rejected([&]{ElectricVofField::prepareElectrostaticSnapshot(*mesh,publication,charge,material,bad);},
            "E_SELECTED both-excluded face accepted as conductor");++guards;
        bad=bc;bad[internal.front()]={Type::Interior,0};
        rejected([&]{ElectricVofField::prepareElectrostaticSnapshot(*mesh,publication,charge,material,bad);},
            "E_SELECTED selected/solid face accepted as interior");++guards;
        bad=bc;for(auto& b:bad)if(b.type!=Type::Interior)b={Type::HomogeneousNeumann,0};
        rejected([&]{ElectricVofField::prepareElectrostaticSnapshot(*mesh,publication,charge,material,bad);},
            "E_SELECTED unanchored component accepted");++guards;
        const auto fullDomain=domainFor(*mesh,rz);const auto fullStocks=phaseStocks(*mesh,fullDomain);
        const auto foreign=Publication::prepareLocalPlicInventory(*mesh,fullDomain,end,fullStocks,{});
        rejected([&]{ElectricVofField::prepareOhmic(*mesh,foreign,charge,end,material,bc);},
            "E_SELECTED foreign material mask accepted for selected source charge");++guards;
        rejected([&]{field.validateFor(*mesh,publication,charge,end);},"E_SELECTED wrong complete clock accepted");++guards;
        check(snapshot(publication)==saved&&charge.integratedCharge()==q,
            "E_SELECTED accepted phase/charge changed after successful/rejected preparations");
        ++scenes;
    }
    check(scenes==4&&guards==24&&internalFaces>0&&excludedFaces>0&&axisFaces>0,"E_SELECTED required XY/RZ actual conductor witnesses missing");
    std::cout<<"ELECTRIC_VOF_SELECTED_DOMAIN_PASS,scenes,"<<scenes<<",guards,"<<guards<<",internal_faces,"<<internalFaces
        <<",excluded_faces,"<<excludedFaces<<",axis_faces,"<<axisFaces<<",potential_error,"<<maximumPotentialError
        <<",integrated_H_error,"<<maximumFieldError<<",affine_ohmic_absQ,"<<maximumCharge<<",Gauss_gate,1e-10,linear_gate,1e-11\n";
}
} // namespace ElectricVofSelectedDomainContract
namespace ElectricVofSelectedChainContract {
using namespace ElectricVofContract;
using Slab=I::GeometricVofAdvectorSlab;
using Transport=I::GeometricVofTransport;
using Mode=Slab::Boundary::Mode;
using ChargeTransport=LAMP::Electric::ElectricVofChargeTransport;
using ChargeBoundary=LAMP::Electric::ElectricVofChargeBoundary;
using Advected=LAMP::Electric::ElectricVofAdvectedCharge;
std::vector<Stock> pureStocks(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,bool gas) {
    std::vector<Stock> stocks(mesh.cellIndexNum());for(int ci:domain.cellIndexByRow()){
        const auto v=domain.cellVolumes()[domain.cellRowByIndex()[ci]];stocks[ci]=gas?Stock{0,v}:Stock{v,0};}
    return stocks;
}
std::vector<double> circulation(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) {
    const auto psi=[](Vector2D p){const long double x=p.x,y=p.y;
        const auto a=x*(x-.5L)*(x-1),b=y*(y-.5L)*(y-1);return 64*a*a*b*b;};
    std::vector<double> velocity(domain.nativeFaceCount());
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);
        if(f.negativeCell<0||f.positiveCell<0)continue;
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if((f.horizontal?a.x:a.y)>(f.horizontal?b.x:b.y))std::swap(a,b);
        const auto integral=(f.horizontal?1.L:-1.L)*(psi(b)-psi(a))*(domain.axisymmetric()?CartesianMeasure::twoPi:1.);
        velocity[row]=static_cast<double>(integral/f.area);
    }
    return velocity;
}
void verify(const QuadtreeMesh& mesh,const Advected& a,bool half,bool open,double reservoir,std::size_t& zeroNetGross) {
    const auto& domain=a.geometry().domain();const auto& old=a.sourceCharge().integratedCharge();const auto& next=a.integratedCharge();
    std::vector<long double> net(old.size()),activity(old.size());long double before=0,after=0,l1=0,physicalBoundary=0,incoming=0,outgoing=0;
    for(const auto& face:a.faces()){const auto actual=domain.face(mesh,face.nativeRow);const auto& accepted=half?a.transport().faces()[face.nativeRow].half:a.transport().faces()[face.nativeRow].full;
        check(face.totalVolume==accepted.total,"E_SELECTED_CHAIN did not consume canonical actual VOF T");
        if(actual.negativeCell<0&&actual.positiveCell<0){check(face.totalVolume==0&&face.integratedCharge==0&&face.lowOrderIntegratedCharge==0&&face.highOrderIntegratedCharge==0,
            "E_SELECTED_CHAIN excluded face transported charge");continue;}
        if(actual.negativeCell>=0){net[actual.negativeCell]+=face.integratedCharge;activity[actual.negativeCell]+=std::abs(face.integratedCharge);}
        if(actual.positiveCell>=0){net[actual.positiveCell]-=face.integratedCharge;activity[actual.positiveCell]+=std::abs(face.integratedCharge);}
        if((actual.negativeCell>=0)!=(actual.positiveCell>=0)){
            const int sign=actual.negativeCell>=0?1:-1;const int ci=sign>0?actual.negativeCell:actual.positiveCell;
            const auto vin=sign>0?accepted.negativeTotal:accepted.positiveTotal,vout=sign>0?accepted.positiveTotal:accepted.negativeTotal;
            const auto oldOut=static_cast<long double>(old[ci])*vout/domain.cellVolumes()[domain.cellRowByIndex()[ci]];
            const auto newIn=static_cast<long double>(reservoir)*vin;const auto expected=sign*(oldOut-newIn);
            bounded(face.integratedCharge-expected,std::abs(expected)+std::abs(face.integratedCharge),1e-12L,
                "E_SELECTED_CHAIN independent signed gross P0 boundary law failed");
            incoming+=newIn;outgoing+=oldOut;physicalBoundary+=sign*static_cast<long double>(face.integratedCharge);
            const bool internal=actual.center.x>0&&actual.center.x<1&&actual.center.y>0&&actual.center.y<1;
            if(internal)check(vin==0&&vout==0&&face.integratedCharge==0,"E_SELECTED_CHAIN internal conductor wall exchanged fluid charge");
            if(accepted.total==0&&vin>0&&vout>0){check(face.integratedCharge!=0,"E_SELECTED_CHAIN zero net volume erased physical charge exchange");++zeroNetGross;}
        }
    }
    for(std::size_t ci=0;ci<old.size();++ci){if(domain.cellRowByIndex()[ci]<0){check(old[ci]==0&&next[ci]==0&&net[ci]==0,
            "E_SELECTED_CHAIN inactive cell contains transported charge");continue;}
        bounded(static_cast<long double>(next[ci])-old[ci]+net[ci],std::abs(old[ci])+std::abs(next[ci])+activity[ci],gate,
            "E_SELECTED_CHAIN independent actual shared-Qt cell ledger failed");
        const double density=next[ci]/domain.cellVolumes()[domain.cellRowByIndex()[ci]];
        check(density>=-.375-1e-12&&density<=(open?reservoir:.375)+1e-12,"E_SELECTED_CHAIN signed density gained a global extremum");
        before+=old[ci];after+=next[ci];l1+=std::abs(old[ci])+std::abs(next[ci])+activity[ci];}
    bounded(after-before+physicalBoundary,l1+std::abs(incoming)+std::abs(outgoing),gate,"E_SELECTED_CHAIN actual global signed-charge ledger failed");
    if(open){bounded(a.report().boundaryIncomingCharge-incoming,l1,1e-12L,"E_SELECTED_CHAIN gross incoming diagnostic failed");
        bounded(a.report().boundaryOutgoingCharge-outgoing,l1,1e-12L,"E_SELECTED_CHAIN gross outgoing diagnostic failed");
        check(a.report().boundaryEnabled&&a.report().relativeBoundaryModelChargeResidual<=gate,"E_SELECTED_CHAIN original open physical charge gate failed");}
    else check(physicalBoundary==0&&!a.report().boundaryEnabled,"E_SELECTED_CHAIN closed mode invented external work/flux");
    check(a.report().cells==domain.cellIndexByRow().size()&&a.report().relativeChargeResidual<=gate,
        "E_SELECTED_CHAIN active count or original charge gate failed");
    if(a.report().limitedLinearFct)check(a.report().reconstructedCells==domain.cellIndexByRow().size()&&a.report().maximumDensityBoundDefect==0,
        "E_SELECTED_CHAIN FCT rebuilt excluded cells or broke signed local bounds");
}
void run() {
    std::size_t stages=0,guards=0,zeroNetGross=0,highOrder=0,maxwellSamples=0,wallFaces=0;
    std::size_t requestedCompletions=0,requestedGuards=0,distinctEndpoint=0,typedOmit=0;
    for(bool rz:{false,true}){
        auto mesh=meshFor();const auto domain=ElectricVofSelectedDomainContract::selectedDomain(*mesh,rz,true);
        const auto stocks=pureStocks(*mesh,domain,rz);const Clock start{0,0,0,Clock::Position::Initial},end{1,0,.002,Clock::Position::End};
        const auto source=Publication::prepareLocalPlicInventory(*mesh,domain,start,stocks,{});const auto saved=snapshot(source);
        std::vector<double> q(mesh->cellIndexNum());for(int ci:domain.cellIndexByRow())q[ci]=domain.cellVolumes()[domain.cellRowByIndex()[ci]]*(ci%2?-.125:.125)*(1+ci%3);
        check(*std::min_element(q.begin(),q.end())<0&&*std::max_element(q.begin(),q.end())>0,"E_SELECTED_CHAIN needs signed source Q");
        const auto charge=ElectricVofCharge::fromIntegrated(*mesh,source,q);const double reservoir=.625;
        std::optional<ElectricVofField> preparedField;
        std::vector<Mode> rules(domain.nativeFaceCount(),Mode::Interior);std::vector<int> internal,excluded;
        std::vector<double> first(rules.size()),last(rules.size());
        std::vector<ElectricVofBoundary> electric(rules.size());
        for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(*mesh,row);
            if(f.negativeCell<0&&f.positiveCell<0){rules[row]=Mode::Closed;excluded.push_back(row);continue;}
            if(f.negativeCell>=0&&f.positiveCell>=0){if(f.horizontal&&f.center.x<.5){first[row]=.125;last[row]=-.125;}continue;}
            const bool wall=f.center.x>0&&f.center.x<1&&f.center.y>0&&f.center.y<1;
            rules[row]=Mode::Closed;
            if(wall){internal.push_back(row);++wallFaces;}
            if(f.horizontal&&(f.center.y==0||f.center.y==1)&&f.center.x<.5){rules[row]=rz?Mode::PureGasReservoir:Mode::PureLiquidReservoir;first[row]=.125;last[row]=-.125;}
            electric[row]=f.area==0?ElectricVofBoundary{Type::HomogeneousNeumann,0}:ElectricVofBoundary{Type::Dirichlet,wall?6.:6.5};
        }
        check(!internal.empty()&&!excluded.empty(),"E_SELECTED_CHAIN missing internal conductor or excluded rows");
        const auto fluid=Slab::Boundary::prepare(*mesh,source,rules);
        const auto openSlab=Slab::prepareBoundary(*mesh,source,fluid,end,first,last);
        const auto openVolume=Transport::prepareLocalPlicBoundarySweep(*mesh,source,openSlab,1000,1.226);
        const auto normals=circulation(*mesh,domain);const auto closedSlab=Slab::prepare(*mesh,source,end,normals,normals);
        const auto closedVolume=Transport::prepareLocalPlicSweep(*mesh,source,closedSlab,1000,1.226);
        std::vector<double> rho(rules.size());for(const auto& f:openVolume.faces())if((f.negativeCell>=0)!=(f.positiveCell>=0)){
            const int sign=f.negativeCell>=0?1:-1;if((sign>0?f.full.negativeTotal:f.full.positiveTotal)>0)rho[f.nativeRow]=reservoir;}
        const auto chargeBoundary=ChargeBoundary::prepare(*mesh,charge,openVolume,rho);
        for(bool open:{false,true})for(auto scheme:{ElectricVofChargeAdvectionScheme::DonorCell,ElectricVofChargeAdvectionScheme::LimitedLinearFct}){
            ElectricVofControls controls;controls.chargeAdvection=scheme;
            const auto moved=open?ChargeTransport::prepareBoundary(*mesh,charge,openVolume,chargeBoundary,controls):ChargeTransport::prepareClosed(*mesh,charge,closedVolume,controls);
            const auto halfOnly=open?ChargeTransport::prepareMidpointBoundary(*mesh,charge,openVolume,chargeBoundary,controls):
                ChargeTransport::prepareMidpointClosed(*mesh,charge,closedVolume,controls);
            controls.chargeAdvection=scheme==ElectricVofChargeAdvectionScheme::DonorCell?
                ElectricVofChargeAdvectionScheme::LimitedLinearFct:ElectricVofChargeAdvectionScheme::DonorCell;
            const auto completed=ChargeTransport::completeEndpoint(*mesh,halfOnly);
            check(completed.midpoint().sameInstance(halfOnly)
                &&completed.midpoint().report().totalSeconds==halfOnly.report().totalSeconds,
                "completion replaced genuine half owner or altered its measured duration");
            ElectricVofChargeParallelContract::payload(moved.midpoint(),completed.midpoint());
            ElectricVofChargeParallelContract::payload(moved.endpoint(),completed.endpoint());
            check(completed.sourceCharge().sameInstance(charge)&&completed.transport().sameInstance(open?openVolume:closedVolume),
                "completion lost original accepted Q/VOF transport owner");++requestedCompletions;
            if(halfOnly.integratedCharge()!=moved.endpoint().integratedCharge())++distinctEndpoint;
            rejected([&]{ChargeTransport::completeEndpoint(*mesh,moved.endpoint());},"completion accepted actual End instead of Midpoint");++requestedGuards;
            auto wrongHalf=halfOnly.timeLayer();++wrongHalf.step;
            rejected([&]{halfOnly.validateFor(*mesh,charge,open?openVolume:closedVolume,wrongHalf);},"requested half accepted wrong fullclock");++requestedGuards;
            rejected([&]{halfOnly.validateFor(*mesh,charge,open?closedVolume:openVolume,halfOnly.timeLayer());},"requested half accepted foreign genuine transport owner");++requestedGuards;
            for(bool half:{true,false}){const auto& a=half?moved.midpoint():moved.endpoint();verify(*mesh,a,half,open,reservoir,zeroNetGross);++stages;highOrder+=a.report().highOrderFaces;
                const ElectricVofMaterials material{2,2,.3,.3};const auto field=ElectricVofField::prepareOhmic(*mesh,a.geometry(),a,material,electric);
                if(!open&&half&&scheme==ElectricVofChargeAdvectionScheme::DonorCell){
                    preparedField.emplace(field);const auto savedPhi=field.potential(),savedQ=field.chargeCandidate().integratedCharge();
                    const auto offered=ElectricVofField::prepareOhmic(*mesh,a.geometry(),a,material,electric,ElectricVofControls{},
                        ElectricVofField::CellGradientRecovery::Recover,&field);
                    ElectricVofPostprocessContract::payload(field,offered);
                    auto changed=electric;bool changedBoundary=false;
                    for(auto& bc:changed)if(bc.type==Type::Dirichlet){bc.value+=.375;
                        if(rz&&!changedBoundary)bc={Type::HomogeneousNeumann,0};changedBoundary=true;}
                    check(changedBoundary,"prepared Field test needs an actual positive-area boundary change");
                    const auto freshBoundary=ElectricVofField::prepareOhmic(*mesh,a.geometry(),a,material,changed);
                    const auto offeredBoundary=ElectricVofField::prepareOhmic(*mesh,a.geometry(),a,material,changed,ElectricVofControls{},
                        ElectricVofField::CellGradientRecovery::Recover,&field);
                    ElectricVofPostprocessContract::payload(freshBoundary,offeredBoundary);
                    check(freshBoundary.potential()!=savedPhi,"prepared Field retained stale boundary values/types");
                    if(!rz){
                        std::fenv_t environment{};check(std::fegetenv(&environment)==0,"prepared Field rounding capture");const int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
                        const unsigned csr=_mm_getcsr();
#endif
                        const auto restore=[&]{std::fesetenv(&environment);
#if defined(_M_X64) || defined(__SSE2__)
                            _mm_setcsr(csr);
#endif
                            errno=error;};
                        try{check(std::fesetround(FE_DOWNWARD)==0,"prepared Field rounding change");
                            std::optional<ElectricVofField> fresh,offeredControl;std::string freshError,offeredError;
                            try{fresh.emplace(ElectricVofField::prepareOhmic(*mesh,a.geometry(),a,material,electric));}
                            catch(const std::exception& e){freshError=e.what();}
                            try{offeredControl.emplace(ElectricVofField::prepareOhmic(*mesh,a.geometry(),a,material,electric,ElectricVofControls{},
                                ElectricVofField::CellGradientRecovery::Recover,&field));}
                            catch(const std::exception& e){offeredError=e.what();}
                            check(fresh.has_value()==offeredControl.has_value()&&freshError==offeredError,"prepared Field changed fresh rounding-control acceptance");
                            if(fresh)ElectricVofPostprocessContract::payload(*fresh,*offeredControl);
                        }catch(...){restore();throw;}restore();
                    }
                    ElectricVofPostprocessContract::bits(field.potential(),savedPhi,"prepared Field offer mutated its old potential");
                    ElectricVofPostprocessContract::bits(field.chargeCandidate().integratedCharge(),savedQ,"prepared Field offer mutated its old charge");
                }
                if(scheme==ElectricVofChargeAdvectionScheme::LimitedLinearFct){
                    const auto omittedField=ElectricVofField::prepareOhmic(*mesh,a.geometry(),a,material,electric,ElectricVofControls{},
                        ElectricVofField::CellGradientRecovery::Omit);
                    ElectricVofPostprocessContract::payload(field,omittedField,false);
                    ElectricVofPostprocessContract::missingCellGradients(omittedField);++typedOmit;
                }
                check(field.report().usedAdvection&&field.report().maximumRelativeGaussResidual<=gate,
                    "E_SELECTED_CHAIN typed Ohmic field/charge consumer failed");
                const auto stress=LAMP::Electric::ElectricVofMaxwellStress::prepare(*mesh,field);
                stress.validateFor(*mesh,field,a.timeLayer());check(stress.report().cells==domain.cellIndexByRow().size()&&
                    stress.report().maximumRelativeSideMeasureDefect<=1e-12,"E_SELECTED_CHAIN Maxwell active side coverage failed");
                for(std::size_t ci=0;ci<q.size();++ci){const auto& c=stress.cells()[ci];
                    if(domain.cellRowByIndex()[ci]<0){check(c.epsilon==0&&c.sideDisplacement==std::array<long double,4>{}&&c.sidePhysicalMeasure==std::array<long double,4>{},
                        "E_SELECTED_CHAIN excluded cell fabricated Maxwell data");continue;}
                    const auto center=mesh->cells(mesh->cellIndexToLabel[ci]).centerCoord();const auto sample=stress.sample(static_cast<int>(ci),center);++maxwellSamples;
                    check(std::isfinite(sample.tensor.xx)&&std::isfinite(sample.tensor.yy)&&std::isfinite(sample.tensor.xy)&&std::isfinite(sample.tensor.hoop)&&
                        sample.tensor.hoop<=0&&sample.electricField.x==sample.displacement.x/c.epsilon&&sample.electricField.y==sample.displacement.y/c.epsilon,
                        "E_SELECTED_CHAIN selected Maxwell sample is not finite/coeval actual RT0");
                    bounded(static_cast<long double>(sample.tensor.xx)+sample.tensor.yy,std::abs(sample.tensor.xx)+std::abs(sample.tensor.yy),1e-12L,
                        "E_SELECTED_CHAIN Maxwell in-plane trace identity failed");
                    if(rz&&c.lower.x==0){const auto axis=stress.sample(static_cast<int>(ci),{0,center.y});check(axis.displacement.x==0&&axis.electricField.x==0,
                        "E_SELECTED_CHAIN actual RZ-axis Maxwell limit failed");}
                }
                const int inactive=static_cast<int>(std::find(domain.cellRowByIndex().begin(),domain.cellRowByIndex().end(),-1)-domain.cellRowByIndex().begin());
                rejected([&]{stress.sample(inactive,mesh->cells(mesh->cellIndexToLabel[inactive]).centerCoord());},"E_SELECTED_CHAIN inactive Maxwell sample accepted");++guards;
            }
            const auto foreignCharge=ElectricVofCharge::fromIntegrated(*mesh,source,q);
            rejected([&]{moved.endpoint().validateFor(*mesh,foreignCharge,open?openVolume:closedVolume,end);},"E_SELECTED_CHAIN foreign same-value Q owner accepted");++guards;
            rejected([&]{halfOnly.validateFor(*mesh,foreignCharge,open?openVolume:closedVolume,halfOnly.timeLayer());},"requested half accepted foreign same-value source Q owner");++requestedGuards;
        }
        auto bad=rho;bad[excluded.front()]=1;rejected([&]{ChargeBoundary::prepare(*mesh,charge,openVolume,bad);},"E_SELECTED_CHAIN excluded reservoir accepted");++guards;
        bad=rho;bad[internal.front()]=1;rejected([&]{ChargeBoundary::prepare(*mesh,charge,openVolume,bad);},"E_SELECTED_CHAIN internal Closed wall reservoir accepted");++guards;
        const auto fullDomain=domainFor(*mesh,rz);const auto fullSource=Publication::prepareLocalPlicInventory(*mesh,fullDomain,start,pureStocks(*mesh,fullDomain,rz),{});
        const auto fullCharge=ElectricVofCharge::fromIntegrated(*mesh,fullSource,std::vector<double>(q.size()));
        if(!rz){
            const ElectricVofMaterials material{2,2,.3,.3};std::vector<ElectricVofBoundary> fullBoundary;
            references(*mesh,fullSource,material,0,fullBoundary);
            const auto foreign=ElectricVofField::prepareElectrostaticSnapshot(*mesh,fullSource,fullCharge,material,fullBoundary);
            check(preparedField&&preparedField->advection(),"prepared Field foreign-domain test lacks its actual typed midpoint");
            const auto& a=*preparedField->advection();
            const auto offered=ElectricVofField::prepareOhmic(*mesh,a.geometry(),a,material,electric,ElectricVofControls{},
                ElectricVofField::CellGradientRecovery::Recover,&foreign);
            ElectricVofPostprocessContract::payload(*preparedField,offered);
        }
        rejected([&]{ChargeTransport::prepareClosed(*mesh,fullCharge,closedVolume);},"E_SELECTED_CHAIN foreign zero-valued active mask accepted");++guards;
        const auto otherDomain=ElectricVofSelectedDomainContract::selectedDomain(*mesh,!rz,true);
        const auto otherSource=Publication::prepareLocalPlicInventory(*mesh,otherDomain,start,pureStocks(*mesh,otherDomain,rz),{});
        const auto otherCharge=ElectricVofCharge::fromIntegrated(*mesh,otherSource,std::vector<double>(q.size()));
        rejected([&]{ChargeTransport::prepareClosed(*mesh,otherCharge,closedVolume);},"E_SELECTED_CHAIN foreign XY/RZ physical measure accepted");++guards;
        rejected([&]{ChargeTransport::prepareClosed(*mesh,charge,openVolume);},"E_SELECTED_CHAIN closed charge entry accepted actual open gross flow");++guards;
        check(snapshot(source)==saved&&q==charge.integratedCharge(),"E_SELECTED_CHAIN success/rejection changed accepted phase/Q");
    }
    check(requestedCompletions==8&&requestedGuards==32&&distinctEndpoint>0&&typedOmit==8,"requested-stage completion/source/law/control/typed-Omit evidence incomplete");
    check(stages==16&&guards==34&&zeroNetGross>0&&highOrder>0&&maxwellSamples>0&&wallFaces>0,
        "E_SELECTED_CHAIN required signed/gross/FCT/Maxwell witnesses incomplete");
    std::cout<<"ELECTRIC_VOF_SELECTED_CHAIN_PASS,half_full_stages,"<<stages<<",guards,"<<guards<<",zero_net_gross_charge_faces,"<<zeroNetGross
        <<",high_order_faces,"<<highOrder<<",Maxwell_samples,"<<maxwellSamples<<",internal_closed_faces,"<<wallFaces
        <<",XY_RZ,actual_2to1,signed_Q,P0_FCT,typed_Ohmic_Maxwell,original_gates"
        <<",requested_completions,"<<requestedCompletions<<",requested_guards,"<<requestedGuards
        <<",distinct_half_full_Q,"<<distinctEndpoint<<",typed_Omit,"<<typedOmit<<"\n";
}
} // namespace ElectricVofSelectedChainContract

namespace ElectricVofSelectedAmrContract {
using namespace ElectricVofContract;
using Transfer=I::GeometricVofMeshTransfer;
using Remap=LAMP::Electric::ElectricVofChargeMeshTransfer;
using Exact=LAMP::Geometry::Detail::FiniteVolume::Dyadic;
using Rectangle=LAMP::Geometry::PlanarCellCut::Rectangle;
bool fluid(const Cell& cell) {const auto p=cell.centerCoord();return !(p.x>1&&p.y<3);}
NativeLeafDomain selected(const QuadtreeMesh& mesh,bool rz,bool exchanged=false,int omit=-1) {
    FaceSegmentGraph graph;graph.build(mesh,[exchanged,omit](const QuadtreeMesh&,const Cell& cell){const auto p=cell.centerCoord();
        return (exchanged?!(p.x>1&&p.y>3):fluid(cell))&&cell.index()!=omit;});
    return NativeLeafDomain::build(mesh,graph,rz);
}
std::unique_ptr<QuadtreeMesh> nozzle() {
    DictionaryMesh d;d.Nx=d.Ny=2;d.minCoord={0,0};d.maxCoord={2,6};d.maxCellNum=512;d.maxInitialIterations=0;d.maxLevel_local={2};
    return std::make_unique<QuadtreeMesh>(d);
}
Rectangle bounds(const QuadtreeMesh& mesh,int ci) {
    const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);const auto c=cell.centerCoord(),h=cell.size();
    return {{c.x-h.x/2,c.y-h.y/2},{c.x+h.x/2,c.y+h.y/2}};
}
long double volume(const Rectangle& b,bool rz) {
    if(b.upper.x<=b.lower.x||b.upper.y<=b.lower.y)return 0;
    return (static_cast<long double>(b.upper.x)-b.lower.x)*(static_cast<long double>(b.upper.y)-b.lower.y)*
        (rz?static_cast<long double>(CartesianMeasure::twoPi)*(b.lower.x+b.upper.x)/2:1);
}
std::vector<Stock> stocks(const QuadtreeMesh& mesh,const NativeLeafDomain& d) {
    std::vector<Stock> out(mesh.cellIndexNum());
    for(int ci:d.cellIndexByRow()){auto liquid=bounds(mesh,ci),gas=liquid;liquid.upper.y=std::min(liquid.upper.y,4.234375);
        gas.lower.y=std::max(gas.lower.y,4.234375);out[ci]={volume(liquid,d.axisymmetric()),volume(gas,d.axisymmetric())};}
    return out;
}
void adapt(QuadtreeMesh& mesh,int level) {
    CellField<double> aim(mesh,0.),gradient(mesh,0.);
    for(int ci=0;ci<mesh.cellIndexNum();++ci)aim.value[ci]=fluid(mesh.cells(mesh.cellIndexToLabel[ci]))?level:0;
    check(mesh.mayAdaptToTargets(aim),"E_SELECTED_AMR actual refine/coarsen preflight failed");
    std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
    std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;std::vector<CellField<double>*> cells{&aim};
    const auto before=mesh.cellIndexNum();mesh.adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&aim);
    check(before!=mesh.cellIndexNum(),"E_SELECTED_AMR actual tree did not change");
}
void verify(const QuadtreeMesh& mesh,const Publication& source,const ElectricVofCharge& q,
    const QuadtreeMesh& target,const Transfer& t,const Remap& result,bool tiny,std::size_t& solids,std::size_t& active,std::size_t& pairs) {
    result.validateFor(mesh,source,q,target,t);const auto& candidate=result.chargeCandidate();const auto& next=candidate.integratedCharge();
    const auto clock=q.timeLayer();check(clock.position==Clock::Position::End&&sameClock(clock,candidate.timeLayer())&&
        sameClock(clock,t.candidate().timeLayer()),"E_SELECTED_AMR changed actual complete End clock");
    check(result.sourceCharge().sameInstance(q)&&result.sourceGeometry().identity()==source.identity()&&
        result.interfaceTransfer().candidate().identity()==t.candidate().identity(),"E_SELECTED_AMR receipt lost actual owners");
    const auto& overlap=t.nativeOverlaps();std::vector<long double> coveredSource(mesh.cellIndexNum()),coveredTarget(target.cellIndexNum());
    std::size_t intersections=0;
    // Independent all-pairs geometric intersections certify the selected active domain.
    for(int si:source.domain().cellIndexByRow())for(int ti:t.candidate().domain().cellIndexByRow()){
        const auto a=bounds(mesh,si),b=bounds(target,ti);
        const Rectangle intersection{{std::max(a.lower.x,b.lower.x),std::max(a.lower.y,b.lower.y)},
            {std::min(a.upper.x,b.upper.x),std::min(a.upper.y,b.upper.y)}};
        const auto cap=volume(intersection,source.domain().axisymmetric());if(cap==0)continue;std::size_t matches=0;
        for(const auto& piece:overlap.pieces())if(piece.sourceCell==si&&piece.targetCell==ti){
            bounded(piece.capacity-cap,cap,1e-12L,"E_SELECTED_AMR independent actual overlap measure");++matches;}
        check(matches==1,"E_SELECTED_AMR selected intersection lacks unique real overlap");coveredSource[si]+=cap;coveredTarget[ti]+=cap;++intersections;}
    check(intersections==overlap.pieces().size(),"E_SELECTED_AMR overlap contains nonphysical/solid pair");pairs+=intersections;
    const auto inspect=[&](const QuadtreeMesh& actual,const Publication& p,const auto& values,const auto& groups,const auto& covered){
        for(int ci=0;ci<actual.cellIndexNum();++ci){
            if(p.domain().cellRowByIndex()[ci]<0){check(groups[ci].empty()&&values[ci]==0&&covered[ci]==0&&
                p.plicCells()[ci].inventory[0]==0&&p.plicCells()[ci].inventory[1]==0,
                "E_SELECTED_AMR inactive cell has overlap/phase/charge stock");++solids;continue;}
            const auto cap=volume(bounds(actual,ci),p.domain().axisymmetric());check(!groups[ci].empty()&&std::isfinite(values[ci]),
                "E_SELECTED_AMR active cell missing coverage or finite Q");
            bounded(covered[ci]-cap,cap,1e-12L,"E_SELECTED_AMR independent per-active-leaf coverage");++active;}
    };
    inspect(mesh,source,q.integratedCharge(),overlap.bySource(),coveredSource);
    inspect(target,t.candidate(),next,overlap.byTarget(),coveredTarget);
    Exact before,after;bool positive=false,negative=false;
    for(double v:q.integratedCharge())before.add(v);
    for(double v:next){after.add(v);positive|=v>0;negative|=v<0;}
    check(before.compare(Exact{})==0&&after.compare(before)==0&&positive&&negative,
        "E_SELECTED_AMR signed zero-net stock or a charge sign was lost");
    const auto& r=result.report();check(r.relativeChargeResidual<=gate&&r.chargeResidual==0&&r.roundingResidual==0&&
        r.pieces==overlap.pieces().size(),"E_SELECTED_AMR original signed exact/global ledger gate");
    check(t.report().relativeGlobalPhaseResidual[0]<=1e-12L&&t.report().relativeGlobalPhaseResidual[1]<=1e-12L&&
        t.report().maximumRelativeCoverageResidual<=1e-12L,"E_SELECTED_AMR original independent Interface phase/coverage gates");
    if(tiny)ElectricVofMovingContract::tinyRemapOracle(result);
}
void run() {
    std::size_t remaps=0,adaptations=0,guards=0,solids=0,active=0,pairs=0,refinements=0,coarsenings=0;
    for(bool rz:{false,true}){
        auto root=nozzle();const auto d=selected(*root,rz);const Clock end{7,.2,.25,Clock::Position::End};const auto inventory=stocks(*root,d);
        check(root->cellIndexNum()==4&&d.cellIndexByRow().size()==3,"E_SELECTED_AMR actual selected L-shape/solid fixture missing");
        const auto source=Publication::prepareLocalPlicInventory(*root,d,end,inventory,{});const auto saved=snapshot(source);
        for(bool tiny:{false,true}){
            std::vector<double> initial(root->cellIndexNum());const auto unit=tiny?7*std::numeric_limits<double>::denorm_min():.125;
            initial[d.cellIndexByRow()[0]]=unit;initial[d.cellIndexByRow()[1]]=-unit;
            const auto charge=ElectricVofCharge::fromIntegrated(*root,source,initial);auto current=source;auto currentQ=charge;auto mesh=root->clone();
            const auto unchanged=Transfer::prepare(*root,source,*mesh,selected(*mesh,rz),end);
            const auto same=Remap::prepare(*root,source,charge,*mesh,unchanged);verify(*root,source,charge,*mesh,unchanged,same,tiny,solids,active,pairs);++remaps;
            check(std::memcmp(same.chargeCandidate().integratedCharge().data(),initial.data(),initial.size()*sizeof(double))==0,
                "E_SELECTED_AMR unchanged transfer altered Q bits");current=unchanged.candidate();currentQ=same.chargeCandidate();
            bool graded=false;std::size_t up=0,down=0;
            for(int level:{1,0,1,2,1,0}){
                auto target=mesh->clone();adapt(*target,level);++adaptations;const auto td=selected(*target,rz);graded|=!td.uniformLeaves();
                const auto oldStock=snapshot(current);const auto oldQ=currentQ.integratedCharge();
                const auto transfer=Transfer::prepare(*mesh,current,*target,td,end);const auto moved=Remap::prepare(*mesh,current,currentQ,*target,transfer);
                verify(*mesh,current,currentQ,*target,transfer,moved,tiny,solids,active,pairs);++remaps;
                up+=moved.report().refinedDonors>0;down+=moved.report().coarsenedTargets>0;
                check(snapshot(current)==oldStock&&currentQ.integratedCharge()==oldQ,"E_SELECTED_AMR candidate mutated predecessor phase/Q");
                current=transfer.candidate();currentQ=moved.chargeCandidate();mesh=std::move(target);
            }
            check(up>=3&&down>=3&&graded,"E_SELECTED_AMR repeated actual selected refinement/coarsening/2:1 missing");refinements+=up;coarsenings+=down;
            for(int ci=0;ci<root->cellIndexNum();++ci){const auto point=root->cells(root->cellIndexToLabel[ci]).centerCoord();
                const auto finalIndex=mesh->cells(mesh->findCellLabel(point,0)).index();
                check(currentQ.integratedCharge()[finalIndex]==initial[ci],"E_SELECTED_AMR repeated signed charge roundtrip changed root stock");}
            check(snapshot(source)==saved&&charge.integratedCharge()==initial,"E_SELECTED_AMR repeated transfer changed original accepted state");
        }
        auto fine=root->clone();adapt(*fine,1);const auto fd=selected(*fine,rz);const auto t=Transfer::prepare(*root,source,*fine,fd,end);
        std::vector<double> q(root->cellIndexNum());q[d.cellIndexByRow()[0]]=.125;q[d.cellIndexByRow()[1]]=-.125;
        const auto charge=ElectricVofCharge::fromIntegrated(*root,source,q);const auto result=Remap::prepare(*root,source,charge,*fine,t);
        const auto foreignQ=ElectricVofCharge::fromIntegrated(*root,source,q);
        rejected([&]{result.validateFor(*root,source,foreignQ,*fine,t);},"E_SELECTED_AMR foreign same-value source Q owner accepted");++guards;
        const auto foreignSource=Publication::prepareLocalPlicInventory(*root,d,end,inventory,{});
        rejected([&]{Remap::prepare(*root,foreignSource,charge,*fine,t);},"E_SELECTED_AMR foreign source geometry owner accepted");++guards;
        const auto foreignTransfer=Transfer::prepare(*root,source,*fine,fd,end);
        rejected([&]{result.validateFor(*root,source,charge,*fine,foreignTransfer);},"E_SELECTED_AMR foreign same-layout target transfer owner accepted");++guards;
        rejected([&]{result.validateFor(*root,source,charge,*root,t);},"E_SELECTED_AMR foreign target native layout accepted");++guards;
        rejected([&]{Remap::prepare(*fine,source,charge,*fine,t);},"E_SELECTED_AMR foreign source native layout accepted");++guards;
        auto wrongClock=end;wrongClock.begin=.21;const auto wrongGeometry=Publication::prepareLocalPlicInventory(*root,d,wrongClock,inventory,{});
        const auto wrongQ=ElectricVofCharge::fromIntegrated(*root,wrongGeometry,q);
        rejected([&]{Remap::prepare(*root,source,wrongQ,*fine,t);},"E_SELECTED_AMR same-time wrong complete End interval accepted");++guards;
        const Clock initial{0,end.end,end.end,Clock::Position::Initial};const auto initialGeometry=Publication::prepareLocalPlicInventory(*root,d,initial,inventory,{});
        const auto initialQ=ElectricVofCharge::fromIntegrated(*root,initialGeometry,q);
        rejected([&]{Remap::prepare(*root,initialGeometry,initialQ,*fine,t);},"E_SELECTED_AMR Initial charge accepted as actual End");++guards;
        const auto fullDomain=domainFor(*root,rz);const auto fullGeometry=Publication::prepareLocalPlicInventory(*root,fullDomain,end,stocks(*root,fullDomain),{});
        const auto fullZero=ElectricVofCharge::fromIntegrated(*root,fullGeometry,std::vector<double>(q.size()));
        rejected([&]{Remap::prepare(*root,source,fullZero,*fine,t);},"E_SELECTED_AMR changed source Q mask bypassed with zero charge");++guards;
        int solid=-1;for(int ci=0;ci<root->cellIndexNum();++ci)if(d.cellRowByIndex()[ci]<0)solid=ci;check(solid>=0,"E_SELECTED_AMR missing solid witness");
        auto bad=q;bad[solid]=std::numeric_limits<double>::denorm_min();
        rejected([&]{ElectricVofCharge::fromIntegrated(*root,source,bad);},"E_SELECTED_AMR inactive source Q stock accepted");++guards;
        const auto fullBad=ElectricVofCharge::fromIntegrated(*root,fullGeometry,bad);
        rejected([&]{Remap::prepare(*root,source,fullBad,*fine,t);},"E_SELECTED_AMR foreign full-mask source leaked solid Q");++guards;
        const auto swapped=selected(*root,rz,true);long double actualVolume=0,swappedVolume=0;
        for(double v:d.cellVolumes())actualVolume+=v;for(double v:swapped.cellVolumes())swappedVolume+=v;
        bounded(swappedVolume-actualVolume,actualVolume,1e-12L,"E_SELECTED_AMR changed-mask equal-volume witness differs");
        check(swapped.cellRowByIndex()!=d.cellRowByIndex(),"E_SELECTED_AMR negative mask unchanged");
        rejected([&]{Transfer::prepare(*root,source,*root,swapped,end);},"E_SELECTED_AMR same-volume exchanged physical mask accepted");++guards;
        rejected([&]{Transfer::prepare(*root,source,*fine,selected(*fine,rz,false,fd.cellIndexByRow().front()),end);},"E_SELECTED_AMR missing active fine child accepted");++guards;
        const auto otherMeasure=selected(*root,!rz);const auto otherGeometry=Publication::prepareLocalPlicInventory(*root,otherMeasure,end,stocks(*root,otherMeasure),{});
        const auto otherQ=ElectricVofCharge::fromIntegrated(*root,otherGeometry,q);
        rejected([&]{Remap::prepare(*root,source,otherQ,*fine,t);},"E_SELECTED_AMR foreign source Q physical measure accepted");++guards;
        ElectricVofControls relaxed;relaxed.relativeConservationTolerance=1e-8;
        rejected([&]{Remap::prepare(*root,source,charge,*fine,t,relaxed);},"E_SELECTED_AMR relaxed original electric gate accepted");++guards;
        check(snapshot(source)==saved&&charge.integratedCharge()==q,"E_SELECTED_AMR negative controls mutated accepted phase/Q");
        std::cout<<"ELECTRIC_VOF_SELECTED_AMR_METRIC,"<<(rz?"RZ":"XY")<<",actual_same_End_step,"<<end.step<<",End_time,"<<end.time()
            <<",source_fluid_leaves,"<<d.cellIndexByRow().size()<<",solid_leaves,"<<root->cellIndexNum()-d.cellIndexByRow().size()<<'\n';
    }
    check(remaps==28&&adaptations==24&&refinements==12&&coarsenings==12&&guards==28&&solids>0&&active>0&&pairs>0,
        "E_SELECTED_AMR bounded signed/zero-net/tiny/owner coverage incomplete");
    std::cout<<"ELECTRIC_VOF_SELECTED_AMR_PASS,remaps,"<<remaps<<",actual_adaptations,"<<adaptations<<",refinements,"<<refinements
        <<",coarsenings,"<<coarsenings<<",guards,"<<guards<<",solid_checks,"<<solids<<",active_checks,"<<active<<",independent_overlaps,"<<pairs
        <<",XY_RZ,selected_L_shape,complete_same_End,signed_zero_net,odd_subnormal,independent_integer_CDF,original_gates,no_full_nozzle_acceptance\n";
}
} // namespace ElectricVofSelectedAmrContract

namespace ElectricVofSecondOrderRemapContract {
using namespace ElectricVofContract;
using namespace ElectricVofMovingContract;
namespace Exact=LAMP::Geometry::Detail::FiniteVolume;
using Scheme=LAMP::Electric::ElectricVofChargeAdvectionScheme;
using Box=LAMP::Geometry::PlanarCellCut::Rectangle;
std::size_t scenes=0,affinePieces=0,bipolarDonors=0,finiteDonors=0,overflowSums=0;
ElectricVofControls controls(){ElectricVofControls c;c.chargeAdvection=Scheme::LimitedLinearFct;return c;}
Box boxFor(const QuadtreeMesh& mesh,int ci){const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);const auto x=cell.centerCoord(),h=cell.size();
    return {{x.x-h.x/2,x.y-h.y/2},{x.x+h.x/2,x.y+h.y/2}};}
// Independent analytic primitive of (a+b*r+c*z) dV; no production centroid,
// polynomial reconstruction, overlap capacity or source partition is reused.
long double analytic(const Box& box,bool rz,long double a,long double b,long double c){
    const long double x0=box.lower.x,x1=box.upper.x,z0=box.lower.y,z1=box.upper.y,zm=(z0+z1)/2;
    if(rz)return CartesianMeasure::twoPi*(z1-z0)*((a+c*zm)*(x1*x1-x0*x0)/2+b*(x1*x1*x1-x0*x0*x0)/3);
    return (x1-x0)*(z1-z0)*(a+b*(x0+x1)/2+c*zm);
}
void ledger(const ChargeRemap& remap){const auto& overlap=remap.interfaceTransfer().nativeOverlaps();
    const auto& old=remap.sourceCharge().integratedCharge();const auto& next=remap.chargeCandidate().integratedCharge();
    Exact::Dyadic before,after;for(double q:old)before.add(q);for(double q:next)after.add(q);after.add(before,-1);
    check(Exact::roundDouble(after).value==remap.report().chargeResidual&&remap.report().chargeResidual==remap.report().roundingResidual,
        "REMAP_P1 global difference is not declared final rounding");
    check(remap.report().relativeChargeResidual<=gate&&remap.report().maximumDonorRelativeChargeResidual<=gate,
        "REMAP_P1 relaxed original charge gate");
    for(std::size_t si=0;si<old.size();++si){const auto& group=overlap.bySource()[si];if(group.size()<2)continue;
        Exact::Dyadic own;for(auto pi:group)own.add(next[overlap.pieces()[pi].targetCell]);own.add(old[si],-1);
        check(own.compare(Exact::Dyadic{})==0,"REMAP_P1 lost a donor's exact signed stock");
    }
    check(remap.report().prolongationRoundingResidual==0,"REMAP_P1 weakened exact donor conservation into a rounding allowance");++scenes;
}
void run(){
    for(bool rz:{false,true}){
        auto mesh=ElectricVofSecondOrderChargeContract::uniformMesh(8);const auto domain=domainFor(*mesh,rz);
        const Clock end{1,.2,.25,Clock::Position::End};const auto stocks=ElectricVofSecondOrderChargeContract::gasStocks(domain);
        const auto geometry=Publication::prepareLocalPlicInventory(*mesh,domain,end,stocks,{});const auto geometryBefore=snapshot(geometry);
        auto fine=mesh->clone();DictionaryMesh::RefineBox refineBox;refineBox.min={.375,.375};refineBox.max={.625,.625};refineBox.level=1;
        fine->applyRefineBoxes({refineBox});check(fine->cellIndexNum()>mesh->cellIndexNum(),"REMAP_P1 no actual refined leaves");
        const auto transfer=VofTransfer::prepare(*mesh,geometry,*fine,domainFor(*fine,rz),end);const auto& overlap=transfer.nativeOverlaps();
        const auto backTransfer=VofTransfer::prepare(*fine,transfer.candidate(),*mesh,domain,end);
        for(int polynomial=0;polynomial<2;++polynomial){const long double a=polynomial?-.4375L:.5L,b=polynomial?0:.375L,c=polynomial?1:.25L;
            std::vector<double> q(mesh->cellIndexNum());for(int ci=0;ci<mesh->cellIndexNum();++ci)q[ci]=static_cast<double>(analytic(boxFor(*mesh,ci),rz,a,b,c));
            const auto source=ElectricVofCharge::fromIntegrated(*mesh,geometry,q);const auto high=ChargeRemap::prepare(*mesh,geometry,source,*fine,transfer,controls());
            const auto low=ChargeRemap::prepare(*mesh,geometry,source,*fine,transfer);ledger(high);ledger(low);
            check(high.report().limitedLinearRemap&&high.report().highOrderRefinedDonors>0&&!low.report().limitedLinearRemap,
                "REMAP_P1 scheme did not select actual high-order refined donors");
            long double highError=0,lowError=0,scale=0;
            for(const auto& p:overlap.pieces())if(p.kind==LAMP::Mesh::NativeLeafOverlapPiece::Kind::Refined){const auto exact=analytic(p.rectangle,rz,a,b,c);
                highError+=std::abs(static_cast<long double>(high.chargeCandidate().integratedCharge()[p.targetCell])-exact);
                lowError+=std::abs(static_cast<long double>(low.chargeCandidate().integratedCharge()[p.targetCell])-exact);scale+=std::abs(exact);++affinePieces;}
            check(highError<=1e-12L*scale&&lowError>1e-4L*scale,"REMAP_P1 actual affine overlap integral is not more accurate than P0");
            if(polynomial)for(std::size_t si=0;si<q.size();++si)if(q[si]==0&&overlap.bySource()[si].size()>1){bool positive=false,negative=false;
                for(auto pi:overlap.bySource()[si]){const auto value=high.chargeCandidate().integratedCharge()[overlap.pieces()[pi].targetCell];positive=positive||value>0;negative=negative||value<0;}
                check(positive&&negative,"REMAP_P1 erased zero-net donor's physical bipolar children");++bipolarDonors;}
            check(q==source.integratedCharge(),"REMAP_P1 changed accepted signed source");
            std::cout<<"ELECTRIC_VOF_REMAP_P1_AFFINE,"<<(rz?"RZ":"XY")<<",polynomial,"<<polynomial<<",L1_P1,"<<highError/scale<<",L1_P0,"<<lowError/scale<<'\n';
        }
        // Sharp signed data: local density bounds are tested independently of
        // the production linear polynomial and its exact inventory budgets.
        std::vector<double> jump(mesh->cellIndexNum()),lo(jump.size()),hi(jump.size());
        for(std::size_t ci=0;ci<jump.size();++ci){const auto x=mesh->cells(mesh->cellIndexToLabel[ci]).centerCoord();
            jump[ci]=domain.cellVolumes()[domain.cellRowByIndex()[ci]]*(x.x+x.y<1?-.75:1.25);lo[ci]=hi[ci]=jump[ci]/domain.cellVolumes()[domain.cellRowByIndex()[ci]];}
        for(const auto f:domain.faces(*mesh))if(f.negativeCell>=0&&f.positiveCell>=0){const int i=f.negativeCell,j=f.positiveCell;
            const auto a=jump[i]/domain.cellVolumes()[domain.cellRowByIndex()[i]],b=jump[j]/domain.cellVolumes()[domain.cellRowByIndex()[j]];
            lo[i]=std::min(lo[i],b);hi[i]=std::max(hi[i],b);lo[j]=std::min(lo[j],a);hi[j]=std::max(hi[j],a);}
        const auto jumpQ=ElectricVofCharge::fromIntegrated(*mesh,geometry,jump);const auto jumped=ChargeRemap::prepare(*mesh,geometry,jumpQ,*fine,transfer,controls());ledger(jumped);
        for(const auto& p:overlap.pieces()){const auto rho=jumped.chargeCandidate().integratedCharge()[p.targetCell]/transfer.candidate().domain().cellVolumes()[transfer.candidate().domain().cellRowByIndex()[p.targetCell]];
            check(rho>=lo[p.sourceCell]-1e-12&&rho<=hi[p.sourceCell]+1e-12,"REMAP_P1 introduced signed density extremum");}
        // Explicit odd binary64 quanta keep the qualified original CDF.
        std::vector<double> tiny(jump.size());std::size_t chosen=0;const auto unit=std::numeric_limits<double>::denorm_min();
        int zeroDonor=-1,unchanged=-1;
        for(std::size_t si=0;si<tiny.size();++si){if(overlap.bySource()[si].size()>1){tiny[si]=(chosen%2?-1:1)*(chosen%3==0?1:chosen%3==1?3:7)*unit;++chosen;
                const auto center=mesh->cells(mesh->cellIndexToLabel[si]).centerCoord();if(center.y==.4375)zeroDonor=static_cast<int>(si);}
            else unchanged=static_cast<int>(si);}
        check(chosen>0&&zeroDonor>=0&&unchanged>=0,"REMAP_P1 tiny/bipolar/unchanged fixture absent");
        // A local subnormal donor must not propagate its reconstruction
        // qualification through the entire refined patch via face order.
        std::vector<double> localTiny(jump.size());for(int ci=0;ci<mesh->cellIndexNum();++ci)
            localTiny[ci]=static_cast<double>(analytic(boxFor(*mesh,ci),rz,.5L,.375L,.25L));
        localTiny[zeroDonor]=unit;const auto localTinyQ=ElectricVofCharge::fromIntegrated(*mesh,geometry,localTiny);
        const auto localTinyFine=ChargeRemap::prepare(*mesh,geometry,localTinyQ,*fine,transfer,controls());ledger(localTinyFine);
        check(localTinyFine.report().finiteEligibilityP0Donors==1&&
            localTinyFine.report().highOrderRefinedDonors+1==localTinyFine.report().refinedDonors,
            "REMAP_P1 propagated local tiny eligibility into normal refined donors");
        // Oblique bipolar reconstruction exercises nonsymmetric RZ weights
        // and finite signed quantization beyond pairwise z-only symmetry.
        const auto volumeCenter=[&](int ci){const auto box=boxFor(*mesh,ci);const long double x0=box.lower.x,x1=box.upper.x;
            return std::array<long double,2>{rz?2*(x1*x1*x1-x0*x0*x0)/(3*(x1*x1-x0*x0)):(x0+x1)/2,
                (static_cast<long double>(box.lower.y)+box.upper.y)/2};};
        const auto center=volumeCenter(zeroDonor);std::vector<double> oblique(jump.size());
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto point=volumeCenter(ci);
            oblique[ci]=static_cast<double>(domain.cellVolumes()[domain.cellRowByIndex()[ci]]*(.73L*(point[0]-center[0])+.29L*(point[1]-center[1])));}
        oblique[zeroDonor]=0;const auto obliqueQ=ElectricVofCharge::fromIntegrated(*mesh,geometry,oblique);
        const auto obliqueFine=ChargeRemap::prepare(*mesh,geometry,obliqueQ,*fine,transfer,controls());ledger(obliqueFine);
        bool positiveOblique=false,negativeOblique=false;for(auto pi:overlap.bySource()[zeroDonor]){
            const auto value=obliqueFine.chargeCandidate().integratedCharge()[overlap.pieces()[pi].targetCell];
            positiveOblique=positiveOblique||value>0;negativeOblique=negativeOblique||value<0;}
        check(positiveOblique&&negativeOblique,"REMAP_P1 erased oblique zero-net donor's bipolar structure");++bipolarDonors;
        const auto tinyQ=ElectricVofCharge::fromIntegrated(*mesh,geometry,tiny);const auto tinyFine=ChargeRemap::prepare(*mesh,geometry,tinyQ,*fine,transfer,controls());
        ledger(tinyFine);tinyRemapOracle(tinyFine);finiteDonors+=tinyFine.report().finiteEligibilityP0Donors;
        check(tinyFine.report().finiteEligibilityP0Donors>0,"REMAP_P1 subnormal finite eligibility not reported");
        const auto tinyBack=ChargeRemap::prepare(*fine,transfer.candidate(),tinyFine.chargeCandidate(),*mesh,backTransfer,controls());
        ledger(tinyBack);check(tinyBack.chargeCandidate().integratedCharge()==tiny,"REMAP_P1 subnormal refine/coarsen cycle lost signed stock");
        // Tiny net charge beside order-one signed gradients cannot generally
        // coexist with binary64 bipolar children and an exact source quantum.
        std::vector<double> cancellation(tiny.size());for(int ci=0;ci<mesh->cellIndexNum();++ci)cancellation[ci]=static_cast<double>(analytic(boxFor(*mesh,ci),rz,-.4375L,0,1));
        cancellation[zeroDonor]=1e-300;cancellation[unchanged]=-0.;
        const auto cancellationQ=ElectricVofCharge::fromIntegrated(*mesh,geometry,cancellation);
        const auto cancellationFine=ChargeRemap::prepare(*mesh,geometry,cancellationQ,*fine,transfer,controls());ledger(cancellationFine);
        check(cancellationFine.report().finiteEligibilityP0Donors>0,"REMAP_P1 strong cancellation did not report finite eligibility");finiteDonors+=cancellationFine.report().finiteEligibilityP0Donors;
        const auto unchangedPiece=overlap.bySource()[unchanged].front();const auto unchangedValue=cancellationFine.chargeCandidate().integratedCharge()[overlap.pieces()[unchangedPiece].targetCell];
        check(std::signbit(unchangedValue)&&unchangedValue==0,"REMAP_P1 unchanged negative zero changed bits");
        // Full signed exact restriction must cancel before any native sum,
        // even when a same-sign prefix would overflow.
        int parent=-1;for(std::size_t ti=0;ti<backTransfer.nativeOverlaps().byTarget().size();++ti)
            if(backTransfer.nativeOverlaps().byTarget()[ti].size()==4){parent=static_cast<int>(ti);break;}
        check(parent>=0,"REMAP_P1 restriction lacks four actual children");const auto& children=backTransfer.nativeOverlaps().byTarget()[parent];
        for(bool residual:{false,true}){std::vector<double> charge(fine->cellIndexNum());const double maximum=std::numeric_limits<double>::max();
            const std::array<double,4> values=residual?std::array<double,4>{maximum,-maximum,unit,-0.}:std::array<double,4>{maximum,maximum,-maximum,-maximum};
            for(std::size_t k=0;k<4;++k)charge[backTransfer.nativeOverlaps().pieces()[children[k]].sourceCell]=values[k];
            const auto source=ElectricVofCharge::fromIntegrated(*fine,transfer.candidate(),charge);
            const auto restricted=ChargeRemap::prepare(*fine,transfer.candidate(),source,*mesh,backTransfer,controls());ledger(restricted);
            check(restricted.chargeCandidate().integratedCharge()[parent]==(residual?unit:0)&&restricted.report().roundingResidual==0,
                "REMAP_P1 exact restriction lost overflow-cancelling prefix or tiny remainder");++overflowSums;
        }
        const auto foreign=ElectricVofCharge::fromIntegrated(*mesh,geometry,jump);
        rejected([&]{jumped.validateFor(*mesh,geometry,foreign,*fine,transfer);},"REMAP_P1 foreign source receipt accepted");
        check(snapshot(geometry)==geometryBefore&&jumpQ.integratedCharge()==jump,"REMAP_P1 success/rejection mutated accepted source");
    }
    check(scenes>=20&&affinePieces>0&&bipolarDonors>0&&finiteDonors>0&&overflowSums==4,"REMAP_P1 required witnesses incomplete");
    std::cout<<"ELECTRIC_VOF_REMAP_P1_PASS,scenes,"<<scenes<<",analytic_pieces,"<<affinePieces<<",zero_net_bipolar_donors,"<<bipolarDonors
        <<",finite_eligibility_P0_donors,"<<finiteDonors<<",overflow_cancellation_sums,"<<overflowSums
        <<",XY_RZ,actual_overlap,original_donor_and_global_gates,no_floor_no_clip_no_redistribution\n";
}
} // namespace ElectricVofSecondOrderRemapContract


namespace ElectricVofWeakFieldGaussContract {
using Publication=LAMP::Interface::GeometricVofPublication;
using Clock=LAMP::TimeLayer;
namespace Exact=LAMP::Geometry::Detail::FiniteVolume;
void run(ElectricVofControls::Backend backend=ElectricVofControls::Backend::Gmres){
    ElectricVofControls controls;controls.backend=backend;
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={5,16};input.Nx=40;input.Ny=80;
    input.maxInitialIterations=0;input.maxCellNum=4096;input.maxLevel_local={0};
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;
    graph.build(mesh,[](const QuadtreeMesh&,const Cell& cell){const auto x=cell.centerCoord();return x.x<1||x.y>3;});
    const auto domain=NativeLeafDomain::build(mesh,graph,true);
    // Initial pure phases use the Interface owner to issue its actual
    // capacity and independent stocks; Domain's rounded metric is not stock.
    std::vector<double> alpha(mesh.cellIndexNum());
    for(int ci:domain.cellIndexByRow())alpha[ci]=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord().y<3?1:0;
    const Clock clock{1,0,5e-6,Clock::Position::End};
    const auto p=Publication::prepareLocalPlic(mesh,domain,clock,alpha,{});
    std::vector<ElectricVofBoundary> boundary(domain.nativeFaceCount());std::size_t conductor=0,axis=0;
    using Type=ElectricVofBoundary::Type;
    for(int row=0;row<domain.nativeFaceCount();++row){const auto face=domain.face(mesh,row);
        if((face.negativeCell>=0)==(face.positiveCell>=0))continue;
        if(face.area==0){boundary[row]={Type::HomogeneousNeumann,0};++axis;}
        else if((!face.horizontal&&face.center.x==1&&face.center.y<3)||
            (face.horizontal&&face.center.y==3&&face.center.x>1)){
            boundary[row]={Type::Dirichlet,6};++conductor;
        }else if(face.horizontal&&face.center.y==16)boundary[row]={Type::Dirichlet,0};
        else boundary[row]={Type::HomogeneousNeumann,0};
    }
    check(conductor>0&&axis>0&&domain.cellIndexByRow().size()==2720,"WEAK_GAUSS selected electrode/axis fixture missing");
    const ElectricVofMaterials material{4.8972277725663718,.4897227772566372,12825.409328871216,0};
    double peak=0,positiveDuration=0;std::size_t weak=0,nonzeroCorrections=0,subUlpCorrections=0,checkedFaces=0;
    for(int scene=0;scene<3;++scene){
        const bool ohmic=scene==2;
        std::vector<double> q(mesh.cellIndexNum());
        if(scene>0)for(int ci:domain.cellIndexByRow())if(mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord().y<1)
            q[ci]=(ci%2?-1:1)*1e-16;
        const auto before=q;const auto source=ElectricVofCharge::fromIntegrated(mesh,p,q);
        const Clock actualClock=ohmic?Clock{2,5e-6,1.075e-5,Clock::Position::Midpoint}:clock;
        const auto stage=ohmic?Publication::prepareLocalPlic(mesh,domain,actualClock,alpha,{}):p;
        const auto field=ohmic?ElectricVofField::prepareOhmic(mesh,stage,source,actualClock,material,boundary,controls):
            ElectricVofField::prepareElectrostaticSnapshot(mesh,p,source,material,boundary,controls);
        field.validateFor(mesh,stage,source,actualClock);
        check(source.integratedCharge()==before,"WEAK_GAUSS rewrote accepted signed Q");
        if(ohmic){positiveDuration=field.effectiveDt();check(positiveDuration==actualClock.time()-clock.time()&&positiveDuration>0&&
            !field.isElectrostaticSnapshot()&&!field.chargeCandidate().sameInstance(source),"WEAK_GAUSS actual Ohmic duration/owner missing");}
        else check(field.chargeCandidate().sameInstance(source)&&field.effectiveDt()==0&&field.isElectrostaticSnapshot(),
            "WEAK_GAUSS snapshot advanced signed Q");
        const auto& terms=field.potentialTerms();check(terms.size()==q.size(),"WEAK_GAUSS public potential terms lost dense native indexing");
        for(std::size_t ci=0;ci<terms.size();++ci){Exact::Dyadic sum;
            for(double component:terms[ci])sum.add(component);
            check(Exact::roundDouble(sum).value==field.potential()[ci],"WEAK_GAUSS displayed potential is not nearest same-owner term sum");
            if(domain.cellRowByIndex()[ci]<0)check(terms[ci]==std::array<double,3>{0,0,0},"WEAK_GAUSS inactive potential terms not zero");
            for(int correction=1;correction<3;++correction)if(terms[ci][correction]!=0){++nonzeroCorrections;
                Exact::Dyadic own(terms[ci][0]);own.add(terms[ci][correction]);
                if(Exact::roundDouble(own).value==terms[ci][0])++subUlpCorrections;}
        }
        std::vector<Exact::Dyadic> divergence(q.size()),gross(q.size()),transferred(q.size());
        for(const auto& face:field.faces()){
            Exact::Dyadic h(face.gradientConstant);
            for(const auto& coefficient:face.gradientCoefficients)
                for(double component:terms[coefficient.nativeCell])h.add(Exact::Dyadic::product({coefficient.value,component}));
            check(Exact::roundDouble(h).value==face.integratedField,"WEAK_GAUSS H does not consume the public same-owner potential terms");++checkedFaces;
            for(int ci:{face.negativeCell,face.positiveCell})if(ci>=0){const int sign=ci==face.negativeCell?1:-1;
                divergence[ci].add(face.displacement,sign);gross[ci].add(std::abs(face.displacement));
                transferred[ci].add(Exact::Dyadic::product({field.effectiveDt(),face.current}),sign);}
        }
        const auto& next=field.chargeCandidate().integratedCharge();
        for(int ci:domain.cellIndexByRow()){
            auto expected=Exact::Dyadic(q[ci]);expected.add(transferred[ci],-1);
            const double expectedQ=transferred[ci].compare(Exact::Dyadic{})==0?q[ci]:Exact::roundDouble(expected).value;
            check(next[ci]==expectedQ,"WEAK_GAUSS independent actual signed Ohmic update changed");
            auto residual=divergence[ci];residual.add(next[ci],-1);gross[ci].add(std::abs(next[ci]));
            check(residual.absolute().compare(gross[ci].multiplied(1e-10))<=0,"WEAK_GAUSS actual shared-face local physical bound failed");
            const auto scale=gross[ci].rounded();const double relative=scale==0?0:static_cast<double>(residual.absolute().rounded()/scale);peak=std::max(peak,relative);
            const auto x=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord();
            if(x.y<.5&&x.x<.5){check(field.potential()[ci]>5.999&&field.potential()[ci]<6,
                "WEAK_GAUSS fixture lacks weak physical field below common six-volt/code potential");++weak;}
        }
        ElectricVofContract::backendReceipt(field,backend);
        check(field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=1e-10,
            "WEAK_GAUSS retained linear/Gauss gates changed");
    }
    check(weak>0&&positiveDuration>0&&nonzeroCorrections>0&&subUlpCorrections>0,"WEAK_GAUSS required actual refinement/duration witnesses absent");
    std::cout<<"ELECTRIC_VOF_WEAK_FIELD_GAUSS_PASS,scenes,3,snapshots,2,positive_duration,"<<positiveDuration
        <<",active_cells,2720,weak_cells,"<<weak<<",nonzero_corrections,"<<nonzeroCorrections<<",sub_ulp_corrections,"<<subUlpCorrections
        <<",same_potential_faces,"<<checkedFaces<<",maximum_relative_Gauss,"<<peak
        <<",same_fixed_operator,2_defect_passes,original_physical_bound,1e-10,no_accepted_Q_change\n";
}
} // namespace ElectricVofWeakFieldGaussContract


namespace ElectricVofOhmicQuantumContract {
namespace Exact=LAMP::Geometry::Detail::FiniteVolume;
using D=Exact::Dyadic;
using Clock=LAMP::TimeLayer;
using Publication=LAMP::Interface::GeometricVofPublication;
constexpr double eta=std::numeric_limits<double>::denorm_min();
bool equal(const D& a,const D& b){return a.compare(b)==0;}
double fromBits(std::uint64_t bits){double value=0;std::memcpy(&value,&bits,sizeof value);return value;}
D quanta(const D& value){return value.multiplied(std::scalbn(1.,1022)).multiplied(std::scalbn(1.,52));}
double relative(const D& residual,const D& scale){
    if(equal(scale,D{}))return 0;
    auto a=residual.absolute(),b=scale;
    const unsigned count=scale.compare(D(eta))<0?2u:scale.compare(D(std::numeric_limits<double>::min()))<0?1u:0u;
    for(unsigned k=0;k<count;++k){a=a.multiplied(std::scalbn(1.,1022));b=b.multiplied(std::scalbn(1.,1022));}
    int exponent=0;std::frexp(b.rounded(),&exponent);const auto factor=std::scalbn(1.,1-exponent);
    a=a.multiplied(factor);b=b.multiplied(factor);
    return static_cast<double>(a.rounded()/b.rounded());
}
void arithmetic(){
    const auto quarterEta=D::product({eta,.25});
    check(relative(quarterEta,D(.25))==eta,"OHMIC_ETA representable eta ratio hidden by denominator normalization");
    const double independentRatio=static_cast<double>(quarterEta.multiplied(std::scalbn(1.,1022)).rounded())/
        std::scalbn(1e-300,1022);
    check(relative(quarterEta,D(1e-300))==independentRatio&&independentRatio>1e-24,
        "OHMIC_ETA normal scale hid sub-eta raw relative");
    for(int numerator:{-7,-5,-3,-1,0,1,3,5,7}){
        const auto raw=D::product({eta,static_cast<double>(numerator),.5});
        const int expected=numerator==0?0:(numerator<0?-1:1)*((std::abs(numerator)+1)/2/2*2);
        check(equal(raw.roundedBinary64Quantum(),D::product({eta,static_cast<double>(expected)})),"OHMIC_ETA nearest-even half quantum");
    }
    const auto wide=D::product({1+std::scalbn(1.,-52),1+std::scalbn(1.,-52)});
    check(equal(wide.roundedBinary64Quantum(),wide)&&!equal(wide,D(Exact::roundDouble(wide).value)),"OHMIC_ETA lost eta-integer 106-bit product");
    auto huge=D::product({std::numeric_limits<double>::max(),2});
    check(equal(huge.roundedBinary64Quantum(),huge),"OHMIC_ETA truncated super-double transfer");
    auto canceled=huge.roundedBinary64Quantum();canceled.add(huge,-1);canceled.add(eta);
    check(equal(canceled,D(eta)),"OHMIC_ETA exact super-double cancellation changed");
    // A fractional eta may move an otherwise NORMAL target across a binary64
    // midpoint. This is a declared finite law, never a normal-bitwise promise.
    D rawTarget(1);rawTarget.add(std::scalbn(1.,-53));rawTarget.add(D::product({eta,.25}));
    D finiteTarget(1);finiteTarget.add(std::scalbn(1.,-53));finiteTarget.add(D::product({eta,.25}).roundedBinary64Quantum());
    check(Exact::roundDouble(rawTarget).value==std::nextafter(1.,2.)&&Exact::roundDouble(finiteTarget).value==1.,
        "OHMIC_ETA normal target RN-tie witness missing");
    // Actual attempt49 native1575, independently retained raw IEEE operands.
    const double dt=fromBits(UINT64_C(0x3f1041a473928ee0));
    const auto a=D::product({dt,fromBits(UINT64_C(0x79480))});
    const auto b=D::product({dt,fromBits(UINT64_C(0x509de13))});
    check(equal(a.roundedBinary64Quantum(),D::product({eta,31}))&&equal(b.roundedBinary64Quantum(),D::product({eta,5242})),
        "OHMIC_ETA actual witness per-face integer oracle");
    auto raw=a;raw.add(b);auto applied=a.roundedBinary64Quantum();applied.add(b.roundedBinary64Quantum());
    check(equal(applied,D::product({eta,5273}))&&Exact::roundDouble(raw).value==fromBits(UINT64_C(0x1499)),
        "OHMIC_ETA actual RN or finite shared target changed");
    check(raw.compare(D::product({eta,5272}))>0&&raw.compare(D::product({eta,5273}))<0,"OHMIC_ETA actual adjacent representable interval");
    for(double count:{5272.,5273.}){auto residual=D::product({eta,count});residual.add(raw,-1);
        auto scale=D::product({eta,count});scale.add(raw);
        check(residual.absolute().compare(scale.multiplied(1e-10))>0,"OHMIC_ETA original infeasible raw gate lost");}
    auto rawError=applied;rawError.add(raw,-1);auto rawScale=applied;rawScale.add(raw);
    check(relative(rawError,rawScale)>3.6e-6&&relative(rawError,rawScale)<3.7e-6,"OHMIC_ETA sub-eta raw relative was rounded to zero");
    // Shared internal incidence cancels delta. Boundary delta, raw residual,
    // and ordinary inventory sigma remain independent exact ledgers.
    const double old[3]={eta,-2*eta,1};
    struct Face{int neg,pos;D raw;};
    const std::array<Face,3> faces{{{0,1,a},{1,2,b},{2,-1,D::product({eta,.5})}}};
    D rawNet[3],finiteNet[3],sumDelta,boundaryRaw,boundaryFinite,boundaryDelta;
    for(const auto& face:faces){const auto t=face.raw.roundedBinary64Quantum();auto delta=t;delta.add(face.raw,-1);
        for(int cell:{face.neg,face.pos})if(cell>=0){const int sign=cell==face.neg?1:-1;rawNet[cell].add(face.raw,sign);finiteNet[cell].add(t,sign);sumDelta.add(delta,sign);}
        if(face.pos<0){boundaryRaw.add(face.raw);boundaryFinite.add(t);boundaryDelta.add(delta);}}
    check(equal(sumDelta,boundaryDelta),"OHMIC_ETA internal quantum defects did not cancel exactly");
    auto rawResidual=boundaryRaw,finiteResidual=boundaryFinite;D sigma;
    for(int i=0;i<3;++i){auto target=D(old[i]);target.add(finiteNet[i],-1);const double stored=Exact::roundDouble(target).value;
        auto rounding=D(stored);rounding.add(target,-1);sigma.add(rounding);
        rawResidual.add(stored);rawResidual.add(old[i],-1);finiteResidual.add(stored);finiteResidual.add(old[i],-1);}
    auto expectedRaw=sigma;expectedRaw.add(boundaryDelta,-1);
    check(equal(finiteResidual,sigma)&&equal(rawResidual,expectedRaw),"OHMIC_ETA raw/finite/delta/sigma global identity");
}
void verify(const QuadtreeMesh& mesh,const ElectricVofField& field,bool zeroK,bool snapshot,bool expectQuantized=true,bool requireSubEtaWitness=true){
    const auto& receipt=field.ohmicTransferReceipt();const auto& report=field.report();
    check(receipt.rule==ElectricVofOhmicTransferReceipt::Rule::MinimumBinary64QuantumNearestEven&&receipt.quantumExponent==-1074&&
        receipt.duration==field.effectiveDt(),"OHMIC_ETA typed same-owner transfer receipt");
    const auto& old=field.implicitInputCharge();const auto& next=field.chargeCandidate().integratedCharge();
    std::vector<D> rawNet(old.size()),finiteNet(old.size()),rawGross(old.size()),finiteGross(old.size()),divD(old.size()),grossD(old.size());
    D boundaryRaw,boundaryFinite,boundaryDelta,sigma,rawGlobalScale,finiteGlobalScale;
    std::size_t quantized=0,nonzeroRawFaces=0,nonzeroFiniteFaces=0;double deltaMax=0,rawMax=0,finiteMax=0;
    for(const auto& face:field.faces()){
        const auto raw=D::product({receipt.duration,face.current});const auto applied=raw.roundedBinary64Quantum();auto delta=applied;delta.add(raw,-1);
        if(!equal(raw,D{}))++nonzeroRawFaces;if(!equal(applied,D{}))++nonzeroFiniteFaces;
        if(!equal(delta,D{}))++quantized;deltaMax=std::max(deltaMax,static_cast<double>(quanta(delta.absolute()).rounded()));
        check(delta.absolute().compare(D::product({eta,.5L}))<=0,"OHMIC_ETA independent nearest-eta error bound");
        if(zeroK||snapshot)check(equal(raw,D{})&&equal(applied,D{})&&equal(delta,D{}),"OHMIC_ETA zero K/snapshot produced transfer");
        for(int cell:{face.negativeCell,face.positiveCell})if(cell>=0){const int sign=cell==face.negativeCell?1:-1;
            rawNet[cell].add(raw,sign);finiteNet[cell].add(applied,sign);rawGross[cell].add(raw.absolute());finiteGross[cell].add(applied.absolute());
            divD[cell].add(face.displacement,sign);grossD[cell].add(std::abs(face.displacement));}
        if((face.negativeCell>=0)!=(face.positiveCell>=0)){const int sign=face.negativeCell>=0?1:-1;
            boundaryRaw.add(raw,sign);boundaryFinite.add(applied,sign);boundaryDelta.add(delta,sign);
            rawGlobalScale.add(raw.absolute());finiteGlobalScale.add(applied.absolute());}
    }
    auto rawGlobal=boundaryRaw,finiteGlobal=boundaryFinite;
    for(std::size_t i=0;i<old.size();++i){auto target=D(old[i]);target.add(finiteNet[i],-1);
        const double expected=equal(finiteNet[i],D{})?old[i]:Exact::roundDouble(target).value;
        check(std::memcmp(&expected,&next[i],sizeof(double))==0,"OHMIC_ETA stored Q differs from actual receipt including signed zero");
        auto error=D(next[i]);error.add(target,-1);sigma.add(error);
        auto finiteResidual=D(next[i]);finiteResidual.add(old[i],-1);finiteResidual.add(finiteNet[i]);
        auto rawResidual=D(next[i]);rawResidual.add(old[i],-1);rawResidual.add(rawNet[i]);
        auto fscale=finiteGross[i],rscale=rawGross[i];fscale.add(std::abs(old[i]));fscale.add(std::abs(next[i]));rscale.add(std::abs(old[i]));rscale.add(std::abs(next[i]));
        check(finiteResidual.absolute().compare(fscale.multiplied(1e-10))<=0,"OHMIC_ETA original finite local tolerance");
        finiteMax=std::max(finiteMax,relative(finiteResidual,fscale));rawMax=std::max(rawMax,relative(rawResidual,rscale));
        divD[i].add(next[i],-1);grossD[i].add(std::abs(next[i]));
        check(divD[i].absolute().compare(grossD[i].multiplied(1e-10))<=0,"OHMIC_ETA original actual physical Gauss gate");
        rawGlobal.add(next[i]);rawGlobal.add(old[i],-1);finiteGlobal.add(next[i]);finiteGlobal.add(old[i],-1);
        for(auto* scale:{&rawGlobalScale,&finiteGlobalScale}){scale->add(std::abs(old[i]));scale->add(std::abs(next[i]));}
    }
    auto expectedRaw=sigma;expectedRaw.add(boundaryDelta,-1);
    check(equal(finiteGlobal,sigma)&&equal(rawGlobal,expectedRaw),
        "OHMIC_ETA independent actual-factor raw/finite/sigma global identities");
    check(relative(rawGlobal,rawGlobalScale)<=1e-10&&relative(finiteGlobal,finiteGlobalScale)<=1e-10&&report.maximumRelativeGaussResidual<=1e-10,
        "OHMIC_ETA independent global balances and retained physical Gauss gate");
    if(!zeroK&&!snapshot&&requireSubEtaWitness){
        if(expectQuantized)check(quantized>0&&rawMax>1e-10,"OHMIC_ETA actual nonzero subnormal K/raw defect witness absent");
        else if(!(nonzeroFiniteFaces==0&&quantized>0&&rawMax>0)){
            std::ostringstream context;context<<std::setprecision(17)
                <<"OHMIC_ETA normal inventory/sub-eta raw residual hidden"
                <<" nonzero_raw_faces="<<nonzeroRawFaces<<" nonzero_finite_faces="<<nonzeroFiniteFaces
                <<" quantized_delta_faces="<<quantized<<" raw_max_relative="<<rawMax<<" finite_max_relative="<<finiteMax
                <<" max_abs_delta_eta="<<deltaMax<<" duration="<<receipt.duration
                <<" finite_faces_zero_condition="<<(nonzeroFiniteFaces==0)<<" quantized_positive_condition="<<(quantized>0)<<" raw_positive_condition="<<(rawMax>0);
            throw std::runtime_error(context.str());
        }
    }
    field.validateFor(mesh,field.geometry(),field.sourceCharge(),field.candidateTime());
}
void run(ElectricVofControls::Backend backend=ElectricVofControls::Backend::Gmres){
    ElectricVofControls controls;controls.backend=backend;
    arithmetic();
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=8;
    input.maxInitialIterations=0;input.maxCellNum=512;input.maxLevel_local={0};QuadtreeMesh mesh(input);
    FaceSegmentGraph graph;graph.build(mesh,[](const QuadtreeMesh&,const Cell&){return true;});
    const auto domain=NativeLeafDomain::build(mesh,graph,true);std::vector<double> alpha(mesh.cellIndexNum()),q(mesh.cellIndexNum());
    for(int ci:domain.cellIndexByRow())alpha[ci]=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord().x<.5?1:0;
    const Clock initial{0,0,0,Clock::Position::Initial},end{1,0,fromBits(UINT64_C(0x3f1041a473928ee0)),Clock::Position::End};
    const auto p=Publication::prepareLocalPlic(mesh,domain,initial,alpha,{}),stage=Publication::prepareLocalPlic(mesh,domain,end,alpha,{});
    q[0]=eta;q[7]=.001;q[15]=-.001;q[63]=-0.;const auto source=ElectricVofCharge::fromIntegrated(mesh,p,q);
    std::vector<ElectricVofBoundary> bc(domain.nativeFaceCount());
    for(int row=0;row<domain.nativeFaceCount();++row){const auto face=domain.face(mesh,row);
        if((face.negativeCell>=0)==(face.positiveCell>=0))continue;
        bc[row]=face.area==0?ElectricVofBoundary{ElectricVofBoundary::Type::HomogeneousNeumann,0}:
            ElectricVofBoundary{ElectricVofBoundary::Type::Dirichlet,.17+face.center.y};}
    const ElectricVofMaterials material{1,1,fromBits(UINT64_C(0x4555555)),0},zero{1,1,0,0};
#if defined(_OPENMP)
    const int oldThreads=omp_get_max_threads(),oldDynamic=omp_get_dynamic();omp_set_dynamic(0);omp_set_num_threads(1);
#endif
    const auto serial=ElectricVofField::prepareOhmic(mesh,stage,source,end,material,bc,controls);ElectricVofContract::backendReceipt(serial,backend);verify(mesh,serial,false,false);
    const auto snapshot=ElectricVofField::prepareElectrostaticSnapshot(mesh,p,source,material,bc,controls);ElectricVofContract::backendReceipt(snapshot,backend);verify(mesh,snapshot,false,true);
    check(snapshot.chargeCandidate().sameInstance(source),"OHMIC_ETA snapshot changed charge owner");
    const auto insulated=ElectricVofField::prepareOhmic(mesh,stage,source,end,zero,bc,controls);ElectricVofContract::backendReceipt(insulated,backend);verify(mesh,insulated,true,false);
    // All inventory scales are NORMAL; every transfer is fractional eta.
    // The raw residual is nonzero despite an unavailable absolute projection.
    std::vector<double> normalQ(q.size());for(std::size_t i=0;i<normalQ.size();++i)normalQ[i]=(i%2?1:-1)*1e-300;
    const auto normalSource=ElectricVofCharge::fromIntegrated(mesh,p,normalQ);auto normalBc=bc;
    for(int row=0;row<domain.nativeFaceCount();++row)if(normalBc[row].type==ElectricVofBoundary::Type::Dirichlet)
        normalBc[row].value+=.37*domain.face(mesh,row).center.x;
    const ElectricVofMaterials normalMaterial{1,1,500*eta,0};
    const auto normalInventory=ElectricVofField::prepareOhmic(mesh,stage,normalSource,end,normalMaterial,normalBc,controls);
    ElectricVofContract::backendReceipt(normalInventory,backend);
    verify(mesh,normalInventory,false,false,false);
#if defined(_OPENMP)
    omp_set_num_threads(4);const auto parallel=ElectricVofField::prepareOhmic(mesh,stage,source,end,material,bc,controls);ElectricVofContract::backendReceipt(parallel,backend);verify(mesh,parallel,false,false);
    ElectricVofPostprocessContract::payload(serial,parallel);
    omp_set_num_threads(oldThreads);omp_set_dynamic(oldDynamic);
#endif
    check(source.integratedCharge()==q,"OHMIC_ETA candidate rewrote source Q");
    std::cout<<"ELECTRIC_VOF_OHMIC_QUANTUM_PASS,eta_exp,-1074,actual_witness_native,1575,finite_charge,5273eta,old_raw_gate,infeasible,normal_RN_tie,declared_change,scenes,4"
        <<",production_ledger_reports,retired,independent_raw_finite_quantum_sigma_global,verified,raw_global_gate,1e-10,physical_Gauss_gate,1e-10\n";
}
} // namespace ElectricVofOhmicQuantumContract

namespace ElectricVofBackendContract {
using namespace ElectricVofContract;
// Producer-port witness on the SAME existing 160-leaf XY/RZ 2:1 mesh.
// This true geometric TPFA operator is deliberately not claimed to reproduce
// Field's hidden LS correction matrix; the original Field physical tests below
// remain the independent PDE/Gauss/finite-Q qualification.
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
void dumpBackendReceipt(const char* label,const ElectricVofPostprocessContract::Receipt& state) {
    std::cerr<<','<<label<<"_errno,"<<state.error<<','<<label<<"_csr,0x"<<std::hex<<state.csr
        <<','<<label<<"_native,";
    const auto* bytes=reinterpret_cast<const unsigned char*>(&state.native);
    for(std::size_t i=0;i<sizeof state.native;++i)std::cerr<<std::setw(2)<<std::setfill('0')<<static_cast<unsigned>(bytes[i]);
    std::cerr<<std::dec<<std::setfill(' ');
}
void dumpBackendDifference(const ElectricVofPostprocessContract::Receipt& one,const ElectricVofPostprocessContract::Receipt& other) {
    std::cerr<<",exit_csr_xor,0x"<<std::hex<<(one.csr^other.csr)<<std::dec<<",exit_native_diff,";
    const auto* a=reinterpret_cast<const unsigned char*>(&one.native);
    const auto* b=reinterpret_cast<const unsigned char*>(&other.native);
    for(std::size_t i=0;i<sizeof one.native;++i)if(a[i]!=b[i])
        std::cerr<<i<<':'<<static_cast<unsigned>(a[i])<<'/'<<static_cast<unsigned>(b[i])<<';';
}
void integratedRows(const QuadtreeMesh& mesh,bool rz,const std::vector<double>& signedQ) {
    using Solver=PolarizationPoissonFVSolverCellCentered;
    using Type=Solver::FittedBoundaryType;
    using ElectricVofPostprocessContract::sameBits;
    using ElectricVofPostprocessContract::bits;
    using ElectricVofPostprocessContract::Runtime;
    using ElectricVofPostprocessContract::Receipt;
    Runtime restore;omp_set_dynamic(0);
    FaceSegmentGraph graph;graph.build(mesh);
    const auto& segments=graph.segments();
    std::vector<double> kappa(segments.size());
    std::vector<Solver::FittedBoundary> boundaries(segments.size());
    std::vector<Solver::SegmentFluxCorrection> corrections(segments.size());
    std::vector<Solver::IntegratedFluxFactor> factors(segments.size());
    for(std::size_t sid=0;sid<segments.size();++sid){
        const auto& face=segments[sid];auto& boundary=boundaries[sid];auto& f=factors[sid];
        f.primary=1.25;f.secondary=.125;f.weight=.015625;
        f.secondaryTransfer=Solver::IntegratedFluxFactor::SecondaryTransfer::MinimumBinary64QuantumNearestEven;
        kappa[sid]=static_cast<double>(static_cast<long double>(f.primary)+static_cast<long double>(f.weight)*f.secondary);
        const double area=rz?CartesianMeasure::face(face.Af,face.faceCenter.x,true):face.Af;
        if(face.L>=0){
            boundary.type=Type::Interior;
            const double c=area/face.df;
            f.affine.coefficients={{mesh.cells(face.K).index(),c},{mesh.cells(face.L).index(),-c}};
            std::sort(f.affine.coefficients.begin(),f.affine.coefficients.end(),[](const auto& a,const auto& b){return a.cellIndex<b.cellIndex;});
        }else if(area==0){
            boundary.type=Type::HomogeneousNeumann; // Actual r=0 face, exact H=0.
        }else{
            boundary.type=Type::Dirichlet;boundary.value=affine(face.faceCenter,rz);
            const auto center=mesh.cells(face.K).centerCoord();
            const long double distance=face.isHorizontal
                ?static_cast<long double>(face.n.y)*(static_cast<long double>(face.faceCenter.y)-center.y)
                :static_cast<long double>(face.n.x)*(static_cast<long double>(face.faceCenter.x)-center.x);
            f.affine.coefficients={{mesh.cells(face.K).index(),static_cast<double>(static_cast<long double>(area)/distance)}};
            f.affine.constant=static_cast<double>(-(static_cast<long double>(area)/distance)*boundary.value);
        }
    }
    std::vector<Type> types;std::vector<std::vector<int>> correctionColumns(segments.size()),affineColumns(segments.size());
    for(std::size_t sid=0;sid<segments.size();++sid){types.push_back(boundaries[sid].type);
        for(const auto& c:corrections[sid].coefficients)correctionColumns[sid].push_back(c.cellIndex);
        for(const auto& c:factors[sid].affine.coefficients)affineColumns[sid].push_back(c.cellIndex);}
    auto ownedGraph=std::make_unique<FaceSegmentGraph>();ownedGraph->build(mesh);const auto* actualGraph=ownedGraph.get();
    const auto prepared=Solver::prepareIntegratedStructure(mesh,std::move(ownedGraph),types,correctionColumns,affineColumns,rz);
    check(prepared&&!ownedGraph&& &prepared->graph()==actualGraph,"v80 structure must own the actual submitted graph");
    CellField<double> rhs;rhs.initialize(mesh,0.);
    check(signedQ.size()==static_cast<std::size_t>(rhs.size()),"F10 native signed-Q shape");
    for(int i=0;i<rhs.size();++i)rhs.value[i]=signedQ[static_cast<std::size_t>(i)];
    const auto originalQ=signedQ;
    struct Result {CellField<double> phi;Solver::SolveReport report;std::vector<std::array<double,3>> terms;};
    const auto oldExpandedReference=[&](const Result& result) {
        namespace Exact=LAMP::Geometry::Detail::FiniteVolume;
        const int n=rhs.size();check(result.report.N==n,"V56 TPFA reference native row map changed");
        std::vector<Eigen::Triplet<double>> entries;Eigen::VectorXd b(n);
        for(int i=0;i<n;++i)b[i]=rhs.value[i];
        for(std::size_t sid=0;sid<segments.size();++sid){const auto& f=segments[sid];
            const int row=mesh.cells(f.K).index();const auto& bc=boundaries[sid];
            const double area=rz?CartesianMeasure::face(f.Af,f.faceCenter.x,true):f.Af;
            if(area==0||bc.type==Type::HomogeneousNeumann)continue;
            if(bc.type==Type::Dirichlet){const double distance=(f.faceCenter-mesh.cells(f.K).centerCoord()).dotProduct(f.n);
                const double c=kappa[sid]*area/distance;entries.emplace_back(row,row,c);b[row]+=c*bc.value;
            }else{const double c=kappa[sid]*area/f.df;entries.emplace_back(row,row,c);entries.emplace_back(row,mesh.cells(f.L).index(),-c);}
        }
        Eigen::SparseMatrix<double> a(n,n);a.setFromTriplets(entries.begin(),entries.end());a.makeCompressed();
        std::vector<Exact::Dyadic> exact(static_cast<std::size_t>(n));
        for(int i=0;i<n;++i)exact[i].add(b[i],-1);
        // Independent old CSC traversal, not the new row-routing implementation.
        for(int col=0;col<a.outerSize();++col)for(Eigen::SparseMatrix<double>::InnerIterator e(a,col);e;++e)
            for(double term:result.terms[e.col()])exact[e.row()].add(Exact::Dyadic::product({e.value(),term}));
        Eigen::VectorXd expanded(n);for(int i=0;i<n;++i)expanded[i]=Exact::roundDouble(exact[i]).value;
        const Eigen::Map<const Eigen::VectorXd> phi(result.phi.value.get(),n);const Eigen::VectorXd projected=a*phi-b;
        const double norm=b.norm();const double er=norm>0?expanded.norm()/norm:expanded.norm();
        const double pr=norm>0?projected.norm()/norm:projected.norm();
        check(expanded.allFinite()&&std::isfinite(er)&&er<=1e-11&&sameBits(result.report.finalResidual,std::max(pr,er)),
            "V56 independent old full-A three-term RN/norm report or original relTol differs");
        const auto& receipt=result.report;check(receipt.integratedExpandedResidualAvailable&&
            receipt.integratedExpandedResidualRows==static_cast<std::size_t>(n),"V56 producer full-A receipt row scope");
    };
    std::size_t pairs=0,fallbackPairs=0,storageSwapPairs=0,storageSwapGeometryGuards=0;
    for(bool superlu:{false,true}){
        Runtime backendRestore; // Keep artificial negative-input flags within this backend fixture.
        Solver::Options options;options.axisymmetric=rz;options.useSuperLU=superlu;options.requireSuperLU=superlu;
        options.useGMRES=!superlu;options.useCG=false;options.useIC=false;options.relTol=1e-11;
        std::fenv_t entry{};check(std::fegetenv(&entry)==0,"F10 common FP capture");const auto entryCsr=_mm_getcsr();
        const auto reset=[&]{check(std::fesetenv(&entry)==0,"F10 common FP restore");_mm_setcsr(entryCsr);errno=173;};
        const auto solve=[&](Receipt* callEntry=nullptr,Receipt* callExit=nullptr){Result out;out.phi.initialize(mesh,0.);Solver solver;
            reset();if(callEntry)callEntry->capture();
            const bool ok=solver.solveFittedIntegratedBoundary(mesh,graph,rhs,out.phi,kappa,boundaries,corrections,
                options,&out.report,&factors,&out.terms);
            if(callExit)callExit->capture();
            if(!ok)throw std::runtime_error(std::string("F10 strict integrated producer failed: ")+out.report.message);
            return out;};
        Receipt serialEntry,serialBodyExit,parallelEntry,parallelBodyExit;
        omp_set_num_threads(1);reset();const auto serial=solve(&serialEntry,&serialBodyExit);Receipt sf;sf.capture();
        omp_set_num_threads(4);unsigned observedTeam=0;
#pragma omp parallel num_threads(4)
        {
#pragma omp master
            {observedTeam=static_cast<unsigned>(omp_get_num_threads());}
        }
        reset();const auto parallel=solve(&parallelEntry,&parallelBodyExit);Receipt pf;pf.capture();
        oldExpandedReference(serial);oldExpandedReference(parallel);
        if(!sf.same(pf)){
            std::cerr<<"F10_GENERIC_SUCCESS_DIAGNOSTIC,scene,"<<(rz?"RZ":"XY")
                <<",backend,"<<(superlu?"SuperLU":"GMRES")<<",entry_equal,"<<serialEntry.same(parallelEntry)
                <<",body_exit_equal,"<<serialBodyExit.same(parallelBodyExit);
            dumpBackendReceipt("entry1",serialEntry);dumpBackendReceipt("entry4",parallelEntry);
            dumpBackendReceipt("body_exit1",serialBodyExit);dumpBackendReceipt("body_exit4",parallelBodyExit);
            dumpBackendReceipt("exit1",sf);dumpBackendReceipt("exit4",pf);dumpBackendDifference(sf,pf);std::cerr<<'\n';
        }
        check(serialEntry.same(parallelEntry),"F10 comparison supplied different actual solver entries");
        check(sf.same(pf),"F10 phase execution changed complete FP/errno");
        check(serial.phi.size()==parallel.phi.size()&&std::memcmp(serial.phi.value.get(),parallel.phi.value.get(),
            static_cast<std::size_t>(serial.phi.size())*sizeof(double))==0,"F10 native phi bits/zero signs");
        for(int i=0;i<serial.phi.size();++i)check(serial.phi.label(i)==parallel.phi.label(i),"F10 native phi label");
        bits(serial.terms,parallel.terms,"F10 three actual potential-term bits/zero signs");
        const auto& a=serial.report;const auto& b=parallel.report;
        check(a.success&&b.success&&a.usedSuperLU==superlu&&b.usedSuperLU==superlu&&
            a.usedGMRES==!superlu&&b.usedGMRES==!superlu,"F10 actual backend not exclusive");
        for(auto member:{&Solver::SolveReport::finalResidual,&Solver::SolveReport::cgResidual})
            check(sameBits(a.*member,b.*member),"F10 original numeric residual report bits");
        check(a.iterations==b.iterations&&a.N==b.N&&a.nnz==b.nnz&&a.gmresIterations==b.gmresIterations&&
            a.gmresRestarts==b.gmresRestarts&&a.usedCG==b.usedCG&&a.cgAttempted==b.cgAttempted&&
            a.cgFallbackToDirect==b.cgFallbackToDirect&&a.cgIterations==b.cgIterations&&a.cgPatternReused==b.cgPatternReused&&
            a.cgTemporalGuessUsed==b.cgTemporalGuessUsed&&a.directPatternReused==b.directPatternReused&&
            a.directMatrixFrozen==b.directMatrixFrozen&&a.fittedFactorizationReused==b.fittedFactorizationReused&&
            a.fittedRhsOnlySolve==b.fittedRhsOnlySolve&&a.message==b.message,
            "F10 original non-timing/non-worker report changed");
        check(a.integratedOrderedRowsAvailable&&b.integratedOrderedRowsAvailable&&a.integratedDefectPhases>0&&
            a.integratedDefectPhases<=2&&a.integratedProjectionPhases>0&&a.integratedProjectionPhases<=a.integratedDefectPhases&&
            a.integratedDefectPhases==b.integratedDefectPhases&&a.integratedProjectionPhases==b.integratedProjectionPhases,
            "F10 actual two-phase receipt missing");
        check(a.integratedParallelDefectRows==0&&a.integratedParallelProjectionRows==0&&
            a.integratedRowsMaximumWorkers==0&&a.integratedRowsFallbackPhases==0,"F10 one-worker path must remain literal");
        if(observedTeam>=2){
            check(b.integratedParallelDefectRows==static_cast<std::size_t>(b.N)*b.integratedDefectPhases&&
                b.integratedParallelProjectionRows==static_cast<std::size_t>(b.N)*b.integratedProjectionPhases&&
                b.integratedRowsMaximumWorkers>=2&&b.integratedRowsMaximumWorkers<=4&&b.integratedRowsFallbackPhases==0,
                "F10 eligible producer did not actually consume parallel rows");
        }else{
            // Root can run this SAME selector with OMP_THREAD_LIMIT=1. The
            // requested team remains four; actual worker shortage must join,
            // restore phase entry and replay literal WITHOUT repeating solveDelta.
            check(observedTeam==1&&b.integratedParallelDefectRows==0&&b.integratedParallelProjectionRows==0&&
                b.integratedRowsMaximumWorkers==0&&b.integratedRowsFallbackPhases==b.integratedDefectPhases+b.integratedProjectionPhases,
                "F10 reduced-team phase fallback/replay receipt missing");++fallbackPairs;
        }
        check(std::memcmp(rhs.value.get(),originalQ.data(),originalQ.size()*sizeof(double))==0,"F10 producer changed independent signed Q");
        ++pairs;
        const auto compareResult=[&](const Result& left,const Result& right){
            check(left.phi.size()==right.phi.size()&&std::memcmp(left.phi.value.get(),right.phi.value.get(),
                static_cast<std::size_t>(left.phi.size())*sizeof(double))==0,"v80 generic/prepared phi bits or zero signs changed");
            bits(left.terms,right.terms,"v80 generic/prepared actual three-term expansion changed");
            check(left.report.success&&right.report.success&&left.report.N==right.report.N&&left.report.nnz==right.report.nnz
                &&left.report.usedSuperLU==right.report.usedSuperLU&&left.report.usedGMRES==right.report.usedGMRES
                &&left.report.gmresIterations==right.report.gmresIterations&&left.report.gmresRestarts==right.report.gmresRestarts
                &&sameBits(left.report.finalResidual,right.report.finalResidual),"v80 original actual-A/backend residual changed");
        };
        const auto preparedSolve=[&]{Result out;out.phi.initialize(mesh,0.);Solver solver;reset();
            const bool ok=solver.solvePreparedIntegratedBoundary(mesh,*prepared,rhs,out.phi,kappa,boundaries,corrections,
                options,&out.report,&factors,&out.terms);
            if(!ok)throw std::runtime_error(std::string("v80 prepared integrated producer failed: ")+out.report.message);
            return out;};
        const auto originalKappa=kappa;const auto originalBoundary=boundaries;const auto originalFactors=factors;
        for(int variant=0;variant<3;++variant){
            if(variant){for(std::size_t sid=0;sid<segments.size();++sid){auto& f=factors[sid];const auto& face=segments[sid];
                    f.primary=1.25+.03125*(variant+1)*(1+face.faceCenter.x+face.faceCenter.y/7);
                    f.secondary=variant==2?0:.0625*(1+face.faceCenter.x/3);f.weight=variant==1?.03125:.0078125;
                    kappa[sid]=static_cast<double>(static_cast<long double>(f.primary)+static_cast<long double>(f.weight)*f.secondary);
                    if(boundaries[sid].type==Type::Dirichlet){boundaries[sid].value=affine(face.faceCenter,rz)+variant*.0625;
                        f.affine.constant=static_cast<double>(-static_cast<long double>(f.affine.coefficients.front().value)*boundaries[sid].value);}}
                for(int i=0;i<rhs.size();++i)rhs.value[i]=signedQ[static_cast<std::size_t>(i)]*(1+variant*.125)+(i%3-1)*.0001;
            }
            // Cold generic versus prepared compares physical payload only.
            omp_set_num_threads(1);reset();const auto generic=solve();
            reset();const auto preparedSerial=preparedSolve();Receipt one;one.capture();
            omp_set_num_threads(8);reset();const auto preparedParallel=preparedSolve();Receipt eight;eight.capture();
            compareResult(generic,preparedSerial);compareResult(preparedSerial,preparedParallel);
            oldExpandedReference(preparedSerial);oldExpandedReference(preparedParallel);
            check(one.same(eight),"v80 prepared numeric1/8 native/rawCSR/errno changed");
            if(variant)check(std::memcmp(serial.phi.value.get(),preparedSerial.phi.value.get(),
                static_cast<std::size_t>(serial.phi.size())*sizeof(double))!=0,"v80 changed material/BC/Q/duration reused old numeric solution");
        }
        kappa=originalKappa;boundaries=originalBoundary;factors=originalFactors;
        for(int i=0;i<rhs.size();++i)rhs.value[i]=signedQ[static_cast<std::size_t>(i)];
        {
            // AMR prepares the target while it is owned by a candidate mesh,
            // then commits its storage into the accepted mesh wrapper. The
            // prepared structure follows that exact storage, not its address.
            Runtime storageSwapRestore;
            auto candidate=mesh.clone(),accepted=mesh.clone();
            auto& oldVertex=accepted->verts(accepted->vertIndexToLabel[0]);
            oldVertex.setCoord(oldVertex.coord()+Vector2D{0,.0009765625});
            auto targetGraph=std::make_unique<FaceSegmentGraph>();targetGraph->build(*candidate);
            const auto targetPrepared=Solver::prepareIntegratedStructure(*candidate,std::move(targetGraph),types,
                correctionColumns,affineColumns,rz);
            check(candidate.get()!=accepted.get(),"prepared storage-swap fixture reused the candidate mesh wrapper");
            accepted->swapStorage(*candidate);
            FaceSegmentGraph freshGraph;freshGraph.build(*accepted);
            CellField<double> targetRhs;targetRhs.initialize(*accepted,0.);
            for(int i=0;i<targetRhs.size();++i)targetRhs.value[i]=signedQ[static_cast<std::size_t>(i)];
            const auto targetSolve=[&](bool reuse){Result out;out.phi.initialize(*accepted,0.);Solver solver;reset();
                const bool ok=reuse
                    ?solver.solvePreparedIntegratedBoundary(*accepted,*targetPrepared,targetRhs,out.phi,kappa,boundaries,corrections,
                        options,&out.report,&factors,&out.terms)
                    :solver.solveFittedIntegratedBoundary(*accepted,freshGraph,targetRhs,out.phi,kappa,boundaries,corrections,
                        options,&out.report,&factors,&out.terms);
                if(!ok)throw std::runtime_error(std::string("prepared storage-swap producer failed: ")+out.report.message);
                return out;};
            std::optional<Result> firstSwapped;Receipt firstSwappedExit;
            for(int threads:{1,8}){omp_set_num_threads(threads);
                const auto fresh=targetSolve(false),reused=targetSolve(true);Receipt exit;exit.capture();
                compareResult(serial,fresh);compareResult(fresh,reused);oldExpandedReference(reused);
                if(firstSwapped){compareResult(*firstSwapped,reused);
                    check(firstSwappedExit.same(exit),"prepared storage-swap numeric1/8 native/rawCSR/errno changed");}
                else{firstSwapped.emplace(std::move(reused));firstSwappedExit=exit;}
            }
            ++storageSwapPairs;
            auto& changedVertex=accepted->verts(accepted->vertIndexToLabel[0]);const auto coordinate=changedVertex.coord();
            changedVertex.setCoord(coordinate+Vector2D{0,.0009765625});
            std::string firstError;Receipt firstEntry,firstExit;
            for(int threads:{1,8}){omp_set_num_threads(threads);CellField<double> phi;phi.initialize(*accepted,7.25);
                auto* storage=phi.value.get();const std::vector<double> saved(static_cast<std::size_t>(phi.size()),7.25);
                std::vector<std::array<double,3>> terms(static_cast<std::size_t>(phi.size()),{{.25,-.5,.75}});const auto savedTerms=terms;
                Solver solver;Solver::SolveReport report;bool ok=false;std::string error;reset();Receipt callEntry;callEntry.capture();
                try{ok=solver.solvePreparedIntegratedBoundary(*accepted,*targetPrepared,targetRhs,phi,kappa,boundaries,corrections,
                    options,&report,&factors,&terms);error=report.message;}
                catch(const std::exception& e){error=std::string(typeid(e).name())+":"+e.what();}
                Receipt exit;exit.capture();
                check(!ok&&error.find("[NativeMeshIdentity]")!=std::string::npos
                    &&error.find("mesh geometry/topology changed since identity construction")!=std::string::npos,
                    "prepared storage-swap structure accepted a changed live vertex or lost its exact identity guard");
                check(phi.value.get()==storage&&std::memcmp(phi.value.get(),saved.data(),saved.size()*sizeof(double))==0,
                    "prepared storage-swap geometry failure changed phi owner/payload");
                bits(terms,savedTerms,"prepared storage-swap geometry failure changed potential terms");
                if(threads==1){firstError=error;firstEntry=callEntry;firstExit=exit;}
                else check(firstError==error&&firstEntry.same(callEntry)&&firstExit.same(exit),
                    "prepared storage-swap geometry failure changed first error/FP/errno between1/8");
                ++storageSwapGeometryGuards;
            }
            changedVertex.setCoord(coordinate);
        }
        std::size_t interior=0;while(interior<segments.size()&&segments[interior].L<0)++interior;
        check(interior<segments.size(),"v80 bad-input fixture lacks actual reversed interior faces");
        for(int kind=0;kind<5;++kind){auto badKappa=kappa;auto badBoundary=boundaries;auto badFactors=factors;
            if(kind==0){badKappa[interior]=-1;badBoundary[interior].value=std::numeric_limits<double>::quiet_NaN();}
            if(kind==1)badFactors[interior].affine.coefficients.front().cellIndex=rhs.size()+3;
            if(kind==2)badBoundary[interior].type=Type::Dirichlet;
            if(kind==3)badFactors[interior].affine.coefficients.front().value*=1.125; // Same schema, broken reverse action.
            if(kind==4){for(std::size_t sid=0;sid<segments.size();++sid){badKappa[sid]=std::numeric_limits<double>::denorm_min();
                    badFactors[sid].primary=badKappa[sid];badFactors[sid].secondary=0;badFactors[sid].weight=0;}}
            std::string firstError;Receipt firstEntry,firstExit;
            for(int threads:{1,8}){omp_set_num_threads(threads);CellField<double> phi;phi.initialize(mesh,7.25);
                auto* storage=phi.value.get();std::vector<double> saved(static_cast<std::size_t>(phi.size()),7.25);
                std::vector<std::array<double,3>> terms(static_cast<std::size_t>(phi.size()),{{.25,-.5,.75}});const auto savedTerms=terms;
                Solver solver;Solver::SolveReport report;bool ok=false;std::string error;
                // Field initialization runs its own OpenMP setup; compare the solve from identical entries.
                reset();Receipt callEntry;callEntry.capture();
                try{ok=solver.solvePreparedIntegratedBoundary(mesh,*prepared,rhs,phi,badKappa,badBoundary,corrections,
                    options,&report,&badFactors,&terms);error=report.message;}
                catch(const std::exception& e){error=std::string(typeid(e).name())+":"+e.what();}
                Receipt exit;exit.capture();check(!ok&&!error.empty(),"v80 prepared invalid dynamic input accepted");
                check(phi.value.get()==storage&&std::memcmp(phi.value.get(),saved.data(),saved.size()*sizeof(double))==0,
                    "v80 prepared failure changed original phi owner/payload");bits(terms,savedTerms,"v80 prepared failure changed potential terms");
                if(threads==1){firstError=error;firstEntry=callEntry;firstExit=exit;}else{
                    check(firstEntry.same(callEntry),"v80 prepared negative test supplied different solve entry environments");
                    const bool same=firstError==error&&firstExit.same(exit);
                    if(!same){
                        std::cerr<<"V80_PREPARED_FAILURE_DIAGNOSTIC,scene,"<<(rz?"RZ":"XY")
                            <<",backend,"<<(superlu?"SuperLU":"GMRES")<<",kind,"<<kind
                            <<",error1,"<<std::quoted(firstError)<<",error8,"<<std::quoted(error)
                            <<",entry_equal,"<<firstEntry.same(callEntry);
                        dumpBackendReceipt("entry1",firstEntry);dumpBackendReceipt("entry8",callEntry);
                        dumpBackendReceipt("exit1",firstExit);dumpBackendReceipt("exit8",exit);dumpBackendDifference(firstExit,exit);
                        std::cerr<<'\n';
                    }
                    check(same,"v80 prepared failure changed original first error/FP/errno between1/8");
                }
            }
        }
    }
    std::cout<<"ELECTRIC_INTEGRATED_ORDERED_ROWS_PASS,"<<(rz?"RZ":"XY")<<",same_existing_true2to1_mesh,1,backend_pairs,"<<pairs
        <<",fallback_pairs,"<<fallbackPairs<<",storage_swap_backend_pairs,"<<storageSwapPairs<<",storage_swap_geometry_guards,"<<storageSwapGeometryGuards
        <<",native_phi_terms_old_report_bits,exact,FP_errno,exact,old_CSC_fullA_reference,independent_TPFA,hidden_Field_LS_operator_claim,0\n";
}
#endif
void run() {
    check(ElectricVofControls{}.backend==Backend::Gmres,"E_BACKEND changed the direct module default");
    ElectricVofControls controls;controls.backend=Backend::SuperLU;
    std::size_t defaults=0,explicitDefaults=0,ohmic=0,snapshots=0,guards=0;
    for(bool rz:{false,true}){
        auto mesh=meshFor();const auto domain=domainFor(*mesh,rz);
        check(!domain.uniformLeaves(),"E_BACKEND requires the existing actual 2:1 affine fixture");
        const auto stocks=phaseStocks(*mesh,domain);
        const Clock acceptedClock{0,.2,.2,Clock::Position::Initial},geometryClock{1,.2,.22,Clock::Position::Midpoint};
        const Clock candidate{1,.2,.22,rz?Clock::Position::Midpoint:Clock::Position::End};
        const auto accepted=Publication::prepareLocalPlicInventory(*mesh,domain,acceptedClock,stocks,{});
        const auto geometry=Publication::prepareLocalPlicInventory(*mesh,domain,geometryClock,stocks,{});
        const auto acceptedBefore=snapshot(accepted),geometryBefore=snapshot(geometry);
        std::vector<double> arbitrary(mesh->cellIndexNum());bool positive=false,negative=false;
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto c=mesh->cells(mesh->cellIndexToLabel[ci]).centerCoord();
            arbitrary[ci]=domain.cellVolumes()[domain.cellRowByIndex()[ci]]*(c.x-.5+.19*(c.y-.5));
            positive=positive||arbitrary[ci]>0;negative=negative||arbitrary[ci]<0;}
        check(positive&&negative,"E_BACKEND lost the existing independent signed-Q fixture");
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        integratedRows(*mesh,rz,arbitrary);
#endif
        for(int kind=0;kind<3;++kind){
            ElectricVofMaterials material{2.,5.,.7,.03};if(kind==2)material.conductivityLiquid=material.conductivityGas=0;
            std::vector<ElectricVofContract::Boundary> boundaries;const auto reference=references(*mesh,geometry,material,candidate.time()-acceptedClock.time(),boundaries);
            auto q=arbitrary;
            if(kind==0){std::vector<long double> analytic(q.size());
                for(const auto& face:reference)scatter(analytic,face.face,static_cast<long double>(face.a)*face.affineH);
                for(std::size_t ci=0;ci<q.size();++ci)q[ci]=static_cast<double>(analytic[ci]);}
            const auto charge=ElectricVofCharge::fromIntegrated(*mesh,accepted,q);
            if(kind==0){
                const auto defaultField=ElectricVofField::prepareOhmic(*mesh,geometry,charge,candidate,material,boundaries);
                verify(*mesh,geometry,charge,defaultField,candidate,material,reference,true,false);++defaults;
                ElectricVofControls gmres;gmres.backend=Backend::Gmres;
                const auto explicitField=ElectricVofField::prepareOhmic(*mesh,geometry,charge,candidate,material,boundaries,gmres);
                verify(*mesh,geometry,charge,explicitField,candidate,material,reference,true,false,Backend::Gmres);
                // Existing byte comparisons include signed zero, all potential
                // terms/face H,D,J/Q and every physical numeric report field.
                ElectricVofPostprocessContract::payload(defaultField,explicitField);
                const auto& x=defaultField.report();const auto& y=explicitField.report();
                check(x.usedGMRES==y.usedGMRES&&x.linearIterations==y.linearIterations&&x.linearRestarts==y.linearRestarts&&
                    x.linearMaximumWorkers==y.linearMaximumWorkers&&x.postprocessParallelUsed==y.postprocessParallelUsed&&
                    x.postprocessParallelFallback==y.postprocessParallelFallback&&x.postprocessWorkers==y.postprocessWorkers&&
                    x.postprocessScratchBytes==y.postprocessScratchBytes,"E_BACKEND default/explicit GMRES non-timing execution report differs");
                ++explicitDefaults;
            }
            const auto field=ElectricVofField::prepareOhmic(*mesh,geometry,charge,candidate,material,boundaries,controls);
            verify(*mesh,geometry,charge,field,candidate,material,reference,kind==0,kind==2,Backend::SuperLU);++ohmic;
            if(kind==1){
                const auto passive=ElectricVofField::prepareElectrostaticSnapshot(*mesh,accepted,charge,material,boundaries,controls);
                ElectricVofSnapshotContract::verifySnapshot(*mesh,accepted,charge,passive,material,Backend::SuperLU);++snapshots;
                auto wrong=candidate;++wrong.step;
                check(wrong.time()==candidate.time(),"E_BACKEND clock guard must retain physical time");
                rejects([&]{field.validateFor(*mesh,geometry,charge,wrong);});++guards;
                rejects([&]{ElectricVofField::prepareOhmic(*mesh,geometry,charge,wrong,material,boundaries,controls);});++guards;
                auto invalid=controls;invalid.backend=static_cast<Backend>(99);
                rejects([&]{ElectricVofField::prepareOhmic(*mesh,geometry,charge,candidate,material,boundaries,invalid);});++guards;
                rejects([&]{ElectricVofField::prepareElectrostaticSnapshot(*mesh,accepted,charge,material,boundaries,invalid);});++guards;
            }
            check(std::memcmp(q.data(),charge.integratedCharge().data(),q.size()*sizeof(double))==0&&
                snapshot(accepted)==acceptedBefore&&snapshot(geometry)==geometryBefore,
                "E_BACKEND successful/rejected fields changed accepted Q or either VOF owner");
        }
    }
    check(defaults==2&&explicitDefaults==2&&ohmic==6&&snapshots==2&&guards==8,"E_BACKEND required direct default/explicit/fullclock witnesses missing");
    // These are the original independent fixtures and thresholds, with only
    // the explicitly selected backend changed. They reject backend failure.
    ElectricVofSelectedDomainContract::run(Backend::SuperLU);
    ElectricVofWeakFieldGaussContract::run(Backend::SuperLU);
    ElectricVofOhmicQuantumContract::run(Backend::SuperLU);
    std::cout<<"ELECTRIC_VOF_BACKEND_PASS,default_GMRES,"<<defaults<<",explicit_GMRES_payload_pairs,"<<explicitDefaults<<",explicit_SuperLU_Ohmic,"<<ohmic
        <<",explicit_SuperLU_snapshot,"<<snapshots<<",enum_fullclock_guards,"<<guards
        <<",selected_original_scenes,4,weak_Gauss_original_scenes,3,quantum_original_scenes,4"
        <<",XY_RZ_true_2to1,independent_affine_signed_Q_zero_K_snapshot_selected_weak_finite_receipts,original_gates_unchanged\n";
}
} // namespace ElectricVofBackendContract

int main(int argc,char** argv) try {
    if(argc==2&&std::string(argv[1])=="--electric-vof-backend-only"){ElectricVofBackendContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-ohmic-quantum-only"){ElectricVofOhmicQuantumContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-weak-field-gauss-only"){ElectricVofWeakFieldGaussContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-charge-remap-second-order-only"){ElectricVofSecondOrderRemapContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-selected-amr-only"){ElectricVofSelectedAmrContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-selected-chain-only"){ElectricVofSelectedChainContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-selected-domain-only"){ElectricVofSelectedDomainContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-charge-second-order-only"){ElectricVofSecondOrderChargeContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-charge-parallel-only"){ElectricVofChargeParallelContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-postprocess-parallel-only"){ElectricVofPostprocessContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-validation-batch-only"){ElectricVofValidationBatchContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-hho-frame-only"){ElectricContract::localFrame();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-counterflow-only"){ElectricVofCounterflowContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-dielectric-inclusion-only"){ElectricVofDielectricInclusionContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-boundary-only"){ElectricVofBoundaryChargeContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-layer-relaxation-only"){ElectricVofLayerRelaxationContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-snapshot-only"){ElectricVofSnapshotContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-maxwell-rt0-only"){ElectricVofMaxwellContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-gradient-only"){ElectricVofGradientContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-binary64-only"){ElectricVofBinary64Contract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-moving-only"){ElectricVofMovingContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-vof-ohmic-only"){ElectricVofContract::run();return 0;}
    if(argc>2||(argc==2&&std::string(argv[1])!="--electric-p0-only"&&std::string(argv[1])!="--electric-p1-only"
        &&std::string(argv[1])!="--electric-p2-only"&&std::string(argv[1])!="--electric-p3-only"&&std::string(argv[1])!="--electric-p4-only"&&std::string(argv[1])!="--electric-p5-only"))
        throw std::invalid_argument("--electric-p0-only | --electric-p1-only | --electric-p2-only | --electric-p3-only | --electric-p4-only | --electric-p5-only");
    if(argc==2&&std::string(argv[1])=="--electric-p0-only"){ElectricContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-p1-only"){ElectricPhaseContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-p2-only"){ElectricPhaseContract::run(true);return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-p3-only"){ElectricPhaseContract::run(true,true);return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-p4-only"){ElectricPhaseContract::run(true,true,true);return 0;}
    if(argc==2&&std::string(argv[1])=="--electric-p5-only"){ElectricPhaseContract::run(true,true,true,true);return 0;}
    ElectricContract::run();ElectricPhaseContract::run(true,true,true,true);
    std::cout<<"Electric phase, charge, current, moving-stage and Maxwell contracts passed\n";
    return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
