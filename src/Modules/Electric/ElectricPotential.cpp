#include "Modules/Electric/ElectricPotential.h"
#include "Modules/Electric/ElectricPhaseGeometry.h"
#include "Modules/Electric/SurfaceChargeField.h"

#include "Solvers/Poisson/ElectricLocalOperator.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
#include <cstring>
#include <type_traits>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Electric
{
namespace
{
namespace Electric = LAMP::ElectricDiscretization;
using TimingClock = std::chrono::steady_clock;
double elapsedSeconds(TimingClock::time_point begin)
{
    return std::chrono::duration<double>(TimingClock::now()-begin).count();
}
constexpr std::array<double,8> gaussPoints{{
    -.96028985649753623168,-.79666647741362673959,
    -.52553240991632898582,-.18343464249564980494,
     .18343464249564980494, .52553240991632898582,
     .79666647741362673959, .96028985649753623168}};
constexpr std::array<double,8> gaussWeights{{
    .10122853629037625915,.22238103445337447054,
    .31370664587788728734,.36268378337836198297,
    .36268378337836198297,.31370664587788728734,
    .22238103445337447054,.10122853629037625915}};
void require(bool valid,const char* reason)
{
    if(!valid) throw std::invalid_argument(std::string("[ElectricPotential] ")+reason);
}
int boundarySide(Boundary boundary)
{
    switch(boundary) {
    case Boundary::west: return 0;
    case Boundary::south: return 1;
    case Boundary::east: return 2;
    case Boundary::north: return 3;
    case Boundary::non: return -1;
    default: throw std::invalid_argument("[ElectricPotential] invalid face boundary");
    }
}
Vector2D sideNormal(int side)
{
    switch(side) {
    case 0:return {-1,0};
    case 1:return {0,-1};
    case 2:return {1,0};
    case 3:return {0,1};
    default:throw std::invalid_argument("[ElectricPotential] invalid cell side");
    }
}
double coordinate(double lower,double upper,double q)
{
    const double width=upper-lower;
    const double result=q>0?upper-width*(1-q)/2:lower+width*(1+q)/2;
    require(std::isfinite(result) && result>lower && result<upper,
        "Gauss point has no represented interior coordinate");
    return result;
}
struct FaceData
{
    int edgeLabel=-1,boundary=-1,offset=-1;
    bool horizontal=false,axis=false,prescribed=false;
    std::array<Vector2D,2> bounds;
    std::vector<Electric::FaceQuadraturePoint> nodes;
    std::array<double,3> coefficients{{0,0,0}};
};
struct CellData
{
    int label=-1;
    Vector2D center;
    double scale=0;
    std::array<Vector2D,2> bounds;
    std::array<int,4> faceRows{{-1,-1,-1,-1}};
    std::array<int,4> localFace{{-1,-1,-1,-1}};
    std::array<double,6> cell{{0,0,0,0,0,0}};
    std::array<double,10> potential{{0,0,0,0,0,0,0,0,0,0}};
    std::array<std::array<double,3>,4> residual{};
};
struct PreparedCell
{
    Electric::LocalOperator op;
    Electric::CellVector cellLoad;
    Eigen::VectorXd faceLoad,condensedLoad;
    std::vector<int> faceRows;
    PreparedCell(Electric::LocalOperator prepared,Electric::CellVector load,
        Eigen::VectorXd faceRhs,std::vector<int> rows)
        :op(std::move(prepared)),cellLoad(std::move(load)),
         faceLoad(std::move(faceRhs)),faceRows(std::move(rows))
    {condensedLoad=op.condensedRhs(cellLoad,faceLoad);}
};
double normalized(double residual,double scale)
{
    require(std::isfinite(residual) && residual>=0 && std::isfinite(scale) && scale>=0,
        "nonfinite residual or residual scale");
    if(residual==0)return 0;
    require(scale>0,"nonzero residual with zero equation scale");
    const double value=residual/scale;
    require(std::isfinite(value),"unrepresentable normalized residual");
    return value;
}
template<std::size_t N,class Derived>
std::array<double,N> toArray(const Eigen::MatrixBase<Derived>& vector)
{
    require(vector.size()==static_cast<Eigen::Index>(N) && vector.allFinite(),
        "invalid recovered polynomial coefficients");
    std::array<double,N> result{};
    for(std::size_t j=0;j<N;++j)result[j]=vector[static_cast<Eigen::Index>(j)];
    return result;
}
void requirePoint(const CellData& cell,Vector2D point)
{
    require(std::isfinite(point.x) && std::isfinite(point.y)
        && point.x>=cell.bounds[0].x && point.x<=cell.bounds[1].x
        && point.y>=cell.bounds[0].y && point.y<=cell.bounds[1].y,
        "physical polynomial query must stay in the requested closed native cell");
}
void validateControls(const ElectricCartesianControls& controls)
{
    require(std::isfinite(controls.relativeLinearTolerance) && controls.relativeLinearTolerance>0
        && controls.relativeLinearTolerance<1 && std::isfinite(controls.relativeSymmetryTolerance)
        && controls.relativeSymmetryTolerance>0 && controls.relativeSymmetryTolerance<=1e-12
        && std::isfinite(controls.relativeConservationTolerance)
        && controls.relativeConservationTolerance>0 && controls.relativeConservationTolerance<=1e-10,
        "invalid linear, symmetry or conservation tolerance");
}
void validateBoundaries(const std::array<ElectricCartesianBoundary,4>& boundaries)
{
    for(const auto& bc:boundaries) {
        require(bc.type==ElectricCartesianBoundary::Type::Dirichlet
            || bc.type==ElectricCartesianBoundary::Type::Neumann,"invalid typed boundary law");
        require(static_cast<bool>(bc.value),"every physical boundary requires an explicit value provider");
    }
}
bool prepareFaceUnknown(FaceData& face,const std::array<ElectricCartesianBoundary,4>& boundaries,
    bool axisymmetric,ElectricCartesianReport& report)
{
    face.prescribed=face.boundary>=0
        && boundaries[static_cast<std::size_t>(face.boundary)].type==ElectricCartesianBoundary::Type::Dirichlet;
    if(face.prescribed) {
        Eigen::Matrix3d mass=Eigen::Matrix3d::Zero();
        Electric::FacePolynomial rhs=Electric::FacePolynomial::Zero();
        const auto& provider=boundaries[static_cast<std::size_t>(face.boundary)].value;
        for(const auto& node:face.nodes) {
            const auto basis=Electric::faceBasis(node.coordinate);
            const double value=provider(node.point);
            const double weight=node.weight*Electric::physicalWeight(node.point,axisymmetric);
            require(std::isfinite(value) && std::isfinite(weight) && weight>0,
                "nonfinite Dirichlet data or physical face weight");
            mass.noalias()+=weight*basis*basis.transpose();rhs.noalias()+=weight*value*basis;
        }
        Eigen::LLT<Eigen::Matrix3d> factor(mass);
        require(factor.info()==Eigen::Success,"Dirichlet weighted face projection is not positive definite");
        face.coefficients=toArray<3>(factor.solve(rhs));
        report.prescribedSkeletonDofs+=3;
    } else {
        require(report.freeSkeletonDofs<=std::numeric_limits<int>::max()-3,"skeleton index overflow");
        face.offset=report.freeSkeletonDofs;report.freeSkeletonDofs+=3;
    }
    return face.prescribed;
}
struct RecoveredElement
{
    std::array<double,6> cell{};
    std::array<double,10> potential{};
    std::vector<std::array<double,3>> residual,physicalPrimal,physicalScale;
};
using SkeletonFactor = Eigen::SimplicialLDLT<ElectricCartesianResult::Matrix>;
struct SkeletonSolution
{
    ElectricCartesianResult::Matrix matrix;
    std::vector<double> rhs,solution,fixedRhs;
    std::shared_ptr<const SkeletonFactor> factor;
    std::vector<RecoveredElement> elements;
    std::vector<std::array<double,3>> balance,balanceScale;
    double matrixNorm=0;
};
void recoverCondensedSkeleton(std::vector<FaceData>& faces,
    const std::vector<PreparedCell>& local,
    ElectricCartesianReport& report,const Eigen::VectorXd& solution,
    SkeletonSolution& result,const std::vector<std::array<double,3>>& sharedLoads)
{
    const auto recoveryBegin=TimingClock::now();
    require(solution.allFinite(),"nonfinite recovered skeleton solution");
    report.maxRelativeCellResidual=0;report.maxRelativeFreeFaceBalance=0;
    for(auto& face:faces)if(!face.axis && !face.prescribed)
        for(int j=0;j<3;++j)face.coefficients[static_cast<std::size_t>(j)]=solution[face.offset+j];
    result.balance.assign(faces.size(),{});result.balanceScale.assign(faces.size(),{});
    result.elements.clear();
    result.elements.reserve(local.size());
    for(const auto& prepared:local) {
        Eigen::VectorXd coefficients(3*static_cast<int>(prepared.faceRows.size()));
        for(std::size_t f=0;f<prepared.faceRows.size();++f) {
            const auto& face=faces[static_cast<std::size_t>(prepared.faceRows[f])];
            for(int j=0;j<3;++j)coefficients[3*static_cast<Eigen::Index>(f)+j]=face.coefficients[static_cast<std::size_t>(j)];
        }
        const auto cell=prepared.op.recoverCell(coefficients,prepared.cellLoad);
        RecoveredElement record;record.cell=toArray<6>(cell);
        record.potential=toArray<10>(prepared.op.reconstruct(cell,coefficients));
        const auto& stiffness=prepared.op.stiffnessMatrix();
        const auto kcc=stiffness.topLeftCorner(6,6);
        const auto kcf=stiffness.topRightCorner(6,coefficients.size());
        const Eigen::VectorXd cellResidual=kcc*cell+kcf*coefficients-prepared.cellLoad;
        const double cellScale=prepared.cellLoad.norm()+kcc.norm()*cell.norm()+kcf.norm()*coefficients.norm();
        report.maxRelativeCellResidual=std::max(report.maxRelativeCellResidual,
            normalized(cellResidual.norm(),cellScale));
        const auto residual=prepared.op.faceResidual(cell,coefficients,prepared.faceLoad);
        const Eigen::VectorXd cellScaleByRow=stiffness.bottomLeftCorner(coefficients.size(),6).cwiseAbs()*cell.cwiseAbs();
        const Eigen::VectorXd faceScaleByRow=stiffness.bottomRightCorner(coefficients.size(),coefficients.size()).cwiseAbs()*coefficients.cwiseAbs();
        record.residual.resize(prepared.faceRows.size());
        record.physicalPrimal.resize(prepared.faceRows.size());
        record.physicalScale.resize(prepared.faceRows.size());
        for(std::size_t f=0;f<prepared.faceRows.size();++f) {
            const auto row=static_cast<std::size_t>(prepared.faceRows[f]);
            const auto& face=faces[row];
            for(int j=0;j<3;++j) {
                const auto localRow=3*static_cast<Eigen::Index>(f)+j;
                record.residual[f][static_cast<std::size_t>(j)]=residual[localRow];
                // r subtracts physical Neumann data; restore it to recover epsilon grad(V).n.
                record.physicalPrimal[f][static_cast<std::size_t>(j)]
                    =residual[localRow]+prepared.faceLoad[localRow];
                record.physicalScale[f][static_cast<std::size_t>(j)]=cellScaleByRow[localRow]
                    +faceScaleByRow[localRow]+std::abs(prepared.faceLoad[localRow]);
                if(!face.prescribed) {
                    result.balance[row][static_cast<std::size_t>(j)]+=residual[localRow];
                    result.balanceScale[row][static_cast<std::size_t>(j)]+=cellScaleByRow[localRow]
                        +faceScaleByRow[localRow]+std::abs(prepared.faceLoad[localRow]);
                }
            }
        }
        result.elements.push_back(std::move(record));
    }
    for(std::size_t face=0;face<faces.size();++face)if(!faces[face].prescribed && !faces[face].axis)
        for(int j=0;j<3;++j) {
            const double load=sharedLoads.empty()?0:sharedLoads[face][static_cast<std::size_t>(j)];
            report.maxRelativeFreeFaceBalance=std::max(report.maxRelativeFreeFaceBalance,
                normalized(std::abs(result.balance[face][static_cast<std::size_t>(j)]-load),
                    result.balanceScale[face][static_cast<std::size_t>(j)]+std::abs(load)));
        }
    report.timing.recoverySeconds+=elapsedSeconds(recoveryBegin);
}
struct PreparedPhaseBackend
{
    std::vector<PreparedCell> local;
    std::vector<std::vector<std::size_t>> incidence;
    std::vector<double> fixedRhs;
    std::shared_ptr<const SkeletonFactor> factor;
    ElectricCartesianControls controls;
    double matrixNorm=0;
};
// One cell-only condensation, SPD solve and numerical-flux recovery for both
// native Cartesian cells and real Hermite phase elements. It consumes prepared
// geometry/operator records and never discovers, fits or refines an interface.
SkeletonSolution assembleCondensedSkeleton(std::vector<FaceData>& faces,
    const std::vector<PreparedCell>& local,const ElectricCartesianControls& controls,
    ElectricCartesianReport& report,
    const std::vector<std::array<double,3>>& sharedLoads = {})
{
    const auto assemblyBegin=TimingClock::now();
    validateControls(controls);
    require(std::any_of(faces.begin(),faces.end(),[](const FaceData& face){return face.prescribed;}),
        "a positive-measure Dirichlet face must anchor the skeleton");
    const int size=report.freeSkeletonDofs;
    Eigen::VectorXd rhs=Eigen::VectorXd::Zero(size),fixedRhs=Eigen::VectorXd::Zero(size);
    require(sharedLoads.empty() || sharedLoads.size()==faces.size(),"shared load has foreign face layout");
    if(!sharedLoads.empty())for(std::size_t f=0;f<faces.size();++f)
        for(int j=0;j<3;++j) {
            const double value=sharedLoads[f][static_cast<std::size_t>(j)];
            require(std::isfinite(value),"nonfinite shared interface load");
            if(faces[f].prescribed || faces[f].axis)
                require(value==0,"a prescribed/axis trace cannot carry a shared load");
            else rhs[faces[f].offset+j]+=value; // Once per shared trace, never per incident cell.
        }
    std::vector<Eigen::Triplet<double>> triplets;
    std::size_t tripletCount=0;
    for(const auto& prepared:local) {
        std::size_t freeFaces=0;
        for(int row:prepared.faceRows) {
            require(row>=0 && static_cast<std::size_t>(row)<faces.size(),
                "local element has a foreign shared face");
            if(!faces[static_cast<std::size_t>(row)].prescribed)++freeFaces;
        }
        tripletCount+=9*freeFaces*freeFaces;
    }
    triplets.reserve(tripletCount);
    for(const auto& prepared:local) {
        const auto& schur=prepared.op.condensedMatrix();
        const auto& load=prepared.condensedLoad;
        for(std::size_t a=0;a<prepared.faceRows.size();++a) {
            require(prepared.faceRows[a]>=0 && prepared.faceRows[a]<static_cast<int>(faces.size()),
                "local element has a foreign shared face");
            const auto& fa=faces[static_cast<std::size_t>(prepared.faceRows[a])];
            require(!fa.axis,"axis face cannot enter a local skeleton");
            if(fa.prescribed)continue;
            for(int ia=0;ia<3;++ia) {
                const int globalA=fa.offset+ia,localA=3*static_cast<int>(a)+ia;
                require(globalA>=0 && globalA<size,"invalid free skeleton row");
                rhs[globalA]+=load[localA];fixedRhs[globalA]+=load[localA];
                for(std::size_t b=0;b<prepared.faceRows.size();++b) {
                    require(prepared.faceRows[b]>=0 && prepared.faceRows[b]<static_cast<int>(faces.size()),
                        "local element has a foreign shared face");
                    const auto& fb=faces[static_cast<std::size_t>(prepared.faceRows[b])];
                    for(int ib=0;ib<3;++ib) {
                        const double value=schur(localA,3*static_cast<int>(b)+ib);
                        require(std::isfinite(value),"nonfinite local condensed entry");
                        if(fb.prescribed) {
                            rhs[globalA]-=value*fb.coefficients[static_cast<std::size_t>(ib)];
                            fixedRhs[globalA]-=value*fb.coefficients[static_cast<std::size_t>(ib)];
                        }
                        else if(value!=0)triplets.emplace_back(globalA,fb.offset+ib,value);
                    }
                }
            }
        }
    }
    require(rhs.allFinite(),"nonfinite assembled skeleton load");
    SkeletonSolution result;
    result.matrix.resize(size,size);result.matrix.setFromTriplets(triplets.begin(),triplets.end());
    result.matrix.makeCompressed();
    const ElectricCartesianResult::Matrix transpose=result.matrix.transpose();
    const ElectricCartesianResult::Matrix difference=result.matrix-transpose;
    const double matrixNorm=result.matrix.norm();result.matrixNorm=matrixNorm;
    report.relativeSymmetryDefect=normalized(difference.norm(),matrixNorm);
    require(report.relativeSymmetryDefect<=controls.relativeSymmetryTolerance,
        "actual assembled skeleton violates symmetry");
    Eigen::VectorXd solution=Eigen::VectorXd::Zero(size);
    report.timing.assemblySeconds+=elapsedSeconds(assemblyBegin);
    if(size>0) {
        const auto factorBegin=TimingClock::now();
        require(matrixNorm>0 && std::isfinite(matrixNorm),"empty or nonfinite free skeleton operator");
        auto factor=std::make_shared<SkeletonFactor>();
        factor->compute(result.matrix);
        require(factor->info()==Eigen::Success,"anchored skeleton SPD factorization failed");
        const auto pivots=factor->vectorD();
        require(pivots.allFinite() && (pivots.array()>0).all(),"skeleton has a nonpositive LDLT pivot");
        report.minimumSkeletonPivot=pivots.minCoeff();
        report.timing.factorSeconds+=elapsedSeconds(factorBegin);
        const auto solveBegin=TimingClock::now();
        solution=factor->solve(rhs);
        require(factor->info()==Eigen::Success && solution.allFinite(),"skeleton SPD solve failed");
        const Eigen::VectorXd residual=result.matrix*solution-rhs;
        report.relativeLinearResidual=normalized(residual.norm(),rhs.norm()+matrixNorm*solution.norm());
        require(report.relativeLinearResidual<=controls.relativeLinearTolerance,
            "skeleton residual exceeds its numerical solve tolerance");
        report.timing.solveSeconds+=elapsedSeconds(solveBegin);
        result.factor=std::move(factor);
    }
    if(size>0) {
        result.fixedRhs.assign(fixedRhs.data(),fixedRhs.data()+fixedRhs.size());
        result.rhs.assign(rhs.data(),rhs.data()+rhs.size());
        result.solution.assign(solution.data(),solution.data()+solution.size());
    }
    recoverCondensedSkeleton(faces,local,report,solution,result,sharedLoads);
    return result;
}

}
struct ElectricCartesianResult::Data
{
    explicit Data(const QuadtreeMesh& mesh,bool rz)
        :identity(mesh),axisymmetric(rz),verticalFaces(mesh.vEdgeIndexNum()),
         cells(mesh.cellIndexNum()),faces(mesh.vEdgeIndexNum()+mesh.hEdgeIndexNum()){}
    NativeMeshIdentity identity;
    bool axisymmetric;
    int verticalFaces;
    ElectricCartesianReport report;
    std::vector<CellData> cells;
    std::vector<FaceData> faces;
    Matrix skeleton;
    std::vector<double> rhs,solution;
    const CellData& cell(int index) const {
        require(index>=0 && index<static_cast<int>(cells.size()),"invalid native cell index");
        return cells[static_cast<std::size_t>(index)];
    }
    const FaceData& face(bool horizontal,int index) const {
        const int count=horizontal?static_cast<int>(faces.size())-verticalFaces:verticalFaces;
        require(index>=0 && index<count,"invalid native edge index");
        return faces[static_cast<std::size_t>((horizontal?verticalFaces:0)+index)];
    }
};
ElectricCartesianResult::ElectricCartesianResult(std::shared_ptr<const Data> data):data_(std::move(data)){}
const ElectricCartesianResult::Data& ElectricCartesianResult::checked() const
{
    if(!data_)throw std::logic_error("[ElectricPotential] moved-from result");
    return *data_;
}
const ElectricCartesianReport& ElectricCartesianResult::report() const{return checked().report;}
bool ElectricCartesianResult::axisymmetric() const{return checked().axisymmetric;}
void ElectricCartesianResult::validateMesh(const QuadtreeMesh& mesh) const{checked().identity.validate(mesh);}
std::size_t ElectricCartesianResult::nativeCellCount() const{return checked().cells.size();}
Vector2D ElectricCartesianResult::cellCenter(int index) const{return checked().cell(index).center;}
double ElectricCartesianResult::cellScale(int index) const{return checked().cell(index).scale;}
std::array<double,6> ElectricCartesianResult::cellPolynomialP2(int index) const{return checked().cell(index).cell;}
std::array<double,10> ElectricCartesianResult::reconstructedPotentialP3(int index) const{return checked().cell(index).potential;}
double ElectricCartesianResult::potentialInCell(int index,Vector2D point) const
{
    const auto& c=checked().cell(index);requirePoint(c,point);
    const auto basis=Electric::p3Basis(point,c.center,c.scale);
    long double value=0;
    for(int j=0;j<10;++j)value+=static_cast<long double>(basis[j])*c.potential[j];
    const double result=static_cast<double>(value);
    require(std::isfinite(result),"nonfinite retained potential query");return result;
}
Vector2D ElectricCartesianResult::gradientInCell(int index,Vector2D point) const
{
    const auto& c=checked().cell(index);requirePoint(c,point);
    const auto basis=Electric::p3Gradient(point,c.center,c.scale);
    long double dx=0,dy=0;
    for(int j=0;j<10;++j) {
        dx+=static_cast<long double>(basis(0,j))*c.potential[j];
        dy+=static_cast<long double>(basis(1,j))*c.potential[j];
    }
    const Vector2D result{static_cast<double>(dx),static_cast<double>(dy)};
    require(std::isfinite(result.x) && std::isfinite(result.y),"nonfinite retained gradient query");return result;
}
std::array<Vector2D,2> ElectricCartesianResult::nativeFaceBounds(bool horizontal,int index) const
{return checked().face(horizontal,index).bounds;}
std::array<double,3> ElectricCartesianResult::nativeFaceTrace(bool horizontal,int index) const
{
    const auto& f=checked().face(horizontal,index);
    require(!f.axis,"zero-area axis has no trace polynomial");return f.coefficients;
}
std::array<double,3> ElectricCartesianResult::outwardFaceResidual(int index,int side) const
{
    const auto& c=checked().cell(index);
    require(side>=0 && side<4 && c.localFace[static_cast<std::size_t>(side)]>=0,
        "invalid side or zero-area axis has no face row");
    return c.residual[static_cast<std::size_t>(side)];
}
const ElectricCartesianResult::Matrix& ElectricCartesianResult::skeletonMatrix() const{return checked().skeleton;}
const std::vector<double>& ElectricCartesianResult::skeletonRhs() const{return checked().rhs;}
const std::vector<double>& ElectricCartesianResult::skeletonSolution() const{return checked().solution;}

ElectricCartesianResult solveElectricPotential(const QuadtreeMesh& mesh,
    const ElectricCartesianProblem& problem,const ElectricCartesianControls& controls)
{
    const auto boundaryBegin=TimingClock::now();
    validateControls(controls);
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    require(problem.epsilonByCell.size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "permittivity array must use the actual native cell index space");
    for(double epsilon:problem.epsilonByCell)
        require(std::isfinite(epsilon) && epsilon>0,"permittivity must be finite and positive");
    require(!problem.axisymmetric || mesh.minCoord.x>=0,"RZ domain has a negative radius");
    validateBoundaries(problem.boundaries);
    auto data=std::make_shared<ElectricCartesianResult::Data>(mesh,problem.axisymmetric);
    data->report.cells=data->cells.size();
    bool anchored=false;
    for(bool horizontal:{false,true}) {
        const auto& labels=horizontal?grid.horizontal():grid.vertical();
        for(int label:labels) {
            const auto& edge=horizontal?mesh.hEdges(label):mesh.vEdges(label);
            const int row=(horizontal?data->verticalFaces:0)+edge.index();
            auto& face=data->faces[static_cast<std::size_t>(row)];
            face.edgeLabel=label;face.horizontal=horizontal;face.boundary=boundarySide(edge.boundary());
            face.bounds={{mesh.verts(edge.edgeToVert[0]).coord(),mesh.verts(edge.edgeToVert[1]).coord()}};
            const auto a=face.bounds[0],b=face.bounds[1];
            const double length=horizontal?b.x-a.x:b.y-a.y;
            require(std::isfinite(length) && length>0,"invalid actual face interval");
            face.axis=problem.axisymmetric && !horizontal && a.x==0 && b.x==0;
            if(face.axis) {
                require(face.boundary==0 && problem.boundaries[0].type==ElectricCartesianBoundary::Type::Neumann,
                    "RZ axis must be explicit homogeneous Neumann");
                for(double q:gaussPoints) {
                    const double value=problem.boundaries[0].value({0,coordinate(a.y,b.y,q)});
                    require(std::isfinite(value) && value==0,"RZ axis Neumann provider is not identically zero at its rule");
                }
                ++data->report.axisFaces;continue;
            }
            ++data->report.physicalFaces;
            for(std::size_t q=0;q<gaussPoints.size();++q) {
                const Vector2D point=horizontal?Vector2D{coordinate(a.x,b.x,gaussPoints[q]),a.y}
                    :Vector2D{a.x,coordinate(a.y,b.y,gaussPoints[q])};
                const double weight=length*gaussWeights[q]/2;
                require(std::isfinite(weight) && weight>0,"unrepresentable physical face rule");
                face.nodes.push_back({point,{0,0},weight,gaussPoints[q]});
            }
            anchored=prepareFaceUnknown(face,problem.boundaries,problem.axisymmetric,data->report)||anchored;
        }
    }
    require(anchored,"a positive-measure Dirichlet face must anchor the skeleton");
    data->report.timing.boundarySeconds=elapsedSeconds(boundaryBegin);
    std::vector<PreparedCell> local;local.reserve(data->cells.size());
    // local storage follows the grid labels, while every result and material
    // record uses the actual native index. No identity is inferred from order.
    std::vector<int> localCellIndices;localCellIndices.reserve(data->cells.size());
    for(int label:grid.cells()) {
        const auto localBegin=TimingClock::now();
        const auto& cell=mesh.cells(label);const int native=cell.index();
        auto& record=data->cells[static_cast<std::size_t>(native)];record.label=label;
        record.bounds={{mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()}};
        const auto lower=record.bounds[0],upper=record.bounds[1];
        const double hx=upper.x-lower.x,hy=upper.y-lower.y;
        record.center=cell.centerCoord();record.scale=std::hypot(hx,hy);
        require(std::isfinite(record.scale) && hx>0 && hy>0 && record.scale>0,"unresolved native cell geometry");
        Electric::LocalGeometry geometry;geometry.center=record.center;geometry.scale=record.scale;
        geometry.axisymmetric=problem.axisymmetric;
        for(std::size_t j=0;j<gaussPoints.size();++j)for(std::size_t i=0;i<gaussPoints.size();++i) {
            const double weight=(hx*gaussWeights[i]/2)*(hy*gaussWeights[j]/2);
            require(std::isfinite(weight) && weight>0,"unrepresentable native volume rule");
            geometry.volume.push_back({{coordinate(lower.x,upper.x,gaussPoints[i]),
                coordinate(lower.y,upper.y,gaussPoints[j])},weight});
        }
        record.faceRows={{mesh.vEdges(cell.cellToVEdge[0][0]).index(),
            data->verticalFaces+mesh.hEdges(cell.cellToHEdge[0][0]).index(),
            mesh.vEdges(cell.cellToVEdge[1][0]).index(),
            data->verticalFaces+mesh.hEdges(cell.cellToHEdge[1][0]).index()}};
        std::vector<int> rows;
        for(int side=0;side<4;++side) {
            const int row=record.faceRows[static_cast<std::size_t>(side)];
            const auto& face=data->faces[static_cast<std::size_t>(row)];
            if(face.axis)continue;
            Electric::FaceQuadrature quadrature;
            quadrature.length=face.horizontal?face.bounds[1].x-face.bounds[0].x:face.bounds[1].y-face.bounds[0].y;
            quadrature.nodes=face.nodes;
            for(auto& node:quadrature.nodes)node.outwardNormal=sideNormal(side);
            record.localFace[static_cast<std::size_t>(side)]=static_cast<int>(rows.size());
            geometry.faces.push_back(std::move(quadrature));rows.push_back(row);
        }
        Electric::LocalOperator op(std::move(geometry),problem.epsilonByCell[static_cast<std::size_t>(native)]);
        data->report.timing.localOperatorsSeconds+=elapsedSeconds(localBegin);
        const auto loadBegin=TimingClock::now();
        const Electric::CellVector cellLoad=problem.source?op.loadCell(problem.source):Electric::CellVector::Zero().eval();
        Eigen::VectorXd faceLoad=Eigen::VectorXd::Zero(3*static_cast<int>(rows.size()));
        for(std::size_t j=0;j<rows.size();++j) {
            const auto& face=data->faces[static_cast<std::size_t>(rows[j])];
            if(face.boundary>=0 && !face.prescribed)
                faceLoad.segment<3>(3*static_cast<Eigen::Index>(j))=op.loadFace(j,
                    problem.boundaries[static_cast<std::size_t>(face.boundary)].value);
        }
        local.emplace_back(std::move(op),cellLoad,std::move(faceLoad),std::move(rows));
        localCellIndices.push_back(native);
        data->report.timing.localLoadsSeconds+=elapsedSeconds(loadBegin);
    }
    auto solved=assembleCondensedSkeleton(data->faces,local,controls,data->report);
    data->skeleton=std::move(solved.matrix);
    data->rhs=std::move(solved.rhs);data->solution=std::move(solved.solution);
    for(std::size_t ci=0;ci<local.size();++ci) {
        auto& record=data->cells[static_cast<std::size_t>(localCellIndices[ci])];
        const auto& recovered=solved.elements[ci];
        record.cell=recovered.cell;record.potential=recovered.potential;
        for(int side=0;side<4;++side) {
            const int face=record.localFace[static_cast<std::size_t>(side)];
            if(face>=0)record.residual[static_cast<std::size_t>(side)]=recovered.residual[static_cast<std::size_t>(face)];
        }
    }
    data->identity.validate(mesh);
    return ElectricCartesianResult(std::move(data));
}

namespace
{
bool sameClock(const Interface::HermiteGalsTimeLayer& a,const Interface::HermiteGalsTimeLayer& b)
{
    return a.step==b.step && a.begin==b.begin && a.end==b.end && a.position==b.position;
}
bool sameGraph(const LevelSetMappedInterface::Graph& a,const LevelSetMappedInterface::Graph& b)
{
    return a.nativeCellIndex==b.nativeCellIndex && a.cellI==b.cellI && a.cellJ==b.cellJ
        && a.coordinateIsX==b.coordinateIsX
        && a.begin.side==b.begin.side && a.begin.rootOrdinal==b.begin.rootOrdinal
        && a.end.side==b.end.side && a.end.rootOrdinal==b.end.rootOrdinal
        && a.localLower.x==b.localLower.x && a.localLower.y==b.localLower.y
        && a.localUpper.x==b.localUpper.x && a.localUpper.y==b.localUpper.y
        && a.controlledIntegration==b.controlledIntegration && !a.nativeSide && !b.nativeSide;
}
using ChargeMoments = Eigen::Matrix<double,2,3>;
void addChargeMoments(std::array<long double,6>& sums,const SurfaceChargeDensityStencil& stencil,
    double s,double zeta,double physicalWeight)
{
    require(std::isfinite(physicalWeight) && physicalWeight>0
        && std::isfinite(zeta) && zeta>-1 && zeta<1,"invalid charge moment integration node");
    const auto weights=stencil.densityWeightsAtArclength(s);
    const auto basis=Electric::faceBasis(zeta);
    for(int mode=0;mode<2;++mode)for(int col=0;col<3;++col)
        sums[static_cast<std::size_t>(3*mode+col)]+=static_cast<long double>(physicalWeight)
            *basis[mode+1]*weights[static_cast<std::size_t>(col)];
}
ChargeMoments chargeMoments(const std::array<long double,6>& sums)
{
    ChargeMoments result;
    for(int mode=0;mode<2;++mode)for(int col=0;col<3;++col)
        result(mode,col)=static_cast<double>(sums[static_cast<std::size_t>(3*mode+col)]);
    require(result.allFinite(),"unrepresentable surface charge moment coefficients");
    return result;
}
ChargeMoments retainedChargeMoments(const ElectricPhaseSkeletonFace& face,
    const SurfaceChargeDensityStencil& stencil,bool rz)
{
    std::array<long double,6> sums{};
    for(const auto& node:face.nodes)
        addChargeMoments(sums,stencil,(node.coordinate+1)*face.length/2,node.coordinate,
            node.weight*Electric::physicalWeight(node.point,rz));
    return chargeMoments(sums);
}
struct RefinedChargeMoments { ChargeMoments moments; double length=0; };
RefinedChargeMoments refinedChargeMoments(const ElectricPhaseGeometry& geometry,int component,int graph,
    const SurfaceChargeDensityStencil& stencil,double traceLength,int refinement)
{
    require(refinement>=1 && refinement<=12,"charge quadrature refinement exceeds its bounded rule");
    const auto& publication=geometry.publication();
    const auto& chain=publication.topology().curves[static_cast<std::size_t>(component)];
    const int intervals=1<<refinement;
    long double complete=0;
    std::array<long double,6> sums{};
    for(int interval=0;interval<intervals;++interval) {
        const double begin=static_cast<double>(interval)/intervals;
        const double end=static_cast<double>(interval+1)/intervals;
        long double intervalLength=0;
        for(const auto& q:publication.quadrature(static_cast<std::size_t>(component),
                static_cast<std::size_t>(graph),begin,end)) {
            require(std::isfinite(q.weight) && q.weight>0,"invalid true-arc charge quadrature weight");
            const double fraction=chain.graphFraction(graph,q.sample.frame.cellChartPoint);
            require(std::isfinite(fraction) && fraction>begin && fraction<end,
                "unresolved charge quadrature sample in its original graph chart");
            long double prefix=0;
            for(const auto& part:publication.quadrature(static_cast<std::size_t>(component),
                    static_cast<std::size_t>(graph),begin,fraction)) {
                require(std::isfinite(part.weight) && part.weight>0,"invalid true-arc charge prefix");
                prefix+=part.weight;
            }
            const double arc=static_cast<double>(complete+prefix);
            // The trace retains the prepared true length; no density/arc stretching is applied.
            const double zeta=static_cast<double>(2*(complete+prefix)/traceLength-1);
            addChargeMoments(sums,stencil,arc,zeta,
                q.weight*Electric::physicalWeight(q.sample.frame.point,geometry.axisymmetric()));
            intervalLength+=q.weight;
        }
        require(intervalLength>0 && std::isfinite(intervalLength),"unresolved charge integration subarc");
        complete+=intervalLength;
    }
    const double length=static_cast<double>(complete);
    require(std::isfinite(length) && length>0,"unrepresentable true-arc charge panel measure");
    return {chargeMoments(sums),length};
}
}
struct ElectricSurfaceChargeMomentOperator::Data
{
    explicit Data(ElectricPhaseGeometry source,Interface::HermiteGalsTimeLayer layer)
        :geometry(std::move(source)),timeLayer(layer){}
    struct Panel { int component=-1,graph=-1,face=-1; };
    ElectricPhaseGeometry geometry;
    Interface::HermiteGalsTimeLayer timeLayer;
    ElectricSurfaceChargeMomentReport report;
    std::vector<Panel> panels;
    std::vector<const SurfaceChargeField*> fields;
    // Moving P_Q retains the actual immutable component geometry/Q owners.
    // All loads, potential/current copies and solver factors retain this P_Q.
    std::shared_ptr<const void> fieldsLifetime;
    Matrix projection,injection;
    const Panel& panel(std::size_t index) const
    {
        require(index<panels.size(),"invalid bound charge panel index");
        return panels[index];
    }
};
ElectricSurfaceChargeMomentOperator::ElectricSurfaceChargeMomentOperator(std::shared_ptr<const Data> data)
    :data_(std::move(data)){}
const ElectricSurfaceChargeMomentOperator::Data& ElectricSurfaceChargeMomentOperator::checked() const
{
    if(!data_)throw std::logic_error("[ElectricSurfaceChargeMomentOperator] moved-from operator");
    return *data_;
}
const ElectricPhaseGeometry& ElectricSurfaceChargeMomentOperator::geometry() const{return checked().geometry;}
Interface::HermiteGalsTimeLayer ElectricSurfaceChargeMomentOperator::timeLayer() const{return checked().timeLayer;}
const ElectricSurfaceChargeMomentReport& ElectricSurfaceChargeMomentOperator::report() const{return checked().report;}
std::size_t ElectricSurfaceChargeMomentOperator::panelCount() const{return checked().panels.size();}
int ElectricSurfaceChargeMomentOperator::interfaceFace(std::size_t index) const{return checked().panel(index).face;}
int ElectricSurfaceChargeMomentOperator::component(std::size_t index) const{return checked().panel(index).component;}
int ElectricSurfaceChargeMomentOperator::graph(std::size_t index) const{return checked().panel(index).graph;}
const ElectricSurfaceChargeMomentOperator::Matrix& ElectricSurfaceChargeMomentOperator::projectionMatrix() const
{return checked().projection;}
const ElectricSurfaceChargeMomentOperator::Matrix& ElectricSurfaceChargeMomentOperator::interfaceInjectionMatrix() const
{return checked().injection;}
ElectricSurfaceChargeMomentOperator ElectricSurfaceChargeMomentOperator::prepare(const ElectricPhaseGeometry& geometry,
    const std::vector<const SurfaceChargeField*>& fields,
    const Interface::HermiteGalsTimeLayer& coevalChargeLayer,
    const ElectricSurfaceChargeMomentControls& controls)
{
    return prepareRetained(geometry,fields,coevalChargeLayer,controls,{});
}
ElectricSurfaceChargeMomentOperator ElectricSurfaceChargeMomentOperator::prepareRetained(const ElectricPhaseGeometry& geometry,
    const std::vector<const SurfaceChargeField*>& fields,
    const Interface::HermiteGalsTimeLayer& coevalChargeLayer,
    const ElectricSurfaceChargeMomentControls& controls,std::shared_ptr<const void> fieldsLifetime)
{
    require(std::isfinite(controls.relativeQuadratureTolerance)
        && controls.relativeQuadratureTolerance>0 && controls.relativeQuadratureTolerance<=1e-9
        && controls.maximumAdditionalRefinements>=1 && controls.maximumAdditionalRefinements<=8,
        "invalid surface charge quadrature controls");
    const auto layer=geometry.timeLayer();
    (void)layer.time();(void)coevalChargeLayer.time();
    require(sameClock(layer,coevalChargeLayer),"declared charge clock differs from publication full clock");
    const auto& publication=geometry.publication();
    const auto& curves=publication.topology().curves;
    const auto& snapshot=publication.state().geometry().snapshot();
    require(fields.size()==curves.size(),"charge fields must cover publication components in their actual order");
    auto data=std::make_shared<Data>(geometry,coevalChargeLayer);
    data->fields=fields;
    data->fieldsLifetime=std::move(fieldsLifetime);
    std::vector<std::size_t> offsets(fields.size()+1,0);
    for(std::size_t c=0;c<fields.size();++c) {
        require(fields[c]!=nullptr,"null component charge owner");
        const auto& panels=fields[c]->geometry();
        require(panels.usesNativeLeaves() && panels.axisymmetric()==geometry.axisymmetric(),
            "charge owner has foreign panel layout or physical measure");
        const auto& native=panels.nativeChain();const auto& curve=curves[c];
        require(native.usesSource(snapshot) && curve.usesSource(snapshot),
            "charge owner is not bound to the publication's actual scalar pointer");
        require(native.closed()==curve.closed() && native.axisymmetric()==curve.axisymmetric()
            && native.beginEndpointKind()==curve.beginEndpointKind()
            && native.endEndpointKind()==curve.endEndpointKind()
            && native.negativeRegion()==curve.negativeRegion()
            && native.graphs().size()==curve.graphs().size()
            && panels.panels().size()==curve.graphs().size()
            && fields[c]->transport().panelCharge().size()==curve.graphs().size(),
            "charge component topology differs from the publication");
        for(std::size_t g=0;g<curve.graphs().size();++g) {
            require(sameGraph(native.graphs()[g],curve.graphs()[g]),
                "charge panel graph identity or oriented endpoint roots differ from the publication");
            const auto& panel=panels.panels()[g];
            require(panel.cellIndex==curve.graphs()[g].nativeCellIndex
                && panel.coordinateIsX==curve.graphs()[g].coordinateIsX
                && std::isfinite(panel.measure) && panel.measure>0,
                "charge material panel does not own its original graph");
        }
        offsets[c+1]=offsets[c]+panels.panels().size();
    }
    const std::size_t count=offsets.back();
    require(count>0 && count<=static_cast<std::size_t>(std::numeric_limits<int>::max()/3)
        && geometry.faces().size()<=static_cast<std::size_t>(std::numeric_limits<int>::max()/3),
        "surface charge moment index overflow or empty panel set");
    data->panels.resize(count);
    std::vector<bool> assigned(count,false);
    for(std::size_t face=0;face<geometry.faces().size();++face) {
        const auto& retained=geometry.faces()[face];
        if(retained.kind!=ElectricPhaseSkeletonFace::Kind::Interface)continue;
        require(retained.component>=0 && static_cast<std::size_t>(retained.component)<fields.size(),
            "interface trace has foreign charge component");
        const auto c=static_cast<std::size_t>(retained.component);
        require(retained.graph>=0 && static_cast<std::size_t>(retained.graph)<offsets[c+1]-offsets[c],
            "interface trace has foreign charge graph");
        const auto flat=offsets[c]+static_cast<std::size_t>(retained.graph);
        require(!assigned[flat],"a material charge panel maps to multiple interface traces");
        assigned[flat]=true;
        data->panels[flat]={retained.component,retained.graph,static_cast<int>(face)};
    }
    require(std::all_of(assigned.begin(),assigned.end(),[](bool value){return value;}),
        "a material charge panel has no actual shared interface trace");
    std::vector<Eigen::Triplet<double>> projection,injection;
    const int baseRefinement=geometry.report().acceptedRefinement;
    require(baseRefinement>=0 && baseRefinement+controls.maximumAdditionalRefinements<=12,
        "prepared geometry leaves no bounded charge quadrature refinement");
    struct PanelPreparation {
        ChargeMoments moments;
        std::array<std::size_t,3> neighbor{};
        int additionalRefinements=0;
        double momentDifference=0,lengthDifference=0;
    };
    std::vector<PanelPreparation> preparations(count);
    std::vector<std::exception_ptr> failures(count);
    // Each panel reads immutable geometry and its density descriptor. Keep
    // every panel's original quadrature and accumulation order; only panels
    // execute concurrently. Matrix entries and diagnostics merge below in
    // the original flat-panel order, and no exception crosses the OpenMP region.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(count>1)
#endif
    for(std::ptrdiff_t flatIndex=0;flatIndex<static_cast<std::ptrdiff_t>(count);++flatIndex) {
        const auto flat=static_cast<std::size_t>(flatIndex);
        try {
            const auto& binding=data->panels[flat];
            const auto c=static_cast<std::size_t>(binding.component);
            const auto g=static_cast<std::size_t>(binding.graph);
            const auto& panel=fields[c]->geometry().panels()[g];
            const auto& face=geometry.faces()[static_cast<std::size_t>(binding.face)];
            const auto stencil=fields[c]->densityStencil(g);
            require(stencil.panelMeasure==panel.measure && std::isfinite(face.length) && face.length>0,
                "charge density descriptor differs from its original physical panel");
            auto& ready=preparations[flat];ready.neighbor=stencil.neighbor;
            ready.lengthDifference=normalized(std::abs(panel.measure-face.length),panel.measure+face.length);
            require(ready.lengthDifference<=controls.relativeQuadratureTolerance,
                "material panel and electric trace true arclength measures differ");
            auto coarse=retainedChargeMoments(face,stencil,geometry.axisymmetric());
            bool converged=false;
            for(int additional=1;additional<=controls.maximumAdditionalRefinements;++additional) {
                const auto fine=refinedChargeMoments(geometry,binding.component,binding.graph,stencil,
                    face.length,baseRefinement+additional);
                const double lengthError=normalized(std::abs(fine.length-face.length),fine.length+face.length);
                require(lengthError<=controls.relativeQuadratureTolerance,
                    "refined charge integral has a different true trace length");
                ready.lengthDifference=std::max(ready.lengthDifference,lengthError);
                const double error=normalized((fine.moments-coarse).norm(),fine.moments.norm()+coarse.norm());
                if(error<=controls.relativeQuadratureTolerance) {
                    ready.moments=fine.moments;converged=true;
                    ready.additionalRefinements=additional;ready.momentDifference=error;
                    break;
                }
                coarse=fine.moments;
            }
            require(converged,"true-arc higher charge moments did not reach their quadrature tolerance");
            for(std::size_t col=0;col<3;++col)
                require(ready.neighbor[col]<offsets[c+1]-offsets[c],
                    "charge density stencil references a foreign component panel");
        } catch(...) {failures[flat]=std::current_exception();}
    }
    projection.reserve(7*count);injection.reserve(3*count);
    for(std::size_t flat=0;flat<count;++flat) {
        if(failures[flat])std::rethrow_exception(failures[flat]);
        const auto& binding=data->panels[flat];
        const auto c=static_cast<std::size_t>(binding.component);
        const auto& ready=preparations[flat];
        data->report.additionalRefinements=std::max(data->report.additionalRefinements,ready.additionalRefinements);
        data->report.maxRelativeMomentDifference=std::max(data->report.maxRelativeMomentDifference,ready.momentDifference);
        data->report.maxRelativePanelLengthDifference=std::max(
            data->report.maxRelativePanelLengthDifference,ready.lengthDifference);
        const int row=3*static_cast<int>(flat);
        projection.emplace_back(row,static_cast<int>(flat),1.0);
        for(int mode=0;mode<2;++mode)for(int col=0;col<3;++col) {
            const int column=static_cast<int>(offsets[c]+ready.neighbor[static_cast<std::size_t>(col)]);
            if(ready.moments(mode,col)!=0)projection.emplace_back(row+mode+1,column,ready.moments(mode,col));
        }
        for(int mode=0;mode<3;++mode)injection.emplace_back(3*binding.face+mode,row+mode,1.0);
    }
    data->projection.resize(3*static_cast<int>(count),static_cast<int>(count));
    data->projection.setFromTriplets(projection.begin(),projection.end());data->projection.makeCompressed();
    data->injection.resize(3*static_cast<int>(geometry.faces().size()),3*static_cast<int>(count));
    data->injection.setFromTriplets(injection.begin(),injection.end());data->injection.makeCompressed();
    data->report.panels=count;
    return ElectricSurfaceChargeMomentOperator(std::move(data));
}
struct ElectricPanelChargeLoad::Data
{
    explicit Data(ElectricSurfaceChargeMomentOperator source):op(std::move(source)){}
    ElectricSurfaceChargeMomentOperator op;
    std::vector<double> charge;
    std::vector<std::array<double,3>> moments,faceMoments;
};
ElectricPanelChargeLoad::ElectricPanelChargeLoad(std::shared_ptr<const Data> data):data_(std::move(data)){}
const ElectricPanelChargeLoad::Data& ElectricPanelChargeLoad::checked() const
{
    if(!data_)throw std::logic_error("[ElectricPanelChargeLoad] moved-from load");
    return *data_;
}
const ElectricSurfaceChargeMomentOperator& ElectricPanelChargeLoad::momentOperator() const{return checked().op;}
const std::vector<double>& ElectricPanelChargeLoad::panelCharge() const{return checked().charge;}
const std::vector<std::array<double,3>>& ElectricPanelChargeLoad::moments() const{return checked().moments;}
const std::vector<std::array<double,3>>& ElectricPanelChargeLoad::faceMoments() const{return checked().faceMoments;}
void ElectricPanelChargeLoad::validateGeometry(const ElectricPhaseGeometry& candidate) const
{
    const auto& bound=checked().op.geometry();
    require(bound.geometryIdentity()==candidate.geometryIdentity()
        && sameClock(bound.timeLayer(),candidate.timeLayer())
        && &bound.publication().state().geometry().snapshot()
            ==&candidate.publication().state().geometry().snapshot()
        && &bound.faces()==&candidate.faces() && &bound.elements()==&candidate.elements(),
        "charge moments require their exact immutable prepared geometry/basis owner");
}
ElectricPanelChargeLoad ElectricSurfaceChargeMomentOperator::apply(const std::vector<double>& charge) const
{
    const auto& owner=checked();
    require(charge.size()==owner.panels.size(),"trial charge has foreign material panel layout");
    for(double value:charge)require(std::isfinite(value),"nonfinite trial integrated charge");
    auto data=std::make_shared<ElectricPanelChargeLoad::Data>(*this);
    data->charge=charge;data->moments.resize(charge.size());
    data->faceMoments.resize(owner.geometry.faces().size());
    const Eigen::Map<const Eigen::VectorXd> q(charge.data(),static_cast<Eigen::Index>(charge.size()));
    const Eigen::VectorXd projected=owner.projection*q;
    require(projected.allFinite(),"unrepresentable trial charge moments");
    for(std::size_t i=0;i<charge.size();++i) {
        data->moments[i]={{charge[i],projected[3*static_cast<Eigen::Index>(i)+1],
            projected[3*static_cast<Eigen::Index>(i)+2]}}; // eta0 copies Q bit-for-bit.
        data->faceMoments[static_cast<std::size_t>(owner.panels[i].face)]=data->moments[i];
    }
    return ElectricPanelChargeLoad(std::move(data));
}

struct ElectricPhaseResult::Data
{
    explicit Data(ElectricPhaseGeometry source):geometry(std::move(source)),solutionLayer(geometry.timeLayer()){}
    ElectricPhaseGeometry geometry;
    Interface::HermiteGalsTimeLayer solutionLayer;
    std::shared_ptr<const PreparedPhaseBackend> prepared;
    ElectricPhaseReport report;
    std::vector<FaceData> faces;
    std::vector<RecoveredElement> elements;
    std::shared_ptr<const ElectricPanelChargeLoad> chargeLoad;
    struct InterfaceDisplacement {
        bool available=false;
        std::array<std::array<double,3>,2> moments{},polynomial{};
        std::shared_ptr<const Eigen::LLT<Eigen::Matrix3d>> factor;
    };
    std::vector<InterfaceDisplacement> displacement;
    std::array<double,2> epsilonByPhase{};
    // Constant Gamma residual rows: r_phase = R_phase*u + a_phase.
    // Row identity is retained independently of the caller's material Q order.
    std::vector<int> interfaceRow;
    std::array<Matrix,2> residualMap;
    std::array<std::vector<double>,2> residualOffset;
    std::vector<double> bulkLoad0; // Actual source moments; no bulk charge state.
    std::array<double,2> exteriorDisplacement{},currentEquationScale{};
    double relativeLinearTolerance=0,relativeConservationTolerance=0;
    std::shared_ptr<const Matrix> skeleton;
    std::vector<double> rhs,solution;

    void finishRecovery(const std::vector<PreparedCell>& local,
        const std::vector<std::vector<std::size_t>>& incidence,SkeletonSolution& solved,
        const ElectricPanelChargeLoad* charge,const ElectricCartesianControls& controls,bool prepareAffine)
    {
        const auto gaussAffineBegin=TimingClock::now();
        this->report.maxRelativeInterfaceBalance=0;
        this->report.maxRelativeLocalGaussResidual=0;
        const std::vector<std::array<double,3>> noCharge;
        const auto& sharedLoads=charge?charge->faceMoments():noCharge;
        this->displacement.resize(geometry.faces().size());
        int interfaceRows=static_cast<int>(this->residualMap[0].rows());
        if(prepareAffine) {
            this->interfaceRow.assign(geometry.faces().size(),-1);interfaceRows=0;
            for(std::size_t face=0;face<geometry.faces().size();++face)
                if(geometry.faces()[face].kind==ElectricPhaseSkeletonFace::Kind::Interface)
                    this->interfaceRow[face]=interfaceRows++;
        }
        std::array<std::vector<Eigen::Triplet<double>>,2> residualEntries;
        if(prepareAffine)for(auto& offset:this->residualOffset)offset.assign(static_cast<std::size_t>(interfaceRows),0);
        for(std::size_t face=0;face<geometry.faces().size();++face) {
            if(geometry.faces()[face].kind!=ElectricPhaseSkeletonFace::Kind::Interface)continue;
            const auto eta=charge?sharedLoads[face]:std::array<double,3>{{0,0,0}};
            for(int moment=0;moment<3;++moment) {
                const auto j=static_cast<std::size_t>(moment);
                this->report.maxRelativeInterfaceBalance=std::max(this->report.maxRelativeInterfaceBalance,
                    normalized(std::abs(solved.balance[face][j]-eta[j]),
                        solved.balanceScale[face][j]+std::abs(eta[j])));
            }
            auto& trace=this->displacement[face];trace.available=true;
            if(prepareAffine) {
                Eigen::Matrix3d mass=Eigen::Matrix3d::Zero();
                for(const auto& node:geometry.faces()[face].nodes) {
                    const auto basis=Electric::faceBasis(node.coordinate);
                    const double weight=node.weight*Electric::physicalWeight(node.point,geometry.axisymmetric());
                    require(std::isfinite(weight) && weight>0,"invalid displacement face mass weight");
                    mass.noalias()+=weight*basis*basis.transpose();
                }
                auto factor=std::make_shared<Eigen::LLT<Eigen::Matrix3d>>(mass);
                require(factor->info()==Eigen::Success,"interface displacement mass is not positive definite");
                trace.factor=std::move(factor);
            }
            require(static_cast<bool>(trace.factor),"missing prepared displacement face factor");
            for(std::size_t owner:incidence[face]) {
                const auto& element=geometry.elements()[owner];
                const auto found=std::find(element.faces.begin(),element.faces.end(),static_cast<int>(face));
                require(found!=element.faces.end(),"interface incidence lost its original local face row");
                const std::size_t localFace=static_cast<std::size_t>(found-element.faces.begin());
                const auto& raw=solved.elements[owner].residual[localFace];
                const std::size_t phase=element.liquid?0:1;
                Electric::FacePolynomial moments;
                for(int mode=0;mode<3;++mode)
                    moments[mode]=(element.liquid?-1.:1.)*raw[static_cast<std::size_t>(mode)];
                trace.moments[phase]=toArray<3>(moments);
                trace.polynomial[phase]=toArray<3>(trace.factor->solve(moments));
                if(prepareAffine) {
                    const auto& prepared=local[owner];
                    require(prepared.op.permittivity()==this->epsilonByPhase[phase],
                        "retained material permittivity differs from the actual local electric operator");
                    const auto& schur=prepared.op.condensedMatrix();
                    const Eigen::Index localRow=3*static_cast<Eigen::Index>(localFace);
                    const int currentRow=this->interfaceRow[face];
                    long double fixed=-prepared.condensedLoad[localRow];
                    for(std::size_t f=0;f<prepared.faceRows.size();++f) {
                        const auto& shared=this->faces[static_cast<std::size_t>(prepared.faceRows[f])];
                        for(int mode=0;mode<3;++mode) {
                            const double coefficient=schur(localRow,3*static_cast<Eigen::Index>(f)+mode);
                            require(std::isfinite(coefficient),"nonfinite condensed current residual entry");
                            if(shared.prescribed)
                                fixed+=static_cast<long double>(coefficient)*shared.coefficients[static_cast<std::size_t>(mode)];
                            else if(coefficient!=0)residualEntries[phase].emplace_back(
                                currentRow,shared.offset+mode,coefficient);
                        }
                    }
                    const double fixedValue=static_cast<double>(fixed);
                    require(std::isfinite(fixedValue),"unrepresentable condensed current source/Dirichlet offset");
                    this->residualOffset[phase][static_cast<std::size_t>(currentRow)]=fixedValue;
                }
            }
        }
        if(prepareAffine) {
            for(std::size_t phase=0;phase<2;++phase) {
                this->residualMap[phase].resize(interfaceRows,this->report.skeleton.freeSkeletonDofs);
                this->residualMap[phase].setFromTriplets(residualEntries[phase].begin(),residualEntries[phase].end());
                this->residualMap[phase].makeCompressed();
            }
        }
        require(this->report.maxRelativeInterfaceBalance<=controls.relativeLinearTolerance
            && this->report.maxRelativeInterfaceBalance<=controls.relativeConservationTolerance,
            "raw interface residual sum minus charge moments exceeds its solve/conservation tolerance");
        // Green's constant test recovers physical D.out from raw rows. External
        // Neumann data are restored before the sign change, never added to Gamma.
        long double exterior=0,bulk=0,globalScale=0;
        std::array<long double,2> exteriorByPhase{},scaleByPhase{};
        this->bulkLoad0.resize(local.size());
        for(std::size_t element=0;element<local.size();++element) {
            const std::size_t phase=geometry.elements()[element].liquid?0:1;
            this->bulkLoad0[element]=local[element].cellLoad[0];
            long double outward=0,equationScale=std::abs(local[element].cellLoad[0]);
            const auto& recovered=solved.elements[element];
            for(std::size_t face=0;face<local[element].faceRows.size();++face) {
                const double displacement=-recovered.physicalPrimal[face][0];
                outward+=displacement;
                equationScale+=recovered.physicalScale[face][0];
                const auto shared=static_cast<std::size_t>(local[element].faceRows[face]);
                if(geometry.faces()[shared].boundarySide>=0) {
                    exterior+=displacement;exteriorByPhase[phase]+=displacement;
                }
            }
            const long double cellResidual=outward-local[element].cellLoad[0];
            this->report.maxRelativeLocalGaussResidual=std::max(this->report.maxRelativeLocalGaussResidual,
                normalized(static_cast<double>(std::abs(cellResidual)),static_cast<double>(equationScale)));
            bulk+=local[element].cellLoad[0];globalScale+=equationScale;
            scaleByPhase[phase]+=equationScale;
        }
        for(std::size_t phase=0;phase<2;++phase) {
            this->exteriorDisplacement[phase]=static_cast<double>(exteriorByPhase[phase]);
            this->currentEquationScale[phase]=static_cast<double>(scaleByPhase[phase]);
            require(std::isfinite(this->exteriorDisplacement[phase])
                && std::isfinite(this->currentEquationScale[phase]),"unrepresentable retained phase current balance data");
        }
        long double surface=0;
        if(charge)for(double q:charge->panelCharge()) {surface+=q;globalScale+=std::abs(q);}
        const double gauss=static_cast<double>(exterior-surface-bulk);
        require(std::isfinite(gauss),"nonfinite global displacement Gauss balance");
        this->report.globalGaussResidual=gauss;
        this->report.relativeGlobalGaussResidual=normalized(std::abs(gauss),static_cast<double>(globalScale));
        require(this->report.maxRelativeLocalGaussResidual<=controls.relativeLinearTolerance
            && this->report.relativeGlobalGaussResidual<=controls.relativeLinearTolerance
            && this->report.maxRelativeLocalGaussResidual<=controls.relativeConservationTolerance
            && this->report.relativeGlobalGaussResidual<=controls.relativeConservationTolerance,
            "raw physical displacement Gauss balance exceeds its solve/conservation tolerance");
        this->report.skeleton.timing.gaussAffineSeconds+=elapsedSeconds(gaussAffineBegin);
    }
    const ElectricPhaseElement& element(std::size_t index) const
    {
        require(index<geometry.elements().size(),"invalid phase element index");
        return geometry.elements()[index];
    }
    const RecoveredElement& recovered(std::size_t index) const
    {
        require(index<elements.size(),"invalid recovered phase element index");
        return elements[index];
    }
    void requirePoint(std::size_t index,Vector2D point) const
    {
        const auto& owner=element(index);
        const auto& source=geometry.publication().state().geometry().snapshot();
        const auto bounds=source.nativeCellBounds(owner.nativeCell);
        require(std::isfinite(point.x) && std::isfinite(point.y)
            && point.x>=bounds[0].x && point.x<=bounds[1].x
            && point.y>=bounds[0].y && point.y<=bounds[1].y,
            "phase polynomial query leaves its actual closed native cell");
        // The caller explicitly selects this phase element. A same-source
        // interface display point can have a residual of either sign; retain
        // both one-sided polynomial traces instead of relocating its material.
    }
};
ElectricPhaseResult::ElectricPhaseResult(std::shared_ptr<const Data> data):data_(std::move(data)){}
const ElectricPhaseResult::Data& ElectricPhaseResult::checked() const
{
    if(!data_)throw std::logic_error("[ElectricPhaseElectrostatics] moved-from result");
    return *data_;
}
const ElectricPhaseGeometry& ElectricPhaseResult::geometry() const{return checked().geometry;}
const ElectricPhaseReport& ElectricPhaseResult::report() const{return checked().report;}
const std::array<double,2>& ElectricPhaseResult::epsilonByPhase() const{return checked().epsilonByPhase;}
Interface::HermiteGalsTimeLayer ElectricPhaseResult::timeLayer() const{return checked().solutionLayer;}
void ElectricPhaseResult::validateSource(const QuadtreeMesh& mesh,
    const Interface::HermiteGalsPublication& publication) const
{checked().geometry.validateSource(mesh,publication);}
bool ElectricPhaseResult::axisymmetric() const{return checked().geometry.axisymmetric();}
std::size_t ElectricPhaseResult::phaseElementCount() const{return checked().elements.size();}
int ElectricPhaseResult::nativeCellIndex(std::size_t index) const{return checked().element(index).nativeCell;}
bool ElectricPhaseResult::liquid(std::size_t index) const{return checked().element(index).liquid;}
std::vector<std::size_t> ElectricPhaseResult::phaseElements(int native,bool negative) const
{
    const auto& data=checked();
    require(native>=0 && static_cast<std::size_t>(native)
        <data.geometry.publication().state().geometry().snapshot().nativeCellCount(),
        "invalid phase native cell index");
    std::vector<std::size_t> result;
    for(std::size_t index=0;index<data.geometry.elements().size();++index) {
        const auto& element=data.geometry.elements()[index];
        if(element.nativeCell==native && element.liquid==negative)result.push_back(index);
    }
    return result;
}
Vector2D ElectricPhaseResult::elementCenter(std::size_t index) const{return checked().element(index).geometry.center;}
double ElectricPhaseResult::elementScale(std::size_t index) const{return checked().element(index).geometry.scale;}
std::array<double,6> ElectricPhaseResult::cellPolynomialP2(std::size_t index) const{return checked().recovered(index).cell;}
std::array<double,10> ElectricPhaseResult::reconstructedPotentialP3(std::size_t index) const{return checked().recovered(index).potential;}
double ElectricPhaseResult::potentialInElement(std::size_t index,Vector2D point) const
{
    const auto& data=checked();data.requirePoint(index,point);
    const auto& element=data.element(index);
    const auto basis=Electric::p3Basis(point,element.geometry.center,element.geometry.scale);
    const auto& coefficients=data.recovered(index).potential;
    long double value=0;
    for(int j=0;j<10;++j)value+=static_cast<long double>(basis[j])*coefficients[static_cast<std::size_t>(j)];
    const double result=static_cast<double>(value);
    require(std::isfinite(result),"nonfinite retained phase potential query");return result;
}
Vector2D ElectricPhaseResult::gradientInElement(std::size_t index,Vector2D point) const
{
    const auto& data=checked();data.requirePoint(index,point);
    const auto& element=data.element(index);
    const auto basis=Electric::p3Gradient(point,element.geometry.center,element.geometry.scale);
    const auto& coefficients=data.recovered(index).potential;
    long double dx=0,dy=0;
    for(int j=0;j<10;++j) {
        dx+=static_cast<long double>(basis(0,j))*coefficients[static_cast<std::size_t>(j)];
        dy+=static_cast<long double>(basis(1,j))*coefficients[static_cast<std::size_t>(j)];
    }
    const Vector2D result{static_cast<double>(dx),static_cast<double>(dy)};
    require(std::isfinite(result.x) && std::isfinite(result.y),"nonfinite retained phase gradient query");return result;
}
std::array<double,3> ElectricPhaseResult::skeletalFaceTrace(std::size_t face) const
{
    const auto& data=checked();
    require(face<data.faces.size(),"invalid phase skeletal face index");
    return data.faces[face].coefficients;
}
std::array<double,3> ElectricPhaseResult::outwardFaceResidual(std::size_t element,std::size_t face) const
{
    const auto& residual=checked().recovered(element).residual;
    require(face<residual.size(),"invalid local phase element face index");
    return residual[face];
}
std::array<double,3> ElectricPhaseResult::physicalOutwardDisplacementMoments(
    std::size_t element,std::size_t face) const
{
    const auto& rows=checked().recovered(element).physicalPrimal;
    require(face<rows.size(),"invalid physical displacement local face index");
    auto result=rows[face];for(double& value:result)value=-value;return result;
}
std::array<double,3> ElectricPhaseResult::normalDisplacementMoments(std::size_t face,bool negative) const
{
    const auto& data=checked();
    require(face<data.displacement.size() && data.displacement[face].available,
        "normal displacement query requires an actual interface trace");
    return data.displacement[face].moments[negative?0:1];
}
std::array<double,3> ElectricPhaseResult::normalDisplacementPolynomialP2(std::size_t face,bool negative) const
{
    const auto& data=checked();
    require(face<data.displacement.size() && data.displacement[face].available,
        "normal displacement polynomial query requires an actual interface trace");
    return data.displacement[face].polynomial[negative?0:1];
}
std::array<double,3> ElectricPhaseResult::interfaceChargeMoments(std::size_t face) const
{
    const auto& data=checked();
    require(face<data.geometry.faces().size()
        && data.geometry.faces()[face].kind==ElectricPhaseSkeletonFace::Kind::Interface,
        "surface charge query requires an actual interface trace");
    return data.chargeLoad?data.chargeLoad->faceMoments()[face]:std::array<double,3>{{0,0,0}};
}
const ElectricPanelChargeLoad* ElectricPhaseResult::surfaceChargeLoad() const{return checked().chargeLoad.get();}
const ElectricPhaseResult::Matrix& ElectricPhaseResult::skeletonMatrix() const{return *checked().skeleton;}
const std::vector<double>& ElectricPhaseResult::skeletonRhs() const{return checked().rhs;}
const std::vector<double>& ElectricPhaseResult::skeletonSolution() const{return checked().solution;}

struct ElectricPanelCurrentOperator::Data
{
    Data(ElectricPhaseResult field,ElectricSurfaceChargeMomentOperator binding)
        :source(std::move(field)),momentOperator(std::move(binding)){}
    ElectricPhaseResult source;
    ElectricSurfaceChargeMomentOperator momentOperator;
    std::array<double,2> epsilon{},conductivity{},ratio{};
    ElectricPanelCurrentReport report;
    std::shared_ptr<const Matrix> matrix; // Fixed affine functional shared by all coeval recoveries.
    std::vector<double> offset,current;
};
ElectricPanelCurrentOperator::ElectricPanelCurrentOperator(std::shared_ptr<const Data> data):data_(std::move(data)){}
const ElectricPanelCurrentOperator::Data& ElectricPanelCurrentOperator::checked() const
{
    if(!data_)throw std::logic_error("[ElectricPanelCurrentOperator] moved-from current operator");
    return *data_;
}
const ElectricSurfaceChargeMomentOperator& ElectricPanelCurrentOperator::momentOperator() const
{return checked().momentOperator;}
const ElectricPhaseGeometry& ElectricPanelCurrentOperator::geometry() const{return checked().momentOperator.geometry();}
Interface::HermiteGalsTimeLayer ElectricPanelCurrentOperator::timeLayer() const{return checked().source.timeLayer();}
const ElectricPanelCurrentReport& ElectricPanelCurrentOperator::report() const{return checked().report;}
const std::array<double,2>& ElectricPanelCurrentOperator::epsilonByPhase() const{return checked().epsilon;}
const std::array<double,2>& ElectricPanelCurrentOperator::conductivityByPhase() const{return checked().conductivity;}
const ElectricPanelCurrentOperator::Matrix& ElectricPanelCurrentOperator::matrix() const{return *checked().matrix;}
const std::vector<double>& ElectricPanelCurrentOperator::offset() const{return checked().offset;}
const std::vector<double>& ElectricPanelCurrentOperator::panelCurrent() const{return checked().current;}
std::vector<double> ElectricPanelCurrentOperator::apply(const std::vector<double>& freeSkeleton) const
{
    const auto& data=checked();
    require(freeSkeleton.size()==static_cast<std::size_t>(data.matrix->cols()),
        "current functional has foreign free skeleton ordering/size");
    for(double value:freeSkeleton)require(std::isfinite(value),"nonfinite trial free skeleton potential");
    const Eigen::Map<const Eigen::VectorXd> u(freeSkeleton.data(),static_cast<Eigen::Index>(freeSkeleton.size()));
    const Eigen::Map<const Eigen::VectorXd> fixed(data.offset.data(),static_cast<Eigen::Index>(data.offset.size()));
    const Eigen::VectorXd result=(*data.matrix)*u+fixed;
    require(result.allFinite(),"unrepresentable affine panel current");
    return std::vector<double>(result.data(),result.data()+result.size());
}
std::array<double,3> ElectricPanelCurrentOperator::normalCurrentMoments(std::size_t panel,bool negative) const
{
    const auto& data=checked();
    auto moments=data.source.normalDisplacementMoments(
        static_cast<std::size_t>(data.momentOperator.interfaceFace(panel)),negative);
    const double ratio=data.ratio[negative?0:1];
    for(double& value:moments) {
        value*=ratio;require(std::isfinite(value),"unrepresentable one-sided normal current moment");
    }
    return moments;
}
std::array<double,3> ElectricPanelCurrentOperator::normalCurrentPolynomialP2(std::size_t panel,bool negative) const
{
    const auto& data=checked();
    auto polynomial=data.source.normalDisplacementPolynomialP2(
        static_cast<std::size_t>(data.momentOperator.interfaceFace(panel)),negative);
    const double ratio=data.ratio[negative?0:1];
    for(double& value:polynomial) {
        value*=ratio;require(std::isfinite(value),"unrepresentable one-sided normal current polynomial");
    }
    return polynomial;
}
ElectricPanelCurrentOperator ElectricPhaseResult::preparePanelCurrent(
    const ElectricSurfaceChargeMomentOperator& momentOperator,const std::array<double,2>& conductivity) const
{
    return preparePanelCurrentBound(momentOperator,conductivity,nullptr);
}
ElectricPanelCurrentOperator ElectricPhaseResult::preparePanelCurrentBound(
    const ElectricSurfaceChargeMomentOperator& momentOperator,const std::array<double,2>& conductivity,
    const ElectricPanelCurrentOperator* retained) const
{
    const auto& field=checked();
    const auto& bound=momentOperator.geometry();
    require(bound.geometryIdentity()==field.geometry.geometryIdentity()
        && sameClock(bound.timeLayer(),field.geometry.timeLayer())
        && &bound.publication().state().geometry().snapshot()
            ==&field.geometry.publication().state().geometry().snapshot()
        && &bound.faces()==&field.geometry.faces() && &bound.elements()==&field.geometry.elements(),
        "current functional requires the exact electric solve geometry/basis/source/clock owner");
    auto data=retained?std::make_shared<ElectricPanelCurrentOperator::Data>(retained->checked())
        :std::make_shared<ElectricPanelCurrentOperator::Data>(*this,momentOperator);
    if(retained) {
        require(&retained->momentOperator().geometry().faces()==&bound.faces()
            && retained->conductivityByPhase()==conductivity
            && retained->epsilonByPhase()==field.epsilonByPhase,
            "retained current operator differs from the fixed geometry/material preparation");
        data->source=*this;data->report={};
    }
    data->epsilon=field.epsilonByPhase;data->conductivity=conductivity;
    for(std::size_t phase=0;phase<2;++phase) {
        const double sigma=conductivity[phase],epsilon=data->epsilon[phase];
        require(std::isfinite(sigma) && sigma>=0 && std::isfinite(epsilon) && epsilon>0,
            "conductivity must be finite nonnegative with the actual positive solve permittivity");
        data->ratio[phase]=sigma/epsilon;
        require(std::isfinite(data->ratio[phase]) && (sigma==0 || data->ratio[phase]>0),
            "positive conductivity/permittivity ratio overflows or underflows its representation");
    }
    const auto count=momentOperator.panelCount();
    require(count==static_cast<std::size_t>(field.residualMap[0].rows())
        && field.residualMap[0].rows()==field.residualMap[1].rows(),
        "current functional has foreign material panel coverage");
    std::vector<int> rowToPanel(count,-1);
    for(std::size_t panel=0;panel<count;++panel) {
        const int face=momentOperator.interfaceFace(panel);
        require(face>=0 && static_cast<std::size_t>(face)<field.interfaceRow.size(),
            "current material panel has foreign interface trace");
        const int row=field.interfaceRow[static_cast<std::size_t>(face)];
        require(row>=0 && static_cast<std::size_t>(row)<count
            && rowToPanel[static_cast<std::size_t>(row)]<0,"current panel repeats or loses a true Gamma row");
        rowToPanel[static_cast<std::size_t>(row)]=static_cast<int>(panel);
    }
    if(!retained) {
    auto matrix=std::make_shared<Matrix>();
    std::vector<Eigen::Triplet<double>> entries;
    entries.reserve(static_cast<std::size_t>(field.residualMap[0].nonZeros()+field.residualMap[1].nonZeros()));
    for(std::size_t phase=0;phase<2;++phase) {
        const auto& map=field.residualMap[phase];
        require(map.cols()==static_cast<Eigen::Index>(field.solution.size()),
            "current residual map differs from the solved free skeleton ordering");
        for(int col=0;col<map.outerSize();++col)for(Matrix::InnerIterator entry(map,col);entry;++entry) {
            const double coefficient=-data->ratio[phase]*entry.value();
            require(std::isfinite(coefficient),"unrepresentable Ohmic current matrix coefficient");
            if(coefficient!=0)entries.emplace_back(
                rowToPanel[static_cast<std::size_t>(entry.row())],static_cast<int>(entry.col()),coefficient);
        }
    }
    matrix->resize(static_cast<int>(count),static_cast<int>(field.solution.size()));
    matrix->setFromTriplets(entries.begin(),entries.end());matrix->makeCompressed();
    for(int col=0;col<matrix->outerSize();++col)
        for(Matrix::InnerIterator entry(*matrix,col);entry;++entry)
            require(std::isfinite(entry.value()),"nonfinite summed Ohmic current matrix coefficient");
    data->matrix=std::move(matrix);
    }
    data->offset.resize(count);data->current.resize(count);
    long double total=0;
    for(std::size_t row=0;row<count;++row) {
        const auto panel=static_cast<std::size_t>(rowToPanel[row]);
        const auto face=static_cast<std::size_t>(momentOperator.interfaceFace(panel));
        const long double fixed=-static_cast<long double>(data->ratio[0])*field.residualOffset[0][row]
            -static_cast<long double>(data->ratio[1])*field.residualOffset[1][row];
        // Dminus=-rminus, Dplus=rplus in the one canonical interface normal.
        const long double current=static_cast<long double>(data->ratio[0])*field.displacement[face].moments[0][0]
            -static_cast<long double>(data->ratio[1])*field.displacement[face].moments[1][0];
        if(!retained)data->offset[panel]=static_cast<double>(fixed);
        data->current[panel]=static_cast<double>(current);
        require(std::isfinite(data->offset[panel]) && std::isfinite(data->current[panel]),
            "unrepresentable original residual panel current or affine offset");
        total+=data->current[panel];
    }
    ElectricPanelCurrentOperator result(data);
    const auto affine=result.apply(field.solution);
    std::vector<long double> rawScale(count,0);
    for(std::size_t element=0;element<field.elements.size();++element) {
        const auto phase=field.geometry.elements()[element].liquid?0:1;
        const auto& faces=field.geometry.elements()[element].faces;
        for(std::size_t localFace=0;localFace<faces.size();++localFace) {
            const int row=field.interfaceRow[static_cast<std::size_t>(faces[localFace])];
            if(row>=0)rawScale[static_cast<std::size_t>(rowToPanel[static_cast<std::size_t>(row)])]
                +=static_cast<long double>(data->ratio[phase])*field.elements[element].physicalScale[localFace][0];
        }
    }
    const Eigen::Map<const Eigen::VectorXd> u(field.solution.data(),static_cast<Eigen::Index>(field.solution.size()));
    const Eigen::VectorXd matrixScale=data->matrix->cwiseAbs()*u.cwiseAbs();
    for(std::size_t panel=0;panel<count;++panel) {
        const long double scale=rawScale[panel]+matrixScale[static_cast<Eigen::Index>(panel)]
            +std::abs(data->offset[panel]);
        data->report.maxRelativeAffineResidual=std::max(data->report.maxRelativeAffineResidual,
            normalized(std::abs(affine[panel]-data->current[panel]),static_cast<double>(scale)));
    }
    long double exterior=0,bulk=0,scale=0;
    for(std::size_t phase=0;phase<2;++phase) {
        exterior+=static_cast<long double>(data->ratio[phase])*field.exteriorDisplacement[phase];
        scale+=static_cast<long double>(data->ratio[phase])*field.currentEquationScale[phase];
    }
    for(std::size_t element=0;element<field.bulkLoad0.size();++element)
        bulk+=static_cast<long double>(data->ratio[field.geometry.elements()[element].liquid?0:1])*field.bulkLoad0[element];
    const double balance=static_cast<double>(total+exterior-bulk);
    data->report.panels=count;data->report.totalPanelCurrent=static_cast<double>(total);
    data->report.exteriorCurrent=static_cast<double>(exterior);data->report.bulkCurrent=static_cast<double>(bulk);
    data->report.globalCurrentBalance=balance;
    require(std::isfinite(data->report.totalPanelCurrent) && std::isfinite(data->report.exteriorCurrent)
        && std::isfinite(data->report.bulkCurrent) && std::isfinite(balance),
        "unrepresentable total physical current balance");
    data->report.relativeGlobalCurrentBalance=normalized(std::abs(balance),static_cast<double>(scale));
    require(data->report.maxRelativeAffineResidual<=field.relativeLinearTolerance
        && data->report.maxRelativeAffineResidual<=field.relativeConservationTolerance
        && data->report.relativeGlobalCurrentBalance<=field.relativeLinearTolerance
        && data->report.relativeGlobalCurrentBalance<=field.relativeConservationTolerance,
        "original Ohmic current affine or physical balance exceeds the solve/conservation tolerance");
    return result;
}

ElectricPhaseResult ElectricPhaseResult::solvePrepared(const ElectricPhaseGeometry& geometry,
    const ElectricPhaseProblem& problem,const ElectricPanelChargeLoad* charge,
    const ElectricCartesianControls& controls)
{
    const auto boundaryBegin=TimingClock::now();
    validateControls(controls);validateBoundaries(problem.boundaries);
    if(charge)charge->validateGeometry(geometry);
    for(double epsilon:problem.epsilonByPhase)
        require(std::isfinite(epsilon) && epsilon>0,"phase permittivity must be finite and positive");
    const auto& source=geometry.publication().state().geometry().snapshot();
    require(!geometry.elements().empty() && !geometry.faces().empty(),"phase geometry is empty");
    auto data=std::make_shared<ElectricPhaseResult::Data>(geometry);
    if(charge)data->chargeLoad=std::make_shared<const ElectricPanelChargeLoad>(*charge);
    data->epsilonByPhase=problem.epsilonByPhase;
    data->relativeLinearTolerance=controls.relativeLinearTolerance;
    data->relativeConservationTolerance=controls.relativeConservationTolerance;
    auto& report=data->report.skeleton;
    report.cells=geometry.elements().size();
    report.physicalFaces=geometry.faces().size();
    data->report.interfaceFaces=geometry.report().interfaceFaces;
    // The adapter omits r=0 faces entirely. The original physical boundary
    // law still has to be explicit homogeneous Neumann, exactly as in P0.
    if(geometry.axisymmetric() && source.minCoord().x==0) {
        require(problem.boundaries[0].type==ElectricCartesianBoundary::Type::Neumann,
            "RZ axis must be explicit homogeneous Neumann");
        const auto lower=source.minCoord(),upper=source.maxCoord();
        for(double q:gaussPoints) {
            const double value=problem.boundaries[0].value({0,coordinate(lower.y,upper.y,q)});
            require(std::isfinite(value) && value==0,
                "RZ axis Neumann provider is not identically zero at its rule");
        }
        for(std::size_t cell=0;cell<source.nativeCellCount();++cell)
            if(source.nativeCellBounds(static_cast<int>(cell))[0].x==0)++report.axisFaces;
    }
    data->faces.reserve(geometry.faces().size());
    for(const auto& retained:geometry.faces()) {
        FaceData face;face.boundary=retained.boundarySide;face.nodes=retained.nodes;
        require(face.boundary>=-1 && face.boundary<4,"invalid retained phase boundary side");
        require(!face.nodes.empty(),"empty retained phase subface");
        if(retained.kind==ElectricPhaseSkeletonFace::Kind::Interface)
            require(face.boundary==-1,"an interface trace cannot carry an exterior boundary law");
        prepareFaceUnknown(face,problem.boundaries,geometry.axisymmetric(),report);
        // The authoritative geometry already retains these quadrature nodes.
        // After boundary projection FaceData needs only coefficients/ordering.
        std::vector<Electric::FaceQuadraturePoint>().swap(face.nodes);
        data->faces.push_back(std::move(face));
    }
    require(std::any_of(data->faces.begin(),data->faces.end(),
        [](const FaceData& face){return face.prescribed;}),
        "a positive-measure Dirichlet face must anchor the phase skeleton");
    report.timing.boundarySeconds=elapsedSeconds(boundaryBegin);
    std::vector<PreparedCell> local;local.reserve(geometry.elements().size());
    std::vector<std::vector<std::size_t>> incidence(geometry.faces().size());
    for(std::size_t index=0;index<geometry.elements().size();++index) {
        const auto localBegin=TimingClock::now();
        const auto& element=geometry.elements()[index];
        require(element.nativeCell>=0 && static_cast<std::size_t>(element.nativeCell)<source.nativeCellCount()
            && element.geometry.axisymmetric==geometry.axisymmetric(),
            "phase element has foreign native identity or physical measure");
        require(element.faces.size()==element.geometry.faces.size(),"phase local face layout differs from its shared connectivity");
        for(std::size_t face=0;face<element.faces.size();++face) {
            const int shared=element.faces[face];
            require(shared>=0 && static_cast<std::size_t>(shared)<geometry.faces().size(),
                "phase element has foreign shared face identity");
            require(std::find(element.faces.begin(),element.faces.begin()+static_cast<std::ptrdiff_t>(face),shared)
                    ==element.faces.begin()+static_cast<std::ptrdiff_t>(face),
                "phase element repeats one physical shared face");
            incidence[static_cast<std::size_t>(shared)].push_back(index);
        }
        Electric::LocalOperator op(element.geometry,problem.epsilonByPhase[element.liquid?0:1]);
        report.timing.localOperatorsSeconds+=elapsedSeconds(localBegin);
        const auto loadBegin=TimingClock::now();
        const Electric::ScalarFunction& provider=problem.sourceByPhase[element.liquid?0:1];
        const Electric::CellVector cellLoad=provider?op.loadCell(provider):Electric::CellVector::Zero().eval();
        Eigen::VectorXd faceLoad=Eigen::VectorXd::Zero(3*static_cast<int>(element.faces.size()));
        for(std::size_t face=0;face<element.faces.size();++face) {
            const auto& retained=data->faces[static_cast<std::size_t>(element.faces[face])];
            if(retained.boundary>=0 && !retained.prescribed)
                faceLoad.segment<3>(3*static_cast<Eigen::Index>(face))=op.loadFace(face,
                    problem.boundaries[static_cast<std::size_t>(retained.boundary)].value);
        }
        local.emplace_back(std::move(op),cellLoad,std::move(faceLoad),element.faces);
        report.timing.localLoadsSeconds+=elapsedSeconds(loadBegin);
    }
    for(std::size_t face=0;face<geometry.faces().size();++face) {
        const auto& retained=geometry.faces()[face];const auto& owners=incidence[face];
        require(owners.size()==static_cast<std::size_t>(retained.boundarySide>=0?1:2),
            "physical phase subface does not have the required one/two incident elements");
        if(retained.kind==ElectricPhaseSkeletonFace::Kind::Interface) {
            require(owners.size()==2 && geometry.elements()[owners[0]].liquid!=geometry.elements()[owners[1]].liquid,
                "interface trace must have one true element from each phase");
        } else {
            for(std::size_t owner:owners)require(geometry.elements()[owner].liquid==retained.liquid,
                "Cartesian phase subface belongs to a different incident phase");
        }
    }
    const std::vector<std::array<double,3>> noCharge;
    const auto& sharedLoads=charge?charge->faceMoments():noCharge;
    auto solved=assembleCondensedSkeleton(data->faces,local,controls,report,sharedLoads);
    data->finishRecovery(local,incidence,solved,charge,controls,true);
    auto backend=std::make_shared<PreparedPhaseBackend>();
    backend->local=std::move(local);backend->incidence=std::move(incidence);
    backend->fixedRhs=std::move(solved.fixedRhs);backend->factor=std::move(solved.factor);
    backend->controls=controls;backend->matrixNorm=solved.matrixNorm;data->prepared=std::move(backend);
    data->elements=std::move(solved.elements);
    data->skeleton=std::make_shared<const Matrix>(std::move(solved.matrix));
    data->rhs=std::move(solved.rhs);data->solution=std::move(solved.solution);
    return ElectricPhaseResult(std::move(data));
}

ElectricPhaseResult ElectricPhaseResult::recoverPrepared(const ElectricPhaseResult& base,
    const ElectricPanelChargeLoad& charge,const std::vector<double>& freeSkeleton,
    const Interface::HermiteGalsTimeLayer& solutionLayer)
{
    const auto& retained=base.checked();
    charge.validateGeometry(retained.geometry);(void)solutionLayer.time();
    require(static_cast<bool>(retained.prepared),"fixed-geometry recovery has no original prepared backend");
    const auto& backend=*retained.prepared;
    const auto size=static_cast<Eigen::Index>(retained.solution.size());
    require(freeSkeleton.size()==retained.solution.size(),"fixed-geometry recovery has foreign skeleton layout");
    const Eigen::Map<const Eigen::VectorXd> u(freeSkeleton.data(),size);
    require(u.allFinite(),"nonfinite final fixed-geometry electric solution");
    // Copy only retained preparation; old recovered fields, RHS and solution
    // are replaced below and need not be allocated/copied for each candidate.
    auto data=std::make_shared<Data>(retained.geometry);
    data->prepared=retained.prepared;data->report=retained.report;data->faces=retained.faces;
    data->displacement=retained.displacement;data->epsilonByPhase=retained.epsilonByPhase;
    data->interfaceRow=retained.interfaceRow;data->residualMap=retained.residualMap;
    data->residualOffset=retained.residualOffset;
    data->relativeLinearTolerance=retained.relativeLinearTolerance;
    data->relativeConservationTolerance=retained.relativeConservationTolerance;data->skeleton=retained.skeleton;
    data->report.skeleton.timing={}; // Do not inherit initial preparation times.
    data->solutionLayer=solutionLayer;
    data->chargeLoad=std::make_shared<const ElectricPanelChargeLoad>(charge);
    Eigen::VectorXd rhs=Eigen::Map<const Eigen::VectorXd>(backend.fixedRhs.data(),size);
    for(std::size_t f=0;f<data->faces.size();++f) {
        const auto& face=data->faces[f];
        if(!face.prescribed && !face.axis)
            for(int mode=0;mode<3;++mode)rhs[face.offset+mode]+=charge.faceMoments()[f][static_cast<std::size_t>(mode)];
    }
    require(rhs.allFinite(),"nonfinite final charge-loaded skeleton RHS");
    auto& report=data->report.skeleton;
    report.relativeLinearResidual=normalized(((*data->skeleton)*u-rhs).norm(),
        rhs.norm()+backend.matrixNorm*u.norm());
    require(report.relativeLinearResidual<=backend.controls.relativeLinearTolerance,
        "final coeval skeleton residual exceeds its numerical solve tolerance");
    SkeletonSolution solved;
    recoverCondensedSkeleton(data->faces,backend.local,report,u,solved,charge.faceMoments());
    data->finishRecovery(backend.local,backend.incidence,solved,&charge,backend.controls,false);
    data->elements=std::move(solved.elements);data->solution=freeSkeleton;
    data->rhs.assign(rhs.data(),rhs.data()+rhs.size());
    return ElectricPhaseResult(std::move(data));
}

ElectricPhaseResult solveElectricPhaseElectrostatics(const ElectricPhaseGeometry& geometry,
    const ElectricPhaseProblem& problem,const ElectricCartesianControls& controls)
{
    return ElectricPhaseResult::solvePrepared(geometry,problem,nullptr,controls);
}
ElectricPhaseResult solveElectricPhaseElectrostatics(const ElectricPhaseGeometry& geometry,
    const ElectricPhaseProblem& problem,const ElectricPanelChargeLoad& charge,const ElectricCartesianControls& controls)
{
    return ElectricPhaseResult::solvePrepared(geometry,problem,&charge,controls);
}

namespace
{
struct AcceptedChargeLevel
{
    Interface::HermiteGalsTimeLayer layer,representationLayer;
    Interface::HermiteGalsIdentity representationIdentity;
    std::vector<double> charge;
};
bool sameChargeBits(const std::vector<double>& a,const std::vector<double>& b)
{
    return a.size()==b.size() && (a.empty()
        || std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0);
}
void validateChargeTimeControls(const ElectricChargeTimeControls& controls)
{
    require((controls.method==ElectricChargeTimeMethod::BackwardEuler
            || controls.method==ElectricChargeTimeMethod::VariableBdf2)
        && std::isfinite(controls.relativeSchurTolerance)
        && controls.relativeSchurTolerance>0 && controls.relativeSchurTolerance<1
        && std::isfinite(controls.relativeChargeBalanceTolerance)
        && controls.relativeChargeBalanceTolerance>0 && controls.relativeChargeBalanceTolerance<=1e-10
        && controls.maximumIterations>0 && controls.restartLength>0
        && std::isfinite(controls.maximumStepRatio) && controls.maximumStepRatio>0
        && controls.maximumStepRatio<=2.0,
        "charge time controls: invalid method, numerical/physical tolerance, iteration or stable step-ratio bound");
}
double validateChargeTarget(const Interface::HermiteGalsTimeLayer& accepted,
    const Interface::HermiteGalsTimeLayer& target,bool midpoint=false)
{
    using Position=Interface::HermiteGalsTimeLayer::Position;
    const double time=accepted.time(),targetTime=target.time();
    const auto position=midpoint?Position::Midpoint:Position::End;
    require((accepted.position==Position::Initial || accepted.position==Position::End)
        && target.position==position
        && accepted.step<std::numeric_limits<std::int64_t>::max()
        && target.step==accepted.step+1 && target.begin==time,
        midpoint
            ?"charge stage preparation: target must be the next Midpoint layer beginning at the actual accepted time"
            :"charge time advance: target must be the next End layer beginning at the actual accepted time");
    // The publication's full clock still spans the Hydro interval. Only the
    // physical increment to its actual stage time enters the charge equation.
    const double dt=midpoint?targetTime-time:target.end-target.begin;
    require(std::isfinite(dt) && dt>0,"charge time advance: nonpositive/unrepresentable time increment");
    return dt;
}
template<class Action>
Eigen::VectorXd solveChargeGmres(const Eigen::VectorXd& rhs,const Eigen::VectorXd& initial,
    const Action& action,const ElectricChargeTimeControls& controls,ElectricChargeTimeReport& report)
{
    require(rhs.size()==initial.size() && rhs.allFinite() && initial.allFinite(),
        "charge Schur solve: nonfinite RHS/initial guess or foreign panel layout");
    Eigen::VectorXd x=initial;
    Eigen::VectorXd ax=action(x),residual=rhs-ax;
    auto checkResidual=[&]() {
        require(ax.allFinite() && residual.allFinite(),"charge Schur solve: nonfinite true residual");
        report.relativeSchurResidual=normalized(residual.norm(),rhs.norm()+ax.norm());
        return report.relativeSchurResidual<=controls.relativeSchurTolerance;
    };
    if(checkResidual())return x;
    const int count=static_cast<int>(rhs.size());
    require(count>0,"charge Schur solve: empty nonzero charge equation");
    while(report.iterations<controls.maximumIterations) {
        const int window=std::min({controls.restartLength,count,controls.maximumIterations-report.iterations});
        Eigen::MatrixXd basis=Eigen::MatrixXd::Zero(count,window+1);
        Eigen::MatrixXd hessenberg=Eigen::MatrixXd::Zero(window+1,window);
        Eigen::VectorXd cosine=Eigen::VectorXd::Zero(window),sine=cosine;
        Eigen::VectorXd g=Eigen::VectorXd::Zero(window+1);
        const double beta=residual.norm();
        require(std::isfinite(beta) && beta>0,"charge Schur solve: unresolved nonzero restart residual");
        basis.col(0)=residual/beta;g[0]=beta;
        int used=0;bool breakdown=false;
        for(int column=0;column<window;++column) {
            Eigen::VectorXd w=action(basis.col(column));
            require(w.size()==count && w.allFinite(),"charge Schur solve: nonfinite matrix-free action");
            // Two-pass modified Gram-Schmidt controls loss of orthogonality;
            // only the small Krylov Hessenberg is dense, never the charge map M.
            for(int pass=0;pass<2;++pass)for(int row=0;row<=column;++row) {
                const double projection=basis.col(row).dot(w);
                hessenberg(row,column)+=projection;w.noalias()-=projection*basis.col(row);
            }
            const double nextNorm=w.norm();
            require(std::isfinite(nextNorm),"charge Schur solve: nonfinite Arnoldi norm");
            hessenberg(column+1,column)=nextNorm;
            if(nextNorm>0)basis.col(column+1)=w/nextNorm;
            else breakdown=true;
            for(int row=0;row<column;++row) {
                const double upper=hessenberg(row,column),lower=hessenberg(row+1,column);
                hessenberg(row,column)=cosine[row]*upper+sine[row]*lower;
                hessenberg(row+1,column)=-sine[row]*upper+cosine[row]*lower;
            }
            const double diagonal=hessenberg(column,column),lower=hessenberg(column+1,column);
            const double rotation=std::hypot(diagonal,lower);
            require(std::isfinite(rotation) && rotation>0,
                "charge Schur solve: singular Arnoldi column (no fallback or regularization)");
            cosine[column]=diagonal/rotation;sine[column]=lower/rotation;
            hessenberg(column,column)=rotation;hessenberg(column+1,column)=0;
            g[column+1]=-sine[column]*g[column];g[column]=cosine[column]*g[column];
            ++report.iterations;used=column+1;
            if(breakdown || std::abs(g[column+1])<=controls.relativeSchurTolerance*rhs.norm())break;
        }
        const Eigen::VectorXd y=hessenberg.topLeftCorner(used,used)
            .template triangularView<Eigen::Upper>().solve(g.head(used));
        require(y.allFinite(),"charge Schur solve: nonfinite Krylov triangular solution");
        x.noalias()+=basis.leftCols(used)*y;
        require(x.allFinite(),"charge Schur solve: nonfinite candidate Q");
        ax=action(x);residual=rhs-ax;
        if(checkResidual())return x;
        require(!breakdown,"charge Schur solve: Arnoldi breakdown failed the true residual gate");
    }
    throw std::invalid_argument("[ElectricFixedGeometryChargeSolver] charge Schur solve exhausted maximumIterations");
}
}

struct ElectricFixedGeometryChargeSolver::Data
{
    Data(ElectricPhaseResult source,ElectricSurfaceChargeMomentOperator binding,ElectricPanelCurrentOperator ohmic)
        :base(std::move(source)),moments(std::move(binding)),current(std::move(ohmic)){}
    ElectricPhaseResult base;
    ElectricSurfaceChargeMomentOperator moments;
    ElectricPanelCurrentOperator current;
    std::shared_ptr<const PreparedPhaseBackend> backend;
    ElectricPhaseResult::Matrix injection; // Actual free B_Gamma P_Q, sparse.
    Eigen::VectorXd fixedCurrent; // I0, including the d-driven field.
    double matrixNorm=0;
    Eigen::VectorXd electricSolve(const Eigen::VectorXd& rhs,int* solves=nullptr,double* seconds=nullptr) const
    {
        const auto solveBegin=TimingClock::now();
        require(rhs.size()==base.skeletonMatrix().rows() && rhs.allFinite(),
            "charge electric solve: nonfinite RHS or foreign free trace layout");
        if(rhs.size()==0)return Eigen::VectorXd{};
        require(static_cast<bool>(backend) && static_cast<bool>(backend->factor),
            "charge electric solve: original SPD factor is unavailable");
        if(solves)++*solves;
        Eigen::VectorXd result=backend->factor->solve(rhs);
        require(result.allFinite(),"charge electric solve: nonfinite original SPD solve result");
        const double relative=normalized((base.skeletonMatrix()*result-rhs).norm(),
            rhs.norm()+matrixNorm*result.norm());
        require(relative<=backend->controls.relativeLinearTolerance,
            "charge electric solve: true inner residual exceeds its numerical solve tolerance");
        if(seconds)*seconds+=elapsedSeconds(solveBegin);
        return result;
    }
    Eigen::VectorXd linearCurrent(const Eigen::VectorXd& charge,int* solves=nullptr,double* seconds=nullptr) const
    {
        require(charge.size()==injection.cols() && charge.allFinite(),
            "charge Schur action: foreign panel layout or nonfinite trial Q");
        Eigen::VectorXd rhs=injection*charge;
        const auto u=electricSolve(rhs,solves,seconds);
        Eigen::VectorXd result=current.matrix()*u;
        require(result.allFinite(),"charge Schur action: nonfinite same-field linear current");
        return result;
    }
};
ElectricFixedGeometryChargeSolver::ElectricFixedGeometryChargeSolver(std::shared_ptr<const Data> data)
    :data_(std::move(data)){}
const ElectricFixedGeometryChargeSolver::Data& ElectricFixedGeometryChargeSolver::checked() const
{
    if(!data_)throw std::logic_error("[ElectricFixedGeometryChargeSolver] moved-from solver");
    return *data_;
}
Interface::HermiteGalsTimeLayer ElectricFixedGeometryChargeSolver::geometrySourceTimeLayer() const
{return checked().base.geometry().timeLayer();}
const ElectricSurfaceChargeMomentOperator& ElectricFixedGeometryChargeSolver::momentOperator() const
{return checked().moments;}
const ElectricPanelCurrentOperator& ElectricFixedGeometryChargeSolver::initialCurrent() const
{return checked().current;}
ElectricFixedGeometryChargeSolver ElectricFixedGeometryChargeSolver::prepare(const ElectricPhaseResult& base,
    const ElectricSurfaceChargeMomentOperator& moments,const std::array<double,2>& conductivity)
{
    const auto& field=base.checked();
    require(static_cast<bool>(field.prepared),"charge preparation: base has no original prepared electric backend");
    auto current=base.preparePanelCurrent(moments,conductivity);
    auto data=std::make_shared<Data>(base,moments,std::move(current));
    data->backend=field.prepared;data->matrixNorm=field.prepared->matrixNorm;
    const auto& full=moments.interfaceInjectionMatrix();
    require(full.rows()==static_cast<Eigen::Index>(3*field.faces.size()),
        "charge preparation: B_Gamma has foreign full skeleton layout");
    std::vector<Eigen::Triplet<double>> entries;
    for(int column=0;column<full.outerSize();++column)
        for(ElectricPhaseResult::Matrix::InnerIterator entry(full,column);entry;++entry) {
            const auto f=static_cast<std::size_t>(entry.row()/3);
            const int mode=static_cast<int>(entry.row()%3);
            const auto& face=field.faces[f];
            require(!face.prescribed && !face.axis && face.offset>=0,
                "charge preparation: B_Gamma touches a prescribed/axis trace");
            entries.emplace_back(face.offset+mode,static_cast<int>(entry.col()),entry.value());
        }
    ElectricPhaseResult::Matrix freeInjection(field.skeleton->rows(),full.cols());
    freeInjection.setFromTriplets(entries.begin(),entries.end());freeInjection.makeCompressed();
    data->injection=freeInjection*moments.projectionMatrix();data->injection.makeCompressed();
    for(int column=0;column<data->injection.outerSize();++column)
        for(ElectricPhaseResult::Matrix::InnerIterator entry(data->injection,column);entry;++entry)
            require(std::isfinite(entry.value()),"charge preparation: nonfinite sparse charge injection");
    const auto size=static_cast<Eigen::Index>(data->backend->fixedRhs.size());
    const Eigen::Map<const Eigen::VectorXd> d(data->backend->fixedRhs.data(),size);
    const auto baseline=data->electricSolve(d);
    data->fixedCurrent=data->current.matrix()*baseline
        +Eigen::Map<const Eigen::VectorXd>(data->current.offset().data(),
            static_cast<Eigen::Index>(data->current.offset().size()));
    require(data->fixedCurrent.allFinite(),"charge preparation: nonfinite physical affine I0");
    return ElectricFixedGeometryChargeSolver(std::move(data));
}

struct ElectricSurfaceChargeTimeState::Data
{
    Data(ElectricFixedGeometryChargeSolver preparation,std::vector<SurfaceChargeField*> owners)
        :solver(std::move(preparation)),fields(std::move(owners)){}
    std::shared_ptr<const std::vector<std::shared_ptr<SurfaceChargeField>>> ownedFields;
    ElectricFixedGeometryChargeSolver solver;
    std::vector<SurfaceChargeField*> fields;
    std::vector<std::size_t> offsets;
    std::shared_ptr<const void> ownerToken=std::make_shared<const int>(0);
    std::shared_ptr<const AcceptedChargeLevel> accepted,previous;
    std::vector<SurfaceChargeTransport> prepareFieldCharge(const std::vector<double>& charge) const
    {
        require(accepted && charge.size()==accepted->charge.size(),
            "charge Field preparation: foreign integrated panel coverage");
        static_assert(std::is_nothrow_move_assignable_v<SurfaceChargeTransport>,
            "atomic Field publication requires noexcept transport move assignment");
        std::vector<SurfaceChargeTransport> staged;staged.reserve(fields.size());
        for(std::size_t component=0;component<fields.size();++component) {
            auto* field=fields[component];
            std::vector<double> q(charge.begin()+static_cast<std::ptrdiff_t>(offsets[component]),
                charge.begin()+static_cast<std::ptrdiff_t>(offsets[component+1]));
            staged.push_back(SurfaceChargeTransport::fromIntegratedCharge(
                field->charge_.panelMeasure(),q,field->panels_.closed()));
            const auto source=field->source(staged.back());
            for(double value:source.cellIntegratedCharge)
                require(std::isfinite(value),"charge Field preparation: reconstructed component source is nonfinite");
        }
        return staged;
    }
    void verifyFields() const
    {
        require(static_cast<bool>(accepted) && static_cast<bool>(ownerToken),
            "charge accepted state: missing accepted receipt or writable ownership");
        require(sameClock(accepted->representationLayer,solver.geometrySourceTimeLayer())
            && accepted->representationIdentity==solver.momentOperator().geometry().geometryIdentity()
            && (!previous || (sameClock(previous->representationLayer,accepted->representationLayer)
                && previous->representationIdentity==accepted->representationIdentity)),
            "charge accepted state: history representation differs from the actual geometry binding");
        for(std::size_t component=0;component<fields.size();++component) {
            require(fields[component]->chargeTimeOwner_.lock().get()==ownerToken.get(),
                "charge accepted state: Field clock/history ownership binding changed");
            const auto& q=fields[component]->transport().panelCharge();
            require(q.size()==offsets[component+1]-offsets[component]
                && (q.empty() || std::memcmp(q.data(),accepted->charge.data()+offsets[component],
                    q.size()*sizeof(double))==0),
                "charge accepted state: Field Q changed outside its bound accepted clock/history");
        }
    }
};
ElectricSurfaceChargeTimeState::ElectricSurfaceChargeTimeState(std::unique_ptr<Data> data):data_(std::move(data)){}
ElectricSurfaceChargeTimeState::~ElectricSurfaceChargeTimeState()=default;
ElectricSurfaceChargeTimeState::ElectricSurfaceChargeTimeState(ElectricSurfaceChargeTimeState&&) noexcept=default;
ElectricSurfaceChargeTimeState& ElectricSurfaceChargeTimeState::operator=(ElectricSurfaceChargeTimeState&&) noexcept=default;
const ElectricSurfaceChargeTimeState::Data& ElectricSurfaceChargeTimeState::checked() const
{
    if(!data_)throw std::logic_error("[ElectricSurfaceChargeTimeState] moved-from accepted state");
    return *data_;
}
ElectricSurfaceChargeTimeState ElectricSurfaceChargeTimeState::initialize(const ElectricFixedGeometryChargeSolver& solver,
    const std::vector<SurfaceChargeField*>& fields,const Interface::HermiteGalsTimeLayer& initial)
{
    const auto& preparation=solver.checked();(void)initial.time();
    require(sameClock(initial,solver.geometrySourceTimeLayer()),
        "charge initialize: initial physical clock must exactly match the fixed geometry source publication");
    const auto& bound=preparation.moments.checked();
    require(fields.size()==bound.fields.size(),"charge initialize: component coverage differs from actual P_Q owners");
    auto data=std::make_unique<Data>(solver,fields);
    data->offsets.push_back(0);
    auto receipt=std::make_shared<AcceptedChargeLevel>();receipt->layer=initial;
    receipt->representationLayer=solver.geometrySourceTimeLayer();
    receipt->representationIdentity=solver.momentOperator().geometry().geometryIdentity();
    for(std::size_t component=0;component<fields.size();++component) {
        require(fields[component] && fields[component]==bound.fields[component],
            "charge initialize: Field must be the actual P_Q component/source/layout owner");
        require(fields[component]->chargeTimeOwner_.expired(),
            "charge initialize: Field already has a live accepted clock/history owner");
        for(std::size_t other=0;other<component;++other)
            require(fields[other]!=fields[component],"charge initialize: duplicate mutable Field owner");
        const auto& q=fields[component]->transport().panelCharge();
        const auto source=fields[component]->source(fields[component]->charge_);
        for(double value:source.cellIntegratedCharge)
            require(std::isfinite(value),"charge initialize: nonfinite reconstructed component source");
        receipt->charge.insert(receipt->charge.end(),q.begin(),q.end());
        data->offsets.push_back(receipt->charge.size());
    }
    require(receipt->charge.size()==preparation.moments.panelCount(),
        "charge initialize: Field panel count differs from actual native P_Q ordering");
    data->accepted=std::move(receipt);
    // Weak markers do not extend this state or its accepted-clock ownership.
    // All allocations and component checks completed before any marker changes.
    for(auto* field:fields)field->chargeTimeOwner_=data->ownerToken;
    return ElectricSurfaceChargeTimeState(std::move(data));
}
Interface::HermiteGalsTimeLayer ElectricSurfaceChargeTimeState::timeLayer() const{return checked().accepted->layer;}
const Interface::HermiteGalsTimeLayer* ElectricSurfaceChargeTimeState::previousTimeLayer() const
{const auto& data=checked();return data.previous?&data.previous->layer:nullptr;}
const std::vector<double>& ElectricSurfaceChargeTimeState::integratedCharge() const{return checked().accepted->charge;}
const std::vector<double>* ElectricSurfaceChargeTimeState::previousIntegratedCharge() const
{const auto& data=checked();return data.previous?&data.previous->charge:nullptr;}

const ElectricFixedGeometryChargeSolver& ElectricSurfaceChargeTimeState::solver() const
{return checked().solver;}
Interface::HermiteGalsTimeLayer ElectricSurfaceChargeTimeState::geometrySourceTimeLayer() const
{return checked().solver.geometrySourceTimeLayer();}
std::size_t ElectricSurfaceChargeTimeState::componentCount() const
{return checked().fields.size();}
const SurfaceChargeField& ElectricSurfaceChargeTimeState::field(std::size_t component) const
{
    const auto& data=checked();data.verifyFields();
    require(component<data.fields.size(),"charge Field query: invalid component");
    return *data.fields[component];
}

struct ElectricFixedGeometryChargeStep::Data
{
    Data(std::vector<double> q,ElectricPhaseResult field,ElectricPanelCurrentOperator ohmic,
        Interface::HermiteGalsTimeLayer physical,Interface::HermiteGalsTimeLayer source)
        :charge(std::move(q)),potential(std::move(field)),current(std::move(ohmic)),
         solutionLayer(physical),geometryLayer(source){}
    std::vector<double> charge;
    ElectricPhaseResult potential;
    ElectricPanelCurrentOperator current;
    Interface::HermiteGalsTimeLayer solutionLayer,geometryLayer;
    ElectricChargeTimeReport report;
    std::shared_ptr<const void> solverToken;
    std::shared_ptr<const AcceptedChargeLevel> sourceReceipt;
};
ElectricFixedGeometryChargeStep::ElectricFixedGeometryChargeStep(std::shared_ptr<const Data> data)
    :data_(std::move(data)){}
const ElectricFixedGeometryChargeStep::Data& ElectricFixedGeometryChargeStep::checked() const
{
    if(!data_)throw std::logic_error("[ElectricFixedGeometryChargeStep] moved-from candidate");
    return *data_;
}
const std::vector<double>& ElectricFixedGeometryChargeStep::integratedCharge() const{return checked().charge;}
const ElectricPhaseResult& ElectricFixedGeometryChargeStep::potential() const{return checked().potential;}
const ElectricPanelCurrentOperator& ElectricFixedGeometryChargeStep::current() const{return checked().current;}
Interface::HermiteGalsTimeLayer ElectricFixedGeometryChargeStep::timeLayer() const{return checked().solutionLayer;}
Interface::HermiteGalsTimeLayer ElectricFixedGeometryChargeStep::geometrySourceTimeLayer() const{return checked().geometryLayer;}
const ElectricChargeTimeReport& ElectricFixedGeometryChargeStep::report() const{return checked().report;}

ElectricFixedGeometryChargeStep ElectricFixedGeometryChargeSolver::advance(const ElectricSurfaceChargeTimeState& state,
    const Interface::HermiteGalsTimeLayer& target,const ElectricChargeTimeControls& controls) const
{
    return advanceCandidate(state,target,controls,false);
}
ElectricFixedGeometryChargeStep ElectricFixedGeometryChargeSolver::advanceCandidate(
    const ElectricSurfaceChargeTimeState& state,const Interface::HermiteGalsTimeLayer& target,
    const ElectricChargeTimeControls& controls,bool midpoint) const
{
    const auto totalBegin=TimingClock::now();
    const auto& prepared=checked();const auto& accepted=state.checked();
    validateChargeTimeControls(controls);
    require(accepted.solver.data_.get()==data_.get(),
        "charge time advance: accepted state belongs to another immutable preparation");
    accepted.verifyFields();
    ElectricChargeTimeReport report;report.dt=validateChargeTarget(accepted.accepted->layer,target,midpoint);
    report.a0=1;report.a1=-1;report.temporalOrder=1;
    if(accepted.previous) {
        report.previousDt=accepted.accepted->layer.end-accepted.accepted->layer.begin;
        require(std::isfinite(report.previousDt) && report.previousDt>0,
            "charge time advance: accepted history has invalid previous increment");
        require(accepted.accepted->layer.begin==accepted.previous->layer.time()
            && accepted.accepted->layer.step==accepted.previous->layer.step+1,
            "charge time advance: actual accepted history clocks are not consecutive");
        if(controls.method==ElectricChargeTimeMethod::VariableBdf2) {
            const double ratio=report.dt/report.previousDt;
            // Compare increments with the endpoint arithmetic error, rather
            // than clipping their actual ratio. A tiny dt at a large clock
            // must not turn an unresolved ratio into a relaxed stability bound.
            const double intervalRoundoff=8*std::numeric_limits<double>::epsilon()
                *(std::abs(midpoint?target.time():target.end)+std::abs(target.begin)
                    +controls.maximumStepRatio*(std::abs(accepted.accepted->layer.end)
                        +std::abs(accepted.accepted->layer.begin)));
            const double ratioExcess=ratio-controls.maximumStepRatio;
            const bool representedBoundary=ratioExcess>0
                && ratioExcess<=std::sqrt(std::numeric_limits<double>::epsilon())
                && report.dt-controls.maximumStepRatio*report.previousDt<=intervalRoundoff;
            require(std::isfinite(ratio) && ratio>0
                && (ratio<=controls.maximumStepRatio || representedBoundary),
                "charge time advance: BDF2 increment ratio exceeds its stable declared bound");
            report.a0=(1+2*ratio)/(1+ratio);report.a1=-(1+ratio);
            report.a2=ratio*ratio/(1+ratio);report.temporalOrder=2;
        }
    }
    const auto size=static_cast<Eigen::Index>(accepted.accepted->charge.size());
    const Eigen::Map<const Eigen::VectorXd> qn(accepted.accepted->charge.data(),size);
    Eigen::VectorXd rhs=-report.a1*qn+report.dt*prepared.fixedCurrent;
    if(report.a2!=0) {
        require(accepted.previous->charge.size()==accepted.accepted->charge.size(),
            "charge time advance: previous integrated history has foreign material panel coverage");
        rhs.noalias()-=report.a2*Eigen::Map<const Eigen::VectorXd>(accepted.previous->charge.data(),size);
    }
    const auto action=[&](const Eigen::VectorXd& q) -> Eigen::VectorXd {
        return report.a0*q-report.dt*prepared.linearCurrent(q,&report.electricSolves,
            &report.timing.schurElectricSolveSeconds);
    };
    const auto schurBegin=TimingClock::now();
    const auto q=solveChargeGmres(rhs,qn,action,controls,report);
    report.timing.schurSeconds=elapsedSeconds(schurBegin);
    const auto momentLoadBegin=TimingClock::now();
    std::vector<double> candidate(q.data(),q.data()+q.size());
    const auto load=prepared.moments.apply(candidate);
    report.timing.momentLoadSeconds=elapsedSeconds(momentLoadBegin);
    const auto finalSolveBegin=TimingClock::now();
    const Eigen::VectorXd electricRhs=Eigen::Map<const Eigen::VectorXd>(
        prepared.backend->fixedRhs.data(),static_cast<Eigen::Index>(prepared.backend->fixedRhs.size()))
        +prepared.injection*q;
    const auto u=prepared.electricSolve(electricRhs,&report.electricSolves);
    report.timing.finalElectricSolveSeconds=elapsedSeconds(finalSolveBegin);
    const auto recoveryBegin=TimingClock::now();
    std::vector<double> skeleton(u.data(),u.data()+u.size());
    auto potential=ElectricPhaseResult::recoverPrepared(prepared.base,load,skeleton,target);
    report.timing.recoverySeconds=elapsedSeconds(recoveryBegin);
    const auto currentBegin=TimingClock::now();
    auto current=potential.preparePanelCurrentBound(prepared.moments,prepared.current.conductivityByPhase(),
        &prepared.current);
    report.timing.currentSeconds=elapsedSeconds(currentBegin);
    const auto balanceBegin=TimingClock::now();
    const auto& raw=current.panelCurrent();
    long double totalQ=0,totalQn=0,totalPrevious=0,totalI=0,totalScale=0;
    for(std::size_t panel=0;panel<candidate.size();++panel) {
        const long double old=accepted.accepted->charge[panel];
        const long double previous=report.a2!=0?accepted.previous->charge[panel]:0;
        const long double t0=static_cast<long double>(report.a0)*candidate[panel];
        const long double t1=static_cast<long double>(report.a1)*old;
        const long double t2=static_cast<long double>(report.a2)*previous;
        const long double ti=static_cast<long double>(report.dt)*raw[panel];
        const long double residual=t0+t1+t2-ti;
        const long double scale=std::abs(t0)+std::abs(t1)+std::abs(t2)+std::abs(ti);
        report.maxRelativePanelEquationResidual=std::max(report.maxRelativePanelEquationResidual,
            normalized(static_cast<double>(std::abs(residual)),static_cast<double>(scale)));
        totalQ+=candidate[panel];totalQn+=old;totalPrevious+=previous;totalI+=raw[panel];totalScale+=scale;
    }
    const long double integrated=report.a0*totalQ+report.a1*totalQn+report.a2*totalPrevious-report.dt*totalI;
    report.integratedEquationResidual=static_cast<double>(integrated);
    report.relativeTotalEquationResidual=normalized(std::abs(report.integratedEquationResidual),
        static_cast<double>(totalScale));
    require(report.maxRelativePanelEquationResidual<=controls.relativeChargeBalanceTolerance
        && report.relativeTotalEquationResidual<=controls.relativeChargeBalanceTolerance,
        "charge time advance: original panel/total current equation exceeds its independent conservation tolerance");
    report.timing.balanceSeconds=elapsedSeconds(balanceBegin);
    auto result=std::make_shared<ElectricFixedGeometryChargeStep::Data>(std::move(candidate),
        std::move(potential),std::move(current),target,geometrySourceTimeLayer());
    report.timing.totalSeconds=elapsedSeconds(totalBegin);
    result->report=report;result->solverToken=data_;result->sourceReceipt=accepted.accepted;
    return ElectricFixedGeometryChargeStep(std::move(result));
}
void ElectricSurfaceChargeTimeState::commit(const ElectricFixedGeometryChargeStep& step)
{
    const auto& currentState=checked();const auto& candidate=step.checked();
    require(candidate.solutionLayer.position==Interface::HermiteGalsTimeLayer::Position::End,
        "charge commit: only an End candidate can be accepted; Midpoint stages are read-only");
    require(candidate.solverToken.get()==currentState.solver.data_.get()
        && candidate.sourceReceipt.get()==currentState.accepted.get(),
        "charge commit: stale/cross-state candidate or foreign immutable preparation");
    currentState.verifyFields();
    (void)validateChargeTarget(currentState.accepted->layer,candidate.solutionLayer);
    require(sameClock(candidate.potential.timeLayer(),candidate.solutionLayer)
        && sameClock(candidate.current.timeLayer(),candidate.solutionLayer)
        && sameClock(candidate.geometryLayer,currentState.solver.geometrySourceTimeLayer())
        && candidate.potential.surfaceChargeLoad()!=nullptr
        && sameChargeBits(candidate.charge,candidate.potential.surfaceChargeLoad()->panelCharge()),
        "charge commit: candidate Q/V/current clock or bound charge load is inconsistent");
    auto staged=currentState.prepareFieldCharge(candidate.charge);
    auto receipt=std::make_shared<AcceptedChargeLevel>();
    receipt->layer=candidate.solutionLayer;receipt->charge=candidate.charge;
    receipt->representationLayer=candidate.geometryLayer;
    receipt->representationIdentity=candidate.potential.geometry().geometryIdentity();
    // All throwing work is complete. Original component owners publish Q and
    // immutable accepted/history receipts together; no trial ever edits Field.
    for(std::size_t component=0;component<currentState.fields.size();++component)
        currentState.fields[component]->charge_=std::move(staged[component]);
    data_->previous=data_->accepted;data_->accepted=std::move(receipt);
}


struct ElectricMovingGeometryChargeStep::Data
{
    Data(std::shared_ptr<ElectricSurfaceChargeTimeState> state,ElectricFixedGeometryChargeStep candidate)
        :staged(std::move(state)),final(std::move(candidate)){}
    std::shared_ptr<ElectricSurfaceChargeTimeState> staged;
    std::shared_ptr<const std::vector<std::shared_ptr<SurfaceChargeField>>> retainedFields;
    ElectricFixedGeometryChargeStep final;
    // Weak provenance never prolongs the old mutable Field lease.
    std::weak_ptr<const void> sourceOwnerToken,sourceSolverToken;
    Interface::HermiteGalsTimeLayer sourceGeometryLayer;
    std::shared_ptr<const AcceptedChargeLevel> sourceReceipt,mappedAccepted,mappedPrevious;
    std::vector<SurfaceChargeFieldStep> remapReports;
    std::vector<std::vector<SurfaceChargeStepReport>> historyReports;
};
ElectricMovingGeometryChargeStep::ElectricMovingGeometryChargeStep(std::shared_ptr<const Data> data)
    :data_(std::move(data)){}
const ElectricMovingGeometryChargeStep::Data& ElectricMovingGeometryChargeStep::checked() const
{
    if(!data_)throw std::logic_error("[ElectricMovingGeometryChargeStep] moved-from candidate");
    return *data_;
}
const ElectricFixedGeometryChargeStep& ElectricMovingGeometryChargeStep::step() const
{return checked().final;}
const std::vector<double>& ElectricMovingGeometryChargeStep::mappedIntegratedCharge() const
{return checked().mappedAccepted->charge;}
const std::vector<double>* ElectricMovingGeometryChargeStep::mappedPreviousIntegratedCharge() const
{const auto& data=checked();return data.mappedPrevious?&data.mappedPrevious->charge:nullptr;}
Interface::HermiteGalsTimeLayer ElectricMovingGeometryChargeStep::sourceTimeLayer() const
{return checked().mappedAccepted->layer;}
Interface::HermiteGalsTimeLayer ElectricMovingGeometryChargeStep::sourceGeometryTimeLayer() const
{return checked().sourceGeometryLayer;}
const Interface::HermiteGalsTimeLayer* ElectricMovingGeometryChargeStep::previousSourceTimeLayer() const
{const auto& data=checked();return data.mappedPrevious?&data.mappedPrevious->layer:nullptr;}
Interface::HermiteGalsTimeLayer ElectricMovingGeometryChargeStep::geometrySourceTimeLayer() const
{return checked().final.geometrySourceTimeLayer();}
std::size_t ElectricMovingGeometryChargeStep::componentCount() const
{return checked().retainedFields->size();}
const SurfaceChargeField& ElectricMovingGeometryChargeStep::field(std::size_t component) const
{
    const auto& data=checked();
    require(component<data.retainedFields->size(),"moving charge Field query: invalid component");
    return *(*data.retainedFields)[component];
}
const std::vector<SurfaceChargeFieldStep>& ElectricMovingGeometryChargeStep::remapReports() const
{return checked().remapReports;}
const std::vector<std::vector<SurfaceChargeStepReport>>& ElectricMovingGeometryChargeStep::historyRemapReports() const
{return checked().historyReports;}

ElectricMovingGeometryChargeStep ElectricSurfaceChargeTimeState::prepareMoving(const QuadtreeMesh& mesh,
    const ElectricPhaseGeometry& geometry,const std::vector<SurfaceChargeMaterialTransfer>& transfers,
    const Interface::HermiteGalsTimeLayer& target,const ElectricPhaseProblem& problem,
    const std::array<double,2>& conductivity,const ElectricSurfaceChargeMomentControls& momentControls,
    const ElectricCartesianControls& electricControls,const ElectricChargeTimeControls& timeControls) const
{
    return prepareMovingCandidate(mesh,geometry,transfers,target,problem,conductivity,
        momentControls,electricControls,timeControls,false);
}
ElectricMovingGeometryChargeStep ElectricSurfaceChargeTimeState::prepareMovingStage(const QuadtreeMesh& mesh,
    const ElectricPhaseGeometry& geometry,const std::vector<SurfaceChargeMaterialTransfer>& transfers,
    const Interface::HermiteGalsTimeLayer& target,const ElectricPhaseProblem& problem,
    const std::array<double,2>& conductivity,const ElectricSurfaceChargeMomentControls& momentControls,
    const ElectricCartesianControls& electricControls,const ElectricChargeTimeControls& timeControls) const
{
    return prepareMovingCandidate(mesh,geometry,transfers,target,problem,conductivity,
        momentControls,electricControls,timeControls,true);
}
ElectricMovingGeometryChargeStep ElectricSurfaceChargeTimeState::prepareMovingCandidate(const QuadtreeMesh& mesh,
    const ElectricPhaseGeometry& geometry,const std::vector<SurfaceChargeMaterialTransfer>& transfers,
    const Interface::HermiteGalsTimeLayer& target,const ElectricPhaseProblem& problem,
    const std::array<double,2>& conductivity,const ElectricSurfaceChargeMomentControls& momentControls,
    const ElectricCartesianControls& electricControls,const ElectricChargeTimeControls& timeControls,
    bool midpoint) const
{
    const auto& source=checked();source.verifyFields();
    validateChargeTimeControls(timeControls);
    const double dt=validateChargeTarget(source.accepted->layer,target,midpoint);
    require(sameClock(geometry.timeLayer(),target),
        "moving charge preparation: target geometry must have the actual complete target candidate clock");
    geometry.validateMesh(mesh);
    if(source.previous)(void)validateChargeTarget(source.previous->layer,source.accepted->layer);
    const auto& sourceMoments=source.solver.momentOperator().checked();
    require(sourceMoments.fields.size()==source.fields.size()
        && transfers.size()==source.fields.size(),
        "moving charge preparation: transfers do not cover actual P_Q component owners");
    auto owners=std::make_shared<std::vector<std::shared_ptr<SurfaceChargeField>>>();
    owners->reserve(source.fields.size());
    std::vector<SurfaceChargeField*> targetFields;
    std::vector<const SurfaceChargeField*> targetReads;
    std::vector<std::size_t> offsets{0};
    std::vector<SurfaceChargeFieldStep> reports;
    std::vector<std::vector<SurfaceChargeStepReport>> historyReports;
    auto mapped=std::make_shared<AcceptedChargeLevel>();
    mapped->layer=source.accepted->layer;mapped->representationLayer=geometry.timeLayer();
    mapped->representationIdentity=geometry.geometryIdentity();
    std::shared_ptr<AcceptedChargeLevel> previous;
    if(source.previous) {
        previous=std::make_shared<AcceptedChargeLevel>();
        previous->layer=source.previous->layer;previous->representationLayer=geometry.timeLayer();
        previous->representationIdentity=geometry.geometryIdentity();
    }
    for(std::size_t component=0;component<source.fields.size();++component) {
        const auto& transfer=transfers[component];
        require(sourceMoments.fields[component]==source.fields[component],
            "moving charge preparation: source Field differs from the actual P_Q owner");
        require(transfer.step()==target.step && transfer.sourceTime()==source.accepted->layer.time()
            && transfer.targetTime()==target.time() && transfer.dt()==dt,
            "moving charge preparation: native transfer differs from the actual physical interval");
        std::vector<std::vector<double>> history;
        if(previous)history.emplace_back(
            source.previous->charge.begin()+static_cast<std::ptrdiff_t>(source.offsets[component]),
            source.previous->charge.begin()+static_cast<std::ptrdiff_t>(source.offsets[component+1]));
        auto result=source.fields[component]->remapPrepared(mesh,transfer,history);
        require(static_cast<bool>(result.field) && result.history.size()==history.size()
            && result.historyReports.size()==history.size(),
            "moving charge preparation: native remap did not return the actual complete history");
        std::shared_ptr<SurfaceChargeField> field(std::move(result.field));
        const auto& q=field->transport().panelCharge();
        mapped->charge.insert(mapped->charge.end(),q.begin(),q.end());
        offsets.push_back(mapped->charge.size());
        if(previous) {
            require(result.history[0].size()==q.size(),
                "moving charge preparation: mapped previous history has foreign target panels");
            previous->charge.insert(previous->charge.end(),result.history[0].begin(),result.history[0].end());
        }
        targetFields.push_back(field.get());targetReads.push_back(field.get());
        owners->push_back(std::move(field));
        reports.push_back(std::move(result.report));
        historyReports.push_back(std::move(result.historyReports));
    }
    // This preparation checks real scalar/curve pointers, oriented native panel
    // order and weighted moments against the actual target publication, not its
    // display points. The retained owners outlive every P_Q/result/solver copy.
    auto moments=ElectricSurfaceChargeMomentOperator::prepareRetained(
        geometry,targetReads,target,momentControls,owners);
    const auto load=moments.apply(mapped->charge);
    auto base=solveElectricPhaseElectrostatics(geometry,problem,load,electricControls);
    auto solver=ElectricFixedGeometryChargeSolver::prepare(base,moments,conductivity);
    auto data=std::make_unique<Data>(std::move(solver),std::move(targetFields));
    data->ownedFields=owners;data->offsets=std::move(offsets);
    data->accepted=mapped;data->previous=previous;
    for(auto* field:data->fields) {
        require(field->chargeTimeOwner_.expired(),"moving charge preparation: target Field already has a writable owner");
        field->chargeTimeOwner_=data->ownerToken;
    }
    auto staged=std::shared_ptr<ElectricSurfaceChargeTimeState>(new ElectricSurfaceChargeTimeState(std::move(data)));
    auto final=staged->solver().advanceCandidate(*staged,target,timeControls,midpoint);
    if(midpoint) {
        // Publish only these privately owned target Fields. No stage accepted
        // receipt is created and the two physical histories remain unchanged.
        auto charge=staged->data_->prepareFieldCharge(final.integratedCharge());
        for(std::size_t component=0;component<staged->data_->fields.size();++component)
            staged->data_->fields[component]->charge_=std::move(charge[component]);
    } else staged->commit(final);
    auto candidate=std::make_shared<ElectricMovingGeometryChargeStep::Data>(
        midpoint?nullptr:staged,std::move(final));
    candidate->retainedFields=owners;
    candidate->sourceOwnerToken=source.ownerToken;candidate->sourceSolverToken=source.solver.data_;
    candidate->sourceGeometryLayer=source.accepted->representationLayer;
    candidate->sourceReceipt=source.accepted;candidate->mappedAccepted=std::move(mapped);
    candidate->mappedPrevious=std::move(previous);
    candidate->remapReports=std::move(reports);candidate->historyReports=std::move(historyReports);
    // A retained immutable candidate owns Field lifetime, never a writable
    // lease. The successful outer commit installs a fresh accepted owner.
    for(auto* field:staged->data_->fields)field->chargeTimeOwner_.reset();
    staged->data_->ownerToken.reset();
    return ElectricMovingGeometryChargeStep(std::move(candidate));
}
void ElectricSurfaceChargeTimeState::commit(const ElectricMovingGeometryChargeStep& step)
{
    const auto& source=checked();const auto& candidate=step.checked();
    require(candidate.staged
        && candidate.final.timeLayer().position==Interface::HermiteGalsTimeLayer::Position::End,
        "moving charge commit: only an End candidate can be accepted; Midpoint stages are read-only");
    require(candidate.sourceOwnerToken.lock().get()==source.ownerToken.get()
        && candidate.sourceSolverToken.lock().get()==source.solver.data_.get()
        && candidate.sourceReceipt.get()==source.accepted.get(),
        "moving charge commit: stale/cross-state candidate or foreign source preparation");
    source.verifyFields();
    const auto& final=candidate.final.checked();
    (void)validateChargeTarget(source.accepted->layer,final.solutionLayer);
    const auto& target=candidate.staged->checked();
    require(!target.ownerToken && target.accepted && target.previous
        && target.previous.get()==candidate.mappedAccepted.get()
        && sameClock(target.accepted->layer,final.solutionLayer)
        && sameClock(target.accepted->representationLayer,final.geometryLayer)
        && target.accepted->representationIdentity==target.solver.momentOperator().geometry().geometryIdentity()
        && target.previous->representationIdentity==target.accepted->representationIdentity
        && sameClock(target.previous->layer,source.accepted->layer)
        && sameClock(target.previous->representationLayer,final.geometryLayer)
        && sameClock(final.potential.timeLayer(),final.solutionLayer)
        && sameClock(final.current.timeLayer(),final.solutionLayer)
        && sameClock(final.geometryLayer,target.solver.geometrySourceTimeLayer())
        && final.solverToken.get()==target.solver.data_.get()
        && final.sourceReceipt.get()==candidate.mappedAccepted.get()
        && sameChargeBits(target.accepted->charge,final.charge),
        "moving charge commit: target Q/field/history/clock preparation is inconsistent");
    require(target.ownedFields && target.ownedFields.get()==candidate.retainedFields.get(),
        "moving charge commit: missing actual target component lifetime");
    for(std::size_t component=0;component<target.fields.size();++component) {
        const auto* field=target.fields[component];
        const auto& q=field->transport().panelCharge();
        require(field==(*target.ownedFields)[component].get() && field->chargeTimeOwner_.expired()
            && q.size()==target.offsets[component+1]-target.offsets[component]
            && (q.empty() || std::memcmp(q.data(),target.accepted->charge.data()+target.offsets[component],
                q.size()*sizeof(double))==0),
            "moving charge commit: target Field was edited or acquired by another writable owner");
    }
    auto replacement=std::make_unique<Data>(target);
    replacement->ownerToken=std::make_shared<const int>(0);
    // All checks and allocations precede this nonthrowing publication. Previous
    // is mapped Qn at its original physical clock; old source leases expire even
    // if callers retain candidates, potentials or current operators.
    for(auto* field:replacement->fields)field->chargeTimeOwner_=replacement->ownerToken;
    data_.swap(replacement);
}

}
