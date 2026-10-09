#include "Modules/Hydrodynamics/Detail/HydroNativeRZRowPreparation.h"
#include "Modules/Hydrodynamics/HydroFittedMomentum.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroFittedProjection.h"
#include <Eigen/SuperLUSupport>
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/IO/DictionaryMesh.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <cstring>
#include <optional>
#include <type_traits>
#include <limits>
#include <algorithm>

namespace {
void check(bool value,const char* why) { if(!value)throw std::runtime_error(why); }
template<class F> void rejects(F action) {
    try { action(); } catch(const std::invalid_argument&) { return; }
    throw std::runtime_error("invalid projection input accepted");
}
void near(double value,double expected,const char* why) {
    check(std::isfinite(value) && std::abs(value-expected)<=4e-12*std::max({1.,std::abs(value),std::abs(expected)}),why);
}
void run(const std::filesystem::path& root,bool refined,bool varying)
{
    DictionaryMesh input;input.Nx=8;input.Ny=10;input.maxCellNum=1024;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    if(refined) {
        DictionaryMesh::RefineBox box;box.min={.625,.3};box.max={.875,.7};box.level=1;mesh.applyRefineBoxes({box});
        check(mesh.cellIndexNum()>80,"real refinement required");
    }
    const auto density=[&](Vector2D p){return varying?1.2+.4*p.x+.3*p.y:1.2;};
    const auto viscosity=[&](Vector2D p){return varying?.04+.01*p.x+.005*p.y:.04;};
    CellField<double> rho(mesh,0.),mu(mesh,0.),fx(mesh,0.),fy(mesh,0.);
    VertField<double> vrho(mesh,0.),vmu(mesh,0.);
    for(int i=0;i<mesh.cellIndexNum();++i) {
        const auto p=mesh.cells(mesh.cellIndexToLabel[i]).centerCoord();
        rho.value[i]=density(p);mu.value[i]=viscosity(p);fx.value[i]=.015*p.x;fy.value[i]=.02+.01*p.y;
    }
    for(int i=0;i<mesh.vertIndexNum();++i) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[i]).coord();vrho.value[i]=density(p);vmu.value[i]=viscosity(p);
    }
    const int nr=mesh.vEdgeIndexNum(),count=nr+mesh.hEdgeIndexNum();const double dt=.003;
    Eigen::VectorXd current(count),previous(count);
    for(int i=0;i<count;++i) {
        const auto& edge=i<nr?mesh.vEdges(mesh.vEdgeIndexToLabel[i]):mesh.hEdges(mesh.hEdgeIndexToLabel[i-nr]);
        const auto p=edge.midpointCoord();
        current[i]=i<nr?.02*p.x*(1.-p.x)*(.2+p.y):.03*p.y*(.5+p.x);
        previous[i]=i<nr?.017*p.x*(1.-p.x)*(.3+p.y):.025*p.y*(.4+p.x);
        if(edge.boundary()!=Boundary::non && !(i>=nr && p.y==1.))current[i]=previous[i]=0.;
    }
    FaceSegmentGraph graph;HydroFittedDomain domain;HydroFittedMomentumSystem fitted;
    if(!refined) {
        graph.build(mesh);std::vector<HydroFittedDomain::Boundary> boundary(graph.segments().size());
        for(std::size_t s=0;s<boundary.size();++s)if(graph.segments()[s].L<0)
            boundary[s].type=graph.segments()[s].side==FaceSide::North?HydroFittedDomain::BoundaryType::PressureOutlet:
                graph.segments()[s].side==FaceSide::West?HydroFittedDomain::BoundaryType::Symmetry:HydroFittedDomain::BoundaryType::NoSlipWall;
        domain=HydroFittedDomain::build(mesh,graph,boundary,true);
        fitted=assembleHydroFittedMomentum(mesh,domain,dt,rho,mu,vrho,vmu,current,previous,fx,fy);
    }
    const auto name=std::string(refined?"refined_":"uniform_")+(varying?"variable":"constant");
    std::ofstream out(root/(name+".csv"));out<<std::setprecision(17);
    int accepted=0,rejected=0,coarseFine=0,radialRows=0,axialRows=0;
    for(int native=0;native<count;++native) {
        const bool radial=native<nr;
        const auto& edge=radial?mesh.vEdges(mesh.vEdgeIndexToLabel[native]):mesh.hEdges(mesh.hEdgeIndexToLabel[native-nr]);
        HydroNativeRZRow::Row row;
        const bool valid=HydroNativeRZRowDetail::regularRow(mesh,edge,radial,rho,mu,vmu,row);
        const bool hanging=edge.masterCell>=0 && edge.slaveCell>=0 && mesh.cells(edge.masterCell).level()!=mesh.cells(edge.slaveCell).level();
        if(hanging) {++coarseFine;check(!valid,"coarse/fine row entered regular arithmetic");}
        if(!valid) {++rejected;out<<native<<",rejected\n";continue;}
        ++accepted;if(radial)++radialRows;else ++axialRows;
        const auto p=edge.midpointCoord();const auto h=mesh.cells(edge.masterCell).size();const auto& g=row.metrics;
        // Independent annular geometry and linear-material midpoint references.
        near(g.volume,Pi*(std::pow(p.x+.5*h.x,2)-std::pow(p.x-.5*h.x,2))*h.y,"annular dual volume");
        near(g.areaEast,2.*Pi*(p.x+.5*h.x)*h.y,"outer cylindrical area");
        near(g.areaWest,2.*Pi*(p.x-.5*h.x)*h.y,"inner cylindrical area");
        near(g.areaNS,Pi*(std::pow(p.x+.5*h.x,2)-std::pow(p.x-.5*h.x,2)),"annular cap area");
        near(g.shearLowDistance,radial?h.y:h.x,"low native shear spacing");
        near(g.shearHighDistance,radial?h.y:h.x,"high native shear spacing");
        near(row.materials.rho,density(p),"face density interpolation");
        const auto& start=mesh.verts(edge.edgeToVert[0]);const auto& end=mesh.verts(edge.edgeToVert[1]);
        if(radial) {row.materials.rhoSouth=vrho.value[start.index()];row.materials.rhoNorth=vrho.value[end.index()];}
        else {row.materials.rhoWest=vrho.value[start.index()];row.materials.rhoEast=vrho.value[end.index()];}
        const auto now=HydroNativeRZRow::values(row.indices,[&](int i){return current[i];});
        const auto pre=HydroNativeRZRow::values(row.indices,[&](int i){return previous[i];});
        const int master=mesh.cells(edge.masterCell).index(),slave=mesh.cells(edge.slaveCell).index();
        double rhs=radial?HydroNativeRZRow::radialRhs(row,now,pre,fx.value[master],fx.value[slave],dt):
            HydroNativeRZRow::axialRhs(row,now,pre,fy.value[master],fy.value[slave],dt);
        std::map<int,double> coefficients;
        const auto emit=[&](int r,int c,double a){check(r==native,"native row numbering");coefficients[c]+=a;};
        if(radial)HydroNativeRZRow::appendRadialMatrix(row,dt,emit);else HydroNativeRZRow::appendAxialMatrix(row,dt,emit);
        if(radial) {
            const double hoop=.5*(mu.value[master]+mu.value[slave])*g.volume/(p.x*p.x);
            coefficients[native]+=hoop;rhs-=hoop*current[native];
        }
        out<<native<<','<<g.volume<<','<<g.areaWest<<','<<g.areaEast<<','<<g.areaNS<<','<<g.dx<<','<<g.dy<<','<<g.shearLowDistance<<','<<g.shearHighDistance<<','<<row.materials.rho<<','<<rhs;
        for(const auto& entry:coefficients)out<<','<<entry.first<<','<<entry.second;
        out<<'\n';
        if(!refined) {
            const int r=domain.velocityRowByNative()[native];check(r>=0,"shared interior row must be active");
            near(fitted.rhs[r],rhs,"fitted/native CN/AB2 RHS including hoop");
            near(fitted.densityByVelocityRow[r],row.materials.rho,"fitted/native projection density");
            for(int c=0;c<static_cast<int>(domain.nativeByVelocityRow().size());++c)
                near(fitted.matrix.coeff(r,c),coefficients[domain.nativeByVelocityRow()[c]],"fitted/native full regular momentum row");
        }
    }
    check(radialRows>0 && axialRows>0 && rejected>0,"both regular orientations and excluded rows required");
    if(refined)check(coarseFine>0,"actual coarse/fine exclusions required");
    if(!refined) {
        check(fitted.regularNativeRows==static_cast<std::size_t>(accepted) && fitted.boundaryPieceRows>0,"fitted routing changed");
        std::ofstream matrix(root/(name+"_fitted.csv"));matrix<<std::setprecision(17);
        for(int c=0;c<fitted.matrix.outerSize();++c)for(Eigen::SparseMatrix<double>::InnerIterator a(fitted.matrix,c);a;++a)matrix<<a.row()<<','<<a.col()<<','<<a.value()<<'\n';
        for(int r=0;r<fitted.rhs.size();++r)matrix<<"rhs,"<<r<<','<<fitted.rhs[r]<<','<<fitted.densityByVelocityRow[r]<<'\n';
    }
    std::cout<<name<<": cells="<<mesh.cellIndexNum()<<" regular="<<accepted<<" excluded="<<rejected<<" coarseFine="<<coarseFine<<'\n';
}
void projectionRows(const std::filesystem::path& root,bool refined,bool varying,double scale)
{
    DictionaryMesh input;input.Nx=8;input.Ny=10;input.maxCellNum=1024;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    if(refined) {
        DictionaryMesh::RefineBox box;box.min={.375,.3};box.max={.625,.7};box.level=1;
        mesh.applyRefineBoxes({box});
    }
    FaceSegmentGraph graph;graph.build(mesh);
    std::vector<HydroFittedDomain::Boundary> boundary(graph.segments().size());
    for(std::size_t i=0;i<boundary.size();++i)if(graph.segments()[i].L<0) {
        const auto side=graph.segments()[i].side;
        auto& b=boundary[i];
        b.type=side==FaceSide::North?HydroFittedDomain::BoundaryType::PressureOutlet:
            side==FaceSide::South?HydroFittedDomain::BoundaryType::PressureReservoir:
            side==FaceSide::West?HydroFittedDomain::BoundaryType::Symmetry:HydroFittedDomain::BoundaryType::NoSlipWall;
        if(side==FaceSide::South){b.incomingDensity=1.;b.incomingViscosity=.01;}
    }
    const auto domain=HydroFittedDomain::build(mesh,graph,boundary,true);
    check(domain.uniformLeaves()!=refined,"projection geometry strategy");
    const int count=static_cast<int>(domain.nativeByVelocityRow().size());
    std::vector<double> density(count);
    for(int i=0;i<count;++i) {
        const auto& face=domain.faces()[domain.nativeByVelocityRow()[i]];
        density[i]=varying?1.+.2*face.center.x+.1*face.center.y:1.;
    }
    const auto name=std::string("projection_")+(refined?"refined_":"uniform_")+(varying?"variable_":"constant_")+(scale==1.?"transfer":"step");
    std::ofstream out(root/(name+".csv"));out<<std::setprecision(17);
    Eigen::SuperLU<Eigen::SparseMatrix<double>> solver;
    Eigen::SparseMatrix<double> frozen;
    std::optional<HydroFittedProjection::System> projection;
    const double* coefficients=nullptr;
    const int* inner=nullptr;
    const int* outer=nullptr;
    for(int history=0;history<2;++history) {
        Eigen::VectorXd velocity(count);
        for(int i=0;i<count;++i)velocity[i]=.02*std::sin((i+1.)*(history+1.)*.17);
        for(int pass=0;pass<3;++pass) {
            if(!projection)projection.emplace(HydroFittedProjection::assemble(domain,scale,density,velocity));
            else projection->updateRightHandSide(velocity);
            const auto& system=projection->equations();
            if(history==0 && pass==0) {
                frozen=system.matrix;solver.compute(frozen);
                coefficients=system.matrix.valuePtr();inner=system.matrix.innerIndexPtr();outer=system.matrix.outerIndexPtr();
                check(solver.info()==Eigen::Success,"projection factorization");
                for(int c=0;c<frozen.outerSize();++c)for(Eigen::SparseMatrix<double>::InnerIterator a(frozen,c);a;++a)
                    out<<"matrix,"<<a.row()<<','<<a.col()<<','<<a.value()<<'\n';
                // Face masses are owned by this preparation, not borrowed
                // from the caller's material input.
                std::fill(density.begin(),density.end(),std::numeric_limits<double>::quiet_NaN());
            }
            check(system.matrix.nonZeros()==frozen.nonZeros()
                && system.matrix.valuePtr()==coefficients && system.matrix.innerIndexPtr()==inner && system.matrix.outerIndexPtr()==outer
                && std::memcmp(inner,frozen.innerIndexPtr(),frozen.nonZeros()*sizeof(int))==0
                && std::memcmp(outer,frozen.outerIndexPtr(),(frozen.outerSize()+1)*sizeof(int))==0
                && std::memcmp(system.matrix.valuePtr(),frozen.valuePtr(),frozen.nonZeros()*sizeof(double))==0,
                "projection history changed the matrix");
            const Eigen::VectorXd multiplier=solver.solve(system.rhs);
            check(solver.info()==Eigen::Success && multiplier.allFinite(),"projection backsolve");
            const double residual=(system.matrix*multiplier-system.rhs).stableNorm();
            check(residual<=1e-11*system.rhs.stableNorm(),"projection full linear residual");
            velocity=projection->correctVelocity(velocity,multiplier);
            const Eigen::VectorXd flux=domain.divergence()*velocity;
            for(int i=0;i<multiplier.size();++i)
                out<<history<<','<<pass<<",cell,"<<i<<','<<system.rhs[i]<<','<<multiplier[i]<<','<<flux[i]<<'\n';
            for(int i=0;i<velocity.size();++i)out<<history<<','<<pass<<",velocity,"<<i<<','<<velocity[i]<<'\n';
        }
    }
    const auto saved=projection->equations().rhs;
    rejects([&]{projection->updateRightHandSide(Eigen::VectorXd::Zero(count-1));});
    Eigen::VectorXd invalid=Eigen::VectorXd::Zero(count);
    invalid[0]=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{projection->updateRightHandSide(invalid);});
    Eigen::VectorXd area=Eigen::VectorXd::Zero(domain.divergence().rows());
    for(int c=0;c<domain.divergence().outerSize();++c)
        for(Eigen::SparseMatrix<double>::InnerIterator a(domain.divergence(),c);a;++a)area[a.row()]+=std::abs(a.value());
    Eigen::Index cell=0;check(area.maxCoeff(&cell)>1.,"overflow probe needs summed physical face area greater than one");
    invalid.setZero();
    for(int c=0;c<domain.divergence().outerSize();++c)
        for(Eigen::SparseMatrix<double>::InnerIterator a(domain.divergence(),c);a;++a)
            if(a.row()==cell)invalid[a.col()]=std::copysign(std::numeric_limits<double>::max(),a.value());
    rejects([&]{projection->updateRightHandSide(invalid);});
    check(std::memcmp(saved.data(),projection->equations().rhs.data(),saved.size()*sizeof(double))==0,
        "rejected predictor changed pressure rhs");
    const Eigen::VectorXd zero=Eigen::VectorXd::Zero(count);
    rejects([&]{projection->correctVelocity(zero,Eigen::VectorXd::Zero(saved.size()-1));});
    Eigen::VectorXd badPressure=Eigen::VectorXd::Zero(saved.size());badPressure[0]=std::numeric_limits<double>::infinity();
    rejects([&]{projection->correctVelocity(zero,badPressure);});
    rejects([&]{HydroFittedProjection::assemble(domain,scale,density,zero);});
    std::fill(density.begin(),density.end(),1.);
    rejects([&]{HydroFittedProjection::assemble(domain,0.,density,zero);});
    check(projection->equations().matrix.valuePtr()==coefficients
        && std::memcmp(coefficients,frozen.valuePtr(),frozen.nonZeros()*sizeof(double))==0,"rejection changed pressure matrix");
    out.flush();check(bool(out),"projection row output");
    std::cout<<name<<": six projections of two histories passed\n";
}
}
int main(int argc,char** argv) try {
    static_assert(!std::is_copy_constructible_v<HydroFittedProjection::System>);
    static_assert(!std::is_default_constructible_v<HydroFittedProjection::System>);
    static_assert(std::is_same_v<decltype(std::declval<const HydroFittedProjection::System&>().equations()),const HydroPressureSystem&>);
    const auto root=argc>1?std::filesystem::path(argv[1]):std::filesystem::current_path()/("rz_rows_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(root);
    for(bool refined:{false,true})for(bool varying:{false,true})run(root,refined,varying);
    for(bool refined:{false,true})for(bool varying:{false,true})for(double scale:{1.,.01})
        projectionRows(root,refined,varying,scale);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
