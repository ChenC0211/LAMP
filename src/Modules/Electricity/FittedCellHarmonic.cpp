#include "Modules/Electricity/FittedCellHarmonic.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace Quantum::Electricity {
namespace {
using Solver=PolarizationPoissonFVSolverCellCentered;
using BC=Solver::FittedBoundaryType;
void require(bool valid,const char* reason) {
    if(!valid) throw std::invalid_argument(std::string("[FittedCellHarmonic] ")+reason);
}
double ownerDistance(const QuadtreeMesh& mesh,const FaceSegment& face) {
    const auto center=mesh.cells(face.K).centerCoord();
    return std::abs(face.isHorizontal?face.faceCenter.y-center.y:face.faceCenter.x-center.x);
}
struct Problem {
    FaceSegmentGraph graph;
    std::vector<Solver::FittedBoundary> boundaries;
    NativeLeafDomain geometry;
    Problem(const QuadtreeMesh& mesh,FaceSegmentGraph g,std::vector<Solver::FittedBoundary> b,bool axis)
        : graph(std::move(g)),boundaries(std::move(b)),geometry(NativeLeafDomain::build(mesh,graph,axis)) {
        require(boundaries.size()==graph.segments().size(),"boundary/graph sizes differ");
        for(std::size_t j=0;j<boundaries.size();++j) {
            const auto& f=graph.segments()[j]; const auto& bc=boundaries[j];
            const bool interior=f.L>=0 && geometry.cellRowByIndex()[mesh.cells(f.L).index()]>=0;
            require(std::isfinite(bc.value) && (interior?bc.type==BC::Interior && bc.value==0.0
                :bc.type==BC::Dirichlet || bc.type==BC::HomogeneousNeumann),"invalid fixed electric boundary");
        }
    }
};
std::vector<double> coefficients(const QuadtreeMesh& mesh,const Problem& problem,const std::vector<double>& epsilon) {
    require(epsilon.size()==static_cast<std::size_t>(mesh.cellIndexNum()),"cell material size differs from mesh");
    for(int c:problem.geometry.cellIndexByRow())
        require(std::isfinite(epsilon[c]) && epsilon[c]>0.0,"active permittivity must be finite and positive");
    const auto& segments=problem.graph.segments();
    std::vector<double> result(segments.size());
    for(std::size_t j=0;j<segments.size();++j) {
        const auto& face=segments[j]; const int k=mesh.cells(face.K).index();
        const int l=face.L<0?-1:mesh.cells(face.L).index();
        const double epsilonK=epsilon[k];
        if(l>=0 && problem.geometry.cellRowByIndex()[l]>=0) {
            const auto neighborCenter=mesh.cells(face.L).centerCoord();
            const double dK=ownerDistance(mesh,face);
            const double dL=std::abs(face.isHorizontal?neighborCenter.y-face.faceCenter.y:neighborCenter.x-face.faceCenter.x);
            result[j]=(dK+dL)/(dK/epsilonK+dL/epsilon[l]);
        } else result[j]=epsilonK;
        require(std::isfinite(result[j]) && result[j]>0.0,"invalid distance-weighted face permittivity");
    }
    return result;
}
void validateSource(const QuadtreeMesh& mesh,const Problem& problem,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& source) {
    if(!source) return;
    require(source->radialAxisSymmetry()==problem.geometry.axisymmetric(),"scalar/electric coordinate measures differ");
    require(source->minCoord().x==mesh.minCoord.x && source->minCoord().y==mesh.minCoord.y
        && source->maxCoord().x==mesh.maxCoord.x && source->maxCoord().y==mesh.maxCoord.y,"scalar/electric bounds differ");
    if(source->usesNativeLeaves()) {
        require(source->nativeCellCount()==static_cast<std::size_t>(mesh.cellIndexNum()),"scalar native count differs");
        for(int c=0;c<mesh.cellIndexNum();++c) {
            const auto bounds=source->nativeCellBounds(c); const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
            const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
            require(bounds[0].x==lo.x && bounds[0].y==lo.y && bounds[1].x==hi.x && bounds[1].y==hi.y
                && source->nativeCellActive(c)==(problem.geometry.cellRowByIndex()[c]>=0),"scalar leaf/active domain differs");
        }
    } else {
        const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
        require(source->nx()==grid.nx() && source->ny()==grid.ny()
            && source->cellSize().x==grid.spacing().x && source->cellSize().y==grid.spacing().y,"scalar lattice differs");
        for(int j=0;j<grid.ny();++j) for(int i=0;i<grid.nx();++i) {
            const int c=mesh.cells(grid.cells()[i+grid.nx()*j]).index();
            require(source->cellActive(i,j)==(problem.geometry.cellRowByIndex()[c]>=0),"scalar active domain differs");
        }
    }
}
Solver::SolveReport solve(const QuadtreeMesh& mesh,const Problem& problem,CellField<double>& potential,
    const std::vector<double>& epsilon,const FittedCellHarmonicControls& controls) {
    require(std::isfinite(controls.linearTolerance) && controls.linearTolerance>0.0 && controls.linearTolerance<1.0
        && std::isfinite(controls.potentialBoundsTolerance) && controls.potentialBoundsTolerance>0.0
        && std::isfinite(controls.potentialScale) && controls.potentialScale>0.0
        && std::isfinite(controls.lowerPotential) && std::isfinite(controls.upperPotential)
        && controls.lowerPotential<=controls.upperPotential,"invalid electric solve or bound controls");
    CellField<double> rhs(mesh,0.0); Solver solver; Solver::Options options;
    options.axisymmetric=problem.geometry.axisymmetric(); options.requireSuperLU=true; options.useSuperLU=true;
    options.useCG=options.useIC=false; options.relTol=controls.linearTolerance;
    Solver::SolveReport report;
    if(!solver.solveFittedBoundary(mesh,problem.graph,rhs,potential,epsilon,problem.boundaries,options,&report))
        throw std::runtime_error("[FittedCellHarmonic] electric solve failed: "+report.message);
    const double tolerance=controls.potentialBoundsTolerance*std::max({1.0,std::abs(controls.lowerPotential),std::abs(controls.upperPotential)});
    require(std::isfinite(tolerance),"potential bound tolerance overflow");
    for(int k=0;k<mesh.cellIndexNum();++k) {
        const double voltage=potential.value[k]*controls.potentialScale;
        require(std::isfinite(voltage) && voltage>=controls.lowerPotential-tolerance && voltage<=controls.upperPotential+tolerance,
            "potential is nonfinite or violates the Dirichlet maximum principle");
    }
    return report;
}
}
struct FittedCellHarmonicState::Data {
    std::shared_ptr<const Problem> problem;
    FittedElectricPotential potential;
    std::vector<double> epsilonCell,epsilonSegment;
    Solver::SolveReport report;
    Stage stage=Stage::Prepared;
    double time=0.0;
};
struct FittedCellHarmonicAccess {
    static const auto& read(const FittedCellHarmonicState& state) {
        if(!state.data_) throw std::logic_error("cell-harmonic state has not been prepared");
        return *state.data_;
    }
    static auto create() {return std::make_shared<FittedCellHarmonicState::Data>();}
    static FittedCellHarmonicState publish(std::shared_ptr<FittedCellHarmonicState::Data> data) {
        FittedCellHarmonicState result;result.data_=std::move(data);return result;
    }
    static auto working(const QuadtreeMesh& mesh,const FittedPotentialInput* seed) {
        if(!seed) return std::make_shared<CellField<double>>(mesh,0.0);
        require(seed->size()==mesh.cellIndexNum(),"seed storage differs from mesh");
        for(int k=0;k<mesh.cellIndexNum();++k)
            require(seed->label(k)==mesh.cellIndexToLabel[k] && std::isfinite((*seed)[k]),"seed labels/values are invalid");
        return std::make_shared<CellField<double>>(*seed->field_);
    }
    static auto freeze(std::shared_ptr<CellField<double>> field,std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source) {
        return FittedElectricPotential(std::move(field),std::move(source));
    }
};
FittedCellHarmonicState::Stage FittedCellHarmonicState::stage() const {return FittedCellHarmonicAccess::read(*this).stage;}
const FittedElectricPotential& FittedCellHarmonicState::potential() const {return FittedCellHarmonicAccess::read(*this).potential;}
const NativeLeafDomain& FittedCellHarmonicState::geometry() const {return FittedCellHarmonicAccess::read(*this).problem->geometry;}
const FaceSegmentGraph& FittedCellHarmonicState::graph() const {return FittedCellHarmonicAccess::read(*this).problem->graph;}
const std::vector<Solver::FittedBoundary>& FittedCellHarmonicState::boundaries() const {return FittedCellHarmonicAccess::read(*this).problem->boundaries;}
const std::vector<double>& FittedCellHarmonicState::epsilonByCell() const {return FittedCellHarmonicAccess::read(*this).epsilonCell;}
const std::vector<double>& FittedCellHarmonicState::epsilonBySegment() const {return FittedCellHarmonicAccess::read(*this).epsilonSegment;}
const Solver::SolveReport& FittedCellHarmonicState::report() const {return FittedCellHarmonicAccess::read(*this).report;}
double FittedCellHarmonicState::time() const {return FittedCellHarmonicAccess::read(*this).time;}
void FittedCellHarmonicState::validateMesh(const QuadtreeMesh& mesh) const {geometry().validateMeshIdentity(mesh);}
std::vector<FittedElectricFaceValue> FittedCellHarmonicState::faceValues(const QuadtreeMesh& mesh) const {
    validateMesh(mesh); require(stage()==Stage::Solved && report().success,"face diagnostics require a solved state");
    constexpr double twoPi=6.283185307179586476925286766559;
    std::vector<FittedElectricFaceValue> result;result.reserve(graph().segments().size());
    for(std::size_t j=0;j<graph().segments().size();++j) {
        const auto& segment=graph().segments()[j]; const auto& bc=boundaries()[j];
        const double owner=potential()[mesh.cells(segment.K).index()];
        const double distance=bc.type==BC::Interior?segment.df:ownerDistance(mesh,segment);
        const double drop=bc.type==BC::Interior?owner-potential()[mesh.cells(segment.L).index()]
            :bc.type==BC::Dirichlet?owner-bc.value:0.0;
        const double charge=geometry().axisymmetric()
            ?epsilonBySegment()[j]*twoPi*segment.faceCenter.x*segment.Af*drop/distance
            :epsilonBySegment()[j]*segment.Af*drop/distance;
        require(distance>0.0 && std::isfinite(drop) && std::isfinite(charge),"invalid solved electric face diagnostic");
        result.push_back({drop,distance,charge});
    }
    return result;
}
FittedCellHarmonicState prepareFittedCellHarmonic(const QuadtreeMesh& mesh,FaceSegmentGraph graph,
    std::vector<Solver::FittedBoundary> boundaries,std::vector<double> epsilon,bool axis,double inactivePotential,
    const FittedPotentialInput* seed,const FittedCellHarmonicControls* controls) {
    require(std::isfinite(inactivePotential),"inactive potential must be finite");
    auto data=FittedCellHarmonicAccess::create();
    data->problem=std::make_shared<const Problem>(mesh,std::move(graph),std::move(boundaries),axis);
    data->epsilonCell=std::move(epsilon);data->epsilonSegment=coefficients(mesh,*data->problem,data->epsilonCell);
    auto field=FittedCellHarmonicAccess::working(mesh,seed);
    for(int k=0;k<mesh.cellIndexNum();++k) if(data->problem->geometry.cellRowByIndex()[k]<0) field->value[k]=inactivePotential;
    if(controls) {data->report=solve(mesh,*data->problem,*field,data->epsilonSegment,*controls);data->stage=FittedCellHarmonicState::Stage::Solved;}
    data->potential=FittedCellHarmonicAccess::freeze(std::move(field),nullptr);
    return FittedCellHarmonicAccess::publish(std::move(data));
}
FittedCellHarmonicState updateFittedCellHarmonic(const QuadtreeMesh& mesh,const FittedCellHarmonicState& previous,
    std::vector<double> epsilon,const FittedCellHarmonicControls& controls,double time,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source) {
    previous.validateMesh(mesh);
    require(std::isfinite(time) && time>=previous.time(),"electric material time reversed or nonfinite");
    auto data=FittedCellHarmonicAccess::create();data->problem=FittedCellHarmonicAccess::read(previous).problem;
    validateSource(mesh,*data->problem,source);
    data->epsilonCell=std::move(epsilon);data->epsilonSegment=coefficients(mesh,*data->problem,data->epsilonCell);
    const FittedPotentialInput seed(previous.potential());auto field=FittedCellHarmonicAccess::working(mesh,&seed);
    data->report=solve(mesh,*data->problem,*field,data->epsilonSegment,controls);
    data->potential=FittedCellHarmonicAccess::freeze(std::move(field),std::move(source));
    data->time=time;data->stage=FittedCellHarmonicState::Stage::Solved;
    return FittedCellHarmonicAccess::publish(std::move(data));
}
}
