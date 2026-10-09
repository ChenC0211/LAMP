#include "Modules/Hydrodynamics/HydroFlowDiagnostics.h"
#include "Core/Field/FieldChecks.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Modules/Hydrodynamics/HydroConservativePressure.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace HydroFlowDiagnostics
{
using FieldChecks::checkFinite;

NativeConservationDiagnostics nativeConservation(const QuadtreeMesh& mesh,const Hydrodynamics& hydro,const NativeLeafDomain& domain)
{
    domain.validateMeshIdentity(mesh);const auto state=hydro.state();
    const auto duals=HydroFittedDetail::MomentumDuals::build(mesh,domain);
    long double total=0,energy=0;
    for(int c:domain.cellIndexByRow())total+=state.density[c]*static_cast<long double>(domain.cellVolumes()[domain.cellRowByIndex()[c]]);
    for(const auto face:domain.faces(mesh)){
        const auto& dual=duals.volumes()[face.nativeRow];long double mass=0;
        for(int c:{face.negativeCell,face.positiveCell})if(c>=0){
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
            const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
            const double x0=std::max(a.x,dual.lower.x),x1=std::min(b.x,dual.upper.x),
                y0=std::max(a.y,dual.lower.y),y1=std::min(b.y,dual.upper.y);
            mass+=state.density[c]*static_cast<long double>(CartesianMeasure::cell((x1-x0)*(y1-y0),.5*(x0+x1),domain.axisymmetric()));
        }
        const double u=face.horizontal?state.velocityY[face.nativeRow-mesh.vEdgeIndexNum()]:state.velocityX[face.nativeRow];
        energy+=.5L*mass*u*u;
    }
    if(!std::isfinite(total) || !std::isfinite(energy) || total<=0 || energy<0)
        throw std::runtime_error("nonfinite native conservation diagnostic");
    return {static_cast<double>(total),static_cast<double>(energy)};
}

template<class Field>
void checkFinite(const HydroReadOnlyField<Field>& field, const std::string& name, int count)
{
    if (field.size() != count || !field.hasStorage())
        throw std::runtime_error(name + " storage does not match the fixed mesh");
    for (int i = 0; i < field.size(); ++i)
        if (!std::isfinite(field[i]))
            throw std::runtime_error(name + " contains a non-finite or invalid value at index " + std::to_string(i));
}

double axisVelocity(const QuadtreeMesh& mesh, const Hydrodynamics& hydro, bool rz)
{
    const auto state = hydro.state();
    double maximum = 0.0;
    if (rz)
        for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
            if (mesh.vEdges(mesh.vEdgeIndexToLabel[i]).midpointCoord().x == 0.0)
                maximum = std::max(maximum, std::abs(state.velocityX[i]));
    return maximum;
}

void collectForcePressureImbalance(const QuadtreeMesh& mesh, const Hydrodynamics& hydro,
    const VEdgeField<double>& faceForceX, const HEdgeField<double>& faceForceY,
    bool rz, HydroInterfaceForce::Report& diagnostic)
{
    const auto state = hydro.state();
    checkFinite(faceForceX, "active total face force x", mesh.vEdgeIndexNum());
    checkFinite(faceForceY, "active total face force y", mesh.hEdgeIndexNum());
    checkFinite(state.externalForceX, "active total cell force x", mesh.cellIndexNum());
    checkFinite(state.externalForceY, "active total cell force y", mesh.cellIndexNum());
    checkFinite(state.physicalPressure, "post-solve physical pressure", mesh.cellIndexNum());
    long double squared = 0.0L, volume = 0.0L;
    auto collect = [&](const Edge& edge, bool horizontal, double faceForce)
    {
        if (edge.masterCell == -1 || edge.slaveCell == -1) return;
        const auto& master = mesh.cells(edge.masterCell);
        const auto& slave = mesh.cells(edge.slaveCell);
        const auto& cellForce = horizontal ? state.externalForceY : state.externalForceX;
        const double gradient = (state.physicalPressure[master.index()] - state.physicalPressure[slave.index()])
            / (horizontal ? master.size().y : master.size().x);
        const double residual = faceForce + 0.5 * (cellForce[master.index()] + cellForce[slave.index()]) - gradient;
        const double measure = master.area() * (rz ? 2.0 * Pi * edge.midpointCoord().x : 1.0);
        if (!std::isfinite(residual) || !std::isfinite(measure) || measure <= 0.0)
            throw std::runtime_error("non-finite force-pressure imbalance or invalid staggered dual volume");
        diagnostic.maxForcePressureImbalance = std::max(diagnostic.maxForcePressureImbalance, std::abs(residual));
        squared += static_cast<long double>(residual) * residual * measure;
        volume += measure;
    };
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
        collect(mesh.vEdges(mesh.vEdgeIndexToLabel[i]), false, faceForceX.value[i]);
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
        collect(mesh.hEdges(mesh.hEdgeIndexToLabel[i]), true, faceForceY.value[i]);
    if (!std::isfinite(squared) || !std::isfinite(volume) || volume <= 0.0L)
        throw std::runtime_error("force-pressure imbalance norm overflow or empty face domain");
    diagnostic.l2ForcePressureImbalance = static_cast<double>(std::sqrt(squared / volume));
    if (!std::isfinite(diagnostic.l2ForcePressureImbalance))
        throw std::runtime_error("non-finite force-pressure imbalance L2 norm");
}

void collectForcePressureImbalance(const QuadtreeMesh& mesh,const Hydrodynamics& hydro,
    const VEdgeField<double>& faceForceX,const HEdgeField<double>& faceForceY,const NativeLeafDomain& domain,
    const Quantum::Hydrodynamics::HydroConservativePressure& pressure,HydroInterfaceForce::Report& diagnostic)
{
    const auto state=hydro.state();pressure.validate(mesh,domain);
    checkFinite(faceForceX,"active native face force x",mesh.vEdgeIndexNum());
    checkFinite(faceForceY,"active native face force y",mesh.hEdgeIndexNum());
    checkFinite(state.externalForceX,"active native cell force x",mesh.cellIndexNum());
    checkFinite(state.externalForceY,"active native cell force y",mesh.cellIndexNum());
    checkFinite(state.physicalPressure,"native physical pressure",mesh.cellIndexNum());
    std::vector<double> values(mesh.cellIndexNum());
    for(int c=0;c<mesh.cellIndexNum();++c)values[c]=state.physicalPressure[c];
    const auto load=pressure.apply(mesh,domain,values,std::vector<double>(domain.nativeFaceCount(),0));
    long double squared=0,volume=0;
    for(const auto face:domain.faces(mesh)) {
        if(face.outwardSign)continue;
        const auto& dual=pressure.duals().volumes()[face.nativeRow];
        long double body=0;
        for(int owner:{face.negativeCell,face.positiveCell}) {
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[owner]);
            const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
            const Vector2D lo{std::max(a.x,dual.lower.x),std::max(a.y,dual.lower.y)},
                hi{std::min(b.x,dual.upper.x),std::min(b.y,dual.upper.y)};
            body+=CartesianMeasure::cell((hi.x-lo.x)*(hi.y-lo.y),(hi.x+lo.x)/2,domain.axisymmetric())
                *(face.horizontal?state.externalForceY[owner]:state.externalForceX[owner]);
        }
        const double force=face.horizontal?faceForceY.value[face.nativeRow-mesh.vEdgeIndexNum()]:faceForceX.value[face.nativeRow];
        const double residual=force+static_cast<double>((body-load[face.nativeRow])/dual.measure);
        if(!std::isfinite(residual) || !std::isfinite(dual.measure) || dual.measure<=0)
            throw std::runtime_error("invalid native force-pressure imbalance");
        diagnostic.maxForcePressureImbalance=std::max(diagnostic.maxForcePressureImbalance,std::abs(residual));
        squared+=static_cast<long double>(residual)*residual*dual.measure;volume+=dual.measure;
    }
    if(!std::isfinite(squared) || !std::isfinite(volume) || volume<=0)
        throw std::runtime_error("invalid native force-pressure norm");
    diagnostic.l2ForcePressureImbalance=static_cast<double>(std::sqrt(squared/volume));
}

GeometryMaterialDiagnostics collectDiagnostics(const QuadtreeMesh& mesh, const Hydrodynamics& hydro,
    const ImmersedGeometry<QuadtreeMesh>* geometry, int interfaceIndex,
    const std::vector<double>& initialLevelSet, bool rz, double interfaceThickness,
    int boundaryCutCells)
{
    const auto state = hydro.state();
    GeometryMaterialDiagnostics result;
    result.rhoMin = result.rhoMax = state.density[0];
    result.muMin = result.muMax = state.viscosity[0];
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        result.rhoMin = std::min(result.rhoMin, state.density[i]);
        result.rhoMax = std::max(result.rhoMax, state.density[i]);
        result.muMin = std::min(result.muMin, state.viscosity[i]);
        result.muMax = std::max(result.muMax, state.viscosity[i]);
    }
    if (!geometry) return result;
    const auto& levelSet = geometry->levelSets[interfaceIndex];
    result.levelSetMin = result.levelSetMax = levelSet.value[0];
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        result.interfaceMaxChange = std::max(result.interfaceMaxChange, std::abs(levelSet.value[i] - initialLevelSet[i]));
        result.levelSetMin = std::min(result.levelSetMin, levelSet.value[i]);
        result.levelSetMax = std::max(result.levelSetMax, levelSet.value[i]);
    }
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        double fraction = 0.0;
        for (int corner = 0; corner < 4; ++corner)
        {
            const int vi = mesh.verts(cell.cellToVert[corner]).index();
            fraction += 0.125 * (1.0 - std::tanh(levelSet.value[vi] / (2.0 * interfaceThickness)));
        }
        const double volume = cell.area() * (rz ? 2.0 * std::acos(-1.0) * cell.centerCoord().x : 1.0);
        result.negativeVolume += fraction * volume;
    }
    if (!std::isfinite(result.negativeVolume) || !std::isfinite(result.interfaceMaxChange))
        throw std::runtime_error("geometry/material volume diagnostics are non-finite");
    result.cutCells = boundaryCutCells >= 0 ? boundaryCutCells : geometry->cutInfoDiagnostics().cutCellCount;
    return result;
}

BoundaryFlowDiagnostics boundaryFlow(const QuadtreeMesh& mesh, const Hydrodynamics& hydro,
    bool axisymmetric, double tolerance, bool constantMaterial)
{
    const auto state = hydro.state();
    long double inflow = 0.0L, outflow = 0.0L, net = 0.0L, massNet = 0.0L, massScale = 0.0L;
    const auto add = [&](const Edge& edge, double velocity) {
        if (edge.masterCell != -1 && edge.slaveCell != -1) return;
        const int owner = edge.masterCell == -1 ? edge.slaveCell : edge.masterCell;
        if (owner < 0) throw std::runtime_error("boundary flow face has no fluid owner");
        const double density = state.density[mesh.cells(owner).index()];
        const double area = edge.length() * (axisymmetric ? 2.0 * Pi * edge.midpointCoord().x : 1.0);
        const long double flux = static_cast<long double>(area) * velocity * (edge.masterCell == -1 ? 1.0L : -1.0L);
        if (!std::isfinite(flux) || !std::isfinite(density) || density <= 0.0 || !std::isfinite(area) || area < 0.0)
            throw std::runtime_error("boundary flow has invalid velocity, density or physical measure");
        if (flux < 0.0L) inflow -= flux;
        else outflow += flux;
        net += flux;
        massNet += density * flux;
        massScale += density * std::abs(flux);
    };
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
        add(mesh.vEdges(mesh.vEdgeIndexToLabel[i]), state.velocityX[i]);
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
        add(mesh.hEdges(mesh.hEdgeIndexToLabel[i]), state.velocityY[i]);
    if (!std::isfinite(inflow) || !std::isfinite(outflow) || !std::isfinite(net) || !std::isfinite(massNet)
        || std::abs(net) > 16.0L*tolerance*(inflow+outflow)
        || (constantMaterial && std::abs(massNet) > 16.0L*tolerance*massScale))
        throw std::runtime_error("accepted open-boundary flow does not close volume flux or constant-material mass flux");
    return {static_cast<double>(inflow), static_cast<double>(outflow),
        static_cast<double>(net), static_cast<double>(massNet)};
}
}
