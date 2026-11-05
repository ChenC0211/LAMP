#pragma once

#include "Modules/Hydrodynamics/HydroTwoPhaseStep.h"
#include "Modules/Electric/ElectricPotential.h"
#include "Modules/Electric/ElectricPhaseGeometry.h"
#include "Modules/Electric/ElectricMaxwellTrace.h"
#include "Modules/Electric/SurfaceChargeField.h"
#include <array>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace LAMP::Coupling
{
// The Electrohydrodynamics App owns one synchronous accepted transaction. The interface
// publication and master clock belong to hydro_; no second geometry is stored.
class ElectrohydrodynamicState
{
public:
    static ElectrohydrodynamicState create(const QuadtreeMesh& mesh,
        Hydrodynamics::HydroTwoPhaseState hydro,Electric::ElectricSurfaceChargeTimeState&& electric)
    {
        validate(mesh,hydro,electric);
        return ElectrohydrodynamicState(std::move(hydro),std::move(electric));
    }
    ElectrohydrodynamicState(ElectrohydrodynamicState&&) noexcept=default;
    ElectrohydrodynamicState& operator=(ElectrohydrodynamicState&&) noexcept=default;
    ElectrohydrodynamicState(const ElectrohydrodynamicState&)=delete;
    ElectrohydrodynamicState& operator=(const ElectrohydrodynamicState&)=delete;
    const Hydrodynamics::HydroTwoPhaseState& hydro() const noexcept { return hydro_; }
    const Electric::ElectricSurfaceChargeTimeState& electric() const noexcept { return electric_; }
    const Interface::HermiteGalsPublication& publication() const noexcept { return hydro_.publication(); }
    Interface::HermiteGalsTimeLayer timeLayer() const { return hydro_.timeLayer(); }
    void validateFor(const QuadtreeMesh& mesh) const { validate(mesh,hydro_,electric_); }
private:
    static bool sameClock(Interface::HermiteGalsTimeLayer a,Interface::HermiteGalsTimeLayer b)
    { return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position; }
    static void require(bool ok,const char* why)
    { if(!ok)throw std::invalid_argument(std::string("[ElectrohydrodynamicState] ")+why); }
    static void validate(const QuadtreeMesh& mesh,const Hydrodynamics::HydroTwoPhaseState& hydro,
        const Electric::ElectricSurfaceChargeTimeState& electric)
    {
        hydro.validateFor(mesh);
        const auto& geometry=electric.solver().momentOperator().geometry();
        geometry.validateSource(mesh,hydro.publication());
        require(geometry.axisymmetric(),"joint moving charge currently requires native open RZ");
        require(sameClock(hydro.timeLayer(),electric.timeLayer())
            &&sameClock(hydro.timeLayer(),electric.geometrySourceTimeLayer()),
            "accepted electric/Hydro physical or geometry full clock differs");
        const auto& curves=hydro.publication().topology().curves;
        require(electric.componentCount()==curves.size(),"accepted component coverage differs");
        for(std::size_t c=0;c<curves.size();++c)
            require(!curves[c].closed()&&!electric.field(c).geometry().closed(),
                "joint material transfer currently requires open native components");
    }
    ElectrohydrodynamicState(Hydrodynamics::HydroTwoPhaseState hydro,
        Electric::ElectricSurfaceChargeTimeState electric)
        :hydro_(std::move(hydro)),electric_(std::move(electric)) {}
    Hydrodynamics::HydroTwoPhaseState hydro_;
    Electric::ElectricSurfaceChargeTimeState electric_;
    friend class ElectrohydrodynamicStep;
};

struct ElectrohydrodynamicStepConfig
{
    Hydrodynamics::HydroTwoPhaseStepConfig hydro;
    Electric::ElectricPhaseProblem electricProblem;
    std::array<double,2> conductivityByPhase{{0,0}};
    Electric::ElectricPhaseGeometryConfig geometry;
    Electric::ElectricSurfaceChargeMomentControls moments;
    Electric::ElectricCartesianControls electric;
    // A normwise Schur stop needs margin for the independent panelwise gate
    // when signed charge makes individual panel scales small.
    Electric::ElectricChargeTimeControls charge=[] {
        Electric::ElectricChargeTimeControls controls;
        controls.relativeSchurTolerance=1e-13;
        return controls;
    }();
    // Passive caller-thread observer, outside the accepted publication boundary.
    // A failure rejects preparation without editing any accepted owner.
    std::function<void(const Electric::ElectricMovingGeometryChargeStep&)> onElectricCandidate;
};

// Wires existing kernels; it contains no interface transport, charge equation,
// Maxwell tensor or force discretization. All stages use their actual Hydro
// transport/source. End is prepared from the original accepted electric history.
class ElectrohydrodynamicStep
{
public:
    static ElectrohydrodynamicStep prepare(const QuadtreeMesh& mesh,
        const ElectrohydrodynamicState& accepted,Interface::HermiteGalsTimeLayer end,
        const ElectrohydrodynamicStepConfig& config)
    {
        accepted.validateFor(mesh);
        require(!config.hydro.interfaceLoad,"joint config cannot replace its authenticated stage wiring");
        require(config.hydro.timeScheme==Hydrodynamics::HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint
            ||config.hydro.timeScheme==Hydrodynamics::HydroTwoPhaseTimeScheme::ExtrapolatedMidpoint,
            "joint electric/Hydro step requires the AC time route");
        require(config.hydro.capillaryModel==Hydrodynamics::HydroCapillaryModel::BalancedCsf,
            "Sharp external electric traction is not implemented; select the explicit balancedCsf route");
        require(end.step>=0&&end.step<=std::numeric_limits<int>::max(),
            "native charge transfer step exceeds its supported integer range");
        std::vector<Electric::ElectricMovingGeometryChargeStep> stages;
        auto hydroConfig=config.hydro;
        hydroConfig.interfaceLoad=[&](const QuadtreeMesh& actualMesh,
            const Hydrodynamics::HydroMaterialTransport& transport) {
            const auto& publication=transport.midpointRate().geometry().publication();
            auto candidate=prepareElectric(actualMesh,accepted,transport,publication,config,true);
            auto trace=Electric::ElectricMaxwellTrace::prepare(candidate.step().potential());
            auto view=trace.surfaceVectorView();
            stages.push_back(std::move(candidate));
            return view;
        };
        auto hydro=Hydrodynamics::HydroTwoPhaseStep::prepare(mesh,accepted.hydro(),end,hydroConfig);
        require(!stages.empty()&&stages.size()==static_cast<std::size_t>(hydro.report().iterations),
            "electric candidates do not cover every actual AC stage");
        const auto& load=hydro.capillaryLoad();
        load.externalSurfaceVector().validateFor(mesh,
            hydro.materialTransport().midpointRate().geometry().publication(),
            hydro.materialTransport().timeLayer());
        stages.back().step().potential().geometry().validateSource(mesh,load.phaseGeometry().publication());
        auto electricEnd=prepareElectric(mesh,accepted,hydro.materialTransport(),
            hydro.candidateState().publication(),config,false);
        return ElectrohydrodynamicStep(std::make_shared<const Data>(
            std::move(hydro),std::move(stages),std::move(electricEnd)));
    }
    const Hydrodynamics::HydroTwoPhaseStep& hydroStep() const noexcept { return data_->hydro; }
    const std::vector<Electric::ElectricMovingGeometryChargeStep>& electricStages() const noexcept
    { return data_->stages; }
    const Electric::ElectricMovingGeometryChargeStep& electricEnd() const noexcept { return data_->end; }
    void commit(const QuadtreeMesh& mesh,ElectrohydrodynamicState& accepted) const
    {
        // Complete every Hydro/geometry check on a LOCAL handle first. This
        // publishes nothing to accepted and can freely throw or allocate.
        auto nextHydro=accepted.hydro_;
        data_->hydro.commit(mesh,nextHydro);
        data_->end.step().potential().geometry().validateSource(mesh,nextHydro.publication());
        require(sameClock(data_->end.step().timeLayer(),nextHydro.timeLayer()),
            "End electric/Hydro full clock differs at joint acceptance");
        static_assert(std::is_nothrow_move_assignable_v<Hydrodynamics::HydroTwoPhaseState>,
            "joint publication requires noexcept Hydro handle assignment");
        // Original electric commit completes all checks/allocations before its own
        // nonthrowing field-lease/data swap. No observer or IO occurs here.
        accepted.electric_.commit(data_->end);
        accepted.hydro_=std::move(nextHydro);
    }
private:
    static void require(bool ok,const char* why)
    { if(!ok)throw std::invalid_argument(std::string("[ElectrohydrodynamicStep] ")+why); }
    static bool sameClock(Interface::HermiteGalsTimeLayer a,Interface::HermiteGalsTimeLayer b)
    { return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position; }
    static Electric::ElectricMovingGeometryChargeStep prepareElectric(const QuadtreeMesh& mesh,
        const ElectrohydrodynamicState& accepted,const Hydrodynamics::HydroMaterialTransport& transport,
        const Interface::HermiteGalsPublication& publication,const ElectrohydrodynamicStepConfig& config,
        bool midpoint)
    {
        const auto clock=publication.state().timeLayer();
        const auto& source=transport.characteristicSource();
        source.validateFor(mesh,source.previous().identity(),source.arrival().identity());
        require(source.previous().geometryIdentity()==accepted.hydro().geometryIdentity()
            &&sameClock(source.previous().timeLayer(),accepted.timeLayer())
            &&sameClock(source.arrival().timeLayer(),transport.endGeometry().timeLayer()),
            "material trajectory lost the actual accepted and endpoint owners");
        require(clock.step==source.arrival().timeLayer().step
            &&clock.begin==source.arrival().timeLayer().begin&&clock.end==source.arrival().timeLayer().end,
            "target electric publication differs from the actual complete Hydro interval");
        const auto& curves=publication.topology().curves;
        require(curves.size()==accepted.electric().componentCount(),"material transfer changes component topology");
        std::vector<Electric::SurfaceChargeMaterialTransfer> transfers;
        transfers.reserve(curves.size());
        for(std::size_t c=0;c<curves.size();++c){
            const auto panels=LevelSetInterfacePanels(curves[c]);
            transfers.push_back(accepted.electric().field(c).prepareMaterialTransfer(mesh,panels,
                [&](Vector2D p){
                    const auto mapped=source.traceForward(p,accepted.timeLayer().time(),clock.time(),
                        config.hydro.spaceTime.characteristic);
                    require(!mapped.boundary()&&mapped.time()==clock.time(),
                        "material trajectory hit a physical boundary before its actual target time");
                    return mapped.point();
                },static_cast<int>(clock.step),accepted.timeLayer().time(),clock.time()));
        }
        auto geometry=Electric::ElectricPhaseGeometry::build(mesh,publication,config.geometry);
        auto candidate=midpoint
            ?accepted.electric().prepareMovingStage(mesh,geometry,transfers,clock,config.electricProblem,
                config.conductivityByPhase,config.moments,config.electric,config.charge)
            :accepted.electric().prepareMoving(mesh,geometry,transfers,clock,config.electricProblem,
                config.conductivityByPhase,config.moments,config.electric,config.charge);
        if(config.onElectricCandidate)config.onElectricCandidate(candidate);
        return candidate;
    }
    struct Data
    {
        Data(Hydrodynamics::HydroTwoPhaseStep hydro,
            std::vector<Electric::ElectricMovingGeometryChargeStep> stages,
            Electric::ElectricMovingGeometryChargeStep end)
            :hydro(std::move(hydro)),stages(std::move(stages)),end(std::move(end)) {}
        Hydrodynamics::HydroTwoPhaseStep hydro;
        std::vector<Electric::ElectricMovingGeometryChargeStep> stages;
        Electric::ElectricMovingGeometryChargeStep end;
    };
    explicit ElectrohydrodynamicStep(std::shared_ptr<const Data> data):data_(std::move(data)) {}
    std::shared_ptr<const Data> data_;
};
}
