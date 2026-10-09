#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include <array>
#include <stdexcept>
#include <type_traits>
#include <utility>

struct HydroStepCandidate::Control
{
    enum class Stage { Rejected, Pending, Accepted, Abandoned };
    HydroSolveReport report;
    Hydrodynamics* owner = nullptr;
    const QuadtreeMesh* mesh = nullptr;
    Stage stage = Stage::Rejected;
};

HydroStepCandidate::HydroStepCandidate(std::shared_ptr<Control> control) noexcept
    : control_(std::move(control)) {}
HydroStepCandidate::~HydroStepCandidate() { abandon(); }
HydroStepCandidate::HydroStepCandidate(HydroStepCandidate&&) noexcept = default;
HydroStepCandidate& HydroStepCandidate::operator=(HydroStepCandidate&& other) noexcept
{
    if (this != &other)
    {
        abandon();
        control_ = std::move(other.control_);
    }
    return *this;
}

const HydroSolveReport& HydroStepCandidate::report() const
{
    if (!control_) throw std::logic_error("moved-from Hydro candidate has no report");
    return control_->report;
}
bool HydroStepCandidate::pending() const noexcept
{ return control_ && control_->stage == Control::Stage::Pending; }
bool HydroStepCandidate::accepted() const noexcept
{ return control_ && control_->stage == Control::Stage::Accepted; }

void HydroStepCandidate::accept(const QuadtreeMesh& mesh)
{
    if (!pending() || !control_->owner)
        throw std::logic_error("Hydro candidate is not pending or its owner no longer exists");
    auto& owner = *control_->owner;
    if (owner.pendingStep_ != control_ || !owner.stepProgress_.pendingAcceptance || owner.stepProgress_.stopped)
        throw std::logic_error("Hydro candidate belongs to an invalidated substep");
    if (control_->mesh != &mesh)
        throw std::invalid_argument("Hydro candidate acceptance requires its original mesh");
    owner.validateStateGeometry(mesh, {*owner.stepProgress_.geometry, owner.fittedDomain_.get()});
    // All potentially throwing checks precede these scalar/ownership changes.
    owner.stepProgress_.acceptedSteps = owner.stepProgress_.completedSteps;
    owner.stepProgress_.acceptedEndTime = owner.stepProgress_.endTime;
    owner.stepProgress_.pendingAcceptance = false;
    control_->stage = Control::Stage::Accepted;
    control_->owner = nullptr;
    owner.pendingStep_.reset();
}

void HydroStepCandidate::abandon() noexcept
{
    if (pending() && control_->owner) control_->owner->discardPendingStep();
}

void Hydrodynamics::discardPendingStep() noexcept
{
    if (!pendingStep_) return;
    pendingStep_->owner = nullptr;
    pendingStep_->stage = HydroStepCandidate::Control::Stage::Abandoned;
    stepProgress_.pendingAcceptance = false;
    stepProgress_.stopped = true;
    pendingStep_.reset();
}

HydroStepCandidate Hydrodynamics::prepareStep(QuadtreeMesh& mesh, const HydroStepInput& input)
{
    // Allocate only the small lifetime/report record before the numerical
    // entry; no new copy of mesh, fields, histories or factorization caches.
    auto control = std::make_shared<HydroStepCandidate::Control>();
    static_assert(std::is_nothrow_move_assignable_v<HydroSolveReport>);
    control->report = advanceStep(mesh, input, true);
    if (control->report.success)
    {
        control->owner = this;
        control->mesh = &mesh;
        control->stage = HydroStepCandidate::Control::Stage::Pending;
        pendingStep_ = control;
    }
    return HydroStepCandidate(std::move(control));
}

struct HydroStepTransaction::Control
{
    Hydrodynamics* owner=nullptr;
    const QuadtreeMesh* mesh=nullptr;
    std::shared_ptr<const NativeMeshIdentity> identity;
    std::array<Field<double>,26> values;
    HydroStepProgress progress;
    Hydrodynamics::StepEntry entry;
    bool initialized,conservative;
    std::optional<HydroMaterialState> materials;
    std::shared_ptr<const NativeMeshIdentity> stepIdentity;
    std::shared_ptr<const HydroFittedDomain> fitted;
    std::shared_ptr<const HydroFittedMomentumLayout> layout;
    Hydrodynamics::PressureMode pressureMode;
    std::vector<Hydrodynamics::BoundaryFaceData> boundaryFaces;
    std::vector<HydroBoundaryPatch> boundaryPatches;
    bool boundaryFixed;
    double maxDivergence,maxSpeed,iterationCoefficient,fixedDt,thickness;
    HydroSolveOptions options;
    HydroSolveReport report;

    static std::array<Field<double>*,26> fields(Hydrodynamics& h) {
        return {&h.vertRho,&h.vertMu,&h.p,&h.rho,&h.mu,&h.externalFx,&h.externalFy,
            &h.directFx,&h.directFy,&h.evapUx,&h.evapUy,&h.ux,&h.uy,&h.uxPre,&h.uyPre,
            &h.uxMid,&h.uyMid,&h.uMag,&h.divU,&h.dxUx,&h.dyUx,&h.dxUy,&h.dyUy,
            &h.magGradU,&h.hyperTanLS,&h.pressurePotential_};
    }
    Control(Hydrodynamics& h,const QuadtreeMesh& m)
        :owner(&h),mesh(&m),identity(h.fittedDomain_?h.fittedDomain_->geometry().meshIdentity():h.stepMeshIdentity_),
        progress(h.stepProgress_),entry(h.stepEntry_),initialized(h.initialStateSet_),conservative(h.conservativeTransportActive_),
        materials(h.materialState_),stepIdentity(h.stepMeshIdentity_),fitted(h.fittedDomain_),layout(h.fittedMomentumLayout_),
        pressureMode(h.pressureMode_),boundaryFaces(h.boundaryFaces_),boundaryPatches(h.frozenBoundaryPatches_),
        boundaryFixed(h.boundaryConfigurationFixed_),maxDivergence(h.maxDivURes_),maxSpeed(h.uMagMax_),
        iterationCoefficient(h.iterCoeff_),fixedDt(h.fixedDt_),thickness(h.hyperTanLS.thickness()),
        options(h.activeOptions_),report(h.activeReport_)
    {
        if(!identity) throw std::logic_error("Hydro rollback requires an initialized mesh identity");
        identity->validate(m);
        const auto source=fields(h);
        for(std::size_t i=0;i<source.size();++i) if(source[i]->size()) values[i]=*source[i];
    }
    void restore() noexcept {
        auto& h=*owner;
        h.discardPendingStep();
        const auto target=fields(h);
        for(std::size_t i=0;i<target.size();++i) target[i]->swap(values[i]);
        h.stepProgress_=progress;h.stepEntry_=entry;h.initialStateSet_=initialized;h.conservativeTransportActive_=conservative;
        h.materialState_.swap(materials);h.stepMeshIdentity_.swap(stepIdentity);h.fittedDomain_.swap(fitted);h.fittedMomentumLayout_.swap(layout);
        h.pressureMode_=pressureMode;h.boundaryFaces_.swap(boundaryFaces);h.frozenBoundaryPatches_.swap(boundaryPatches);
        h.boundaryConfigurationFixed_=boundaryFixed;h.maxDivURes_=maxDivergence;h.uMagMax_=maxSpeed;
        h.iterCoeff_=iterationCoefficient;h.fixedDt_=fixedDt;h.hyperTanLS.setThickness(thickness);
        static_assert(std::is_nothrow_move_assignable_v<HydroSolveOptions> && std::is_nothrow_move_assignable_v<HydroSolveReport>);
        h.activeOptions_=std::move(options);h.activeReport_=std::move(report);
        // Factors may have been mutated in place. Reassemble them on retry;
        // never restore a pointer to a factor from the rejected material layer.
        h.momentumLinearCache_.reset();h.pressureLinearCache_.reset();h.coupledLinearCache_.reset();
        // Geometric interpolation scratch is rebuilt by the existing kernel.
        h.supCells_.reset();h.supSonCells_.reset();
    }
};

HydroStepTransaction Hydrodynamics::beginStepTransaction(const QuadtreeMesh& mesh)
{
    if(stepTransaction_ || pendingStep_ || !initialStateSet_ || !stepProgress_.geometry ||
        stepProgress_.stopped || stepProgress_.pendingAcceptance || stepEntry_==StepEntry::Legacy ||
        stepProgress_.completedSteps!=stepProgress_.acceptedSteps || stepProgress_.endTime!=stepProgress_.acceptedEndTime)
        throw std::logic_error("Hydro rollback requires a ready accepted managed state and no open transaction");
    validateStateGeometry(mesh,{*stepProgress_.geometry,fittedDomain_.get()});
    auto control=std::make_shared<HydroStepTransaction::Control>(*this,mesh);
    stepTransaction_=control;
    return HydroStepTransaction(std::move(control));
}
void Hydrodynamics::validateStepTransactionInput() const
{
    if(stepTransaction_ && stepProgress_.completedSteps!=stepTransaction_->progress.completedSteps)
        throw std::logic_error("commit or roll back the current Hydro transaction before another step");
}
void Hydrodynamics::discardStepTransaction() noexcept
{
    if(stepTransaction_) {stepTransaction_->owner=nullptr;stepTransaction_.reset();}
}
HydroStepTransaction::HydroStepTransaction(std::shared_ptr<Control> control) noexcept :control_(std::move(control)) {}
HydroStepTransaction::~HydroStepTransaction() {abandon();}
HydroStepTransaction::HydroStepTransaction(HydroStepTransaction&&) noexcept = default;
HydroStepTransaction& HydroStepTransaction::operator=(HydroStepTransaction&& other) noexcept
{
    if(this!=&other){abandon();control_=std::move(other.control_);}return *this;
}
bool HydroStepTransaction::pending() const noexcept {return control_ && control_->owner;}
void HydroStepTransaction::rollback()
{
    if(!pending()) throw std::logic_error("Hydro transaction has no live pending owner");
    control_->identity->validate(*control_->mesh);
    auto& h=*control_->owner;
    if(h.stepTransaction_!=control_) throw std::logic_error("Hydro transaction was replaced");
    control_->restore();h.discardStepTransaction();control_.reset();
}
void HydroStepTransaction::commit()
{
    if(!pending()) throw std::logic_error("Hydro transaction has no live pending owner");
    auto& h=*control_->owner;
    control_->identity->validate(*control_->mesh);
    if(h.stepTransaction_!=control_ || h.pendingStep_ || h.stepProgress_.stopped || h.stepProgress_.pendingAcceptance ||
        h.stepProgress_.acceptedSteps!=control_->progress.acceptedSteps+1 ||
        h.stepProgress_.completedSteps!=h.stepProgress_.acceptedSteps)
        throw std::logic_error("Hydro transaction commit requires one accepted physical step");
    h.discardStepTransaction();control_.reset();
}
void HydroStepTransaction::abandon() noexcept
{
    if(!pending()) return;
    try {rollback();}
    catch(...) {
        // An externally changed mesh cannot receive the old arrays. Retain a
        // stopped state rather than silently pretending restoration succeeded.
        auto& h=*control_->owner;h.discardPendingStep();h.stepProgress_.stopped=true;
        h.discardStepTransaction();control_.reset();
    }
}
