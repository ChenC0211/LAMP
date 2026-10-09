#pragma once

#include <string>

#include "Modules/Plasma/PlasmaRunInputs.h"
#include "Modules/Plasma/PicCaseInput.h"
#include "Modules/Plasma/PicStepSession.h"

// ----------------------------------------------------------------------------
// ExplicitPICModule
//
// Explicit PIC initialization and macro-step engine. The caller drives the
// step loop; each step retains its particle, field, MCC, AMR and output order.
// ----------------------------------------------------------------------------

class ExplicitPICModule
{
public:
    explicit ExplicitPICModule(const PlasmaRunInputs& inputs);

    std::string name() const { return "ExplicitPIC"; }

    static constexpr PicIntegrator integrator = PicIntegrator::Explicit;
    // Initialize once, lend a synchronous step session to the caller, then
    // perform the original final output only after successful completion.
    int withSession(PicCaseInput input, const PicSessionDriver& driver);

private:
    PlasmaRunInputs cfg_;
};
