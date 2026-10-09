#pragma once

#include <string>

#include "Modules/Plasma/PlasmaRunInputs.h"
#include "Modules/Plasma/PicCaseInput.h"
#include "Modules/Plasma/PicStepSession.h"

// ----------------------------------------------------------------------------
// ImplicitPICModule
//
// Semi-implicit FV--PIC macro-step with polarization closure (PaperB).
// This module is designed to reuse ExplicitPIC sub-systems (deposit/gather/push,
// MCC, AMR infrastructure) with its own macro-step physics. The caller controls
// repetition through the same scoped step session as ExplicitPIC.
// ----------------------------------------------------------------------------

class ImplicitPICModule
{
public:
    explicit ImplicitPICModule(const PlasmaRunInputs& inputs);

    std::string name() const { return "ImplicitPIC"; }

    static constexpr PicIntegrator integrator = PicIntegrator::Implicit;
    // Initialize once, lend a synchronous step session to the caller, then
    // perform the original final output only after successful completion.
    int withSession(PicCaseInput input, const PicSessionDriver& driver);

private:
    PlasmaRunInputs cfg_;
};
