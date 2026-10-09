#pragma once

#include "Solvers/Poisson/PoissonBC.h"

// Owns the boundary values actually used by one field solve or implicit run.
// Construction preserves the supplied initial values. Only update() evaluates
// driven Dirichlet data; the driver controls its original physical time.
class PicFieldBoundaryState
{
public:
    PicFieldBoundaryState(const PicDomainBC& domain, const PicEmbeddedBC& embedded,
        const PicEmbeddedBCMap& embeddedById);
    PicFieldBoundaryState(const PicFieldBoundaryState&) = delete;
    PicFieldBoundaryState& operator=(const PicFieldBoundaryState&) = delete;

    void update(double timeSI);
    // References/pointers remain valid across update(), but observe new values.
    const PicDomainBC& domain() const { return domain_; }
    const PicEmbeddedBC& embedded() const { return embedded_; }
    const PicEmbeddedBCMap& embeddedById() const { return embeddedById_; }

    // Read the current implicit solve's boundary value, without evaluating a
    // different time. Disabled entries retain the original default fallback.
    bool embeddedVoltage(int boundaryId, double& valueCode) const;

    // The original waveform and diagnostic phase policy: disabled/constant,
    // zero amplitude or nonpositive frequency yields dc and phase zero;
    // unrecognized enabled types retain the existing sine fallback.
    static double evaluate(const PicBoundarySide::TimeFunction& function, double timeSI,
        double* phaseRadians = nullptr, double* phaseDegrees = nullptr);
    static double voltage(const PicBoundarySide& side, double timeSI);
    static bool embeddedVoltage(const PicEmbeddedBC& fallback, const PicEmbeddedBCMap& byId,
        int boundaryId, double timeSI, double& valueCode);

private:
    PicDomainBC domain_;
    PicEmbeddedBC embedded_;
    PicEmbeddedBCMap embeddedById_;
};
