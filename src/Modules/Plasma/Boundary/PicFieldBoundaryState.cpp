#include "Modules/Plasma/Boundary/PicFieldBoundaryState.h"

#include <cctype>
#include <cmath>

namespace {
const PicEmbeddedBC* electrode(const PicEmbeddedBC& fallback, const PicEmbeddedBCMap& byId, int id)
{
    const PicEmbeddedBC* result = nullptr;
    const auto it = byId.find(id);
    if (it != byId.end() && it->second.enabled)
        result = &it->second;
    else if (fallback.enabled)
        result = &fallback;
    return result && result->type == PicBCType::Dirichlet ? result : nullptr;
}
}

PicFieldBoundaryState::PicFieldBoundaryState(const PicDomainBC& domain, const PicEmbeddedBC& embedded,
    const PicEmbeddedBCMap& embeddedById)
    : domain_(domain), embedded_(embedded), embeddedById_(embeddedById) {}

double PicFieldBoundaryState::evaluate(const PicBoundarySide::TimeFunction& tf, double t_si,
    double* phaseRadians, double* phaseDegrees)
{
    if (phaseRadians) *phaseRadians = 0.0;
    if (phaseDegrees) *phaseDegrees = 0.0;
    if (!tf.enabled)
        return tf.dc;
    std::string type = tf.type;
    for (char& c : type)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (type == "constant" || type == "const")
        return tf.dc;
    const double ac = tf.ac;
    const double f = tf.freq_hz;
    if (!(ac != 0.0) || !(f > 0.0))
        return tf.dc;
    const double t = t_si - tf.t0_s;
    constexpr double kPi = 3.141592653589793238462643383279502884;
    const double arg = 2.0 * kPi * f * t + tf.phase_rad;
    if (phaseRadians || phaseDegrees)
    {
        double phase = std::fmod(arg, 2.0 * kPi);
        if (phase < 0.0)
            phase += 2.0 * kPi;
        if (phaseRadians) *phaseRadians = phase;
        if (phaseDegrees) *phaseDegrees = phase * 180.0 / kPi;
    }
    if (type == "cosine" || type == "cos")
        return tf.dc + ac * std::cos(arg);
    return tf.dc + ac * std::sin(arg);
}

double PicFieldBoundaryState::voltage(const PicBoundarySide& side, double timeSI)
{
    if (side.type != PicBCType::Dirichlet) return 0.0;
    return side.timeDirichlet.enabled ? evaluate(side.timeDirichlet, timeSI) : side.value;
}

void PicFieldBoundaryState::update(double timeSI)
{
    for (auto* side : {&domain_.left, &domain_.right, &domain_.bottom, &domain_.top})
        if (side->type == PicBCType::Dirichlet && side->timeDirichlet.enabled)
            side->value = evaluate(side->timeDirichlet, timeSI);
    const auto updateEmbedded = [&](PicEmbeddedBC& boundary) {
        if (boundary.enabled && boundary.type == PicBCType::Dirichlet && boundary.timeDirichlet.enabled)
            boundary.g = evaluate(boundary.timeDirichlet, timeSI);
    };
    updateEmbedded(embedded_);
    for (auto& entry : embeddedById_) updateEmbedded(entry.second);
}

bool PicFieldBoundaryState::embeddedVoltage(int boundaryId, double& valueCode) const
{
    const auto* boundary = electrode(embedded_, embeddedById_, boundaryId);
    if (!boundary) return false;
    valueCode = boundary->g;
    return true;
}

bool PicFieldBoundaryState::embeddedVoltage(const PicEmbeddedBC& fallback, const PicEmbeddedBCMap& byId,
    int boundaryId, double timeSI, double& valueCode)
{
    const auto* boundary = electrode(fallback, byId, boundaryId);
    if (!boundary) return false;
    valueCode = boundary->timeDirichlet.enabled ? evaluate(boundary->timeDirichlet, timeSI) : boundary->g;
    return true;
}
