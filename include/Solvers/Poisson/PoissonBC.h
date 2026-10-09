#pragma once

#include <string>
#include <unordered_map>

// -----------------------------------------------------------------------------
// Poisson boundary-condition PODs for FV/EB solvers.
//
// - Keep this header lightweight (no Core/Mesh heavy includes).
// - Provide compatibility aliases (Pic* names) for existing ExplicitPIC code.
// -----------------------------------------------------------------------------

namespace Solvers
{
namespace Poisson
{

enum class BCType
{
    Dirichlet,
    Neumann,
    Robin,
    Periodic
};

struct BoundarySide
{
    BCType type = BCType::Dirichlet;
    double value = 0.0; // for Dirichlet/Neumann: phi or dphi/dn

    // Optional time-varying Dirichlet value:
    //   phi(t) = dc + ac * sin(2*pi*freq*(t - t0) + phase)
    // Only applied when type == Dirichlet.
    struct TimeFunction
    {
        bool enabled = false;
        std::string type = "sine"; // sine | cosine | constant
        double dc = 0.0;
        double ac = 0.0;
        double freq_hz = 0.0;
        double phase_rad = 0.0;
        double t0_s = 0.0;
    };
    TimeFunction timeDirichlet;
};

struct DomainBC
{
    BoundarySide left;
    BoundarySide right;
    BoundarySide bottom;
    BoundarySide top;
};

// Generic Robin form on embedded boundary:
//   a*phi + b*dphi/dn = g
struct EmbeddedBC
{
    bool enabled = false;
    BCType type = BCType::Robin;
    double a = 1.0;
    double b = 0.0;
    double g = 0.0;

    // Optional time-varying Dirichlet value for embedded boundaries:
    //   g(t) = dc + ac * sin(2*pi*freq*(t - t0) + phase)
    // Only applied when type == Dirichlet.
    BoundarySide::TimeFunction timeDirichlet;
};

} // namespace Poisson
} // namespace Solvers

// -----------------------------------------------------------------------------
// Compatibility aliases (to minimize churn in ExplicitPIC for now)
// -----------------------------------------------------------------------------
using PicBCType = Solvers::Poisson::BCType;
using PicBoundarySide = Solvers::Poisson::BoundarySide;
using PicDomainBC = Solvers::Poisson::DomainBC;
using PicEmbeddedBC = Solvers::Poisson::EmbeddedBC;
using PicEmbeddedBCMap = std::unordered_map<int, PicEmbeddedBC>;
