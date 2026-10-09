#include "Modules/Electricity/SurfaceChargeRelaxation.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Quantum::Electricity
{
namespace
{
void require(bool condition, const char* message)
{
    if (!condition) throw std::invalid_argument(message);
}

double finite(double value, const char* message)
{
    require(std::isfinite(value), message);
    return value;
}
}

SurfaceChargeRelaxationResult relaxSurfaceCharge(
    const SurfaceChargeTransport& history, const std::vector<double>* previousCharge,
    double dt, double previousDt,
    const std::function<std::vector<double>(const std::vector<double>&)>& evaluateCurrent,
    const SurfaceChargeRelaxationOptions& options)
{
    SurfaceChargeRelaxationResult result;
    auto& report = result.report;
    report.stage = "input";
    try
    {
        require(std::isfinite(dt) && dt > 0.0, "charge dt must be finite and positive");
        require(std::isfinite(previousDt) && (previousCharge ? previousDt > 0.0 : previousDt == 0.0),
            "previous charge and positive previous dt must be supplied together");
        require(static_cast<bool>(evaluateCurrent), "charge relaxation requires an actual normal-current evaluator");
        require(std::isfinite(options.underRelaxation) && options.underRelaxation > 0.0 &&
            options.underRelaxation <= 1.0 && options.maxIterations > 0,
            "under-relaxation must be in (0,1] and iteration limit positive");
        require(std::isfinite(options.densityChangeTolerance) && options.densityChangeTolerance > 0.0 &&
            std::isfinite(options.relativeEquationTolerance) && options.relativeEquationTolerance > 0.0 &&
            options.relativeEquationTolerance < 1.0, "invalid charge iteration change/equation tolerances");
        const auto& old = history.panelCharge();
        const auto& area = history.panelMeasure();
        if (previousCharge)
            (void)SurfaceChargeTransport::fromIntegratedCharge(area,*previousCharge,history.closed());
        // Same coefficients as the validated density update, now acting on Q.
        // Use the step ratio to avoid overflowing 2*dt or dt+previousDt.
        const double ratio = previousCharge ? finite(dt/previousDt,"charge step ratio overflow") : 0.0;
        require(!previousCharge || ratio > 0.0, "charge step ratio underflow");
        report.temporalOrder = previousCharge ? 2 : 1;
        report.alpha = previousCharge ? 1.0 + ratio/(1.0+ratio) : 1.0;
        report.beta = previousCharge ? -ratio/(1.0+ratio) : 0.0;
        const double historyWeight = finite(report.beta*ratio,"charge history weight overflow");
        auto trial = old;
        std::vector<double> next(trial.size());
        for (int iteration = 0; iteration < options.maxIterations; ++iteration)
        {
            report.iterations = iteration+1;
            report.stage = "normal-current";
            auto current = evaluateCurrent(trial);
            require(current.size() == trial.size(), "normal current must cover every material panel");
            report.stage = "charge-equation";
            report.maxDensityChange = report.maxEquationDensityResidual = report.equationDensityScale = 0.0;
            long double integratedResidual = 0.0L;
            for (std::size_t i = 0; i < trial.size(); ++i)
            {
                finite(current[i],"nonfinite normal current");
                const double historyTerm = previousCharge
                    ? finite(historyWeight*(old[i]-(*previousCharge)[i]),"charge history term overflow") : 0.0;
                const double increment = finite(dt*current[i],"charge current increment overflow");
                require(current[i] == 0.0 || increment != 0.0, "nonzero charge current increment underflow");
                const double residual = finite(report.alpha*(trial[i]-old[i])+historyTerm-increment,
                    "charge equation residual overflow");
                // Algebraically the old q_calc=b/a followed by omega mixing.
                // Correct around the current iterate, without a duplicate
                // area-expansion term or a second material advection update.
                const double correction = finite(-residual/report.alpha,"charge correction overflow");
                next[i] = finite(trial[i]+options.underRelaxation*correction,"relaxed charge overflow");
                const double change = finite(std::abs(next[i]-trial[i])/area[i],"charge density change overflow");
                const double densityResidual = finite(std::abs(residual)/area[i],"charge density residual overflow");
                const double densityScale = finite((report.alpha*std::max(std::abs(trial[i]),std::abs(old[i]))
                    +std::abs(historyTerm)+std::abs(increment))/area[i],"charge equation scale overflow");
                report.maxDensityChange = std::max(report.maxDensityChange,change);
                report.maxEquationDensityResidual = std::max(report.maxEquationDensityResidual,densityResidual);
                report.equationDensityScale = std::max(report.equationDensityScale,densityScale);
                integratedResidual += residual;
            }
            report.integratedEquationResidual = finite(static_cast<double>(integratedResidual),
                "integrated charge equation residual overflow");
            report.relativeEquationResidual = report.equationDensityScale > 0.0
                ? report.maxEquationDensityResidual/report.equationDensityScale : 0.0;
            finite(report.relativeEquationResidual,"relative charge equation residual overflow");
            if (report.maxDensityChange <= options.densityChangeTolerance &&
                report.maxEquationDensityResidual <= options.relativeEquationTolerance*report.equationDensityScale)
            {
                // Current and the electric callback's last field correspond to
                // trial. Do not publish next with a one-iteration-old E field.
                report.stage = "complete";
                report.reason = "charge change and full material time equation converged at the same current/charge iterate";
                report.success = true;
                result.integratedCharge = std::move(trial);
                result.integratedNormalCurrent = std::move(current);
                return result;
            }
            // Reuse the existing signed-integral/density validation; never
            // silently clamp negative charge or an unrepresentable density.
            (void)SurfaceChargeTransport::fromIntegratedCharge(area,next,history.closed());
            trial.swap(next);
        }
        report.stage = "iteration-limit";
        report.reason = "surface-charge under-relaxation did not satisfy both convergence gates before maxIterations";
    }
    catch (const std::exception& error)
    {
        report.reason = error.what();
    }
    return result;
}
}
