#include "Modules/Electricity/DielectricElectrostatics.h"
#include <cmath>
#include <stdexcept>
namespace {
void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}
}
void validateDielectricConfig(const DielectricConfig& config)
{
    require(!config.nativeAxisNeckHarmonics || (config.axisymmetric && config.interfaceDiscretization=="sharpHermiteP3"),
        "native axis-neck harmonics require sharp axisymmetric electrostatics");
    const auto validTolerance = [](double value) {
        return std::isfinite(value) && value > 0.0 && value < 1.0;
    };
    if (config.boundaryFluxTolerance)
    {
        require(config.interfaceDiscretization == "diffuseTanh",
            "boundaryFlux is only consumed by diffuseTanh electrostatics");
        require(validTolerance(*config.boundaryFluxTolerance), "boundaryFlux must be finite and in (0,1)");
    }
    if (config.chargeClosureTolerance)
    {
        require(config.model == "leakyDielectric", "chargeClosure requires leakyDielectric");
        require(validTolerance(*config.chargeClosureTolerance), "chargeClosure must be finite and in (0,1)");
    }
    if (config.sharpFluxChecks)
    {
        require(config.interfaceDiscretization == "sharpHermiteP3",
            "explicit sharpFluxChecks are not consumed by diffuse electrostatics");
        config.sharpFluxChecks->validate();
    }
    require(config.model == "perfectDielectric" || config.model == "leakyDielectric",
        "model must be perfectDielectric or leakyDielectric");
    require(config.interfaceDiscretization == "diffuseTanh" ||
            config.interfaceDiscretization == "sharpHermiteP3",
            "interfaceDiscretization must be diffuseTanh or sharpHermiteP3");
    require(config.electrodeAxis == "x" || config.electrodeAxis == "y", "electrodeAxis must be x or y");
    require(!config.axisymmetric || config.interfaceDiscretization == "sharpHermiteP3",
        "RZ electrostatics requires sharpHermiteP3; diffuse RZ electrostatics is unsupported");
    require(std::isfinite(config.permittivityNegative) && config.permittivityNegative > 0.0 &&
            std::isfinite(config.permittivityPositive) && config.permittivityPositive > 0.0,
            "both permittivities must be finite and positive in consistent code units");
    require(std::isfinite(config.lowerPotential) && std::isfinite(config.upperPotential) &&
            std::isfinite(config.upperPotential - config.lowerPotential), "electrode potentials and their difference must be finite");
    require(std::isfinite(config.linearResidualTolerance) && config.linearResidualTolerance > 0.0 &&
            config.linearResidualTolerance <= 1e-2, "linearResidualTolerance must be in (0,0.01]");
    if (config.model == "leakyDielectric")
    {
        require(config.interfaceDiscretization == "sharpHermiteP3", "leakyDielectric requires sharpHermiteP3");
        require(std::isfinite(config.conductivityNegative) && config.conductivityNegative >= 0.0 &&
            std::isfinite(config.conductivityPositive) && config.conductivityPositive >= 0.0 &&
            std::isfinite(config.initialSurfaceChargeDensity), "invalid conductivity/initial surface charge");
    }
}
