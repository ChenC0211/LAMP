#pragma once
#include "Hydrodynamics/IO/HydrodynamicsInput.h"
#include "Electrohydrodynamics/ElectrohydrodynamicStep.h"
#include <memory>
#include <optional>
#include <string>
#include <vector>

// Joint input belongs only to its standalone App. Common fluid data has no
// charge or electric configuration; this owner carries the complete provenance.
struct ElectrohydrodynamicsInput
{
    HydrodynamicsInput fluid;
    LAMP::Coupling::ElectrohydrodynamicStepConfig electric;
    DictionaryInput::ConfigurationTrace trace;
    using Document=AppIO::CaseDocument;
    std::vector<Document> documents;
    double initialSurfaceChargeDensity(Vector2D) const;
private:
    std::shared_ptr<const LAMP::ScalarExpression> chargeExpression_;
    friend ElectrohydrodynamicsInput readElectrohydrodynamicsInput(const std::string&,std::optional<int>);
};
ElectrohydrodynamicsInput readElectrohydrodynamicsInput(const std::string& systemDir="system",
    std::optional<int> threadsOverride={});
