#pragma once

#include "Core/Base/ScalarExpression.h"
#include "Core/IO/ConfigurationTrace.h"
#include "Modules/Electric/ElectricPotential.h"
#include "Modules/Electric/ElectricMaxwellTrace.h"
#include "Modules/Electric/ElectricPhaseGeometry.h"
#include "Modules/Interface/HermiteGalsPublication.h"
#include "RunControls.h"
#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

// Formal fixed-geometry sharp-interface electrostatics data. Dictionary paths belong only to this
// reader; the App consumes typed controls and immutable expression providers.
struct ElectricCaseInput
{
    RunControls run;
    bool radial = false;
    Vector2D domainMin{0,0}, domainMax{1,1};
    int nx = 0, ny = 0, negativePhase = 1, positivePhase = 2;
    int parallelThreads = 2;
    LAMP::Electric::ElectricPhaseProblem electric;
    std::array<double,2> conductivity{{0,0}};
    LAMP::Electric::ElectricPhaseGeometryConfig geometry;
    LAMP::Interface::HermiteGalsPublicationConfig publication;
    LAMP::Electric::ElectricCartesianControls linear;
    LAMP::Electric::ElectricChargeTimeControls charge;
    LAMP::Electric::ElectricSurfaceChargeMomentControls moments;
    LAMP::Electric::ElectricMaxwellControls maxwell;
    DictionaryInput::ConfigurationTrace trace;
    struct Document { std::string file, sourceText, effectiveText; };
    std::vector<Document> documents;
    LAMP::ScalarExpression::Jet initialJet(Vector2D) const;
    double initialChargeDensity(Vector2D) const;
private:
    std::shared_ptr<const LAMP::ScalarExpression> initialPhi_, initialCharge_;
    friend ElectricCaseInput readElectricCaseInput(std::optional<int>);
};
ElectricCaseInput readElectricCaseInput(std::optional<int> threadsOverride = {});
