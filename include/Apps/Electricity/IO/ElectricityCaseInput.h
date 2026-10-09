#pragma once
#include "Modules/Electricity/IO/DictionaryDielectric.h"
#include "Runtime/RunControls.h"

struct ElectricityCaseInput {
    RunControls run;
    DielectricDomainInput domain;
};
ElectricityCaseInput readElectricityCaseInput();
