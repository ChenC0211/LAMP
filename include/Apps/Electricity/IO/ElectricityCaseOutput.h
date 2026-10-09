#pragma once
#include "Electricity/IO/ElectricityCaseInput.h"
#include "Modules/Electricity/DielectricEvolution.h"
#include "Core/IO/TecplotWriter.h"
#include <memory>

class ElectricityCaseOutput {
public:
    ElectricityCaseOutput(const ElectricityCaseInput&, int finalStep);
    void write(const QuadtreeMesh&,
        const Quantum::Electricity::DielectricEvaluation&,
        int step, const Quantum::Electricity::DielectricChargeUpdate* update=nullptr);
private:
    RunControls run_;
    const int finalStep_;
    bool radial_=false;
    std::unique_ptr<TecplotWriter::History> history_;
};
