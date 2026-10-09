#pragma once

#include <string>
class QuadtreeMesh;
class UnitsContext;
class LevelSetInterfacePanels;
enum class InterfaceMotionStage;
struct ElectrosprayInterfaceState;
struct InterfaceMotionReport;

// Three short console observations per selected step, no duplicate field
// files. Actual panel representatives, not a curvature error/acceptance gate.
void writeElectrosprayInterfaceCurvature(int step, InterfaceMotionStage stage,
    const LevelSetInterfacePanels& panels, const UnitsContext& units);

// Candidate geometry diagnostic at the original pre-commit output stage. These
// files alone do NOT mark an accepted frame. The driver writes the accepted-step
// histories and Tecplot manifest later. No volume/charge acceptance is changed.
void writeElectrosprayInterfaceMotion(const QuadtreeMesh& mesh,
    const ElectrosprayInterfaceState& state, const InterfaceMotionReport& report,
    const UnitsContext& units, const std::string& prefix);
