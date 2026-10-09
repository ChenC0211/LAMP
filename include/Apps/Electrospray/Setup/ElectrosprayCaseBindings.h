#pragma once
class DictionaryGeometry;
class DictionaryMesh;
struct FittedFlowCaseConfig;
struct FittedInterfaceCaseConfig;
struct FittedElectricCaseConfig;

// Read-only dictionary identity/capability checks, without building a mesh.
// Exact leaf alignment and numerical geometry quality remain initialization checks.
void validateElectrosprayCaseBindings(const DictionaryGeometry& geometry,
    const DictionaryMesh& mesh, const FittedFlowCaseConfig& flow,
    const FittedInterfaceCaseConfig& interfaceConfig, const FittedElectricCaseConfig& electric);
