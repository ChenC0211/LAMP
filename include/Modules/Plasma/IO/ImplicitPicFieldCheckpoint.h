#pragma once
#include <filesystem>
#include "Core/Base/Vector2D.h"
class QuadtreeMesh;
template<typename T> class CellField;
namespace ImplicitPicFieldCheckpoint {
void write(const QuadtreeMesh& checkpointMesh,
    const std::filesystem::path& dir,
    const CellField<double>& phiField,
    const CellField<double>& rhsField,
    const CellField<double>& rhoField,
    const CellField<double>& qField,
    const CellField<Vector2D>& ecellField,
    const CellField<Vector2D>& ecellNewField,
    const CellField<Vector2D>& ehalfField);
bool read(const QuadtreeMesh& checkpointMesh,
    const std::filesystem::path& dir,
    CellField<double>& phiField,
    CellField<double>& rhsField,
    CellField<double>& rhoField,
    CellField<double>& qField,
    CellField<Vector2D>& ecellField,
    CellField<Vector2D>& ecellNewField,
    CellField<Vector2D>& ehalfField);
}
