#include "Core/Geometry/FittedDomainGeometry.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Field/CellField.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace FittedDomainGeometry
{
// Evaluate complete phase expressions on interpolated cell-center level sets;
// averaging vertex phase IDs would assign the curved liquid cap incorrectly.
std::vector<int> cellPhases(const QuadtreeMesh& mesh,
    const ImmersedGeometry<QuadtreeMesh>& geometry)
{
    std::vector<CellField<double>> cellLevelSets;
    for (const auto& field : geometry.levelSets)
        cellLevelSets.push_back(*field.toCellField(mesh));
    std::vector<int> result(mesh.cellIndexNum(), -1);
    const int originalCount = static_cast<int>(cellLevelSets.size());
    for (std::size_t p = 0; p < geometry.phases.size(); ++p)
    {
        const auto operations = geometry.phases[p].operation();
        const auto a = geometry.phases[p].levelSetIndexA();
        const auto b = geometry.phases[p].levelSetIndexB();
        if (operations.empty() || a.size() != operations.size() || b.size() != operations.size())
            throw std::runtime_error("invalid cell phase expression");
        for (std::size_t j = 0; j < operations.size(); ++j)
        {
            if (a[j].value < 0 || b[j].value < 0
                || a[j].value >= originalCount + static_cast<int>(j)
                || b[j].value >= originalCount + static_cast<int>(j)
                || (operations[j] != LevelSetOperation::non && operations[j] != LevelSetOperation::min
                    && operations[j] != LevelSetOperation::max)
                || (operations[j] == LevelSetOperation::non
                    && (a[j].value != b[j].value || a[j].isNegative != b[j].isNegative)))
                throw std::runtime_error("invalid cell phase operation/reference");
        }
        std::vector<double> values(originalCount + operations.size());
        for (int i = 0; i < mesh.cellIndexNum(); ++i)
        {
            for (int j = 0; j < originalCount; ++j) values[j] = cellLevelSets[j].value[i];
            for (std::size_t j = 0; j < operations.size(); ++j)
            {
                const double va = (a[j].isNegative ? -1.0 : 1.0) * values[a[j].value];
                const double vb = (b[j].isNegative ? -1.0 : 1.0) * values[b[j].value];
                if (!std::isfinite(va) || !std::isfinite(vb))
                    throw std::runtime_error("non-finite cell phase operand");
                values[originalCount + j] = operations[j] == LevelSetOperation::non ? va
                    : operations[j] == LevelSetOperation::min ? std::min(va, vb) : std::max(va, vb);
            }
            if (values.back() < 0.0)
            {
                if (result[i] != -1) throw std::runtime_error("overlapping physical phases at cell center");
                result[i] = static_cast<int>(p) + 1;
            }
        }
    }
    return result;
}

namespace
{
bool coversFace(const Line& line, const FaceSegment& segment, double tolerance)
{
    if (line.type() != LineType::straight) return false;
    const bool horizontal = segment.isHorizontal;
    const double normal0 = horizontal ? line.points[0].y : line.points[0].x;
    const double normal1 = horizontal ? line.points[1].y : line.points[1].x;
    const double faceNormal = horizontal ? segment.faceCenter.y : segment.faceCenter.x;
    const double tangent0 = horizontal ? line.points[0].x : line.points[0].y;
    const double tangent1 = horizontal ? line.points[1].x : line.points[1].y;
    const double faceTangent = horizontal ? segment.faceCenter.x : segment.faceCenter.y;
    return std::abs(normal0 - faceNormal) <= tolerance
        && std::abs(normal1 - faceNormal) <= tolerance
        && faceTangent - 0.5 * segment.Af >= std::min(tangent0, tangent1) - tolerance
        && faceTangent + 0.5 * segment.Af <= std::max(tangent0, tangent1) + tolerance;
}

}

const DictionaryGeometry::BoundaryEntity& bindFace(const DictionaryGeometry& definition,
    const FaceSegment& segment, int ownerPhase, int excludedPhase, double tolerance)
{
    const DictionaryGeometry::BoundaryEntity* found = nullptr;
    for (const auto& boundary : definition.boundaries)
    {
        if (boundary.ignored) continue;
        if (segment.L < 0)
        {
            if (boundary.kind != GeometricBoundary::domain || boundary.adjacentPhase + 1 != ownerPhase) continue;
        }
        else
        {
            if (boundary.kind != GeometricBoundary::fixed
                || !((boundary.phase + 1 == excludedPhase && boundary.opposedPhase + 1 == ownerPhase)
                    || (boundary.opposedPhase + 1 == excludedPhase && boundary.phase + 1 == ownerPhase))) continue;
        }
        for (int lineId : boundary.lineIds)
        {
            if (lineId < 0 || lineId >= static_cast<int>(definition.lines.size()))
                throw std::runtime_error("invalid physical boundary line ID");
            if (!coversFace(definition.lines[lineId], segment, tolerance)) continue;
            if (found) throw std::runtime_error("ambiguous fitted electric boundary coverage");
            found = &boundary;
        }
    }
    if (!found) throw std::runtime_error("electric boundary is missing, cuts a cell, or is not face fitted: cell="
        + std::to_string(segment.K) + " face=(" + std::to_string(segment.faceCenter.x)
        + "," + std::to_string(segment.faceCenter.y) + ")");
    return *found;
}

}
