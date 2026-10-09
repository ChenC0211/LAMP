#include "Modules/Hydrodynamics/HydroNodalTraction.h"
#include "Modules/Interface/Compatibility/UniformInterfaceContour.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Field/CellField.h"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace NodalTractionDetail
{
[[noreturn]] inline void fail(const std::string& stage, const std::string& reason)
{
    throw std::runtime_error("[SharpInterface:Hydro " + stage + "] " + reason);
}

template<typename Mesh>
bool validCell(const Mesh& mesh, int label)
{
    if (label < 0 || label >= mesh.cellLabelCapacity()) return false;
    const int index = mesh.cells(label).index();
    return index >= 0 && index < mesh.cellIndexNum() && mesh.cellIndexToLabel[index] == label;
}

template<typename Mesh>
double volume(const Mesh& mesh, int index, bool axisymmetric)
{
    const int label = mesh.cellIndexToLabel[index];
    if (!validCell(mesh, label)) fail("forcing", "invalid cell index-to-label mapping");
    const auto& cell = mesh.cells(label);
    const double result = cell.area() * (axisymmetric ? 2.0 * Pi * cell.centerCoord().x : 1.0);
    if (!std::isfinite(result) || result <= 0.0)
        fail("forcing", "cell area/axisymmetric volume must be finite and positive");
    return result;
}
}


namespace HydroNodalTraction
{
Report spread(QuadtreeMesh& mesh, const UniformInterfaceContour& contour,
    const Field<double>& tractionX, const Field<double>& tractionY,
    CellField<double>& fx, CellField<double>& fy)
{
    using namespace NodalTractionDetail;
    contour.validateMesh(mesh);
    const bool axisymmetric = contour.axisymmetric();
    const auto& nodes = contour.nodes();
    if (tractionX.size() != static_cast<int>(nodes.size()) || tractionY.size() != static_cast<int>(nodes.size())
        || !tractionX.value || !tractionY.value)
        fail("forcing", "traction fields do not match the rebuilt interface nodes");
    if (fx.size() != mesh.cellIndexNum() || fy.size() != mesh.cellIndexNum() || !fx.value || !fy.value)
        fail("forcing", "output force fields must match the current mesh cells");
    if (fx.value.get() == fy.value.get())
        fail("forcing", "x/y output force fields must not alias");
    std::vector<double> volumes(mesh.cellIndexNum()), deltaX(mesh.cellIndexNum(), 0.0), deltaY(mesh.cellIndexNum(), 0.0);
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        volumes[i] = volume(mesh, i, axisymmetric);
        if (!std::isfinite(fx.value[i]) || !std::isfinite(fy.value[i]))
            fail("forcing", "existing output force field contains non-finite values");
    }
    long double inputX = 0.0L, inputY = 0.0L;
    for (const auto& node : nodes)
    {
        const int index = node.index();
        if (index < 0 || index >= tractionX.size()) fail("forcing", "invalid traction field index");
        const double tx = tractionX.value[index], ty = tractionY.value[index];
        const double measure = node.ds() * (axisymmetric ? 2.0 * Pi * node.coord().x : 1.0);
        if (!std::isfinite(tx) || !std::isfinite(ty) || !std::isfinite(measure) || measure <= 0.0)
            fail("forcing", "traction and positive interface measure must be finite");
        const double forceX = tx * measure, forceY = ty * measure;
        if (!std::isfinite(forceX) || !std::isfinite(forceY)) fail("forcing", "integrated node force overflow");
        inputX += forceX; inputY += forceY;
        double weightSum = 0.0;
        for (int j = 0; j < 25; ++j)
        {
            const double weight = node.cellStencilCoeff[j];
            if (!std::isfinite(weight) || weight < 0.0) fail("forcing", "invalid spreading coefficient");
            if (weight == 0.0) continue;
            const int cell = node.cellStencil[j];
            if (cell < 0 || cell >= mesh.cellIndexNum()) fail("forcing", "invalid nonzero spreading support index");
            weightSum += weight;
            deltaX[cell] += forceX * weight / volumes[cell];
            deltaY[cell] += forceY * weight / volumes[cell];
        }
        if (!std::isfinite(weightSum) || std::abs(weightSum - 1.0) > 1e-12)
            fail("forcing", "spreading weights do not sum to one; rebuild the interface");
    }
    long double outputX = 0.0L, outputY = 0.0L;
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        if (!std::isfinite(deltaX[i]) || !std::isfinite(deltaY[i])
            || !std::isfinite(fx.value[i] + deltaX[i]) || !std::isfinite(fy.value[i] + deltaY[i]))
            fail("forcing", "spread/accumulated force field overflow");
        outputX += static_cast<long double>(deltaX[i]) * volumes[i];
        outputY += static_cast<long double>(deltaY[i]) * volumes[i];
    }
    Report report;
    report.forceX = static_cast<double>(inputX); report.forceY = static_cast<double>(inputY);
    report.spreadForceX = static_cast<double>(outputX); report.spreadForceY = static_cast<double>(outputY);
    if (!std::isfinite(report.forceX) || !std::isfinite(report.forceY)
        || !std::isfinite(report.spreadForceX) || !std::isfinite(report.spreadForceY))
        fail("forcing", "integrated surface force overflow");
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        fx.value[i] += deltaX[i];
        fy.value[i] += deltaY[i];
    }
    return report;
}

}
