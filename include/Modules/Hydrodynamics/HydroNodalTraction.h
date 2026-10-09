#pragma once

class QuadtreeMesh;
class UniformInterfaceContour;
template<typename T> class Field;
template<typename T> class CellField;

// Existing normalized Peskin quadrature for a fixed uniform closed contour.
// Traction components are supplied in code units and indexed by contour node.
// Accumulate into caller fields only after the complete load is validated.
namespace HydroNodalTraction
{
struct Report
{
    double forceX = 0.0, forceY = 0.0;
    double spreadForceX = 0.0, spreadForceY = 0.0;
};
Report spread(QuadtreeMesh& mesh, const UniformInterfaceContour& contour,
    const Field<double>& tractionX, const Field<double>& tractionY,
    CellField<double>& forceX, CellField<double>& forceY);
}
