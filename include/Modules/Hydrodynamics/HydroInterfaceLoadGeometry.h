#pragma once

#include <memory>
#include <vector>
#include "Core/Base/Enumeration.h"
#include "Core/Geometry/InterfaceGeometry.h"
#include "Core/Geometry/InterfaceTraceGeometry.h"

class QuadtreeMesh;
template<typename T> class VertField;

// Immutable geometry for a single fixed-uniform interface-load state.
// Geometric phase fraction retains the original SW--NE PL definition. The
// load indicator defaults to that fraction; an explicit compact C2 layer can
// instead regularize the same Hermite scalar. Neither is transported VOF.
// Optional axisymmetric RZ (x=r,y=z) currently requires Regularized and an
// interface/band separated from positive-area walls unless explicit constrained
// normal-velocity intervals and boundary contacts are supplied. At r=0 the common scalar
// uses even radial symmetry; smooth orthogonal axis contact is permitted.
// The retained PL fraction and dualArea remain MERIDIONAL two-dimensional
// quantities: positiveFraction is not an r-weighted RZ volume fraction.
// No material, traction coefficient, electric solver, clock or shape lives here.
class HydroInterfaceLoadGeometry
{
public:
    enum class Layer { GeometricPL, Regularized };
    struct NormalVelocityBoundary
    {
        Boundary side = Boundary::non;
        double start = 0.0, end = 0.0; // Entire native boundary-face intervals.
    };
    struct Options
    {
        Layer layer = Layer::GeometricPL;
        double halfWidth = 0.0; // Regularized: positive, in phi units; PL: zero.
        bool includeTangentialLoads = false; // Also retain the rotated G(H) components.
        bool axisymmetric = false; // RZ: r>=0, Regularized only; even scalar at r=0.
        // Opt-in for a resolved open interface at a prescribed-normal-velocity
        // boundary. The caller must match these intervals to its actual
        // constrained velocity rows; pressure/traction outlets are NOT covered.
        // Does not provide mechanical wall reaction or a contact-angle law.
        // Currently regularized normal loads only; no tangential zero extension.
        std::vector<NormalVelocityBoundary> normalVelocityBoundaries;
        std::vector<Vector2D> boundaryContacts;
        // Explicit complete equal-level leaf layout, retaining actual native
        // labels/indices. Not hanging-node AMR or active-fluid scalar support.
        CartesianGridScope gridScope = CartesianGridScope::RootOnly;
        // Explicit scalar selection; never a face-only curvature replacement.
        // Cubic retains the original separately reconstructed D3. Quartic
        // uses D1/D2/D3 of its same selected cell and requires 7/8-node data.
        LevelSetHermiteGeometry::Reconstruction reconstruction = LevelSetHermiteGeometry::Reconstruction::Cubic;
    };
    struct Face
    {
        bool horizontal = false; // horizontal edge carries the y component.
        int edgeIndex = -1;
        // Low/high refer to Cartesian position (slave/master), NOT phase.
        int lowCellIndex = -1, highCellIndex = -1;
        Vector2D point;
        double gradient = 0.0, dualArea = 0.0;
        LevelSetHermiteGeometry::Frame frame; // Always meridional curvature, including in RZ.
        double curvatureTangentialDerivative = 0.0; // Meridional d_s kappa for trace reconstruction.
        // Four-face average of the OTHER pressure-gradient component:
        // x face: Gy(H); y face: Gx(H). No phase-normal division or new delta.
        // R(GH) contributes -transverse on x faces, +transverse on y faces.
        double transverseGradient = 0.0;
        // Planar: kappa_m and dualArea. RZ: kappa_m+n_r/r at the projected
        // interface point and 2*pi*r_face*dualArea at the native load face.
        double surfaceCurvature = 0.0;
        double dualMeasure = 0.0;
        bool boundaryEndpointProjection = false;
        double endpointKktResidual = 0.0;
    };

    HydroInterfaceLoadGeometry(const QuadtreeMesh& mesh, const VertField<double>& phi);
    HydroInterfaceLoadGeometry(const QuadtreeMesh& mesh, const VertField<double>& phi, Options options);
    ~HydroInterfaceLoadGeometry();
    HydroInterfaceLoadGeometry(HydroInterfaceLoadGeometry&&) noexcept;
    HydroInterfaceLoadGeometry& operator=(HydroInterfaceLoadGeometry&&) noexcept;
    HydroInterfaceLoadGeometry(const HydroInterfaceLoadGeometry&) = delete;
    HydroInterfaceLoadGeometry& operator=(const HydroInterfaceLoadGeometry&) = delete;

    const InterfaceGeometry& geometry() const;
    const LevelSetHermiteGeometry::Snapshot& snapshot() const;
    const std::vector<double>& positiveFraction() const; // Meridional geometric PL, NOT RZ volume fraction; cell solver order.
    const std::vector<double>& loadIndicator() const; // H used by the active pressure-face G(H).
    const Options& options() const;
    // x faces then y. Default: nonzero primary G(H); tangential option:
    // nonzero primary OR transverse G(H). Boundary-normal rows stay excluded.
    const std::vector<Face>& faces() const;
    const std::shared_ptr<const InterfaceTraceGeometry>& traceGeometry() const;
    // Reject another mesh, changed consumed mesh metadata, changed phi bits,
    // or moved-from state. Never updates/rebuilds this retained snapshot.
    void validateSource(const QuadtreeMesh& mesh, const VertField<double>& phi) const;

private:
    struct Impl;
    std::unique_ptr<const Impl> impl_;
    std::shared_ptr<const InterfaceTraceGeometry> traceGeometry_;
    const Impl& checked() const;
};
