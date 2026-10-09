#include "Modules/Interface/Compatibility/Detail/InterfaceCorrectionRegularization.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <cmath>
#include <stdexcept>

namespace InterfaceVolumeCorrection
{
QuarticRegularization geometricRegularization(const QuadtreeMesh& mesh,
    const InterfaceState& state,double strength)
{
    if (!std::isfinite(strength) || strength<0)
        throw std::invalid_argument("geometric regularization strength must be finite and nonnegative");
    QuarticRegularization rows;
    if (strength==0) return rows;
    std::vector<LevelSetInterfacePanels::ArcSample> samples;
    long double length=0,gradientSquare=0;
    for (const auto& surface:state.topology().surfaces()) {
        const auto& chain=surface.geometry.nativeChain();
        for (int k=0;k<static_cast<int>(chain.graphs().size());++k)
            for (auto sample:chain.quadrature(k)) {
                if (!(sample.gradientMagnitude>0) || !std::isfinite(sample.gradientMagnitude))
                    throw std::runtime_error("geometric regularization encountered a critical interface");
                length+=sample.weight;
                gradientSquare+=sample.weight*static_cast<long double>(sample.gradientMagnitude)*sample.gradientMagnitude;
                samples.push_back(std::move(sample));
            }
    }
    if (!(length>0) || !std::isfinite(gradientSquare))
        throw std::runtime_error("invalid interface scale for geometric regularization");
    const long double reference=std::sqrt(gradientSquare/length);
    rows.reserve(3*samples.size());
    for (const auto& q:samples) {
        const auto& f=q.frame;
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[f.nativeCellIndex]);
        const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
        const long double h=std::sqrt((static_cast<long double>(b.x)-a.x)*(static_cast<long double>(b.y)-a.y));
        const long double scale=reference/q.gradientMagnitude*std::sqrt(strength*q.weight/h);
        const long double nx=f.normal.x,ny=f.normal.y,tx=-ny,ty=nx;
        const auto chart=LevelSetHermiteGeometry::cellChartPoint(f);
        // phi, h*d_s(phi), h^2*d_ss(phi) increments along the current curve.
        // d_ss(delta phi)=t.H(delta phi).t-kappa*n.grad(delta phi).
        // These are a regularity penalty, not a replacement curvature model.
        rows.push_back(state.geometry()->nativeQuarticDifferentialSensitivity(f.nativeCellIndex,chart,
            {{scale,0,0,0,0,0}}));
        rows.push_back(state.geometry()->nativeQuarticDifferentialSensitivity(f.nativeCellIndex,chart,
            {{0,scale*h*tx,scale*h*ty,0,0,0}}));
        rows.push_back(state.geometry()->nativeQuarticDifferentialSensitivity(f.nativeCellIndex,chart,
            {{0,-scale*h*h*f.curvature*nx,-scale*h*h*f.curvature*ny,
                scale*h*h*tx*tx,2*scale*h*h*tx*ty,scale*h*h*ty*ty}}));
    }
    return rows;
}
}
