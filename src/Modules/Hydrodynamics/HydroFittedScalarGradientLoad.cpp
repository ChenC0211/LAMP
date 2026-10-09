#include "Modules/Hydrodynamics/HydroFittedScalarGradientLoad.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroNativePressureStencil.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
void require(bool ok, const char* message)
{
    if (!ok) throw std::invalid_argument(std::string("[HydroFittedScalarGradientLoad] ") + message);
}
}

HydroFittedScalarGradientLoad::HydroFittedScalarGradientLoad(const QuadtreeMesh& mesh)
    : forceX(mesh, 0.0), forceY(mesh, 0.0) {}

HydroFittedScalarGradientLoad buildFittedScalarGradientLoad(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const NativeLeafDomain& sourceGeometry,
    const std::vector<double>& scalarByCell, const std::vector<double>& coefficientByNative,
    double scale)
{
    domain.validateMeshIdentity(mesh);
    sourceGeometry.validateMeshIdentity(mesh);
    require(domain.axisymmetric() == sourceGeometry.axisymmetric()
        && domain.cellRowByIndex() == sourceGeometry.cellRowByIndex()
        && domain.cellIndexByRow() == sourceGeometry.cellIndexByRow(),
        "source and receiving selected cells or coordinate measures disagree");
    require(scalarByCell.size() == static_cast<std::size_t>(mesh.cellIndexNum())
        && coefficientByNative.size() == domain.faces().size(), "source dimensions disagree");
    require(std::isfinite(scale) && scale > 0.0, "positive finite conversion scale required");
    Eigen::VectorXd scalar(domain.cellIndexByRow().size());
    for (int k : domain.cellIndexByRow()) {
        require(std::isfinite(scalarByCell[k]), "non-finite active scalar");
        scalar[domain.cellRowByIndex()[k]] = scalarByCell[k];
    }
    const Eigen::VectorXd gradientLoad = domain.pressureGradientLoad() * scalar;
    require(gradientLoad.allFinite(), "non-finite pressure-compatible scalar gradient");
    HydroFittedScalarGradientLoad result(mesh);
    double identityScale = 0.0;
    long double integralX = 0.0L, integralY = 0.0L;
    std::map<std::pair<int,int>,supCoeff> centers;
    const auto childCenter = [&](const Cell& cell,int child) -> const supCoeff& {
        const auto key=std::make_pair(cell.index(),child);
        const auto found=centers.find(key);
        if(found!=centers.end())return found->second;
        auto stencil=HydroNativePressureStencil::childCenter(mesh,cell,child);
        for(int column:stencil.ind)
            require(column>=0 && column<mesh.cellIndexNum() && domain.cellRowByIndex()[column]>=0,
                "native scalar pressure interpolation consumes an inactive column");
        return centers.emplace(key,std::move(stencil)).first->second;
    };
    for (const auto& face : domain.faces()) {
        const double coefficient = coefficientByNative[face.nativeRow];
        require(std::isfinite(coefficient), "non-finite load coefficient");
        if (face.negativeCell < 0 || face.positiveCell < 0) {
            require(coefficient == 0.0, "boundary coefficient requires a separate traction closure");
            continue;
        }
        if (face.velocityRow < 0) continue;
        const int low = face.negativeCell, high = face.positiveCell;
        const double integrated = gradientLoad[face.velocityRow];
        double independent = face.area * (scalarByCell[high] - scalarByCell[low]);
        const bool coarseFine = mesh.cells(mesh.cellIndexToLabel[low]).level()
            != mesh.cells(mesh.cellIndexToLabel[high]).level();
        if(coarseFine) {
            long double sum=0.0L,absoluteSum=0.0L;
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            HydroNativePressureStencil::visitGradient(mesh,edge,face.horizontal,childCenter,
                [&](int column,double gradient) {
                    const long double term=static_cast<long double>(face.dualVolume*gradient)*scalarByCell[column];
                    sum+=term;absoluteSum+=std::abs(term);
                });
            independent=static_cast<double>(sum);
            const double arithmeticScale=static_cast<double>(absoluteSum);
            require(std::isfinite(independent) && std::isfinite(arithmeticScale),
                "unrepresentable independent native scalar pressure load");
            // Multi-column interpolation may cancel a constant only to roundoff.
            // Scale the SAME arithmetic identity gate by the terms being
            // summed; do not replace or zero the actual matrix product.
            identityScale=std::max(identityScale,arithmeticScale);
        }
        result.maxGradientIdentityError = std::max(result.maxGradientIdentityError, std::abs(integrated - independent));
        identityScale = std::max({identityScale, std::abs(integrated), std::abs(independent)});
        require(coarseFine || scalarByCell[low] != scalarByCell[high] || integrated == 0.0,
            "uniform scalar has a spurious pressure gradient");
        if (integrated == 0.0) continue;
        require(std::isfinite(face.dualVolume) && face.dualVolume > 0.0, "invalid momentum dual volume");
        const double gradient = integrated / face.dualVolume;
        const double density = coefficient * gradient;
        const double force = scale * density;
        require(std::isfinite(force) && (density == 0.0 || force != 0.0), "load conversion overflow or underflow");
        if (face.horizontal) result.forceY.value[face.nativeRow - mesh.vEdgeIndexNum()] = force;
        else result.forceX.value[face.nativeRow] = force;
        result.faces.push_back({face.nativeRow, low, high, gradient, face.dualVolume, force});
        (face.horizontal ? integralY : integralX) += static_cast<long double>(force) * face.dualVolume;
        result.maxAbsComponent = std::max(result.maxAbsComponent, std::abs(force));
    }
    require(result.maxGradientIdentityError <= 512.0 * std::numeric_limits<double>::epsilon() * identityScale,
        "scalar gradient disagrees with native pressure-gradient load");
    result.integralX = static_cast<double>(integralX); result.integralY = static_cast<double>(integralY);
    require(std::isfinite(result.integralX) && std::isfinite(result.integralY), "face-load integral overflow");
    return result;
}
