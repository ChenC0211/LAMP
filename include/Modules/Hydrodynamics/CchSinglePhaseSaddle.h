#pragma once

#include "Modules/Hydrodynamics/CchPressureLoad.h"
#include <memory>
#include <vector>

namespace Quantum::Hydrodynamics
{
// Caller-owned physics prepares A_u; this contract copies and validates it.
class CchSaddleVelocityBlock
{
public:
    using Matrix = CchPressureLoad::Matrix;
    static CchSaddleVelocityBlock build(CchDofLayout layout, const Matrix& integralVelocity,int tangentialCount=0);
    // Explicit coercive extension, retaining the certified SPD part. Supplied
    // transport must satisfy the original 1e-12 skew tolerance; no projection.
    CchSaddleVelocityBlock withSkew(const Matrix& integralSkew) const;
    const CchDofLayout& layout() const noexcept;
    const Matrix& matrix() const noexcept;
    double symmetryDefect() const noexcept;
    int tangentialCount() const noexcept;
private:
    struct Data;
    explicit CchSaddleVelocityBlock(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
};

class CchSaddleRhs
{
public:
    CchSaddleRhs(CchDofLayout layout, const std::vector<double>& integralMomentum,
        const std::vector<double>& prescribedNormalVelocity,const std::vector<double>& tangentialMomentum={});
    const std::vector<double>& integralMomentum() const noexcept { return momentum_; }
    const std::vector<double>& prescribedNormalVelocity() const noexcept { return prescribed_; }
    const std::vector<double>& tangentialMomentum() const noexcept {return tangential_;}
    void validateFor(const CchDofLayout& expected, const QuadtreeMesh& mesh) const;
private:
    CchDofLayout layout_;
    std::vector<double> momentum_, prescribed_,tangential_;
};

struct CchSaddleConfig
{
    double algebraicRelativeTolerance=1e-11, momentumRelativeTolerance=1e-11;
    double continuityRelativeTolerance=1e-11, pressureMeanRelativeTolerance=1e-11;
    double boundaryFluxRelativeTolerance=1e-12;
};
struct CchSaddleReport
{
    double algebraicRelativeResidual=0, freeMomentumDefect=0, maxCellContinuityDefect=0;
    double maxPressureMeanDefect=0, maxBoundaryFluxDefect=0, maxAbsGaugeMultiplier=0;
    int freeVelocityCount=0, pressureGaugeCount=0,freeTangentialCount=0;
    int linearRefinements=0;
};
class CchSaddleResult
{
public:
    const CchNativeValues& values() const noexcept { return values_; }
    const CchSaddleReport& report() const noexcept { return report_; }
    const std::vector<double>& tangentialVelocity() const noexcept {return tangential_;}
private:
    friend class CchSinglePhaseSaddle;
    CchSaddleResult(CchNativeValues values, CchSaddleReport report,std::vector<double> tangential={});
    CchNativeValues values_;
    CchSaddleReport report_;
    std::vector<double> tangential_;
};

// Closed-domain monolithic reference solve. Boundary velocities are imposed;
// each pressure component has volume-weighted mean zero, never a cell pin.
class CchSinglePhaseSaddle
{
public:
    using Matrix = Eigen::SparseMatrix<double, Eigen::ColMajor, Eigen::Index>;
    static CchSinglePhaseSaddle buildClosed(CchPressureLoad pressure,
        CchSaddleVelocityBlock velocity, CchSaddleConfig config={});
    const CchPressureLoad& pressureLoad() const noexcept;
    const CchSaddleVelocityBlock& velocityBlock() const noexcept;
    const CchDofLayout& layout() const noexcept;
    const Matrix& matrix() const noexcept;
    const std::vector<int>& freeVelocityDofs() const noexcept;
    const std::vector<int>& pressureComponents() const noexcept;
    int pressureGaugeCount() const noexcept;
    CchSaddleResult solve(const QuadtreeMesh& mesh, const CchSaddleRhs& rhs) const;
private:
    struct Data;
    explicit CchSinglePhaseSaddle(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
};
}
