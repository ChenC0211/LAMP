#pragma once

#include "Core/Field/CellField.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/DepositQK.h"
#include <cstdint>
#include <vector>

class QuadtreeMesh;
class FaceSegmentGraph;
class NormalizationContext;
class ParticleType;
class PicStepParticleBalance;

// Owns the orbit-averaged deposition stage and its polarization coefficients.
// The mesh outlives this object; caches/particles are borrowed only per call.
// Q/N are SI amounts, rho is normalized, mu is F/m, kappa is dimensionless.
class ImplicitPicPolarization
{
public:
    struct DepositStats { std::uint64_t particles = 0, entries = 0; };
    struct NeutralityCorrection { bool applied = false; double totalQ = 0.0, meanRhoSI = 0.0; };

    ImplicitPicPolarization(QuadtreeMesh& mesh, std::size_t speciesCount, std::size_t segmentCount);
    ImplicitPicPolarization(const ImplicitPicPolarization&) = delete;
    ImplicitPicPolarization& operator=(const ImplicitPicPolarization&) = delete;

    DepositStats deposit(const std::vector<ParticleType>& species,
        const std::vector<double>& macroWeight, const std::vector<CouplingWeightsCache>& cachesBar,
        int depositFillOmpMinParticles, PicStepParticleBalance& chargeBalance);
    void assembleCells(const std::vector<ParticleType>& species, const NormalizationContext& norm, double dt_si);
    NeutralityCorrection removeMeanCharge(const NormalizationContext& norm);
    double totalCharge() const;
    void assembleFaces(const FaceSegmentGraph& faceGraph);

    // Kept separate to retain the driver's AMR policy-field reset ordering.
    void resetDepositedFields();
    void resetCoefficients(std::size_t segmentCount);
    // Optional AMR reset: clear particle arrays, reinitialize deposit workspaces.
    void resetScratch();

    const CellField<double>& chargeAmount() const { return Qstar; }
    const CellField<double>& chargeDensity() const { return rhoStar; }
    const std::vector<CellField<double>>& numberAmounts() const { return Nstar; }
    const std::vector<double>& susceptibilitySI() const { return muCell_Fpm; }
    const std::vector<double>& cellCoefficients() const { return kappaByCell; }
    const std::vector<double>& faceCoefficients() const { return kappaBySeg; }

private:
    friend class ImplicitPicRefinement;
    friend class ImplicitPicFieldSolve;
    void requireCurrentSizes() const;
    const QuadtreeMesh& mesh;
    CellField<double> Qstar, rhoStar;
    std::vector<CellField<double>> Nstar;
    std::vector<double> muCell_Fpm, kappaByCell, kappaBySeg;
    std::vector<std::vector<double>> depositQpScratch, depositWpScratch;
    std::vector<DepositQK::Workspace> depositDualWorkspaces;
};
