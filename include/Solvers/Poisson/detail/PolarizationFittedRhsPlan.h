#pragma once

// Private data of an already validated fitted operator. This is an owned
// native-index assembly plan, not another mesh and not a borrowed live graph.
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Core/Mesh/QuadtreeMesh.h"
#include "Solvers/Poisson/PolarizationPoissonFVSolverCellCentered.h"

namespace PolarizationFittedBoundaryDetail
{
class RhsPlan
{
    using Boundary = PolarizationPoissonFVSolverCellCentered::FittedBoundary;
    using Type = PolarizationPoissonFVSolverCellCentered::FittedBoundaryType;
    struct FaceUpdate
    {
        int row = -1;
        bool zeroFlux = false, dirichlet = false;
        double conductance = 0.0, boundaryValue = 0.0;
    };
    std::vector<int> labels_, cells_;
    std::vector<double> volumes_;
    std::vector<FaceUpdate> faces_;
    std::vector<std::pair<std::size_t, std::size_t>> reversePairs_;
    bool corrected_ = false;

    static void require(bool condition, const char* message)
    {
        if (!condition) throw std::runtime_error(std::string("prepared fitted RHS: ") + message);
    }

public:
    // Called only after the original full graph/coefficient validators. The
    // plan becomes const and reusable only after a successful full solve.
    RhsPlan(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
        const std::vector<Boundary>& boundaries, const std::vector<int>& cells,
        const std::vector<int>& rows, const std::vector<double>& volumes, bool corrected)
        : labels_(mesh.cellIndexToLabel, mesh.cellIndexToLabel + mesh.cellIndexNum()),
          cells_(cells), volumes_(volumes), faces_(graph.segments().size()), corrected_(corrected)
    {
        const auto& segments = graph.segments();
        const auto& offsets = graph.cellOffsets();
        for (std::size_t sid = 0; sid < segments.size(); ++sid)
        {
            const auto& face = segments[sid];
            faces_[sid].row = rows[mesh.cells(face.K).index()];
            require(faces_[sid].row >= 0, "missing active owner row");
            faces_[sid].zeroFlux = boundaries[sid].type == Type::HomogeneousNeumann;
            if (!corrected || boundaries[sid].type != Type::Interior) continue;
            const int neighbor = mesh.cells(face.L).index();
            std::size_t reverse = segments.size();
            for (std::size_t other = offsets[neighbor]; other < offsets[neighbor + 1]; ++other)
                if (segments[other].isHorizontal == face.isHorizontal &&
                    segments[other].edgeLabel == face.edgeLabel)
                {
                    require(reverse == segments.size(), "duplicate reverse face");
                    reverse = other;
                }
            require(reverse < segments.size(), "missing reverse face");
            if (sid < reverse) reversePairs_.emplace_back(sid, reverse);
        }
    }

    // Capture the actual conductance where the shared assembler computes it;
    // do not re-derive it from rounded geometry or pre-sum boundary RHS terms.
    void recordDirichlet(std::size_t sid, double conductance, double value)
    {
        auto& face = faces_[sid];
        face.dirichlet = true;
        face.conductance = conductance;
        face.boundaryValue = value;
    }

    Eigen::VectorXd assemble(const CellField<double>& rhs, const CellField<double>& phi,
        const std::vector<double>* constants) const
    {
        require(rhs.size() == static_cast<int>(labels_.size()) && rhs.value &&
                phi.size() == rhs.size() && phi.value, "field size/storage does not match the prepared operator");
        require(&rhs != &phi && rhs.value.get() != phi.value.get(), "rhs/phi aliasing is unsupported");
        require(corrected_ == (constants != nullptr), "correction mode does not match the prepared operator");
        for (int ci = 0; ci < rhs.size(); ++ci)
        {
            require(rhs.label(ci) == labels_[ci] && phi.label(ci) == labels_[ci],
                "field labels do not match the prepared native identity");
            require(std::isfinite(rhs.value[ci]) && std::isfinite(phi.value[ci]),
                "non-finite rhs or initial phi");
        }
        if (constants)
        {
            require(constants->size() == faces_.size(), "correction constant count does not cover every face");
            for (std::size_t sid = 0; sid < faces_.size(); ++sid)
            {
                require(std::isfinite((*constants)[sid]), "non-finite correction constant");
                require(!faces_[sid].zeroFlux || (*constants)[sid] == 0.0,
                    "homogeneous Neumann correction constant must be exactly zero");
            }
            for (const auto& pair : reversePairs_)
                require((*constants)[pair.second] == -(*constants)[pair.first],
                    "interior correction constants must be exactly antisymmetric");
        }
        Eigen::VectorXd b(cells_.size());
        for (std::size_t row = 0; row < cells_.size(); ++row)
            b[row] = rhs.value[cells_[row]] * volumes_[row];
        // Keep the original per-segment operation order. Combining all static
        // boundary terms first would change floating-point charge histories.
        for (std::size_t sid = 0; sid < faces_.size(); ++sid)
        {
            const auto& face = faces_[sid];
            if (constants) b[face.row] -= (*constants)[sid];
            if (face.dirichlet) b[face.row] += face.conductance * face.boundaryValue;
        }
        require(b.allFinite(), "assembled RHS contains non-finite entries");
        return b;
    }
};
}
