#pragma once

#include <Eigen/Sparse>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

// One assembled pressure/projection equation. The selected row strategy owns
// its gauge, signs and scaling; this object adds no normalization or solver.
// Advance consumes the completed system before publishing the pressure potential.
struct HydroPressureSystem
{
    Eigen::SparseMatrix<double> matrix;
    Eigen::VectorXd rhs;
};

// One complete-domain boundary projection. Native XY/RZ preparation supplies
// face metrics in velocity order; equations, correction and flux checks then
// consume the same owned incidence/constraints, without rereading live geometry.
// Initial compatibility uses mass scale 1; a time step supplies its original dt.
class HydroBoundaryProjection
{
public:
    struct Face
    {
        int high = -1, low = -1; // low < 0: high is the boundary owner
        double area = 0., inverseMass = 0., sign = 0., fixedVelocity = 0.;
        bool fixed = false;
    };
    struct Flux { Eigen::VectorXd integrated, absolute; };
    template<class ReadVelocity>
    HydroBoundaryProjection(int cells, std::vector<Face> faces, ReadVelocity velocity)
        : faces_(std::move(faces))
    {
        if (cells <= 0 || faces_.empty()) throw std::invalid_argument("empty boundary projection");
        std::vector<Eigen::Triplet<double>> entries;
        equations_.rhs = Eigen::VectorXd::Zero(cells);
        auto& rhs = equations_.rhs;
        for (std::size_t row = 0; row < faces_.size(); ++row)
        {
            const auto& f = faces_[row];
            if (f.high < 0 || f.high >= cells || f.low < -1 || f.low >= cells || f.low == f.high
                || (f.low >= 0 && f.fixed)) throw std::invalid_argument("invalid boundary projection incidence");
            const double predicted = velocity(row);
            if (f.low < 0)
            {
                rhs[f.high] -= f.sign * f.area * predicted;
                if (f.fixed) continue;
                entries.emplace_back(f.high, f.high, f.area * f.area * f.inverseMass);
            }
            else
            {
                const double coefficient = f.area * f.area * f.inverseMass;
                rhs[f.high] += f.area * predicted; rhs[f.low] -= f.area * predicted;
                entries.emplace_back(f.high, f.high, coefficient); entries.emplace_back(f.low, f.low, coefficient);
                entries.emplace_back(f.high, f.low, -coefficient); entries.emplace_back(f.low, f.high, -coefficient);
            }
        }
        equations_.matrix.resize(cells, cells);
        equations_.matrix.setFromTriplets(entries.begin(), entries.end());
    }
    HydroBoundaryProjection(const HydroBoundaryProjection&) = delete;
    HydroBoundaryProjection(HydroBoundaryProjection&&) = default;
    const HydroPressureSystem& equations() const { return equations_; }
    const Face& face(std::size_t row) const { return faces_.at(row); }

    // Initial fields already contain prescribed values; preserve those entries.
    // Step correction imposes them even if the momentum predictor differs.
    // Read/write may address the same field: each face only reads its own value.
    template<class ReadVelocity, class WriteVelocity>
    void apply(const Eigen::VectorXd& multiplier, ReadVelocity velocity, WriteVelocity write,
        bool preserveFixed, const char* nonfiniteReason) const
    {
        if (multiplier.size() != equations_.rhs.size()) throw std::invalid_argument("boundary multiplier size mismatch");
        for (std::size_t row = 0; row < faces_.size(); ++row)
        {
            const auto& f = faces_[row];
            if (preserveFixed && f.fixed) continue;
            const double load = f.low < 0 ? -f.sign * f.area * multiplier[f.high]
                : f.area * (multiplier[f.high] - multiplier[f.low]);
            const double value = f.fixed ? f.fixedVelocity : velocity(row) - f.inverseMass * load;
            if (!std::isfinite(value)) throw std::runtime_error(nonfiniteReason);
            write(row, value);
        }
    }
    template<class ReadVelocity>
    Flux fluxes(ReadVelocity velocity) const
    {
        Flux result{Eigen::VectorXd::Zero(equations_.rhs.size()), Eigen::VectorXd::Zero(equations_.rhs.size())};
        for (std::size_t row = 0; row < faces_.size(); ++row)
        {
            const auto& f = faces_[row];
            const double q = f.area * velocity(row);
            if (f.low < 0)
            {
                result.integrated[f.high] += f.sign * q;
                result.absolute[f.high] += std::abs(q);
            }
            else
            {
                result.integrated[f.high] -= q; result.integrated[f.low] += q;
                result.absolute[f.high] += std::abs(q); result.absolute[f.low] += std::abs(q);
            }
        }
        return result;
    }
private:
    std::vector<Face> faces_;
    HydroPressureSystem equations_;
};
