#include "Modules/Interface/Compatibility/InterfaceVolumeCorrection.h"
#include "Modules/Interface/InterfaceCellMeasure.h"

#include <Eigen/Sparse>
#ifdef EIGEN_SUPERLU_SUPPORT
#include <Eigen/SuperLUSupport>
#endif
#include <cmath>
#include <algorithm>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <tuple>
#include <sstream>
#include <iomanip>

namespace InterfaceVolumeCorrection
{
namespace
{
void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(std::string("[InterfaceVolumeCorrection] ")+reason);
}

double represented(long double value, const char* reason)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result) && (value == 0.0L || result != 0.0), reason);
    return result;
}

struct Row
{
    std::map<int, long double> coefficients;
    long double norm = 0.0L, rhs = 0.0L;
};
struct LinearStep : LinearReport
{
    std::vector<double> increments;
};
// One original row-normalized minimum-norm/SuperLU kernel. Scalar-specific
// assembly owns column identities; scales map phi-equivalent unknowns to
// physical returned data. Residuals check those returned physical increments.
LinearStep solveRows(std::map<int,Row> byCell,
    const std::vector<long double>& negativeVolumeByCell,
    const std::vector<long double>& targetNegativeVolumeByCell, bool incrementInput,
    const std::vector<double>& scales, double linearTolerance,
    const std::vector<Row>& regularization = {})
{
    std::map<int, int> columnByDof; // Scalar-specific input column order is retained.
    for (auto &item : byCell)
    {
        const int cell = item.first;
        auto& row = item.second;
        require(std::isfinite(negativeVolumeByCell[cell])
            && std::isfinite(targetNegativeVolumeByCell[cell]), "nonfinite constrained cell volume");
        for (const auto& entry : row.coefficients)
            if (entry.second != 0.0L)
            {
                row.norm = std::hypot(row.norm, entry.second);
                columnByDof.emplace(entry.first, 0);
            }
        require(std::isfinite(row.norm) && row.norm > 0.0L, "zero or nonfinite free sensitivity row");
        const long double defect = incrementInput ? targetNegativeVolumeByCell[cell]
            : targetNegativeVolumeByCell[cell]-negativeVolumeByCell[cell];
        require(std::isfinite(defect), "nonfinite cell volume defect");
        row.rhs = defect/row.norm;
        require(std::isfinite(row.rhs) && (defect == 0.0L || row.rhs != 0.0L),
            "unrepresentable normalized volume defect");
    }
    for (const auto& row:regularization) for (const auto& entry:row.coefficients)
        if (entry.second!=0) columnByDof.emplace(entry.first,0);
    require(byCell.size() <= columnByDof.size(), "more volume constraints than free nodal DOFs: rank deficient");
    require(byCell.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
        "volume row count exceeds sparse index range");

    LinearStep result;
    result.rows = static_cast<int>(byCell.size());
    result.columns = static_cast<int>(columnByDof.size());
    for (const auto& item : byCell)
        if (result.largestScaledDefectCell < 0 || std::abs(item.second.rhs) > result.largestScaledDefect)
        {
            result.largestScaledDefectCell = item.first;
            result.largestScaledDefect = std::abs(item.second.rhs);
            result.correspondingRowNorm = item.second.norm;
        }
    result.increments.assign(scales.size(),0.0);
    int column = 0;
    for (auto& item : columnByDof) item.second = column++;
    std::vector<Eigen::Triplet<double>> entries;
    Eigen::VectorXd rhs(result.rows);
    int rowIndex = 0;
    for (const auto& item : byCell)
    {
        const auto& row = item.second;
        rhs[rowIndex] = represented(row.rhs, "normalized volume defect cannot be represented by SuperLU");
        for (const auto& entry : row.coefficients)
            if (entry.second != 0.0L)
            {
                const long double normalized = entry.second/row.norm;
                require(normalized != 0.0L, "nonzero sensitivity underflowed during row normalization");
                entries.emplace_back(rowIndex, columnByDof.at(entry.first),
                    represented(normalized, "normalized sensitivity cannot be represented by SuperLU"));
            }
        ++rowIndex;
    }
    Eigen::SparseMatrix<double> matrix(result.rows, result.columns);
    matrix.setFromTriplets(entries.begin(), entries.end());
    matrix.makeCompressed();
#ifndef EIGEN_SUPERLU_SUPPORT
    throw std::runtime_error("[InterfaceVolumeCorrection] SuperLU is required; SparseLU fallback is disabled");
#else
    Eigen::VectorXd correction;
    {
        // Equality-constrained quadratic correction. The constraints remain A;
        // an empty penalty gives the SAME Euclidean minimum-norm problem.
        // No Gram squaring, diagonal pin, or relaxed volume equation.
        std::vector<Eigen::Triplet<double>> penalties;
        int r=0;
        for (const auto& row:regularization) {
            for (const auto& entry:row.coefficients) if (entry.second!=0)
                penalties.emplace_back(r,columnByDof.at(entry.first),
                    represented(entry.second,"unrepresentable regularization coefficient"));
            ++r;
        }
        Eigen::SparseMatrix<double> penalty(r,result.columns);
        penalty.setFromTriplets(penalties.begin(),penalties.end());
        Eigen::SparseMatrix<double> hessian=penalty.transpose()*penalty;
        for (int c=0;c<result.columns;++c) hessian.coeffRef(c,c)+=1;
        double objectiveScale=1;
        for (int c=0;c<result.columns;++c) objectiveScale=std::max(objectiveScale,hessian.coeff(c,c));
        require(std::isfinite(objectiveScale),"nonfinite regularization Hessian");
        std::vector<Eigen::Triplet<double>> blocks;
        for (int k=0;k<hessian.outerSize();++k)
            for (Eigen::SparseMatrix<double>::InnerIterator entry(hessian,k);entry;++entry) {
                require(std::isfinite(entry.value()),"nonfinite regularization Hessian coefficient");
                blocks.emplace_back(static_cast<int>(entry.row()),static_cast<int>(entry.col()),entry.value()/objectiveScale);
            }
        for (int k=0;k<matrix.outerSize();++k)
            for (Eigen::SparseMatrix<double>::InnerIterator entry(matrix,k);entry;++entry) {
                blocks.emplace_back(static_cast<int>(entry.col()),result.columns+static_cast<int>(entry.row()),entry.value());
                blocks.emplace_back(result.columns+static_cast<int>(entry.row()),static_cast<int>(entry.col()),entry.value());
            }
        const int size=result.columns+result.rows;
        Eigen::SparseMatrix<double> system(size,size);
        system.setFromTriplets(blocks.begin(),blocks.end()); system.makeCompressed();
        Eigen::VectorXd systemRhs=Eigen::VectorXd::Zero(size);
        systemRhs.tail(result.rows)=rhs;
        Eigen::SuperLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(system);
        require(solver.info()==Eigen::Success,"SuperLU constrained regularization factorization failed");
        Eigen::VectorXd solution=solver.solve(systemRhs);
        require(solver.info()==Eigen::Success && solution.allFinite(),"SuperLU constrained regularization solve failed");
        Eigen::VectorXd error=system*solution-systemRhs;
        const double rhsNorm=systemRhs.stableNorm();
        // The two KKT blocks have different natural scales. The stationarity
        // RHS is zero; normalizing it by the small volume defect is meaningless.
        // Check stationarity relative to its two objective/dual terms AND all
        // constraints relative to b. The physical returned-increment residual
        // is checked independently below, with its original tolerance.
        const auto relativeResidual=[&](const Eigen::VectorXd& x,const Eigen::VectorXd& e) {
            const double objectiveNorm=(hessian*x.head(result.columns)).stableNorm()/objectiveScale;
            const double dualNorm=(matrix.transpose()*x.tail(result.rows)).stableNorm();
            const double stationarityScale=std::max(objectiveNorm,dualNorm);
            const double stationarityError=e.head(result.columns).stableNorm();
            const double constraintError=e.tail(result.rows).stableNorm();
            require(std::isfinite(stationarityScale) && std::isfinite(stationarityError)
                && std::isfinite(constraintError),"nonfinite KKT block residual");
            return std::max(stationarityScale>0?stationarityError/stationarityScale:stationarityError,
                rhsNorm>0?constraintError/rhsNorm:constraintError);
        };
        double errorNorm=relativeResidual(solution,error);
        // Reuse the SAME factorization to refine the SAME equality system.
        // Retain only improving iterates; no diagonal shift or tolerance change.
        for (int refinement=0;refinement<4 && errorNorm>linearTolerance;++refinement) {
            const Eigen::VectorXd increment=solver.solve(-error);
            require(solver.info()==Eigen::Success && increment.allFinite(),"SuperLU KKT refinement solve failed");
            Eigen::VectorXd candidate=solution+increment;
            Eigen::VectorXd candidateError=system*candidate-systemRhs;
            if (!candidate.allFinite() || !candidateError.allFinite()) break;
            const double candidateNorm=relativeResidual(candidate,candidateError);
            if (!(candidateNorm<errorNorm)) break;
            solution=std::move(candidate);error=std::move(candidateError);errorNorm=candidateNorm;
        }
        require(error.allFinite() && std::isfinite(errorNorm) && std::isfinite(rhsNorm),"nonfinite complete KKT residual");
        result.normalEquationResidual=errorNorm;
        if (!std::isfinite(result.normalEquationResidual) || result.normalEquationResidual>linearTolerance) {
            std::ostringstream reason;
            reason<<std::setprecision(17)<<"[InterfaceVolumeCorrection] complete normalized KKT residual exceeds tolerance"
                <<"; rows="<<result.rows<<", columns="<<result.columns<<", residual="<<result.normalEquationResidual
                <<", tolerance="<<linearTolerance<<", rhsNorm="<<rhsNorm
                <<", correctionNorm="<<solution.head(result.columns).stableNorm()
                <<", multiplierNorm="<<solution.tail(result.rows).stableNorm()
                <<", largest scaled defect cell="<<result.largestScaledDefectCell;
            throw std::runtime_error(reason.str());
        }
        correction=solution.head(result.columns);
    }
    require(correction.allFinite(), "nonfinite constrained nodal correction");
    for (const auto& item : columnByDof)
        result.increments[item.first] = scales[item.first]==1.0 ? correction[item.second]
            : represented(static_cast<long double>(correction[item.second])*scales[item.first],
                "physical data correction cannot be represented");

    // Check the ORIGINAL assembled sensitivities after row normalization,
    // including every constraint; do not accept only the Gram solve residual.
    long double residualNorm = 0.0L, originalRhsNorm = 0.0L;
    for (const auto& item : byCell)
    {
        const auto& row = item.second;
        long double action = 0.0L;
        for (const auto& entry : row.coefficients)
        {
            // Check the returned physical gradient increment, including the
            // terminal scale multiplication, not only the raw SuperLU vector.
            const long double amplitude = scales[entry.first]==1.0 ? result.increments[entry.first]
                : static_cast<long double>(result.increments[entry.first])/scales[entry.first];
            action += (entry.second/row.norm)*amplitude;
        }
        require(std::isfinite(action), "nonfinite full normalized sensitivity action");
        residualNorm = std::hypot(residualNorm, action-row.rhs);
        originalRhsNorm = std::hypot(originalRhsNorm, row.rhs);
    }
    require(std::isfinite(residualNorm) && std::isfinite(originalRhsNorm), "nonfinite full normalized residual norm");
    result.linearResidual = represented(originalRhsNorm > 0.0L ? residualNorm/originalRhsNorm : residualNorm,
        "unrepresentable full normalized sensitivity residual");
    require(result.linearResidual <= linearTolerance, "complete normalized sensitivity residual exceeds tolerance");
    result.usedSuperLU = true;
    return result;
#endif
}

Result solveImpl(const InterfaceCellMeasure::Measure &measure,
                 const std::vector<long double> &targetNegativeVolumeByCell,
                 const std::vector<std::uint8_t> &fixedValue, double linearTolerance, bool incrementInput,
                 const std::vector<GradientDof> &freeGradients,
                 const InterfacePhaseFlux::FluxSensitivity *fluxSensitivity = nullptr,
                 long double halfDt = 0)
{
    require(std::isfinite(linearTolerance) && linearTolerance > 0.0,
            "linear tolerance must be positive and finite");
    require(!measure.arcs.empty() && !fixedValue.empty(),
            "free arcs and native vertex mask must be nonempty");
    require(targetNegativeVolumeByCell.size() == measure.negativeVolumeByCell.size(),
            "target and measured native-cell counts disagree");
    require(fixedValue.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
            "native vertex count exceeds sparse index range");
    for (const auto fixed : fixedValue)
        require(fixed <= 1, "fixed-value mask must be binary");
    require(freeGradients.size() <=
                static_cast<std::size_t>(std::numeric_limits<int>::max()) - fixedValue.size(),
            "appended gradient count exceeds sparse index range");
    const int nativeValues = static_cast<int>(fixedValue.size());
    std::map<std::pair<int, unsigned>, std::size_t> gradientByIdentity;
    for (std::size_t index = 0; index < freeGradients.size(); ++index)
    {
        const auto &dof = freeGradients[index];
        require(dof.vertexIndex >= 0 && dof.vertexIndex < nativeValues &&
                    (dof.component == 1 || dof.component == 2) && std::isfinite(dof.scale) && dof.scale > 0.0,
                "invalid appended first-derivative DOF");
        require(gradientByIdentity.emplace(std::make_pair(dof.vertexIndex, dof.component), index).second,
                "duplicate appended first-derivative identity");
    }

    std::map<int, Row> byCell;
    const auto assembleVertex = [&](Row &row, int vertex, const long double *sensitivity) {
        require(vertex >= 0 && static_cast<std::size_t>(vertex) < fixedValue.size(),
                "arc source vertex is outside the native mask");
        const long double coefficient = sensitivity[0];
        require(std::isfinite(coefficient), "nonfinite nodal value sensitivity");
        if (!fixedValue[vertex] && coefficient != 0.0L)
        {
            auto &accumulated = row.coefficients[vertex];
            accumulated += coefficient;
            require(std::isfinite(accumulated), "nonfinite aggregated nodal sensitivity");
        }
        for (unsigned component = 1; component <= 2; ++component)
        {
            const auto found = gradientByIdentity.find(std::make_pair(vertex, component));
            if (found == gradientByIdentity.end())
                continue;
            const auto index = found->second;
            const long double physical = sensitivity[component];
            const long double scaled = physical * freeGradients[index].scale;
            require(std::isfinite(physical) && std::isfinite(scaled) && (physical == 0.0L || scaled != 0.0L),
                    "nonfinite or underflowed appended gradient sensitivity");
            if (scaled == 0.0L)
                continue;
            auto &accumulated = row.coefficients[nativeValues + static_cast<int>(index)];
            accumulated += scaled;
            require(std::isfinite(accumulated), "nonfinite aggregated gradient sensitivity");
        }
    };
    for (const auto &arc : measure.arcs)
    {
        require(arc.cellIndex >= 0 &&
                    static_cast<std::size_t>(arc.cellIndex) < measure.negativeVolumeByCell.size(),
                "arc has invalid native control-volume owner");
        auto &row = byCell[arc.cellIndex];
        if (arc.nativeSourceCell >= 0)
        {
            require(static_cast<std::size_t>(arc.nativeSourceCell) < measure.negativeVolumeByCell.size() &&
                        !arc.independentSensitivity.empty(),
                    "native arc lacks its independent scalar sensitivity");
            int previous = -1;
            for (const auto &entry : arc.independentSensitivity)
            {
                require(entry.vertexIndex > previous,
                        "native sensitivity vertices must be unique and sorted");
                for (auto value : entry.coefficients)
                    require(std::isfinite(value), "nonfinite independent native sensitivity");
                assembleVertex(row, entry.vertexIndex, entry.coefficients.data());
                previous = entry.vertexIndex;
            }
        }
        else
        {
            require(arc.nativeSourceCell == -1 && arc.independentSensitivity.empty(),
                    "mixed uniform/native arc sensitivity addressing");
            for (int corner = 0; corner < 4; ++corner)
            {
                const int vertex = arc.vertexIndices[corner];
                for (int previous = 0; previous < corner; ++previous)
                    require(vertex != arc.vertexIndices[previous],
                            "arc source cell has repeated corner vertices");
                assembleVertex(row, vertex, arc.sensitivity.data() + 4 * corner);
            }
        }
    }

    if (fluxSensitivity)
    {
        require(std::isfinite(halfDt) && halfDt > 0, "invalid flux-variation time weight");
        require(fluxSensitivity->outwardFluxByCell.size() == measure.negativeVolumeByCell.size(),
                "phase flux variation has a different native-cell count");
        for (std::size_t cell = 0; cell < fluxSensitivity->outwardFluxByCell.size(); ++cell)
        {
            int previous = -1;
            for (const auto &entry : fluxSensitivity->outwardFluxByCell[cell])
            {
                const int vertex = entry.vertexIndex;
                require(vertex > previous && vertex < nativeValues,
                        "flux sensitivity vertices must be valid, unique and sorted");
                previous = vertex;
                std::array<long double, 4> weighted{};
                for (int k = 0; k < 4; ++k)
                {
                    weighted[k] = halfDt * entry.coefficients[k];
                    require(std::isfinite(entry.coefficients[k]) && std::isfinite(weighted[k]) &&
                                (entry.coefficients[k] == 0 || weighted[k] != 0),
                            "nonfinite or underflowed time-weighted flux variation");
                }
                bool contributes = !fixedValue[vertex] && weighted[0] != 0;
                for (unsigned k = 1; k <= 2; ++k)
                    contributes =
                        contributes || (weighted[k] != 0 && gradientByIdentity.count({vertex, k}) != 0);
                // A pure-cell flux row can enter the balance system even if
                // it has no free volume arc. Fixed-only variations are not
                // unknowns; the nonlinear gates still check every such cell.
                if (contributes)
                    assembleVertex(byCell[static_cast<int>(cell)], vertex, weighted.data());
            }
        }
    }

    std::vector<double> scales(fixedValue.size()+freeGradients.size(),1.0);
    for(std::size_t i=0;i<freeGradients.size();++i) scales[fixedValue.size()+i]=freeGradients[i].scale;
    auto step=solveRows(std::move(byCell),measure.negativeVolumeByCell,targetNegativeVolumeByCell,
        incrementInput,scales,linearTolerance);
    Result result;
    static_cast<LinearReport&>(result)=step;
    result.deltaValue.assign(step.increments.begin(),step.increments.begin()+nativeValues);
    result.deltaGradient.assign(step.increments.begin()+nativeValues,step.increments.end());
    return result;
}
}

Result solve(const InterfaceCellMeasure::Measure& measure,
    const std::vector<long double>& targetNegativeVolumeByCell,
    const std::vector<std::uint8_t>& fixedValue, double linearTolerance)
{
    return solveImpl(measure,targetNegativeVolumeByCell,fixedValue,linearTolerance,false,{});
}

Result solveDefect(const InterfaceCellMeasure::Measure& measure,
    const std::vector<long double>& deltaVolumeByCell,
    const std::vector<std::uint8_t>& fixedValue, double linearTolerance,
    const std::vector<GradientDof>& freeGradients)
{
    return solveImpl(measure,deltaVolumeByCell,fixedValue,linearTolerance,true,freeGradients);
}

Result solveBalanceDefect(const InterfaceCellMeasure::Measure &measure,
                          const InterfacePhaseFlux::FluxSensitivity &fluxSensitivity, long double halfDt,
                          const std::vector<long double> &negativeResidualByCell,
                          const std::vector<std::uint8_t> &fixedValue, double linearTolerance,
                          const std::vector<GradientDof> &freeGradients)
{
    return solveImpl(measure, negativeResidualByCell, fixedValue, linearTolerance, true, freeGradients,
                     &fluxSensitivity, halfDt);
}

QuarticResult solveQuarticBalanceDefect(const InterfaceCellMeasure::QuarticMeasure& measure,
    const InterfacePhaseFlux::QuarticFluxSensitivity& fluxSensitivity,long double halfDt,
    const std::vector<long double>& negativeResidualByCell,
    const std::vector<QuarticDof>& freeData,double linearTolerance,
    const QuarticRegularization& regularization)
{
    using Sensitivity=Quantum::Geometry::ActiveQuarticHermiteSensitivity;
    using Key=std::tuple<int,int,int>;
    require(std::isfinite(linearTolerance) && linearTolerance>0,
        "linear tolerance must be positive and finite");
    require(!measure.arcs.empty() && !freeData.empty(),"quartic free arcs and free data must be nonempty");
    require(negativeResidualByCell.size()==measure.negativeVolumeByCell.size()
        && fluxSensitivity.outwardFluxByCell.size()==measure.negativeVolumeByCell.size(),
        "quartic volume/flux/residual native-cell counts disagree");
    require(std::isfinite(halfDt) && halfDt>0,"invalid quartic flux-variation time weight");
    require(freeData.size()<=static_cast<std::size_t>(std::numeric_limits<int>::max()),
        "quartic data count exceeds sparse index range");
    const auto key=[](const auto& data) {
        const bool valid=data.kind==Sensitivity::Kind::Vertex ? data.component>=0 && data.component<4
            : data.kind==Sensitivity::Kind::EdgeMidpoint ? data.component>=0 && data.component<2
            : data.kind==Sensitivity::Kind::CellCenter && data.component==0;
        require(data.index>=0 && valid,"invalid independent quartic data identity");
        return Key{static_cast<int>(data.kind),data.index,data.component};
    };
    std::map<Key,int> freeByIdentity;
    std::vector<double> scales;
    for(std::size_t i=0;i<freeData.size();++i)
    {
        require(std::isfinite(freeData[i].scale) && freeData[i].scale>0,"invalid quartic correction scale");
        require(freeByIdentity.emplace(key(freeData[i]),static_cast<int>(i)).second,
            "duplicate free quartic data identity");
        scales.push_back(freeData[i].scale);
    }
    std::map<int,Row> byCell;
    const auto assemble=[&](int cell,const std::vector<Sensitivity>& entries,long double factor,bool volumeArc) {
        require(cell>=0 && static_cast<std::size_t>(cell)<measure.negativeVolumeByCell.size(),
            "invalid quartic control-volume owner");
        if(volumeArc) byCell[cell]; // Retain every free-arc row, including rank failures.
        Key previous{-1,-1,-1};
        for(const auto& entry:entries)
        {
            const Key identity=key(entry);
            require(previous<identity,"quartic sensitivity identities must be unique and sorted");
            previous=identity;
            const long double weighted=factor*entry.coefficient;
            require(std::isfinite(entry.coefficient) && std::isfinite(weighted)
                && (entry.coefficient==0 || weighted!=0),"nonfinite or underflowed quartic sensitivity");
            const auto found=freeByIdentity.find(identity);
            if(found==freeByIdentity.end() || weighted==0) continue;
            const int column=found->second;
            const long double scaled=weighted*scales[column];
            require(std::isfinite(scaled) && scaled!=0,"unrepresentable scaled quartic sensitivity");
            auto& accumulated=byCell[cell].coefficients[column];
            accumulated+=scaled;
            require(std::isfinite(accumulated),"nonfinite aggregated quartic sensitivity");
        }
    };
    for(const auto& arc:measure.arcs)
    {
        require(arc.nativeSourceCell>=0
            && static_cast<std::size_t>(arc.nativeSourceCell)<measure.negativeVolumeByCell.size()
            && !arc.independentSensitivity.empty(),"quartic arc lacks its native independent sensitivity");
        assemble(arc.cellIndex,arc.independentSensitivity,1.0L,true);
    }
    for(std::size_t cell=0;cell<fluxSensitivity.outwardFluxByCell.size();++cell)
        assemble(static_cast<int>(cell),fluxSensitivity.outwardFluxByCell[cell],halfDt,false);
    std::vector<Row> penalty;
    for (const auto& entries:regularization) {
        Row row;
        for (const auto& entry:entries) {
            require(std::isfinite(entry.coefficient),"nonfinite geometric regularization coefficient");
            const auto found=freeByIdentity.find(key(entry));
            if (found!=freeByIdentity.end()) row.coefficients[found->second]+=entry.coefficient*scales[found->second];
        }
        if (!row.coefficients.empty()) penalty.push_back(std::move(row));
    }
    auto step=solveRows(std::move(byCell),measure.negativeVolumeByCell,negativeResidualByCell,
        true,scales,linearTolerance,penalty);
    QuarticResult result;
    static_cast<LinearReport&>(result)=step;
    result.deltaData=std::move(step.increments);
    return result;
}
} // namespace InterfaceVolumeCorrection
