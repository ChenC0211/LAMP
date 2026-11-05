#include "Core/Geometry/FittedDomainGeometry.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Field/CellField.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <limits>

namespace FittedDomainGeometry
{
namespace
{
// Exact byte signature, not a collision-prone hash. The owned sampler's
// source geometry/CSG must match the definition later used for face binding.
// Boundary associations deliberately remain independently checked by bindFace.
struct DefinitionKey
{
    std::string value;
    template<class T> void add(T item)
    { value.append(reinterpret_cast<const char*>(&item),sizeof(item)); }
    void add(const std::string& text)
    { add(text.size()); value.append(text); }
    void point(Vector2D p) { add(p.x); add(p.y); }
    void references(const std::vector<SignInt>& refs)
    {
        add(refs.size());
        for(auto ref:refs){add(ref.value);add(static_cast<int>(ref.isNegative));}
    }
};
std::string initialDefinitionKey(const DictionaryGeometry& definition)
{
    DefinitionKey key;
    key.add(static_cast<int>(definition.hasDomain));
    if(definition.hasDomain){key.point(definition.domainMin);key.point(definition.domainMax);}
    key.add(definition.lines.size());
    for(const auto& line:definition.lines)
    {
        key.add(static_cast<int>(line.type()));key.point(line.points[0]);key.point(line.points[1]);
        if(line.type()==LineType::arc)
        {key.point(line.arcCenter());key.add(line.radius());key.add(line.radian());}
        if(line.type()==LineType::spline)
        {const auto points=line.featurePoints();key.add(points.size());for(auto p:points)key.point(p);}
    }
    key.add(definition.levelSetName.size());for(const auto& name:definition.levelSetName)key.add(name);
    key.add(definition.innerPoints.size());for(auto p:definition.innerPoints)key.point(p);
    key.add(definition.moving.size());for(bool moving:definition.moving)key.add(static_cast<int>(moving));
    key.add(definition.linesOnInterface.size());for(const auto& refs:definition.linesOnInterface)key.references(refs);
    key.add(definition.signClosureLines.size());for(const auto& refs:definition.signClosureLines)key.references(refs);
    key.add(definition.phaseName.size());for(const auto& name:definition.phaseName)key.add(name);
    key.add(definition.levelSetOperation.size());
    for(const auto& operations:definition.levelSetOperation)
    {key.add(operations.size());for(auto operation:operations)key.add(static_cast<int>(operation));}
    key.add(definition.levelSetIndexA.size());for(const auto& refs:definition.levelSetIndexA)key.references(refs);
    key.add(definition.levelSetIndexB.size());for(const auto& refs:definition.levelSetIndexB)key.references(refs);
    key.add(definition.boundaryLines.size());for(const auto& refs:definition.boundaryLines)key.references(refs);
    return std::move(key.value);
}
struct Expression
{
    std::vector<LevelSetOperation> operations;
    std::vector<SignInt> a, b;
};
struct Bounds
{
    Vector2D lower, upper;
    bool intersectsOpen(Vector2D lo, Vector2D hi) const
    {
        return upper.x > lo.x && lower.x < hi.x && upper.y > lo.y && lower.y < hi.y;
    }
};
Bounds bounds(Vector2D a, Vector2D b)
{
    ExactCurveDistance::requireFinitePoint(a);
    ExactCurveDistance::requireFinitePoint(b);
    return {{std::min(a.x,b.x),std::min(a.y,b.y)},
            {std::max(a.x,b.x),std::max(a.y,b.y)}};
}
void include(Bounds& box, Vector2D p)
{
    ExactCurveDistance::requireFinitePoint(p);
    box.lower.x=std::min(box.lower.x,p.x); box.lower.y=std::min(box.lower.y,p.y);
    box.upper.x=std::max(box.upper.x,p.x); box.upper.y=std::max(box.upper.y,p.y);
}
struct CurveEnvelope { Bounds box; Vector2D start,end; };
CurveEnvelope envelope(const Line& line, SignInt reference)
{
    CurveEnvelope result{bounds(line.points[0],line.points[1]),line.points[0],line.points[1]};
    if (line.type()==LineType::arc)
    {
        const ExactCurveDistance::CircularArc arc(line.arcCenter(),line.points[0],
            line.radius(),line.radian());
        result={bounds(arc.start(),arc.end()),arc.start(),arc.end()};
        const double angle=std::atan2(line.points[0].y-line.arcCenter().y,
            line.points[0].x-line.arcCenter().x);
        const double sweep=line.radian()*ExactCurveDistance::pi;
        if(sweep!=0.0) for(int k=-6;k<=6;++k)
        {
            const double t=(k*0.5*ExactCurveDistance::pi-angle)/sweep;
            if(t>0.0 && t<1.0) include(result.box,arc.point(t));
        }
    }
    else if(line.type()==LineType::spline)
    {
        const auto controls=line.featurePoints();
        if(controls.size()<4) throw std::invalid_argument("initial CSG spline needs four control points");
        // The existing sampled cubic B-spline and its chords lie in this hull.
        for(auto p:controls) include(result.box,p);
        result.start=controls.front(); result.end=controls.back();
    }
    else if(line.type()!=LineType::straight)
        throw std::invalid_argument("unsupported initial CSG primitive");
    if(reference.isNegative) std::swap(result.start,result.end);
    return result;
}
struct SignRange { int lower,upper; };
SignRange signedRange(SignRange value, bool negative)
{
    return negative ? SignRange{-value.upper,-value.lower} : value;
}
SignRange evaluate(const Expression& expression, std::vector<SignRange> values)
{
    for(std::size_t i=0;i<expression.operations.size();++i)
    {
        const auto a=signedRange(values[expression.a[i].value],expression.a[i].isNegative);
        const auto b=signedRange(values[expression.b[i].value],expression.b[i].isNegative);
        const auto operation=expression.operations[i];
        values.push_back(operation==LevelSetOperation::non ? a : operation==LevelSetOperation::min
            ? SignRange{std::min(a.lower,b.lower),std::min(a.upper,b.upper)}
            : SignRange{std::max(a.lower,b.lower),std::max(a.upper,b.upper)});
    }
    return values.back();
}
double evaluate(const Expression& expression, std::vector<double> values)
{
    for(std::size_t i=0;i<expression.operations.size();++i)
    {
        const double a=(expression.a[i].isNegative ? -1.0:1.0)*values[expression.a[i].value];
        const double b=(expression.b[i].isNegative ? -1.0:1.0)*values[expression.b[i].value];
        if(!std::isfinite(a)||!std::isfinite(b)) throw std::runtime_error("non-finite initial CSG operand");
        values.push_back(expression.operations[i]==LevelSetOperation::non ? a
            : expression.operations[i]==LevelSetOperation::min ? std::min(a,b):std::max(a,b));
    }
    return values.back();
}
}
struct InitialPhaseSamplerData
{
    std::vector<std::function<ImmersedGeometry<QuadtreeMesh>::InitialCurveSample(Vector2D)>> primitives;
    std::vector<Expression> expressions;
    std::string definitionKey;
    std::vector<std::vector<Bounds>> signBounds;
    Vector2D lower,upper;
    std::vector<double> sample(Vector2D point) const
    {
        ExactCurveDistance::requireFinitePoint(point);
        std::vector<double> values;
        for(const auto& primitive:primitives)
        {
            const auto value=primitive(point);
            if(!std::isfinite(value.unsignedDistance)||value.unsignedDistance<0.0)
                throw std::runtime_error("non-finite initial CSG primitive");
            values.push_back(value.negative ? -value.unsignedDistance:value.unsignedDistance);
        }
        return values;
    }
};
InitialPhaseSampler::InitialPhaseSampler(QuadtreeMesh& mesh, const DictionaryGeometry& definition)
{
    const auto n=definition.levelSetName.size(), p=definition.phaseName.size();
    if(!n||!p||definition.innerPoints.size()!=n||definition.moving.size()!=n
        ||definition.linesOnInterface.size()!=n
        ||(!definition.signClosureLines.empty()&&definition.signClosureLines.size()!=n)
        ||definition.levelSetOperation.size()!=p||definition.levelSetIndexA.size()!=p
        ||definition.levelSetIndexB.size()!=p||definition.boundaryLines.size()!=p)
        throw std::invalid_argument("invalid initial CSG definition sizes");
    auto data=std::make_shared<InitialPhaseSamplerData>();
    data->lower=mesh.minCoord; data->upper=mesh.maxCoord;
    data->definitionKey=initialDefinitionKey(definition);
    for(std::size_t phase=0;phase<p;++phase)
    {
        Expression expression{definition.levelSetOperation[phase],definition.levelSetIndexA[phase],
            definition.levelSetIndexB[phase]};
        if(expression.operations.empty()||expression.a.size()!=expression.operations.size()
            ||expression.b.size()!=expression.operations.size())
            throw std::invalid_argument("invalid initial CSG expression size");
        for(std::size_t j=0;j<expression.operations.size();++j)
        {
            const auto a=expression.a[j], b=expression.b[j];
            const auto operation=expression.operations[j];
            if(a.value<0||b.value<0||static_cast<std::size_t>(a.value)>=n+j
                ||static_cast<std::size_t>(b.value)>=n+j
                ||(operation!=LevelSetOperation::non&&operation!=LevelSetOperation::min&&operation!=LevelSetOperation::max)
                ||(operation==LevelSetOperation::non&&(a.value!=b.value||a.isNegative!=b.isNegative)))
                throw std::invalid_argument("invalid initial CSG operation/reference");
        }
        data->expressions.push_back(std::move(expression));
    }
    data->signBounds.resize(n);
    for(std::size_t ls=0;ls<n;++ls)
    {
        auto references=definition.linesOnInterface[ls];
        if(!definition.signClosureLines.empty()) references.insert(references.end(),
            definition.signClosureLines[ls].begin(),definition.signClosureLines[ls].end());
        if(references.empty()) throw std::invalid_argument("empty initial CSG sign path");
        std::vector<CurveEnvelope> path;
        for(auto reference:references)
        {
            if(reference.value<0||reference.value>=static_cast<int>(definition.lines.size()))
                throw std::invalid_argument("invalid initial CSG zero-based line reference");
            path.push_back(envelope(definition.lines[reference.value],reference));
            data->signBounds[ls].push_back(path.back().box);
        }
        // initialCurveParity always includes these connectors, including the
        // legacy implicit closure. They are sign-only, never physical faces.
        for(std::size_t i=0;i<path.size();++i)
            data->signBounds[ls].push_back(bounds(path[i].end,path[(i+1)%path.size()].start));
    }
    DictionaryGeometry initialDefinition=definition;
    ImmersedGeometry<QuadtreeMesh> seed(initialDefinition,mesh);
    for(std::size_t ls=0;ls<n;++ls)
        data->primitives.push_back(seed.initialCurveSampler(mesh,definition.linesOnInterface[ls],static_cast<int>(ls)));
    data_=std::move(data);
}
int InitialPhaseSampler::phaseCount() const { return static_cast<int>(data_->expressions.size()); }
double InitialPhaseSampler::levelSetValue(int levelSetIndex, Vector2D point) const
{
    if(levelSetIndex<0||levelSetIndex>=static_cast<int>(data_->primitives.size()))
        throw std::invalid_argument("initial CSG level-set index outside geometry");
    ExactCurveDistance::requireFinitePoint(point);
    const auto value=data_->primitives[levelSetIndex](point);
    if(!std::isfinite(value.unsignedDistance)||value.unsignedDistance<0.0)
        throw std::runtime_error("non-finite initial CSG primitive");
    return value.negative ? -value.unsignedDistance:value.unsignedDistance;
}
double InitialPhaseSampler::phaseValue(int phaseId, Vector2D point) const
{
    if(phaseId<1||phaseId>phaseCount()) throw std::invalid_argument("initial CSG phase ID outside geometry");
    return evaluate(data_->expressions[phaseId-1],data_->sample(point));
}
int InitialPhaseSampler::phaseAt(Vector2D point) const
{
    const auto values=data_->sample(point);
    int phase=-1;
    for(std::size_t p=0;p<data_->expressions.size();++p)
        if(evaluate(data_->expressions[p],values)<0.0)
        {
            if(phase!=-1) throw std::runtime_error("overlapping initial CSG physical phases");
            phase=static_cast<int>(p)+1;
        }
    return phase;
}
FixedSolidSelection selectFixedSolids(const QuadtreeMesh& mesh, const DictionaryGeometry& definition,
    const InitialPhaseSampler& sampler, const std::vector<int>& excludedPhaseIds, double tolerance)
{
    if(!(tolerance>=0.0)||!std::isfinite(tolerance)||!definition.periodicPairs.empty()
        ||mesh.minCoord.x!=sampler.data_->lower.x||mesh.minCoord.y!=sampler.data_->lower.y
        ||mesh.maxCoord.x!=sampler.data_->upper.x||mesh.maxCoord.y!=sampler.data_->upper.y
        ||definition.phaseName.size()!=sampler.data_->expressions.size())
        throw std::invalid_argument("invalid fixed-solid domain/tolerance");
    if(initialDefinitionKey(definition)!=sampler.data_->definitionKey)
        throw std::invalid_argument("fixed-solid sampler source geometry/CSG has changed");
    std::vector<bool> selected(sampler.phaseCount()+1,false);
    for(int phase:excludedPhaseIds)
    {
        if(phase<1||phase>sampler.phaseCount()||selected[phase])
            throw std::invalid_argument("invalid or duplicate excluded solid phase ID");
        selected[phase]=true;
    }
    if(excludedPhaseIds.empty()) throw std::invalid_argument("fixed-solid selection requires explicit excluded phase IDs");
    FixedSolidSelection result;
    result.excludedPhaseByCell.assign(mesh.cellIndexNum(),0);
    std::vector<int> centerPhase(mesh.cellIndexNum(),-1);
    for(int i=0;i<mesh.cellIndexNum();++i)
    {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[i]);
        if(cell.index()!=i) throw std::runtime_error("invalid fixed-solid native leaf index");
        const auto center=cell.centerCoord(), size=cell.size();
        if(!(size.x>0.0)||!(size.y>0.0)||!std::isfinite(size.x)||!std::isfinite(size.y))
            throw std::runtime_error("invalid fixed-solid native leaf extent");
        const auto lo=mesh.verts(cell.cellToVert[0]).coord();
        const auto hi=mesh.verts(cell.cellToVert[2]).coord();
        ExactCurveDistance::requireFinitePoint(lo); ExactCurveDistance::requireFinitePoint(hi);
        if(!(lo.x<hi.x)||!(lo.y<hi.y)) throw std::runtime_error("invalid fixed-solid native vertex extent");
        const auto values=sampler.data_->sample(center);
        std::vector<SignRange> signs;
        for(std::size_t ls=0;ls<values.size();++ls)
        {
            bool possibleChange=false;
            for(const auto& box:sampler.data_->signBounds[ls]) possibleChange|=box.intersectsOpen(lo,hi);
            signs.push_back(possibleChange ? SignRange{-1,1}
                : values[ls]<0.0 ? SignRange{-1,-1}:values[ls]>0.0 ? SignRange{1,1}:SignRange{0,0});
        }
        for(int phase:excludedPhaseIds)
        {
            const auto sign=evaluate(sampler.data_->expressions[phase-1],signs);
            if(sign.lower<0 && sign.upper>=0)
                throw std::runtime_error("selected solid cuts a leaf or has unsupported whole-leaf geometry: cell="+std::to_string(cell.label()));
            if(sign.upper<0)
            {
                if(result.excludedPhaseByCell[i]) throw std::runtime_error("overlapping selected solid phases");
                result.excludedPhaseByCell[i]=phase;
            }
        }
        centerPhase[i]=sampler.phaseAt(center);
        if(result.excludedPhaseByCell[i]&&centerPhase[i]!=result.excludedPhaseByCell[i])
            throw std::runtime_error("fixed-solid phase certificate disagrees with initial CSG");
    }
    FaceSegmentGraph graph; graph.build(mesh);
    result.boundaryIdsBySegment.assign(graph.segments().size(),0);
    for(std::size_t i=0;i<graph.segments().size();++i)
    {
        const auto& segment=graph.segments()[i];
        const int owner=mesh.cells(segment.K).index();
        if(owner<0||owner>=mesh.cellIndexNum()) throw std::runtime_error("invalid fitted face owner leaf");
        if(result.excludedPhaseByCell[owner]) continue;
        const int other=segment.L<0 ? -1:mesh.cells(segment.L).index();
        if(segment.L>=0&&(other<0||other>=mesh.cellIndexNum()))
            throw std::runtime_error("invalid fitted face neighbor leaf");
        if(other>=0&&!result.excludedPhaseByCell[other]) continue;
        if(centerPhase[owner]<1) throw std::runtime_error("uncovered fitted boundary owner phase");
        const int excluded=other<0 ? -1:result.excludedPhaseByCell[other];
        const auto& boundary=bindFace(definition,segment,centerPhase[owner],excluded,tolerance);
        if(boundary.id<1) throw std::runtime_error("invalid fitted boundary ID");
        result.boundaryIdsBySegment[i]=boundary.id;
    }
    return result;
}

// Evaluate complete phase expressions on interpolated cell-center level sets;
// averaging vertex phase IDs would assign the curved liquid cap incorrectly.
std::vector<int> cellPhases(const QuadtreeMesh& mesh,
    const ImmersedGeometry<QuadtreeMesh>& geometry)
{
    std::vector<CellField<double>> cellLevelSets;
    for (const auto& field : geometry.levelSets)
        cellLevelSets.push_back(*field.toCellField(mesh));
    std::vector<int> result(mesh.cellIndexNum(), -1);
    const int originalCount = static_cast<int>(cellLevelSets.size());
    for (std::size_t p = 0; p < geometry.phases.size(); ++p)
    {
        const auto operations = geometry.phases[p].operation();
        const auto a = geometry.phases[p].levelSetIndexA();
        const auto b = geometry.phases[p].levelSetIndexB();
        if (operations.empty() || a.size() != operations.size() || b.size() != operations.size())
            throw std::runtime_error("invalid cell phase expression");
        for (std::size_t j = 0; j < operations.size(); ++j)
        {
            if (a[j].value < 0 || b[j].value < 0
                || a[j].value >= originalCount + static_cast<int>(j)
                || b[j].value >= originalCount + static_cast<int>(j)
                || (operations[j] != LevelSetOperation::non && operations[j] != LevelSetOperation::min
                    && operations[j] != LevelSetOperation::max)
                || (operations[j] == LevelSetOperation::non
                    && (a[j].value != b[j].value || a[j].isNegative != b[j].isNegative)))
                throw std::runtime_error("invalid cell phase operation/reference");
        }
        std::vector<double> values(originalCount + operations.size());
        for (int i = 0; i < mesh.cellIndexNum(); ++i)
        {
            for (int j = 0; j < originalCount; ++j) values[j] = cellLevelSets[j].value[i];
            for (std::size_t j = 0; j < operations.size(); ++j)
            {
                const double va = (a[j].isNegative ? -1.0 : 1.0) * values[a[j].value];
                const double vb = (b[j].isNegative ? -1.0 : 1.0) * values[b[j].value];
                if (!std::isfinite(va) || !std::isfinite(vb))
                    throw std::runtime_error("non-finite cell phase operand");
                values[originalCount + j] = operations[j] == LevelSetOperation::non ? va
                    : operations[j] == LevelSetOperation::min ? std::min(va, vb) : std::max(va, vb);
            }
            if (values.back() < 0.0)
            {
                if (result[i] != -1) throw std::runtime_error("overlapping physical phases at cell center");
                result[i] = static_cast<int>(p) + 1;
            }
        }
    }
    return result;
}

namespace
{
bool coversFace(const Line& line, const FaceSegment& segment, double tolerance)
{
    if (line.type() != LineType::straight) return false;
    const bool horizontal = segment.isHorizontal;
    const double normal0 = horizontal ? line.points[0].y : line.points[0].x;
    const double normal1 = horizontal ? line.points[1].y : line.points[1].x;
    const double faceNormal = horizontal ? segment.faceCenter.y : segment.faceCenter.x;
    const double tangent0 = horizontal ? line.points[0].x : line.points[0].y;
    const double tangent1 = horizontal ? line.points[1].x : line.points[1].y;
    const double faceTangent = horizontal ? segment.faceCenter.x : segment.faceCenter.y;
    return std::abs(normal0 - faceNormal) <= tolerance
        && std::abs(normal1 - faceNormal) <= tolerance
        && faceTangent - 0.5 * segment.Af >= std::min(tangent0, tangent1) - tolerance
        && faceTangent + 0.5 * segment.Af <= std::max(tangent0, tangent1) + tolerance;
}

}

const DictionaryGeometry::BoundaryEntity& bindFace(const DictionaryGeometry& definition,
    const FaceSegment& segment, int ownerPhase, int excludedPhase, double tolerance)
{
    const DictionaryGeometry::BoundaryEntity* found = nullptr;
    for (const auto& boundary : definition.boundaries)
    {
        if (boundary.ignored) continue;
        if (segment.L < 0)
        {
            if (boundary.kind != GeometricBoundary::domain || boundary.adjacentPhase + 1 != ownerPhase) continue;
        }
        else
        {
            if (boundary.kind != GeometricBoundary::fixed
                || !((boundary.phase + 1 == excludedPhase && boundary.opposedPhase + 1 == ownerPhase)
                    || (boundary.opposedPhase + 1 == excludedPhase && boundary.phase + 1 == ownerPhase))) continue;
        }
        for (int lineId : boundary.lineIds)
        {
            if (lineId < 0 || lineId >= static_cast<int>(definition.lines.size()))
                throw std::runtime_error("invalid physical boundary line ID");
            if (!coversFace(definition.lines[lineId], segment, tolerance)) continue;
            if (found) throw std::runtime_error("ambiguous fitted electric boundary coverage");
            found = &boundary;
        }
    }
    if (!found) throw std::runtime_error("electric boundary is missing, cuts a cell, or is not face fitted: cell="
        + std::to_string(segment.K) + " face=(" + std::to_string(segment.faceCenter.x)
        + "," + std::to_string(segment.faceCenter.y) + ")");
    return *found;
}

}
