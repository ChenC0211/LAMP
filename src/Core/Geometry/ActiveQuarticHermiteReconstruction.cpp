#include "Core/Geometry/ActiveQuarticHermiteReconstruction.h"
#include "Core/Geometry/ExactPointKey.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Geometry/LevelSetQueryReconstruction.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <exception>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace ActiveQuarticHermiteReconstruction
{
struct Prepared::Impl
{
    using Stencil=LevelSetPolynomial::QueryReconstruction::HermiteDataStencil;
    enum class Kind { Node, HorizontalEdge, VerticalEdge, Center };
    struct Datum
    {
        Stencil stencil;
        Kind kind;
        std::size_t target, nodeValue=0;
        std::vector<std::size_t> samples;
    };
    LAMP::Geometry::ActiveQuarticHermitePatch layout;
    std::vector<Vector2D> points;
    std::vector<Datum> data;
    Report report;
    bool radial=false;
};
struct StencilCache::Impl
{
    std::vector<std::uint64_t> identity;
    std::shared_ptr<const Prepared::Impl> prepared;
};
StencilCache::StencilCache() : impl_(std::make_unique<Impl>()) {}
StencilCache::~StencilCache() = default;
Prepared::Prepared(std::shared_ptr<const Impl> impl) : impl_(std::move(impl)) {}
const std::vector<Vector2D>& Prepared::samplePoints() const
{
    if (!impl_) throw std::invalid_argument("[ActiveQuarticHermiteReconstruction] moved-from batch");
    return impl_->points;
}
Prepared::AffectedData Prepared::affectedBy(const std::vector<std::uint8_t>& marked) const
{
    if (!impl_ || marked.size()!=impl_->points.size()
        || std::any_of(marked.begin(),marked.end(),[](auto x){return x>1;}))
        throw std::invalid_argument("[ActiveQuarticHermiteReconstruction] invalid sample influence mask");
    AffectedData result;
    result.nodes.assign(impl_->layout.nodes.size(),0);
    result.edgeMidpoints.assign(impl_->layout.edgeMidpoints.size(),0);
    for (const auto& datum:impl_->data) {
        if (datum.kind==Impl::Kind::Center) continue;
        const bool hit=std::any_of(datum.samples.begin(),datum.samples.end(),[&](auto i){return marked[i]!=0;})
            || (datum.kind==Impl::Kind::Node && marked[datum.nodeValue]);
        if (datum.kind==Impl::Kind::Node) result.nodes[datum.target]=hit;
        else result.edgeMidpoints[datum.target]=hit;
    }
    return result;
}

namespace
{
using Patch=LAMP::Geometry::ActiveQuarticHermitePatch;
using EdgeKey=std::array<int,2>;
void require(bool valid,const char* reason)
{if(!valid)throw std::runtime_error(std::string("[ActiveQuarticHermiteReconstruction] ")+reason);}
struct Side
{
    Vector2D anchor,spacing;
    bool horizontal=false,dependent=false;
};
std::uint64_t bits(double value)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t), "binary64 sample keys required");
    std::uint64_t result;
    std::memcpy(&result, &value, sizeof(value));
    return result;
}
// Compare the exact data consumed by spacing, support admission and native
// dependent-side ownership. Addresses/counts alone are not a mesh identity:
// adapt/cache reload can replace labels, coordinates or connectivity in place.
std::vector<std::uint64_t> stencilIdentity(const QuadtreeMesh& mesh,
    const std::vector<std::uint8_t>& active, const Options& options)
{
    std::vector<std::uint64_t> key;
    key.reserve(8 + 6*mesh.vertIndexNum() + 9*mesh.cellIndexNum()
        + 6*(mesh.hEdgeIndexNum()+mesh.vEdgeIndexNum()));
    const auto integer = [&](int value) { key.push_back(static_cast<std::uint64_t>(value)); };
    const auto point = [&](Vector2D value) { key.push_back(bits(value.x)); key.push_back(bits(value.y)); };
    point(mesh.minCoord); point(mesh.maxCoord);
    integer(mesh.vertIndexNum()); integer(mesh.cellIndexNum());
    integer(mesh.hEdgeIndexNum()); integer(mesh.vEdgeIndexNum());
    for (int k=0; k<mesh.vertIndexNum(); ++k)
    {
        const int label=mesh.vertIndexToLabel[k];
        const auto& vertex=mesh.verts(label);
        integer(label); integer(vertex.index()); integer(vertex.supCell); point(vertex.coord());
    }
    for (int c=0; c<mesh.cellIndexNum(); ++c)
    {
        const int label=mesh.cellIndexToLabel[c];
        const auto& cell=mesh.cells(label);
        integer(label); integer(cell.index()); integer(cell.level());
        for (int vertex:cell.cellToVert) integer(vertex);
    }
    for (bool horizontal:{true,false})
    {
        const int count=horizontal?mesh.hEdgeIndexNum():mesh.vEdgeIndexNum();
        for (int k=0; k<count; ++k)
        {
            const int label=horizontal?mesh.hEdgeIndexToLabel[k]:mesh.vEdgeIndexToLabel[k];
            const auto& edge=horizontal?mesh.hEdges(label):mesh.vEdges(label);
            integer(label); integer(edge.index()); integer(edge.masterCell); integer(edge.slaveCell);
            for (int vertex:edge.edgeToVert) integer(vertex);
        }
    }
    integer(options.radialAxisSymmetry);
    for (const auto* mask:{&active,&options.oneSidedSamples,&options.extensionSideSamples})
    {
        key.push_back(mask->size());
        for (auto value:*mask) key.push_back(value);
    }
    return key;
}
// All mesh/chart borrows end with preparation. The published batch owns only
// geometry, sample indices and output identities, never a sampled field.
class PlanBuilder
{
  public:
    PlanBuilder(const QuadtreeMesh& mesh,const std::vector<std::uint8_t>& active,
        const Options& options, Prepared::Impl& plan)
        : mesh_(mesh),active_(active),options_(options),locator_(mesh),
          reconstruction_(options.radialAxisSymmetry),bounds_(mesh.cellIndexNum()),plan_(plan),report_(plan.report)
    {
        indices_.reserve(static_cast<std::size_t>(mesh.vertIndexNum())*4);
        for(int c=0;c<mesh.cellIndexNum();++c)
        {
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
            bounds_[c]={{mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()}};
        }
    }
    const std::array<Vector2D,2>& bounds(int c) const {return bounds_.at(c);}
    std::size_t pointIndex(Vector2D point)
    {
        ++report_.scalarRequests;
        require(std::isfinite(point.x) && std::isfinite(point.y),"nonfinite sample coordinate");
        const auto key=LAMP::Geometry::exactPointKey(point);
        const auto found=indices_.find(key);
        if(found!=indices_.end())return found->second;
        const auto incident=locator_.incident(point);
        bool active=false;
        for(int k=0;k<incident.count;++k)active=active || active_[incident.indices[k]];
        require(active,"scalar query leaves active native domain");
        ++report_.scalarEvaluations;
        const auto result=plan_.points.size();
        plan_.points.push_back(point);
        indices_.emplace(key,result);
        return result;
    }
    template<std::size_t N> Prepared::Impl::Datum& datum(Vector2D anchor,Vector2D spacing,
        std::array<int,2> half,const std::array<int,N>& vertices, Prepared::Impl::Kind kind, std::size_t target)
    {
        const auto contains=[&](const auto& mask)
        {return !mask.empty() && std::all_of(vertices.begin(),vertices.end(),[&](int k){return mask[k]!=0;});};
        const auto* mask=contains(options_.oneSidedSamples)?&options_.oneSidedSamples:
            contains(options_.extensionSideSamples)?&options_.extensionSideSamples:nullptr;
        try
        {
            ++report_.stencilRequests;
            const std::array<std::uint64_t,7> key{{bits(anchor.x),bits(anchor.y),bits(spacing.x),bits(spacing.y),
                static_cast<std::uint64_t>(half[0]),static_cast<std::uint64_t>(half[1]),
                mask==&options_.oneSidedSamples?1u:mask?2u:0u}};
            auto found=stencils_.find(key);
            if (found==stencils_.end())
            {
                auto stencil=reconstruction_.prepareHermiteDataStencil(anchor,spacing,half,
                    [&](Vector2D lo,Vector2D hi){return admitted(lo,hi,mask);});
                found=stencils_.emplace(key,std::move(stencil)).first;
                ++report_.stencilBuilds;
            }
            plan_.data.push_back({found->second,kind,target});
            auto& data=plan_.data.back();
            reconstruction_.visitHermiteDataSamples(data.stencil,
                [&](Vector2D point){data.samples.push_back(pointIndex(point));});
            return data;
        }
        catch(const std::exception& error)
        {
            std::ostringstream message;
            message<<std::setprecision(17)<<"[ActiveQuarticHermiteReconstruction] datum vertices=";
            for(int vertex:vertices)message<<vertex<<',';
            message<<" anchor=("<<anchor.x<<','<<anchor.y<<") spacing=("<<spacing.x<<','<<spacing.y
                <<") half=("<<half[0]<<','<<half[1]<<") chart="
                <<(mask==&options_.oneSidedSamples?"one-sided":mask?"extension":"active")<<": "<<error.what();
            throw std::runtime_error(message.str());
        }
    }
  private:
    const QuadtreeMesh& mesh_;
    const std::vector<std::uint8_t>& active_;
    const Options& options_;
    CartesianLeafLocator locator_;
    LevelSetPolynomial::QueryReconstruction reconstruction_;
    std::vector<std::array<Vector2D,2>> bounds_;
    // Exact-key lookup only; hash-table iteration never determines ordering.
    std::unordered_map<LAMP::Geometry::ExactPointKey,std::size_t,LAMP::Geometry::ExactPointKeyHash> indices_;
    std::map<std::array<std::uint64_t,7>,Prepared::Impl::Stencil> stencils_;
    Prepared::Impl& plan_;
    Report& report_;
    bool admitted(Vector2D lo,Vector2D hi,const std::vector<std::uint8_t>* mask) const
    {
        if(lo.x<mesh_.minCoord.x || lo.y<mesh_.minCoord.y || hi.x>mesh_.maxCoord.x
            || hi.y>mesh_.maxCoord.y)return false;
        bool covered=false;
        for(int c:locator_.intersectingBox(lo,hi))
        {
            const auto& b=bounds_[c];
            if(std::max(lo.x,b[0].x)>=std::min(hi.x,b[1].x)
                || std::max(lo.y,b[0].y)>=std::min(hi.y,b[1].y))continue;
            if(!active_[c])return false;
            if(mask)for(int label:mesh_.cells(mesh_.cellIndexToLabel[c]).cellToVert)
                if(!(*mask)[mesh_.verts(label).index()])return false;
            covered=true;
        }
        return covered;
    }
};
}

Prepared prepare(const QuadtreeMesh& mesh,const std::vector<std::uint8_t>& activeCells,
    const Options& options, StencilCache* cache, Report* report)
{
    LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
    require(activeCells.size()==static_cast<std::size_t>(mesh.cellIndexNum()),"active mask must use full native cell indices");
    require(!options.radialAxisSymmetry || mesh.minCoord.x==0,"radial extension requires actual r=0 boundary");
    for(auto bit:activeCells)require(bit<=1,"active mask must be binary");
    StencilCache localCache;
    auto& workspace=*(cache?cache:&localCache)->impl_;
    auto identity=stencilIdentity(mesh,activeCells,options);
    if (workspace.prepared && workspace.identity==identity)
    {
        if (report)
        {
            *report=workspace.prepared->report;
            report->reusedTopology=true;
            report->stencilBuilds=0;
        }
        return Prepared(workspace.prepared);
    }
    // Failed preparation must not leave an old identity marked reusable.
    workspace.prepared.reset();
    workspace.identity.clear();
    auto plan=std::make_shared<Prepared::Impl>();
    plan->radial=options.radialAxisSymmetry;
    const double infinity=std::numeric_limits<double>::infinity(),nan=std::numeric_limits<double>::quiet_NaN();
    auto& result=plan->layout;
    result.activeCells=activeCells;
    result.nodes.assign(mesh.vertIndexNum(),{nan,nan,nan,nan});
    result.cellCenters.assign(mesh.cellIndexNum(),nan);
    std::vector<Vector2D> spacing(mesh.vertIndexNum(),{infinity,infinity});
    std::vector<std::uint8_t> used(mesh.vertIndexNum(),0);
    std::map<EdgeKey,Side> sides;
    PlanBuilder sample(mesh,activeCells,options,*plan);
    constexpr int sideCorners[4][2]={{0,1},{3,2},{0,3},{1,2}};
    for(int c=0;c<mesh.cellIndexNum();++c)
    {
        if(!activeCells[c])continue;
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        const auto& bounds=sample.bounds(c);
        const Vector2D h=bounds[1]-bounds[0];
        for(int label:cell.cellToVert)
        {
            const int k=mesh.verts(label).index();used[k]=1;
            spacing[k].x=std::min(spacing[k].x,h.x);spacing[k].y=std::min(spacing[k].y,h.y);
        }
        for(int s=0;s<4;++s)
        {
            const auto& a=mesh.verts(cell.cellToVert[sideCorners[s][0]]);
            const auto& b=mesh.verts(cell.cellToVert[sideCorners[s][1]]);
            auto inserted=sides.emplace(EdgeKey{{a.index(),b.index()}},Side{a.coord(),h,s<2,false});
            auto& side=inserted.first->second;
            side.spacing.x=std::min(side.spacing.x,h.x);side.spacing.y=std::min(side.spacing.y,h.y);
        }
    }
    require(!sides.empty(),"empty active domain");
    for(const auto* mask:{&options.oneSidedSamples,&options.extensionSideSamples})
    {
        require(mask->empty() || mask->size()==used.size(),"sample chart requires full native vertex mask");
        for(std::size_t k=0;k<mask->size();++k)
            require((*mask)[k]<=1 && (!(*mask)[k] || used[k]),"sample chart includes invalid/inactive vertex");
    }
    // These are data-ownership exclusions, not another hanging interpolation.
    // The scalar owner will supply the exact coarse trace and its transpose.
    for(bool horizontal:{true,false})
    {
        const int count=horizontal?mesh.hEdgeIndexNum():mesh.vEdgeIndexNum();
        for(int i=0;i<count;++i)
        {
            const auto& edge=horizontal?mesh.hEdges(mesh.hEdgeIndexToLabel[i]):mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
            if(edge.masterCell<0 || edge.slaveCell<0)continue;
            const auto& a=mesh.cells(edge.masterCell);const auto& b=mesh.cells(edge.slaveCell);
            if(a.level()==b.level() || !activeCells[a.index()] || !activeCells[b.index()])continue;
            const EdgeKey key{{mesh.verts(edge.edgeToVert[0]).index(),mesh.verts(edge.edgeToVert[1]).index()}};
            const auto found=sides.find(key);
            require(found!=sides.end(),"active coarse/fine subedge lacks native full-side identity");
            found->second.dependent=true;
            // The authoritative coarse side also touches these fine cells.
            // Keep its full tangential span, but use the finest actual
            // incident normal spacing for its normal-derivative stencil.
            const auto& coarse=a.level()<b.level()?a:b;
            const auto& fine=a.level()<b.level()?b:a;
            const auto& coarseBounds=sample.bounds(coarse.index());
            const auto point=mesh.verts(edge.edgeToVert[0]).coord();
            const int side=horizontal?(point.y==coarseBounds[0].y?0:1):
                (point.x==coarseBounds[0].x?2:3);
            const EdgeKey coarseKey{{mesh.verts(coarse.cellToVert[sideCorners[side][0]]).index(),
                mesh.verts(coarse.cellToVert[sideCorners[side][1]]).index()}};
            auto& coarseSide=sides.at(coarseKey);
            const auto& fineBounds=sample.bounds(fine.index());
            if(horizontal)coarseSide.spacing.y=std::min(coarseSide.spacing.y,fineBounds[1].y-fineBounds[0].y);
            else coarseSide.spacing.x=std::min(coarseSide.spacing.x,fineBounds[1].x-fineBounds[0].x);
        }
    }
    for(int k=0;k<mesh.vertIndexNum();++k)
    {
        const auto& vertex=mesh.verts(mesh.vertIndexToLabel[k]);
        if(!used[k] || vertex.supCell>=0)continue;
        const auto point=vertex.coord();
        auto& datum=sample.datum(point,spacing[k],{{0,0}},std::array<int,1>{{k}},Prepared::Impl::Kind::Node,k);
        datum.nodeValue=sample.pointIndex(point);
    }
    for(const auto& entry:sides)
    {
        const auto& side=entry.second;
        Patch::EdgeMidpoint data{entry.first,nan,nan};
        if(!side.dependent)
        {
            sample.datum(side.anchor,side.spacing,side.horizontal?std::array<int,2>{{1,0}}:
                std::array<int,2>{{0,1}},entry.first,side.horizontal?Prepared::Impl::Kind::HorizontalEdge:
                Prepared::Impl::Kind::VerticalEdge,result.edgeMidpoints.size());
        }
        result.edgeMidpoints.push_back(data);
    }
    for(int c=0;c<mesh.cellIndexNum();++c)
    {
        if(!activeCells[c])continue;
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        std::array<int,4> corners;
        for(int k=0;k<4;++k)corners[k]=mesh.verts(cell.cellToVert[k]).index();
        const auto& bounds=sample.bounds(c);
        sample.datum(bounds[0],bounds[1]-bounds[0],{{1,1}},corners,Prepared::Impl::Kind::Center,c);
    }
    workspace.identity=std::move(identity);
    workspace.prepared=plan;
    if (report) *report=plan->report;
    return Prepared(std::move(plan));
}

Patch reconstruct(const Prepared& prepared, const std::vector<double>& values)
{
    LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
    require(prepared.impl_ && values.size()==prepared.impl_->points.size(),"prepared scalar array size mismatch");
    for (double value:values) require(std::isfinite(value),"nonfinite supplied scalar");
    const auto& plan=*prepared.impl_;
    const LevelSetPolynomial::QueryReconstruction reconstruction(plan.radial);
    Patch result=plan.layout;
    const auto evaluate=[&](std::size_t index) {
        LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
        const auto& datum=plan.data[index];
        const auto jet=reconstruction.evaluateHermiteDataSamples(datum.stencil,
            [&](std::size_t i){return values[datum.samples[i]];});
        switch (datum.kind)
        {
        case Prepared::Impl::Kind::Node:
            result.nodes[datum.target]={values[datum.nodeValue],jet.derivatives.dx,jet.derivatives.dy,jet.derivatives.dxy};
            break;
        case Prepared::Impl::Kind::Center: result.cellCenters[datum.target]=jet.value; break;
        default:
            auto& edge=result.edgeMidpoints[datum.target];
            edge.value=jet.value;
            edge.normalDerivative=datum.kind==Prepared::Impl::Kind::HorizontalEdge?jet.derivatives.dy:jet.derivatives.dx;
        }
    };
    // Independent output slots, immutable geometry/current samples. No user
    // callback or shared reduction in this region; retain ordered failures.
    if (!plan.data.empty()) evaluate(0);
    std::vector<std::exception_ptr> errors(plan.data.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 32) if(plan.data.size()>64)
#endif
    for (std::ptrdiff_t i=1; i<static_cast<std::ptrdiff_t>(plan.data.size()); ++i)
    {
        try {evaluate(static_cast<std::size_t>(i));}
        catch (...) {errors[i]=std::current_exception();}
    }
    for (std::size_t i=1; i<errors.size(); ++i)
        if (errors[i])
        {
            try {std::rethrow_exception(errors[i]);}
            catch (const std::exception& error) {
                throw std::runtime_error("[ActiveQuarticHermiteReconstruction] prepared datum="
                    +std::to_string(i)+", native target="+std::to_string(plan.data[i].target)+": "+error.what());
            }
        }
    return result;
}

Patch reconstruct(const QuadtreeMesh& mesh,const std::vector<std::uint8_t>& activeCells,
    const std::function<double(Vector2D)>& scalar,const Options& options, StencilCache* cache, Report* report)
{
    require(static_cast<bool>(scalar),"missing scalar sampler");
    Report current;
    const auto prepared=prepare(mesh,activeCells,options,cache,&current);
    std::vector<double> values;
    values.reserve(prepared.samplePoints().size());
    for (Vector2D point:prepared.samplePoints())
    {
        LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
        const double value=scalar(point);
        require(std::isfinite(value),"nonfinite supplied scalar");
        values.push_back(value);
    }
    auto result=reconstruct(prepared,values);
    if (report) *report=current;
    return result;
}
}
