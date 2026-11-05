#include "Core/Geometry/ActiveQuarticHermiteSkeleton.h"

#include "Core/Geometry/Detail/LevelSetHermiteQuarticLeafData.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Geometry
{
namespace
{
using Scalar=Detail::LevelSetHermiteScalar::QuarticLeafData;
ActiveQuarticHermitePatch skeletonPatch(const ActiveQuarticHermiteSkeletonInput& input)
{
    ActiveQuarticHermitePatch patch;
    patch.activeCells=input.activeCells;patch.nodes=input.nodes;patch.edgeMidpoints=input.edgeMidpoints;
    return patch; // no fabricated centers or cell geometry
}
template<class F> auto query(const std::string& stage,F&& action) -> decltype(action())
{
    try {Detail::ExactExpansion::validateEnvironment();return action();}
    catch(const std::exception& error) {throw std::runtime_error("[ActiveQuarticHermiteSkeleton] "+stage+": "+error.what());}
}
std::string edgeName(ActiveQuarticHermiteSkeleton::EdgeKey edge)
{return "edge=("+std::to_string(edge[0])+","+std::to_string(edge[1])+")";}
}
struct ActiveQuarticHermiteSkeleton::Impl
{
    NativeMeshIdentity mesh;
    Scalar scalar;
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    Impl(const QuadtreeMesh& native,const ActiveQuarticHermiteSkeletonInput& input,bool radial)
        : mesh(native),scalar(native,skeletonPatch(input),radial,nullptr,Scalar::Build::Skeleton)
    {
        for(int v:scalar.skeletonVertices())vertices.push_back({v,scalar.vertexMasterSide(v)});
        for(const auto& e:scalar.skeletonEdges())edges.push_back({e,scalar.edgeParentSide(e)});
    }
};
ActiveQuarticHermiteSkeleton::ActiveQuarticHermiteSkeleton(const QuadtreeMesh& mesh,
    const ActiveQuarticHermiteSkeletonInput& input,bool radialAxisSymmetry)
    : impl_(query("construct",[&]{return std::make_shared<const Impl>(mesh,input,radialAxisSymmetry);})) {}
const ActiveQuarticHermiteSkeleton::Impl& ActiveQuarticHermiteSkeleton::implementation() const
{
    if(!impl_)throw std::runtime_error("[ActiveQuarticHermiteSkeleton] unbound moved-from owner");
    return *impl_;
}
const std::vector<ActiveQuarticHermiteSkeleton::Vertex>& ActiveQuarticHermiteSkeleton::vertices() const
{return implementation().vertices;}
const std::vector<ActiveQuarticHermiteSkeleton::Edge>& ActiveQuarticHermiteSkeleton::edges() const
{return implementation().edges;}
std::size_t ActiveQuarticHermiteSkeleton::independentDofCount() const
{return implementation().scalar.independentDofCount();}
void ActiveQuarticHermiteSkeleton::validateMesh(const QuadtreeMesh& mesh) const
{query("mesh binding",[&]{implementation().mesh.validate(mesh);});}
ActiveQuarticHermiteSkeleton::Jet ActiveQuarticHermiteSkeleton::nodeJet(int vertex) const
{return query("vertex="+std::to_string(vertex),[&]{return implementation().scalar.nodeJet(vertex);});}
std::vector<ActiveQuarticHermiteSensitivity> ActiveQuarticHermiteSkeleton::nodeSensitivity(int vertex,int component) const
{return query("vertex="+std::to_string(vertex)+" sensitivity component="+std::to_string(component),
    [&]{return implementation().scalar.pullbackNodeSensitivity(vertex,component);});}
std::array<double,2> ActiveQuarticHermiteSkeleton::edgeMidpoint(EdgeKey edge) const
{return query(edgeName(edge)+" midpoint",[&]{return implementation().scalar.skeletonMidpoint(edge);});}
ActiveQuarticHermiteSkeleton::Jet ActiveQuarticHermiteSkeleton::edgeJet(EdgeKey edge,double coordinate) const
{return query(edgeName(edge)+" trace",[&]{return implementation().scalar.skeletonTrace(edge,coordinate);});}
} // namespace LAMP::Geometry
