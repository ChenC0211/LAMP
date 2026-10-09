#include "Modules/Interface/HermiteGalsMixedClosure.h"

#include "Core/Geometry/CompactQuarticDerivative.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Interface
{
namespace
{
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
}
HermiteGalsMixedClosureResult::HermiteGalsMixedClosureResult(HermiteGalsIdentity primary,HermiteGalsTimeLayer time,
    std::vector<double> mixed,HermiteGalsMixedClosureReport report)
    : primary_(primary),time_(time),mixed_(std::move(mixed)),report_(report) {}

HermiteGalsMixedClosureResult closeHermiteGalsMixedUniform(const QuadtreeMesh& mesh,const HermiteGalsCandidate& candidate)
{
    HermiteGalsIdentity identity;int vertex=-1;const char* stage="input";
    try{
        identity=candidate.identity();candidate.validateMesh(mesh);
        require(!candidate.hasGeometry(),"mixed closure requires an unpublished primary revision");
        const auto& primary=candidate.primary();
        const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
        const int nx=grid.nx(),ny=grid.ny();
        require(nx>=2 && ny>=2,"compact Q4 closure needs at least two cells per direction");
        require(primary.vertices.size()==static_cast<std::size_t>(mesh.vertIndexNum()),"complete uniform vertex skeleton required");
        require(primary.edges.size()==static_cast<std::size_t>(mesh.hEdgeIndexNum()+mesh.vEdgeIndexNum()),
            "complete uniform edge skeleton required");
        const bool radial=primary.measure==HermiteGalsMeasure::AxisymmetricVolume;
        require(!radial || grid.minCoord().x>=0,"RZ closure requires a nonnegative radial domain");
        const auto idAt=[&](int i,int j){return mesh.verts(grid.vertices()[j*(nx+1)+i]).index();};
        const auto coordinate=[&](int id){return mesh.verts(mesh.vertIndexToLabel[id]).coord();};
        std::map<std::array<int,2>,double> normals;
        for(const auto& edge:primary.edges)normals.emplace(edge.vertices,edge.dn);
        const auto normal=[&](int a,int b){const auto found=normals.find({{a,b}});
            require(found!=normals.end(),"missing native uniform master edge");return found->second;};
        // Validate every expected side, not just count/orientation. No query
        // silently substitutes a primary edge spanning multiple leaf cells.
        for(int j=0;j<=ny;++j)for(int i=0;i<nx;++i)normal(idAt(i,j),idAt(i+1,j));
        for(int j=0;j<ny;++j)for(int i=0;i<=nx;++i)normal(idAt(i,j),idAt(i,j+1));
        const bool axis=radial && grid.minCoord().x==0;
        if(axis){
            for(int j=0;j<=ny;++j){vertex=idAt(0,j);require(primary.vertices[vertex].dx==0,"nonzero radial first jet on axis");}
            for(int j=0;j<ny;++j)require(normal(idAt(0,j),idAt(0,j+1))==0,"nonzero radial edge derivative on axis");
        }
        std::vector<double> mixed(primary.vertices.size());HermiteGalsMixedClosureReport report;
        report.vertices=mixed.size();
        for(int j=0;j<=ny;++j)for(int i=0;i<=nx;++i){
            vertex=idAt(i,j);stage="directional closure";
            const auto rowCoord=[&](int s){const double a=coordinate(idAt(s/2,j)).x;
                return s%2?a+(coordinate(idAt(s/2+1,j)).x-a)*.5:a;};
            const auto columnCoord=[&](int s){const double a=coordinate(idAt(i,s/2)).y;
                return s%2?a+(coordinate(idAt(i,s/2+1)).y-a)*.5:a;};
            const auto rowValue=[&](int s){return s%2?normal(idAt(s/2,j),idAt(s/2+1,j)):primary.vertices[idAt(s/2,j)].dy;};
            const auto columnValue=[&](int s){return s%2?normal(idAt(i,s/2),idAt(i,s/2+1)):primary.vertices[idAt(i,s/2)].dx;};
            // Even scalar -> even phi_y across r=0, so its centered odd
            // derivative is exactly zero. Axis phi_x data were checked above.
            const double dx=axis && i==0?0:Quantum::Geometry::compactQuarticDerivativeAtNode(2*i,2*nx+1,rowCoord,rowValue,require);
            const double dy=Quantum::Geometry::compactQuarticDerivativeAtNode(2*j,2*ny+1,columnCoord,columnValue,require);
            const double value=.5*dx+.5*dy;
            require(std::isfinite(value),"nonfinite mixed closure result");
            mixed[vertex]=value;
            const double scale=std::max({1.0,std::abs(dx),std::abs(dy)});
            const double disagreement=2*std::abs(.5*dx/scale-.5*dy/scale);
            report.maxRelativeDirectionalDisagreement=std::max(report.maxRelativeDirectionalDisagreement,disagreement);
        }
        stage="final mesh binding";vertex=-1;candidate.validateMesh(mesh);
        return HermiteGalsMixedClosureResult(identity,candidate.timeLayer(),std::move(mixed),report);
    }catch(const std::exception& error){
        throw std::runtime_error("[Interface/HermiteGalsMixedClosure primary="+std::to_string(identity.value())
            +" vertex="+std::to_string(vertex)+"] "+stage+": "+error.what());
    }
}
} // namespace Quantum::Interface
