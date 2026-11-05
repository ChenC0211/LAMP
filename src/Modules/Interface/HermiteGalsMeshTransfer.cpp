#include "Modules/Interface/HermiteGalsMeshTransfer.h"

#include "Core/Geometry/ActiveQuarticHermiteSkeleton.h"
#include "Core/Geometry/CellChartPoint.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace LAMP::Interface
{
namespace
{
using Snapshot = LevelSetHermiteGeometry::Snapshot;
using Key = std::array<int,2>;
using Chart = LAMP::Geometry::CellChartPoint;
void require(bool ok,const char* why) { if(!ok)throw std::invalid_argument(why); }
bool same(Vector2D a,Vector2D b) { return a.x==b.x && a.y==b.y; }
bool contains(const std::array<Vector2D,2>& a,const std::array<Vector2D,2>& b)
{return a[0].x<=b[0].x && a[0].y<=b[0].y && a[1].x>=b[1].x && a[1].y>=b[1].y;}
std::vector<std::array<int,3>> topologySignature(const HermiteGalsPublication& source)
{
    std::vector<std::array<int,3>> result;
    for(const auto& c:source.topology().curves){const auto a=static_cast<int>(c.beginEndpointKind()),b=static_cast<int>(c.endEndpointKind());
        result.push_back({c.closed()?1:0,std::min(a,b),std::max(a,b)});}
    std::sort(result.begin(),result.end());return result;
}
}
struct HermiteGalsMeshTransfer::Data
{
    Data(HermiteGalsIdentity old,HermiteGalsPublication next,HermiteGalsMeshTransferReport diagnostic)
        :source(old),publication(std::move(next)),report(std::move(diagnostic)){}
    HermiteGalsIdentity source;
    HermiteGalsPublication publication;
    HermiteGalsMeshTransferReport report;
};
HermiteGalsMeshTransfer::HermiteGalsMeshTransfer(std::shared_ptr<const Data> data):data_(std::move(data)){}
const HermiteGalsMeshTransfer::Data& HermiteGalsMeshTransfer::checked() const
{require(bool(data_),"[Interface/HermiteGalsMeshTransfer] moved-from result");return *data_;}
HermiteGalsIdentity HermiteGalsMeshTransfer::sourceIdentity() const {return checked().source;}
const HermiteGalsPublication& HermiteGalsMeshTransfer::publication() const {return checked().publication;}
const HermiteGalsMeshTransferReport& HermiteGalsMeshTransfer::report() const {return checked().report;}

HermiteGalsMeshTransfer transferHermiteGalsFirstJetMesh(const QuadtreeMesh& sourceMesh,
    const HermiteGalsPublication& source,const QuadtreeMesh& targetMesh,HermiteGalsMeshTransferKind kind,
    const HermiteGalsPublicationConfig& config)
{return HermiteGalsMeshTransfer::prepare(sourceMesh,source,targetMesh,kind,config);}

HermiteGalsMeshTransfer HermiteGalsMeshTransfer::prepare(const QuadtreeMesh& sourceMesh,
    const HermiteGalsPublication& source,const QuadtreeMesh& targetMesh,HermiteGalsMeshTransferKind kind,
    const HermiteGalsPublicationConfig& config)
{
    HermiteGalsIdentity identity;const char* stage="input";
    try{
        identity=source.state().identity();source.validateMesh(sourceMesh);
        require(kind==HermiteGalsMeshTransferKind::ExactRefinement || kind==HermiteGalsMeshTransferKind::FirstJetProjection,
            "unknown transfer kind");
        const bool exact=kind==HermiteGalsMeshTransferKind::ExactRefinement;
        const auto& old=source.state();const auto& input=old.primary();const auto& scalar=old.geometry().snapshot();
        require(input.transport==HermiteGalsTransport::FirstJet && input.phaseVolumes.empty(),"requires FirstJet source without inventory");
        require(input.geometryCells.size()==static_cast<std::size_t>(sourceMesh.cellIndexNum()),"requires full scalar support");
        require(scalar.usesNativeLeaves() && scalar.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,
            "requires native Q4 source");
        require(same(sourceMesh.minCoord,targetMesh.minCoord) && same(sourceMesh.maxCoord,targetMesh.maxCoord),"mesh transfer changes domain");
        const bool radial=input.measure==HermiteGalsMeasure::AxisymmetricVolume;
        HermiteGalsMeshTransferReport report;report.kind=kind;
        const auto point=[&](int v){return targetMesh.verts(targetMesh.vertIndexToLabel[v]).coord();};
        const auto box=[&](int c){const auto& cell=targetMesh.cells(targetMesh.cellIndexToLabel[c]);
            return std::array<Vector2D,2>{targetMesh.verts(cell.cellToVert[0]).coord(),targetMesh.verts(cell.cellToVert[2]).coord()};};
        stage="nested target coverage";
        for(int c=0;c<targetMesh.cellIndexNum();++c){const auto b=box(c);int overlaps=0,parent=-1;
            for(int p:scalar.nativeCellCandidates(b[0],b[1])){const auto a=scalar.nativeCellBounds(p);
                if(std::max(a[0].x,b[0].x)>=std::min(a[1].x,b[1].x) || std::max(a[0].y,b[0].y)>=std::min(a[1].y,b[1].y))continue;
                require(contains(a,b) || contains(b,a),"requires nested source/target leaves");++overlaps;
                if(contains(a,b))parent=p;
            }
            require(overlaps>0,"target leaf has no source coverage");
            if(parent>=0){require(overlaps==1,"target leaf has ambiguous source parent");const auto a=scalar.nativeCellBounds(parent);
                if(same(a[0],b[0]) && same(a[1],b[1]))++report.unchangedCells;else ++report.refinedCells;
            }else{require(!exact,"ExactRefinement cannot coarsen or span source leaves");++report.coarsenedCells;}
        }
        stage="target master topology";
        LAMP::Geometry::ActiveQuarticHermiteSkeletonInput topology;
        topology.activeCells.assign(targetMesh.cellIndexNum(),1);topology.nodes.resize(targetMesh.vertIndexNum());
        std::map<Key,int> edgeOwner;std::vector<int> vertexOwner(targetMesh.vertIndexNum(),-1);
        for(int c=0;c<targetMesh.cellIndexNum();++c){const auto& cell=targetMesh.cells(targetMesh.cellIndexToLabel[c]);
            for(int k=0;k<4;++k){const int v=targetMesh.verts(cell.cellToVert[k]).index();if(vertexOwner[v]<0)vertexOwner[v]=c;
                Key key{v,targetMesh.verts(cell.cellToVert[(k+1)%4]).index()};const auto a=point(key[0]),b=point(key[1]);
                if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);edgeOwner.emplace(key,c);}
        }
        for(const auto& e:edgeOwner)topology.edgeMidpoints.push_back({e.first,0,0});
        // Zero data serve only to ask the existing constraint engine which
        // target entries are independent. This scalar is never published.
        const LAMP::Geometry::ActiveQuarticHermiteSkeleton owners(targetMesh,topology,radial);
        std::optional<Snapshot> restricted;
        if(exact){stage="exact scalar restriction";restricted.emplace(Snapshot::refineNativeQuartic(targetMesh,scalar));}
        const auto sample=[&](int cell,Chart chart,Vector2D display){
            return exact?restricted->scalarJetAtNativeCellChart(cell,chart):scalar.scalarJetNative(display);};
        stage="first-jet transfer";
        HermiteGalsPrimaryState primary;primary.measure=input.measure;primary.transport=HermiteGalsTransport::FirstJet;
        std::vector<double> mixedByVertex(targetMesh.vertIndexNum(),0),centers;
        const Key none{-1,-1};
        for(const auto& v:owners.vertices())if(v.masterSide==none){const int c=vertexOwner.at(v.index);require(c>=0,"master vertex lacks owning leaf");
            const auto p=point(v.index);const auto b=box(c);const Chart chart{{0,p.x==b[1].x},{0,p.y==b[1].y}};
            const auto j=sample(c,chart,p);primary.vertices.push_back({v.index,j.value,j.dx,j.dy});mixedByVertex[v.index]=j.dxy;
        }
        std::sort(primary.vertices.begin(),primary.vertices.end(),[](const auto& a,const auto& b){return a.vertex<b.vertex;});
        for(const auto& e:owners.edges())if(e.parentSide==none){const int c=edgeOwner.at(e.vertices);
            const auto a=point(e.vertices[0]),b=point(e.vertices[1]);const auto bounds=box(c);const bool horizontal=a.y==b.y;
            const Chart chart=horizontal?Chart{{.5,false},{0,a.y==bounds[1].y}}:Chart{{0,a.x==bounds[1].x},{.5,false}};
            const auto j=sample(c,chart,a+(b-a)*.5);primary.edges.push_back({e.vertices,j.value,horizontal?j.dy:j.dx});
        }
        std::sort(primary.edges.begin(),primary.edges.end(),[](const auto& a,const auto& b){return a.vertices<b.vertices;});
        for(int c=0;c<targetMesh.cellIndexNum();++c){primary.geometryCells.push_back(c);
            centers.push_back(sample(c,{{.5,false},{.5,false}},targetMesh.cells(targetMesh.cellIndexToLabel[c]).centerCoord()).value);}
        HermiteGalsCandidate candidate(targetMesh,old.timeLayer(),primary);
        stage="unique target publication";
        auto publication=[&] {
            if(exact){HermiteGalsClosure closure;closure.cellBubbleCoefficients=centers;
                for(const auto& v:primary.vertices)closure.vertexMixedDerivatives.push_back(mixedByVertex[v.vertex]);
                return HermiteGalsPublication::publishScalar(targetMesh,std::move(candidate),closure,std::move(*restricted),config);
            }
            const auto mixed=closeHermiteGalsMixedAmr(targetMesh,candidate,config.volume.symmetryPlanes);
            return publishHermiteGalsFirstJetGeometry(targetMesh,std::move(candidate),mixed,{candidate.identity(),centers},config);
        }();
        stage="topology preservation";
        require(topologySignature(source)==topologySignature(publication),"mesh projection changed components or endpoint kinds");
        report.sourceVolume=source.volumes().report().totalRawLiquid;report.targetVolume=publication.volumes().report().totalRawLiquid;
        const auto capacity=source.volumes().report().totalCapacity;
        require(std::isfinite(capacity) && capacity>0,"invalid source domain measure");
        report.relativeDomainVolumeChange=(report.targetVolume-report.sourceVolume)/capacity;
        require(std::isfinite(report.relativeDomainVolumeChange),"nonfinite mesh-transfer volume diagnostic");
        return HermiteGalsMeshTransfer(std::make_shared<const Data>(identity,std::move(publication),report));
    }catch(const std::exception& e){throw std::invalid_argument("[Interface/HermiteGalsMeshTransfer source="+std::to_string(identity.value())
        +" stage="+stage+"] "+e.what());}
}
}
