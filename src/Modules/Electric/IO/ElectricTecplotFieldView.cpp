#include "Modules/Electric/IO/ElectricTecplotFieldView.h"

#include "Modules/Electric/ElectricPotential.h"
#include "Modules/Electric/ElectricPhaseGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Electric
{
namespace
{
using Layer=Interface::HermiteGalsTimeLayer;
using Piece=std::vector<Vector2D>;
using Pieces=std::vector<Piece>;
using Triangle=std::array<std::size_t,3>;

void require(bool value,const std::string& reason)
{
    if(!value)throw std::runtime_error("[Electric/TecplotFieldView] "+reason);
}
std::string number(double value)
{
    std::ostringstream out;out.imbue(std::locale::classic());
    out<<std::setprecision(17)<<value;return out.str();
}
const char* position(Layer::Position value)
{
    switch(value){
    case Layer::Position::Initial:return "Initial";
    case Layer::Position::Midpoint:return "Midpoint";
    case Layer::Position::End:return "End";
    }
    throw std::runtime_error("[Electric/TecplotFieldView] unavailable clock position");
}
long double cross(Vector2D a,Vector2D b,Vector2D c)
{
    return (static_cast<long double>(b.x)-a.x)*(static_cast<long double>(c.y)-a.y)
        -(static_cast<long double>(b.y)-a.y)*(static_cast<long double>(c.x)-a.x);
}
bool identical(Vector2D a,Vector2D b){return a.x==b.x && a.y==b.y;}
double distance(Vector2D a,Vector2D b){return std::hypot(a.x-b.x,a.y-b.y);}
Vector2D midpoint(Vector2D a,Vector2D b)
{
    return {a.x+(b.x-a.x)/2,a.y+(b.y-a.y)/2};
}
void orient(Piece& points,Vector2D outward)
{
    require(points.size()>=2,"boundary display piece has fewer than two points");
    // Geometry normals retain the actual phase-element orientation.
    const auto direction=points.back()-points.front();
    double projection=-direction.x*outward.y+direction.y*outward.x;
    if(projection==0){
        for(std::size_t i=1;i<points.size();++i){
            const auto segment=points[i]-points[i-1];
            projection=-segment.x*outward.y+segment.y*outward.x;
            if(projection!=0)break;
        }
    }
    require(std::isfinite(projection) && projection!=0,"unresolved display boundary orientation");
    if(projection<0)std::reverse(points.begin(),points.end());
}
Piece cartesian(const LevelSetMappedSegment& cut,const LevelSetMappedSegment::Range& range,
    Vector2D outward)
{
    const auto a=cut.begin(),b=cut.end();
    const auto point=[&](long double fraction){
        const long double offset=range.offset+range.length*fraction;
        if(range.fromEnd)
            return Vector2D{static_cast<double>(static_cast<long double>(b.x)
                -(static_cast<long double>(b.x)-a.x)*offset),
                static_cast<double>(static_cast<long double>(b.y)
                -(static_cast<long double>(b.y)-a.y)*offset)};
        return Vector2D{static_cast<double>(static_cast<long double>(a.x)
            +(static_cast<long double>(b.x)-a.x)*offset),
            static_cast<double>(static_cast<long double>(a.y)
            +(static_cast<long double>(b.y)-a.y)*offset)};
    };
    Piece result{point(0),point(1)};orient(result,outward);return result;
}
Piece interfacePiece(const Interface::HermiteGalsPublication& publication,
    const ElectricPhaseSkeletonFace& face,Vector2D outward)
{
    require(face.component>=0 && face.graph>=0 && !face.nodes.empty(),
        "interface display face has no retained graph or nodes");
    const auto& chain=publication.topology().curves.at(static_cast<std::size_t>(face.component));
    const auto first=chain.frameAt(face.graph,0),last=chain.frameAt(face.graph,1);
    Piece result;result.reserve(face.nodes.size()+2);result.push_back(first.point);
    // The electric geometry retains nodes in the original graph-fraction order.
    // They are already certified same-source samples, not a new root search.
    for(const auto& node:face.nodes)
        if(!identical(result.back(),node.point))result.push_back(node.point);
    if(!identical(result.back(),last.point))result.push_back(last.point);
    orient(result,outward);return result;
}
Piece closeLoop(Pieces pieces,double joinBudget,const std::string& context)
{
    require(!pieces.empty(),context+": missing display boundary");
    Piece polygon=std::move(pieces.back());pieces.pop_back();
    while(!pieces.empty()){
        std::size_t found=pieces.size();
        for(std::size_t p=0;p<pieces.size();++p){
            if(distance(polygon.back(),pieces[p].front())>joinBudget)continue;
            require(found==pieces.size(),context+": ambiguous display endpoint connection");
            found=p;
        }
        if(found==pieces.size()){
            if(distance(polygon.back(),polygon.front())<=joinBudget)
                require(false,context+": a phase element has several boundary loops; hole display is unsupported");
            require(false,context+": original boundary cannot close within floating display precision");
        }
        auto next=std::move(pieces[found]);pieces.erase(pieces.begin()+static_cast<std::ptrdiff_t>(found));
        // Adjacent same-source pieces can embed the same endpoint one ULP
        // apart. Only display incidence is shared; no retained root, geometry,
        // integral, electric DOF or numerical tolerance is changed.
        polygon.insert(polygon.end(),next.begin()+1,next.end());
    }
    require(distance(polygon.back(),polygon.front())<=joinBudget,context+": display boundary is open");
    polygon.pop_back();
    bool changed=true;
    while(changed && polygon.size()>3){
        changed=false;
        for(std::size_t i=0;i<polygon.size();++i){
            const auto before=polygon[(i+polygon.size()-1)%polygon.size()];
            const auto here=polygon[i],after=polygon[(i+1)%polygon.size()];
            const long double dot=(static_cast<long double>(here.x)-before.x)*(after.x-here.x)
                +(static_cast<long double>(here.y)-before.y)*(after.y-here.y);
            if(identical(here,before) || (cross(before,here,after)==0 && dot>=0)){
                polygon.erase(polygon.begin()+static_cast<std::ptrdiff_t>(i));changed=true;break;
            }
        }
    }
    require(polygon.size()>=3,context+": degenerate display polygon");
    long double area=0;
    for(std::size_t i=1;i+1<polygon.size();++i)area+=cross(polygon.front(),polygon[i],polygon[i+1]);
    require(std::isfinite(area) && area!=0,context+": zero/nonfinite display polygon area");
    if(area<0)std::reverse(polygon.begin(),polygon.end());
    return polygon;
}
bool inside(Vector2D point,Vector2D a,Vector2D b,Vector2D c)
{
    return cross(a,b,point)>=0 && cross(b,c,point)>=0 && cross(c,a,point)>=0;
}
std::vector<Triangle> triangulate(const Piece& polygon,const std::string& context)
{
    std::vector<std::size_t> remaining(polygon.size());std::iota(remaining.begin(),remaining.end(),0);
    std::vector<Triangle> result;result.reserve(polygon.size()-2);
    while(remaining.size()>3){
        bool clipped=false;
        for(std::size_t i=0;i<remaining.size();++i){
            const auto a=remaining[(i+remaining.size()-1)%remaining.size()];
            const auto b=remaining[i],c=remaining[(i+1)%remaining.size()];
            if(cross(polygon[a],polygon[b],polygon[c])<=0)continue;
            bool occupied=false;
            for(const auto node:remaining){
                if(node==a || node==b || node==c)continue;
                if(inside(polygon[node],polygon[a],polygon[b],polygon[c])){occupied=true;break;}
            }
            if(occupied)continue;
            result.push_back({a,b,c});remaining.erase(remaining.begin()+static_cast<std::ptrdiff_t>(i));
            clipped=true;break;
        }
        require(clipped,context+": display polygon is not a resolved simple boundary");
    }
    require(cross(polygon[remaining[0]],polygon[remaining[1]],polygon[remaining[2]])>0,
        context+": last display triangle is degenerate");
    result.push_back({remaining[0],remaining[1],remaining[2]});return result;
}
TecplotWriter::Zone zone(const ElectricPhaseResult& field,bool liquid,int phaseId)
{
    TecplotWriter::Zone result;
    result.title=liquid?"liquid electric field":"gas electric field";
    result.element=TecplotWriter::Element::Triangle;
    const bool radial=field.axisymmetric();
    result.coordinateNames=radial?std::array<std::string,2>{{"r [code]","z [code]"}}
        :std::array<std::string,2>{{"x [code]","y [code]"}};
    result.fields={{"V [code]",TecplotWriter::Location::Nodal,{}},
        {radial?"Er [code]":"Ex [code]",TecplotWriter::Location::Nodal,{}},
        {radial?"Ez [code]":"Ey [code]",TecplotWriter::Location::Nodal,{}},
        {"E_mag [code]",TecplotWriter::Location::Nodal,{}},
        {"phase_id",TecplotWriter::Location::Nodal,{}}};
    const auto clock=field.timeLayer(),source=field.geometry().timeLayer();
    result.solutionTime=clock.time();result.strandId=liquid?1:2;
    result.auxiliary={{"Units","code"},{"Method","HHO"},
        {"CoordinateSystem",radial?"axisymmetricRZ":"cartesianXY"},
        {"AcceptedStep",std::to_string(clock.step)},
        {"ElectricClockPosition",position(clock.position)},
        {"ElectricClockBegin",number(clock.begin)},{"ElectricClockEnd",number(clock.end)},
        {"GeometryIdentity",std::to_string(field.geometry().geometryIdentity().value())},
        {"GeometrySourceStep",std::to_string(source.step)},
        {"GeometrySourceClockPosition",position(source.position)},
        {"GeometrySourceClockBegin",number(source.begin)},{"GeometrySourceClockEnd",number(source.end)},
        {"PhaseId",std::to_string(phaseId)},{"PhaseSign",liquid?"negative phi":"positive phi"},
        {"Sampling","element-local one-sided reconstructed P3; E=-grad(V); no cell or phase averaging"},
        {"DisplayGeometry","piecewise-linear chords through original electric graph quadrature nodes and endpoints; Cartesian ranges unchanged"},
        {"DisplayRefinement","each phase polygon triangulated then each triangle split once at edge midpoints"},
        {"DisplayMeaning","visualization only; not phase quadrature or a conservation measure"}};
    return result;
}
void append(TecplotWriter::Zone& target,const ElectricPhaseResult& field,std::size_t element,
    const Piece& polygon,const std::vector<Triangle>& triangles,int phaseId)
{
    for(const auto& triangle:triangles){
        const auto a=polygon[triangle[0]],b=polygon[triangle[1]],c=polygon[triangle[2]];
        const std::array<Vector2D,6> points{{a,b,c,midpoint(a,b),midpoint(b,c),midpoint(c,a)}};
        require(target.points.size()<=static_cast<std::size_t>(std::numeric_limits<int>::max())-points.size(),
            "display node count exceeds Tecplot FE index range");
        const int first=static_cast<int>(target.points.size());
        for(const auto point:points){
            const auto gradient=field.gradientInElement(element,point);
            const double potential=field.potentialInElement(element,point);
            target.points.push_back(point);
            target.fields[0].values.push_back(potential);
            target.fields[1].values.push_back(-gradient.x);
            target.fields[2].values.push_back(-gradient.y);
            target.fields[3].values.push_back(std::hypot(gradient.x,gradient.y));
            target.fields[4].values.push_back(phaseId);
        }
        target.connectivity.push_back({first,first+3,first+5});
        target.connectivity.push_back({first+3,first+1,first+4});
        target.connectivity.push_back({first+5,first+4,first+2});
        target.connectivity.push_back({first+3,first+4,first+5});
    }
}
}
std::vector<TecplotWriter::Zone> buildElectricTecplotFieldView(const QuadtreeMesh& mesh,
    const ElectricPhaseResult& field,const std::array<int,2>& phaseIds)
{
    require(phaseIds[0]>0 && phaseIds[1]>0 && phaseIds[0]!=phaseIds[1],"phase IDs must be distinct positive integers");
    const auto& geometry=field.geometry();const auto& publication=geometry.publication();
    field.validateSource(mesh,publication);
    const auto& snapshot=publication.state().geometry().snapshot();
    std::array<TecplotWriter::Zone,2> zones{{zone(field,true,phaseIds[0]),zone(field,false,phaseIds[1])}};
    std::vector<std::array<Pieces,2>> axis(static_cast<std::size_t>(mesh.cellIndexNum()));
    if(field.axisymmetric()){
        // r=0 is real meridional display geometry even though its zero physical
        // measure correctly removes it from the electric skeleton.
        for(const auto face:publication.domain().faces(mesh)){
            if(face.horizontal || face.center.x!=0)continue;
            const int cell=face.positiveCell>=0?face.positiveCell:face.negativeCell;
            require(cell>=0 && cell<mesh.cellIndexNum(),"axis display boundary has no native cell");
            const Vector2D outward{face.positiveCell>=0?-1.:1.,0};
            for(int half=0;half<2;++half){
                const auto& cut=publication.facePartitions().halfFace(face.nativeRow,half);
                for(const auto& range:cut.ranges())
                    axis[static_cast<std::size_t>(cell)][range.negative?0:1].push_back(cartesian(cut,range,outward));
            }
        }
    }
    const auto& elements=geometry.elements();const auto& faces=geometry.faces();
    for(std::size_t e=0;e<elements.size();++e){
        const auto& element=elements[e];Pieces boundary=axis[static_cast<std::size_t>(element.nativeCell)][element.liquid?0:1];
        require(element.faces.size()==element.geometry.faces.size(),"phase display incidence lost local face ordering");
        for(std::size_t local=0;local<element.faces.size();++local){
            const auto& face=faces.at(static_cast<std::size_t>(element.faces[local]));
            const auto& nodes=element.geometry.faces[local].nodes;
            require(!nodes.empty(),"phase display local face has no oriented node");
            const auto outward=nodes.front().outwardNormal;
            if(face.kind==ElectricPhaseSkeletonFace::Kind::Interface)
                boundary.push_back(interfacePiece(publication,face,outward));
            else{
                const auto& cut=publication.facePartitions().halfFace(face.nativeRow,face.half);
                const auto& range=cut.ranges().at(static_cast<std::size_t>(face.range));
                require(range.negative==element.liquid,"display range belongs to the other phase");
                boundary.push_back(cartesian(cut,range,outward));
            }
        }
        const auto bounds=snapshot.nativeCellBounds(element.nativeCell);
        const double scale=std::max({1.,std::abs(bounds[0].x),std::abs(bounds[0].y),
            std::abs(bounds[1].x),std::abs(bounds[1].y)});
        const double joinBudget=64*std::numeric_limits<double>::epsilon()*scale;
        const std::string context="geometry="+std::to_string(geometry.geometryIdentity().value())
            +" cell="+std::to_string(element.nativeCell)+" phase="+std::to_string(phaseIds[element.liquid?0:1]);
        const auto polygon=closeLoop(std::move(boundary),joinBudget,context);
        for(const auto point:polygon)
            require(point.x>=bounds[0].x && point.x<=bounds[1].x
                && point.y>=bounds[0].y && point.y<=bounds[1].y,context+": display point leaves original native cell");
        append(zones[element.liquid?0:1],field,e,polygon,triangulate(polygon,context),phaseIds[element.liquid?0:1]);
    }
    std::vector<TecplotWriter::Zone> result;
    for(auto& value:zones)if(!value.points.empty())result.push_back(std::move(value));
    require(!result.empty(),"no actual phase elements to display");return result;
}
}
