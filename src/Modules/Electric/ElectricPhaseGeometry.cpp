#include "Modules/Electric/ElectricPhaseGeometry.h"

#include "Modules/Interface/HermiteGalsPhaseQuadrature.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <exception>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <string>
#include <utility>

namespace LAMP::Electric
{
namespace
{
namespace Electric=ElectricDiscretization;
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point start)
{return std::chrono::duration<double>(Clock::now()-start).count();}
constexpr std::array<double,8> gx{{-.96028985649753623168,-.79666647741362673959,
    -.52553240991632898582,-.18343464249564980494,.18343464249564980494,
    .52553240991632898582,.79666647741362673959,.96028985649753623168}};
constexpr std::array<double,8> gw{{.10122853629037625915,.22238103445337447054,
    .31370664587788728734,.36268378337836198297,.36268378337836198297,
    .31370664587788728734,.22238103445337447054,.10122853629037625915}};
void require(bool valid,const char* reason)
{
    if(!valid)throw std::invalid_argument(std::string("[ElectricPhaseGeometry] ")+reason);
}
void require(bool valid,const std::string& reason)
{
    if(!valid)throw std::invalid_argument("[ElectricPhaseGeometry] "+reason);
}
bool sameClock(const Interface::HermiteGalsTimeLayer& a,const Interface::HermiteGalsTimeLayer& b)
{return a.step==b.step && a.begin==b.begin && a.end==b.end && a.position==b.position;}
// Compensated sums remain effective on MSVC, where long double is double.
struct Sum
{
    double total=0,correction=0;
    void add(double value)
    {
        const double next=total+value;
        correction+=std::abs(total)>=std::abs(value)?(total-next)+value:(value-next)+total;
        total=next;
    }
    double value() const{return total+correction;}
};
using Moments=std::array<double,21>; // Ordinary dA monomials through P5.
struct BasisOrigin {Vector2D center;double scale;};
struct Rule
{
    std::vector<ElectricPhaseElement> elements;
    std::vector<ElectricPhaseSkeletonFace> faces;
    std::vector<Moments> moments;
    std::vector<Electric::LocalOperator> operators;
    ElectricPhaseGeometryReport report;
};
int elementIndex(const std::vector<std::array<int,2>>& indices,int cell,bool liquid)
{
    require(cell>=0 && cell<static_cast<int>(indices.size()),"invalid native face/graph incidence");
    const int result=indices[static_cast<std::size_t>(cell)][liquid?0:1];
    require(result>=0,"face points to an absent phase element");
    return result;
}
void appendFace(Rule& rule,int element,int face,double sign)
{
    auto& owner=rule.elements[static_cast<std::size_t>(element)];
    const auto& source=rule.faces[static_cast<std::size_t>(face)];
    Electric::FaceQuadrature q;q.length=source.length;q.nodes=source.nodes;
    for(auto& node:q.nodes) {
        node.outwardNormal.x*=sign;node.outwardNormal.y*=sign;
    }
    owner.faces.push_back(face);owner.geometry.faces.push_back(std::move(q));
}
double difference(const Eigen::MatrixXd& a,const Eigen::MatrixXd& b)
{
    require(a.rows()==b.rows() && a.cols()==b.cols() && a.allFinite() && b.allFinite(),
        "quadrature refinement changed the local matrix layout or produced nonfinite entries");
    const double numerator=(a-b).norm(),scale=a.norm()+b.norm();
    require(std::isfinite(numerator) && std::isfinite(scale),"matrix norm overflow during refinement");
    if(numerator==0)return 0;
    require(scale>0,"nonzero matrix difference has zero scale");return numerator/scale;
}
void validateTopology(const Interface::HermiteGalsPublication& publication,int cellCount)
{
    const auto& volumes=publication.volumes().cells();
    require(volumes.size()==static_cast<std::size_t>(cellCount),"P1 requires the complete active native domain");
    std::vector<std::vector<std::pair<int,int>>> incidences(static_cast<std::size_t>(cellCount));
    for(std::size_t c=0;c<publication.topology().curves.size();++c) {
        const auto& graphs=publication.topology().curves[c].graphs();
        for(std::size_t g=0;g<graphs.size();++g) {
            const auto& graph=graphs[g];
            require(!graph.nativeSide && graph.nativeCellIndex>=0 && graph.nativeCellIndex<cellCount,
                "P1 requires ordinary noncoincident native interface graphs");
            incidences[static_cast<std::size_t>(graph.nativeCellIndex)].push_back(
                {static_cast<int>(c),static_cast<int>(g)});
        }
    }
    std::vector<bool> seen(static_cast<std::size_t>(cellCount),false);
    for(const auto& cell:volumes) {
        require(cell.cell>=0 && cell.cell<cellCount && !seen[static_cast<std::size_t>(cell.cell)],
            "duplicate or invalid native volume cell");
        seen[static_cast<std::size_t>(cell.cell)]=true;
        const auto& list=incidences[static_cast<std::size_t>(cell.cell)];
        require(cell.certifiedPhase==-1 || cell.certifiedPhase==0 || cell.certifiedPhase==1,
            "invalid certified phase");
        if(cell.certifiedPhase!=0) {
            require(list.empty(),"pure phase cell unexpectedly owns a free interface graph");
            continue;
        }
        require(!list.empty(),"cut cell has no regular interface graph");
        const int component=list.front().first;
        const auto& chain=publication.topology().curves[static_cast<std::size_t>(component)];
        std::vector<bool> belongs(chain.graphs().size(),false);
        for(const auto& item:list) {
            require(item.first==component,"several interface components occupy one P1 cell");
            belongs[static_cast<std::size_t>(item.second)]=true;
        }
        int runs=0;
        for(std::size_t g=0;g<belongs.size();++g) {
            const bool previous=g>0?belongs[g-1]:chain.closed()?belongs.back():false;
            if(belongs[g] && !previous)++runs;
        }
        if(runs==0 && std::all_of(belongs.begin(),belongs.end(),[](bool x){return x;}))runs=1;
        require(runs==1,"one phase has disconnected interface graph runs in a P1 cell");
    }
    for(bool present:seen)require(present,"missing native volume cell");
}
Rule makeRule(const QuadtreeMesh& mesh,const Interface::HermiteGalsPublication& publication,
    const ElectricPhaseGeometryConfig& config,int refinement,ElectricPhaseGeometryTimings& timings,
    const std::vector<BasisOrigin>* fixedBasis=nullptr)
{
    Rule rule;
    const auto& snapshot=publication.state().geometry().snapshot();
    const bool rz=publication.domain().axisymmetric();
    const int cellCount=mesh.cellIndexNum();
    std::vector<std::array<int,2>> indices(static_cast<std::size_t>(cellCount),{{-1,-1}});
    for(int cell=0;cell<cellCount;++cell) {
        const auto& source=publication.volumes().cells()[static_cast<std::size_t>(cell)];
        // Native publication's volume order is checked independently below;
        // use the source cell identity rather than treating its ordinal as one.
        const auto bounds=snapshot.nativeCellBounds(source.cell);
        const auto& native=mesh.cells(mesh.cellIndexToLabel[source.cell]);
        for(bool liquid:{true,false}) {
            if(source.certifiedPhase!=0 && (source.certifiedPhase<0)!=liquid)continue;
            ElectricPhaseElement element;element.nativeCell=source.cell;element.liquid=liquid;
            element.geometry.center=native.centerCoord();
            element.geometry.scale=std::hypot(bounds[1].x-bounds[0].x,bounds[1].y-bounds[0].y);
            element.geometry.axisymmetric=rz;
            require(std::isfinite(element.geometry.scale) && element.geometry.scale>0,"unresolved native cell scale");
            indices[static_cast<std::size_t>(source.cell)][liquid?0:1]=static_cast<int>(rule.elements.size());
            rule.elements.push_back(std::move(element));rule.moments.push_back(Moments{});
        }
    }
    const auto phaseStarted=Clock::now();
    Interface::HermiteGalsPhaseQuadratureConfig volumeConfig;
    const double tolerance=1e-8*std::pow(.1,refinement);
    volumeConfig.relativeMomentTolerance=tolerance;volumeConfig.relativeVolumeTolerance=tolerance;
    volumeConfig.maximumDepth=config.maximumDepth;
    volumeConfig.maximumSliceQueriesPerCell=config.maximumSliceQueriesPerCell;
    volumeConfig.maximumPoints=config.maximumPoints;
    const auto quadrature=Interface::HermiteGalsPhaseQuadrature::build(mesh,publication,volumeConfig);
    require(quadrature.geometryIdentity()==publication.state().identity()
        && sameClock(quadrature.timeLayer(),publication.state().timeLayer()),"volume cubature source or clock differs");
    for(const auto& point:quadrature.points()) {
        const int element=elementIndex(indices,point.cell,point.liquid);
        const double measure=rz?CartesianMeasure::radialWeight(point.point.x):1.;
        require(std::isfinite(measure) && measure>0 && std::isfinite(point.physicalWeight)
            && point.physicalWeight>0,"positive dA cannot be recovered from a zero/nonfinite physical weight");
        const double area=static_cast<double>(static_cast<long double>(point.physicalWeight)/measure);
        require(std::isfinite(area) && area>0,"unrepresentable positive ordinary phase weight");
        auto& target=rule.elements[static_cast<std::size_t>(element)];
        target.geometry.volume.push_back({point.point,area});
    }
    rule.report.volumePoints=quadrature.points().size();
    require(!fixedBasis || fixedBasis->size()==rule.elements.size(),
        "frozen phase basis origins changed layout");
    for(std::size_t e=0;e<rule.elements.size();++e) {
        auto& geometry=rule.elements[e].geometry;
        if(fixedBasis) {
            geometry.center=(*fixedBasis)[e].center;geometry.scale=(*fixedBasis)[e].scale;
        } else {
            // Ordinary centroid and RMS coordinate scale condition the basis.
            // Freeze both from level zero so refinement compares the same DOFs;
            // no physical point, measure, face length or polynomial space changes.
            Sum area,x,y;
            const Vector2D origin=geometry.center;
            for(const auto& node:geometry.volume) {
                area.add(node.weight);
                x.add(node.weight*(node.point.x-origin.x));
                y.add(node.weight*(node.point.y-origin.y));
            }
            require(std::isfinite(area.value()) && area.value()>0,
                "missing positive phase area for the polynomial basis origin");
            geometry.center={origin.x+x.value()/area.value(),origin.y+y.value()/area.value()};
            Sum secondMoment;
            for(const auto& node:geometry.volume) {
                const double dx=node.point.x-geometry.center.x,dy=node.point.y-geometry.center.y;
                secondMoment.add(node.weight*(dx*dx+dy*dy));
            }
            geometry.scale=std::sqrt(secondMoment.value()/(2*area.value()));
        }
        require(std::isfinite(geometry.scale) && geometry.scale>0,
            "unresolved ordinary phase RMS polynomial scale");
        const auto bounds=snapshot.nativeCellBounds(rule.elements[e].nativeCell);
        require(std::isfinite(geometry.center.x) && std::isfinite(geometry.center.y)
            && geometry.center.x>=bounds[0].x && geometry.center.x<=bounds[1].x
            && geometry.center.y>=bounds[0].y && geometry.center.y<=bounds[1].y,
            "phase polynomial basis origin lies outside its native cell");
        std::array<Sum,21> sums;
        for(const auto& node:geometry.volume) {
            const double x=(node.point.x-geometry.center.x)/geometry.scale;
            const double y=(node.point.y-geometry.center.y)/geometry.scale;
            std::array<double,6> xp{{1,0,0,0,0,0}},yp{{1,0,0,0,0,0}};
            for(int d=1;d<=5;++d){xp[d]=xp[d-1]*x;yp[d]=yp[d-1]*y;}
            int moment=0;
            for(int degree=0;degree<=5;++degree)for(int a=0;a<=degree;++a)
                sums[static_cast<std::size_t>(moment++)].add(node.weight*xp[a]*yp[degree-a]);
        }
        for(std::size_t j=0;j<sums.size();++j)rule.moments[e][j]=sums[j].value();
        const long double capacity=(static_cast<long double>(bounds[1].x)-bounds[0].x)
            *(static_cast<long double>(bounds[1].y)-bounds[0].y);
        const long double area=rule.moments[e][0];
        require(std::isfinite(area) && area>0 && std::isfinite(capacity) && capacity>0,
            "missing or invalid positive phase element area");
        const double relative=static_cast<double>(area/capacity);
        require(std::isfinite(relative) && relative>=config.minimumRelativePhaseArea,
            "ordinary phase area requires small-cut agglomeration before P1");
        rule.report.minimumRelativePhaseArea=std::min(rule.report.minimumRelativePhaseArea,relative);
    }
    timings.phaseQuadrature_s+=elapsed(phaseStarted);
    const auto cartesianStarted=Clock::now();
    const auto& partitions=publication.facePartitions();
    require(&partitions.source()==&snapshot,"retained face partition uses another scalar");
    const auto& domain=publication.domain();
    for(const auto nativeFace:domain.faces(mesh)) {
        // Axis is natural and has no skeleton polynomial; all nonaxis
        // subfaces retain the actual publication's two half-face identities.
        if(rz && !nativeFace.horizontal && nativeFace.center.x==0)continue;
        require(nativeFace.negativeCell>=0 || nativeFace.positiveCell>=0,"inactive native face");
        const int boundary=nativeFace.negativeCell<0?(nativeFace.horizontal?1:0)
            :nativeFace.positiveCell<0?(nativeFace.horizontal?3:2):-1;
        for(int half=0;half<2;++half) {
            const auto& trace=partitions.traceHalfFace(nativeFace.nativeRow,half);
            require(trace.usesSource(snapshot) && trace.coincidentRanges().empty(),
                "positive-length coincident phase face is outside P1");
            const auto& cut=partitions.halfFace(nativeFace.nativeRow,half);
            const auto a=cut.begin(),b=cut.end();
            const long double width=nativeFace.horizontal?static_cast<long double>(b.x)-a.x
                :static_cast<long double>(b.y)-a.y;
            require(width>0 && std::isfinite(width),"unordered retained half-face");
            const auto& ranges=cut.ranges();
            for(std::size_t r=0;r<ranges.size();++r) {
                const auto& range=ranges[r];
                require(std::isfinite(range.offset) && std::isfinite(range.length)
                    && range.offset>=0 && range.length>0,"invalid retained phase range");
                ElectricPhaseSkeletonFace face;face.kind=ElectricPhaseSkeletonFace::Kind::Cartesian;
                face.boundarySide=boundary;face.liquid=range.negative;face.nativeRow=nativeFace.nativeRow;
                face.half=half;face.range=static_cast<int>(r);
                face.length=static_cast<double>(width*range.length);
                require(std::isfinite(face.length) && face.length>0,"unrepresentable retained phase face length");
                const Vector2D normal=nativeFace.horizontal?Vector2D{0,1}:Vector2D{1,0};
                for(std::size_t q=0;q<gx.size();++q) {
                    const long double fraction=(1.L+gx[q])/2;
                    const long double distance=range.offset+range.length*(range.fromEnd?1-fraction:fraction);
                    const long double position=range.fromEnd
                        ?(nativeFace.horizontal?b.x:b.y)-width*distance
                        :(nativeFace.horizontal?a.x:a.y)+width*distance;
                    const double coordinate=static_cast<double>(position);
                    const Vector2D point=nativeFace.horizontal?Vector2D{coordinate,a.y}:Vector2D{a.x,coordinate};
                    require(std::isfinite(coordinate) && coordinate>(nativeFace.horizontal?a.x:a.y)
                        && coordinate<(nativeFace.horizontal?b.x:b.y),
                        "retained phase face has no resolved interior Gauss coordinate");
                    const double weight=static_cast<double>(width*range.length*gw[q]/2);
                    require(std::isfinite(weight) && weight>0,"unrepresentable positive phase face weight");
                    face.nodes.push_back({point,normal,weight,gx[q]});
                }
                const int id=static_cast<int>(rule.faces.size());rule.faces.push_back(std::move(face));
                if(nativeFace.negativeCell>=0)appendFace(rule,elementIndex(indices,nativeFace.negativeCell,range.negative),id,1.);
                if(nativeFace.positiveCell>=0)appendFace(rule,elementIndex(indices,nativeFace.positiveCell,range.negative),id,-1.);
                ++rule.report.cartesianSubfaces;
            }
        }
    }
    timings.cartesianFaces_s+=elapsed(cartesianStarted);
    const auto interfaceStarted=Clock::now();
    const int pieces=1<<refinement;
    std::vector<std::pair<std::size_t,std::size_t>> graphRows;
    for(std::size_t c=0;c<publication.topology().curves.size();++c)
        for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g)
            graphRows.emplace_back(c,g);
    std::vector<ElectricPhaseSkeletonFace> interfaceFaces(graphRows.size());
    std::vector<std::exception_ptr> graphFailures(graphRows.size());
    const auto graphCount=static_cast<std::ptrdiff_t>(graphRows.size());
    // Each retained graph queries immutable source/cuts with local root scratch.
    // Keep every graph's arithmetic serial; only independent graphs run together.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) if(graphCount>=8)
#endif
    for(std::ptrdiff_t item=0;item<graphCount;++item) {
        const auto index=static_cast<std::size_t>(item);
        try {
            const auto c=graphRows[index].first,g=graphRows[index].second;
            const auto& chain=publication.topology().curves[c];
            ElectricPhaseSkeletonFace face;face.kind=ElectricPhaseSkeletonFace::Kind::Interface;
            face.component=static_cast<int>(c);face.graph=static_cast<int>(g);
            struct Node {Interface::HermiteGalsSurfaceQuadrature sample;long double arclength;};
            std::vector<Node> nodes;long double complete=0;
            for(int piece=0;piece<pieces;++piece) {
                const double begin=static_cast<double>(piece)/pieces,end=static_cast<double>(piece+1)/pieces;
                const auto samples=publication.quadrature(c,g,begin,end);
                long double intervalLength=0;
                for(const auto& q:samples) {
                    require(q.sample.geometryIdentity==publication.state().identity()
                        && sameClock(q.sample.timeLayer,publication.state().timeLayer()),
                        "interface quadrature source or clock differs");
                    require(std::isfinite(q.weight) && q.weight>0,"nonpositive same-curve quadrature weight");
                    const double fraction=chain.graphFraction(static_cast<int>(g),q.sample.frame.cellChartPoint);
                    require(std::isfinite(fraction) && fraction>begin && fraction<end,
                        "interface Gauss fraction is unresolved in its retained graph chart");
                    long double partial=0;
                    for(const auto& prefix:publication.quadrature(c,g,begin,fraction)) {
                        require(std::isfinite(prefix.weight) && prefix.weight>0,"invalid same-curve arclength prefix");
                        partial+=prefix.weight;
                    }
                    nodes.push_back({q,complete+partial});intervalLength+=q.weight;
                }
                require(intervalLength>0 && std::isfinite(intervalLength),"unresolved true interface arclength");
                complete+=intervalLength;
            }
            face.length=static_cast<double>(complete);
            require(std::isfinite(face.length) && face.length>0,"unrepresentable shared interface length");
            for(const auto& node:nodes) {
                const double zeta=static_cast<double>(2*node.arclength/complete-1);
                require(std::isfinite(zeta) && zeta>-1 && zeta<1,
                    "normalized interface arclength requires further resolution");
                face.nodes.push_back({node.sample.sample.frame.point,node.sample.sample.frame.normal,
                    node.sample.weight,zeta});
            }
            interfaceFaces[index]=std::move(face);
        } catch(...) {graphFailures[index]=std::current_exception();}
    }
    // Original graph order defines face IDs, accumulated work limits and the
    // first reported failure. No concurrent append or cross-graph reduction.
    for(std::size_t index=0;index<graphRows.size();++index) {
        if(graphFailures[index])std::rethrow_exception(graphFailures[index]);
        auto& face=interfaceFaces[index];
        require(face.nodes.size()<=config.maximumPoints-rule.report.surfacePoints,"interface point work limit");
        rule.report.surfacePoints+=face.nodes.size();
        const auto& row=graphRows[index];
        const int nativeCell=publication.topology().curves[row.first].graphs()[row.second].nativeCellIndex;
        const int id=static_cast<int>(rule.faces.size());rule.faces.push_back(std::move(face));
        appendFace(rule,elementIndex(indices,nativeCell,true),id,1.);
        appendFace(rule,elementIndex(indices,nativeCell,false),id,-1.);
        ++rule.report.interfaceFaces;
    }
    timings.interfaceFaces_s+=elapsed(interfaceStarted);
    rule.report.elements=rule.elements.size();
    rule.report.acceptedRefinement=refinement;
    const auto operatorsStarted=Clock::now();
    rule.operators.reserve(rule.elements.size());
    for(const auto& element:rule.elements) {
        require(!element.geometry.volume.empty() && !element.geometry.faces.empty(),
            "phase element lacks volume or boundary support");
        rule.operators.emplace_back(element.geometry,1.);
        const auto& op=rule.operators.back();
        Eigen::VectorXd constant=Eigen::VectorXd::Zero(op.stiffnessMatrix().cols());
        constant[0]=1;
        for(std::size_t f=0;f<op.faceCount();++f)constant[6+3*static_cast<Eigen::Index>(f)]=1;
        Electric::PotentialVector expected=Electric::PotentialVector::Zero();expected[0]=1;
        const Electric::PotentialVector reconstructed=op.reconstructionMatrix()*constant;
        const double reconstruction=(reconstructed-expected).norm();
        const double scale=op.stiffnessMatrix().norm()*constant.norm();
        const double defect=(op.stiffnessMatrix()*constant).norm();
        require(std::isfinite(reconstruction) && std::isfinite(scale) && scale>0 && std::isfinite(defect),
            "nonfinite constant Green or energy-kernel diagnostic");
        rule.report.maxConstantReconstructionDefect=std::max(rule.report.maxConstantReconstructionDefect,reconstruction);
        rule.report.maxRelativeConstantEnergyDefect=std::max(rule.report.maxRelativeConstantEnergyDefect,defect/scale);
    }
    timings.localOperators_s+=elapsed(operatorsStarted);
    return rule;
}
void compare(Rule& fine,const Rule& coarse)
{
    require(fine.elements.size()==coarse.elements.size() && fine.faces.size()==coarse.faces.size(),
        "quadrature refinement changed phase element or skeleton identities");
    for(std::size_t e=0;e<fine.elements.size();++e) {
        const auto& element=fine.elements[e];const auto& before=coarse.elements[e];
        require(element.nativeCell==before.nativeCell && element.liquid==before.liquid && element.faces==before.faces,
            "quadrature refinement changed native phase incidence");
        const auto& a=fine.operators[e];const auto& b=coarse.operators[e];
        fine.report.maxRelativeReconstructionDifference=std::max(fine.report.maxRelativeReconstructionDifference,
            difference(a.reconstructionMatrix(),b.reconstructionMatrix()));
        fine.report.maxRelativeGradientGramDifference=std::max(fine.report.maxRelativeGradientGramDifference,
            difference(a.gradientGramMatrix(),b.gradientGramMatrix()));
        fine.report.maxRelativeStiffnessDifference=std::max(fine.report.maxRelativeStiffnessDifference,
            difference(a.stiffnessMatrix(),b.stiffnessMatrix()));
        fine.report.maxRelativeCondensedDifference=std::max(fine.report.maxRelativeCondensedDifference,
            difference(a.condensedMatrix(),b.condensedMatrix()));
        const long double scale=fine.moments[e][0]+coarse.moments[e][0];
        require(scale>0 && std::isfinite(scale),"invalid ordinary P5 comparison scale");
        for(std::size_t j=0;j<fine.moments[e].size();++j) {
            const double relative=static_cast<double>(std::abs(fine.moments[e][j]-coarse.moments[e][j])/scale);
            require(std::isfinite(relative),"nonfinite ordinary P5 quadrature difference");
            fine.report.maxRelativeOrdinaryP5MomentDifference=std::max(fine.report.maxRelativeOrdinaryP5MomentDifference,relative);
        }
    }
}
// A final rejected rule reports one element, without per-step diagnostic work.
std::string constantDiagnosis(const Rule& rule)
{
    std::size_t worst=0;double largest=-1;
    for(std::size_t e=0;e<rule.operators.size();++e) {
        const auto& op=rule.operators[e];
        Eigen::VectorXd constant=Eigen::VectorXd::Zero(op.stiffnessMatrix().cols());
        constant[0]=1;
        for(std::size_t f=0;f<op.faceCount();++f)constant[6+3*static_cast<Eigen::Index>(f)]=1;
        auto defect=(op.reconstructionMatrix()*constant).eval();defect[0]-=1;
        if(defect.norm()>largest){largest=defect.norm();worst=e;}
    }
    const auto& element=rule.elements[worst];
    std::ostringstream out;out.precision(17);
    out<<" worst_cell="<<element.nativeCell<<" phase="<<(element.liquid?"liquid":"gas")
        <<" basis_center=("<<element.geometry.center.x<<","<<element.geometry.center.y<<")"
        <<" dA="<<rule.moments[worst][0]<<" basis_scale="<<element.geometry.scale;
    return out.str();
}
bool resolved(const ElectricPhaseGeometryReport& report,double tolerance)
{
    return report.maxRelativeReconstructionDifference<=tolerance
        && report.maxRelativeGradientGramDifference<=tolerance
        && report.maxRelativeStiffnessDifference<=tolerance
        && report.maxRelativeCondensedDifference<=tolerance
        && report.maxRelativeOrdinaryP5MomentDifference<=tolerance
        && report.maxConstantReconstructionDefect<=1e-11
        && report.maxRelativeConstantEnergyDefect<=1e-11;
}
}
struct ElectricPhaseGeometry::Data
{
    explicit Data(Interface::HermiteGalsPublication source):publication(std::move(source)){}
    Interface::HermiteGalsPublication publication;
    std::vector<ElectricPhaseElement> elements;
    std::vector<ElectricPhaseSkeletonFace> faces;
    ElectricPhaseGeometryReport report;
};
ElectricPhaseGeometry::ElectricPhaseGeometry(std::shared_ptr<const Data> data):data_(std::move(data)){}
const ElectricPhaseGeometry::Data& ElectricPhaseGeometry::checked() const
{
    if(!data_)throw std::logic_error("[ElectricPhaseGeometry] moved-from geometry");return *data_;
}
ElectricPhaseGeometry ElectricPhaseGeometry::build(const QuadtreeMesh& mesh,Interface::HermiteGalsPublication publication,
    const ElectricPhaseGeometryConfig& config)
{
    const auto totalStarted=Clock::now();
    ElectricPhaseGeometryTimings timings;
    require(std::isfinite(config.relativeMatrixTolerance) && config.relativeMatrixTolerance>0
        && config.relativeMatrixTolerance<=1e-6 && config.maximumRefinements>=1 && config.maximumRefinements<=8
        && std::isfinite(config.minimumRelativePhaseArea) && config.minimumRelativePhaseArea>0
        && config.minimumRelativePhaseArea<.5 && config.maximumPoints>0 && config.maximumSliceQueriesPerCell>0
        && config.maximumDepth>=0 && config.maximumDepth<=40,"invalid P1 quadrature or small-cut controls");
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    publication.validateMesh(mesh);
    require(publication.state().primary().transport==Interface::HermiteGalsTransport::FirstJet,
        "P1 requires the unique FirstJet publication");
    const auto& snapshot=publication.state().geometry().snapshot();
    require(snapshot.usesNativeLeaves() && !snapshot.usesC2Macro()
        && snapshot.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,
        "P1 requires the current authoritative native Q4 geometry");
    require(publication.domain().cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum())
        && snapshot.nativeCellCount()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "P1 requires complete native leaf support");
    for(int c=0;c<mesh.cellIndexNum();++c)require(snapshot.nativeCellActive(c),"inactive native cell in P1");
    require(!publication.domain().axisymmetric() || (snapshot.radialAxisSymmetry() && mesh.minCoord.x==0),
        "RZ P1 requires the actual even r=0 source");
    validateTopology(publication,mesh.cellIndexNum());
    auto coarse=makeRule(mesh,publication,config,0,timings);
    std::vector<BasisOrigin> fixedBasis;fixedBasis.reserve(coarse.elements.size());
    for(const auto& element:coarse.elements)
        fixedBasis.push_back({element.geometry.center,element.geometry.scale});
    for(int level=1;level<=config.maximumRefinements;++level) {
        auto fine=makeRule(mesh,publication,config,level,timings,&fixedBasis);
        const auto comparisonStarted=Clock::now();compare(fine,coarse);
        timings.matrixComparison_s+=elapsed(comparisonStarted);
        if(resolved(fine.report,config.relativeMatrixTolerance)) {
            publication.validateMesh(mesh);
            auto data=std::make_shared<Data>(std::move(publication));
            data->elements=std::move(fine.elements);data->faces=std::move(fine.faces);data->report=fine.report;
            timings.total_s=elapsed(totalStarted);
            timings.remaining_s=std::max(0.,timings.total_s-timings.phaseQuadrature_s
                -timings.cartesianFaces_s-timings.interfaceFaces_s-timings.localOperators_s
                -timings.matrixComparison_s);
            data->report.timings=timings;
            return ElectricPhaseGeometry(std::move(data));
        }
        coarse=std::move(fine);
    }
    const auto& r=coarse.report;
    std::ostringstream message;message.precision(17);
    message<<"[ElectricPhaseGeometry] actual R/G/K/S, ordinary P5 or constant Green quadrature unresolved:"
        <<" refinement="<<r.acceptedRefinement<<" R="<<r.maxRelativeReconstructionDifference
        <<" G="<<r.maxRelativeGradientGramDifference<<" K="<<r.maxRelativeStiffnessDifference
        <<" S="<<r.maxRelativeCondensedDifference<<" P5="<<r.maxRelativeOrdinaryP5MomentDifference
        <<" constant="<<r.maxConstantReconstructionDefect<<" kernel="<<r.maxRelativeConstantEnergyDefect
        <<constantDiagnosis(coarse);
    throw std::runtime_error(message.str());
}
const Interface::HermiteGalsPublication& ElectricPhaseGeometry::publication() const{return checked().publication;}
Interface::HermiteGalsIdentity ElectricPhaseGeometry::geometryIdentity() const{return publication().state().identity();}
Interface::HermiteGalsTimeLayer ElectricPhaseGeometry::timeLayer() const{return publication().state().timeLayer();}
bool ElectricPhaseGeometry::axisymmetric() const{return publication().domain().axisymmetric();}
const std::vector<ElectricPhaseElement>& ElectricPhaseGeometry::elements() const{return checked().elements;}
const std::vector<ElectricPhaseSkeletonFace>& ElectricPhaseGeometry::faces() const{return checked().faces;}
const ElectricPhaseGeometryReport& ElectricPhaseGeometry::report() const{return checked().report;}
void ElectricPhaseGeometry::validateMesh(const QuadtreeMesh& mesh) const{publication().validateMesh(mesh);}
void ElectricPhaseGeometry::validateSource(const QuadtreeMesh& mesh,const Interface::HermiteGalsPublication& source) const
{
    validateMesh(mesh);source.validateMesh(mesh);
    require(geometryIdentity()==source.state().identity() && sameClock(timeLayer(),source.state().timeLayer())
        && &publication().state().geometry().snapshot()==&source.state().geometry().snapshot(),
        "geometry primary identity, full clock or scalar owner differs");
}
}
