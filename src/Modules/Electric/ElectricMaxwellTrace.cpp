#include "Modules/Electric/ElectricMaxwellTrace.h"

#include "Modules/Electric/ElectricPotential.h"
#include "Modules/Electric/ElectricPhaseGeometry.h"
#include "Solvers/Poisson/ElectricLocalOperator.h"
#include <Eigen/Cholesky>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace LAMP::Electric
{
namespace
{
namespace Electric = ElectricDiscretization;
void require(bool valid,const char* reason,const std::string& context = {})
{
    if(!valid)throw std::invalid_argument("[ElectricMaxwellTrace] "+context+reason);
}
bool sameClock(const Interface::HermiteGalsTimeLayer& a,const Interface::HermiteGalsTimeLayer& b)
{return a.step==b.step && a.begin==b.begin && a.end==b.end && a.position==b.position;}
// Compensation is needed on MSVC too, where long double has double precision.
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
double relative(double residual,double scale,const std::string& context)
{
    require(std::isfinite(residual) && residual>=0 && std::isfinite(scale) && scale>=0,
        "nonfinite arithmetic residual family",context);
    if(scale==0) {
        require(residual==0,"nonzero residual has zero physical family scale",context);
        return 0;
    }
    const double result=residual/scale;
    require(std::isfinite(result),"unrepresentable relative residual",context);
    return result;
}
double cubicBasis(double zeta){return .5*(5*zeta*zeta*zeta-3*zeta);}
double valueP2(const std::array<double,3>& coefficients,const Electric::FacePolynomial& basis)
{
    Sum value;
    for(int j=0;j<3;++j)value.add(coefficients[static_cast<std::size_t>(j)]*basis[j]);
    return value.value();
}
double valueP3(const std::array<double,4>& coefficients,const Electric::FacePolynomial& basis,double p3)
{
    Sum value;
    for(int j=0;j<3;++j)value.add(coefficients[static_cast<std::size_t>(j)]*basis[j]);
    value.add(coefficients[3]*p3);return value.value();
}
template<std::size_t N>
bool finite(const std::array<double,N>& values)
{
    return std::all_of(values.begin(),values.end(),[](double value){return std::isfinite(value);});
}
using ChartKey=std::tuple<double,bool,double,bool>;
ChartKey chartKey(const LAMP::Geometry::CellChartPoint& chart)
{return {chart.x.offset,chart.x.fromUpper,chart.y.offset,chart.y.fromUpper};}
struct CompletionNode
{
    Electric::FacePolynomial basis;
    double p3=0,weight=0,sharedPotential=0,averagePotential=0,psi=0;
};
}
struct ElectricMaxwellTrace::Data
{
    explicit Data(ElectricPhaseResult field):source(std::move(field)){}
    struct Face
    {
        bool available=false;
        std::array<double,4> potential{};
        std::array<double,3> tangentialField{};
        std::array<std::array<double,3>,2> normalDisplacement{};
        std::array<double,3> charge{};
        // Replay of the accepted electric geometry's original piece-wise GL8
        // arclength arithmetic. Keys retain both endpoint-relative chart
        // coordinates; neither rounded points nor rounded fractions identify nodes.
        std::vector<long double> arclengthPrefix;
        std::map<ChartKey,double> coordinates;
        std::vector<ElectricMaxwellSample> samples;
        Vector2D integrated;
    };
    ElectricPhaseResult source;
    ElectricMaxwellReport report;
    std::vector<std::size_t> interface;
    std::vector<Face> faces;
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> graphFaces;
    const Face& face(std::size_t index) const
    {
        require(index<faces.size() && faces[index].available,
            "query requires an original interface face",
            "geometry="+std::to_string(source.geometry().geometryIdentity().value())+
                " face="+std::to_string(index)+" ");
        return faces[index];
    }
};
ElectricMaxwellTrace::ElectricMaxwellTrace(std::shared_ptr<const Data> data):data_(std::move(data)){}
const ElectricMaxwellTrace::Data& ElectricMaxwellTrace::checked() const
{
    if(!data_)throw std::logic_error("[ElectricMaxwellTrace] moved-from trace");
    return *data_;
}
ElectricMaxwellTrace ElectricMaxwellTrace::prepare(const ElectricPhaseResult& source,const ElectricMaxwellControls& controls)
{
    require(std::isfinite(controls.relativeProjectionTolerance)
        && controls.relativeProjectionTolerance>0 && controls.relativeProjectionTolerance<=1e-8,
        "projection tolerance must be in (0,1e-8]");
    require(std::isfinite(controls.relativeChargeJumpTolerance)
        && controls.relativeChargeJumpTolerance>0 && controls.relativeChargeJumpTolerance<=1e-10,
        "physical charge-jump tolerance must be in (0,1e-10]");
    auto data=std::make_shared<Data>(source);
    const auto& geometry=source.geometry();
    const auto& faces=geometry.faces();
    const auto& elements=geometry.elements();
    const auto& epsilon=source.epsilonByPhase();
    data->report.maxRelativeChargeJumpResidual=source.report().maxRelativeInterfaceBalance;
    require(std::isfinite(data->report.maxRelativeChargeJumpResidual)
        && data->report.maxRelativeChargeJumpResidual>=0
        && data->report.maxRelativeChargeJumpResidual<=controls.relativeChargeJumpTolerance,
        "original equilibrated charge-jump receipt exceeds its physical tolerance");
    for(double value:epsilon)
        require(std::isfinite(value) && value>0,"actual solve permittivity must be positive");
    (void)source.timeLayer().time();
    (void)geometry.timeLayer().time();
    data->faces.resize(faces.size());
    const std::size_t missing=std::numeric_limits<std::size_t>::max();
    std::vector<std::array<std::size_t,2>> owners(faces.size(),{{missing,missing}});
    for(std::size_t element=0;element<elements.size();++element) {
        const auto& entry=elements[element];
        for(int id:entry.faces) {
            require(id>=0 && static_cast<std::size_t>(id)<faces.size(),"invalid original element-face incidence");
            const auto index=static_cast<std::size_t>(id);
            if(faces[index].kind!=ElectricPhaseSkeletonFace::Kind::Interface)continue;
            auto& owner=owners[index][entry.liquid?0:1];
            require(owner==missing,"interface has duplicate same-phase incidence");
            owner=element;
        }
    }
    for(std::size_t id=0;id<faces.size();++id) {
        const auto& face=faces[id];
        if(face.kind!=ElectricPhaseSkeletonFace::Kind::Interface)continue;
        const std::string context="geometry="+std::to_string(geometry.geometryIdentity().value())+
            " face="+std::to_string(id)+" ";
        require(owners[id][0]!=missing && owners[id][1]!=missing,
            "interface lacks original two-sided element incidence",context);
        require(std::isfinite(face.length) && face.length>0 && !face.nodes.empty(),
            "interface has unresolved positive arclength",context);
        const auto shared=source.skeletalFaceTrace(id);
        const auto minus=source.normalDisplacementPolynomialP2(id,true);
        const auto plus=source.normalDisplacementPolynomialP2(id,false);
        const auto eta=source.interfaceChargeMoments(id);
        require(finite(shared) && finite(minus) && finite(plus) && finite(eta),
            "nonfinite original electric trace or charge moment",context);
        std::array<std::array<Sum,3>,3> massSums;
        std::array<Sum,3> crossSums;
        std::vector<CompletionNode> nodes;nodes.reserve(face.nodes.size());
        for(const auto& node:face.nodes) {
            require(std::isfinite(node.point.x) && std::isfinite(node.point.y)
                && std::isfinite(node.coordinate) && node.coordinate>=-1 && node.coordinate<=1,
                "nonfinite same-source point or arclength coordinate",context);
            const double normalNorm=std::hypot(node.outwardNormal.x,node.outwardNormal.y);
            require(std::isfinite(normalNorm) && std::abs(normalNorm-1)<=1e-12,
                "same-source interface normal must be unit",context);
            require(std::isfinite(node.weight) && node.weight>0,
                "ordinary interface quadrature weight must be positive",context);
            CompletionNode q;
            q.basis=Electric::faceBasis(node.coordinate);q.p3=cubicBasis(node.coordinate);
            q.weight=node.weight*Electric::physicalWeight(node.point,geometry.axisymmetric());
            require(std::isfinite(q.weight) && q.weight>0,
                "physical interface quadrature weight must be positive",context);
            q.sharedPotential=valueP2(shared,q.basis);
            const double pMinus=source.potentialInElement(owners[id][0],node.point);
            const double pPlus=source.potentialInElement(owners[id][1],node.point);
            q.averagePotential=.5*pMinus+.5*pPlus; // Fixed symmetric omega, PDF §14.
            require(std::isfinite(q.sharedPotential) && std::isfinite(q.averagePotential),
                "unrepresentable same-field cubic trace samples",context);
            for(int j=0;j<3;++j) {
                crossSums[static_cast<std::size_t>(j)].add(q.weight*q.basis[j]*q.p3);
                for(int k=0;k<3;++k)
                    massSums[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)]
                        .add(q.weight*q.basis[j]*q.basis[k]);
            }
            nodes.push_back(std::move(q));
        }
        Eigen::Matrix3d mass;
        Electric::FacePolynomial cross,chargeMoments;
        for(int j=0;j<3;++j) {
            cross[j]=crossSums[static_cast<std::size_t>(j)].value();
            chargeMoments[j]=eta[static_cast<std::size_t>(j)];
            for(int k=0;k<3;++k)
                mass(j,k)=massSums[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)].value();
        }
        require(mass.allFinite() && cross.allFinite(),"nonfinite weighted interface mass",context);
        Eigen::LLT<Eigen::Matrix3d> factor(mass);
        require(factor.info()==Eigen::Success,"weighted interface mass is not positive definite",context);
        const auto diagonal=factor.matrixL().toDenseMatrix().diagonal().eval();
        require(diagonal.allFinite() && (diagonal.array()>0).all(),
            "weighted interface mass has invalid positive pivots",context);
        const Electric::FacePolynomial alpha=factor.solve(cross),charge=factor.solve(chargeMoments);
        require(alpha.allFinite() && charge.allFinite(),"weighted trace projection failed",context);
        // The original solve already accepts D_n^+-D_n^- in its full recovery
        // arithmetic scale. Preserve that receipt: do not add a roundoff /
        // roundoff physical-flux gate for a legitimate constant potential.
        std::array<Sum,3> chargeProjectionResidual;
        Sum chargeProjectionScale;
        for(int j=0;j<3;++j) {
            for(int k=0;k<3;++k) {
                const double term=mass(j,k)*charge[k];
                chargeProjectionResidual[static_cast<std::size_t>(j)].add(term);
                chargeProjectionScale.add(std::abs(term));
            }
            chargeProjectionResidual[static_cast<std::size_t>(j)].add(-chargeMoments[j]);
            chargeProjectionScale.add(std::abs(chargeMoments[j]));
        }
        double chargeProjectionError=0;
        for(int j=0;j<3;++j)
            chargeProjectionError=std::max(chargeProjectionError,relative(
                std::abs(chargeProjectionResidual[static_cast<std::size_t>(j)].value()),
                chargeProjectionScale.value(),context));
        Sum psiNorm,fit;
        for(auto& q:nodes) {
            q.psi=q.p3-q.basis.dot(alpha);
            psiNorm.add(q.weight*q.psi*q.psi);
            fit.add(q.weight*q.psi*(q.averagePotential-q.sharedPotential));
        }
        const double denominator=psiNorm.value(),numerator=fit.value();
        require(std::isfinite(denominator) && denominator>0 && std::isfinite(numerator),
            "cubic complementary mode has unresolved positive norm",context);
        const double completion=numerator/denominator;
        require(std::isfinite(completion),"unrepresentable cubic completion",context);
        auto& target=data->faces[id];target.available=true;
        for(int j=0;j<3;++j)
            target.potential[static_cast<std::size_t>(j)]=shared[static_cast<std::size_t>(j)]-completion*alpha[j];
        target.potential[3]=completion;
        const double derivativeScale=-2/face.length;
        target.tangentialField={{
            derivativeScale*(target.potential[1]+target.potential[3]),
            derivativeScale*3*target.potential[2],
            derivativeScale*5*target.potential[3]}};
        require(finite(target.potential) && finite(target.tangentialField),
            "unrepresentable completed potential or arclength derivative",context);
        std::array<Sum,3> projectionResidual;
        Sum projectionScale,optimalityResidual,optimalityScale,integratedX,integratedY;
        target.samples.reserve(face.nodes.size());
        const std::array<double,3> chargeCoefficients{{charge[0],charge[1],charge[2]}};
        target.normalDisplacement={{minus,plus}};
        target.charge=chargeCoefficients;
        for(std::size_t index=0;index<nodes.size();++index) {
            const auto& q=nodes[index];const auto& original=face.nodes[index];
            const double completed=valueP3(target.potential,q.basis,q.p3);
            const double dnMinus=valueP2(minus,q.basis),dnPlus=valueP2(plus,q.basis);
            const double density=valueP2(chargeCoefficients,q.basis);
            const double et=valueP2(target.tangentialField,q.basis);
            require(std::isfinite(completed) && std::isfinite(dnMinus) && std::isfinite(dnPlus)
                && std::isfinite(density) && std::isfinite(et),
                "unrepresentable completed electric field sample",context);
            for(int j=0;j<3;++j) {
                projectionResidual[static_cast<std::size_t>(j)].add(
                    q.weight*q.basis[j]*(completed-q.sharedPotential));
                // A full arithmetic trace family avoids dividing zero modes
                // by their own roundoff. No unit floor or altered equation.
                projectionScale.add(q.weight*std::abs(q.basis[j])*
                    (std::abs(completed)+std::abs(q.sharedPotential)));
            }
            optimalityResidual.add(q.weight*q.psi*(completed-q.averagePotential));
            optimalityScale.add(q.weight*std::abs(q.psi)*
                (std::abs(completed)+std::abs(q.averagePotential)));
            ElectricMaxwellSample sample;
            sample.point=original.point;sample.normal=original.outwardNormal;
            sample.tangent={-sample.normal.y,sample.normal.x};
            sample.coordinate=original.coordinate;sample.physicalWeight=q.weight;
            sample.normalDisplacement={{dnMinus,dnPlus}};
            sample.tangentialField=et;sample.chargeDensity=density;
            const double plusStress=dnPlus*dnPlus/epsilon[1];
            const double minusStress=dnMinus*dnMinus/epsilon[0];
            const double tangentStress=(epsilon[1]-epsilon[0])*et*et;
            require(std::isfinite(plusStress) && std::isfinite(minusStress) && std::isfinite(tangentStress),
                "Maxwell stress product overflows",context);
            sample.normalTraction=.5*(plusStress-minusStress-tangentStress);
            sample.tangentialTraction=density*et;
            require(std::isfinite(sample.normalTraction) && std::isfinite(sample.tangentialTraction),
                "Maxwell traction overflows",context);
            sample.tractionJump={
                sample.normalTraction*sample.normal.x+sample.tangentialTraction*sample.tangent.x,
                sample.normalTraction*sample.normal.y+sample.tangentialTraction*sample.tangent.y};
            require(std::isfinite(sample.tractionJump.x) && std::isfinite(sample.tractionJump.y),
                "Maxwell traction vector overflows",context);
            const double x=q.weight*sample.tractionJump.x,y=q.weight*sample.tractionJump.y;
            require(std::isfinite(x) && std::isfinite(y),"Maxwell traction integral product overflows",context);
            integratedX.add(x);integratedY.add(y);
            target.samples.push_back(std::move(sample));
        }
        double projectionError=chargeProjectionError;
        for(int j=0;j<3;++j) {
            projectionError=std::max(projectionError,relative(
                std::abs(projectionResidual[static_cast<std::size_t>(j)].value()),projectionScale.value(),context));
        }
        const double optimalityError=relative(
            std::abs(optimalityResidual.value()),optimalityScale.value(),context);
        require(projectionError<=controls.relativeProjectionTolerance,
            "completed cubic violates the original shared P2 projection",context);
        require(optimalityError<=controls.relativeProjectionTolerance,
            "completed cubic violates weighted one-mode optimality",context);
        target.integrated={integratedX.value(),integratedY.value()};
        require(std::isfinite(target.integrated.x) && std::isfinite(target.integrated.y),
            "Maxwell traction integral accumulation overflows",context);
        data->report.maxRelativeProjectionResidual=std::max(
            data->report.maxRelativeProjectionResidual,projectionError);
        data->report.maxRelativeOptimalityResidual=std::max(
            data->report.maxRelativeOptimalityResidual,optimalityError);
        // The electric geometry already accepted this precise GL8 rule.
        // Recover its piece prefix once for subsequent same-field point
        // queries, retaining every original node's accepted coordinate.
        require(face.component>=0 && face.graph>=0,"interface lacks original graph identity",context);
        const auto component=static_cast<std::size_t>(face.component),graph=static_cast<std::size_t>(face.graph);
        const auto& publication=geometry.publication();
        const auto& chain=publication.topology().curves.at(component);
        require(graph<chain.graphs().size(),"interface graph identity is out of range",context);
        const int refinement=geometry.report().acceptedRefinement;
        require(refinement>=0 && refinement<=8,"invalid accepted geometry refinement",context);
        const int pieces=1<<refinement;
        target.arclengthPrefix.reserve(static_cast<std::size_t>(pieces)+1);
        target.arclengthPrefix.push_back(0);
        long double complete=0;std::size_t originalIndex=0;
        for(int piece=0;piece<pieces;++piece) {
            const double begin=static_cast<double>(piece)/pieces,end=static_cast<double>(piece+1)/pieces;
            long double intervalLength=0;
            for(const auto& q:publication.quadrature(component,graph,begin,end)) {
                require(originalIndex<face.nodes.size(),"point-query replay exceeds original quadrature",context);
                const auto& original=face.nodes[originalIndex++];
                require(q.sample.geometryIdentity==geometry.geometryIdentity()
                    && sameClock(q.sample.timeLayer,geometry.timeLayer())
                    && q.sample.frame.point.x==original.point.x && q.sample.frame.point.y==original.point.y
                    && q.sample.frame.normal.x==original.outwardNormal.x
                    && q.sample.frame.normal.y==original.outwardNormal.y && q.weight==original.weight,
                    "point-query replay differs from original accepted quadrature",context);
                const double fraction=chain.graphFraction(static_cast<int>(graph),q.sample.frame.cellChartPoint);
                require(std::isfinite(fraction) && fraction>begin && fraction<end,
                    "point-query replay has unresolved graph fraction",context);
                target.coordinates.emplace(chartKey(q.sample.frame.cellChartPoint),original.coordinate);
                require(std::isfinite(q.weight) && q.weight>0,"invalid point-query arclength weight",context);
                intervalLength+=q.weight;
            }
            require(std::isfinite(intervalLength) && intervalLength>0,
                "point-query replay has unresolved positive piece length",context);
            complete+=intervalLength;target.arclengthPrefix.push_back(complete);
        }
        require(originalIndex==face.nodes.size() && static_cast<double>(complete)==face.length,
            "point-query replay lost original accepted arclength",context);
        require(data->graphFaces.emplace(std::make_pair(component,graph),id).second,
            "duplicate original electric interface graph",context);
        data->report.quadratureNodes+=target.samples.size();
        data->interface.push_back(id);
    }
    data->report.interfaceFaces=data->interface.size();
    return ElectricMaxwellTrace(std::move(data));
}
const ElectricPhaseGeometry& ElectricMaxwellTrace::geometry() const{return checked().source.geometry();}
Interface::HermiteGalsTimeLayer ElectricMaxwellTrace::timeLayer() const{return checked().source.timeLayer();}
Interface::HermiteGalsTimeLayer ElectricMaxwellTrace::geometrySourceTimeLayer() const
{return checked().source.geometry().timeLayer();}
const ElectricMaxwellReport& ElectricMaxwellTrace::report() const{return checked().report;}
const std::vector<std::size_t>& ElectricMaxwellTrace::interfaceFaces() const{return checked().interface;}
std::array<double,4> ElectricMaxwellTrace::potentialPolynomialP3(std::size_t face) const
{return checked().face(face).potential;}
std::array<double,3> ElectricMaxwellTrace::tangentialFieldPolynomialP2(std::size_t face) const
{return checked().face(face).tangentialField;}
const std::vector<ElectricMaxwellSample>& ElectricMaxwellTrace::quadrature(std::size_t face) const
{return checked().face(face).samples;}
ElectricMaxwellSample ElectricMaxwellTrace::sample(std::size_t component,std::size_t graph,
    const Interface::HermiteGalsSurfaceSample& query) const
{
    const auto& d=checked();const auto& geometry=d.source.geometry();
    const std::string context="geometry="+std::to_string(geometry.geometryIdentity().value())+
        " component="+std::to_string(component)+" graph="+std::to_string(graph)+" ";
    require(query.geometryIdentity==geometry.geometryIdentity()
        && sameClock(query.timeLayer,geometry.timeLayer()),"point-query geometry/full clock differs",context);
    const auto found=d.graphFaces.find(std::make_pair(component,graph));
    require(found!=d.graphFaces.end(),"point query requires an original electric interface graph",context);
    const auto& target=d.face(found->second);const auto& publication=geometry.publication();
    const auto& chain=publication.topology().curves.at(component);
    const auto& originalGraph=chain.graphs().at(graph);
    const auto& frame=query.frame;
    require(frame.nativeCellIndex==originalGraph.nativeCellIndex && frame.hasCellChartPoint
        && !frame.hasCellLocalPoint,"point-query lacks the original native graph chart",context);
    const auto chart=LevelSetHermiteGeometry::cellChartPoint(frame);
    const auto& scalar=publication.state().geometry().snapshot();
    require(chain.usesSource(scalar),"point-query graph uses another scalar",context);
    require(chain.graphAt(frame.nativeCellIndex,chart)==static_cast<int>(graph),
        "point-query chart belongs to another actual graph",context);
    // NativeGraphChain::frameAtInterval obtains exactly this same chart jet.
    // Authenticate it directly rather than round-tripping a displayed global
    // point (or re-solving a root from a rounded graph fraction).
    const auto actual=scalar.jetAtNativeCellChart(frame.nativeCellIndex,chart);
    require(frame.point.x==actual.point.x && frame.point.y==actual.point.y
        && frame.normal.x==actual.normal.x && frame.normal.y==actual.normal.y
        && frame.curvature==actual.curvature && frame.levelSetResidual==actual.levelSetResidual,
        "point-query frame differs from its actual source scalar jet",context);
    const double fraction=chain.graphFraction(static_cast<int>(graph),chart);
    require(std::isfinite(fraction) && fraction>=0 && fraction<=1,
        "point-query original graph fraction is unresolved",context);
    // Certify the retained graph query with the publication's original
    // interface/normal/RZ gates. Its round-trip chart is not substituted for
    // the caller's certified chart or used to identify a quadrature node.
    (void)publication.sample(component,graph,fraction);
    double coordinate=0;
    const auto cached=target.coordinates.find(chartKey(chart));
    if(cached!=target.coordinates.end())coordinate=cached->second;
    else {
        const int pieces=static_cast<int>(target.arclengthPrefix.size())-1;
        require(pieces>0,"point-query lacks accepted arclength pieces",context);
        const int piece=fraction==1?pieces-1:static_cast<int>(fraction*pieces);
        require(piece>=0 && piece<pieces,"point-query arclength piece is out of range",context);
        const double begin=static_cast<double>(piece)/pieces;
        long double partial=0;
        // Exact endpoints/piece boundaries have no partial arc to integrate.
        if(fraction>begin)
            for(const auto& prefix:publication.quadrature(component,graph,begin,fraction)) {
                require(std::isfinite(prefix.weight) && prefix.weight>0,
                    "invalid same-curve point-query arclength prefix",context);
                partial+=prefix.weight;
            }
        const long double complete=target.arclengthPrefix.back();
        coordinate=static_cast<double>(2*(target.arclengthPrefix[static_cast<std::size_t>(piece)]+partial)/complete-1);
    }
    require(std::isfinite(coordinate)
        && (fraction==0?coordinate==-1:fraction==1?coordinate==1:coordinate>-1 && coordinate<1),
        "point-query normalized original arclength is unresolved",context);
    const auto basis=Electric::faceBasis(coordinate);
    ElectricMaxwellSample result;
    result.point=frame.point;result.normal=frame.normal;result.tangent={-frame.normal.y,frame.normal.x};
    result.coordinate=coordinate;result.physicalWeight=0;
    result.normalDisplacement={{valueP2(target.normalDisplacement[0],basis),
        valueP2(target.normalDisplacement[1],basis)}};
    result.chargeDensity=valueP2(target.charge,basis);
    result.tangentialField=valueP2(target.tangentialField,basis);
    require(finite(result.normalDisplacement) && std::isfinite(result.chargeDensity)
        && std::isfinite(result.tangentialField),"unrepresentable same-field point query",context);
    const auto& epsilon=d.source.epsilonByPhase();
    const double dnMinus=result.normalDisplacement[0],dnPlus=result.normalDisplacement[1],et=result.tangentialField;
    const double plusStress=dnPlus*dnPlus/epsilon[1],minusStress=dnMinus*dnMinus/epsilon[0];
    const double tangentStress=(epsilon[1]-epsilon[0])*et*et;
    require(std::isfinite(plusStress) && std::isfinite(minusStress) && std::isfinite(tangentStress),
        "point-query Maxwell stress product overflows",context);
    result.normalTraction=.5*(plusStress-minusStress-tangentStress);
    result.tangentialTraction=result.chargeDensity*et;
    require(std::isfinite(result.normalTraction) && std::isfinite(result.tangentialTraction),
        "point-query Maxwell traction overflows",context);
    result.tractionJump={
        result.normalTraction*result.normal.x+result.tangentialTraction*result.tangent.x,
        result.normalTraction*result.normal.y+result.tangentialTraction*result.tangent.y};
    require(std::isfinite(result.tractionJump.x) && std::isfinite(result.tractionJump.y),
        "point-query Maxwell traction vector overflows",context);
    return result;
}
Interface::HermiteGalsSurfaceVectorView ElectricMaxwellTrace::surfaceVectorView() const
{
    const auto& d=checked();
    require(sameClock(d.source.timeLayer(),d.source.geometry().timeLayer()),
        "surface load view requires coeval electric/geometry full clocks");
    return Interface::HermiteGalsSurfaceVectorView(d.source.geometry().publication(),
        [trace=*this](std::size_t component,std::size_t graph,const Interface::HermiteGalsSurfaceSample& q) {
            return trace.sample(component,graph,q).tractionJump;
        });
}
Vector2D ElectricMaxwellTrace::integratedTractionJump(std::size_t face) const
{return checked().face(face).integrated;}
void ElectricMaxwellTrace::validateSource(const QuadtreeMesh& mesh,
    const Interface::HermiteGalsPublication& publication,const Interface::HermiteGalsTimeLayer& electricLayer) const
{
    const auto& source=checked().source;
    source.validateSource(mesh,publication);
    (void)electricLayer.time();
    require(sameClock(source.timeLayer(),electricLayer),"actual electric clock differs",
        "geometry="+std::to_string(source.geometry().geometryIdentity().value())+" ");
}
}
