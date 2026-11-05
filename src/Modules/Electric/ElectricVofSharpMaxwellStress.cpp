#include "Modules/Electric/ElectricVofSharpMaxwellStress.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Electric
{
namespace
{
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
void require(bool ok,const char* reason)
{if(!ok)throw std::invalid_argument(std::string("[ElectricVofSharpMaxwellStress] ")+reason);}
double rounded(const FV::Dyadic& value)
{const auto x=FV::roundDouble(value).value;require(std::isfinite(x),"unrepresentable derived tensor or traction");return x;}
ElectricVofSharpMaxwellTensor tensor(Vector2D e,double epsilon)
{
    require(std::isfinite(e.x)&&std::isfinite(e.y)&&std::isfinite(epsilon)&&epsilon>0,"invalid electric sample/material");
    auto xx=FV::Dyadic::product({.5,epsilon,e.x,e.x});xx.add(FV::Dyadic::product({-.5,epsilon,e.y,e.y}));
    auto yy=FV::Dyadic::product({.5,epsilon,e.y,e.y});yy.add(FV::Dyadic::product({-.5,epsilon,e.x,e.x}));
    auto hoop=FV::Dyadic::product({-.5,epsilon,e.x,e.x});hoop.add(FV::Dyadic::product({-.5,epsilon,e.y,e.y}));
    return {rounded(xx),rounded(FV::Dyadic::product({epsilon,e.x,e.y})),rounded(yy),rounded(hoop)};
}
FV::Dyadic traction(const ElectricVofSharpMaxwellTensor& a,Vector2D v,Vector2D n,FV::Dyadic* absoluteTerms=nullptr)
{
    FV::Dyadic result;
    const auto add=[&](std::initializer_list<long double> factors){const auto term=FV::Dyadic::product(factors);
        result.add(term);if(absoluteTerms)absoluteTerms->add(term.absolute());};
    add({v.x,a.xx,n.x});add({v.x,a.xy,n.y});add({v.y,a.xy,n.x});add({v.y,a.yy,n.y});return result;
}
double relativeResidual(const FV::Dyadic& residual,const FV::Dyadic& scale)
{
    if(scale.compare(FV::Dyadic{})==0){require(residual.compare(FV::Dyadic{})==0,"nonzero stress defect at zero physical scale");return 0;}
    const auto ratio=residual.absolute().rounded()/scale.rounded();require(std::isfinite(ratio),"unrepresentable stress residual");return ratio;
}
}
struct ElectricVofSharpMaxwellStress::Data
{
    ElectricVofSharpField field;
    explicit Data(ElectricVofSharpField f):field(std::move(f)){}
};
ElectricVofSharpMaxwellStress::ElectricVofSharpMaxwellStress(std::shared_ptr<const Data> d):data_(std::move(d)){}
ElectricVofSharpMaxwellStress ElectricVofSharpMaxwellStress::prepare(const QuadtreeMesh& mesh,const ElectricVofSharpField& field)
{
    field.geometry().validateFor(mesh,field.geometry().source());
    field.chargeCandidate().validateFor(mesh,field.geometry(),field.candidateTime());
    return ElectricVofSharpMaxwellStress(std::make_shared<Data>(field));
}
const ElectricVofSharpField& ElectricVofSharpMaxwellStress::field() const
{require(bool(data_),"empty stress owner");return data_->field;}
bool ElectricVofSharpMaxwellStress::axisymmetric() const{return field().geometry().source().domain().axisymmetric();}
ElectricVofSharpMaxwellSample ElectricVofSharpMaxwellStress::sample(int cell,unsigned phase,Vector2D point) const
{return sampleAt(cell,phase,{{point.x,point.y},{0,0}});}
ElectricVofSharpMaxwellSample ElectricVofSharpMaxwellStress::sampleAt(int cell,unsigned phase,
    Geometry::PlanarCellCut::AnchoredPoint point) const
{
    require(phase<2,"invalid phase");const auto e=field().sampleElectricFieldAt(cell,phase,point);
    return {e,tensor(e,field().materials().epsilon[phase])};
}
ElectricVofSharpMaxwellInterface ElectricVofSharpMaxwellStress::interfaceSample(std::size_t arc,double parameter) const
{
    ElectricVofSharpMaxwellInterface out;out.field=field().interfaceField(arc,parameter);
    const auto& e=out.field;const auto& eps=field().materials().epsilon;
    out.tensor={tensor(e.electricField[0],eps[0]),tensor(e.electricField[1],eps[1])};
    auto normal=traction(out.tensor[1],e.normal,e.normal);normal.add(traction(out.tensor[0],e.normal,e.normal),-1);
    FV::Dyadic tangentScale;
    auto tangent=traction(out.tensor[1],e.tangent,e.normal,&tangentScale);
    tangent.add(traction(out.tensor[0],e.tangent,e.normal,&tangentScale),-1);
    FV::Dyadic expectedNormal,normalScale;
    for(unsigned phase=0;phase<2;++phase){const double sign=phase?1:-1;
        auto n=FV::Dyadic::product({sign,.5,eps[phase],e.normalElectricField[phase],e.normalElectricField[phase]});
        auto t=FV::Dyadic::product({-sign,.5,eps[phase],e.tangentialElectricField,e.tangentialElectricField});
        expectedNormal.add(n);expectedNormal.add(t);normalScale.add(n.absolute());normalScale.add(t.absolute());}
    auto expectedTangent=FV::Dyadic::product({e.surfaceChargeDensity,e.tangentialElectricField});
    tangentScale.add(expectedTangent.absolute());
    for(unsigned phase=0;phase<2;++phase)tangentScale.add(FV::Dyadic::product({eps[phase],e.normalElectricField[phase],e.tangentialElectricField}).absolute());
    normalScale.add(normal.absolute());auto nr=normal;nr.add(expectedNormal,-1);auto tr=tangent;tr.add(expectedTangent,-1);
    out.normalJump=rounded(normal);out.tangentialJump=rounded(tangent);
    out.expectedNormalJump=rounded(expectedNormal);out.expectedTangentialJump=rounded(expectedTangent);
    out.relativeNormalIdentityResidual=relativeResidual(nr,normalScale);
    out.relativeTangentialJumpResidual=relativeResidual(tr,tangentScale);
    require(out.relativeNormalIdentityResidual<=1e-12,"normal tensor identity violates original stress gate");
    require(out.relativeTangentialJumpResidual<=1e-10,"charge/tangential traction jump violates original electric gate");return out;
}
bool ElectricVofSharpMaxwellStress::sameInstance(const ElectricVofSharpMaxwellStress& other)const noexcept{return data_==other.data_;}
void ElectricVofSharpMaxwellStress::validateFor(const QuadtreeMesh& mesh,const ElectricVofSharpField& source)const
{
    require(field().sameInstance(source),"foreign field owner");source.geometry().validateFor(mesh,source.geometry().source());
    source.chargeCandidate().validateFor(mesh,source.geometry(),source.candidateTime());
}
}
