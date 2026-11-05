#include "Modules/Electric/ElectricVofRefinement.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Electric
{
namespace
{
void require(bool ok,const char* why)
{ if(!ok)throw std::invalid_argument(std::string("[ElectricVofRefinement] ")+why); }
struct CellSample
{
    long double density=0,x=0,y=0,hx=0,hy=0;
};
}
ElectricVofRefinementIndicators::ElectricVofRefinementIndicators(
    const Interface::GeometricVofPublication& geometry,const ElectricVofCharge& charge,
    std::vector<long double> jump,long double maximum,std::size_t faces)
    :geometry_(geometry),charge_(charge),densityJump_(std::move(jump)),
      maximumDensityJump_(maximum),interiorFaces_(faces){}
void ElectricVofRefinementIndicators::validateFor(const QuadtreeMesh& mesh,
    const Interface::GeometricVofPublication& geometry,const ElectricVofCharge& charge) const
{
    require(geometry_.identity()==geometry.identity()&&charge_.sameInstance(charge),
        "indicator requires its actual geometry and charge owners");
    charge.validateFor(mesh,geometry,geometry.timeLayer());
    require(geometry.timeLayer().position==LAMP::TimeLayer::Position::End,
        "indicator requires the complete End clock");
}
ElectricVofRefinementIndicators ElectricVofRefinement::sample(const QuadtreeMesh& mesh,
    const Interface::GeometricVofPublication& geometry,const ElectricVofCharge& charge)
{
    charge.validateFor(mesh,geometry,geometry.timeLayer());
    require(geometry.timeLayer().position==LAMP::TimeLayer::Position::End,
        "indicator requires the complete End clock");
    const auto& domain=geometry.domain();
    const auto count=static_cast<std::size_t>(mesh.cellIndexNum());
    const auto& rows=domain.cellRowByIndex();const auto& indices=domain.cellIndexByRow();
    const auto& volume=domain.cellVolumes();const auto& q=charge.integratedCharge();
    require(rows.size()==count&&q.size()==count&&volume.size()==indices.size(),
        "charge, selected maps and actual volumes have different native coverage");
    std::vector<CellSample> cells(count);
    std::vector<long double> jump(count);
    for(std::size_t native=0;native<count;++native){
        const int row=rows[native];
        if(row<0){require(q[native]==0,"excluded native cell carries nonzero charge");continue;}
        require(static_cast<std::size_t>(row)<indices.size()&&indices[row]==static_cast<int>(native),
            "selected native-cell map is not invertible");
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
        const auto center=cell.centerCoord(),size=cell.size();auto& value=cells[native];
        value.hx=size.x;value.hy=size.y;value.x=center.x;value.y=center.y;
        require(std::isfinite(value.hx)&&value.hx>0&&std::isfinite(value.hy)&&value.hy>0&&
            std::isfinite(volume[row])&&volume[row]>0&&std::isfinite(q[native]),
            "nonfinite charge, native size or actual selected volume");
        if(domain.axisymmetric()){
            require(value.x>0,"RZ leaf volume centroid requires positive cell-centre radius");
            value.x+=value.hx*(value.hx/value.x/12);
        }
        value.density=static_cast<long double>(q[native])/volume[row];
        require(std::isfinite(value.x)&&std::isfinite(value.y)&&std::isfinite(value.density),
            "native volume centroid or physical charge density is not finite");
    }
    long double maximum=0;std::size_t faces=0;
    for(const auto face:domain.faces(mesh)){
        if(face.negativeCell<0||face.positiveCell<0)continue;
        const auto a=static_cast<std::size_t>(face.negativeCell),b=static_cast<std::size_t>(face.positiveCell);
        require(a<count&&b<count&&rows[a]>=0&&rows[b]>=0,
            "actual internal subface lacks two selected native cells");
        const auto& left=cells[a];const auto& right=cells[b];
        const long double distance=std::abs(face.horizontal?right.y-left.y:right.x-left.x);
        const long double difference=std::abs(right.density-left.density);
        require(std::isfinite(distance)&&distance>0&&std::isfinite(difference),
            "nonfinite charge jump or nonpositive normal volume-centroid distance");
        const long double leftJump=((face.horizontal?left.hy:left.hx)/distance)*difference;
        const long double rightJump=((face.horizontal?right.hy:right.hx)/distance)*difference;
        require(std::isfinite(leftJump)&&std::isfinite(rightJump),
            "physical charge-density jump indicator is not finite");
        jump[a]=std::max(jump[a],leftJump);jump[b]=std::max(jump[b],rightJump);
        maximum=std::max({maximum,leftJump,rightJump});++faces;
    }
    return ElectricVofRefinementIndicators(geometry,charge,std::move(jump),maximum,faces);
}
}
