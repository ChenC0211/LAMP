#include "Interface/Setup/InterfaceRunSetup.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>

using Config=InterfaceTransportConfig;
InterfacePrescribedData::InterfacePrescribedData(const Config& c)
    :initial_(c.initialPhi,c.constants,c.radial?LAMP::ScalarExpression::Coordinates::Meridional:LAMP::ScalarExpression::Coordinates::Planar),
     origin_(c.velocityOrigin),scale_(c.velocityScale)
{
    c.validate();using Expression=LAMP::ScalarExpression;
    if(c.prescribedBoundary)boundary_=std::make_shared<const Expression>(c.boundaryPhi,c.constants,
        c.radial?Expression::Coordinates::MeridionalTime:Expression::Coordinates::PlanarTime);
    if(c.velocity==Config::Velocity::Polynomial){
        if(!(std::isfinite(scale_.x) && std::isfinite(scale_.y) && scale_.x>0 && scale_.y>0))throw std::invalid_argument("invalid polynomial scales");
        for(int k=0;k<16;++k){coefficients_[k][0]=std::make_shared<const Expression>(c.coefficients[k].x,c.constants,Expression::Coordinates::Time);
            coefficients_[k][1]=std::make_shared<const Expression>(c.coefficients[k].y,c.constants,Expression::Coordinates::Time);}
    }
}
LAMP::ScalarExpression::Jet InterfacePrescribedData::initial(Vector2D p) const{return initial_.sample(p.x,p.y);}
LAMP::ScalarExpression::Jet InterfacePrescribedData::boundary(Vector2D p,double time) const
{if(!boundary_)throw std::logic_error("no prescribed boundary data");return boundary_->sample(p.x,p.y,time);}
LAMP::Geometry::PolynomialVelocityData InterfacePrescribedData::polynomial(double previous,double arrival) const
{
    LAMP::Geometry::PolynomialVelocityData data;data.origin=origin_;data.scale=scale_;
    for(int k=0;k<16;++k){if(!coefficients_[k][0])throw std::logic_error("no polynomial velocity data");
        data.previous[k]={coefficients_[k][0]->sample(previous).value,coefficients_[k][1]->sample(previous).value};
        data.arrival[k]={coefficients_[k][0]->sample(arrival).value,coefficients_[k][1]->sample(arrival).value};}
    return data;
}
std::unique_ptr<QuadtreeMesh> makeInterfaceTransportMesh(const Config& c,bool amr)
{
    DictionaryMesh input;input.minCoord=c.domainMin;input.maxCoord=c.domainMax;input.Nx=input.Ny=c.resolution;input.maxInitialIterations=0;
    input.maxCellNum=std::max(16384,16*c.resolution*c.resolution);
    auto mesh=std::make_unique<QuadtreeMesh>(input);
    if(amr){DictionaryMesh::RefineBox box;box.level=c.refinementLevel;box.min=c.refineMin;box.max=c.refineMax;mesh->applyRefineBoxes({box});}
    return mesh;
}
LevelSetBoundaryTrace::Domain interfaceTransportTraceDomain(const QuadtreeMesh& mesh,bool radial)
{
    Vector2D h=mesh.maxCoord-mesh.minCoord;for(int c=0;c<mesh.cellIndexNum();++c){const auto size=mesh.cells(mesh.cellIndexToLabel[c]).size();h.x=std::min(h.x,size.x);h.y=std::min(h.y,size.y);}
    return {mesh.minCoord,mesh.maxCoord,h,radial};
}
