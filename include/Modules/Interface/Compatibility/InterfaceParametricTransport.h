#pragma once
// S08_FROZEN_BASELINE_BEGIN
// Frozen experimental baseline: parametricGals curve transport and fluxes.
// Retained only for explicit A/B comparisons and legacy case regression.
// Do not add conservativeHermiteGals V2 features here or use this path as
// an implicit fallback. Keep V2 state and accepted geometry independently owned.
// Scope and source identities: docs/task/S08_Frozen_Baselines.json.
// S08_FROZEN_BASELINE_END
#include "Modules/Interface/Compatibility/InterfaceParametricProjection.h"
#include "Core/Base/Vector2D.h"
#include "Core/Geometry/UnitPolynomialRoots.h"

// Explicit parametric GALS producer. The flow and native flux callbacks must
// describe the same velocity over a step. Frozen and time-linear fluxes use
// the same space-time integration; no Hydro implementation is
// owned here. Complete rectangular XY/RZ native leaves and resolved regular
// space-time face charts, including temporal returns. Nonregular contacts
// use the same Green integral on shrinking material-time rectangles, with
// their unresolved contribution charged to the common flux error budget.
// Unresolved events beyond that budget are rejected. Fixed mixed leaves use
// their actual shared faces; dynamic mesh transfer is a separate operation.
namespace InterfaceParametricTransport {
using Candidate=InterfaceParametricProjection::Candidate;
using Curve=InterfaceParametricProjection::Curve;
using Source=InterfaceParametricProjection::Source;
struct MapJet {
    Vector2D point,dx,dy,timeDerivative,displacement;
    double coordinateResolution=0;
    MapJet()=default;
    // Retain the integrated increment; subtracting two absolute positions
    // cannot resolve a face-relative residual near a tangency.
    MapJet(Vector2D p,Vector2D x,Vector2D y,Vector2D t,Vector2D move,double resolution=0)
        :point(p),dx(x),dy(y),timeDerivative(t),displacement(move),coordinateResolution(resolution){}
};
// A fixed, read-only map for this step. Independent calls can execute
// concurrently and must return the same jet for identical arguments.
using Flow=std::function<MapJet(Vector2D,double)>;
struct NativeFlux {
    // Same physical grid and immutable native leaf identity as the inventory.
    InterfaceParametricProjection::Grid grid;
    int retainedMomentDegree=0;
    bool retainedVelocityProducts=false;
    // Internal weighted functional: the retained polynomial is the base q,
    // while this nonnegative multiplier (bounded by one) weights its work.
    // It does not change the trajectories or the inventory's zeroth flux.
    std::function<double(bool,int,int,double)> functionalWeight;
    int productFunctional=-1;
    int planeCount(bool horizontal) const {return grid.leaves?int(grid.leaves->planes[horizontal].size()):(horizontal?grid.ny:grid.nx)+1;}
    int stripCount(bool horizontal,int line) const {return grid.leaves?int(grid.leaves->planes[horizontal].at(line).faces.size()):(horizontal?grid.nx:grid.ny);}
    int faceOffset(bool horizontal,int line) const {return grid.leaves?grid.leaves->planes[horizontal].at(line).offset:line*stripCount(horizontal,line);}
    int faceCount(bool horizontal) const {const int last=planeCount(horizontal)-1;return faceOffset(horizontal,last)+stripCount(horizontal,last);}
    std::size_t cellCount() const {return grid.size();}
    bool exteriorPlane(bool horizontal,int line) const {return line==0||line==planeCount(horizontal)-1;}
    double planeCoordinate(bool horizontal,int line) const {
        return grid.leaves?grid.leaves->planes[horizontal].at(line).coordinate:
            (horizontal?grid.bottom+grid.height*line/grid.ny:grid.left+grid.width*line/grid.nx);
    }
    InterfaceParametricProjection::Grid::Face face(bool horizontal,int line,int strip) const {
        if(grid.leaves)return grid.leaves->planes[horizontal].at(line).faces.at(strip);
        const double begin=horizontal?grid.left:grid.bottom,h=horizontal?grid.width/grid.nx:grid.height/grid.ny;
        return {begin+strip*h,begin+(strip+1)*h,
            line>0?(horizontal?strip*grid.ny+line-1:(line-1)*grid.ny+strip):-1,
            line<(horizontal?grid.ny:grid.nx)?(horizontal?strip*grid.ny+line:line*grid.ny+strip):-1,-1};
    }
    // Use a neighboring physical cell capacity for the common plane budget.
    // Horizontal planes share a budget no larger than the smallest RZ cell.
    double planeCapacity(bool horizontal,int line) const {
        if(grid.leaves)return grid.leaves->planes[horizontal].at(line).capacity;
        return faceCapacity(horizontal,line,0);
    }
    double faceCapacity(bool horizontal,int line,int strip) const {
        double result=std::numeric_limits<double>::infinity();const auto f=face(horizontal,line,strip);
        for(int slot:{f.negative,f.positive})if(slot>=0)result=std::min(result,grid.capacity(slot));return result;
    }
    // +coordinate normal; native line/strip indices. Positions are physical.
    // Read-only callbacks for the same step; independent faces may query
    // concurrently, just as the trajectory callback above does.
    std::function<double(bool,int,int,double,double)> integral;
    std::function<double(bool,int,int,double)> density;
    // Bound of the physical |density| over the entire native face.
    // Used to account for unresolved vertex/tangency intervals in flux units.
    std::function<double(bool,int,int)> densityMagnitudeBound;
    // Optional retained density on each native HALF face, with respect to
    // its unit parameter (thus already including the half-face measure).
    // Required to publish gross directional volumes, never fitted to samples.
    std::function<std::vector<Quantum::Geometry::UnitPolynomialTerm>(bool,int,int,int)> polynomial;
};
// Endpoint traces of q(s,theta)=(1-theta)q_begin(s)+theta*q_end(s).
// theta is normalized elapsed time on the SAME interval as Flow. The
// frozen overloads below explicitly supply two identical endpoint traces.
struct TimeLinearFlux {
    NativeFlux begin,end;
    double spatialIntegral(bool horizontal,int line,int strip,double a,double b,double theta) const {
        const double first=begin.integral(horizontal,line,strip,a,b),last=end.integral(horizontal,line,strip,a,b);
        return first+(last-first)*theta;
    }
    // Normalized time integral; retain signed limits for returning crossings.
    double timeIntegral(bool horizontal,int line,int strip,double a,double b,double lo,double hi) const {
        return (hi-lo)*spatialIntegral(horizontal,line,strip,a,b,.5*(lo+hi));
    }
};
struct SourceReport {
    double positionApproximation=0,tangentApproximation=0,sourceJoin=0,sourceTangentJoin=0;
    std::size_t maximumSourcePieces=0;int subdivisions=0;
};
struct FluxReport {
    int crossingIntervals=0,vertexEvents=0,tangentTurns=0,quadratureRefinements=0;
    int spaceTimePatches=0,timeCharts=0,materialCharts=0,chartSubdivisions=0;
    int boundedVertexIntervals=0;
    int boundedContactPatches=0,contactEdgeRefinements=0;
    long long mapQueries=0,mapCacheHits=0;
    double crossingResidual=0,quadratureError=0,exteriorFlux=0;
    double sheetPositionError=0;
    double vertexIntervalError=0;
    double contactIntegralError=0;
    double directionalQuadratureError=0,directionalBoundAdjustment=0,momentQuadratureError=0;
    // Sum of face integration errors / their own adjacent-cell capacities.
    // Physical-volume fields above remain conservative absolute bounds.
    double quadratureFractionError=0,directionalFractionError=0;
    double totalFluxFractionError=0;
};
// Surface measure transported across the physical box, per unit constant
// surface tension. Positive advection is outward. endpointWork is the time
// integral of the conormal term already present in weak curvature traction.
// Their difference belongs in the open-domain surface-energy balance.
struct SurfaceFlux {
    double advection=0,endpointWork=0,quadratureError=0,criticalError=0,budget=0;
    int timeCharts=0,materialCharts=0,refinements=0,criticalPatches=0;
};
struct FaceVolumes {
    static constexpr int productDegree=20;
    using Products=std::array<double,productDegree+1>;
    // B_k^20(xi), times (r_min/r)^2 on off-axis horizontal RZ faces.
    // Axis-touching traces use their analytically regular polynomial form.
    // All weights lie in [0,1]; the original zeroth moment remains separate.
    static double productWeight(int k,double xi,bool horizontal,bool rz,double lower,double upper){
        if(k<0||k>productDegree||xi<0||xi>1||!(upper>lower))throw std::invalid_argument("invalid face product chart");
        double coefficient=1;for(int j=1;j<=k;j++)coefficient*=double(productDegree-j+1)/j;
        double value=coefficient*std::pow(xi,k)*std::pow(1-xi,productDegree-k);
        if(horizontal&&rz&&lower>0)value*=std::pow(lower/(lower+(upper-lower)*xi),2);return value;
    }
    static double productWeightBound(int k,double a,double b,bool horizontal,bool rz,double lower,double upper){
        const double peak=std::clamp(double(k)/productDegree,a,b);
        double value=productWeight(k,peak,false,false,0,1);
        if(horizontal&&rz&&lower>0)value*=std::pow(lower/(lower+(upper-lower)*a),2);
        return value==0?0:std::nextafter(value*(1+64*std::numeric_limits<double>::epsilon()),std::numeric_limits<double>::infinity());
    }
    // Coefficients representing the constant function in that SAME basis.
    static Products constantProduct(bool horizontal,bool rz,double lower,double upper){Products out{};
        const double b=horizontal&&rz&&lower>0?(upper-lower)/lower:0;
        for(int k=0;k<=productDegree;k++)out[k]=1+2*b*k/productDegree+b*b*k*(k-1)/(productDegree*(productDegree-1));return out;
    }
    struct Face {
        double total=0,liquid=0;
        double totalPositive=0,totalNegative=0,liquidPositive=0,liquidNegative=0;
        // Integral of xi*q and xi^2*q on the SAME face-time wet region,
        // xi in [0,1] along the entire face. total/liquid above remain the
        // unique zeroth moments used to advance the inventory.
        std::array<double,2> totalHigher{},liquidHigher{};
        // Positive-direction moments on physical exterior faces, from the
        // same trace and wet region. Negative = positive - signed.
        std::array<double,2> totalPositiveHigher{},liquidPositiveHigher{};
        Products totalProducts{},liquidProducts{},totalPositiveProducts{},liquidPositiveProducts{};
    };
    // +coordinate orientation, NativeFlux::faceOffset+strip, vertical then horizontal.
    std::array<std::vector<Face>,2> faces;
    bool hasDirections=false;
    bool hasExteriorDirectionMoments=false;
    int momentDegree=0;
    bool hasVelocityProducts=false;
    SurfaceFlux surface;
};
Source transportSource(const Curve&,const InterfaceParametricProjection::Grid&,double dt,const Flow&,SourceReport&);
// Returns a new ledger; on any failure the given inventory remains unchanged.
std::vector<double> advanceInventory(const Curve&,const Source&,const Flow&,const NativeFlux&,
    double dt,const std::vector<double>& inventory,FluxReport&);
std::vector<double> advanceInventory(const Curve&,const Source&,const Flow&,const TimeLinearFlux&,
    double dt,const std::vector<double>& inventory,FluxReport&,FaceVolumes* retained=nullptr);
// No negative phase: zero wet region, while the same native traces still
// supply total directional, higher and complete velocity-product fluxes.
void backgroundFlux(const TimeLinearFlux&,double dt,FluxReport&,FaceVolumes&);
struct Step {Candidate candidate;SourceReport transport;FluxReport flux;double inventoryBalance=0;FaceVolumes faceVolumes;};
using Observer=std::function<void(const Source&,const Candidate&)>;
Step prepare(const Candidate& old,double dt,const Flow&,const NativeFlux&,
    const InterfaceParametricProjection::Options&,const Observer& observer={});
Step prepare(const Candidate& old,double dt,const Flow&,const TimeLinearFlux&,
    const InterfaceParametricProjection::Options&,const Observer& observer={});
} // namespace InterfaceParametricTransport
