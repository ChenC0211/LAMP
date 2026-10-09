#include "Modules/Hydrodynamics/HydroCompatibleMomentBoundary.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianMeasure.h"

namespace Quantum::Hydrodynamics {
void HydroCompatibleMomentBoundary::validate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space) const{
    space.validate(mesh,domain);
    if(!identity_||identity_!=domain.meshIdentity()||faces_.size()!=std::size_t(domain.nativeFaceCount()))
        throw std::invalid_argument("moment boundary belongs to another native domain");
}
HydroCompatibleMomentBoundary HydroCompatibleMomentBoundary::build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentViscosity& stress,const std::vector<HydroBoundaryPatch>& patches){
    space.validate(mesh,domain);stress.validate(mesh,domain,space);
    validateHydroBoundaryPatches(patches,mesh.minCoord,mesh.maxCoord,
        {(mesh.maxCoord.x-mesh.minCoord.x)/mesh.Nx,(mesh.maxCoord.y-mesh.minCoord.y)/mesh.Ny},domain.axisymmetric(),true);
    HydroCompatibleMomentBoundary result;result.identity_=domain.meshIdentity();result.generalTrace_=!patches.empty();result.faces_.resize(domain.nativeFaceCount());
    for(const auto& patch:patches){result.open_=result.open_||patch.type!=HydroBoundaryType::FreeSlip;
        result.pressureOutlet_=result.pressureOutlet_||patch.type==HydroBoundaryType::PressureOutlet;}
    result.fixed_.assign(stress.matrix().rows(),false);result.values_=result.traction_=Eigen::VectorXd::Zero(stress.matrix().rows());
    int outlets=0;
    for(const auto face:domain.faces(mesh))if(face.outwardSign&&face.area>0){
        auto& saved=result.faces_[face.nativeRow];const int normal=space.faces()[face.nativeRow].first,tangent=stress.firstTangent(face.nativeRow);
        const auto side=face.horizontal?(face.outwardSign<0?HydroBoundarySide::South:HydroBoundarySide::North)
            :(face.outwardSign<0?HydroBoundarySide::West:HydroBoundarySide::East);
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const double lo=face.horizontal?a.x:a.y,hi=face.horizontal?b.x:b.y;
        for(const auto& patch:patches)if(patch.side==side&&patch.start<hi&&patch.end>lo){
            if(patch.start>lo||patch.end<hi)throw std::invalid_argument("moment patch cuts through a native boundary face");
            saved.kind=patch.type==HydroBoundaryType::NormalVelocity?Kind::NormalVelocity:
                (patch.type==HydroBoundaryType::FreeSlip?Kind::FreeSlip:Kind::PressureOutlet);
            if(saved.kind==Kind::NormalVelocity)saved.outwardVelocity=patch.value;
            else if(saved.kind==Kind::PressureOutlet)saved.pressure=patch.value;
            saved.incomingPhase=patch.incomingPhase;saved.incomingVelocity=patch.incomingVelocity;
            if(saved.kind==Kind::NormalVelocity&&!saved.incomingVelocity){Vector2D v{0,0};if(face.horizontal)v.y=face.outwardSign*patch.value;else v.x=face.outwardSign*patch.value;saved.incomingVelocity=v;}
        }
        if(saved.kind==Kind::PressureOutlet){
            // Constant full traction pairs with the SAME weighted face mean.
            result.traction_[normal]=-face.outwardSign*saved.pressure*face.area;++outlets;
        }else{
            result.fixed_[normal]=result.fixed_[normal+1]=true;
            result.values_[normal]=face.outwardSign*saved.outwardVelocity;
            const auto& f=space.faces()[face.nativeRow];
            result.normalReference_.push_back({face.horizontal,face.horizontal?a.y:a.x,lo,hi,f.area,f.tangentCenter,f.tangentScale,result.values_[normal],saved.kind});
        }
        if(saved.kind==Kind::Wall)result.fixed_[tangent]=result.fixed_[tangent+1]=true;
    }
    if(result.open_&&outlets==0)throw std::invalid_argument("moment patches contain no positive-area pressure outlet");
    result.assembleNormalConstraints(mesh,domain,space);
    return result;
}
void HydroCompatibleMomentBoundary::assembleNormalConstraints(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space){
    std::vector<Eigen::Triplet<double>> entries;normalValues_=Eigen::VectorXd::Zero(2*normalReference_.size());nestedNormal_=true;
    for(int k=0;k<int(normalReference_.size());k++){const auto& reference=normalReference_[k];double covered=0;normalValues_[2*k]=reference.value;
        for(const auto face:domain.faces(mesh))if(face.outwardSign&&face.area>0&&face.horizontal==reference.horizontal){
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            if((face.horizontal?a.y:a.x)!=reference.plane)continue;
            const double lo=face.horizontal?a.x:a.y,hi=face.horizontal?b.x:b.y;
            if(hi<=reference.lo||lo>=reference.hi)continue;
            if(faces_[face.nativeRow].kind!=reference.kind)throw std::invalid_argument("normal boundary physical kind changed");
            if(lo<reference.lo||hi>reference.hi){
                // A coarse face can cover several original weak data
                // intervals. Its first two moments cannot express those
                // subinterval tests: the complete field integrates them.
                nestedNormal_=false;const double a=std::max(lo,reference.lo),b=std::min(hi,reference.hi);
                covered+=(b-a)*(domain.axisymmetric()?CartesianMeasure::twoPi*(face.horizontal?.5*(a+b):reference.plane):1.);continue;
            }
            const auto& f=space.faces()[face.nativeRow];const double ratio=f.area/reference.area;
            entries.emplace_back(2*k,f.first,ratio);
            entries.emplace_back(2*k+1,f.first,ratio*(f.tangentCenter-reference.center)/reference.scale);
            entries.emplace_back(2*k+1,f.first+1,ratio*f.tangentScale/reference.scale);
            covered+=f.area;
        }
        if(std::abs(covered-reference.area)>256*std::numeric_limits<double>::epsilon()*reference.area)
            throw std::invalid_argument("normal boundary reference is not covered by target faces");
    }
    normalConstraints_.resize(normalValues_.size(),space.mass().rows());normalConstraints_.setFromTriplets(entries.begin(),entries.end());
}
HydroCompatibleMomentBoundary HydroCompatibleMomentBoundary::inheritNormalConstraints(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentBoundary& source) const{
    validate(mesh,domain,space);if(open_!=source.open_||pressureOutlet_!=source.pressureOutlet_||generalTrace_!=source.generalTrace_)
        throw std::invalid_argument("normal boundary kind changed during refinement");
    auto result=*this;result.normalReference_=source.normalReference_;result.assembleNormalConstraints(mesh,domain,space);return result;
}
}
