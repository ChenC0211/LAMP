#include "Modules/Hydrodynamics/HydroConservativePressure.h"
#include "Core/Base/Interpolation.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

namespace Quantum::Hydrodynamics
{
namespace {
void require(bool ok,const char* message)
{ if(!ok)throw std::invalid_argument(std::string("[HydroConservativePressure] ")+message); }
const Cell& cellAt(const QuadtreeMesh& mesh,int i) { return mesh.cells(mesh.cellIndexToLabel[i]); }
using Trace=std::map<int,double>;
struct Reconstruction { Vector2D center; std::array<Trace,2> gradient; };
Reconstruction reconstruct(const QuadtreeMesh& mesh,int owner)
{
    const auto& cell=cellAt(mesh,owner);const auto center=cell.centerCoord(),h=cell.size();
    std::vector<int> cells{owner};
    for(const auto& side:cell.neiCell)for(int label:side)if(label>=0)cells.push_back(mesh.cells(label).index());
    std::sort(cells.begin(),cells.end());cells.erase(std::unique(cells.begin(),cells.end()),cells.end());
    require(cells.size()>=3,"cell pressure reconstruction needs two-dimensional native support");
    std::vector<Vector2D> points;
    for(int c:cells) { const auto p=cellAt(mesh,c).centerCoord();points.push_back({(p.x-center.x)/h.x,(p.y-center.y)/h.y}); }
    long double rank=0;
    for(std::size_t i=1;i<points.size();++i)for(std::size_t j=i+1;j<points.size();++j)
        rank=std::max(rank,std::abs((static_cast<long double>(points[i].x)-points[0].x)*(points[j].y-points[0].y)
            -(static_cast<long double>(points[j].x)-points[0].x)*(points[i].y-points[0].y)));
    require(rank>0,"singular native pressure reconstruction");
    const auto value=Interpolation::leastSquareInterpolationCoeff(points,{0,0});
    const std::array<std::vector<double>,2> shifted{{Interpolation::leastSquareInterpolationCoeff(points,{1,0}),
        Interpolation::leastSquareInterpolationCoeff(points,{0,1})}};
    Reconstruction result;result.center=center;
    for(int axis=0;axis<2;++axis) {
        std::array<long double,2> moment{};
        for(std::size_t i=0;i<cells.size();++i) {
            const double w=shifted[axis][i]-value[i];require(std::isfinite(w),"nonfinite pressure trace");
            moment[0]+=w*points[i].x;moment[1]+=w*points[i].y;
            if(cells[i]!=owner) {
                const double derivative=w/(axis?h.y:h.x);
                result.gradient[axis][cells[i]]+=derivative;
                result.gradient[axis][owner]-=derivative;
            }
        }
        for(int k=0;k<2;++k)require(std::abs(moment[k]-(k==axis?1.L:0.L))<1e-10L,
            "native pressure trace lost affine consistency");
    }
    return result;
}
Vector2D weightedPoint(Vector2D a,Vector2D b,bool rz)
{
    Vector2D result{(a.x+b.x)/2,(a.y+b.y)/2};
    if(rz && a.x!=b.x)result.x=(2./3.)*(a.x*a.x+a.x*b.x+b.x*b.x)/(a.x+b.x);
    return result;
}
}

HydroConservativePressure HydroConservativePressure::topology(const QuadtreeMesh& mesh,const NativeLeafDomain& domain)
{
    domain.validateMeshIdentity(mesh);
    require(domain.cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum()),"requires the complete native fluid domain");
    HydroConservativePressure result;result.identity_=domain.meshIdentity();result.axisymmetric_=domain.axisymmetric();
    result.duals_=HydroFittedDetail::MomentumDuals::build(mesh,domain);
    const int count=domain.nativeFaceCount();
    result.boundaryCoefficient_.assign(count,0);result.boundaryPoints_.resize(count);result.anchors_.resize(count);
    std::vector<Eigen::Triplet<double>> incidence;
    for(const auto face:domain.faces(mesh)) {
        result.anchors_[face.nativeRow]=std::max(face.negativeCell,face.positiveCell);
        if(face.negativeCell>=0)incidence.emplace_back(face.negativeCell,face.nativeRow,face.area);
        if(face.positiveCell>=0)incidence.emplace_back(face.positiveCell,face.nativeRow,-face.area);
    }
    result.divergence_.resize(mesh.cellIndexNum(),count);result.divergence_.setFromTriplets(incidence.begin(),incidence.end());
    return result;
}

HydroConservativePressure HydroConservativePressure::buildDivergenceAdjoint(const QuadtreeMesh& mesh,const NativeLeafDomain& domain)
{
    require(domain.uniformLeaves(),"adjoint pressure requires a compatible AMR Hodge before mixed-leaf use");
    auto result=topology(mesh,domain);result.divergenceAdjoint_=true;
    result.gradientLoad_=-result.divergence_.transpose();
    for(const auto face:domain.faces(mesh))if(face.outwardSign!=0) {
        result.boundaryCoefficient_[face.nativeRow]=face.outwardSign*face.area;
        result.boundaryPoints_[face.nativeRow]=face.center;
    }
    return result;
}

HydroConservativePressure HydroConservativePressure::build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain)
{
    auto result=topology(mesh,domain);const int count=domain.nativeFaceCount();
    std::vector<Reconstruction> reconstruction;
    for(int c=0;c<mesh.cellIndexNum();++c)reconstruction.push_back(reconstruct(mesh,c));
    std::vector<Trace> rows(count);
    const auto add=[&](int row,int owner,Vector2D point,double measure) {
        const auto& fit=reconstruction[owner];rows[row][owner]+=measure;
        const Vector2D offset=point-fit.center;
        for(int axis=0;axis<2;++axis)for(const auto& term:fit.gradient[axis])
            rows[row][term.first]+=measure*(axis?offset.y:offset.x)*term.second;
    };
    for(const auto& segment:result.duals_.segments()) {
        const int row=std::max(segment.negativeNative,segment.positiveNative);
        const int component=row<mesh.vEdgeIndexNum()?0:1;
        if(segment.normalAxis!=component)continue;
        const auto point=weightedPoint(segment.begin,segment.end,result.axisymmetric_);
        if(segment.negativeNative>=0 && segment.positiveNative>=0) {
            require(segment.ownerCell>=0,"normal interior pressure segment lacks a real owner");
            add(segment.negativeNative,segment.ownerCell,point,segment.measure);
            add(segment.positiveNative,segment.ownerCell,point,-segment.measure);
        } else {
            require(segment.primalNative==row,"pressure exterior is not its native boundary face");
            result.boundaryCoefficient_[row]+=(segment.negativeNative>=0?1:-1)*segment.measure;
            result.boundaryPoints_[row]=point;
        }
    }
    for(const auto face:domain.faces(mesh)) {
        result.anchors_[face.nativeRow]=std::max(face.negativeCell,face.positiveCell);
        if(!result.axisymmetric_ || face.horizontal)continue;
        const auto& volume=result.duals_.volumes()[face.nativeRow];
        // On each owner half-piece, p/r dV = 2*pi*p dr dz. Midpoint
        // quadrature is exact for the same affine pressure trace used above.
        for(int owner:{face.negativeCell,face.positiveCell})if(owner>=0) {
            const auto& cell=cellAt(mesh,owner);
            const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
            const Vector2D a{std::max(lo.x,volume.lower.x),std::max(lo.y,volume.lower.y)},
                b{std::min(hi.x,volume.upper.x),std::min(hi.y,volume.upper.y)};
            require(b.x>a.x && b.y>a.y,"empty radial pressure half-piece");
            add(face.nativeRow,owner,{(a.x+b.x)/2,(a.y+b.y)/2},-CartesianMeasure::twoPi*(b.x-a.x)*(b.y-a.y));
        }
    }
    std::vector<Eigen::Triplet<double>> entries;
    for(int row=0;row<count;++row) {
        auto& coefficients=rows[row];require(!coefficients.empty(),"empty pressure load row");
        // Restore only roundoff in the constant identity C*1+b=0. The
        // largest coefficient absorbs it, avoiding a tiny new distant link.
        auto largest=coefficients.begin();long double sum=result.boundaryCoefficient_[row],magnitude=std::abs(sum);
        for(auto it=coefficients.begin();it!=coefficients.end();++it) {
            require(std::isfinite(it->second),"pressure coefficient overflow");
            sum+=it->second;magnitude+=std::abs(it->second);
            if(std::abs(it->second)>std::abs(largest->second))largest=it;
        }
        require(std::abs(sum)<=64*std::numeric_limits<double>::epsilon()*magnitude,"pressure constant mode does not close");
        largest->second-=static_cast<double>(sum);
        for(const auto& entry:coefficients)if(entry.second!=0)entries.emplace_back(row,entry.first,entry.second);
    }
    result.gradientLoad_.resize(count,mesh.cellIndexNum());result.gradientLoad_.setFromTriplets(entries.begin(),entries.end());
    return result;
}

void HydroConservativePressure::validate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) const
{
    domain.validateMeshIdentity(mesh);
    require(identity_ && identity_==domain.meshIdentity() && axisymmetric_==domain.axisymmetric(),"pressure operator/domain source differs");
}
std::vector<double> HydroConservativePressure::apply(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const std::vector<double>& scalar,const std::vector<double>& boundaryScalar) const
{
    validate(mesh,domain);
    require(scalar.size()==static_cast<std::size_t>(gradientLoad_.cols()) && boundaryScalar.size()==anchors_.size(),"pressure input size mismatch");
    for(double p:scalar)require(std::isfinite(p),"nonfinite cell pressure");
    for(double p:boundaryScalar)require(std::isfinite(p),"nonfinite boundary pressure");
    std::vector<double> result(anchors_.size());
    for(int row=0;row<gradientLoad_.rows();++row) {
        const double anchor=scalar[anchors_[row]];
        long double load=static_cast<long double>(boundaryCoefficient_[row])*(static_cast<long double>(boundaryScalar[row])-anchor);
        for(Matrix::InnerIterator it(gradientLoad_,row);it;++it)
            load+=static_cast<long double>(it.value())*(static_cast<long double>(scalar[it.col()])-anchor);
        result[row]=static_cast<double>(load);require(std::isfinite(load) && std::isfinite(result[row]),"pressure load overflow");
    }
    return result;
}
}
