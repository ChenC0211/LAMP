#include "Modules/Hydrodynamics/CchEnergyViscosity.h"

#include "Core/Mesh/QuadtreeMesh.h"
#include <Eigen/QR>
#include <Eigen/SVD>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Hydrodynamics
{
namespace
{
using Matrix=CchStrainOperator::Matrix;
using Axis=CchDofLayout::NormalAxis;
void require(bool ok,const std::string& reason)
{
    if (!ok) throw std::invalid_argument("[CchEnergyViscosity] "+reason);
}
void validate(const QuadtreeMesh& mesh,const CchNativeValues& values,const CchPhysicalDivergence& metrics)
{
    values.validateFor(metrics.layout(),mesh);
    for (int u : metrics.axisVelocityDofs())
        require(values.normalVelocity()[u]==0,"nonzero radial axis velocity");
}
std::vector<double> multiply(const Matrix& matrix,const std::vector<double>& input)
{
    const Eigen::Map<const Eigen::VectorXd> u(input.data(),static_cast<Eigen::Index>(input.size()));
    const Eigen::VectorXd result=matrix*u;
    require(result.allFinite(),"nonfinite operator result");
    return {result.data(),result.data()+result.size()};
}
struct LocalSample {int dof;double sign;Vector2D point;};
struct Fit {std::vector<LocalSample> samples;Eigen::MatrixXd coefficients;};
Eigen::Matrix<double,1,6> basis(double x,double y)
{
    Eigen::Matrix<double,1,6> result;result<<1,x,y,x*x,x*y,y*y;return result;
}
Fit fit(const CchPhysicalDivergence& metrics,const std::vector<NativeLeafFace>& faces,
    const std::vector<int>& support,const std::vector<int>& nativeFaces,Axis axis,
    Vector2D center,Vector2D width,bool reflected,int cell)
{
    Fit result;
    for (int u : support) if (metrics.layout().normalAxis(u)==axis) {
        result.samples.push_back({u,1,faces[u].center});
        if (reflected && faces[u].center.x>0)
            result.samples.push_back({u,metrics.reflectAxisNormalVelocity(axis,1),{-faces[u].center.x,faces[u].center.y}});
    }
    Eigen::MatrixXd a(result.samples.size(),6);
    for (Eigen::Index i=0;i<a.rows();++i) {
        const auto p=result.samples[static_cast<std::size_t>(i)].point;
        a.row(i)=basis((p.x-center.x)/width.x,(p.y-center.y)/width.y);
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(a);
    qr.setThreshold(1e-12);
    require(qr.rank()==6,"rank-deficient local P2 support at cell "+std::to_string(cell));
    const Eigen::MatrixXd identity=Eigen::MatrixXd::Identity(a.rows(),a.rows());
    if (reflected) {
        // Enforce regular parity algebraically, not by cancellation of large
        // 1/r coefficients. Odd radial P2: r,r*z. Even axial P2: 1,z,r*r,z*z.
        const std::vector<int> columns=axis==Axis::First ? std::vector<int>{1,4} : std::vector<int>{0,2,3,5};
        Eigen::MatrixXd regular(a.rows(),columns.size());
        for (Eigen::Index j=0;j<regular.cols();++j) regular.col(j)=a.col(columns[static_cast<std::size_t>(j)]);
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> parity(regular);parity.setThreshold(1e-12);
        require(parity.rank()==regular.cols(),"rank-deficient axis-regular strain support");
        const Eigen::MatrixXd coefficients=parity.solve(identity);
        result.coefficients=Eigen::MatrixXd::Zero(6,a.rows());
        for (Eigen::Index j=0;j<regular.cols();++j) result.coefficients.row(columns[static_cast<std::size_t>(j)])=coefficients.row(j);
    }
    else {
        // Interpolate this cell's actual normal faces. On a regular cell this
        // preserves the central normal difference, rather than filtering it.
        std::vector<Eigen::Index> constraints;
        for (int u : nativeFaces) if (metrics.layout().normalAxis(u)==axis)
            for (Eigen::Index j=0;j<a.rows();++j)
                if (result.samples[static_cast<std::size_t>(j)].dof==u) {constraints.push_back(j);break;}
        Eigen::MatrixXd c(constraints.size(),6),selector=Eigen::MatrixXd::Zero(constraints.size(),a.rows());
        for (Eigen::Index i=0;i<c.rows();++i) {c.row(i)=a.row(constraints[static_cast<std::size_t>(i)]);selector(i,constraints[static_cast<std::size_t>(i)])=1;}
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(c,Eigen::ComputeThinU|Eigen::ComputeFullV);
        svd.setThreshold(1e-12);
        require(svd.rank()==c.rows(),"dependent native face constraints at cell "+std::to_string(cell));
        const Eigen::MatrixXd particular=svd.solve(selector);
        const Eigen::MatrixXd nullspace=svd.matrixV().rightCols(6-c.rows());
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> reduced(a*nullspace);reduced.setThreshold(1e-12);
        require(reduced.rank()==nullspace.cols(),"rank-deficient constrained strain support");
        result.coefficients=particular+nullspace*reduced.solve(identity-a*particular);
    }
    require(result.coefficients.allFinite(),"nonfinite local reconstruction");
    return result;
}
}
struct CchStrainOperator::Data
{
    explicit Data(CchPhysicalDivergence input) : metrics(std::move(input)) {}
    CchPhysicalDivergence metrics;
    std::vector<CchStrainPoint> points;
    Matrix matrix;
    int components=3;
};
CchStrainOperator::CchStrainOperator(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
const CchPhysicalDivergence& CchStrainOperator::metrics() const noexcept {return data_->metrics;}
const CchDofLayout& CchStrainOperator::layout() const noexcept {return metrics().layout();}
const std::vector<CchStrainPoint>& CchStrainOperator::points() const noexcept {return data_->points;}
int CchStrainOperator::components() const noexcept {return data_->components;}
const Matrix& CchStrainOperator::matrix() const noexcept {return data_->matrix;}
std::vector<CchStrainPoint> CchStrainOperator::volumeQuadrature(const QuadtreeMesh& mesh,
    const CchPhysicalDivergence& metrics)
{
    const auto& layout=metrics.layout();layout.validateMesh(mesh);
    std::vector<CchStrainPoint> result;result.reserve(static_cast<std::size_t>(layout.pressureCount())*4);
    const double g=1/std::sqrt(3.0);
    for (int p=0;p<layout.pressureCount();++p) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[layout.nativeCell(p)]);
        const auto c=cell.centerCoord(),w=cell.size();
        for (double x : {-g,g}) for (double y : {-g,g}) {
            const Vector2D point{c.x+x*w.x/2,c.y+y*w.y/2};
            const double weight=metrics.cellVolumes()[p]/4*(layout.domain().axisymmetric() ? point.x/c.x : 1);
            require(std::isfinite(weight) && weight>0,"invalid physical volume quadrature");
            result.push_back({p,point,weight});
        }
    }
    return result;
}
CchStrainOperator CchStrainOperator::build(const QuadtreeMesh& mesh,CchPhysicalDivergence metrics,
    const std::vector<CchStrainPoint>& points)
{
    metrics.layout().validateMesh(mesh);
    require(!points.empty(),"empty strain quadrature");
    auto data=std::make_shared<Data>(std::move(metrics));data->points=points;
    const auto& layout=data->metrics.layout();const bool rz=layout.domain().axisymmetric();
    data->components=rz ? 4 : 3;
    require(points.size()<=static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max()/data->components),"too many strain points");
    std::vector<Vector2D> centers(layout.pressureCount()),widths(layout.pressureCount());
    std::vector<std::vector<std::size_t>> pointRows(layout.pressureCount());
    for (int p=0;p<layout.pressureCount();++p) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[layout.nativeCell(p)]);
        centers[p]=cell.centerCoord();widths[p]=cell.size();
    }
    for (std::size_t q=0;q<points.size();++q) {
        const auto& sample=points[q];const int p=sample.cellDof;
        require(p>=0 && p<layout.pressureCount(),"invalid quadrature cell");
        require(std::isfinite(sample.point.x) && std::isfinite(sample.point.y)
            && std::isfinite(sample.physicalWeight) && sample.physicalWeight>=0,"invalid point/physical weight");
        const auto c=centers[p],w=widths[p];
        require(sample.point.x>=c.x-w.x/2 && sample.point.x<=c.x+w.x/2
            && sample.point.y>=c.y-w.y/2 && sample.point.y<=c.y+w.y/2,"quadrature point outside native cell");
        pointRows[p].push_back(q);
    }
    std::vector<NativeLeafFace> faces;faces.reserve(layout.velocityCount());
    std::vector<std::vector<int>> cellFaces(layout.pressureCount()),neighbors(layout.pressureCount());
    for (int u=0;u<layout.velocityCount();++u) {
        const auto f=layout.domain().face(mesh,layout.nativeFace(u));faces.push_back(f);
        const int a=f.negativeCell<0 ? -1 : layout.pressureDof(f.negativeCell);
        const int b=f.positiveCell<0 ? -1 : layout.pressureDof(f.positiveCell);
        if (a>=0) cellFaces[a].push_back(u);
        if (b>=0) cellFaces[b].push_back(u);
        if (a>=0 && b>=0) {neighbors[a].push_back(b);neighbors[b].push_back(a);}
    }
    std::vector<Eigen::Triplet<double,Eigen::Index>> entries;
    for (int p=0;p<layout.pressureCount();++p) if (!pointRows[p].empty()) {
        std::vector<int> patch{p},frontier{p};
        for (int ring=0;ring<2;++ring) {
            std::vector<int> next;
            for (int c : frontier) for (int n : neighbors[c])
                if (std::find(patch.begin(),patch.end(),n)==patch.end()) {patch.push_back(n);next.push_back(n);}
            frontier=std::move(next);
        }
        std::vector<int> support;
        for (int c : patch) support.insert(support.end(),cellFaces[c].begin(),cellFaces[c].end());
        std::sort(support.begin(),support.end());support.erase(std::unique(support.begin(),support.end()),support.end());
        const bool reflected=rz && std::any_of(support.begin(),support.end(),[&](int u) {return faces[u].area==0;});
        const auto c=reflected ? Vector2D{0,centers[p].y} : centers[p];const auto w=widths[p];
        const std::array<Fit,2> fits{fit(data->metrics,faces,support,cellFaces[p],Axis::First,c,w,reflected,p),
            fit(data->metrics,faces,support,cellFaces[p],Axis::Second,c,w,reflected,p)};
        for (std::size_t q : pointRows[p]) {
            const auto point=points[q].point;const double x=(point.x-c.x)/w.x,y=(point.y-c.y)/w.y;
            Eigen::Matrix<double,1,6> dx,dy;
            dx<<0,1/w.x,0,2*x/w.x,y/w.x,0;
            dy<<0,0,1/w.y,0,x/w.y,2*y/w.y;
            const auto value=basis(x,y);
            for (int axis=0;axis<2;++axis) for (int component=0;component<data->components;++component) {
                Eigen::Matrix<double,1,6> functional=Eigen::Matrix<double,1,6>::Zero();
                if (component==0 && axis==0) functional=dx;
                if (component==1 && axis==1) functional=dy;
                if (rz && component==2 && axis==0) {
                    require(point.x>0 || reflected,"hoop axis limit requires retained axis");
                    functional=point.x==0 ? dx : value/point.x;
                }
                if (component==data->components-1) functional=(axis==0 ? dy : dx)/std::sqrt(2.0);
                const Eigen::RowVectorXd coefficients=functional*fits[axis].coefficients;
                require(coefficients.allFinite(),"nonfinite strain coefficients");
                for (Eigen::Index j=0;j<coefficients.size();++j) if (coefficients[j]!=0)
                    entries.emplace_back(static_cast<Eigen::Index>(q)*data->components+component,
                        fits[axis].samples[static_cast<std::size_t>(j)].dof,
                        coefficients[j]*fits[axis].samples[static_cast<std::size_t>(j)].sign);
            }
        }
    }
    data->matrix.resize(static_cast<Eigen::Index>(points.size())*data->components,layout.velocityCount());
    data->matrix.setFromTriplets(entries.begin(),entries.end());data->matrix.prune(0.0);data->matrix.makeCompressed();
    for (Eigen::Index i=0;i<data->matrix.nonZeros();++i) require(std::isfinite(data->matrix.valuePtr()[i]),"nonfinite assembled strain");
    return CchStrainOperator(std::move(data));
}
std::vector<double> CchStrainOperator::apply(const QuadtreeMesh& mesh,const CchNativeValues& values) const
{
    validate(mesh,values,metrics());return multiply(matrix(),values.normalVelocity());
}
struct CchEnergyViscosity::Data
{
    explicit Data(CchStrainOperator input) : strain(std::move(input)) {}
    CchStrainOperator strain;
    std::vector<double> viscosity,weights;
    Matrix matrix;
};
CchEnergyViscosity::CchEnergyViscosity(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
const CchStrainOperator& CchEnergyViscosity::strain() const noexcept {return data_->strain;}
const std::vector<double>& CchEnergyViscosity::viscosity() const noexcept {return data_->viscosity;}
const std::vector<double>& CchEnergyViscosity::rowWeights() const noexcept {return data_->weights;}
const Matrix& CchEnergyViscosity::matrix() const noexcept {return data_->matrix;}
CchEnergyViscosity CchEnergyViscosity::build(CchStrainOperator strain,const std::vector<double>& viscosity)
{
    require(viscosity.size()==strain.points().size(),"viscosity dimensions do not match quadrature");
    auto data=std::make_shared<Data>(std::move(strain));data->viscosity=viscosity;
    for (std::size_t q=0;q<viscosity.size();++q) {
        require(std::isfinite(viscosity[q]) && viscosity[q]>=0,"negative/nonfinite viscosity");
        const long double product=2*static_cast<long double>(viscosity[q])*data->strain.points()[q].physicalWeight;
        const double weight=static_cast<double>(product);
        require(std::isfinite(weight),"nonfinite viscosity quadrature weight");
        for (int k=0;k<data->strain.components();++k) data->weights.push_back(weight);
    }
    const auto& e=data->strain.matrix();data->matrix.resize(e.cols(),e.cols());
    if (std::any_of(data->weights.begin(),data->weights.end(),[](double x) {return x!=0;})) {
        Matrix weighted=e;
        for (Eigen::Index row=0;row<weighted.outerSize();++row)
            for (Matrix::InnerIterator entry(weighted,row);entry;++entry) entry.valueRef()*=data->weights[static_cast<std::size_t>(row)];
        data->matrix=e.transpose()*weighted;
        data->matrix.prune(0.0);
    }
    data->matrix.makeCompressed();
    for (Eigen::Index i=0;i<data->matrix.nonZeros();++i) require(std::isfinite(data->matrix.valuePtr()[i]),"nonfinite assembled viscosity");
    return CchEnergyViscosity(std::move(data));
}
std::vector<double> CchEnergyViscosity::apply(const QuadtreeMesh& mesh,const CchNativeValues& values) const
{
    validate(mesh,values,strain().metrics());return multiply(matrix(),values.normalVelocity());
}
double CchEnergyViscosity::dissipation(const QuadtreeMesh& mesh,const CchNativeValues& values) const
{
    if (std::none_of(rowWeights().begin(),rowWeights().end(),[](double x) {return x!=0;})) {
        validate(mesh,values,strain().metrics());return 0;
    }
    const auto e=strain().apply(mesh,values);long double total=0;
    for (std::size_t i=0;i<e.size();++i) total+=static_cast<long double>(rowWeights()[i])*e[i]*e[i];
    const double result=static_cast<double>(total);require(std::isfinite(result),"nonfinite viscous dissipation");return result;
}
}

