#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineNative.h"
#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineAssembly.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentMassCurrent.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentComplement.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentBoundary.h"
#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineBasis.h"
#include "Core/Geometry/ParametricHermiteCurve.h"
#include <algorithm>
#include <limits>
#include <exception>
#include <stdexcept>

namespace Quantum::Hydrodynamics {
HydroCompatibleMomentSpline::NativeLayout::Matrix HydroCompatibleMomentSpline::NativeLayout::localCoordinateMap(
    const HydroCompatibleMomentComplement* complement,std::vector<int>& offsets) const{
    constexpr int terms=81;const int nq=int(polynomial.cols()),ne=complement?int(complement->basisMap().cols()):0;
    offsets.assign(cells.size()+1,0);
    for(int slot=0;slot<int(cells.size());slot++)offsets[slot+1]=offsets[slot]+terms+(complement?int(complement->cells().at(cells[slot].index).basis.size()):0);
    Matrix transform;
    if(!complement)transform=polynomial;
    else{
        std::vector<Eigen::Triplet<double>> entries;
        for(int col=0;col<polynomial.outerSize();col++)for(Matrix::InnerIterator it(polynomial,col);it;++it)
            entries.emplace_back(offsets[it.row()/terms]+it.row()%terms,col,it.value());
        const Eigen::SparseMatrix<double,Eigen::RowMajor> map=complement->basisMap();
        for(int slot=0;slot<int(cells.size());slot++){
            const auto& rows=complement->cells().at(cells[slot].index).basis;
            for(int k=0;k<int(rows.size());k++)for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(map,rows[k]);it;++it)
                entries.emplace_back(offsets[slot]+terms+k,nq+it.col(),it.value());
        }
        transform.resize(offsets.back(),nq+ne);transform.setFromTriplets(entries.begin(),entries.end());
    }
    return transform;
}
HydroCompatibleMomentSpline::NativeVolumeTransport HydroCompatibleMomentSpline::NativeLayout::volumeTransport(
    const HydroCompatibleMomentMassCurrent& current,const HydroCompatibleMomentComplement* complement) const{
    constexpr int terms=81;NativeVolumeTransport result;
    result.solenoidalCoordinates=int(polynomial.cols());result.complementCoordinates=complement?int(complement->basisMap().cols()):0;
    const auto assembly=coordinateAssembly(complement);
    std::lock_guard<std::mutex> lock(mutex);
    const auto& offsets=assembly->offsets;
    std::vector<Eigen::MatrixXd> blocks(cells.size());std::vector<double> geometricScales(cells.size());
    std::vector<std::exception_ptr> errors(cells.size());
    // Each cell owns its quadrature reduction. The final sparse assembly is
    // deterministic and uses the same q / physical-complement coordinate map.
#pragma omp parallel for schedule(dynamic)
    for(int slot=0;slot<int(cells.size());slot++){try{
        const int count=offsets[slot+1]-offsets[slot];const auto& cell=current.cells().at(cells[slot].index);
        if(complement&&cell.partition!=complement->cells().at(cell.index).partition)
            throw std::invalid_argument("volume current quadrature did not split the physical complement support");
        auto& block=blocks[slot];block=Eigen::MatrixXd::Zero(count,count);
        Eigen::Matrix<double,6,Eigen::Dynamic> jets(6,count);
        // Batch the SAME weak-form quadrature into dense matrix products.
        // Storage is bounded independently of the number of cut panels.
        constexpr int batch=12;
        Eigen::Matrix<double,2*batch,Eigen::Dynamic> values(2*batch,count),derivatives(2*batch,count);
        int used=0;
        const auto addCurrent=[&](){
            if(used)block.noalias()-=derivatives.topRows(2*used).transpose()*values.topRows(2*used);
            used=0;
        };
        for(const auto& p:cell.current){
            if(p.weight.x==0&&p.weight.y==0)continue;
            jets.leftCols(terms)=cellJets(slot,p.point);if(complement)jets.rightCols(count-terms)=complement->cellJets(cell.index,p.point);
            // Direct conservative weak form, including every q/eta cross
            // term. No assembled skew projection or metric correction.
            values.middleRows<2>(2*used)=jets.topRows<2>();
            derivatives.middleRows<2>(2*used)=p.weight.x*jets.middleRows<2>(2)+p.weight.y*jets.bottomRows<2>();
            if(++used==batch)addCurrent();
        }
        addCurrent();
        // phi_0=1, while every other scalar basis function has zero mean
        // under d^2 dx dy. Hence the coefficient of the constant in f is
        // <d^2 f>/<d^2>. This determines the entire observed constant-mode
        // defect; it is not a fit to a measured kinetic residual.
        auto& geometry=assembly->geometricGram[slot];
        const double residual=cell.afterMoments[0]-cell.beforeMoments[0]+cell.boundaryMoments[0];
        geometricScales[slot]=residual;
        if(residual!=0&&geometry.size()==0){
        Eigen::MatrixXd freshGeometry=Eigen::MatrixXd::Zero(count,count);
        const double beta=cell.rational?cell.size.x/cell.lower.x:0;
        const double norm=cell.size.x*cell.size.y*(1+beta+beta*beta/3);
        const MomentSplineBasis::GaussTwelve rule;
        const auto& xs=cell.partition[0];const auto& ys=cell.partition[1];
        const auto addGeometry=[&](){
            if(used)freshGeometry.noalias()+=values.topRows(2*used).transpose()*values.topRows(2*used);
            used=0;
        };
        for(std::size_t x=1;x<xs.size();x++)for(std::size_t y=1;y<ys.size();y++)
        for(int i=0;i<12;i++)for(int j=0;j<12;j++){
            const double hx=xs[x]-xs[x-1],hy=ys[y]-ys[y-1];
            const Vector2D p{xs[x-1]+.5*hx*(1+rule.x[i]),ys[y-1]+.5*hy*(1+rule.x[j])};
            const double d=1+beta*(p.x-cell.lower.x)/cell.size.x;
            const double weight=.25*hx*hy*rule.w[i]*rule.w[j]*d*d/norm;
            jets.leftCols(terms)=cellJets(slot,p);if(complement)jets.rightCols(count-terms)=complement->cellJets(cell.index,p);
            values.middleRows<2>(2*used)=std::sqrt(weight)*jets.topRows<2>();
            if(++used==batch)addGeometry();
        }
        addGeometry();
        if(!freshGeometry.allFinite())throw std::invalid_argument("nonfinite constant-mode observation metric");
        geometry=std::move(freshGeometry);
        }
        if(!block.allFinite()||!geometry.allFinite())throw std::invalid_argument("nonfinite complete volume momentum integral");
    }catch(...){errors[slot]=std::current_exception();}}
    for(const auto& error:errors)if(error)std::rethrow_exception(error);
    result.volume=assembly->assemble([&](int slot)->const Eigen::MatrixXd&{return blocks[slot];});
    result.geometricMetric=assembly->assemble([&](int slot)->const Eigen::MatrixXd&{return assembly->geometricGram[slot];},&geometricScales);return result;
}
HydroCompatibleMomentSpline::NativeMomentumTransport HydroCompatibleMomentSpline::NativeLayout::momentumTransport(
    const HydroCompatibleMomentMassCurrent& current,const HydroCompatibleMomentBoundary& boundary,
    const HydroCompatibleMomentComplement* complement,const Eigen::VectorXd& coordinates) const{
    constexpr int terms=81,degree=10;using namespace MomentSplineBasis::Bernstein;
    const auto require=[](bool ok,const char* message){if(!ok)throw std::invalid_argument(message);};
    NativeMomentumTransport result;static_cast<NativeVolumeTransport&>(result)=volumeTransport(current,complement);
    const int nq=result.solenoidalCoordinates,n=result.volume.rows();
    result.boundaryLoad=Eigen::VectorXd::Zero(n);
    const auto& mass=current.massFlux();const auto& phases=current.facePhase();
    require(coordinates.size()==nq&&mass.products.size()==faces.size(),"complete face momentum source differs");
    const auto assembly=coordinateAssembly(complement);const auto& offsets=assembly->offsets;const auto& transform=assembly->rowTransform;
    std::vector<Eigen::Triplet<double>> faceEntries,metricEntries;
    const auto& rule=Geometry::Parametric::gaussEight();
    const auto gram=[&](const Eigen::VectorXd& weights){Eigen::Matrix<double,11,11> h;
        require(weights.size()==21&&weights.allFinite(),"invalid shared face work functional");
        for(int i=0;i<=degree;i++)for(int j=0;j<=degree;j++)h(i,j)=double(choose(degree,i)*choose(degree,j)/choose(2*degree,i+j))*weights[i+j];return h;};
    for(int row=0;row<int(faces.size());row++){
        const auto& f=faces[row];if(!f.boundary)continue;
        const auto& c=cells[f.cell];const auto& port=boundary.faces().at(row);const auto& phase=phases.at(row);
        if(HydroCompatibleMomentBoundary::impermeable(port.kind)){
            require(mass.products[row].isZero(0)&&phase.total==0,"closed complete-velocity boundary carries mass");continue;}
        const double normal=f.horizontal?f.chart.begin.y:f.chart.begin.x;
        const double sign=normal==(f.horizontal?c.lower.y:c.lower.x)?-1.:1.;
        const double totalIn=sign>0?phase.totalNegative:phase.totalPositive,liquidIn=sign>0?phase.liquidNegative:phase.liquidPositive;
        const double volume=c.size.x*c.size.y*(rz?2*std::acos(-1.)*(c.lower.x+.5*c.size.x):1.);
        if(current.tracksInterface()){
            require(port.incomingPhase==1||port.incomingPhase==-1,"complete momentum port has no declared incoming phase");
            require(std::abs(liquidIn-(port.incomingPhase<0?totalIn:0.))<=1e-12*volume,"complete momentum incoming phase differs from actual geometry");
        }
        require(totalIn==0||bool(port.incomingVelocity),"complete momentum backflow needs prescribed incoming velocity");
        if(totalIn>0&&port.kind==HydroCompatibleMomentBoundary::Kind::PressureOutlet)++result.pressureBackflowFaces;
        if(phase.totalPositive>0&&phase.totalNegative>0)++result.bidirectionalFaces;
        const Eigen::VectorXd out=sign>0?mass.exteriorPositiveProducts[row]:Eigen::VectorXd(mass.exteriorPositiveProducts[row]-mass.products[row]);
        const Eigen::VectorXd in=out-sign*mass.products[row];
        const double directionScale=std::max(out.cwiseAbs().maxCoeff(),in.cwiseAbs().maxCoeff());
        require(out.minCoeff()>=-512*std::numeric_limits<double>::epsilon()*directionScale
            &&in.minCoeff()>=-512*std::numeric_limits<double>::epsilon()*directionScale,"shared directional face mass is negative");
        const auto ho=gram(out),hi=gram(in);
        const int count=offsets[f.cell+1]-offsets[f.cell];
        std::vector<int> active;for(int k=offsets[f.cell];k<offsets[f.cell+1];k++)
            for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(transform,k);it;++it)active.push_back(int(it.col()));
        std::sort(active.begin(),active.end());active.erase(std::unique(active.begin(),active.end()),active.end());
        Eigen::MatrixXd map=Eigen::MatrixXd::Zero(count,active.size());
        for(int k=0;k<count;k++)for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(transform,offsets[f.cell]+k);it;++it)
            map(k,std::lower_bound(active.begin(),active.end(),it.col())-active.begin())=it.value();
        const auto native=faceTrace(row);NativeFaceTrace local;local.x.resize(11,count);local.y.resize(11,count);
        local.x.leftCols(terms)=native.x;local.y.leftCols(terms)=native.y;
        if(complement){const auto extra=complement->faceTrace(c.index,f.chart,f.horizontal);
            local.x.rightCols(count-terms)=extra.x;local.y.rightCols(count-terms)=extra.y;}
        const Eigen::MatrixXd x=local.x*map,y=local.y*map;
        // Verify algebraic traces against the actual physical basis, not
        // against another polynomial interpolation of the same samples.
        const double beta=rz&&f.horizontal&&f.chart.begin.x>0?(f.chart.end.x-f.chart.begin.x)/f.chart.begin.x:0;
        for(int g=0;g<8;g++){const double t=.5*(1+rule.x[g]),d=1+beta*t;const auto p=f.chart.begin+(f.chart.end-f.chart.begin)*t;
            Eigen::Matrix<double,1,11> b;for(int k=0;k<=degree;k++)b[k]=double(choose(degree,k))*std::pow(t,k)*std::pow(1-t,degree-k)/d;
            Eigen::MatrixXd actual(2,count);actual.leftCols(terms)=cellJets(f.cell,p).topRows<2>();
            if(complement)actual.rightCols(count-terms)=complement->cellJets(c.index,p).topRows(2);
            const Eigen::MatrixXd expected=actual*map;Eigen::MatrixXd represented(2,active.size());represented.row(0)=b*x;represented.row(1)=b*y;
            result.traceDefect=std::max(result.traceDefect,(represented-expected).norm()/(1+expected.norm()));
        }
        Eigen::VectorXd d(11);for(int k=0;k<=degree;k++)d[k]=1+beta*k/degree;
        const Eigen::MatrixXd a=x.transpose()*ho*x+y.transpose()*ho*y;
        const Eigen::MatrixXd metric=a+x.transpose()*hi*x+y.transpose()*hi*y;
        Eigen::Vector2d exterior=Eigen::Vector2d::Zero();if(port.incomingVelocity)exterior<<port.incomingVelocity->x,port.incomingVelocity->y;
        const Eigen::VectorXd supplied=exterior[0]*x.transpose()*hi*d+exterior[1]*y.transpose()*hi*d;
        Eigen::VectorXd q=Eigen::VectorXd::Zero(active.size());for(int k=0;k<q.size();k++)if(active[k]<nq)q[k]=coordinates[active[k]];
        const Eigen::VectorXd ux=x*q,uy=y*q,dx=ux-exterior[0]*d,dy=uy-exterior[1]*d;
        result.boundaryKineticFlux+=.5*(ux.dot(ho*ux)+uy.dot(ho*uy)-exterior.squaredNorm()*d.dot(hi*d));
        result.boundaryInflowDissipation+=.5*(dx.dot(hi*dx)+dy.dot(hi*dy));
        result.boundaryVolume+=sign*phase.total;
        for(int i=0;i<int(active.size());i++){result.boundaryLoad[active[i]]+=supplied[i];
            for(int j=0;j<int(active.size());j++){if(a(i,j)!=0)faceEntries.emplace_back(active[i],active[j],a(i,j));
                if(metric(i,j)!=0)metricEntries.emplace_back(active[i],active[j],metric(i,j));}}
    }
    require(result.traceDefect<2e-10,"complete face momentum trace differs from physical velocity basis");
    Matrix face(n,n);face.setFromTriplets(faceEntries.begin(),faceEntries.end());result.central=result.volume+face;
    result.boundaryMetric.resize(n,n);result.boundaryMetric.setFromTriplets(metricEntries.begin(),metricEntries.end());
    return result;
}
}
