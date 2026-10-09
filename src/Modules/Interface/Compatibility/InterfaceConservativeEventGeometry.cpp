#include "Core/Geometry/ActiveHermitePatch.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Modules/Interface/Compatibility/InterfaceState.h"
#include "Modules/Interface/InterfaceCellMeasure.h"
#include "Modules/Interface/Compatibility/InterfaceVolumeCorrection.h"
#include "Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <memory>
#include <stdexcept>
#include <limits>
#include <map>
#include <Eigen/SuperLUSupport>

// Local constrained Newton synchronization for an already declared interior
// event. Independent FV inventories never enter this file as writable data.
// Inequalities protect pure Bernstein cells and admitted graph/axis traces;
// actual native discovery remains the final topology certificate.
namespace {
using Vector=Eigen::VectorXd;using Matrix=Eigen::MatrixXd;
Vector boundedMinimum(const Matrix& a,const Vector& b,const Matrix& g,const Vector& bound)
{
    const int n=static_cast<int>(a.cols()),m=static_cast<int>(a.rows()),p=static_cast<int>(g.rows());
    if(n<=0||n>2048||m<=0||m>n||p<=0||g.cols()!=n||b.size()!=m||bound.size()!=p
        ||!a.allFinite()||!b.allFinite()||!g.allFinite()||!bound.allFinite())
        throw std::runtime_error("invalid or nonlocal bounded volume system");
    Vector x=Vector::Zero(n),y=Vector::Zero(m),z=Vector::Ones(p),s=bound.cwiseMax(1.0);
    const auto positiveStep=[](const Vector& value,const Vector& delta) {
        double alpha=1;for(int k=0;k<value.size();++k)if(delta[k]<0)alpha=std::min(alpha,-value[k]/delta[k]);return alpha;
    };
    for(int iteration=0;iteration<60;++iteration) {
        const Vector rd=x+a.transpose()*y+g.transpose()*z,rp=a*x-b,ri=g*x+s-bound;
        const double mu=s.dot(z)/p;
        if(std::max({rd.lpNorm<Eigen::Infinity>(),rp.lpNorm<Eigen::Infinity>(),ri.lpNorm<Eigen::Infinity>(),mu})<1e-10) {
            if((g*x-bound).maxCoeff()>1e-9 || (a*x-b).lpNorm<Eigen::Infinity>()>1e-9)throw std::runtime_error("bounded linear residual");
            return x;
        }
        const Vector ratio=z.cwiseQuotient(s);
        const Matrix hessian=Matrix::Identity(n,n)+g.transpose()*ratio.asDiagonal()*g;
        if(!hessian.allFinite())throw std::runtime_error("nonfinite bounded volume Hessian");
        Matrix kkt=Matrix::Zero(n+m,n+m);kkt.topLeftCorner(n,n)=hessian;kkt.topRightCorner(n,m)=a.transpose();kkt.bottomLeftCorner(m,n)=a;
        Eigen::SparseMatrix<double> sparse=kkt.sparseView();Eigen::SuperLU<Eigen::SparseMatrix<double>> solver;solver.compute(sparse);
        if(solver.info()!=Eigen::Success)throw std::runtime_error("bounded KKT factorization failed");
        const auto direction=[&](const Vector& cent,Vector& dx,Vector& dy,Vector& ds,Vector& dz) {
            Vector rhs(n+m);rhs.head(n)=-rd+g.transpose()*((cent-z.cwiseProduct(ri)).cwiseQuotient(s));rhs.tail(m)=-rp;
            Vector solution=solver.solve(rhs);
            if(solver.info()!=Eigen::Success||!solution.allFinite())throw std::runtime_error("bounded KKT solve failed");
            const double scale=std::max(1.0,rhs.lpNorm<Eigen::Infinity>()+sparse.norm()*solution.norm());
            if(!std::isfinite(scale))throw std::runtime_error("unrepresentable bounded KKT scale");
            const double backward=(sparse*solution-rhs).lpNorm<Eigen::Infinity>()/scale;
            if(!std::isfinite(backward)||backward>1e-10)throw std::runtime_error("bounded KKT backward residual failed");
            dx=solution.head(n);dy=solution.tail(m);ds=-ri-g*dx;dz=(-cent-z.cwiseProduct(ds)).cwiseQuotient(s);
        };
        Vector dx,dy,ds,dz;direction(s.cwiseProduct(z),dx,dy,ds,dz);
        const double primal=positiveStep(s,ds),dual=positiveStep(z,dz);
        const double affine=(s+primal*ds).dot(z+dual*dz)/p;
        const double sigma=std::pow(std::max(0.0,affine/mu),3);
        Vector cent=s.cwiseProduct(z)+ds.cwiseProduct(dz)-Vector::Constant(p,sigma*mu);
        direction(cent,dx,dy,ds,dz);
        const double ap=std::min(1.,.995*positiveStep(s,ds)),ad=std::min(1.,.995*positiveStep(z,dz));
        x+=ap*dx;s+=ap*ds;y+=ad*dy;z+=ad*dz;
        if(!x.allFinite()||s.minCoeff()<=0||z.minCoeff()<=0)throw std::runtime_error("bounded KKT left positive cone");
    }
    throw std::runtime_error("bounded KKT iteration budget exhausted");
}
InterfaceVolumeCorrection::Result phaseBoundedStep(const QuadtreeMesh& mesh,
    const InterfaceCellMeasure::Measure& measure,const std::vector<long double>& defects,
    const std::vector<std::uint8_t>& fixed,const std::vector<InterfaceVolumeCorrection::GradientDof>& gradients,
    const Quantum::Geometry::ActiveHermitePatch& patch,const Quantum::Geometry::ActiveHermitePatch& reference,
    const std::vector<int>& phase,double h)
{
    std::map<std::pair<int,int>,int> columns;
    for(int v=0;v<static_cast<int>(fixed.size());++v)if(!fixed[v])columns[{v,0}]=static_cast<int>(columns.size());
    for(const auto& q:gradients)columns[{q.vertexIndex,static_cast<int>(q.component)}]=static_cast<int>(columns.size());
    const int n=static_cast<int>(columns.size());std::map<int,Vector> rows;
    for(const auto& arc:measure.arcs) {
        auto it=rows.emplace(arc.cellIndex,Vector::Zero(n)).first;
        for(const auto& q:arc.independentSensitivity)for(int k=0;k<3;++k) {
            const auto col=columns.find({q.vertexIndex,k});
            if(col!=columns.end())it->second[col->second]+=static_cast<double>(q.coefficients[k])*(k==0?h:1);
        }
    }
    Matrix a(rows.size(),n);Vector b(rows.size());int index=0;
    for(const auto& row:rows){const double norm=row.second.norm();if(!(norm>0))throw std::runtime_error("empty bounded volume row");a.row(index)=row.second/norm;b[index++]=static_cast<double>(defects[row.first])/norm;}
    std::vector<Vector> inequalities;std::vector<double> bounds;
    for(int c=0;c<mesh.cellIndexNum();++c)if(phase[c]!=0) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto size=cell.size();
        for(int corner=0;corner<4;++corner) {
            const int v=mesh.verts(cell.cellToVert[corner]).index();if(fixed[v])continue;
            const auto& node=patch.nodes[v];const int sx=(corner==0||corner==3)?1:-1,sy=corner<2?1:-1,sign=phase[c];
            for(int control=0;control<4;++control) {
                const double cx=(control&1)?sx*size.x/3:0,cy=(control&2)?sy*size.y/3:0;
                const double value=sign*(node.value+cx*node.dx+cy*node.dy+cx*cy*node.dxy)/h;
                Vector row=Vector::Zero(n);for(int k=0;k<3;++k) {
                    const auto col=columns.find({v,k});if(col!=columns.end())row[col->second]=-sign*(k==0?1:(k==1?cx:cy)/h);
                }
                if(row.squaredNorm()!=0){inequalities.push_back(std::move(row));bounds.push_back(value);}
            }
        }
    }
    // Preserve the source's nodal phase during this zero-time correction.
    // Physical motion belongs to the preceding advection, not this solve.
    for(int v=0;v<static_cast<int>(fixed.size());++v)if(!fixed[v]&&patch.nodes[v].value!=0) {
        const double sign=patch.nodes[v].value>0?1:-1;Vector row=Vector::Zero(n);
        row[columns.at({v,0})]=-sign;inequalities.push_back(std::move(row));bounds.push_back(std::abs(patch.nodes[v].value)/h);
    }
    
    for(const auto& item:rows) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[item.first]);const auto size=cell.size();
        std::array<double,16> controls;std::array<Vector,16> variations;
        for(int j=0;j<4;++j)for(int i=0;i<4;++i) {
            const bool right=i>=2,top=j>=2;const int corner=top?(right?2:3):(right?1:0);
            const int v=mesh.verts(cell.cellToVert[corner]).index();const auto& node=patch.nodes[v];
            const double cx=(i==1||i==2)?(right?-1:1)*size.x/3:0,cy=(j==1||j==2)?(top?-1:1)*size.y/3:0;
            const int k=4*j+i;controls[k]=(node.value+cx*node.dx+cy*node.dy+cx*cy*node.dxy)/h;
            variations[k]=Vector::Zero(n);
            for(int derivative=0;derivative<3;++derivative) {const auto col=columns.find({v,derivative});if(col!=columns.end())variations[k][col->second]=derivative==0?1:(derivative==1?cx:cy)/h;}
        }
        for(int axis=0;axis<2;++axis) {
            int sign=0;bool valid=true;
            for(int j=0;j<(axis?3:4);++j)for(int i=0;i<(axis?4:3);++i) {
                const int lo=4*j+i,hi=lo+(axis?4:1);const double d=controls[hi]-controls[lo];
                if(d!=0){const int next=d>0?1:-1;if(sign&&next!=sign)valid=false;sign=next;}
            }
            if(!valid||!sign)continue;
            for(int j=0;j<(axis?3:4);++j)for(int i=0;i<(axis?4:3);++i) {
                const int lo=4*j+i,hi=lo+(axis?4:1);Vector row=-sign*(variations[hi]-variations[lo]);
                if(row.squaredNorm()!=0){inequalities.push_back(std::move(row));bounds.push_back(sign*(controls[hi]-controls[lo]));}
            }
        }
    }
    
    for(const auto& item:rows) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[item.first]);
        if(mesh.verts(cell.cellToVert[0]).coord().x!=0)continue;
        const int bottom=mesh.verts(cell.cellToVert[0]).index(),top=mesh.verts(cell.cellToVert[3]).index();
        const auto& a0=reference.nodes[bottom];const auto& b0=reference.nodes[top];
        if(a0.value*b0.value>=0)continue;
        const double height=cell.size().y;
        const auto controls=[&](const auto& a,const auto& b){return std::array<double,4>{a.value/h,(a.value+height*a.dy/3)/h,(b.value-height*b.dy/3)/h,b.value/h};};
        const auto initial=controls(a0,b0),current=controls(patch.nodes[bottom],patch.nodes[top]);
        std::array<Vector,4> v;for(auto& q:v)q=Vector::Zero(n);
        for(int end=0;end<2;++end) {
            const int id=end?top:bottom;const auto value=columns.find({id,0}),dy=columns.find({id,2});
            const int outer=end?3:0,inner=end?2:1;
            if(value!=columns.end()){v[outer][value->second]=1;v[inner][value->second]=1;}
            if(dy!=columns.end())v[inner][dy->second]=(end?-1:1)*height/(3*h);
        }
        // A regular cap need not be monotone over its entire axis cell.
        // Partition the ORIGINAL cubic into certified one-phase pieces and
        // monotone crossing pieces. Protect both with linear Bernstein
        // constraints, including the original phase at crossing endpoints.
        // This keeps a remote extremum from becoming a second axis root.
        struct AxisPiece {std::array<double,4> initial,current;std::array<Vector,4> variation;int depth=0;};
        std::vector<AxisPiece> pending{{initial,current,v,0}};
        const auto protect=[&](const Vector& variation,double value,int phaseSign) {
            if(variation.squaredNorm()!=0){inequalities.push_back(-phaseSign*variation);bounds.push_back(phaseSign*value);}
            else if(phaseSign*value<0)throw std::runtime_error("fixed axis data violate their original phase chart");
        };
        while(!pending.empty()) {
            auto q=std::move(pending.back());pending.pop_back();
            const auto range=std::minmax_element(q.initial.begin(),q.initial.end());
            const int phaseSign=*range.first>=0&&*range.second>0?1:*range.second<=0&&*range.first<0?-1:0;
            if(phaseSign){for(int k=0;k<4;++k)protect(q.variation[k],q.current[k],phaseSign);continue;}
            std::array<double,3> derivative;for(int k=0;k<3;++k)derivative[k]=q.initial[k+1]-q.initial[k];
            const auto derivativeRange=std::minmax_element(derivative.begin(),derivative.end());
            const int direction=*derivativeRange.first>=0&&*derivativeRange.second>0?1:
                *derivativeRange.second<=0&&*derivativeRange.first<0?-1:0;
            if(direction) {
                for(int k=0;k<3;++k)protect(q.variation[k+1]-q.variation[k],q.current[k+1]-q.current[k],direction);
                for(int k:{0,3})if(q.initial[k]!=0)protect(q.variation[k],q.current[k],q.initial[k]>0?1:-1);
                continue;
            }
            if(q.depth==16)throw std::runtime_error("event inventory correction cannot partition a regular axis trace");
            AxisPiece lo,hi;lo.depth=hi.depth=q.depth+1;
            const auto split=[](const auto& a,auto& l,auto& r) {
                l[0]=a[0];l[1]=(a[0]+a[1])*.5;l[2]=(a[0]+2*a[1]+a[2])*.25;
                l[3]=(a[0]+3*a[1]+3*a[2]+a[3])*.125;
                r[0]=l[3];r[1]=(a[1]+2*a[2]+a[3])*.25;r[2]=(a[2]+a[3])*.5;r[3]=a[3];
            };
            split(q.initial,lo.initial,hi.initial);split(q.current,lo.current,hi.current);split(q.variation,lo.variation,hi.variation);
            pending.push_back(std::move(hi));pending.push_back(std::move(lo));
        }
    }
    
    Matrix g(inequalities.size(),n);Vector bound(bounds.size());for(int k=0;k<bound.size();++k){g.row(k)=inequalities[k];bound[k]=bounds[k];}
    const Vector solution=boundedMinimum(a,b,g,bound);
    InterfaceVolumeCorrection::Result result;result.rows=static_cast<int>(a.rows());result.columns=n;result.usedSuperLU=true;
    result.deltaValue.assign(fixed.size(),0);result.deltaGradient.assign(gradients.size(),0);
    for(int v=0;v<static_cast<int>(fixed.size());++v)if(!fixed[v])result.deltaValue[v]=h*solution[columns.at({v,0})];
    for(std::size_t k=0;k<gradients.size();++k)result.deltaGradient[k]=solution[columns.at({gradients[k].vertexIndex,static_cast<int>(gradients[k].component)})];
    return result;
}
}

void InterfaceConservativeGeometry::correctEventInventory(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceState& previous,
    const InterfaceConservativePhase& inventory,Quantum::Geometry::ActiveHermitePatch& patch,
    InterfaceConservativeGeometryReport& report)
{
    using Patch=Quantum::Geometry::ActiveHermitePatch;
    const auto require=[](bool ok,const char* message) {
        if(!ok)throw std::runtime_error(std::string("[InterfaceEventInventory] ")+message);
    };
    require(domain.uniformLeaves()&&domain.axisymmetric()&&mesh.minCoord.x==0,
        "interior-event correction currently requires uniform RZ leaves");
    require(previous.topology().liquids().size()==2&&previous.topology().surfaces().size()==2,
        "interior-event correction requires the declared two children");
    require(patch.cornerSectors.empty(),"event inventory correction has no contact-sector transfer");
    const double h=std::min(domain.spacing().x,domain.spacing().y);
    constexpr long double cutFractionTolerance=1e-6L;
    const auto reference=patch;auto work=patch;
    std::vector<std::uint8_t> base(mesh.cellIndexNum(),0),support=base,affected=base;
    std::vector<std::uint8_t> fixed(mesh.vertIndexNum(),1);
    for(int c:previous.eventCurvatureCells())base.at(c)=support.at(c)=1;
    // The same two-owner support used by event surface-force moments. Shared
    // corrections also affect its incident neighbours, accounted below.
    for(const auto face:domain.faces(mesh))if(face.negativeCell>=0&&face.positiveCell>=0
        &&(base[face.negativeCell]||base[face.positiveCell]))support[face.negativeCell]=support[face.positiveCell]=1;
    for(int c:domain.cellIndexByRow())if(support[c])for(int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert) {
        const auto& vertex=mesh.verts(label);const auto p=vertex.coord();
        require(vertex.supCell<0,"event correction cannot free a hanging node");
        // Preserve the entire original outer-boundary nodal trace. The
        // reflection condition is not an additional correction unknown.
        if(p.y==mesh.minCoord.y||p.y==mesh.maxCoord.y||p.x==mesh.maxCoord.x)continue;
        fixed[vertex.index()]=0;
    }
    std::vector<InterfaceVolumeCorrection::GradientDof> gradients;
    for(int v=0;v<mesh.vertIndexNum();++v)if(!fixed[v]) {
        if(mesh.verts(mesh.vertIndexToLabel[v]).coord().x!=0)gradients.push_back({v,1,1/h});
        gradients.push_back({v,2,1/h});
    }
    for(int c:domain.cellIndexByRow())for(int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert)
        if(!fixed[mesh.verts(label).index()]){affected[c]=1;break;}
    report.eventInventorySupportCells=std::count(affected.begin(),affected.end(),std::uint8_t(1));
    const auto geometry=[&](const Patch& data) {
        auto state=std::make_unique<InterfaceState>(mesh);
        InterfaceGeometryReconstruction::initialize(mesh,data,true,{},*state,previous.topology().symmetryPlanes());
        require(state->topology().liquids().size()==2&&state->topology().surfaces().size()==2,
            "local correction changed the two-child topology");
        return state;
    };
    auto state=geometry(work);const auto originalPhase=state->topology().cellPhase();
    for(int c:domain.cellIndexByRow())if(affected[c]&&originalPhase[c]==0)++report.eventInventoryCutCells;
    const auto residual=[&](const InterfaceCellMeasure::Measure& measure,bool record) {
        long double cut=0,full=0,pure=0;
        for(int c:domain.cellIndexByRow())if(affected[c]) {
            const auto error=std::abs(measure.negativeVolumeByCell[c]-inventory.liquidVolumes()[c])/inventory.capacities()[c];
            require(std::isfinite(error),"nonfinite local geometric/inventory discrepancy");
            full=std::max(full,error);
            if(originalPhase[c]==0)cut=std::max(cut,error);else pure=std::max(pure,error);
        }
        if(record){report.eventInventoryFinalCutError=static_cast<double>(cut);
            report.eventInventorySupportError=static_cast<double>(full);report.eventInventoryPureTailError=static_cast<double>(pure);}
        return cut;
    };
    const auto evaluate=[&](const InterfaceState& geometryState) {
        InterfacePhaseGeometry::MappedFacePartitions partitions(mesh,domain,geometryState.geometry());
        return InterfaceCellMeasure::evaluateRefined(mesh,domain,partitions,geometryState);
    };
    auto measure=evaluate(*state);
    for(int iteration=0;;++iteration) {
        const auto error=residual(measure,true);
        if(iteration==0)report.eventInventoryInitialCutError=static_cast<double>(error);
        if(error<=cutFractionTolerance)break;
        require(iteration<12,"local cut-volume iteration budget exhausted");
        measure.arcs.erase(std::remove_if(measure.arcs.begin(),measure.arcs.end(),
            [&](const auto& arc){return !affected[arc.cellIndex];}),measure.arcs.end());
        std::vector<long double> defects(inventory.liquidVolumes().size());
        for(int c:domain.cellIndexByRow())defects[c]=inventory.liquidVolumes()[c]-measure.negativeVolumeByCell[c];
        const auto delta=phaseBoundedStep(mesh,measure,defects,fixed,gradients,work,reference,originalPhase,h);
        double limit=0;for(double v:delta.deltaValue)limit=std::max(limit,std::abs(v));
        for(double v:delta.deltaGradient)limit=std::max(limit,h*std::abs(v));
        double alpha=limit>h/8?h/(8*limit):1;bool accepted=false;std::string rejection;
        for(int backtrack=0;backtrack<8;++backtrack,alpha*=.5)try {
            auto trial=work;
            for(int v=0;v<mesh.vertIndexNum();++v)trial.nodes[v].value+=alpha*delta.deltaValue[v];
            for(std::size_t k=0;k<gradients.size();++k) {
                auto& node=trial.nodes[gradients[k].vertexIndex];
                (gradients[k].component==1?node.dx:node.dy)+=alpha*delta.deltaGradient[k];
            }
            auto next=geometry(trial);
            // Preserving each original pure cell also retains both resolved
            // liquid identity anchors; final publication checks real lineage.
            for(int c:domain.cellIndexByRow())if(originalPhase[c]!=0)
                require(next->topology().cellPhase()[c]==originalPhase[c],"correction invaded an originally pure cell");
            auto nextMeasure=evaluate(*next);
            require(residual(nextMeasure,false)<error,"local cut-volume residual did not decrease");
            // This owned measure includes the accepted geometry's volume and
            // native sensitivities. Carry it into the next Newton iteration;
            // recomputing the same immutable state's full domain duplicates
            // the most expensive part of this local correction. Rejected
            // trials never replace the current state or its measure.
            measure=std::move(nextMeasure);
            work=std::move(trial);state=std::move(next);accepted=true;
            ++report.eventInventoryIterations;break;
        }catch(const std::exception& failure){rejection=failure.what();++report.eventInventoryBacktracks;}
        if(!accepted)throw std::runtime_error("[InterfaceEventInventory] local correction line search exhausted: "+rejection);
    }
    if(!report.eventInventoryIterations)return;
    // Keep the existing sparse saved-jet payload truthful after this last
    // correction. Source 6 denotes final event-volume data; old before/mixed
    // entries remain observations of the earlier reconstruction stages.
    std::vector<unsigned char> recorded(mesh.vertIndexNum(),0);
    for(auto& jet:report.nodalJets) {
        recorded[jet.vertex]=1;const auto& node=work.nodes[jet.vertex];
        if(jet.value!=node.value||jet.pure!=std::array<double,3>{node.dx,node.dy,node.dxy}) {
            jet.scalarCorrection+=node.value-jet.value;jet.source=6;
            jet.value=node.value;jet.pure={node.dx,node.dy,node.dxy};
        }
    }
    for(int v=0;v<mesh.vertIndexNum();++v)if(!recorded[v]) {
        const auto& node=work.nodes[v];
        if(node.dx==0&&node.dy==0&&node.dxy==0)continue;
        InterfaceConservativeGeometryReport::NodalJet jet;jet.vertex=v;jet.source=6;jet.value=node.value;
        jet.before=jet.mixed=jet.pure={node.dx,node.dy,node.dxy};
        jet.rejectedValues.fill(std::numeric_limits<double>::quiet_NaN());report.nodalJets.push_back(jet);
    }
    for(auto& jet:report.axisJets){jet.value=work.nodes[jet.vertex].value;jet.dyPure=work.nodes[jet.vertex].dy;}
    patch=std::move(work);
}
