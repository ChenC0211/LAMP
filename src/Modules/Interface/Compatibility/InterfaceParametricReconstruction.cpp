#include "Modules/Interface/Compatibility/InterfaceParametricProjection.h"
#include <Eigen/Dense>
#include <numeric>
#include <sstream>
#include <cstdlib>

namespace InterfaceParametricProjection {
namespace {
using namespace Quantum::Geometry::Parametric;
using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;
using Row=Eigen::RowVectorXd;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
struct Error {double squared=0,energy=0;std::vector<double> panels;};
Error error(const Curve& curve,const Source& source) {
    Error out;const auto& rule=gaussEight();
    for(const auto& arc:curve.arcs) {
        double local=0;const auto knots=source.partition(arc.a,arc.b);
        for(std::size_t i=1;i<knots.size();++i) {
            const auto exact=restrictedArc(source.arcs[source.index(knots[i-1])],knots[i-1],knots[i]);
            const auto fitted=restrictedArc(arc,knots[i-1],knots[i]);
            for(int q=0;q<8;++q){const double t=.5*(1+rule.x[q]),w=.5*(knots[i]-knots[i-1])*rule.w[q];
                const Point dd=exact.sampleFraction(t,2),defect=fitted.sampleFraction(t,2)-dd;
                local+=w*dot(defect,defect);out.energy+=w*dot(dd,dd);}
        }
        out.panels.push_back(local);out.squared+=local;
    }
    require(out.energy>0 && std::isfinite(out.squared),"nonfinite global reconstruction norm");return out;
}
void refineApproximation(std::vector<double>& knots,const Curve& curve,const Error& e) {
    std::vector<std::size_t> order(e.panels.size());std::iota(order.begin(),order.end(),0);
    std::sort(order.begin(),order.end(),[&](auto a,auto b){return e.panels[a]>e.panels[b];});
    double marked=0;for(auto index:order){const auto& arc=curve.arcs[index];const double middle=.5*(arc.a+arc.b);
        require(middle>arc.a && middle<arc.b,"continuous Hermite refinement exhausted material resolution");knots.push_back(middle);marked+=e.panels[index];if(marked>=.5*e.squared)break;}
    std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());
}
Curve prediction(const Source& source,const std::vector<double>& knots) {
    // The affine origin must belong to the SAME continuous space as its
    // correction. A broken panelwise H2 fit plus a C1 correction remains
    // broken, and its missing join areas disagree with phase-face transport.
    // Unique material nodes own both endpoint jets; volume constraints may
    // move them together in the global solve below.
    Curve out;out.boundary=source.boundary;for(std::size_t i=1;i<knots.size();++i){const double a=knots[i-1],b=knots[i],end=b==1&&!source.boundary.open()?0:b,h=b-a;
        const Point p=source.sample(a),delta=source.sample(end)-p,d0=source.sample(a,1)*h,d1=source.sample(end,1)*h;
        out.arcs.push_back({a,b,{{p,d0,delta*3-d0*2-d1,delta*(-2)+d0+d1,{}}}});}
    return out;
}
using ScalarPolynomial=std::array<double,5>;
const std::array<ScalarPolynomial,5> cardinal{{{{1,0,-3,2,0}},{{0,0,3,-2,0}},{{0,1,-2,1,0}},{{0,0,-1,1,0}},{{0,0,1,-2,1}}}};
std::array<int,5> coarseColumns(int i,int n,bool open=false){const int next=open?i+1:(i+1)%n;return {{3*i,3*next,3*i+1,3*next+1,3*i+2}};}
Matrix coarseCoordinates(int nodes) {
    const int d=3*nodes;const double h=1./nodes;Matrix stiffness=Matrix::Zero(d,d);Row mean=Row::Zero(d);
    for(int cell=0;cell<nodes;++cell){const auto ids=coarseColumns(cell,nodes);
        for(int i=0;i<5;++i){for(int k=0;k<5;++k)mean[ids[i]]+=h*cardinal[i][k]/(k+1);
            for(int j=0;j<5;++j)for(int k=2;k<5;++k)for(int l=2;l<5;++l)
                stiffness(ids[i],ids[j])+=cardinal[i][k]*cardinal[j][l]*k*(k-1)*l*(l-1)/((k+l-3)*h*h*h);}}
    // Remove the constant mode before Cholesky. No small mean penalty is
    // added to a matrix with an increasingly large finest-interval scale.
    Eigen::LLT<Matrix> factor(stiffness.bottomRightCorner(d-1,d-1));
    require(factor.info()==Eigen::Success,"coarse periodic Hermite metric is unresolved");
    Matrix transform=Matrix::Zero(d,d);for(int i=0;i<nodes;++i)transform(3*i,0)=1;
    transform.bottomRightCorner(d-1,d-1)=factor.matrixU().solve(Matrix::Identity(d-1,d-1));
    const Row average=mean*transform;
    for(int i=0;i<nodes;++i)for(int j=1;j<d;++j)transform(3*i,j)-=average[j];
    return transform;
}
Matrix boundaryCoordinates(int nodes,int component,const CurveBoundary& boundary) {
    // Homogeneous endpoint conditions in the same C1 correction space:
    // axis -> r=0 and z'=0; fixed rim -> both coordinate increments zero.
    // Only the two-axis axial coordinate retains a constant null mode.
    const int raw=3*nodes+2,d=3*nodes;const double h=1./nodes;
    Matrix stiffness=Matrix::Zero(raw,raw);Row mean=Row::Zero(raw);
    for(int cell=0;cell<nodes;cell++){const auto ids=coarseColumns(cell,nodes,true);
        for(int i=0;i<5;i++){for(int k=0;k<5;k++)mean[ids[i]]+=h*cardinal[i][k]/(k+1);
            for(int j=0;j<5;j++)for(int k=2;k<5;k++)for(int l=2;l<5;l++)
                stiffness(ids[i],ids[j])+=cardinal[i][k]*cardinal[j][l]*k*(k-1)*l*(l-1)/((k+l-3)*h*h*h);}}
    const int firstFixed=boundary.axis(0)&&component==1?1:0,lastFixed=3*nodes+(boundary.axis(1)&&component==1?1:0);
    std::vector<int> free;for(int i=0;i<raw;i++)if(i!=firstFixed&&i!=lastFixed)free.push_back(i);
    Matrix reduced(d,d);for(int i=0;i<d;i++)for(int j=0;j<d;j++)reduced(i,j)=stiffness(free[i],free[j]);
    const int first=boundary.constantMode(component)?1:0;Eigen::LLT<Matrix> factor(reduced.bottomRightCorner(d-first,d-first));
    require(factor.info()==Eigen::Success,"boundary Hermite metric is unresolved");
    Matrix result=Matrix::Zero(raw,d);const Matrix normalized=factor.matrixU().solve(Matrix::Identity(d-first,d-first));
    for(int i=first;i<d;i++)result.row(free[i]).tail(d-first)=normalized.row(i-first);
    if(first){for(int i=0;i<=nodes;i++)result(3*i,0)=1;const Row average=mean*result;
        for(int i=0;i<=nodes;i++)for(int j=1;j<d;j++)result(3*i,j)-=average[j];}
    return result;
}
using Detail=std::array<Arc,2>;
double detailEnergy(const Detail& a,const Detail& b) {
    double out=0;for(int side=0;side<2;++side)for(int q=0;q<5;++q)
        out+=.25*gw[q]*dot(a[side].sampleFraction(.5*(1+gx[q]),2),b[side].sampleFraction(.5*(1+gx[q]),2));
    return out;
}
const std::array<Detail,3>& detailCoordinates() {
    // Four midpoint/child-bubble functions minus the parent bubble leave
    // three new modes. Their end values and slopes vanish. H2 orthogonality
    // to the parent bubble implies zero mean by two integrations by parts,
    // hence orthogonality to every ancestor quartic, not only that bubble.
    static const auto modes=[](){
        std::array<Detail,4> seeds;for(auto& mode:seeds){mode[0].a=0;mode[0].b=.5;mode[1].a=.5;mode[1].b=1;}
        const auto set=[](Arc& arc,ScalarPolynomial values){for(int k=0;k<5;++k)arc.p[k].x=values[k];};
        set(seeds[0][0],{{0,0,3,-2,0}});set(seeds[0][1],{{1,0,-3,2,0}});
        set(seeds[1][0],{{0,0,-.5,.5,0}});set(seeds[1][1],{{0,.5,-1,.5,0}});
        set(seeds[2][0],{{0,0,1,-2,1}});set(seeds[3][1],{{0,0,1,-2,1}});
        Arc bubble;bubble.a=0;bubble.b=1;set(bubble,{{0,0,1,-2,1}});
        const Detail parent{{restrictedArc(bubble,0,.5),restrictedArc(bubble,.5,1)}};
        const auto subtract=[](Detail& a,const Detail& b,double scale){for(int side=0;side<2;++side)for(int k=0;k<5;++k)a[side].p[k]+=b[side].p[k]*(-scale);};
        std::array<Detail,3> out;int accepted=0;
        for(auto candidate:seeds){for(int pass=0;pass<2;++pass){subtract(candidate,parent,detailEnergy(candidate,parent)/detailEnergy(parent,parent));
                for(int j=0;j<accepted;++j)subtract(candidate,out[j],detailEnergy(candidate,out[j]));}
            const double energy=detailEnergy(candidate,candidate);if(energy<1e-20)continue;
            for(auto& arc:candidate)for(auto& p:arc.p)p=p/std::sqrt(energy);out[accepted++]=candidate;if(accepted==3)break;}
        require(accepted==3,"quartic Hermite split does not have three energy details");return out;
    }();return modes;
}
struct Space {
    int nodes,dofs;CurveBoundary boundary;std::vector<double> knots;std::vector<Matrix> shapes,radialShapes;
    // A leaf sees the coarse modes and its ancestor details only. Keep the
    // same global H2 coordinates without storing or evaluating zero columns.
    std::vector<std::vector<int>> columns;
    Space(std::vector<double> partition,int original,CurveBoundary ends):nodes(int(partition.size())-1),dofs(3*nodes),boundary(std::move(ends)),knots(std::move(partition)) {
        const Matrix coarse=boundary.open()?boundaryCoordinates(original,1,boundary):coarseCoordinates(original);const int initial=3*original;
        const Matrix radial=boundary.open()?boundaryCoordinates(original,0,boundary):Matrix();
        shapes.assign(nodes,Matrix::Zero(5,initial));
        columns.assign(nodes,std::vector<int>(initial));for(auto& ids:columns)std::iota(ids.begin(),ids.end(),0);
        if(boundary.open())radialShapes.assign(nodes,Matrix::Zero(5,initial));
        for(int i=0;i<nodes;++i){const int owner=std::min(original-1,int(.5*(knots[i]+knots[i+1])*original));const auto ids=coarseColumns(owner,original,boundary.open());
            for(int j=0;j<5;++j){Arc phi;phi.a=double(owner)/original;phi.b=double(owner+1)/original;for(int k=0;k<5;++k)phi.p[k].x=cardinal[j][k];
                const auto restricted=restrictedArc(phi,knots[i],knots[i+1]);
                for(int k=0;k<5;++k){shapes[i].row(k).head(initial)+=restricted.p[k].x*coarse.row(ids[j]);
                    if(boundary.open())radialShapes[i].row(k).head(initial)+=restricted.p[k].x*radial.row(ids[j]);}}
            if(boundary.constantMode(1)){shapes[i].col(0).setZero();shapes[i](0,0)=1;}
        }
        int added=initial;const auto& details=detailCoordinates();
        std::function<void(double,double)> split=[&](double a,double b){const auto next=std::upper_bound(knots.begin(),knots.end(),a);
            if(next==knots.end() || *next>=b)return;const double middle=.5*(a+b),scale=std::pow(b-a,1.5);
            // Second-derivative energy scales as amplitude^2 / length^3.
            require(std::binary_search(knots.begin(),knots.end(),middle),"material hierarchy is not nested dyadic refinement");const int first=added;added+=3;
            for(int i=index(a);i<nodes && knots[i]<b;++i){const int side=knots[i]<middle?0:1;
                const int local=int(columns[i].size());for(int j=0;j<3;j++)columns[i].push_back(first+j);
                shapes[i].conservativeResize(5,local+3);if(boundary.open())radialShapes[i].conservativeResize(5,local+3);
                for(int j=0;j<3;++j){auto phi=details[j][side];phi.a=side?middle:a;phi.b=side?b:middle;for(auto& p:phi.p)p=p*scale;
                    const auto restricted=restrictedArc(phi,knots[i],knots[i+1]);for(int k=0;k<5;++k){shapes[i](k,local+j)=restricted.p[k].x;
                        if(boundary.open())radialShapes[i](k,local+j)=restricted.p[k].x;}}}
            split(a,middle);split(middle,b);
        };
        for(int i=0;i<original;++i)split(double(i)/original,double(i+1)/original);
        require(added==dofs,"material hierarchy lost a Hermite degree of freedom");
    }
    int index(double s)const{return std::min(nodes-1,int(std::upper_bound(knots.begin(),knots.end(),s)-knots.begin())-1);}
    std::size_t storedCoefficients()const {
        std::size_t count=0;for(const auto& a:shapes)count+=a.size();for(const auto& a:radialShapes)count+=a.size();return count;
    }
    Row localBasis(int i,double t,int derivative=0,int component=1)const {
        Row out=Row::Zero(columns[i].size());const double h=knots[i+1]-knots[i];
        const auto& values=boundary.open()&&component==0?radialShapes:shapes;
        for(int k=4;k>=derivative;--k){double factor=1;for(int j=0;j<derivative;++j)factor*=(k-j)/h;out=out*t+factor*values[i].row(k);}return out;
    }
    template<class Dense>
    void addBasis(Dense& target,int offset,int i,double t,int derivative,int component,double weight)const {
        const auto local=localBasis(i,t,derivative,component);
        for(int k=0;k<int(columns[i].size());k++)target[offset+columns[i][k]]+=weight*local[k];
    }
    Curve curve(const Vector& coefficients)const {
        Curve out;out.boundary=boundary.homogeneous();for(int i=0;i<nodes;++i){Arc arc;arc.a=knots[i];arc.b=knots[i+1];
            Vector localX(columns[i].size()),localY(columns[i].size());for(int k=0;k<int(columns[i].size());k++){localX[k]=coefficients[columns[i][k]];localY[k]=coefficients[dofs+columns[i][k]];}
            const Vector x=(boundary.open()?radialShapes[i]:shapes[i])*localX,y=shapes[i]*localY;
            for(int k=0;k<5;++k)arc.p[k]={x[k],y[k]};out.arcs.push_back(arc);}return out;
    }
    Vector prolong(const Curve& old)const {
        Vector out=Vector::Zero(2*dofs);
        for(int i=0;i<nodes;++i){const auto arc=restrictedArc(old.arcs[old.index(knots[i])],knots[i],knots[i+1]);
            for(int q=0;q<5;++q){const double t=.5*(1+gx[q]),w=.5*(arc.b-arc.a)*gw[q];const auto p=arc.sampleFraction(t),dd=arc.sampleFraction(t,2);
                addBasis(out,0,i,t,2,0,w*dd.x);addBasis(out,dofs,i,t,2,1,w*dd.y);if(boundary.constantMode(0))out[0]+=w*p.x;if(boundary.constantMode(1))out[dofs]+=w*p.y;}}
        return out;
    }
};
Curve compose(const Curve& base,const Curve& correction) {
    std::vector<double> knots;for(const auto& a:base.arcs)knots.push_back(a.a);for(const auto& a:correction.arcs)knots.push_back(a.a);knots.push_back(1);
    std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());Curve out;out.boundary=base.boundary;
    for(std::size_t i=1;i<knots.size();++i){auto arc=restrictedArc(base.arcs[base.index(knots[i-1])],knots[i-1],knots[i]);
        const auto delta=restrictedArc(correction.arcs[correction.index(knots[i-1])],knots[i-1],knots[i]);
        for(int j=0;j<5;++j)arc.p[j]+=delta.p[j];out.arcs.push_back(arc);}return out;
}
struct Objective {Matrix mean;Vector gradient,meanTarget;};
Objective objective(const Curve& base,const Source& source,const Space& space) {
    const int d=space.dofs,means=int(space.boundary.constantMode(0))+int(space.boundary.constantMode(1));
    Objective out;out.gradient=Vector::Zero(2*d);out.mean=Matrix::Zero(means,2*d);out.meanTarget=Vector::Zero(means);
    int meanRows[2]={-1,-1},row=0;for(int c=0;c<2;c++)if(space.boundary.constantMode(c)){meanRows[c]=row;out.mean(row++,c*d)=1;}
    std::vector<double> knots=space.knots;for(const auto& a:base.arcs)knots.push_back(a.a);for(const auto& a:source.arcs)knots.push_back(a.a);
    std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());const auto& rule=gaussEight();
    for(std::size_t k=1;k<knots.size();++k){const double a=knots[k-1],b=knots[k];const int owner=space.index(a);const double h=space.knots[owner+1]-space.knots[owner];
        const auto exact=restrictedArc(source.arcs[source.index(a)],a,b);
        const auto predicted=restrictedArc(base.arcs[base.index(a)],a,b);
        for(int q=0;q<8;++q){const double t=.5*(1+rule.x[q]),weight=.5*(b-a)*rule.w[q],u=(a-space.knots[owner])/h+(b-a)/h*t;
            const auto defect=predicted.sampleFraction(t,2)-exact.sampleFraction(t,2),position=exact.sampleFraction(t)-predicted.sampleFraction(t);
            if(meanRows[0]>=0)out.meanTarget[meanRows[0]]+=weight*position.x;if(meanRows[1]>=0)out.meanTarget[meanRows[1]]+=weight*position.y;
            space.addBasis(out.gradient,0,owner,u,2,0,weight*defect.x);space.addBasis(out.gradient,d,owner,u,2,1,weight*defect.y);
        }
    }
    // The H2 null modes, when allowed by the physical endpoint conditions,
    // are fixed by exact material-mean constraints. All other coordinates
    // are H2 orthonormal, with no extra mean condition at a fixed rim.
    for(int c=0;c<2;c++)if(meanRows[c]>=0)out.gradient[c*d]=-out.meanTarget[meanRows[c]];return out;
}
using Intervals=std::vector<std::pair<double,double>>;
struct ConstraintSpaceDefect:std::runtime_error {
    int cell;
    Intervals support;
    ConstraintSpaceDefect(int owner,Intervals intervals):std::runtime_error("volume constraint has no resolved independent direction in this material space"),cell(owner),support(std::move(intervals)){}
};
struct Measures {std::vector<double> volume;std::map<int,Row> jacobian;std::map<int,Intervals> support;double total=0;std::size_t incidences=0;double incidenceWidth=0;};
Measures measure(const Curve& curve,const Grid& grid,const Space& space) {
    Measures out;out.volume.resize(grid.size());PhaseMeasure phase(curve);
    for(std::size_t k=0;k<grid.size();k++){const auto c=grid.rectangle(k);
        const double v=grid.axisymmetric?phase.rectanglePhaseMoments(c.left,c.right,c.bottom,c.top,true).values[0]:phase.simpleRectangle(c.left,c.right,c.bottom,c.top);out.volume[k]=v;out.total+=v;}
    out.incidences=phase.endpointIncidences;out.incidenceWidth=phase.maximumEndpointIncidenceWidth;
    const std::array<std::vector<double>,2> lines{grid.lines(0),grid.lines(1)};
    const DomainRectangle domain{grid.left,grid.left+grid.width,grid.bottom,grid.bottom+grid.height};
    for(const auto& arc:curve.arcs){const auto wet=rectangleIntervals(arc,domain);std::vector<double> cuts;
        for(const auto& ab:wet){cuts.push_back(ab[0]);cuts.push_back(ab[1]);}if(cuts.empty())continue;
        for(int axis=0;axis<2;++axis)for(std::size_t line=1;line+1<lines[axis].size();line++){
            auto p=arc.coordinate(axis);p[0]-=lines[axis][line];const auto r=roots(p,0,1);cuts.insert(cuts.end(),r.begin(),r.end());}
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());const int owner=space.index(arc.a);const double h=space.knots[owner+1]-space.knots[owner];
        for(std::size_t k=1;k<cuts.size();++k){const double a=cuts[k-1],b=cuts[k],middle=.5*(a+b);
            if(std::none_of(wet.begin(),wet.end(),[&](const auto& ab){return middle>=ab[0]&&middle<=ab[1];}))continue;
            const auto p=arc.sampleFraction(middle);
            const int slot=grid.locate(p);require(slot>=0,"domain-restricted volume derivative has no native cell");
            auto found=out.jacobian.try_emplace(slot,Row::Zero(2*space.dofs));auto& row=found.first->second;
            auto& support=out.support[slot];if(support.empty() || support.back().first!=arc.a)support.emplace_back(arc.a,arc.b);
            const auto accumulate=[&](double node,double w){const double t=.5*(a+b+(b-a)*node),weight=.5*(b-a)*(arc.b-arc.a)*w*(grid.axisymmetric?2*std::acos(-1.)*arc.sampleFraction(t).x:1.);
                const double u=(arc.a-space.knots[owner])/h+(arc.b-arc.a)/h*t;const auto tangent=arc.sampleFraction(t,1);
                space.addBasis(row,0,owner,u,0,0,weight*tangent.y);space.addBasis(row,space.dofs,owner,u,0,1,-weight*tangent.x);};
            if(grid.axisymmetric){const auto& rule=gaussEight();for(int q=0;q<8;q++)accumulate(rule.x[q],rule.w[q]);}
            else for(int q=0;q<5;q++)accumulate(gx[q],gw[q]);
        }
    }
    return out;
}
struct Constraint {Row row;double rhs=0;int cell=-1;};
struct QpResult {Vector value;std::size_t active=0;double stationarity=0;};
// The nonlinear merit uses the original FV bounds in volume-fraction units.
// Its linear prediction uses the SAME analytic volume rows as the QP. This
// separates a bad full Newton step from insufficient material directions.
double volumeViolation(const Measures& state,const Vector& coefficients,const Objective& objective,
    const Grid& grid,const std::vector<double>& inventory,const Options& options,const Vector* increment=nullptr) {
    double violation=0,total=state.total;
    for(std::size_t k=0;k<inventory.size();k++) {
        double volume=state.volume[k];const auto row=state.jacobian.find(int(k));
        if(increment&&row!=state.jacobian.end()){const double delta=row->second.dot(*increment);volume+=delta;total+=delta;}
        const double capacity=grid.capacity(k),radius=options.inventoryFractionUncertainty*capacity;
        const double lower=std::clamp(inventory[k]-radius,0.,capacity),upper=std::clamp(inventory[k]+radius,0.,capacity);
        violation+=std::max({0.,lower-volume,volume-upper})/capacity;
    }
    violation+=std::abs(total-std::accumulate(inventory.begin(),inventory.end(),0.))/grid.totalCapacity();
    violation+=(objective.mean*coefficients-objective.meanTarget).lpNorm<1>();
    return violation;
}
double dotRoundoff(const Row& row,const Vector& value) {
    const double u=std::numeric_limits<double>::epsilon(),nu=(row.size()+2)*u;
    return nu/(1-nu)*row.cwiseAbs().dot(value.cwiseAbs());
}
// One orthogonal factorization owns the primal face and its dual response.
// Reconstructing the face after a constraint enters/leaves avoids accumulating
// small primal errors through many dual steps. Such errors are not evidence
// that the geometric approximation space needs more degrees of freedom.
struct WorkingFace {
    Matrix rows,q,r;
    Vector rhs,value,multiplier;
    WorkingFace(Matrix constraints,Vector values,const Vector& free):rows(std::move(constraints)),rhs(std::move(values)) {
        const int n=int(free.size()),m=int(rows.rows());
        require(m<=n,"volume working face has too many independent constraints");
        if(!m){q=Matrix(n,0);r=Matrix(0,0);value=free;multiplier=Vector(0);return;}
        Eigen::HouseholderQR<Matrix> factor(rows.transpose());
        q=factor.householderQ()*Matrix::Identity(n,m);
        r=factor.matrixQR().topLeftCorner(m,m).triangularView<Eigen::Upper>();
        require(r.diagonal().array().abs().minCoeff()>0,"volume working face lost an independent constraint");
        Vector normalCoordinates=r.transpose().triangularView<Eigen::Lower>().solve(rows*free-rhs);
        value=free-q*normalCoordinates;
        // Refine the original constraints, not their normal equations.
        // Corrections remain in the same row space and preserve stationarity.
        for(int pass=0;pass<2;pass++){
            const Vector correction=r.transpose().triangularView<Eigen::Lower>().solve(rows*value-rhs);
            normalCoordinates+=correction;value-=q*correction;
        }
        multiplier=r.triangularView<Eigen::Upper>().solve(normalCoordinates);
        require(value.allFinite()&&multiplier.allFinite(),"nonfinite volume working-face solution");
    }
    void response(const Row& row,Vector& normal,Vector& dual)const {
        Vector coordinates=q.transpose()*row.transpose();normal=row.transpose()-q*coordinates;
        // Reorthogonalization resolves the independent part without subtracting
        // two large products involving ill-conditioned Lagrange multipliers.
        const Vector correction=q.transpose()*normal;coordinates+=correction;normal-=q*correction;
        dual=r.triangularView<Eigen::Upper>().solve(coordinates);
    }
};
QpResult project(const Measures& state,const Vector& current,const Objective& objective,const Grid& grid,
    const std::vector<double>& inventory,const Options& options) {
    const double total=std::accumulate(inventory.begin(),inventory.end(),0.);Row sum=Row::Zero(current.size());for(const auto& pair:state.jacobian)sum+=pair.second;
    std::vector<Constraint> equal,box;
    // A component wholly outside (or covering) the domain has no moving
    // interface in any cell. Its volume equality is already constant; retain
    // its inventory and test feasibility without dividing by a zero row.
    if(sum.squaredNorm()>0)equal.push_back({sum,total-state.total+sum.dot(current),-1});
    else require(std::abs(total-state.total)<=options.measureFractionReserve*grid.totalCapacity(),"constant domain volume disagrees with independent inventory");
    for(int i=0;i<objective.mean.rows();i++)equal.push_back({objective.mean.row(i),objective.meanTarget[i],-1});
    const auto normalize=[&](Constraint& c){const double scale=c.row.norm();require(scale>0 && std::isfinite(scale),"global volume constraint has no resolved response");c.row/=scale;c.rhs/=scale;};
    for(auto& c:equal)normalize(c);
    for(std::size_t k=0;k<inventory.size();++k){const double capacity=grid.capacity(k),radius=options.inventoryFractionUncertainty*capacity,reserve=options.measureFractionReserve*capacity;
        const double lower=std::clamp(inventory[k]-radius,0.,capacity),upper=std::clamp(inventory[k]+radius,0.,capacity);const auto found=state.jacobian.find(int(k));
        if(found==state.jacobian.end()||found->second.squaredNorm()==0){
            if(!(state.volume[k]>=lower-reserve && state.volume[k]<=upper+reserve)){std::ostringstream message;message<<std::setprecision(17)<<"global projection has unresolved pure-cell support: cell="<<k<<" geometry_fraction="<<state.volume[k]/capacity<<" inventory_fraction="<<inventory[k]/capacity;throw std::runtime_error(message.str());}continue;}
        const auto& row=found->second;const double shift=row.dot(current);
        // Subtract the neighboring physical volumes before adding the small
        // coefficient response. Otherwise their cancellation discards that
        // response even when the final bound is well resolved.
        box.push_back({row,(upper-state.volume[k])+shift,int(k)});normalize(box.back());
        box.push_back({-row,(state.volume[k]-lower)-shift,int(k)});normalize(box.back());
    }
    const Vector free=-objective.gradient;std::vector<int> active;
    const auto constraints=[&](){Matrix c(equal.size()+active.size(),current.size());Vector rhs(c.rows());
        for(int i=0;i<int(equal.size());++i){c.row(i)=equal[i].row;rhs[i]=equal[i].rhs;}
        for(int i=0;i<int(active.size());++i){c.row(equal.size()+i)=box[active[i]].row;rhs[equal.size()+i]=box[active[i]].rhs;}
        return WorkingFace(std::move(c),std::move(rhs),free);};
    WorkingFace face=constraints();Vector value=face.value,lambda=face.multiplier;
    // Dual feasible active set: a new violated bound is approached only until
    // an old multiplier reaches zero. Dropping that bound and continuing the
    // SAME dual step avoids cycling between full incompatible equality solves.
    int changes=0;
    for(int iteration=0;iteration<16*int((box.size()+1)*(box.size()+1));++iteration){
        int add=-1;double violation=0;
        for(int i=0;i<int(box.size());++i){if(std::find(active.begin(),active.end(),i)!=active.end())continue;
            const double actual=box[i].row.dot(value),v=actual-box[i].rhs;
            // A cancellation-prone dot product is resolved by the magnitudes
            // of its summands, not by its possibly vanishing final value.
            // This arithmetic bound does not spend the geometric error budget.
            const double tolerance=dotRoundoff(box[i].row,value)+4*std::numeric_limits<double>::epsilon()*std::abs(box[i].rhs);
            if(v>std::max(violation,tolerance)){violation=v;add=i;}}
        if(add<0){const double stationarity=(value+objective.gradient+face.rows.transpose()*lambda).norm();
            if(std::getenv("QUANTUM_HERMITE_RECONSTRUCTION_DIAGNOSTICS"))std::cout<<std::setprecision(17)
                <<"HERMITE_VOLUME_SOLVE changes="<<changes<<" active="<<active.size()<<" current_norm="<<current.norm()
                <<" value_norm="<<value.norm()<<" stationarity="<<stationarity<<std::endl;
            return {value,active.size(),stationarity};}
        if(std::getenv("QUANTUM_HERMITE_RECONSTRUCTION_DIAGNOSTICS")&&(iteration<8||(iteration&(iteration-1))==0)){
            const int cell=box[add].cell;std::cout<<std::setprecision(17)<<"HERMITE_VOLUME_ADD iteration="<<iteration<<" cell="<<cell
                <<" violation="<<violation<<" geometry_fraction="<<state.volume[cell]/grid.capacity(cell)
                <<" inventory_fraction="<<inventory[cell]/grid.capacity(cell)<<" response="<<state.jacobian.at(cell).norm()
                <<" value_norm="<<value.norm()<<" active="<<active.size()<<std::endl;}
        double pending=0;
        for(;;){require(++changes<16*int((box.size()+1)*(box.size()+1)),"global Hermite dual active set exhausted its iteration budget");
            Vector normal,r;face.response(box[add].row,normal,r);
            const double rankResolution=64*std::numeric_limits<double>::epsilon()*std::max(1.,face.rows.norm()*r.norm());
            const bool independent=normal.norm()>rankResolution;
            value=face.value;if(independent)value-=pending*normal;lambda=face.multiplier-pending*r;
            const double curvature=normal.squaredNorm(),distance=box[add].row.dot(value)-box[add].rhs;
            double full=independent?std::max(0.,distance)/curvature:std::numeric_limits<double>::infinity(),partial=std::numeric_limits<double>::infinity();int drop=-1;
            for(int i=0;i<int(active.size());++i)if(r[equal.size()+i]>0){const double limit=std::max(0.,lambda[equal.size()+i])/r[equal.size()+i];if(limit<partial){partial=limit;drop=i;}}
            // A numerically dependent violated constraint with no dual release
            // identifies insufficient local material directions. Enrich that
            // constraint's actual curve support within the SAME nested space;
            // do not invert its roundoff-size compliance or change its bound.
            if(!independent && !std::isfinite(partial)){
                if(std::getenv("QUANTUM_HERMITE_RECONSTRUCTION_DIAGNOSTICS")){
                    const int cell=box[add].cell;const double capacity=grid.capacity(cell);
                    std::cout<<std::setprecision(17)<<"HERMITE_VOLUME_DEPENDENCY cell="<<cell<<" bound="<<(add%2?"lower":"upper")
                        <<" geometry_fraction="<<state.volume[cell]/capacity<<" inventory_fraction="<<inventory[cell]/capacity
                        <<" unknowns="<<current.size()<<" active="<<active.size()<<" equalities="<<equal.size()<<" distance="<<distance
                        <<" linear_fraction_violation="<<distance*state.jacobian.at(cell).norm()/capacity
                        <<" normal="<<normal.norm()<<" rank_resolution="<<rankResolution
                        <<" face_residual="<<(face.rows*value-face.rhs).lpNorm<Eigen::Infinity>()
                        <<" dependent_rhs_defect="<<r.dot(face.rhs)-box[add].rhs<<std::endl;
                }
                throw ConstraintSpaceDefect(box[add].cell,state.support.at(box[add].cell));
            }
            const double step=std::min(full,partial);
            if(!std::isfinite(step)){std::ostringstream message;message<<std::setprecision(17)<<"global volume active-set step is nonfinite: cell="<<box[add].cell<<" unknowns="<<current.size()<<" active="<<active.size()<<" distance="<<distance<<" curvature="<<curvature<<" normal="<<normal.norm()<<" multiplier_response="<<(active.empty()?0.:r.tail(active.size()).maxCoeff());throw std::runtime_error(message.str());}
            pending+=step;
            if(full<=partial){active.push_back(add);face=constraints();value=face.value;lambda=face.multiplier;break;}
            require(drop>=0,"global dual bound has no owner");active.erase(active.begin()+drop);
            face=constraints();
        }
    }
    throw std::runtime_error("global Hermite volume active set did not converge");
}
double maximumDisplacement(const Curve& correction) {
    double out=0;for(const auto& arc:correction.arcs)for(double t:{0.,.25,.5,.75,1.})out=std::max(out,norm(arc.sampleFraction(t)));return out;
}
void auditVolumeDerivative(const Curve& curve,const Grid& grid,const Space& space,const Measures& state,int cell) {
    const auto& row=state.jacobian.at(cell);const double analytic=row.norm();if(!(analytic>0))return;
    const Curve direction=space.curve(row.transpose()/analytic);double bound=0;
    for(const auto& arc:direction.arcs){double value=0;for(const auto& p:arc.p)value+=norm(p);bound=std::max(bound,value);}
    if(!(bound>0))return;const auto box=grid.rectangle(cell);const double h=std::min(box.right-box.left,box.top-box.bottom);
    const auto volume=[&](double step){auto displacement=direction;for(auto& arc:displacement.arcs)for(auto& p:arc.p)p=p*step;
        const auto moved=compose(curve,displacement);PhaseMeasure phase(moved);
        return grid.axisymmetric?phase.rectanglePhaseMoments(box.left,box.right,box.bottom,box.top,true).values[0]
            :phase.simpleRectangle(box.left,box.right,box.bottom,box.top);};
    for(double relative:{1e-3,1e-4,1e-5}){const double step=relative*h/bound,finite=(volume(step)-volume(-step))/(2*step);
        std::cout<<std::setprecision(17)<<"HERMITE_VOLUME_DERIVATIVE cell="<<cell<<" relative_displacement="<<relative
            <<" analytic="<<analytic<<" finite_difference="<<finite<<" relative_difference="<<std::abs(finite-analytic)/analytic<<std::endl;}
}
void validateGeometryProblem(const Source& source,const Grid& grid,const Options& options){
    require(grid.nx>0 && grid.ny>0 && grid.width>0 && grid.height>0 && std::isfinite(grid.left) && std::isfinite(grid.bottom) && std::isfinite(grid.width) && std::isfinite(grid.height),"invalid global reconstruction grid");
    require(options.materialSegments>0 && options.materialSegments<=65536 && options.relativeH2>0 && std::isfinite(options.relativeH2) &&
        options.inventoryFractionUncertainty>0 && options.measureFractionReserve>0 && std::isfinite(options.inventoryFractionTolerance) && options.inventoryFractionTolerance>options.inventoryFractionUncertainty+options.measureFractionReserve,"invalid global reconstruction budget");
    require(options.maximumMaterialSegments>=options.materialSegments&&options.maximumMaterialSegments<=65536,"invalid material-space resource ceiling");
    require(!source.arcs.empty() && source.arcs.front().a==0 && source.arcs.back().b==1,"global reconstruction requires the complete material interval");
    source.validateBoundary(grid.left,grid.left+grid.width,grid.bottom,grid.bottom+grid.height,grid.axisymmetric);
    for(std::size_t i=0;i<source.arcs.size();++i){const auto& arc=source.arcs[i];require(arc.b>arc.a && (i==0 || source.arcs[i-1].b==arc.a),"noncontiguous source interval");for(const auto& p:arc.p)require(std::isfinite(p.x) && std::isfinite(p.y),"nonfinite source coefficient");}
    require((!grid.axisymmetric||grid.left>=0)&&(!(source.boundary.axis(0)||source.boundary.axis(1))||(grid.axisymmetric&&grid.left==0)),"inconsistent RZ material boundary and projection domain");
}
std::vector<double> coarseKnots(const Options& options){
    std::vector<double> knots;for(int i=0;i<=options.materialSegments;++i)knots.push_back(double(i)/options.materialSegments);return knots;
}
void retainKnots(std::vector<double>& knots,const Options& options){
    for(double s:options.minimumMaterialKnots){require(std::isfinite(s)&&s>=0&&s<=1,"invalid retained material knot");knots.push_back(s);}
    std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());
}
void checkSpaceSize(const std::vector<double>& knots,const Options& options,const char* cause,double error){
    if(knots.size()>std::size_t(options.maximumMaterialSegments)+1){std::ostringstream message;
        message<<std::setprecision(17)<<"global Hermite correction exhausted material-space budget: requested="<<knots.size()-1
            <<" maximum="<<options.maximumMaterialSegments<<" cause="<<cause<<" last_h2="<<error;
        throw MaterialResolutionLimit(message.str());}
}
struct FreeApproximation {Objective metric;Vector coefficients;Curve curve;Error defect;};
FreeApproximation freeApproximation(const Curve& base,const Source& source,const Space& space){
    FreeApproximation out;out.metric=objective(base,source,space);out.coefficients=-out.metric.gradient;
    out.curve=compose(base,space.curve(out.coefficients));out.defect=error(out.curve,source);return out;
}
}
std::vector<double> resolveMaterialSpace(const Source& source,const Grid& grid,const Options& options){
    validateGeometryProblem(source,grid,options);auto knots=coarseKnots(options);const Curve base=prediction(source,knots);retainKnots(knots,options);
    double relative=0;const char* cause="retained_space";
    for(;;){checkSpaceSize(knots,options,cause,relative);Space space(knots,options.materialSegments,source.boundary);
        const auto fit=freeApproximation(base,source,space);relative=std::sqrt(fit.defect.squared/fit.defect.energy);
        if(fit.defect.squared<=options.relativeH2*options.relativeH2*fit.defect.energy)return knots;
        refineApproximation(knots,fit.curve,fit.defect);cause="free_H2";
    }
}
Candidate prepare(const Source& source,const Grid& grid,const std::vector<double>& inventory,const Options& options) {
    validateGeometryProblem(source,grid,options);require(inventory.size()==grid.size(),"invalid global reconstruction inventory size");
    const double budget=options.inventoryFractionUncertainty+options.measureFractionReserve,total=std::accumulate(inventory.begin(),inventory.end(),0.);
    for(std::size_t k=0;k<inventory.size();++k){const double capacity=grid.capacity(k);if(!(std::isfinite(inventory[k]) && inventory[k]>=-budget*capacity && inventory[k]<=capacity*(1+budget))){std::ostringstream message;message<<std::setprecision(17)<<"independent inventory exceeds the declared physical uncertainty budget: cell="<<k<<" fraction="<<inventory[k]/capacity<<" fraction_budget="<<budget;throw std::runtime_error(message.str());}}
    auto knots=coarseKnots(options);
    const Curve base=prediction(source,knots);Curve previousCorrection;std::size_t refinements=0;
    // Keep the affine origin unchanged: enriching the admissible space must
    // not redefine the objective. Space checks the nested dyadic hierarchy.
    retainKnots(knots,options);
    const char* refinementCause="retained_space";double lastRelativeError=0;
    for(;;){checkSpaceSize(knots,options,refinementCause,lastRelativeError);
        Space space(knots,options.materialSegments,source.boundary);const auto fit=freeApproximation(base,source,space);
        const auto& metric=fit.metric;
        // Adding volume constraints cannot decrease the minimum H2 error.
        // An underresolved unconstrained space therefore needs enrichment
        // before any expensive nonlinear volume solve can yield acceptance.
        const auto& free=fit.coefficients;const auto& freeCurve=fit.curve;const auto& lowerBound=fit.defect;
        const bool diagnostics=std::getenv("QUANTUM_HERMITE_RECONSTRUCTION_DIAGNOSTICS")!=nullptr;
        lastRelativeError=std::sqrt(lowerBound.squared/lowerBound.energy);
        if(diagnostics)std::cout<<std::setprecision(17)<<"HERMITE_SPACE arcs="<<space.nodes<<" free_h2="<<lastRelativeError<<" source_energy="<<lowerBound.energy
            <<" basis_coefficients="<<space.storedCoefficients()<<" dense_equivalent="<<std::size_t(5)*space.dofs*space.nodes*(source.boundary.open()?2:1)<<std::endl;
        if(lowerBound.squared>options.relativeH2*options.relativeH2*lowerBound.energy){refinementCause="free_H2";refineApproximation(knots,freeCurve,lowerBound);++refinements;continue;}
        // Start the constrained solve at the minimizer of its own geometric
        // objective, not at the lower-order Hermite affine origin. The latter
        // can miss a resolved cut cell before any volume constraint is applied.
        Vector coefficients=previousCorrection.arcs.empty()?free:space.prolong(previousCorrection);Curve curve=compose(base,space.curve(coefficients));
        bool feasible=false,enrich=false;Report report;double meritPenalty=1;
        for(int iteration=0;iteration<12;++iteration){auto state=measure(curve,grid,space);QpResult solved;
            // A feasible unconstrained minimizer is already the solution.
            // Check with the FV uncertainty alone: no integration reserve is
            // spent to avoid an otherwise necessary geometric correction.
            // In particular, clipping a rounded negative zero to a hard zero
            // inside the QP must not move an already resolved stationary arc.
            bool freeFeasible=iteration==0&&previousCorrection.arcs.empty()
                &&std::abs(state.total-total)<=options.inventoryFractionUncertainty*grid.totalCapacity();
            for(std::size_t k=0;freeFeasible&&k<inventory.size();k++)
                freeFeasible=std::abs(state.volume[k]-inventory[k])<=options.inventoryFractionUncertainty*grid.capacity(k);
            try{if(freeFeasible)solved={coefficients,0,0};else solved=project(state,coefficients,metric,grid,inventory,options);}
            catch(const ConstraintSpaceDefect& defect){if(diagnostics){std::cout<<"HERMITE_SPACE volume_support="<<defect.support.size()<<std::endl;
                    if(refinements<8)auditVolumeDerivative(curve,grid,space,state,defect.cell);}
                for(const auto& interval:defect.support){const double middle=.5*(interval.first+interval.second);
                    require(middle>interval.first && middle<interval.second,"volume-space refinement exhausted material resolution");knots.push_back(middle);}
                previousCorrection=space.curve(coefficients);refinementCause="volume_support";enrich=true;break;}
            Measures actual;
            if(freeFeasible)actual=std::move(state);
            else {
                const Curve modelCurve=compose(base,space.curve(solved.value));
                const auto modelError=error(modelCurve,source);
                if(modelError.squared>options.relativeH2*options.relativeH2*modelError.energy) {
                    // This is the minimum-H2 solution of the current volume
                    // linearization. Do not chase nonlinear feasibility with
                    // a correction outside the declared geometric accuracy.
                    // Enrich the same nested space using that objective's
                    // existing panel indicator; retain the last accepted iterate.
                    if(diagnostics)std::cout<<std::setprecision(17)<<"HERMITE_SPACE linearized_h2="
                        <<std::sqrt(modelError.squared/modelError.energy)<<" budget="<<options.relativeH2<<std::endl;
                    refineApproximation(knots,modelCurve,modelError);previousCorrection=space.curve(coefficients);
                    refinementCause="linearized_H2";enrich=true;break;
                }
                const Vector direction=solved.value-coefficients;
                const double violation=volumeViolation(state,coefficients,metric,grid,inventory,options);
                const double predicted=volumeViolation(state,solved.value,metric,grid,inventory,options,&direction);
                const double decrease=violation-predicted;
                require(decrease>0,"linearized volume projection does not reduce the original constraint violation");
                const double slope=(coefficients-free).dot(direction),curvature=direction.squaredNorm();
                meritPenalty=std::max(meritPenalty,(slope+curvature)/(.9*decrease));
                const double meritSlope=slope-meritPenalty*decrease;
                const double objectiveBefore=.5*(coefficients-free).squaredNorm();
                bool accepted=false;double step=1;
                for(int search=0;search<32;search++,step*=.5) {
                    const Vector trial=coefficients+step*direction;
                    Curve trialCurve=compose(base,space.curve(trial));auto trialState=measure(trialCurve,grid,space);
                    const double trialViolation=volumeViolation(trialState,trial,metric,grid,inventory,options);
                    const double change=.5*(trial-free).squaredNorm()-objectiveBefore+meritPenalty*(trialViolation-violation);
                    if(change<=1e-4*step*meritSlope) {
                        if(diagnostics)std::cout<<std::setprecision(17)<<"HERMITE_VOLUME_STEP iteration="<<iteration
                            <<" scale="<<step<<" violation="<<violation<<" predicted="<<predicted<<" actual="<<trialViolation
                            <<" penalty="<<meritPenalty<<" step_norm="<<step*direction.norm()<<std::endl;
                        coefficients=trial;curve=std::move(trialCurve);actual=std::move(trialState);accepted=true;break;
                    }
                }
                require(accepted,"nonlinear volume merit line search did not resolve a descending step");
            }
            previousCorrection=space.curve(coefficients);
            report={};report.arcs=curve.arcs.size();report.refinements=refinements;report.activeCellConstraints=solved.active;report.stationarityResidual=solved.stationarity;
            report.basisCoefficients=space.storedCoefficients();report.denseBasisCoefficients=std::size_t(5)*space.dofs*space.nodes*(source.boundary.open()?2:1);
            report.maximumGeometryDisplacement=maximumDisplacement(previousCorrection);report.totalGeometryInventoryDifference=actual.total-total;
            report.endpointIncidences=actual.incidences;report.maximumEndpointIncidenceWidth=actual.incidenceWidth;
            for(std::size_t k=0;k<inventory.size();++k)report.maximumInventoryFractionError=std::max(report.maximumInventoryFractionError,std::abs(actual.volume[k]-inventory[k])/grid.capacity(k));
            feasible=report.maximumInventoryFractionError<=budget && std::abs(report.totalGeometryInventoryDifference)<=options.measureFractionReserve*grid.totalCapacity();

            if(feasible)break;
        }
        if(enrich){std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());++refinements;continue;}
        require(feasible,"global Hermite nonlinear volume projection did not meet its original measure budget");
        double coordinateEnergy=0;for(int c=0;c<2;c++){const int first=int(space.boundary.constantMode(c));coordinateEnergy+=coefficients.segment(c*space.dofs+first,space.dofs-first).squaredNorm();}
        const double representedEnergy=previousCorrection.energy();
        require(std::abs(coordinateEnergy-representedEnergy)<=1e-8*std::max({1e-30,coordinateEnergy,representedEnergy}),"hierarchical coordinates changed the material H2 objective");
        const auto e=error(curve,source);report.relativeH2=std::sqrt(e.squared/e.energy);lastRelativeError=report.relativeH2;
        if(diagnostics)std::cout<<"HERMITE_SPACE constrained_h2="<<report.relativeH2<<" fraction="<<report.maximumInventoryFractionError<<" displacement="<<report.maximumGeometryDisplacement<<std::endl;
        if(report.relativeH2<=options.relativeH2){
            Intervals unresolvedCharts;double sourceLongitudinal=1;
            for(std::size_t i=0;i<curve.arcs.size();++i){const auto& a=curve.arcs[i];const auto& b=curve.arcs[(i+1)%curve.arcs.size()];
                if(!curve.boundary.open()||i+1<curve.arcs.size()){
                    report.positionJoin=std::max(report.positionJoin,norm(a.sampleFraction(1)-b.p[0]));report.tangentJoin=std::max(report.tangentJoin,norm(a.sampleFraction(1,1)-b.sampleFraction(0,1)));}
                const auto chord=a.increment(a.a,a.b);require(dot(chord,chord)>0,"global projected arc has coincident endpoints");std::array<double,5> longitudinal{};
                for(int j=1;j<5;++j)longitudinal[j-1]=j*dot(a.p[j],chord)/dot(chord,chord);const double minimum=range(longitudinal)[0];
                report.minimumLongitudinalDerivative=std::min(report.minimumLongitudinalDerivative,minimum);
                if(minimum<=0){unresolvedCharts.emplace_back(a.a,a.b);const auto cuts=source.partition(a.a,a.b);
                    for(std::size_t k=1;k<cuts.size();k++){const auto part=restrictedArc(source.arcs[source.index(cuts[k-1])],cuts[k-1],cuts[k]);std::array<double,9> tangent{};
                        for(int j=1;j<9;j++)tangent[j-1]=j*dot(part.p[j],chord)/dot(chord,chord)*(a.b-a.a)/(part.b-part.a);
                        sourceLongitudinal=std::min(sourceLongitudinal,range(tangent)[0]);}}
            }
            double coordinateScale=1;for(const auto& arc:curve.arcs)for(const auto& p:arc.p)coordinateScale=std::max(coordinateScale,norm(p));
            if(!(report.positionJoin<=64*std::numeric_limits<double>::epsilon()*coordinateScale&&report.tangentJoin<2e-8)){
                std::ostringstream message;message<<std::setprecision(17)<<"global reconstructed material first jet is unresolved: arcs="<<curve.arcs.size()<<" open="<<curve.boundary.open()
                    <<" position="<<report.positionJoin<<" tangent="<<report.tangentJoin<<" longitudinal="<<report.minimumLongitudinalDerivative<<" h2="<<report.relativeH2;
                throw std::runtime_error(message.str());}
            if(!unresolvedCharts.empty()){
                // H2 accuracy is an integral requirement; it does not imply
                // that one material arc is a regular graph over its chord.
                // Resolve that chart in the SAME nested approximation space,
                // leaving the source, inventory and every budget unchanged.
                std::cout<<std::setprecision(17)<<"HERMITE_CHART_REFINEMENT arcs="<<curve.arcs.size()<<" charts="<<unresolvedCharts.size()
                    <<" projected_longitudinal="<<report.minimumLongitudinalDerivative<<" source_longitudinal="<<sourceLongitudinal<<" h2="<<report.relativeH2<<std::endl;
                for(const auto& interval:unresolvedCharts){const double middle=.5*(interval.first+interval.second);
                    require(middle>interval.first&&middle<interval.second,"regular material chart exhausted parameter resolution");knots.push_back(middle);}
                std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());refinementCause="regular_chart";++refinements;continue;
            }
            if(curve.boundary.axis(0))sampleSurfaceRZ(curve.arcs.front(),0,true);if(curve.boundary.axis(1))sampleSurfaceRZ(curve.arcs.back(),1,true);
            Candidate out;out.grid_=grid;out.inventory_=inventory;out.boundary_=std::make_shared<const Boundary>(std::vector<Boundary::Component>{{1,std::make_shared<const Curve>(std::move(curve))}});out.report_=report;return out;
        }
        // All old polynomials remain admissible under knot insertion. Refine
        // the largest contributors to the SAME curvature objective, not tiny
        // grid-cut intervals or one fixed-endpoint volume bubble at a time.
        refinementCause="constrained_H2";refineApproximation(knots,curve,e);++refinements;
    }
}
} // namespace InterfaceParametricProjection
