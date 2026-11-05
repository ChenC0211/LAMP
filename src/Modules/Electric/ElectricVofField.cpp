#include "Modules/Electric/ElectricVofField.h"
#include "Modules/Electric/ElectricVofChargeTransport.h"
#include "Solvers/Poisson/PolarizationPoissonFVSolverCellCentered.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <cfenv>
#include <cerrno>
#include <cstring>
#include <cstdint>
#include <memory>
#include <optional>
#include <exception>
#include <sstream>
#include <iomanip>
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#include <float.h>
#include <xmmintrin.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#endif

namespace LAMP::Electric
{
namespace
{
namespace FV=Geometry::Detail::FiniteVolume;
using Clock=LAMP::TimeLayer;
using Vof=Interface::GeometricVofPublication;
using Solver=PolarizationPoissonFVSolverCellCentered;
using Coefficients=std::map<int,long double>;
using Timer=std::chrono::steady_clock;
void require(bool ok,const char* reason)
{if(!ok)throw std::invalid_argument(std::string("[ElectricVofField] ")+reason);}
bool sameClock(Clock a,Clock b)
{return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
double elapsed(Timer::time_point begin){return std::chrono::duration<double>(Timer::now()-begin).count();}
double finite(long double value,const char* reason)
{
    const double result=static_cast<double>(value);
    require(std::isfinite(value)&&std::isfinite(result),reason);return result;
}
std::uint64_t chargeIdentity()
{
    static std::atomic<std::uint64_t> next{1};auto value=next.load(std::memory_order_relaxed);
    do{require(value!=std::numeric_limits<std::uint64_t>::max(),"charge identity space exhausted");}
    while(!next.compare_exchange_weak(value,value+1,std::memory_order_relaxed));return value;
}
struct CellGradient {Coefficients x,y;};
void add(Coefficients& to,const Coefficients& from,long double factor)
{for(const auto& entry:from)to[entry.first]+=factor*entry.second;}
std::vector<CellGradient> geometricGradients(const QuadtreeMesh& mesh,const FaceSegmentGraph& graph,
    const NativeLeafDomain& domain,ElectricVofReport& report)
{
    const auto n=static_cast<std::size_t>(mesh.cellIndexNum());
    std::vector<std::set<int>> neighbors(n);
    for(const auto& segment:graph.segments())if(segment.L>=0){
        const int k=mesh.cells(segment.K).index(),l=mesh.cells(segment.L).index();
        require(k>=0&&l>=0&&static_cast<std::size_t>(k)<n&&static_cast<std::size_t>(l)<n,"gradient native incidence outside complete mesh");
        if(domain.cellRowByIndex()[l]>=0)neighbors[k].insert(l);
    }
    std::vector<CellGradient> result(n);
    for(std::size_t i=0;i<n;++i){
        if(domain.cellRowByIndex()[i]<0)continue;
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[i]);const auto center=cell.centerCoord();
        const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        const long double scale=std::min(static_cast<long double>(hi.x)-lo.x,static_cast<long double>(hi.y)-lo.y);
        require(std::isfinite(scale)&&scale>0,"invalid actual gradient cell scale");
        long double xx=0,xy=0,yy=0;
        for(const int j:neighbors[i]){const auto point=mesh.cells(mesh.cellIndexToLabel[j]).centerCoord();
            const long double x=(static_cast<long double>(point.x)-center.x)/scale,y=(static_cast<long double>(point.y)-center.y)/scale;
            const auto distance=x*x+y*y;require(std::isfinite(distance)&&distance>0,"invalid native geometric gradient distance");
            const auto w=1/distance;xx+=w*x*x;xy+=w*x*y;yy+=w*y*y;
        }
        const auto det=xx*yy-xy*xy;
        require(std::isfinite(det)&&det>64*std::numeric_limits<double>::epsilon()*(xx+yy)*(xx+yy),
            "geometric native gradient lacks two independent directions");
        long double selfX=0,selfY=0;
        for(const int j:neighbors[i]){const auto point=mesh.cells(mesh.cellIndexToLabel[j]).centerCoord();
            const long double x=(static_cast<long double>(point.x)-center.x)/scale,y=(static_cast<long double>(point.y)-center.y)/scale;
            const auto w=1/(x*x+y*y),gx=w*(yy*x-xy*y)/(det*scale),gy=w*(xx*y-xy*x)/(det*scale);
            require(std::isfinite(gx)&&std::isfinite(gy),"nonfinite coefficient-independent gradient");
            result[i].x[j]=gx;result[i].y[j]=gy;selfX-=gx;selfY-=gy;
        }
        result[i].x[static_cast<int>(i)]=selfX;result[i].y[static_cast<int>(i)]=selfY;
        long double xX=0,xY=0,yX=0,yY=0;
        for(const auto& e:result[i].x){const auto p=mesh.cells(mesh.cellIndexToLabel[e.first]).centerCoord();xX+=e.second*(p.x-center.x);xY+=e.second*(p.y-center.y);}
        for(const auto& e:result[i].y){const auto p=mesh.cells(mesh.cellIndexToLabel[e.first]).centerCoord();yX+=e.second*(p.x-center.x);yY+=e.second*(p.y-center.y);}
        const auto defect=std::max({std::abs(xX-1),std::abs(xY),std::abs(yX),std::abs(yY-1)});
        report.maximumAffineGradientDefect=std::max(report.maximumAffineGradientDefect,finite(defect,"nonfinite affine gradient certificate"));
        require(defect<=1e-12L,"coefficient-independent native gradient failed its affine moments");
    }
    return result;
}
double harmonic(double a,double b,long double da,long double db)
{
    require(a>=0&&b>=0&&std::isfinite(a)&&std::isfinite(b)&&da>0&&db>0,"invalid face material resistance input");
    if(a==0||b==0)return 0;
    const long double least=std::min(a,b);
    const auto denominator=da*(least/a)+db*(least/b);
    require(std::isfinite(denominator)&&denominator>0,"unrepresentable face material resistance");
    const auto value=finite(least*((da+db)/denominator),"unrepresentable unique face material");
    require(value>0,"positive face material below supported scalar range");return value;
}
double mixture(const std::array<long double,2>& stock,long double capacity,double a,double b)
{
    // Constitutive P0 averaging only; authoritative independent stocks and
    // their finite closure certificate are never rewritten or normalized.
    const auto value=finite((static_cast<long double>(a)*stock[0]+static_cast<long double>(b)*stock[1])/capacity,
        "unrepresentable raw-stock material coefficient");
    require(value>=0,"negative raw-stock material coefficient");
    if((a>0&&stock[0]>0)||(b>0&&stock[1]>0))require(value>0,"positive material coefficient below supported scalar range");
    return value;
}
void validateDenseCharge(const QuadtreeMesh& mesh,const Vof& witness,const std::vector<double>& q)
{
    const auto& rows=witness.domain().cellRowByIndex();
    require(q.size()==static_cast<std::size_t>(mesh.cellIndexNum())&&rows.size()==q.size(),
        "integrated charge must cover dense native cells of the actual selected domain");
    for(std::size_t i=0;i<q.size();++i){
        require(std::isfinite(q[i]),"nonfinite integrated charge input");
        require(rows[i]>=0||q[i]==0,"inactive native cell cannot hold integrated charge");
    }
}
struct FaceOperator {Coefficients full,correction;long double constant=0;};

// Pure geometric stock: no mesh pointer, VOF/charge/clock or accepted-state
// owner. A live publication validation precedes every use of this identity.
struct FieldGeometryControl
{
    int rounding=-1;
    unsigned native=0,csr=0;
    bool valid=false;
    static FieldGeometryControl read() noexcept
    {
        FieldGeometryControl result;result.rounding=std::fegetround();
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        result.valid=result.rounding!=-1&&_controlfp_s(&result.native,0,0)==0;
        result.csr=_mm_getcsr()&~0x3fu;
#endif
        return result;
    }
    bool same(const FieldGeometryControl& other) const noexcept
    {return valid&&other.valid&&rounding==other.rounding&&native==other.native&&csr==other.csr;}
};
struct PreparedElectricFace
{
    NativeLeafFace native;
    long double negativeDistance=0,positiveDistance=0,dirichletCoefficient=0;
    Coefficients correction;
    std::vector<ElectricVofGradientCoefficient> full;
};
struct PreparedElectricSegment {int nativeRow=-1;double sign=0;};
struct FieldGeometryPreparation
{
    std::shared_ptr<const NativeMeshIdentity> identity;
    std::vector<int> selectedRows,selectedCells;
    bool axisymmetric=false;
    FieldGeometryControl control;
    std::vector<ElectricVofBoundary::Type> boundaryTypes;
    std::shared_ptr<const Solver::PreparedIntegratedStructure> integratedStructure;
    std::vector<CellGradient> gradients;
    std::vector<PreparedElectricFace> faces;
    std::vector<PreparedElectricSegment> segments;
    std::vector<std::size_t> incidenceOffsets,incidentFaces;
    std::size_t nonorthogonalFaces=0,axisFaces=0;
    double maximumAffineGradientDefect=0;
    bool matches(const NativeLeafDomain& domain,const std::vector<ElectricVofBoundary>& boundary,
        const FieldGeometryControl& current) const
    {
        if(identity!=domain.meshIdentity()||axisymmetric!=domain.axisymmetric()
            ||selectedRows!=domain.cellRowByIndex()||selectedCells!=domain.cellIndexByRow()
            ||boundaryTypes.size()!=boundary.size()||!control.same(current))return false;
        for(std::size_t row=0;row<boundary.size();++row)if(boundaryTypes[row]!=boundary[row].type)return false;
        return true;
    }
};
void validateBoundary(const NativeLeafFace& f,const ElectricVofBoundary& bc,bool axisymmetric)
{
    require(std::isfinite(bc.value),"nonfinite static electric boundary value");
    if(f.negativeCell<0&&f.positiveCell<0){
        require(bc.type==ElectricVofBoundary::Type::Interior&&bc.value==0,
            "both-excluded native face requires an Interior zero sentinel, not a physical boundary");return;
    }
    const bool interior=f.negativeCell>=0&&f.positiveCell>=0;
    require(interior?(bc.type==ElectricVofBoundary::Type::Interior&&bc.value==0):
        (bc.type==ElectricVofBoundary::Type::Dirichlet||bc.type==ElectricVofBoundary::Type::HomogeneousNeumann),"incorrect interior/exterior electric boundary kind");
    if(bc.type==ElectricVofBoundary::Type::HomogeneousNeumann)require(bc.value==0,"only homogeneous electric Neumann is supported");
    if(!interior&&f.area==0)require(axisymmetric&&f.center.x==0&&bc.type==ElectricVofBoundary::Type::HomogeneousNeumann&&bc.value==0,
        "zero-area true axis requires homogeneous Neumann");
}
std::shared_ptr<const FieldGeometryPreparation> prepareFieldGeometry(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const std::vector<ElectricVofBoundary>& boundary,const FieldGeometryControl& control)
{
    auto prepared=std::make_shared<FieldGeometryPreparation>();prepared->identity=domain.meshIdentity();
    prepared->selectedRows=domain.cellRowByIndex();prepared->selectedCells=domain.cellIndexByRow();
    prepared->axisymmetric=domain.axisymmetric();prepared->control=control;
    const auto n=static_cast<std::size_t>(mesh.cellIndexNum());
    auto graph=std::make_unique<FaceSegmentGraph>();
    if(domain.cellIndexByRow().size()==n)graph->build(mesh);
    else{
        // CellActiveFn is opaque even when this caller supplies a pure mask.
        // Its complete call invalidates every enclosing live-mesh proof.
        {NativeMeshIdentity::ValidationBatch::Pause pause;
            graph->build(mesh,[&domain](const QuadtreeMesh&,const Cell& cell){return domain.cellRowByIndex()[cell.index()]>=0;});}
        domain.validateMeshIdentity(mesh);
    }
    ElectricVofReport report;prepared->gradients=geometricGradients(mesh,*graph,domain,report);
    prepared->maximumAffineGradientDefect=report.maximumAffineGradientDefect;
    std::map<std::pair<bool,int>,int> nativeByEdge;prepared->faces.resize(boundary.size());prepared->boundaryTypes.reserve(boundary.size());
    for(int row=0;row<domain.nativeFaceCount();++row){auto& face=prepared->faces[row];face.native=domain.face(mesh,row);
        const auto& f=face.native;const auto& bc=boundary[row];prepared->boundaryTypes.push_back(bc.type);
        validateBoundary(f,bc,domain.axisymmetric());
        if(f.negativeCell<0&&f.positiveCell<0)continue;
        require(nativeByEdge.emplace(std::make_pair(f.horizontal,f.edgeLabel),row).second,"duplicate actual electric native edge");
        require(f.area>=0&&std::isfinite(f.area),"invalid actual physical electric face measure");
        FaceOperator op;const int axis=f.horizontal?1:0;
        const auto coord=[&](Vector2D p){return axis?p.y:p.x;};const auto tangent=[&](Vector2D p){return axis?p.x:p.y;};
        if(f.negativeCell>=0&&f.positiveCell>=0){
            const auto a=mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).centerCoord(),b=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).centerCoord();
            const long double da=static_cast<long double>(coord(f.center))-coord(a),db=static_cast<long double>(coord(b))-coord(f.center),distance=da+db;
            require(da>0&&db>0&&std::isfinite(distance),"invalid actual face half distances");
            face.negativeDistance=da;face.positiveDistance=db;
            const auto c=static_cast<long double>(f.area)/distance;
            op.full[f.negativeCell]+=c;op.full[f.positiveCell]-=c;
            const long double shift=static_cast<long double>(tangent(b))-tangent(a);
            if(shift!=0){++prepared->nonorthogonalFaces;const auto factor=c*shift/2;
                add(op.correction,axis?prepared->gradients[f.negativeCell].x:prepared->gradients[f.negativeCell].y,factor);
                add(op.correction,axis?prepared->gradients[f.positiveCell].x:prepared->gradients[f.positiveCell].y,factor);}
        }else{const int cell=f.negativeCell>=0?f.negativeCell:f.positiveCell;require(cell>=0,"electric exterior face has no actual owner");
            if(f.area==0)++prepared->axisFaces;
            else if(bc.type==ElectricVofBoundary::Type::Dirichlet){const long double sign=f.negativeCell>=0?1:-1;
                const auto center=mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord();const long double distance=sign*(static_cast<long double>(coord(f.center))-coord(center));
                require(std::isfinite(distance)&&distance>0,"invalid actual Dirichlet face distance");
                const auto c=sign*f.area/distance;op.full[cell]+=c;face.dirichletCoefficient=c;
                const long double shift=static_cast<long double>(tangent(f.center))-tangent(center);
                if(shift!=0)add(op.correction,axis?prepared->gradients[cell].x:prepared->gradients[cell].y,c*shift);
            }
        }
        add(op.full,op.correction,1);
        face.correction=std::move(op.correction);face.full.reserve(op.full.size());
        for(const auto& entry:op.full)face.full.push_back({entry.first,finite(entry.second,"unrepresentable shared gradient coefficient")});
    }
    prepared->segments.reserve(graph->segments().size());std::vector<unsigned> orientedCount(boundary.size());
    for(const auto& segment:graph->segments()){
        const auto found=nativeByEdge.find({segment.isHorizontal,segment.edgeLabel});require(found!=nativeByEdge.end(),"solver segment lacks actual native electric face");
        const int row=found->second;const auto& face=prepared->faces[row].native;const int k=mesh.cells(segment.K).index();
        require(k==face.negativeCell||k==face.positiveCell,"solver graph/native owner incidence differs");
        prepared->segments.push_back({row,k==face.negativeCell?1.:-1.});++orientedCount[row];
    }
    prepared->incidenceOffsets.resize(n+1);
    for(std::size_t row=0;row<prepared->faces.size();++row){const auto& face=prepared->faces[row].native;
        const unsigned expected=static_cast<unsigned>(face.negativeCell>=0)+static_cast<unsigned>(face.positiveCell>=0);
        require(orientedCount[row]==expected,"electric graph does not cover one shared physical face exactly");
        for(const int c:{face.negativeCell,face.positiveCell})if(c>=0)++prepared->incidenceOffsets[static_cast<std::size_t>(c)+1];
    }
    for(std::size_t i=0;i<n;++i)prepared->incidenceOffsets[i+1]+=prepared->incidenceOffsets[i];
    auto cursor=prepared->incidenceOffsets;prepared->incidentFaces.resize(prepared->incidenceOffsets.back());
    for(std::size_t row=0;row<prepared->faces.size();++row){const auto& face=prepared->faces[row].native;
        for(const int c:{face.negativeCell,face.positiveCell})if(c>=0)prepared->incidentFaces[cursor[static_cast<std::size_t>(c)]++]=row;
    }
    std::vector<Solver::FittedBoundaryType> types;types.reserve(prepared->segments.size());
    std::vector<std::vector<int>> correctionColumns(prepared->segments.size()),affineColumns(prepared->segments.size());
    for(std::size_t sid=0;sid<prepared->segments.size();++sid){const auto row=prepared->segments[sid].nativeRow;
        const auto type=boundary[row].type;
        types.push_back(type==ElectricVofBoundary::Type::Interior?Solver::FittedBoundaryType::Interior:
            type==ElectricVofBoundary::Type::Dirichlet?Solver::FittedBoundaryType::Dirichlet:Solver::FittedBoundaryType::HomogeneousNeumann);
        for(const auto& entry:prepared->faces[row].correction)correctionColumns[sid].push_back(entry.first);
        for(const auto& entry:prepared->faces[row].full)affineColumns[sid].push_back(entry.nativeCell);
    }
    // FaceSegmentGraph declares a destructor and has no move constructor.
    // Transfer its unique owner, rather than copying a nominal rvalue graph.
    prepared->integratedStructure=Solver::prepareIntegratedStructure(mesh,std::move(graph),types,
        correctionColumns,affineColumns,domain.axisymmetric());
    // A change to FP controls during construction cannot certify future reuse.
    if(!prepared->control.same(FieldGeometryControl::read()))prepared->control.valid=false;
    return prepared;
}
template<class T>using FieldPostArray=std::unique_ptr<T[]>;
template<class T>FieldPostArray<T> fieldPostArray(std::size_t count)
{return FieldPostArray<T>(new T[count]{});}
struct FieldPostBudget
{
    static constexpr std::size_t limit=64u*1024u*1024u;
    std::size_t used=0;
    void add(std::size_t count,std::size_t width)
    {if(used>limit||count>(limit-used)/width)throw std::length_error("optional electric postprocess scratch cap");
        used+=count*width;}
};
struct FieldPostFace
{
    double h=0,d=0,j=0;
};
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
struct FieldPostFloatingState
{
    std::fenv_t native{};unsigned csr=0;bool valid=false;
    bool capture() noexcept {csr=_mm_getcsr();valid=std::fegetenv(&native)==0;return valid;}
    void restore() const noexcept
    {if(valid&&std::fesetenv(&native)!=0)std::terminate();_mm_setcsr(csr);}
    bool enterFresh(const FieldPostFloatingState& from) const noexcept
    {auto clean=from.native;clean._Fe_stat=0;if(std::fesetenv(&clean)!=0)return false;
        _mm_setcsr(from.csr&~0x3fu);return true;}
};
bool fieldPostSameControl(const FieldPostFloatingState& a,const FieldPostFloatingState& b) noexcept
{
    auto x=a.native,y=b.native;x._Fe_stat=0;y._Fe_stat=0;
    return std::memcmp(&x,&y,sizeof x)==0&&(a.csr&~0x3fu)==(b.csr&~0x3fu);
}
bool fieldPostQualified(FieldPostFloatingState& caller) noexcept
{
#if defined(_M_FP_STRICT) && _M_FP_STRICT && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    if(omp_in_parallel()||omp_get_max_threads()<2||!caller.capture())return false;
    unsigned control=0;
    if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM
        ||(control&_MCW_RC)!=_RC_NEAR||std::fegetround()!=FE_TONEAREST
        ||(caller.csr&0x1f80u)!=0x1f80u||(caller.csr&0xe040u)!=0)return false;
#ifdef _EM_AMBIGUOUS
    if(control&_EM_AMBIGUOUS)return false;
#endif
    return sizeof(double)==8&&sizeof(long double)==8
        &&std::numeric_limits<double>::is_iec559&&std::numeric_limits<double>::digits==53
        &&std::numeric_limits<long double>::digits==53;
#else
    (void)caller;return false;
#endif
}
// Declare before all optional storage; arm before the candidate identity is
// created. Thus no new capture can fail after the no-fail publication point,
// and deletion/unlocking cannot leak allocator FP/errno into the caller.
struct FieldPostCleanup
{
    FieldPostFloatingState finalState;int error=0;bool armed=false;
    bool arm() noexcept
    {error=errno;armed=finalState.capture();errno=error;return armed;}
    ~FieldPostCleanup() noexcept
    {if(armed){finalState.restore();errno=error;}}
};
template<class Kernel>
bool fieldPostTeam(int team,Kernel&& kernel,unsigned& actualWorkers)
{
    FieldPostFloatingState prepared;const int parentErrno=errno;
    if(!prepared.capture())return false;
    std::atomic<int> failed{0};unsigned long raisedNative=0;unsigned raisedCsr=0;
#pragma omp parallel num_threads(team)
    {
        FieldPostFloatingState worker;const int workerErrno=errno;
        const bool captured=worker.capture(),ready=captured&&worker.enterFresh(prepared);
        if(!ready)failed.store(1,std::memory_order_relaxed);
        if(ready)errno=parentErrno;
#pragma omp single
        {actualWorkers=std::max(actualWorkers,static_cast<unsigned>(omp_get_num_threads()));}
        kernel(ready,failed);
        if(ready){
            const int numericalErrno=errno;FieldPostFloatingState after;
            const bool got=after.capture();
            if(!got||numericalErrno!=parentErrno||!fieldPostSameControl(after,prepared))
                failed.store(1,std::memory_order_relaxed);
#pragma omp critical(lamp_electric_post_floating_flags)
            {if(got){raisedNative|=after.native._Fe_stat;raisedCsr|=after.csr&0x3fu;}}
        }
        if(captured)worker.restore();errno=workerErrno;
    }
    if(failed.load(std::memory_order_relaxed))return false;
    auto merged=prepared.native;merged._Fe_stat|=raisedNative;
    if(std::fesetenv(&merged)!=0)std::terminate();
    _mm_setcsr(prepared.csr|raisedCsr);errno=parentErrno;return true;
}
template<class Kernel>
bool fieldPostRows(std::size_t count,int team,Kernel&& kernel,unsigned& actualWorkers)
{
    return fieldPostTeam(team,[&](bool ready,std::atomic<int>& failed){
#pragma omp for schedule(static) nowait
        for(std::ptrdiff_t row=0;row<static_cast<std::ptrdiff_t>(count);++row){
            if(!ready)continue;
            try{kernel(static_cast<std::size_t>(row));}catch(...){failed.store(1,std::memory_order_relaxed);}
        }
    },actualWorkers);
}
#endif
#if defined(_MSC_VER)
#pragma float_control(push)
#pragma float_control(precise,on)
#pragma float_control(except,on)
#pragma fenv_access(on)
#endif
std::optional<Timer::time_point> fieldPostTimingStart(bool& available) noexcept
{
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    FieldPostFloatingState current;const int error=errno;
    if(!current.capture()){errno=error;available=false;return {};}
    unsigned control=0;const bool masked=_controlfp_s(&control,0,0)==0&&(control&_MCW_EM)==_MCW_EM&&(current.csr&0x1f80u)==0x1f80u;
    if(!masked){current.restore();errno=error;available=false;return {};}
    const auto began=Timer::now();current.restore();errno=error;return began;
#else
    available=false;return {};
#endif
}
void fieldPostTimingFinish(const std::optional<Timer::time_point>& began,double& seconds,bool& available) noexcept
{
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    if(!began){available=false;return;}
    FieldPostFloatingState current;const int error=errno;
    if(!current.capture()){errno=error;available=false;return;}
    unsigned control=0;const bool masked=_controlfp_s(&control,0,0)==0&&(control&_MCW_EM)==_MCW_EM&&(current.csr&0x1f80u)==0x1f80u;
    if(!masked){current.restore();errno=error;available=false;return;}
    seconds+=std::chrono::duration<double>(Timer::now()-*began).count();
    current.restore();errno=error;
#else
    (void)began;(void)seconds;available=false;
#endif
}
struct FieldPostStageTiming
{
    double& seconds;bool& available;std::optional<Timer::time_point> began;
    FieldPostStageTiming(double& s,bool& a):seconds(s),available(a),began(fieldPostTimingStart(a)){}
    ~FieldPostStageTiming() noexcept{fieldPostTimingFinish(began,seconds,available);}
};
#if defined(_MSC_VER)
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif
// Pure bit resource bounds; no floating-point arithmetic or new numerical gate.
bool fieldPostRange(double value,int& least,int& greatest) noexcept
{
    std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof bits);bits&=UINT64_C(0x7fffffffffffffff);
    const auto fraction=bits&UINT64_C(0x000fffffffffffff);const int exponent=static_cast<int>(bits>>52);
    if(!exponent&&!fraction)return false;
    least=std::min(least,exponent?exponent-1023-52:-1074);
    greatest=std::max(greatest,exponent?exponent-1023:-1023);return true;
}
}
struct ElectricVofCharge::Data
{
    Vof meshWitness;
    Clock clock;
    std::vector<double> charge;
    std::uint64_t identity;
    Data(const Vof& witness,Clock t,std::vector<double> q):meshWitness(witness),clock(t),charge(std::move(q)),identity(chargeIdentity()){}
};
ElectricVofCharge::ElectricVofCharge(std::shared_ptr<const Data> data):data_(std::move(data)){}
ElectricVofCharge ElectricVofCharge::fromIntegrated(const QuadtreeMesh& mesh,const Vof& vof,const std::vector<double>& q)
{
    vof.validateFor(mesh,vof.identity(),vof.timeLayer());const auto clock=vof.timeLayer();clock.time();
    require(vof.geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic,"Eulerian charge requires actual LocalPlic mesh witness");
    require(clock.position!=Clock::Position::Midpoint,"initial charge input requires accepted Initial/End source");
    require(vof.cells().size()==q.size(),"integrated charge/VOF dense native size mismatch");
    validateDenseCharge(mesh,vof,q);
    return ElectricVofCharge(std::make_shared<Data>(vof,clock,q));
}
ElectricVofCharge ElectricVofCharge::fromPrepared(const QuadtreeMesh& mesh,const Vof& witness,Clock clock,std::vector<double> q)
{
    witness.validateMesh(mesh);clock.time();
    require(witness.geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic&&q.size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "prepared charge requires actual LocalPlic native mesh");
    validateDenseCharge(mesh,witness,q);
    return ElectricVofCharge(std::make_shared<Data>(witness,clock,std::move(q)));
}
std::uint64_t ElectricVofCharge::identity() const{require(bool(data_),"consumed charge owner");return data_->identity;}
Clock ElectricVofCharge::timeLayer() const{require(bool(data_),"consumed charge owner");return data_->clock;}
const std::vector<double>& ElectricVofCharge::integratedCharge() const{require(bool(data_),"consumed charge owner");return data_->charge;}
void ElectricVofCharge::validateFor(const QuadtreeMesh& mesh,Clock clock) const
{require(bool(data_),"consumed charge owner");clock.time();require(sameClock(clock,data_->clock),"full charge clock mismatch");data_->meshWitness.validateMesh(mesh);}
void ElectricVofCharge::validateFor(const QuadtreeMesh& mesh,const Vof& geometry,Clock clock) const
{
    validateFor(mesh,clock);
    geometry.validateFor(mesh,geometry.identity(),clock);
    validateSelectedDomain(geometry.domain());
    require(data_->meshWitness.identity()==geometry.identity()&&sameClock(data_->meshWitness.timeLayer(),clock),
        "charge requires its actual geometry owner and complete geometry clock");
}
void ElectricVofCharge::validateForInBatch(const QuadtreeMesh& mesh,Clock clock,
    const Detail::ElectricVofValidationBatch& batch) const
{
    require(bool(data_),"consumed charge owner");clock.time();require(sameClock(clock,data_->clock),"full charge clock mismatch");
    if(!batch.sameBinding(mesh,data_->meshWitness.domain()))data_->meshWitness.validateMesh(mesh);
}
void ElectricVofCharge::validateSelectedDomain(const NativeLeafDomain& domain) const
{
    require(bool(data_),"consumed charge owner");const auto& actual=data_->meshWitness.domain();
    require(actual.axisymmetric()==domain.axisymmetric()&&actual.cellRowByIndex()==domain.cellRowByIndex(),
        "electric charge source and consumer select different actual cells or physical measure");
}
void Detail::ElectricVofValidationBatch::rememberValidated(const QuadtreeMesh& mesh,const NativeLeafDomain& domain)
{if(!identity_){mesh_=&mesh;identity_=domain.meshIdentity();}}
bool Detail::ElectricVofValidationBatch::sameBinding(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) const
{return mesh_==&mesh&&identity_&&domain.meshIdentity()==identity_;}
bool ElectricVofCharge::sameInstance(const ElectricVofCharge& other) const noexcept{return data_&&data_==other.data_;}
struct ElectricVofField::Data
{
    Vof geometry;
    ElectricVofCharge source,candidate;
    std::shared_ptr<const ElectricVofAdvectedCharge> advection;
    std::vector<double> implicitInput;
    Clock candidateClock;
    double dt=0;
    std::vector<double> phi,epsilon,conductivity;
    std::vector<std::array<double,3>> potentialTerms;
    std::shared_ptr<const FieldGeometryPreparation> preparedGeometry;
    std::vector<Vector2D> potentialGradients;
    CellGradientRecovery cellGradientRecovery=CellGradientRecovery::Recover;
    std::vector<ElectricVofFace> faces;
    ElectricVofOhmicTransferReceipt ohmicReceipt;
    ElectricVofReport report;
    Data(const Vof& g,const ElectricVofCharge& q,Clock t,const std::vector<double>& input,const ElectricVofAdvectedCharge* adv):
        geometry(g),source(q),candidate(q),advection(adv?std::make_shared<ElectricVofAdvectedCharge>(*adv):nullptr),
        implicitInput(input),candidateClock(t){}
};
ElectricVofField::ElectricVofField(std::shared_ptr<const Data> data):data_(std::move(data)){}
ElectricVofField ElectricVofField::prepareElectrostaticSnapshot(const QuadtreeMesh& mesh,const Vof& geometry,
    const ElectricVofCharge& source,const ElectricVofMaterials& material,
    const std::vector<ElectricVofBoundary>& boundary,const ElectricVofControls& controls,CellGradientRecovery recovery,const ElectricVofField* preparedGeometrySource)
{
    NativeMeshIdentity::ValidationBatch meshValidation(mesh);
    return prepareImpl(mesh,geometry,source,geometry.timeLayer(),source.integratedCharge(),nullptr,material,boundary,controls,true,nullptr,recovery,preparedGeometrySource);
}
ElectricVofField ElectricVofField::prepareOhmic(const QuadtreeMesh& mesh,const Vof& geometry,
    const ElectricVofCharge& source,Clock candidateClock,const ElectricVofMaterials& material,
    const std::vector<ElectricVofBoundary>& boundary,const ElectricVofControls& controls,CellGradientRecovery recovery,const ElectricVofField* preparedGeometrySource)
{
    NativeMeshIdentity::ValidationBatch meshValidation(mesh);
    return prepareImpl(mesh,geometry,source,candidateClock,source.integratedCharge(),nullptr,material,boundary,controls,false,nullptr,recovery,preparedGeometrySource);
}
ElectricVofField ElectricVofField::prepareOhmic(const QuadtreeMesh& mesh,const Vof& geometry,
    const ElectricVofAdvectedCharge& advection,const ElectricVofMaterials& material,
    const std::vector<ElectricVofBoundary>& boundary,const ElectricVofControls& controls,CellGradientRecovery recovery,const ElectricVofField* preparedGeometrySource)
{
    NativeMeshIdentity::ValidationBatch meshValidation(mesh);
    Detail::ElectricVofValidationBatch batch;
    advection.validateForInBatch(mesh,advection.sourceCharge(),advection.transport(),advection.timeLayer(),batch);
    require(geometry.identity()==advection.transport().midpoint().identity()||geometry.identity()==advection.transport().endpoint().identity(),
        "Ohmic advection entry requires actual retained transport geometry");
    return prepareImpl(mesh,geometry,advection.sourceCharge(),advection.timeLayer(),advection.integratedCharge(),&advection,material,boundary,controls,false,&batch,recovery,preparedGeometrySource);
}
ElectricVofField ElectricVofField::prepareImpl(const QuadtreeMesh& mesh,const Vof& geometry,
    const ElectricVofCharge& source,Clock candidateClock,const std::vector<double>& input,const ElectricVofAdvectedCharge* advection,
    const ElectricVofMaterials& material,const std::vector<ElectricVofBoundary>& boundary,const ElectricVofControls& controls,
    bool electrostaticSnapshot,Detail::ElectricVofValidationBatch* offeredBatch,CellGradientRecovery recovery,const ElectricVofField* preparedGeometrySource)
{
    NativeMeshIdentity::ValidationBatch meshValidation(mesh);
    const auto started=Timer::now();geometry.validateFor(mesh,geometry.identity(),geometry.timeLayer());
    Detail::ElectricVofValidationBatch localBatch;auto& batch=offeredBatch?*offeredBatch:localBatch;
    batch.rememberValidated(mesh,geometry.domain());
    require(geometry.geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic,"electric field requires genuine LocalPlic material source");
    const auto geometryClock=geometry.timeLayer(),old=source.timeLayer();old.time();candidateClock.time();geometryClock.time();
    source.validateForInBatch(mesh,old,batch);
    if(electrostaticSnapshot){
        require((geometryClock.position==Clock::Position::Initial||geometryClock.position==Clock::Position::End)&&
            sameClock(geometryClock,old)&&sameClock(candidateClock,old),
            "electrostatic snapshot requires one complete accepted Initial/End geometry/charge/field clock");
        require(source.data_->meshWitness.identity()==geometry.identity()&&
            sameClock(source.data_->meshWitness.timeLayer(),old),
            "electrostatic snapshot requires the actual accepted charge geometry owner/fullclock");
        require(advection==nullptr,"electrostatic snapshot cannot relabel an advected receipt as accepted charge");
    }else{
    require(old.position!=Clock::Position::Midpoint&&old.step<std::numeric_limits<std::int64_t>::max()&&candidateClock.step==old.step+1&&
        candidateClock.begin==old.time()&&candidateClock.end>candidateClock.begin&&
        (candidateClock.position==Clock::Position::Midpoint||candidateClock.position==Clock::Position::End),
        "Ohmic candidate requires accepted Begin and a complete next Midpoint/End interval");
    require(geometryClock.step==candidateClock.step&&geometryClock.begin==candidateClock.begin&&geometryClock.end==candidateClock.end&&
        (geometryClock.position==Clock::Position::Midpoint||geometryClock.position==Clock::Position::End),
        "frozen material geometry must retain the same complete candidate interval");
    }
    const auto duration=electrostaticSnapshot?0:candidateClock.time()-old.time();
    if(!electrostaticSnapshot)require(std::isfinite(duration)&&duration>0,"invalid actual Ohmic stage duration");
    require(std::isfinite(controls.relativeLinearTolerance)&&controls.relativeLinearTolerance>0&&controls.relativeLinearTolerance<=1e-11&&
        std::isfinite(controls.relativeConservationTolerance)&&controls.relativeConservationTolerance>0&&controls.relativeConservationTolerance<=1e-10,
        "cannot loosen original electric linear/conservation tolerances");
    require(controls.backend==ElectricVofControls::Backend::Gmres||controls.backend==ElectricVofControls::Backend::SuperLU,
        "unknown strict integrated electric backend");
    require(std::isfinite(material.epsilonLiquid)&&material.epsilonLiquid>0&&std::isfinite(material.epsilonGas)&&material.epsilonGas>0&&
        std::isfinite(material.conductivityLiquid)&&material.conductivityLiquid>=0&&std::isfinite(material.conductivityGas)&&material.conductivityGas>=0,
        "permittivity must be positive and conductivity nonnegative; no gas floor");
    require(input.size()==source.integratedCharge().size(),"implicit charge input/native layout mismatch");
    for(const auto value:input)require(std::isfinite(value),"nonfinite implicit integrated charge input");
    require(recovery==CellGradientRecovery::Recover||recovery==CellGradientRecovery::Omit,"unknown cell potential gradient recovery request");
    const bool recoverCellGradients=recovery==CellGradientRecovery::Recover;
    auto data=std::make_shared<Data>(geometry,source,candidateClock,input,advection);data->cellGradientRecovery=recovery;
    data->dt=duration;data->ohmicReceipt.duration=duration;
    data->report.usedAdvection=advection!=nullptr;
    data->report.electrostaticSnapshot=electrostaticSnapshot;
    const auto n=static_cast<std::size_t>(mesh.cellIndexNum());const auto& domain=geometry.domain();
    require(geometry.plicCells().size()==n&&source.integratedCharge().size()==n&&domain.cellRowByIndex().size()==n,
        "electric/charge source requires the actual dense native selected domain");
    source.validateSelectedDomain(domain);
    validateDenseCharge(mesh,geometry,input);
    require(boundary.size()==static_cast<std::size_t>(domain.nativeFaceCount()),"electric boundary does not cover actual native faces");
    data->report.cells=domain.cellIndexByRow().size();data->report.nativeFaces=boundary.size();data->epsilon.resize(n);data->conductivity.resize(n);
    const auto materialCell=[&](std::size_t i){
        if(domain.cellRowByIndex()[i]<0)return 0.;
        const auto& cell=geometry.plicCells()[i];const auto capacity=geometry.cells()[i].capacity;
        require(std::isfinite(capacity)&&capacity>0,"invalid native electric material capacity");
        data->epsilon[i]=mixture(cell.inventory,capacity,material.epsilonLiquid,material.epsilonGas);
        data->conductivity[i]=mixture(cell.inventory,capacity,material.conductivityLiquid,material.conductivityGas);
        require(data->epsilon[i]>0,"native permittivity must remain strictly positive");
        const auto defect=std::abs((cell.inventory[0]+cell.inventory[1])-capacity)/capacity;
        return finite(defect,"nonfinite material closure report");
    };
    const auto literalMaterial=[&]{data->report.maximumPhaseCapacityDefect=0;
        for(std::size_t i=0;i<n;++i)data->report.maximumPhaseCapacityDefect=std::max(data->report.maximumPhaseCapacityDefect,materialCell(i));};
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    FieldPostFloatingState materialEntry;const int materialErrno=errno;const bool materialQualified=fieldPostQualified(materialEntry);
    if(materialEntry.valid)materialEntry.restore();errno=materialErrno;
    bool materialParallel=false;
    if(materialQualified)try{
        FieldPostCleanup cleanup;
        const int team=static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(omp_get_max_threads()),std::max<std::size_t>(1,n)));
        FieldPostBudget budget;budget.add(static_cast<std::size_t>(team),sizeof(double)+16384u);
        auto maxima=fieldPostArray<double>(static_cast<std::size_t>(team));
        materialEntry.restore();errno=materialErrno;unsigned workers=0;
        if(!fieldPostRows(n,team,[&](std::size_t row){auto& maximum=maxima[static_cast<std::size_t>(omp_get_thread_num())];
            maximum=std::max(maximum,materialCell(row));},workers))throw std::runtime_error("optional electric material worker unavailable");
        for(int worker=0;worker<team;++worker)data->report.maximumPhaseCapacityDefect=std::max(data->report.maximumPhaseCapacityDefect,maxima[static_cast<std::size_t>(worker)]);
        if(!cleanup.arm())throw std::runtime_error("optional electric material FP receipt");materialParallel=true;
    }catch(...){materialEntry.restore();errno=materialErrno;}
    if(!materialParallel)literalMaterial();
#else
    literalMaterial();
#endif
    const auto geometryControl=FieldGeometryControl::read();
    if(preparedGeometrySource&&preparedGeometrySource->data_){const auto& offered=preparedGeometrySource->data_->preparedGeometry;
        if(offered&&offered->matches(domain,boundary,geometryControl))data->preparedGeometry=offered;
    }
    if(!data->preparedGeometry)data->preparedGeometry=prepareFieldGeometry(mesh,domain,boundary,geometryControl);
    const auto& prepared=*data->preparedGeometry;const auto& graph=prepared.integratedStructure->graph();const auto& gradients=prepared.gradients;
    data->report.maximumAffineGradientDefect=prepared.maximumAffineGradientDefect;
    data->report.nonorthogonalFaces=prepared.nonorthogonalFaces;data->report.axisFaces=prepared.axisFaces;
    data->faces.resize(boundary.size());
    const auto faceRow=[&](std::size_t row){const auto& geometric=prepared.faces[row];const auto& f=geometric.native;const auto& bc=boundary[row];
        validateBoundary(f,bc,domain.axisymmetric());auto& face=data->faces[row];face.nativeRow=static_cast<int>(row);
        if(f.negativeCell<0&&f.positiveCell<0)return;
        face.negativeCell=f.negativeCell;face.positiveCell=f.positiveCell;face.physicalMeasure=f.area;
        if(f.negativeCell>=0&&f.positiveCell>=0){
            face.epsilon=harmonic(data->epsilon[f.negativeCell],data->epsilon[f.positiveCell],geometric.negativeDistance,geometric.positiveDistance);
            face.conductivity=harmonic(data->conductivity[f.negativeCell],data->conductivity[f.positiveCell],geometric.negativeDistance,geometric.positiveDistance);
        }else{const int cell=f.negativeCell>=0?f.negativeCell:f.positiveCell;
            face.epsilon=data->epsilon[cell];face.conductivity=data->conductivity[cell];
        }
        face.implicitCoefficient=finite(static_cast<long double>(face.epsilon)+static_cast<long double>(duration)*face.conductivity,
            "unrepresentable per-face epsilon+dt*sigma");
        require(face.implicitCoefficient>0,"implicit face coefficient must remain strictly positive");
        long double constant=0;
        if(f.area>0&&bc.type==ElectricVofBoundary::Type::Dirichlet)constant=-geometric.dirichletCoefficient*bc.value;
        face.gradientConstant=finite(constant,"unrepresentable shared gradient boundary constant");
        face.gradientCoefficients=geometric.full;
    };
    std::vector<double> aBySegment(graph.segments().size());std::vector<Solver::FittedBoundary> fitted(graph.segments().size());
    std::vector<Solver::SegmentFluxCorrection> correction(graph.segments().size());
    std::vector<Solver::IntegratedFluxFactor> integratedFlux(graph.segments().size());
    const auto segmentRow=[&](std::size_t sid){const auto& segment=prepared.segments[sid];
        const int row=segment.nativeRow;const auto& face=data->faces[row];const double sign=segment.sign;const auto& bc=boundary[row];
        aBySegment[sid]=face.implicitCoefficient;
        fitted[sid].type=bc.type==ElectricVofBoundary::Type::Interior?Solver::FittedBoundaryType::Interior:
            bc.type==ElectricVofBoundary::Type::Dirichlet?Solver::FittedBoundaryType::Dirichlet:Solver::FittedBoundaryType::HomogeneousNeumann;
        fitted[sid].value=bc.value;
        for(const auto& entry:prepared.faces[row].correction){
            const auto canonical=finite(static_cast<long double>(face.implicitCoefficient)*entry.second,"unrepresentable corrected implicit coefficient");
            correction[sid].coefficients.push_back({entry.first,sign*canonical});}
        // No independent material correction or post-split re-equilibration:
        // all geometric correction coefficients are scaled by this SAME a_f.
        correction[sid].constant=0;
        auto& factor=integratedFlux[sid];factor.primary=face.epsilon;
        factor.secondary=face.conductivity;factor.weight=duration;
        factor.secondaryTransfer=Solver::IntegratedFluxFactor::SecondaryTransfer::MinimumBinary64QuantumNearestEven;
        factor.affine.constant=sign*face.gradientConstant;
        for(const auto& coefficient:face.gradientCoefficients)
            factor.affine.coefficients.push_back({coefficient.nativeCell,sign*coefficient.value});
    };
    const auto literalNumeric=[&]{for(std::size_t row=0;row<boundary.size();++row)faceRow(row);
        for(std::size_t sid=0;sid<prepared.segments.size();++sid)segmentRow(sid);};
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    FieldPostFloatingState numericEntry;const int numericErrno=errno;const bool numericQualified=fieldPostQualified(numericEntry);
    if(numericEntry.valid)numericEntry.restore();errno=numericErrno;
    bool numericParallel=false;
    if(numericQualified)try{
        FieldPostCleanup cleanup;
        const int team=static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(omp_get_max_threads()),
            std::max<std::size_t>(1,std::max(boundary.size(),prepared.segments.size()))));
        FieldPostBudget budget;budget.add(static_cast<std::size_t>(team),16384u);unsigned workers=0;
        if(!fieldPostTeam(team,[&](bool ready,std::atomic<int>& failed){
            const auto rows=[&](std::size_t count,const auto& kernel){
#pragma omp for schedule(static)
                for(std::ptrdiff_t row=0;row<static_cast<std::ptrdiff_t>(count);++row){
                    if(!ready||failed.load(std::memory_order_relaxed))continue;
                    try{kernel(static_cast<std::size_t>(row));}catch(...){failed.store(1,std::memory_order_relaxed);}
                }
            };
            rows(boundary.size(),faceRow);rows(prepared.segments.size(),segmentRow);
        },workers))throw std::runtime_error("optional electric coefficient worker unavailable");
        if(!cleanup.arm())throw std::runtime_error("optional electric coefficient FP receipt");numericParallel=true;
    }catch(...){
        // No solver or candidate identity has consumed these private rows.
        // Vector payloads append in the shared kernel, so clear before replay.
        for(auto& value:correction){value.coefficients.clear();value.constant=0;}
        for(auto& value:integratedFlux)value.affine.coefficients.clear();
        numericEntry.restore();errno=numericErrno;
    }
    if(!numericParallel)literalNumeric();
#else
    literalNumeric();
#endif
    data->report.geometrySeconds=elapsed(started);
    CellField<double> rhs,phi;rhs.initialize(mesh,0.);phi.initialize(mesh,0.);
    for(std::size_t i=0;i<n;++i)rhs.value[i]=input[i];
    const bool superlu=controls.backend==ElectricVofControls::Backend::SuperLU;
    Solver solver;Solver::Options options;options.useSuperLU=superlu;options.requireSuperLU=superlu;options.useGMRES=!superlu;options.useCG=false;options.useIC=false;
    options.relTol=controls.relativeLinearTolerance;options.axisymmetric=domain.axisymmetric();options.reuseFittedFactorization=false;
    Solver::SolveReport solved;
    const bool success=solver.solvePreparedIntegratedBoundary(mesh,*prepared.integratedStructure,rhs,phi,aBySegment,fitted,correction,options,&solved,&integratedFlux,&data->potentialTerms);
    require(success,solved.message.c_str());
    require(solved.success&&(superlu?(solved.usedSuperLU&&!solved.usedGMRES):(solved.usedGMRES&&!solved.usedSuperLU)),
        "strict integrated field solve did not use the configured backend exclusively");
    require(solved.N==static_cast<int>(domain.cellIndexByRow().size()),"electric solve did not retain the actual compact active unknown count");
    data->report.usedSuperLU=solved.usedSuperLU;data->report.relativeLinearResidual=solved.finalResidual;
    data->report.usedGMRES=solved.usedGMRES;data->report.linearIterations=solved.gmresIterations;
    data->report.linearRestarts=solved.gmresRestarts;data->report.linearMaximumWorkers=solved.gmresMaximumWorkers;
    data->report.assembleSeconds=solved.timeAssemble_s;data->report.factorSeconds=solved.timeFactor_s;data->report.solveSeconds=solved.timeSolve_s;
    data->report.integratedExpandedResidualAvailable=solved.integratedExpandedResidualAvailable;
    data->report.integratedExpandedResidualParallelUsed=solved.integratedExpandedResidualParallelUsed;
    data->report.integratedExpandedResidualParallelFallback=solved.integratedExpandedResidualParallelFallback;
    data->report.integratedExpandedResidualWorkers=solved.integratedExpandedResidualWorkers;
    data->report.integratedExpandedResidualRows=solved.integratedExpandedResidualRows;
    data->report.integratedExpandedResidualScratchBytes=solved.integratedExpandedResidualScratchBytes;
    data->report.integratedExpandedResidualTimingAvailable=solved.integratedExpandedResidualTimingAvailable;
    data->report.integratedExpandedResidualRouteSeconds=solved.integratedExpandedResidualRouteSeconds;
    data->report.integratedExpandedResidualSeconds=solved.integratedExpandedResidualSeconds;
    require(data->potentialTerms.size()==n,"native potential expansion does not cover actual mesh");
    data->phi.resize(n);for(std::size_t i=0;i<n;++i){data->phi[i]=phi.value[i];require(std::isfinite(data->phi[i]),"nonfinite native electric potential");
        FV::Dyadic exact;for(double component:data->potentialTerms[i]){require(std::isfinite(component),"nonfinite native potential term");exact.add(component);}
        require(FV::roundDouble(exact).value==data->phi[i],"potential projection differs from actual solved expansion");}
    bool postTimingAvailable=true;double postFaceSeconds=0,postCellSeconds=0,postFinalizeSeconds=0,postSeconds=0;
    const auto postBegan=fieldPostTimingStart(postTimingAvailable);const auto baseReport=data->report;
    const double gaussTolerance=controls.relativeConservationTolerance;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    FieldPostFloatingState postprocessEntry;const int postprocessErrno=errno;const bool postprocessQualified=fieldPostQualified(postprocessEntry);
    if(postprocessEntry.valid)postprocessEntry.restore();errno=postprocessErrno;
#endif
    const auto postprocess=[&](bool parallel) {
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        FieldPostCleanup cleanup;
#endif
        const auto faceCount=data->faces.size();std::size_t teamCount=1;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        if(parallel)teamCount=std::min<std::size_t>(static_cast<std::size_t>(omp_get_max_threads()),std::max<std::size_t>(1,std::max(n,faceCount)));
#endif
        const int team=static_cast<int>(teamCount);FieldPostBudget budget;
        const auto nextCount=electrostaticSnapshot?0u:n;
        if(parallel){
            require(n<=static_cast<std::size_t>(std::numeric_limits<int>::max())&&faceCount<=static_cast<std::size_t>(std::numeric_limits<int>::max()),"optional electric postprocess index range");
            budget.add(faceCount,sizeof(FieldPostFace));
            if(!electrostaticSnapshot)budget.add(faceCount,sizeof(FV::Dyadic)+128u);
            budget.add(n,sizeof(double)+(recoverCellGradients?sizeof(Vector2D):0u));budget.add(nextCount,sizeof(double));
            budget.add(teamCount+1,1024u*1024u);
            budget.add(1024u,1u); // array/vector metadata and allocation allowance
        }
        auto faceValues=fieldPostArray<FieldPostFace>(faceCount);
        auto transfers=electrostaticSnapshot?FieldPostArray<FV::Dyadic>{}:fieldPostArray<FV::Dyadic>(faceCount);
        auto gaussRelative=fieldPostArray<double>(n);
        const auto& offsets=prepared.incidenceOffsets;const auto& incident=prepared.incidentFaces;
        std::vector<double> next(nextCount);std::vector<Vector2D> recovered(recoverCellGradients?n:0u);
        if(parallel&&(next.capacity()!=nextCount||recovered.capacity()!=(recoverCellGradients?n:0u)))throw std::runtime_error("optional electric vector allocation contract");
        unsigned actualWorkers=0;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        if(parallel){if(errno!=postprocessErrno)throw std::runtime_error("optional electric setup errno");postprocessEntry.restore();errno=postprocessErrno;}
#endif
        const auto rows=[&](std::size_t count,auto&& kernel){
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
            if(parallel){if(!fieldPostRows(count,team,kernel,actualWorkers))throw std::runtime_error("optional electric row worker unavailable");return;}
#else
            (void)team;
#endif
            for(std::size_t row=0;row<count;++row)kernel(row);
        };
        {
            FieldPostStageTiming timing(postFaceSeconds,postTimingAvailable);
            rows(faceCount,[&](std::size_t row){const auto& face=data->faces[row];auto& value=faceValues[row];
                if(face.negativeCell<0&&face.positiveCell<0)return;
                FV::Dyadic h(face.gradientConstant);
                for(const auto& coefficient:face.gradientCoefficients)for(double component:data->potentialTerms[coefficient.nativeCell])h.add(FV::Dyadic::product({coefficient.value,component}));
                value.h=finite(h.rounded(),"unrepresentable native integrated electric field");
                value.d=finite(static_cast<long double>(face.epsilon)*value.h,"unrepresentable shared native displacement");
                value.j=face.conductivity==0?0:finite(static_cast<long double>(face.conductivity)*value.h,"unrepresentable shared native Ohmic current");
                // The finite shared transfer is the law, not a diagnostic projection.
                if(!electrostaticSnapshot)transfers[row]=FV::Dyadic::product({duration,value.j}).roundedBinary64Quantum();
            });
        }
        {
            FieldPostStageTiming timing(postCellSeconds,postTimingAvailable);
            rows(n,[&](std::size_t i){
                FV::Dyadic netDisplacement,integratedCurrent,gaussScale,gaussBudget;
                int least=std::numeric_limits<int>::max(),greatest=std::numeric_limits<int>::min();
                if(parallel)fieldPostRange(input[i],least,greatest);
                for(std::size_t p=offsets[i];p<offsets[i+1];++p){const auto row=incident[p];const auto& face=data->faces[row];const auto& f=faceValues[row];const int sign=static_cast<int>(i)==face.negativeCell?1:-1;
                    netDisplacement.add(f.d,sign);gaussScale.add(std::abs(f.d));
                    if(!electrostaticSnapshot){integratedCurrent.add(transfers[row],sign);
                        if(parallel){int a=std::numeric_limits<int>::max(),b=std::numeric_limits<int>::min(),c=a,d=b;
                            const bool hasDuration=fieldPostRange(duration,a,b),hasCurrent=fieldPostRange(f.j,c,d);
                            if(hasDuration&&hasCurrent){least=std::min(least,a+c);greatest=std::max(greatest,b+d+1);}}}
                }
                gaussBudget=gaussScale.multiplied(controls.relativeConservationTolerance);
                const auto oldQ=input[i];double nextQ=oldQ;
                // No transfer preserves the actual old payload, including signed zero.
                // A snapshot retains its accepted owner and needs no next-Q storage.
                if(!electrostaticSnapshot&&integratedCurrent.compare(FV::Dyadic{})!=0){
                    auto exact=FV::Dyadic(oldQ);exact.add(integratedCurrent,-1);nextQ=exact.roundedDouble();}
                auto residual=netDisplacement;residual.add(nextQ,-1);gaussScale.add(std::abs(nextQ));
                gaussBudget.add(FV::Dyadic::product({gaussTolerance,std::abs(nextQ)}));
                if(residual.absolute().compare(gaussBudget)>0){
                    // Failure-only source context. The exact physical gate and all
                    // represented stocks, fields and acceptance arithmetic are intact.
                    std::ostringstream context;context<<std::setprecision(std::numeric_limits<double>::max_digits10);
                    const auto center=mesh.cells(mesh.cellIndexToLabel[i]).centerCoord();
                    const auto putClock=[&](const char* label,Clock clock){context<<' '<<label<<"=("<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<')';};
                    context<<"[ElectricVofField] raw physical displacement Gauss balance exceeds original electric bound"
                        <<" native_cell="<<i<<" active_row="<<domain.cellRowByIndex()[i]<<" center=("<<center.x<<','<<center.y<<')'
                        <<" geometry="<<geometry.identity().value()<<" charge_geometry="<<source.data_->meshWitness.identity().value()<<" charge="<<source.identity()
                        <<" snapshot="<<electrostaticSnapshot<<" duration="<<duration;
                    putClock("geometry_clock",geometryClock);putClock("source_clock",old);putClock("candidate_clock",candidateClock);
                    context<<" epsilon="<<data->epsilon[i]<<" conductivity="<<data->conductivity[i]<<" phi="<<data->phi[i]
                        <<" phi_terms=("<<data->potentialTerms[i][0]<<','<<data->potentialTerms[i][1]<<','<<data->potentialTerms[i][2]<<')'
                        <<" input_Q="<<oldQ<<" next_Q="<<nextQ<<" net_D="<<netDisplacement.rounded()
                        <<" residual="<<residual.rounded()<<" scale="<<gaussScale.rounded()<<" budget="<<gaussBudget.rounded()
                        <<" tolerance="<<gaussTolerance<<" solve_relative="<<solved.finalResidual;
                    for(std::size_t row=0;row<faceCount;++row){auto face=data->faces[row];
                        if(face.negativeCell!=static_cast<int>(i)&&face.positiveCell!=static_cast<int>(i))continue;
                        face.integratedField=faceValues[row].h;face.displacement=faceValues[row].d;face.current=faceValues[row].j;
                        context<<" face=("<<face.nativeRow<<",neg="<<face.negativeCell<<",pos="<<face.positiveCell
                            <<",area="<<face.physicalMeasure<<",epsilon="<<face.epsilon<<",H="<<face.integratedField
                            <<",D="<<face.displacement<<",J="<<face.current<<",constant="<<face.gradientConstant<<",coeffs=";
                        for(const auto& coefficient:face.gradientCoefficients)context<<coefficient.nativeCell<<':'<<coefficient.value<<'@'<<data->phi[coefficient.nativeCell]<<';';
                        context<<')';
                    }
                    throw std::invalid_argument(context.str());
                }
                const auto scale=gaussScale.rounded();const auto relative=scale==0?0:finite(residual.absolute().rounded()/scale,"nonfinite Gauss report");
                gaussRelative[i]=relative;
                if(parallel&&least!=std::numeric_limits<int>::max()){
                    unsigned carry=0;for(std::size_t count=offsets[i+1]-offsets[i]+4;count;count>>=1)++carry;
                    if(static_cast<long long>(greatest)-least+carry+4>1024)throw std::length_error("optional electric cell exact span");
                }
                if(!electrostaticSnapshot)next[i]=nextQ;
            });
        }
        FieldPostStageTiming finalTiming(postFinalizeSeconds,postTimingAvailable);
        auto report=baseReport;
        for(std::size_t i=0;i<n;++i)report.maximumRelativeGaussResidual=std::max(report.maximumRelativeGaussResidual,gaussRelative[i]);
        if(recoverCellGradients)rows(n,[&](std::size_t i){FV::Dyadic gx,gy;
            for(const auto& coefficient:gradients[i].x)for(double component:data->potentialTerms[coefficient.first])gx.add(FV::Dyadic::product({coefficient.second,component}));
            for(const auto& coefficient:gradients[i].y)for(double component:data->potentialTerms[coefficient.first])gy.add(FV::Dyadic::product({coefficient.second,component}));
            const double x=FV::roundDouble(gx).value,y=FV::roundDouble(gy).value;
            require(std::isfinite(x)&&std::isfinite(y),"unrepresentable actual cell LS potential gradient");recovered[i]={x,y};});
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        if(parallel&&!cleanup.arm())throw std::runtime_error("optional electric final FP receipt");
        if(!parallel)cleanup.arm();
#endif
        ElectricVofCharge candidate=source;
        if(!electrostaticSnapshot)candidate=ElectricVofCharge(std::make_shared<ElectricVofCharge::Data>(geometry,candidateClock,std::move(next)));
        // All physical checks precede identity creation; publication cannot throw.
        static_assert(noexcept(data->candidate=std::move(candidate)),"candidate publication must not throw");
        static_assert(noexcept(data->potentialGradients=std::move(recovered)),"gradient publication must not throw");
        static_assert(noexcept(data->report=report),"report publication must not throw");
        for(std::size_t row=0;row<faceCount;++row){auto& face=data->faces[row];const auto& f=faceValues[row];face.integratedField=f.h;face.displacement=f.d;face.current=f.j;}
        data->candidate=std::move(candidate);
        if(recoverCellGradients)data->potentialGradients=std::move(recovered);
        report.postprocessParallelUsed=parallel&&actualWorkers>1;report.postprocessWorkers=parallel?actualWorkers:0;report.postprocessScratchBytes=parallel?budget.used:0;data->report=report;
    };
    bool postprocessComplete=false;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    if(postprocessQualified){try{postprocess(true);postprocessComplete=true;}catch(...){postprocessEntry.restore();errno=postprocessErrno;}}
#endif
    if(!postprocessComplete){postprocess(false);
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        data->report.postprocessParallelFallback=postprocessQualified;
#endif
    }
    fieldPostTimingFinish(postBegan,postSeconds,postTimingAvailable);
    data->report.postprocessTimingAvailable=postTimingAvailable;data->report.postprocessFaceSeconds=postFaceSeconds;data->report.postprocessCellSeconds=postCellSeconds;
    data->report.postprocessFinalizeSeconds=postFinalizeSeconds;data->report.postprocessSeconds=postSeconds;
    data->report.totalSeconds=elapsed(started);return ElectricVofField(std::move(data));
}
const Vof& ElectricVofField::geometry() const{require(bool(data_),"consumed electric field");return data_->geometry;}
const ElectricVofCharge& ElectricVofField::sourceCharge() const{require(bool(data_),"consumed electric field");return data_->source;}
bool ElectricVofField::isElectrostaticSnapshot() const{require(bool(data_),"consumed electric field");return data_->report.electrostaticSnapshot;}
const ElectricVofAdvectedCharge* ElectricVofField::advection() const{require(bool(data_),"consumed electric field");return data_->advection.get();}
const std::vector<double>& ElectricVofField::implicitInputCharge() const{require(bool(data_),"consumed electric field");return data_->implicitInput;}
const ElectricVofCharge& ElectricVofField::chargeCandidate() const{require(bool(data_),"consumed electric field");return data_->candidate;}
Clock ElectricVofField::geometryTime() const{return geometry().timeLayer();}
Clock ElectricVofField::chargeSourceTime() const{return sourceCharge().timeLayer();}
Clock ElectricVofField::candidateTime() const{require(bool(data_),"consumed electric field");return data_->candidateClock;}
double ElectricVofField::effectiveDt() const{require(bool(data_),"consumed electric field");return data_->dt;}
const std::vector<double>& ElectricVofField::potential() const{require(bool(data_),"consumed electric field");return data_->phi;}
const std::vector<std::array<double,3>>& ElectricVofField::potentialTerms() const{require(bool(data_),"consumed electric field");return data_->potentialTerms;}
const std::vector<Vector2D>& ElectricVofField::cellPotentialGradients() const{require(bool(data_),"consumed electric field");
    require(data_->cellGradientRecovery==CellGradientRecovery::Recover,"cell potential gradients were not requested");return data_->potentialGradients;}
const std::vector<double>& ElectricVofField::epsilonByCell() const{require(bool(data_),"consumed electric field");return data_->epsilon;}
const std::vector<double>& ElectricVofField::conductivityByCell() const{require(bool(data_),"consumed electric field");return data_->conductivity;}
const std::vector<ElectricVofFace>& ElectricVofField::faces() const{require(bool(data_),"consumed electric field");return data_->faces;}
const ElectricVofOhmicTransferReceipt& ElectricVofField::ohmicTransferReceipt() const{require(bool(data_),"consumed electric field");return data_->ohmicReceipt;}
const ElectricVofReport& ElectricVofField::report() const{require(bool(data_),"consumed electric field");return data_->report;}
void ElectricVofField::validateFor(const QuadtreeMesh& mesh,const Vof& geometry,const ElectricVofCharge& source,Clock candidate) const
{
    NativeMeshIdentity::ValidationBatch meshValidation(mesh);
    require(bool(data_),"consumed electric field");geometry.validateFor(mesh,data_->geometry.identity(),data_->geometry.timeLayer());
    Detail::ElectricVofValidationBatch batch;batch.rememberValidated(mesh,geometry.domain());
    source.validateForInBatch(mesh,data_->source.timeLayer(),batch);require(source.sameInstance(data_->source),"foreign same-value source charge owner");
    candidate.time();require(sameClock(candidate,data_->candidateClock),"complete electric candidate clock mismatch");
    if(data_->advection)data_->advection->validateForInBatch(mesh,source,data_->advection->transport(),candidate,batch);
    data_->candidate.validateForInBatch(mesh,candidate,batch);
}
}
