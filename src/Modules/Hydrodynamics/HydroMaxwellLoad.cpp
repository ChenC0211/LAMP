#include "Modules/Hydrodynamics/HydroMaxwellLoad.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <cfenv>
#include <chrono>
#include <cerrno>
#include <cstring>
#include <cstdint>
#include <limits>
#include <exception>
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#include <omp.h>
#endif
#if defined(_MSC_VER) && defined(_M_X64)
#include <float.h>
#include <xmmintrin.h>
#endif

namespace LAMP::Hydrodynamics
{
namespace
{
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Matrix=HydroStrainOperator::Matrix;
using Clock=LAMP::TimeLayer;
using MaxwellWallClock=std::chrono::steady_clock;
struct MaxwellTimingState
{
    std::fenv_t native{};unsigned csr=0;int error=errno;bool captured=false;
    MaxwellTimingState() noexcept
    {
#if defined(_MSC_VER) && defined(_M_X64)
        csr=_mm_getcsr();
#endif
        captured=std::fegetenv(&native)==0;errno=error;
    }
    ~MaxwellTimingState() noexcept
    {
        if(captured&&std::fesetenv(&native)!=0)std::terminate();
#if defined(_MSC_VER) && defined(_M_X64)
        _mm_setcsr(csr);
#endif
        errno=error;
    }
};
struct MaxwellStamp {MaxwellWallClock::time_point value{};bool valid=false;};
MaxwellStamp maxwellNow() noexcept
{
    const MaxwellTimingState state;if(!state.captured)return {};
    return {MaxwellWallClock::now(),true};
}
double maxwellSeconds(MaxwellStamp a,MaxwellStamp b,bool& valid) noexcept
{
    const MaxwellTimingState state;
    valid=a.valid&&b.valid&&state.captured;if(!valid)return 0;
    return std::chrono::duration<double>(b.value-a.value).count();
}
bool maxwellNormalOrZero(double value) noexcept
{
    std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof(bits));bits&=UINT64_C(0x7fffffffffffffff);
    return bits==0||(bits>=UINT64_C(0x0010000000000000)&&bits<UINT64_C(0x7ff0000000000000));
}
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
// UCRT's FE_ALL_EXCEPT omits denormal. Retain its entire native state and
// raw MXCSR, rather than synthesizing flags with feraiseexcept.
struct MaxwellFloatingState
{
    std::fenv_t native{};unsigned csr=0;bool valid=false;
    bool capture() noexcept
    {csr=_mm_getcsr();valid=std::fegetenv(&native)==0;return valid;}
    void restore() const noexcept
    {if(valid&&std::fesetenv(&native)!=0)std::terminate();_mm_setcsr(csr);}
    bool enterFresh(const MaxwellFloatingState& source) const noexcept
    {
        auto clean=source.native;clean._Fe_stat=0;
        if(std::fesetenv(&clean)!=0)return false;
        _mm_setcsr(source.csr&~0x3fu);return true;
    }
};
bool maxwellParallelQualified(MaxwellFloatingState& caller) noexcept
{
    if(omp_in_parallel()||omp_get_max_threads()<=1||!caller.capture())return false;
    unsigned int control=0;
    return _controlfp_s(&control,0,0)==0&&(control&_MCW_EM)==_MCW_EM
        &&(caller.csr&0x1f80u)==0x1f80u;
}
#endif
void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[HydroMaxwellLoad] ")+why);}
bool sameClock(Clock a,Clock b)
{return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
void qualifiedBase(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,const HydroPhaseGeometry& phase,Clock clock)
{
    clock.time();require(clock.position==Clock::Position::Midpoint,"Maxwell force requires a complete Midpoint clock");
    require(phase.model()==HydroPhaseGeometry::Model::GeometricVof,"Maxwell force requires genuine VOF phase geometry");
    require(sameClock(phase.timeLayer(),clock),"phase/force complete clock mismatch");
    phase.validateVofFor(mesh,phase.vofIdentity(),clock);pressure.layout().validateMesh(mesh);
    require(pressure.layout().sameIdentity(phase.layout()),"pressure/phase actual native layout mismatch");
    require(&pressure.divergence().matrix()==&phase.metrics().divergence().matrix(),"pressure must retain the actual phase divergence owner");
}
void qualified(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,const HydroPhaseGeometry& phase,
    const Electric::ElectricVofMaxwellStress& stress,Clock clock)
{
    qualifiedBase(mesh,pressure,phase,clock);stress.validateFor(mesh,stress.field(),clock);
    require(stress.field().geometry().identity()==phase.vofPublication().identity(),"Maxwell stress has a different actual VOF owner");
    require(stress.axisymmetric()==phase.layout().domain().axisymmetric(),"stress/phase metric mismatch");
}
void qualifiedSharp(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,const HydroPhaseGeometry& phase,
    const Electric::ElectricVofSharpMaxwellStress& stress,Clock clock)
{
    qualifiedBase(mesh,pressure,phase,clock);const auto& field=stress.field();const auto& smooth=field.geometry();
    stress.validateFor(mesh,field);smooth.validateFor(mesh,phase.vofPublication());
    require(smooth.sourceIdentity()==phase.vofIdentity()&&sameClock(smooth.sourceTime(),clock)
        &&sameClock(field.candidateTime(),clock),"sharp stress lacks the actual coeval Midpoint geometry/field");
    require(field.advection()!=nullptr&&field.effectiveDt()==clock.time()-clock.begin,
        "sharp Midpoint force requires the genuine stage-advection Ohmic receipt");
    require(field.sourceCharge().timeLayer().position!=Clock::Position::Midpoint
        &&field.sourceCharge().timeLayer().time()==clock.begin,"sharp force lost its accepted interval-begin charge");
    field.implicitInputCharge().validateFor(mesh,smooth,clock);field.chargeCandidate().validateFor(mesh,smooth,clock);
    require(stress.axisymmetric()==phase.layout().domain().axisymmetric(),"sharp stress/phase metric mismatch");
}
template<class Tensor>std::array<double,4> packed(const Tensor& t,bool rz)
{
    const auto shear=FV::Dyadic::product({std::sqrt(2.),t.xy}).roundedDouble();
    return {t.xx,t.yy,rz?t.hoop:shear,shear};
}
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
using MaxwellProducts=FV::Binary64ExactProducts;
constexpr std::size_t maximumCompiledTerms=8000000;
constexpr std::size_t maximumCompiledBytes=256u*1024u*1024u;
bool maxwellCompiledControls() noexcept
{
#if !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    const int saved=errno;unsigned control=0;const bool read=_controlfp_s(&control,0,0)==0;errno=saved;
    const unsigned csr=_mm_getcsr();
    if(!read||std::fegetround()!=FE_TONEAREST||(control&_MCW_RC)!=_RC_NEAR
        ||(control&_MCW_EM)!=_MCW_EM||(csr&0x1f80u)!=0x1f80u||(csr&0xe040u)!=0)return false;
#ifdef _EM_AMBIGUOUS
    if(control&_EM_AMBIGUOUS)return false;
#endif
    return sizeof(double)==8&&std::numeric_limits<double>::is_iec559&&std::numeric_limits<double>::digits==53;
#else
    return false;
#endif
}
struct MaxwellCompiledColumns
{
    std::vector<std::size_t> outer;
    std::vector<std::uint32_t> rows;
    std::vector<MaxwellProducts::Factor> factors;
    MaxwellProducts::Result sum(int column,const std::vector<MaxwellProducts::Value>& values) const
    {
        const auto first=outer[static_cast<std::size_t>(column)],count=outer[static_cast<std::size_t>(column)+1]-first;
        return MaxwellProducts::sum(count?factors.data()+first:nullptr,count?rows.data()+first:nullptr,
            count,values.data(),values.size());
    }
};
template<class Weight>
MaxwellCompiledColumns compileMaxwellColumns(const Matrix& source,const std::vector<unsigned char>& active,Weight weight)
{
    require(source.isCompressed()&&source.rows()>=0&&source.cols()>=0
        &&static_cast<unsigned long long>(source.rows())<=std::numeric_limits<std::uint32_t>::max(),
        "compiled Maxwell row range unsupported");
    MaxwellCompiledColumns result;const auto n=static_cast<std::size_t>(source.cols());result.outer.resize(n+1,0);
    for(Eigen::Index row=0;row<source.outerSize();++row){
        if(!active.empty()&&!active[static_cast<std::size_t>(row)])continue;
        for(Matrix::InnerIterator e(source,row);e;++e){
            require(e.col()>=0&&e.col()<source.cols(),"compiled Maxwell column out of range");
            ++result.outer[static_cast<std::size_t>(e.col())+1];}}
    for(std::size_t u=0;u<n;++u){
        require(result.outer[u+1]<=maximumCompiledTerms-result.outer[u],"compiled Maxwell entry cap exceeded");
        result.outer[u+1]+=result.outer[u];}
    result.rows.resize(result.outer.back());result.factors.resize(result.outer.back());auto next=result.outer;
    for(Eigen::Index row=0;row<source.outerSize();++row){
        if(!active.empty()&&!active[static_cast<std::size_t>(row)])continue;
        const double actualWeight=weight(row);require(maxwellNormalOrZero(actualWeight),"compiled Maxwell requires normal geometry weight");
        for(Matrix::InnerIterator e(source,row);e;++e){
            require(maxwellNormalOrZero(e.value()),"compiled Maxwell requires normal geometry coefficient");
            const auto at=next[static_cast<std::size_t>(e.col())]++;
            result.rows[at]=static_cast<std::uint32_t>(row);
            result.factors[at]=MaxwellProducts::compile(actualWeight,e.value());}}
    return result;
}
struct MaxwellAssemblyGeometry
{
    struct Side {std::size_t point;int negativeNative,positiveNative;bool vertical;};
    struct Exterior {std::size_t point;int nativeCell;double sign;bool vertical;};
    // Pure geometry only. The pairing is part of the key: A raw/side Green
    // columns must never alias B1's partitioned cell-mean strain columns.
    HydroStrainOperator geometry;
    bool green=false;
    MaxwellCompiledColumns volumeColumns,negativeColumns,positiveColumns,boundaryColumns;
    std::vector<int> volumeCells;
    std::vector<Side> sides;
    std::vector<Exterior> exterior;
    std::size_t axisPoints=0,cellSidePoints=0,retainedBytes=0;
    MaxwellAssemblyGeometry(const HydroStrainOperator& source,bool useGreen):geometry(source),green(useGreen)
    {
        const auto& layout=geometry.layout();const auto& points=geometry.points();
        const auto& e=green?geometry.stressVolumeVelocityMatrix():geometry.matrix();
        const auto& b=geometry.stressFaceVelocityMatrix();const auto& faces=geometry.stressFacePoints();
        require(e.nonZeros()>=0&&static_cast<unsigned long long>(e.nonZeros())<=maximumCompiledTerms
            &&points.size()<=maximumCompiledBytes/sizeof(int),"compiled Maxwell geometry cap exceeded");
        std::vector<unsigned char> activeBoundary(static_cast<std::size_t>(b.rows()),0),
            activeNegative(activeBoundary.size(),0),activePositive(activeBoundary.size(),0);
        volumeCells.reserve(points.size());for(const auto& point:points)volumeCells.push_back(layout.nativeCell(point.cellDof));
        for(std::size_t q=0;q<faces.size();++q){const auto& face=faces[q];
            const bool internal=face.negativeCellDof>=0&&face.positiveCellDof>=0;
            if(face.physicalWeight==0){if(!internal)++axisPoints;continue;}
            const bool vertical=layout.normalAxis(face.velocityDof)==HydroDofLayout::NormalAxis::First;
            if(green){const int negative=face.negativeCellDof<0?-1:layout.nativeCell(face.negativeCellDof),
                    positive=face.positiveCellDof<0?-1:layout.nativeCell(face.positiveCellDof);
                require(negative>=0||positive>=0,"Maxwell side quadrature has no actual owner");
                sides.push_back({q,negative,positive,vertical});
                if(negative>=0){activeNegative[q*2]=activeNegative[q*2+1]=1;++cellSidePoints;}
                if(positive>=0){activePositive[q*2]=activePositive[q*2+1]=1;++cellSidePoints;}}
            if(internal)continue;
            const int p=face.negativeCellDof>=0?face.negativeCellDof:face.positiveCellDof;
            require(p>=0,"exterior Maxwell quadrature has no actual owner");
            exterior.push_back({q,layout.nativeCell(p),face.negativeCellDof>=0?1.:-1.,vertical});
            activeBoundary[q*2]=activeBoundary[q*2+1]=1;}
        std::size_t terms=static_cast<std::size_t>(e.nonZeros());
        const auto count=[&](const std::vector<unsigned char>& active){
            for(Eigen::Index row=0;row<b.outerSize();++row)if(active[static_cast<std::size_t>(row)])
                for(Matrix::InnerIterator entry(b,row);entry;++entry){
                    require(terms<maximumCompiledTerms,"compiled Maxwell entry cap exceeded");++terms;}};
        count(activeBoundary);if(green){count(activeNegative);count(activePositive);}
        std::size_t estimated=sizeof(*this);
        const auto account=[&](std::size_t count,std::size_t width){
            require(width==0||count<=(maximumCompiledBytes-estimated)/width,"compiled Maxwell byte cap exceeded");estimated+=count*width;};
        account(terms,sizeof(std::uint32_t)+sizeof(MaxwellProducts::Factor));
        account((green?4:2)*(static_cast<std::size_t>(layout.velocityCount())+1),sizeof(std::size_t));
        account(volumeCells.capacity(),sizeof(int));account(exterior.capacity(),sizeof(Exterior));account(sides.capacity(),sizeof(Side));
        const int components=green?2:geometry.components();
        volumeColumns=compileMaxwellColumns(e,{},[&](Eigen::Index row){return (green?1.:-1.)*points[static_cast<std::size_t>(row/components)].physicalWeight;});
        boundaryColumns=compileMaxwellColumns(b,activeBoundary,[&](Eigen::Index row){return faces[static_cast<std::size_t>(row/2)].physicalWeight;});
        if(green){
            negativeColumns=compileMaxwellColumns(b,activeNegative,[&](Eigen::Index row){return -faces[static_cast<std::size_t>(row/2)].physicalWeight;});
            positiveColumns=compileMaxwellColumns(b,activePositive,[&](Eigen::Index row){return faces[static_cast<std::size_t>(row/2)].physicalWeight;});}
        require(volumeColumns.rows.size()+negativeColumns.rows.size()+positiveColumns.rows.size()+boundaryColumns.rows.size()<=maximumCompiledTerms,
            "compiled Maxwell entry cap exceeded");
        estimated=sizeof(*this);
        for(const auto* columns:{&volumeColumns,&negativeColumns,&positiveColumns,&boundaryColumns}){
            account(columns->outer.capacity(),sizeof(std::size_t));account(columns->rows.capacity(),sizeof(std::uint32_t));
            account(columns->factors.capacity(),sizeof(MaxwellProducts::Factor));}
        account(volumeCells.capacity(),sizeof(int));account(exterior.capacity(),sizeof(Exterior));account(sides.capacity(),sizeof(Side));retainedBytes=estimated;
    }
};
std::shared_ptr<const MaxwellAssemblyGeometry> maxwellAssemblyGeometry(const HydroStrainOperator& geometry,bool green,bool& hit)
{
    static std::mutex mutex;static std::shared_ptr<const MaxwellAssemblyGeometry> last;
    std::lock_guard<std::mutex> lock(mutex);hit=false;
    if(last&&last->green==green
        &&(green?last->geometry.isStressCompatible()==geometry.isStressCompatible():&last->geometry.matrix()==&geometry.matrix())
        &&&last->geometry.stressVolumeVelocityMatrix()==&geometry.stressVolumeVelocityMatrix()
        &&&last->geometry.stressFaceVelocityMatrix()==&geometry.stressFaceVelocityMatrix()){hit=true;return last;}
    const auto prepared=std::make_shared<MaxwellAssemblyGeometry>(geometry,green);
    last=prepared;return prepared; // Publish complete, cap-checked entries only.
}
#endif
struct SharpQuadrature
{
    std::vector<HydroStrainPoint> volume;
    std::vector<HydroStressFacePoint> faces;
    std::vector<unsigned> volumePhases,facePhases;
};
SharpQuadrature sharpQuadrature(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,
    const Interface::VofSmoothGeometry& smooth)
{
    SharpQuadrature result;const auto& layout=phase.layout();const bool rz=layout.domain().axisymmetric();
    for(const auto& cell:smooth.cells()){
        const int p=layout.pressureDof(cell.nativeCell);require(p>=0,"smooth volume cell is outside the actual Hydro layout");
        for(unsigned s=0;s<2;++s)for(const auto& q:cell.phaseQuadrature[s]){
            const double weight=static_cast<double>(q.weight);const Vector2D point{static_cast<double>(q.point.x),static_cast<double>(q.point.y)};
            require(std::isfinite(weight)&&weight>0,"unrepresentable actual smooth phase volume weight");
            result.volume.push_back({p,point,weight});result.volumePhases.push_back(s);}}
    std::map<int,Geometry::PlanarCellCut::Point> crossings;
    for(const auto& node:smooth.nodes())require(crossings.emplace(node.nativeRow,node.point).second,"multiple smooth crossings of one native face");
    constexpr std::array<long double,4> nodes{{-.861136311594052575223946488893L,-.339981043584856264802665759103L,
        .339981043584856264802665759103L,.861136311594052575223946488893L}};
    constexpr std::array<long double,4> weights{{.347854845137453857373063949222L,.652145154862546142626936050778L,
        .652145154862546142626936050778L,.347854845137453857373063949222L}};
    for(int u=0;u<layout.velocityCount();++u){const int row=layout.nativeFace(u);const auto f=layout.domain().face(mesh,row);
        const long double low=(f.horizontal?f.center.x:f.center.y)-f.tangentSpacing/2;
        const long double high=low+f.tangentSpacing;std::vector<long double> cuts{low,high};const auto found=crossings.find(row);
        if(found!=crossings.end()){const auto p=found->second;
            require(f.horizontal?p.y==f.center.y:p.x==f.center.x,"smooth crossing is not on its actual native face");
            const auto value=f.horizontal?p.x:p.y;require(value>low&&value<high,"smooth crossing is not transverse inside the actual face");cuts.push_back(value);}
        std::sort(cuts.begin(),cuts.end());
        const int negative=f.negativeCell<0?-1:layout.pressureDof(f.negativeCell),positive=f.positiveCell<0?-1:layout.pressureDof(f.positiveCell);
        require(negative>=0||positive>=0,"native face has no actual Hydro cell");
        for(std::size_t interval=0;interval+1<cuts.size();++interval){const auto center=(cuts[interval]+cuts[interval+1])/2;
            const auto radius=(cuts[interval+1]-cuts[interval])/2;require(radius>0,"empty sharp native-face phase interval");
            for(unsigned g=0;g<4;++g){const auto coordinate=center+radius*nodes[g];Vector2D point=f.center;
                if(f.horizontal)point.x=static_cast<double>(coordinate);else point.y=static_cast<double>(coordinate);
                const Geometry::PlanarCellCut::Point query{point.x,point.y};unsigned s;
                if(f.negativeCell>=0){s=smooth.phaseAt(f.negativeCell,query);
                    if(f.positiveCell>=0)require(s==smooth.phaseAt(f.positiveCell,query),"shared smooth native face has inconsistent phase traces");}
                else s=smooth.phaseAt(f.positiveCell,query);
                const auto physical=radius*weights[g]*(rz?CartesianMeasure::twoPi*point.x:1.L);
                const double weight=static_cast<double>(physical);
                require(std::isfinite(weight)&&(f.area==0?weight==0:weight>0),"unrepresentable sharp native-face physical weight");
                result.faces.push_back({u,row,negative,positive,point,weight});result.facePhases.push_back(s);}}
    }
    return result;
}
}
struct HydroMaxwellLoad::Data
{
    using Stress=std::variant<Electric::ElectricVofMaxwellStress,Electric::ElectricVofSharpMaxwellStress>;
    HydroPressureLoad pressure;HydroPhaseGeometry phase;Stress stress;
    HydroStrainOperator geometry;
    // The sole retained force is one RNE of the complete exact row.
    std::vector<double> rhs;
    HydroMaxwellLoadReport report;
    bool validationTimingAvailable=false;
    template<class S>Data(HydroPressureLoad p,HydroPhaseGeometry g,S s,HydroStrainOperator e):
        pressure(std::move(p)),phase(std::move(g)),stress(std::move(s)),geometry(std::move(e)){}
};
HydroMaxwellLoad::HydroMaxwellLoad(std::shared_ptr<const Data> data):data_(std::move(data)){}
HydroMaxwellLoad HydroMaxwellLoad::prepare(const QuadtreeMesh& mesh,HydroPressureLoad pressure,
    HydroPhaseGeometry phase,Electric::ElectricVofMaxwellStress stress)
{
    const auto validationStarted=maxwellNow();
    qualified(mesh,pressure,phase,stress,phase.timeLayer());
    const auto validationFinished=maxwellNow();
    auto geometry=HydroStrainOperator::buildStressCompatible(mesh,phase.metrics());
    auto data=std::make_shared<Data>(std::move(pressure),std::move(phase),std::move(stress),std::move(geometry));
    data->report.validationSeconds=maxwellSeconds(validationStarted,validationFinished,data->validationTimingAvailable);
    return assemble(std::move(data),{},{});
}
HydroMaxwellLoad HydroMaxwellLoad::prepareRt0Green(const QuadtreeMesh& mesh,HydroPressureLoad pressure,
    HydroPhaseGeometry phase,Electric::ElectricVofMaxwellStress stress)
{
    const auto validationStarted=maxwellNow();
    qualified(mesh,pressure,phase,stress,phase.timeLayer());
    const auto validationFinished=maxwellNow();
    auto geometry=HydroStrainOperator::buildRt0Green(mesh,phase.metrics());
    auto data=std::make_shared<Data>(std::move(pressure),std::move(phase),std::move(stress),std::move(geometry));
    data->report.validationSeconds=maxwellSeconds(validationStarted,validationFinished,data->validationTimingAvailable);
    return assemble(std::move(data),{},{});
}
HydroMaxwellLoad HydroMaxwellLoad::prepare(const QuadtreeMesh& mesh,HydroPressureLoad pressure,
    HydroPhaseGeometry phase,Electric::ElectricVofMaxwellStress stress,HydroStrainOperator geometry)
{
    const auto validationStarted=maxwellNow();
    qualified(mesh,pressure,phase,stress,phase.timeLayer());
    require(geometry.isRt0GreenCompatible(),"prepared Maxwell geometry requires actual RT0 Green capability");
    geometry.layout().validateMesh(mesh);
    require(geometry.layout().sameIdentity(phase.layout()),"prepared Maxwell geometry has a different actual native layout");
    require(&geometry.metrics().divergence().matrix()==&phase.metrics().divergence().matrix(),
        "prepared Maxwell geometry must retain the actual phase divergence owner");
    const auto validationFinished=maxwellNow();
    auto data=std::make_shared<Data>(std::move(pressure),std::move(phase),std::move(stress),std::move(geometry));
    data->report.validationSeconds=maxwellSeconds(validationStarted,validationFinished,data->validationTimingAvailable);
    return assemble(std::move(data),{},{});
}
HydroMaxwellLoad HydroMaxwellLoad::prepareSharp(const QuadtreeMesh& mesh,HydroPressureLoad pressure,
    HydroPhaseGeometry phase,Electric::ElectricVofSharpMaxwellStress stress)
{
    const auto validationStarted=maxwellNow();
    qualifiedSharp(mesh,pressure,phase,stress,phase.timeLayer());
    const auto validationFinished=maxwellNow();
    const auto quadrature=sharpQuadrature(mesh,phase,stress.field().geometry());
    auto geometry=HydroStrainOperator::buildStressCompatible(mesh,phase.metrics(),quadrature.volume,quadrature.faces);
    auto data=std::make_shared<Data>(std::move(pressure),std::move(phase),std::move(stress),std::move(geometry));
    data->report.validationSeconds=maxwellSeconds(validationStarted,validationFinished,data->validationTimingAvailable);
    return assemble(std::move(data),quadrature.volumePhases,quadrature.facePhases);
}
HydroMaxwellLoad HydroMaxwellLoad::assemble(std::shared_ptr<Data> data,
    const std::vector<unsigned>& volumePhases,const std::vector<unsigned>& facePhases)
{
    const auto& layout=data->geometry.layout();const bool rz=layout.domain().axisymmetric();
    const auto& points=data->geometry.points();const auto& faces=data->geometry.stressFacePoints();
    const auto* sharp=std::get_if<Electric::ElectricVofSharpMaxwellStress>(&data->stress);
    const auto* mixed=std::get_if<Electric::ElectricVofMaxwellStress>(&data->stress);
    const bool green=mixed!=nullptr;
    const auto& e=green?data->geometry.stressVolumeVelocityMatrix():data->geometry.matrix();
    const auto& b=data->geometry.stressFaceVelocityMatrix();
    const int components=green?2:data->geometry.components(),nu=layout.velocityCount();
    require(!sharp||(volumePhases.size()==points.size()&&facePhases.size()==faces.size()),"sharp quadrature phase count mismatch");
    require(!green||(e.rows()==static_cast<Eigen::Index>(points.size())*2&&e.cols()==nu),
        "RT0 Green pairing requires the actual raw volume velocity dimensions");
    std::vector<FV::Dyadic> force(nu);
    std::vector<std::array<double,4>> tensor(points.size());
    std::vector<std::array<double,2>> negativeTraction(faces.size()),positiveTraction(faces.size()),traction(faces.size());
    data->report.volumePoints=points.size();
    data->report.pairing=green?HydroMaxwellPairing::Rt0Green:HydroMaxwellPairing::CellMeanStrain;
    const auto volumeSample=[&](std::size_t q,int cell){
        if(green){const auto sample=mixed->sample(cell,points[q].point);
            return std::array<double,4>{sample.tensorDivergence.x,sample.tensorDivergence.y,0,0};}
        return packed(sharp->sample(cell,volumePhases[q],points[q].point).tensor,rz);
    };
    const auto sideSample=[&](int cell,std::size_t q,bool vertical){const auto t=mixed->sample(cell,faces[q].point).tensor;
        return std::array<double,2>{vertical?t.xx:t.xy,vertical?t.xy:t.yy};};
    // A uses the represented tensor Green identity, not a cell-center qE
    // force replacement. Keep both actual side tractions independently exact;
    // do not round their jump before native contraction. B1 is unchanged.
    const auto serialAssembly=[&](){
        for(std::size_t q=0;q<points.size();++q){const auto& point=points[q];
            tensor[q]=volumeSample(q,layout.nativeCell(point.cellDof));
            for(int k=0;k<components;++k){const auto weighted=FV::Dyadic::product({(green?1.:-1.)*point.physicalWeight,tensor[q][k]});
                for(Matrix::InnerIterator entry(e,static_cast<Eigen::Index>(q)*components+k);entry;++entry){
                    force[entry.col()].add(weighted.multiplied(entry.value()));}}
        }
        for(std::size_t q=0;q<faces.size();++q){const auto& face=faces[q];
            const bool internal=face.negativeCellDof>=0&&face.positiveCellDof>=0;
            if(face.physicalWeight==0){if(!internal)++data->report.axisPoints;continue;}
            const bool vertical=layout.normalAxis(face.velocityDof)==HydroDofLayout::NormalAxis::First;
            if(green){
                if(face.negativeCellDof>=0)negativeTraction[q]=sideSample(layout.nativeCell(face.negativeCellDof),q,vertical);
                if(face.positiveCellDof>=0)positiveTraction[q]=sideSample(layout.nativeCell(face.positiveCellDof),q,vertical);
                for(int side=0;side<2;++side){if((side==0?face.negativeCellDof:face.positiveCellDof)<0)continue;
                    ++data->report.cellSidePoints;const auto& t=side==0?negativeTraction[q]:positiveTraction[q];
                    for(int axis=0;axis<2;++axis){const auto weighted=FV::Dyadic::product({(side==0?-1.:1.)*face.physicalWeight,t[axis]});
                        for(Matrix::InnerIterator entry(b,static_cast<Eigen::Index>(q)*2+axis);entry;++entry){
                            force[entry.col()].add(weighted.multiplied(entry.value()));}}}
            }
            if(internal)continue;
            ++data->report.boundaryPoints;const int p=face.negativeCellDof>=0?face.negativeCellDof:face.positiveCellDof;
            require(p>=0,"exterior Maxwell quadrature has no actual owner");const double sign=face.negativeCellDof>=0?1.:-1.;
            if(green){const auto& t=face.negativeCellDof>=0?negativeTraction[q]:positiveTraction[q];traction[q]={sign*t[0],sign*t[1]};}
            else {const auto t=sharp->sample(layout.nativeCell(p),facePhases[q],face.point).tensor;
                traction[q]={sign*(vertical?t.xx:t.xy),sign*(vertical?t.xy:t.yy)};}
            for(int axis=0;axis<2;++axis){const auto weighted=FV::Dyadic::product({face.physicalWeight,traction[q][axis]});
                for(Matrix::InnerIterator entry(b,static_cast<Eigen::Index>(q)*2+axis);entry;++entry){
                    force[entry.col()].add(weighted.multiplied(entry.value()));}}
        }
    };
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    MaxwellFloatingState caller;
    if(maxwellParallelQualified(caller)&&maxwellCompiledControls()){
        try{
            const auto planStarted=maxwellNow();bool planHit=false;
            const auto plan=maxwellAssemblyGeometry(data->geometry,green,planHit);const auto planFinished=maxwellNow();
            const auto& volumeColumns=plan->volumeColumns;const auto& boundaryColumns=plan->boundaryColumns;
            data->report.axisPoints=plan->axisPoints;data->report.boundaryPoints=plan->exterior.size();data->report.cellSidePoints=plan->cellSidePoints;
            std::vector<MaxwellProducts::Value> weightedVolume(e.rows()),weightedBoundary(b.rows()),
                weightedNegative(green?b.rows():0),weightedPositive(green?b.rows():0);
            const auto faceSamples=green?plan->sides.size():plan->exterior.size();
            std::vector<std::exception_ptr> sampleFailures(points.size()+faceSamples),failures(nu);
            MaxwellFloatingState prepared;
            if(!prepared.capture())throw std::runtime_error("Maxwell parallel environment capture failed");
            unsigned long raisedNative=0;unsigned raisedCsr=0;int setupFailure=0;
            const auto numericStarted=maxwellNow();
#pragma omp parallel
            {
                MaxwellFloatingState worker;const bool captured=worker.capture();const bool ready=captured&&worker.enterFresh(prepared);
                if(!ready){
#pragma omp atomic write
                    setupFailure=1;
                }
#pragma omp for schedule(static)
                for(std::ptrdiff_t index=0;index<static_cast<std::ptrdiff_t>(points.size());++index){
                    if(!ready)continue;const auto q=static_cast<std::size_t>(index);
                    try{tensor[q]=volumeSample(q,plan->volumeCells[q]);
                        for(int k=0;k<components;++k){require(maxwellNormalOrZero(tensor[q][k]),"compiled Maxwell requires normal represented tensor/divergence");
                            weightedVolume[static_cast<Eigen::Index>(q)*components+k]=MaxwellProducts::decode(tensor[q][k]);}
                    }catch(...){sampleFailures[q]=std::current_exception();}
                }
#pragma omp for schedule(static)
                for(std::ptrdiff_t index=0;index<static_cast<std::ptrdiff_t>(faceSamples);++index){
                    if(!ready)continue;const auto i=static_cast<std::size_t>(index);
                    try{
                        if(green){const auto& item=plan->sides[i];const auto q=item.point;
                            const bool exterior=item.negativeNative<0||item.positiveNative<0;
                            if(item.negativeNative>=0)negativeTraction[q]=sideSample(item.negativeNative,q,item.vertical);
                            if(item.positiveNative>=0)positiveTraction[q]=sideSample(item.positiveNative,q,item.vertical);
                            for(int axis=0;axis<2;++axis){const auto row=static_cast<Eigen::Index>(q)*2+axis;
                                if(item.negativeNative>=0){require(maxwellNormalOrZero(negativeTraction[q][axis]),"compiled Maxwell requires normal represented side traction");
                                    weightedNegative[row]=MaxwellProducts::decode(negativeTraction[q][axis]);}
                                if(item.positiveNative>=0){require(maxwellNormalOrZero(positiveTraction[q][axis]),"compiled Maxwell requires normal represented side traction");
                                    weightedPositive[row]=MaxwellProducts::decode(positiveTraction[q][axis]);}
                                if(exterior){traction[q][axis]=item.negativeNative>=0?negativeTraction[q][axis]:-positiveTraction[q][axis];
                                    weightedBoundary[row]=MaxwellProducts::decode(traction[q][axis]);}}
                        }else{const auto& item=plan->exterior[i];const auto q=item.point;
                            const auto t=sharp->sample(item.nativeCell,facePhases[q],faces[q].point).tensor;
                            traction[q]={item.sign*(item.vertical?t.xx:t.xy),item.sign*(item.vertical?t.xy:t.yy)};
                            for(int axis=0;axis<2;++axis){const auto row=static_cast<Eigen::Index>(q)*2+axis;
                                require(maxwellNormalOrZero(traction[q][axis]),"compiled Maxwell requires normal represented traction");
                                weightedBoundary[row]=MaxwellProducts::decode(traction[q][axis]);}}
                    }catch(...){sampleFailures[points.size()+i]=std::current_exception();}
                }
#pragma omp for schedule(static) nowait
                for(int u=0;u<nu;++u){if(!ready)continue;
                    try{auto v=volumeColumns.sum(u,weightedVolume),bnd=boundaryColumns.sum(u,weightedBoundary);
                        if(green){const auto negative=plan->negativeColumns.sum(u,weightedNegative),positive=plan->positiveColumns.sum(u,weightedPositive);
                            v.sum.add(negative.sum);v.sum.add(positive.sum);}
                        v.sum.add(bnd.sum);force[u]=std::move(v.sum);
                    }catch(...){failures[u]=std::current_exception();}
                }
                if(ready){std::fenv_t after{};const bool got=std::fegetenv(&after)==0;const unsigned flags=_mm_getcsr()&0x3fu;
                    if(!got){
#pragma omp atomic write
                        setupFailure=1;
                    }
#pragma omp critical(lamp_maxwell_floating_flags)
                    {if(got)raisedNative|=after._Fe_stat;raisedCsr|=flags;}}
                if(captured)worker.restore();
            }
            auto merged=prepared.native;merged._Fe_stat|=raisedNative;
            if(std::fesetenv(&merged)!=0)std::terminate();_mm_setcsr(prepared.csr|raisedCsr);
            if(setupFailure)throw std::runtime_error("Maxwell parallel worker environment setup failed");
            for(const auto& failure:sampleFailures)if(failure)std::rethrow_exception(failure);
            for(const auto& failure:failures)if(failure)std::rethrow_exception(failure);
            const auto numericFinished=maxwellNow();bool planTiming=false,numericTiming=false;
            data->report.compiledPlanSeconds=maxwellSeconds(planStarted,planFinished,planTiming);
            data->report.samplingAndContractionSeconds=maxwellSeconds(numericStarted,numericFinished,numericTiming);
            data->report.compiledContraction=true;data->report.compiledGeometryHit=planHit;
            data->report.compiledVolumeTerms=volumeColumns.rows.size()+plan->negativeColumns.rows.size()+plan->positiveColumns.rows.size();
            data->report.compiledBoundaryTerms=boundaryColumns.rows.size();data->report.compiledRetainedBytes=plan->retainedBytes;
            data->report.compiledTimingAvailable=data->validationTimingAvailable&&planTiming&&numericTiming;
        }catch(...){
            // Same literal fallback and original caller environment contract.
            for(int u=0;u<nu;++u)force[u]={};
            data->report.axisPoints=0;data->report.boundaryPoints=0;data->report.cellSidePoints=0;
            data->report.compiledContraction=false;data->report.compiledGeometryHit=false;data->report.compiledTimingAvailable=false;
            data->report.compiledVolumeTerms=0;data->report.compiledBoundaryTerms=0;data->report.compiledRetainedBytes=0;
            data->report.compiledPlanSeconds=0;data->report.samplingAndContractionSeconds=0;caller.restore();serialAssembly();
        }
    }else serialAssembly();
#else
    serialAssembly();
#endif
    const auto rneStarted=data->report.compiledContraction?maxwellNow():MaxwellStamp{};
    data->rhs.resize(nu);
    for(int u=0;u<nu;++u)data->rhs[u]=force[u].roundedDouble();
    if(data->report.compiledContraction){bool finalTiming=false;
        data->report.finalRneSeconds=maxwellSeconds(rneStarted,maxwellNow(),finalTiming);
        data->report.compiledTimingAvailable=data->report.compiledTimingAvailable&&finalTiming;}
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    // Retired exact stock is temporary: release it under its CURRENT state.
    const int finalErrno=errno;MaxwellFloatingState finalState;const bool finalCaptured=finalState.capture();
    std::vector<FV::Dyadic>().swap(force);
    if(finalCaptured)finalState.restore();errno=finalErrno;
#endif
    return HydroMaxwellLoad(std::move(data));
}
const HydroPressureLoad& HydroMaxwellLoad::pressureLoad() const{require(bool(data_),"consumed load owner");return data_->pressure;}
const HydroPhaseGeometry& HydroMaxwellLoad::phaseGeometry() const{require(bool(data_),"consumed load owner");return data_->phase;}
const Electric::ElectricVofMaxwellStress& HydroMaxwellLoad::stress() const
{require(bool(data_)&&std::holds_alternative<Electric::ElectricVofMaxwellStress>(data_->stress),"load does not own an A Maxwell stress");return std::get<Electric::ElectricVofMaxwellStress>(data_->stress);}
bool HydroMaxwellLoad::isSharp() const noexcept{return data_&&std::holds_alternative<Electric::ElectricVofSharpMaxwellStress>(data_->stress);}
const Electric::ElectricVofSharpMaxwellStress& HydroMaxwellLoad::sharpStress() const
{require(isSharp(),"load does not own a B1 sharp Maxwell stress");return std::get<Electric::ElectricVofSharpMaxwellStress>(data_->stress);}
const HydroStrainOperator& HydroMaxwellLoad::geometry() const{require(bool(data_),"consumed load owner");return data_->geometry;}
const std::vector<double>& HydroMaxwellLoad::momentumRhsLoad() const{require(bool(data_),"consumed load owner");return data_->rhs;}
const HydroMaxwellLoadReport& HydroMaxwellLoad::report() const{require(bool(data_),"consumed load owner");return data_->report;}
void HydroMaxwellLoad::validateFor(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,const HydroPhaseGeometry& phase,
    const Electric::ElectricVofMaxwellStress& stress,Clock clock) const
{
    require(bool(data_)&&!isSharp(),"load does not own an A Maxwell stress");qualified(mesh,pressure,phase,stress,clock);
    require(&pressure.matrix()==&data_->pressure.matrix(),"different actual pressure load owner");
    require(&phase.cellVolumes()==&data_->phase.cellVolumes(),"different actual phase geometry owner");
    require(&stress.cells()==&this->stress().cells(),"different actual Maxwell stress owner");
    require(sameClock(data_->phase.timeLayer(),clock),"different complete load clock");
}
void HydroMaxwellLoad::validateSharpFor(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,const HydroPhaseGeometry& phase,
    const Electric::ElectricVofSharpMaxwellStress& stress,Clock clock) const
{
    require(isSharp(),"load does not own a B1 sharp Maxwell stress");qualifiedSharp(mesh,pressure,phase,stress,clock);
    require(&pressure.matrix()==&data_->pressure.matrix(),"different actual pressure load owner");
    require(&phase.cellVolumes()==&data_->phase.cellVolumes(),"different actual phase geometry owner");
    require(stress.sameInstance(sharpStress()),"different actual sharp Maxwell stress owner");
    require(sameClock(data_->phase.timeLayer(),clock),"different complete load clock");
}
}
