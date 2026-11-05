#include "Modules/Electric/ElectricVofChargeTransport.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <atomic>
#include <cfenv>
#include <cerrno>
#include <cstring>
#include <memory>
#include <exception>
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#include <omp.h>
#include <float.h>
#include <xmmintrin.h>
#endif

namespace LAMP::Electric
{
namespace
{
namespace FV=Geometry::Detail::FiniteVolume;
using Clock=LAMP::TimeLayer;
using Vof=Interface::GeometricVofPublication;
using Transport=Interface::GeometricVofTransport;
using Transfer=Interface::GeometricVofMeshTransfer;
using Timer=std::chrono::steady_clock;
void require(bool ok,const char* reason)
{if(!ok)throw std::invalid_argument(std::string("[ElectricVofChargeTransport] ")+reason);}
bool sameClock(Clock a,Clock b)
{return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
double elapsed(Timer::time_point started){return std::chrono::duration<double>(Timer::now()-started).count();}
void controlsValid(const ElectricVofControls& c)
{
    require(std::isfinite(c.relativeLinearTolerance)&&c.relativeLinearTolerance>0&&c.relativeLinearTolerance<=1e-11&&
        std::isfinite(c.relativeConservationTolerance)&&c.relativeConservationTolerance>0&&c.relativeConservationTolerance<=1e-10,
        "cannot relax original electric linear/conservation gates");
    require(c.chargeAdvection==ElectricVofChargeAdvectionScheme::DonorCell||c.chargeAdvection==ElectricVofChargeAdvectionScheme::LimitedLinearFct,
        "unsupported signed charge advection scheme");
}
double represented(const FV::Dyadic& exact)
{const double value=FV::roundDouble(exact).value;require(std::isfinite(value),"unrepresentable final signed charge sum");return value;}
struct Ledger
{
    FV::Dyadic residual,scale,budget;
    void term(double q,int sign,double tolerance)
    {residual.add(q,sign);scale.add(std::abs(q));budget.add(FV::Dyadic::product({tolerance,std::abs(q)}));}
    void activity(double q,double tolerance)
    {scale.add(std::abs(q));budget.add(FV::Dyadic::product({tolerance,std::abs(q)}));}
    void validate() const
    {require(residual.absolute().compare(budget)<=0,"physical signed charge ledger exceeds original bound");}
    double relative() const
    {if(scale.compare(FV::Dyadic{})==0){require(residual.compare(FV::Dyadic{})==0,"nonzero zero-scale charge ledger");return 0;}
        auto s=scale,r=residual.absolute();const FV::Dyadic maximum(std::numeric_limits<double>::max());
        while(s.compare(maximum)>0){s=s.multiplied(.5L);r=r.multiplied(.5L);}
        const double value=represented(s);require(value>0,"unrepresentable charge ledger scale");return represented(r)/value;}
};
// Defined below as the private Electric synchronous validation helper.
struct Stage
{
    std::vector<double> q;
    std::vector<ElectricVofChargeFace> faces;
    ElectricVofChargeAdvectionReport report;
    std::unique_ptr<Ledger> boundaryModelExtra;
};
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
struct ChargeAdvanceFloatingState
{
    std::fenv_t native{};unsigned csr=0;bool valid=false;
    bool capture() noexcept {csr=_mm_getcsr();valid=std::fegetenv(&native)==0;return valid;}
    void restore() const noexcept
    {if(valid&&std::fesetenv(&native)!=0)std::terminate();_mm_setcsr(csr);}
    bool enterFresh(const ChargeAdvanceFloatingState& from) const noexcept
    {auto clean=from.native;clean._Fe_stat=0;if(std::fesetenv(&clean)!=0)return false;
        _mm_setcsr(from.csr&~0x3fu);return true;}
};
bool chargeAdvanceSameControl(const ChargeAdvanceFloatingState& a,const ChargeAdvanceFloatingState& b) noexcept
{
    auto x=a.native,y=b.native;x._Fe_stat=0;y._Fe_stat=0;
    return std::memcmp(&x,&y,sizeof x)==0&&(a.csr&~0x3fu)==(b.csr&~0x3fu);
}
bool chargeAdvanceQualified(ChargeAdvanceFloatingState& caller) noexcept
{
#if defined(_M_FP_STRICT) && _M_FP_STRICT && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    if(!caller.capture())return false;
    if(omp_in_parallel()||omp_get_max_threads()<2)return false;
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
// Declare before optional storage and arm immediately before private Stage
// publication. Cleanup cannot leak allocator FP/errno into the caller.
struct ChargeAdvanceCleanup
{
    ChargeAdvanceFloatingState finalState;int error=0;bool armed=false;
    bool arm() noexcept
    {error=errno;armed=finalState.capture();errno=error;return armed;}
    ~ChargeAdvanceCleanup() noexcept
    {if(armed){finalState.restore();errno=error;}}
};
struct ChargeAdvanceBudget
{
    static constexpr std::size_t limit=64u*1024u*1024u;
    std::size_t used=0;
    void add(std::size_t count,std::size_t width)
    {if(used>limit||count>(limit-used)/width)throw std::length_error("optional charge advance scratch cap");
        used+=count*width;}
};
template<class Kernel>
bool chargeAdvanceRows(std::size_t count,int team,Kernel&& kernel,unsigned& actualWorkers)
{
    ChargeAdvanceFloatingState prepared;const int parentErrno=errno;
    if(!prepared.capture())return false;
    std::atomic<int> failed{0};unsigned long raisedNative=0;unsigned raisedCsr=0;
#pragma omp parallel num_threads(team)
    {
        ChargeAdvanceFloatingState worker;const int workerErrno=errno;
        const bool captured=worker.capture(),ready=captured&&worker.enterFresh(prepared);
        if(!ready)failed.store(1,std::memory_order_relaxed);
        if(ready)errno=parentErrno;
#pragma omp single
        {actualWorkers=std::max(actualWorkers,static_cast<unsigned>(omp_get_num_threads()));}
#pragma omp for schedule(static) nowait
        for(std::ptrdiff_t row=0;row<static_cast<std::ptrdiff_t>(count);++row){
            if(!ready)continue;
            try{kernel(static_cast<std::size_t>(row));}
            catch(...){failed.store(1,std::memory_order_relaxed);}
        }
        if(ready){
            const int numericalErrno=errno;ChargeAdvanceFloatingState after;
            const bool got=after.capture();
            if(!got||numericalErrno!=parentErrno||!chargeAdvanceSameControl(after,prepared))
                failed.store(1,std::memory_order_relaxed);
#pragma omp critical(lamp_charge_advance_floating_flags)
            {if(got){raisedNative|=after.native._Fe_stat;raisedCsr|=after.csr&0x3fu;}}
        }
        if(captured)worker.restore();errno=workerErrno;
    }
    if(failed.load(std::memory_order_relaxed))return false;
    auto merged=prepared.native;merged._Fe_stat|=raisedNative;
    if(std::fesetenv(&merged)!=0)std::terminate();
    _mm_setcsr(prepared.csr|raisedCsr);errno=parentErrno;return true;
}
#endif
void requireInactiveFaceZero(const Interface::GeometricVofTransportFace& f)
{
    require(f.half.total==0&&f.half.positiveTotal==0&&f.half.negativeTotal==0&&
        f.full.total==0&&f.full.positiveTotal==0&&f.full.negativeTotal==0,
        "both-excluded charge face cannot carry accepted net or gross volume");
}
Stage advanceLiteral(const QuadtreeMesh& mesh,const ElectricVofCharge& source,const Transport& t,bool half,const ElectricVofControls& controls,
    const ElectricVofChargeBoundary* boundary=nullptr)
{
    const auto started=Timer::now();const auto& domain=t.source().domain();const auto n=source.integratedCharge().size();
    Stage result;result.report.cells=domain.cellIndexByRow().size();result.report.nativeFaces=t.faces().size();
    require(t.faces().size()==static_cast<std::size_t>(domain.nativeFaceCount()),"actual charge transport lacks native faces");
    result.faces.resize(t.faces().size());std::vector<std::vector<std::size_t>> outgoing(n);std::vector<unsigned> faceSeen(result.faces.size());
    FV::Dyadic boundaryProductError,rawIncomingMagnitude,boundaryNetError;
    struct Counterflow {int sign=0;long double incomingVolume=0,outgoingVolume=0;double incomingCanonical=0,outgoingCanonical=0;};
    std::map<std::size_t,Counterflow> counterflow;
    result.report.boundaryEnabled=boundary!=nullptr;
    for(const auto& f:t.faces()){
        require(f.nativeRow>=0&&static_cast<std::size_t>(f.nativeRow)<result.faces.size(),"actual charge face row outside native layout");
        require(++faceSeen[f.nativeRow]==1,"duplicate actual charge face row");const auto actual=domain.face(mesh,f.nativeRow);
        require(f.negativeCell==actual.negativeCell&&f.positiveCell==actual.positiveCell,"charge face incidence differs from actual native geometry");
        const auto volume=half?f.half.total:f.full.total;require(std::isfinite(volume),"nonfinite accepted total volume flux");
        auto& face=result.faces[f.nativeRow];face={f.nativeRow,f.negativeCell,f.positiveCell,volume,0};
        if(f.negativeCell<0&&f.positiveCell<0){requireInactiveFaceZero(f);continue;}
        if(f.negativeCell<0||f.positiveCell<0){
            if(!boundary){const auto& accepted=half?f.half:f.full;
                require(volume==0&&accepted.positiveTotal==0&&accepted.negativeTotal==0,
                    "closed charge entry rejects nonzero exterior total flow");continue;}
            const auto* fluid=t.advector().boundary();require(fluid!=nullptr,"open charge receipt lacks actual F1 boundary owner");
            const auto mode=fluid->rules().at(f.nativeRow);
            if(mode!=Interface::GeometricVofAdvectorSlab::Boundary::Mode::Closed)++result.report.boundaryReservoirFaces;
            const int inside=f.negativeCell>=0?f.negativeCell:f.positiveCell;
            require(inside>=0&&static_cast<std::size_t>(inside)<n,"charge exterior face lacks one actual inside donor");
            const int sign=f.negativeCell>=0?1:-1;const auto outward=sign>0?volume:-volume;
            const auto& accepted=half?f.half:f.full;
            if(accepted.positiveTotal>0&&accepted.negativeTotal>0){
                auto& both=counterflow[static_cast<std::size_t>(f.nativeRow)];both.sign=sign;
                both.incomingVolume=sign>0?accepted.negativeTotal:accepted.positiveTotal;
                both.outgoingVolume=sign>0?accepted.positiveTotal:accepted.negativeTotal;
                require(std::isfinite(both.incomingVolume)&&std::isfinite(both.outgoingVolume),"nonfinite accepted counterflow gross volume");
                const auto exact=FV::Dyadic::product({boundary->inflowChargeDensities().at(f.nativeRow),-sign*both.incomingVolume});
                const auto rounded=FV::roundDouble(exact);require(std::isfinite(rounded.value),"unrepresentable signed counterflow reservoir charge");
                both.incomingCanonical=rounded.value;boundaryProductError.add(rounded.error,sign);rawIncomingMagnitude.add(exact.absolute());
                result.report.maximumBoundaryProductRoundingError=std::max(result.report.maximumBoundaryProductRoundingError,represented(rounded.error.absolute()));
                outgoing[inside].push_back(static_cast<std::size_t>(f.nativeRow));
                continue;
            }
            if(outward<0){
                const auto exact=FV::Dyadic::product({boundary->inflowChargeDensities().at(f.nativeRow),volume});
                const auto rounded=FV::roundDouble(exact);require(std::isfinite(rounded.value),"unrepresentable signed reservoir charge flux");
                face.integratedCharge=rounded.value;boundaryProductError.add(rounded.error,sign);rawIncomingMagnitude.add(exact.absolute());
                result.report.maximumBoundaryProductRoundingError=std::max(result.report.maximumBoundaryProductRoundingError,represented(rounded.error.absolute()));
            }else if(outward>0)outgoing[inside].push_back(static_cast<std::size_t>(f.nativeRow));
            continue;
        }
        require(static_cast<std::size_t>(f.negativeCell)<n&&static_cast<std::size_t>(f.positiveCell)<n,"charge face cell outside native layout");
        if(volume!=0)outgoing[volume>0?f.negativeCell:f.positiveCell].push_back(static_cast<std::size_t>(f.nativeRow));
    }
    for(std::size_t i=0;i<n;++i){auto& rows=outgoing[i];if(rows.empty())continue;++result.report.donors;
        std::sort(rows.begin(),rows.end());const int row=domain.cellRowByIndex()[i];require(row>=0,"source charge lies outside complete domain");
        const double capacity=domain.cellVolumes()[row];require(std::isfinite(capacity)&&capacity>0,"invalid actual charge donor capacity");
        std::vector<FV::Dyadic> weights;FV::Dyadic out;weights.reserve(rows.size()+1);
        for(const auto fid:rows){const auto both=counterflow.find(fid);
            FV::Dyadic weight(both==counterflow.end()?std::abs(result.faces[fid].totalVolume):both->second.outgoingVolume);
            weights.push_back(weight);out.add(weight);}
        require(out.compare(FV::Dyadic(capacity))<=0,"charge donor outgoing actual volume exceeds physical capacity");
        FV::Dyadic retained(capacity);retained.add(out,-1);weights.push_back(retained);
        const auto allocation=FV::partitionDouble(std::abs(source.integratedCharge()[i]),weights);
        require(allocation.values.size()==weights.size()&&allocation.sumError.compare(FV::Dyadic{})==0,"binary64 donor partition did not conserve source");
        const bool negative=std::signbit(source.integratedCharge()[i]);FV::Dyadic own;
        for(std::size_t k=0;k<allocation.values.size();++k){own.add(allocation.exactAllocations[k]);
            if(weights[k].compare(FV::Dyadic{})>0&&allocation.values[k]==0&&source.integratedCharge()[i]!=0)++result.report.quantizedZeroPieces;}
        require(own.compare(FV::Dyadic(std::abs(source.integratedCharge()[i])))==0,"own donor charge partition mismatch");
        for(std::size_t k=0;k<rows.size();++k){const auto fid=rows[k];const double amount=std::copysign(allocation.values[k],negative?-1.:1.);
            const auto both=counterflow.find(fid);
            if(both!=counterflow.end())both->second.outgoingCanonical=both->second.sign>0?amount:-amount;
            else result.faces[fid].integratedCharge=result.faces[fid].totalVolume>0?amount:-amount;}
        result.report.maximumOutflowCfl=std::max(result.report.maximumOutflowCfl,static_cast<double>(out.rounded()/capacity));
    }
    FV::Dyadic counterOutgoingMagnitude,counterOutgoingBudget;
    for(const auto& entry:counterflow){const auto fid=entry.first;const auto& both=entry.second;
        FV::Dyadic combined(both.outgoingCanonical);combined.add(both.incomingCanonical);
        const auto rounded=FV::roundDouble(combined);require(std::isfinite(rounded.value),"unrepresentable canonical counterflow charge flux");
        result.faces[fid].integratedCharge=rounded.value;boundaryNetError.add(rounded.error,both.sign);
        result.report.maximumBoundaryNetRoundingError=std::max(result.report.maximumBoundaryNetRoundingError,represented(rounded.error.absolute()));
        counterOutgoingMagnitude.add(std::abs(both.outgoingCanonical));
        counterOutgoingBudget.add(FV::Dyadic::product({controls.relativeConservationTolerance,std::abs(both.outgoingCanonical)}));
    }
    std::vector<FV::Dyadic> net(n);Ledger global;FV::Dyadic before,after,rounding;
    FV::Dyadic incoming,outgoingCharge,boundaryNet;
    for(const auto& f:result.faces){const auto qt=f.integratedCharge;if(qt!=0)++result.report.nonzeroFaceCharges;
        if(f.negativeCell>=0)net[f.negativeCell].add(qt);if(f.positiveCell>=0)net[f.positiveCell].add(qt,-1);
        global.activity(qt,controls.relativeConservationTolerance);
        if(boundary&&((f.negativeCell>=0)!=(f.positiveCell>=0))){
            const int sign=f.negativeCell>=0?1:-1;global.residual.add(qt,sign);boundaryNet.add(qt,sign);
            const auto outward=sign>0?f.totalVolume:-f.totalVolume;const auto both=counterflow.find(static_cast<std::size_t>(f.nativeRow));
            if(both!=counterflow.end()){incoming.add(both->second.incomingCanonical,-sign);outgoingCharge.add(both->second.outgoingCanonical,sign);}
            else if(outward<0)incoming.add(qt,-sign);else if(outward>0)outgoingCharge.add(qt,sign);
        }}
    result.q.resize(n);
    for(std::size_t i=0;i<n;++i){const auto old=source.integratedCharge()[i];FV::Dyadic exact(old);exact.add(net[i],-1);
        FV::RoundedDouble rounded;if(net[i].compare(FV::Dyadic{})==0)rounded.value=old;else rounded=FV::roundDouble(exact);
        require(std::isfinite(rounded.value),"unrepresentable advected integrated charge");result.q[i]=rounded.value;rounding.add(rounded.error);
        result.report.maximumCellRoundingError=std::max(result.report.maximumCellRoundingError,represented(rounded.error.absolute()));
        before.add(old);after.add(rounded.value);global.term(old,-1,controls.relativeConservationTolerance);global.term(rounded.value,1,controls.relativeConservationTolerance);}
    require(global.residual.compare(rounding)==0,"charge ledger is not explained by final target rounding");global.validate();
    result.report.sourceTotalCharge=represented(before);result.report.candidateTotalCharge=represented(after);
    result.report.chargeResidual=represented(global.residual);result.report.relativeChargeResidual=global.relative();
    if(boundary){
        // Distinguish native shared-Qt conservation from the declared inlet
        // rho*T model. No rounded Report total participates in either gate.
        auto model=global;model.residual.add(boundaryProductError,-1);model.residual.add(boundaryNetError,-1);
        model.scale.add(rawIncomingMagnitude);model.scale.add(counterOutgoingMagnitude);model.budget.add(counterOutgoingBudget);
        // Multiply the exact inlet activity by tolerance without converting it
        // through a native subtotal (the helper has exact add/product inputs).
        for(const auto& f:result.faces)if((f.negativeCell>=0)!=(f.positiveCell>=0)){
            const int sign=f.negativeCell>=0?1:-1;const auto outward=sign>0?f.totalVolume:-f.totalVolume;
            const auto both=counterflow.find(static_cast<std::size_t>(f.nativeRow));
            if(both!=counterflow.end())model.budget.add(FV::Dyadic::product({controls.relativeConservationTolerance,
                std::abs(boundary->inflowChargeDensities().at(f.nativeRow)),both->second.incomingVolume}));
            else if(outward<0)model.budget.add(FV::Dyadic::product({controls.relativeConservationTolerance,
                std::abs(boundary->inflowChargeDensities().at(f.nativeRow)),std::abs(f.totalVolume)}));
        }
        FV::Dyadic expected=rounding;expected.add(boundaryProductError,-1);expected.add(boundaryNetError,-1);
        require(model.residual.compare(expected)==0,"declared inlet charge difference is not explained by product/cell RNE");model.validate();
        result.report.boundaryIncomingCharge=represented(incoming);result.report.boundaryOutgoingCharge=represented(outgoingCharge);
        result.report.boundaryNetCharge=represented(boundaryNet);result.report.boundaryProductRoundingResidual=represented(boundaryProductError);
        result.report.boundaryNetRoundingResidual=represented(boundaryNetError);
        result.report.boundaryModelChargeResidual=represented(model.residual);result.report.relativeBoundaryModelChargeResidual=model.relative();
        if(controls.chargeAdvection==ElectricVofChargeAdvectionScheme::LimitedLinearFct){
            result.boundaryModelExtra=std::make_unique<Ledger>(model);
            result.boundaryModelExtra->residual.add(global.residual,-1);
            result.boundaryModelExtra->scale.add(global.scale,-1);
            result.boundaryModelExtra->budget.add(global.budget,-1);
        }
    }
    result.report.roundingResidual=represented(rounding);result.report.totalSeconds=elapsed(started);return result;
}
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
Stage advanceParallel(const QuadtreeMesh& mesh,const ElectricVofCharge& source,const Transport& t,bool half,const ElectricVofControls& controls,
    const ElectricVofChargeBoundary* boundary)
{
    const auto started=Timer::now();const auto& domain=t.source().domain();const auto n=source.integratedCharge().size();
    ChargeAdvanceCleanup cleanup;
    Stage result;result.report.cells=domain.cellIndexByRow().size();result.report.nativeFaces=t.faces().size();
    require(t.faces().size()==static_cast<std::size_t>(domain.nativeFaceCount()),"actual charge transport lacks native faces");
    result.faces.resize(t.faces().size());std::vector<std::vector<std::size_t>> outgoing(n);std::vector<unsigned> faceSeen(result.faces.size());
    FV::Dyadic boundaryProductError,rawIncomingMagnitude,boundaryNetError;
    struct Counterflow {int sign=0;long double incomingVolume=0,outgoingVolume=0;double incomingCanonical=0,outgoingCanonical=0;};
    std::map<std::size_t,Counterflow> counterflow;
    result.report.boundaryEnabled=boundary!=nullptr;
    for(const auto& f:t.faces()){
        require(f.nativeRow>=0&&static_cast<std::size_t>(f.nativeRow)<result.faces.size(),"actual charge face row outside native layout");
        require(++faceSeen[f.nativeRow]==1,"duplicate actual charge face row");const auto actual=domain.face(mesh,f.nativeRow);
        require(f.negativeCell==actual.negativeCell&&f.positiveCell==actual.positiveCell,"charge face incidence differs from actual native geometry");
        const auto volume=half?f.half.total:f.full.total;require(std::isfinite(volume),"nonfinite accepted total volume flux");
        auto& face=result.faces[f.nativeRow];face={f.nativeRow,f.negativeCell,f.positiveCell,volume,0};
        if(f.negativeCell<0&&f.positiveCell<0){requireInactiveFaceZero(f);continue;}
        if(f.negativeCell<0||f.positiveCell<0){
            if(!boundary){const auto& accepted=half?f.half:f.full;
                require(volume==0&&accepted.positiveTotal==0&&accepted.negativeTotal==0,
                    "closed charge entry rejects nonzero exterior total flow");continue;}
            const auto* fluid=t.advector().boundary();require(fluid!=nullptr,"open charge receipt lacks actual F1 boundary owner");
            const auto mode=fluid->rules().at(f.nativeRow);
            if(mode!=Interface::GeometricVofAdvectorSlab::Boundary::Mode::Closed)++result.report.boundaryReservoirFaces;
            const int inside=f.negativeCell>=0?f.negativeCell:f.positiveCell;
            require(inside>=0&&static_cast<std::size_t>(inside)<n,"charge exterior face lacks one actual inside donor");
            const int sign=f.negativeCell>=0?1:-1;const auto outward=sign>0?volume:-volume;
            const auto& accepted=half?f.half:f.full;
            if(accepted.positiveTotal>0&&accepted.negativeTotal>0){
                auto& both=counterflow[static_cast<std::size_t>(f.nativeRow)];both.sign=sign;
                both.incomingVolume=sign>0?accepted.negativeTotal:accepted.positiveTotal;
                both.outgoingVolume=sign>0?accepted.positiveTotal:accepted.negativeTotal;
                require(std::isfinite(both.incomingVolume)&&std::isfinite(both.outgoingVolume),"nonfinite accepted counterflow gross volume");
                const auto exact=FV::Dyadic::product({boundary->inflowChargeDensities().at(f.nativeRow),-sign*both.incomingVolume});
                const auto rounded=FV::roundDouble(exact);require(std::isfinite(rounded.value),"unrepresentable signed counterflow reservoir charge");
                both.incomingCanonical=rounded.value;boundaryProductError.add(rounded.error,sign);rawIncomingMagnitude.add(exact.absolute());
                result.report.maximumBoundaryProductRoundingError=std::max(result.report.maximumBoundaryProductRoundingError,represented(rounded.error.absolute()));
                outgoing[inside].push_back(static_cast<std::size_t>(f.nativeRow));
                continue;
            }
            if(outward<0){
                const auto exact=FV::Dyadic::product({boundary->inflowChargeDensities().at(f.nativeRow),volume});
                const auto rounded=FV::roundDouble(exact);require(std::isfinite(rounded.value),"unrepresentable signed reservoir charge flux");
                face.integratedCharge=rounded.value;boundaryProductError.add(rounded.error,sign);rawIncomingMagnitude.add(exact.absolute());
                result.report.maximumBoundaryProductRoundingError=std::max(result.report.maximumBoundaryProductRoundingError,represented(rounded.error.absolute()));
            }else if(outward>0)outgoing[inside].push_back(static_cast<std::size_t>(f.nativeRow));
            continue;
        }
        require(static_cast<std::size_t>(f.negativeCell)<n&&static_cast<std::size_t>(f.positiveCell)<n,"charge face cell outside native layout");
        if(volume!=0)outgoing[volume>0?f.negativeCell:f.positiveCell].push_back(static_cast<std::size_t>(f.nativeRow));
    }
    const int routingErrno=errno;ChargeAdvanceFloatingState routingState;
    if(!routingState.capture()||errno!=routingErrno)throw std::runtime_error("optional charge setup FP receipt");
    // Only the two independent row loops are reordered. The actual face
    // prepass above and every counterflow/scatter/global loop below stay literal.
    struct DonorResult {std::size_t quantized=0;double cfl=0;bool active=false;};
    ChargeAdvanceBudget budget;std::size_t maxOutgoing=0;
    for(const auto& rows:outgoing)maxOutgoing=std::max(maxOutgoing,rows.size());
    const auto teamCount=std::min<std::size_t>(static_cast<std::size_t>(omp_get_max_threads()),std::max<std::size_t>(1,n));
    const int team=static_cast<int>(teamCount);
    budget.add(n,sizeof(DonorResult)+sizeof(double));budget.add(result.faces.size(),sizeof(Counterflow*));
    // Binary64 CDF inputs/retained differences span <4096 exact bits. Reserve
    // full row vectors, allocator growth and division temporaries for each
    // actual maximum donor width, not a guessed 2:1 face-degree limit.
    if(maxOutgoing==std::numeric_limits<std::size_t>::max())throw std::length_error("optional charge donor width");
    budget.add(teamCount+1,1024u*1024u);
    budget.add(teamCount,(maxOutgoing+1)>ChargeAdvanceBudget::limit/8192u?ChargeAdvanceBudget::limit:8192u*(maxOutgoing+1));
    std::unique_ptr<DonorResult[]> donorReport(new DonorResult[n]{});
    std::unique_ptr<double[]> targetError(new double[n]{});
    std::unique_ptr<Counterflow*[]> counterByFace(new Counterflow*[result.faces.size()]{});
    for(auto& entry:counterflow)counterByFace[entry.first]=&entry.second;
    if(errno!=routingErrno)throw std::runtime_error("optional charge setup errno");
    routingState.restore();errno=routingErrno;
    unsigned donorWorkers=0,targetWorkers=0;
    if(!chargeAdvanceRows(n,team,[&](std::size_t i){auto& rows=outgoing[i];if(rows.empty())return;donorReport[i].active=true;
        std::sort(rows.begin(),rows.end());const int row=domain.cellRowByIndex()[i];require(row>=0,"source charge lies outside complete domain");
        const double capacity=domain.cellVolumes()[row];require(std::isfinite(capacity)&&capacity>0,"invalid actual charge donor capacity");
        std::vector<FV::Dyadic> weights;FV::Dyadic out;weights.reserve(rows.size()+1);
        for(const auto fid:rows){const auto* both=counterByFace[fid];
            FV::Dyadic weight(!both?std::abs(result.faces[fid].totalVolume):both->outgoingVolume);
            weights.push_back(weight);out.add(weight);}
        require(out.compare(FV::Dyadic(capacity))<=0,"charge donor outgoing actual volume exceeds physical capacity");
        FV::Dyadic retained(capacity);retained.add(out,-1);weights.push_back(retained);
        const auto allocation=FV::partitionDouble(std::abs(source.integratedCharge()[i]),weights);
        require(allocation.values.size()==weights.size()&&allocation.sumError.compare(FV::Dyadic{})==0,"binary64 donor partition did not conserve source");
        const bool negative=std::signbit(source.integratedCharge()[i]);FV::Dyadic own;
        for(std::size_t k=0;k<allocation.values.size();++k){own.add(allocation.exactAllocations[k]);
            if(weights[k].compare(FV::Dyadic{})>0&&allocation.values[k]==0&&source.integratedCharge()[i]!=0)++donorReport[i].quantized;}
        require(own.compare(FV::Dyadic(std::abs(source.integratedCharge()[i])))==0,"own donor charge partition mismatch");
        for(std::size_t k=0;k<rows.size();++k){const auto fid=rows[k];const double amount=std::copysign(allocation.values[k],negative?-1.:1.);
            auto* both=counterByFace[fid];
            if(both)both->outgoingCanonical=both->sign>0?amount:-amount;
            else result.faces[fid].integratedCharge=result.faces[fid].totalVolume>0?amount:-amount;}
        donorReport[i].cfl=static_cast<double>(out.rounded()/capacity);
    },donorWorkers))throw std::runtime_error("optional charge donor worker unavailable");
    for(std::size_t i=0;i<n;++i)if(donorReport[i].active){++result.report.donors;
        result.report.quantizedZeroPieces+=donorReport[i].quantized;
        result.report.maximumOutflowCfl=std::max(result.report.maximumOutflowCfl,donorReport[i].cfl);}
    FV::Dyadic counterOutgoingMagnitude,counterOutgoingBudget;
    for(const auto& entry:counterflow){const auto fid=entry.first;const auto& both=entry.second;
        FV::Dyadic combined(both.outgoingCanonical);combined.add(both.incomingCanonical);
        const auto rounded=FV::roundDouble(combined);require(std::isfinite(rounded.value),"unrepresentable canonical counterflow charge flux");
        result.faces[fid].integratedCharge=rounded.value;boundaryNetError.add(rounded.error,both.sign);
        result.report.maximumBoundaryNetRoundingError=std::max(result.report.maximumBoundaryNetRoundingError,represented(rounded.error.absolute()));
        counterOutgoingMagnitude.add(std::abs(both.outgoingCanonical));
        counterOutgoingBudget.add(FV::Dyadic::product({controls.relativeConservationTolerance,std::abs(both.outgoingCanonical)}));
    }
    std::vector<FV::Dyadic> net(n);Ledger global;FV::Dyadic before,after,rounding;
    FV::Dyadic incoming,outgoingCharge,boundaryNet;
    for(const auto& f:result.faces){const auto qt=f.integratedCharge;if(qt!=0)++result.report.nonzeroFaceCharges;
        if(f.negativeCell>=0)net[f.negativeCell].add(qt);if(f.positiveCell>=0)net[f.positiveCell].add(qt,-1);
        global.activity(qt,controls.relativeConservationTolerance);
        if(boundary&&((f.negativeCell>=0)!=(f.positiveCell>=0))){
            const int sign=f.negativeCell>=0?1:-1;global.residual.add(qt,sign);boundaryNet.add(qt,sign);
            const auto outward=sign>0?f.totalVolume:-f.totalVolume;const auto both=counterflow.find(static_cast<std::size_t>(f.nativeRow));
            if(both!=counterflow.end()){incoming.add(both->second.incomingCanonical,-sign);outgoingCharge.add(both->second.outgoingCanonical,sign);}
            else if(outward<0)incoming.add(qt,-sign);else if(outward>0)outgoingCharge.add(qt,sign);
        }}
    result.q.resize(n);
    if(!chargeAdvanceRows(n,team,[&](std::size_t i){const auto old=source.integratedCharge()[i];FV::Dyadic exact(old);exact.add(net[i],-1);
        FV::RoundedDouble rounded;if(net[i].compare(FV::Dyadic{})==0)rounded.value=old;else rounded=FV::roundDouble(exact);
        require(std::isfinite(rounded.value),"unrepresentable advected integrated charge");result.q[i]=rounded.value;
        targetError[i]=represented(rounded.error.absolute());
        // net is dead after this target closure. Reuse its PRIVATE existing
        // slot for the exact RN error rather than retain a second Dyadic array.
        net[i]=std::move(rounded.error);
    },targetWorkers))throw std::runtime_error("optional charge target worker unavailable");
    for(std::size_t i=0;i<n;++i){const auto old=source.integratedCharge()[i];rounding.add(net[i]);
        result.report.maximumCellRoundingError=std::max(result.report.maximumCellRoundingError,targetError[i]);
        before.add(old);after.add(result.q[i]);global.term(old,-1,controls.relativeConservationTolerance);
        global.term(result.q[i],1,controls.relativeConservationTolerance);}
    require(global.residual.compare(rounding)==0,"charge ledger is not explained by final target rounding");global.validate();
    result.report.sourceTotalCharge=represented(before);result.report.candidateTotalCharge=represented(after);
    result.report.chargeResidual=represented(global.residual);result.report.relativeChargeResidual=global.relative();
    if(boundary){
        // Distinguish native shared-Qt conservation from the declared inlet
        // rho*T model. No rounded Report total participates in either gate.
        auto model=global;model.residual.add(boundaryProductError,-1);model.residual.add(boundaryNetError,-1);
        model.scale.add(rawIncomingMagnitude);model.scale.add(counterOutgoingMagnitude);model.budget.add(counterOutgoingBudget);
        // Multiply the exact inlet activity by tolerance without converting it
        // through a native subtotal (the helper has exact add/product inputs).
        for(const auto& f:result.faces)if((f.negativeCell>=0)!=(f.positiveCell>=0)){
            const int sign=f.negativeCell>=0?1:-1;const auto outward=sign>0?f.totalVolume:-f.totalVolume;
            const auto both=counterflow.find(static_cast<std::size_t>(f.nativeRow));
            if(both!=counterflow.end())model.budget.add(FV::Dyadic::product({controls.relativeConservationTolerance,
                std::abs(boundary->inflowChargeDensities().at(f.nativeRow)),both->second.incomingVolume}));
            else if(outward<0)model.budget.add(FV::Dyadic::product({controls.relativeConservationTolerance,
                std::abs(boundary->inflowChargeDensities().at(f.nativeRow)),std::abs(f.totalVolume)}));
        }
        FV::Dyadic expected=rounding;expected.add(boundaryProductError,-1);expected.add(boundaryNetError,-1);
        require(model.residual.compare(expected)==0,"declared inlet charge difference is not explained by product/cell RNE");model.validate();
        result.report.boundaryIncomingCharge=represented(incoming);result.report.boundaryOutgoingCharge=represented(outgoingCharge);
        result.report.boundaryNetCharge=represented(boundaryNet);result.report.boundaryProductRoundingResidual=represented(boundaryProductError);
        result.report.boundaryNetRoundingResidual=represented(boundaryNetError);
        result.report.boundaryModelChargeResidual=represented(model.residual);result.report.relativeBoundaryModelChargeResidual=model.relative();
        if(controls.chargeAdvection==ElectricVofChargeAdvectionScheme::LimitedLinearFct){
            result.boundaryModelExtra=std::make_unique<Ledger>(model);
            result.boundaryModelExtra->residual.add(global.residual,-1);
            result.boundaryModelExtra->scale.add(global.scale,-1);
            result.boundaryModelExtra->budget.add(global.budget,-1);
        }
    }
    result.report.roundingResidual=represented(rounding);
    result.report.donorParallelUsed=donorWorkers>1&&result.report.donors>0;
    result.report.targetParallelUsed=targetWorkers>1&&n>0;
    result.report.parallelWorkers=std::max(donorWorkers,targetWorkers);result.report.parallelScratchBytes=budget.used;
    result.report.totalSeconds=elapsed(started);
    if(!cleanup.arm())throw std::runtime_error("optional charge final FP receipt");
    static_assert(noexcept(Stage(std::move(result))),"charge Stage publication must not throw");
    return result;
}
#endif
Stage advanceLow(const QuadtreeMesh& mesh,const ElectricVofCharge& source,const Transport& t,bool half,const ElectricVofControls& controls,
    const ElectricVofChargeBoundary* boundary=nullptr)
{
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    const int error=errno;ChargeAdvanceFloatingState entry;const bool qualified=chargeAdvanceQualified(entry);
    if(entry.valid)entry.restore();errno=error;
    if(qualified)try{return advanceParallel(mesh,source,t,half,controls,boundary);}
    catch(...){
        // Whole private attempt has unwound. No receipt/Q identity is published.
        // The literal function reproduces its original first error; no PDE exists.
        entry.restore();errno=error;
        auto replay=advanceLiteral(mesh,source,t,half,controls,boundary);
        replay.report.parallelFallback=true;return replay;
    }
#endif
    return advanceLiteral(mesh,source,t,half,controls,boundary);
}
// Signed unsplit P1 candidate, then a shared-face multidimensional FCT.
// Physical boundaries deliberately retain the authenticated P0 gross model.
double fctRatio(const FV::Dyadic& available,const FV::Dyadic& requested)
{
    const FV::Dyadic zero;
    require(available.compare(zero)>=0&&requested.compare(zero)>=0,"negative FCT budget/activity");
    if(requested.compare(zero)==0||available.compare(requested)>=0)return 1;
    if(available.compare(zero)==0)return 0;
    auto p=requested,q=available;const FV::Dyadic largest(std::numeric_limits<double>::max());
    while(p.compare(largest)>0){p=p.multiplied(.5L);q=q.multiplied(.5L);}
    double numerator=represented(q),denominator=represented(p);
    if(FV::Dyadic(numerator).compare(q)>0)numerator=std::nextafter(numerator,0.);
    if(FV::Dyadic(denominator).compare(p)<0)denominator=std::nextafter(denominator,std::numeric_limits<double>::infinity());
    require(std::isfinite(numerator)&&std::isfinite(denominator)&&denominator>0,"invalid FCT ratio representation");
    double ratio=std::min(1.,numerator/denominator);
    while(requested.multiplied(ratio).compare(available)>0)ratio=std::nextafter(ratio,0.);
    require(std::isfinite(ratio)&&ratio>=0&&ratio<=1,"invalid shared FCT coefficient");return ratio;
}
Stage advanceFct(const QuadtreeMesh& mesh,const ElectricVofCharge& source,const Transport& t,bool half,const ElectricVofControls& controls,
    const ElectricVofChargeBoundary* boundary,bool threaded,bool replay=false)
{
    const auto started=Timer::now();
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    ChargeAdvanceCleanup cleanup;
#endif
    auto result=replay?advanceLiteral(mesh,source,t,half,controls,boundary):advanceLow(mesh,source,t,half,controls,boundary);
    for(auto& face:result.faces)face.lowOrderIntegratedCharge=face.highOrderIntegratedCharge=face.integratedCharge;
    const auto& domain=t.source().domain();const auto& old=source.integratedCharge();const auto n=old.size();
    struct Cell {
        double volume=0,rho=0,lowRho=0,minimum=0,maximum=0,slopeMinimum=0,slopeMaximum=0;
        long double x=0,y=0,hx=0,hy=0,a=0,b=0,c=0,bx=0,by=0,gx=0,gy=0;
        long double displacementX=0,displacementY=0,dilation=0;
    };
    std::vector<Cell> cells(n);std::vector<Vector2D> points(result.faces.size());
    // Each pair changes role only after every cell has consumed its old role:
    // low net/exact -> positive/negative activity -> final net/RNE error.
    std::vector<FV::Dyadic> net(n),scratch(n),anti(result.faces.size()),upper(n),lower(n);
    std::vector<double> nativeLower(n),nativeUpper(n),rPlus(n,1.),rMinus(n,1.);
    std::vector<std::size_t> offsets,adjacent;
    struct Totals {
        Ledger global;FV::Dyadic before,after,rounding;
        std::size_t highOrderFaces=0,limitedFaces=0,nonzeroFaceCharges=0;
        double minimumFaceLimiter=1,maximumCellRoundingError=0,maximumDensityBoundDefect=0;
    } totals;
    std::vector<Totals> chunks;
    const bool rz=domain.axisymmetric();
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    std::atomic<int> failed{0};unsigned actualWorkers=0;std::size_t scratchBytes=0;int team=1;
    if(threaded){
        const int setupErrno=errno;ChargeAdvanceFloatingState setup;
        if(!setup.capture())throw std::runtime_error("optional FCT setup FP receipt");
        const auto count=std::min<std::size_t>(static_cast<std::size_t>(omp_get_max_threads()),std::max<std::size_t>(1,n));
        team=static_cast<int>(count);ChargeAdvanceBudget budget;
        budget.add(n+1,2*sizeof(std::size_t));budget.add(t.faces().size(),2*sizeof(std::size_t));
        budget.add(count+1,1024u*1024u);budget.add(count,sizeof(Totals));scratchBytes=budget.used;
        chunks.resize(count);offsets.resize(n+1);
        for(const auto& f:t.faces())for(int ci:{f.negativeCell,f.positiveCell})if(ci>=0)++offsets[static_cast<std::size_t>(ci)+1];
        std::partial_sum(offsets.begin(),offsets.end(),offsets.begin());adjacent.resize(offsets.back());
        auto cursor=offsets;
        // Keep each cell's original transport-face order: WLS is rounded math.
        for(std::size_t fi=0;fi<t.faces().size();++fi){const auto& f=t.faces()[fi];
            for(int ci:{f.negativeCell,f.positiveCell})if(ci>=0)adjacent[cursor[ci]++]=fi;}
        if(errno!=setupErrno)throw std::runtime_error("optional FCT setup errno");
        setup.restore();errno=setupErrno;
    }
#else
    require(!threaded,"threaded FCT is unavailable");
#endif
    const auto rows=[&](std::size_t count,const auto& kernel){
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        if(threaded){
#pragma omp for schedule(static)
            for(std::ptrdiff_t row=0;row<static_cast<std::ptrdiff_t>(count);++row){
                if(failed.load(std::memory_order_relaxed))continue;
                try{kernel(static_cast<std::size_t>(row));}catch(...){failed.store(1,std::memory_order_relaxed);}
            }
            return;
        }
#endif
        for(std::size_t row=0;row<count;++row)kernel(row);
    };
    const auto sample=[&](int ci,long double x,long double y,double rho,double lowRho){auto& c=cells[ci];
        const auto dx=(x-c.x)/c.hx,dy=(y-c.y)/c.hy;const auto square=dx*dx+dy*dy;
        require(std::isfinite(square)&&square>0,"invalid native FCT reconstruction sample");
        const auto w=1/square,difference=static_cast<long double>(rho)-c.rho;
        c.a+=w*dx*dx;c.b+=w*dx*dy;c.c+=w*dy*dy;c.bx+=w*dx*difference;c.by+=w*dy*difference;
        c.minimum=std::min({c.minimum,rho,lowRho});c.maximum=std::max({c.maximum,rho,lowRho});
        c.slopeMinimum=std::min(c.slopeMinimum,rho);c.slopeMaximum=std::max(c.slopeMaximum,rho);
    };
    const auto facePoint=[&](std::size_t fi){const auto& f=t.faces()[fi];
        if(f.negativeCell<0&&f.positiveCell<0){requireInactiveFaceZero(f);return;}
        const auto actual=domain.face(mesh,f.nativeRow);auto point=actual.center;
        if(rz&&actual.horizontal){const long double r=point.x,h=actual.tangentSpacing;
            require(r>0,"FCT horizontal RZ face has no weighted centroid");point.x=static_cast<double>(r+h*(h/r/12));}
        points[f.nativeRow]=point;
    };
    const auto cellFace=[&](int ci,const auto& f,int sign){if(ci<0)return;auto& c=cells[ci];
        const auto point=points[f.nativeRow];const auto volume=result.faces[f.nativeRow].totalVolume;
        const auto swept=sign*volume;c.displacementX+=(swept/c.volume)*((static_cast<long double>(point.x)-c.x)/c.hx);
        c.displacementY+=(swept/c.volume)*((static_cast<long double>(point.y)-c.y)/c.hy);c.dilation+=swept/c.volume;
        net[ci].add(result.faces[f.nativeRow].integratedCharge,sign);
    };
    const auto faceSample=[&](int ci,const auto& f){
        if(f.negativeCell>=0&&f.positiveCell>=0){const int other=ci==f.negativeCell?f.positiveCell:f.negativeCell;const auto& c=cells[other];
            sample(ci,c.x,c.y,c.rho,c.lowRho);
        }else if(boundary){const auto& accepted=half?f.half:f.full;const int sign=f.negativeCell>=0?1:-1;
            const auto incoming=sign>0?accepted.negativeTotal:accepted.positiveTotal;
            if(incoming>0){const double rho=boundary->inflowChargeDensities().at(f.nativeRow);const auto point=points[f.nativeRow];sample(ci,point.x,point.y,rho,rho);}}
    };
    auto& global=totals.global;auto& before=totals.before;auto& after=totals.after;auto& rounding=totals.rounding;
    const auto foldFace=[&](std::size_t fi,Totals& sum){const auto& f=result.faces[fi];
        if(f.negativeCell>=0&&f.positiveCell>=0&&f.totalVolume!=0)++sum.highOrderFaces;
        if(anti[f.nativeRow].compare(FV::Dyadic{})!=0){if(f.limiter<1)++sum.limitedFaces;
            sum.minimumFaceLimiter=std::min(sum.minimumFaceLimiter,f.limiter);}
        const auto value=f.integratedCharge;if(value!=0)++sum.nonzeroFaceCharges;
        if(!threaded){if(f.negativeCell>=0)net[f.negativeCell].add(value);if(f.positiveCell>=0)net[f.positiveCell].add(value,-1);}
        sum.global.activity(value,controls.relativeConservationTolerance);
        if((f.negativeCell>=0)!=(f.positiveCell>=0))sum.global.residual.add(value,f.negativeCell>=0?1:-1);
    };
    const auto finishCell=[&](std::size_t ci,Totals& sum){
        if(domain.cellRowByIndex()[ci]<0)return;
        sum.rounding.add(scratch[ci]);sum.before.add(old[ci]);sum.after.add(result.q[ci]);
        sum.maximumCellRoundingError=std::max(sum.maximumCellRoundingError,represented(scratch[ci].absolute()));
        sum.global.term(old[ci],-1,controls.relativeConservationTolerance);sum.global.term(result.q[ci],1,controls.relativeConservationTolerance);
        const auto density=result.q[ci]/cells[ci].volume;require(std::isfinite(density),"nonfinite limited charge density");
        sum.maximumDensityBoundDefect=std::max({sum.maximumDensityBoundDefect,
            cells[ci].minimum-density,density-cells[ci].maximum});
    };
    const auto targetCell=[&](std::size_t ci,Totals& sum){
        if(threaded)for(auto k=offsets[ci];k<offsets[ci+1];++k){const auto& f=t.faces()[adjacent[k]];
            net[ci].add(result.faces[f.nativeRow].integratedCharge,f.negativeCell==static_cast<int>(ci)?1:-1);}
        if(domain.cellRowByIndex()[ci]<0){require(old[ci]==0&&net[ci].compare(FV::Dyadic{})==0,"inactive FCT target received charge");result.q[ci]=old[ci];return;}
        FV::Dyadic exact(old[ci]);exact.add(net[ci],-1);FV::RoundedDouble rounded;
        if(net[ci].compare(FV::Dyadic{})==0)rounded.value=old[ci];else rounded=FV::roundDouble(exact);
        require(std::isfinite(rounded.value)&&rounded.value>=nativeLower[ci]&&rounded.value<=nativeUpper[ci],"FCT target exceeds native signed local bounds");
        result.q[ci]=rounded.value;scratch[ci]=std::move(rounded.error);finishCell(ci,sum);
    };
    const auto pipeline=[&]{
        rows(n,[&](std::size_t ci){
            if(domain.cellRowByIndex()[ci]<0){require(old[ci]==0&&result.q[ci]==0,"inactive FCT cell cannot contain charge");return;}
            auto& c=cells[ci];const auto& native=mesh.cells(mesh.cellIndexToLabel[ci]);
            const auto center=native.centerCoord(),h=native.size();const auto row=domain.cellRowByIndex().at(ci);
            require(row>=0,"FCT source outside complete native domain");c.volume=domain.cellVolumes().at(row);
            c.hx=h.x;c.hy=h.y;c.x=center.x;c.y=center.y;
            if(rz){require(c.x>0,"FCT RZ cell centroid is not positive");c.x+=c.hx*(c.hx/c.x/12);}
            c.rho=old[ci]/c.volume;c.lowRho=result.q[ci]/c.volume;
            require(std::isfinite(c.rho)&&std::isfinite(c.lowRho),"unrepresentable signed Q/V for linear charge advection");
            c.minimum=std::min(c.rho,c.lowRho);c.maximum=std::max(c.rho,c.lowRho);c.slopeMinimum=c.slopeMaximum=c.rho;
        });
        if(threaded){
            rows(t.faces().size(),facePoint);
            rows(n,[&](std::size_t ci){for(auto k=offsets[ci];k<offsets[ci+1];++k){const auto& f=t.faces()[adjacent[k]];
                cellFace(static_cast<int>(ci),f,f.negativeCell==static_cast<int>(ci)?1:-1);faceSample(static_cast<int>(ci),f);}});
        }else for(std::size_t fi=0;fi<t.faces().size();++fi){const auto& f=t.faces()[fi];facePoint(fi);
            if(f.negativeCell<0&&f.positiveCell<0)continue;
            cellFace(f.negativeCell,f,1);cellFace(f.positiveCell,f,-1);
            if(f.negativeCell>=0)faceSample(f.negativeCell,f);if(f.positiveCell>=0)faceSample(f.positiveCell,f);
        }
        rows(n,[&](std::size_t ci){if(domain.cellRowByIndex()[ci]<0)return;
            auto& c=cells[ci];const auto determinant=c.a*c.c-c.b*c.b;
            if(c.a>0&&c.c>0&&determinant>64*std::numeric_limits<long double>::epsilon()*c.a*c.c){
                c.gx=(c.c*c.bx-c.b*c.by)/determinant;c.gy=(c.a*c.by-c.b*c.bx)/determinant;
            }else{c.gx=c.a>0?c.bx/c.a:0;c.gy=c.c>0?c.by/c.c:0;}
            require(std::isfinite(c.gx)&&std::isfinite(c.gy),"unrepresentable charge reconstruction gradient");
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);const auto center=cell.centerCoord(),h=cell.size();long double limiter=1;
            for(const long double x:{static_cast<long double>(center.x)-h.x/2,static_cast<long double>(center.x)+h.x/2})
                for(const long double y:{static_cast<long double>(center.y)-h.y/2,static_cast<long double>(center.y)+h.y/2}){
                    const auto delta=c.gx*(x-c.x)/c.hx+c.gy*(y-c.y)/c.hy;
                    if(delta>0)limiter=std::min(limiter,(static_cast<long double>(c.slopeMaximum)-c.rho)/delta);
                    else if(delta<0)limiter=std::min(limiter,(static_cast<long double>(c.slopeMinimum)-c.rho)/delta);}
            limiter=std::clamp(limiter,0.L,1.L);c.gx*=limiter;c.gy*=limiter;
            scratch[ci]=FV::Dyadic(old[ci]);scratch[ci].add(net[ci],-1);
            const auto densityLower=FV::Dyadic::product({c.volume,c.minimum}),densityUpper=FV::Dyadic::product({c.volume,c.maximum});
            auto lowerValue=represented(densityLower),upperValue=represented(densityUpper);
            if(FV::Dyadic(lowerValue).compare(densityLower)<0)lowerValue=std::nextafter(lowerValue,std::numeric_limits<double>::infinity());
            if(FV::Dyadic(upperValue).compare(densityUpper)>0)upperValue=std::nextafter(upperValue,-std::numeric_limits<double>::infinity());
            nativeLower[ci]=std::min(result.q[ci],lowerValue);nativeUpper[ci]=std::max(result.q[ci],upperValue);
            auto exactLower=FV::Dyadic(nativeLower[ci]),exactUpper=FV::Dyadic(nativeUpper[ci]);
            if(scratch[ci].compare(exactLower)<0)exactLower=scratch[ci];if(scratch[ci].compare(exactUpper)>0)exactUpper=scratch[ci];
            upper[ci]=exactUpper;upper[ci].add(scratch[ci],-1);lower[ci]=scratch[ci];lower[ci].add(exactLower,-1);
            net[ci]=FV::Dyadic{};scratch[ci]=FV::Dyadic{};
        });
        rows(result.faces.size(),[&](std::size_t fi){auto& f=result.faces[fi];
            if(f.negativeCell<0||f.positiveCell<0||f.totalVolume==0)return;
            const int donor=f.totalVolume>0?f.negativeCell:f.positiveCell;const auto& c=cells[donor];const auto p=points[f.nativeRow];
            const long double predicted=c.rho+c.gx*((static_cast<long double>(p.x)-c.x)/c.hx-c.displacementX/2)
                +c.gy*((static_cast<long double>(p.y)-c.y)/c.hy-c.displacementY/2)-c.rho*c.dilation/2;
            require(std::isfinite(predicted),"unrepresentable multidimensional charge face prediction");
            f.highOrderIntegratedCharge=represented(FV::Dyadic::product({f.totalVolume,predicted}));
            auto& a=anti[f.nativeRow];a=FV::Dyadic(f.highOrderIntegratedCharge);a.add(f.lowOrderIntegratedCharge,-1);
            if(!threaded){const int sign=a.compare(FV::Dyadic{});
                if(sign>0){scratch[f.negativeCell].add(a);net[f.positiveCell].add(a);}
                else if(sign<0){net[f.negativeCell].add(a,-1);scratch[f.positiveCell].add(a,-1);}}
        });
        rows(n,[&](std::size_t ci){if(domain.cellRowByIndex()[ci]<0)return;
            if(threaded)for(auto k=offsets[ci];k<offsets[ci+1];++k){const auto& f=t.faces()[adjacent[k]];
                const auto& a=anti[f.nativeRow];const int orientation=f.negativeCell==static_cast<int>(ci)?1:-1;
                const int sign=orientation*a.compare(FV::Dyadic{});
                if(sign>0)scratch[ci].add(a,orientation);else if(sign<0)net[ci].add(a,-orientation);}
            rPlus[ci]=fctRatio(upper[ci],net[ci]);rMinus[ci]=fctRatio(lower[ci],scratch[ci]);
            net[ci]=FV::Dyadic{};scratch[ci]=FV::Dyadic{};
        });
        rows(result.faces.size(),[&](std::size_t fi){auto& f=result.faces[fi];const auto& a=anti[f.nativeRow];const int sign=a.compare(FV::Dyadic{});if(!sign)return;
            f.limiter=sign>0?std::min(rMinus[f.negativeCell],rPlus[f.positiveCell]):std::min(rPlus[f.negativeCell],rMinus[f.positiveCell]);
            const auto allowed=a.multiplied(f.limiter);FV::Dyadic exact(f.lowOrderIntegratedCharge);exact.add(allowed);
            auto value=represented(exact);FV::Dyadic delta(value);delta.add(f.lowOrderIntegratedCharge,-1);
            if(delta.absolute().compare(allowed.absolute())>0){value=std::nextafter(value,f.lowOrderIntegratedCharge);
                delta=FV::Dyadic(value);delta.add(f.lowOrderIntegratedCharge,-1);}
            require(delta.absolute().compare(allowed.absolute())<=0&&(delta.compare(FV::Dyadic{})==0||delta.compare(FV::Dyadic{})==sign),
                "represented FCT face exceeds its signed shared budget");f.integratedCharge=value;
        });
        if(threaded)rows(chunks.size(),[&](std::size_t chunk){auto& sum=chunks[chunk];
            // Fixed contiguous ranges preserve each scalar expression; only
            // exact integer subtotals cross chunks, never rounded double sums.
            const auto foldRange=[&](std::size_t count,const auto& kernel){
                const auto width=count/chunks.size(),extra=count%chunks.size();
                const auto first=chunk*width+std::min(chunk,extra),last=first+width+(chunk<extra);
                for(auto row=first;row<last;++row)kernel(row,sum);
            };
            foldRange(result.faces.size(),foldFace);foldRange(n,targetCell);
        });
        else{
            for(std::size_t fi=0;fi<result.faces.size();++fi)foldFace(fi,totals);
            rows(n,[&](std::size_t ci){targetCell(ci,totals);});
        }
    };
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    if(threaded){
        const int parentErrno=errno;ChargeAdvanceFloatingState prepared;
        if(!prepared.capture())throw std::runtime_error("optional FCT worker FP receipt");
        unsigned long raisedNative=0;unsigned raisedCsr=0;
#pragma omp parallel num_threads(team)
        {
            ChargeAdvanceFloatingState worker;const int workerErrno=errno;
            const bool captured=worker.capture(),ready=captured&&worker.enterFresh(prepared);
            if(!ready)failed.store(1,std::memory_order_relaxed);if(ready)errno=parentErrno;
#pragma omp single
            {actualWorkers=static_cast<unsigned>(omp_get_num_threads());}
            pipeline();
            if(ready){const int numericalErrno=errno;ChargeAdvanceFloatingState after;
                const bool got=after.capture();
                if(!got||numericalErrno!=parentErrno||!chargeAdvanceSameControl(after,prepared))failed.store(1,std::memory_order_relaxed);
#pragma omp critical(lamp_charge_advance_floating_flags)
                {if(got){raisedNative|=after.native._Fe_stat;raisedCsr|=after.csr&0x3fu;}}
            }
            if(captured)worker.restore();errno=workerErrno;
        }
        if(failed.load(std::memory_order_relaxed))throw std::runtime_error("optional FCT worker unavailable");
        auto merged=prepared.native;merged._Fe_stat|=raisedNative;
        if(std::fesetenv(&merged)!=0)std::terminate();_mm_setcsr(prepared.csr|raisedCsr);errno=parentErrno;
        result.report.targetParallelUsed=result.report.targetParallelUsed||actualWorkers>1;
        result.report.parallelWorkers=std::max(result.report.parallelWorkers,static_cast<std::size_t>(actualWorkers));
        result.report.parallelScratchBytes=std::max(result.report.parallelScratchBytes,scratchBytes);
    }else
#endif
        pipeline();
    for(const auto& sum:chunks){
        global.residual.add(sum.global.residual);global.scale.add(sum.global.scale);global.budget.add(sum.global.budget);
        before.add(sum.before);after.add(sum.after);rounding.add(sum.rounding);
        totals.highOrderFaces+=sum.highOrderFaces;totals.limitedFaces+=sum.limitedFaces;totals.nonzeroFaceCharges+=sum.nonzeroFaceCharges;
        totals.minimumFaceLimiter=std::min(totals.minimumFaceLimiter,sum.minimumFaceLimiter);
        totals.maximumCellRoundingError=std::max(totals.maximumCellRoundingError,sum.maximumCellRoundingError);
        totals.maximumDensityBoundDefect=std::max(totals.maximumDensityBoundDefect,sum.maximumDensityBoundDefect);
    }
    result.report.limitedLinearFct=true;result.report.reconstructedCells=domain.cellIndexByRow().size();
    result.report.highOrderFaces=totals.highOrderFaces;result.report.limitedFaces=totals.limitedFaces;result.report.nonzeroFaceCharges=totals.nonzeroFaceCharges;
    result.report.minimumFaceLimiter=totals.minimumFaceLimiter;result.report.maximumCellRoundingError=totals.maximumCellRoundingError;
    result.report.maximumDensityBoundDefect=totals.maximumDensityBoundDefect;
    require(global.residual.compare(rounding)==0,"FCT shared charge ledger differs from final target RNE");global.validate();
    result.report.sourceTotalCharge=represented(before);result.report.candidateTotalCharge=represented(after);
    result.report.chargeResidual=represented(global.residual);result.report.relativeChargeResidual=global.relative();
    result.report.roundingResidual=represented(rounding);
    if(boundary){require(bool(result.boundaryModelExtra),"FCT lost the authenticated boundary model ledger");auto model=global;
        model.residual.add(result.boundaryModelExtra->residual);model.scale.add(result.boundaryModelExtra->scale);model.budget.add(result.boundaryModelExtra->budget);
        model.validate();result.report.boundaryModelChargeResidual=represented(model.residual);result.report.relativeBoundaryModelChargeResidual=model.relative();}
    result.report.totalSeconds=elapsed(started);
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    if(threaded&&!cleanup.arm())throw std::runtime_error("optional FCT final FP receipt");
#endif
    return result;
}
Stage advance(const QuadtreeMesh& mesh,const ElectricVofCharge& source,const Transport& t,bool half,const ElectricVofControls& controls,
    const ElectricVofChargeBoundary* boundary=nullptr)
{
    if(controls.chargeAdvection==ElectricVofChargeAdvectionScheme::DonorCell){
        auto result=advanceLow(mesh,source,t,half,controls,boundary);
        for(auto& face:result.faces)face.lowOrderIntegratedCharge=face.highOrderIntegratedCharge=face.integratedCharge;return result;
    }
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    const int error=errno;ChargeAdvanceFloatingState entry;
    const auto workers=static_cast<std::size_t>(std::max(1,omp_get_max_threads()));
    // Each worker needs a complete batch of cells, not a few isolated faces.
    const bool qualified=source.integratedCharge().size()/workers>=16&&chargeAdvanceQualified(entry);
    if(entry.valid)entry.restore();errno=error;
    if(qualified)try{return advanceFct(mesh,source,t,half,controls,boundary,true);}
    catch(...){entry.restore();errno=error;
        auto result=advanceFct(mesh,source,t,half,controls,boundary,false,true);result.report.parallelFallback=true;return result;}
#endif
    return advanceFct(mesh,source,t,half,controls,boundary,false);
}

void remapValid(const QuadtreeMesh& sourceMesh,const Vof& geometry,const ElectricVofCharge& source,
    const QuadtreeMesh& targetMesh,const Transfer& transfer)
{
    const auto clock=source.timeLayer();clock.time();require(clock.position==Clock::Position::End,"charge AMR requires actual accepted End");
    source.validateFor(sourceMesh,clock);geometry.validateFor(sourceMesh,geometry.identity(),clock);
    require(geometry.geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic&&transfer.report().localPlic,
        "charge AMR requires actual LocalPlic transfer");
    require(transfer.sourceIdentity()==geometry.identity(),"foreign source VOF transfer owner");
    transfer.candidate().validateFor(targetMesh,transfer.candidate().identity(),clock);transfer.nativeOverlaps().validateFor(sourceMesh,targetMesh);
}

// Same fixed native-face WLS and corner limiter as charge advection, without
// its time prediction. The mean is located at the physical volume centroid.
struct RemapLinearCell
{
    long double x=0,y=0,hx=0,hy=0,a=0,b=0,c=0,bx=0,by=0,gx=0,gy=0;
    double rho=0,minimum=0,maximum=0;
    double slopeLimiter=1;
    bool needed=false,densityEligible=false,eligible=true;
};
Vector2D remapCentroid(const Geometry::PlanarCellCut::Rectangle& box,bool rz)
{
    const long double x=(static_cast<long double>(box.lower.x)+box.upper.x)/2;
    const long double y=(static_cast<long double>(box.lower.y)+box.upper.y)/2;
    const long double h=static_cast<long double>(box.upper.x)-box.lower.x;
    return {static_cast<double>(rz?x+h*(h/x/12):x),static_cast<double>(y)};
}
std::vector<RemapLinearCell> remapReconstruction(const QuadtreeMesh& mesh,const Vof& geometry,
    const Mesh::NativeLeafOverlap& overlap,const std::vector<double>& q)
{
    const auto& domain=geometry.domain();const bool rz=domain.axisymmetric();std::vector<RemapLinearCell> cells(q.size());
    for(std::size_t ci=0;ci<q.size();++ci){auto& c=cells[ci];
        if(domain.cellRowByIndex()[ci]<0){c.eligible=false;continue;}
        const auto& box=overlap.sourceRectangles()[ci];const auto point=remapCentroid(box,rz);
        c.x=point.x;c.y=point.y;c.hx=static_cast<long double>(box.upper.x)-box.lower.x;c.hy=static_cast<long double>(box.upper.y)-box.lower.y;
        c.rho=q[ci]/domain.cellVolumes()[domain.cellRowByIndex()[ci]];c.minimum=c.maximum=c.rho;
        // This is a declared finite-arithmetic eligibility branch, not a floor:
        // all such signed stocks retain the original exact binary64 partition.
        c.densityEligible=std::isfinite(c.rho)&&(c.rho!=0||q[ci]==0);
        c.eligible=c.densityEligible&&std::fpclassify(q[ci])!=FP_SUBNORMAL;
        for(const auto pi:overlap.bySource()[ci])c.needed=c.needed||overlap.pieces()[pi].kind==Mesh::NativeLeafOverlapPiece::Kind::Refined;
    }
    const auto sample=[&](int ci,int other){auto& c=cells[ci];const auto& n=cells[other];if(!c.needed||!c.eligible)return;
        if(!n.densityEligible){c.eligible=false;return;}const auto dx=(n.x-c.x)/c.hx,dy=(n.y-c.y)/c.hy;
        const auto square=dx*dx+dy*dy,difference=static_cast<long double>(n.rho)-c.rho;
        if(!std::isfinite(square)||square<=0||!std::isfinite(difference)){c.eligible=false;return;}
        const auto w=1/square;c.a+=w*dx*dx;c.b+=w*dx*dy;c.c+=w*dy*dy;c.bx+=w*dx*difference;c.by+=w*dy*difference;
        c.minimum=std::min(c.minimum,n.rho);c.maximum=std::max(c.maximum,n.rho);
    };
    for(const auto f:domain.faces(mesh))if(f.negativeCell>=0&&f.positiveCell>=0){sample(f.negativeCell,f.positiveCell);sample(f.positiveCell,f.negativeCell);}
    for(std::size_t ci=0;ci<cells.size();++ci){auto& c=cells[ci];if(!c.needed||!c.eligible)continue;
        const auto determinant=c.a*c.c-c.b*c.b;
        if(c.a>0&&c.c>0&&determinant>64*std::numeric_limits<long double>::epsilon()*c.a*c.c){
            c.gx=(c.c*c.bx-c.b*c.by)/determinant;c.gy=(c.a*c.by-c.b*c.bx)/determinant;
        }else{c.gx=c.a>0?c.bx/c.a:0;c.gy=c.c>0?c.by/c.c:0;}
        if(!std::isfinite(c.gx)||!std::isfinite(c.gy)){c.eligible=false;continue;}
        const auto& box=overlap.sourceRectangles()[ci];long double limiter=1;
        for(const long double x:{static_cast<long double>(box.lower.x),static_cast<long double>(box.upper.x)})
            for(const long double y:{static_cast<long double>(box.lower.y),static_cast<long double>(box.upper.y)}){
                const auto delta=c.gx*((x-c.x)/c.hx)+c.gy*((y-c.y)/c.hy);
                if(!std::isfinite(delta)){c.eligible=false;continue;}
                if(delta>0)limiter=std::min(limiter,(static_cast<long double>(c.maximum)-c.rho)/delta);
                else if(delta<0)limiter=std::min(limiter,(static_cast<long double>(c.minimum)-c.rho)/delta);
            }
        limiter=std::clamp(limiter,0.L,1.L);c.slopeLimiter=static_cast<double>(limiter);c.gx*=limiter;c.gy*=limiter;
    }
    return cells;
}
// Inward finite inventory bounds contain the exact P0 allocation, including
// its source-quantum discrepancy from the continuous density bound.
double remapBound(const FV::Dyadic& exact,bool lower)
{
    const double maximum=std::numeric_limits<double>::max();
    if(exact.compare(FV::Dyadic(maximum))>0)return maximum;
    if(exact.compare(FV::Dyadic(-maximum))<0)return -maximum;
    auto value=represented(exact);
    if(lower&&FV::Dyadic(value).compare(exact)<0)value=std::nextafter(value,maximum);
    if(!lower&&FV::Dyadic(value).compare(exact)>0)value=std::nextafter(value,-maximum);
    return value;
}
// Independent RNE need not cancel for arbitrary zero-mean signed corrections.
// Quantize their ACTIVITY toward zero on one
// common binary64 lattice, then contract the excess sign in one bounded pass.
// Each correction stays between zero and its bounded candidate; both signed sums
// become exactly equal. No charge is sent to another donor or added to Q.
void remapMeanZeroQuantize(std::vector<double>& values,int minimumQuantum)
{
    double largest=0;for(double value:values)largest=std::max(largest,std::abs(value));
    if(largest==0)return;int exponent=0;std::frexp(largest,&exponent);
    const int quantum=std::max({exponent-std::numeric_limits<double>::digits,-1074,minimumQuantum});
    FV::Dyadic positive,negative;
    for(auto& value:values){const auto units=std::floor(std::scalbn(std::abs(value),-quantum));
        value=std::copysign(std::scalbn(units,quantum),value);
        if(value>0)positive.add(value);else negative.add(-value);}
    auto excess=positive;excess.add(negative,-1);const int sign=excess.compare(FV::Dyadic{});excess=excess.absolute();
    for(auto& value:values){if(excess.compare(FV::Dyadic{})==0)break;if((value>0?1:value<0?-1:0)!=sign)continue;
        const FV::Dyadic magnitude(std::abs(value));const auto take=excess.compare(magnitude)<0?excess:magnitude;
        auto remaining=magnitude;remaining.add(take,-1);const auto rounded=FV::roundDouble(remaining);
        require(rounded.error.compare(FV::Dyadic{})==0,"AMR common-lattice contraction is not representable");
        value=std::copysign(rounded.value,value);excess.add(take,-1);
    }
    require(excess.compare(FV::Dyadic{})==0,"AMR signed activity could not balance");FV::Dyadic own;
    for(double value:values)own.add(value);
    require(own.compare(FV::Dyadic{})==0,"AMR correction quantization lost its exact zero mean");
}
bool remapLinearDonor(const RemapLinearCell& c,double source,const std::vector<std::size_t>& group,
    const Mesh::NativeLeafOverlap& overlap,const std::vector<FV::Dyadic>& weights,const std::vector<double>& baseline,
    bool rz,std::vector<double>& values,FV::Dyadic& error,double& limiter,double& quantizationMaximum,FV::Dyadic& quantizationL1)
{
    if(!c.eligible)return false;const auto count=group.size();std::vector<FV::Dyadic> delta(count);FV::Dyadic sum;
    std::size_t pivot=0;
    for(std::size_t k=0;k<count;++k){if(weights[k].compare(weights[pivot])>0)pivot=k;
        const auto point=remapCentroid(overlap.pieces()[group[k]].rectangle,rz);
        const auto change=c.gx*((static_cast<long double>(point.x)-c.x)/c.hx)+c.gy*((static_cast<long double>(point.y)-c.y)/c.hy);
        if(!std::isfinite(change))return false;delta[k]=weights[k].multiplied(change);sum.add(delta[k]);
    }
    // Correct only the finite evaluation of the zero-mean linear integral.
    // The largest geometric piece owns the deterministic exact remainder.
    delta[pivot].add(sum,-1);limiter=1;const FV::Dyadic zero;
    std::vector<double> lower(count),upper(count);
    for(std::size_t k=0;k<count;++k){lower[k]=std::min(baseline[k],remapBound(weights[k].multiplied(c.minimum),true));
        upper[k]=std::max(baseline[k],remapBound(weights[k].multiplied(c.maximum),false));
        const int sign=delta[k].compare(zero);FV::Dyadic budget;
        if(sign>0){budget=FV::Dyadic(upper[k]);budget.add(baseline[k],-1);}
        else if(sign<0){budget=FV::Dyadic(baseline[k]);budget.add(lower[k],-1);}
        if(sign)limiter=std::min(limiter,fctRatio(budget,delta[k].absolute()));
    }
    std::vector<double> correction(count);const FV::Dyadic maximum(std::numeric_limits<double>::max());
    for(std::size_t k=0;k<count;++k){delta[k]=delta[k].multiplied(limiter);if(delta[k].absolute().compare(maximum)>0)return false;
        correction[k]=represented(delta[k]);
        if(FV::Dyadic(correction[k]).absolute().compare(delta[k].absolute())>0)correction[k]=std::nextafter(correction[k],0.);
    }
    int exponent=0;std::frexp(std::abs(source),&exponent);
    const int sourceQuantum=source==0?-1074:std::max(exponent-std::numeric_limits<double>::digits,-1074);
    remapMeanZeroQuantize(correction,sourceQuantum);values.resize(count);FV::Dyadic own;
    // One finite-arithmetic qualification decision, never an iterative
    // accuracy repair. Every accepted child and the complete donor sum retain
    // the original exact binary64 conservation contract.
    for(std::size_t k=0;k<count;++k){FV::Dyadic exact(baseline[k]);exact.add(correction[k]);const auto rounded=FV::roundDouble(exact);
        if(rounded.error.compare(zero)!=0)return false;values[k]=rounded.value;own.add(values[k]);
        require(values[k]>=lower[k]&&values[k]<=upper[k],"linear AMR charge exceeds signed donor bounds");
        auto adjustment=delta[k];adjustment.add(correction[k],-1);quantizationL1.add(adjustment.absolute());
        quantizationMaximum=std::max(quantizationMaximum,represented(adjustment.absolute()));
    }
    own.add(source,-1);require(own.compare(zero)==0,"linear AMR donor lost its exact source charge");
    error=std::move(own);limiter*=c.slopeLimiter;
    return true;
}
}
void Detail::ElectricVofValidationBatch::validateTransport(const QuadtreeMesh& mesh,const ElectricVofCharge& source,const Transport& t)
{
    // Only the exact immutable Transport and the same original arguments have
    // completed this public validator earlier in this callback-free stack.
    if(!(mesh_==&mesh&&transport_&&t.sameInstance(*transport_))){
        t.validateFor(mesh,t.source(),t.advector(),t.liquidDensity(),t.gasDensity());
        rememberValidated(mesh,t.source().domain());transport_=&t;
    }
    require(t.source().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic,"charge transport requires actual LocalPlic transport");
    source.validateForInBatch(mesh,t.sourceTime(),*this);
    source.validateSelectedDomain(t.source().domain());
    require(source.timeLayer().position!=Clock::Position::Midpoint,"charge advection source is not accepted Initial/End");
    require(source.integratedCharge().size()==static_cast<std::size_t>(mesh.cellIndexNum()),"charge/native layout mismatch");
    const auto end=t.endTime();require(end.position==Clock::Position::End&&source.timeLayer().step<std::numeric_limits<std::int64_t>::max()&&end.step==source.timeLayer().step+1&&end.begin==source.timeLayer().time(),
        "actual transport has no next accepted-source interval");
    require(sameClock(t.midpoint().timeLayer(),Clock{end.step,end.begin,end.end,Clock::Position::Midpoint})&&
        sameClock(t.endpoint().timeLayer(),end),"actual transport stage clocks disagree");
}
struct ElectricVofChargeBoundary::Data
{
    ElectricVofCharge source;Transport transport;std::vector<double> density;bool closed=false;
    Data(const ElectricVofCharge& s,const Transport& t,const std::vector<double>& rho,bool c):source(s),transport(t),density(rho),closed(c){}
};
ElectricVofChargeBoundary::ElectricVofChargeBoundary(std::shared_ptr<const Data> data):data_(std::move(data)){}
ElectricVofChargeBoundary ElectricVofChargeBoundary::prepare(const QuadtreeMesh& mesh,const ElectricVofCharge& source,
    const Transport& t,const std::vector<double>& rho)
{
    Detail::ElectricVofValidationBatch batch;batch.validateTransport(mesh,source,t);const auto* fluid=t.advector().boundary();
    if(fluid)fluid->validateFor(mesh,t.source());const bool closed=!fluid||fluid->allClosed();
    require(rho.size()==t.faces().size(),"inlet charge density must cover actual native faces");
    // The actual F1 owner certifies endpoint normals and its bounded pure-phase
    // counterflow law. Constant-per-face reservoir rho consumes each actual
    // incoming gross volume, never the net volume after temporal cancellation.
    const auto& domain=t.source().domain();std::vector<unsigned> seen(rho.size());
    for(const auto& f:t.faces()){
        require(f.nativeRow>=0&&static_cast<std::size_t>(f.nativeRow)<rho.size()&&++seen[f.nativeRow]==1,"invalid/duplicate actual inlet charge row");
        const auto actual=domain.face(mesh,f.nativeRow);const auto value=rho[f.nativeRow];require(std::isfinite(value),"nonfinite signed inlet charge density");
        require(f.negativeCell==actual.negativeCell&&f.positiveCell==actual.positiveCell,"inlet charge incidence differs from native face");
        if(f.negativeCell<0&&f.positiveCell<0){
            requireInactiveFaceZero(f);require(value==0,"both-excluded charge face cannot carry reservoir density");
            if(fluid)require(fluid->rules().at(f.nativeRow)==Interface::GeometricVofAdvectorSlab::Boundary::Mode::Closed,
                "both-excluded charge face cannot be a reservoir");
            continue;
        }
        const bool exterior=(f.negativeCell>=0)!=(f.positiveCell>=0);const int sign=f.negativeCell>=0?1:-1;
        const auto outward=sign>0?f.full.total:-f.full.total;
        bool inflow=false;
        if(exterior&&!closed){const auto mode=fluid->rules().at(f.nativeRow);
            if(actual.area==0)require(domain.axisymmetric()&&actual.center.x==0&&mode==Interface::GeometricVofAdvectorSlab::Boundary::Mode::Closed&&
                f.full.total==0&&f.half.total==0&&value==0,"true axis cannot carry reservoir charge");
            if(mode==Interface::GeometricVofAdvectorSlab::Boundary::Mode::Closed)require(f.full.total==0&&f.half.total==0,"closed inlet charge face has nonzero actual T");
            else inflow=(sign>0?f.full.negativeTotal:f.full.positiveTotal)>0;
        }
        if(!inflow)require(value==0,"charge density must be zero on interior/closed/outflow/unused rows");
    }
    return ElectricVofChargeBoundary(std::make_shared<Data>(source,t,rho,closed));
}
const ElectricVofCharge& ElectricVofChargeBoundary::sourceCharge()const{require(bool(data_),"consumed inlet charge boundary");return data_->source;}
const Transport& ElectricVofChargeBoundary::transport()const{require(bool(data_),"consumed inlet charge boundary");return data_->transport;}
const std::vector<double>& ElectricVofChargeBoundary::inflowChargeDensities()const{require(bool(data_),"consumed inlet charge boundary");return data_->density;}
bool ElectricVofChargeBoundary::allClosed()const{require(bool(data_),"consumed inlet charge boundary");return data_->closed;}
bool ElectricVofChargeBoundary::sameInstance(const ElectricVofChargeBoundary& other)const noexcept{return data_&&data_==other.data_;}
void ElectricVofChargeBoundary::validateFor(const QuadtreeMesh& mesh,const ElectricVofCharge& source,const Transport& t)const
{
    Detail::ElectricVofValidationBatch batch;validateForInBatch(mesh,source,t,batch);
}
void ElectricVofChargeBoundary::validateForInBatch(const QuadtreeMesh& mesh,const ElectricVofCharge& source,const Transport& t,
    Detail::ElectricVofValidationBatch& batch)const
{
    require(bool(data_),"consumed inlet charge boundary");batch.validateTransport(mesh,source,t);
    require(source.sameInstance(data_->source)&&t.sameInstance(data_->transport),"foreign source/actual transport inlet charge owner");
    require(sameClock(t.sourceTime(),data_->transport.sourceTime())&&sameClock(t.endTime(),data_->transport.endTime()),"complete inlet charge boundary interval mismatch");
}
struct ElectricVofAdvectedCharge::Data
{
    ElectricVofCharge source;Transport transport;Vof geometry;Clock clock;
    std::shared_ptr<const ElectricVofChargeBoundary> boundary;
    const ElectricVofControls controls;
    std::vector<double> q;std::vector<ElectricVofChargeFace> faces;ElectricVofChargeAdvectionReport report;
    Data(const ElectricVofCharge& s,const Transport& t,bool half,Stage result,const ElectricVofControls& c,
        const ElectricVofChargeBoundary* b=nullptr):source(s),transport(t),geometry(half?t.midpoint():t.endpoint()),clock(geometry.timeLayer()),
        boundary(b?std::make_shared<ElectricVofChargeBoundary>(*b):nullptr),controls(c),
        q(std::move(result.q)),faces(std::move(result.faces)),report(result.report){}
};
ElectricVofAdvectedCharge::ElectricVofAdvectedCharge(std::shared_ptr<const Data> d):data_(std::move(d)){}
const ElectricVofCharge& ElectricVofAdvectedCharge::sourceCharge()const{require(bool(data_),"consumed charge advection receipt");return data_->source;}
const Transport& ElectricVofAdvectedCharge::transport()const{require(bool(data_),"consumed charge advection receipt");return data_->transport;}
const Vof& ElectricVofAdvectedCharge::geometry()const{require(bool(data_),"consumed charge advection receipt");return data_->geometry;}
Clock ElectricVofAdvectedCharge::timeLayer()const{require(bool(data_),"consumed charge advection receipt");return data_->clock;}
const std::vector<double>& ElectricVofAdvectedCharge::integratedCharge()const{require(bool(data_),"consumed charge advection receipt");return data_->q;}
const std::vector<ElectricVofChargeFace>& ElectricVofAdvectedCharge::faces()const{require(bool(data_),"consumed charge advection receipt");return data_->faces;}
const ElectricVofChargeAdvectionReport& ElectricVofAdvectedCharge::report()const{require(bool(data_),"consumed charge advection receipt");return data_->report;}
const ElectricVofChargeBoundary* ElectricVofAdvectedCharge::boundary()const{require(bool(data_),"consumed charge advection receipt");return data_->boundary.get();}
bool ElectricVofAdvectedCharge::sameInstance(const ElectricVofAdvectedCharge& other)const noexcept{return data_&&data_==other.data_;}
void ElectricVofAdvectedCharge::validateFor(const QuadtreeMesh& mesh,const ElectricVofCharge& source,const Transport& t,Clock stage)const
{
    Detail::ElectricVofValidationBatch batch;validateForInBatch(mesh,source,t,stage,batch);
}
void ElectricVofAdvectedCharge::validateForInBatch(const QuadtreeMesh& mesh,const ElectricVofCharge& source,const Transport& t,Clock stage,
    Detail::ElectricVofValidationBatch& batch)const
{
    require(bool(data_),"consumed charge advection receipt");batch.validateTransport(mesh,source,t);
    require(source.sameInstance(data_->source)&&t.sameInstance(data_->transport),"foreign charge/transport receipt owner");
    if(data_->boundary)data_->boundary->validateForInBatch(mesh,source,t,batch);
    require(sameClock(stage,data_->clock),"complete advected charge stage clock mismatch");
    const auto& actual=stage.position==Clock::Position::Midpoint?t.midpoint():t.endpoint();
    require(actual.identity()==data_->geometry.identity(),"foreign actual stage VOF geometry");geometry().validateFor(mesh,actual.identity(),stage);
}
ElectricVofChargeTransport::ElectricVofChargeTransport(ElectricVofAdvectedCharge half,ElectricVofAdvectedCharge full):midpoint_(std::move(half)),endpoint_(std::move(full)){}
ElectricVofChargeTransport ElectricVofChargeTransport::prepareClosed(const QuadtreeMesh& mesh,const ElectricVofCharge& source,const Transport& t,const ElectricVofControls& controls)
{
    controlsValid(controls);Detail::ElectricVofValidationBatch batch;batch.validateTransport(mesh,source,t);
    auto half=advance(mesh,source,t,true,controls);auto full=advance(mesh,source,t,false,controls);
    return ElectricVofChargeTransport(ElectricVofAdvectedCharge(std::make_shared<ElectricVofAdvectedCharge::Data>(source,t,true,std::move(half),controls)),
        ElectricVofAdvectedCharge(std::make_shared<ElectricVofAdvectedCharge::Data>(source,t,false,std::move(full),controls)));
}
ElectricVofChargeTransport ElectricVofChargeTransport::prepareBoundary(const QuadtreeMesh& mesh,const ElectricVofCharge& source,
    const Transport& t,const ElectricVofChargeBoundary& boundary,const ElectricVofControls& controls)
{
    controlsValid(controls);boundary.validateFor(mesh,source,t);
    if(boundary.allClosed())return prepareClosed(mesh,source,t,controls);
    auto half=advance(mesh,source,t,true,controls,&boundary);auto full=advance(mesh,source,t,false,controls,&boundary);
    return ElectricVofChargeTransport(ElectricVofAdvectedCharge(std::make_shared<ElectricVofAdvectedCharge::Data>(source,t,true,std::move(half),controls,&boundary)),
        ElectricVofAdvectedCharge(std::make_shared<ElectricVofAdvectedCharge::Data>(source,t,false,std::move(full),controls,&boundary)));
}
ElectricVofAdvectedCharge ElectricVofChargeTransport::prepareMidpointClosed(const QuadtreeMesh& mesh,const ElectricVofCharge& source,
    const Transport& t,const ElectricVofControls& controls)
{
    controlsValid(controls);Detail::ElectricVofValidationBatch batch;batch.validateTransport(mesh,source,t);
    auto half=advance(mesh,source,t,true,controls);
    return ElectricVofAdvectedCharge(std::make_shared<ElectricVofAdvectedCharge::Data>(source,t,true,std::move(half),controls));
}
ElectricVofAdvectedCharge ElectricVofChargeTransport::prepareMidpointBoundary(const QuadtreeMesh& mesh,const ElectricVofCharge& source,
    const Transport& t,const ElectricVofChargeBoundary& boundary,const ElectricVofControls& controls)
{
    controlsValid(controls);boundary.validateFor(mesh,source,t);
    if(boundary.allClosed())return prepareMidpointClosed(mesh,source,t,controls);
    auto half=advance(mesh,source,t,true,controls,&boundary);
    return ElectricVofAdvectedCharge(std::make_shared<ElectricVofAdvectedCharge::Data>(source,t,true,std::move(half),controls,&boundary));
}
ElectricVofChargeTransport ElectricVofChargeTransport::completeEndpoint(const QuadtreeMesh& mesh,const ElectricVofAdvectedCharge& half)
{
    require(half.timeLayer().position==Clock::Position::Midpoint,"endpoint completion requires an actual Midpoint charge receipt");
    const auto& source=half.sourceCharge();const auto& t=half.transport();const auto& controls=half.data_->controls;
    controlsValid(controls);half.validateFor(mesh,source,t,t.midpoint().timeLayer());
    auto full=advance(mesh,source,t,false,controls,half.boundary());
    return ElectricVofChargeTransport(half,ElectricVofAdvectedCharge(std::make_shared<ElectricVofAdvectedCharge::Data>(
        source,t,false,std::move(full),controls,half.boundary())));
}
const ElectricVofCharge& ElectricVofChargeTransport::sourceCharge()const{return midpoint_.sourceCharge();}
const Transport& ElectricVofChargeTransport::transport()const{return midpoint_.transport();}
struct ElectricVofChargeMeshTransfer::Data
{
    Vof geometry;ElectricVofCharge source,candidate;Transfer transfer;ElectricVofChargeRemapReport report;
    Data(const Vof& g,const ElectricVofCharge& s,const Transfer& t,ElectricVofCharge q,ElectricVofChargeRemapReport r):
        geometry(g),source(s),candidate(std::move(q)),transfer(t),report(r){}
};
ElectricVofChargeMeshTransfer::ElectricVofChargeMeshTransfer(std::shared_ptr<const Data> d):data_(std::move(d)){}
ElectricVofChargeMeshTransfer ElectricVofChargeMeshTransfer::prepare(const QuadtreeMesh& sourceMesh,const Vof& geometry,const ElectricVofCharge& source,
    const QuadtreeMesh& targetMesh,const Transfer& transfer,const ElectricVofControls& controls)
{
    const auto started=Timer::now();controlsValid(controls);remapValid(sourceMesh,geometry,source,targetMesh,transfer);
    source.validateSelectedDomain(geometry.domain());
    const auto& sourceActive=geometry.domain().cellRowByIndex();const auto& targetActive=transfer.candidate().domain().cellRowByIndex();
    const auto& overlap=transfer.nativeOverlaps();const auto& pieces=overlap.pieces();const auto& old=source.integratedCharge();
    require(overlap.bySource().size()==old.size()&&sourceActive.size()==old.size()&&
        overlap.byTarget().size()==static_cast<std::size_t>(targetMesh.cellIndexNum())&&targetActive.size()==overlap.byTarget().size(),
        "actual charge overlap layouts disagree");
    std::vector<double> contributions(pieces.size());std::vector<unsigned> seen(pieces.size());ElectricVofChargeRemapReport report;report.pieces=pieces.size();
    const bool rz=geometry.domain().axisymmetric();
    const bool linear=controls.chargeAdvection==ElectricVofChargeAdvectionScheme::LimitedLinearFct&&
        std::any_of(pieces.begin(),pieces.end(),[](const auto& p){return p.kind==Mesh::NativeLeafOverlapPiece::Kind::Refined;});
    const auto reconstruction=linear?remapReconstruction(sourceMesh,geometry,overlap,old):std::vector<RemapLinearCell>{};
    report.limitedLinearRemap=linear;FV::Dyadic prolongationRounding,quantizationL1;
    for(std::size_t si=0;si<old.size();++si){auto group=overlap.bySource()[si];
        if(sourceActive[si]<0){require(group.empty()&&old[si]==0,"inactive source charge must have zero stock and no overlap");continue;}
        require(!group.empty(),"source charge lacks actual overlap coverage");
        std::sort(group.begin(),group.end(),[&](std::size_t a,std::size_t b){const auto& x=pieces.at(a).rectangle;const auto& y=pieces.at(b).rectangle;
            return std::tie(x.lower.x,x.lower.y,x.upper.x,x.upper.y)<std::tie(y.lower.x,y.lower.y,y.upper.x,y.upper.y);});
        std::vector<FV::Dyadic> weights;FV::Dyadic whole;bool refined=false;
        for(const auto pi:group){require(pi<pieces.size()&&++seen[pi]==1&&pieces[pi].sourceCell==static_cast<int>(si),"invalid/duplicate source charge overlap");
            const auto w=FV::rectangleCapacityExact(pieces[pi].rectangle,geometry.domain().axisymmetric());weights.push_back(w);whole.add(w);
            refined=refined||pieces[pi].kind==Mesh::NativeLeafOverlapPiece::Kind::Refined;}
        require(whole.compare(FV::rectangleCapacityExact(overlap.sourceRectangles().at(si),geometry.domain().axisymmetric()))==0,
            "actual charge rectangles do not exactly partition source geometry");
        if(refined)++report.refinedDonors;
        if(group.size()==1&&pieces[group[0]].kind==Mesh::NativeLeafOverlapPiece::Kind::Unchanged){contributions[group[0]]=old[si];continue;}
        const auto allocation=FV::partitionDouble(std::abs(old[si]),weights);require(allocation.sumError.compare(FV::Dyadic{})==0,"AMR binary64 source partition lost charge");
        std::vector<double> baseline(group.size());for(std::size_t k=0;k<group.size();++k)baseline[k]=std::copysign(allocation.values[k],old[si]);
        FV::Dyadic ownError,donorQuantizationL1;std::vector<double> high;double limiter=1,quantizationMaximum=0;bool highOrder=false;
        if(linear&&refined){highOrder=remapLinearDonor(reconstruction[si],old[si],group,overlap,weights,baseline,rz,
                high,ownError,limiter,quantizationMaximum,donorQuantizationL1);
            if(highOrder){++report.highOrderRefinedDonors;if(limiter<1)++report.limitedRefinedDonors;
                quantizationL1.add(donorQuantizationL1);report.maximumLinearQuantizationAdjustment=
                    std::max(report.maximumLinearQuantizationAdjustment,quantizationMaximum);
                if(old[si]==0&&quantizationMaximum>0){++report.zeroNetQuantizedDonors;report.maximumZeroNetQuantizationAdjustment=
                    std::max(report.maximumZeroNetQuantizationAdjustment,quantizationMaximum);}
                report.minimumDonorLimiter=std::min(report.minimumDonorLimiter,limiter);prolongationRounding.add(ownError);
                if(old[si]!=0)report.maximumDonorRelativeChargeResidual=std::max(report.maximumDonorRelativeChargeResidual,
                    represented(ownError.absolute())/std::abs(old[si]));
            }else ++report.finiteEligibilityP0Donors;}
        FV::Dyadic own;for(std::size_t k=0;k<group.size();++k){contributions[group[k]]=highOrder?high[k]:baseline[k];own.add(contributions[group[k]]);
            if(old[si]!=0&&contributions[group[k]]==0)++report.quantizedZeroPieces;}
        own.add(old[si],-1);require(own.compare(ownError)==0,"signed AMR source charge partition mismatch");
    }
    std::vector<double> next(overlap.byTarget().size());FV::Dyadic rounding=prolongationRounding,restrictionRounding,before,after;Ledger global;
    for(std::size_t ti=0;ti<next.size();++ti){const auto& group=overlap.byTarget()[ti];
        if(targetActive[ti]<0){require(group.empty(),"inactive target charge has physical overlap");continue;}
        require(!group.empty(),"target charge lacks actual overlap coverage");
        FV::Dyadic exact;bool unchanged=group.size()==1,coarsened=false;
        for(const auto pi:group){require(pi<pieces.size()&&seen[pi]==1&&pieces[pi].targetCell==static_cast<int>(ti),"invalid target charge overlap");exact.add(contributions[pi]);
            unchanged=unchanged&&pieces[pi].kind==Mesh::NativeLeafOverlapPiece::Kind::Unchanged;coarsened=coarsened||pieces[pi].kind==Mesh::NativeLeafOverlapPiece::Kind::Coarsened;}
        FV::RoundedDouble rounded;if(unchanged){rounded.value=contributions[group[0]];++report.unchangedTargets;}else rounded=FV::roundDouble(exact);
        if(coarsened)++report.coarsenedTargets;require(std::isfinite(rounded.value),"unrepresentable target AMR charge");next[ti]=rounded.value;
        rounding.add(rounded.error);restrictionRounding.add(rounded.error);
        report.maximumCellRoundingError=std::max(report.maximumCellRoundingError,represented(rounded.error.absolute()));
        after.add(rounded.value);global.term(rounded.value,1,controls.relativeConservationTolerance);
    }
    for(const auto q:old){before.add(q);global.term(q,-1,controls.relativeConservationTolerance);}
    require(global.residual.compare(rounding)==0,"AMR charge global difference is not prolongation plus restriction rounding");global.validate();
    report.sourceTotalCharge=represented(before);report.candidateTotalCharge=represented(after);report.chargeResidual=represented(global.residual);
    report.relativeChargeResidual=global.relative();report.roundingResidual=represented(rounding);
    report.prolongationRoundingResidual=represented(prolongationRounding);report.restrictionRoundingResidual=represented(restrictionRounding);
    report.linearQuantizationL1Adjustment=represented(quantizationL1);
    report.totalSeconds=elapsed(started);
    auto candidate=ElectricVofCharge::fromPrepared(targetMesh,transfer.candidate(),source.timeLayer(),std::move(next));
    return ElectricVofChargeMeshTransfer(std::make_shared<Data>(geometry,source,transfer,std::move(candidate),report));
}
const Vof& ElectricVofChargeMeshTransfer::sourceGeometry()const{require(bool(data_),"consumed charge AMR receipt");return data_->geometry;}
const ElectricVofCharge& ElectricVofChargeMeshTransfer::sourceCharge()const{require(bool(data_),"consumed charge AMR receipt");return data_->source;}
const Transfer& ElectricVofChargeMeshTransfer::interfaceTransfer()const{require(bool(data_),"consumed charge AMR receipt");return data_->transfer;}
const ElectricVofCharge& ElectricVofChargeMeshTransfer::chargeCandidate()const{require(bool(data_),"consumed charge AMR receipt");return data_->candidate;}
const ElectricVofChargeRemapReport& ElectricVofChargeMeshTransfer::report()const{require(bool(data_),"consumed charge AMR receipt");return data_->report;}
void ElectricVofChargeMeshTransfer::validateFor(const QuadtreeMesh& sourceMesh,const Vof& geometry,const ElectricVofCharge& source,
    const QuadtreeMesh& targetMesh,const Transfer& transfer)const
{
    require(bool(data_),"consumed charge AMR receipt");remapValid(sourceMesh,geometry,source,targetMesh,transfer);
    source.validateSelectedDomain(geometry.domain());
    require(source.sameInstance(data_->source)&&geometry.identity()==data_->geometry.identity()&&
        transfer.sourceIdentity()==data_->transfer.sourceIdentity()&&transfer.candidate().identity()==data_->transfer.candidate().identity(),"foreign charge AMR receipt owner");
    data_->candidate.validateFor(targetMesh,source.timeLayer());
}
}
