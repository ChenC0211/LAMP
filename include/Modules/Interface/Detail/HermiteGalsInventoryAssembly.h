#pragma once

#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Modules/Interface/Detail/HermiteGalsVolumeRoundoff.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace LAMP::Interface::Detail
{
// One private numerical owner for both the explicit old rate API and the
// interval API. unitVolume is old dt for rates, exactly one for integrated
// volumes. Keep the old multiplication/subtraction order in the rate branch.
struct InventoryFace
{
    int negativeCell=-1,positiveCell=-1;
    long double totalPositive=0,totalNegative=0,liquidPositive=0,liquidNegative=0;
    long double rawPositive=0,rawNegative=0,liquidSigned=0,liquidVolume=0;
};
struct InventoryFaceBounds
{
    long double lowerPositive=0,upperPositive=0,lowerNegative=0,upperNegative=0;
};
struct InventorySum
{
    long double value=0,correction=0;
    void add(long double x){const auto a=x-correction,b=value+a;correction=(b-value)-a;value=b;}
};
class HermiteInventoryAssembly
{
public:
    std::vector<long double> capacity,before;
    HermiteInventoryAssembly(std::size_t count,const std::vector<int>& cells,const std::vector<double>& measures,
        const std::vector<HermiteGalsCellVolume>& source,HermiteGalsIdentity identity)
        :capacity(count,0),before(count,0),cells_(cells),id_(identity)
    {
        require(cells.size()==measures.size() && cells.size()==source.size(),"advance: previous inventory layout mismatch");
        for(std::size_t i=0;i<cells.size();++i){const int cell=cells[i];capacity[cell]=measures[i];
            require(source[i].cell==cell,"advance: previous inventory cell mismatch");before[cell]=source[i].liquid;
            require(std::isfinite(before[cell]) && before[cell]>=0 && before[cell]<=capacity[cell],"advance: old inventory outside capacity");
            totalCapacity_.add(capacity[cell]);totalBefore_.add(before[cell]);}
        require(std::isfinite(totalCapacity_.value) && totalCapacity_.value>0 && std::isfinite(totalBefore_.value),
            "advance: nonfinite capacity/inventory total");
    }
    HermiteGalsInventoryAssemblyReport finish(std::vector<InventoryFace>& faces,long double unitVolume,
        double divergenceTolerance,double conservationTolerance,std::vector<HermiteGalsCellVolume>& volumes,
        const std::vector<InventoryFaceBounds>* bounds=nullptr) const
    {
        const auto count=capacity.size();const auto& cells=cells_;
        const auto roundoffBudget=hermiteInventoryRoundoffBudget;
        if(bounds){
            require(unitVolume==1 && bounds->size()==faces.size(),"advance: invalid integrated certificate layout/units");
            for(std::size_t k=0;k<faces.size();++k){const auto& f=faces[k];const auto& b=(*bounds)[k];
                const auto valid=[](long double lo,long double hi,long double value,long double raw,long double total){
                    return std::isfinite(lo) && std::isfinite(hi) && lo>=0 && lo<=value && value<=hi
                        && lo<=raw && raw<=hi && hi<=total;};
                require(valid(b.lowerPositive,b.upperPositive,f.liquidPositive,f.rawPositive,f.totalPositive)
                    && valid(b.lowerNegative,b.upperNegative,f.liquidNegative,f.rawNegative,f.totalNegative),
                    "advance: invalid certified phase bounds");
            }
        }
        std::vector<long double> exchange(count,0),divergence(count,0),outflow(count,0);
        HermiteGalsInventoryAssemblyReport report;report.activeFaces=faces.size();
        for(const auto& f:faces){
            const auto totalVolume=unitVolume*(f.totalPositive-f.totalNegative);
            if(f.negativeCell>=0){divergence[f.negativeCell]+=totalVolume;outflow[f.negativeCell]+=unitVolume*f.totalPositive;}
            if(f.positiveCell>=0){divergence[f.positiveCell]-=totalVolume;outflow[f.positiveCell]+=unitVolume*f.totalNegative;}
            report.internalFaces+=f.negativeCell>=0 && f.positiveCell>=0;
            report.counterflowFaces+=f.totalPositive>0 && f.totalNegative>0;
        }
        long double maxCfl=0,maxDivergence=0;
        for(int cell:cells){
            require(std::isfinite(outflow[cell]) && outflow[cell]<=capacity[cell],"advance: gross outflow CFL exceeds one cell="+std::to_string(cell));
            const auto residual=std::abs(divergence[cell])/capacity[cell];
            require(std::isfinite(residual) && residual<=divergenceTolerance,"advance: total flux is not divergence free cell="+std::to_string(cell));
            maxCfl=std::max(maxCfl,outflow[cell]/capacity[cell]);maxDivergence=std::max(maxDivergence,residual);}
        const auto assemble=[&]{std::fill(exchange.begin(),exchange.end(),0);
            for(auto& f:faces){f.liquidSigned=f.liquidPositive-f.liquidNegative;f.liquidVolume=unitVolume*f.liquidSigned;
                require(std::isfinite(f.liquidVolume),"advance: nonfinite shared phase volume");
                if(f.negativeCell>=0)exchange[f.negativeCell]-=f.liquidVolume;
                if(f.positiveCell>=0)exchange[f.positiveCell]+=f.liquidVolume;}};
        assemble();std::size_t sweeps=0;
        for(;;){bool bounded=true;
            for(int cell:cells){const auto value=before[cell]+exchange[cell];
                require(std::isfinite(value) && value>=-roundoffBudget*capacity[cell] && value<=capacity[cell]*(1+roundoffBudget),
                    "advance: inventory violation exceeds roundoff budget cell="+std::to_string(cell));
                bounded=bounded && value>=0 && value<=capacity[cell];}
            if(bounded)break;
            require(sweeps<128,"advance: roundoff flux repair work limit exhausted");++sweeps;
            // A just-filled cell can be one ulp above capacity even when its
            // interface-side outgoing liquid is below total. Shrinking only
            // incoming liquid pushes this defect upstream through full cells
            // until a prescribed boundary cannot change. Use the AVAILABLE
            // internal gross component and the neighbour's inventory room.
            // Every increment remains within the original per-face budget;
            // neither total nor boundary flux is altered.
            bool increased=false;
            for(int cell:cells) {
                const auto value=before[cell]+exchange[cell];
                const bool excess=value>capacity[cell],deficit=value<0;
                if(!excess && !deficit)continue;
                for(auto& f:faces) {
                    if(f.negativeCell<0 || f.positiveCell<0 || (f.negativeCell!=cell && f.positiveCell!=cell))continue;
                    const int neighbour=f.negativeCell==cell?f.positiveCell:f.negativeCell;
                    const auto current=before[cell]+exchange[cell],other=before[neighbour]+exchange[neighbour];
                    const auto defect=excess?current-capacity[cell]:-current;
                    if(defect<=0)break;
                    const auto room=excess?capacity[neighbour]-other:other;
                    if(room<=0)continue;
                    const bool positive=excess?(f.negativeCell==cell):(f.positiveCell==cell);
                    auto& liquid=positive?f.liquidPositive:f.liquidNegative;
                    const auto total=positive?f.totalPositive:f.totalNegative;
                    const auto maximum=bounds?(positive?(*bounds)[static_cast<std::size_t>(&f-faces.data())].upperPositive
                        :(*bounds)[static_cast<std::size_t>(&f-faces.data())].upperNegative):total;
                    const auto raw=positive?f.rawPositive:f.rawNegative;
                    const auto budget=roundoffBudget*std::min(capacity[cell],capacity[neighbour]);
                    const auto margin=std::min(unitVolume*(maximum-liquid),budget-unitVolume*std::abs(liquid-raw));
                    const auto amount=std::min({defect,room,margin});
                    if(amount<=0)continue;
                    auto next=liquid+amount/unitVolume;
                    if(next==liquid)next=std::nextafter(liquid,std::numeric_limits<long double>::infinity());
                    const auto change=unitVolume*(next-liquid);
                    if(!std::isfinite(next) || next>maximum || change<=0 || change>room ||
                        unitVolume*std::abs(next-raw)>budget)continue;
                    liquid=next;increased=true;assemble();
                }
            }
            if(increased)continue;
            if(bounds){
                // Consume only the reducible part above the certificate's
                // lower bound. The unqualified rate/integral expressions
                // below remain unchanged for their existing consumers.
                std::vector<long double> incoming(count,0),outgoing(count,0);
                for(std::size_t k=0;k<faces.size();++k){const auto& f=faces[k];const auto& b=(*bounds)[k];
                    if(f.negativeCell<0 || f.positiveCell<0)continue;
                    const auto positive=f.liquidPositive-b.lowerPositive,negative=f.liquidNegative-b.lowerNegative;
                    outgoing[f.negativeCell]+=positive;incoming[f.positiveCell]+=positive;
                    outgoing[f.positiveCell]+=negative;incoming[f.negativeCell]+=negative;
                }
                const auto fraction=[&](long double defect,long double available){if(defect<=0)return 1.L;
                    require(available>=defect && available>0,"advance: no certified flux capacity for roundoff repair");
                    const auto f=1-defect/available;return f<=0?0.L:std::nextafter(f,0.L);};
                for(int cell:cells){const auto value=before[cell]+exchange[cell];
                    incoming[cell]=fraction(value-capacity[cell],incoming[cell]);outgoing[cell]=fraction(-value,outgoing[cell]);}
                const auto reduce=[](long double liquid,long double lower,long double keep){
                    if(keep==1 || liquid==lower)return liquid;
                    auto next=std::max(lower,std::min(liquid,lower+(liquid-lower)*keep));
                    if(next==liquid)next=std::nextafter(liquid,lower);
                    return next;
                };
                bool changed=false;
                for(std::size_t k=0;k<faces.size();++k){auto& f=faces[k];const auto& b=(*bounds)[k];
                    if(f.negativeCell<0 || f.positiveCell<0)continue;
                    const auto positive=reduce(f.liquidPositive,b.lowerPositive,std::min(outgoing[f.negativeCell],incoming[f.positiveCell])),
                        negative=reduce(f.liquidNegative,b.lowerNegative,std::min(outgoing[f.positiveCell],incoming[f.negativeCell]));
                    require(positive>=b.lowerPositive && positive<=b.upperPositive && negative>=b.lowerNegative && negative<=b.upperNegative,
                        "advance: certified reduction escaped phase bounds");
                    const auto correction=std::max(std::abs(f.rawPositive-positive),std::abs(f.rawNegative-negative));
                    require(correction<=roundoffBudget*std::min(capacity[f.negativeCell],capacity[f.positiveCell]),
                        "advance: certified flux correction exceeds roundoff budget");
                    changed=changed || positive!=f.liquidPositive || negative!=f.liquidNegative;
                    f.liquidPositive=positive;f.liquidNegative=negative;
                }
                require(changed,"advance: certified roundoff repair has no representable progress");
                assemble();continue;
            }
            std::vector<long double> incoming(count,0),outgoing(count,0);
            for(const auto& f:faces)if(f.negativeCell>=0 && f.positiveCell>=0){
                outgoing[f.negativeCell]+=unitVolume*f.liquidPositive;incoming[f.positiveCell]+=unitVolume*f.liquidPositive;
                outgoing[f.positiveCell]+=unitVolume*f.liquidNegative;incoming[f.negativeCell]+=unitVolume*f.liquidNegative;}
            const auto reduction=[&](long double defect,long double available){if(defect<=0)return 1.L;
                require(available>=defect && available>0,"advance: no interior flux capacity for roundoff repair");
                const auto fraction=1-defect/available;return fraction<=0?0.L:std::nextafter(fraction,0.L);};
            for(int cell:cells){const auto value=before[cell]+exchange[cell];
                incoming[cell]=reduction(value-capacity[cell],incoming[cell]);outgoing[cell]=reduction(-value,outgoing[cell]);}
            for(auto& f:faces)if(f.negativeCell>=0 && f.positiveCell>=0){
                f.liquidPositive*=std::min(outgoing[f.negativeCell],incoming[f.positiveCell]);
                f.liquidNegative*=std::min(outgoing[f.positiveCell],incoming[f.negativeCell]);
                const auto correction=unitVolume*std::max(std::abs(f.rawPositive-f.liquidPositive),std::abs(f.rawNegative-f.liquidNegative));
                require(correction<=roundoffBudget*std::min(capacity[f.negativeCell],capacity[f.positiveCell]),"advance: shared flux correction exceeds roundoff budget");}
            assemble();
        }
        InventorySum after,boundary;std::size_t repaired=0;long double maxCorrection=0;
        for(const auto& f:faces){
            if(f.negativeCell<0)boundary.add(-f.liquidVolume);else if(f.positiveCell<0)boundary.add(f.liquidVolume);
            if(f.rawPositive!=f.liquidPositive || f.rawNegative!=f.liquidNegative){++repaired;
                maxCorrection=std::max(maxCorrection,unitVolume*std::max(std::abs(f.rawPositive-f.liquidPositive),std::abs(f.rawNegative-f.liquidNegative))
                    /std::min(capacity[f.negativeCell],capacity[f.positiveCell]));}}
        for(int cell:cells){const double value=static_cast<double>(before[cell]+exchange[cell]);
            require(std::isfinite(value) && value>=0 && value<=capacity[cell],"advance: final inventory outside strict bounds");
            volumes.push_back({cell,value});after.add(value);}
        const auto closure=after.value-totalBefore_.value+boundary.value;
        require(std::isfinite(closure) && std::abs(closure)<=conservationTolerance*totalCapacity_.value,"advance: global phase-volume ledger does not close");
        report.roundoffRepairSweeps=sweeps;report.repairedFaces=repaired;
        report.liquidBefore=totalBefore_.value;report.liquidAfter=after.value;report.boundaryOutward=boundary.value;report.closure=closure;
        report.maxOutflowCfl=maxCfl;report.maxRelativeDivergence=maxDivergence;report.maxRelativeFluxCorrection=maxCorrection;
        return report;
    }
private:
    const std::vector<int>& cells_;
    HermiteGalsIdentity id_;
    InventorySum totalCapacity_,totalBefore_;
    void require(bool ok,const std::string& message) const
    {if(!ok)throw std::invalid_argument("[Interface/HermiteGalsPhaseFlux geometry="+std::to_string(id_.value())+"] "+message);}
};
}
