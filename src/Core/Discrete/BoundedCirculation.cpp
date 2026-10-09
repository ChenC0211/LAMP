#include "Core/Discrete/BoundedCirculation.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace Quantum::Discrete
{
namespace
{
namespace Exact=Quantum::Geometry::Detail::ExactExpansion;
using Number=Exact::Expansion;
void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[BoundedCirculation] ")+why);}
struct Arc {std::size_t from,to;Number lower,upper,flow;};
struct Link {std::size_t arc,to;bool forward;};
}
BoundedCirculationResult solveBoundedCirculation(std::size_t nodes,
    const std::vector<BoundedCirculationArc>& input)
{
    Exact::validateEnvironment();
    require(nodes>0 && nodes<=std::numeric_limits<std::size_t>::max()-2,"invalid node count");
    const auto source=nodes,sink=nodes+1,missing=std::numeric_limits<std::size_t>::max();
    std::vector<Arc> arcs;arcs.reserve(input.size());
    std::vector<std::vector<Link>> graph(nodes+2);
    std::vector<Number> balance(nodes);
    const auto append=[&](std::size_t from,std::size_t to,Number lo,Number hi,Number flow){
        const auto index=arcs.size();arcs.push_back({from,to,std::move(lo),std::move(hi),std::move(flow)});
        graph[from].push_back({index,to,true});graph[to].push_back({index,from,false});
    };
    for(const auto& a:input){
        require(a.from<nodes && a.to<nodes && a.from!=a.to,"invalid arc endpoints");
        require(std::isfinite(a.lower) && std::isfinite(a.upper) && std::isfinite(a.initial)
            && a.lower>=0 && a.lower<=a.initial && a.initial<=a.upper,"invalid bounded arc");
        const auto flow=Exact::fromDouble(a.initial);
        balance[a.from]=Exact::subtract(balance[a.from],flow);balance[a.to]=Exact::add(balance[a.to],flow);
        append(a.from,a.to,Exact::fromDouble(a.lower),Exact::fromDouble(a.upper),flow);
    }
    Number required;
    for(std::size_t node=0;node<nodes;++node){const int sign=Exact::sign(balance[node]);
        if(sign>0){append(source,node,{},balance[node],{});required=Exact::add(required,balance[node]);}
        else if(sign<0)append(node,sink,{},Exact::subtract({},balance[node]),{});
    }
    // A conservative Edmonds-Karp bound, with checked integer arithmetic.
    const auto max=std::numeric_limits<std::size_t>::max();const auto edges=arcs.size();
    require(edges==0 || (edges<=max/edges && nodes+2<=max/(edges*edges)),"network work bound overflow");
    const auto workLimit=edges? (nodes+2)*edges*edges:0;
    BoundedCirculationResult result;Number sent;
    const auto residual=[&](const Link& link){const auto& a=arcs[link.arc];
        return link.forward?Exact::subtract(a.upper,a.flow):Exact::subtract(a.flow,a.lower);};
    while(Exact::sign(Exact::subtract(required,sent))>0){
        std::vector<std::size_t> parent(nodes+2,missing);std::vector<Link> path(nodes+2);std::vector<std::size_t> queue{source};
        parent[source]=source;
        for(std::size_t head=0;head<queue.size() && parent[sink]==missing;++head){const auto node=queue[head];
            for(const auto link:graph[node])if(parent[link.to]==missing && Exact::sign(residual(link))>0){
                parent[link.to]=node;path[link.to]=link;queue.push_back(link.to);if(link.to==sink)break;}
        }
        require(parent[sink]!=missing,"infeasible bounded circulation");
        auto amount=Exact::subtract(required,sent);
        for(auto node=sink;node!=source;node=parent[node]){const auto capacity=residual(path[node]);
            if(Exact::sign(Exact::subtract(capacity,amount))<0)amount=capacity;}
        require(Exact::sign(amount)>0 && result.augmentations<workLimit,"exact augmentation failed to progress");
        for(auto node=sink;node!=source;node=parent[node]){const auto link=path[node];auto& flow=arcs[link.arc].flow;
            flow=link.forward?Exact::add(flow,amount):Exact::subtract(flow,amount);}
        sent=Exact::add(sent,amount);++result.augmentations;
    }
    std::fill(balance.begin(),balance.end(),Number{});
    std::vector<Number> rounded(nodes);
    for(std::size_t k=0;k<input.size();++k){const auto& a=arcs[k];
        require(Exact::sign(Exact::subtract(a.flow,a.lower))>=0 && Exact::sign(Exact::subtract(a.upper,a.flow))>=0,
            "exact selected flow escaped bounds");
        balance[a.from]=Exact::subtract(balance[a.from],a.flow);balance[a.to]=Exact::add(balance[a.to],a.flow);
        const auto enclosure=Exact::enclosure(a.flow);const double value=Exact::toDouble(a.flow);
        require(std::isfinite(value) && value>=enclosure.lo && value<=enclosure.hi
            && value>=input[k].lower && value<=input[k].upper,"rounded selected flow escaped bounds");
        result.flows.push_back(value);const auto f=Exact::fromDouble(value);
        rounded[a.from]=Exact::subtract(rounded[a.from],f);rounded[a.to]=Exact::add(rounded[a.to],f);
    }
    for(std::size_t node=0;node<nodes;++node){require(Exact::sign(balance[node])==0,"exact circulation ledger did not close");
        const auto b=Exact::enclosure(rounded[node]);result.roundedBalance.push_back({b.lo,b.hi});}
    return result;
}
}
