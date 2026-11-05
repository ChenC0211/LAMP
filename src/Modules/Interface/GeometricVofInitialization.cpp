#include "Modules/Interface/GeometricVofInitialization.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include <array>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace LAMP::Interface
{
namespace
{
namespace Cut=Geometry::PlanarCellCut;
void require(bool ok,const std::string& message)
{if(!ok)throw std::invalid_argument("[GeometricVofInitialization] "+message);}
// The represented binary64 sites are the actual sites both queried and cut.
// Never cut at hidden ideal coordinates while querying rounded world points.
std::vector<double> coordinates(double lower,double upper,unsigned count)
{
    std::vector<double> result(count+1);result.front()=lower;result.back()=upper;
    for(unsigned i=1;i<count;++i)result[i]=static_cast<double>(static_cast<long double>(lower)
        +(static_cast<long double>(upper)-lower)*i/count);
    for(unsigned i=0;i<count;++i)
        require(std::isfinite(result[i])&&std::isfinite(result[i+1])&&result[i]<result[i+1],
            "fixed subdivision has unrepresentable/invalid query sites");
    return result;
}
double midpoint(double lower,double upper)
{
    const auto result=static_cast<double>(static_cast<long double>(lower)
        +(static_cast<long double>(upper)-lower)/2);
    require(std::isfinite(result)&&lower<result&&result<upper,
        "fixed fan has no representable interior centre");
    return result;
}
}
GeometricVofPublication GeometricVofInitialization::prepare(const QuadtreeMesh& mesh,bool rz,
    const GeometricVofInitialValueQuery& query,const GeometricVofInitializationConfig& config)
{
    FaceSegmentGraph graph;graph.build(mesh);
    return prepare(mesh,NativeLeafDomain::build(mesh,graph,rz),query,config);
}
GeometricVofPublication GeometricVofInitialization::prepare(const QuadtreeMesh& mesh,NativeLeafDomain domain,
    const GeometricVofInitialValueQuery& query,const GeometricVofInitializationConfig& config)
{
    require(static_cast<bool>(query),"empty initial scalar query");
    const unsigned n=config.subcellsPerAxis;
    require(n>=1&&n<=32,"subcellsPerAxis must be in [1,32]");
    const auto tolerance=config.geometry.relativeVolumeTolerance;
    require(std::isfinite(tolerance)&&tolerance>0&&tolerance<=1e-12L,
        "volume tolerance must be finite, positive and no looser than 1e-12");
    domain.validateMeshIdentity(mesh);const bool rz=domain.axisymmetric();
    require(domain.cellRowByIndex().size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "selected domain map must cover native leaf indices");
    std::vector<std::array<long double,2>> inventories(static_cast<std::size_t>(mesh.cellIndexNum()));
    for(int native=0;native<mesh.cellIndexNum();++native){
        if(domain.cellRowByIndex()[native]<0)continue;
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
        const Cut::Rectangle rectangle{mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()};
        const auto capacity=Cut::capacity(rectangle,rz);
        const Cut::Point origin{rectangle.lower.x,rectangle.lower.y};
        const auto x=coordinates(rectangle.lower.x,rectangle.upper.x,n),
            y=coordinates(rectangle.lower.y,rectangle.upper.y,n);
        bool anyNegative=false,anyPositive=false,anyZero=false;
        const auto sample=[&](double a,double b){const auto value=query({a,b});
            require(std::isfinite(value),"nonfinite initial scalar query cell="+std::to_string(native));
            anyNegative=anyNegative||value<0;anyPositive=anyPositive||value>0;anyZero=anyZero||value==0;
            return static_cast<long double>(value);};
        const auto index=[&](unsigned i,unsigned j){return static_cast<std::size_t>(j)*(n+1)+i;};
        std::vector<long double> corner(static_cast<std::size_t>(n+1)*(n+1));
        for(unsigned j=0;j<=n;++j)for(unsigned i=0;i<=n;++i)corner[index(i,j)]=sample(x[i],y[j]);
        std::vector<long double> centres(static_cast<std::size_t>(n)*n);
        for(unsigned j=0;j<n;++j)for(unsigned i=0;i<n;++i)
            centres[static_cast<std::size_t>(j)*n+i]=sample(midpoint(x[i],x[i+1]),midpoint(y[j],y[j+1]));
        auto& stocks=inventories[static_cast<std::size_t>(native)];
        // A strictly same-sign sampled PL cell has an exactly known phase.
        // Use its independent Core rectangle metric rather than summing full
        // triangles, which could put a pure RZ stock just above capacity.
        if(!anyZero&&(!anyNegative||!anyPositive)){
            stocks[anyNegative?0:1]=capacity;continue;
        }
        std::array<long double,2> corrections{};
        const auto add=[&](unsigned phase,long double value){
            require(std::isfinite(value)&&value>=0,"invalid independent triangle phase measure");
            const auto adjusted=value-corrections[phase],next=stocks[phase]+adjusted;
            corrections[phase]=(next-stocks[phase])-adjusted;stocks[phase]=next;
        };
        for(unsigned j=0;j<n;++j)for(unsigned i=0;i<n;++i){
            const auto cx=midpoint(x[i],x[i+1]),cy=midpoint(y[j],y[j+1]);
            const Cut::Point centre{static_cast<long double>(cx)-origin.x,static_cast<long double>(cy)-origin.y};
            const auto phiCentre=centres[static_cast<std::size_t>(j)*n+i];
            const std::array<Cut::Point,4> sites{{
                {static_cast<long double>(x[i])-origin.x,static_cast<long double>(y[j])-origin.y},
                {static_cast<long double>(x[i+1])-origin.x,static_cast<long double>(y[j])-origin.y},
                {static_cast<long double>(x[i+1])-origin.x,static_cast<long double>(y[j+1])-origin.y},
                {static_cast<long double>(x[i])-origin.x,static_cast<long double>(y[j+1])-origin.y}}};
            const std::array<long double,4> phi{{corner[index(i,j)],corner[index(i+1,j)],
                corner[index(i+1,j+1)],corner[index(i,j+1)]}};
            for(unsigned k=0;k<4;++k){
                const std::array<Cut::Point,3> triangle{{centre,sites[k],sites[(k+1)%4]}};
                const std::array<long double,3> values{{phiCentre,phi[k],phi[(k+1)%4]}};
                bool negative=false,positive=false;
                for(const auto value:values){negative=negative||value<0;positive=positive||value>0;}
                require(negative||positive,"positive-area all-zero scalar triangle has no phase cell="+std::to_string(native));
                const auto measure=[&](unsigned phase,const std::array<long double,3>& cutValues){
                    try{return Cut::negativeMeasure(triangle,cutValues,origin,rz);}
                    catch(const std::exception& error){
                        std::ostringstream detail;detail.precision(std::numeric_limits<long double>::max_digits10);
                        detail<<"[GeometricVofInitialization] "<<error.what()<<" [native="<<native<<",sub_i="<<i<<",sub_j="<<j
                            <<",fan="<<k<<",phase="<<phase<<",axisymmetric="<<rz<<",capacity="<<capacity
                            <<",origin=("<<origin.x<<','<<origin.y<<"),triangle=";
                        for(unsigned q=0;q<3;++q)detail<<" ("<<triangle[q].x<<','<<triangle[q].y<<";phi="<<values[q]<<";cut="<<cutValues[q]<<')';
                        detail<<']';throw std::runtime_error(detail.str());
                    }
                };
                // Same-sign triangles have a known empty opposite phase;
                // Core measures the full occupied triangle without root work.
                if(!positive)add(0,measure(0,values));
                else if(!negative){const std::array<long double,3> opposite{{-values[0],-values[1],-values[2]}};
                    add(1,measure(1,opposite));}
                else{
                    add(0,measure(0,values));
                    const std::array<long double,3> opposite{{-values[0],-values[1],-values[2]}};
                    add(1,measure(1,opposite));
                }
            }
        }
        // A same-sign PL fan can touch zero on its boundary without owning
        // any positive-area opposite phase. The loop above has already rejected
        // every all-zero positive-area triangle. Its union is therefore exactly
        // the occupied rectangle, just as in the strict-sign fast path; use the
        // independent rectangle measure instead of a rounded full-triangle sum.
        if(!anyNegative||!anyPositive){
            stocks={};stocks[anyNegative?0:1]=capacity;
        }
        for(const auto stock:stocks)require(std::isfinite(stock)&&stock>=0&&stock<=capacity,
            "independent phase inventory outside physical capacity cell="+std::to_string(native));
        require(std::abs(stocks[0]+stocks[1]-capacity)<=tolerance*capacity,
            "independent phase inventory closure exceeds tolerance cell="+std::to_string(native));
    }
    return GeometricVofPublication::prepareLocalPlicInventory(mesh,std::move(domain),{},inventories,config.geometry);
}
}
