#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafOverlap.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Modules/Hydrodynamics/Detail/HydroVirtualFaceProlongation.h"
#include <Eigen/SVD>
#include <algorithm>
#include <atomic>
#include <cfenv>
#include <cerrno>
#include <tuple>
#if defined(_MSC_VER)
#include <float.h>
#endif
#include <cstdint>
#include <cstring>
#include <exception>
#include <memory>
#include <new>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <stdexcept>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif

namespace LAMP::Hydrodynamics::Detail
{
// An entry contains mathematical coefficients and local ordinal maps only.
// No prepared owner, mesh pointer, native index, phase, time or field survives.
struct RegularOperatorTemplateKey
{
    std::vector<std::uint64_t> words;
    bool operator<(const RegularOperatorTemplateKey& other) const noexcept {return words<other.words;}
};
struct RegularOperatorTemplateEntry
{
    bool terminal=false,hasPort=false;
    int rank=0;
    HydroRegularRow row; // dimensionless; support is local ordinal only
    HydroRegularNormalPort port; // support is canonical union ordinal only
    std::size_t accountedBytes(const RegularOperatorTemplateKey& key) const noexcept
    {
        std::size_t bytes=sizeof(*this)+sizeof(key)+512+key.words.capacity()*sizeof(std::uint64_t);
        bytes+=(row.columns.capacity()+row.support.capacity()+port.support.capacity())*sizeof(int);
        for(const auto& w:row.weights)bytes+=w.capacity()*sizeof(double);
        for(const auto& w:port.weights)bytes+=w.capacity()*sizeof(double);
        return bytes;
    }
};
struct HydroRegularOperatorTemplateBank
{
    static constexpr std::size_t entryLimit=32768,byteLimit=64U*1024U*1024U;
    std::mutex mutex;
    std::map<RegularOperatorTemplateKey,std::shared_ptr<const RegularOperatorTemplateEntry>> entries;
    std::size_t retainedBytes=0,rowEntries=0;
    HydroRegularTemplateStatistics work;
    std::shared_ptr<const RegularOperatorTemplateEntry> find(const RegularOperatorTemplateKey& key,bool& room)
    {
        std::lock_guard<std::mutex> lock(mutex);
        const auto found=entries.find(key);room=entries.size()<entryLimit&&retainedBytes<byteLimit;
        return found==entries.end()?nullptr:found->second;
    }
    bool publish(RegularOperatorTemplateKey key,std::shared_ptr<const RegularOperatorTemplateEntry> value)
    {
        const auto bytes=value->accountedBytes(key);
        std::lock_guard<std::mutex> lock(mutex);
        if(entries.size()>=entryLimit||bytes>byteLimit-retainedBytes||entries.find(key)!=entries.end())return false;
        const bool terminal=value->terminal;
        const auto inserted=entries.emplace(std::move(key),std::move(value));
        if(inserted.second){retainedBytes+=bytes;if(terminal)++rowEntries;}
        return inserted.second;
    }
    void accumulate(const HydroRegularPreparationReport& r)
    {
        std::lock_guard<std::mutex> lock(mutex);++work.preparations;
        work.lookups+=r.templateLookups;work.rowQueries+=r.templateRowQueries;
        work.rowHits+=r.templateRowHits;work.rowBuilds+=r.templateRowBuilds;work.routingHits+=r.templateRoutingHits;
        work.portQueries+=r.templatePortQueries;work.portHits+=r.templatePortHits;work.portBuilds+=r.templatePortBuilds;
        work.fallbacks+=r.templateFallbacks;work.actualSvdCalls+=r.svdCount;work.publishedEntries+=r.templatePublishedEntries;
    }
    HydroRegularTemplateStatistics statistics()
    {
        std::lock_guard<std::mutex> lock(mutex);auto result=work;
        result.entries=entries.size();result.rowEntries=rowEntries;result.routingEntries=entries.size()-rowEntries;
        result.retainedBytes=retainedBytes;result.entryLimit=entryLimit;result.byteLimit=byteLimit;return result;
    }
    void report(HydroRegularPreparationReport& out)
    {
        std::lock_guard<std::mutex> lock(mutex);
        out.templateEntryCount=entries.size();out.templateRowEntries=rowEntries;
        out.templateRoutingEntries=entries.size()-rowEntries;out.templateRetainedBytes=retainedBytes;
        out.templateEntryLimit=entryLimit;out.templateByteLimit=byteLimit;
    }
};


namespace VirtualFaceProlongationImplementation
{
void require(bool ok,const char*why){if(!ok)throw std::invalid_argument(std::string("[HydroVirtualFaceProlongation] ")+why);}
using Accumulator=std::map<int,long double>;
using Term=std::pair<const HydroVirtualRow*,long double>;
void validate(const HydroVirtualRow&row,int count)
{for(const auto&w:row)require(w.velocityDof>=0&&w.velocityDof<count&&std::isfinite(w.weight),"invalid original native weight");}
HydroVirtualRow compact(const Accumulator&source,long double scale=1)
{
    HydroVirtualRow result;result.reserve(source.size());
    for(const auto&w:source){const auto value=w.second*scale;const auto native=static_cast<double>(value);
        require(std::isfinite(value)&&std::isfinite(native),"unrepresentable derived native weight");
        require(value==0||native!=0,"derived native weight underflow");if(native!=0)result.push_back({w.first,native});}
    return result;
}
HydroVirtualRow sum(std::initializer_list<Term> terms)
{
    Accumulator result;
    for(const auto&term:terms){require(std::isfinite(term.second),"nonfinite linear map coefficient");
        for(const auto&w:*term.first)result[w.velocityDof]+=term.second*static_cast<long double>(w.weight);}
    return compact(result);
}
HydroVirtualFaceProlongation prepare(int nx,int ny,Vector2D spacing,
    int count,const std::vector<HydroVirtualCellInput>&cells,const std::vector<HydroVirtualFaceInput>&faces)
{
    require(nx>0&&ny>0&&count>0&&!cells.empty()&&faces.size()==static_cast<std::size_t>(count),"incomplete rectangular native metadata");
    require(std::isfinite(spacing.x)&&std::isfinite(spacing.y)&&spacing.x>0&&spacing.x==spacing.y,"only positive square XY spacing is supported");
    const auto cellCount=static_cast<std::size_t>(nx)*ny;
    require(cellCount<=static_cast<std::size_t>(std::numeric_limits<int>::max())&&cells.size()<=cellCount,"logical grid exceeds index scope");
    HydroVirtualFaceProlongation out;out.nx=nx;out.ny=ny;
    out.cellPressure.assign(cellCount,-1);out.rows[0].resize(static_cast<std::size_t>(nx+1)*ny);out.rows[1].resize(static_cast<std::size_t>(ny+1)*nx);
    std::array<std::vector<bool>,2> assigned{{std::vector<bool>(out.rows[0].size()),std::vector<bool>(out.rows[1].size())}};
    const auto index=[&](int axis,int x,int y){require((axis==0||axis==1)&&x>=0&&y>=0&&x<(axis==0?nx+1:nx)&&y<(axis==0?ny:ny+1),"virtual face index exceeds rectangle");
        return static_cast<std::size_t>(y)*(axis==0?nx+1:nx)+x;};
    const auto put=[&](int axis,int x,int y,HydroVirtualRow row){const auto at=index(axis,x,y);require(!assigned[axis][at],"overlapping actual or internal virtual faces");
        require(!row.empty(),"derived virtual face has no original source");assigned[axis][at]=true;out.rows[axis][at]=std::move(row);};
    for(std::size_t p=0;p<cells.size();++p){const auto&c=cells[p];const int w=c.upper[0]-c.lower[0],h=c.upper[1]-c.lower[1];
        require((w==1||w==2)&&w==h&&c.lower[0]>=0&&c.lower[1]>=0&&c.upper[0]<=nx&&c.upper[1]<=ny,"actual leaf is not a supported1x1/2x2 rectangle");
        for(int a=0;a<2;++a){validate(c.value[a],count);for(int b=0;b<2;++b)validate(c.gradient[a][b],count);}
        for(int y=c.lower[1];y<c.upper[1];++y)for(int x=c.lower[0];x<c.upper[0];++x){auto&owner=out.cellPressure[static_cast<std::size_t>(y)*nx+x];
            require(owner<0,"actual leaf interiors overlap");owner=static_cast<int>(p);}}
    for(int p:out.cellPressure)require(p>=0,"gap in complete native rectangle");
    std::vector<bool> dofSeen(count);std::vector<Accumulator> divergence(cells.size());
    for(const auto&f:faces){require(f.velocityDof>=0&&f.velocityDof<count&&!dofSeen[f.velocityDof]&&(f.normalAxis==0||f.normalAxis==1),"invalid or duplicate native face");
        dofSeen[f.velocityDof]=true;validate(f.tangentGradient,count);const int n=f.normalAxis,t=1-n,length=f.upper[t]-f.lower[t];
        require((f.negativePressure>=0||f.lower[n]==0)&&(f.positivePressure>=0||f.lower[n]==(n==0?nx:ny)),
            "exterior native face is not on the physical rectangle boundary");
        require(f.lower[n]==f.upper[n]&&(length==1||length==2)&&f.negativePressure!=f.positivePressure&&f.negativePressure>=-1&&f.positivePressure>=-1&&
            (f.negativePressure>=0||f.positivePressure>=0),"invalid actual face geometry/incidence");
        for(int side=0;side<2;++side){const int p=side==0?f.negativePressure:f.positivePressure;if(p<0)continue;
            require(static_cast<std::size_t>(p)<cells.size(),"native face pressure index outside original cells");const auto&c=cells[p];
            require(f.lower[n]==(side==0?c.upper[n]:c.lower[n])&&f.lower[t]>=c.lower[t]&&f.upper[t]<=c.upper[t],"native face is not on its declared actual cell side");
            divergence[p][f.velocityDof]+=(side==0?1.L:-1.L)*length;}
        const HydroVirtualRow unit{{f.velocityDof,1}};
        for(int k=0;k<length;++k){auto logical=f.lower;logical[t]+=k;
            const auto row=length==1?unit:sum({{&unit,1},{&f.tangentGradient,(k==0?-.5L:.5L)*spacing.x}});
            put(n,logical[0],logical[1],row);}}
    for(std::size_t p=0;p<cells.size();++p){const auto&c=cells[p];if(c.upper[0]-c.lower[0]==1)continue;const int x=c.lower[0],y=c.lower[1];
        const auto known=[&](int axis,int a,int b)->const HydroVirtualRow&{const auto at=index(axis,a,b);require(assigned[axis][at],"coarse boundary lacks an actual subface");return out.rows[axis][at];};
        // q0=(xb,xt,yl,yr); retain z.q0 with z=(1,-1,-1,1).
        const auto xb0=sum({{&c.value[0],1},{&c.gradient[0][1],-.5L*spacing.y}});
        const auto xt0=sum({{&c.value[0],1},{&c.gradient[0][1], .5L*spacing.y}});
        const auto yl0=sum({{&c.value[1],1},{&c.gradient[1][0],-.5L*spacing.x}});
        const auto yr0=sum({{&c.value[1],1},{&c.gradient[1][0], .5L*spacing.x}});
        const auto circulation=sum({{&xb0,1},{&xt0,-1},{&yl0,-1},{&yr0,1}});
        // divergence[] is Dnative/h; d=Dnative/(4h). No rounded area ratio.
        const auto d=compact(divergence[p],.25L);
        const auto b0=sum({{&d,1},{&known(0,x,y),1},{&known(1,x,y),1}});
        const auto b1=sum({{&d,1},{&known(0,x+2,y),-1},{&known(1,x+1,y),1}});
        const auto b2=sum({{&d,1},{&known(0,x,y+1),1},{&known(1,x,y+2),-1}});
        const auto xb=sum({{&circulation,.25L},{&b2,.25L},{&b0,.5L},{&b1,-.25L}});
        const auto yl=sum({{&b0,1},{&xb,-1}}),yr=sum({{&b1,1},{&xb,1}}),xt=sum({{&b2,1},{&yl,1}});
        put(0,x+1,y,xb);put(0,x+1,y+1,xt);put(1,x,y+1,yl);put(1,x+1,y+1,yr);}
    for(int axis=0;axis<2;++axis)for(bool present:assigned[axis])require(present,"gap in virtual face coverage");
    return out;
}
}
HydroVirtualFaceProlongation prepareHydroVirtualFaceProlongation(int nx,int ny,Vector2D spacing,
    int nativeVelocityCount,const std::vector<HydroVirtualCellInput>& cells,
    const std::vector<HydroVirtualFaceInput>& faces)
{
    return VirtualFaceProlongationImplementation::prepare(nx,ny,spacing,nativeVelocityCount,cells,faces);
}
namespace
{
// One prepare owns this cache. Only the main observation SVD is reusable:
// its complete matrix, options, threshold and stable full floating environment
// must match. Rows, native IDs, weights and all later inverses remain fresh.
class ExactMainSvdCache
{
    using Svd=Eigen::JacobiSVD<Eigen::MatrixXd>;
    static constexpr std::size_t entryLimit=4096,byteLimit=32U*1024U*1024U;
    struct FloatingState {
        std::array<unsigned char,sizeof(std::fenv_t)> environment{};
        unsigned csr=0;
        bool capture() noexcept
        {
            std::fenv_t raw{};if(std::fegetenv(&raw)!=0)return false;
            std::memcpy(environment.data(),&raw,sizeof(raw));
#if defined(_M_X64) || defined(__SSE2__)
            csr=_mm_getcsr();
#endif
            return true;
        }
        bool operator==(const FloatingState& other) const noexcept
        {return csr==other.csr&&environment==other.environment;}
        void restore() const noexcept
        {
            std::fenv_t raw{};std::memcpy(&raw,environment.data(),sizeof(raw));
            if(std::fesetenv(&raw)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
            _mm_setcsr(csr);
#endif
        }
    };
    // Preserve only the optional cache bookkeeping's environment. The real
    // compute/setThreshold/rank and every downstream operation are outside it.
    struct CacheEnvironment {
        FloatingState state;bool valid=false;
        CacheEnvironment() noexcept:valid(state.capture()){}
        ~CacheEnvironment() noexcept {if(valid)state.restore();}
    };
    struct Key {
        std::vector<std::uint64_t> words;
        bool operator<(const Key& other) const noexcept {return words<other.words;}
    };
    std::map<Key,std::shared_ptr<const Svd>> entries_;
#if defined(_MSC_VER) && defined(_M_X64)
    bool enabled_=sizeof(double)==8&&std::numeric_limits<double>::digits==53;
#else
    bool enabled_=false; // Other floating-environment ABIs retain literal compute.
#endif
    static std::size_t storageUpperBound(const Eigen::MatrixXd& matrix) noexcept
    {
        if(matrix.rows()<=0||matrix.cols()<=0||matrix.cols()>6||
            matrix.size()>static_cast<Eigen::Index>(byteLimit/64))return byteLimit+1;
        const auto m=static_cast<std::size_t>(matrix.rows()),n=static_cast<std::size_t>(matrix.cols());
        // Eigen 3.4.1 thin U/full V, work/scaled matrices, both possible QR
        // preconditioners, their vectors/permutations, key and allocator slack.
        // This deliberately overcounts inactive buffers retained from an earlier
        // smaller ring. The observation support only grows within one row.
        return sizeof(Svd)+sizeof(Key)+sizeof(FloatingState)+512+
            8*(6*m*n+8*n*n+16*(m+n)+256);
    }
    static bool makeKey(const Eigen::MatrixXd& matrix,unsigned options,std::size_t thresholdScale,
        const FloatingState& state,Key& key,std::uint64_t& thresholdBits)
    {
        if(options!=(Eigen::ComputeThinU|Eigen::ComputeFullV)||thresholdScale==0||
            thresholdScale>static_cast<std::size_t>(std::numeric_limits<int>::max())||
            storageUpperBound(matrix)>byteLimit)return false;
        // The original threshold is 128*2^-52*integer = 2^-45*integer.
        // Native counts fit int, so this value is exact. Encode it with integer
        // operations to leave its original floating expression after compute.
        unsigned top=0;for(auto value=thresholdScale;value>>=1;)++top;
        thresholdBits=(std::uint64_t(1023+static_cast<int>(top)-45)<<52)|
            ((std::uint64_t(thresholdScale)<<(52-top))&UINT64_C(0x000fffffffffffff));
        constexpr std::size_t environmentWords=(sizeof(std::fenv_t)+7)/8;
        key.words.reserve(5+environmentWords+static_cast<std::size_t>(matrix.size()));
        key.words.push_back(static_cast<std::uint64_t>(matrix.rows()));
        key.words.push_back(static_cast<std::uint64_t>(matrix.cols()));
        key.words.push_back(options);key.words.push_back(thresholdBits);key.words.push_back(state.csr);
        for(std::size_t i=0;i<environmentWords;++i){std::uint64_t word=0;
            const auto offset=8*i,count=std::min<std::size_t>(8,state.environment.size()-offset);
            std::memcpy(&word,state.environment.data()+offset,count);key.words.push_back(word);}
        for(Eigen::Index i=0;i<matrix.size();++i){std::uint64_t bits=0;
            std::memcpy(&bits,matrix.data()+i,sizeof(bits));
            if((bits&UINT64_C(0x7ff0000000000000))==UINT64_C(0x7ff0000000000000))return false;
            key.words.push_back(bits);}
        return true;
    }
public:
    std::size_t requests=0,hits=0,misses=0,stored=0,bytes=0;
    class Selection
    {
        ExactMainSvdCache& cache_;
        Svd literal_;std::shared_ptr<Svd> working_;std::shared_ptr<const Svd> borrowed_;
        const Svd* selected_=&literal_;
        const Eigen::MatrixXd* matrix_=nullptr;unsigned options_=0;
        std::size_t* computeCount_=nullptr;
        std::optional<Key> pending_;FloatingState before_;
        std::uint64_t thresholdBits_=0;
        bool owned_=false,stable_=false;
    public:
        explicit Selection(ExactMainSvdCache& cache):cache_(cache){}
        void compute(const Eigen::MatrixXd& matrix,unsigned options,std::size_t thresholdScale,
            std::size_t& computeCount)
        {
            matrix_=&matrix;options_=options;computeCount_=&computeCount;
            owned_=false;stable_=false;++cache_.requests;
            {
                CacheEnvironment environment;
                borrowed_.reset();pending_.reset();
                if(cache_.enabled_&&environment.valid){
                    try{
                        Key key;
                        if(makeKey(matrix,options,thresholdScale,environment.state,key,thresholdBits_)){
                            const auto found=cache_.entries_.find(key);
                            if(found!=cache_.entries_.end()){
                                borrowed_=found->second;selected_=borrowed_.get();++cache_.hits;return;
                            }
                            const auto cost=storageUpperBound(matrix);
                            if(cache_.stored<entryLimit&&cost<=byteLimit-cache_.bytes){
                                pending_=std::move(key);before_=environment.state;
                                if(!working_||working_.use_count()!=1)working_=std::make_shared<Svd>();
                                owned_=true;
                            }
                        }
                    }catch(const std::bad_alloc&){
                        // An optional cache allocation cannot reject an original
                        // fit. Keep every existing entry and use literal compute.
                        cache_.enabled_=false;pending_.reset();owned_=false;
                    }
                }
            }
            ++cache_.misses;
            Svd& actual=owned_?*working_:literal_;selected_=&actual;
            actual.compute(matrix,options);++computeCount;
        }
        void setThreshold(double threshold)
        {
            std::uint64_t bits=0;if constexpr(sizeof(double)==sizeof(bits))std::memcpy(&bits,&threshold,sizeof(bits));
            if(borrowed_){
                if(bits==thresholdBits_)return;
                // Defensive literal path if the caller ever changes the original
                // threshold expression without changing this private key encoder.
                cache_.enabled_=false;borrowed_.reset();pending_.reset();
                literal_.compute(*matrix_,options_);++*computeCount_;selected_=&literal_;
            }
            Svd& actual=owned_?*working_:literal_;actual.setThreshold(threshold);
            if(owned_&&pending_){FloatingState after;
                stable_=after.capture()&&after==before_&&bits==thresholdBits_&&actual.info()==Eigen::Success;}
        }
        void retain()
        {
            // Publish only the final selected ring; a published object is never
            // mutated or copied. Later rows only borrow const decomposition data.
            if(!owned_||!stable_||!pending_||!cache_.enabled_)return;
            CacheEnvironment environment;if(!environment.valid)return;
            try{
                const auto cost=storageUpperBound(*matrix_);
                if(cache_.stored<entryLimit&&cost<=byteLimit-cache_.bytes){
                    const auto inserted=cache_.entries_.emplace(std::move(*pending_),working_);
                    if(inserted.second){++cache_.stored;cache_.bytes+=cost;}
                }
            }catch(const std::bad_alloc&){cache_.enabled_=false;}
        }
        auto rank() const {return selected_->rank();}
        const auto& matrixU() const {return selected_->matrixU();}
        const auto& matrixV() const {return selected_->matrixV();}
        const auto& singularValues() const {return selected_->singularValues();}
    };
};

void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[HydroRegularAdvectorGeometry] ")+why);}
void fp()
{
    require(std::fegetround()==FE_TONEAREST,"round-to-nearest required; mode not changed");
#if defined(_M_X64) || defined(__SSE2__)
    require((_mm_getcsr() & ((1U<<15)|(1U<<6)))==0,"gradual underflow required; mode not changed");
#endif
}
double coord(Vector2D p,int d){return d==0?p.x:p.y;}
constexpr int exponents[10][2]={{0,0},{1,0},{0,1},{2,0},{1,1},{0,2},{3,0},{2,1},{1,2},{0,3}};
constexpr double basisPower[5][5]={{1,0,-11,18,-8},{0,0,-5,14,-8},{0,1,-4,5,-2},{0,0,1,-3,2},{0,0,16,-32,16}};
std::array<double,5> basis(double s,bool derivative=false)
{
    std::array<double,5> out{};
    for(int i=0;i<5;++i)for(int k=4;k>=(derivative?1:0);--k)
        out[i]=out[i]*s+basisPower[i][k]*(derivative?k:1);
    return out;
}
constexpr double gauss[5]={.0469100770306680036,.230765344947158454,.5,.769234655052841546,.953089922969331996};
constexpr double gaussWeight[5]={.118463442528094544,.239314335249683234,.284444444444444444,.239314335249683234,.118463442528094544};
std::array<double,10> mean(const HydroRegularFace& f,Vector2D p,Vector2D h,bool rz,int degree=2)
{
    std::array<double,10> out{};const double radius=(f.lower.x+f.upper.x)/2;
    for(int q=0;q<5;++q){const auto x=f.lower+(f.upper-f.lower)*gauss[q];
        const double a=(x.x-p.x)/h.x,b=(x.y-p.y)/h.y;
        const std::array<double,10> v{1,a,b,a*a,a*b,b*b,a*a*a,a*a*b,a*b*b,b*b*b};
        const double w=gaussWeight[q]*(rz?x.x/radius:1);
        for(int k=0;k<(degree==2?6:10);++k)out[k]+=w*v[k];}
    return out;
}
std::vector<int> columns(Vector2D p,int component,bool rz,bool wall,int degree=2)
{
    std::vector<int> out;
    for(int k=0;k<(degree==2?6:10);++k){if(wall && exponents[k][component]==0)continue;
        if(rz && p.x==0 && exponents[k][0]%2!=(component==0?1:0))continue;
        out.push_back(k);}
    return out;
}
// Orthonormal coefficient directions for the declared lower physical space,
// intersected with the same axis/wall conditions as the full P2 columns.
Eigen::MatrixXd prioritySpace(const HydroRegularRow& row,bool rz,int degree=2)
{
    std::vector<Eigen::VectorXd> vectors;
    const auto add=[&](int a,double x,int b=-1,double y=0.) {
        Eigen::VectorXd v=Eigen::VectorXd::Zero(row.columns.size());
        for(auto term:{std::pair<int,double>{a,x},{b,y}})if(term.second!=0){
            const auto it=std::find(row.columns.begin(),row.columns.end(),term.first);
            if(it==row.columns.end())return; // This entire physical direction is incompatible.
            v[it-row.columns.begin()]=term.second;
        }
        if(v.norm()>0)vectors.push_back(v.normalized());
    };
    if(!rz){add(0,1);if(degree>=1){add(1,1);add(2,1);}}
    else if(row.component==0){if(degree>=1)add(0,row.point.x,1,row.scale.x);} // physical r.
    else {
        add(0,1);if(degree>=1)add(2,1);
        // Remove the constant component of r^2; this spans the same space
        // with 1 and avoids a nearly parallel pair far from the axis.
        if(degree>=2&&!row.zeroNormal)add(1,2*row.point.x,3,row.scale.x);
    }
    Eigen::MatrixXd out(row.columns.size(),vectors.size());
    for(std::size_t k=0;k<vectors.size();++k)out.col(k)=vectors[k];return out;
}
// Resolve observable lower-degree directions first. A missing transverse
// slope must not contaminate a constant whose observations are available.
// The threshold belongs to the original lower observation matrix: scaling a
// roundoff-only projected residual by its own largest value invents rank.
Eigen::MatrixXd identifiableInverse(const HydroRegularRow& row,bool rz,
    const Eigen::MatrixXd& observations,double threshold,std::size_t& svdCount)
{
    const int count=static_cast<int>(observations.rows());
    Eigen::MatrixXd inverse=Eigen::MatrixXd::Zero(observations.cols(),count);
    Eigen::MatrixXd seen(count,0);int previous=0;
    for(int degree=0;degree<=2;++degree){
        const Eigen::MatrixXd space=prioritySpace(row,rz,degree);
        const int added=static_cast<int>(space.cols())-previous;
        if(added==0)continue;
        const Eigen::MatrixXd directions=space.rightCols(added),values=observations*directions;
        previous=static_cast<int>(space.cols());
        const Eigen::MatrixXd residual=values-seen*(seen.transpose()*values);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(residual,Eigen::ComputeThinU|Eigen::ComputeThinV);++svdCount;
        int rank=0;while(rank<svd.singularValues().size()&&svd.singularValues()[rank]>threshold)++rank;
        const Eigen::MatrixXd coefficients=directions-inverse*values;
        for(int k=0;k<rank;++k)inverse.noalias()+=coefficients*svd.matrixV().col(k)*(svd.matrixU().col(k).transpose()/svd.singularValues()[k]);
        const int old=static_cast<int>(seen.cols());seen.conservativeResize(count,old+rank);
        seen.rightCols(rank)=svd.matrixU().leftCols(rank);
    }
    require(seen.cols()==row.priorityRank,"inconsistent identifiable lower-space rank");
    return inverse;
}
// The new path has an explicit numerical contract: analytic normalized line
// moments, geometric canonical observation order and the original P2 space.
// It does not claim bit identity with world-coordinate Gaussian sampling.
struct OperatorTemplateEnvironment
{
    std::fenv_t environment{};int error=errno;unsigned csr=0;bool captured=false;
    OperatorTemplateEnvironment() noexcept
    {
        captured=std::fegetenv(&environment)==0;
#if defined(_M_X64) || defined(__SSE2__)
        csr=_mm_getcsr();
#endif
    }
    ~OperatorTemplateEnvironment() noexcept
    {
        if(captured&&std::fesetenv(&environment)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
        if(captured)_mm_setcsr(csr);
#endif
        errno=error;
    }
    bool eligible() const noexcept
    {
#if defined(_MSC_VER) && defined(_M_X64) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
        unsigned control=0;
        return captured&&sizeof(double)==8&&std::numeric_limits<double>::digits==53&&
            std::fegetround()==FE_TONEAREST&&_controlfp_s(&control,0,0)==0&&
            (control&_MCW_EM)==_MCW_EM&&(control&_MCW_RC)==_RC_NEAR&&
            (csr&0x1f80U)==0x1f80U&&(csr&((1U<<15)|(1U<<6)))==0;
#else
        return false; // Other ABIs keep the original literal preparation.
#endif
    }
};
// Optional exact routing of collinear edge containment. The original quadratic
// search below remains the complete fallback. All comparisons here are integer;
// the ONLY floating operations are the same Vector2D subtractions, after BOTH
// components have an exact-zero/zero-operand/Sterbenz certificate. Normal-or-zero
// inputs and results make the omitted repetitions FP-neutral (including DAZ's
// denormal-operand status). No numerical reconstruction is cached or changed.
struct ExactEdgeParentRecord
{
    std::uint64_t lower=0,upper=0,length=0;
    std::size_t ordinal=0;int edge=-1;
};
std::uint64_t edgeParentBits(double value) noexcept
{std::uint64_t word=0;std::memcpy(&word,&value,sizeof word);return word;}
bool edgeParentNormalOrZero(std::uint64_t word) noexcept
{
    const auto magnitude=word&UINT64_C(0x7fffffffffffffff);
    const auto exponent=magnitude&UINT64_C(0x7ff0000000000000);
    return magnitude==0||(exponent!=0&&exponent!=UINT64_C(0x7ff0000000000000));
}
std::uint64_t edgeParentOrder(std::uint64_t word) noexcept
{
    if((word&UINT64_C(0x7fffffffffffffff))==0)word=0; // Numeric -0 == +0.
    return (word>>63)!=0?~word:(word^UINT64_C(0x8000000000000000));
}
bool edgeParentExactDifference(double upper,double lower) noexcept
{
    const auto a=edgeParentBits(upper),b=edgeParentBits(lower);
    if(!edgeParentNormalOrZero(a)||!edgeParentNormalOrZero(b))return false;
    auto x=a&UINT64_C(0x7fffffffffffffff),y=b&UINT64_C(0x7fffffffffffffff);
    if(x==0||y==0||a==b)return true;
    if((a>>63)!=(b>>63))return false;
    if(x>y)std::swap(x,y);
    // For normal finite magnitudes, adding one exponent bit doubles exactly.
    // If x already has the largest exponent, real 2*x exceeds every finite y.
    return (x>>52)==2046||y<=x+UINT64_C(0x0010000000000000);
}
bool exactEdgeParents(std::vector<HydroRegularEdge>& edges,
    const std::vector<std::pair<std::pair<int,double>,std::vector<int>>>& groups) noexcept
{
    OperatorTemplateEnvironment restore;
    if(!restore.eligible()||std::numeric_limits<double>::radix!=2||
        std::numeric_limits<double>::min_exponent!=-1021||std::numeric_limits<double>::max_exponent!=1024)return false;
    try{
        constexpr std::size_t byteLimit=64U*1024U*1024U;
        const auto count=edges.size();if(count==0||count>static_cast<std::size_t>(std::numeric_limits<int>::max()))return false;
        std::size_t largest=0,total=0;
        for(const auto& group:groups){if(group.second.size()>count-total)return false;total+=group.second.size();largest=std::max(largest,group.second.size());}
        if(total!=count||largest==0)return false;
        // Exact-capacity arrays, not geometrically growing vectors. Account the
        // simultaneous five payloads plus conservative fixed control reserve.
        std::size_t bytes=1024;
        const auto account=[&](std::size_t n,std::size_t width){if(n>(byteLimit-bytes)/width)return false;bytes+=n*width;return true;};
        if(!account(count,sizeof(int))||!account(count,sizeof(unsigned char))||
            !account(largest,sizeof(ExactEdgeParentRecord))||!account(largest,sizeof(std::uint64_t))||
            !account(largest+1,sizeof(std::size_t)))return false;
        auto parents=std::make_unique<int[]>(count);
        auto seen=std::make_unique<unsigned char[]>(count);
        auto records=std::make_unique<ExactEdgeParentRecord[]>(largest);
        auto upperKeys=std::make_unique<std::uint64_t[]>(largest);
        auto tree=std::make_unique<std::size_t[]>(largest+1);
        for(std::size_t i=0;i<count;++i)parents[i]=edges[i].parent;
        const auto absent=std::numeric_limits<std::size_t>::max();
        for(const auto& group:groups){
            const int normal=group.first.first;if(normal<0||normal>1||!edgeParentNormalOrZero(edgeParentBits(group.first.second)))return false;
            const int tangent=1-normal;const auto plane=edgeParentOrder(edgeParentBits(group.first.second));
            const auto size=group.second.size();if(size==0)continue;
            for(std::size_t ordinal=0;ordinal<size;++ordinal){const int index=group.second[ordinal];
                if(index<0||static_cast<std::size_t>(index)>=count||seen[index])return false;
                seen[index]=1;const auto& edge=edges[index];
                if(edge.normalAxis!=normal||!edgeParentExactDifference(edge.upper.x,edge.lower.x)||
                    !edgeParentExactDifference(edge.upper.y,edge.lower.y))return false;
                if(edgeParentOrder(edgeParentBits(coord(edge.lower,normal)))!=plane||
                    edgeParentOrder(edgeParentBits(coord(edge.upper,normal)))!=plane)return false;
                const auto difference=edge.upper-edge.lower; // Original two-component expression.
                if(!edgeParentNormalOrZero(edgeParentBits(difference.x))||!edgeParentNormalOrZero(edgeParentBits(difference.y)))return false;
                const auto length=edgeParentBits(coord(difference,tangent));
                if((length>>63)!=0||(length&UINT64_C(0x7ff0000000000000))==0)return false;
                const auto lower=edgeParentOrder(edgeParentBits(coord(edge.lower,tangent))),upper=edgeParentOrder(edgeParentBits(coord(edge.upper,tangent)));
                if(lower>=upper)return false;
                records[ordinal]={lower,upper,length,ordinal,index};upperKeys[ordinal]=upper;
            }
            std::sort(records.get(),records.get()+size,[](const auto& a,const auto& b){return a.lower<b.lower||(a.lower==b.lower&&a.ordinal<b.ordinal);});
            std::sort(upperKeys.get(),upperKeys.get()+size);
            const auto ranks=static_cast<std::size_t>(std::unique(upperKeys.get(),upperKeys.get()+size)-upperKeys.get());
            std::fill_n(tree.get(),ranks+1,absent);
            const auto better=[&](std::size_t a,std::size_t b){
                if(a==absent)return b;if(b==absent)return a;
                return records[a].length>records[b].length||(records[a].length==records[b].length&&records[a].ordinal<records[b].ordinal)?a:b;
            };
            for(std::size_t begin=0;begin<size;){std::size_t end=begin+1;while(end<size&&records[end].lower==records[begin].lower)++end;
                // Equal lower endpoints are inclusive: publish the entire batch
                // before querying. Reversed upper ranks implement >= inclusively.
                for(std::size_t i=begin;i<end;++i){const auto rank=ranks-static_cast<std::size_t>(std::lower_bound(upperKeys.get(),upperKeys.get()+ranks,records[i].upper)-upperKeys.get());
                    for(auto k=rank;k<=ranks;k+=k&(~k+1))tree[k]=better(tree[k],i);}
                for(std::size_t i=begin;i<end;++i){const auto rank=ranks-static_cast<std::size_t>(std::lower_bound(upperKeys.get(),upperKeys.get()+ranks,records[i].upper)-upperKeys.get());
                    auto winner=absent;for(auto k=rank;k!=0;k-=k&(~k+1))winner=better(winner,tree[k]);
                    if(winner!=absent&&records[winner].length>records[i].length)parents[records[i].edge]=records[winner].edge;}
                begin=end;
            }
        }
        // No mutation precedes this nonthrowing commit. Any failed proof,
        // allocation or resource cap above replays the ORIGINAL complete loop.
        for(std::size_t i=0;i<count;++i)edges[i].parent=parents[i];
        return true;
    }catch(...){return false;}
}
struct RegularRowDependency
{
    bool complete=false;int portSpan=-1;
    std::vector<int> visitedCells;
    std::fenv_t sourceEnvironment{};unsigned sourceCsr=0;int sourceErrno=0;
    void capture(const std::vector<int>& visited,bool retainedChain,int span,
        const OperatorTemplateEnvironment& environment) noexcept
    {
        OperatorTemplateEnvironment restore;
        if(!restore.captured||!retainedChain||!environment.eligible())return;
        try{
            visitedCells.assign(visited.begin(),visited.end());portSpan=span;
            sourceEnvironment=environment.environment;sourceCsr=environment.csr;sourceErrno=environment.error;complete=true;
        }catch(...){visitedCells.clear();complete=false;}
    }
};
struct OperatorTemplateFace
{
    Vector2D lower{0,0},upper{0,0};
    int axis=0;bool negative=false,positive=false,external=false;
    auto tuple() const {return std::tie(lower.x,lower.y,upper.x,upper.y,axis,negative,positive,external);}
    bool operator<(const OperatorTemplateFace& other) const {return tuple()<other.tuple();}
    bool operator==(const OperatorTemplateFace& other) const {return tuple()==other.tuple();}
};
void templateWord(RegularOperatorTemplateKey& key,std::uint64_t word){key.words.push_back(word);}
void templateNumber(RegularOperatorTemplateKey& key,double value)
{
    if(value==0)value=0; // Canonical mathematical zero, never a tolerance bin.
    std::uint64_t word=0;std::memcpy(&word,&value,sizeof(word));templateWord(key,word);
}
void templateFaceWords(RegularOperatorTemplateKey& key,const OperatorTemplateFace& f)
{
    templateWord(key,static_cast<std::uint64_t>(f.axis));
    templateWord(key,static_cast<unsigned>(f.negative)|(static_cast<unsigned>(f.positive)<<1)|(static_cast<unsigned>(f.external)<<2));
    for(double value:{f.lower.x,f.lower.y,f.upper.x,f.upper.y})templateNumber(key,value);
}
bool normalizedTemplateFace(const HydroRegularFace& f,Vector2D point,Vector2D h,OperatorTemplateFace& out)
{
    out.axis=f.normalAxis;out.negative=f.negativePressure>=0;out.positive=f.positivePressure>=0;out.external=f.external;
    out.lower={(f.lower.x-point.x)/h.x,(f.lower.y-point.y)/h.y};
    out.upper={(f.upper.x-point.x)/h.x,(f.upper.y-point.y)/h.y};
    for(double value:{out.lower.x,out.lower.y,out.upper.x,out.upper.y})if(!std::isfinite(value))return false;
    return out.axis==0?out.lower.x==out.upper.x&&out.upper.y>out.lower.y:
        out.axis==1&&out.lower.y==out.upper.y&&out.upper.x>out.lower.x;
}
double templatePower(double value,int power)
{double result=1;for(int k=0;k<power;++k)result*=value;return result;}
// (1/(b-a))*integral_a^b x^k dx, in factored form without cancellation.
double templateIntervalMean(double a,double b,int power)
{
    double value=0;for(int k=0;k<=power;++k)value+=templatePower(a,k)*templatePower(b,power-k);
    return value/(power+1);
}
double templateMoment(const OperatorTemplateFace& f,int radial,int axial,bool rz,double rho)
{
    if(f.axis==0)return templatePower(f.lower.x,radial)*templateIntervalMean(f.lower.y,f.upper.y,axial);
    const double ordinary=templateIntervalMean(f.lower.x,f.upper.x,radial);
    const double moment=rz?(rho*ordinary+templateIntervalMean(f.lower.x,f.upper.x,radial+1))/
        (rho+templateIntervalMean(f.lower.x,f.upper.x,1)):ordinary;
    return moment*templatePower(f.lower.y,axial);
}
bool finiteTemplateEntry(const RegularOperatorTemplateEntry& entry,std::size_t observations,
    std::size_t unionCount,const std::vector<int>& columns,bool port)
{
    const auto& r=entry.row;
    if(!entry.terminal||entry.rank!=static_cast<int>(columns.size())||r.rank!=entry.rank||
        r.p2Rank!=r.rank||r.p2ColumnCount!=static_cast<int>(columns.size())||r.columns!=columns||
        r.support.size()!=observations||r.priorityRank!=r.priorityDimension||entry.hasPort!=port)return false;
    for(double value:{r.coefficientNorm,r.deficit,r.p2Deficit,r.priorityDeficit,r.nullResidual})if(!std::isfinite(value))return false;
    for(const auto& weights:r.weights){if(weights.size()!=observations)return false;for(double w:weights)if(!std::isfinite(w))return false;}
    if(port){if(entry.port.weights[0].size()!=entry.port.support.size()||entry.port.weights[1].size()!=entry.port.support.size())return false;
        for(int i:entry.port.support)if(i<0||static_cast<std::size_t>(i)>=unionCount)return false;
        for(const auto& weights:entry.port.weights)for(double w:weights)if(!std::isfinite(w))return false;}
    return true;
}
using OperatorTemplateLine=std::pair<int,double>;
using OperatorTemplateLines=std::vector<std::pair<OperatorTemplateLine,std::vector<int>>>;
using OperatorTemplateSpans=std::vector<std::pair<OperatorTemplateLine,std::vector<std::pair<int,int>>>>;
template<class Groups,class Key>
auto findGroup(Groups& groups,const Key& key)
{
    const auto found=std::lower_bound(groups.begin(),groups.end(),key,
        [](const auto& group,const auto& value){return group.first<value;});
    return found!=groups.end()&&!(key<found->first)?found:groups.end();
}
template<class Groups,class Key>
const auto& groupValues(const Groups& groups,const Key& key)
{
    const auto found=findGroup(groups,key);require(found!=groups.end(),"missing ordered geometry group");return found->second;
}
void uniqueIndices(std::vector<int>& values)
{std::sort(values.begin(),values.end());values.erase(std::unique(values.begin(),values.end()),values.end());}
void mergeIndices(std::vector<int>& target,const std::vector<int>& added)
{
    std::vector<int> merged;merged.reserve(target.size()+added.size());
    std::set_union(target.begin(),target.end(),added.begin(),added.end(),std::back_inserter(merged));target.swap(merged);
}
void expandSupport(std::vector<int>& support,const std::vector<std::vector<int>>& neighbors)
{
    const auto before=support;for(int p:before)support.insert(support.end(),neighbors[p].begin(),neighbors[p].end());uniqueIndices(support);
}
OperatorTemplateLines orderedLines(std::vector<std::pair<OperatorTemplateLine,int>> routes)
{
    std::sort(routes.begin(),routes.end(),[](const auto& a,const auto& b){
        return a.first<b.first||(!(b.first<a.first)&&a.second<b.second);});
    OperatorTemplateLines result;
    for(const auto& route:routes){if(result.empty()||result.back().first<route.first)result.push_back({route.first,{}});
        result.back().second.push_back(route.second);}
    return result;
}
// The actual line-selection/tie rule is deliberately shared with the old law:
// nearest centers, <= chooses the lower side, then the same axis power mask.
bool templatePortObservations(const std::vector<HydroRegularFace>& faces,
    const OperatorTemplateLines& lines,const OperatorTemplateSpans& segments,
    Vector2D point,int component,bool rz,std::vector<int>& observations,int* selectedSpan=nullptr)
{
    if(rz&&component==0&&point.x==0)return true;
    const int tangent=1-component;const double target=coord(point,tangent);
    const OperatorTemplateLine key{component,coord(point,component)};const auto found=findGroup(lines,key);
    const auto spansFound=findGroup(segments,key);if(found==lines.end()||spansFound==segments.end())return false;
    const auto& ids=found->second;const auto& spans=spansFound->second;
    const auto span=std::lower_bound(spans.begin(),spans.end(),target,[&](const auto& s,double x){return coord(faces[ids[s.second-1]].upper,tangent)<x;});
    if(span==spans.end()||coord(faces[ids[span->first]].lower,tangent)>target)return false;
    if(selectedSpan)*selectedSpan=static_cast<int>(span-spans.begin());
    const auto position=std::lower_bound(ids.begin()+span->first,ids.begin()+span->second,target,
        [&](int u,double x){return coord(faces[u].lower+faces[u].upper,tangent)/2<x;});
    int left=static_cast<int>(position-ids.begin())-1,right=left+1;
    while(observations.size()<4&&(left>=span->first||right<span->second)){
        if(right==span->second||(left>=span->first&&target-coord(faces[ids[left]].lower+faces[ids[left]].upper,tangent)/2<=
            coord(faces[ids[right]].lower+faces[ids[right]].upper,tangent)/2-target))observations.push_back(ids[left--]);
        else observations.push_back(ids[right++]);}
    if(rz&&component==1&&point.x==0&&observations.size()>2)observations.resize(2);
    return !observations.empty();
}
template<class GeometryData,class Datum>
bool preparedOperatorTemplate(GeometryData& d,const Datum& datum,Vector2D h,int component,int mode,bool rz,
    const std::vector<std::vector<int>>& neighbors,const OperatorTemplateLines& lines,const OperatorTemplateSpans& segments,
    HydroRegularOperatorTemplateBank& bank,HydroRegularRow& result,std::optional<HydroRegularNormalPort>& resultPort,
    RegularRowDependency* dependency=nullptr)
{
    ++d.report.templateRowQueries;
    OperatorTemplateEnvironment environment;
    if(!environment.eligible())return false;
    if(!std::isfinite(h.x)||!std::isfinite(h.y)||h.x<=0||h.y<=0)return false;
    const bool axis=rz&&datum.point.x==0,hasPort=mode==0&&datum.kind==HydroRegularRow::Kind::Vertex;
    const double rho=rz?datum.point.x/h.x:0,aspect=h.y/h.x;
    if(!std::isfinite(rho)||!std::isfinite(aspect)||aspect<=0||(rz&&!axis&&rho==0))return false;
    std::vector<int> observations;int selectedSpan=-1;bool retainedChain=true;
    if(hasPort){++d.report.templatePortQueries;
        if(!templatePortObservations(d.faces,lines,segments,datum.point,component,rz,observations,&selectedSpan))return false;}
    std::vector<OperatorTemplateFace> portFaces(observations.size());
    for(std::size_t k=0;k<observations.size();++k)
        if(!normalizedTemplateFace(d.faces[observations[k]],datum.point,h,portFaces[k]))return false;
    const auto cols=columns(datum.point,component,rz,mode==1,2);
    std::vector<int> support=datum.seeds;
    try{
        for(int ring=1;ring<=4;++ring){
            expandSupport(support,neighbors);
            std::vector<int> selected;for(int p:support)for(int u:d.cells[p].faces)
                if(d.faces[u].normalAxis==component&&d.faces[u].area>0)selected.push_back(u);uniqueIndices(selected);
            if(selected.empty())return false;
            bool exhausted=true;for(int p:support)for(int q:neighbors[p])if(!std::binary_search(support.begin(),support.end(),q))exhausted=false;
            std::vector<std::pair<OperatorTemplateFace,int>> local;
            local.reserve(selected.size());
            for(int u:selected){OperatorTemplateFace f;if(!normalizedTemplateFace(d.faces[u],datum.point,h,f))return false;local.emplace_back(f,u);}
            std::sort(local.begin(),local.end(),[](const auto& a,const auto& b){return a.first<b.first;});
            for(std::size_t i=1;i<local.size();++i)if(local[i-1].first==local[i].first)return false;
            std::vector<int> native;native.reserve(local.size());for(const auto& f:local)native.push_back(f.second);
            std::vector<std::pair<OperatorTemplateFace,int>> united=local;
            for(std::size_t i=0;i<observations.size();++i)if(!std::binary_search(selected.begin(),selected.end(),observations[i]))united.emplace_back(portFaces[i],observations[i]);
            std::sort(united.begin(),united.end(),[](const auto& a,const auto& b){return a.first<b.first;});
            for(std::size_t i=1;i<united.size();++i)if(united[i-1].first==united[i].first)return false;
            RegularOperatorTemplateKey key;key.words.reserve(16+6*(local.size()+portFaces.size()));
            templateWord(key,1); // Normalized-moment/operator format version.
            for(int value:{static_cast<int>(datum.kind),component,mode,static_cast<int>(rz),static_cast<int>(axis),ring,static_cast<int>(exhausted)})templateWord(key,value);
            templateNumber(key,rho);templateNumber(key,aspect);templateWord(key,cols.size());for(int col:cols)templateWord(key,col);
            templateWord(key,local.size());for(const auto& f:local)templateFaceWords(key,f.first);
            templateWord(key,portFaces.size());for(const auto& f:portFaces)templateFaceWords(key,f);
            bool room=false;++d.report.templateLookups;const auto retained=bank.find(key,room);
            if(retained&&!retained->terminal){
                if(retained->rank<0||retained->rank>=static_cast<int>(cols.size())||ring==4||exhausted)return false;
                ++d.report.templateRoutingHits;continue;
            }
            std::shared_ptr<const RegularOperatorTemplateEntry> completed=retained;
            if(!completed){
                if(!room)return false; // The bounded library never grows beyond its declared scope.
                HydroRegularRow row;row.kind=datum.kind;row.component=component;row.zeroNormal=mode==1;
                row.point={rho,0};row.scale={1,1};row.columns=cols;row.degree=2;row.rings=ring;
                row.support.resize(local.size());for(std::size_t i=0;i<local.size();++i)row.support[i]=static_cast<int>(i);
                Eigen::MatrixXd m(local.size(),cols.size());
                for(std::size_t i=0;i<local.size();++i)for(std::size_t j=0;j<cols.size();++j){
                    m(i,j)=templateMoment(local[i].first,exponents[cols[j]][0],exponents[cols[j]][1],rz,rho);
                    if(!std::isfinite(m(i,j)))return false;}
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(m,Eigen::ComputeThinU|Eigen::ComputeFullV);++d.report.svdCount;
                svd.setThreshold(128*std::numeric_limits<double>::epsilon()*std::max<std::size_t>(local.size(),6));
                if(svd.info()!=Eigen::Success)return false;row.rank=static_cast<int>(svd.rank());
                if(row.rank<static_cast<int>(cols.size())){
                    if(ring==4||exhausted)return false; // Uncovered deficient hierarchy uses the original path.
                    auto route=std::make_shared<RegularOperatorTemplateEntry>();route->rank=row.rank;
                    if(bank.publish(std::move(key),std::move(route)))++d.report.templatePublishedEntries;
                    continue;
                }
                Eigen::MatrixXd inverse=Eigen::MatrixXd::Zero(cols.size(),local.size());
                for(int k=0;k<row.rank;++k)inverse.noalias()+=svd.matrixV().col(k)*(svd.matrixU().col(k).transpose()/svd.singularValues()[k]);
                const Eigen::MatrixXd lower=prioritySpace(row,rz),observed=m*lower;
                row.priorityDimension=static_cast<int>(lower.cols());row.priorityRank=row.priorityDimension;
                row.p2Rank=row.rank;row.p2ColumnCount=static_cast<int>(cols.size());
                row.p2Deficit=(Eigen::MatrixXd::Identity(cols.size(),cols.size())-inverse*m).norm();
                row.priorityDeficit=(inverse*observed-lower).norm();row.coefficientNorm=inverse.norm();row.deficit=row.p2Deficit;
                constexpr int outputs[4]={0,1,2,4};
                for(int k=0;k<4;++k){auto& weights=row.weights[k];weights.assign(local.size(),0);const auto col=std::find(cols.begin(),cols.end(),outputs[k]);
                    if(col!=cols.end())for(std::size_t i=0;i<local.size();++i)weights[i]=inverse(col-cols.begin(),i);}
                auto entry=std::make_shared<RegularOperatorTemplateEntry>();entry->terminal=true;entry->rank=row.rank;entry->row=std::move(row);entry->hasPort=hasPort;
                if(hasPort){auto& port=entry->port;port.component=component;
                    if(!observations.empty()){
                        const int tangent=1-component;const bool evenAxis=rz&&component==1&&axis;
                        const std::vector<int> powers=evenAxis?std::vector<int>{0,2}:std::vector<int>{0,1,2,3};
                        const int count=static_cast<int>(observations.size()),degree=static_cast<int>(powers.size());
                        Eigen::MatrixXd moments(count,degree);
                        for(int i=0;i<count;++i)for(int k=0;k<degree;++k){
                            moments(i,k)=templateMoment(portFaces[i],tangent==0?powers[k]:0,tangent==1?powers[k]:0,rz,rho);
                            if(!std::isfinite(moments(i,k)))return false;}
                        Eigen::JacobiSVD<Eigen::MatrixXd> portSvd(moments.leftCols(count),Eigen::ComputeThinU|Eigen::ComputeThinV);++d.report.svdCount;
                        portSvd.setThreshold(128*std::numeric_limits<double>::epsilon()*std::max(count,6));
                        if(portSvd.info()!=Eigen::Success||portSvd.rank()!=count)return false;port.observedModes=count;
                        const Eigen::MatrixXd rightInverse=portSvd.solve(Eigen::MatrixXd::Identity(count,count));
                        std::map<int,std::array<double,2>> sparse;
                        const auto ordinal=[&](int u){const auto found=std::find_if(united.begin(),united.end(),[&](const auto& f){return f.second==u;});
                            return static_cast<int>(found-united.begin());};
                        for(int i=0;i<count;++i){auto& w=sparse[ordinal(observations[i])];w[0]+=rightInverse(0,i);
                            if(!evenAxis&&count>1)w[1]+=rightInverse(1,i);}
                        for(int k=count;k<degree;++k){const int power=powers[k];if(power>2)continue;
                            const int polynomial=power==0?0:power==1?1+tangent:power==2?(tangent==0?3:5):(tangent==0?6:9);
                            const auto col=std::find(cols.begin(),cols.end(),polynomial);if(col==cols.end())continue;
                            const Eigen::VectorXd correction=-rightInverse*moments.col(k);
                            for(std::size_t i=0;i<local.size();++i){const double seed=inverse(col-cols.begin(),i);auto& w=sparse[ordinal(native[i])];
                                w[0]+=correction[0]*seed;if(!evenAxis)w[1]+=(count>1?correction[1]:(power==1?1.:0.))*seed;}}
                        for(const auto& item:sparse){port.support.push_back(item.first);for(int k=0;k<2;++k)port.weights[k].push_back(item.second[k]);}
                    }
                }
                if(!finiteTemplateEntry(*entry,local.size(),united.size(),cols,hasPort))return false;
                ++d.report.templateRowBuilds;if(hasPort)++d.report.templatePortBuilds;
                completed=entry;
                if(bank.publish(std::move(key),std::move(entry)))++d.report.templatePublishedEntries;
            }else{
                if(!finiteTemplateEntry(*completed,local.size(),united.size(),cols,hasPort))return false;
                ++d.report.templateRowHits;if(hasPort)++d.report.templatePortHits;
            }
            result=completed->row;result.kind=datum.kind;result.index=datum.index;result.point=datum.point;result.scale=h;result.support=std::move(native);
            const double scales[4]={1,h.x,h.y,h.x*h.y};
            for(int k=0;k<4;++k)for(double& w:result.weights[k]){w/=scales[k];if(!std::isfinite(w))return false;}
            if(hasPort){resultPort=completed->port;resultPort->vertex=datum.index;
                for(int& u:resultPort->support)u=united[u].second;
                for(double& w:resultPort->weights[1]){w/=coord(h,1-component);if(!std::isfinite(w))return false;}}
            if(dependency)dependency->capture(support,retainedChain,selectedSpan,environment);
            return true;
        }
    }catch(const std::bad_alloc&){return false;}
    return false;
}

struct OperatorTemplateAccounting
{
    HydroRegularOperatorTemplateBank& bank;HydroRegularPreparationReport& report;
    ~OperatorTemplateAccounting() noexcept {bank.accumulate(report);}
};
// Compact C3 cardinal interpolation: four quartic pieces share the
// complete third jet. Endpoint derivatives use the same <=5-node
// Lagrange operator; H=-sum L' is consequently C2, with exact edge means.
struct CompatibleTerm {int index=0;std::array<double,5> powers{};};
struct CompatibleSegment {std::vector<CompatibleTerm> nodes,edges;};
using CompatibleInterval=std::array<CompatibleSegment,4>;
constexpr double macroBasis[4][8][5]={
    {
        {1.,0.,0.,0.,-1./24},
        {0.,0.,0.,0.,1./24},
        {0.,1./4,0.,0.,-5./192},
        {0.,0.,0.,0.,-1./64},
        {0.,0.,1./32,0.,-35./4608},
        {0.,0.,0.,0.,11./4608},
        {0.,0.,0.,1./384,-25./18432},
        {0.,0.,0.,0.,-1./6144}
    },
    {
        {23./24,-1./6,-1./4,-1./6,1./8},
        {1./24,1./6,1./4,1./6,-1./8},
        {43./192,7./48,-5./32,-5./48,13./192},
        {-1./64,-1./16,-3./32,-1./16,11./192},
        {109./4608,37./1152,-11./768,-35./1152,23./1536},
        {11./4608,11./1152,11./768,11./1152,-5./512},
        {23./18432,11./4608,-1./3072,-13./4608,23./18432},
        {-1./6144,-1./1536,-1./1024,-1./1536,13./18432}
    },
    {
        {1./2,-2./3,0.,1./3,-1./8},
        {1./2,2./3,0.,-1./3,1./8},
        {17./96,-5./24,-1./16,1./6,-11./192},
        {-17./96,-5./24,1./16,1./6,-13./192},
        {5./192,-1./36,-1./64,17./576,-5./512},
        {5./192,1./36,-1./64,-17./576,23./1536},
        {1./576,-1./576,-1./768,5./2304,-13./18432},
        {-1./576,-1./576,1./768,5./2304,-23./18432}
    },
    {
        {1./24,-1./6,1./4,-1./6,1./24},
        {23./24,1./6,-1./4,1./6,-1./24},
        {1./64,-1./16,3./32,-1./16,1./64},
        {-43./192,7./48,5./32,-5./48,5./192},
        {11./4608,-11./1152,11./768,-11./1152,11./4608},
        {109./4608,-37./1152,-11./768,35./1152,-35./4608},
        {1./6144,-1./1536,1./1024,-1./1536,1./6144},
        {-23./18432,11./4608,1./3072,-13./4608,25./18432}
    }
};
std::map<int,std::array<double,3>> compatibleDerivatives(int node,int intervals)
{
    const int count=std::min(5,intervals+1),begin=std::clamp(node-2,0,intervals+1-count);
    std::map<int,std::array<double,3>> result;
    for(int k=0;k<count;++k){
        std::array<long double,5> polynomial{};polynomial[0]=1;int degree=0;long double denominator=1;
        for(int m=0;m<count;++m)if(m!=k){
            const long double offset=begin+m-node;
            for(int power=degree+1;power>0;--power)polynomial[power]=polynomial[power-1]-offset*polynomial[power];
            polynomial[0]*=-offset;++degree;denominator*=k-m;
        }
        result[begin+k]={static_cast<double>(polynomial[1]/denominator),
            static_cast<double>(2*polynomial[2]/denominator),static_cast<double>(6*polynomial[3]/denominator)};
    }
    return result;
}
std::vector<CompatibleInterval> compatibleIntervals(int count)
{
    std::vector<std::map<int,std::array<double,3>>> derivatives(count+1);
    for(int i=0;i<=count;++i)derivatives[i]=compatibleDerivatives(i,count);
    std::vector<CompatibleInterval> result(count);
    for(int cell=0;cell<count;++cell){
        std::map<int,std::array<double,8>> datum;
        datum[cell][0]=1;datum[cell+1][1]=1;
        for(int end=0;end<2;++end)for(const auto& item:derivatives[cell+end]){
            datum[item.first][2+end]=item.second[0];datum[item.first][4+end]=item.second[1];datum[item.first][6+end]=item.second[2];
        }
        for(int half=0;half<4;++half){
            auto& out=result[cell][half];
            for(const auto& item:datum){CompatibleTerm term;term.index=item.first;
                for(int k=0;k<5;++k){long double value=0;for(int j=0;j<8;++j)value+=item.second[j]*macroBasis[half][j][k];term.powers[k]=static_cast<double>(value);}
                out.nodes.push_back(term);
            }
            // H*E_j=-H*sum_{i<=j} L_i'. In the quarter-cell chart
            // H/(H/4)=4. Partition unity makes all other E rows exact zero.
            std::array<long double,5> prefix{};
            for(std::size_t i=0;i+1<out.nodes.size();++i){
                for(int k=0;k<5;++k)prefix[k]+=out.nodes[i].powers[k];
                for(int j=out.nodes[i].index;j<out.nodes[i+1].index;++j){CompatibleTerm term;term.index=j;
                    for(int k=0;k<4;++k)term.powers[k]=static_cast<double>(-4*(k+1)*prefix[k+1]);
                    out.edges.push_back(term);
                }
            }
        }
    }
    return result;
}
struct CompatiblePatch {int index=0;Vector2D lower,upper;};
int compatibleSegment(double point,double lower,double width)
{
    int segment=0;const double quarter=width/4;
    for(int k=1;k<4;++k)if(point>=lower+k*quarter)segment=k;
    return segment;
}
CompatiblePatch compatiblePatch(const HydroRegularCell& cell,int pressure,Vector2D point)
{
    const auto h=(cell.upper-cell.lower)*.25;
    const int x=compatibleSegment(point.x,cell.lower.x,4*h.x),y=compatibleSegment(point.y,cell.lower.y,4*h.y);
    const Vector2D lower=cell.lower+Vector2D{x*h.x,y*h.y},
        upper{x==3?cell.upper.x:cell.lower.x+(x+1)*h.x,
              y==3?cell.upper.y:cell.lower.y+(y+1)*h.y};
    return {16*pressure+4*y+x,lower,upper};
}
template<class Geometry>
const HydroRegularCell& tensorCell(const Geometry& d,int index)
{return d.virtualProlongation?d.virtualCells[index]:d.cells[index];}
template<class Geometry>
int tensorCellAt(const Geometry& d,int pressure,Vector2D point)
{
    if(!d.virtualProlongation)return pressure;
    std::array<int,2> logical{};
    for(int axis=0;axis<2;++axis){const auto& lines=d.virtualLines[axis];const auto x=coord(point,axis);
        require(x>=lines.front()&&x<=lines.back(),"point outside derived virtual grid");
        auto upper=std::upper_bound(lines.begin(),lines.end(),x);if(upper==lines.end())--upper;
        require(upper!=lines.begin(),"invalid virtual cell lookup");logical[axis]=static_cast<int>(upper-lines.begin())-1;}
    return logical[1]*d.nx+logical[0];
}
template<class Geometry>
CompatiblePatch preparedCompatiblePatch(const Geometry& d,int pressure,Vector2D point)
{
    const int index=tensorCellAt(d,pressure,point);
    return compatiblePatch(tensorCell(d,index),index,point);
}
template<class Geometry>
const HydroVirtualRow& virtualRow(const Geometry& d,int slot)
{
    const auto& v=*d.virtualProlongation;const auto xCount=v.rows[0].size();
    return static_cast<std::size_t>(slot)<xCount?v.rows[0][slot]:v.rows[1][slot-xCount];
}
template<class Geometry>
std::optional<std::array<double,5>> compatibleTrace(const Geometry& d,const HydroRegularCoefficients& coefficients,int pressure,int normal,Vector2D point)
{
    const auto& cell=tensorCell(d,pressure);const int tangent=1-normal;
    if(coord(point,normal)!=coord(cell.lower,normal)&&coord(point,normal)!=coord(cell.upper,normal))return {};
    const auto logical=d.logicalCells[pressure];
    const int node=logical[normal]+(coord(point,normal)==coord(cell.upper,normal)?1:0);
    const int half=compatibleSegment(coord(point,tangent),coord(cell.lower,tangent),coord(cell.upper-cell.lower,tangent));
    std::array<long double,5> sums{};
    for(const auto& term:d.compatibleBasis[tangent][logical[tangent]][half].edges){
        const int u=d.compatibleVelocity[normal][normal==0?term.index*(d.nx+1)+node:node*d.nx+term.index];
        for(int k=0;k<5;++k)sums[k]+=coefficients.compatibleNativeMeans[u]*term.powers[k];
    }
    std::array<double,5> powers{};for(int k=0;k<5;++k)powers[k]=static_cast<double>(sums[k]);return powers;
}
template<class Geometry>
HydroAdvectorCellJet compatibleJet(const Geometry& d,const HydroRegularCoefficients& coefficients,const CompatiblePatch& piece,Vector2D point)
{
    const auto h=piece.upper-piece.lower;const double x=(point.x-piece.lower.x)/h.x,y=(point.y-piece.lower.y)/h.y;
    std::array<long double,5> bx{1,x,x*x,x*x*x,x*x*x*x},by{1,y,y*y,y*y*y,y*y*y*y};
    Vector2D value;LevelSetBoundaryTrace::Jacobian derivative;
    for(int c=0;c<2;++c){long double v=0,dx=0,dy=0;const auto& powers=coefficients.cellDofs[piece.index][c];
        for(int j=0;j<5;++j)for(int i=0;i<5;++i){const double a=powers[5*j+i];v+=a*bx[i]*by[j];
            if(i)dx+=i*a*bx[i-1]*by[j]/h.x;if(j)dy+=j*a*bx[i]*by[j-1]/h.y;}
        (c==0?value.x:value.y)=static_cast<double>(v);
        if(c==0){derivative.xx=static_cast<double>(dx);derivative.xy=static_cast<double>(dy);}
        else{derivative.yx=static_cast<double>(dx);derivative.yy=static_cast<double>(dy);}
    }
    // At a native nodal coordinate L is cardinal exactly. Restrict the
    // SAME tensor expression using that identity rather than cancellation
    // among rounded power coefficients (especially an impermeable trace).
    for(int c=0;c<2;++c)if(const auto trace=compatibleTrace(d,coefficients,piece.index/16,c,point)){
        const auto& cell=tensorCell(d,piece.index/16);const int t=1-c;
        const double width=coord(cell.upper-cell.lower,t)/4,
            lower=coord(cell.lower,t)+compatibleSegment(coord(point,t),coord(cell.lower,t),4*width)*width,
            z=(coord(point,t)-lower)/width;
        long double v=0,derivativeT=0;
        for(int k=4;k>=0;--k){derivativeT=derivativeT*z+v;v=v*z+(*trace)[k];}
        if(c==0){value.x=static_cast<double>(v);derivative.xy=static_cast<double>(derivativeT/width);}
        else{value.y=static_cast<double>(v);derivative.yx=static_cast<double>(derivativeT/width);}
    }
    require(std::isfinite(value.x)&&std::isfinite(value.y)&&std::isfinite(derivative.xx)&&std::isfinite(derivative.xy)
        &&std::isfinite(derivative.yx)&&std::isfinite(derivative.yy),"nonfinite compatible jet");
    return {{value,derivative},piece.lower,piece.upper};
}
struct Datum
{
    HydroRegularRow::Kind kind;
    int index=0;
    Vector2D point{0,0};
    std::vector<int> seeds;
    std::array<std::vector<int>,2> wallFaces;
    std::array<int,2> base{},wall{};
    std::array<int,2> normalPort{{-1,-1}};
};
void sortUnique(std::vector<double>& v)
{std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());}
// Same ordered tensor contributions, with only map lookup and the immutable
// wall-mask filter pre-resolved. No point, velocity, phase or clock is retained.
struct VirtualJetScatterPlan
{
    struct Weight {std::size_t slot=0;double weight=0;};
    struct Row {int virtualVelocity=-1;std::vector<Weight> weights;};
    std::vector<int> nativeVelocity;
    std::vector<Row> rows;
};
struct VirtualJetScatterSlot
{
    std::once_flag ready;
    VirtualJetScatterPlan plan;
};
}
// One completed topology/row stock, shared by current evaluation and remapping.
// It contains no Regular/State/Momentum Data owner or live mesh pointer.
struct RegularSourcePreparation
{
    static constexpr std::size_t metadataLimit=128u*1024u*1024u;
    bool radial=false;
    const HydroRegularOperatorTemplateBank* bank=nullptr;
    std::vector<HydroRegularRow> rows;
    std::vector<HydroRegularNormalPort> ports;
    std::vector<HydroRegularCell> cells;
    std::vector<HydroRegularFace> faces;
    std::vector<Datum> datums;
    std::vector<RegularRowDependency> dependencies;
    std::vector<std::vector<int>> neighbors;
    OperatorTemplateLines lines;
    OperatorTemplateSpans spans;
    using DatumKey=std::tuple<int,std::uint64_t,std::uint64_t>;
    std::vector<std::pair<DatumKey,std::size_t>> datumByPoint;
    std::size_t metadataBytes=0;
    static std::uint64_t bits(double x) noexcept {std::uint64_t b=0;std::memcpy(&b,&x,sizeof b);return b;}
    static DatumKey key(const Datum& datum) noexcept
    {return {static_cast<int>(datum.kind),bits(datum.point.x),bits(datum.point.y)};}
};
struct HydroRegularRemapPreparation
{
    std::shared_ptr<const RegularSourcePreparation> source;
    std::vector<int> sourceNativeByTarget;
};
struct HydroRegularAdvectorGeometry::Data
{
    explicit Data(HydroMomentumGeometry g):geometry(std::move(g)),stock(std::make_shared<RegularSourcePreparation>()),
        cells(stock->cells),faces(stock->faces),rows(stock->rows),normalPorts(stock->ports),rowDependencies(stock->dependencies),datums(stock->datums){}
    HydroMomentumGeometry geometry;
    std::shared_ptr<RegularSourcePreparation> stock; // Mutable only during this unpublished preparation.
    std::vector<Vector2D> vertices;
    std::vector<int> vertexParent,usedVertices;
    std::vector<HydroRegularCell>& cells;
    std::vector<HydroRegularEdge> edges;
    std::vector<HydroRegularFace>& faces;
    std::vector<HydroRegularRow>& rows;
    std::vector<HydroRegularNormalPort>& normalPorts;
    std::shared_ptr<const RegularSourcePreparation> completePreparation;
    std::vector<RegularRowDependency>& rowDependencies;
    bool dependencyRecording=true;
    std::vector<HydroRegularRow>& rowView() noexcept {return rows;}
    const std::vector<HydroRegularRow>& rowView() const noexcept {return completePreparation?completePreparation->rows:rows;}
    std::vector<HydroRegularNormalPort>& portView() noexcept {return normalPorts;}
    const std::vector<HydroRegularNormalPort>& portView() const noexcept {return completePreparation?completePreparation->ports:normalPorts;}
    void recordDependency(RegularRowDependency dependency) noexcept
    {
        OperatorTemplateEnvironment restore;
        if(!restore.captured||!dependencyRecording)return;
        try{rowDependencies.push_back(std::move(dependency));}
        catch(...){dependencyRecording=false;rowDependencies.clear();}
    }
    std::vector<Datum>& datums;
    std::vector<std::pair<std::pair<HydroRegularRow::Kind,int>,std::size_t>> datumIndex;
    std::size_t datum(HydroRegularRow::Kind kind,int index) const
    {return groupValues(datumIndex,std::make_pair(kind,index));}
    HydroRegularPreparationReport report;
    bool compatible=false;int nx=0,ny=0;
    std::vector<std::array<int,2>> logicalCells;
    std::array<std::vector<int>,2> compatibleVelocity;
    std::array<std::vector<CompatibleInterval>,2> compatibleBasis;
    std::optional<HydroVirtualFaceProlongation> virtualProlongation;
    std::vector<HydroRegularCell> virtualCells;
    std::array<std::vector<double>,2> virtualLines;
    // Data is already bound to the original exact NativeMeshIdentity/layout.
    // The key also retains the actual tensor cell and complete native mask;
    // a mesh address alone never selects this cache. New AMR data starts empty.
    mutable std::mutex virtualJetScatterMutex;
    mutable std::map<std::pair<int,std::vector<int>>,std::shared_ptr<VirtualJetScatterSlot>> virtualJetScatter;
    mutable std::atomic<std::size_t> virtualJetScatterPrepared{0};
};
template<class GeometryData>
void retainRegularSourcePreparation(GeometryData& d,std::vector<std::vector<int>>& neighbors,
    OperatorTemplateLines& lines,OperatorTemplateSpans& spans,
    const HydroRegularOperatorTemplateBank* bank) noexcept
{
    OperatorTemplateEnvironment environment;
    if(!environment.eligible()||d.compatible||!d.dependencyRecording||d.rowDependencies.size()!=d.rows.size())return;
    if(std::none_of(d.rowDependencies.begin(),d.rowDependencies.end(),[](const auto& r){return r.complete;}))return;
    try{
        std::vector<std::pair<RegularSourcePreparation::DatumKey,std::size_t>> index;
        index.reserve(d.datums.size());
        for(std::size_t k=0;k<d.datums.size();++k)index.emplace_back(RegularSourcePreparation::key(d.datums[k]),k);
        std::sort(index.begin(),index.end(),[](const auto& a,const auto& b){return a.first<b.first;});
        for(std::size_t k=1;k<index.size();++k)if(index[k-1].first==index[k].first)return;
        std::size_t bytes=sizeof(RegularSourcePreparation)+512u;
        const auto add=[&](std::size_t count,std::size_t width){
            if(count>(RegularSourcePreparation::metadataLimit-bytes)/width)return false;bytes+=count*width;return true;};
        if(!add(d.cells.capacity(),sizeof(HydroRegularCell))||!add(d.faces.capacity(),sizeof(HydroRegularFace))
            ||!add(d.datums.capacity(),sizeof(Datum))||!add(neighbors.capacity(),sizeof(std::vector<int>))
            ||!add(d.rowDependencies.capacity(),sizeof(RegularRowDependency))
            ||!add(index.capacity(),sizeof(index[0]))||!add(lines.capacity(),sizeof(lines[0]))||!add(spans.capacity(),sizeof(spans[0])))return;
        for(const auto& cell:d.cells){if(!add(cell.faces.capacity(),sizeof(int)))return;
            for(const auto& knots:cell.supportKnots)if(!add(knots.capacity(),sizeof(double)))return;}
        for(const auto& datum:d.datums){if(!add(datum.seeds.capacity(),sizeof(int)))return;
            for(const auto& faces:datum.wallFaces)if(!add(faces.capacity(),sizeof(int)))return;}
        for(const auto& row:d.rowDependencies)if(!add(row.visitedCells.capacity(),sizeof(int)))return;
        for(const auto& row:neighbors)if(!add(row.capacity(),sizeof(int)))return;
        for(const auto& line:lines)if(!add(line.second.capacity(),sizeof(int)))return;
        for(const auto& line:spans)if(!add(line.second.capacity(),sizeof(std::pair<int,int>)))return;
        // Current evaluation and future remapping share this one completed stock.
        // No old Regular/State/Momentum owner or full topology copy is retained.
        auto& completed=*d.stock;completed.radial=d.geometry.layout().domain().axisymmetric();completed.bank=bank;
        completed.neighbors=std::move(neighbors);completed.lines=std::move(lines);completed.spans=std::move(spans);
        completed.datumByPoint=std::move(index);completed.metadataBytes=bytes;d.completePreparation=d.stock;
    }catch(...){return;} // Optional source retention; current numerical rows remain intact.
}
struct RegularRemapContext
{
    std::shared_ptr<const HydroRegularRemapPreparation> request;
    std::vector<int> pressureMap,faceMap;
    std::vector<unsigned char> cleanCell;
    std::vector<std::pair<std::pair<OperatorTemplateLine,int>,bool>> cleanSpan;
    HydroRegularPreparationReuseReport report;
    bool ready=false;
    const OperatorTemplateLines* targetLines=nullptr;
    const OperatorTemplateSpans* targetSpans=nullptr;
    template<class GeometryData>
    RegularRemapContext(const std::shared_ptr<const HydroRegularRemapPreparation>& offer,const GeometryData& d,
        const std::vector<std::vector<int>>& neighbors,const OperatorTemplateLines& lines,const OperatorTemplateSpans& spans,
        const HydroRegularOperatorTemplateBank* currentBank,bool requested):request(offer),targetLines(&lines),targetSpans(&spans)
    {
        report.requested=requested;report.targetDatums=d.datums.size();OperatorTemplateEnvironment restore;
        if(!restore.eligible()||!request||!request->source||d.compatible)return;
        const auto& source=*request->source;
        if(source.bank!=currentBank||source.radial!=d.geometry.layout().domain().axisymmetric())return;
        try{
            pressureMap.assign(source.cells.size(),-1);faceMap.assign(source.faces.size(),-1);cleanCell.assign(source.cells.size(),0);
            std::vector<std::pair<int,int>> pressureByNative;pressureByNative.reserve(source.cells.size());
            for(std::size_t p=0;p<source.cells.size();++p)pressureByNative.emplace_back(source.cells[p].nativeCell,static_cast<int>(p));
            std::sort(pressureByNative.begin(),pressureByNative.end());
            for(std::size_t p=0;p<d.cells.size();++p){const int native=d.cells[p].nativeCell;
                if(native<0||static_cast<std::size_t>(native)>=request->sourceNativeByTarget.size())continue;
                const auto found=findGroup(pressureByNative,request->sourceNativeByTarget[native]);if(found==pressureByNative.end())continue;
                const int old=found->second;if(pressureMap[old]!=-1){pressureMap[old]=-2;continue;}
                const auto& a=source.cells[old];const auto& b=d.cells[p];
                if(!samePoint(a.lower,b.lower)||!samePoint(a.upper,b.upper))continue;pressureMap[old]=static_cast<int>(p);
            }
            using FaceKey=std::tuple<int,std::uint64_t,std::uint64_t,std::uint64_t,std::uint64_t>;
            const auto faceKey=[](const HydroRegularFace& f){return FaceKey{f.normalAxis,
                RegularSourcePreparation::bits(f.lower.x),RegularSourcePreparation::bits(f.lower.y),
                RegularSourcePreparation::bits(f.upper.x),RegularSourcePreparation::bits(f.upper.y)};};
            std::vector<std::pair<FaceKey,int>> targetFaces;targetFaces.reserve(d.faces.size());
            for(std::size_t u=0;u<d.faces.size();++u)targetFaces.emplace_back(faceKey(d.faces[u]),static_cast<int>(u));
            std::sort(targetFaces.begin(),targetFaces.end());
            for(std::size_t u=1;u<targetFaces.size();++u)if(targetFaces[u-1].first==targetFaces[u].first)return;
            for(std::size_t u=0;u<source.faces.size();++u){const auto& a=source.faces[u];const auto found=findGroup(targetFaces,faceKey(a));
                if(found==targetFaces.end())continue;const auto& b=d.faces[found->second];
                const auto mapped=[&](int p){return p<0?-1:(static_cast<std::size_t>(p)<pressureMap.size()?pressureMap[p]:-2);};
                if((a.negativePressure>=0&&mapped(a.negativePressure)<0)||(a.positivePressure>=0&&mapped(a.positivePressure)<0)
                    ||mapped(a.negativePressure)!=b.negativePressure||mapped(a.positivePressure)!=b.positivePressure
                    ||a.external!=b.external||std::memcmp(&a.area,&b.area,sizeof(a.area))!=0
                    ||std::memcmp(a.supportWidth.data(),b.supportWidth.data(),sizeof(a.supportWidth))!=0)continue;
                faceMap[u]=found->second;
            }
            for(std::size_t p=0;p<source.cells.size();++p){const int target=pressureMap[p];if(target<0)continue;
                std::vector<int> faces,adjacent;bool complete=true;
                for(int u:source.cells[p].faces){if(u<0||static_cast<std::size_t>(u)>=faceMap.size()||faceMap[u]<0){complete=false;break;}faces.push_back(faceMap[u]);}
                if(!complete)continue;std::sort(faces.begin(),faces.end());
                if(faces!=d.cells[target].faces)continue;
                for(int q:source.neighbors[p]){if(q<0||static_cast<std::size_t>(q)>=pressureMap.size()||pressureMap[q]<0){complete=false;break;}adjacent.push_back(pressureMap[q]);}
                if(!complete)continue;std::sort(adjacent.begin(),adjacent.end());
                if(adjacent.size()!=neighbors[target].size()||!std::equal(adjacent.begin(),adjacent.end(),neighbors[target].begin()))continue;
                cleanCell[p]=1;
            }
            for(const auto& item:source.spans){const auto& sourceIds=groupValues(source.lines,item.first);
                const auto targetLine=findGroup(lines,item.first);const auto targetSpans=findGroup(spans,item.first);
                for(std::size_t k=0;k<item.second.size();++k){bool good=false;const auto span=item.second[k];
                    if(targetLine!=lines.end()&&targetSpans!=spans.end()&&span.first<span.second){
                        const int first=sourceIds[span.first],last=sourceIds[span.second-1];
                        if(faceMap[first]>=0&&faceMap[last]>=0)for(const auto& candidate:targetSpans->second){
                            if(candidate.second-candidate.first!=span.second-span.first)continue;
                            if(targetLine->second[candidate.first]!=faceMap[first]||targetLine->second[candidate.second-1]!=faceMap[last])continue;
                            good=true;for(int i=0;i<span.second-span.first;++i){const int old=sourceIds[span.first+i];
                                if(faceMap[old]<0||targetLine->second[candidate.first+i]!=faceMap[old]){good=false;break;}}
                            if(good)break;
                        }
                    }
                    cleanSpan.emplace_back(std::make_pair(item.first,static_cast<int>(k)),good);
                }
            }
            report.enabled=true;report.retainedMetadataBytes=source.metadataBytes;ready=true;
        }catch(...){request.reset();pressureMap.clear();faceMap.clear();cleanCell.clear();cleanSpan.clear();ready=false;}
    }
    static bool samePoint(Vector2D a,Vector2D b) noexcept
    {return RegularSourcePreparation::bits(a.x)==RegularSourcePreparation::bits(b.x)
        &&RegularSourcePreparation::bits(a.y)==RegularSourcePreparation::bits(b.y);}
    static bool sameControl(const OperatorTemplateEnvironment& current,const RegularRowDependency& source) noexcept
    {
#if defined(_MSC_VER) && defined(_M_X64)
        auto a=current.environment,b=source.sourceEnvironment;a._Fe_stat=0;b._Fe_stat=0;
        return std::memcmp(&a,&b,sizeof a)==0&&(current.csr&~0x3fu)==(source.sourceCsr&~0x3fu);
#else
        return false;
#endif
    }
    bool row(const Datum& datum,Vector2D h,int component,int mode,const std::vector<int>& cols,
        HydroRegularRow& result,std::optional<HydroRegularNormalPort>& resultPort,RegularRowDependency& dependency) noexcept
    {
        if(!ready)return false;OperatorTemplateEnvironment environment;
        if(!environment.eligible()){++report.fpRejectedRows;return false;}
        try{
            const auto& source=*request->source;const auto found=findGroup(source.datumByPoint,RegularSourcePreparation::key(datum));
            if(found==source.datumByPoint.end()){++report.stencilRejectedRows;return false;}
            const auto& priorDatum=source.datums[found->second];const int index=(mode==0?priorDatum.base:priorDatum.wall)[component];
            if(index<0||static_cast<std::size_t>(index)>=source.rows.size()||static_cast<std::size_t>(index)>=source.dependencies.size())return false;
            const auto& prior=source.rows[index];const auto& certificate=source.dependencies[index];
            if(!certificate.complete){++report.stencilRejectedRows;return false;}
            if(!sameControl(environment,certificate)){++report.fpRejectedRows;return false;}
            if(prior.component!=component||prior.zeroNormal!=(mode==1)||prior.columns!=cols||!samePoint(prior.point,datum.point)||!samePoint(prior.scale,h)){
                ++report.columnRejectedRows;return false;}
            std::vector<int> seeds,walls;
            for(int p:priorDatum.seeds){if(p<0||static_cast<std::size_t>(p)>=pressureMap.size()||pressureMap[p]<0){++report.stencilRejectedRows;return false;}seeds.push_back(pressureMap[p]);}
            std::sort(seeds.begin(),seeds.end());
            if(seeds.size()!=datum.seeds.size()||!std::equal(seeds.begin(),seeds.end(),datum.seeds.begin())){++report.stencilRejectedRows;return false;}
            for(int p:certificate.visitedCells)if(p<0||static_cast<std::size_t>(p)>=cleanCell.size()||!cleanCell[p]){++report.stencilRejectedRows;return false;}
            for(int u:priorDatum.wallFaces[component]){if(u<0||static_cast<std::size_t>(u)>=faceMap.size()||faceMap[u]<0){++report.stencilRejectedRows;return false;}walls.push_back(faceMap[u]);}
            if(walls!=datum.wallFaces[component]){++report.stencilRejectedRows;return false;}
            const bool hasPort=mode==0&&datum.kind==HydroRegularRow::Kind::Vertex;
            if(hasPort&&!(source.radial&&component==0&&datum.point.x==0)){
                const OperatorTemplateLine key{component,coord(datum.point,component)};
                const auto span=findGroup(cleanSpan,std::make_pair(key,certificate.portSpan));
                if(span==cleanSpan.end()||!span->second){++report.portSpanRejectedRows;return false;}
            }
            // Only canonical normalized-helper rows are complete. Equal actual
            // face geometry, point/scale and full selected dependency closure
            // reproduce that helper's geometric order, even if native labels
            // permute. Literal sorted-native fallback rows are never eligible.
            auto completed=prior;completed.index=datum.index;
            for(int& u:completed.support){if(u<0||static_cast<std::size_t>(u)>=faceMap.size()||faceMap[u]<0){++report.stencilRejectedRows;return false;}u=faceMap[u];}
            std::optional<HydroRegularNormalPort> port;
            if(hasPort){const int oldPort=priorDatum.normalPort[component];
                if(oldPort<0||static_cast<std::size_t>(oldPort)>=source.ports.size())return false;
                port=source.ports[oldPort];port->vertex=datum.index;
                for(int& u:port->support){if(u<0||static_cast<std::size_t>(u)>=faceMap.size()||faceMap[u]<0){++report.portSpanRejectedRows;return false;}u=faceMap[u];}
            }
            // A remapped completed row may itself become the next independent
            // source receipt. Its dependency labels are all fresh target labels.
            auto nextDependency=certificate;nextDependency.sourceEnvironment=environment.environment;
            nextDependency.sourceCsr=environment.csr;nextDependency.sourceErrno=environment.error;
            for(int& p:nextDependency.visitedCells)p=pressureMap[p];
            std::sort(nextDependency.visitedCells.begin(),nextDependency.visitedCells.end());
            if(hasPort&&!(source.radial&&component==0&&datum.point.x==0)){
                // Span indices need not survive insertion of unrelated spans.
                const OperatorTemplateLine key{component,coord(datum.point,component)};
                const auto& oldIds=groupValues(source.lines,key);const auto oldSpan=groupValues(source.spans,key)[certificate.portSpan];
                const int first=faceMap[oldIds[oldSpan.first]],last=faceMap[oldIds[oldSpan.second-1]];
                const auto& ids=groupValues(*targetLines,key);const auto& spans=groupValues(*targetSpans,key);nextDependency.portSpan=-1;
                for(std::size_t k=0;k<spans.size();++k)if(ids[spans[k].first]==first&&ids[spans[k].second-1]==last){nextDependency.portSpan=static_cast<int>(k);break;}
                if(nextDependency.portSpan<0)return false;
            }
            result=std::move(completed);resultPort=std::move(port);dependency=std::move(nextDependency);
            ++report.reusedRows;if(hasPort)++report.reusedPorts;return true;
        }catch(...){return false;} // Optional resource failure retains literal helper.
    }
};
std::shared_ptr<const HydroRegularRemapPreparation> HydroRegularAdvectorGeometry::makeRemapPreparation(
    const QuadtreeMesh& sourceMesh,const HydroRegularAdvectorGeometry& source,
    const QuadtreeMesh& targetMesh,const Mesh::NativeLeafOverlap& overlaps)
{
    source.validateFor(sourceMesh,source.momentumGeometry());overlaps.validateFor(sourceMesh,targetMesh);
    OperatorTemplateEnvironment restore;
    if(!restore.eligible()||source.usesCompatibleTensor()||!source.data_->completePreparation)return {};
    try{
        auto request=std::make_shared<HydroRegularRemapPreparation>();request->source=source.data_->completePreparation;
        request->sourceNativeByTarget.assign(overlaps.byTarget().size(),-1);
        const auto sameRectangle=[](const auto& a,const auto& b){return RegularRemapContext::samePoint(a.lower,b.lower)
            &&RegularRemapContext::samePoint(a.upper,b.upper);};
        for(std::size_t target=0;target<overlaps.byTarget().size();++target){const auto& ids=overlaps.byTarget()[target];if(ids.size()!=1)continue;
            const auto& piece=overlaps.pieces()[ids.front()];if(piece.kind!=Mesh::NativeLeafOverlapPiece::Kind::Unchanged
                ||piece.targetCell!=static_cast<int>(target)||piece.sourceCell<0
                ||static_cast<std::size_t>(piece.sourceCell)>=overlaps.bySource().size()
                ||overlaps.bySource()[piece.sourceCell].size()!=1
                ||overlaps.bySource()[piece.sourceCell].front()!=ids.front()
                ||!sameRectangle(piece.rectangle,overlaps.sourceRectangles()[piece.sourceCell])
                ||!sameRectangle(piece.rectangle,overlaps.targetRectangles()[target]))continue;
            request->sourceNativeByTarget[target]=piece.sourceCell;
        }
        return request;
    }catch(...){return {};}
}
HydroRegularAdvectorGeometry::HydroRegularAdvectorGeometry(std::shared_ptr<const Data> d):data_(std::move(d)){}
const HydroMomentumGeometry& HydroRegularAdvectorGeometry::momentumGeometry() const noexcept{return data_->geometry;}
const std::vector<HydroRegularCell>& HydroRegularAdvectorGeometry::cells() const noexcept{return data_->cells;}
const std::vector<HydroRegularEdge>& HydroRegularAdvectorGeometry::edges() const noexcept{return data_->edges;}
const std::vector<HydroRegularFace>& HydroRegularAdvectorGeometry::faces() const noexcept{return data_->faces;}
const std::vector<HydroRegularRow>& HydroRegularAdvectorGeometry::rows() const noexcept{return data_->rowView();}
const std::vector<HydroRegularNormalPort>& HydroRegularAdvectorGeometry::normalPorts() const noexcept{return data_->portView();}
const HydroRegularPreparationReport& HydroRegularAdvectorGeometry::report() const noexcept{return data_->report;}
bool HydroRegularAdvectorGeometry::usesCompatibleTensor() const noexcept{return data_->compatible;}
bool HydroRegularAdvectorGeometry::usesVirtualCompatibleTensor() const noexcept{return data_->virtualProlongation.has_value();}
std::size_t HydroRegularAdvectorGeometry::cachedVirtualJetPlans() const noexcept
{return data_->virtualJetScatterPrepared.load(std::memory_order_relaxed);}
bool HydroRegularAdvectorGeometry::sameIdentity(const HydroRegularAdvectorGeometry& other) const noexcept{return data_==other.data_;}
void HydroRegularAdvectorGeometry::validateFor(const QuadtreeMesh& mesh,const HydroMomentumGeometry& expected) const
{fp();momentumGeometry().validateMesh(mesh);require(momentumGeometry().sameIdentity(expected),"different momentum geometry instance");}

HydroRegularAdvectorGeometry HydroRegularAdvectorGeometry::prepare(const QuadtreeMesh& mesh,HydroMomentumGeometry geometry,
    const HydroAdvectorBoundary* boundary)
{
    fp(); // The private caller has checked this live mesh and boundary, without a callback.
    bool remapRequested=false;
    const auto remapRequest=[&]{OperatorTemplateEnvironment restore;return geometry.takeRegularRemapPreparation(remapRequested);}();
    ExactMainSvdCache mainSvdCache;
    auto d=std::make_shared<Data>(std::move(geometry));
    const auto templateBank=d->geometry.regularTemplates();
    OperatorTemplateAccounting templateAccounting{*templateBank,d->report};
    const auto& l=d->geometry.layout();const bool rz=l.domain().axisymmetric();
    const bool compatible=!rz&&l.domain().uniformLeaves()&&l.pressureCount()==mesh.cellIndexNum()
        &&(!boundary||boundary->completeBoundaryLines());
    constexpr int degree=2; // actual general AMR/RZ reconstruction only
    d->vertices.resize(mesh.vertIndexNum());d->vertexParent.assign(mesh.vertIndexNum(),-1);
    d->cells.resize(l.pressureCount());d->faces.resize(l.velocityCount());
    using VertexUse=std::pair<int,int>;
    using EdgeUse=std::pair<std::array<int,2>,int>;
    std::vector<VertexUse> vertexCells;
    std::vector<EdgeUse> edgeCells;
    require(d->cells.size()<=vertexCells.max_size()/4&&d->cells.size()<=edgeCells.max_size()/4,
        "regular topology incidence overflow");
    vertexCells.reserve(4*d->cells.size());edgeCells.reserve(4*d->cells.size());
    constexpr int corners[4][2]={{0,1},{3,2},{0,3},{1,2}};
    for(int p=0;p<l.pressureCount();++p){auto& c=d->cells[p];c.nativeCell=l.nativeCell(p);
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c.nativeCell]);
        for(int k=0;k<4;++k){const auto& v=mesh.verts(cell.cellToVert[k]);c.vertices[k]=v.index();d->vertices[v.index()]=v.coord();vertexCells.emplace_back(v.index(),p);}
        c.lower=d->vertices[c.vertices[0]];c.upper=d->vertices[c.vertices[2]];
        for(const auto& k:corners)edgeCells.push_back({{c.vertices[k[0]],c.vertices[k[1]]},p});
    }
    std::sort(vertexCells.begin(),vertexCells.end());vertexCells.erase(std::unique(vertexCells.begin(),vertexCells.end()),vertexCells.end());
    std::sort(edgeCells.begin(),edgeCells.end());edgeCells.erase(std::unique(edgeCells.begin(),edgeCells.end()),edgeCells.end());
    std::vector<std::vector<int>> seeds;
    std::vector<std::pair<OperatorTemplateLine,int>> edgeRoutes;
    for(std::size_t begin=0;begin<edgeCells.size();){std::size_t end=begin+1;
        while(end<edgeCells.size()&&edgeCells[end].first==edgeCells[begin].first)++end;
        require(d->edges.size()<static_cast<std::size_t>(std::numeric_limits<int>::max()),"regular edge index overflow");
        HydroRegularEdge e;e.vertices=edgeCells[begin].first;e.lower=d->vertices[e.vertices[0]];e.upper=d->vertices[e.vertices[1]];
        e.normalAxis=e.lower.y==e.upper.y?1:0;const int i=static_cast<int>(d->edges.size());
        d->edges.push_back(e);seeds.emplace_back();seeds.back().reserve(end-begin);
        for(auto k=begin;k<end;++k)seeds.back().push_back(edgeCells[k].second);
        edgeRoutes.push_back({{e.normalAxis,coord(e.lower,e.normalAxis)},i});begin=end;
    }
    for(auto& c:d->cells)for(int k=0;k<4;++k){const std::array<int,2> key{c.vertices[corners[k][0]],c.vertices[corners[k][1]]};
        const auto found=std::lower_bound(d->edges.begin(),d->edges.end(),key,
            [](const auto& edge,const auto& value){return edge.vertices<value;});
        require(found!=d->edges.end()&&found->vertices==key,"missing regular cell edge");c.edges[k]=static_cast<int>(found-d->edges.begin());}
    std::vector<EdgeUse>().swap(edgeCells);
    const auto groups=orderedLines(std::move(edgeRoutes));
    if(!exactEdgeParents(d->edges,groups)){
    for(const auto& group:groups)for(int i:group.second){auto& e=d->edges[i];const int t=1-e.normalAxis;double longest=coord(e.upper-e.lower,t);
        for(int j:group.second){const auto& other=d->edges[j];const double length=coord(other.upper-other.lower,t);
            if(length>longest && coord(other.lower,t)<=coord(e.lower,t) && coord(other.upper,t)>=coord(e.upper,t)){e.parent=j;longest=length;}}}
    }
    for(std::size_t i=0;i<d->edges.size();++i){const auto& e=d->edges[i];if(e.parent<0)continue;const auto& master=d->edges[e.parent];
        mergeIndices(seeds[e.parent],seeds[i]);const int t=1-e.normalAxis;
        for(int v:e.vertices)if(coord(master.lower,t)<coord(d->vertices[v],t) && coord(d->vertices[v],t)<coord(master.upper,t)){
            require(d->vertexParent[v]<0 || d->vertexParent[v]==e.parent,"incompatible active master sides");d->vertexParent[v]=e.parent;}}
    std::vector<std::vector<int>> neighbors(d->cells.size());
    for(int u=0;u<l.velocityCount();++u){auto& f=d->faces[u];const auto actual=l.domain().face(mesh,l.nativeFace(u));
        const auto& edge=actual.horizontal?mesh.hEdges(actual.edgeLabel):mesh.vEdges(actual.edgeLabel);
        f.normalAxis=actual.horizontal?1:0;f.lower=mesh.verts(edge.edgeToVert[0]).coord();f.upper=mesh.verts(edge.edgeToVert[1]).coord();
        if(coord(f.lower,1-f.normalAxis)>coord(f.upper,1-f.normalAxis))std::swap(f.lower,f.upper);
        f.negativePressure=actual.negativeCell<0?-1:l.pressureDof(actual.negativeCell);f.positivePressure=actual.positiveCell<0?-1:l.pressureDof(actual.positiveCell);
        f.area=d->geometry.metrics().faceMeasures()[u];f.external=actual.outwardSign!=0;
        const int ids[2]={f.negativePressure,f.positivePressure};
        for(int side=0;side<2;++side)if(ids[side]>=0){auto& c=d->cells[ids[side]];c.faces.push_back(u);f.supportWidth[side]=coord(c.upper-c.lower,f.normalAxis)/2;
            c.supportKnots[f.normalAxis].push_back(coord(c.lower+c.upper,f.normalAxis)/2);
            for(int k=0;k<2;++k){c.supportKnots[k].push_back(coord(f.lower,k));c.supportKnots[k].push_back(coord(f.upper,k));}}
        if(ids[0]>=0 && ids[1]>=0){neighbors[ids[0]].push_back(ids[1]);neighbors[ids[1]].push_back(ids[0]);}
        if(f.area>0){const auto& cell=d->cells[ids[0]>=0?ids[0]:ids[1]];const auto h=cell.upper-cell.lower;
            for(int q=0;q<5;++q){const auto x=f.lower+(f.upper-f.lower)*gauss[q];const auto bx=basis((x.x-cell.lower.x)/h.x),by=basis((x.y-cell.lower.y)/h.y);
                const double w=gaussWeight[q]*(rz?x.x/((f.lower.x+f.upper.x)/2):1);
                for(int y=0;y<5;++y)for(int xindex=0;xindex<5;++xindex)f.meanWeights[5*y+xindex]+=w*bx[xindex]*by[y];}}
    }
    for(auto& adjacent:neighbors)uniqueIndices(adjacent);
    for(auto& cell:d->cells)for(auto& knots:cell.supportKnots)sortUnique(knots);
    // Select a derived chart from metadata, never as a failed-fit fallback.
    bool virtualCandidate=!compatible&&!rz&&l.pressureCount()==mesh.cellIndexNum()&&
        (!boundary||boundary->completeBoundaryLines());
    double virtualSpacing=std::numeric_limits<double>::infinity();int virtualNx=0,virtualNy=0;
    const Vector2D virtualOrigin=mesh.minCoord;
    for(const auto& cell:d->cells)virtualSpacing=std::min(virtualSpacing,std::min(cell.upper.x-cell.lower.x,cell.upper.y-cell.lower.y));
    if(virtualCandidate){const auto ext=mesh.maxCoord-mesh.minCoord;const double nx=ext.x/virtualSpacing,ny=ext.y/virtualSpacing;
        virtualCandidate=std::isfinite(nx)&&std::isfinite(ny)&&nx>=2&&ny>=2&&
            nx<std::numeric_limits<int>::max()&&ny<std::numeric_limits<int>::max();
        if(virtualCandidate){virtualNx=static_cast<int>(std::llround(nx));virtualNy=static_cast<int>(std::llround(ny));
            virtualCandidate=virtualOrigin.x+virtualNx*virtualSpacing==mesh.maxCoord.x&&
                virtualOrigin.y+virtualNy*virtualSpacing==mesh.maxCoord.y&&
                static_cast<std::size_t>(virtualNx)*virtualNy<=4*d->cells.size();}
        bool fine=false,coarse=false;
        for(const auto& cell:d->cells){const auto size=cell.upper-cell.lower;
            fine=fine||(size.x==virtualSpacing&&size.y==virtualSpacing);coarse=coarse||(size.x==2*virtualSpacing&&size.y==2*virtualSpacing);
            virtualCandidate=virtualCandidate&&((size.x==virtualSpacing&&size.y==virtualSpacing)||
                (size.x==2*virtualSpacing&&size.y==2*virtualSpacing));
            for(int axis=0;axis<2;++axis)for(double x:{coord(cell.lower,axis),coord(cell.upper,axis)}){
                const auto index=std::llround((x-coord(virtualOrigin,axis))/virtualSpacing);
                virtualCandidate=virtualCandidate&&coord(virtualOrigin,axis)+index*virtualSpacing==x;}}
        virtualCandidate=virtualCandidate&&fine&&coarse;
        std::vector<std::pair<OperatorTemplateLine,int>> boundaryRoutes;
        for(std::size_t u=0;u<d->faces.size();++u){const auto& face=d->faces[u];if(face.external)boundaryRoutes.push_back({{face.normalAxis,coord(face.lower,face.normalAxis)},static_cast<int>(u)});}
        const auto boundaryCounts=orderedLines(std::move(boundaryRoutes));
        for(const auto& count:boundaryCounts)virtualCandidate=virtualCandidate&&count.second.size()>=2;
    }
    if(compatible){
        const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
        d->compatible=true;d->nx=grid.nx();d->ny=grid.ny();d->logicalCells.resize(d->cells.size());
        for(auto& cell:d->cells)for(int axis=0;axis<2;++axis){
            for(int quarter=1;quarter<4;++quarter)
                cell.supportKnots[axis].push_back(coord(cell.lower,axis)+quarter*(coord(cell.upper-cell.lower,axis)/4));
            sortUnique(cell.supportKnots[axis]);
        }
        const auto lower=grid.minCoord(),spacing=grid.spacing();
        d->compatibleVelocity[0].assign((d->nx+1)*d->ny,-1);
        d->compatibleVelocity[1].assign(d->nx*(d->ny+1),-1);
        for(std::size_t p=0;p<d->cells.size();++p){
            const auto& c=d->cells[p];const int x=static_cast<int>(std::llround((c.lower.x-lower.x)/spacing.x)),
                y=static_cast<int>(std::llround((c.lower.y-lower.y)/spacing.y));
            require(x>=0&&x<d->nx&&y>=0&&y<d->ny,"invalid compatible cell map");d->logicalCells[p]={x,y};
        }
        for(std::size_t u=0;u<d->faces.size();++u){const auto& f=d->faces[u];
            const int x=static_cast<int>(std::llround((f.lower.x-lower.x)/spacing.x)),
                y=static_cast<int>(std::llround((f.lower.y-lower.y)/spacing.y));
            const int index=f.normalAxis==0?y*(d->nx+1)+x:y*d->nx+x;
            require(index>=0&&static_cast<std::size_t>(index)<d->compatibleVelocity[f.normalAxis].size(),"invalid compatible face map");
            require(d->compatibleVelocity[f.normalAxis][index]<0,"duplicate compatible face map");
            d->compatibleVelocity[f.normalAxis][index]=static_cast<int>(u);
        }
        for(const auto& ids:d->compatibleVelocity)for(int u:ids)require(u>=0,"incomplete compatible face map");
        d->compatibleBasis={compatibleIntervals(d->nx),compatibleIntervals(d->ny)};
        for(const auto& axis:d->compatibleBasis)for(const auto& interval:axis)for(const auto& half:interval)
            d->report.weights+=5*(half.nodes.size()+half.edges.size());
        templateBank->report(d->report);
        return HydroRegularAdvectorGeometry(std::move(d));
    }
    // Coplanar physical subfaces, partitioned at gaps in the selected domain.
    // Preparation alone owns the moment solves; apply only multiplies weights.
    using LineKey=std::pair<int,double>;
    std::vector<std::pair<LineKey,int>> faceRoutes;faceRoutes.reserve(d->faces.size());
    for(std::size_t u=0;u<d->faces.size();++u){const auto& f=d->faces[u];
        if(f.area>0)faceRoutes.push_back({{f.normalAxis,coord(f.lower,f.normalAxis)},static_cast<int>(u)});}
    auto lines=orderedLines(std::move(faceRoutes));OperatorTemplateSpans segments;segments.reserve(lines.size());
    for(auto& item:lines){auto& ids=item.second;const int t=1-item.first.first;segments.push_back({item.first,{}});
        std::sort(ids.begin(),ids.end(),[&](int a,int b){return coord(d->faces[a].lower,t)<coord(d->faces[b].lower,t);});
        int begin=0;for(int i=1;i<=static_cast<int>(ids.size());++i)
            if(i==static_cast<int>(ids.size())||coord(d->faces[ids[i-1]].upper,t)!=coord(d->faces[ids[i]].lower,t)){
                segments.back().second.emplace_back(begin,i);begin=i;}}
    for(std::size_t begin=0;begin<vertexCells.size();){std::size_t end=begin+1;const int vertex=vertexCells[begin].first;
        while(end<vertexCells.size()&&vertexCells[end].first==vertex)++end;d->usedVertices.push_back(vertex);
        if(d->vertexParent[vertex]<0){std::vector<int> incident;incident.reserve(end-begin);
            for(auto k=begin;k<end;++k)incident.push_back(vertexCells[k].second);
            d->datums.push_back({HydroRegularRow::Kind::Vertex,vertex,d->vertices[vertex],std::move(incident)});}begin=end;}
    std::vector<VertexUse>().swap(vertexCells);
    for(std::size_t e=0;e<d->edges.size();++e)if(d->edges[e].parent<0)d->datums.push_back({HydroRegularRow::Kind::Edge,static_cast<int>(e),(d->edges[e].lower+d->edges[e].upper)*.5,seeds[e]});
    for(std::size_t p=0;p<d->cells.size();++p)d->datums.push_back({HydroRegularRow::Kind::Center,static_cast<int>(p),(d->cells[p].lower+d->cells[p].upper)*.5,{static_cast<int>(p)}});
    // Late boundary binding selects a precomputed finite mode, never a fit.
    // Same exact axis/plane predicates, indexed once instead of E*F scans.
    // Insert in native-u order: every retained wall list keeps its old order.
    std::vector<std::pair<LineKey,int>> exteriorRoutes;
    for(std::size_t u=0;u<d->faces.size();++u){const auto& f=d->faces[u];
        if(f.external&&f.area>0)exteriorRoutes.push_back({{f.normalAxis,coord(f.lower,f.normalAxis)},static_cast<int>(u)});}
    const auto exteriorFacesByLine=orderedLines(std::move(exteriorRoutes));
    std::vector<std::pair<int,std::vector<int>>> wallEdges;
    std::vector<std::pair<LineKey,int>> wallEdgeRoutes;
    for(std::size_t e=0;e<d->edges.size();++e){const auto& edge=d->edges[e];if(edge.parent>=0)continue;std::vector<int> walls;const int a=edge.normalAxis,t=1-a;
        const auto found=findGroup(exteriorFacesByLine,LineKey{a,coord(edge.lower,a)});
        if(found!=exteriorFacesByLine.end())for(int u:found->second){const auto& f=d->faces[u];
            if(std::min(coord(f.upper,t),coord(edge.upper,t))>std::max(coord(f.lower,t),coord(edge.lower,t)))walls.push_back(u);}
        if(!walls.empty()){
            wallEdgeRoutes.push_back({{a,coord(edge.lower,a)},static_cast<int>(wallEdges.size())});
            wallEdges.emplace_back(static_cast<int>(e),std::move(walls));}}
    const auto wallEdgesByLine=orderedLines(std::move(wallEdgeRoutes));
    RegularRemapContext remap(virtualCandidate?nullptr:remapRequest,*d,neighbors,lines,segments,templateBank.get(),remapRequested);
    for(auto& datum:d->datums){Vector2D h{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
        for(int p:datum.seeds){const auto w=d->cells[p].upper-d->cells[p].lower;h.x=std::min(h.x,w.x);h.y=std::min(h.y,w.y);}
        for(int a=0;a<2;++a){const int t=1-a;const auto found=findGroup(wallEdgesByLine,LineKey{a,coord(datum.point,a)});
            if(found==wallEdgesByLine.end())continue;
            // Within each axis, records remain in the original global edge order.
            for(std::size_t index:found->second){const auto& item=wallEdges[index];const auto& e=d->edges[item.first];
                if(coord(datum.point,t)>=coord(e.lower,t) && coord(datum.point,t)<=coord(e.upper,t))
                    datum.wallFaces[a].insert(datum.wallFaces[a].end(),item.second.begin(),item.second.end());}}
        for(int component=0;component<2;++component)for(int mode=0;mode<2;++mode){const auto cols=columns(datum.point,component,rz,mode==1,degree);
            if(mode==1 && (datum.wallFaces[component].empty() || cols==columns(datum.point,component,rz,false,degree))){datum.wall[component]=datum.base[component];continue;}
            HydroRegularRow templateRow;std::optional<HydroRegularNormalPort> templatePort;RegularRowDependency dependency;
            const bool reused=remap.row(datum,h,component,mode,cols,templateRow,templatePort,dependency);
            if(!reused){++remap.report.recomputedRows;if(mode==0&&datum.kind==HydroRegularRow::Kind::Vertex)++remap.report.recomputedPorts;}
            if(reused||preparedOperatorTemplate(*d,datum,h,component,mode,rz,neighbors,lines,segments,*templateBank,templateRow,templatePort,&dependency)){
                if(templatePort){datum.normalPort[component]=static_cast<int>(d->portView().size());
                    ++d->report.normalPorts;d->report.normalPortWeights+=2*templatePort->support.size();
                    d->report.weights+=2*templatePort->support.size();d->portView().push_back(std::move(*templatePort));}
                d->report.weights+=4*templateRow.support.size();d->report.maxSupport=std::max(d->report.maxSupport,templateRow.support.size());
                d->report.maxCoefficientNorm=std::max(d->report.maxCoefficientNorm,templateRow.coefficientNorm);
                (mode==0?datum.base:datum.wall)[component]=static_cast<int>(d->rowView().size());d->rowView().push_back(std::move(templateRow));d->recordDependency(std::move(dependency));continue;
            }
            ++d->report.templateFallbacks;
            HydroRegularRow row;row.kind=datum.kind;row.index=datum.index;row.component=component;row.zeroNormal=mode==1;row.point=datum.point;row.scale=h;row.columns=cols;row.degree=degree;
            std::vector<int> support=datum.seeds;Eigen::MatrixXd m,inverse;ExactMainSvdCache::Selection svd(mainSvdCache);
            for(int ring=1;ring<=4;++ring){expandSupport(support,neighbors);
                std::vector<int> faces;for(int p:support)for(int u:d->cells[p].faces)if(d->faces[u].normalAxis==component && d->faces[u].area>0)faces.push_back(u);uniqueIndices(faces);
                require(!faces.empty(),"no native component observations");row.support.assign(faces.begin(),faces.end());m.resize(faces.size(),cols.size());
                for(std::size_t i=0;i<row.support.size();++i){const auto moments=mean(d->faces[row.support[i]],datum.point,h,rz,degree);for(std::size_t j=0;j<cols.size();++j)m(i,j)=moments[cols[j]];}
                svd.compute(m,Eigen::ComputeThinU|Eigen::ComputeFullV,
                    std::max<std::size_t>(row.support.size(),degree==2?6:10),d->report.svdCount);
                svd.setThreshold(128*std::numeric_limits<double>::epsilon()*std::max<std::size_t>(row.support.size(),degree==2?6:10));row.rank=static_cast<int>(svd.rank());row.rings=ring;
                bool exhausted=true;for(int p:support)for(int q:neighbors[p])if(!std::binary_search(support.begin(),support.end(),q))exhausted=false;
                if(row.rank==static_cast<int>(cols.size()) || ring==4 || exhausted)break;
            }
            svd.retain();
            inverse=Eigen::MatrixXd::Zero(cols.size(),row.support.size());
            for(int k=0;k<row.rank;++k)inverse.noalias()+=svd.matrixV().col(k)*(svd.matrixU().col(k).transpose()/svd.singularValues()[k]);

            const Eigen::MatrixXd lower=prioritySpace(row,rz),observed=m*lower;
            row.priorityDimension=static_cast<int>(lower.cols());row.priorityRank=row.priorityDimension;
            if(row.rank<static_cast<int>(cols.size()) && lower.cols()>0){
                Eigen::JacobiSVD<Eigen::MatrixXd> low(observed,Eigen::ComputeThinU|Eigen::ComputeThinV);++d->report.svdCount;
                low.setThreshold(128*std::numeric_limits<double>::epsilon()*std::max<std::size_t>(row.support.size(),6));row.priorityRank=static_cast<int>(low.rank());
                if(row.priorityRank==row.priorityDimension){
                    Eigen::MatrixXd right=Eigen::MatrixXd::Zero(lower.cols(),row.support.size());
                    for(int k=0;k<row.priorityRank;++k)right.noalias()+=low.matrixV().col(k)*(low.matrixU().col(k).transpose()/low.singularValues()[k]);
                    // Preserve the already accepted fully identifiable weights.
                    const Eigen::MatrixXd correction=(lower-inverse*observed)*right;
                    inverse+=correction;
                }else{
                    const double threshold=low.threshold()*low.singularValues()[0];
                    const Eigen::MatrixXd hierarchy=identifiableInverse(row,rz,m,threshold,d->report.svdCount);
                    // Observation-null correction preserves the full least-squares
                    // fit while reproducing the identifiable degree hierarchy.
                    const Eigen::MatrixXd correction=hierarchy-inverse*(m*hierarchy);
                    inverse+=correction;
                }
            }
            row.p2Rank=row.rank;row.p2ColumnCount=static_cast<int>(cols.size());
            row.p2Deficit=(Eigen::MatrixXd::Identity(cols.size(),cols.size())-inverse*m).norm();
        
            row.priorityDeficit=(inverse*observed-lower).norm();
            row.coefficientNorm=inverse.norm();row.deficit=(Eigen::MatrixXd::Identity(cols.size(),cols.size())-inverse*m).norm();
            if(row.rank<static_cast<int>(cols.size())){++d->report.rankDeficient;const auto null=svd.matrixV().col(cols.size()-1);for(std::size_t j=0;j<cols.size();++j)row.nullPolynomial[cols[j]]=null[j];row.nullResidual=(m*null).cwiseAbs().maxCoeff();}
            constexpr int outputs[4]={0,1,2,4};const double scales[4]={1,h.x,h.y,h.x*h.y};
            for(int k=0;k<4;++k){auto& w=row.weights[k];w.assign(row.support.size(),0);const auto it=std::find(cols.begin(),cols.end(),outputs[k]);
                if(it!=cols.end())for(std::size_t i=0;i<w.size();++i){w[i]=inverse(it-cols.begin(),i)/scales[k];require(std::isfinite(w[i]),"nonfinite prepared weight");}}
            if(mode==0&&datum.kind==HydroRegularRow::Kind::Vertex){
                HydroRegularNormalPort port;port.vertex=datum.index;port.component=component;
                const int tangent=1-component;const double target=coord(datum.point,tangent),scale=coord(h,tangent);
                // Axis radial normal velocity and its axial derivative vanish.
                if(!(rz&&component==0&&datum.point.x==0)){
                    const LineKey key{component,coord(datum.point,component)};const auto found=findGroup(lines,key);
                    require(found!=lines.end(),"missing normal-port coordinate line");const auto& ids=found->second;
                    const auto& spans=groupValues(segments,key);
                    const auto span=std::lower_bound(spans.begin(),spans.end(),target,[&](const auto& s,double x){return coord(d->faces[ids[s.second-1]].upper,tangent)<x;});
                    require(span!=spans.end()&&coord(d->faces[ids[span->first]].lower,tangent)<=target,"normal port outside contiguous native span");
                    const auto position=std::lower_bound(ids.begin()+span->first,ids.begin()+span->second,target,[&](int u,double x){return coord(d->faces[u].lower+d->faces[u].upper,tangent)/2<x;});
                    int left=static_cast<int>(position-ids.begin())-1,right=left+1;
                    std::vector<int> observations;
                    while(observations.size()<4&&(left>=span->first||right<span->second)){
                        if(right==span->second||(left>=span->first&&target-coord(d->faces[ids[left]].lower+d->faces[ids[left]].upper,tangent)/2<=coord(d->faces[ids[right]].lower+d->faces[ids[right]].upper,tangent)/2-target))observations.push_back(ids[left--]);
                        else observations.push_back(ids[right++]);}
                    const bool evenAxis=rz&&component==1&&datum.point.x==0;
                    const std::vector<int> powers=evenAxis?std::vector<int>{0,2}:std::vector<int>{0,1,2,3};
                    if(observations.size()>powers.size())observations.resize(powers.size());
                    const int count=static_cast<int>(observations.size()),degree=static_cast<int>(powers.size());
                    Eigen::MatrixXd moments=Eigen::MatrixXd::Zero(count,degree);
                    for(int i=0;i<count;++i){const auto& f=d->faces[observations[i]];
                        for(int q=0;q<5;++q){const auto x=f.lower+(f.upper-f.lower)*gauss[q];const double z=(coord(x,tangent)-target)/scale;
                            const double weight=gaussWeight[q]*(rz?x.x/((f.lower.x+f.upper.x)/2):1);
                            for(int k=0;k<degree;++k)moments(i,k)+=weight*std::pow(z,powers[k]);}}
                    Eigen::JacobiSVD<Eigen::MatrixXd> portSvd(moments.leftCols(count),Eigen::ComputeThinU|Eigen::ComputeThinV);++d->report.svdCount;
                    portSvd.setThreshold(128*std::numeric_limits<double>::epsilon()*std::max(count,6));
                    require(portSvd.rank()==count,"deficient native line moments");port.observedModes=count;
                    const Eigen::MatrixXd rightInverse=portSvd.solve(Eigen::MatrixXd::Identity(count,count));
                    // Missing line modes use the same selected-degree seed. Observable
                    // line moments remain exact; P2 has no cubic coefficient.
                    std::map<int,std::array<double,2>> sparse;
                    for(int i=0;i<count;++i){auto& w=sparse[observations[i]];w[0]+=rightInverse(0,i);
                        if(!evenAxis&&count>1)w[1]+=rightInverse(1,i)/scale;}
                    for(int k=count;k<degree;++k){const int power=powers[k];if(power>row.degree)continue;
                        const int polynomial=power==0?0:power==1?1+tangent:power==2?(tangent==0?3:5):(tangent==0?6:9);
                        const auto col=std::find(cols.begin(),cols.end(),polynomial);if(col==cols.end())continue;
                        const Eigen::VectorXd correction=-rightInverse*moments.col(k);
                        for(std::size_t i=0;i<row.support.size();++i){const double seed=inverse(col-cols.begin(),i);auto& w=sparse[row.support[i]];
                            w[0]+=correction[0]*seed;
                            if(!evenAxis)w[1]+=(count>1?correction[1]:(power==1?1.:0.))*seed/scale;}}
                    for(const auto& item:sparse){port.support.push_back(item.first);for(int k=0;k<2;++k){require(std::isfinite(item.second[k]),"nonfinite native port weight");port.weights[k].push_back(item.second[k]);}}
                }
                datum.normalPort[component]=static_cast<int>(d->portView().size());
                ++d->report.normalPorts;d->report.normalPortWeights+=2*port.support.size();
                d->report.weights+=2*port.support.size();d->portView().push_back(std::move(port));
            }
            d->report.weights+=4*row.support.size();d->report.maxSupport=std::max(d->report.maxSupport,row.support.size());d->report.maxCoefficientNorm=std::max(d->report.maxCoefficientNorm,row.coefficientNorm);
            (mode==0?datum.base:datum.wall)[component]=static_cast<int>(d->rowView().size());d->rowView().push_back(std::move(row));d->recordDependency(RegularRowDependency{});
        }
    }
    d->datumIndex.reserve(d->datums.size());
    for(std::size_t i=0;i<d->datums.size();++i)d->datumIndex.emplace_back(std::make_pair(d->datums[i].kind,d->datums[i].index),i);
    std::sort(d->datumIndex.begin(),d->datumIndex.end());
    if(virtualCandidate){
        const auto centerRow=[&](int pressure,int component)->const HydroRegularRow&{
            const auto& datum=d->datums[d->datum(HydroRegularRow::Kind::Center,pressure)];
            return d->rowView()[datum.base[component]];};
        const auto weights=[](const HydroRegularRow& row,int kind){HydroVirtualRow result;
            for(std::size_t k=0;k<row.support.size();++k)if(row.weights[kind][k]!=0)
                result.push_back({row.support[k],row.weights[kind][k]});return result;};
        const auto index=[&](double coordinate,int axis){return static_cast<int>(std::llround((coordinate-coord(virtualOrigin,axis))/virtualSpacing));};
        std::vector<HydroVirtualCellInput> inputs(d->cells.size());
        for(std::size_t p=0;p<d->cells.size();++p){const auto& cell=d->cells[p];auto& input=inputs[p];
            for(int axis=0;axis<2;++axis){input.lower[axis]=index(coord(cell.lower,axis),axis);input.upper[axis]=index(coord(cell.upper,axis),axis);}
            for(int c=0;c<2;++c){const auto& row=centerRow(static_cast<int>(p),c);input.value[c]=weights(row,0);
                for(int axis=0;axis<2;++axis)input.gradient[c][axis]=weights(row,1+axis);}}
        std::vector<std::pair<OperatorTemplateLine,int>> boundaryRoutes;
        for(std::size_t u=0;u<d->faces.size();++u){const auto& f=d->faces[u];if(f.external)
            boundaryRoutes.push_back({{f.normalAxis,coord(f.lower,f.normalAxis)},static_cast<int>(u)});}
        auto boundaryLines=orderedLines(std::move(boundaryRoutes));
        for(auto& line:boundaryLines){const int tangent=1-line.first.first;
            std::sort(line.second.begin(),line.second.end(),[&](int a,int b){return coord(d->faces[a].lower,tangent)<coord(d->faces[b].lower,tangent);});}
        std::vector<HydroVirtualFaceInput> faceInputs;faceInputs.reserve(d->faces.size());
        for(std::size_t u=0;u<d->faces.size();++u){const auto& f=d->faces[u];HydroVirtualFaceInput input;
            input.velocityDof=static_cast<int>(u);input.normalAxis=f.normalAxis;input.negativePressure=f.negativePressure;input.positivePressure=f.positivePressure;
            for(int axis=0;axis<2;++axis){input.lower[axis]=index(coord(f.lower,axis),axis);input.upper[axis]=index(coord(f.upper,axis),axis);}
            const int tangent=1-f.normalAxis;
            if(input.upper[tangent]-input.lower[tangent]==2){
                if(f.external){const auto& ids=groupValues(boundaryLines,OperatorTemplateLine{f.normalAxis,coord(f.lower,f.normalAxis)});
                    const auto found=std::find(ids.begin(),ids.end(),static_cast<int>(u));require(found!=ids.end(),"missing virtual boundary line incidence");
                    const auto offset=static_cast<std::size_t>(found-ids.begin());const int a=offset==0?ids[0]:ids[offset-1],b=offset==0?ids[1]:ids[offset];
                    const double xa=coord((d->faces[a].lower+d->faces[a].upper)*.5,tangent),xb=coord((d->faces[b].lower+d->faces[b].upper)*.5,tangent);
                    require(xb>xa,"invalid virtual boundary derivative span");input.tangentGradient={{a,-1/(xb-xa)},{b,1/(xb-xa)}};
                }else {const int owner=f.negativePressure>=0?f.negativePressure:f.positivePressure;
                    input.tangentGradient=weights(centerRow(owner,f.normalAxis),1+tangent);}}
            faceInputs.push_back(std::move(input));}
        d->virtualProlongation=prepareHydroVirtualFaceProlongation(virtualNx,virtualNy,{virtualSpacing,virtualSpacing},l.velocityCount(),inputs,faceInputs);
        d->compatible=true;d->nx=virtualNx;d->ny=virtualNy;
        for(int axis=0;axis<2;++axis){auto& lines=d->virtualLines[axis];const int count=axis==0?d->nx:d->ny;
            lines.reserve(count+1);for(int k=0;k<=count;++k)lines.push_back(coord(virtualOrigin,axis)+k*virtualSpacing);}
        d->virtualCells.resize(static_cast<std::size_t>(d->nx)*d->ny);d->logicalCells.resize(d->virtualCells.size());
        for(int y=0;y<d->ny;++y)for(int x=0;x<d->nx;++x){const auto p=static_cast<std::size_t>(y)*d->nx+x;auto& cell=d->virtualCells[p];
            cell.lower={d->virtualLines[0][x],d->virtualLines[1][y]};cell.upper={d->virtualLines[0][x+1],d->virtualLines[1][y+1]};
            cell.nativeCell=d->cells[d->virtualProlongation->cellPressure[p]].nativeCell;d->logicalCells[p]={x,y};}
        d->compatibleVelocity[0].resize(d->virtualProlongation->rows[0].size());d->compatibleVelocity[1].resize(d->virtualProlongation->rows[1].size());
        int slot=0;for(auto& axis:d->compatibleVelocity)for(int& u:axis)u=slot++;
        d->compatibleBasis={compatibleIntervals(d->nx),compatibleIntervals(d->ny)};
        for(auto& cell:d->cells)for(int axis=0;axis<2;++axis){auto& knots=cell.supportKnots[axis];knots.clear();
            const int quarters=4*static_cast<int>(std::llround(coord(cell.upper-cell.lower,axis)/virtualSpacing));
            for(int k=0;k<=quarters;++k)knots.push_back(k==quarters?coord(cell.upper,axis):coord(cell.lower,axis)+k*(virtualSpacing/4));}
        d->report.virtualCells=d->virtualCells.size();d->report.virtualFaces=static_cast<std::size_t>(slot);
        for(const auto& axis:d->virtualProlongation->rows)for(const auto& row:axis)d->report.prolongationWeights+=row.size();
        d->report.weights+=d->report.prolongationWeights;
    }
    templateBank->report(d->report);
    retainRegularSourcePreparation(*d,neighbors,lines,segments,templateBank.get());
    if(remapRequested){OperatorTemplateEnvironment restore;d->geometry.publishRegularPreparationReuseReport(remap.report);}
    return HydroRegularAdvectorGeometry(std::move(d));
}

HydroRegularCoefficients HydroRegularAdvectorGeometry::apply(const QuadtreeMesh& mesh,const HydroNativeValues& values,const HydroAdvectorBoundary* boundary) const
{
    fp();const auto& d=*data_;values.validateFor(d.geometry.layout(),mesh);if(boundary)boundary->validateFor(mesh,d.geometry);
    const auto& velocity=values.normalVelocity();std::vector<bool> wall(d.faces.size(),false);
    if(boundary)for(int u:boundary->velocityDofs()){require(velocity[u]==0,"nonzero prescribed wall velocity");wall[u]=true;}
    for(std::size_t u=0;u<d.faces.size();++u)if(d.faces[u].area==0)require(velocity[u]==0,"nonzero axis velocity");
    if(d.compatible){
        require(!boundary||boundary->completeBoundaryLines(),"partial impermeable boundary requires boundary-aware source preparation");
        std::vector<double> prolonged;
        if(d.virtualProlongation){prolonged.resize(d.virtualProlongation->rows[0].size()+d.virtualProlongation->rows[1].size());
            const auto count=static_cast<int>(prolonged.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if(count>=256)
#endif
            for(int u=0;u<count;++u){long double value=0;for(const auto& w:virtualRow(d,u))value+=static_cast<long double>(w.weight)*velocity[w.velocityDof];
                prolonged[u]=static_cast<double>(value);}
            for(double value:prolonged)require(std::isfinite(value),"nonfinite prolonged virtual mean");}
        const auto& tensorVelocity=d.virtualProlongation?prolonged:velocity;
        const auto tensorCount=d.virtualProlongation?d.virtualCells.size():d.cells.size();
        HydroRegularCoefficients out;out.geometry_=data_;out.cellDofs.resize(16*tensorCount);out.faceCorrections.assign(d.faces.size(),0);out.compatibleNativeMeans=tensorVelocity;
        const auto count=static_cast<int>(tensorCount);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if(count>=256)
#endif
        for(int p=0;p<count;++p){const auto logical=d.logicalCells[p];
            for(int hy=0;hy<4;++hy)for(int hx=0;hx<4;++hx){auto& target=out.cellDofs[16*p+4*hy+hx];
                const auto& bx=d.compatibleBasis[0][logical[0]][hx];const auto& by=d.compatibleBasis[1][logical[1]][hy];
                for(int c=0;c<2;++c){const auto& xx=c==0?bx.nodes:bx.edges;const auto& yy=c==0?by.edges:by.nodes;
                    std::array<long double,25> polynomial{};
                    for(const auto& y:yy){std::array<long double,5> row{};
                        for(const auto& x:xx){const int u=d.compatibleVelocity[c][c==0?y.index*(d.nx+1)+x.index:y.index*d.nx+x.index];
                            for(int k=0;k<5;++k)row[k]+=tensorVelocity[u]*x.powers[k];}
                        for(int j=0;j<5;++j)for(int i=0;i<5;++i)polynomial[5*j+i]+=row[i]*y.powers[j];
                    }
                    for(int k=0;k<25;++k)target[c][k]=static_cast<double>(polynomial[k]);
                }
            }
        }
        for(const auto& cell:out.cellDofs)for(const auto& component:cell)for(double value:component)
            require(std::isfinite(value),"nonfinite compatible polynomial");
        return out;
    }
    using Jet=std::array<std::array<double,4>,2>;
    std::vector<Jet> nodes(d.vertices.size());
    std::vector<std::array<std::array<double,2>,2>> mid(d.edges.size());
    std::vector<std::array<double,2>> centers(d.cells.size());
    for(const auto& datum:d.datums)for(int c=0;c<2;++c){bool prescribed=false;for(int f:datum.wallFaces[c])prescribed=prescribed||wall[f];
        const auto& row=d.rowView()[prescribed?datum.wall[c]:datum.base[c]];std::array<double,4> value{};
        for(int k=0;k<4;++k)for(std::size_t i=0;i<row.support.size();++i)value[k]+=row.weights[k][i]*velocity[row.support[i]];
        if(datum.kind==HydroRegularRow::Kind::Vertex){
            const auto& port=d.portView()[datum.normalPort[c]];value[0]=0;value[2-c]=0;
            if(!prescribed)for(std::size_t i=0;i<port.support.size();++i){value[0]+=port.weights[0][i]*velocity[port.support[i]];value[2-c]+=port.weights[1][i]*velocity[port.support[i]];}
            nodes[datum.index][c]=value;}
        else if(datum.kind==HydroRegularRow::Kind::Edge)mid[datum.index][c]={value[0],value[1+d.edges[datum.index].normalAxis]};
        else centers[datum.index][c]=value[0];}
    std::vector<unsigned char> nodeState(d.vertices.size()),edgeState(d.edges.size());
    std::function<void(int)> node,edge;
    std::function<Jet(int,double)> edgeJet=[&](int e,double s){edge(e);const auto& ge=d.edges[e];const int a=ge.vertices[0],b=ge.vertices[1],n=ge.normalAxis,t=1-n;
        node(a);node(b);const double h=coord(ge.upper-ge.lower,t);const auto bv=basis(s),bd=basis(s,true);Jet result{};
        for(int c=0;c<2;++c){const std::array<double,5> v{nodes[a][c][0],nodes[b][c][0],h*nodes[a][c][1+t],h*nodes[b][c][1+t],mid[e][c][0]};
            const std::array<double,5> dn{nodes[a][c][1+n],nodes[b][c][1+n],h*nodes[a][c][3],h*nodes[b][c][3],mid[e][c][1]};
            for(int k=0;k<5;++k){result[c][0]+=bv[k]*v[k];result[c][1+t]+=bd[k]*v[k]/h;result[c][1+n]+=bv[k]*dn[k];result[c][3]+=bd[k]*dn[k]/h;}}
        return result;};
    node=[&](int v){if(nodeState[v]==2)return;require(nodeState[v]==0,"cyclic active node constraint");nodeState[v]=1;const int e=d.vertexParent[v];
        if(e>=0){const auto& ge=d.edges[e];const int t=1-ge.normalAxis;nodes[v]=edgeJet(e,(coord(d.vertices[v],t)-coord(ge.lower,t))/coord(ge.upper-ge.lower,t));}nodeState[v]=2;};
    edge=[&](int e){if(edgeState[e]==2)return;require(edgeState[e]==0,"cyclic active edge constraint");edgeState[e]=1;const auto& ge=d.edges[e];
        if(ge.parent>=0){const auto& master=d.edges[ge.parent];const int t=1-ge.normalAxis;const auto jet=edgeJet(ge.parent,(coord((ge.lower+ge.upper)*.5,t)-coord(master.lower,t))/coord(master.upper-master.lower,t));
            for(int c=0;c<2;++c)mid[e][c]={jet[c][0],jet[c][1+ge.normalAxis]};}edgeState[e]=2;};
    for(int v:d.usedVertices)node(v);for(std::size_t e=0;e<d.edges.size();++e)edge(static_cast<int>(e));
    HydroRegularCoefficients out;out.geometry_=data_;out.cellDofs.resize(d.cells.size());out.faceCorrections.resize(d.faces.size());
    constexpr int corners[4][2]={{0,0},{1,0},{1,1},{0,1}},mids[4][2]={{4,0},{4,1},{0,4},{1,4}};
    for(std::size_t p=0;p<d.cells.size();++p){const auto& cell=d.cells[p];const auto h=cell.upper-cell.lower;
        for(int c=0;c<2;++c){auto& raw=out.cellDofs[p][c];for(int k=0;k<4;++k){const int x=corners[k][0],y=corners[k][1];const auto& j=nodes[cell.vertices[k]][c];
                raw[5*y+x]=j[0];raw[5*y+x+2]=h.x*j[1];raw[5*(y+2)+x]=h.y*j[2];raw[5*(y+2)+x+2]=h.x*h.y*j[3];}
            for(int k=0;k<4;++k){const int x=mids[k][0],y=mids[k][1];const auto& z=mid[cell.edges[k]][c];raw[5*y+x]=z[0];if(k<2)raw[5*(y+2)+x]=z[1]*h.y;else raw[5*y+x+2]=z[1]*h.x;}
            raw[24]=centers[p][c];for(double v:raw)require(std::isfinite(v),"nonfinite Q4 datum");}}
    for(std::size_t u=0;u<d.faces.size();++u){const auto& f=d.faces[u];if(f.area==0)continue;const int p=f.negativePressure>=0?f.negativePressure:f.positivePressure;
        double meanValue=0;for(int k=0;k<25;++k)meanValue+=f.meanWeights[k]*out.cellDofs[p][f.normalAxis][k];out.faceCorrections[u]=15./8*(velocity[u]-meanValue);
        require(std::isfinite(out.faceCorrections[u]),"nonfinite native mean correction");}
    return out;
}

namespace
{
void inside(const HydroRegularCell& cell,Vector2D point)
{
    require(std::isfinite(point.x)&&std::isfinite(point.y)&&point.x>=cell.lower.x&&point.x<=cell.upper.x
        &&point.y>=cell.lower.y&&point.y<=cell.upper.y,"point outside regular cell");
}
double supportWidth(const HydroRegularFace& face,int pressure)
{
    require(pressure==face.negativePressure||pressure==face.positivePressure,"missing regular face support");
    const auto width=face.supportWidth[pressure==face.positivePressure?1:0];
    require(std::isfinite(width)&&width>0,"invalid regular face support width");return width;
}
// Exact algebraic restriction q(origin+scale*t), never a sampled fit.
std::array<long double,5> compose(const double* power,long double origin,long double scale)
{
    std::array<long double,5> out{};
    for(int k=4;k>=0;--k){for(int j=4;j>0;--j)out[j]=origin*out[j]+scale*out[j-1];out[0]=origin*out[0]+power[k];}
    return out;
}
}
HydroRegularCellStencil HydroRegularAdvectorGeometry::compatibleCellStencil(int pressure,const std::vector<int>* boundaryDofs) const
{
    const auto& d=*data_;HydroRegularCellStencil out;out.pressure=pressure;out.geometry_=data_;
    if(d.virtualProlongation){if(boundaryDofs)out.excludedVelocityDofs_=*boundaryDofs;return out;}
    if(boundaryDofs){const auto logical=d.logicalCells[pressure];
        const auto& bx=d.compatibleBasis[0][logical[0]][0];const auto& by=d.compatibleBasis[1][logical[1]][0];
        for(int c=0;c<2;++c){const auto& xx=c==0?bx.nodes:bx.edges;const auto& yy=c==0?by.edges:by.nodes;
            for(const auto& y:yy)for(const auto& x:xx){const int u=d.compatibleVelocity[c][c==0?y.index*(d.nx+1)+x.index:y.index*d.nx+x.index];
                if(std::binary_search(boundaryDofs->begin(),boundaryDofs->end(),u))out.excludedVelocityDofs_.push_back(u);}}
        std::sort(out.excludedVelocityDofs_.begin(),out.excludedVelocityDofs_.end());
        out.excludedVelocityDofs_.erase(std::unique(out.excludedVelocityDofs_.begin(),out.excludedVelocityDofs_.end()),out.excludedVelocityDofs_.end());
    }
    return out;
}
std::vector<HydroRegularCellStencil> HydroRegularAdvectorGeometry::cellStencils(const QuadtreeMesh& mesh,
    const std::vector<int>& pressures,const HydroAdvectorBoundary* boundary) const
{
    fp();const auto& d=*data_;d.geometry.validateMesh(mesh);if(boundary)boundary->validateFor(mesh,d.geometry);
    require(d.compatible,"compatible stencil batch requires the prepared compatible tensor source");
    require(!boundary||boundary->completeBoundaryLines(),"partial impermeable boundary requires boundary-aware source preparation");
    std::vector<HydroRegularCellStencil> result;result.reserve(pressures.size());
    for(int p:pressures){require(p>=0&&static_cast<std::size_t>(p)<d.cells.size(),"stencil pressure index out of range");
        result.push_back(compatibleCellStencil(p,boundary?&boundary->velocityDofs():nullptr));}
    return result;
}
HydroRegularCellStencil HydroRegularAdvectorGeometry::cellStencil(const QuadtreeMesh& mesh,int pressure,const HydroAdvectorBoundary* boundary) const
{
    fp();const auto& d=*data_;d.geometry.validateMesh(mesh);if(boundary)boundary->validateFor(mesh,d.geometry);
    return cellStencilPrepared(pressure,boundary);
}
HydroRegularCellStencil HydroRegularAdvectorGeometry::cellStencilPrepared(int pressure,const HydroAdvectorBoundary* boundary) const
{
    fp();const auto& d=*data_;
    require(pressure>=0&&static_cast<std::size_t>(pressure)<d.cells.size(),"stencil pressure index out of range");
    if(d.compatible){
        require(!boundary||boundary->completeBoundaryLines(),"partial impermeable boundary requires boundary-aware source preparation");
        return compatibleCellStencil(pressure,boundary?&boundary->velocityDofs():nullptr);
    }
    using Sparse=std::map<int,double>;
    using Jet=std::array<std::array<Sparse,4>,2>;
    using Mid=std::array<std::array<Sparse,2>,2>;
    using Raw=std::array<std::array<Sparse,25>,2>;
    const auto add=[](Sparse& into,const Sparse& from,double factor){if(factor!=0)for(const auto& w:from)into[w.first]+=factor*w.second;};
    std::set<int> wall;if(boundary)wall.insert(boundary->velocityDofs().begin(),boundary->velocityDofs().end());
    const auto allowed=[&](int u){return d.faces[u].area>0&&wall.count(u)==0;};
    const auto datumJet=[&](HydroRegularRow::Kind kind,int index){Jet jet{};const auto& datum=d.datums[d.datum(kind,index)];
        for(int c=0;c<2;++c){bool prescribed=false;for(int f:datum.wallFaces[c])prescribed=prescribed||wall.count(f)!=0;
            const auto& row=d.rowView()[prescribed?datum.wall[c]:datum.base[c]];
            for(int k=0;k<4;++k)for(std::size_t i=0;i<row.support.size();++i)if(allowed(row.support[i]))jet[c][k][row.support[i]]+=row.weights[k][i];
            if(kind==HydroRegularRow::Kind::Vertex){const auto& port=d.portView()[datum.normalPort[c]];jet[c][0].clear();jet[c][2-c].clear();
                if(!prescribed)for(std::size_t i=0;i<port.support.size();++i)if(allowed(port.support[i])){
                    jet[c][0][port.support[i]]+=port.weights[0][i];jet[c][2-c][port.support[i]]+=port.weights[1][i];}}
        }return jet;};
    std::map<int,Jet> nodes;std::map<int,Mid> mids;std::map<int,Raw> rawCells;
    std::set<int> pendingNodes,pendingEdges;
    std::function<const Jet&(int)> node;std::function<const Mid&(int)> edge;
    const auto edgeJet=[&](int e,double s){const auto& mid=edge(e);const auto& ge=d.edges[e];
        const auto& a=node(ge.vertices[0]);const auto& b=node(ge.vertices[1]);const int n=ge.normalAxis,t=1-n;
        const double h=coord(ge.upper-ge.lower,t);const auto bv=basis(s),bd=basis(s,true);Jet result{};
        for(int c=0;c<2;++c){const Sparse* values[5]={&a[c][0],&b[c][0],&a[c][1+t],&b[c][1+t],&mid[c][0]};
            const Sparse* normal[5]={&a[c][1+n],&b[c][1+n],&a[c][3],&b[c][3],&mid[c][1]};
            for(int k=0;k<5;++k){const double scale=(k==2||k==3)?h:1;
                add(result[c][0],*values[k],bv[k]*scale);add(result[c][1+t],*values[k],bd[k]*scale/h);
                add(result[c][1+n],*normal[k],bv[k]*scale);add(result[c][3],*normal[k],bd[k]*scale/h);}}
        return result;};
    node=[&](int v)->const Jet&{const auto found=nodes.find(v);if(found!=nodes.end())return found->second;
        require(pendingNodes.insert(v).second,"cyclic local node constraint");Jet jet{};const int e=d.vertexParent[v];
        if(e<0)jet=datumJet(HydroRegularRow::Kind::Vertex,v);
        else {const auto& ge=d.edges[e];const int t=1-ge.normalAxis;jet=edgeJet(e,(coord(d.vertices[v],t)-coord(ge.lower,t))/coord(ge.upper-ge.lower,t));}
        pendingNodes.erase(v);return nodes.emplace(v,std::move(jet)).first->second;};
    edge=[&](int e)->const Mid&{const auto found=mids.find(e);if(found!=mids.end())return found->second;
        require(pendingEdges.insert(e).second,"cyclic local edge constraint");const auto& ge=d.edges[e];Mid mid{};
        if(ge.parent<0){const auto jet=datumJet(HydroRegularRow::Kind::Edge,e);for(int c=0;c<2;++c)mid[c]={jet[c][0],jet[c][1+ge.normalAxis]};}
        else {const auto& master=d.edges[ge.parent];const int t=1-ge.normalAxis;
            const auto jet=edgeJet(ge.parent,(coord((ge.lower+ge.upper)*.5,t)-coord(master.lower,t))/coord(master.upper-master.lower,t));
            for(int c=0;c<2;++c)mid[c]={jet[c][0],jet[c][1+ge.normalAxis]};}
        pendingEdges.erase(e);return mids.emplace(e,std::move(mid)).first->second;};
    const auto rawCell=[&](int p)->const Raw&{const auto found=rawCells.find(p);if(found!=rawCells.end())return found->second;
        const auto& cell=d.cells[p];const auto h=cell.upper-cell.lower;Raw raw{};
        constexpr int corners[4][2]={{0,0},{1,0},{1,1},{0,1}},positions[4][2]={{4,0},{4,1},{0,4},{1,4}};
        for(int c=0;c<2;++c){for(int k=0;k<4;++k){const int x=corners[k][0],y=corners[k][1];const auto& j=node(cell.vertices[k])[c];
                add(raw[c][5*y+x],j[0],1);add(raw[c][5*y+x+2],j[1],h.x);add(raw[c][5*(y+2)+x],j[2],h.y);add(raw[c][5*(y+2)+x+2],j[3],h.x*h.y);}
            for(int k=0;k<4;++k){const int x=positions[k][0],y=positions[k][1];const auto& z=edge(cell.edges[k])[c];
                add(raw[c][5*y+x],z[0],1);if(k<2)add(raw[c][5*(y+2)+x],z[1],h.y);else add(raw[c][5*y+x+2],z[1],h.x);}
            raw[c][24]=datumJet(HydroRegularRow::Kind::Center,p)[c][0];}
        return rawCells.emplace(p,std::move(raw)).first->second;};
    const auto pack=[](const Sparse& values){std::vector<HydroRegularScalarWeight> result;for(const auto& w:values)if(w.second!=0){require(std::isfinite(w.second),"nonfinite local source stencil");result.push_back({w.first,w.second});}return result;};
    HydroRegularCellStencil out;out.pressure=pressure;out.geometry_=data_;const auto& raw=rawCell(pressure);
    for(int c=0;c<2;++c)for(int k=0;k<25;++k)out.cellDofs[c][k]=pack(raw[c][k]);
    for(int u:d.cells[pressure].faces){const auto& f=d.faces[u];Sparse correction;
        if(f.area>0){const int p=f.negativePressure>=0?f.negativePressure:f.positivePressure;const auto& owner=rawCell(p);
            if(allowed(u))correction[u]=15./8;for(int k=0;k<25;++k)add(correction,owner[f.normalAxis][k],-15./8*f.meanWeights[k]);}
        out.faceCorrections.push_back(pack(correction));}
    return out;
}
std::vector<HydroRegularJetWeight> HydroRegularAdvectorGeometry::cellJetWeights(const HydroRegularCellStencil& stencil,Vector2D point) const
{
    fp();require(stencil.geometry_==data_,"different local stencil geometry identity");
    require(stencil.pressure>=0&&static_cast<std::size_t>(stencil.pressure)<cells().size(),"stencil pressure index out of range");
    const auto& cell=cells()[stencil.pressure];inside(cell,point);
    if(data_->compatible){
        const auto& d=*data_;const auto patch=preparedCompatiblePatch(d,stencil.pressure,point);const int tensorPressure=patch.index/16;
        const auto logical=d.logicalCells[tensorPressure];const auto& sourceCell=tensorCell(d,tensorPressure);
        const int hx=(patch.index-16*tensorPressure)%4,hy=(patch.index-16*tensorPressure)/4;
        const auto h=patch.upper-patch.lower;const double sx=(point.x-patch.lower.x)/h.x,sy=(point.y-patch.lower.y)/h.y;
        const auto evaluate=[](const std::array<double,5>& a,double z){std::array<double,2> result{};
            for(int k=4;k>=0;--k){result[1]=result[1]*z+result[0];result[0]=result[0]*z+a[k];}return result;};
        const auto& bx=d.compatibleBasis[0][logical[0]][hx];const auto& by=d.compatibleBasis[1][logical[1]][hy];
        std::shared_ptr<VirtualJetScatterSlot> scatter;
        if(d.virtualProlongation){
            {std::lock_guard<std::mutex> lock(d.virtualJetScatterMutex);
                auto& slot=d.virtualJetScatter[{tensorPressure,stencil.excludedVelocityDofs_}];
                if(!slot)slot=std::make_shared<VirtualJetScatterSlot>();scatter=slot;}
            std::call_once(scatter->ready,[&]{
                VirtualJetScatterPlan plan;std::set<int> nativeVelocity;
                // compatibleIntervals retains the same ordered node/edge
                // indices in all sixteen quarters; only their powers differ.
                for(int c=0;c<2;++c){const auto& xx=c==0?bx.nodes:bx.edges;const auto& yy=c==0?by.edges:by.nodes;
                    for(const auto& y:yy)for(const auto& x:xx){
                        const int u=d.compatibleVelocity[c][c==0?y.index*(d.nx+1)+x.index:y.index*d.nx+x.index];
                        for(const auto& w:virtualRow(d,u))
                            if(!std::binary_search(stencil.excludedVelocityDofs_.begin(),stencil.excludedVelocityDofs_.end(),w.velocityDof))nativeVelocity.insert(w.velocityDof);
                    }}
                plan.nativeVelocity.assign(nativeVelocity.begin(),nativeVelocity.end());
                for(int c=0;c<2;++c){const auto& xx=c==0?bx.nodes:bx.edges;const auto& yy=c==0?by.edges:by.nodes;
                    for(const auto& y:yy)for(const auto& x:xx){
                        const int u=d.compatibleVelocity[c][c==0?y.index*(d.nx+1)+x.index:y.index*d.nx+x.index];
                        VirtualJetScatterPlan::Row row;row.virtualVelocity=u;
                        for(const auto& w:virtualRow(d,u)){
                            if(std::binary_search(stencil.excludedVelocityDofs_.begin(),stencil.excludedVelocityDofs_.end(),w.velocityDof))continue;
                            const auto slot=std::lower_bound(plan.nativeVelocity.begin(),plan.nativeVelocity.end(),w.velocityDof)-plan.nativeVelocity.begin();
                            row.weights.push_back({static_cast<std::size_t>(slot),w.weight});}
                        plan.rows.push_back(std::move(row));
                    }}
                scatter->plan=std::move(plan);d.virtualJetScatterPrepared.fetch_add(1,std::memory_order_relaxed);
            });
        }
        std::vector<HydroRegularJetWeight> out;
        std::vector<std::array<long double,6>> prolongedWeights(scatter?scatter->plan.nativeVelocity.size():0);
        std::size_t tensorRow=0;
        for(int c=0;c<2;++c){const auto& xx=c==0?bx.nodes:bx.edges;const auto& yy=c==0?by.edges:by.nodes;
            for(const auto& y:yy){auto vy=evaluate(y.powers,sy);
                if(c==1&&point.y==sourceCell.lower.y)vy[0]=y.index==logical[1]?1.:0.;
                if(c==1&&point.y==sourceCell.upper.y)vy[0]=y.index==logical[1]+1?1.:0.;
                for(const auto& x:xx){auto vx=evaluate(x.powers,sx);
                if(c==0&&point.x==sourceCell.lower.x)vx[0]=x.index==logical[0]?1.:0.;
                if(c==0&&point.x==sourceCell.upper.x)vx[0]=x.index==logical[0]+1?1.:0.;
                const int u=d.compatibleVelocity[c][c==0?y.index*(d.nx+1)+x.index:y.index*d.nx+x.index];
                if(!d.virtualProlongation&&std::binary_search(stencil.excludedVelocityDofs_.begin(),stencil.excludedVelocityDofs_.end(),u))continue;
                HydroRegularJetWeight row;row.velocityDof=u;row.jet[c]=vx[0]*vy[0];
                row.jet[2+2*c]=vx[1]*vy[0]/h.x;row.jet[3+2*c]=vx[0]*vy[1]/h.y;
                if(d.virtualProlongation){
                    const auto& prepared=scatter->plan.rows[tensorRow++];
                    require(prepared.virtualVelocity==u,"different cached virtual tensor support order");
                    for(const auto& w:prepared.weights){auto& target=prolongedWeights[w.slot];
                        for(int k=0;k<6;++k)target[k]+=static_cast<long double>(w.weight)*row.jet[k];}}
                else out.push_back(row);
            }}
        }
        if(d.virtualProlongation)for(std::size_t slot=0;slot<prolongedWeights.size();++slot){HydroRegularJetWeight row;row.velocityDof=scatter->plan.nativeVelocity[slot];bool nonzero=false;
            for(int k=0;k<6;++k){row.jet[k]=static_cast<double>(prolongedWeights[slot][k]);require(std::isfinite(row.jet[k]),"nonfinite prolonged symbolic weight");nonzero=nonzero||row.jet[k]!=0;}
            if(nonzero)out.push_back(row);}
        return out;
    }
    require(stencil.faceCorrections.size()==cell.faces.size(),"different local stencil face inventory");
    const auto h=cell.upper-cell.lower;const auto bx=basis((point.x-cell.lower.x)/h.x),by=basis((point.y-cell.lower.y)/h.y);
    const auto dx=basis((point.x-cell.lower.x)/h.x,true),dy=basis((point.y-cell.lower.y)/h.y,true);
    std::map<int,std::array<long double,6>> weights;
    const auto add=[&](const auto& list,int c,double value,double x,double y){for(const auto& w:list){auto& jet=weights[w.velocityDof];
        jet[c]+=static_cast<long double>(w.weight)*value;jet[2+2*c]+=static_cast<long double>(w.weight)*x;jet[3+2*c]+=static_cast<long double>(w.weight)*y;}};
    for(int c=0;c<2;++c)for(int y=0;y<5;++y)for(int x=0;x<5;++x)add(stencil.cellDofs[c][5*y+x],c,bx[x]*by[y],dx[x]*by[y]/h.x,bx[x]*dy[y]/h.y);
    for(std::size_t i=0;i<cell.faces.size();++i){const auto& face=faces()[cell.faces[i]];const int n=face.normalAxis,t=1-n;
        if(face.area==0)continue;const auto width=supportWidth(face,stencil.pressure),length=coord(face.upper-face.lower,t);
        const auto s=(coord(point,t)-coord(face.lower,t))/length,xi=(coord(point,n)-coord(face.lower,n))/width;if(s<0||s>1||xi<-1||xi>1)continue;
        const auto g=16*s*s*(1-s)*(1-s),b=(1-xi*xi)*(1-xi*xi),gt=32*s*(1-s)*(1-2*s)/length,bn=-4*xi*(1-xi*xi)/width;
        add(stencil.faceCorrections[i],n,b*g,n==0?bn*g:b*gt,n==0?b*gt:bn*g);}
    std::vector<HydroRegularJetWeight> out;for(const auto& w:weights){HydroRegularJetWeight row;row.velocityDof=w.first;bool nonzero=false;
        for(int k=0;k<6;++k){row.jet[k]=static_cast<double>(w.second[k]);require(std::isfinite(row.jet[k]),"nonfinite local jet weight");nonzero=nonzero||row.jet[k]!=0;}
        if(nonzero)out.push_back(row);}return out;
}
void HydroRegularAdvectorGeometry::validateCoefficients(const HydroRegularCoefficients& coefficients,int pressure) const
{
    fp();require(coefficients.geometry_==data_,"different regular coefficient geometry identity");
    const auto tensorCount=data_->virtualProlongation?data_->virtualCells.size():cells().size();
    const auto tensorFaces=data_->virtualProlongation?data_->virtualProlongation->rows[0].size()+data_->virtualProlongation->rows[1].size():faces().size();
    require(coefficients.cellDofs.size()==(data_->compatible?16:1)*tensorCount&&coefficients.faceCorrections.size()==faces().size(),"different regular coefficient inventory");
    require(coefficients.compatibleNativeMeans.size()==(data_->compatible?tensorFaces:0),"different compatible native coefficient inventory");
    require(pressure>=0&&static_cast<std::size_t>(pressure)<cells().size(),"regular pressure index out of range");
}
HydroAdvectorCellJet HydroRegularAdvectorGeometry::cellJet(const HydroRegularCoefficients& coefficients,int pressure,Vector2D point) const
{
    validateCoefficients(coefficients,pressure);const auto& cell=cells()[pressure];inside(cell,point);
    if(data_->compatible)return compatibleJet(*data_,coefficients,preparedCompatiblePatch(*data_,pressure,point),point);
    const auto h=cell.upper-cell.lower;const auto bx=basis((point.x-cell.lower.x)/h.x),by=basis((point.y-cell.lower.y)/h.y);
    const auto dx=basis((point.x-cell.lower.x)/h.x,true),dy=basis((point.y-cell.lower.y)/h.y,true);
    std::array<long double,2> value{},jx{},jy{};
    for(int component=0;component<2;++component)for(int y=0;y<5;++y)for(int x=0;x<5;++x){const auto datum=coefficients.cellDofs[pressure][component][5*y+x];
        require(std::isfinite(datum),"nonfinite regular Q4 datum");value[component]+=datum*bx[x]*by[y];
        jx[component]+=datum*dx[x]*by[y]/h.x;jy[component]+=datum*bx[x]*dy[y]/h.y;}
    for(int u:cell.faces){const auto& face=faces()[u];const int n=face.normalAxis,t=1-n;const auto delta=coefficients.faceCorrections[u];
        require(std::isfinite(delta),"nonfinite regular face correction");if(delta==0)continue;
        const auto width=supportWidth(face,pressure),length=coord(face.upper-face.lower,t);
        const auto s=(coord(point,t)-coord(face.lower,t))/length,xi=(coord(point,n)-coord(face.lower,n))/width;
        if(s<0||s>1||xi<-1||xi>1)continue;
        const auto g=16*s*s*(1-s)*(1-s),b=(1-xi*xi)*(1-xi*xi);
        const auto gt=32*s*(1-s)*(1-2*s)/length,bn=-4*xi*(1-xi*xi)/width;
        value[n]+=delta*b*g;(n==0?jx:jy)[n]+=delta*bn*g;(t==0?jx:jy)[n]+=delta*b*gt;
    }
    Vector2D lower,upper;
    for(int axis=0;axis<2;++axis){const auto& knots=cell.supportKnots[axis];const auto coordinate=coord(point,axis);
        require(knots.size()>=2&&knots.front()==coord(cell.lower,axis)&&knots.back()==coord(cell.upper,axis),"incomplete regular support partition");
        auto index=static_cast<std::size_t>(std::upper_bound(knots.begin(),knots.end(),coordinate)-knots.begin());if(index==knots.size())--index;
        (axis==0?lower.x:lower.y)=knots[index-1];(axis==0?upper.x:upper.y)=knots[index];
    }
    for(int c=0;c<2;++c)require(std::isfinite(static_cast<double>(value[c]))&&std::isfinite(static_cast<double>(jx[c]))&&std::isfinite(static_cast<double>(jy[c])),"nonfinite regular jet");
    return {{{static_cast<double>(value[0]),static_cast<double>(value[1])},
        {static_cast<double>(jx[0]),static_cast<double>(jy[0]),static_cast<double>(jx[1]),static_cast<double>(jy[1])}},lower,upper};
}
HydroAdvectorCellJetPair HydroRegularAdvectorGeometry::pieceJetPair(
    const HydroRegularCoefficients& previous,const HydroRegularCoefficients& arrival,
    const HydroAdvectorPiece& piece,Vector2D point) const
{
    validateCoefficients(previous,piece.pressure);validateCoefficients(arrival,piece.pressure);
    const auto& cell=cells()[piece.pressure];inside(cell,point);
    require(point.x>=piece.lower.x&&point.x<=piece.upper.x&&point.y>=piece.lower.y&&point.y<=piece.upper.y,
        "point outside selected piece");
    if(data_->compatible){auto patch=preparedCompatiblePatch(*data_,piece.pressure,(piece.lower+piece.upper)*.5);
        auto a=compatibleJet(*data_,previous,patch,point),b=compatibleJet(*data_,arrival,patch,point);
        a.lower=b.lower=piece.lower;a.upper=b.upper=piece.upper;return {a,b};}
    const auto h=cell.upper-cell.lower;const auto bx=basis((point.x-cell.lower.x)/h.x),by=basis((point.y-cell.lower.y)/h.y);
    const auto dx=basis((point.x-cell.lower.x)/h.x,true),dy=basis((point.y-cell.lower.y)/h.y,true);
    const std::array<const HydroRegularCoefficients*,2> fields{{&previous,&arrival}};
    std::array<std::array<long double,2>,2> value{},jx{},jy{};
    for(int endpoint=0;endpoint<2;++endpoint)
        for(int component=0;component<2;++component)for(int y=0;y<5;++y)for(int x=0;x<5;++x){
            const auto datum=fields[endpoint]->cellDofs[piece.pressure][component][5*y+x];
            require(std::isfinite(datum),"nonfinite regular Q4 datum");
            // Do not pre-multiply basis products: each original double
            // multiply/divide and long-double accumulation stays unchanged.
            value[endpoint][component]+=datum*bx[x]*by[y];
            jx[endpoint][component]+=datum*dx[x]*by[y]/h.x;jy[endpoint][component]+=datum*bx[x]*dy[y]/h.y;
        }
    for(int u:cell.faces){const auto& face=faces()[u];const int n=face.normalAxis,t=1-n;
        const std::array<double,2> delta{{previous.faceCorrections[u],arrival.faceCorrections[u]}};
        for(double correction:delta)require(std::isfinite(correction),"nonfinite regular face correction");
        if(delta[0]==0&&delta[1]==0)continue;
        const auto width=supportWidth(face,piece.pressure),length=coord(face.upper-face.lower,t);
        const auto s=(coord(point,t)-coord(face.lower,t))/length,xi=(coord(point,n)-coord(face.lower,n))/width;
        if(s<0||s>1||xi<-1||xi>1)continue;
        const auto g=16*s*s*(1-s)*(1-s),b=(1-xi*xi)*(1-xi*xi);
        const auto gt=32*s*(1-s)*(1-2*s)/length,bn=-4*xi*(1-xi*xi)/width;
        for(int endpoint=0;endpoint<2;++endpoint){if(delta[endpoint]==0)continue;
            value[endpoint][n]+=delta[endpoint]*b*g;
            (n==0?jx:jy)[endpoint][n]+=delta[endpoint]*bn*g;
            (t==0?jx:jy)[endpoint][n]+=delta[endpoint]*b*gt;
        }
    }
    const auto pack=[&](int endpoint){
        for(int c=0;c<2;++c)require(std::isfinite(static_cast<double>(value[endpoint][c]))
            &&std::isfinite(static_cast<double>(jx[endpoint][c]))&&std::isfinite(static_cast<double>(jy[endpoint][c])),
            "nonfinite regular jet");
        return HydroAdvectorCellJet{{{static_cast<double>(value[endpoint][0]),static_cast<double>(value[endpoint][1])},
            {static_cast<double>(jx[endpoint][0]),static_cast<double>(jy[endpoint][0]),
             static_cast<double>(jx[endpoint][1]),static_cast<double>(jy[endpoint][1])}},piece.lower,piece.upper};
    };
    // The reconstruction already selected this piece; no cell partition
    // lookup is needed to produce bounds that would immediately be replaced.
    const HydroAdvectorCellJetPair result{pack(0),pack(1)};
    fp();return result;
}
namespace
{
// Degree sixteen is the actual Q4xQ4 composition with a quadratic curve.
// Spatial derivative times a quadratic P has the same degree; theta moments
// are integrated directly below, so no degree-seventeen storage is needed.
struct CurvePower
{
    std::array<long double,17> a{};
    int degree=0;
};
void trim(CurvePower& p)
{while(p.degree>0&&p.a[p.degree]==0)--p.degree;}
CurvePower constantPower(long double a)
{CurvePower p;p.a[0]=a;return p;}
CurvePower addPower(CurvePower p,const CurvePower& q)
{
    p.degree=std::max(p.degree,q.degree);
    for(int i=0;i<=q.degree;++i)p.a[i]+=q.a[i];
    trim(p);return p;
}
CurvePower timesQuadratic(const CurvePower& p,const std::array<long double,3>& q)
{
    int degree=2;while(degree>0&&q[degree]==0)--degree;
    require(p.degree+degree<=16,"curve polynomial degree exceeds actual Q4");
    CurvePower out;out.degree=p.degree+degree;
    for(int i=0;i<=p.degree;++i)for(int j=0;j<=degree;++j)out.a[i+j]+=p.a[i]*q[j];
    trim(out);return out;
}
std::array<long double,3> shiftedQuadratic(long double a,long double b,long double c,long double lo,long double width)
{return {a+lo*(b+lo*c),width*(b+2*lo*c),width*width*c};}
std::array<long double,2> curveMoments(const CurvePower& p,long double lo,long double width)
{
    std::array<long double,2> out{};
    for(int i=0;i<=p.degree;++i){
        out[0]+=p.a[i]/(i+1);
        out[1]+=p.a[i]*(lo/(i+1)+width/(i+2));
    }
    for(auto& x:out)x*=width;
    return out;
}
}
std::vector<HydroRegularPiecePolynomial> HydroRegularAdvectorGeometry::piecePolynomials(
    const HydroRegularCoefficients& coefficients,const std::vector<HydroAdvectorPiece>& pieces) const
{
    fp();std::vector<HydroRegularPiecePolynomial> result;result.reserve(pieces.size());
    for(const auto& piece:pieces){
        validateCoefficients(coefficients,piece.pressure);
        const auto& cell=cells()[piece.pressure];inside(cell,piece.lower);inside(cell,piece.upper);
        const auto h=cell.upper-cell.lower,w=piece.upper-piece.lower;
        require(w.x>0&&w.y>0,"invalid polynomial piece bounds");
        for(int axis=0;axis<2;++axis)for(double knot:cell.supportKnots[axis])
            require(!(coord(piece.lower,axis)<knot&&knot<coord(piece.upper,axis)),"piece crosses regular support knot");
        HydroRegularPiecePolynomial out;out.pressure=piece.pressure;out.lower=piece.lower;out.upper=piece.upper;out.geometry_=data_;
        if(data_->compatible){
            const auto patch=preparedCompatiblePatch(*data_,piece.pressure,(piece.lower+piece.upper)*.5);const auto span=patch.upper-patch.lower;
            require(piece.lower.x>=patch.lower.x&&piece.upper.x<=patch.upper.x&&piece.lower.y>=patch.lower.y&&piece.upper.y<=patch.upper.y,"compatible polynomial crosses quarter knot");
            if(piece.lower.x==patch.lower.x&&piece.lower.y==patch.lower.y&&piece.upper.x==patch.upper.x&&piece.upper.y==patch.upper.y){
                for(int c=0;c<2;++c)for(int k=0;k<25;++k)out.powers[c][k]=coefficients.cellDofs[patch.index][c][k];
                result.push_back(std::move(out));continue;
            }
            std::array<std::array<long double,5>,5> bx{},by{};
            for(int k=0;k<5;++k){double power[5]={};power[k]=1;
                bx[k]=compose(power,(static_cast<long double>(piece.lower.x)-patch.lower.x)/span.x,static_cast<long double>(w.x)/span.x);
                by[k]=compose(power,(static_cast<long double>(piece.lower.y)-patch.lower.y)/span.y,static_cast<long double>(w.y)/span.y);}
            for(int c=0;c<2;++c)for(int j=0;j<5;++j)for(int i=0;i<5;++i){const double value=coefficients.cellDofs[patch.index][c][5*j+i];
                for(int y=0;y<5;++y)for(int x=0;x<5;++x)out.powers[c][5*y+x]+=value*bx[i][x]*by[j][y];}
            result.push_back(std::move(out));continue;
        }
        std::array<std::array<long double,5>,5> bx{},by{};
        for(int i=0;i<5;++i){
            bx[i]=compose(basisPower[i],(static_cast<long double>(piece.lower.x)-cell.lower.x)/h.x,static_cast<long double>(w.x)/h.x);
            by[i]=compose(basisPower[i],(static_cast<long double>(piece.lower.y)-cell.lower.y)/h.y,static_cast<long double>(w.y)/h.y);
        }
        for(int component=0;component<2;++component){
            std::array<std::array<long double,5>,5> rows{};
            for(int y=0;y<5;++y)for(int x=0;x<5;++x){
                const auto datum=coefficients.cellDofs[piece.pressure][component][5*y+x];
                require(std::isfinite(datum),"nonfinite regular Q4 datum");
                for(int k=0;k<5;++k)rows[y][k]+=datum*bx[x][k];
            }
            for(int y=0;y<5;++y)for(int x=0;x<5;++x)for(int k=0;k<5;++k)
                out.powers[component][5*y+x]+=rows[k][x]*by[k][y];
        }
        const auto middle=piece.lower+w*.5;
        constexpr double bubblePower[5]={1,0,-2,0,1};
        for(int u:cell.faces){
            const auto& face=faces()[u];const int n=face.normalAxis,t=1-n;
            const auto delta=coefficients.faceCorrections[u];
            require(std::isfinite(delta),"nonfinite regular face correction");if(delta==0)continue;
            const auto width=supportWidth(face,piece.pressure),length=coord(face.upper-face.lower,t);
            require(length>0,"invalid regular face tangent support");
            const auto xi=(coord(middle,n)-coord(face.lower,n))/width;
            if(xi<-1||xi>1||coord(middle,t)<coord(face.lower,t)||coord(middle,t)>coord(face.upper,t))continue;
            const auto bn=compose(bubblePower,(static_cast<long double>(coord(piece.lower,n))-coord(face.lower,n))/width,static_cast<long double>(coord(w,n))/width);
            const auto bt=compose(basisPower[4],(static_cast<long double>(coord(piece.lower,t))-coord(face.lower,t))/length,static_cast<long double>(coord(w,t))/length);
            for(int y=0;y<5;++y)for(int x=0;x<5;++x)
                out.powers[n][5*y+x]+=delta*(n==0?bn[x]*bt[y]:bt[x]*bn[y]);
        }
        for(const auto& component:out.powers)for(auto value:component)
            require(std::isfinite(value)&&std::isfinite(static_cast<double>(value)),"nonfinite piece polynomial");
        result.push_back(std::move(out));
    }
    fp();return result;
}
HydroAdvectorCurveIntegral HydroRegularAdvectorGeometry::curveIntegral(const HydroRegularPiecePolynomial& piece,
    Vector2D r0,Vector2D r1,Vector2D r2,LevelSetBoundaryTrace::Jacobian P0,
    LevelSetBoundaryTrace::Jacobian P1,LevelSetBoundaryTrace::Jacobian P2,double lo,double hi) const
{
    fp();require(piece.geometry_==data_&&piece.pressure>=0&&static_cast<std::size_t>(piece.pressure)<cells().size(),"different curve polynomial geometry");
    require(std::isfinite(lo)&&std::isfinite(hi)&&lo>=0&&hi>lo&&hi<=1,"invalid curve span");
    const auto finitePoint=[](Vector2D p){return std::isfinite(p.x)&&std::isfinite(p.y);};
    const auto finiteMatrix=[](const auto& p){return std::isfinite(p.xx)&&std::isfinite(p.xy)&&std::isfinite(p.yx)&&std::isfinite(p.yy);};
    require(finitePoint(r0)&&finitePoint(r1)&&finitePoint(r2)&&finiteMatrix(P0)&&finiteMatrix(P1)&&finiteMatrix(P2),"nonfinite curve data");
    const long double left=lo,width=static_cast<long double>(hi)-lo;
    const auto h=piece.upper-piece.lower;
    std::array<std::array<long double,3>,2> chart{};
    for(int axis=0;axis<2;++axis){
        const long double a=coord(r0,axis),b=coord(r1,axis),c=coord(r2,axis);
        const auto check=[&](long double theta){
            const auto x=a+theta*(b+theta*c);const auto lower=coord(piece.lower,axis),upper=coord(piece.upper,axis);
            require(std::isfinite(x)&&x>=lower&&x<=upper,"curve outside actual piece");
        };
        check(lo);check(hi);if(c!=0){const auto stationary=-b/(2*c);if(stationary>lo&&stationary<hi)check(stationary);}
        chart[axis]=shiftedQuadratic((a-coord(piece.lower,axis))/coord(h,axis),b/coord(h,axis),c/coord(h,axis),left,width);
    }
    const std::array<double,4> p0{P0.xx,P0.xy,P0.yx,P0.yy},p1{P1.xx,P1.xy,P1.yx,P1.yy},p2{P2.xx,P2.xy,P2.yx,P2.yy};
    std::array<std::array<long double,3>,4> arrival{};
    for(int k=0;k<4;++k)arrival[k]=shiftedQuadratic(p0[k],p1[k],p2[k],left,width);
    HydroAdvectorCurveIntegral result;
    for(int component=0;component<2;++component){
        CurvePower value,dx,dy;
        for(int y=4;y>=0;--y){
            CurvePower row=constantPower(piece.powers[component][5*y+4]),rowDx;
            for(int x=3;x>=0;--x){
                rowDx=addPower(timesQuadratic(rowDx,chart[0]),row);
                row=addPower(timesQuadratic(row,chart[0]),constantPower(piece.powers[component][5*y+x]));
            }
            dy=addPower(timesQuadratic(dy,chart[1]),value);
            dx=addPower(timesQuadratic(dx,chart[1]),rowDx);
            value=addPower(timesQuadratic(value,chart[1]),row);
        }
        const auto moments=curveMoments(value,left,width);
        result.V0[component]=moments[0];result.V1[component]=moments[1];
        for(int k=0;k<=dx.degree;++k)dx.a[k]/=h.x;
        for(int k=0;k<=dy.degree;++k)dy.a[k]/=h.y;
        for(int column=0;column<2;++column){
            const auto derivative=addPower(timesQuadratic(dx,arrival[column]),timesQuadratic(dy,arrival[2+column]));
            const auto k=curveMoments(derivative,left,width);const int index=2*component+column;
            result.K0[index]=k[0];result.K1[index]=k[1];
        }
    }
    for(const auto& values:{result.V0,result.V1})for(auto v:values)
        require(std::isfinite(v)&&std::isfinite(static_cast<double>(v)),"nonfinite curve integral");
    for(const auto& values:{result.K0,result.K1})for(auto v:values)
        require(std::isfinite(v)&&std::isfinite(static_cast<double>(v)),"nonfinite curve integral");
    fp();return result;
}
HydroAdvectorCurveIntegralPair HydroRegularAdvectorGeometry::curveIntegralPair(
    const HydroRegularPiecePolynomial& previous,const HydroRegularPiecePolynomial& arrivalPiece,
    Vector2D r0,Vector2D r1,Vector2D r2,LevelSetBoundaryTrace::Jacobian P0,
    LevelSetBoundaryTrace::Jacobian P1,LevelSetBoundaryTrace::Jacobian P2,double lo,double hi) const
{
    fp();
    require(previous.geometry_==data_&&arrivalPiece.geometry_==data_
        &&previous.pressure>=0&&static_cast<std::size_t>(previous.pressure)<cells().size()
        &&arrivalPiece.pressure>=0&&static_cast<std::size_t>(arrivalPiece.pressure)<cells().size(),
        "different curve polynomial geometry");
    require(previous.pressure==arrivalPiece.pressure
        &&previous.lower.x==arrivalPiece.lower.x&&previous.lower.y==arrivalPiece.lower.y
        &&previous.upper.x==arrivalPiece.upper.x&&previous.upper.y==arrivalPiece.upper.y,
        "different curve polynomial piece");
    require(std::isfinite(lo)&&std::isfinite(hi)&&lo>=0&&hi>lo&&hi<=1,"invalid curve span");
    const auto finitePoint=[](Vector2D p){return std::isfinite(p.x)&&std::isfinite(p.y);};
    const auto finiteMatrix=[](const auto& p){return std::isfinite(p.xx)&&std::isfinite(p.xy)
        &&std::isfinite(p.yx)&&std::isfinite(p.yy);};
    require(finitePoint(r0)&&finitePoint(r1)&&finitePoint(r2)
        &&finiteMatrix(P0)&&finiteMatrix(P1)&&finiteMatrix(P2),"nonfinite curve data");
    const long double left=lo,width=static_cast<long double>(hi)-lo;
    const auto h=previous.upper-previous.lower;
    require(finitePoint(h)&&h.x>0&&h.y>0,"invalid polynomial piece bounds");
    std::array<std::array<long double,3>,2> chart{};
    for(int axis=0;axis<2;++axis){
        const long double a=coord(r0,axis),b=coord(r1,axis),c=coord(r2,axis);
        const auto check=[&](long double theta){
            const auto x=a+theta*(b+theta*c);
            require(std::isfinite(x)&&x>=coord(previous.lower,axis)&&x<=coord(previous.upper,axis),
                "curve outside actual piece");
        };
        check(lo);check(hi);
        if(c!=0){const auto stationary=-b/(2*c);if(stationary>lo&&stationary<hi)check(stationary);}
        chart[axis]=shiftedQuadratic((a-coord(previous.lower,axis))/coord(h,axis),
            b/coord(h,axis),c/coord(h,axis),left,width);
        for(auto value:chart[axis])require(std::isfinite(value),"nonfinite curve data");
    }
    const std::array<double,4> p0{P0.xx,P0.xy,P0.yx,P0.yy},
        p1{P1.xx,P1.xy,P1.yx,P1.yy},p2{P2.xx,P2.xy,P2.yx,P2.yy};
    std::array<std::array<long double,3>,4> arrival{};
    for(int k=0;k<4;++k){
        arrival[k]=shiftedQuadratic(p0[k],p1[k],p2[k],left,width);
        for(auto value:arrival[k])require(std::isfinite(value),"nonfinite curve data");
    }
    // Only exact zero coefficients remove a monomial. Include the two
    // first-derivative neighbours needed by the same physical Jacobian.
    std::array<bool,25> active{},needed{};
    int maxX=0,maxY=0;
    for(int y=0;y<5;++y)for(int x=0;x<5;++x){
        const int index=5*y+x;
        for(int component=0;component<2;++component)
            active[index]=active[index]||previous.powers[component][index]!=0
                ||arrivalPiece.powers[component][index]!=0;
        if(active[index]){
            needed[index]=true;if(x>0)needed[index-1]=true;if(y>0)needed[index-5]=true;
            maxX=std::max(maxX,x);maxY=std::max(maxY,y);
        }
    }
    std::array<std::array<CurvePower,5>,2> axisPowers{};
    for(int axis=0;axis<2;++axis){
        axisPowers[axis][0]=constantPower(1);
        const int maximum=axis==0?maxX:maxY;
        for(int degree=1;degree<=maximum;++degree)
            axisPowers[axis][degree]=timesQuadratic(axisPowers[axis][degree-1],chart[axis]);
    }
    // M[k][i,j] = width * integral_0^1 s^k X(s)^i Y(s)^j ds.
    // A basis has degree <=16; k<=3 needs only denominators 1..20.
    std::array<long double,20> measures{};
    for(int n=0;n<20;++n)measures[n]=width/(n+1);
    std::array<std::array<long double,4>,25> moments{};
    for(int y=0;y<=maxY;++y)for(int x=0;x<=maxX;++x){
        const int index=5*y+x;if(!needed[index])continue;
        const auto& bx=axisPowers[0][x];const auto& by=axisPowers[1][y];
        require(bx.degree+by.degree<=16,"curve polynomial degree exceeds actual Q4");
        CurvePower basis;basis.degree=bx.degree+by.degree;
        for(int i=0;i<=bx.degree;++i)for(int j=0;j<=by.degree;++j)
            basis.a[i+j]+=bx.a[i]*by.a[j];
        for(int i=0;i<=basis.degree;++i)for(int k=0;k<4;++k)
            moments[index][k]+=basis.a[i]*measures[i+k];
    }
    HydroAdvectorCurveIntegralPair result;
    for(int y=4;y>=0;--y)for(int x=4;x>=0;--x){
        const int index=5*y+x;if(!active[index])continue;
        const auto& m=moments[index];
        const std::array<long double,2> value{m[0],left*m[0]+width*m[1]};
        std::array<long double,2> derivative0{},derivative1{};
        for(int column=0;column<2;++column)for(int k=0;k<3;++k){
            if(x>0){
                const auto& dx=moments[index-1];const long double factor=static_cast<long double>(x)/h.x;
                derivative0[column]+=factor*arrival[column][k]*dx[k];
                derivative1[column]+=factor*arrival[column][k]*(left*dx[k]+width*dx[k+1]);
            }
            if(y>0){
                const auto& dy=moments[index-5];const long double factor=static_cast<long double>(y)/h.y;
                derivative0[column]+=factor*arrival[2+column][k]*dy[k];
                derivative1[column]+=factor*arrival[2+column][k]*(left*dy[k]+width*dy[k+1]);
            }
        }
        for(int endpoint=0;endpoint<2;++endpoint){
            const auto& piece=endpoint==0?previous:arrivalPiece;
            auto& integral=endpoint==0?result.previous:result.arrival;
            for(int component=0;component<2;++component){
                const auto coefficient=piece.powers[component][index];
                integral.V0[component]+=coefficient*value[0];
                integral.V1[component]+=coefficient*value[1];
                for(int column=0;column<2;++column){
                    const int entry=2*component+column;
                    integral.K0[entry]+=coefficient*derivative0[column];
                    integral.K1[entry]+=coefficient*derivative1[column];
                }
            }
        }
    }
    for(const auto* integral:{&result.previous,&result.arrival}){
        for(const auto& values:{integral->V0,integral->V1})for(auto value:values)
            require(std::isfinite(value)&&std::isfinite(static_cast<double>(value)),"nonfinite curve integral");
        for(const auto& values:{integral->K0,integral->K1})for(auto value:values)
            require(std::isfinite(value)&&std::isfinite(static_cast<double>(value)),"nonfinite curve integral");
    }
    fp();return result;
}
HydroNormalPolynomial HydroRegularAdvectorGeometry::normalPolynomial(const HydroRegularCoefficients& coefficients,int pressure,
    int normal,Vector2D begin,Vector2D end) const
{
    validateCoefficients(coefficients,pressure);require(normal==0||normal==1,"invalid regular normal axis");
    const auto& cell=cells()[pressure];inside(cell,begin);inside(cell,end);const int tangent=1-normal;
    require(coord(begin,normal)==coord(end,normal)&&coord(end,tangent)>coord(begin,tangent),"invalid canonical regular normal interval");
    for(double knot:cell.supportKnots[tangent])require(!(coord(begin,tangent)<knot&&knot<coord(end,tangent)),"normal interval crosses regular support knot");
    HydroNormalPolynomial out;out.length=static_cast<long double>(coord(end,tangent))-coord(begin,tangent);
    out.radialMeasure=momentumGeometry().layout().domain().axisymmetric();out.radialTangent=out.radialMeasure&&normal==1;out.radiusBegin=begin.x;out.radiusEnd=end.x;
    if(data_->compatible){
        const auto patch=preparedCompatiblePatch(*data_,pressure,(begin+end)*.5);const auto h=patch.upper-patch.lower;
        const double z=(coord(begin,normal)-coord(patch.lower,normal))/coord(h,normal);
        std::array<long double,5> normalPowers{1,z,z*z,z*z*z,z*z*z*z};double tangentPowers[5]={};
        if(const auto trace=compatibleTrace(*data_,coefficients,patch.index/16,normal,(begin+end)*.5)){
            for(int j=0;j<5;++j)tangentPowers[j]=(*trace)[j];
        }else for(int j=0;j<5;++j){long double value=0;for(int i=0;i<5;++i)
            value+=coefficients.cellDofs[patch.index][normal][normal==0?5*j+i:5*i+j]*normalPowers[i];tangentPowers[j]=static_cast<double>(value);}
        for(int reverse=0;reverse<2;++reverse){const auto origin=reverse?end:begin;const long double sign=reverse?-1:1;
            auto& powers=reverse?out.endCoefficients:out.beginCoefficients;
            powers=compose(tangentPowers,(static_cast<long double>(coord(origin,tangent))-coord(patch.lower,tangent))/coord(h,tangent),sign/coord(h,tangent));
        }
        return out;
    }
    const auto h=cell.upper-cell.lower;const auto normalBasis=basis((coord(begin,normal)-coord(cell.lower,normal))/coord(h,normal));
    for(int reverse=0;reverse<2;++reverse){const auto origin=reverse?end:begin;auto& result=reverse?out.endCoefficients:out.beginCoefficients;
        const long double sign=reverse?-1:1,z=(static_cast<long double>(coord(origin,tangent))-coord(cell.lower,tangent))/coord(h,tangent);
        for(int j=0;j<5;++j){const auto restriction=compose(basisPower[j],z,sign/coord(h,tangent));long double weight=0;
            for(int i=0;i<5;++i){const auto datum=coefficients.cellDofs[pressure][normal][normal==0?5*j+i:5*i+j];
                require(std::isfinite(datum),"nonfinite regular Q4 datum");weight+=datum*normalBasis[i];}
            for(int k=0;k<5;++k)result[k]+=weight*restriction[k];}
        for(int u:cell.faces){const auto& face=faces()[u];if(face.normalAxis!=normal)continue;
            const auto delta=coefficients.faceCorrections[u];require(std::isfinite(delta),"nonfinite regular face correction");if(delta==0)continue;
            const auto width=supportWidth(face,pressure),length=coord(face.upper-face.lower,tangent);
            const auto xi=(coord(origin,normal)-coord(face.lower,normal))/width;
            const auto mid=(static_cast<long double>(coord(begin,tangent))+coord(end,tangent))/2;
            if(xi<-1||xi>1||mid<coord(face.lower,tangent)||mid>coord(face.upper,tangent))continue;
            const auto restriction=compose(basisPower[4],(static_cast<long double>(coord(origin,tangent))-coord(face.lower,tangent))/length,sign/length);
            const long double weight=delta*(1-xi*xi)*(1-xi*xi);for(int k=0;k<5;++k)result[k]+=weight*restriction[k];
        }
        for(auto value:result)require(std::isfinite(value)&&std::isfinite(static_cast<double>(value)),"nonfinite regular normal polynomial");
    }
    return out;
}
}

namespace LAMP::Hydrodynamics
{
std::shared_ptr<Detail::HydroRegularOperatorTemplateBank> HydroMomentumGeometry::makeRegularTemplates()
{return std::make_shared<Detail::HydroRegularOperatorTemplateBank>();}
HydroRegularTemplateStatistics HydroMomentumGeometry::regularTemplateStatistics() const
{return regularTemplates()->statistics();}
}
