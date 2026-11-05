#include "Modules/Electric/ElectricVofRefinement.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{
using namespace LAMP::Electric;
using Publication=LAMP::Interface::GeometricVofPublication;
using Clock=LAMP::TimeLayer;
constexpr long double tolerance=4e-14L;
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class F>void rejects(F&& action,const char* why)
{bool caught=false;try{action();}catch(const std::invalid_argument&){caught=true;}check(caught,why);}
void close(long double actual,long double expected,const char* why)
{check(expected==0?actual==0:std::abs(actual-expected)<=tolerance*std::abs(expected),why);}
std::unique_ptr<QuadtreeMesh> mesh(bool graded,double width=1)
{
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={width,1};input.Nx=4;input.Ny=4;
    input.maxInitialIterations=0;input.maxCellNum=4096;input.maxLevel_local={1};
    auto result=std::make_unique<QuadtreeMesh>(input);
    if(graded){DictionaryMesh::RefineBox box;box.min={width/2,0};box.max=input.maxCoord;box.level=1;
        result->applyRefineBoxes({box});}
    return result;
}
NativeLeafDomain domain(const QuadtreeMesh& grid,bool rz,bool selected=false)
{
    FaceSegmentGraph graph;graph.build(grid,[selected](const QuadtreeMesh&,const Cell& cell){
        const auto p=cell.centerCoord();return !selected||p.x<.5||p.y>.5;});
    return NativeLeafDomain::build(grid,graph,rz);
}
Publication geometry(const QuadtreeMesh& grid,const NativeLeafDomain& d,
    Clock clock={3,.1,.2,Clock::Position::End})
{
    std::vector<std::array<long double,2>> stocks(grid.cellIndexNum());
    for(int native:d.cellIndexByRow())stocks[native]={d.cellVolumes()[d.cellRowByIndex()[native]],0};
    return Publication::prepareLocalPlicInventory(grid,d,clock,stocks,{});
}
// Independent moment quotient: integral(r^2 dr)/integral(r dr).
long double radialCentroid(const QuadtreeMesh& grid,int native,bool rz)
{
    const auto& cell=grid.cells(grid.cellIndexToLabel[native]);
    const long double center=cell.centerCoord().x,h=cell.size().x;
    if(!rz)return center;
    const long double a=center-h/2,b=center+h/2;
    return ((b*b*b-a*a*a)/3)/((b*b-a*a)/2);
}
template<class F>std::vector<double> integrated(const QuadtreeMesh& grid,const NativeLeafDomain& d,F&& density)
{
    std::vector<double> q(grid.cellIndexNum());
    for(int native:d.cellIndexByRow())q[native]=static_cast<double>(
        density(native)*d.cellVolumes()[d.cellRowByIndex()[native]]);
    return q;
}
void unchanged(const QuadtreeMesh& grid,const Publication& p,const ElectricVofCharge& charge,
    const std::vector<double>& before,const std::vector<std::array<long double,2>>& stocks)
{
    charge.validateFor(grid,p,p.timeLayer());check(charge.integratedCharge()==before,"indicator changed integrated signed Q");
    for(std::size_t ci=0;ci<stocks.size();++ci)
        check(p.plicCells()[ci].inventory==stocks[ci],"indicator changed independent phase stocks");
}
std::vector<std::array<long double,2>> snapshot(const Publication& p)
{std::vector<std::array<long double,2>> out;for(const auto& c:p.plicCells())out.push_back(c.inventory);return out;}
void constantsAndLinear(std::size_t& scenes,std::size_t& subfaces)
{
    for(bool rz:{false,true})for(bool graded:{false,true})for(bool selected:{false,true}){
        auto grid=mesh(graded);const auto d=domain(*grid,rz,selected);const auto p=geometry(*grid,d);
        const auto savedStocks=snapshot(p);
        for(long double density:{0.L,1.L,-8.L}){
            const auto q=integrated(*grid,d,[&](int){return density;});
            const auto charge=ElectricVofCharge::fromIntegrated(*grid,p,q);
            const auto values=ElectricVofRefinement::sample(*grid,p,charge);values.validateFor(*grid,p,charge);
            check(values.densityJump().size()==static_cast<std::size_t>(grid->cellIndexNum()),"indicator lost native indexing");
            close(values.maximumDensityJump(),0,"constant signed physical density has a false gradient");
            for(auto value:values.densityJump())close(value,0,"constant/inactive native indicator is not zero");
            unchanged(*grid,p,charge,q,savedStocks);++scenes;
        }
        // An affine radial/x density straddles zero. Exact volume centroids
        // make eta=h_x on every cell with a radial/x neighbour, including the
        // first RZ cell and a genuine split coarse/fine interface.
        const auto q=integrated(*grid,d,[&](int ci){return radialCentroid(*grid,ci,rz)-.5L;});
        const auto charge=ElectricVofCharge::fromIntegrated(*grid,p,q);
        const auto values=ElectricVofRefinement::sample(*grid,p,charge);
        std::vector<bool> radialNeighbour(grid->cellIndexNum());
        std::size_t internal=0,split=0;
        for(const auto face:d.faces(*grid))if(face.negativeCell>=0&&face.positiveCell>=0){
            ++internal;if(!face.horizontal){radialNeighbour[face.negativeCell]=radialNeighbour[face.positiveCell]=true;}
            const auto& a=grid->cells(grid->cellIndexToLabel[face.negativeCell]);
            const auto& b=grid->cells(grid->cellIndexToLabel[face.positiveCell]);
            if(a.level()!=b.level())++split;
        }
        check(values.interiorFaces()==internal&&internal>0,"actual native subfaces were skipped or duplicated");
        check(!graded||split>0,"graded fixture lacks actual coarse/fine subfaces");
        long double maximum=0;
        for(int ci:d.cellIndexByRow())if(radialNeighbour[ci]){
            const long double expected=grid->cells(grid->cellIndexToLabel[ci]).size().x;
            close(values.densityJump()[ci],expected,"affine density lost true normal scale or RZ volume centroid");
            maximum=std::max(maximum,expected);
        }
        close(values.maximumDensityJump(),maximum,"affine native maximum differs from physical scale");
        for(int ci=0;ci<grid->cellIndexNum();++ci)if(d.cellRowByIndex()[ci]<0)
            close(values.densityJump()[ci],0,"excluded cell received an implicit exterior charge gradient");
        unchanged(*grid,p,charge,q,savedStocks);subfaces+=internal;++scenes;
    }
}
void signedLayer(std::size_t& scenes)
{
    for(bool rz:{false,true})for(bool graded:{false,true}){
        auto grid=mesh(graded);const auto d=domain(*grid,rz);const auto p=geometry(*grid,d);
        const auto q=integrated(*grid,d,[&](int ci){return grid->cells(grid->cellIndexToLabel[ci]).centerCoord().x<.5?-1.L:1.L;});
        const auto charge=ElectricVofCharge::fromIntegrated(*grid,p,q);
        const auto values=ElectricVofRefinement::sample(*grid,p,charge);
        std::vector<long double> expected(grid->cellIndexNum());
        bool actualSplit=false;
        for(const auto face:d.faces(*grid))if(face.negativeCell>=0&&face.positiveCell>=0&&!face.horizontal){
            const int a=face.negativeCell,b=face.positiveCell;
            const auto& ca=grid->cells(grid->cellIndexToLabel[a]);const auto& cb=grid->cells(grid->cellIndexToLabel[b]);
            if((ca.centerCoord().x<.5)==(cb.centerCoord().x<.5))continue;
            const long double distance=radialCentroid(*grid,b,rz)-radialCentroid(*grid,a,rz);
            expected[a]=2*ca.size().x/distance;expected[b]=2*cb.size().x/distance;
            actualSplit|=ca.level()!=cb.level();
        }
        check(!graded||actualSplit,"signed layer does not cross actual coarse/fine faces");
        for(std::size_t ci=0;ci<expected.size();++ci)
            close(values.densityJump()[ci],expected[ci],"signed layer was cancelled or refined with the wrong native scale");
        check(values.maximumDensityJump()>1,"opposite-sign charge layer vanished under net-charge cancellation");++scenes;
    }
}
void axialRectangles(std::size_t& scenes)
{
    // Rectangular native cells distinguish normal h_y from h_x. The signed
    // affine axial/y density has derivative 2 and a true zero crossing.
    for(bool rz:{false,true}){
        auto grid=mesh(false,2);const auto d=domain(*grid,rz);const auto p=geometry(*grid,d);
        const auto q=integrated(*grid,d,[&](int ci){
            return 2*static_cast<long double>(grid->cells(grid->cellIndexToLabel[ci]).centerCoord().y)-1;
        });
        const auto charge=ElectricVofCharge::fromIntegrated(*grid,p,q);
        const auto values=ElectricVofRefinement::sample(*grid,p,charge);
        for(int ci:d.cellIndexByRow())close(values.densityJump()[ci],
            2*grid->cells(grid->cellIndexToLabel[ci]).size().y,
            "axial/y jump used tangential spacing instead of the true normal scale");
        ++scenes;
    }
}
void guards(std::size_t& count)
{
    auto grid=mesh(false);const auto d=domain(*grid,false);const auto p=geometry(*grid,d);
    const std::vector<double> q(grid->cellIndexNum());const auto charge=ElectricVofCharge::fromIntegrated(*grid,p,q);
    const auto values=ElectricVofRefinement::sample(*grid,p,charge);
    const auto other=geometry(*grid,d);const auto otherCharge=ElectricVofCharge::fromIntegrated(*grid,p,q);
    rejects([&]{ElectricVofRefinement::sample(*grid,other,charge);},"same mesh/clock foreign geometry owner accepted");++count;
    rejects([&]{values.validateFor(*grid,p,otherCharge);},"foreign equal-valued charge owner accepted");++count;
    const auto rz=geometry(*grid,domain(*grid,true));
    rejects([&]{ElectricVofRefinement::sample(*grid,rz,charge);},"different physical volume measure accepted");++count;
    const auto selected=geometry(*grid,domain(*grid,false,true));
    rejects([&]{ElectricVofRefinement::sample(*grid,selected,charge);},"different selected native domain accepted");++count;
    const auto wrongClock=geometry(*grid,d,{3,.05,.2,Clock::Position::End});
    rejects([&]{ElectricVofRefinement::sample(*grid,wrongClock,charge);},"equal-time different complete clock accepted");++count;
    const auto initial=geometry(*grid,d,{0,0,0,Clock::Position::Initial});
    const auto initialQ=ElectricVofCharge::fromIntegrated(*grid,initial,q);
    rejects([&]{ElectricVofRefinement::sample(*grid,initial,initialQ);},"indicator accepted Initial instead of End");++count;
    auto foreign=mesh(false,2);
    rejects([&]{ElectricVofRefinement::sample(*foreign,p,charge);},"equal-count different native mesh accepted");++count;
    auto huge=q;huge[0]=std::numeric_limits<double>::max();
    const auto hugeQ=ElectricVofCharge::fromIntegrated(*grid,p,huge);
    // On implementations where long double has binary64 range (MSVC), this
    // density is unrepresentable. Extended-range platforms can observe it.
    if(std::numeric_limits<long double>::max_exponent==std::numeric_limits<double>::max_exponent){
        rejects([&]{ElectricVofRefinement::sample(*grid,p,hugeQ);},"nonfinite physical density was accepted");++count;
    }
    auto tiny=q;tiny[0]=std::numeric_limits<double>::denorm_min();
    const auto tinyQ=ElectricVofCharge::fromIntegrated(*grid,p,tiny);
    check(ElectricVofRefinement::sample(*grid,p,tinyQ).maximumDensityJump()>0,
        "signed subnormal charge layer was erased by a density floor");
    values.validateFor(*grid,p,charge);
}
}
int main()
{
    try{
        std::size_t scenes=0,subfaces=0,rejections=0;
        constantsAndLinear(scenes,subfaces);signedLayer(scenes);axialRectangles(scenes);guards(rejections);
        std::cout<<"PASS electric_vof_refinement scenes="<<scenes<<" native_internal_subfaces="<<subfaces
            <<" rejection_guards="<<rejections<<" XY_RZ signed_density true_volume_centroids passive_End\n";
        return 0;
    }catch(const std::exception& error){std::cerr<<"FAIL electric_vof_refinement: "<<error.what()<<'\n';return 1;}
}
