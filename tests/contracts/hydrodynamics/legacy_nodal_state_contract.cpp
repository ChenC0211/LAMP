#include "Modules/Hydrodynamics/LegacyNodalState.h"
#include "Modules/Hydrodynamics/LegacyNodalForcing.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>

using NodalState = LegacyNodalState;
auto& geometryOf(NodalState& state) { return state.geometry(); }
void setEvaporation(NodalState& state,int label,double value) { state.setEvaporationVelocity(0,label,value); }
double evaporation(NodalState& state,int label,Time layer=Time::current) { return state.evaporationVelocity(0,layer)[label]; }
void setWallEvaporation(NodalState& state,double value) { state.setNoSlipEvaporationVelocity(0,value); }
double wallEvaporation(NodalState& state) { return state.noSlipEvaporationVelocity(0); }
void temporal(NodalState& state,const InterfaceNode& node) { state.temporalOutput("node",{},{},0,node,evaporation(state,node.label()),1.25,3,true); }
namespace {
void check(bool value,const char* why) { if(!value)throw std::runtime_error(why); }
void near(double actual,double expected,const char* why) { check(std::isfinite(actual) && std::abs(actual-expected)<1e-11,why); }
DictionaryGeometry circle()
{
    DictionaryGeometry input;
    input.points={{.73,.5},{.5,.73},{.27,.5},{.5,.27}};
    for(int i=0;i<4;++i) {
        Line line;line.points[0]=input.points[i];line.points[1]=input.points[(i+1)%4];
        line.setType(LineType::arc);line.setArcCenter({.5,.5});line.setRadius(.23);line.setRadian(.5);
        line.setGeoBoundary(GeometricBoundary::moving);line.setImmersed(true);input.lines.push_back(line);
    }
    input.levelSetName={"liquid"};input.innerPoints={{.5,.5}};input.moving={true};
    input.linesOnInterface={{{0,false},{1,false},{2,false},{3,false}}};return input;
}
}
int main(int argc,char** argv) try
{
    const auto root=argc>1?std::filesystem::path(argv[1]):std::filesystem::current_path()/
        ("legacy_nodal_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(root);std::filesystem::create_directories(root/"roundtrip");
    DictionaryMesh input;input.Nx=input.Ny=16;input.maxCellNum=1024;input.maxInitialIterations=0;
    input.ensureLevelSetCount(1);QuadtreeMesh mesh(input);
    auto dictionary=circle();ImmersedGeometry<QuadtreeMesh> immersed(dictionary,mesh);
    immersed.initializeLevelSets(mesh,0.,false);immersed.refreshLevelSetDerivedData(mesh,0.,false);
    immersed.checkCutInfo(mesh,true);
    NodalState state;state.initializeGeometry(mesh,immersed);
    auto& geometry=geometryOf(state);
    const auto& labels=geometry.nodeIndexToLabel[0];check(labels.size()>8,"actual closed interface nodes");
    LegacyNodalForcing::initializeCapillary(labels,geometry.allNodes(0),.125,state.Fxs[0],state.Fys[0]);
    for(int label:labels)setEvaporation(state,label,.5+.03125*label);
    VertField<double> ux(mesh,2.),uy(mesh,-.5),preUx(mesh,1.),preUy(mesh,.25);
    state.updateGeometry(mesh,immersed,ux,uy,preUx,preUy,.02,.01,2);
    check(geometry.preNodeIndexToLabel[0].size()==labels.size(),"first geometry history retained");
    for(int i=0;i<state.Uxs[0].size();++i) { near(state.Uxs[0].value[i],3.,"half-step extrapolated x velocity");near(state.Uys[0].value[i],-1.25,"half-step extrapolated y velocity"); }
    for(int label:labels) {
        near(evaporation(state,label),0.,"moving geometry resets current evaporation");
        near(evaporation(state,label,Time::previous),.5+.03125*label,"first evaporation history");
        setEvaporation(state,label,1.5+.03125*label);
    }
    state.preUxs=state.Uxs;state.preUys=state.Uys;state.prePreUxs=state.Uxs;state.prePreUys=state.Uys;
    state.updateGeometry(mesh,immersed,ux,uy,preUx,preUy,.02,.01,3);
    check(geometry.prePreNodeIndexToLabel[0].size()==labels.size(),"second geometry history retained");
    for(int label:labels) {
        near(evaporation(state,label),0.,"second moving geometry reset");
        near(evaporation(state,label,Time::previous),1.5+.03125*label,"previous evaporation history");
        near(evaporation(state,label,Time::prePrevious),.5+.03125*label,"second previous evaporation history");
        setEvaporation(state,label,2.5+.03125*label);
    }
    for(int i=0;i<state.Uxs[0].size();++i) {
        state.preUxs[0].value[i]=1.;state.preUys[0].value[i]=.25;
        state.prePreUxs[0].value[i]=-1.5;state.prePreUys[0].value[i]=.75;
    }
    std::ofstream output(root/"nodes.csv");output<<std::setprecision(17)<<"label,x,y,ds,nx,ny,curvature,ux,uy,fx,fy\n";
    for(int label:labels) {
        const auto& node=geometry.interfaceNode(0,label);const int i=node.index();
        output<<label<<','<<node.coord().x<<','<<node.coord().y<<','<<node.ds()<<','<<node.normal().x<<','<<node.normal().y<<','<<node.curvature()<<','<<state.Uxs[0].value[i]<<','<<state.Uys[0].value[i]<<','<<state.Fxs[0].value[i]<<','<<state.Fys[0].value[i]<<'\n';
    }
    const auto& node=geometry.interfaceNode(0,labels[labels.size()/2]);
    Field<double> scalar(static_cast<int>(labels.size()),2.75);
    for(Time layer:{Time::current,Time::previous,Time::prePrevious}) {
        near(state.interfaceInterpolation(node.coord(),scalar,0,false,false,layer,node.normal()),2.75,"constant scalar node interpolation");
        near(state.interfaceInterpolation(node.coord(),scalar,0,true,false,layer,node.normal()),0.,"constant tangential derivative");
        const double vx=layer==Time::current?3.:(layer==Time::previous?1.:-1.5);
        const double vy=layer==Time::current?-1.25:(layer==Time::previous?.25:.75);
        near(state.interfaceInterpolation(node.coord(),scalar,0,false,true,layer,node.normal()),-vx*node.normal().y+vy*node.normal().x,"tangential velocity history");
        Field<double> polynomial(static_cast<int>(labels.size()),0.);
        for(int label:labels) {
            const auto& related=geometry.interfaceNode(0,label);
            const auto offset=related.coord()-node.coord();
            const double eta=-offset.x*node.normal().y+offset.y*node.normal().x;
            polynomial.value[related.index()]=.5*eta*eta+1.25*eta+2.;
        }
        near(state.interfaceInterpolation(node.coord(),polynomial,0,false,false,layer,node.normal()),2.,"quadratic node interpolation");
        near(state.interfaceInterpolation(node.coord(),polynomial,0,true,false,layer,node.normal()),1.25,"quadratic tangential derivative");
    }
    // Exercise the legacy cache payload; initialization leaves this optional set empty.
    geometry.noSlipNodeNum=1;geometry.noSlipNodes().push_back(node);
    state.noSlipUx.initialize(1,1.5);state.noSlipUy.initialize(1,-2.);
    state.noSlipFx.initialize(1,3.5);state.noSlipFy.initialize(1,-4.);
    setWallEvaporation(state,4.25);
    // Preserve the legacy temporal writer's velocity slot while checking its node payload.
    state.Uxs.resize(2);state.Uys.resize(2);state.Uxs[1]=state.Uxs[0];state.Uys[1]=state.Uys[0];
    std::filesystem::create_directories(root/"output");
    const auto previousDirectory=std::filesystem::current_path();std::filesystem::current_path(root);
    state.output("nodes",0);temporal(state,node);std::filesystem::current_path(previousDirectory);
    std::ofstream(root/"node-size.txt")<<sizeof(InterfaceNode)<<'\n';
    const auto prefix=root.generic_string()+"/";
    state.writeCache(prefix);
    const auto verifyCache=[&](NodalState& restored) {
    auto& restoredGeometry=geometryOf(restored);
    check(restoredGeometry.nodeIndexToLabel==geometry.nodeIndexToLabel && restoredGeometry.prePreNodeIndexToLabel==geometry.prePreNodeIndexToLabel,"SIC2 label/history roundtrip");
    for(int label:labels) {
        const auto& a=geometry.interfaceNode(0,label);const auto& b=restoredGeometry.interfaceNode(0,label);
        near(b.coord().x,a.coord().x,"SIC2 coordinate");near(b.curvature(),a.curvature(),"SIC2 curvature");
        near(restored.Uxs[0].value[a.index()],state.Uxs[0].value[a.index()],"SIC2 velocity");
        near(restored.Uys[0].value[a.index()],state.Uys[0].value[a.index()],"SIC2 y velocity");
        near(restored.preUxs[0].value[a.index()],1.,"SIC2 previous x velocity");
        near(restored.preUys[0].value[a.index()],.25,"SIC2 previous y velocity");
        near(restored.prePreUxs[0].value[a.index()],-1.5,"SIC2 second previous x velocity");
        near(restored.prePreUys[0].value[a.index()],.75,"SIC2 second previous y velocity");
        near(restored.Fxs[0].value[a.index()],state.Fxs[0].value[a.index()],"SIC2 force");
        near(restored.Fys[0].value[a.index()],state.Fys[0].value[a.index()],"SIC2 y force");
        near(evaporation(restored,label),2.5+.03125*label,"SIC2 current evaporation");
        near(evaporation(restored,label,Time::previous),1.5+.03125*label,"SIC2 previous evaporation");
        near(evaporation(restored,label,Time::prePrevious),.5+.03125*label,"SIC2 second previous evaporation");
    }
    check(restoredGeometry.noSlipNodeNum==1 && restoredGeometry.noSlipNodes().size()==1,"SIC2 no-slip nodes");
    near(restored.noSlipUx.value[0],1.5,"SIC2 no-slip x velocity");near(restored.noSlipUy.value[0],-2.,"SIC2 no-slip y velocity");
    near(restored.noSlipFx.value[0],3.5,"SIC2 no-slip x force");near(restored.noSlipFy.value[0],-4.,"SIC2 no-slip y force");
    near(wallEvaporation(restored),4.25,"SIC2 wall evaporation");
    };
    NodalState restored;restored.readCache(argc>2?argv[2]:prefix);verifyCache(restored);
    restored.writeCache((root/"roundtrip").generic_string()+"/");
    NodalState legacy;legacy.readCache(QUANTUM_SIC2_FIXTURE_DIR);verifyCache(legacy);
    std::filesystem::create_directories(root/"fixture_roundtrip");
    legacy.writeCache((root/"fixture_roundtrip").generic_string()+"/");
    const auto bytes=[](const std::filesystem::path& path) {
        std::ifstream file(path,std::ios::binary);
        return std::vector<char>(std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>());
    };
    check(bytes(std::filesystem::path(QUANTUM_SIC2_FIXTURE_DIR)/"interface_cache.bin")==bytes(root/"fixture_roundtrip/interface_cache.bin"),"frozen legacy SIC2 byte roundtrip");
    auto fixedDictionary=circle();fixedDictionary.moving[0]=false;
    ImmersedGeometry<QuadtreeMesh> fixedGeometry(fixedDictionary,mesh);
    fixedGeometry.initializeLevelSets(mesh,0.,false);fixedGeometry.refreshLevelSetDerivedData(mesh,0.,false);fixedGeometry.checkCutInfo(mesh,true);
    NodalState fixed;fixed.initializeGeometry(mesh,fixedGeometry);
    const int fixedLabel=geometryOf(fixed).nodeIndexToLabel[0][0];setEvaporation(fixed,fixedLabel,7.25);
    for(int loop=2;loop<=3;++loop)fixed.updateGeometry(mesh,fixedGeometry,ux,uy,preUx,preUy,.02,.01,loop);
    for(Time layer:{Time::current,Time::previous,Time::prePrevious})near(evaporation(fixed,fixedLabel,layer),7.25,"fixed geometry retains evaporation at all layers");
    std::cout<<"Legacy nodal geometry, velocity/history, interpolation and SIC2 checks passed: "<<root.string()<<'\n';return 0;
}
catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
