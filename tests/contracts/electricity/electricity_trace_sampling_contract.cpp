#include "Core/IO/DictionaryMesh.h"
#include "Modules/Electricity/Detail/DielectricTraceSampling.h"
#include "Modules/Electricity/Detail/DielectricCubicTrace.h"
#include <iostream>

namespace {
void check(bool ok,const char* message) { if(!ok)throw std::runtime_error(message); }
}

void nativeElectricSamplingContract()
{
    namespace Sampling=Quantum::Electricity::Detail::DielectricTraceSampling;
    namespace P3=Quantum::Electricity::Detail::DielectricCubicTrace;
    using Snapshot=LevelSetHermiteGeometry::Snapshot;
    DictionaryMesh input;input.Nx=64;input.Ny=84;input.maxCellNum=12000;input.maxInitialIterations=0;
    input.maxCoord={4,5.2359877559829888};
    QuadtreeMesh mesh(input);FaceSegmentGraph faces;faces.build(mesh);
    Quantum::Geometry::ActiveHermitePatch patch;
    patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());
    // An attached narrowing column, with exactly one liquid-center row in
    // the uppermost interior horizontal face's original 6x5 stencil.
    // Quadratic phi is represented exactly apart from evaluation roundoff.
    for(int v=0;v<mesh.vertIndexNum();++v) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();
        const double z=p.y-input.maxCoord.y;
        patch.nodes[v]={p.x*p.x-.028*.028-.0078*z*z,2*p.x,-.0156*z,0};
    }
    auto source=std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(mesh,patch,true));
    Sampling::Lookup lookup(mesh,source);
    std::size_t checkedFaces=0;const FaceSegment* criticalFace=nullptr;
    for(const auto& face:faces.segments()) {
        if(face.L<0 || !face.isHorizontal || face.n.y<0)continue;
        const auto owner=Quantum::Mesh::cellAddress(mesh,mesh.cells(face.K));
        if(owner.i!=0 || owner.j<76)continue;
        const int row=static_cast<int>(owner.j)+1;
        const auto window=lookup.faceWindow(mesh,face);
        std::vector<int> expected;
        for(int j=row-3;j<=std::min(input.Ny-1,row+2);++j)
            for(int i=0;i<=2;++i)expected.push_back(mesh.cells(j*input.Nx+i).index());
        check(window.indices==expected,"native electric face stencil lost/added a nominal 6x5 member");
        ++checkedFaces;
        if(row!=83)continue;
        criticalFace=&face;
        std::size_t negative=0,oldNegative=0,oldCount=0;
        const auto bounds=source->nativeCellBounds(mesh.cells(face.K).index());
        const long double hx=static_cast<long double>(bounds[1].x)-bounds[0].x;
        const long double hy=static_cast<long double>(bounds[1].y)-bounds[0].y;
        for(const auto& sample:window.samples) {
            negative+=sample.phase==P3::Phase::Negative;
            if(std::abs(static_cast<long double>(sample.point.x)-face.faceCenter.x)<=2*hx
                && std::abs(static_cast<long double>(sample.point.y)-face.faceCenter.y)<=2.5L*hy) {
                ++oldCount;oldNegative+=sample.phase==P3::Phase::Negative;
            }
        }
        check(window.indices.size()==12 && negative==1 && oldNegative==0 && oldCount<12,
            "critical neck control does not expose the rounded-width phase loss");
        const double z=face.faceCenter.y-input.maxCoord.y;
        const double r=std::sqrt(.028*.028+.0078*z*z);
        const auto frame=source->jetAtNativeCellChart(mesh.cells(face.K).index(),
            Quantum::Geometry::cellChartFromLocal({r/(bounds[1].x-bounds[0].x),1}));
        const double ks=source->curvatureTangentialDerivativeAtNativeCellChart(frame.nativeCellIndex,frame.cellChartPoint);
        const P3::Reconstruction trace(frame.point,frame.normal,frame.curvature,ks,2,1,
            window.samples,std::max(mesh.cellSize.x,mesh.cellSize.y),true);
        const P3::Jump jump{-.2*frame.normal.y,-.2*frame.curvature*frame.normal.x,
            -.2*(ks*frame.normal.x-frame.curvature*frame.curvature*frame.normal.y)};
        const auto potentials=Sampling::evaluate(window,[](Vector2D p){return 3-.2*p.y;});
        const auto value=trace.evaluate(potentials,jump);
        check(std::abs(value.tangentialElectricField-.2*frame.normal.x)<1e-9
            && std::abs(value.normalDisplacement-.4*frame.normal.y)<1e-9,
            "restored native stencil failed the two-phase harmonic electric trace");
        const auto local=lookup.window(frame);
        std::vector<int> expectedLocal;
        for(int j=80;j<=83;++j)for(int i=0;i<=4;++i)
            expectedLocal.push_back(mesh.cells(j*input.Nx+i).index());
        check(local.indices==expectedLocal,"native axis-local 5x5 stencil lost nominal members");
        std::cout<<"Native dyadic electric stencil: nominal="<<window.indices.size()
            <<" negative="<<negative<<" old_rounded_count="<<oldCount<<" old_negative="<<oldNegative
            <<" P3_condition="<<trace.conditionEstimate()<<" local_count="<<local.indices.size()<<'\n';
    }
    check(criticalFace && checkedFaces==7,"native sampling control missed its physical faces");
    for(double widening:{.005,0.}) {
        // The first neck has resolved liquid farther into the domain, but no
        // liquid center in the boundary-truncated box. The constant subcell
        // cylinder has no liquid centers anywhere and MUST still reject.
        for(int v=0;v<mesh.vertIndexNum();++v) {
            const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();const double z=p.y-input.maxCoord.y;
            patch.nodes[v]={p.x*p.x-.025*.025-widening*z*z,2*p.x,-2*widening*z,0};
        }
        auto thin=std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(mesh,patch,true));
        Sampling::Lookup clipped(mesh,thin),oneSided(mesh,thin,Sampling::NativeBoundaryWindow::OneSided);
        const auto old=clipped.faceWindow(mesh,*criticalFace),full=oneSided.faceWindow(mesh,*criticalFace);
        const auto negative=[](const Sampling::Window& w){return std::count_if(w.samples.begin(),w.samples.end(),
            [](const auto& s){return s.phase==P3::Phase::Negative;});};
        std::vector<int> expected;
        for(int row=78;row<=83;++row)for(int column=0;column<5;++column)
            expected.push_back(mesh.cells(row*input.Nx+column).index());
        check(full.indices==expected && old.samples.size()==12 && negative(old)==0,
            "one-sided electric footprint is not the prescribed in-domain 6x5 box");
        const int owner=mesh.cells(criticalFace->K).index();const auto bounds=thin->nativeCellBounds(owner);
        const double z=criticalFace->faceCenter.y-input.maxCoord.y;
        const double r=std::sqrt(.025*.025+widening*z*z);
        const auto frame=thin->jetAtNativeCellChart(owner,
            Quantum::Geometry::cellChartFromLocal({r/(bounds[1].x-bounds[0].x),1}));
        const double ks=thin->curvatureTangentialDerivativeAtNativeCellChart(frame.nativeCellIndex,frame.cellChartPoint);
        const auto fit=[&](const Sampling::Window& w){return P3::Reconstruction(frame.point,frame.normal,
            frame.curvature,ks,2,1,w.samples,std::max(mesh.cellSize.x,mesh.cellSize.y),true);};
        const auto missing=[&](const Sampling::Window& w){
            try {fit(w);}catch(const std::exception& error){return std::string(error.what()).find("both phases")!=std::string::npos;}
            return false;
        };
        check(missing(old),"clipped electric neck control unexpectedly has two-phase support");
        if(widening==0) {
            check(negative(full)==0 && missing(full),"one-sided stencil fabricated an unresolved liquid potential");
            std::cout<<"Native one-sided electric stencil: entirely subcell cylinder still rejects missing phase\n";
            continue;
        }
        check(negative(full)==2,"one-sided neck did not retain its two real liquid centers");
        const auto trace=fit(full);
        constexpr double a=.1;
        const auto potential=[&](Vector2D p){const double z=p.y-input.maxCoord.y;return 3+a*(z*z-.5*p.x*p.x);};
        const Vector2D n=frame.normal,t{-n.y,n.x},electric{a*frame.point.x,-2*a*(frame.point.y-input.maxCoord.y)};
        const auto dot=[](Vector2D u,Vector2D v){return u.x*v.x+u.y*v.y;};
        const auto jac=[a](Vector2D v){return Vector2D{a*v.x,-2*a*v.y};};
        const double en=dot(electric,n),et=dot(electric,t),k=frame.curvature;
        const P3::Jump jump{-en,-(dot(jac(t),n)+k*et),
            -(-k*dot(jac(n),n)+2*k*dot(jac(t),t)+ks*et-k*k*en)};
        const auto value=trace.evaluate(Sampling::evaluate(full,potential),jump);
        const double error=std::max({std::abs(value.interfacePotential-potential(frame.point)),
            std::abs(value.tangentialElectricField-et),std::abs(value.normalDisplacement-2*en)});
        check(error<1e-9,"one-sided boundary electric field lost analytic RZ harmonic transmission");
        const auto local=oneSided.window(frame);
        check(local.indices.size()==25 && negative(local)==1,"one-sided local trace lost fixed 5x5 support");
        const auto localValue=fit(local).evaluate(Sampling::evaluate(local,potential),jump);
        check(std::abs(localValue.tangentialElectricField-et)<1e-9
            && std::abs(localValue.normalDisplacement-2*en)<1e-9,
            "one-sided local electric reconstruction lost the shared analytic trace");
        std::cout<<"Native one-sided electric stencil: face="<<full.indices.size()<<" negative="<<negative(full)
            <<" local="<<local.indices.size()<<" harmonic_error="<<error<<" condition="<<trace.conditionEstimate()<<'\n';
    }
}
