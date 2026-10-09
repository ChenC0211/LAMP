#include "Modules/Interface/Compatibility/Detail/InterfaceParametricWorldSheet.h"
#include <optional>

namespace InterfaceParametricTransport::WorldSheet {
namespace {
double magnitude(Bounds b){return std::max(std::abs(b.lo),std::abs(b.hi));}
double speed(Bounds x,Bounds y){return std::hypot(magnitude(x),magnitude(y));}
struct Panel {
    std::optional<Chart> chart;
    Sheet critical;
    std::array<double,2> bound{},value{};
    double error=0;int depth=0,axis=0,side=0;
};
}
SurfaceFlux measureSurfaceFlux(const std::vector<Sheet>& sheets,const std::vector<SheetBounds>& bounds,
    Sampler& sampler,const InterfaceParametricProjection::Grid& grid){
    SurfaceFlux result;
    const double pi=std::acos(-1.),radius=grid.left+grid.width;
    // One absolute surface-measure budget per transported component. The
    // combined report sums both estimates and budgets across components.
    result.budget=1e-10*std::max(grid.width,grid.height)*(grid.axisymmetric?2*pi*radius:1.);
    std::vector<Panel> panels;FluxReport charts;
    const auto plane=[&](int axis,int side){return axis?grid.bottom+side*grid.height:grid.left+side*grid.width;};
    const auto limits=[&](int axis){return axis?std::array<double,2>{grid.left,radius}:std::array<double,2>{grid.bottom,grid.bottom+grid.height};};
    const auto weightBound=[&](int axis,int side){return grid.axisymmetric?2*pi*(axis?radius:plane(axis,side)):1.;};
    const auto evaluate=[&](Panel p){
        const auto& chart=*p.chart;const auto extent=chart.envelope(chart.a,chart.b).z;const auto interval=limits(p.axis);
        p.value={};p.error=0;
        if(extent.hi<interval[0]||extent.lo>interval[1])return p;
        const auto integrand=[&](double q){
            const auto x=chart.curve(q);std::array<double,2> value{};
            if(x.z<interval[0]||x.z>interval[1])return value;
            const auto a=x.materialJet,b=x.timeJet;const double length=norm(a),an=p.axis?a.y:a.x,bn=p.axis?b.y:b.x;
            require(length>0&&std::isfinite(length),"degenerate material tangent in surface transport");
            const double w=(p.side?1.:-1.)*(grid.axisymmetric?2*pi*(p.axis?x.z:plane(p.axis,p.side)):1.);
            if(chart.independentTime){
                require(an!=0,"singular time chart in surface transport");
                value={w*length*bn/std::abs(an),w*std::copysign(1.,an)*dot(a,b)/length};
            }else{
                require(bn!=0,"singular material chart in surface transport");
                value={w*std::copysign(length,bn),w*std::copysign(1.,an)*dot(a,b)/length*std::abs(an/bn)};
            }
            require(std::isfinite(value[0])&&std::isfinite(value[1]),"nonfinite surface transport functional");return value;
        };
        const auto integrate=[&](double lo,double hi){std::array<long double,2> sum{};const auto& rule=gaussEight();
            for(int k=0;k<8;k++){const auto v=integrand(.5*(lo+hi+(hi-lo)*rule.x[k]));
                for(int j=0;j<2;j++)sum[j]+=.5L*(hi-lo)*rule.w[k]*v[j];}return sum;};
        const double mid=.5*(chart.a+chart.b);const auto coarse=integrate(chart.a,chart.b),left=integrate(chart.a,mid),right=integrate(mid,chart.b);
        for(int j=0;j<2;j++){p.value[j]=double(left[j]+right[j]);p.error+=double(std::abs(left[j]+right[j]-coarse[j]));}
        // A crossing of a physical corner is a cut in the integration domain,
        // not an event that a few quadrature points may silently miss.
        if(extent.lo<interval[0]||extent.hi>interval[1])
            p.error=std::max(p.error,2*(chart.b-chart.a)*(p.bound[0]+p.bound[1]));
        return p;
    };
    const auto append=[&](const Sheet& sheet,const SheetBounds& box,int axis,int side,int depth){
        const double coordinate=plane(axis,side),weight=weightBound(axis,side);if(weight==0)return;
        const auto f=axis?box.y:box.x,z=axis?box.x:box.y;const auto interval=limits(axis);
        const double pad=sheet.positionEnclosure+64*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(coordinate),std::abs(interval[0]),std::abs(interval[1])});
        if(f.lo>coordinate+pad||f.hi<coordinate-pad||z.hi<interval[0]-pad||z.lo>interval[1]+pad)return;
        const auto fs=axis?box.ys:box.xs,ft=axis?box.yt:box.xt;
        const double ds=speed(box.xs,box.ys),dt=speed(box.xt,box.yt);
        if(dt==0)return; // exactly stationary material: both functionals vanish
        if(!strict(fs)&&!strict(ft)){
            // The retained sheet has material degree eight and time degree
            // four. Coarea gives a numerical critical-patch estimate without
            // dividing by a vanishing crossing derivative. Repeated crossings
            // are counted with their polynomial multiplicity. As for the
            // sheet's jet enclosures, this is not an arbitrary-Flow certificate.
            Panel p;p.critical=sheet;p.axis=axis;p.side=side;p.depth=depth;
            p.bound={4*weight*ds,8*weight*dt};p.error=p.bound[0]+p.bound[1];panels.push_back(std::move(p));return;
        }
        for(auto chart:regularCharts(sheet,box,sampler,axis,coordinate,charts)){
            Panel p;p.axis=axis;p.side=side;p.depth=depth;
            const auto denominator=chart.independentTime?fs:ft;
            const double minimum=std::min(std::abs(denominator.lo),std::abs(denominator.hi));
            p.bound=chart.independentTime?std::array<double,2>{weight*ds*magnitude(ft)/minimum,weight*dt}:
                std::array<double,2>{weight*ds,weight*dt*magnitude(fs)/minimum};
            p.chart=std::move(chart);panels.push_back(evaluate(std::move(p)));
        }
    };
    for(int axis=0;axis<2;axis++)for(int side=0;side<2;side++)
        for(std::size_t k=0;k<sheets.size();k++)append(sheets[k],bounds[k],axis,side,0);
    for(int iteration=0;!panels.empty();iteration++){
        double error=0;std::size_t worst=0;
        for(std::size_t k=0;k<panels.size();k++){error+=panels[k].error;if(panels[k].error>panels[worst].error)worst=k;}
        if(error<=result.budget)break;
        const auto old=panels[worst];
        require(iteration<100000&&old.depth<52,"surface transport exhausted its quadrature/chart budget");
        panels[worst]=std::move(panels.back());panels.pop_back();++result.refinements;
        if(old.chart){
            const double mid=.5*(old.chart->a+old.chart->b);
            require(mid>old.chart->a&&mid<old.chart->b,"surface transport exhausted crossing resolution");
            auto first=old,last=old;first.chart->b=last.chart->a=mid;++first.depth;++last.depth;
            panels.push_back(evaluate(std::move(first)));panels.push_back(evaluate(std::move(last)));
        }else{
            const auto& s=old.critical;
            const auto box=enclose(s);const auto fs=old.axis?box.ys:box.xs,ft=old.axis?box.yt:box.xt;
            const bool material=magnitude(fs)>=magnitude(ft);
            const double mid=material?.5*(s.a+s.b):.5*(s.c+s.d);
            require(material?(mid>s.a&&mid<s.b):(mid>s.c&&mid<s.d),"surface critical patch exhausted material/time resolution");
            for(int half=0;half<2;half++){
                auto child=material?s.restrict(.5*half,.5*(half+1),0,1):s.restrict(0,1,.5*half,.5*(half+1));
                // Re-evaluate the SAME chart against the actual map locally.
                // A parent jet-error estimate must not impose a false floor
                // along an exactly pinned material endpoint.
                estimateSheet(child,sampler);append(child,enclose(child),old.axis,old.side,old.depth+1);
            }
        }
    }
    for(const auto& p:panels){result.advection+=p.value[0];result.endpointWork+=p.value[1];result.quadratureError+=p.error;
        if(!p.chart){++result.criticalPatches;result.criticalError+=p.error;}}
    result.timeCharts=charts.timeCharts;result.materialCharts=charts.materialCharts;
    require(std::isfinite(result.advection)&&std::isfinite(result.endpointWork)&&result.quadratureError<=result.budget,"invalid integrated surface transport");return result;
}
}
