// Offline readers/observations for ordinary B2/B3 accepted output. No solver,
// restart, filtering or geometric reconstruction is performed here.
using System;
using System.IO;
using System.Linq;
using System.Collections.Generic;
using System.Globalization;
using System.Text.RegularExpressions;

public static class HydroStandardAnalysis {
    static readonly CultureInfo CI=CultureInfo.InvariantCulture;
    static double N(string s) { double v=double.Parse(s,CI); Check(!double.IsNaN(v)&&!double.IsInfinity(v),"Nonfinite observation input"); return v; }
    static void Check(bool ok,string message) { if(!ok)throw new InvalidDataException(message); }
    public sealed class Field {
        public double Time,PressureTime;
        public bool PressureAvailable;
        public double[][] Values;
        public string[] Names;
        public int[][] Cells;
        public double[] Column(string name) { int i=Array.IndexOf(Names,name);Check(i>=0,"Missing field "+name);return Values[i]; }
    }
    public static Field ReadField(string path) {
        string text=File.ReadAllText(path),header=text.Substring(0,text.IndexOf("ZONE ",StringComparison.Ordinal));
        var f=new Field();f.Names=Regex.Matches(header,"\"([^\"]+)\"").Cast<Match>().Select(m=>m.Groups[1].Value).ToArray();
        Match count=Regex.Match(text,@"NODES=(\d+), ELEMENTS=(\d+), ZONETYPE=FEQUADRILATERAL, DATAPACKING=BLOCK");
        Check(count.Success&&f.Names.Length>=2,"Expected one native FE/BLOCK quadrilateral field");
        int nodes=int.Parse(count.Groups[1].Value),cells=int.Parse(count.Groups[2].Value);
        Check(nodes==4*cells,"Expected independent native cell vertices");
        f.Time=N(Regex.Match(text,@"SOLUTIONTIME=([^,\s]+)").Groups[1].Value);
        f.PressureTime=N(Regex.Match(text,"AUXDATA PressureTime=\"([^\"]+)\"").Groups[1].Value);
        f.PressureAvailable=Regex.IsMatch(text,"AUXDATA PressureAvailable=\"1\"");
        Check(text.Contains("nested 5x5 and 9x9 complete-field samples"),"B2/B3 requires complete velocity observations");
        Match location=Regex.Match(text,@"VARLOCATION=\(\[([\s\d,]+)\]=CELLCENTERED\)");
        Check(location.Success,"Missing cell field location");
        int[] ids=location.Groups[1].Value.Split(new[]{',','\r','\n',' ','\t'},StringSplitOptions.RemoveEmptyEntries).Select(int.Parse).ToArray();
        Check(ids.SequenceEqual(Enumerable.Range(3,f.Names.Length-2)),"Unexpected mixed nodal/cell field layout");
        string tail=text.Substring(location.Index+location.Length);
        tail=Regex.Replace(tail,@"(?m)^AUXDATA[^\r\n]*\r?$","");
        string[] tokens=tail.Split((char[])null,StringSplitOptions.RemoveEmptyEntries);
        Check(tokens.Length==2*nodes+(f.Names.Length-2)*cells+4*cells,"Truncated or extra field data");
        int at=0;f.Values=new double[f.Names.Length][];
        for(int k=0;k<f.Names.Length;k++){f.Values[k]=new double[k<2?nodes:cells];for(int j=0;j<f.Values[k].Length;j++)f.Values[k][j]=N(tokens[at++]);}
        f.Cells=new int[cells][];
        for(int k=0;k<cells;k++){f.Cells[k]=new int[4];for(int j=0;j<4;j++){int id=int.Parse(tokens[at++])-1;Check(id>=0&&id<nodes,"Bad field connectivity");f.Cells[k][j]=id;}}
        return f;
    }
    public sealed class StaticFrame {
        public double Time,PressureTime,MeanInside=double.NaN,MeanOutside=double.NaN,PressureJump=double.NaN;
        public double InnerVolume,OuterVolume,SpeedCenter,Speed5,Speed9,MaximumCellWidth;
        public int InsideCells,OutsideCells;
        public bool PressureAvailable;
    }
    public static StaticFrame Static(Field f,bool rz,double radius) {
        var r=new StaticFrame();r.Time=f.Time;r.PressureTime=f.PressureTime;r.PressureAvailable=f.PressureAvailable;
        r.SpeedCenter=f.Column("speed").Max();r.Speed5=f.Column("speed_sampled_max_5").Max();r.Speed9=f.Column("speed_sampled_max_9").Max();
        Check(r.Speed9>=r.Speed5&&r.Speed5>=r.SpeedCenter*(1-1e-12),"Inconsistent nested speed samples");
        double[] x=f.Column("x"),y=f.Column("y"),p=f.Column("p_volume_mean"),phase=f.Column("liquid_fraction_geometry");
        double pin=0,pout=0;
        for(int i=0;i<f.Cells.Length;i++){
            int[] ids=f.Cells[i];double lo=ids.Min(k=>x[k]),hi=ids.Max(k=>x[k]),bottom=ids.Min(k=>y[k]),top=ids.Max(k=>y[k]);
            double h=Math.Max(hi-lo,top-bottom),cx=.5*(lo+hi),cy=.5*(bottom+top),d=Math.Sqrt(cx*cx+cy*cy);
            double v=(hi-lo)*(top-bottom)*(rz?Math.PI*(lo+hi):1);Check(v>0,"Nonpositive physical cell volume");
            r.MaximumCellWidth=Math.Max(r.MaximumCellWidth,h);
            // Fixed analytical reference region, with an endpoint purity check.
            // Pressure itself belongs to the solved midpoint, never t=0.
            if(d<=radius-2*h&&Math.Abs(phase[i]-1)<=1e-10){r.InsideCells++;r.InnerVolume+=v;pin+=v*p[i];}
            if(d>=radius+2*h&&Math.Abs(phase[i])<=1e-10){r.OutsideCells++;r.OuterVolume+=v;pout+=v*p[i];}
        }
        if(r.PressureAvailable&&r.InnerVolume>0&&r.OuterVolume>0){r.MeanInside=pin/r.InnerVolume;r.MeanOutside=pout/r.OuterVolume;r.PressureJump=r.MeanInside-r.MeanOutside;}
        return r;
    }
    sealed class Arc {
        public double A,B;public double[] X=new double[5],Y=new double[5];
        public double V(double t,bool y,int d=0){double value=0;double[] p=y?Y:X;for(int k=4;k>=d;k--){double factor=1;for(int j=0;j<d;j++)factor*=k-j;value=value*t+factor*p[k];}return value;}
    }
    public sealed class Geometry {
        public double Time,Measure,Volume,CurvatureRelativeL2,RadialDistanceRelativeL2,QuadratureEstimate;
        public double MeanRadius=double.NaN,Mode2Cos=double.NaN,Mode2Sin=double.NaN,MirrorRelativeL2=double.NaN,PolarSamplingDifference=double.NaN;
        public bool PolarAvailable;public string PolarReason="not_requested";
        public int Step,Arcs;public double[] PolarRadii;
    }
    static readonly double[] G={-.9602898564975363,-.7966664774136267,-.525532409916329,-.1834346424956498,.1834346424956498,.525532409916329,.7966664774136267,.9602898564975363};
    static readonly double[] W={.1012285362903763,.2223810344533745,.3137066458778873,.362683783378362,.362683783378362,.3137066458778873,.2223810344533745,.1012285362903763};
    static double[] Rule(Arc a,double lo,double hi,bool rz,double radius) {
        var sum=new double[4];double kref=(rz?2:1)/radius;
        for(int i=0;i<G.Length;i++){
            double t=.5*(lo+hi+(hi-lo)*G[i]),x=a.V(t,false),y=a.V(t,true),dx=a.V(t,false,1),dy=a.V(t,true,1);
            double speed=Math.Sqrt(dx*dx+dy*dy);Check(speed>0,"Nonregular accepted arc");
            double k=(dx*a.V(t,true,2)-dy*a.V(t,false,2))/(speed*speed*speed);
            if(rz){Check(x>0,"Nonpositive radius inside an RZ free arc");k+=dy/(speed*x);}
            double w=.5*(hi-lo)*W[i],measure=w*speed*(rz?2*Math.PI*x:1),d=(Math.Sqrt(x*x+y*y)-radius)/radius;
            sum[0]+=measure;sum[1]+=measure*Math.Pow(k/kref-1,2);sum[2]+=measure*d*d;
            sum[3]+=rz?w*Math.PI*x*x*dy:w*.5*(x*dy-y*dx);
        }
        return sum;
    }
    static double[] Integrate(Arc a,double lo,double hi,bool rz,double radius,double budget,ref double estimate,int depth) {
        double mid=.5*(lo+hi);double[] coarse=Rule(a,lo,hi,rz,radius),left=Rule(a,lo,mid,rz,radius),right=Rule(a,mid,hi,rz,radius);
        double error=0;for(int k=0;k<4;k++){left[k]+=right[k];error=Math.Max(error,Math.Abs(left[k]-coarse[k]));}
        Check(!double.IsNaN(error)&&!double.IsInfinity(error),"Nonfinite geometry integral");
        if(error<=budget){estimate+=error;return left;}
        Check(depth<20&&mid>lo&&mid<hi,"Offline geometry quadrature unresolved");
        left=Integrate(a,lo,mid,rz,radius,budget/2,ref estimate,depth+1);right=Integrate(a,mid,hi,rz,radius,budget/2,ref estimate,depth+1);
        for(int k=0;k<4;k++)left[k]+=right[k];return left;
    }
    static double Eval(double[] p,double t){double v=0;for(int k=p.Length-1;k>=0;k--)v=v*t+p[k];return v;}
    static double[] Derivative(double[] p){if(p.Length==1)return new[]{0.0};var d=new double[p.Length-1];for(int k=1;k<p.Length;k++)d[k-1]=k*p[k];return d;}
    static List<double> Roots(double[] p) {
        int degree=p.Length-1;while(degree>0&&p[degree]==0)degree--;
        var result=new List<double>();if(degree==0)return result;Array.Resize(ref p,degree+1);
        var cuts=Roots(Derivative(p));cuts.Insert(0,0);cuts.Add(1);
        for(int i=1;i<cuts.Count;i++){
            double a=cuts[i-1],b=cuts[i],fa=Eval(p,a),fb=Eval(p,b);if(fa==0&&a>0&&a<1)result.Add(a);
            if((fa<0&&fb>0)||(fa>0&&fb<0)){
                for(int j=0;j<60;j++){double m=.5*(a+b),fm=Eval(p,m);if(m==a||m==b)break;if(fm==0){a=b=m;break;}if((fa<0)==(fm<0)){a=m;fa=fm;}else b=m;}
                result.Add(.5*(a+b));
            }
        }
        result.Sort();return result.Distinct().ToList();
    }
    sealed class AngularPiece {public Arc Arc;public double Lo,Hi,ThetaLo,ThetaHi;}
    static double PositiveAngle(double angle){double two=2*Math.PI;angle%=two;if(angle<0)angle+=two;return angle;}
    static double[] Polar(List<Arc> arcs,int count) {
        var pieces=new List<AngularPiece>();double theta0=Math.Atan2(arcs[0].V(0,true),arcs[0].V(0,false)),theta=theta0;
        foreach(Arc a in arcs){
            // d(theta)/dt has sign cross(X,X'). Check its extrema using
            // derivative roots of the actual polynomial, not point sorting.
            var determinant=new double[8];for(int i=0;i<5;i++)for(int j=1;j<5;j++)determinant[i+j-1]+=j*(a.X[i]*a.Y[j]-a.Y[i]*a.X[j]);
            var critical=Roots(Derivative(determinant));critical.Add(0);critical.Add(1);
            Check(critical.Min(t=>Eval(determinant,t))>0,"Polar mode needs one regular star-shaped curve about the reference center");
            for(int j=0;j<16;j++){
                double lo=j/16.0,hi=(j+1)/16.0;
                double al=Math.Atan2(a.V(lo,true),a.V(lo,false)),ah=Math.Atan2(a.V(hi,true),a.V(hi,false));
                double delta=PositiveAngle(ah-al);Check(delta>0&&delta<Math.PI,"Unresolved polar arc angle");
                pieces.Add(new AngularPiece{Arc=a,Lo=lo,Hi=hi,ThetaLo=theta,ThetaHi=theta+delta});theta+=delta;
            }
        }
        Check(Math.Abs(theta-theta0-2*Math.PI)<1e-8,"Polar curve winding is not one");
        var radii=new double[count];
        for(int j=0;j<count;j++){
            double target=theta0+PositiveAngle(2*Math.PI*j/count-theta0);int low=0,high=pieces.Count-1;
            while(low<high){int mid=(low+high)/2;if(pieces[mid].ThetaHi<target)low=mid+1;else high=mid;}
            var p=pieces[low];double a=p.Lo,b=p.Hi,anchor=Math.Atan2(p.Arc.V(p.Lo,true),p.Arc.V(p.Lo,false));
            for(int iteration=0;iteration<52;iteration++){
                double mid=.5*(a+b),current=p.ThetaLo+PositiveAngle(Math.Atan2(p.Arc.V(mid,true),p.Arc.V(mid,false))-anchor);
                if(current<target)a=mid;else b=mid;
            }
            double t=.5*(a+b),x=p.Arc.V(t,false),y=p.Arc.V(t,true);radii[j]=Math.Sqrt(x*x+y*y);
        }
        return radii;
    }
    static double[] Modes(double[] r,double radius,int stride) {
        int n=r.Length/stride;double mean=0,c=0,s=0,mirror=0;
        for(int j=0;j<n;j++){
            double v=r[j*stride],theta=2*Math.PI*j/n;mean+=v;c+=v*Math.Cos(2*theta);s+=v*Math.Sin(2*theta);
            double dx=v-r[((n-j)%n)*stride],dy=v-r[((n/2-j+n)%n)*stride];mirror+=.5*(dx*dx+dy*dy);
        }
        return new[]{mean/n,2*c/n/radius,2*s/n/radius,Math.Sqrt(mirror/n)/radius};
    }
    public static Geometry ReadGeometry(string path,bool rz,double radius,bool modes) {
        string[] lines=File.ReadAllLines(path);Check(lines.Length>1,"No active interface in requested frame");
        string[] names=lines[0].Split(',');var columns=new Dictionary<string,int>();for(int i=0;i<names.Length;i++)columns.Add(names[i],i);
        var arcs=new List<Arc>();var r=new Geometry();string component=null;
        for(int i=1;i<lines.Length;i++){
            if(string.IsNullOrWhiteSpace(lines[i]))continue;string[] v=lines[i].Split(',');Check(v.Length==names.Length,"Truncated polynomial frame");
            int step=(int)N(v[columns["step"]]);double time=N(v[columns["time"]]);string id=v[columns["component"]];
            if(component==null){component=id;r.Step=step;r.Time=time;}
            Check(id==component&&step==r.Step&&time==r.Time,"Single-component accepted frame required");
            Check((v[columns["axisymmetric"]]=="1")==rz,"Geometry coordinate mismatch");
            Check(v[columns["begin_kind"]]==(rz?"axis":"periodic")&&v[columns["end_kind"]]==(rz?"axis":"periodic"),"Closed drop/axis meridian required");
            var a=new Arc{A=N(v[columns["a"]]),B=N(v[columns["b"]])};Check(a.B>a.A,"Invalid material interval");
            for(int k=0;k<5;k++){a.X[k]=N(v[columns["p"+k+"x"]]);a.Y[k]=N(v[columns["p"+k+"y"]]);}arcs.Add(a);
        }
        arcs=arcs.OrderBy(a=>a.A).ToList();Check(arcs[0].A==0&&arcs[arcs.Count-1].B==1,"Incomplete material interval");
        for(int i=1;i<arcs.Count;i++)Check(arcs[i-1].B==arcs[i].A,"Material interval gap/overlap");
        var sum=new double[4];double error=0;
        foreach(Arc a in arcs)for(int part=0;part<4;part++){
            double[] q=Integrate(a,part/4.0,(part+1)/4.0,rz,radius,1e-11*(a.B-a.A)/4,ref error,0);
            for(int k=0;k<4;k++)sum[k]+=q[k];
        }
        Check(sum[0]>0&&sum[3]>0,"Nonpositive oriented surface/volume");
        r.Arcs=arcs.Count;r.Measure=sum[0];r.Volume=sum[3];r.CurvatureRelativeL2=Math.Sqrt(sum[1]/sum[0]);r.RadialDistanceRelativeL2=Math.Sqrt(sum[2]/sum[0]);r.QuadratureEstimate=error;
        if(modes){try{
            r.PolarRadii=Polar(arcs,1024);double[] fine=Modes(r.PolarRadii,radius,1),coarse=Modes(r.PolarRadii,radius,2);
            r.MeanRadius=fine[0];r.Mode2Cos=fine[1];r.Mode2Sin=fine[2];r.MirrorRelativeL2=fine[3];
            r.PolarSamplingDifference=Math.Max(Math.Abs(fine[0]-coarse[0])/radius,Math.Max(Math.Abs(fine[1]-coarse[1]),Math.Max(Math.Abs(fine[2]-coarse[2]),Math.Abs(fine[3]-coarse[3]))));
            r.PolarAvailable=true;r.PolarReason="ok";
        }catch(InvalidDataException e){r.PolarAvailable=false;r.PolarReason=e.Message;}}
        return r;
    }
    public sealed class Fit {
        public bool Available,AtSearchBoundary;public string Reason;
        public double Omega,Damping,Amplitude,RelativeRms,Begin,End,FrequencyLower,FrequencyUpper,DampingLower,DampingUpper;
    }
    static double[] FitCost(double[] t,double[] y,double omega,double damping) {
        double yy=0,yf=0,ff=0;for(int i=0;i<t.Length;i++){
            double f=Math.Exp(-damping*(t[i]-t[0]))*(1-Math.Cos(2*omega*t[i]));yy+=y[i]*y[i];yf+=y[i]*f;ff+=f*f;
        }
        double amplitude=ff>0?yf/ff:0,error=0;
        for(int i=0;i<t.Length;i++){double f=amplitude*Math.Exp(-damping*(t[i]-t[0]))*(1-Math.Cos(2*omega*t[i]));error+=Math.Pow(y[i]-f,2);}
        return new[]{error,amplitude,yy};
    }
    public sealed class ReferenceMetric {
        public string Name; public bool Available;
        public double AbsoluteRms,RelativeL2,MaximumSampledError,ComputedMinimum,ComputedMinimumTime,
            ComputedMaximum,ComputedMaximumTime,ReferenceMinimum,ReferenceMinimumTime,ReferenceMaximum,ReferenceMaximumTime,
            ComputedLast,ReferenceLast;
    }
    public sealed class ReferenceComparison {
        public double Begin,End; public int Samples;
        public double[][] Aligned; public ReferenceMetric[] Metrics;
    }
    // Original Hysing files: t, area, circularity, centroid_y, mean_u_y.
    // Align only at accepted simulation endpoints inside the supplied common
    // interval. No extrapolation, shifted clock or smoothing. Time norms use
    // trapezoidal quadrature of the sampled squared differences.
    public static ReferenceComparison CompareHysing(string path,double[] time,double[][] value,double begin,double end) {
        Check(time.Length>=2&&value.Length==4,"Incomplete Hysing observations");
        for(int k=0;k<4;k++)Check(value[k].Length==time.Length,"Mismatched Hysing clocks");
        var reference=new List<double[]>();
        foreach(string line in File.ReadLines(path)) {
            if(string.IsNullOrWhiteSpace(line))continue;
            string[] parts=line.Split((char[])null,StringSplitOptions.RemoveEmptyEntries);
            Check(parts.Length==5,"Hysing reference requires five columns");
            var row=new double[5];for(int k=0;k<5;k++)row[k]=N(parts[k]);
            Check(row[1]>=0&&(reference.Count==0||row[0]>reference[reference.Count-1][0]),"Invalid reference area or clock");reference.Add(row);
        }
        Check(reference.Count>=2,"Empty Hysing reference");
        begin=Math.Max(begin,reference[0][0]);end=Math.Min(end,reference[reference.Count-1][0]);
        var aligned=new List<double[]>();int bracket=0;
        for(int i=0;i<time.Length;i++) {
            Check(!double.IsNaN(time[i])&&!double.IsInfinity(time[i])&&(i==0||time[i]>time[i-1]),"Invalid accepted Hysing clock");
            for(int k=0;k<4;k++)Check(!double.IsNaN(value[k][i])&&!double.IsInfinity(value[k][i]),"Nonfinite Hysing quantity");
            if(time[i]<begin||time[i]>end)continue;
            while(bracket+1<reference.Count-1&&reference[bracket+1][0]<time[i])bracket++;
            double fraction=(time[i]-reference[bracket][0])/(reference[bracket+1][0]-reference[bracket][0]);
            Check(fraction>=0&&fraction<=1,"Hysing comparison would extrapolate");
            var row=new double[9];row[0]=time[i];for(int k=0;k<4;k++){
                row[1+k]=value[k][i];row[5+k]=(1-fraction)*reference[bracket][k+1]+fraction*reference[bracket+1][k+1];}aligned.Add(row);
        }
        Check(aligned.Count>=2,"No common interval with two accepted Hysing endpoints");
        var result=new ReferenceComparison {Begin=aligned[0][0],End=aligned[aligned.Count-1][0],Samples=aligned.Count,
            Aligned=aligned.ToArray(),Metrics=new ReferenceMetric[4]};
        string[] names={"area","circularity","centroid_y","mean_u_y"};bool areaAvailable=reference.Exists(row=>row[1]>0);
        for(int k=0;k<4;k++) {
            var metric=new ReferenceMetric {Name=names[k],Available=k!=0||areaAvailable};result.Metrics[k]=metric;if(!metric.Available)continue;
            metric.ComputedMinimum=metric.ReferenceMinimum=double.PositiveInfinity;metric.ComputedMaximum=metric.ReferenceMaximum=double.NegativeInfinity;
            double square=0,norm=0;
            for(int i=0;i<aligned.Count;i++) {
                double[] row=aligned[i];double a=row[1+k],b=row[5+k],error=a-b;
                metric.MaximumSampledError=Math.Max(metric.MaximumSampledError,Math.Abs(error));
                if(a<metric.ComputedMinimum){metric.ComputedMinimum=a;metric.ComputedMinimumTime=row[0];}
                if(a>metric.ComputedMaximum){metric.ComputedMaximum=a;metric.ComputedMaximumTime=row[0];}
                if(b<metric.ReferenceMinimum){metric.ReferenceMinimum=b;metric.ReferenceMinimumTime=row[0];}
                if(b>metric.ReferenceMaximum){metric.ReferenceMaximum=b;metric.ReferenceMaximumTime=row[0];}
                if(i>0){double[] old=aligned[i-1];double w=.5*(row[0]-old[0]),e=old[1+k]-old[5+k];square+=w*(e*e+error*error);norm+=w*(old[5+k]*old[5+k]+b*b);}
                metric.ComputedLast=a;metric.ReferenceLast=b;
            }
            metric.AbsoluteRms=Math.Sqrt(square/(result.End-result.Begin));metric.RelativeL2=norm>0?Math.Sqrt(square/norm):double.NaN;
        }
        return result;
    }
    public static Fit FitEnergy(double[] t,double[] y,double referenceOmega) {
        var r=new Fit{Available=false,Reason="insufficient_full_window"};
        Check(t.Length==y.Length&&t.Length>0,"Fit samples differ");r.Begin=t[0];r.End=t[t.Length-1];
        if(t.Length<16||r.Begin!=0||r.End<4*Math.PI/referenceOmega)return r;
        double gap=0;for(int i=1;i<t.Length;i++){Check(t[i]>t[i-1],"Fit clock is not increasing");gap=Math.Max(gap,t[i]-t[i-1]);}
        if(2*referenceOmega*gap>Math.PI/2){r.Reason="insufficient_temporal_sampling";return r;}
        if(y.Max()<=0){r.Reason="no_kinetic_signal";return r;}
        r.FrequencyLower=.5*referenceOmega;r.FrequencyUpper=1.5*referenceOmega;
        r.DampingLower=-2/r.End;r.DampingUpper=10/r.End;
        double best=double.PositiveInfinity,wbest=0,bbest=0;
        for(int i=0;i<=200;i++)for(int j=0;j<=24;j++){
            double w=r.FrequencyLower+(r.FrequencyUpper-r.FrequencyLower)*i/200,b=r.DampingLower+(r.DampingUpper-r.DampingLower)*j/24;
            double error=FitCost(t,y,w,b)[0];if(error<best){best=error;wbest=w;bbest=b;}
        }
        double dw=(r.FrequencyUpper-r.FrequencyLower)/200,db=(r.DampingUpper-r.DampingLower)/24;
        for(int iteration=0;iteration<100;iteration++){
            double nextW=wbest,nextB=bbest;
            for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++){
                double w=wbest+i*dw,b=bbest+j*db;if(w<r.FrequencyLower||w>r.FrequencyUpper||b<r.DampingLower||b>r.DampingUpper)continue;
                double error=FitCost(t,y,w,b)[0];if(error<best){best=error;nextW=w;nextB=b;}
            }
            if(nextW==wbest&&nextB==bbest){dw*=.5;db*=.5;}else{wbest=nextW;bbest=nextB;}
            if(dw<1e-9*referenceOmega&&db<1e-9/r.End)break;
        }
        double[] final=FitCost(t,y,wbest,bbest);r.Omega=wbest;r.Damping=bbest;r.Amplitude=final[1];r.RelativeRms=Math.Sqrt(final[0]/final[2]);
        r.AtSearchBoundary=wbest-r.FrequencyLower<1e-5*referenceOmega||r.FrequencyUpper-wbest<1e-5*referenceOmega||bbest-r.DampingLower<1e-5/r.End||r.DampingUpper-bbest<1e-5/r.End;
        r.Available=true;r.Reason=r.AtSearchBoundary?"search_boundary_not_accepted":"diagnostic_fit";return r;
    }
}
