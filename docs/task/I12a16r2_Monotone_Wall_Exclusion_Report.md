# I12a16r2 / 原RK2壁面单调性报告

**PASS（本项）。原a16驱动、I12 / Gate I仍FAIL。** 见 [TODO](I12a16r2_Monotone_Wall_Exclusion_TODO.md) 与 [证明](../api/I12a16r2_Monotone_Wall_Exclusion.md)。本项没有改变原tracer、库存装配/reflux、几何闭合或误差预算。

## 结果

只修改1production+1test：原PolynomialVelocitySlab.cpp的noNewWall判据，以及原bounded_spacetime_flux_contract。对target=P−W*d，d>=0，W<=0可严格排除从内部新触及lower wall，W>=0可排除upper wall。IEEE舍入及零/下溢保留这个弱单调关系。原初始无hit、域内midpoint/target、CFL、轴及覆盖条件全部先行；未认证方向仍用原空间界。没有修改返回的数值bounds，也没有以积分E小来忽略pending。

新增18个完整闭参数区域（XY四边、RZ南北边×常速/变空间与时间速度/亚正规速度），包含旧时间端点、壁面起点和内部点。4,050条原actual tracer对照检查branch、B、J及坐标方向；8项真实初始/随后触壁拒绝均通过。原57主积分、r1的28边界积分/1,134trace/48检查完整输出前缀逐字不变。contract用时.6415153s、exit0。

原四个rotation/vortex×uniform/AMR场景都已越过此前内部face1的new-wall-contact拒绝，依次完成face0–3的积分，然后在 **nativeface4、negative3、positive4** 被初始边界hit排除证书拒绝。四者仍为半步、41panels、pending1、depth40：rotation E=5.7006049186674472e-16，vortex E=2.8527199994097243e-16，面预算仍1.5625e-8。

smoke2.9921844s、full2.9644777s，均exit1、completed0/planned4。full保留原30/400步请求。实际mid/end预测、首边界面history认证保持通过；accepted时间0及旧identity/M保持。未执行到半步装配/闭合，不能把局部算子进展当成完整步骤。

## 回归与归档

39targets首次构建47.4745784s、exit0、无warning。26CTest161.5242398s，25PASS+原driver真实FAIL，未设WILL_FAIL。24项旧完整输出、上述扩展contract的整个旧前缀、旧crossing probe的exit1与输出全部保持。

架构仍12targets、383TU、908files、3118includes、5exact legacy edges；81frozen、3PDF及r1 archive SHA不变。首次审计PASS，源变更范围1+1，无新增production文件。归档交付2源文件、TODO/证明/本文/入口、39runtime及原DLL，并核对每项payload SHA。

证据：[验收](../../runs/baselines/I12a16r2_20261002/verification.json)、[驱动](../../runs/baselines/I12a16r2_20261002/driver_verification.json)、[回归](../../runs/baselines/I12a16r2_20261002/regression_comparison.json)、[范围](../../runs/baselines/I12a16r2_20261002/source_audit.json)。

## 下一项

面4位于x=.5、y∈[0,.125]。原rotation/vortex在该闭面上的uy计算包含x−.5，实际为零；旧一般区间normalize把退化中心坐标的减法包成跨零小区间，导致在南壁端点无法排除初始hit。下一项需认证原polynomial求值的精确零/符号事实，不能凭采样或放宽pending规则处理。随后仍须认证真实后续边界事件、完成同四个first steps及原时间/空间/AMR/曲率gate。独立M与几何预算的用户授权分离继续保持，未进入Hydro。
