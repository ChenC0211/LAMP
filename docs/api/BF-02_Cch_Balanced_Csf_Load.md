# BF-02：同一压力 C 的恒曲率 CSF 实际载荷

`CchBalancedCsfLoad` 归属唯一 `lamp_hydrodynamics_cch`，持有不可变 `CchPressureLoad` 和 `CchRegularizedColor`。本接口当前仅支持给定常曲率；不替换旧 `CchInterfaceNormalLoad`，不推进界面、材料、速度或压力。BF-02规定范围已PASS，实际证据见 [报告](../task/BF-02_Report.md)；真实变量曲率和完整动力学仍待后续gate。

```cpp
auto pressure = CchPressureLoad::build(phase.metrics().divergence());
auto color = CchRegularizedColor::build(mesh, phase, colorConfig);
auto csf = CchBalancedCsfLoad::constantCurvature(mesh, pressure, color, sigma, kappa0);
csf.validateFor(mesh, pressure, color);
// Actual RHS of A*u+C*p=rhs: ADD, with no extra minus sign.
const auto& S = csf.momentumRhsLoad();
```

原点查询、epsilon、相几何和求积说明见 [BF-01 API](BF-01_Cch_Regularized_Color.md)。唯一构造是

\[
g=Cc_\epsilon,\qquad S_\sigma=+\sigma\kappa_0 g,
\quad C=-D_I^T.
\]

`colorGradient()` 返回原 C 的积分面载荷，`momentumRhsLoad()` 返回可直接加到既有 saddle RHS 的物理载荷。不能再次乘面面积、除以密度或引入另一套 gradient；AMR 真子面与 RZ 物理面积均由原 C 承担。`phi<0` 为液相，凸滴 kappa0>0 时静态液内压力较高。旧 H12/H13 `balanceLoads()` 的调用负号不适用于此 API。

所有 native 行均保留，包括固定壁和面积零的 RZ 轴面。零 sigma 或零 kappa0 时 RHS 每行精确为0，仍保留真实 Cc 诊断；sigma 必须有限非负，kappa0 必须有限且可正可负。系数 sigma*kappa0 溢出或两个非零输入的乘积下溢成零时明确拒绝，逐面载荷溢出亦拒绝；不引入系数floor。

## 不可变来源

构造检查同一 layout identity、原 divergence matrix owner、合法原 mesh、geometry 和完整 clock，再调用已有 typed `CchNativeValues` / `CchPressureLoad::apply`。对象及其按值副本共享只读存储；`pressureLoad()`、`color()`、`layout()`、`geometryIdentity()`、`timeLayer()` 返回保留来源，sigma/kappa getter返回实际参数。

`validateFor(mesh, actualPressure, identity, clock)` 要求同一个压力 matrix owner，并复用 color 的 mesh/source/full-clock 检查。`validateFor(mesh, actualPressure, actualColor)` 另外要求同一不可变 color owner，故同 publication 下不同 epsilon、甚至重新构造的等值 color 都不会被误当成已绑定输入；共享存储的副本合法。

从同一 divergence 新建 C 可用于新载荷构造；已构造载荷对另一个 C owner 的运行校验仍拒绝。regrid/publication/time改变后应以新的同源 phase/color/C 重建载荷，不可重贴身份。

## 自由行、压力 gauge 与壁反力

对常曲率 beta=sigma*kappa0，每个真实连通压力分量的闭域参考解是

\[
p_K^*=\beta(c_K-\overline c_{V,component}),\quad u=0.
\]

它仅用作测试参考；生产工厂不写入 pressure，测试必须实际解原 saddle。自由行的 `Cp-S` 应接近舍入尺度；固定壁行可有反力

\[
R_{wall}=(Au+Cp-S)_{wall}
=-\beta\overline c_{V,component}(C1_{component})_{wall}.
\]

不能要求壁面 C1=0，也不能把去 gauge 后的完整 Cp-S 当作自由行平衡残差。纯气 c=0 全行载荷为零；纯液 c=1 的自由行精确为零，固定边界仍可能有 C1 载荷。闭域固定速度为零时壁反力功为零。开放/牵引边界及真实变量曲率支撑不是 BF-02 的能力声明。
