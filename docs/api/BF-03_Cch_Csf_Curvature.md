# BF-03：同源真实曲率的局部 CSF 投影

当前实施/验收状态见 [BF-03 TODO](../task/BF-03_TODO.md)。`CchBalancedCsfLoad::capillary` 扩展同一个不可变load owner，不替换BF-02显式constant工厂，也不改变BF-01、旧H12/H13或压力/材料算子。

```cpp
auto csf = CchBalancedCsfLoad::capillary(mesh, pressure, color, sigma);
csf.validateFor(mesh, pressure, color);
const auto& rhs = csf.momentumRhsLoad(); // Actual positive RHS, add directly.
const auto& k = csf.curvatureAverages();
const auto& support = csf.supportReports();
// k[f] represents a projected curvature only when support[f].curvatureDefined.
```

## 唯一来源和固定物理核

从保留color的同一FirstJet publication读取实际 `quadrature(c,g)`；曲率为原 `q.sample.curvature`，权重为原 `q.physicalWeight`。该值已是RZ总曲率/真实旋转表面积权重，不再乘r或2*pi，不重新求导color、插值另一套phi或输入解析κ。

面中心来自真实native face，相邻pressure leaf物理对角线给出

\[
R_f=\epsilon+\max_{K\text{ adjacent }f}\sqrt{h_{x,K}^2+h_{y,K}^2},
\qquad W_f(x)=\max(0,1-|x-x_f|/R_f)\,V_f(x).
\]

`V_f` 要求q所属leaf在当前selected压力分量，且面中心到q的闭线段全在该分量的selected leaf闭包内。相同全局连通分量并不足以跨孔借曲率；实现复用Core exact expansion，将原端点差作为有理ray参数，以精确cross-product比较slab/union/gap，不以普通double商相等桥接正宽hole。对已验证不重叠native leaves，精确面积等于bbox面积时可证明域为凸闭矩形并直接判断；其它域使用精确slab。`selectedSegmentVisible` 是同一纯几何predicate的独立数学端口，不能覆盖factory实际selected域。固定物理桶只加速查询，完整访问实际R范围，未增加ring上限或远场扩支撑。该CPP使用严格浮点/禁LTO编译边界，未修改Core。

该规则并不证明所有 `phi/|grad phi|` 的force band都在到Γ的epsilon邻域内。若实际有效面没有本核内真实正支撑，明确拒绝，不增大R。固定物理epsilon加密逼近同hat连续平均；要验证Sharp一致性，使用epsilon=O(h)。边界/轴/selected截断会影响投影阶数，需依据独立积分报告。

## 正权平均与实际 RHS

\[
\kappa_f=\frac{\sum_qw_qW_f(x_q)\kappa_q}{\sum_qw_qW_f(x_q)},\qquad
g=Cc,\quad S_f=+\sigma\kappa_fg_f.
\]

`positiveCurvatureAverage(weights,curvatures)` 是本factory调用的纯数学kernel，可单独验证常量再现和正权平均。它不是几何/曲率覆盖端口：factory仅消费原publication。weights有限非负、至少一项正，values有限、非空同长度；zero weight忽略。literal常κ直接准确再现；其它输入按最大weight/曲率标度归一后用anchor-delta累加。归一化weight/非零κ/moment下溢时明确unsupported，不能因mixed sign把可表示的非零均值默默变0。没有曲率clip/cap或人工系数floor，不可表示的平均/系数/载荷明确拒绝。

g仅调用原 `CchPressureLoad::apply` 一次，仍是未消元的完整native rows；不乘第二面积/密度，不构造第二gradient。任何有真实支撑的非零g都照算，包括roundoff signal小的面。`relativeGradientRoundoff` 只形成审计标志，不抹掉弱力；color积分误差传播另作为估计报告。

## 缺支撑、壁反力与零项

| status | 曲率是否有定义 | 实际行为 |
| --- | --- | --- |
| Projected | 是 | 正权真实支撑得到κ；S=σ κ g，g为零时S仍精确零 |
| ZeroGradient | 否 | 无支撑且g精确零，S=0 |
| PureBoundaryReaction | 否 | 无支撑、prescribed边界且邻cell的BF-01 transitionVolumeUpper精确为0；g保留C1诊断，variable factory S=0 |
| SigmaZero | 否 | 合法来源/config后sigma=0，准确零RHS，投影未计算 |

无局部支撑的非零free g、或真实wall-band均拒绝；不能从远Γ猜一个κ。纯相没有band时自然得到物理零force，固定纯液C1项不被误当成几何界面。BF-02显式constant工厂仍保留其已验收完整wall stress；当前variable工厂的boundary-only诊断不代表κ=0已投影，也不支持尚未验收的open/traction边界。

`curvatureAverages` 中不可用面以0占位，必须看curvatureDefined/status；`constantCurvature()` 对任何capillary模式拒绝，包括sigma=0。σ0时所有projection元数据未计算，report.projectionSkippedSigmaZero明确为true；来源与配置检查仍执行。

## 查询、诊断和身份

config默认 `maxCandidateQueriesPerFace=262144`、`maxVisibilityCellQueriesPerFace=262144`、`relativeGradientRoundoff=256*machine_epsilon`。两种work上限须正，分别控制surface候选与ray leaf查询，耗尽拒绝；roundoff须有限且在(0,1)，仅审计。没有独立epsilon或可调kernel半径。

逐面报告radius/center/component、candidate/support/visibility成本、totalWeight/minmax/maxDistance、geometry-defined状态，以及boundary/near-boundary/near-axis/AMR标志。roundoff scale来自同一C行的绝对色函数载荷和；colorGradientErrorEstimate传播BF-01积分误差估计，不能冒称真误差上界。near-boundary/AMR标志是截断/可能降阶诊断，不参与曲率过滤。

geometry/full clock、同一D/C/color owner校验沿用 [BF-02 API](BF-02_Cch_Balanced_Csf_Load.md)，并逐q检查原geometry/full clock。按值副本共享存储，regrid/new publication/new color后需重建，不重贴身份。对象没有IO、主时钟、独立几何、材料更新或速度接受状态。
