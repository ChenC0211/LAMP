# BF-01：只读正则化 CSF color

`CchRegularizedColor` 位于 `lamp_hydrodynamics_cch`，由同一 `CchPhaseGeometry` 保留的 Pure FirstJet publication 派生力场。它保留不可变 phase/config、native pressure 顺序的单元 color 和求积报告；按值复制共享只读存储。该对象不推进 color，不修改 `phi`、sharp 相体积、材料质量、`F_rho/R/uR` 或 Interface accepted state。BF-01规定范围已PASS，证据见 [报告](../task/BF-01_Report.md)；此接口尚不构造CSF载荷，不代表后续动力学验收通过。

```cpp
CchRegularizedColorConfig cfg;
cfg.epsilon = physicalBandHalfWidth; // 显式、统一的物理长度
const auto color = CchRegularizedColor::build(mesh, phase, cfg);
color.validateFor(mesh, phase);
const auto& c = color.colors();      // c[p]，不是 native cell index
const auto& audit = color.report();
```

## 定义与核函数

依照 [BF-CSF V1 §3](../roadmap/LAMP_BF_CSF_Method_V1.md)，原 Q4 jet 给出

\[
d_h=\frac{\phi_h}{|\nabla\phi_h|},\qquad
c_\epsilon(x)=H_\epsilon[-d_h(x)],\qquad
c_K=\frac1{V_K}\int_Kc_\epsilon(x)\,dV.
\]

`phi<0` 为 liquid。`phi/|grad phi|` 是具有长度量纲的局部距离近似，通常不是真正 signed distance；此对象不 redistance 或重建另一条界面。正的标量输入缩放在可表示、可辨识条件下应消去；输入溢出或退化仍拒绝。

`heaviside(s,e)` 在 `s<=-e` 为 0、`s>=e` 为 1，内部为 `1/2+s/(2e)+sin(pi*s/e)/(2*pi)`。端点附近采用等价正值级数减少相消，不 clipping。`delta(s,e)` 是它对长度参数的导数：内部 `cos²((pi/2)*(s/e))/e`，带外及端点精确为 0。计算先取 `s/e`，避免巨大但有限的 `e` 在 `2*e` 中溢出；极小 `e` 导致 delta 不可表示时明确拒绝。两个函数均要求有限 `s` 和有限正 `e`。

## 配置与拒绝条件

| 字段 | 当前默认值/条件 | 语义 |
|:---|:---|:---|
| `epsilon` | 默认 0，必须显式给出有限正值 | 所有选中 pressure cells 共用的物理半带宽；不会自动扩大 |
| `relativeTolerance` | `1e-9`，须在 `(0,1)` | 每 cell 的积分误差估计预算为 `relativeTolerance*V_K` |
| `relativeGradientResolution` | `256*machine_epsilon`，须在 `(0,1)` | 带内要求 `|grad phi| > resolution*gradientScale` |
| `maxDepth` | 20，允许 `[0,40]` | 每次二分的最大深度，耗尽则拒绝 |
| `maxQueriesPerCell` | 262144，配置允许至少 320 | 原 Q4 点 jet 与 affine 包络查询的合计上限；两项分别报告 |

`gradientScale` 为原 cell 的 `max(|phiLower|,|phiUpper|)/min(hx,hy)`；它来自原 Q4 包络，须有限且正。带内梯度非有限、低于局部相对分辨阈值或原标量恒零时抛错。带外有确定非零符号的常值场可以直接取得纯相，不要求人为梯度下限。

非 pure 的每个最终接受 panel 必须先证明整个闭矩形上的 `gradientLower > up(relativeGradientResolution*gradientScale)`；`up` 表示向正无穷的外向取整。未能证明时强制细分，几何认证优先于 GL 误差，即使采样误差为零也不能接受。无法在深度或工作上限内完成证明会明确拒绝；这是一项保守的全 panel 条件，可能拒绝带外梯度弱而尚不能被整体认证 pure 的输入，不能将认证失败等同于已证明原几何必然退化。

当前原生分辨条件为 `epsilon >= max(hx,hy)/2`。sharp cut cell 在求积前检查；其他 cell 若宽度不足而包络仍允许力带存在，必须继续细分直到子矩形被证明 pure，或被证明全在带内后拒绝窄带配置。不能凭 GL 没采到窄带放行，无法结束认证时仍按明确工作上限拒绝。这个条件不证明动力学已充分分辨；接口不触发 regrid。错误报告保留 cell index，构造失败不返回部分 color。

上限的配置合法性与实际工作量分开：320 是允许的最小值，但一个无需细分的 mixed root panel 已需 10 次初始化包络、9 次 panel 包络及 `8²+16²` 次点 jet，共 339 次，因此可能在合法小上限下触发工作量拒绝。纯相仍需要包络查询，只能免除点 jet 查询。

## 真实体积、正权求积和误差报告

XY 使用 `dx*dy`，RZ 使用 `2*pi*r*dr*dz`；每个实际 leaf cell（含 2:1 AMR 与 selected 域）独立积分。RZ 子矩形容量为 `pi*(r1²-r0²)*(z1-z0)` 的稳定等价乘积，径向权重只乘一次，包含 `r=0` 的 cell 无 `1/r` 查询。

原 Q4 包络能够判定整个子矩形位于力带外时，直接返回该子矩形真实容量乘精确 0/1。其余先满足上述整 panel 几何认证，再使用张量 GL8/GL16 正权平均，以 companion liquid/gas 权重和真实子矩形容量归一，故常值与 `0<=c<=1` 不依赖 clipping 或材料体积修正。未认证 panel 优先细分，认证完成后选择当前误差最大的活动 panel，均沿物理较长方向二分。即使 GL 点全为 0 或 1，仍以包络给出的 color 跨度补充误差估计，防止只凭采样纯相跳过细带。

`absoluteErrorEstimate` 主要由 GL8/GL16 差值累加，属于自适应求积**误差估计**。它不是所有任意 Q4 输入的严格积分误差上界；通过内部预算仍需解析平均、独立积分和加密核验。`transitionVolumeEstimate` 同样是采样估计；`transitionVolumeLower/Upper` 来自最终活动 panel 的“全带内/可能带内”包络分类，可较宽，不能将二者混称同等精度的数值积分。

原多项式梯度包络复用 `nativeQuarticAffineEnclosure`。对每个坐标方向，使用数学恒等式

\[
P'(x)=\frac{8[P(x+a)-P(x-a)]-[P(x+2a)-P(x-2a)]}{12a},
\quad a=h_{axis}/4,
\]

它对该方向次数不超过 4 的 Q4 多项式精确成立，不是近似有限差分或 refit。平移坐标用 TwoSum 保留舍入偏移，区间组合向外取整并保留原 affine remainder。平移可能超出 cell，查询的是**同一多项式的数学延拓**；并未读取邻接 cell、physical ghost 或墙外状态。这样的延拓只用于包络，实际体积分始终在原 cell 内。

认证所用范数上下界按分量最大值缩放，再以外向取整的商、平方、和、`sqrt` 与回缩乘积形成，避免直接平方的大数溢出。Heaviside 区间包络使用 `H(-1+a)=a/2*(1-sinc(pi*a))` 的递减交错级数：第 16/17 项外向区间和分别给出下/上界，正半轴由对称关系获得，避免将未经定向的 libm sine 值当成认证界。数学包络可与已知 `[0,1]` 相交；实际 point/cell color 数值始终不 clipping。这些几何包络保证与 GL 积分误差估计的性质不同。

`cellReports()[p]` 报告点/包络查询数、活动 panels、深度、误差估计、带体积估计/包络、输入标度和整 cell 纯相认证。`minimumResolvedGradientRatio` 是实际带内点样本的最小 `g/scale`；`minimumCertifiedGradientRatio` 是最终非 pure panels 的最小 `gradientLower/scale`，反映全 panel 认证。二者不能互换：没有对应点样本或非 pure panels 时，相应字段保持正无穷，不能作为实测最小值。总 `report()` 分别累计两类查询，还给出真实域体积、纯相 cell 数，以及 `sum(V_K*|c_K-alpha_sharp,K|)` 和最大 color 差；这些差异是力 color 与 sharp 材料账本的诊断，不用于修正材料。

## 身份、时钟与消费者

`build` 验证原 phase 的 mesh/layout、geometry identity、完整 `{step,begin,end,position}` 时层及 FirstJet 语义。`geometryIdentity()`、`timeLayer()` 和 `phaseGeometry()` 返回保留来源。`validateFor(mesh, expectedPhase)` 同时验证两份 phase、相同 layout identity 和同一 `D_I` matrix owner；只相同尺寸或重新装配的等值 divergence 不足。另一重载 `validateFor(mesh, identity, clock)` 仅提供 phase 原有的 mesh/source/full-clock 验证，不附带外部 `D_I` owner 比较。

regrid、几何 publication 或完整时层改变后应从实际新 phase 重建 color，禁止重贴身份或复用旧数组。selected 域依循原 pressure layout；可保留更广的 FirstJet 标量支撑，但只输出选中 pressure cells 的 color。

`pointColor(p,x)` 是原 native cell 内的只读点查询，用于诊断，不能替代 `colors()[p]` 的体积平均。它验证 pressure index 并沿用有限/梯度条件；接口不传 mesh，因此消费者应先调用对应 `validateFor`，并遵守原 native-cell 查询范围。后续 BF-02 需以现有同源 `C=-D_I^T` 消费此 color；BF-01 不输出压力、动量或 CSF 力。

任务状态和实测证据以 [BF-01 TODO](../task/BF-01_TODO.md) 及后续报告为准，不能由本 API 文档推断 PASS。
