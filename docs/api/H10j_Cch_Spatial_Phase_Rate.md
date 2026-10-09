# H10j / 同源空间相率 API

**当前状态（2026-10-05）：** 本owner与H10l2同一F/真实R保持；[H10o](../task/H10o_Bounded_Final_Native_Convection_Report.md)显式闭合空间算子和v2完整84项通过，下一[H08a](../task/H08a_SpaceTime_Material_Flux_Upgrade_TODO.md)。材料仍仅四硬约束，R/E_sp/G不要求零；本空间owner不改称同源时空积分，H11未启动。

2026-10-04，已实现；[任务与结果](../task/H10j_Spatial_Phase_Rate_Owner_Report.md)。生产入口仍为 `CchPhaseGeometry.h/.cpp`，属于原 `lamp_hydrodynamics_cch`，没有新库、target 或 Interface 状态。

```cpp
const auto reconstruction = CchAdvectorReconstruction::prepare(mesh, momentumGeometry);
const auto advector = reconstruction.evaluate(mesh, actualValues,
    phase.geometryIdentity(), phase.timeLayer());
const auto rates = CchPhaseRate::prepareSpatial(mesh, phase, advector);
```

调用者包含 `CchSpatialAdvector.h` 以使用完整 advector/identity 类型。原 phase header 可独立包含。

工厂核验同一 native layout、原 D_I instance、mesh、geometry revision 和完整 clock，保存原 immutable advector、momentumGeometry 和 FirstJet publication。没有接收裸速度重新拟合的入口，也没有 setters、时间推进、接受事务或失败 fallback。copies 共享同一 storage；advector identity 来自实际 evaluate，caller 声明不能代替 H14 的实际事务验证。

| 只读 port | 意义 |
| --- | --- |
| `geometry()` | 原同源 phase publication/metrics |
| `isSpatial()` | 恒为true，保留兼容查询，无legacy分支 |
| `spatialAdvector()/momentumGeometry()/advectorIdentity()` | 原 actual advector 与原几何 instance |
| `atomRates()` | reconstruction canonical atom 顺序的单份液/气相率 |
| `faceRates()/totalFaceRates()/cellOutwardRates()` | primal native 有向两相率、原 A_f*u_f、原 D_I incidence 求和 |
| `dualFaceRates()` | 同一 momentumGeometry.segments 顺序的实际 atom 求和 |
| `spatialReport()` | 查询/范围/event/refinement数量、原 tangent 缺口、误差估计与 closure |

所有相率为有符号 `long double`，单位物理测度/时间，不是 dt 积分；液/气分别积累，不从 total 减大相恢复小相。XY 为长度测度，RZ 为 2πr 真测度，axis 率严格零。

每个原 atom 从原 owning geometry view 查询一次 NativeQuarticTrace。ordinary 相范围积分同一实际 normal polynomial；正长度 zero 范围消费独立 minus/plus 相范围：`Q_l = ∫(a⁺χ_minus + a⁻χ_plus)dA`。同一零净面上的 actual polynomial 可变号，两相可分别非零。两侧事件不合并 ordinal 或近似全局 root 坐标；原 endpoint-relative offset/length 进入 public Core shifted 积分，opaque root chart 仅用于请求原 trace refinement。

正部分使用现 Core `integratePositiveUnitPolynomialRange` 和 1e−13 范围工作尺度；至多四轮 same-source refinement，每轮重读当前 event inventory，分别比较两相变化。ordinary/zero 的 terminal ranges 和邻轮变化是误差估计，**不是认证 bounds**；工厂报告实际估计，不修正通量以强制 closure。atom phase sum 按 1e−12 原 polynomial gross 尺度核验；有限性、精度预算或原几何查询失败明确拒绝。

selected availability 以 Hydro layout 为准。外部或 selected 零迹缺少正式流入相 port 时拒绝；额外 scalar 支撑、clip 内部迹和空间对称平均不能替代 boundary law。只有相关实际正部分严格零且误差估计严格零，缺少该侧才可接受。prepared 源错配、不同 D_I、改网格/regrid 或错误 FP 环境拒绝。

原 56 个 n4 selected2/3 未识别 tangent 保留，检查有限值、来源和共享；不称二维 affine 精确。full 和 n≥8 selected 的已识别族通过独立相积分；一般非多项式场、极细相或退化拓扑的普遍精度尚未证明。RZ physical mean 与 point 的 O(h²) 差异沿用 H10j0 的原声明。

H10l2已删除原evaluate和optional legacy存储，所有现行消费者只调用prepareSpatial。H08唯一Flux与H09共享同一dual F；现行缺陷/单位见[H09 API](H09b_Cch_Dual_Material.md)。材料四硬约束保持，quarter/严格material commuting/R=0仅诊断，Pure FirstJet独立；最终native与完整H10b仍未通过，当前入口见[接口索引](README.md)。

2026-10-05 [H08a2r7](../task/H08a2r7_Complete_Normal_Polynomial_Report.md)：空间 rate 的唯一法向多项式转换接受完整四次（RZ 密度五次），使用实际 `Detail::normalRatePolynomial` 和原 Core range/positive-range 积分。普通、共面变号正部和极短区间的高次消费已独立核验；当前实际速度仍为旧二次来源，时间积分仍为中点，不能由表示扩展声称时空相积分已通过。
