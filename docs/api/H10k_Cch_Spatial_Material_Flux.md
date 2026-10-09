# H10k / 唯一空间材料通量与区间缺陷 API

**当前状态（2026-10-05）：** 本owner与H10l2同一F/真实R保持；[H10o](../task/H10o_Bounded_Final_Native_Convection_Report.md)显式闭合空间算子和v2完整84项通过，下一[H08a](../task/H08a_SpaceTime_Material_Flux_Upgrade_TODO.md)。材料仍仅四硬约束，R/E_sp/G不要求零；本空间owner不改称同源时空积分，H11未启动。

2026-10-04，已实现；[实际验收](../task/H10k_Material_Spatial_Owner_Migration_Report.md)。唯一生产入口为现 `CchMaterialTransport.h/.cpp`，属于原 `lamp_hydrodynamics_cch`。没有新库、target 或材料状态。

```cpp
const auto phaseRate = CchPhaseRate::prepareSpatial(mesh, midpointPhase, actualAdvector);
const auto flux = CchMaterialFlux::prepare(mesh, phaseRate, rhoLiquid, rhoGas);
flux.validateFor(mesh, actualMomentumGeometry, midpointPhase.geometryIdentity(),
    actualAdvector.identity(), midpointPhase.timeLayer());
const auto material = CchMaterialTransport::prepare(mesh, beginPhase, flux, endPhase);
material.validateFor(mesh, actualMomentumGeometry, material.sources(), material.timeLayer());
```

`actualAdvector.identity()` 指实际 reconstruction evaluate 返回的 identity；调用者也可直接消费 `phaseRate.advectorIdentity()`。实际 advector 的工厂与 port 见[H10j0](H10j0_Cch_Spatial_Advector.md)。

`CchMaterialFlux` 是唯一瞬时 F owner，保留原 immutable spatial phase、actual advector、原 momentumGeometry/layout/D_I、geometry identity 和完整 clock。对每个原 canonical atom 只计算一次 `F=rho_l*Q_l+rho_g*Q_g`；primal/dual 视图只按原 reconstruction 索引补偿求和，不再次查询几何、拟合速度或乘密度。液/气乘积分别有限且非零可表示；有符号相分量的合法相消允许为零。密度有限且严格正，无 floor/真空。

| Flux 只读 port | 来源与单位 |
| --- | --- |
| `phaseRate()/momentumGeometry()/layout()` | 原 prepared 对象及其实际 instance |
| `geometryIdentity()/advectorIdentity()/timeLayer()` | 原瞬时来源与完整时钟，不能由 caller 改标签 |
| `liquidDensity()/gasDensity()` | 同一 F 的两密度 |
| `atomRates()` | canonical atom 顺序，质量/时间 |
| `faceRates()` | 原 native velocity DOF 顺序，同一 atom 求和 |
| `dualFaceRates()` | 原 momentumGeometry.segments 顺序，同一 atom 求和 |

所有率为 signed `long double`，正号沿正坐标法向；RZ 含真实 2πr。Flux 可消费合法 Initial/Midpoint/End，但不包含 dt、R、库存或接受事务。copies 保留同一 storage，不能默认构造或写入数组。

`CchMaterialTransport` 是 Flux 的区间消费者。新 typed prepare 保留传入 Flux 的**同一 atom/primal/dual 数组**，只生成原 `dt*F`、三层 pressure 几何质量和 R。原 prepare 接收 spatial phaseRate 时仅编排同一 Flux 工厂；不会另乘密度。`midpointFlux()/atomRates()/dualFaceRates()/momentumGeometry()` 暴露原 storage。

`CchMaterialSources` 为三层 geometry identities 加实际 advector stamp。spatial `validateFor` 精确比较 stamp、完整 clock、原 layout/D_I/mesh；带 expected momentumGeometry 的重载进一步核验实际 instance。同坐标 fresh 对象或零/伪造 stamp 不代表同源。来源校验不是 H14 物理候选生成/接受证明。

区间时钟、相体积先作差、pressure R 与全域 balance 沿用 [H08](H08_Cch_Material_Transport.md)。材料积分仍是**中点时间近似**；非零 R 不拒绝、不填零、不修改几何。错来源/实例/时钟、regrid、非有限值及不可表示相乘积/非零积分明确拒绝；无 fallback。

H10l2已迁移最后dual/KEP两份契约，删除CchPhaseRate::evaluate及Transport legacy分支。CchDualMaterial::faceRates()直接返回本owner的dualFaceRates()原引用，quarter数值延拓已退役；isSpatial()恒为true，仅保留兼容查询。真实R/E_sp/G及来源定义见[H09 API](H09b_Cch_Dual_Material.md)，最终native空间精度与material KE仍待联合收口，当前入口见[接口索引](README.md)。

独立真实 CSV 覆盖64布局、128 phase factories、384材料准备和三密度；68.6万 atom/94.7万 view/10.8万 cell记录通过。原56个粗 selected 未识别 tangent 不称二维 affine 精确；其他来源/共享/有限性仍核验。极小流速触发 Core 精确乘积尾部下溢时明确拒绝；有限实际测试不能证明任意退化几何精度。规定几何运动与 Hydro MMS 不相同，当前研究是空间通量/诊断身份验证，不能说明 R 已降为小余项或耦合 PDE 通过。Pure FirstJet 独立，无库存/bubble；H10/H11/H08a 顺序保持。
