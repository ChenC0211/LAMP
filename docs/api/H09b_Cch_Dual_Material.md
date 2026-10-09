# H09 / 唯一共享材料通量与非守恒动量一致性

2026-10-04，当前接口由 [H10l2](../task/H10l2_Same_F_Dual_Material_Retry_Report.md)迁移。原H09b PASS及H10l FAIL见各自历史报告；旧quarter延拓、SPD材料交换与legacy相率入口已退役。公开头为 [CchDualMaterial.h](../../include/Modules/Hydrodynamics/CchDualMaterial.h)，唯一实现归原 `lamp_hydrodynamics_cch`。

```cpp
auto dualMaterial = CchDualMaterial::prepare(mesh, material,
    actualBeginMass, actualMidpointMass, actualEndMass);
dualMaterial.validateFor(mesh, actualSources, actualMidpointClock);
const auto& Fdual = dualMaterial.faceRates();
const auto& Idual = dualMaterial.faceIntegrals();
const auto& continuity = dualMaterial.cells();
auto uRdt = dualMaterial.consistencyImpulse(mesh, actualMidpointVelocity,
    actualSources, actualMidpointClock);
auto materialEnergy = dualMaterial.materialKineticChange(mesh,
    actualMidpointVelocity, actualSources, actualMidpointClock);
```

## 唯一来源和存储

prepared对象保留H08唯一空间材料Flux的区间消费者，以及三个H09a真实质量对象；复制共享只读data。三层质量对应实际开始/中点/结束publication identity与完整time layer，须共用同一momentumGeometry instance、native layout、原D_I和两相密度。material的actual advector/geometry instance须与此geometry一致。来源校验不能代替H14实际候选生成/接受事务。

`faceRates()`直接返回 `material().dualFaceRates()` 的**同一vector引用**，顺序为原momentumGeometry的真实segments，正号沿正坐标法向。H09不重新积分、乘密度、按面积分配primal率或求quarter；相交积分和唯一密度加权已由[空间相率](H10j_Cch_Spatial_Phase_Rate.md)与[材料Flux](H10k_Cch_Spatial_Material_Flux.md)完成。内部两侧消费同一值并反号，AMR共同侧是实际子面之和。

`faceIntegrals()`直接消费 H08 的 `dualFaceIntegrals()`，H09 不再自行按 dt 重算。旧 `prepare` 仍给出具名中点近似；H08a 的 `prepareSpaceTime` 给出真正沿同一几何路径的时间积分，其 `faceRates()` 为同一积分的区间平均率。两种入口都用相同 incidence 与真实三层质量核算 R/uR；H08a 的规定积分、一致性与小余项收敛范围已 [PASS](../task/H08a_SpaceTime_Material_Flux_Report.md)。齐次不可穿透条件由调用者在同一actual B上显式绑定[typed law](H10l1_Cch_Advector_Boundary_Closure.md)；H09不凭零native均值推断墙面，也不后置压F。

`material()/beginMass()/midpointMass()/endMass()`返回原对象；`geometry()/layout()/sources()/timeLayer()/deltaTime()`返回同一来源。mesh不借指针；等效mesh与原mesh析构寿命受原identity契约保护。坐标或拓扑变化后旧对象拒绝，必须显式transfer/rebuild并延续物理clock，不能给旧数组重贴标签。

## 真实缺陷与诊断

`cells()`按native velocity DOF保存真实三层几何质量之差及共享率的实际incidence：

\[
\Delta m_d=\rho_l(V_{l,d}^{end}-V_{l,d}^{begin})+\rho_g(V_{g,d}^{end}-V_{g,d}^{begin}),
\quad q_d=\sum_s s_{ds}F_s,
\quad r_d=\Delta m_d+\Delta t q_d,\qquad R_d=r_d/\Delta t.
\]

相容积先作差，避免独立大质量相减。`massChange/outwardIntegral/integratedDefect/defectRate`分别对应Δm、dt·q、r、R；`relativeDefect=abs(r)/(abs(delta_m)+sum(abs(dt*F)))`，全零时严格0。非零R不拒绝、不填零，也不修复材料或界面。

原真实overlapFraction定义T，只用于具名诊断：

\[
E_{sp}=q_d-Tq_p,\quad G=\Delta m_d-T\Delta m_p,
\qquad r_d=Tr_p+G+\Delta t E_{sp}.
\]

| 字段 | 意义 / 单位 |
| --- | --- |
| `outwardRate` | 真实dual q / 质量·时间⁻¹ |
| `mappedPrimalOutwardRate` | Tq_p / 质量·时间⁻¹ |
| `spatialExchangeRate` | E_sp / 质量·时间⁻¹ |
| `mappedPrimalMassChange` | TΔm_p / 质量 |
| `geometryTransferChange` | G / 质量 |
| `mappedPrimalIntegratedDefect` | Tr_p / 质量 |

E_sp/G/R均可非零；不能把T映射质量冒称真实dual质量。`maxAbsSpatialExchange`与`maxRelativeSpatialExchange`报告真实E_sp，后者分母为dual gross加mapped primal gross。`maxDefectDecomposition`核验含dt·E_sp的身份，沿用256·double epsilon和实际gross尺度，零work要求精确零。原`maxQuarterExchange/maxDualExchange`字段删除。

`report().components`按两个速度方向记录三层总几何质量、总质量变化、边界材料积分、signed/L1缺陷、最大绝对R/相对缺陷及总收支身份。内部反号共享使总收支归结边界；检查是诊断自洽，不要求材料严格守恒。

## 动量与材料动能

`consistencyImpulse()`为同一中点速度返回 `u_mid*r_d`。无外力常速度U满足

\[
m_{end}U=m_{begin}U-U\,outwardIntegral+U r_d.
\]

`materialKineticChange()`按方向返回 `sum(0.5*u_mid²*r_d)`；该项来自同一R，生产公式在迁移中保持。它不是一般速度的端点KE差，完整KEP预算和H14时间离散仍需分别验收，不能宣称非守恒材料下无黏动能无条件严格守恒。

返回load数组是派生量，调用者须保留prepared owner并核验实际sources/full clock。速度端口验证相同layout和RZ轴速度为0；裸CchNativeValues不自带物理clock。来源/instance/D_I/密度/clock错配、regrid、非有限计算、不可表示的非零积分或一致性项明确拒绝，无fallback。

本项只接通唯一材料率与R/uR；[H10o](../task/H10o_Bounded_Final_Native_Convection_Report.md)已完成显式闭合最终native空间/KE算子及v2完整研究，H08a 时空入口也已在规定范围通过，H11 局部受控耗散也已在[规定范围通过](../task/H11_Controlled_Dissipation_Report.md)，下一单任务为 H12，本轮未启动；完整 App、全耦合二阶、一般切界面 PDE 和动态效率未验收，具体范围见报告。Pure FirstJet独立，无库存/bubble。复现短契约：`tools/build.cmd -Target quantum_cch_dual_material_contract`，再用标准CTest运行同名测试。
