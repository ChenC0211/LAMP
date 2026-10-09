# H08 / 唯一材料积分与非守恒缺陷

## H08a 时空入口（2026-10-06，规定范围 PASS）

```cpp
CchSpaceTimeConfig quadrature;
quadrature.characteristic.maxStep = resolvedCharacteristicStep;
auto material = CchMaterialTransport::prepareSpaceTime(mesh,
    originalFirstJet, actualEndpointMotion, rhoLiquid, rhoGas,
    quadrature, prescribedOuterFirstJet, publicationConfig);
const auto& phaseIntegrals = material.phaseAtomIntegrals();
const auto& atomIntegrals = material.atomIntegrals();
const auto& primalIntegrals = material.faceIntegrals();
const auto& dualIntegrals = material.dualFaceIntegrals();
```

`prepareSpaceTime` 保留一个原始 Pure FirstJet publication 和一个实际 CCH endpoint source。所有积分时刻的几何直接从该原始 FirstJet 沿同一离散 motion/J 生成；不在三个独立快照间插值，也不从上一积分节点递推。mid/end 质量取真实候选几何。`characteristicSource()` 保留父区间 source 的两端 identity 和完整 clock；内部采样保留 `[begin,t]` 子区间语义，中点和结束保留完整父 clock。

每个时刻调用同源空间相率，再用 Gauss 3/5 正权对及全 atom 累积误差估计细分时间区间。默认相对预算 1e-7、初始四段、最大深度 24、最多 8192 次几何采样；工作量尺度来自真实 atom 测度与两端完整法向多项式系数，不用密度或 R 反推积分预算。`spaceTimeReport()` 给出采样数、最终分段数、最大深度及误差估计。该估计不是严格区间包络；尚未解析的尖锐瞬态仍需加密研究。预算耗尽明确失败，不回退中点。

一次密度加权 `F_atom=rho_l*Q_l+rho_g*Q_g` 后，primal/dual 只按原 canonical atom 索引求和。时空模式的 `atomRates()/faceRates()/dualFaceRates()` 是同一积分除以 dt；`midpointRate()/midpointFlux()` **仅为即时中点诊断**，不能作为区间通量替代。`isSpaceTime()` 区分入口。旧 `prepare` 保留原中点公式及共享率引用，历史 H08 结果仍可比较。

若提供 boundary callback，其语义是整条物理外边界在每个采样时刻的给定 FirstJet，包含切向/出流点和角点。内部点始终由同一 motion pullback，域外所需的内部历史仍明确拒绝。无 callback 时沿用原入流命中规则。当前要求完整 native FirstJet 支持、固定网格和可消费的实际运动；这不是完整 H14/H15 流体时间推进或通用边界求解器。

步积分缺陷 `r_int=Δm+ΣsF` 由真实几何质量差加唯一材料积分测得，缺陷率 `R=r_int/dt`。H09 消费同一 dual 积分并保留真实 `E_sp/G/uR`、材料动能及 KEP/skew。无库存、bubble、R=0 或 quarter commuting 约束。24 项短检查、28 个时空窗口及原 H10b 回归通过，详见 [H08a 报告](../task/H08a_SpaceTime_Material_Flux_Report.md)；通过范围是积分、材料一致性与小余项收敛，不代表完整流体 App、全耦合二阶或效率。

## 原中点入口与历史验收

历史状态：2026-10-03，**H08 PASS**；[任务书](../task/H08_Material_Flux_Continuity_Defect_TODO.md)、[报告](../task/H08_Material_Flux_Continuity_Defect_Report.md)。沿用[完全非守恒材料路线](../task/H07_Pure_Gals_Hydro_Coupling_Design.md)；当时下一项为 H09：同源 dual 材料与动量一致性。本项不推进速度或提交材料状态。

**重设计边界：** H10k已把H08 spatial材料路径迁移至[唯一Flux owner](H10k_Cch_Spatial_Material_Flux.md)。下述三层pressure质量/相差/R/中点积分定义保留；新F逐canonical atom加权后primal/dual仅求和，sources增加actual advector stamp，带expected momentumGeometry重载核验同一instance。H10l2已删除旧face-constant分支，所有消费者使用同一spatial owner；isSpatial()恒为true，仅保留兼容查询。四材料硬约束、独立Pure FirstJet与无库存/bubble保持；quarter/精确交换/R=0仅诊断，不改称时空积分或完整PDE通过。

公开头：[CchMaterialTransport.h](../../include/Modules/Hydrodynamics/CchMaterialTransport.h)，唯一实现 owner 为正式 `lamp_hydrodynamics_cch`，复用[H07a](H07a_Cch_Phase_Geometry.md)相测度/相率与同一 H05 metrics/D_I。

```cpp
auto material = CchMaterialTransport::prepare(mesh,
    beginGeometry, midpointPhaseRate, endGeometry, rhoLiquid, rhoGas);
material.validateFor(mesh, actualGeometrySources, actualMidpointClock);
const auto& massFlux = material.faceRates();
const auto& integratedFlux = material.faceIntegrals();
const auto& continuity = material.cells();
```

## 来源与时间层

三层必须使用同一 immutable layout、同一保留 D_I 和实际 mesh。开始层为 previous Initial/End，step=mid.step−1、time=mid.begin；中点与结束有相同 `{step,begin,end}`，位置分别是 Midpoint/End，且中点时刻必须严格位于可表示区间内部。Initial 支持明确的 restart 时刻。候选几何各有自己的 identity；不能以同一 time() 代替完整时钟比较，也不能将不同布局按相同尺寸直接混用。

`sources()` 返回开始/中点/结束三个 geometry identity和spatial路径的actual advector stamp，`timeLayer()` 返回完整中点时钟。`validateFor` 要求调用者提供当前实际来源，检查三层身份、实际 mesh 和完整候选区间。prepared 对象共享只读 data，保留三层几何及中点相率的寿命；不借 mesh 指针、不拥有主时钟。坐标/拓扑变化后旧对象拒绝使用，须显式 transfer/rebuild。区间内 regrid 不在本接口支持范围。

本接口校验来源绑定，不能证明调用者提供的任意三个 publication 必由同一物理 characteristic 产生。Interface producer 和后续 App/corrector 负责真实候选生成与接受；短契约另有实际 FirstJet producer 接线检查。

## 材料量与单位

两相密度有限且严格正，值保留在同一对象中；无密度 floor，本任务不支持真空。每层几何质量为 `rho_l*V_l+rho_g*V_g`，其中容积来自同一层的 H07a 几何，不新建材料库存。

`faceRates()` 按原生 velocity DOF 保存唯一材料流量 `rho_l*q_l+rho_g*q_g`，正号沿 native 正坐标法向。两侧相率分别消费原 H07a 结果，不能从 `Q−q_l` 的大数相减抹掉小气相或重新插值 `rho_f*u_f`。`faceIntegrals()` 是 `dt*faceRates()`，即冻结方法 §4.2 允许的**中点时间求积近似**，并非精确时空相通量。XY 为单位深度质量/时间和质量；RZ 包含同源 `2*pi*r`。

每个 `CchMaterialCell` 按 pressure DOF 保存三层几何质量、`massChange`、`outwardIntegral`、`integratedDefect`、`defectRate`、`relativeDefect`：

\[
\Delta m=\rho_l(V_l^{end}-V_l^{begin})+\rho_g(V_g^{end}-V_g^{begin}),\quad
r_{int}=\Delta m+\sum_f s_{cf} I_{\rho,f},\quad R=r_{int}/\Delta t.
\]

先对相容积作差，避免直接相减两份大材料质量。外流符号只读消费同一 D_I；共享内部 native face 使用同一份材料积分，2:1 使用真实 fine subface。

相对缺陷为 `abs(r_int)/(abs(delta_m)+sum(abs(I_rho,f)))`，分母和分子均为零时返回0。有量纲缺陷和原始量同时保留；非零材料缺陷本身不触发守恒性拒绝，不修改端点质量或几何。

`report()` 保存三层总几何质量、总质量变化、选中域边界材料积分、signed/L1积分缺陷、最大绝对R/相对cell缺陷，以及总cell缺陷与质量变化+边界积分的代数收支误差。后者门槛为256*double epsilon，它只证明诊断与共享通量的代数一致性，不要求物理缺陷为零。

材料负值/零/NaN/Inf、来源不符、实际mesh变更、非有限计算、不能表示的正密度乘积/非零面材料积分和退化中点时刻明确拒绝，已有prepared对象保持不变。非法来源不能以重贴标签绕过；裸数组只读访问后，消费者仍负责持有同一prepared对象并核验实际候选。

## 验证与范围

309个检查场景、6,780个cell、16,444个face检查及28项拒绝；另有8组单面注入验证实际相邻cell的相反贡献。XY/RZ、uniform/2:1、选中半域、纯相/平面/圆球、2^-40小相和密度比1/1e3/1e6覆盖。非多项式几何材料mass/face测量末档阶4.740–5.414；光滑面材料流量的中点时间积分末档阶1.997–2.000。

实际 FirstJet 平面跨面步骤保留了非零局部缺陷。整面重合时面相率可能有时间跳跃，单次中点求积不能预先保证局部二阶；缺陷如实暴露，完整两相时间精度由 H14/IH 验收。H09 原范围及 H10l2 同一 F 迁移已 PASS，最终 native 对流算子见 H10o 报告；当前 H08a 时空入口的验收范围见上文。完整 App、通用边界与低黏 EHD 仍未完成。

复现：`tools/build.cmd -Target quantum_cch_material_transport_contract`，随后 `ctest --preset local -R "^cch_material_transport_contract$" --output-on-failure`（CMake在PATH时）。
