# I08：唯一 phase-face flux

入口为 `Modules/Interface/HermiteGalsPhaseFlux.h`，namespace `Quantum::Interface`。本接口只准备一次 midpoint FV 库存更新及其唯一通量结果；不接受耦合步骤，不发布新的 scalar。

```cpp
HermiteGalsFluxVelocity velocity(mesh, midpointState, faceVelocityPolynomials);
const auto flux = advanceHermiteGalsPhaseVolumes(mesh, previousState, midpointState, velocity, config);
// flux.faces() 是共享通量；flux.phaseVolumes() 是下一 End 层的 M。
// I10/I11 将 M 与同 revision 的 first-jet candidate 闭合并发布。
```

## 输入与来源

`previousState` 提供真正需要更新的独立库存 M；可以是 Initial/restart 或前一步 End，不能是 Midpoint。`midpointState` 提供本次唯一 Q4 scalar 和 interval。要求 `mid.step=old.step+1`、`mid.begin=old.time()`，全区间 `dt=mid.end-mid.begin`；velocity identity/time 必须等于该 midpoint。两个 State 的 measure、active cells 和库存覆盖一致；每个库存 cell 都有 scalar 支持。midpoint 的 M 不用于替换 old M。

`HermiteGalsFluxVelocity` 深拷贝局部 normal velocity coefficients，并保留 midpoint handle。每个 actual native face 一条，按 `VEdge index` 后接 `HEdge index` 排序，只包含 active faces；漏项、重复、多余和次序错误均拒绝。half 0 是 begin→midpoint，half 1 是 midpoint→end；水平面沿 +x 参数化、法向为 +y，竖直面沿 +y 参数化、法向为 +x。

每半面的系数为 `u_n(t)=sum(a_k t^k), k=0..10`。它是由速度 owner 提供的局部只读 reconstruction，既不是新的 Hydro DOF，也不是全局 C1 高阶速度状态。XY 长度、RZ `2*pi*r` metric 由模块统一加入，调用方不能再乘一次 metric。RZ 含轴完整域及偶对称 scalar 的限制沿用已有 Snapshot。

内部面必须声明 `Interior`；active 域边界须显式为 `Open` 或 `Impermeable`。后者要求 normal polynomial 恒零，零净积分但有正反流不够。r=0 轴必须 Impermeable 且物理径向 normal velocity 为零，不利用面积为零隐藏不规则速度。Open 的流入 phase 来自本次显式提供的 midpoint scalar boundary trace；此处不生成物理边界历史。

## 唯一通量和积分

使用同一 source 的 `MappedFacePartitions`，正反向分别积分：

`Q+ = ∫max(u_n,0)dA`，`Q- = ∫max(-u_n,0)dA`，
`Q_l+ = ∫H(-phi)max(u_n,0)dA`，`Q_l- = ∫H(-phi)max(-u_n,0)dA`。

唯一 entry 保留四个非负 gross rates，`totalRate=Q+-Q-`，`liquidRate=Q_l+-Q_l-`。其中 `liquidVolume=dt*liquidRate` 只生成一次；negative owner 减它，positive owner 加同一个值。外边界贡献按相同方向单独累加。不能用 `liquidRate/totalRate` 代替四个方向约束：counterflow 可以净值零而 gross outflow 非零。

Core 新 `integratePositiveUnitPolynomialRange` 将同一 polynomial 通过 exact expansion 平移/缩放到 range 的局部 [0,1]，再调用原有 positive integrator。原 `offset` 和 `length` 分开保留，即使两个浮点全局端点重合也能积分。不对原 scalar/velocity 再采样拟合，不删除小系数。超过 exact expansion 可表示范围仍明确拒绝。

经证书确认整半面都是液相时，直接共享该半面的完整正反向积分，保证 pure identity。部分湿润面按所有 negative ranges 累加。I12r5 为 native Q4 增加孤立内部精确重根的相测度支持：在原隔离不足时，以 exact-expansion 的零余式证明公因子并核验无重因子商；原 bounded isolator 只隔离所有不同根，相区间符号仍从原 scalar 的精确样本取得。偶重根不翻相，奇重根翻相，近重根不被合并；根包围细化仍指向同一原 scalar 的零集。

这不改变 phase flux 的空间/时间公式，也不赋予相切点普通负相变分。含此类重根认证 chart 的 `nativeQuarticNegativePhaseVariation` 明确拒绝，不用小分母产生伪导数。旧 cubic/uniform 查询合同保持；算术范围或隔离不能认证时仍拒绝。positive-length coincident zero trace 仍不支持，不人为赋予中点相分数。I12r6 将正则孤立接触接入完整 native graph discovery：同一原标量的相位及法向导数包围证书决定单元内 0/1/2 条曲线关联，负相片不跨零测度点连接。此项能力仍要求原有 regular graph 认证，不能放行真正奇异零点。

## 有界性与舍入修复

I12r1 增加空间子集积分的一致性处理。Core 对 full/liquid ranges 独立计算正部积分时，其误差范围内可能出现 `Q_l^±>Q^±`，尤其在速度零点附近。仅在发生这种超界时，使用同一 polynomial 另积分互补气相 `G`，将液相估计 `L` 调整为 `Q^±*L/(L+G)`；`G=0` 时直接使用固定的 `Q^±`。这不改由速度唯一决定的 full flux，也不改变 midpoint 时间取样。

初始超量须同时落在 liquid/full 积分误差估计之和及 `64*eps(double)*(Q_half^++Q_half^-)` 内；完整调整须同时受 liquid/full/gas 误差估计之和、相同 gross 舍入预算约束，且 liquid 原误差加调整不得超过既有 half-face absolute integration budget。最终每面的单方向调整乘 dt 还须不超过 `64*eps*min(adjacent capacities)`，边界使用唯一内部 cell capacity。所有检查失败均拒绝。调整量计入误差报告；两个非负 half-face 的 full/liquid 分量按同一顺序相加，保持 subset inequality。

要求 `0<=Q_l+<=Q+`、`0<=Q_l-<=Q-`；总通量的 `|dt*sum(s Q)|/capacity` 通过 divergence gate。gross-outflow CFL 是实际离开该 cell 的两个方向分量之和乘 dt，再除 capacity，必须 `<=1`。超出这些约束拒绝。

先用 raw high-order flux 组装全部 M。若全部严格有界，通量不改。非舍入级库存越界拒绝；仅 `64*eps(double)*capacity` 内的小 defect 允许通过内部共享 face 调整：

- underfill 缩减该 cell 的内部液相流出；overfill 缩减流入；同面两侧共享相同系数。
- 只缩减对应 liquid gross direction，不改 total flux、prescribed boundary flux 或 cell M。
- 每个被改面的累计单方向修正也受 `64*eps*min(capacity_left,capacity_right)` 限制。
- 每次重新组装共享流量后检查严格 bounds；缺内部修复容量或 128 次 work cap 耗尽均拒绝。

这只是机器舍入修复，不是可处理有限幅值 overshoot 的 limiter，也没有 upwind/FCT 后备路径。上层须为宏观拒绝重新准备满足约束的 dt/中点候选；本接口不静默选择新算法或改时钟。

## Result / Report

结果持有 immutable old/midpoint States，保持旧库存及 scalar 生命周期；`geometry()` 始终是原 midpoint 的同一地址/identity。`phaseVolumes()` 是有限体积所得 End 库存，结果不会复制或修改 first jets、求 bubble 或绑定新 geometry。

`faces()` 的 `integratedLiquidPositiveRate/NegativeRate` 保存独立空间积分原值；`rawLiquidPositiveRate/NegativeRate` 保存空间一致性处理后、库存微修复前的值；`liquidPositiveRate/NegativeRate` 保存最终值。总通量、dt 乘积和积分误差估计同时保留，只有最终通量供后续 Hydro 使用，其余只作诊断。

Report 包含 before/after/boundary/closure、gross CFL、divergence、counterflow faces、库存微修复次数/规模等。I12r1 增加 `subsetReconciliationCount`（发生一致性调整的 half-face/direction 数）、`directionalIdentityCount`（其中 gas 积分及其误差均为零的数）和 `maxRelativeSubsetAdjustment`（最大调整除以所属 half-face total gross rate）。这些量与库存修复的 `repairedFaces/maxRelativeFluxCorrection` 分开；Report 不决定分支。

具名容限默认：`relativeIntegrationTolerance=1e-13`、`relativeDivergenceTolerance=1e-12`、`relativeConservationTolerance=1e-12`。正反向 primitive 的 absolute budget 为该半面 density bound × 相对积分容限，range 再按长度分配。global closure 除以总 capacity；微修复预算独立且不可放宽。正多项式误差估计不包含 scalar 表示、几何根位置或物理 metric 舍入误差，因此仍须独立解析/MMS 验证。

## 当前边界

actual AMR subfaces 已按同一 entry 写入两侧；I09 已实现 coarse full-face 与 fine subflux 聚合，I10/I11 已实现局部 bubble 闭合及同源发布。I08 的解析 midpoint、face-cut 收敛和静止检查不等于 I12 的旋转/vortex/AMR-crossing 全步 benchmark。I12r1 修复空间 gross subset 缺陷；I12r5 处理 native Q4 面内精确相切的相区间与通量，I12r6 接入正则孤立接触的完整拓扑与发布。任意步内跨面的时间通量误差仍存在。尚无 Gate I、CCH 接线、动态 regrid 或低黏/无黏稳定性完成声明。

见 [I08 TODO](../task/I08_Unique_Phase_Flux_TODO.md)、[I08 报告](../task/I08_Unique_Phase_Flux_Report.md)、[I12r1 修复报告](../task/I12r1_Gross_Phase_Flux_Bounds_Report.md)、[I12r5 面相切报告](../task/I12r5_Native_Q4_Tangent_Phase_Report.md) 与 [I12r6 孤立接触报告](../task/I12r6_Native_Isolated_Contact_Report.md)。
