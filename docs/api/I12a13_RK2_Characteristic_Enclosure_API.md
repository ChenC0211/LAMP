# I12a13 / 实际 RK2 特征区间 API

本项消费 [I12a12 polynomial source](I12a12_Polynomial_Velocity_Interval_Source_API.md)，按[同源区间设计](../task/I12_RK2_Characteristic_Enclosure_Design.md)包围旧 Core RK2 的完整单步、无 hit 分支。查询不发布 geometry，不更新独立守恒 M，不生成 face transfer。

`PolynomialVelocitySlab::encloseInteriorRk2(domain, lower, upper, arrivalBegin, arrivalEnd)` 返回 `PolynomialRk2Bounds`。Core 接受退化空间 box 和闭时间区间；`slab.timeBegin() <= arrivalBegin <= arrivalEnd <= slab.timeEnd()`。回溯时长始终从 slab 的旧时刻起算，不能把 arrival 子区间起点当成新的 old state。

| 返回字段 | 含义 |
| --- | --- |
| `arrivalLower/Upper` | 原输入空间 box |
| `midpointLower/Upper` | 原运算顺序下的 RK2 中点包围 |
| `departureLower/Upper` | 同一旧时刻上的离散 B_h 包围 |
| `jacobian[4]` | 原离散链式 Jacobian 的 xx、xy、yx、yy 范围 |
| `maximumL1Travel` | 各实际采样分支保守速度 L1 × 完整回溯时长上界，须严格小于 h/2；不是误差估计 |

先证明输入、midpoint 与 departure 在物理域内，满足原 sample 半 cell 限制、RZ 轴的精确零径向速度/切向导数，且不触发初始或新 wall contact，之后才返回。域/尺度检查复用原 tracer 构造规则。opaque `activeSegment` 拒绝；不靠有限轨迹样本给任意覆盖函数作证明。返回范围可能保守，无法认证时抛出异常；调用者可在后续任务中细分，不能把失败当零通量或裁剪回溯点。

普通正时长按实际 `(u_arrival(H)+u_old(H))*.5` 传播，包含浮点时间插值的实际分支；半步时间下溢到零时额外覆盖 arrival 分支。区间含旧时刻时也包住 B=X、J=I；纯旧时刻查询直接返回恒等结果。RZ 精确轴线可成功，包含轴邻域的一般 box 可能因过宽的负 r 范围而拒绝。

`encloseHermiteGalsCharacteristic(mesh, interval, lower, upper)` 返回不可默认构造的 `HermiteGalsCharacteristicEnclosure`。`bounds()` 与 `interval()` 均为只读引用，`validateMesh(mesh)` 使用原生几何/拓扑身份。Interface 沿用 interval 的正时间宽度约束，Core 另支持单点时间查询。

结果保留 interval、old state / 单一 snapshot、owned polynomial velocity、characteristic 的原 owner；复制共享这些来源，移动后的源对象拒绝访问。独立创建但取值相同的速度 source 仍是不同来源。polynomial 必须显式创建；opaque 或已有 constant source 不因样本相同而自动获得 polynomial 认证。旧 scalar support 必须覆盖所有 native leaves。`NativeMeshIdentity` 比较几何/拓扑数据，允许数据完全相同的另一 mesh 对象，拒绝实际网格变化。

调用层可读错误包括无 owned polynomial、旧标量支持不全、过期 mesh、moved-from；Core 错误包装中保留 previous geometry identity。非法 box/time、浮点 rounding/FTZ、域、轴、CFL 或尚未证明的边界分支均明确拒绝。查询成功或失败都不修改 old geometry/M。

范围证据来自向外舍入算术、共用 RK2 表达及充分分支条件；样本包含性测试用于验证实现。此 API 不证明连续 ODE 轨迹、边界历史、相事件完备、相通量积分、一般速度场的库存有界性或 I12 时间阶。
