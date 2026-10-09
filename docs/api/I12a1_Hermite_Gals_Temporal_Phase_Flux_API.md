# I12a1 / Hermite–GALS 时间区间相体积契约

这是数据/来源契约，不是时间积分器，也不产生新的 accepted state。头文件为 `include/Modules/Interface/HermiteGalsTemporalPhaseFlux.h`，位于既有 `Quantum::Interface` V2 主线；旧 `InterfaceTemporalPhaseFlux.*`、I08 midpoint 数值算子和 I09 消费类型均保留。

## 来源与实际身份

`HermiteGalsTemporalVelocitySource` 持有真实的 previous/arrival `HermiteGalsVelocityView`。构造时传入旧 state、arrival geometry identity/time 和两端 Query，经原 view 构造器校验空 callback、身份和时间。只保存输入，不执行查询。previous 必须为 Initial/End，arrival 的真实采样时间晚于 previous，且其 begin 与旧物理时间一致；arrival 可以是现有 End 或 Midpoint layer，实际覆盖终点始终用 `time()`。

`HermiteGalsTemporalCharacteristicSource` 持有现行 coherent explicit midpoint RK2 的 `LevelSetBoundaryTrace::Domain`，含实际 coverage callback、cellSize、XY/RZ 模式，并保存 `NativeMeshIdentity`。它不替换 tracer、不调用 tracer、不提供另一种特征方法。其 domain 可与 velocity source 的原 endpoint views 直接交给既有 `traceHermiteGalsWithJacobian`。

两个来源都采用现有模块的 `shared_ptr<const Data/Revision>` owner 模式。复制共享真实来源，`sameSource()` 比较不可变 owner 的身份；独立构造视为不同来源，即使碰巧具有相同 geometry/time 或相同函数值。没有数值计数器、地址哈希、调用方自由填写的 token，亦不把 geometry ID 当速度 ID。公共只读 accessor 可追溯到实际持有的 view/domain；不存在与对象脱离的裸来源编号。

这是用户补充“细节与现有程序冲突时以现有程序路线为主”的适配。既有模式可参照 `InterfaceConservativePhase` 的共享 revision，以及 `InterfaceNormalFluxView::scaled` 对同一空间 source 的共享身份。新类型没有扩展这些 legacy 算法。

与现有 velocity/tracer 一样，任意 callback 的 capture 必须由 provider 保持有效、不可变且可重复。本契约拥有 callable 的副本，但无法证明 callback 的外部可变引用不会改变；这不是事件覆盖或函数等价性证书。来源 identity 只用于进程内所有权，不是 restart 序列化 ID。

## 时间区间

`HermiteGalsPhaseFluxIntervalSource(mesh, previousAccepted, velocity, characteristic, begin, end)` 保存同一旧 `HermiteGalsState` handle，绝不重建 geometry。App 决定哪个 Initial/End handle 是 accepted；I01 本身没有全局 accepted 标记，本契约不冒充耦合 commit 检查。Midpoint state 不能充当旧 accepted 输入。

构造要求：

- begin/end 有限、end>begin、dt=end-begin 有限且可表示。
- previous velocity 的 geometry identity 与完整 time-layer metadata 对应旧 state。
- 区间位于旧时刻至 arrival 实际采样时间的范围内；子区间可以从旧时刻之后开始，但仍共享同一旧 source。
- characteristic 与旧 state 的网格及测度相容。

`timeBegin()`、`timeEnd()`、`dt()` 暴露真实区间，不创建 `Position::Midpoint` 或新的伪 geometry identity。`sameInterval()` 比较旧 geometry identity、两个来源 owner 和区间端点。相同来源的 half/full 区间保留相同旧 authority，但 extent 不同。

mesh identity 沿用当前精确 geometry/connectivity 内容校验，没有另造 mesh revision counter。相同内容的 mesh clone 可通过；坐标/拓扑改变后拒绝，恢复完全相同内容后重新相容。`meshIdentity()` 提供只读旧身份，`validateMesh(mesh)` 执行实际校验。构造/组合和 `face(mesh,row)` 都检查 mesh；只读历史 volume/metadata 查询不要求外部 mesh 继续存活。

## 积分体积及派生 rate

`HermiteGalsIntegratedFaceVolumes` 输入包含四个长精度非负积分量：

```text
integratedTotalPositiveVolume
integratedTotalNegativeVolume
integratedLiquidPositiveVolume
integratedLiquidNegativeVolume
```

`HermiteGalsIntegratedFaceTransfer` 复制输入并只读持有；不存在 rate 存储。它由真实 native face 查询得到相邻 cell、方向和既有测度。`PositiveX` 只用于竖直面，`PositiveY` 只用于水平面；不可用任意反向标记将另一份 transfer 冒充同一 native owner。单元外向符号由后续 consumer 使用相邻身份确定。

signed volume 是相应 positive 减 negative；average rate 现场除以 interval.dt()，调用方不能另传 dt。XY 沿用单位深度面积约定；RZ 沿用真实旋转体积。RZ 轴等零测度面只接受零 total/liquid volume。

所有分量必须有限且 >=0，液相分量逐方向 <= total；`-0` 按数值零接受。没有负值宽容、clipping 或舍入修复。这里校验的是调用方给定的代数量，不负责求积误差预算、CFL、库存有界性或 bubble 可达性。派生 rate 或合并量发生溢出时拒绝，不改变原积分量。

## 受检组合

`left.appended(mesh, right)` 返回新只读结果，输入不变。只允许 `[a,c]` 后接 `[c,b]`，端点必须精确相等，不以新时间容差填补间隙或吞掉重叠。要求双方旧 geometry、velocity source、characteristic source、实际 face/方向/相邻 owner 全部一致，且都在传入 mesh 上有效。

四个 integrated components 逐项相加，结果区间为 `[a,b]`。不平均两段 rate、不进行事件判断、求积、reflux 或库存更新。调用方可分段组合，但不能倒序、重复或将不同来源的两段静默相加。浮点加法有原有舍入性质，不声称任意重排都逐 bit 相等。

## 范围

该类型可存储人工或未来算子提供的积分值，构造成功不证明时空积分已经完成。I12a1 不实现事件检测、时间求积或 crossing 修复，也不改字典、输出和任何现行推进。验收与实际执行记录见 [TODO](../task/I12a1_SpaceTime_PhaseFlux_Contract_TODO.md) 和 [报告](../task/I12a1_SpaceTime_PhaseFlux_Contract_Report.md)。
