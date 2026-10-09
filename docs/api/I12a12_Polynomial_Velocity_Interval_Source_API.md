# I12a12 / 同源多项式速度区间源 API

本接口供一般运动区间源使用。只提供速度/一阶空间导数、连续 box 包围和现有 tracer 的同源 endpoint/substage views；不提供相事件证书或一般时空液量积分。实施范围见 [TODO](../task/I12a12_Polynomial_Velocity_Interval_Source_TODO.md) 与[设计](../task/I12_General_Velocity_Interval_Source_Design.md)。

## Core 数据与求值

`Quantum::Geometry::PolynomialVelocityData` 包含 `origin`、正 `scale` 和 `previous/arrival` 两份 16 个 `Vector2D` 系数。索引 `4*i+j` 表示 `xi^i eta^j`，每个 Vector2D 存两个速度分量。`PolynomialVelocitySlab(begin,end,data)` 校验后复制数据；`data()` 只读。

`evaluate(point,time)` 返回 `LevelSetBoundaryTrace::VelocityJet`，其 derivative 顺序为 xx/xy/yx/yy。空间值和导数来自同一 Horner 递推，导数已换回物理坐标。时间采用两 endpoint 的线性插值；old/end 时直接取该 endpoint，归一化时间恰为 .5 时用与现有 Core 一致的平均表达。非线性真实时间场需要明确其重构误差，不能把该插值称为任意时刻的连续真解。

`enclose(lower,upper,begin,end)` 接受闭 box，可退化为面/点或单一时刻。返回两个 value ranges 与四个 derivative ranges，分别含 lower/upper。strict out-of-line 代码复用既有 outward interval arithmetic，并包含实际求值的端点/中点分支。范围可以保守；它不证明 RK2 轨迹保持在 box 中，也不证明未采样的相区间为空。

所有范围要求有限、排序正确且时间处于 slab 内。非有限系数/坐标、非正 scale、无效时间或算术 overflow 明确拒绝；非最近舍入或 FTZ/DAZ 模式拒绝，不尝试切换这些控制模式。空间多项式有全局代数定义；实际物理域、半 cell 位移和 RZ 轴条件仍由原 tracer 检查。

## Interface 来源

`HermiteGalsTemporalVelocitySource::polynomial(old,arrivalIdentity,arrivalTime,data)` 创建 owned slab，并从该同一对象生成实际 old/end callbacks。old time/identity 及实际到达时间继续经原 source constructor 校验。

`polynomialVelocity()` 返回该 owned slab 的 const pointer。source 拷贝保留相同 owner；另行构造同值数据仍是不同 source。moved-from 查询拒绝；原 opaque 和 constant source 返回 null，不能凭采样值获得多项式能力。

`polynomialArrivalView(actualCandidate,actualArrival)` 产生同一来源的真实 substage view。begin、step、end extent 与实际物理时间必须处于原来源允许范围，actual candidate 必须为非零 ID。End 和 Midpoint 的时间语义保持原 `HermiteGalsTimeLayer` 规则；该 view 不重建 geometry、不分配假 geometry ID，并持有 owner 以保证 source 局部变量销毁后的查询有效。

原 constant 工厂、`constantArrivalView()` 和已有平移 producer 均不变。显式 polynomial source 即使系数碰巧为常量，也不自动获得 constant capability；不得把它直接送进只认证常速度的 phase-flux producer。本接口不会修改 M、发布 geometry 或提交 accepted step。

## 消费顺序

1. 构造单一来源，保留同一 previous/arrival geometry/time 及 characteristic owner。
2. first-jet transport 使用其实际 endpoint/substage views，仍调用原 `traceHermiteGalsWithJacobian`。
3. 后续事件算子从同一 `polynomialVelocity()` 取得空间—时间速度范围，并独立证明实际离散 RK2 的覆盖。
4. 只有事件/求积和共享库存验收通过后，才能用于一般运动完整推进。I12a12 本身不满足这一步。
