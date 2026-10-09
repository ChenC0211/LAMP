# I12 / 一般运动的同源速度区间数据

[r23](I12a11r23_Mesh_Scaled_Geometry_Budget_Report.md) 已通过非多项式平移空间收敛。下一步仍按[已授权的区间通量路线](I12_Phase_Flux_Time_Integration_Design_Change.md)推进旋转/涡旋，保持现有 RK2 B/J；本文明确必要的数据表示，不改变 characteristic 积分方法，也不提前批准未实现的事件/库存算法。

## 当前缺口

`HermiteGalsTemporalVelocitySource` 保存两个不透明 endpoint callbacks。只有 owned constant 工厂能从数据本身证明中间时间值并生成同源 substage view。单靠若干点处的 callback 值，不能排除采样之间的速度变号或相事件；给任意 callback 贴一个“多项式”标签也不能解决这个问题。

现有 Core tracer 在两个 endpoint 之间线性插值速度及其同源空间导数，使用 explicit midpoint RK2，并在边界分支中缩短 substep。I12a12 为这个已有时间插值提供明确的 owned 空间多项式数据和包围能力，不替换实际 tracer，不把另一个解析旋转映射冒充 B_h。

## 唯一数据与计算

固定坐标 `xi=(x-origin.x)/scale.x`、`eta=(y-origin.y)/scale.y`，每个速度分量的两个 endpoint 都由 4×4 个 power coefficients 定义。时间插值为 `(1-theta)*u_old + theta*u_end`，theta 由该 slab 的真实 old/end 时间决定。尺度、系数和时间均在构造时复制并校验，后续不可写。

这是一份只读重构输入，不是额外全局速度未知量。后续 Hydro 仍拥有 native face velocity；本任务只为 Interface verification 和未来明确的重构视图建立消费契约。全局 cubic 数据能覆盖现有 rotation 和 differential-rotation 的空间速度；后续 native piecewise reconstruction 仍须明确其自身覆盖和边界契约。

值与空间导数由同一多项式 Horner/自动一阶递推获得，不另拟合梯度。Interface 工厂的 endpoint callbacks 捕获同一不可变 Core slab；substage view 捕获该 owner，并携带真实 candidate ID 与时间。原 opaque/constant 入口不变，opaque source 无法凭采样结果获得新能力。

## 包围的含义

Core 的 out-of-line strict/no-LTO 实现复用既有 outward interval arithmetic，对矩形空间 box 及时间区间计算值/导数范围，包含浮点坐标换算及时间混合的端点/中点分支。任何 overflow、非有限值、非最近舍入或 FTZ/DAZ 环境都明确拒绝；不主动切换调用者的舍入或 FTZ/DAZ 控制模式。

区间范围可以保守，空间/时间进一步细分可降低依赖膨胀。本项只给出速度多项式的包围，不能从它直接推断离散 RK2 的所有中间点、边界分支或旧 Hermite 相符号。下一项必须把这些范围沿同一个实际 RK2 运算传播，处理跨叶、速度符号、纯相/界面和边界覆盖，才能开始一般时空 phase-flux 积分。

包围证明依赖逐运算区间包含性和 strict 浮点前提。测试中的独立解析极值及密集采样只是核验实现；不能将“未采到反例”写成证明。

后续传播必须包含 Core 对实际 endpoint 查询的再次时间混合、RK 中点和边界 substep 运算。实数层面相同的线性时间函数可能因运算顺序产生舍入差，不能把直接 `slab.evaluate(point,time)` 与 Core 内部混合声称为逐位相同，也不能只用直接 slab 点值代替实际离散阶段的范围。

## 物理时间与比较

时间上有真实非线性变化的外部场（例如 cosine amplitude vortex），本源表示的是两个 endpoint 定义的线性时间重构，不能声称中间值等于连续 cosine 真解。后续 benchmark 必须明确这一离散输入，仍以原连续真解检查时间阶。完整 endpoint trace 沿用 Core 已有插值；半步/全步共享同一个线性 slab，避免各 quadrature node 偷换来源。

原 callback/constant 算子与所有历史结果保留用于明确 A/B。I12a12 的验收只证明数据、包围和 RK2 coherence；一般通量有界性、真正的时空事件覆盖、一般时间收敛与 Gate I 都尚未完成。

实施记录（2026-10-02）：[I12a12 报告](I12a12_Polynomial_Velocity_Interval_Source_Report.md)记录算子及 20 项回归 PASS，18 项旧数值输出不变。完整 I12 仍 FAIL。
