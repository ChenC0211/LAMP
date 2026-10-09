# I12a6 / 同源 native Q4 矩形液量 API

`HermiteGalsPhaseRegionSource(mesh, state, config)` 保留原 immutable state/view，使用已有完整 NativeComponents 发现。复制共享 source，移动后旧 handle 拒绝访问。`source()` 返回原 state；`components()` 返回完整曲线分量数。

`measure(mesh, lower, upper)` 返回两个物理 double 端点界定的矩形负/正相测量。宽高必须严格正，坐标有限，RZ 下半径非负；矩形必须完全落在 active scalar 支撑内。精确 mesh clone 可用，任何结构变化仍由原 mesh identity 拒绝。支撑覆盖按端点顺序逐条验证，没有“面积差足够小就算覆盖”的容差。

结果含 sourceIdentity、原 bounds、capacity、liquid/gas、两种 Green 公式的 raw negative/positive volumes、GL8 误差估计、最大相对 ledger 缺陷、叶片/arc/根细化计数。结果是空间测量值，不是时空来源或库存装配。

切割叶片保留原 snapshot 的完整 graph，以目标矩形与补集对整片 source leaf 进行区域归属。固定边界与自由界面均使用目标交集左下角作为 I12a5 势原点。按交集长宽确定主 Green 方向；分别直接测量两相，采用较小相形成稳定互补。纯相证书给出严格 0/capacity。原 M 从不参加求积或被修改。

`relativeLedgerTolerance` 和 `relativeMomentEstimateTolerance` 默认均 1e-10；分别控制两个 Green 公式/容量加和差与浮点 GL8 估计，**均不等于严格积分误差界，更不改变 I12a3 的 64 eps 库存预算**。原曲线 GL8 控制仍由 Core 决定，外部 config 不会静默更改原控制规则。主矩严格越界直接拒绝，没有 clip/fallback。辅助 raw moment 供诊断，不是另一套生产几何。

每叶片最多 256 次受请求的 retained root 细化重试；替换后重做全部相关 fixed/free moments。整个 query 只写局部数据，失败不改变 prepared source 或 accepted state。

`boundaryContacts` 与 `symmetryPlanes` 原样交给原 complete discovery。现有 Core 对 partial active domain 要求显式接触点（即使是纯场）；本项保留该限制。正长度零迹、真正奇异零点、未解析图或达到工作上限仍明确拒绝。没有外推或自动补齐接触条件。

矩形端点按调用者提供的 double 定义。本接口没有证明通过舍入产生的扫掠端点等于精确时空映射，不能直接把它宣称为已完成的一般时间积分。后续扫掠任务须独立处理这一点。

见 [任务](../task/I12a6_Rectangular_Phase_Measure_TODO.md) 与 [报告](../task/I12a6_Rectangular_Phase_Measure_Report.md)。
