# I12p3a / FirstJet 时间步 API

[任务](../task/I12p3a_FirstJet_Step_Producer_TODO.md) / [接口](../../include/Modules/Interface/HermiteGalsFirstJetStep.h)。

```cpp
auto prepared = prepareHermiteGalsFirstJetStep(
    mesh, accepted, traceDomain, previousVelocity, arrivalVelocity,
    prescribedOuterFirstJet, publicationConfig);
const auto& candidate = prepared.publication();
```

输入源必须是同网格、完整 native Q4 支持的 FirstJet publication，无 M，时间位置为 Initial 或 End。目标时钟来自 `arrivalVelocity.timeLayer()`，必须是源之后的下一 step，position 为 Midpoint 或 End、begin 等于源实际时间。previous velocity 的 geometry identity 和完整时间层必须与源一致。域、RZ 轴语义和 cell scale 复用原 first-jet 校验；原 tracer 的半 cell CFL 限制保持。

`arrivalVelocity` 的 geometry identity 表示该速度视图的来源，可以与源以及新发布的 identity 不同。结果同时保留 `sourceIdentity()`、`arrivalVelocityIdentity()`；新 `publication()` 的时间层等于目标时间层。该接口只准备候选，App 决定是否替换 accepted handle，不在模块内部 commit 或改变主时钟。

`HermiteGalsBoundaryFirstJetQuery` 接收外边界到达点、目标时间层和实际 trace 命中的边界 side，返回该到达点/时间的 phi 与完整物理梯度。只在到达点位于实际物理外边界且 trace 确实命中边界时调用。它不是域外延拓或内部解析重置：内部 vertex/edge 命中仍由原 pullback 拒绝，center 命中明确拒绝；缺失/非有限边界数据或 callback 抛错也拒绝。provider 捕获的状态须保持只读、存活和可重复。

模块按原顺序对 vertex/edge 做同一 pullback：phi 从旧 scalar 的 departure 值取得，梯度为该同一映射 J 的转置乘旧梯度，中心用同一 departure 值，再调用原 I06 和纯 publication。没有新增 RK2/J 公式、mixed stencil、库存闭合或 normal/curvature 算法；没有隐式方法切换。

H08a 增加中立 `HermiteGalsResolvedMotion` 重载，绑定源 identity/完整 clock、速度来源标识、目标 clock 和只读 point/J 查询。完整内部结果必须覆盖精确目标时间段且有有限 J；物理边界命中须有 side 且不得伪造内部 J。Interface 不依赖 Hydro；实际 CCH source 的来源、父区间和网格验证由 Hydro owner 保留。原 RK2 重载只作同一中立消费的适配，保留原 clock/domain/CFL 及边界拒绝行为。

中立重载默认 `InflowHistory`，与上文原规则相同；显式 `PrescribedOuterFirstJet` 要求 callback 提供全部外边界点的目标 FirstJet，包括切向、出流及角点（角点使用对应组合 side）。该选项表达给定物理边界数据，不是内部标量解析重置或域外延拓。H08a 提供 boundary 时采用该显式规则。新入口已通过 [H08a 规定范围验收](../task/H08a_SpaceTime_Material_Flux_Report.md)，完整流体 App 与全耦合二阶仍未验收。

`report()` 记录 transportedVertices/Edges、prescribedVertices/Edges 和 tracedCenters。它们是实际处理的 DOF 数量，不是 trace 调用次数或误差证书；例如外边界检查可能发生额外 trace。报告不参与数值分支。

失败保留源 publication 和外部网格；成功返回只读候选与唯一 scalar。配置、时钟、给定速度、真实边界数据及 IO 留在 App/调用方。现有独立基准通过此接口运行，但正式 App 字典/输出仍属后续 I12p3b/I12p3c。
