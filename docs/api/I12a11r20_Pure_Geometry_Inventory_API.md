# I12a11r20 / 纯相几何与库存差 API 补充

本项不增加调用入口。`HermiteGalsBubbleVolumeConstraint::GeometryAccuracy` 的几何容限现在同时作用于 cut 与经过原 scalar 认证的 pure cell：`abs(V_geometry-M)/capacity <= config.relativeVolumeTolerance`。默认仍为 1e-10。

纯相 center 保持预测值，不求根、不制造微滴或改变几何相位。非端点库存与几何纯相并不等价：库存是独立守恒量，其预算内偏差由同一 scalar 上的有界库存积分表示承载。`relativePureVolumeTolerance`、`maxPureRelativeResidual` 和原逐 cell 的 target/geometric volume 如实报告该差额。

超出几何预算时拒绝，错误为 `pure-cell mismatch exceeds geometry volume budget`。原 M 严格容量界限、未解析成分检测、拓扑认证和 cut 修正的 Hessian 预算不变。`InventoryRoundoff` 的 pure 门槛仍为 `min(relativeVolumeTolerance,64*epsilon)`，数值行为保持。

I11 的实现没有修改；发布仍独立比较 raw 几何积分与原 M，并保留 state identity、first jets 与独立库存。仅当各独立容限都满足时成功发布；收紧发布容限导致拒绝时，candidate 保持可用。

实际输运必须消费 [r18 库存区域](I12a11r18_Bounded_Inventory_Region_API.md) 和 [r19 共享面传递](I12a11r19_Bounded_Inventory_Sweep_API.md)，不能把允许差额的结果送入只按几何体积输运的旧路径。I10 本身不认证未来调用者的选择；本项测试验证了真实 I11 publication→r18→r19→原库存→reflux 链路。完整 first-jet 半步/全步驱动尚未切换。

依据[设计补充](../task/I12_Bounded_Inventory_Transport_Design.md)，参见 [TODO](../task/I12a11r20_Pure_Geometry_Inventory_TODO.md) 与[报告](../task/I12a11r20_Pure_Geometry_Inventory_Report.md)。
