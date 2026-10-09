# I12a11r17 / 几何精度预算 bubble API

后续变更：本文记录 r17 原契约；GeometryAccuracy 的纯相规则已由 [r20](I12a11r20_Pure_Geometry_Inventory_API.md) 扩展为同一几何容限，旧 InventoryRoundoff 保持。cut 精度预算与求根规则不变。

入口仍为 `closeHermiteGalsBubble(mesh,candidate,mixed,prediction,config)`，使用同一 I10 内核和只读结果。

显式选择 `config.volumeConstraint=HermiteGalsBubbleVolumeConstraint::GeometryAccuracy`，并按 `candidate.primary().geometryCells` 完整顺序提供 `config.accuracyBudgets`。每项为 `{cell, maximumHessianComponentChange}`；预算是物理坐标中 scalar Hessian 任一分量的最大新增变化，必须有限且非负。零预算固定预测 center，仍允许已有几何误差在请求的体积容限内。

cut cell 使用 `config.relativeVolumeTolerance` 作为相对容量的几何残差门槛，不再暗中压至库存 roundoff。预测已达标则不改 center。否则原 Core 拓扑/位移证书与如下精度区间相交：

`abs(c-c_pred) <= H_budget / max(32/hx^2,32/hy^2,16/(hx*hy))`。

端点向内处理浮点舍入；对称/单向证书、bracket 和最终实际存储值均受同一原始预算约束。原 safeguarded root、单调性、求积工作量和拓扑检查复用。目标在区间中不能达到请求的几何容限时明确拒绝，M/first jets/mixed/identity 均不修改。

`HermiteGalsBubbleCellReport` 增加请求的 `maximumHessianComponentChange` 与实际 `hessianComponentChangeBound`；总报告增加 `volumeConstraint` 和 `maxHessianComponentChangeBound`。原 initial/final volume、target、signed/L1 residual、迭代及修正报告保持。该界控制新增 scalar Hessian 扰动，不是完整曲率误差上界；RZ/曲率与实际收敛仍需独立验收。

本项仅放宽 cut cell：pure 仍使用 `min(relativeVolumeTolerance,64*epsilon)`，未解析成分和拓扑拒绝保持。缺失、乱序、非法预算拒绝。`InventoryRoundoff` 契约要求预算列表为空，保留旧数值行为与默认值用于过渡期显式 A/B；没有失败后自动切换。

完整驱动尚未切换。几何—库存差超过机器精度后，必须先完成共享通量有界性协调，不能直接提高旧 roundoff 修复常数。现有 I11 可以发布满足其独立几何容限的结果；测试中的单次发布不等于完整时间推进或 Gate I 通过。

见[授权变更](../task/I12_Geometry_Inventory_Accuracy_Design_Change.md)、[TODO](../task/I12a11r17_Accuracy_Limited_Bubble_TODO.md) 与[报告](../task/I12a11r17_Accuracy_Limited_Bubble_Report.md)。
