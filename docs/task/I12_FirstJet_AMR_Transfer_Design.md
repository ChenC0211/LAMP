# I12 / FirstJet 零时间 AMR 传递设计

本设计实现[纯路线](I12_Pure_Hermite_Gals_Design_Change.md)的动态网格环节。沿用现有 native Q4、exact restriction、master/hanging 约束与纯 publication；没有库存或体积闭合。

调用方显式选择两个数学操作，失败不得互相切换：

1. `ExactRefinement`：目标每叶必须包含于唯一旧叶，直接使用原 `Snapshot::refineNativeQuartic`。原 Core 已核对每叶 25 个 exact restricted datum 的恒等式。不能将 rounded 查询值重新建成 scalar；唯一返回 snapshot 就是该限制结果。新 primary/closure 是其 master first-jet、edge、mixed 和 center 的只读派生记录。数据读取使用目标原生 cell chart 的 cardinal 位置，避免先舍入显示坐标再反演局部坐标。此零时间 mixed restriction 不是独立 phi_xy 输运；下一物理步仍走原 I06。
2. `FirstJetProjection`：用于明确的粗化或一般嵌套叶重排。目标 master first jet/edge/center 从旧 accepted scalar 采样；原 I06 重算兼容 mixed，原纯 publication 构造唯一 scalar。一般分片 Q4 在粗化后不能精确表示，故明确报告近似投影和体积变化，不宣称 exact equality。全局目标多项式再现与非多项式收敛是独立验收。

两者均保持完整物理域和 XY/RZ 度量。目标 topology 的 master/hanging 归属由既有 Core skeleton 构造，不在 Interface 新写细边约束。source publication 验证旧网格，目标几何验证新网格；new identity、新网格、原时间层，accepted 事务仍由 App/调用方完成。

原 publication 的唯一 snapshot 绑定后检查提取为 private 方法。既有时间推进仍先由自身 closure 建 Q4 再调用它；exact refinement 仅由受检 mesh-transfer producer 将既有 Core 证明的 snapshot 交给它。没有公开“传任意 snapshot 就发布”的 API，也没有第二套 normal/curvature/cut/体积求积。

报告只含几何体积变化、目标 refined/coarsened/unchanged 叶计数和方法标识，不控制数值分支。投影后组件数量和边界端点类型必须保持；这些拓扑检查不是全场误差证书，不能替代后续曲率、空间收敛及动态重网格试验。丢失细结构时拒绝该 coarse candidate，由调用方保留旧状态；本接口不偷偷返回另一个方法的结果。

需要特别验证：exact regrid 后的首次 I06 闭合是否产生可接受的投影误差；不能只验证静态标量限制而不继续实际输运。若这一环节或粗化阶数失败，保留 FAIL 并据证据修复，不能放弃真实动态 AMR。
