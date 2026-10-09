# I12 / 同一特征与旧相的完整覆盖

[I12a13](I12a13_RK2_Characteristic_Enclosure_Report.md) 给出实际完整单步 RK2 的 departure box，保留同一 old/velocity/characteristic source。本项把它连接到已存在的 `Snapshot::nativeQuarticBoxSign` 及 [r18 库存权重](I12_Bounded_Inventory_Transport_Design.md)。本项不改变任何实际回溯或几何。

## 从 B_h box 到相

设 Q 为已认证 departure box，旧 native leaf 为 L_c。枚举 native leaf search 的所有 closed-box 候选，求 `Q ∩ L_c`；保留零宽度的边/点交集。对每个交集调用同一 old snapshot 的严格 Q4 Bernstein 符号检查，不能只看 donor centroid、端点或若干时间节点。若每叶均为同一严格符号，则所有 `B_h(X,t)` 都在该相内；否则结果为 Unresolved。Unresolved 不证明实际 B_h image 确有零点，因为 box 可能保守。

额外核验这些矩形的并集确实覆盖 Q：取所有 x 端点的有序唯一切分，逐个闭 x 条带选出完整覆盖该条带的矩形；按 y 下界排序检查其闭 y 区间并集是否覆盖 Q 的 y 范围。退化 x 使用其唯一坐标作条带，退化 y 同样通过区间包含检查。算法只比较原有坐标，不取中点、累加面积或引入容差；仅在接缝接触的叶也保留用于相与库存范围。

该检查与原 native source 的完整 active scalar 支持共同建立覆盖，不允许将 inactive 或缺失标量当纯气相。Q4 API 已负责 outward chart conversion、authoritative Bernstein controls 及浮点环境要求，本层不复制这套算术。

## 独立库存的含量范围

几何入口使用 `wl=1,wg=0`。显式库存入口仅接受 source identity 与 snapshot 指针都等于 enclosure.old 的 `HermiteGalsInventoryRegionSource`，保留其 owner，并使用其已校验权重。

每叶负相的含量范围为 `[wl,wl]`，正相为 `[wg,wg]`，未确定为 `[wg,wl]`。所有叶的最小下界/最大上界形成整体范围。该操作只选取已有有限数，不做权重重算、舍入修复或 M reset。不同叶权重可跳变，所以跨叶时即使几何纯相，含量也未必为单个常数；0 符号返回不能把该叶排除。

这是 r18 的 `rho = wl H + wg (1-H)` 在同一旧几何上的范围，不是第二份界面。纯气相可以有预算内小量独立液体库存，纯液相含量也可能略低于 1；不能继续用机器精度几何体积硬匹配消掉这些差异。

## 后续消费边界

结果为私有构造、只读、来源保留的覆盖对象，记录每叶交集、严格相状态和含量上下界。它既不能伪造时间节点 geometry，也不自动生成 face transfer。未来 quadrature 可利用纯相/含量范围和速度范围界定误差，在未确定区域继续细分或累计有证明的误差预算；不能仅以粗细求积相近宣告无薄层事件。

本项不认证边界历史、缩步过程、连续 ODE 或完整时空积分。I12 的 static/translation/rotation/vortex/AMR/curvature 原验收矩阵保持。
