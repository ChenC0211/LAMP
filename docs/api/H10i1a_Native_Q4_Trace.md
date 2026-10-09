# H10i1a / 同源 native Q4 只读相迹

2026-10-04；Core 来源接口，消费入口 `LevelSetMappedSegment::nativeQuarticTrace(source,begin,end,normalAxis)`，返回 move-only owning immutable `NativeQuarticTrace`。[实施设计](../task/H10i1a_Native_Quartic_Trace_Port_Design.md)、[验收](../task/H10i1a_Native_Quartic_Trace_Port_Report.md)。本接口不是材料通量、clip volume 或 publication owner。

normal 为正坐标轴 0/1，输入段必须与该轴垂直、非零、被实际 active native Q4 覆盖。输入 traversal 反转不交换物理 minus/plus。普通 constructor 的恒零拒绝保持；只有显式 factory 可查询 typed 零区间。

`ordinaryRanges()` 是普通 value 相区间，保留原 endpoint-relative Range；`coincidentRanges()` 是恒零 chart，interval 没有 negative 标记。每个 chart 的 minus/plus 包含 availability、实际 nativeCellIndex、normalOrder 和单侧相 ranges。它们来自该侧同源 1–4 阶 normal restriction 的首个非恒零多项式；negative 已包含负侧的阶数奇偶号。零判据为所有 exact coefficients 为零，无 epsilon 采样、平均、外推或第二 Q4。

source active 支持与 Hydro selected 支持不同。本 port 只认证前者；inactive/真实外侧显式 unavailable。后续 consumer 须限制 selected 支持并提供真实流入边界相数据。轴的物理测度/材料率由后续 metric owner 处理，本 port 本身不带 XY/RZ 权重。

`events()` 保留 value 或 normal derivative 的 root/chart/native 来源；normalOrder=0 表示 value incidence，normalSide=0。正长度零区间不创造 value 根序列，但查询端点还可以与一个仅接触该点的 native 闭叶相交：这类原 point-only incidence 保留，不能将其当作正长度根。各 side 的 normal 事件独立，Event 索引不表示跨侧根相等或全局 root 顺序。root.chart 属于其 value/side inventory；coincidentChart 与 normalSide 标明归属。

`rootAtInputEndpoint(index,atEnd)` 使用内部 retained identity；`refined(index)` 返回新对象，在同阶 restriction 上重做有界细化。失败保持原对象。所有 terminal offsets/widths 均为分别保留的近似表示，不是 bounds；不能从 rounded global 端点差重构细相、合并事件或证明根相等。需要跨侧 exact overlay 的后续积分必须在 Core retained chart 中完成。

`ordinary()` 仅在没有任何零 chart 时返回原 legacy view，否则拒绝，避免原 bool 相/variation 改义。所有查询检查 FP 环境及 moved-from；对象持有源，`usesSource` 比较实际 snapshot identity。源必须通过所有 alias 保持不可变。

全零二维叶、未覆盖、非轴对齐/错误轴、非 native Q4、未认证近临界根、算术范围或 bounded root/refinement 预算不足明确拒绝。高阶 normal 测量能力不降低原 publication 的正规曲线要求。[H10i1b 闭合矩](H10i1b_Native_Q4_Region_Moments.md)已由真实 RegionSource/dual mass 消费，已出版 dual 内部八组与16盒通过；native-face 八组仍在旧 publication 入口拒绝。下一 [H10i1c](../task/H10i1c_Native_Zero_Edge_Topology_TODO.md) 单独迁移正规 native 零边拓扑，完整来源链继续逐项实现。H08/H09 四硬约束与独立 Pure FirstJet 保持，H10/H11/H08a 顺序不变。
