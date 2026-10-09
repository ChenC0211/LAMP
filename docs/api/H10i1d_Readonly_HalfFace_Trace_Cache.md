# H10i1d / 同源只读半面 trace cache

`InterfacePhaseGeometry::MappedFacePartitions::prepareNativeTrace(mesh, domain, owningSnapshot)` 显式准备 native Q4 cache。[实施前设计](../task/H10i1d_Readonly_HalfFace_Trace_Cache_Design.md)、[refinement 身份核对](../task/H10i1d_Refinement_Chart_Identity_Note.md)。默认 constructor 继续准备 legacy cache，恒零拒绝保持；新 factory 不产生任何 volume 或材料率。

两种 cache 均 move-only、owns 同一个 immutable Snapshot、借用同一个 immutable NativeLeafDomain。调用者保证 domain 存活，变更 mesh/domain 后必须先 `validate`；错 domain 地址、活动集合/XY/RZ/leaf 布局、solid-only slot、非法 nativeRow/half、moved-from 拒绝。typed factory 仅接受 native EnrichedQuartic，沿用原 Core source/FP/工作预算。

每个真实 native face 的两个查询仍是 actual vertex→face.center 与 face.center→vertex，保留原物理端点和 normalAxis（竖直=0、水平=1）。每个 active half 唯一 owns 原 `NativeQuarticTrace`；`traceHalfFace(row,half)` 返回只读原查询，保留 value 与 zero charts、两侧 availability/cell/normalOrder/ranges 和 normal events。完整零多项式的无穷根不变成 finite root inventory；inactive/外部相不外推。

`halfFace(row,half)` 在 typed 模式直接返回同一查询的 `ordinary()`；含正长度 zero support 时明确拒绝，不用对称 half 分相填补。legacy cache 的 `traceHalfFace` / typed event count / refinement 明确拒绝。`rootIncidenceCount` 仅计 value（normalOrder=0）incidences；`traceEventIncidenceCount` 是全部 value/normal events，并仅对 typed cache 可用。这些是事件统计，不是 phase rate、mass 或 causal flux。

`refineTraceEvent(row,half,eventIndex)` 复用原 immutable same-root refinement，成功后替换唯一 typed query，失效该 half 此前返回的 typed/ordinary 引用。派生测度必须丢弃并重新计算；source、actual segment、normalAxis、value/event/grid-piece inventory 和 source cell/order/side/coincidentChart 身份保持。公开 Root 的 chart/parameter/lower/upper 属于新的 child chart，不能要求 chart ID 不变或跨 chart 比较 bounds。原 Core 负责 parent exact bracket 的同根认证。任意失败保留原 query 和计数。

typed ordinary `refineRoot` 核对当前 legacy root 与 typed value event 的对应后，委托同一个 typed refinement；不会留下两个不同步的 enclosure。重合 query 的 legacy refineRoot 拒绝；有 normal coefficient events 的查询用 `refineTraceEvent`，不会混入普通 root ordinal。

普通→zero 转换由相邻真实 half/subedge incidence 表示，不为测试合并出非实际 face。RZ 水平面的两个几何等长 halves 一般不是各占半个物理面积，验收使用各自真实 `pi*(r_hi²-r_lo²)`；竖直面使用 `2*pi*r*dz`。高阶 normal Taylor 与临界 normal coefficient events 是合法 typed trace 数据，不因此被当作 regular interface publication。

H10i1d 只迁移 cache。新 factory 仍未被 FirstJet publication/H07 材料路径采用；旧 native-face publication 拒绝保持。后续 H10i1e/f 分别迁移派生测量/正式消费，再复核 H10i 与 H10j。Pure FirstJet primary/Q4/transport/clock、F/R/uR/KEP不变；无库存/bubble/体积拟合/修复，quarter、精确交换及 R=0 仅诊断。
