# H10i1d / 实施前只读半面 trace cache 设计

2026-10-04；前置 H10i1c 的175 archived payload/61当前交付已独立核验。实施前目录 `runs/checks/H10i1d/20261004_111021_237`，2 production/1既有短契约/1既有research driver，0新工具/target；原 gate 保持。[TODO](H10i1d_Readonly_HalfFace_Trace_Cache_TODO.md)。

现有 `MappedFacePartitions` owns Snapshot、借用 immutable NativeLeafDomain、按 actual nativeRow/half 保留 legacy query，默认构造遇恒零拒绝，成功细化后替换单一 query。新增显式静态 `prepareNativeTrace(mesh,domain,owningSnapshot)` 返回同一 move-only cache 类型；私有构造参数选择正式来源种类，不增加生产字典或调试开关。typed 模式只接受 native EnrichedQuartic，沿用原 layout/activity 校验、actual edge vertices/face.center、首 active query 的 strict FP 检查与既有 OpenMP 安全事务；无 second Q4/root。

typed 模式每个 active actual half-face 唯一 owns `NativeQuarticTrace`。normalAxis 来自实际水平/竖直面，begin/end 仍是原 vertex→midpoint / midpoint→vertex。`traceHalfFace(row,half)` 返回只读原查询，默认 legacy cache 调用它拒绝。`halfFace` 在 typed 模式仅调用同一 query 的 `ordinary()`；positive-length zero/mixed query 保持明确拒绝，不能选 half 相或对称相冒作材料率。solid-only slot、非法索引、错 domain/活动/网格、moved-from 仍拒绝。

`rootIncidenceCount` 在 typed 模式只计 normalOrder=0 的 value events；新增 `traceEventIncidenceCount` 计全部 typed value/normal events，并要求 typed 模式。两种身份不能混用或按 terminal ranges 反造。缓存不生成任何 volume、phase measure 或材料通量。

显式 `refineTraceEvent(row,half,eventIndex)` 复用 immutable `NativeQuarticTrace::refined`。替换前核验同一 Snapshot、actual endpoints、normalAxis、gridPieceCount、value/event incidence、event normalOrder/side/chart/native cell，失败保留旧 query 和全部计数。成功失效此前该 half 的 typed/ordinary 引用，所有派生测度必须重新计算。typed ordinary `refineRoot` 先核验 legacy value root 与同位 typed value event 的原身份，再委托同一 typed refinement；不能同时维护两套不同 enclosure。含 zero query 的 legacy refineRoot 继续拒绝。

固定验收：原 default cache/flux 先保存标准构建观察，再检查 ordinary 新旧逐字段相同；原16场景追加实际半面 typed query，c 的Core/旧trace观察identity次序保持。XY/RZ、uniform/真实2:1/hanging、正负相及正非仿射 Q4，actual half/whole/subface 两侧测度 ≤1e−12 实际物理面积；对反向Core reference、normal阶数/事件、普通→zero混合、selected/外侧 availability、owning lifetime/move、同一来源细化和失败事务分别检验。缓存允许原 Core trace 支持的高阶 normal Taylor/事件，但不把它们当作已认证 regular interface publication。一般全零/unresolved/非轴对齐 trace 保持原拒绝。

T0 标准 Release/17受影响CTest（原16加既有 Hermite–GALS flux）、bclip/cside/FirstJet事务、64布局/157810普通来源及内部affine。相同旧数值输出/源输入保持，原1e−12尺度与严格FP/工作预算不放宽。旧native publication8组仍明确拒绝，完整来源消费和材料通量迁移留给e/f及H10j以后；H11/H08a等待最终空间通过。

材料硬约束仅唯一共享F、内部反号/真实AMR子面和、来源一致、KEP/skew；quarter/压力到dual精确交换及R=0是诊断。Pure FirstJet主状态/geometry/transport/clock保持独立，无库存/bubble/体积拟合/质量修复。
