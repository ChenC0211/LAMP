# H10i1d / retained refinement 身份核对

2026-10-04；首轮短契约 `cache_preflight` 保存了 wrapper 在真实普通 value event 细化时的错误拒绝。检查原 `MappedHermiteSegment::BasicQuery::refineRoot`：同一个原 root 的 exact parent bracket 经原引擎认证后追加 child chart，公开 `Root::chart/parameter/lower/upper` 属于这个新 chart，不能要求 chart ID 不变或跨 chart 比较 terminal lower/upper。

正式稳定身份是同一 owning scalar、actual input segment/normalAxis、原 root/event index及source native cell、typed normalOrder/normalSide/coincidentChart、value/event/grid-piece inventory。`NativeQuarticTrace::refined` 继续认证同一个 exact root，原 shared incidences/phase endpoints统一更新。wrapper保留这些核验，只移除误把新 child chart 当作身份变化的断言；ordinary legacy view 与 typed event 在当前同一 chart 的对应检查仍保留。独立非二进制normal根物理相测度仍按原1e−12 actual scale核验，原Core工作预算及严格FP不改。

该说明补充已保存的[实施前设计](H10i1d_Readonly_HalfFace_Trace_Cache_Design.md)，不改写它。不是放宽同源/same-root要求、扩大数值容限或替换root engine。
