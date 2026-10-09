# H10i1a / NativeQuarticTrace API 与实施前设计

2026-10-04；[任务书](H10i1a_Native_Quartic_Trace_Port_TODO.md)。前置 H10i1 的 108 payload 已核验，上一轮为 progress。本轮只实现只读 Core 来源能力，4 production / 1 既有测试，不迁移几何或材料消费者。

在 `LevelSetMappedSegment` owner 增加显式 static `nativeQuarticTrace(source,begin,end,normalAxis)` factory，返回独立 move-only owning immutable `NativeQuarticTrace`。只接受 native Q4 及 normalAxis 0/1 对应的非零轴对齐段；normal 为正坐标轴，不随 input traversal 反转。结果不包含速度、密度、物理时钟或接受状态。

普通相区间与正长度恒零区间分型。普通区间保留原 Range 的 endpoint-relative offset/length；零区间逐实际 grid/chart 记录其内部侧、native cell、availability、首个非恒零 normal Taylor order 与该系数的相区间。minus/plus 分别对应负/正法向侧，负侧偶奇号由 order 决定；缺外侧保持 unavailable。该 Core 结果只描述源的 active 支持，selected 支持由后续 Hydro consumer 显式限制。不可从 scalar 的外侧数据暗自推定材料边界条件。

所有 restriction 系数 exact 为零才允许零区间。逐 active incident leaf 验证相同零迹；普通共享迹仍用原正 clearing factor proportionality 认证。每侧用同一 `restrictLocal` 的 1–4 阶 normal restriction 选择首条非恒零多项式；二维恒零、未覆盖、算术范围/预算不足明确拒绝。normal Taylor 根复用原 bounded isolation/square-free/dyadic/refinement，引入一个内部 fixed native chart/order layout：所有子区间和 refinement 均重新限制同阶偏导，不能落回 value phi。不存在第二 Q4 或第二 root engine。

内部保留所有 rational grid cuts、normal root charts 和独立两侧事件。公开 Event 标明 value 或 normal order、normal side、coincident chart 及 native root identity。两侧事件不按 rounded global root 合并；未来需要跨侧 overlay 的积分在 Core 持有的 retained charts 上认证，不能从 terminal ranges 反推精确事件相同。孤立 normal 根只作为事件，不制造 free interface 或无限 value 根库存。`rootAtInputEndpoint` 用 retained identity；`refined(eventIndex)` 复制候选并只收紧同源同阶根，失败不影响原对象。

无零区间时 `ordinary()` 提供原 legacy 查询；有零区间明确拒绝该 view，避免原 bool Range/variation 被零迹改义。原 constructor 的输出、zero rejection 和所有 production consumers 保持。源 owning、wrong source、moved-from 与 FP 环境检查采用既有契约。

实施与验证按任务书：先现有 root enclosure 短契约的实际新 port 输出，包括高阶/非仿射零迹、导数事件、小相、原生面/AMR/缺侧/方向/分片与事务；再研究 driver 读真实已发布八组 dual 场景 1,408 段/40 零段，保留原 publication/mass/clip 拒绝记录。原普通来源输出和 affine/KEP 回归保持。只用独立解析多项式/原 legacy 查询作对照，不以新 kernel 自证；不声称裁剪质量、F、最终 PDE 或 H10 已完成。
