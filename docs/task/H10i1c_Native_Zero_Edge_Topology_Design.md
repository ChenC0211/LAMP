# H10i1c / 实施前正规 native 零边拓扑设计

2026-10-04；H10i1b 的 219 个 archived payload 与54个当前交付已核验。证据 `runs/checks/H10i1c/20261004_103851_399`。生产5文件/测试2文件上限、原正规性/精度/来源 gate 保持；仅新增显式 Core discovery，不改变旧默认入口。

`discoverNativeTraceComponents` 只接受同一 owning immutable native Q4。普通 edge、graph、negative-fragment 和 chain 使用原引擎；显式 trace 模式允许原 retained value query 表示全长 zero chart。实际 native subedge 恒零必须在每个 incident leaf 上 exact 认证，normal 一阶 restriction 无根且非零；整个原 leaf 零集的同号 normal-partial 证书排除其他 transverse 零。因为多项式沿一条正长度边恒零，它在该原 leaf 的整条 supporting side 上恒零；每条 normal 截线已有这个零，统一 crossing 方向排除第二个内部零。因此 incident leaf 的二维内部相可认证，不由几个节点或边上一阶导数单独推断。

每条实际 zero subedge 只有一个 canonical geometry owner（实际 incident cell ID 的确定性选择），保留两侧 availability/cell/phase 与真实 mesh edge/vertex 身份。端点使用 actual vertex node 与普通端点 incidence 拼接；不在无限零集上伪造 rootOrdinal。正规方向一致、内部 valence=2、声明/axis/symmetry endpoint valence=1及原≥3 graph规则保留。一般 crossing、全零叶、同边多支或缺完整证书拒绝。

新 Graph 只携带由 discovery 构造的 owning immutable `NativeSideSupport`，public 手工 flag 不能创建它。其 source/cell/actual side/subedge extent/incident phase只读；Core reader再次核验 matching Snapshot 和原 scalar 的完整证书。side graph 的普通 endpoint ordinal 保持未设置，legacy sensitivity明确拒绝。普通 graph 的原内点/normal/owner门槛不改。

attached sampling 在同一 retained rational cell chart 中沿已证实的 zero side 插值，只有 tangent coordinate 是严格内点，normal coordinate保留准确的0/1边身份。frame/normal仍从相同原 Snapshot 求值并执行原正规性/axis检查；不是新拟合或解析问题分支。geometry-only moments 对这个已认证的直线支撑用 exact retained端点及Green primitive计算，普通图继续原受控积分。来源生命周期、反向/partition/AMR实际几何、两侧方向都验收。

negative-fragment inventory把 zero subedge 的一条 interface arc放到实际液侧叶（可不同于geometry owner），另一侧不造负片；普通负相边段仍按原正长度union。内部zero边不当作共享普通negative interval。真实外侧不可用，不延拓；外边zero曲线仅在声明真实端点时允许，内部负片不存在时显式保留 `negativeRegion=-1`，不能伪造外部液库存。新结果的几何外边账本与后续typed派生volume迁移单独注明，不静默改旧入口。

固定 gate：原 H10i native-face八组的实际 Core topology/帧/测度闭式参考≤1e−12真实尺度，dual内部八组原图数据保持；XY/RZ/符号/非仿射/真2:1与hanging、多component、ordinary→side拼接、axis/外边/错源/生命周期/失败事务。标准Release16受影响CTest、H10i1b闭合矩、64普通来源和内部affine保持。旧publication八组仍记录拒绝；不生成材料F、改R/uR/KEP或进入H11/H08a/H10j。投影/remap与完整publication接线未作为本任务已经通过的能力。
