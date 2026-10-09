# H10i1e / 身份校验与原边界范围核对

实施前设计保持原样。本补充记录具体接口与已观察的限制；不修改数学门槛。

`NativeComponents::cellPhase` 是公开 vector，调用者可以修改；因此 `measureNativeTraceVolumes` 除了检查 source/domain，还在同一 Snapshot 上调用原 Core complete discovery，比较全部 phase/component/graph/side 身份与数量。新增 contacts/symmetry 参数用于同样的边界声明。复核采用同步借用的 shared_ptr alias，缓存持有原 allocation 的实际所有权，未复制 Snapshot 或另造标量；验证期间不接受调用者 phase mask。发布测量读取最终 publication 内相同实际 half cache 和 topology；复核产生的临时拓扑只用于身份核对。

两方向账本保存真正的全局 fixed + free axial/radial 浮点结果，和由原完整纯相证书或原 local integration 得到的 phase volume 分开。研究 CSV 的 `axial_liquid` 必须记录前者，不能用证书体积冒充原轴向 raw ledger。改变此观察字段后最终二进制/回归重建复核。

扩展边界 preflight 中，RZ 轴与物理域上下角点的交会被原 Core 拒绝：`LevelSetNativeInterfaceChain` 的 axis endpoint 条件要求 z 严格处于域内部，且 non-axis contacts 要求 r>0；二者在角点无法同时满足。本任务不修改 Core。4 个实际 FirstJet 候选的原拒绝及状态保持被保留为短契约；XY 全边界、RZ 外径边界与 XY/RZ selected 边界仍须实际通过。原 8 组 native-face 的 RZ 轴端点处于域内部，可实际发布。该角点尚未认证，不计为成功；后续如需扩展 production 边界，应先独立 Core boundary identity 任务。

物理/selected 外边界有时仅存在气侧，Core 可以保留零边几何而没有 negativeRegion。此时不得给缺失的外部液体添加 free moment。短契约记录 sideGraphs 与实际 sideLiquidOwners 的差别，气域测量为 0，未凭 source 延拓制造液相。
