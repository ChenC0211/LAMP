# H10i1e / FirstJet trace 派生测量与发布

`InterfaceCellMeasure::discoverNativeTraceVolumes(mesh, domain, owningSnapshot, contacts, symmetry)` 显式调用原 Core complete native trace discovery 与 `MappedFacePartitions::prepareNativeTrace`，返回完整几何测量及两方向账本。仅支持原 native EnrichedQuartic；原 `discoverVolumes` / bubble / sensitivity 入口保持 legacy 语义。

`measureNativeTraceVolumes(mesh, domain, partitions, topology, contacts, symmetry)` 复用同一 immutable cache/topology；所有 source 地址、layout、活动集、complete component/graph/phase/side 来源验证必须通过。公开拓扑的完整性由同一 Snapshot 的原 Core discovery 重新认证，不能把手改 phase vector 或缺少 component 的输入当作合法证书。无 Snapshot 复制、不同标量、供调用者修改的源或库存。

`DiscoveredNativeTraceVolumes` 保留原 `Volumes` 的 axial fixed/free/phase volume、误差估计及统计；`NativeTraceLedger` 另外保存原全局起点下的 radial fixed/free/capacity、实际 assembled axial/radial volume 与闭合诊断。纯内部相的最终 phase volume 来自完整证书和实际 CV capacity，可能因浮点舍入与 global ledger 略有不同；两者分开记录。cut 内部分支采用原 local Green/smaller-phase integration，容限与 strict capacity 不变。

fixed phase ledger 仅包括原 ordinary negative ranges。整段 zero support 没有 half phase；正规 side graph 从原 `nativeMoments` 得到 free contribution，归其真实可用液侧 CV。canonical graph owner 可以是气侧。没有液侧时保留合法几何，但不制造外部液体/体积。shared fixed 面只积分一次并反号写入相邻 CV，真实 AMR subedge 不合并；free side 不重复计算。两方向容量、fixed closure 与 complete fixed+free 账本均显式核验。

`integrateHermiteGalsPhaseVolumes(mesh, state, config)` 根据原 transport 明确分派：FirstJet 使用 native trace 测量；ConservativeVolume 保持原测量。`integrateHermiteGalsNativeTracePhaseVolumes(mesh, state, domain, partitions, topology, config)` 是发布事务的显式复用入口，只接受 FirstJet 与同一 Snapshot。`result.nativeTraceLedger()` 仅在 FirstJet 可用；historical port 明确拒绝，不返回空/伪零诊断。没有初始化 M、几何修复或从诊断回写 primary。

FirstJet 发布 Data 在 final allocation 内先固定 domain，再绑定同一 typed cache/topology，随后测量；体积结果、raw frame/sample、quadrature、half query 全部属于最终返回的原 source/revision/time。所有检查完成后才 consume candidate。失败保留 candidate 和既有 accepted source/clock；成功保留原 phi/grad、mixed closure、预测 center 与 identity。

正规 XY/RZ uniform/AMR native 面可以实际发布；ordinary→side、多分量、外部、selected 支撑有短契约。原 Core 对临界/全零/高阶 normal 退化及 RZ 轴与域角点的边界身份拒绝保持；后者不能计为已支持边界。普通 graph 的 root/FP/controlled 预算未改变。

本接口不是 H07 face material rate、因果相通量或 Hydro 动量质量入口。H10i1f 迁移正式几何消费链，H10j 再构造空间 phase-rate owner；不得把无速度空间迹当作因果材料 F。仅四类材料硬约束保持，quarter/精确交换/R=0 仅诊断。H10/H11/H08a 尚未闭合。
