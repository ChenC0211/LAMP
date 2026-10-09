# H10i1c / 同源正规 native 零边拓扑

显式入口 `LevelSetInterfacePanels::discoverNativeTraceComponents(mesh, owningSnapshot, axisymmetric, boundaryContacts, symmetryPlanes)`，仅接受同源 immutable native Q4。[实施前设计](../task/H10i1c_Native_Zero_Edge_Topology_Design.md)。返回沿原法向定向的 `NativeComponents`，普通根与图继续用原 discovery；默认 `discoverNativeComponents` 保持沿边恒零拒绝。

恒零实际 subedge 必须在所有 active incident leaf 上同时通过完整零集的同号 normal-partial 证书，以及原 normal restriction 的严格一阶非零、无根证明。前者排除 leaf 内其他分支，不能由边上几个采样替代。一般 crossing、高阶退化、全零 leaf、未认证分支及非二价拼接继续拒绝；不增加 AMR 例外、第二个 Q4/root 引擎或 epsilon 零判据。

`Graph::nativeSide` 是 discovery 构造的只读 owning `NativeSideSupport`。它保留同一 Snapshot、canonical actual incident geometry owner、owner W/S/E/N side、normalAxis、真实 edgeLabel 和 vertex IDs、定向物理端点，以及两侧 actual nativeCellIndex/phase/availability。owner 由实际 incident cell ID 确定；每条真实 subedge 只输出一条几何 arc。普通 `Endpoint::side/rootOrdinal` 保持 −1，不能给恒零多项式伪造有限根。读端核验 scalar/cell/source 和完整公开地址，借用支撑到另一图或另一 Snapshot 拒绝。

frame/quadrature/graphFraction/graphAt 使用原 rational cell chart：tangent 是严格内点，normal 保留精确 0/1 side 身份；法向仍求值于同一个原 scalar。普通图的原 strict interior gate 保持。`nativeMoments` 对已证明的直线支撑用原端点的 exact Green primitive；XY/RZ 物理测度不变，返回浮点结果，`errorEstimate=0` 表示没有数值 quadrature 余项，不是任意算术严格区间证明。反向符号、子弧 partition 与实际 hanging geometry 分别验证。

实际液侧 fragment 持有唯一 interface arc，可能不同于 geometry owner。另一侧保持纯相证书；`cellPhase` 表示没有内部相界面的 leaf 内相，不因边附着 graph 置零。负片只通过正长度同相边段连接。active/physical 外侧 unavailable 不外推；真实声明边界曲线若液相只在 unavailable 外侧，`negativeRegion=-1`，不制造外部液库存。axis/symmetry/declaration 和原 valence/最少 graph 数规则保持。

本入口只完成 Core topology、frame/sampling 和 geometry moments。旧 sensitivity 明确拒绝 attached graph；`integrateNativeQuarticRegionMoments` 尚未消费该种 graph。半面 cache、FirstJet 派生 volume/publication 和 H07 消费分别归 H10i1d/e/f。projection/load/remap 未在本任务独立验收，不能凭此 API 宣称完整出版迁移、材料通量或 PDE 通过。Pure FirstJet primary/Q4/transport/clock、材料 F/R/uR/KEP 均不改变。
