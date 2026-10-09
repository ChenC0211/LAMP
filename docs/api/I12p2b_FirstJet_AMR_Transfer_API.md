# I12p2b / FirstJet 动态 AMR API

[设计](../task/I12_FirstJet_AMR_Transfer_Design.md) / [任务](../task/I12p2b_FirstJet_AMR_Transfer_TODO.md) / [接口](../../include/Modules/Interface/HermiteGalsMeshTransfer.h)。

```cpp
auto result = transferHermiteGalsFirstJetMesh(
    oldMesh, accepted, newMesh,
    HermiteGalsMeshTransferKind::ExactRefinement, config);
const auto& next = result.publication();
```

调用方必须显式选择 `ExactRefinement` 或 `FirstJetProjection`。前者只接受目标叶包含于唯一旧叶的同网格/细化，复用 Core 的 exact Q4 restriction，并直接发布该 snapshot；不会从 rounded primary/closure 记录重建“近似相等”的几何。后者允许同域嵌套叶粗化/重排，从旧 scalar 采样目标 master first jet、edge 和 center，调用原 I06 与纯 publication。一般分片 Q4 的粗化是近似操作，不承诺 exact equality 或体积守恒。

源 publication 必须是无 M 的 FirstJet、完整 native Q4 支持，且与 `oldMesh` 相符。两个网格必须同物理域，所有相交叶须嵌套。XY/RZ 度量继承源状态；调用方通过 `config` 提供目标发布所需的边界接触点和对称面，不自动推断新的边界条件。

返回对象只读：`sourceIdentity()` 记录源 revision，`publication()` 有新 identity，保留原 step/begin/end/position，网格切换不推进物理时间。exact 模式中的 mixed/center 是零时间 scalar restriction 的派生记录；下一实际时间步仍由 I06 生成 mixed closure。

`report()` 的叶计数按目标叶统计：同边界叶计为 unchanged，被某旧叶包含的更细叶计为 refined，覆盖多个旧叶的粗叶计为 coarsened。`sourceVolume`、`targetVolume` 为几何测量，`relativeDomainVolumeChange=(targetVolume-sourceVolume)/sourceDomainCapacity`，分母是整个物理域的 XY/RZ 测度，不是液体初始体积。时间驱动中的 `(Vn-V0)/V0` 是另一项诊断，不得混用。

cut、拓扑、法向、曲率和体积共享返回 publication 的唯一 scalar。传递还检查闭合/开口分量及无方向端点类型的多重集合不变；这只能拒绝已检测到的分量/端点丢失，不证明每个分量的位置、形状或误差满足某个上界。过度粗化可能被拒绝，不能静默删除小液滴。

失败不修改源 publication 或外部网格，也不改用另一传递模式。此接口只准备新 publication；何时替换 accepted state/mesh 属于 App 或调用方的同一事务。本模块接口尚不等于 production App 接线或 Hydro 耦合。
