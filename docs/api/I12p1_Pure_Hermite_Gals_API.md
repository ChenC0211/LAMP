# I12p1 / 纯 Hermite–GALS API

[设计变更](../task/I12_Pure_Hermite_Gals_Design_Change.md) / [任务](../task/I12p1_Pure_Hermite_Gals_Publication_TODO.md)。新 API 位于原 [State](../../include/Modules/Interface/HermiteGalsState.h) 与 [Publication](../../include/Modules/Interface/HermiteGalsPublication.h) 中。

`HermiteGalsPrimaryState.transport = HermiteGalsTransport::FirstJet` 显式声明非守恒模式。该状态只保留 first-jet skeleton 和 geometry support，`phaseVolumes` 必须为空；填入零或重算体积都会被拒绝。`ConservativeVolume` 是历史接口的兼容默认值；调用方必须明确选用当前纯路线。

```cpp
const auto mixed = closeHermiteGalsMixedAmr(mesh, candidate);
HermiteGalsCenterPrediction centers{candidate.identity(), tracedCenterValues};
auto geometry = publishHermiteGalsFirstJetGeometry(
    mesh, std::move(candidate), mixed, centers, config);
```

`tracedCenterValues` 按 geometryCells 排序，由同一旧 snapshot 与实际回溯得到。模块不拟合这些值来满足体积；返回 closure 中的中心值与输入逐元素相等。mixed result 与 center prediction 必须属于 candidate 的同一 revision，布局/有限性及 mesh/time 校验保持。失败不消费 candidate，成功保留 identity/time 并消费 candidate。

`geometry.state().geometry()`、`volumes()`、`facePartitions()`、topology、normal 和 curvature 共享同一 scalar。几何测量继续保留真实 XY 面积或 RZ 的 2πr 体积及有效性检查。`volumes().report().totalRawLiquid` 是当前几何体积诊断；不得使用 `initializeCandidate()` 将其注入纯状态。

此模式没有独立库存，因此 publication 的 `maxRelativeInventoryResidual` 及 measurement 的 `totalSourceInventory`、`geometricInventoryDifference`、`maxCellRelativeInventoryDifference` 为 NaN（不适用）。控制台输出 `inventory=NA,flux_ledger=NA`，并记录 `V0,Vn,(Vn-V0)/V0`。它不表示零守恒误差。

当前模块验证入口：

- `quantum_interface_hermite_gals_benchmark --first-jet-contract`
- `quantum_interface_hermite_gals_benchmark --first-jet-smoke`
- `quantum_interface_hermite_gals_benchmark --first-jet-full`

尚未注册新的 production App 字典选项。纯界面通过完整精度 gate 后，仍需另行完成 Hydro 耦合契约。
