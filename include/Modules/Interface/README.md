# Interface 当前结构

`HermiteGalsState` → trace / vertex、edge transport → mixed/AMR closure → `HermiteGalsPublication` 是共用基础。当前 FirstJet 时间步在 `HermiteGalsFirstJetStep`，网格事务在 `HermiteGalsMeshTransfer`。

`InterfaceCellMeasure` 和 `InterfacePhaseGeometry` 是同源只读几何测量/分区服务；FirstJet 的体积是诊断。`Detail` 为内部算术/闭合辅助，不拥有第二份几何。

旧 Hydro/EHD、`Compatibility` 支撑及其独立target已归档退役，不属于当前生产接口。GALS 仍有历史守恒研究 API，其存在不表示当前入口执行库存输运；移除前需检查现有算子消费者。

当前 App 仅保留入口、严格输入、FirstJet 驱动与同源输出。完整研究在 cases/benchmarks，短契约在 tests/contracts/interface。

`GeometricVofPublication` 是新增候选的不可变alpha/真实native-AMR初始化几何owner；实际fan切割、空间trace与Gamma同源，不冒用FirstJet身份。尚无保守转移和时间推进/App接线。见[API](../../../docs/api/V2b1_Geometric_Vof_Publication.md)。生产不保留研究光滑性监测。
