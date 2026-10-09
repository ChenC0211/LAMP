# Interface 当前结构

`HermiteGalsState` → trace / vertex、edge transport → mixed/AMR closure → `HermiteGalsPublication` 是共用基础。当前 FirstJet 时间步在 `HermiteGalsFirstJetStep`，网格事务在 `HermiteGalsMeshTransfer`。

`InterfaceCellMeasure` 和 `InterfacePhaseGeometry` 是同源只读几何测量/分区服务；FirstJet 的体积是诊断。`Detail` 为内部算术/闭合辅助，不拥有第二份几何。

`Compatibility` 维持旧 Hydro/EHD 必需接口，独立 CMake owner 为 `lamp_interface_compatibility`。当前 FirstJet App 不调用它，新代码不得增加依赖。旧独立路线和专属 App/测试/工具已退役归档。GALS 仍有历史守恒研究 API，其存在不表示当前入口执行库存输运；移除前需检查现有算子消费者。

当前 App 仅保留入口、严格输入、FirstJet 驱动与同源输出。完整研究在 cases/benchmarks，短契约在 tests/contracts/interface。
