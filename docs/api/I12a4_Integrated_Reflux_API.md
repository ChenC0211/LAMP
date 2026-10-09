# I12a4 / 区间积分 reflux API

使用新区间入口时，同时包含 `HermiteGalsTemporalPhaseFlux.h` 与 `HermiteGalsReflux.h`；后者对 inventory 来源仅作前置声明，保持旧 I09 公共头的依赖边界。

`assembleHermiteGalsReflux(mesh, integratedInventory)` 重载返回 `HermiteGalsIntegratedRefluxResult`。它保留同一个 I12a3 immutable 结果；`source().faces()/rawFaces()/phaseVolumes()` 与输入共享 owner，不复制库存、不更新 M、不创建几何或重新求积。

`HermiteGalsIntegratedSideTransfer` 按 cell 升序及 W/S/E/N 排列，记录 actual subface 索引、outward sign、四个非负 gross volume 的派生和，以及实际子面 signed total/liquid volume 的顺序和。各量保持全局正 x/正 y 方向；乘 outward sign 才得到单元外向量。signed volume 直接累加源实际面，避免由 rate 或重新分组 gross 差重建守恒交换量。

`averageTotalFlux(sideRow)`、`averageLiquidFlux(sideRow)` 只除以 `source().interval().dt()`，不存平均 rate。非法 side row、不可表示平均 rate、moved-from handle 和陈旧 mesh 明确拒绝。纯积分量可以合法存在而派生 rate 不可表示，这不会迫使库存转换成 rate。

新旧重载共用同一 native incidence、side ordering、2:1 fine-neighbor pairing 和粗细诊断实现。旧类型仍保持 midpoint 语义；新类型始终使用实际 interval，绝不借用某个时间节点 geometry identity。

报告提供 side/coarse/pair 数量及 total/liquid volume 的粗细外向抵消残差。汇总消费的是 I12a3 **最终**共享面量，包含原预算内舍入修复；原始量仍从 `source().rawFaces()` 可追溯。

见 [TODO](../task/I12a4_Integrated_Reflux_TODO.md)。此 API 是 I09 的区间消费，不代表一般曲面时空积分或 I12 完整 driver 已接线。
