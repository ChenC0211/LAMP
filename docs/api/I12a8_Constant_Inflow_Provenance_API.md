# I12a8 / 常相入流历史与区间身份 API

`HermiteGalsConstantInflowSource(mesh, baseExtent, entries)` 绑定一份没有 BC 的真实 base interval，并复制完整、有序的 actual inflow face 条目。每个条目为 nativeFace 与 `HermiteGalsInflowPhase::Gas/Liquid`。相态在整个面和该时间范围内为常量，由调用者显式给定，不从旧 phi 推断。

支持 owned constant、XY 单轴/RZ 轴向速度，完整 active 原域、无 opaque coverage。当前速度下的所有 actual inflow faces 必须各出现一次，出流/零法向/内部面不能混入；零速度使用空列表。此接口不包含空间相分区、时间插值或任意 callback。

source 的副本共享 immutable Data；独立创建的相同数值条目不具备相同身份。`baseExtent()`、`entries()` 是只读借用；`phaseFor(nativeFace)` 只接受已有实际入流面。`validateInterval(mesh, interval)` 检查 mesh、旧 primary、velocity/characteristic owner 与时间覆盖。baseExtent 禁止已有 BC，所有权因此无循环。

`HermiteGalsPhaseFluxIntervalSource` 新增末参数 `const HermiteGalsConstantInflowSource&` 的构造重载；`inflowSource()` 返回随 interval 保持生命期的只读来源指针，无 BC 时为 nullptr。原重载保持无 BC 语义。sameInterval 与 transfer.appended 均包含 BC owner；组合后的 interval 保持同一 source。I12a3 原有的 sameInterval 检查会拒绝整个面集合中的 BC 混源，I09 继续持有同一最终结果。

`integrateHermiteGalsTranslationFace(mesh, interval, nativeFace, regions, config)` 消费完整单轴平移面：实际入流面必须有显式历史，其 liquid 是 0 或 total；其余面使用 I12a7 扫掠。结果含 transfer、prescribedInflow，以及仅对扫掠路径提供的可选 swept report。无 BC 入流、内部扫掠出域、未支持空间/时间源均拒绝。

`LevelSetHermiteGeometry::constantNormalFaceTransferVolume` 为边界与 I12a7 内部 total 共享的 strict Core 度量：对原 double 速度/时间/端点作精确代数，最后得到 XY 面积或含 pi 的 RZ 体积浮点值。它不计算相指示、也不外推 scalar。入流与内部纯相使用同一个时间/几何计算次序。

新入口只形成共享通量；M 仍由 I12a3 独立装配，边界量不参与库存舍入修复。未接 first-jet/完整 driver，不能把纯相面集合合同当作完整推进证明。

任务与验收状态见 [TODO](../task/I12a8_Constant_Inflow_Provenance_TODO.md) / [报告](../task/I12a8_Constant_Inflow_Provenance_Report.md)。
