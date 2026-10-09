# I12a5 / 同源 Q4 局部 Green 势报告

**PASS。** 新 quartic 重载在同一 retained rational 曲线位置内先平移坐标，再计算 Green 势。它为薄扫掠区域测量提供数值基础；尚未实现子区域液量或一般时间积分，I12 整体仍 FAIL。

## 实现

3 个 production 文件、2 个测试相关文件，符合 [TODO](I12a5_Shifted_Green_Moments_TODO.md)。新增有限物理原点参数，XY 使用 x−x0、y−y0，RZ 使用 pi(r−r0)(r+r0)、2pi r(z−z0)。原零原点入口和 nativeMoments 保留旧运算路径；旧根、graph、分区、volume sensitivity 和控制容限未改。见 [API](../api/I12a5_Shifted_Green_Moments_API.md)。

## 验证

- 主程序和 23 个相关测试目标 Release 构建通过，288.73 秒，无 warning。
- **22/22 CTest 通过**，111.23 秒；其中 **21 项已有测试的完整输出逐字一致**，含 architecture 的 12 targets / 379 translation units / 897 files / 3079 includes / 5 legacy edges。
- 新合同 **145 组检查、8 类拒绝**：独立多项式原函数覆盖 XY/RZ、两个 graph 方向、正反相、GL8/受控路径、多个域内外原点、同源分区加和、逐字段零原点一致及全部 formal/independent volume sensitivity。固定边界使用同一势后，轴向和径向两种 Green 公式均给出独立解析体积。
- 16 组大坐标检查将界面偏移置于 2^30 附近且小于全局 double 间距；先形成全局坐标会完全丢失该偏移，新路径保持非零局部矩。包含正反相、两方向、XY/RZ 和分区。
- 非有限原点、异源 cut、空 cut 拒绝；失败后同一源仍可正确测量。源数据没有写入口。
- 81 个冻结/集成参考文件、3 PDF 与前序 I12a4 archive 不变。旧 crossing probe 仍退出 1，完整输出不变；未重跑完整长程 benchmark，也未宣布 I12 通过。

初次测试编译的聚合初始化括号错误，以及分区测试未供应切分后的完整外边界、无轴环域误启用 r=0 symmetry 的 fixture 错误均已修正；生产数值规则没有为这些测试改变，原始日志保留。

证据与逐文件哈希：`runs/baselines/I12a5_20261001/`；运行 `quantum_interface_shifted_green_contract.exe` 可复现新合同。下一项为同源 native Q4 矩形子区域液量测量，之后继续时空积分与 driver 接线。
