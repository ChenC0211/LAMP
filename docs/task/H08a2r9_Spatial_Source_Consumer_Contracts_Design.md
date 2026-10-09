# H08a2r9 / 消费者契约迁移实施前设计

2026-10-05。上一轮为 progress；H08a2r8 的 74 payload、14 当前文件和 1382 源码/输入已复核。本任务检查归 `runs/checks/H08a2r9/20261005_155951_186`。production 0，test ≤2，仅原 phase/momentum contract。

断言迁移清单：

| 旧断言/参考 | 保留要求与拟替代检查 |
| --- | --- |
| phase 用单点中心差分要求 B-div=DI/V | 原 208 basis、缩放/叠加/零场不减；对同一实际 J 做足阶分片物理积分，与独立 incidence 和原 DI 比较 |
| momentum 用独立再拟合旧 tangent/trace/lift 预测同一旧公式 | 独立输入解析均值保留；值/完整 atom/物理 GL3 一致，加原 affine/RZ 闭式准确性，不能只自比较 |
| dual flux=Σ overlap.volume×cell 平均散度 | 与同一实际场在真实 overlap 的物理散度积分比较；旧交换差另列诊断 |
| jet 用 volume-coordinate 二次差分、要求点态常散度 | 实际 piece 中对物理加权速度做完整次数求导；RZ 无轴除零，积分散度仍与原 DI 相容 |
| 原切向跳跃必须大于阈值、旧 rank/存储具体常数 | 原源仍诊断跳跃和 rank；不将“必须不连续”或一次拟合计数当新源要求；实际信息、局部支撑/成本显式检查 |

关键前置：不能把 `tangentUnidentified==0` 的原 affine/RZ 精确族悄悄换成更窄的“完整 P2 可辨识”子集。先以原 `AdvectorReference(kind=0)` 的输入、原 64 布局及原闭式值/J 核验已经完成的正则内核。当前接受场作为同输入对照；新内核仍不发布 identity/clock/材料状态。原值/解析 J 门槛 1e-12 保持，n32 的 J 只补充观察（原 jet 矩阵到 n16）。

前置通过后才实现通用物理积分/求导消费者，并对当前源与正则暂存数据执行同样数学判据和负向扰动。若原精确族出现违反，保留完整矩阵和最坏实际系数见证，任务记 FAIL（精度先决），退役临时入口并恢复原测试；不能删断言掩盖失败，也不能在本任务改生产拟合。

临时前置模式在原 momentum 契约中构建并捕获源码/runtime，结束时移除。正式交付仍只有原标准 CMake/CTest 入口；未经实际来源切换，不重跑完整 H10b 或声称 H08a/H11 通过。

执行结果：64 布局前置完成，6 个原适用布局失败，见[报告](H08a2r9_Spatial_Source_Consumer_Contracts_Report.md)。按上述预定失败规则，未实施断言迁移；原测试恢复，生产未改。
