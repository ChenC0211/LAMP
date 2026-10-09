# I12a5 / 同源 Q4 局部 Green 势 API

在 `LevelSetMappedInterface::integrateQuartic(source, graph, cuts, regions, axisymmetric)` 后增加显式 `Vector2D greenOrigin` 的重载。原入口委托零原点；原 arithmetic path 保留。

这是积分势的常数原点，不是新的 scalar 坐标系或几何来源。所有 cuts 仍必须持有完全相同的 Snapshot。原点可在域外（包括负 r0），但两个分量必须有限。RZ 的实际曲线/region 仍须满足原来的非负半径约束。

| 返回 moment | XY | RZ |
| --- | --- | --- |
| axialMoment | ∫(y−y0) n_y ds | ∫2πr(z−z0) n_z ds |
| radialMoment | ∫(x−x0) n_x ds | ∫π(r−r0)(r+r0) n_r ds |

坐标差在 strict Core 的 retained rational physical point 内形成后才近似。调用方的固定边界积分必须使用同一原点；free moment 单独并不具有平移不变性。两套势的散度都对应原 XY 面积/RZ 体积。

原 25 个 formal sensitivity 和 native independentSensitivity 仍表示原液体体积变化，不能解释为单独 free moment 对数据的导数。graph/cut/region 身份和原细化拒绝规则不变。GL8 控制和误差诊断不变；局部原点并未把浮点误差估计变成严格包围，也没有承诺任意薄区域达到库存舍入门槛。

`nativeMoments` 及其端点分部积分公式未增加此能力。尚未提供 rectangle volume、时间积分或 driver 接线；后续任务将组合固定边界和 free moment。

见 [任务](../task/I12a5_Shifted_Green_Moments_TODO.md) 与 [报告](../task/I12a5_Shifted_Green_Moments_Report.md)。
