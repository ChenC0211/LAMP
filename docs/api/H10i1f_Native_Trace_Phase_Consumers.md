# H10i1f / 原生 trace 的正式几何消费

`CchPhaseGeometry::build(mesh, metrics, publication)` 使用 publication 原 `traceHalfFace`；普通 half 仍消费原 ordinary ranges。每个真实 native/AMR subface 只构造一份物理测度，不合并子面，不复制 Snapshot。

`faceTraces()` 返回只读 `CchPhaseFaceTrace`，与原 velocity DOF 一一对应：

- `kind`：Ordinary、InteriorCoincident、BoundaryCoincident。
- `negativeCell/positiveCell`：Hydro **所选 domain** 的真实 incident cell；负值表示该侧不可用。publication 有额外活动支撑也不能增加 selected 内侧。
- `minus/plus`：各可用侧完整 liquid/gas 物理空间迹；不可用侧为零，不能解释成外部相数据。
- `coincidentHalfFaces`：含正长度 zero range 的实际 half 数量。

`faceMeasures()` 是空间测度：普通面保留严格相区间；重合部分在内部用原两侧空间迹平均，在边界仅用实际内侧。存在普通与重合部分时分别积累。XY 用长度，RZ 用真实 `2*pi*r` 测度；原 128 eps 面容量 gate 保持。报告增加 `coincidentFaces/boundaryCoincidentFaces`。

这里的两侧平均是明确定义的无速度诊断，并非 zero face 的相分配或因果材料通量。原 `CchPhaseRate::evaluate` 对任意重合 face 的**非零实际 normal velocity**明确拒绝，包括极小非零值；精确零速度给出零率。没有默认外部相、容差零速度或旧入口 fallback。普通 face 常值率、原来源/完整时间/轴 gate 不变。H10j 的新 causal owner 单独设计与验收。

`HermiteGalsPhaseRegionSource` 按原 transport 选择 topology：FirstJet 使用 complete native trace discovery；historical ConservativeVolume 保留旧入口。`nativeSideGraphs()` 给出原 complete topology 的 side graph 数量。每个 side graph 的实际可用 incident leaf 必须有同源完整纯内部相证书，graph 保留连通性与 provenance，但不进入旧 regular arc integrator。

实际 rectangle、dual volume 和 pressure/dual overlap 的容量按真实交集计算：完整纯相证书直接给出 capacity/0；普通 cut 仍调用原 typed clip/root/controlled Green integration。外部 native side 不再添加一份区域内部体积。原严格非负/capacity、coverage、误差预算与失败事务保持。

prepared region source 通过 owning geometry view alias 保留原 Snapshot、identity、mesh 和 full time；`source()` 可核验完整来源。`CchMomentumMass` 本轮没有改算法，仍直接消费 publication state 和 RegionSource，质量不是 pressure cell 密度平均或材料通量更新。

支持范围由 [报告](../task/H10i1f_Native_Trace_Phase_Consumers_Report.md)限定。原 critical/allzero/high-order normal 与未认证 RZ 轴/域角点继续拒绝。没有新库存、bubble、体积拟合、质量修复、source perturbation、material F、R/uR 或 convection 系数。

材料仅四类硬约束：唯一共享 F、内部反号/真实 AMR 子面和、来源一致、KEP/skew。quarter/精确交换/R=0 为诊断；本接口的几何 PASS 不构成材料 PDE、原 H10b 或完整流体/EHD PASS。
