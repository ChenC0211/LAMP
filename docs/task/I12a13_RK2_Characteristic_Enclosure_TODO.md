# I12a13 / 同源 RK2 特征时空包围

目标：只解决 owned polynomial source 下，空间 box × arrival-time interval 的实际 RK2 单步内部分支及 B/J 包围事实，为后续相事件完整覆盖提供可用前提。

允许修改：现有 `Core/Geometry/LevelSetBoundaryTrace.h`（只抽取共用代数表达）、`PolynomialVelocitySlab.h/.cpp`（已有 strict TU）、`HermiteGalsTemporalPhaseFlux.h/.cpp`（来源保留结果），共 5 production 文件。新增 `tests/interface_characteristic_enclosure_contract.cpp` 与 `tests/CMakeLists.txt`，共 2 tests。本文、设计/API/报告、README 与本项证据。

禁止修改：原 RK2 运算顺序、CFL/轴/边界事件/缩步规则；旧时间插值、速度/几何/M、bubble、相积分与 reflux；冻结文件、完整 benchmark 参数或原 gate。不构造第二个实际 tracer，不改变旧 source 行为，不把包围失败变成零通量。

数学约束：arrival time 从同一旧时刻回溯，允许闭时间段包含旧时刻的零时长映射。先包围实际 arrival 查询，再包围同一 RK 中点处的 old/arrival 查询及原 midpoint 平均，传播原 B 与链式 J 表达。只有原完整单步、无物理边界 hit 的分支获得充分覆盖证明时才返回成功；不得把连续 ODE 精确解当 B_h。范围可保守，未证明的 box 拒绝而不假称命中或未命中。

实现要求：先写设计。原 double tracer 和新区间计算共用 midpoint/target/J 代数表达，保持原乘加括号与时间混合。复用已有 outward arithmetic 和 strict 编译对象。检查 closed-domain、实际 sample 的半 cell 限制、RZ 轴、无新 wall contact、opaque coverage；任何外推、非法时间/mesh/source/finite scale 拒绝。Interface 结果私有构造且只读，保留完整 interval/velocity/characteristic/old geometry owner，提供输入、midpoint/departure box 与 J bounds；没有虚构的连续时间 geometry ID。

验收：T0 无新增 compiler warning（依照现行 02 架构规范；完整重建触发的未改旧文件告警须单列、核验文件身份，不能称完整构建零告警），旧 trace/vertex/edge、a12、interval/flux/inventory、r23 smoke/跨面等相关输出不变；冻结/hash/范围审计。T1 XY uniform/AMR、RZ interior/axis-line、rotation/cubic/time-reversal/shear/translation/zero source；实际公共 tracer 的 B/J、单步数与无 hit 均在返回范围内。包含闭时间段、End/Midpoint 视图、零时长和退化 face/point box、独立 affine RK2 解析式、box 收缩。域外、中点/终点未覆盖、CFL、轴、opaque source/coverage、过期 mesh、moved-from、owner/lifetime、非法时间与浮点环境均明确检查；失败与查询不修改旧 geometry/M。

失败规则：旧数值行为改变、实际样本逃逸、错误宣布 interior 或源身份丢失，均 FAIL 并保留数据；不得放宽数学/浮点预算。未能证明的边界/substep box 只允许拒绝，完整 boundary history 另立任务；不因此把父 I12 记为通过或把一般相通量写成已完成。

交付：实际源码、算子与相关回归、来源/冻结/runtime 审计、API/报告及 archive。下一项将包围连接到旧 Hermite 相区间及事件覆盖；本项不同时实施时空 phase quadrature。
