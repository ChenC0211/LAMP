# I12a11r23 / 网格尺度几何误差预算验证

目标：只验证将几何—库存差按截断误差尺度接受，能否在保持 r21 Hessian 限制与独立 M 守恒的同时，通过原非多项式空间矩阵；不把固定求根容限当成网格精度指标。

允许修改：`tests/interface_hermite_gals_benchmark.cpp`、`tests/CMakeLists.txt`（0 production + 最多 2 test-related）；本文、独立设计说明、报告、README 与本项证据。仅在现有配置接口上接线，不复制数值内核。

禁止修改：production 算式、Core 几何/拓扑认证、库存权重/共享面/roundoff/reflux、原 first-jet 及 mixed 路线、冻结文件、原 I12 收敛 gate；不改初态、网格、dt 或终点，不调整原 r21/r22 命令行为。

数学约束：M 仍唯一共享更新且严格有界，geometry 仍唯一。Q4 first-jet 的四阶位置精度目标，对规则界面对应相对 cell 容量误差 O(h³)；该尺度关系不等于全局误差证明，仍须实际通过原位置/梯度/曲率/M 验收。当前现有 API 提供单一相对几何容限，本项使用全网格最小叶宽作为保守统一尺度，不把它宣称为成熟的局部 AMR 误差估计器。

实现要求：先形成设计说明。在单位方域验证配置中预先固定 `tau_V=max(1e-10,(h_min/L)^3)`，`L=max(domainWidth,domainHeight)`，系数 1；实际 leaf/cell 尺寸从当前 mesh 读取。I10 GeometryAccuracy、I11 publication 和 r18 inventory region 必须使用完全相同的 tau_V，原矩量/求积容限不变。每叶 Hessian 预算继续 `min(hx,hy)^2`。新命令显式选择该预算；r21/r22 固定 1e-10 对照保持。记录 h_min、L、tau_V、所有原过程诊断及实际误差，不按结果调节系数或阶次。

验收：T0 原 geometry 和 fixed-budget smoke/跨面/architecture 输出不变，无 warning。T1 新预算的 Q2 smoke 和 48 组平面跨面满足原解析门槛；初态和输运 M 不重置。T2 原 12 组非多项式空间运行全部完成到 .04，三档 n=8/16/32 的梯度/曲率逐档递减，末档阶数 >=2/1.5；原 ledger/reflux <=1e-10。同一 XY uniform n=32 的 dt=.005/.0025 都完成到 .04，和同矩阵 dt=.01 比较；空间截断误差底限与一般运动时间阶要求分开，不伪造时间收敛。

失败规则：任何原空间 gate、实际推进或预算同源性断言失败即 FAIL；完整矩阵及回归保留。不得修改预算重跑成通过；若本项不通过，另立有证据的任务。缺明确接口才 BLOCKED。

交付：最小配置接线、设计依据、所有计划运行、同终点 r8/r22 对照、回归/冻结/hash 审计、报告及 archive。即使本项通过，rotation/vortex、一般时间阶、动态 AMR 及完整 Gate I 仍未因此通过。

完成记录（2026-10-02）：**PASS（本项）**。原空间 12/12 完成、四条梯度/曲率序列通过原 gate；两个 dt 控制完成，7 项回归通过，5 项旧输出不变。详细数值、限制及归档见[报告](I12a11r23_Mesh_Scaled_Geometry_Budget_Report.md)。I12 / Gate I 仍 FAIL。
