# I12p1 / 纯 Hermite–GALS 状态、发布及推进

目标：让实际 first-jet 输运直接发布同源 Q4 几何，不依赖独立库存/相通量或体积闭合。

允许修改：HermiteGalsState.h/.cpp、HermiteGalsPublication.h/.cpp、HermiteGalsPhaseVolume.cpp（5 production）；interface_hermite_gals_benchmark.cpp、tests/CMakeLists.txt（2 tests）；规范/设计/报告/证据。

禁止修改：actual tracer、I05/I06 核函数、Core 几何算法、守恒 flux/bubble/reflux、Hydro/Electricity、冻结 PDF/源码、原历史 case/容限。

数学约束：同一 B_h/J_Bh 与唯一 snapshot；phi_xy 派生、AMR master/hanging 约束不变。FirstJet 模式无 M，禁止几何体积回写成 M。Q4 中心值直接采用 actual trace 预测；无体积 root/重标定/clip/reinitialization。V(phi) 漂移明确为诊断，inventory residual 不适用。

实现要求：先按[设计变更](I12_Pure_Hermite_Gals_Design_Change.md)更新规范。新增显式模式及纯 publication 入口；与旧入口共用同一发布/表面检查内核。candidate 失败事务不变；保留历史守恒入口可回归。benchmark 复用已有 vertex/edge/center trace，不复制数值内核。

验收：T0 <=5 production/2 tests、39 targets、旧 26 CTest、81 freeze/3 PDF/r7 archive 保持；旧模式输出保持。T1 XY/RZ uniform/AMR 纯发布，无 M、中心系数原值、同 scalar、失败回滚、错误 mode/identity/布局/非有限输入拒绝；原 Q2 静止/平移精度保留。T2 新 smoke 和原完整 static 4/translation 40/rotation 30/vortex 400 步；记录几何体积漂移及位置/梯度/曲率误差。完整失败原样报告，后续收敛属于 p2，不将本任务冒充完整 I12。

失败规则：无 M 模式暗调体积/旧结果不明改变/几何多 owner/trace 不一致/新契约失败均 FAIL；数值失败不标外部 BLOCKED。

交付：代码、方法变更/规范、公开 API 验证、原算例结果、审计和归档；一个 task 完成后进入下一项。

结果：**PASS**，见[报告](I12p1_Pure_Hermite_Gals_Publication_Report.md)。8 组发布 / 88 guards、12 完整算例 / 1,036 接受步、28 CTest 及 3 个直接受影响的旧契约通过。体积漂移如实记录；完整 I12 / Gate I 仍需收敛与 AMR 验证。
