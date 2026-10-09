# I12a11r17 / 几何精度预算下的局部 bubble 闭合

目标：只改变局部 I10 的 cut cell 接受关系：在显式几何误差与 Hessian 扰动预算下闭合，不再暗中追求库存机器精度。依据[已授权设计变更](I12_Geometry_Inventory_Accuracy_Design_Change.md)。

允许修改：`include/Modules/Interface/HermiteGalsBubble.h`、`src/Modules/Interface/HermiteGalsBubble.cpp`（2 production）；新增 `tests/interface_accuracy_bubble_contract.cpp`、`tests/CMakeLists.txt`（2 test-related）；TODO/Report/API/README 与本项证据。

禁止修改：Core 数值核、库存/通量/reflux、完整驱动、I11、冻结文件/PDF、原 I12 benchmark 门槛、其它模块。

数学约束：保留原 first jets/mixed/M/identity/time，只有 cell center 派生值可改。新契约显式采用原 `relativeVolumeTolerance` 为 cut 几何预算；pure 仍保留原 roundoff 预算及拓扑规则。逐 cell 输入有限非负物理 Hessian 分量变化上限；按设计中的 bubble 全单元界换算 center 区间并与原 topology/displacement 证书取交集。根、区间扩展和最终实际存储值均不得越过该界。目标在精度区间内不可达则拒绝，不通过改 M、放大预算或造气泡接受。

实现要求：在同一 I10 内核中增加显式契约选择与具名逐 cell 精度预算；历史默认保留为显式 A/B，完整驱动本项不切换，无隐式 fallback。新契约的 inactive cut 保留预测 coefficient；active cut 使用同一 safeguarded root、同源求积和 topology 检查。Report 增加契约、逐 cell 和全局 Hessian 扰动界，Report 不驱动分支。

验收：T0 Release 无 warning，architecture/help 与相邻 I10/I11/I07/平移 smoke/库存/reflux 回归通过，原数值输出不变。T1 弱 corner cut 的 128 epsilon 目标差不再触发破坏几何的修正；XY/RZ、uniform/AMR 的实际修正满足独立采样 Hessian/C1 检查与几何预算；零/过小预算、错误顺序/缺失/非有限预算拒绝；M/identity 不变。原保守误差和工作预算不修改。

失败规则：断言/预算/回归失败即 FAIL，保留原日志；具体缺失接口才 BLOCKED。该算子 PASS 不代表完整新推进已通过；共享通量有界性协调是后续切换前提。

交付：最小源码变更、契约/API、实际测试与原输出比较、冻结/hash 审计、报告和 evidence archive。
