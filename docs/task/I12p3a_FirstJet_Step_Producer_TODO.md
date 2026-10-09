# I12p3a / 可复用纯 first-jet 时间步 producer

目标：将测试驱动中已验证的纯时间步组合归属到 Interface 模块，供正式 App 复用；不改变离散方程。

允许修改：新增 HermiteGalsFirstJetStep.h/.cpp（2 production）；interface_hermite_gals_benchmark.cpp、tests/CMakeLists.txt（2 tests）；任务/API/报告/证据。

禁止修改：Core tracer、I03/I04 pullback、I05/I06 closure、publication/mesh-transfer 算式、旧守恒驱动、字典/App/IO、Hydro/Electricity、case/冻结资料和容限。

数学约束：源为完整 FirstJet accepted publication（Initial/End，无 M），输入 velocity views 的来源/时间/domain 按既有 first-jet 验证；目标为下一步 Midpoint 或 End。调用方提供时钟和真实外边界 first jet，producer 不推进主时钟、不 commit accepted handle。顶点/边用原 transport；中心用同一旧 snapshot 与实际 trace；mixed 用原 I06；唯一几何用原纯 publication。禁止修改中心值来拟合体积。

实现要求：只提取组合。显式 boundary query 仅在目标点位于物理外边界且实际 trace 命中边界时使用，返回目标时间 first jet；不能用它填补内部回溯点/中心的边界命中或重置内部解。缺边界数据、非法来源/时间/域、非法几何均拒绝，源保持。结果保留源 identity、arrival velocity identity、新 publication 及 DOF 计数；arrival velocity 的几何 provenance 不冒充新发布 identity。复用原私有输入校验，不复制 RK2/J、pullback 或 stencil。

测试驱动只准备规定速度、外边界参考、配置和时钟，所有纯 full/收敛/crossing/动态 AMR 的真实时间步调用新 API。删除已替代的测试纯步组合；旧 ConservativeVolume 的历史 predict/finish 保留，不作 fallback。

验收：T0 <=2 production/2 tests，39targets无新warning，原35CTest及完整数值输出保持（架构仅允许新文件/include计数变化），81freeze/3PDF/p2b/p3archive不变。T1 新 producer 的 XY/RZ、uniform/AMR、非多项式、midpoint/end、外边界与非法输入/拒绝事务契约；不添加新的数值容限。T2 原纯 full 和动态 full 各12组/1,036步、原dt/网格/速度/终时，与归档完整数值输出一致。p2a/p2b研究由原CTest覆盖，无解析内部重置。

失败规则：源被修改、时间/provenance丢失、内部边界命中被解析填充、改容限/删失败组、方法自动切换或原数值输出变动均FAIL；内部实现问题不记外部BLOCKED。

交付：模块producer/API、测试接线、原输出比较/范围审计/归档与报告；完成后仅进入I12p3b，不提前声明Gate I或Hydro通过。

结果：**PASS**，见[报告](I12p3a_FirstJet_Step_Producer_Report.md)。16场景/80拒绝、36CTest及2组full（24算例、2,072步、36动态regrid）通过；34旧数值输出与两组full完整输出均逐字保持。下一项I12p3b，整体Gate I未通过。
