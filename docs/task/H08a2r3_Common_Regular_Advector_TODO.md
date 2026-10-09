# H08a2r3 / 共用正则advector的可行性

状态：**PASS（构造可行性，未部署）**，见[报告](H08a2r3_Common_Regular_Advector_Report.md)。144结构组、原32族96窗口和sin六窗口通过；临时入口退役，1379源码/输入恢复，4/4回归及2024原确定性行保持。下一[H08a2r4](H08a2r4_Regular_Advector_Preparation_TODO.md)实现同一owner的几何准备。由[H08a2r2失败](H08a2r2_Skeleton_Motion_Selection_Report.md)触发；[实施前候选设计](H08a2r3_Common_Regular_Advector_Design_Change.md)及其实施细化固定本轮预算。生产当前保持H08a2r1，H08a3原消费目标保留，H11未启动。

目标：只验证“同一个C1高阶速度来源，兼容原native面均值，并为完整FirstJet骨架提供同源二阶梯度”能否在实际XY/RZ/AMR与边界上成立。不同时部署材料时间积分或完整CCH App。

允许修改：先在`runs/checks/H08a2r3/<timestamp>`准备审计；最多临时扩展现有`tests/contracts/hydrodynamics/cch_momentum_geometry_contract.cpp`与`cch_characteristic_contract.cpp`两个测试，标准CMake targets构建。**本可行性任务0 production、≤2 test**，关闭时审计完整归档并恢复。先登记旧atom次数、cellDivergence语义和所有消费者；可行后另拆唯一owner部署任务，不能以预算为由留下永久测试速度替代来源。

禁止修改：Interface primary/geometry、库存/bubble/体积拟合、材料修复、四项材料硬约束、压力/黏性/电场/App/字典。不得将解析参考写入重构权重、选择或接收状态。

数学约束：唯一共享F、内部反号/真实AMR子面和、实际几何/advector/identity/fullclock来源一致、KEP/skew与真实R/uR。候选field及其导数必须同源，所有面积分也来自该field。保留native面均值及cell积分散度；明确放弃一般输入下“全局C1且逐cell点态常散度”的过约束，不改D_I或C=−D_I^T。不得用新速度追踪、旧B积分材料。

实现要求：依设计验证共享master的局部P2正则base和实际face cardinal均值校正；C1/可辨识P2/边界/轴/AMR不能只靠样点宣称成立。rank、支持扩展、系数范数、atom次数与复杂度都显式记录。极小selected/disconnected支持的不可辨识矩独立给证书，保留统一最小范数结果及全部C1/面矩检查，不把它计作P2精度通过。同一新field的轨迹与J需使用同一离散运动，并区分与旧B的数值变更。纯零场/仿射成功不能替代非仿射骨架阶。

验收：执行设计中固定的T0/T1/T2。T0冻结、标准Release和实际source/exe绑定；T1全部实际连接、native面矩、Gauss、零场、可辨识P2/RZ正则及selected边界scaled缺陷≤1e−11；T2原32族/96窗口/完整实际DOF的source及有限时间phi/grad收敛、同一map差分与时间误差分离，含原sin与二次剪切反例。value/derivative候选空间阶≥2.5/1.8，有限时间phi/grad≥1.8；任何阈值/缩放/fixture细化须在数值运行前记录，不能事后放宽。

失败规则：同源、矩相容、C1、规定精度、边界或可完成成本任一失败判FAIL；未声明/无证书的rank缺口不能作为通过。保存全部行和反例，恢复临时改动。真正外部依赖阻断才BLOCKED。不得恢复原选支补丁并冒称目标完成，也不得暗增压力到dual精确commuting或R=0约束。

交付：候选数学构造、实际矩阵/轨迹证据、source/exe/恢复冻结、精度与成本及明确未支持范围。PASS仅表示构造可行；后续部署必须按owner预算拆分，并重验H10o最终native算子/原affine/H10b-v2完整28窗口84项，再解决H08a3全DOF消费及H08a时空材料积分。H11保持在这些前置之后。
