# H08a2r7 / 完整法向多项式表示与相率转换

状态：**PASS（完整表示）**，见[报告](H08a2r7_Complete_Normal_Polynomial_Report.md)及[实施前设计](H08a2r7_Complete_Normal_Polynomial_Design.md)。1632高精度复算、9/9回归及原输出保持；下一[H08a2r8](H08a2r8_Regular_Advector_Evaluation_TODO.md)实现正则共同解析内核。输入为 [H08a2r6](H08a2r6_Benchmark_Polynomial_Schema_Report.md)：材料/dual 参考及三处 benchmark schema 已准备好，实际 B/atom/J 仍是原来源。

目标：使同一个 `CchNormalPolynomial` 完整表示和积分四次法向速度，RZ 包括线性半径后的五次被积项；关闭 phase rate 与 phase/momentum 参考的次数截断。该表示准备不等于正则速度已部署。

允许修改：production ≤5，限 `CchSpatialAdvector.h/.cpp`、`CchPhaseGeometry.cpp` 及必要的同模块 Detail 声明；test ≤2，限现有 `cch_phase_geometry_contract.cpp`、`cch_momentum_geometry_contract.cpp`。不增加长期驱动或 case 开关。实施前列出准确文件及转换入口。

禁止修改：实际 B/J 重构、characteristic 选支/证书、FirstJet、接受状态、F/R/uR/KEP、网格与字典、原 gate/容差和历史证据。不得将合成高次速度发布为实际 owner，也不在本项顺便切换正则 B。

数学约束：材料仍仅唯一共享 F、内部反号/AMR 子面和、来源一致、KEP/skew 四硬约束；几何独立，无库存/bubble，不强制 R=0。begin/end chart 的速度与径向测度方向一致，显式独立 width 不从 rounded 全局端点差重建。保留低次已有行为，高次为完整贡献而非截断或采样拟合。

实现要求：数组含五项，`value/integral` 消费完整库存；局部移位后积分以保护极短区间。`ratePolynomial` 按完整系数转换为现有认证多项式项，保留各长度/测度因子，不复制第二套转换。若需直接测试私有转换，可将原函数提取为同模块 Detail 纯函数供生产与原测试共同调用；不得为测试添加运行时数值分支。两份参考使用全部系数及独立闭式或 GL3（四次×线性测度）积分，清点固定三项 work-scale/线性检查/零数组假设。material/dual 和 benchmark 现有通用消费者直接复用。

验收：T0 标准构建、架构/目录及受影响原契约；原二次来源/身份/时钟及数值 gate 保持。T1 0–4 次单项式与混合四次、双端 chart、XY/RZ 两径向方向、正负法向及非零极短区间，对独立闭式积分；验证真正生产 `ratePolynomial` 的高次系数和物理积分，含显式非零 c3/c4，不用复制转换作为参考。T2 phase/momentum 与已升级 material/dual 回归通过；相关现有导出按名称回读完整五项，新增高项在尚未切换的实际二次来源中应为零，原物理列保持。记录实际舍入差，禁止放宽阈值。完整 H10o/H10b-v2 留到实际来源切换后必跑，不能用表示测试替代。

失败规则：高次项丢失、双端/极短区间错误、来源或物理 gate 改变即 FAIL；修复或退役候选，不靠截断恢复。仅真实外部依赖可 BLOCKED。

交付：完整生产表示/原相率转换、两份既有契约升级、逐项独立证据及源码/runtime 冻结。随后同步部署同源正则 B/atom/J 和运动证书，重验原 native/affine/H10b-v2，再返回 H08a3 桥接及 H08a 时空积分；H11 未启动。
