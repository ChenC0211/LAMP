# I12a15r1 / 合并所有已认证候选分割的积分界

目标：只解决一个数学事实——同一父 face×time 上，切向和时间两种完整候选分割均给出有效积分包围，应保留两者的交集，不能仅因选择另一方向继续细化就丢掉已获得的证书。首先用 a15 原 uniform shear 反例检验，在原精度/工作预算下完成认证。

允许修改：`src/Modules/Interface/HermiteGalsTemporalPhaseFlux.cpp`（1 production）；`tests/interface_bounded_spacetime_flux_contract.cpp`（1 test-related）；本 TODO、设计/API/报告、入口和证据。若证据表明仅此数学改动不能解决反例，本项如实 FAIL 并归档，不在本项扩大算法范围。

禁止修改：owned velocity 数据/实际求值、原 RK2 表达与 branch 守卫、Core 区间/符号引擎、phase/inventory source、M/几何闭合与发布、reflux、原 benchmark/gate、任何既有测试容差和 panel/depth 上限。不得增加 trial 次数或把采样用于相/积分认证。

数学约束：父真实四个 gross 积分同时属于父旧范围、切向完整二分的子范围和、时间完整二分的子范围和；三者交集仍包含同一真值。四项 endpoint liquid<=total 的顺序在交集下保持。所有非零通量的最终树叶仍须认证原完整 RK2 分支；未选候选的条件 [0,1] 估计不能免除 selected tree 的 pending 守卫。

实现要求：每次已有两组候选的计算完成后，先求共同父证书交集，再按原工作启发选保留树分支；把共同证书写入父结点并照旧逐层与子和相交。不改变方向选择所比较的候选自身，不增加评价工作，不把未保留分区伪装成树叶。只存在一种可表示切分时使用该候选原证书。

验收：T0 范围 <=1 production+1 test、无新增 compiler warning、81 frozen/3 PDF/a15 archive 身份不变、全部 37 runtime 重建、原 22 项回归及其 21 个数值输出/架构保持。T1 a15 原 57 个主积分场景、8 组 Core 面矩、864 次实际 B/J、24 guards 全部运行；uniform shear 要求 E<=5.859375e-6、4096 panels/40 depth 不变。增加严格 Core 区间交集端点/代表值包含性的非嵌套证书测试；新算子全部既有包含性、同源性、时间组成、重入、边界拒绝继续满足。评价次数仍 <=1+4*(最终叶数-1)，不加 work budget。验收失败时保留同一原反例，不改为 expected success。

失败规则：交集为空/液总顺序损坏、参考逃出证书、pending branch 被接收、原行为改变、原预算仍耗尽，均 FAIL。仅真正依赖不可获得外部输入时才 BLOCKED。I12 的完整 driver/gate 不因本项通过而降级。

交付：源码、完整 T0/T1、同一反例前后误差与工作量、API/设计/报告和 hash archive。完成归档后再进入一般边界历史/库存可行性与原旋转/涡旋全驱动矩阵。
执行状态：**FAIL**。原剪切在 4096 panels 时 E=6.6763159339198168e-6，超过 5.859375e-6；真实保持普通测试失败，未扩展本项算法范围。共同候选证书微小收紧但不足以达标；完整记录见 [报告](I12a15r1_All_Partition_Certificates_Report.md)。