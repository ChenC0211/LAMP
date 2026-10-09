# H08a2r8 / 正则场共同解析求值与 atom 内核

状态：**PASS（共同解析内核）**，见[报告](H08a2r8_Regular_Advector_Evaluation_Report.md)及[实施前设计](H08a2r8_Regular_Advector_Evaluation_Design.md)。144组、9/9回归及原9849行保持；下一[H08a2r9](H08a2r9_Spatial_Source_Consumer_Contracts_TODO.md)迁移空间消费者契约。输入为 [H08a2r4 准备层](H08a2r4_Regular_Advector_Preparation_Report.md)及 [H08a2r7 完整表示](H08a2r7_Complete_Normal_Polynomial_Report.md)。当前实际 B 仍旧，现有正则数据仅含 Q4 tensor DOFs 与真实面矩修正幅值，尚无生产解析求值/完整 atom 提取实现。

目标：用这份既有准备数据和系数实现同一个正则场的物理值、解析 Jacobian、分片完整法向多项式，为原 owner 的原子式来源切换提供可核验的共同内核。不得分别拟合运动、相率、材料速度。

允许修改：production ≤5，限 `Modules/Hydrodynamics` 的既有正则 Detail 实现/声明及必要的同模块内核；test ≤2，优先扩展既有 momentum/phase 契约。准确文件、分片划分及调用关系写入实施前设计。无新 App/长期审计脚本，几何属于原 reconstruction。

禁止修改：本项不提前只切换某个实际 B/atom/J 消费者，不加入第二个可发布速度 owner、运行时双路线开关或 fallback；不动 FirstJet、材料 F/R/uR、KEP、字典及历史证据，不恢复库存/bubble/R=0。后续生产切换必须同步覆盖这些查询和运动证书。

数学约束：与 H08a2r3/r4 的同一 P2→共享 Q4→真实面矩紧支撑修正一致；修正保持 C1、轴正则及选定墙面法向零，native 面均值和 AMR 实际子面和保持。只保留物理积分散度/Gauss 与原 D_I 的一致，不要求点态散度 cell 常量或 pressure→dual 交换。值/导数/atom 必须解析来自同一系数和分片支撑，不能采样后拟合成另一场。

实现要求：逐项使用 `cellDofs`、`faceCorrections`、`supportKnots` 和真实面邻接/支撑宽度。Q4 tensor 与面修正的值/导数解析可复用，atom 同时提供独立 begin/end chart 的全部四次项，含 RZ 测度所需坐标。支撑边缘、cell/skeleton、AMR coarse/fine 交点按同一个闭包评价；缺少支撑或身份不符明确拒绝。`apply` 与查询不得增加 SVD/动态拟合，准备可复用；只添加正式下一来源必要的共同数学内核，不把测试参考提升为 owner。

验收：T0 标准构建、架构/目录及受影响原契约，当前未切换消费者的原 gate 保持。T1 在原 XY/RZ、uniform/AMR、selected/disconnected、轴、绑定/未绑定墙面范围，对独立 tensor/闭式导数及物理积分验证值/J/完整 atom 的一致；内部两侧 C1、真实 native 面矩、Gauss 积分与原 D_I，精度预算沿用 r3/r4，不用点态常散度旧约束替代。T2 P2 可识别再现和缺秩证书原样保留；普通高次、非多项式 native 输入及极短 atom 积分明确覆盖 c3/c4，记录查询/更新成本和 storage。原 96 运动窗口及 H10o/H10b-v2 在实际来源切换后必跑，不能用局部内核通过代替。

失败规则：任何 B/J/atom 来源分裂、C1/面矩/轴不符、高次截断或退化为重新拟合即 FAIL；修正或退役候选，不放宽门槛。仅真实外部依赖可 BLOCKED。

交付：共同解析内核、独立数值与来源证据、未切换范围和下一原子式部署清单。部署必须同时迁移实际 owner 查询、piece/atom overlay、完整次数 C1 运动证书及相关来源参考；依预算分清准备/切换/完整 benchmark 验收，不得将长期保持旧场当作完成。最终仍返回 H08a3/H08a，H11 未启动。
