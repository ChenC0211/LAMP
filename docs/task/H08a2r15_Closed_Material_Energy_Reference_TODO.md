# H08a2r15 / 闭域真实材料动能参考

状态：**PASS（消费者迁移）**，见[报告](H08a2r15_Closed_Material_Energy_Reference_Report.md)。下一唯一任务为[r16非线性运动可信参考](H08a2r16_Resolved_Nonlinear_Motion_Reference_TODO.md)。前置 [r14](H08a2r14_Identifiable_Constant_Subspace_Report.md)已关闭缺秩常数与原 dual affine；[r11 恢复失败](H08a2r11_Common_Regular_Source_Deployment_Resume_Report.md)中原四个 closed streamfunction 的零功参考已由本项迁移。实际 production 来源仍旧，不能声称正则部署已通过。

目标：用原唯一 F 的真实 dual 出流和三层质量独立构造 R，检查现有 R/uR 与 KEP/skew 动能身份。解除“native D_I=0 必然给出所有 dual R=0/对流零功”的消费者过约束，不修改方法去迎合旧零功断言。

预算：**production 0，test≤1**，只改既有 `tests/contracts/hydrodynamics/cch_kep_convection_contract.cpp`。原 streamfunction、n4 XY/RZ × uniform/central-AMR、同一几何/材料/边界及 transported probe 全部保持；不新增长期 target、owner、运行时来源开关或第二 F。先写独立参考设计再实施。

在原静态纯材料闭域 fixture 中，逐真实共享 segment 以 dual 侧符号求 q，独立核验质量差为零、R=q、外边界质量率严格零。按当前 `momentumRates` 记号，对流功应满足 `sum(u_i C_i)=0.5 sum(q_i u_i²)`；不能预设右边为零。继续从实际 `momentumLoads=-C+uR`、质量变化和材料动能项检查完整能量身份。其他具有边界功、变质量的原测试仍验原通用身份，参考不得从 production report 回填。

T0：标准构建/9 项回归、预算外文件逐字保持。T1：全部原 KEP 矩阵及四个 closed 输入，原 1e−12 门槛、const/affine/native moment/密度比/来源与 clock/拒绝不减。加入遗漏或错号 uR、错误假设 q=0 的有效负控，同时保留原 wall/轴/共享零通量断言。旧实际来源满足零功的原数值记录应保持，不为正则来源单设容差或公式。

T2：隔离冻结 r11 四文件候选，链接 r14 准备库，以同一新消费者复跑完整原 KEP 契约；保持实际原四输入，确认非零 q/R 被正确计入，不能清零、补偿 F 或改材料质量。反例/失败即如实记录，禁止以 native D_I 代替 dual q、删除难例或降低门槛。

交付报告说明旧约束为何不成立、同一 F 的独立动能预算、覆盖/负控/代价/来源绑定；检查和临时 probe 归 `runs/checks`，关闭时归档，不部署正则来源。本项后先补齐 r11 非线性跨分片的可信 point/J 参考，再恢复同步部署及原完整 H10b-v2；H11 未启动。Pure FirstJet 独立，无库存/bubble；材料只保留唯一 F、内部反号/AMR 子面守恒、同源及 KEP/skew 四硬约束。
