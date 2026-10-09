# H08a2r10 / 缺秩准备保留可辨识低阶族

状态：**PASS（准备层精度修复）**，见[报告](H08a2r10_Identifiable_Lower_Space_Preparation_Report.md)及[设计](H08a2r10_Identifiable_Lower_Space_Preparation_Design.md)。下一唯一任务返回[H08a2r9](H08a2r9_Spatial_Source_Consumer_Contracts_TODO.md)消费者迁移。输入为 [H08a2r9 精度前置 FAIL](H08a2r9_Spatial_Source_Consumer_Contracts_Report.md)：原 64 布局中，6 个原解析门槛适用布局在正则准备的缺秩最小范数解下失败。此项不切换实际 source。

目标：只解决同一正则准备层在信息不足时的再现空间选择，保留原可辨识 XY affine/RZ 允许解析族，同时保持可辨识布局的完整 P2 准确性与原 C1/真实 native 面矩结构。

允许修改：production ≤2，`CchRegularAdvectorGeometry.h/.cpp`；test ≤1，原 `cch_momentum_geometry_contract.cpp`。实施前给出缺秩子空间、投影/分级求解及真实面均值矩的设计，不增加第二 owner、生产算法开关、长期实验脚本或算例专用判断。若确需其他文件先拆项。

禁止修改：实际 B/J/atom 部署、characteristic、FirstJet、phase/material/dual/KEP、最终 native convection、H10b 数据和原 gate；不迁移 H08a2r9 的消费者断言，不调整原门槛或精确性布局条件。

数学约束：材料仅四硬约束，几何独立。零信息的二次分量不得污染已声明可再现空间。XY 保留两分量 `{1,x,y}`；RZ 必须保留原允许族 `u_r=a r, u_z=b+c r²+d z`，并满足实际轴奇偶及指定 wall 模式。完整物理 P2 可辨识时仍须 P2 精确。RZ 两列域的 r 与 r² 不能同时任意辨识；不要额外强制这个不可能交集，缺信息须有 rank/deficit/null 证书。若原允许族本身缺信息，继续明确报告，不伪造精确性。

实现要求：在真实 face-mean 矩上以可辨识性定义统一的层级/最小改动准备，优先空间及余空间均需独立数值依据；支撑按同一邻接扩展规则，不加 half/disconnected/AMR 特例。保留 master/hanging 共用数据、边界绑定、紧支撑面矩修正及同一解析消费。拟合仅发生于 prepare；apply/query 不新增 SVD、分支拟合或速度状态。给出稳定性、权重范数/规模和代价，不能以任意大权重或扩大到全域隐藏问题。

验收：T0 标准 Release 构建、架构/目录及受影响契约，预算外文件/旧实际来源输出保持。T1 原 64 布局预检不缩小：56 个原适用布局的值/native 中心/J 原 1e-12 门槛通过，原 8 个不适用仍明确；n32 的 J 继续补充观察。六个冻结见证必须直接重算通过，并以解析族基向量及线性组合验证再现关系，不能只拟合一个混合输入。T2 重跑 H08a2r8 的 144 结构组/288 非零更新：独立 Bernstein 值/J/完整 atom、C1、真实面矩/Gauss、tiny、axis/wall、拒绝测试，以及原完整 P2 可辨识布局精度全部保持；对原 n8/16/32 的 source/运动可行性族检查不因缺秩处理损失阶数。不得把非恒定点散度恢复成硬约束，实际来源尚未切换时不将运动原型测试称为生产桥接。

失败规则：任何原精度/覆盖退化、不可辨识交集硬塞、边界或来源破坏、隐含 fallback、超预算即 FAIL，退役候选并记录证据；真实外部依赖才可 BLOCKED。

交付：准备层改动与数学依据、原六个反例的关闭证据、结构/精度/成本及回归。通过后下一项返回 H08a2r9 消费者迁移，再同步部署实际来源/运动证书，重验原 H10o/H10b-v2，返回 H08a3/H08a。H11 未启动。
