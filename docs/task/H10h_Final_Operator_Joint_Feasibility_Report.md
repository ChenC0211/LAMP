# H10h / 最终算子联合可行性与材料输入尺度

2026-10-04，**固定区间候选 FAIL，工程与独立证书核验 PASS；H10 未关闭。** [任务书](H10h_Final_Operator_Joint_Feasibility_TODO.md)、[设计](H10h_Final_Operator_Joint_Feasibility_Design.md)。开始前已核验 H10g 的 158 个 payload、21 份当前交付、1,042 份源文件和 346 份输入。证据冻结于 `runs/baselines/H10h_20261004`，实际观察位于 `runs/checks/H10h/20261004_061552_898`。

只修改既有研究 driver 的被动观察与文档。H10c production/原短契约、Interface、H08/H09、App、字典和冻结 PDF 均未修改。未实现新 QP、修正 target、配对 production API 或新材料通量；未重跑完整 H10b/audit32，未进入 H11/H08a。

## 固定候选为何停止

预先固定最终 K 的常数硬矩、材料一次矩区间、二次余项与行 l1 预算；一次区间为 `B_i=chi*S_i*h_i^3/L^2`，系数 1，`chi=|rho_l-rho_g|/max(rho_l,rho_g)`，`S_i=gross_i+M_i*|V_i|_1/h_i`。真实 V 由原 H07 体积率恢复，不由解析解提供。无调整 H10f 的 8 或本候选的 1。

对无横向周期跳跃的 affine 样本，最终共享 skew 的一次残差必须满足 `sum r=c0`、`sum z*r=c1`。固定区间要求对所有 lambda：

\[
|\lambda_0c_0+\lambda_1c_1|\le\sum_i B_i|\lambda_0+\lambda_1z_i|.
\]

此处不仅检查标量条件：常数行乘子 `alpha_i=2*lambda0*z_i+lambda1*z_i^2` 与归一化一次行乘子 `beta_i=h_i*(lambda0+lambda1*z_i)` 消去**任意共享 pair** 的最终矩阵列，见设计推导。真实内部段/声明端口的逐 pair 数值消去误差为 0；扩大 patch 也不能修复本固定区间的全局不相容。

实际 owner 的 `intervals/intervals32` 合计 48 组、46,680 条 native 行、1,062 个完整二维 zonotope 法向。独立 PowerShell 从原行数据重新求和，无通量/系数/状态求解。最终证书有 **4 组 FAIL**，均为 XY 2:1：

| n | rho_l/rho_g | 最大缺口/支持比 | 判定 |
| --- | --- | --- | --- |
| 8 | 1,000 | 1.38289377268 | 不可行 |
| 8 | 1,000,000 | 1.50189273893 | 不可行 |
| 32 | 1,000 | 1.27922386261 | 不可行 |
| 32 | 1,000,000 | 1.38250840177 | 不可行 |

n8、rho=1e6 的实际 lambda=(-0.46875,1)，缺口 79.2421704846414，支持 52.7615377786830。另 44 组仅通过全局必要条件，不意味着局部 one-ring 可行或 PDE 通过。独立证书与 C++ 原值的最大 gross 归一化差异 1.22e−15；phase/mass/q/R/尺度定义最大差异 2.19e−16。固定候选据此退役，不继续 QP/KKT 实现。

初次临时 PowerShell 预检误用 `Math.Max(0, double)` 的整数重载，漏判 rho=1000 的两组并改变了归一化 violation；已在 `rejected_integer_overload_observer` 保存脚本和原输出，**不作为最终证据**。最终独立脚本显式使用 double，并按实际入/出通量 gross 核验抵消后的 q/R；没有改容差。C++ 观察无该问题。上述四组最终判定同时由实际 owner 和独立复算确认。

## 比 patch 更上游的输入问题

对于非 alias、非固定轴内部行，既有材料一致性机制给出

\[
A=M^{-1}[\operatorname{diag}(q)/2-K],\qquad K^T=-K,
\quad K1=q/2.
\]

因此 `A_ii=q_i/(2*M_i)` 无法由共享 skew 修正改变，且 `||A_i||_1 >= |q_i|/M_i`。`h_i*|q_i|/M_i` 是输入强制的导数行尺度下界，不是某个求解器产生的坏系数。

实际 XY 2:1、n16、rho_l/rho_g=1e6、dof103（component0，坐标 (0.1875,0.53125)）：

| 量 | 实际值 |
| --- | --- |
| 三层液体体积 | 均为 0 |
| midpoint gas/dual 体积与质量（rho_g=1） | 0.00390625 |
| h、质量变化率 | 0.0625、0 |
| 原 advector | (0.4,0) |
| q、R | 均为 12499.9875 |
| `h*abs(q)/M` | 199999.8 |
| 该 box 的纯气体常速度物理外流参考 | 0 |

该行既非边界也非固定轴。rho=1000 的同一行也有 q=12.4875；XY uniform n16 同样出现纯气体尺度 199999.8，RZ 也有对应纯气体大 q，故这不是单独的 AMR 分支问题。纯相参考只是只读真实 box 测度闭合，没有替代 production 通量。

来源检查显示：现 H09 将 primal cut 材料率按面积分配给 quarter 边界，再强制每个 quarter 的外流等于 pressure 材料散度的体积权重。dual 真实质量却由自己的实际相区域积分。体积映射不识别 quarter/dual 的相占据，因而压力区域中的大密度变化可被分配到纯气体 dual 中，形成大 q/R。它代数上满足旧 H09 commuting 契约，却不是局部材料 PDE 准确性的证明。这个输入问题不能继续靠最终 J 的矩约束掩盖。

这证明本固定候选不可行，并定位旧空间材料转移的不足；**不是任意二阶方法不可能的证明，也不能仅由行范数推论全部时间步不稳定**。

## 用户授权与下一任务

用户本轮要求：H08/H09 重设计只以“唯一共享材料通量、内部面反号共享/AMR 子面守恒、来源一致、KEP/skew 结构”为硬约束。已在[耦合设计](H07_Pure_Gals_Hydro_Coupling_Design.md)和[H10i 任务书](H10i_Spatial_Material_Flux_Transfer_Design_TODO.md)固定。旧 quarter 分配和精确 pressure→dual 散度/缺陷交换改为历史诊断；新设计允许非零空间交换差 `E_sp=q_dual-T*q_primal`，真实缺陷分解包含 E_sp。

Pure FirstJet 几何保持完全独立，无库存/bubble、体积拟合或质量修复。精度和 R 收敛仍须验证；H10g/h FAIL 不因解除旧约束改称 PASS。下一单任务 H10i 先设计同一材料 owner 的相分割感知空间求值，之后另行实施并先验证内部 affine，再重跑原 H10b。空间通过之后才 H08a 时空积分；H11 仍等待。

## 工程验证与证据范围

- 标准 Release 构建研究与原契约通过；本轮 CTest **4/4**（KEP 契约、architecture、自检、repository layout）。原契约仍 124 场景、1,356 trace、6,144 load、24 crossing、19 reject、4 regrid；SUMMARY/NATIVE_MOMENTS 与 H10g 完全一致。
- 原 n8/16 affine 回归 **3,760 行**；内部最大误差 1.58492e−13，原 XY AMR n16 dof911 反例误差 7.41629e−14。新材料观察器没有参与生产载荷。
- 同一 48 组的 affine/compatibility/bounds/sources 共 **8 个原 CSV** 与 H10g 逐字一致。去掉本轮被动观察代码后 driver 与 H10g 完全相同；原 full/RK4/donor/gate 数值路径保持。
- 原完整 H10b 最近的 28 窗口/15,116 步/84 gate 仍有 **14 FAIL（10 空间、4 donor）**；本轮未运行，不以观察器 exit0 代替精度 PASS。
- 只变 research driver/README 与文档；原 production、短测试、App、输入和 CCH PDF 哈希保持。manifest 保存平面 SHA256 与独立核验，初始候选和失败证据完整归档。

实际观察可执行 SHA256：`78A000D3381D69C324101AC76B9D84D964A02BD5A31CEDC20EF96C0DD660BB74`；源码 SHA256：`EC89411356AE2ACB2BC10EB5CFC4332DDC0BF556EB516B605C30DA576E88D7A5`。目标仍 active，未宣称完整流体/EHD 达标。
