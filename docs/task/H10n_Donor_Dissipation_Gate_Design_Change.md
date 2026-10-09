# H10n / donor 耗散指标设计变更

2026-10-05，实施前定义。前置 [H10m](H10m_Final_Native_Same_F_Convection_Report.md) 的 12 个实际 donor 窗口与独立 Fourier 模式证明：正确 O(h) 上风耗散在 n16/n32、T=1 的累计损失阶仅约 0.367，累计量的有限时窗饱和与半离散空间阶不能混用。本项改变研究验收指标，不改变数值方法、物理参数、离散性质或任何轨迹。

## 独立能量身份

每个实际 internal segment 保留一份 F；显式周期两端按原 Study.paired 关系组成一条连接，计一次；未配对外段仍保留。令实际 native 外流为 q、现有材料缺陷为 R，原研究质量 M 不变，则 donor 分支的功满足

`P = sum(M_i*u_i*rhs_i) = -D - P_open + C - P_axis`，

`D = 0.5*sum_internal(|F|*(u_i-u_j)^2) + 0.5*sum_periodic(|F|*(u_i-u_j)^2)`，

`C = sum((R_i-q_i/2)*u_i^2)`。

原未配对外端采用 owner trace，故 `P_open=sum(sign*F*u_i^2/2)`；本固定平移研究的这些段为零通量。RZ 固定轴上的真实 u 必须为0，因而轴反力功为0；不得删除这些 native 行或其 half-dual 质量。

观察器从真实共享段及端口独立累加 D、q、P_open，与原 rhs 的 native 功及真实 R 核对。记录 C 和固定材料检查，不把材料项直接丢弃。只在本固定材料/闭合端口场景核验 `-P=D`；一般非守恒材料不因此要求 R=0。身份采用实际功项绝对值之和作 gross 尺度，零尺度要求精确零，门槛仍1e−12。每组初值和原积分终值均核验，并保存逐段/逐native输入供独立重组。

## 新旧 gate 对应

研究验收版本命名 `H10b-dissipation-v2`。四组旧 `*_donor_loss_order` 原公式及 `[0.7,1.4]` 区间保留在 `legacy_gates.csv`，保留 FAIL，但不混入新版本的通过计数。新版本对应四项为 `*_donor_decay_order`，检验原初始样本的 `D/E(0)` 在 n16→32 的空间阶，采用同一事前区间 `[0.7,1.4]`。记录明确版本、映射和上下界，不能称原84gate通过。

其余80项定义保持：空间阶≥1.7、最细L2≤0.05、Linf下降、动能漂移≤1e−6且≤donor实际损失的1%、全部边界/材料/AMR及时间门槛不变。新版本仍有84项；观察器的身份断言是额外检查，不替换任一原中心 gate。实际 donor 累计损失必须仍为正、完整运行同一 T=1 和32n步。

## 实现与验证边界

只改既有研究 driver 的被动诊断/输出和 gate 选择。Study 构造、rhs、exact/error/energy、RK4/integrate，材料准备、F/R/质量、生产算子、四组网格与时间参数保持。`donor-audit` 直接复用原12窗口，原 donor_audit/trajectory/Fourier 数据须与H10m冻结字节一致；新增逐段/native和功表由独立脚本重算。full 与 donor-audit 共用同一 D 诊断和新旧 gate 命名。

本轮重跑12 donor窗口，并以原 n8 XY AMR probe 验证 full 使用的 spatial 路径；不重复已知未通过的完整中心研究，不据此声称中心PDE通过。标准构建、原8项必要回归和最终文档检查必须通过。源码范围0 production/0 test/1 research，无新增target或生产路径。H10n完成后继续最终native有界空间/KE联合设计，内部affine先于完整研究，H11/H08a等待。
