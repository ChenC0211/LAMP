# H10m / 同一 F 下的最终算子与独立 donor 门槛复核

2026-10-05，**H10m FAIL / H10 未关闭**。本轮确认了一个不能靠修改中心算子解决的验收冲突：同一 F、质量、R 和原 donor 积分不变时，四项原累计损失阶 gate 均失败；均匀网格的正确上风 Fourier 模式也给出同样的低阶。当前中心算子的空间问题仍在，本轮没有宣称修复或交付新的有界算子。下一唯一 [H10n](H10n_Donor_Dissipation_Gate_Contract_TODO.md)先明确耗散对照指标的数学定义，随后继续最终 native 空间/KE 联合设计；H11/H08a等待。

## 不改变材料硬约束

[H10l2](H10l2_Same_F_Dual_Material_Retry_Report.md) 的同一 F 迁移保持。材料硬约束仅唯一共享 F、内部反号/真实 AMR 子面和、来源一致、KEP/skew；不恢复 quarter、严格 material commuting、R=0、库存或 bubble。真实 R/uR、材料动能项和独立 Pure FirstJet 保持。H10m 没有改生产、测试、输入、Core、Interface、App 或 PDF。

## donor 反例：有限时窗累计量不等于空间耗散率

按[实施前设计](H10m_Final_Native_Same_F_Convection_Design.md)，新增 `donor-audit` 直接调用原 `Study/rhs/integrate`。四组合 n8/16/32、T=1、32n 步、全部 native/half-dual 原样保留；完成 **12 窗口 / 7,168 步 / 396 时序行**，原四 gate 全部 FAIL，进程正常以 exit 1 结束。观察器没有替换原 full、RK4、donor 或门槛。

| 原组合 | n16 累计损失 | n32 累计损失 | 原损失阶，要求≥0.7 | 初始耗散率阶，仅诊断 |
| --- | ---: | ---: | ---: | ---: |
| XY uniform | 0.90909980 | 0.70481436 | 0.367195 | 0.986072 |
| XY 2:1 | 0.86806822 | 0.64945203 | 0.418585 | 0.981356 |
| RZ uniform | 0.74841828 | 0.58013398 | 0.367459 | 0.986196 |
| RZ 2:1 | 0.71942552 | 0.54038463 | 0.412859 | 0.985507 |

均匀周期横向分量的独立符号为 `lambda=-n*(1-exp(-i*2*pi/n))`。原 RK4 的精确离散振幅是 `P4(lambda/(32*n))^(32*n)`，能量损失为 `1-|振幅|²`。六个 uniform 组合逐 native Fourier 终值最大误差 **1.643e-14**，共 2,800 行；包括实际横向半控制体质量，未剔除边界。独立 Python 复算同一闭式表达式，与保存的实际损失符合原 1e−12 精度。

该模式连续时间累计损失 n16/n32 约为 0.912477/0.707632，阶 0.366789；半离散衰减率 `2*n*(1-cos(2*pi/n))` 的阶为 0.986072。累计量接近饱和会压低所测阶，不能因此把正确上风对照判作没有 O(h) 耗散。全向量结果还包含重复周期 normal native 与 RZ 非零均值，所以必须保留上表实际全场数据，不能用单 Fourier 模式替换实际 gate。

原 donor rhs 不读取中心 weights/couplings/apply。只修改获准的 CchKepConvection 文件，无法改变固定 donor 的这四项失败；因此当前“原 84 gate 全通过”的任务目标无法在既定范围内实现。本轮没有静默放宽 gate，也没有为不可改变的 donor 指标继续试配中心系数。H10n 将先形成显式指标设计变更，旧四项 FAIL 和其他空间/能量门槛必须保留，不能追溯改判。

## 当前算子的空间事实仍须处理

原 n8/16/32 affine 回归四组合共 **15,344 行**，其中 13,600 内部行最大误差 **9.335e-13 ≤ 1e−10**。同样 15,344 行的独立最终矩阵与实际 apply 最大误差 **1.236e-15**。内部、外部和固定轴贡献全量保存，未以内部通过代替整体通过。

| 组合 | n16 全域瞬时 L2 | n32 全域瞬时 L2 | 末档阶，要求≥1.7 | n32 最大无量纲行范数 |
| --- | ---: | ---: | ---: | ---: |
| XY uniform | 0.0548552 | 0.0186483 | 1.556583 | 2 |
| XY 2:1 | 0.062466 | 0.0324136 | 0.946471 | 7.3357 |
| RZ uniform | 0.0620346 | 0.0215769 | 1.523585 | 2.20708 |
| RZ 2:1 | 0.0702788 | 0.035795 | 0.973330 | 8.82787 |

周期半控制体与 AMR 系数增长问题仍在。三档范数增长是当前候选的负面证据，不是所有 one-ring 二阶方案不可能的证明。本轮未实现新最终算子、未新增周期 quotient API 或失败 fallback；生产仍 H10c。

## 同一 F 后的材料目标复核

原真实 moving FirstJet 的 XY/RZ、uniform/2:1、n4/8/16/32、rho1/1e3/1e6 共 **48 组 / 46,680 native 行**重新运行，原 dt=0.1/n 保持。同一 F 下，此族横向 affine 的全局材料功/PDE 相容缺口降至实际 gross 舍入量级：材料功最大 5.551e−16，PDE 两条件最大 6.475e−18；旧 quarter 反例不能继续用于否定新路线。逐 native 与逐 segment 表达式最大差 3.911e-16；CSV 独立重组矩与 RMS 最大尺度差 3.096e-15。

这仍不等于局部点 PDE 一致性。RZ 2:1、rho1000 的几何 centroid 目标对应横向导数 RMS 在 n16/n32 为 0.0382144/0.0275484，Linf 为 2.93704/3.30361；实际横向点 PDE 目标为零。后续应从最终点采样算子联合空间矩/skew，不能把几何一次矩直接强加为精确 PDE，也不能只消除全局缺口。本次没有构造 mixed-material 最终算子；原 selected 来源缺口、tiny/counterflow 与四硬约束的 H10l2 证据保持，不冒称这些场景在新算子下通过。

## 验收、冻结与限制

标准 Release 六 targets 构建通过；必要 CTest **8/8**，包括架构、自测、目录、phase、material、momentum、dual、KEP，总 181.51 秒。研究新增 Fourier 表达式首次构建有一个 complex/scalar 类型错误，已修正；没有数值断言放宽。最终文档/架构检查另 **3/3**。原内部 affine、全行矩阵和材料 observer exit 0；donor 原 gate exit 1 为实际 FAIL，不能混写成全部研究 PASS。

实际 scope 为 **0 production / 0 test / 1 existing research 源文件**；1,044 源成员中 1,043 保持，只有既有 driver 增加观察模式。原 Study/rhs/full/RK4/donor/gates 可逆核验保持。`compiled_runtime.json` 绑定 33 source / 8 binaries；最终研究 SHA256 为 `7B0612C368A4273705412F09973F66CEFE37FDD38D34740571D5847E4F4DE51C`。原观察运行绑定起始冻结 exe，donor 绑定本次 exe；生产库、短测试 binary 和原研究函数保持，完整来源保存在检查目录。

**未重跑原 28 窗口 / 84 gate 完整中心研究**：当前中心空间/有界前置未通过，且独立固定 donor gate 已确认失败。12 个 donor 窗口不是 28 个完整窗口，48 材料诊断也不是 PDE 轨迹。T0 工程与部分恒等式/原内部回归通过，最终算子目标和 T3 未通过，因此本项明确 FAIL；无外部阻塞，goal 保持 active/progress。

检查目录 `runs/checks/H10m/20261005_074800_106`；交付、独立核验脚本、原始 CSV、实际运行退出记录、source/runtime 和 flat manifest 见[冻结基线](../../runs/baselines/H10m_20261005/README.md)。历史报告和前置 archive 不改判。后续 H10n 只处理对照指标，最终空间离散和完整中心研究仍须单独完成；不进入 H11/H08a。
