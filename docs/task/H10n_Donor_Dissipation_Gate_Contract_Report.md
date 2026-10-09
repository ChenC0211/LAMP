# H10n / donor 耗散率验收定义报告

2026-10-05，**H10n PASS / H10 未关闭**。按[事前设计变更](H10n_Donor_Dissipation_Gate_Design_Change.md)交付显式 `H10b-dissipation-v2`：四项donor对照使用初始离散耗散率的空间阶，旧累计损失阶保留原公式/FAIL，中心其余80项定义保持。生产算子、F/R/质量、原网格/时间窗/数值积分没有改变。本项修正独立对照的验收定义，不能当作中心空间精度通过。

## 实际数值与历史保持

四组合n8/16/32完成原 **12窗口 / 7,168步 / 396时序行**；最终执行正常exit0。新旧门槛均采用原事前阶区间 `[0.7,1.4]`，只改变被测物理量。

| 组合 | 初始D/E空间阶 | 新gate | 原累计损失阶 | legacy gate |
| --- | ---: | --- | ---: | --- |
| xy_uniform | 0.986072 | PASS | 0.367195 | FAIL |
| xy_amr | 0.981356 | PASS | 0.418585 | FAIL |
| rz_uniform | 0.986196 | PASS | 0.367459 | FAIL |
| rz_amr | 0.985507 | PASS | 0.412859 | FAIL |

`donor_audit.csv`、`trajectory.csv`、`donor_fourier.csv`与H10m冻结输出逐文件SHA256相同，旧`gates.csv`与本次`legacy_gates.csv`也逐字相同。保留全部half-dual/native和轴自由度，未缩短时窗或改参考。这里证明原保存的完整时序/误差及Fourier数据保持，并结合原数值路径可逆源码核验；没有声称H10m曾保存所有AMR逐步全场快照。

六组均匀横向模式继续用独立Fourier/RK4闭式参考逐native核验，**2,800行最大误差 1.643e-14**。新版本`gates.csv`中四项donor通过，`legacy_gates.csv`原四项失败；`gate_contract.csv`记录版本、四组旧/新名称、上下界及保留80项。未将历史84-gate结果追溯改判。

## 独立逐段能量身份

新观察器直接读取同一实际F、原周期配对和native质量。内部共享面计一次，周期两端按原几何配对计一次，未配对外段保留；从跳跃平方独立构造非负D，另记录原rhs功P、真实材料项C及开口边界功：

`P = -D - P_open + sum((R-q/2)*u²)`。

轴上真实u/rhs均为0，轴反力功为0；一般材料项没有被丢弃。本固定材料研究进一步核验 `-P=D`，不要求其他非守恒材料场景R为零。所有归一化使用实际通量/功项gross，零尺度要求精确零，没有经验floor。

原12窗口的初值/终值共 **24份快照、30,688 native行、92,288独立连接行**，覆盖原实际segment计数 93,632（周期两端在连接表合并一次）。独立Python从逐段F和逐点u重算q、上风载荷、R分解、D、能量和功，再核对summary；不读取production系数作为参考修正。

| 独立核验 | 最大尺度误差 |
| --- | ---: |
| shared segment耗散 | 1.630e-16 |
| native incidence | 8.699e-17 |
| 原donor native rhs | 1.467e-16 |
| 材料缺陷分解 | 8.078e-17 |
| 含材料/边界项的功身份 | 6.039e-17 |
| 本固定场景 -P=D | 6.062e-17 |
| summary独立重组 | 1.253e-14 |

全部在原1e−12内。上述最大值包括额外probe，probe另有496 native行、1,520连接行。材料功、开口边界功和周期反号分别验证；完整值在`actual_verification.json`，没有用总功抵消掩盖单点rhs差异。

## full路径接入与范围

full和donor-audit共用同一独立D诊断。原 n8 XY AMR `probe`完成中心+donor **2窗口/512步**，既有五项边界/材料/AMR gate通过；新诊断核验初值和donor终值。该中心粗网格终值L2仍 **0.6303866**，并不是精度通过。原其余80项的源码定义及网格/时间条件保持，full新版仍84项；本轮没有重跑已知空间前置未通过的28窗口完整中心研究。

scope为 **0 production / 0 test / 1 existing research源文件**。原Study构造、实际材料准备、rhs、exact/error/energy、RK4/integrate，以及除显式四项指标替换外的原driver，可逆还原到H10m源码；1,044源成员中1,043保持，无新target或永久工具。App/PDF与原算例输入保持。变更只在既有研究诊断/输出和gate，生产仍H10c。

## 工程验收与交付

标准Release六targets通过，**8/8必要CTest**通过（174.70秒），最终文档/架构另3项核验。初次新增观察器的一条混合auto声明编译失败，改为明确long double后通过，没有改数值门槛。所有数值运行均正常结束，没有超时重启。source/runtime绑定33份source与8份binary；仅研究exe改变，生产库和短测试binary与H10m相同。最终exe SHA256：`EBC000B8109536DEEF78769395621F4A2CD65384256C089DD49D918892D397BA`。

检查目录 `runs/checks/H10n/20261005_082057_279`；[冻结基线](../../runs/baselines/H10n_20261005/README.md)保存原始CSV、执行/退出记录、独立核验、设计前快照、canonical交付与flat manifest。H10m及此前失败报告/档案保持，不覆盖旧数据。

下一唯一[H10o](H10o_Bounded_Final_Native_Convection_TODO.md)实现最终native有界空间/KE联合算子：真正解决周期半控制体、AMR系数增长及点PDE一致性；先原内部affine，再明确v2的完整84项与单列legacy结果。本项未修复最终算子、未关闭H10，H11/H08a仍等待，goal active/progress。
