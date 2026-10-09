# H10o / 有界最终 native 空间对流算子

状态：**PASS，2026-10-05**；显式闭合算子、原内部反例、独立矩阵和v2完整84项通过，见[报告](H10o_Bounded_Final_Native_Convection_Report.md)。报告明确一般trace/切界面与效率范围，下一单任务H08a，H11未启动。

## 目标与范围

在 H10l2 同一 F、真实 R/三层质量与用户指定的四材料硬约束上，实现真正用于 production apply 的有界最终 native 对流算子，联合 constant/affine 空间一致性、二阶 PDE 余项和 material-KE/skew。不能以仅存在证书、单segment trace 或全局功平衡代替最终误差收敛。

允许现有 `CchKepConvection.h/.cpp` **至多2 production**、既有 `cch_kep_convection_contract.cpp` **1 test**，既有研究 driver/README、实施前设计/API/报告/当前入口。不得增加源码变体、CMake target、永久工具；临时定位归本任务 runs/checks，失败实现关闭时退役并归档。

禁止改已验收 Core/FirstJet、actual advector/boundary law、相积分、唯一 F、三层质量/R、D_I/C、黏性/电场、App/字典/PDF。Pure FirstJet 完全独立，无库存、bubble、体积拟合或质量修复。禁止加黏度、全局upwind、AMR位置/算例特例、失败fallback、隐藏自由度或任意扩大stencil。

## 数学与实现要求

1. 从最终实际矩阵含显式边界端口、固定轴和 uR 推导完整算子。硬常数和 material-KE/skew 联合；有可识别光滑意义的一次矩要求精确或事前设计的二阶余项。一般 sharp cut 不能强压任意精确 affine；任何余项须有真实质量/h误差范数依据，不按失败结果拟合预算。
2. 点 PDE 目标必须来自真实质量和同一 actual B。H10m已证明同一F下48组横向全局缺口接近舍入，但RZ的flux-centroid目标仍有局部误差；不能把几何矩直接当成点 PDE。重新核验必要兼容条件，不能用旧quarter反例判断新源不可行。
3. 在实际 geometry one-ring 上对最终共享 skew pair 做尺度化、秩揭示、约束最小改动。报告原始一次/二次余项、系数l1/算子范数、支撑和准备/更新成本。旧全局精确P2已知有过约束反例，不得原样重试；必须联合有界性与二阶余项，不能先任意拟合trace再负转置配对。
4. 周期端口必须显式声明几何关系；重复native采样及各自真实half-dual质量/误差保留。端口信息不足应显式拒绝，不猜测周期性、不删除boundary行。一般生产边界/接受事务仍归H14/H15；研究数据不得变成production问题特例。
5. 保留mesh/source/clock/instance/copy/lifetime/regrid/nonfinite拒绝，KEP只消费原共享材料F/R，不能另外修正F或偷偷令R=0。材料仅四硬约束；精度与稳定性验收仍必须通过。

## 验收 T0/T1/T2/T3

**T0：** 标准Release、架构/目录/受影响KEP与必要phase/material/momentum/dual回归；2/1预算、原成员/targets/inputs/App/PDF保持，实际source/runtime/退出码和manifest绑定。失败候选不得留在生产成为第二路线。

**T1：** 独立最终矩阵与实际apply核验；constant、material-KE/skew、四材料硬约束与原D_I/C共同保持。XY/RZ uniform/2:1、rho1/1e3/1e6、full/selected、fixed axis、显式boundary、tiny/counterflow与三层几何来源均须有明确实际覆盖。报告最终矩余项与系数界，非零R本身不作失败判据。

**T2：** 原内部affine counterexample n8/16/32首先通过，原1e−10回归及全部native误差来源保留；周期关系不能用于不满足其边界条件的affine样本。内部/coarse–fine/外边界/固定轴误差贡献加回整体。真实强制几何族与coupled制造族分别报告，不能拿目标诊断冒充新算子PDE通过。

**T3：** 前置通过后完整执行原28窗口/15,116步、n8/16/32、T=1和原时间研究。验收使用明确 `H10b-dissipation-v2` 的84 gate：原80项定义保持，新4项donor初始耗散率阶按H10n，旧4项累计损失阶继续另表输出。中心截断/终值空间阶均≥1.7、最细L2≤0.05、Linf下降、原KE/材料功/边界/时间门槛不变；不能称旧84gate已通过。实际进程运行到终态，保存全时序与失败前缀。

## 失败和交付

必要来源、结构、有界性、PDE精度或完整gate未通过即本项FAIL，保留反例和明确下一数学动作，不缩小为“诊断PASS”，不调阈值。若任务边界确需其他owner，应先按真实数学职责拆分，不扩大源文件堆积。仅实际外部依赖才可能BLOCKED。

交付最终production代码、独立短契约/完整研究结果、性能与限制、API/报告及冻结证据。H10整体完成前不进入H11；随后再按用户授权独立启动H08a同源时空材料通量升级。
