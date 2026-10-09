# H10m / 同一材料 F 下的最终 native 对流算子收口

状态：**FAIL，2026-10-05**，见[报告](H10m_Final_Native_Same_F_Convection_Report.md)。固定 donor gate 存在独立验收冲突，当前空间问题仍未修复；未交付新算子、未重跑原完整中心研究。下一唯一[H10n](H10n_Donor_Dissipation_Gate_Contract_TODO.md)明确对照耗散指标，原任务范围/门槛如下保留。前置 [H10l2](H10l2_Same_F_Dual_Material_Retry_Report.md) 完成同一F/真实R迁移；启动先复核冻结档案/current canonical。旧H10c–h的失败历史保留，不从旧quarter输入得出新路线不可行结论。本项直接处理最终native operator及原H10b验收，不再扩展单segment trace补丁。

## 目标与范围

在已验收 actual B、唯一F、真实三层材料质量/R 的固定来源上，联合 constant/affine spatial moment reproduction、二阶PDE一致性与现material-KE/skew结构，构造有界且可复用的native算子。保持Pure FirstJet独立，无库存/bubble；材料四硬约束固定，quarter/精确material commuting/R=0不作为约束。

允许现 `CchKepConvection.h/.cpp`（**至多2 production**）、现 `cch_kep_convection_contract.cpp`（**1 test**），既有 `cases/benchmarks/cch_convection/driver.cpp`及README、实施前设计/API/report/current状态。原CMake/targets，无新源码/tools/variant。临时诊断归 `runs/checks/H10m/<timestamp>`，关闭冻结baselines。若确需额外owner，先按实际数学边界拆任务，不扩散文件。

禁止修改已验收Core/FirstJet、advector/boundary、phase积分、唯一Flux、三层dual几何/质量/R、D_I/C、黏性/电场、App/字典/PDF及原H10b gates。禁止人工黏度、全局upwind、AMR位置/算例特例、trace后补丁、修正质量/通量或失败fallback。

## 实施前推导及实现

- 从实际 `apply` 的最终native矩阵推导含boundary port、固定轴、本体skew和uR的完整恒等式，明确真实mass/DOF采样与点PDE或体积矩目标的关系。旧trace affine exact不能代替这一推导；不能把几何矩目标未经核验当作点导数目标。
- 同一F固定后，constant/material-KE身份先作代数硬约束。可识别光滑affine内部行要求精确或明确二阶余项，与skew联合求解；不能对所有sharp cut行强压任意exact affine方程导致新的过约束。任何非零余项须来自一致目标/已证兼容性并满足预先固定的实际网格尺度预算，不能按失败数据拟合。
- 原segment支撑不足时使用真实one-ring patch作局部约束最小改动，在最终共享skew pair层施加矩条件；采用尺度化/rank-revealing处理，报告可行性、系数l1/算子norm、支持大小及准备/更新代价。无显式全局高阶速度未知量或任意扩散大stencil。不能用仅存在证书替代实际生产apply。
- 用当前真实F重新检查XY/RZ uniform/2:1、三密度及原selected缺口的兼容条件，独立参考不参与接受系数或状态。若固定目标仍不相容，先定位真实来源/边界/采样与不可满足方程，保留反例；不得为PASS放宽门槛或另造F。
- 保留完整来源/clock/instance/copy/lifetime/regrid/nonfinite拒绝。一般production边界/耦合接受事务仍归H14/H15；研究周期端口保持显式作用，不可偷偷删掉boundary native DOFs来改善整体误差。

## 验收 T0/T1/T2/T3

**T0：** 标准Release、architecture/self-test/layout、受影响KEP及phase/material/momentum/dual必要回归；2/1预算、无新source/targets/inputs，App/PDF保持。实际compiled source/exe/lib/DLL、退出记录和flat manifest绑定；失败候选退役，不残留第二路径。

**T1：** 实际最终matrix与apply独立核验，constant、material-KE/skew、原D_I/C及四材料硬约束共同通过。XY/RZ、uniform/2:1、rho1/1e3/1e6，fixed axis、full/selected/显式boundary、tiny/counterflow与三层几何。检查实际一次/二次矩余项及有界系数，不能只看segment trace或隐藏外边界误差。

**T2：** 首先重跑原内部affine counterexample n8/16并补n32，要求精确或经预先设计验收的二阶；保留原1e−10精确回归门槛与逐native真实误差来源，不能把改目标或排除行后的数值冒称原回归。独立分解内部/coarse–fine/外边界/固定轴，行平方贡献加回整体。真实材料强制运动族与coupled制造族分列，非零R不作为失败本身，PDE误差和收敛必须真实测量。

**T3：** T1/T2通过后完整重跑原H10b **28窗口/84 gates**：原空间阶≥1.7、最细L2≤0.05、Linf下降，以及全部时间、KE/材料功、donor对照和边界门槛；参数/步数/时间窗/DOF误差范围不改。跟踪实际运行句柄至终态，保存完整轨迹与失败前缀；观察器exit0、内部affine或局部身份均不能替代完整PDE通过。

失败规则：必要来源/结构/有界性/精度、原完整gate或预算不满足即FAIL并保存counterexample与明确下一数学动作，不能缩小目标为“仅诊断通过”。仅真实外部依赖反复阻断才BLOCKED。原H10b未全部通过前不进入H11。通过本项后再审核H10整体，并依用户授权单独启动H08a同源时空材料通量升级；不将时间积分升级混入本项。
