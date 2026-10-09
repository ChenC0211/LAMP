# H07 / 纯 Hermite–GALS 与 CCH 的耦合边界

现行设计，2026-10-03。用户本轮明确建议采用“完全接受非守恒材料输运”，本设计按该授权固定路线。依据原[纯路线变更](I12_Pure_Hermite_Gals_Design_Change.md)与[数值硬约束](../01_LAMP_Numerical_Invariants.md)。冻结 CCH PDF 保持不变；其 §4–5 的守恒材料库存前提由本设计取代。压力相容、动量输运一致性、黏性非负和电荷守恒不因此取消。

Interface 唯一 FirstJet publication 提供同一 Q4 的几何相体积、真实 native face trace、identity 和时间层，primary 不含 M。现有程序 `phi<0` 为液相；沿用此约定。几何不接受流体库存拟合或体积校正。

H07a 消费该 publication，测量 `V_l/V_c` 和各真实子面的液气面积，绑定 Hydro 的同一 native layout、D_I 和物理测度。当前相率由同一actual compatible B上的causal空间积分产生，H10l2已退役原面常值evaluate；唯一材料F由H08逐atom密度加权，H09直接共享原dual视图。共享内部 face 的外流代数抵消，但这不证明下一时刻 GALS 几何体积满足 FV 质量更新，也不能将 `dt*q_l` 称为已验证的时空相通量。

[H07b](../api/H07b_Cch_Phase_Quadrature_Viscosity.md)已提供位于实际液气区域内的正权体积求积，绑定相同 geometry/time，并接入 H06 的 `E^T H_mu E`。按总体积分数在同一组全单元点混合黏度，或中心/顶点判相，都不能代替分相积分。H07a/H07b 均已验收，H07 关闭；H08已绑定材料时间积分/非零连续性缺陷，[H09a](../api/H09a_Cch_Momentum_Geometry.md)已实现真实dual几何/相积分质量，[H09b](../api/H09b_Cch_Dual_Material.md)已通过共享流量、真实dual缺陷和uR一致性验收，H09关闭。

定义只读几何质量

\[
m_{geo,c}(t)=\rho_l V_{l,c}(\phi(t))+\rho_g[V_c-V_{l,c}(\phi(t))].
\]

纯 GALS 下一层几何质量一般不满足

\[
m_{geo,c}^{n+1}-m_{geo,c}^{n}+\sum_f s_{cf}\int F_{\rho,f}\,dt=0.
\]

生产材料质量直接取 `m_geo`；密度与黏度/Hodge 来自当前唯一几何的分相积分。不存在另一份守恒材料质量状态，不做质量修复、界面体积拟合或通量 redistribution。几何质量的变化就是实际材料变化，不以它覆写一个名义守恒库存。

[H08](../api/H08_Cch_Material_Transport.md)已绑定开始/中点/结束几何及唯一中点时间积分材料流量，并显式计算逐单元材料连续性缺陷

\[
R_c=(m_{geo,c}^{n+1}-m_{geo,c}^{n})/\Delta t+\sum_f s_{cf}\bar F_{\rho,f}.
\]

`R_c` 可以非零，不能填零冒称守恒。H09a 在 momentum dual volume 与实际 cell 的交集中测量同源相容积，直接形成几何质量；不把 cell 平均密度乘交集体积当作真实相积分。H09b 给出相应 dual 缺陷。生产不再要求这些缺陷机器精度为零；仍要求其定义、空间/时间收敛、AMR/regrid 误差和真实运行漂移可观测。具体归一化、界限及接受规则须在对应 TODO 中先定后测，不由本文件发明宽松数值容差。

**历史 H09b 实现（非新设计硬要求）：** 固定 mesh 的流量转移采用唯一 primal 材料流量，满足 `delta_u*R_F = R_m*delta_c`；几何权重 `R_m(f,c)=V_fc/V_c`。真实 cut dual 质量 `m_geo,f` 一般不等于 `R_m*m_geo,c`，故历史 H09b 独立计算

\[
R_{dual}=\Delta m_{geo,dual}/\Delta t+\delta_u\bar F_{rho,dual}
=R_m R_c+(\Delta m_{geo,dual}-R_m\Delta m_{geo,c})/\Delta t.
\]

第二项是明确的几何质量转移差异，不能填零、用平均密度覆盖真实dual质量，或只把 cell 缺陷插值成 dual 缺陷。H09a 已提供真实交集和几何权重；[H09b验收](H09b_Dual_Material_Consistency_Report.md)核验的是上述旧转移关系的代数自洽，不证明任意 cut 材料空间精度。重设计按下述授权允许非零空间交换差，旧报告保持历史原判。regrid 后需要重新绑定几何/布局，不能继续使用固定mesh公式的旧系数；H09b测试了实际Interface transfer后重建与常速度续算，一般速度transfer仍归H15。

动量离散必须处理这一材料变化。连续层示意为

\[
\partial_t(m u)+\operatorname{div}(F_\rho u)=F+uR,
\]

等价于相同材料来源的速度形式，保证无外力常速度不因 `R` 产生虚假加速。`uR` 是非守恒材料输运的显式一致性项，不能暗中修复材料、改变几何或重新插值一套 `rho_f*u_f`。H09b已验证同源积分项 `u_mid*r_dual` 及高密度比常速度材料更新；一般对流/边界/完整midpoint identity仍待H10/H14/H15验证，本文件不是已实现的 NS 方程。

相应动能收支含材料项 `|u|^2 R/2`。必须把它与压力功、黏性/受控耗散和外力功分开核算；不能声称非守恒材料下无黏动能无条件严格守恒，也不能用额外黏度掩盖该项。[H10c最终native联合矩](../api/H10c_Cch_Native_Convection_Moments.md)消费H09同一F/R，在[原配对装配](H10_Compatible_Center_Convection_Design.md)上增加共享skew系数，修复内部affine反例并保持半离散材料动能/显式边界功恒等式；无材料修复或隐藏黏度。[完整复验](H10c_Native_Convection_Moments_Report.md)仍 **FAIL**，[H10d](H10d_Final_Operator_Consistency_Report.md)候选FAIL且已退役，[H10e](H10e_Bounded_Operator_Consistency_Report.md)联合有界候选FAIL并退役，48组真实材料空间矩相容性诊断完成；[H10f](H10f_Compatible_Bounded_Moments_Report.md)定位真实材料几何矩/PDE目标差异；[H10g](H10g_Point_Pde_Material_Moments_Report.md)完成真实质量误差范数及原材料矩来源审计，目标修正FAIL并退役；[H10h](H10h_Final_Operator_Joint_Feasibility_Report.md)固定最终矩区间四组FAIL并退役，定位旧空间转移的纯气体大q/R。用户本轮解除附加材料过约束，详见下述四硬约束；H10i的source缺口经H10i1a–f补齐，H10j0/j/k/l1及H10l2完成actual B/phase/唯一F与真实dual缺陷；H10m历史空间失败保持；H10n明确donor指标后，[H10o](H10o_Bounded_Final_Native_Convection_Report.md)的显式闭合算子完成v2原28窗口/84项验收。下一单任务H08a，H11未启动；一般切界面与外部trace能力边界见报告。变步长corrector和流体联合验收仍待对应任务；历史H10b失败保持。

本授权固定材料路线。H08材料流量/缺陷核算、H09a真实dual几何质量与H09b共享dual流量/真实非零缺陷/uR积分项均已实现并验收，H09 PASS，下一项H10。中点时间积分已验证光滑面率二阶，整面phase rate时间跳跃的实际非零缺陷仍保留，不提前宣称完整sharp-interface时间二阶。密度比1e6的常速度材料一致性已通过；一般高密度比流动、无黏传输及真正的两相时步，必须在对应后续gate通过后才可宣称可用。旧守恒baseline可比较，但其零质量残差不再是新路线的验收目标。

正式只读跨模块依赖 `lamp_hydrodynamics_cch -> quantum_interface` 用于消费已有 publication。Interface 不依赖 Hydro；旧 Hydro、Compatibility 与 App 不进入 CCH 依赖闭包。不复制 snapshot，不穿透 Interface detail 实现。

## 用户指定的H10c后顺序（2026-10-04）

H10c最终native空间离散通过之后，启动[H08a同源时空材料通量升级](H08a_SpaceTime_Material_Flux_Upgrade_TODO.md)。Pure FirstJet/Hermite–GALS几何保持完全独立，Hydro在同一实际步内几何/characteristic上积分共享相通量并构造唯一F_rho，目标是使R成为小的可收敛离散余项。仍不引入库存、bubble、逐cell体积拟合或质量修复；继续既有R/uR、材料动能项和KEP相容对流。当前H08是已通过的中点近似，不能将其历史证据重新称为时空积分；本段只登记后续工作。

## 用户授权解除 H08/H09 过约束（2026-10-04）

H08/H09 重设计的材料耦合**只有四项硬约束**：唯一共享材料通量；内部面反号共享/AMR 子面守恒；来源一致；KEP/skew 结构。不同真实面的 primal/dual 求值必须属于同一 F_rho owner，同一共享面只有一份值。这里的守恒是内部面代数抵消和子面求和，不是非守恒 Interface 几何的逐 cell 质量库存闭合。

旧 quarter 体积/tensor 分配、`q_dual=R_m*q_primal`、真实 dual 质量等于映射质量、旧缺陷无空间交换差和严格 R=0 不再作为新设计的等式约束；可用于比较和诊断。纯相/混相材料精度、实际导数行尺度、R 的空间/时间收敛及原动量研究仍须验收，不能用解除等式豁免 PDE 精度。Pure FirstJet 几何完全独立，无库存/bubble、体积拟合或质量修复；其他模块原硬约束保留。

固定 mesh 上，定义 `E_sp=q_dual-R_m*q_primal`，新缺陷分解为

\[
R_{dual}=R_mR_c+(\Delta m_{dual}-R_m\Delta m_c)/\Delta t+E_{sp}
=\Delta m_{dual}/\Delta t+q_{dual}.
\]

E_sp 可以非零，须和几何转移差、实际 R 分别报告；同层 uR/材料动能项和 KEP/skew 消费实际 dual F/R。新结构并不依赖 E_sp=0。regrid 后按实际源重建，不外推旧映射。

[H10l2](H10l2_Same_F_Dual_Material_Retry_Report.md)完成H08/H09唯一共享空间F迁移并删除旧quarter/legacy。真实R/E_sp/G和材料动能身份经独立三层质量/物理积分验收；Pure FirstJet独立，无库存/bubble。H10o已通过显式闭合最终native空间/KE算子及v2完整研究，旧四donor gate失败保持；下一唯一H08a，按真实材料兼容性设计时空积分余项，不能强迫任意sharp-cut affine或反修几何。H11未启动。clip账本中的真实重合界面矩与材料面反号共享分别核验，不施加旧交换关系。

**基础设施收口（2026-10-05 用户授权）：** r11 是 common regular source/characteristic 的最后一道关口。除非出现可复现且直接阻断 H08a 的数学反例，不再拆 source/Jacobian/低秩精修任务，不新增硬约束、特殊 case，也不放宽原门槛。r11 通过后立即回到 H08a 的真正同源时空相通量积分；必要的 FirstJet 运动消费接入属于该积分任务，不再单列桥接先决任务。目标仅为 R 成为可收敛的小离散余项，保持 Pure FirstJet 独立、无 M_K/bubble、无 R=0 或 quarter commuting 强制。材料硬约束仍仅为唯一共享材料通量、内部面反号共享/AMR 子面守恒、来源一致和 KEP/skew。
