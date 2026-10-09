# LAMP 数值硬约束

**版本：** 2026-09-30  
**状态：** 数值实现的硬门槛；优先于历史 TODO 和算例调参经验。

本文只规定**不能被实现细节破坏的数学关系**。具体公式来源见三份冻结方法文档。

**2026-10-09 BF-CSF方法变更注记：** [BF-01只读color](task/BF-01_Report.md)与[BF-02恒曲率载荷](task/BF-02_Report.md)规定范围已PASS，新方法规范为 [LAMP BF-CSF V1](roadmap/LAMP_BF_CSF_Method_V1.md)。`color_csf=H_eps(-phi/|grad phi|)`仅是同一FirstJet publication的只读力场派生，不能替代sharp相材料质量、成为输运状态或生成第二F_rho。整非pure panel须完成梯度认证；正权真实体积平均/统一epsilon/来源身份保持。BF-02直接返回实际RHS `+sigma*kappa0*C(color_csf)`，零sigma/曲率精确零；自由行/gauge/完整壁反力分别核验。原H12/H13 balanceLoads约定、C=-D_I^T、真实R/uR及mu=0继续保留。真实变量曲率及动力学gate尚未通过，旧H13判定保持，H14/H15不启动。

**2026-10-02 现行 Interface 变更：** 依用户授权，采用[纯 Hermite–GALS 非守恒路线](task/I12_Pure_Hermite_Gals_Design_Change.md)。I-R1/I-R4、I-R5 的 master/hanging 与唯一 geometry 要求保留；I-R2 改为 first-jet skeleton + 实际回溯中心预测，无独立 M，无体积校正；I-R3 及 I-R5 的体积求和/reflux 条款仅适用于历史守恒模式。新路线必须报告几何体积漂移并验证空间/时间/曲率收敛，不能以零库存残差声称守恒。下文保留旧数学关系用于追溯；Hydro H-R3/H07 的耦合前提须在接入前单独重设计。

## 1. Interface：纯 Hermite–GALS

### I-R1 唯一高阶几何

界面位置、法向、曲率、cut geometry 和 phase-face cut 必须来自同一份 accepted Hermite geometry。禁止为了质量、曲率或受力再建立第二份 production 几何。

### I-R2 状态分解固定

production state 为：

\[
(\text{transported first-jet Hermite skeleton})
\]

其中：

- `phi`、`grad(phi)` 直接随 characteristic 推进；
- `phi_xy` 只通过 compact compatibility closure 获得；
- 当前 FirstJet 模式不携带独立 `M_K`；
- Q4 内部值 `c_K` 直接来自同一旧几何的实际回溯中心预测，不做体积校正。

`V_K(phi_h)` 是当前唯一几何的测量值，其时间漂移与解析误差必须记录。历史 ConservativeVolume 模式仍保留独立 M 及[几何—库存精度预算](task/I12_Geometry_Inventory_Accuracy_Design_Change.md)，供显式 A/B；不得把旧库存约束隐式带回纯路线。

不得把 `phi_xy`、全局参数曲线节点或多自由度 KKT correction 重新定义为 production 主状态。

### I-R3 历史守恒模式：phase flux 唯一且共享

本条只约束显式 ConservativeVolume 比较；当前 FirstJet 路线不装配守恒相通量/库存。

每个内部 native face 的液相体积通量 `Q_l,f` 只计算一次，并以相反符号写入相邻 control volume：

\[
M_K^{n+1}=M_K^n-\Delta t\sum_{f\subset\partial K}s_{Kf}Q_{l,f}^{n+1/2}.
\]

2:1 coarse–fine 接口必须满足

\[
Q_{l,\text{coarse}}=\sum_j Q_{l,\text{fine},j}.
\]

**禁止**对更新后的体积分数直接 clipping；极小超界只能通过守恒的 flux-level redistribution 处理。

### I-R4 characteristic coherence

`B_h(x)` 与 `J_Bh(x)` 必须来自**同一个实际离散回溯公式**。`J_Bh` 不得用与 production tracer 不一致的事后 stencil 近似。

### I-R5 AMR skeleton

- refinement：优先 exact polynomial restriction；
- hanging data：由 master trace 代数生成，不独立推进；
- coarsening：纯路线验证 Hermite 几何投影精度及体积漂移；历史守恒模式的 phase volume 仍须严格求和；
- regrid 后只能发布一份新的 authoritative Hermite geometry。

### Interface 验收门槛

- `B_h/J_Bh` 与离散 map 数值微分一致并按设计阶收敛；
- uniform polynomial reproduction：目标多项式在舍入误差量级闭合；
- 纯路线不声称 internal phase-flux cancellation 或 coarse–fine reflux；这些仅是历史守恒模式的回归项；
- 静止界面不得产生系统 drift；
- 平移/旋转/vortex/AMR crossing 必须分别报告 interface、gradient、curvature、phase-volume error；
- 电喷雾 tip 验收必须单独报告 curvature noise 和 AMR sensitivity。

## 2. Hydrodynamics：CCH

**2026-10-03 现行耦合授权：** 用户选择[完全非守恒材料输运](task/H07_Pure_Gals_Hydro_Coupling_Design.md)。H-R3 中机器精度材料/dual-mass closure 改为显式非零缺陷核算与收敛验收；材料质量直接由唯一几何取得，无独立守恒库存或界面修复。唯一材料通量、常速度动量一致性仍必须成立，动能预算单列材料缺陷项；H-R1/H-R2/H-R5/H-R6/H-R7 和电荷/Gauss 关系保留。

### H-R1 原生未知量

pressure `p_c` 在 pressure cell；normal velocity `u_f` 在 native staggered face。AMR coarse face 由实际 fine subfaces 表示。

production Hydro 不得依赖持续存在的全局 C1/Q8 高阶速度未知量来获得局部 AMR 精度。

### H-R2 唯一散度与压力载荷

积分散度：

\[
D_I=B A_f,
\]

压力载荷必须由代码关系生成：

\[
C=-D_I^T.
\]

硬检查：

\[
\eta_{adj}=\frac{|u^TCp+(D_Iu)^Tp|}{\|u\|\|p\|}<10^{-12},
\]

verification case 目标接近 `1e-13`。

禁止独立装配第二套 pressure-gradient stencil；AMR 的高阶精度只能进入 compatible Hodge/reconstruction。

### H-R3 质量—动量唯一通量

\[
F_{\rho,f}=\rho_lQ_{l,f}+\rho_g(Q_f-Q_{l,f}).
\]

`F_rho` 是质量方程与动量对流的唯一质量通量。禁止由插值 `rho_f` 再构造另一套 momentum flux。

每个 momentum dual volume 从同一几何取得材料质量并显式核算自身非零连续性缺陷；AMR/regrid 前后检查来源与误差。动量必须处理该缺陷以保持无外力常速度，不允许用几何材料冒充守恒库存。具体离散一致性项由 H09/H10 验证。

**2026-10-04 用户授权的 H08/H09 重设计边界：** 材料耦合只保留四项硬约束：唯一共享 F_rho、内部面反号共享/AMR 实际子面求和、实际几何/advector/identity/clock 来源一致、KEP/skew 及 R/uR 结构。quarter 体积分配、`q_dual=T*q_primal` 和旧缺陷精确交换不再是新设计硬约束；固定 mesh 的空间交换差 `E_sp=q_dual-T*q_primal` 必须显式报告，真实 `R_dual=Delta_m_dual/dt+q_dual` 保持。精度和缺陷收敛仍是验收要求，无库存/bubble 或反向几何修复；见[耦合设计](task/H07_Pure_Gals_Hydro_Coupling_Design.md)。

### H-R4 低耗散对流

光滑、解析充分区域使用 symmetric/KEP 中心主通量。稳定化只能以显式、非负耗散的局部项加入，且满足“只耗散、不造能”。

`donor/upwind` 只能作为显式 benchmark/fallback-free comparison，不能成为 production 默认，也不能因为密度跳跃自动开启全局一阶耗散。

### H-R5 能量型黏性与无黏极限

\[
K_\mu=E_h^T H_\mu E_h,
\qquad u^TK_\mu u\ge 0.
\]

两相黏度按同一 Hermite cut geometry 分相积分。

当 `mu_l = mu_g = 0` 时，`K_mu` 必须在舍入误差量级为零；禁止设置人工最小黏度。

### H-R6 统一电—毛细法向应力

\[
J_\Gamma=\sigma\kappa+n\cdot[[T_E]]n.
\]

normal capillary 与 normal Maxwell stress 必须先组成 `J_Gamma`，再通过与压力相容的 interface measure 进入 momentum。

当 `J_Gamma=J0` 为常数时，必须存在

\[
p=p_0-J_0\alpha
\]

使 `Cp + F_Gamma,n` 达到线性求解/舍入误差量级。禁止毛细和 Maxwell normal load 各走一套 production deposition。

### H-R7 时间层

最终 production 采用单步 midpoint predictor–corrector。同一轮 corrector 中的 geometry、`Q_l/F_rho`、material Hodge、viscosity、curvature/Maxwell 与 momentum 必须属于同一中点候选。

不得用“上一步界面力 + 下一步速度”的一步滞后来掩盖强耦合问题。

### Hydro 验收门槛

- `D_I/C` adjoint identity：见 H-R2；
- `u^T K_mu u >= -O(eps)*scale`；
- `mu=0` sequence 连续，不因取消黏性突然切换算法；
- constant-velocity/high-density-ratio transport 不产生系统性界面速度误差；
- smooth inviscid benchmark 不出现 donor-cell 型 `O(h)` 系统耗散；
- static drop：先 analytic curvature，再 Hermite curvature；
- AMR static drop crossing：无系统 coarse–fine parasitic-current source；
- full second-order time accuracy 只有在 corrector 收敛后才能声明。

## 3. Electricity：SI-HHO-Q

### E-R1 HHO 能量结构

局部算子由 weighted potential reconstruction + symmetric nonnegative stabilization 构造。第一版固定 `k=2`：

- cell/face unknown：`P2`；
- local reconstructed potential：`P3`。

straight planar/RZ element 上的 `P3` reproduction 必须达到舍入误差量级。

### E-R2 只静态消元 cell DOF

局部 block：

\[
S_T=K_{ff}-K_{fc}K_{cc}^{-1}K_{cf}.
\]

禁止显式形成 `Kcc^{-1}`；使用稳定的局部 factorization/recovery。全局未知量仅保留 Cartesian phase-subface trace 与 Hermite interface trace。

### E-R3 全局 symmetry / SPD 目标

有 Dirichlet anchor 且 small-cut 已正确处理时，全局 skeleton matrix `S` 目标为 SPD。

硬诊断：

\[
\eta_{sym}=\frac{\|S-S^T\|_F}{\|S\|_F}<10^{-12}.
\]

### E-R4 equilibrated flux 唯一

potential、normal displacement、conduction current 和 Maxwell trace 必须从同一 HHO reconstructed field / equilibrated numerical flux 恢复。

禁止为 current 或 Maxwell force 再差分一套 `E`。

### E-R5 Gauss 与 interface jump

普通内部 face 的 flux moments 必须平衡；interface trace 必须满足

\[
D_n^+-D_n^-=q_h
\]

在 `P2(interface)` moment space 中成立。

运行时至少记录：local Gauss、internal face balance、interface charge-jump、global Gauss closure。production 目标相对残差 `<=1e-10`，verification 应更接近线性/舍入误差。

### E-R6 charge 时间离散

fixed geometry 下 `Q -> V -> I` 视为线性 affine map，surface charge 使用 variable-step BDF2–Schur；第一步使用 backward Euler。

不得继续使用欠松弛 Picard 作为 production charge relaxation；取消显式 `Delta t sigma/epsilon` 稳定性硬门槛不等于允许任意大步长保持精度。

### E-R7 small cut

第一版优先使用 cell agglomeration。small-cut 判据和 algebraic conditioning gate 必须显式记录；禁止通过无限增大 stabilization penalty 隐藏几何病态。

## 4. 跨模块硬约束

1. Interface 发布的 accepted geometry identity 在一个耦合候选内不可被其他模块偷偷重建。
2. `alpha/Q_l`、Hydro material state、Electric phase elements、Maxwell traction 必须记录所属 geometry/time layer。
3. rejected candidate 不得污染 accepted state、历史 charge、AMR topology 或输出 frame manifest。
4. 数值求解容限与物理守恒验收容限必须是不同的具名参数。
5. runtime residual 只能证明“离散方程被求解”，不能替代 MMS/解析收敛对“离散方程逼近连续 PDE”的验证。

## 5. 禁止项

出现以下任一情况，任务直接判定 `FAIL` 或 `BLOCKED`：

- 为让算例继续运行而静默切换旧算法；
- 通过提高黏度、加入全局 upwind、降低精度门槛来“修复”无黏不稳定；
- 同一物理量存在两套 production 通量/梯度/界面；
- clipping 守恒变量后只报告修正后的总量；
- 线性残差通过但 physical closure 未检查；
- AMR correction 破坏 `C=-D_I^T`；
- Electricity 重新引入 nonlocal P3 correction 修改全局 Poisson operator；
- 一个任务同时改变离散方程和验收标准。
