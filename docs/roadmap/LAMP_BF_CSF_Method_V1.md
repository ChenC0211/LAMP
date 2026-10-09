---
title: "LAMP-BF-CSF V1：基于 Hermite–GALS 的单流体连续表面力与相容离散方法"
subtitle: "数值方法规范、离散方程、源码映射与 Codex 分阶段实施任务"
date: "2026-10-09"
lang: zh-CN
---

# LAMP-BF-CSF V1：基于 Hermite–GALS 的单流体连续表面力与相容离散方法

**文档定位：数值方法设计提案 / Codex 实施基准（尚未实现、尚未验收）**  
**适用项目：LAMP，Hydrodynamics CCH / Pure FirstJet Interface**  
**目的：以规范的 one-fluid CSF / balanced-force 方法替代 H13 失败的变量曲率候选；保留已经验收的非毛细空间骨架。**

**2026-10-09 实施进度：** [BF-01只读color](../task/BF-01_Report.md)及[BF-02恒曲率载荷](../task/BF-02_Report.md)规定范围已PASS；本文件完整方法及后续动力学仍未验收，BF-03/H14/H15本轮未启动。下文公式与原物理gate保持。

> **给 Codex 的首要说明。** 该文档是新的方法变更提案，不是当前代码已实现功能。请先审查 `docs/roadmap/05_LAMP_Module_Task_Roadmap.md`、`docs/task/H13_Hermite_Capillary_Traction_Report.md` 与现有 CCH API；按第 15 节逐任务实施。不要把原 H13 的 `capillary()` 改名后宣布完成，不要同时启动 H14/H15，不要修改 Pure FirstJet 为体积守恒方法，也不要抹去历史失败记录。各阶段先写可执行断言，再改生产代码；未满足 gate 则报告 `FAIL / BLOCKED`。

## 0. 方法决策与适用边界

### 0.1 需要纠正的数值方法定位

现有 LAMP 使用 **一套原生 MAC 压力–速度自由度**、高阶 Sharp Hermite 几何和分相材料积分；当前 H13 的变量毛细载荷基于 `C*alpha` 的局部应力映射。它有严格的常应力平衡，但并不等同于已经验证的 CSF/balanced-force 两相动力学。项目已记录完整无黏毛细波的表面能/动能失配，且 H14 完整时间推进尚未完成。

**V1 采用单流体两相模型，而不是单相模型**：液体与气体具有不同密度、黏度；用单一速度和压力空间求解一流体不可压缩 Navier–Stokes；界面力作为连续分布的体积力，经与压力共用的离散梯度施加。**Sharp Hermite 负责几何；CSF color 负责连续力；几何材料积分仍可保持相分辨。**

V1 是一个经过文献验证的**方法家族**在 LAMP 网格与接口上的**具体设计提案**，并不保证当前低黏/无黏算例直接通过。尤其，原 H13 的变量曲率离散失败，不应靠仅增大平滑宽度或人工黏度掩盖。

### 0.2 V1 必须保留的代码约束

1. **唯一界面**：`Pure FirstJet Hermite–GALS publication`，`phi < 0` 为液相。不得新增独立相函数输运、PLIC、保存型 `M_K`、bubble 体积修正或 redistance。
2. **唯一流体拓扑**：原生 MAC face normal velocity `u_f` 与 cell-centered pressure `p_c`；2:1 AMR 使用真实 fine subfaces。
3. **积分连续性**：`D_I=B A_f`，`D_I u=0`；**压力载荷只能从** `C=-D_I^T` 构造，不额外建立另一套压力梯度 stencil。
4. **材料与惯性**：保留 H07–H10 的同源 sharp 相几何质量、唯一共享 `F_rho`、实际 dual 质量、真实非零 `R` 及对应 `uR`；不暗改为 `R=0`。
5. **对流、黏性**：保留 KEP/skew + 局部非负耗散；`K_mu=E_h^T H_mu E_h \succeq 0`，`mu=0` 严格消失，禁止 viscosity floor。
6. **几何身份**：新的 color、曲率映射和力必须消费与本次 `CchPhaseGeometry` 相同的 `mesh / layout / geometry identity / full time layer`；regrid 后全部重建。
7. **开发状态**：历史 H12 常应力工厂及 H13 失败候选继续作为 A/B baseline；不能将当前 H13 单步研究驱动当成 H14 production time integrator。

### 0.3 两个不同的“相指示量”必须明确区分

为避免在代码中用同一个变量名隐藏两种物理职责，定义：

- `alpha_sharp[K]`：现有几何液相体积分数，`int_K H(-phi_h) dV / V_K`。**用于真实两相体积、材料惯性与已有质量通量，非独立时间推进。**
- `color_csf[K] = c_K^eps`：由**同一 Hermite 几何**派生的正则化 color（见 §3），**只负责构造 CSF/balanced-force 的有限厚度力层和可视化诊断**。它不是一个新增演化状态，不得用于直接覆写 sharp material mass。

这一选择属于 *sharp material / regularized force one-fluid formulation*。它在连续极限上逼近相同两相界面，但离散有限网格上两种测度不完全一致，所以必须额外审计体积、能量和低密度相的响应。

> **非本版范围**：若未来选择密度、黏度也由 `c_eps` 正则化，必须**同时**重新定义动量质量、对流质量通量、双时间层材料变化与 `R`；不允许仅把 `rho` 由 sharp 改成 `rho(c_eps)` 而继续消费 sharp `F_rho`。本版禁止这种混合实现。

## 1. 连续方程与符号

### 1.1 两相一流体 Navier–Stokes

令 `n=grad(phi_h)/|grad(phi_h)|` 为**液体指向气体**的法向，`kappa=div(n)`；对凸液滴 `kappa>0`。单一速度与压力满足

\[
\nabla\cdot\mathbf u=0,
\]

\[
\partial_t(\rho\mathbf u)+\nabla\cdot(\rho\mathbf u\otimes\mathbf u)
=-\nabla p+\nabla\cdot(2\mu D(\mathbf u))
+\mathbf f_\sigma+\mathbf f_E+\mathbf f_b,
\qquad D(\mathbf u)=\frac{\nabla\mathbf u+\nabla\mathbf u^T}{2}.
\tag{1.1}
\]

在分布意义下，`chi_l=H(-phi_h)`，

\[
\rho=\rho_g+(\rho_l-\rho_g)\chi_l,
\qquad
\mu=\mu_g+(\mu_l-\mu_g)\chi_l.
\tag{1.2}
\]

这里的系数允许强跳跃；单流体方法并不等于把 `rho_l=rho_g` 或 `mu_l=mu_g`。

### 1.2 CSF 的连续力

采用宽度为 `eps` 的液相 color `c_eps \simeq H(-phi)`。对常数 `sigma>=0`，定义

\[
\boxed{\mathbf f_{\sigma,\epsilon}=\sigma\,\kappa_h\,\nabla c_\epsilon.}
\tag{1.3}
\]

由于 `grad(c_eps)\approx -n delta_eps`，此式对凸液滴产生指向液体内部的力，且静态平衡为

\[
\nabla p=\sigma\kappa\nabla c_\epsilon,
\qquad
p\approx p_g+\sigma\kappa c_\epsilon
\tag{1.4}
\]

（常数曲率时严格成立）。**这个符号选择与源码中 `CchInterfaceNormalLoad::balanceLoads()` 的负号约定不同，进入鞍点右端时必须显式转换**（见 §5.3）。

连续界面能量 `E_sigma= sigma |Gamma|`，完整能量平衡涉及 `f_sigma · u` 与面积变化率；不能只依据式 (1.4) 宣称变量曲率下动能–表面能离散相容。

### 1.3 RZ 轴对称

二维子午面坐标 `(r,z)`，物理体积权重

\[
dV=2\pi r\,dr\,dz,
\qquad dA_\Gamma=2\pi r\,ds,
\quad
\kappa_{RZ}=\kappa_m+\frac{n_r}{r}.
\tag{1.5}
\]

轴上采用 Hermite publication 的正则极限，不额外计算 `1/r` 的奇异差分；所有颜色积分和压力梯度都消费 Core 原生真实 RZ 几何测度。

## 2. 离散空间与共同算子

网格为 native quadtree 2:1 AMR。压力单元集合 `C`，真实原生面集合 `F`。

\[
(D_Iu)_K=\sum_{f\subset\partial K}s_{Kf}A_fu_f,
\qquad D_I=BA_f,
\tag{2.1}
\]

\[
\boxed{C=-D_I^T,\qquad (Cp)_f=A_f(p_R-p_L)}
\tag{2.2}
\]

（式中内部面从 `L` 指向 `R`，自由边界需要沿原 `C` 的反力/牵引约定处理）。

因此

\[
u^TCp=-(D_Iu)^Tp,
\qquad D_Iu=0\Rightarrow u^TCp=0.
\tag{2.3}
\]

动量加速度是 `H_{rho,u}^{-1} C p`，不是不带质量权重的 `C p`。本版不修改 AMR Hodge / pressure acceleration 的后续 H15 工作；不能将 `C=-D_I^T` 误说成 H15 高阶 AMR 精度已经通过。

## 3. CSF 正则化 color 的唯一构造

### 3.1 为什么不直接复用旧 `alpha_sharp`

若取 `c=alpha_sharp` 且采用相同的 `kappa_f`，

\[
S_{\sigma,f}=\sigma\kappa_f(C\alpha)_f
\tag{3.1}
\]

在数值上与 H13 的某类对角 `J_f C alpha`（转成物理 RHS 后）属于同一结构，**不会因为重命名为 CSF 自动消除原失败**。V1 的有效方法变更包括：明确有限厚度 color、同源 CSF 梯度、曲率在相同力带内的投影、独立能量/谱诊断；原 H13 保留为单独 A/B baseline。

### 3.2 确定正则化函数

使用有界、单调、连续可微的 regularized Heaviside：

\[
H_\epsilon(s)=
\begin{cases}
0,&s\le-\epsilon,\\
\displaystyle\frac12+\frac{s}{2\epsilon}
+\frac{1}{2\pi}\sin\!\left(\frac{\pi s}{\epsilon}\right),&|s|<\epsilon,\\
1,&s\ge\epsilon.
\end{cases}
\tag{3.2}
\]

在过渡带内相应的 delta kernel 为

\[
\delta_\epsilon(s)=H_\epsilon'(s)=\begin{cases}
\dfrac{1+\cos(\pi s/\epsilon)}{2\epsilon}, & |s|<\epsilon,\\
0, & |s|\ge\epsilon.
\end{cases}
\tag{3.3}
\]

**数值输入应是具有长度量纲的界面符号距离 `d_h`，而不是未经缩放的任意 level-set 标量。** 由于 Pure FirstJet 不保证 `phi` 始终为 signed distance，V1 建议先采用

\[
d_h(\mathbf x)=\frac{\phi_h(\mathbf x)}{|\nabla\phi_h(\mathbf x)|}
\tag{3.4}
\]

作为局部一阶符号距离近似，仅在正则化力带内使用。要求 `|grad phi|` 大于预先给定的**相对于输入局部标度**的可辨识阈值；不足时报告 unsupported/degenerate geometry，而不是夹一个无量纲 `1e-12` 或让 `d=0`。若验证发现切割尖端或 AMR 上该近似不够精确，再从同一 Hermite 曲线局部投影构造真正的 signed distance；此升级不能悄悄改变 `phi` 或独立创建另一条界面。

对于位置 `x`，

\[
c_\epsilon(x)=H_\epsilon[-d_h(x)].
\tag{3.5}
\]

`d_h` 与曲率均由同一接受时间层的 Hermite publication 求值。

### 3.3 单元 color 是物理体积平均，不是中心点函数值

\[
\boxed{
 c_K^\epsilon=
\frac{1}{V_K}\int_KH_\epsilon[-d_h(x)]\,dV.}
\tag{3.6}
\]

XY 使用 `dx dy`；RZ 使用 `2*pi*r dr dz`。严格使用正权、有界的单元积分，允许沿 `d_h=±epsilon` 划分子区间并作自适应 Gauss 求积；禁止把 `H_eps(-phi(center))` 直接当成全单元 color。

代数与几何要求：

- `0 <= c_K^eps <= 1`；纯气/纯液且过渡带不相交时**精确返回 0/1**；
- 固定同一网格与 publication，多次构造逐位稳定；
- uniform/2:1 AMR 的各 cell 独立以**真实体积积分**定义 color；regrid 后按同一几何重积分；
- 不独立输运、clipping、体积修补或改变 `alpha_sharp`；
- 存储 `||c_eps-alpha_sharp||`、力带面积/体积、色函数对 `epsilon` 的敏感性。

### 3.4 力带宽度（显式数值参数）

定义正常向分辨长度 `h_perp` 为当前界面邻域的网格长度尺度。首版**先在 uniform XY/RZ 验证**，建议实验初值 `epsilon=1.5 h`，并同时测试 `epsilon/h=1.0, 1.5, 2.0`。这些是**扫描参数，不是已验证的最优值**。

对于 2:1 AMR，不能在两个相邻 cell 内独立任意更换 `epsilon_K` 而不报告 color/force 的变化。首版建议对每个连续界面片段及其 force band 使用一致的物理长度 `epsilon`（参考邻域最细空间尺度），在该值下构造全局同一 `c_eps`；随后才研究分辨率适配的连续 `epsilon(x)` 及其影响。若力带未被足够的真实 cell / 面分辨，则**拒绝该配置或触发 AMR 加密请求**，而不是自动引入任意扩大宽度。

注意 `epsilon` 的收敛路径：固定物理 `epsilon` 加密网格得到的是**正则化方程**的数值收敛；若要验证对原 Sharp Interface 问题的收敛，需令 `epsilon = O(h)` 并检查毛细动力学。

## 4. 曲率的同源投影：从 Hermite 界面到 native force face

### 4.1 不重建另一套曲率

使用接受的 FirstJet publication 的几何积分点 `q`、物理面积权重 `w_q`、法向 `n_q` 和曲率 `kappa_q`；不从 `c_eps` 反算低阶曲率，不做 level-set redistance，不采用单元中心二阶差分曲率覆盖 Hermite 结果。

为每个 native 速度面 `f`，构造一个局部、非负、紧支撑权函数 `W_f(x) >= 0`，其支持覆盖 `f` 相邻 pressure leaves 加上整个 `epsilon` 力带所需的邻域。`W_f` 应使用物理空间尺度定义且在 XY/RZ/2:1 AMR 上相容（第一版可采用面邻域的一/二环正权帽函数）。

\[
\boxed{
\kappa_f=
\frac{\sum_{q\in\Gamma_h}w_qW_f(x_q)\kappa_q}
{\sum_{q\in\Gamma_h}w_qW_f(x_q)}
}
\tag{4.1}
\]

其中 RZ 的 `w_q` 已包含 `2*pi*r_q ds`，不能再次乘 `r`。只有在 `|(C c_eps)_f|` 达到由真实面通量/力尺度设定的可辨识值时，才要求有非零的曲率支撑分母；若 `C c_eps=0`，该面可直接令毛细贡献零。**若力带有效而缺少曲率支撑，必须报告失败，不能从远方界面猜一个曲率。**

定义并监视：`supportCount`、`totalWeight`、`kappaMin/Max`、近壁接触线和 AMR 处的降阶标志。只允许固定的、可测的局部正权投影；不要让 Codex 调一个新的全局滤波参数以压低伪流。

### 4.2 曲率常量再现是代数硬条件

若 `kappa_q=kappa0` 对所有可用支撑点成立，则必须有

\[
\boxed{\kappa_f=\kappa_0}
\tag{4.2}
\]

到舍入误差，与位置、两相密度、黏度、RZ 测度、AMR 粗细接口和 `epsilon` 无关。

首版只要求局部正权平均的常量再现及非退化性；二维曲界面上的至少一阶/二阶曲率恢复可作为后续改进，但必须先完成完整真实 Hermite 曲率的网格收敛检查。

## 5. Balanced-force 离散：最关键的动量右端

### 5.1 不使用第二套界面力梯度

定义 `g_csf = C c_eps`（**积分梯度载荷**，量纲为物理面面积乘无量纲 color 跳跃）。

在内部真实面 `f` 上：

\[
(g_{\rm csf})_f=A_f(c_R^\epsilon-c_L^\epsilon).
\tag{5.1}
\]

\[
\boxed{(S_{\sigma})_f=\sigma\kappa_f (g_{\rm csf})_f.}
\tag{5.2}
\]

`S_sigma` 是**作用于动量方程 RHS 的实际毛细驱动载荷**，不是旧 API 所称的 `balanceLoads()`。它必须使用原 `CchPressureLoad` 的 `apply`；不得在 AMR、axis、ghost cell 中再建压力梯度/CSF gradient stencil。不要额外乘 `A_f/d_f` 或 `1/rho`：后者由动量质量/Hodge 中统一处理，否则会产生双重度量除法。

该格式离散了 `sigma kappa grad c_eps`，具有真正的**连续力层**：`c_eps` 在原 Hermite 界面周围以有限厚度平滑过渡；梯度对所有 face 来自一张 cell-centered color。它的连续力平衡不是在两相上强制显式压力跳跃，而是在细化极限趋于正确分布应力。

### 5.2 常曲率静平衡必须代数精确

若 `kappa_f=kappa0`，

\[
S_{\sigma}=\sigma\kappa_0 Cc^\epsilon.
\tag{5.3}
\]

因此对给定边界/压力 gauge 允许的解，取

\[
\boxed{p^*=\sigma\kappa_0 c^\epsilon+p_0,}
\qquad C p^*=S_\sigma
\tag{5.4}
\]

在**自由动量行**应严格成立。闭合域中采用真体积加权均值 gauge；由给定外边界的常数压力/反力项由既有边界约定负责，不能强迫 `C 1=0` 在固定边界全部行上成立。

不要把已知解析 `p^*` 直接写入 solver 作为静滴特殊分支。必须实际求解既有鞍点系统，再检测速度、压力和自由行平衡残差。

### 5.3 与现有 H12/H13 API 的负号映射

现有 H12 / H13 将负向法向应力的平衡项记为 `F_balance`，并在构造 `A u + C p = rhs` 时使用 `rhs -= F_balance`。

本设计选择直接返回物理动量右端

\[
\boxed{\texttt{momentumRhsLoad} = S_\sigma=+\sigma\kappa_f(Cc^\epsilon)_f.}
\tag{5.5}
\]

对应源码中 FirstJet `phi<0`、法向为液→气、凸滴 `kappa>0`。在提供 legacy 对照接口时，可定义 `legacyBalance = -S_sigma`，但**生产只允许一个符号入口**。写一个凸圆滴解析 sign contract：`p_l>p_g`、无外力 `u=0`、`A u+C p-S_sigma=0`；错误符号必须在这一测试中被拒绝。

### 5.4 为什么这依然不足以保证变量曲率稳定

当 `kappa_f` 非常数时，一般不存在一个标量 `p` 使 `C p = diag(kappa_f) C c_eps` 在所有行同时成立。该非梯度部分是驱动毛细波的物理作用；但它也可能因为曲率转移、界面运动速度以及离散力不相容而错误注入能量。

所以 **V1 的代数硬保证是常曲率平衡，不是全变量曲率的无条件离散能量守恒**。新格式通过 `epsilon` 控制的连续力带和正权曲率同源投影建立一个规范 A/B 基线；能否解决 H13，必须由完整毛细波、不同 `epsilon` 和完整谱实验决定。若仅改变 `epsilon` 而性能未改善，必须返回 FAIL，不得调 `sigma` 或最低黏度以满足指标。

### 5.5 CSS / 变分增强仅作为后续研究

更深层的毛细能量形式为

\[
E_\sigma(q)=\sigma A_\Gamma(q),
\qquad F_\sigma^{\rm var}=-\mathcal A(q)^T\nabla_qE_\sigma(q),
\tag{5.6}
\]

其中 `q` 是真实 FirstJet 几何状态，`dot q = A(q) u` **必须包括实际特征回溯、Jacobian、mixed closure 和已接受几何发布的完整方向导数**，不能把它简化为某个点上的 `u·n`。如果这一力与离散面积变化之间形成准确的功配对，才可推导相应的能量关系。

**本版不把式 (5.6) 作为实现任务**。H13 历史实验证明，局部功拟合、压力 lifting、Gram 反馈等候选可能改善某个频率，却导致负刚度模态或远壁面力。将来研究能量相容变分法时，应独立开设计变更并保留基线，不能偷偷混入当前 BF-CSF 工厂。

## 6. 两相材料、惯性、对流与非守恒几何余量

### 6.1 保留 sharp 物性：不把 color 当新材料库存

当前 `alpha_sharp` 继续决定

\[
m_K=\rho_l V_{l,K}^{\rm geom}+\rho_g(V_K-V_{l,K}^{\rm geom}),
\tag{6.1}
\]

\[
M_i=\rho_l|\omega_i^*\cap\Omega_l|
+\rho_g|\omega_i^*\cap\Omega_g|.
\tag{6.2}
\]

上述量由当前 **相分辨 Hermite 几何积分**生成。不是 `rho(c_eps)` 或 `rho(center)`。动量 Hodge 及相分辨黏性继续消费现有 H07b/H09a 的 phase quadrature；这避免为了构造 CSF 而重写材料骨架。

### 6.2 同一 native 面质量通量继续唯一共享

\[
F_{\rho,f}=\rho_lQ_{l,f}+\rho_gQ_{g,f},
\qquad Q_{l,f}+Q_{g,f}=A_fu_f.
\tag{6.3}
\]

其时空积分、dual 面转移及 AMR 共享应由既有 H08a/H09/H10 接口完成。**禁止用** `rho(c_eps) u_f A_f` **生成第二套 convection mass flux**。

`c_eps` 只用于 CSF，故本方法**不能自动修复 Pure FirstJet 的材料体积漂移**。当前显式余量继续记录：

\[
R_K=
\frac{m_K^{n+1}-m_K^n}{\Delta t}
+\sum_{f\subset\partial K}s_{Kf}\bar F_{\rho,f}.
\tag{6.4}
\]

同理，动量 dual 上保留实际 `R_i^dual`，不以 primal 的简单插值冒充它。

### 6.3 `uR` 与能量记录不得删除

在动量守恒形式中保留

\[
\partial_t(Mu)+\mathcal N(F_\rho,u)
=\cdots+uR.
\tag{6.5}
\]

正确匹配时，均匀速度即使 `R != 0` 也不因不一致质量输运产生伪加速度。动能账本中单独列出

\[
P_R=\frac12\sum_i |u_i|^2R_i^{dual}
\tag{6.6}
\]

及其时间中点版本、物理边界质量通量项。不得将 `R` 并入“毛细力耗散”或“压力功残差”，也不得将 `R=0` 写为 V1 验收硬约束。

### 6.4 动量对流、黏性不可改变

保持目前 H10/H11 的共享 `F_rho` KEP/skew 对流，局部耗散 `D_num>=0`；保持

\[
K_\mu=E_h^TH_\mu E_h,\quad u^T K_\mu u\ge 0,
\quad \mu_l=\mu_g=0\Rightarrow K_\mu=0.
\tag{6.7}
\]

如果引入新的色函数或曲率数据导致对流/黏性数组身份不匹配，应拒绝并修正 source identity，不应 silently fallback 到 cell-average 密度、donor-cell 或人工黏性。

## 7. 单体压力–速度系统与时间推进接口

### 7.1 半离散目标方程

在同一 dual mass / actual convection / current geometry 来源下，目标为

\[
\boxed{
\frac{d}{dt}(H_{\rho,u}u)
+\mathcal N_{\mathrm{KEP/lim}}(F_\rho,u)
+K_\mu u+C p
=S_\sigma+S_E+F_b+uR_{\rm dual},
\qquad D_Iu=0.
}
\tag{7.1}
\]

`S_sigma` 由 §5 定义。该式是预期的**物理/数学接口**，并非声称当前仓库已经有完整 H14 production advance。

### 7.2 建议 midpoint / Picard 框架（待 H14 单独完成）

对于冻结对流和材料等中点数据，一步可以写成

\[
\frac{H_{\rho,u}^{n+1}u^{n+1}
-H_{\rho,u}^n u^n}{\Delta t}
+N^{n+1/2}+K_\mu^{n+1/2}\bar u
+C p^{n+1/2}
=S_\sigma^{n+1/2}+F_b^{n+1/2}
+\bar uR^{n+1/2},
\tag{7.2}
\]

\[
D_Iu^{n+1}=0,
\qquad\bar u=(u^{n+1}+u^n)/2.
\tag{7.3}
\]

若对流在 Picard 外迭代冻结，形成的最终**参考**单体系统仍为

\[
\boxed{
\begin{bmatrix}
A_u&C_f\\C_f^T&0
\end{bmatrix}
\begin{bmatrix}u_f^{n+1}\\p^{n+1/2}\end{bmatrix}
=\begin{bmatrix}r_u+S_\sigma^{n+1/2}\\0\end{bmatrix},
\quad A_u=\frac{H_{\rho,u}^{n+1}}{\Delta t}+\frac12K_\mu^{n+1/2}+K_{bc}.
}
\tag{7.4}
\]

式 (7.4) 的 RHS 必须同时包含正确的 `H^n u^n/dt`、已冻结对流、`uR`、体力/边界项；**不能为了保留对称矩阵而把真正的 material/convection 残差丢掉**。若将对流或毛细力做全隐式线性化，矩阵可能不再对称；不可继续无条件声称可使用 MINRES。

### 7.3 一个研究时间步的冻结/更新顺序

1. 取得接受的 `geometry^n`、`u^n`、`p^n` 与完整 identity；从 H07/H09 准备 sharp 几何质量。
2. 用现有 actual `CchAdvectorReconstruction` 做**同源**半步 Hermite 预测，得到 `geometry^(n+1/2)`；若该接口未准备好，限于**研究驱动**，不得伪称 H14。
3. 对同一中点 publication 独立派生 `color_csf`，计算 `C color` 和 `kappa_f`，形成 `S_sigma`；**所有对象同一个中点 clock**。
4. 使用既有唯一时空材料 `F_rho`、dual mass 与 `R_dual`，冻结 KEP/有限耗散项，装配单体系统。
5. 解 `u^(n+1),p^(n+1/2)`；以 `u_mid` 更新 Hermite 几何，必要时进行有限次数 Picard 校正，绝不另建一张“仅用于 CSF 的几何”。
6. 对最终实际接受的 begin/mid/end publication 重新计算并记录 `E_sigma`、`P_sigma`、`R_dual`、物理体积变化、压力功、黏性耗散、动能。超过残差/接受门槛则拒绝整个时间步而不是接受部分状态。

这只是未来 H14 的拟议算法顺序。**第一版 Codex 任务只实现静态 prepared capillary force 与其隔离研究驱动，不得一次性合并完整时序开发。**

### 7.4 毛细时间步的限制

CSF 是强界面力，即使对流和黏性稳定，显式/半显式毛细处理一般仍有尺度

\[
\Delta t_\sigma=O\!\left(\sqrt{\frac{(\rho_l+\rho_g)h_{\min}^3}{\sigma}}\right).
\tag{7.5}
\]

第一版必须进行 `dt` 减半和显式毛细限制扫描。**不能声称**仅改成 CSF 就可以解除表面张力刚性；未来全隐式界面 Jacobian 是独立研究任务。对固定网格低黏/无黏的研究驱动要报告 CFL、capillary Courant 和每一阶段候选几何身份。

## 8. AMR、轴对称、接触线与边界载荷

### 8.1 真实子面而非 virtual coarse face

对于 2:1 AMR，相邻 pressure leaves 的共同几何面可拆为 `f1,f2`。本方法按每个**真实 native subface**计算

\[
(S_\sigma)_{f_j}=\sigma\kappa_{f_j}A_{f_j}
(c_{R_j}^\epsilon-c_L^\epsilon),\quad j=1,2.
\tag{8.1}
\]

不能先在虚拟粗面上计算力，再人为按面积均分给两个细面。必须保持原有面方向/自由度身份及 `C=-D_I^T`，`c_K` 由**各自真实 cell**体积平均产生。常 `kappa0` 时，AMR 的每个自由行都应满足与压力同源的 `C` 平衡。

连续力带与不连续 AMR 尺度之间会有额外的截断误差，因此必须独立进行界面水平移动/旋转扫掠测试和 regrid 前后同源 identity 测试，不能只测固定静滴中心位置。

### 8.2 RZ 轴的正确处理

RZ 的 `c_K` 用 `2*pi*r` 积分，原生径向/轴向 `A_f` 及 `C` 均使用真实旋转体测度。轴面 `r=0` 的径向原生面积为零，不人工设一个无限大的 `kappa/r` 或 unit normal 纠正。

如果 Hermite 界面与轴交汇，使用出版几何的正则总曲率；正则化距离在轴附近的计算必须满足反射对称 `d(-r,z)=d(r,z)`，确保 `c_eps` 是偶函数。对于尖端曲率不能分辨的区域，**拒绝并要求加密/重新设计**，不能通过 cap `kappa` 或改小 `sigma` 制造伪稳定。

### 8.3 固壁、固定接触线和外边界

- **固定壁面**：沿用 H04/H12 的速度消元、压力 gauge 和 reaction 结构。固定 wall 处 `C c_eps` 可能有非零载荷，属于壁面 reaction；不能把它填给自由速度行，也不能假定外壁 `C1=0`。
- **pinned contact**：不得改变 FirstJet 的固定接触位置或任意强加接触角；力带碰到壁面时，要记录消失于求解子空间的壁反力。若需要动态接触角、润湿能和接触线功，属于另一个物理模型，V1 不纳入。
- **open/traction**：H15 尚未完成的开放边界不因新增 CSF 而自动支持。只对目前已明确验收的闭合或受控 Dirichlet 研究域声明可用。
- **selected/disconnected domain**：一个 publication 可以拥有更大的界面图谱，color/curvature 只在当前 Hydro selected physical domain 积分；不要借孔洞、断开分量之外的界面补足力支撑。

## 9. 电喷雾中 Maxwell 牵引的后续兼容接口

V1 首先验证 `sigma` 毛细力且 `E=0`，不把未知电学误差混入 H13。将来接入电学模块时，连续一流体形式可以在 Maxwell 张量分片光滑、体电荷模型明确的条件下考虑

\[
\mathbf f_E=\nabla\cdot T_E
\quad\text{（包含相应的界面分布应力）}.
\tag{9.1}
\]

或者对无体电力的相内模型，由两侧 Maxwell traction 差构造 CSF 力层；**二者只能选一种物理计入方式，不能重复施加。**

为预留实现接口，建议将界面电牵引分为法向和切向：

\[
t_{E,n}=\mathbf n\cdot[\![T_E]\!]\mathbf n,
\qquad
\mathbf t_{E,t}=(I-\mathbf n\mathbf n^T)[\![T_E]\!]\mathbf n.
\tag{9.2}
\]

法向电牵引可按与 CSF 相同的 color gradient 及 face interpolation 进入 `S_{E,n}`，**但其相对于 `Cc` 的符号必须通过解析均匀场压力跳跃测试决定**；不能直接把 `t_{E,n}` 数值加到 `sigma*kappa` 而忽略当前法向/跳跃定义。切向项在同一 `delta_eps` 力带中投影到 MAC 速度空间：

\[
(S_{E,t})_i\approx
\int_\Omega \delta_\epsilon(d_h)\,
\mathbf t_{E,t}^{\rm ext}(x)\cdot\mathcal R_i(x)\,dV.
\tag{9.3}
\]

其中 `mathcal R_i` 是与实际 native 速度值和运动学接口一致的重构测试函数，`t_{E,t}^{ext}` 是由同一 Hermite 界面 quadrature 获得的局部延拓；归一化和单位检查须保证力/功量纲正确。该式仅为未来设计接口，**不是已验收的 EHD traction deposition**。

## 10. 能量、压力功和几何守恒诊断

### 10.1 严格代数项与非严格项分开

已建立的代数关系：

\[
u^TCp=-(D_Iu)^Tp,
\quad u^TK_\mu u\ge0,
\quad \mathcal D_{\rm num}\ge0,
\tag{10.1}
\]

新 V1 应满足的常曲率 BF：

\[
C(\sigma\kappa_0 c_\epsilon)-S_\sigma=0.
\tag{10.2}
\]

但完整变量曲率毛细功

\[
P_\sigma=u^TS_\sigma
\tag{10.3}
\]

与实际 Hermite 面积能的变化率

\[
\dot E_\sigma=\sigma\frac{d}{dt}A_\Gamma(\phi_h)
\tag{10.4}
\]

**通常不具备代数恒等式** `P_sigma+dot E_sigma=0`；其误差需要单独记录。

### 10.2 完整收支残差的报告方式

在无外部功、封闭域的一个研究时步，定义

\[
K^n=\frac12(u^n)^TH_\rho^nu^n,
\qquad E_\sigma^n=\sigma A_\Gamma(\phi_h^n).
\tag{10.5}
\]

一个通用诊断余额为

\[
\boxed{
\mathcal R_E=
\frac{(K^{n+1}+E_\sigma^{n+1})-(K^n+E_\sigma^n)}{\Delta t}
+\mathcal D_\mu+\mathcal D_{\rm num}
-P_R-P_{\rm boundary}-P_E-P_b.
}
\tag{10.6}
\]

这里 `P_R` 的离散计算必须与真实 H09/H10 moment identity 对齐；`D_mu`、`D_num`、`P_boundary`、`P_E` 采用同一时间积分和符号。`R_E` 可能包含时间截断误差、CSF 功不相容、Hermite 面积更新误差、材料非守恒和求解器误差；不得提前把它解释成单独的“毛细误差”。

建议同时输出**独立的毛细–表面能缺陷**：

\[
\mathcal R_{\sigma}=P_\sigma+
\frac{E_\sigma^{n+1}-E_\sigma^n}{\Delta t}.
\tag{10.7}
\]

由于 `P_sigma` 应使用中点速度与中点力，其与面积增量所代表的时间位置须一致。`R_sigma` 若不随 `dt,h` 收敛，是阻止 BF-CSF V1 进入 production 的证据，而不是通过“总动能看起来稳定”可以忽略的误差。

### 10.3 压力、毛细力与惯性的分开导出

所有静态/动态算例至少记录：

- `||C p - S_sigma||`（静态）、`||C p - S_sigma + A u - r_other||`（动态）；
- `max_K |D_Iu| / local_flux_scale`、真体积压力均值；
- `K, E_sigma, D_mu, D_num, P_R, P_sigma, P_boundary`；
- `A_Gamma`, `V_l`, `sum_K R_K`, `max|R_dual|`, `||c_eps-alpha_sharp||`；
- `epsilon/h`, `sigma`, `rho_l/rho_g`, `mu_l/mu_g`, `max|u|`, 归一化 `Ca_spurious`；
- `kappa` 的局部极值/矩统计、`support_count`, AMR regrid 次数与时间层身份。

## 11. 与 LAMP 当前源码的具体映射

本段只是实施映射，不构成允许单次任务修改全部文件的授权。

| 项目已有对象 | V1 使用方式 | 禁止事项 |
|:---|:---|:---|
| `Interface` Pure FirstJet publication | 几何的唯一真值，提供 `phi / grad / curvature / surface quadrature / identity / clock` | 不新增第二条界面或存储型 color 输运 |
| `CchPhaseGeometry` | 保留 `alpha_sharp`、原 phase volumes 和相区域 | 不用 `c_eps` 替换其质量业务 |
| `CchPressureLoad` | `C=-D_I^T` 的唯一入口；`C(c_eps)` | 不另建 CSF 梯度 stencil |
| `CchRegularAdvectorGeometry` / `CchSpatialAdvector` | 后续时间积分消费同一 actual velocity/Jacobian | 不独立插值另一套界面速度 |
| `CchMaterialTransport` / `CchDualMaterial` | 唯一质量通量和真实 `R/R_dual` | 不以正则化 color 重建第二套质量流率 |
| `CchMomentumGeometry` | 原动量 dual mass/metric | 不凭 color 修正密度质量 |
| `CchKepConvection` | 原 KEP/skew 与 controlled dissipation | 不回退 donor/添加 viscosity floor |
| `CchEnergyViscosity` / `CchPhaseViscosity` | `E^T H_mu E` 原分相黏性 | 不用几何力代替应力张量 |
| `CchSinglePhaseSaddle` | reference saddle solve + 真实 gauge、闭域边界 | 不加第二套压力梯度或放松残差门槛 |
| `CchInterfaceNormalLoad` | H12/H13 历史基线，只做 A/B | 不在其现有 `capillary()` 内静默替换旧方法 |

### 11.1 建议新增职责（类型名供 Codex 审核后确定）

1. `CchRegularizedColor`：从**已发布同源几何**派生 `c_eps[K]`，同时保留 `epsilon`、几何身份、积分报告；**不持有自主演化状态**。
2. `CchBalancedCsfLoad`：持有同一 `CchPressureLoad`、`CchRegularizedColor`、curvature/source identity；输出 native `momentumRhsLoad`，并暴露 `C(c)`、`kappa_f`、support metrics。源/边界身份错误时拒绝。
3. 研究 driver：与冻结 H13/现有 H12 比较，**不创建生产 App 或时间步 acceptance state**。

若代码布局规则要求复用特定 namespace 或文件名，请先查询 `docs/README.md` 及现有 CMake 约束再确定；不要直接按上述建议名字构造绕过现有 owner 的捷径。

### 11.2 最小 C++ 端口示意（非可直接编译代码）

```cpp
// All objects must refer to ONE native mesh/metrics/FirstJet publication/time.
auto pressure = CchPressureLoad::build(phase.metrics().divergence());
auto color = CchRegularizedColor::build(mesh, phase, epsConfig);
auto csf = CchBalancedCsfLoad::build(
    mesh, pressure, phase, color, surfaceTension, curvatureConfig);

csf.validateFor(mesh, pressure,
    phase.geometryIdentity(), phase.timeLayer());

// The load is already the physical RHS of A*u + C*p = rhs.
const auto& S_sigma = csf.momentumRhsLoad();
for (int i = 0; i < freeDofCount; ++i)
    rhs[i] += S_sigma[i];

// Hard checks in test and debug builds:
// C = -D_I.transpose();
// sigma==0 -> S_sigma exactly zero;
// constant kappa0 -> S_sigma == pressure.apply(sigma*kappa0*color);
// c in [0,1], no hidden clipping/second F_rho;
// full identity / RZ / AMR constraints.
```

注意 `CchPressureLoad::apply()` 的真实签名和 `CchNativeValues` 绑定，必须从现有 H12 工厂代码复制**接口用法**，而不是由示意代码推测可直接调用向量 overload。Codex 每个 task 需要用实际头文件/CTests 证明 identity 一致。

## 12. 离散算法摘要

**Algorithm 1：LAMP BF-CSF 的一次 capillary evaluation（只读，非时步）**

**输入**：native mesh、同一 `CchPhaseGeometry`/`CchPressureLoad`、接受的 FirstJet snapshot、`sigma>=0`、`epsilon` 配置。

1. 核验 `mesh + native layout + pressure operator owner + geometry publication + clock`。
2. 在同一 Q4 几何上计算有量纲的局部 signed-distance proxy `d_h`，不改 `phi_h`；退化区域早拒绝。
3. 对每个 selected pressure leaf，用真实 XY/RZ 物理体积求积得到 `color_csf[K]=V_K^{-1}\int H_eps(-d_h) dV`，并检查 boundedness/纯相常值。
4. 从原 FirstJet surface quadrature，以真实 `w_q` 和局部正权 `W_f` 形成面曲率 `kappa_f`。所有有效 CSF 面必须有对应的真实曲率支撑。
5. 只调用**唯一** `CchPressureLoad` 计算 `g_csf=C color_csf`。
6. 形成物理 RHS `S_sigma[f]=sigma*kappa_f*g_csf[f]`；`sigma=0` 时准确零，constant-curvature regression 要逐自由行比较 `C(sigma*kappa0*c)`。
7. 输出 `S_sigma` 与不影响解的派生诊断；对象不更新 material、face `F_rho`、pressure、velocity 或 Interface accepted state。

**Algorithm 2：研究型时间推进**：完全按照 §7.3，先重用 H13 已存在的真实 `native/source/Hermite` 运动研究路径；仅当稳定通过静态、短窗、全周期及谱验收后，才另开 H14 生产集成任务。

## 13. 验证矩阵与验收指标：先算子，再静态，最后动态

所有阈值应先从**已有 H12/H13 的冻结 gate**继承；下面新列的目标若旧 gate 未明确，则作为**待预先冻结的目标**，不可算完以后按结果调整。验收必须同时保存原 H13 A/B，不能只报新方案较好的一组。

| Gate | 物理/数学测试 | 实际必须输出 | 通过条件 |
|:---|:---|:---|:---|
| BF-V0 | 现有基础回归 | H01–H12 中相关 CTest、原失败结果索引 | 数值及来源不变，架构门槛不退步 |
| BF-V1 | color kernel / 几何积分 | `H_eps`/`delta_eps`、XY/RZ 常值、圆/球、`c in [0,1]`、pure exact | 有界/纯相常值，GL 收敛，真测度一致 |
| BF-V2 | 相同来源与压力伴随 | `mesh/layout/clock`、随机 `u,p`、AMR 子面 | `u^T C p+(D_Iu)^T p` 接近舍入尺度；错误来源拒绝 |
| BF-V3 | 零 sigma、常曲率 | `S_sigma=0`、`S_sigma - C(sigma*kappa0*c)` | `sigma=0` 精确零；free-row 常应力平衡相对误差接近舍入 |
| BF-V4 | 单相/纯相与 planar interface | 远离界面速度、支撑面、零曲率、相流率 | 纯相 CSF 零，无缺失支撑/壁面泄漏 |
| BF-V5 | 非圆曲线曲率映射 | Ellipse / sinusoid / varying kappa，XY/RZ & AMR | 不采用解析曲率替代真实 Q4；曲率与载荷误差随 `h` 收敛 |
| BF-V6 | 72 种静滴几何 A/B | XY/RZ × uniform/AMR × 三分辨率 × 三中心偏移 × 两几何 | 原 H13 静态全门槛，不能保留原七组例外后仍判 PASS |
| BF-V7 | 低黏/无黏静滴 | `rho_l/rho_g=1,1e3,1e6`；`mu_l` 从有限到 0 | 伪流与加密趋势、鞍点残差均满足原验收；不能由 `mu_floor` 稳定 |
| BF-V8 | 无黏毛细波（完整周期） | `omega`、完整初值幅值、`K,E_sigma`、`R_sigma,R_E`、全部谱 | 保留原频率末档目标 `>=1.8` 与原能量收敛；不产生未解析增长模态 |
| BF-V9 | 时间步与 epsilon 扫描 | `dt, dt/2`, `epsilon/h=1.0,1.5,2.0` | 区分时间/空间误差；不得只挑一个最佳宽度 |
| BF-V10 | 高密度比平移/剪切 | `rho` up to `1e6`、`F_rho/R/uR`、常速度 | 不破坏原 H08–H11 所有质量动量一致性 |
| BF-V11 | AMR 运动、轴极限、contact | moving interface, regrid, `r=0`, pinned endpoint | 可重复、物理权重正确、完整身份、无额外边界伪功 |
| BF-V12 | 电喷雾前的 Maxwell 制造应力 | 定场电牵引正/切向、静态合应力 | 在毛细单独 gate PASS 后另做，不纳入本版首批实现 |

### 13.1 静滴测试：不接受只看 `max |u|` 的选择性报告

对于圆/球静滴，先**指定解析常曲率** `kappa0` 只测试代数 BF；然后使用真实 Hermite publication 的 `kappa_q`，比较：

1. 单元 color 与 `epsilon` 的响应；
2. 恒应力代数平衡及 fixed-wall reaction；
3. `||C p-S_sigma||`、实际 `max|u|`、`Ca_spurious`；
4. `rho_l/rho_g=1,1e3,1e6`、`mu_l=0` 与非零黏度；
5. XY/RZ、uniform/实际 AMR、不同 subcell interface offset。

特别强调：解析常曲率的机器平衡只验算了 `C c_eps` 的同源性；真实曲率驱动的伪流和变量几何仍是另一项验收。

### 13.2 动态波的完整物理检查

先使用两相、无黏、弱振幅、已知色散关系的毛细波。对无限深层可用

\[
\omega^2=\frac{\sigma k^3}{\rho_l+\rho_g};
\tag{13.1}
\]

有限深度必须使用与实际边界及两相深度对应的理论色散公式，不得拿无限深频率评估有壁面的有限深度算例。

测试必须至少包含：完整一个周期、三档以上网格、可变 `dt` 检查、真 `E_sigma`、真实 `H_rho` 动能、`P_R`、所有模态/完整初值与壁面反力功；**不得只报告首次过零频率**。如果能量误差在 `dt` 减半后不变或随 `h` 增大，则明确判定为尚未解决的空间/耦合问题。

### 13.3 建议的工程可执行阈值类别

- **代数恒等式**：以 `1e-12` 量级的双精度相对残差为**起始建议**；若已有 H12/H13 gate 更严格，则继承更严格者，禁止改为宽松阈值。
- **空间阶**：真实光滑标准算例以原 H13 `>=1.8` 目标作二阶门槛；不把某一组的末档假阳性当成全范围通过。
- **物理能量**：要求时空加密下误差逐渐下降，并完整报告量级；不事先宣称 CSF 能量机器守恒。
- **强界面/超高密度比**：允许研究阶段 FAIL，但必须保留全部原始 `rho=1e6,mu=0` 控制；不得选择性删除。

## 14. 理论性质与未解决风险

### 14.1 可由离散结构直接保证

1. 只读派生 color，不改变 Hermite accepted-state。
2. `sigma=0` 毛细力精确零。
3. `C=-D_I^T` 的压力功恒等式继续成立。
4. `kappa=constant` 时 `C(sigma*kappa*c)-S_sigma=0` 的自由行代数平衡。
5. `K_mu` 半正定及 `mu=0` 算子消失，前提是原对象身份/测试不变。
6. 原内部共享面通量的成对抵消；**不是** Hermite 几何相体积的逐 cell 严格守恒。

### 14.2 只能作为设计目标，必须数值验证

1. BF-CSF 在 XY/RZ、AMR 上真实变量曲率力的空间阶。
2. 高密度比静滴伪流随网格加密的收敛。
3. 无黏毛细波频率、完整初值响应和动能–表面能交换。
4. `epsilon=O(h)` 下对 Sharp 两相连续问题的收敛。
5. 运动学体积漂移与非守恒 `R` 的大小及收敛。
6. 接触线/尖端区域 CSF 展宽导致的几何与 EHD 应力偏差。
7. 真正的二阶完整时间推进（属于 H14）。

### 14.3 可预期的失效机制与应对

| 风险 | 原因 | 正确处理（不是数值补丁） |
|:---|:---|:---|
| 重新出现 H13 毛细波能量失配 | `S_sigma` 与真实 FirstJet 运动仍可能不共轭 | 使用完整周期及实际源方向导数诊断；必要时单独研究变分/界面功配对 |
| 重构位置导致 color 异常 | `phi` 不是 signed-distance 或梯度退化 | 用 `d_h` 标度/同源最近点距离；独立网格收敛与拒绝退化场 |
| `epsilon` 被粗网格欠分辨 | 狭窄力带积分和 AMR coarse–fine 跳变 | force-band 加密/正权分区积分/宽度敏感性；不靠调滤波 |
| `c_eps` 与 sharp `alpha` 不同 | 正则化力与真实分相材料测度不同 | 双账本和明确一致性误差；不得冒充质量守恒 |
| 极小液滴/尖端被过度展宽 | 有限厚度 CSF 的分辨率限制 | 要求 `epsilon << curvature radius` 与足够单元分辨，超界拒绝 |
| 接触线能量残余 | 界面端点力、壁反力和接触模型缺失 | 单独壁面/endpoint 功预算，固定 contact 语义不变 |
| Maxwell 正/切向力重复计入 | 既用体 Maxwell 散度又加界面 traction | EHD 阶段统一选择力的物理表示，并做解析制造校验 |

## 15. 可直接交给 Codex 的分阶段任务

### 15.0 严格执行规则

每次**只交给 Codex 一个 task ID**，每个任务保留：目标、允许修改、禁止修改、数学约束、实现、验收、失败规则、交付。遵守当前项目默认的 `<=5` 个 production 文件、`<=2` 个测试文件预算；超过先拆任务。状态只允许 `PASS / FAIL / BLOCKED`。**BF-00 是新方法选择的立法/文档变更，不得把原 H13 的失败改成通过。**

#### BF-00：方法切换确认与冻结基线（先做，0 生产文件）

**目标**：确定新 BF-CSF 与原 H13 的并行 A/B 关系；冻结输入、原动态失败与验收门槛。

- **允许修改**：现行数值方法索引、roadmap 中的**方法变更注记**及本设计/新任务文件。
- **禁止修改**：`src/`、`include/`、Interface 及 H13 已有结果与 FAIL 判定；不启动 H14。
- **数学约束**：承认 `C=-D_I^T`、Pure FirstJet、唯一 `F_rho`/真实 `R/uR`；明确 `c_eps` 仅是派生力场。
- **验收**：记录仓库身份/版本、现有 H12 与 H13 原报告、当前受影响 CTests、原 long-wave FAIL 和全套 gate，确认 `sigma>0` 物理 RHS 符号。
- **失败规则**：任何引用文件缺失/版本不明记 `BLOCKED`，不得以模型推测代替源码。
- **交付**：`BF-00_Report.md`、新 task 入口、现有状态索引引用；0 数值变化。

#### BF-01：Regularized color 的唯一只读派生

**目标**：从当前 FirstJet publication 生成 `color_csf[K]` 与 `epsilon`/求积报告。

- **允许修改**：一个 color `.h/.cpp` 类对与至多两个对应测试/注册文件；确需 CMake 配合先在预算内核对。
- **禁止修改**：Interface 产出、`CchPhaseGeometry` 材料质量、`CchMaterialTransport`、对流、鞍点和旧 H13。
- **数学约束**：式 (3.2)–(3.6)、真物理体积、boundedness、纯相精确常值、同一 geometry/clock；不独立推进 `c`。
- **实现**：实现 `d_h` 的有限/退化审计、正权体积分、RZ/XY 和 selected 域身份绑定。
- **验收**：解析 `H_eps/delta_eps`、平界面网格位置扫掠、圆/球、RZ near-axis、uniform/AMR 只读结果；旧质量/压力等测试不变。
- **失败规则**：不允许通过 clipping、强行以 `phi` 作 signed distance 或直接 cell-center 取值放宽误差。
- **交付**：真实 `c` 数据、求积收敛、测试日志与 `PASS/FAIL/BLOCKED`。

#### BF-02：C 共用的恒曲率 CSF 载荷

**目标**：实现 `g=C c_eps` 与 `S_sigma=sigma*kappa0*g`，建立准确符号和 BF 代数门槛。

- **允许修改**：一组新 CSF load `.h/.cpp` 和 1–2 份独立合同/注册测试。
- **禁止修改**：现有 H12/H13 工厂、`C`、`D_I`、`H_rho`、`F_rho`、`R`、KEP 或黏性。
- **数学约束**：`S_sigma` 是实际 RHS；`sigma=0` 精确零；对自由面 `C(sigma*kappa0*c)-S_sigma=0`；压力 gauge / wall reaction 正确。
- **实现**：copy H12 的 owner/identity 检查结构，不复制其负号定义。支持 XY/RZ、真 AMR 子面、压力分量和 selected 域。
- **验收**：随机强界面对齐、同一 convex drop `p_l>p_g`、静态实际鞍点、误配来源拒绝。
- **失败规则**：出现错误符号、不相容边界行或重建第二套 gradient 即 FAIL。
- **交付**：新的常曲率工厂/测试、完整 free-row residual、旧 H12 回归。

#### BF-03：真实 Hermite 曲率投影与一般 CSF 力

**目标**：将 §4 的 `kappa_q` 正权投影到 native force faces，形成真正的 variable-curvature V1 候选。

- **允许修改**：仅 BF-02 的新 CSF load 类及相关测试；超预算分任务，不改旧 H13。
- **禁止修改**：原 actual source 的插值/Jacobian 规则、FirstJet、质量/KEP/黏性、原 baseline。
- **数学约束**：常 `kappa` 精确再现、局部真实支撑、同一 `C c_eps`、真实 RZ/AMR 物理权重、零支撑拒绝。
- **验收**：真实 Hermite 圆/球/椭圆/正弦曲线、两种 epsilon、uniform/AMR、曲率与载荷独立积分复核。
- **失败规则**：不可拟合频率、默认过滤局部高曲率、扩远程支撑绕过缺失，或偷换解析曲率。
- **交付**：真实曲率误差/积分记录/载荷误差及 PASS/FAIL。

#### BF-04：静滴全范围与旧 H13 A/B

**目标**：在**实际** `CchSinglePhaseSaddle` 和 `CchPhaseViscosity` 上比较 CSF、旧 H12 常跳跃和原 H13 候选。

- **允许修改**：研究 driver/算例输入/基准测试与报告，不修改压力–速度/材料生产模块。
- **禁止修改**：旧 H13 数据、H14/H15、FirstJet、`R/uR`、boundary 契约。
- **数学约束**：实际 `Cp-S_sigma` 静态 residual，真实 `rho/mu`、正交力/壁 reaction，真实 null-space gauge。
- **验收**：完整原 H13 静态矩阵（72 geometries × densities × viscosities），特别保留原七组例外及 RZ/AMR。
- **失败规则**：某些旧失败不改善时不得以解析曲率、人工黏度或更宽 `epsilon` 替换正式门槛。
- **交付**：A/B 表、全部输入/输出、原始失败与新状态，不能提前宣称生产能力。

#### BF-05：无黏毛细波与能量–谱完整 gate

**目标**：只更换冻结旧 H13 研究驱动中的毛细 RHS 为 BF-03；保留原真实 native/source/FirstJet、质量/连续性算法和完全相同的初值。

- **允许修改**：隔离研究驱动、测试参数/结果/新报告；如需修改非毛细生产代码则先暂停并另立设计。
- **禁止修改**：现有 H14、实际 source/Jacobian、材料/动量算子、旧 gate 与已保存失败数据。
- **数学约束**：两相真实 `rho,mu`、`D_Iu=0`、真实 `E_sigma`、`P_sigma`、`R/uR`、完整壁面响应、实际共时几何。
- **验收**：完整一个周期、三至四网格、至少两 `dt`、所有设定 `epsilon/h`、密度比序列；末档频率原目标、完整初值全部谱、能量残差随加密下降、低黏极限无人工 floor。
- **失败规则**：即使频率达二阶，若能量/谱仍 FAIL 则**新方案 FAIL**；禁止仅报告最优 `epsilon` 或使用对照振子替代真实运动。
- **交付**：逐时刻物理账本、全谱/全初值、旧方案 A/B，是否允许进入 BF-06 的明确决定。

#### BF-06：AMR/RZ/接触线与 H14 接入许可（仅 BF-05 全 PASS 后）

**目标**：在真 AMR regrid、轴线、固定接触点及长时间实际界面运动中确认 BF-CSF 没有引入新的不相容性；决定是否冻结为 H13 的替代方法。

- **允许修改**：本模块必要的局部 AMR/边界力带支持和对应测试；任何 H15 职责先拆单任务。
- **禁止修改**：修改 Pure FirstJet 的物理语义、重写 pressure-gradient stencil、隐瞒 `R`、绕过已有 wall reaction。
- **验收**：XY/RZ 动态 AMR、尖端/轴线、pinned contact、断连域、全部报告和身份恢复；通过后出具给 H14 的**明确输入/输出/中点时序契约**。
- **失败规则**：AMR/H15 所需的高阶 metric 或 boundary 超出本任务时，显式标记 `BLOCKED` 并给出依赖项，不把研究驱动误称为 production。
- **交付**：方法冻结版、独立全 gate 报告、下一 H14 所需接口和遗留风险。

## 16. Codex 的最后验收清单（每次 review 直接核对）

- [ ] 用了一张真正同源的 FirstJet `phi_h`，而非偷建另一张界面。
- [ ] `color_csf` 只读；`alpha_sharp`、真 dual mass、`F_rho`、`R/uR` 不变。
- [ ] 不再将旧 `J_f C alpha` 改名冒充新增 CSF；确实实现了 `H_eps` 正则化物理体积 color。
- [ ] `C=-D_I^T` 保留，CSF 通过 **同一个 C** 计算 `C color`。
- [ ] 物理 RHS 符号正确：静态凸液滴 `p_l>p_g`，`C p=S_sigma`。
- [ ] `sigma=0` 精确零，常 `kappa` 每个自由行代数平衡。
- [ ] 每个支持面曲率来自同一 Hermite publication，RZ 不重复乘 `2*pi*r`。
- [ ] AMR 使用真实 fine subfaces，重建后不复用旧身份/clock 数据。
- [ ] `mu=0` 仍无黏；没有 curvature cap、minimum viscosity 或隐藏的 donor。
- [ ] 动态物理数据：`K+E_sigma`、`P_sigma`、`P_R` 和所有真实边界功均独立报告。
- [ ] 完整毛细波而非仅首次过零；历史 H13 FAIL 不被覆盖或伪称通过。
- [ ] 每个 task 以 `PASS / FAIL / BLOCKED` 结束，有可复算输入/源码/结果身份。

## 17. 参考文献与方法来源定位

以下是**外部公开文献提供的方法依据**；LAMP 的特定 `epsilon`、局部曲率投影、源码接口和任务分解是本文件的工程推演，**不是论文原文已有实现**。

1. Brackbill, J. U., Kothe, D. B., & Zemach, C., *A continuum method for modeling surface tension*, **J. Comput. Phys.** 100 (1992) 335–354. DOI: [10.1016/0021-9991(92)90240-Y](https://doi.org/10.1016/0021-9991(92)90240-Y). **用途**：连续有限厚度 color 与 CSF 基本模型。
2. François, M. M., et al., *A balanced-force algorithm for continuous and sharp interfacial surface tension models within a volume tracking framework*, **J. Comput. Phys.** 213 (2006) 141–173. DOI: [10.1016/j.jcp.2005.08.004](https://doi.org/10.1016/j.jcp.2005.08.004). **用途**：压力–表面张力离散平衡及连续/Sharp 表面力比较。
3. Popinet, S., *An accurate adaptive solver for surface-tension-driven interfacial flows*, **J. Comput. Phys.** 228 (2009) 5838–5866. DOI: [10.1016/j.jcp.2009.04.042](https://doi.org/10.1016/j.jcp.2009.04.042). **用途**：几何相分数 + balanced-force CSF + AMR；静滴和毛细波验证范式。
4. Nangia, N., Griffith, B. E., Patankar, N. A., & Bhalla, A. P. S., *A robust incompressible Navier–Stokes solver for high density ratio multiphase flows*, **J. Comput. Phys.** 390 (2019) 548–594. DOI: [10.1016/j.jcp.2019.03.042](https://doi.org/10.1016/j.jcp.2019.03.042). **用途**：高密度比、交错网格、质量–动量一致对流、monolithic saddle 与 well-balanced 力。
5. Valle, N., Trias, F. X., & Castro, J., *An energy-preserving level set method for multiphase flows*, **J. Comput. Phys.** 400 (2020) 108991. DOI: [10.1016/j.jcp.2019.108991](https://doi.org/10.1016/j.jcp.2019.108991). **用途**：表面能与动能的变分/离散交换；对 V1 §5.5 未来能源结构有参考意义。
6. Zeng, Y., et al., *A consistent adaptive level set framework for incompressible two-phase flows with high density ratios and high Reynolds numbers*, **J. Comput. Phys.** 478 (2023) 111971. DOI: [10.1016/j.jcp.2023.111971](https://doi.org/10.1016/j.jcp.2023.111971). **用途**：强材料对比时同源质量–动量输运与 AMR 收敛，**其守恒密度策略不直接等于当前 Pure FirstJet 非守恒路线**。
7. LAMP 当前项目 `docs/api/H03_Cch_Pressure_Load.md`、`H06_Cch_Energy_Viscosity.md`、`H07a_Cch_Phase_Geometry.md`、`H08a2r1_Cch_Characteristic_Source.md`、`H10k_Cch_Spatial_Material_Flux.md`、`H12_Cch_Constant_Normal_Stress.md`、`H13_Cch_Hermite_Capillary_Traction.md`，以及 `docs/task/H13_Hermite_Capillary_Traction_Report.md`。**用途**：现有源码行为、接口身份、已通过/未通过的验收边界。

---

**最终方法原则（供 Codex 保留在方法文档首页）**：

> 只在**同一 Hermite–GALS 几何**上派生正则化 CSF color；真正单流体两相惯性、唯一材料质量通量、真实 `R/uR`、KEP/skew、能量型黏性与兼容 MAC saddle 全部保留；毛细力使用与压力共用的 `C=-D_I^T` 形成 **balanced-force continuum surface force**。常曲率压力平衡是代数要求，变量曲率低黏/无黏动力学是必须经静滴、完整毛细波、能量、谱、AMR 和 RZ 验证的研究目标。新方案验收失败时，必须保留完整证据，而不是回退到单相近似、人工数值黏度、任意曲率滤波或伪造体积守恒。
