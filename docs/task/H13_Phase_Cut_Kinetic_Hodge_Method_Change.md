# H13 / 同源分相速度与局部物理动能Hodge联合方法变更

状态：2026-10-09 / 081136，**用户明确回复“授权”，已批准本联合方案**。现有H13唯一任务继续FAIL；授权允许实施，不代表候选已实现或验收，不启动H14/H15。

**当前修订授权（2026-10-09）：** 用户明确授权加入有初值和演化方程的独立切向模态，并将其动能、弱力、黏性和材料能量一并耦合。新增q属于同一source的组合速度状态；此明确授权替代下文“所有辅助系数仅由U静态确定、不携带独立速度状态”的相关限制。原native U/p、D/C、FirstJet、唯一F/R/uR及物理门槛保持；不再就本范围内正常实现重复申请。§8原单格反例的推论已由实际邻域rank复核纠正，见§9。

授权登记：本次用户答复直接对应上一轮H13 phase-aware唯一source、局部物理动能Hodge及对应弱力/材料/黏性能量配对提案，批准至多24生产文件、10份已有短契约和K_phys公开定义变更。原K_diag与历史FAIL继续单列，原全部频率、收敛、效率门槛保持；两项已验证的小维护补充纳入预算。无需再次申请本范围内的正常实现或修错。

## 1. 为什么需要改变冻结方法

[条件下界证据](H13_Phase_Cut_Kinetic_Metric_Certificate.md)完成36组原算子响应及36组真值投影。原native相质量没有算错，但跨相切向面均值的平方遗漏物理速度方差；半格n64动能范数漏约4.85%，渐近是一阶。在准确同源normal与伴随弱力的假设下，原对角M给出一阶mobility下界。因此继续只改四文件normal端点/载荷不能关闭这个算子族的缺陷。证据不等于一般新方法已被证明成功。

2026-10-08授权明确保留M和原native动能，并限四production文件/两测试；本方案改变这些离散性质，现已依据总纲§2.1获得额外方法授权。此前KEP两行算术精修、characteristic两组可信参考同步现纳入联合预算，不另开任务。

## 2. 唯一owner与局部物理场

全局未知仍为native正常面平均U及pressure p；Pure FirstJet的accepted状态、相体积来源及真实密度保持。phase-aware prepared source统一拥有几何、native稀疏basis和所有派生查询，不另设accepted速度或材料库存。

XY物理场按各相重构 `R_phys(U)`：以原native面平均的H(div) lift起步，使各cell散度为 `D*U/V`；在其上加入局部streamfunction的curl补空间，不改散度。补空间包括分相富集 `curl(phi_+*W)`，`phi_+=max(phi,0)` 只读当前唯一几何。它允许无黏切向跳跃，同时保持界面normal连续。共享face端点/均值约束、master/hanging约束与同相jet连续性一次性消元；不逐case切换公式，不使用解析压力/频率确定系数。prepared系数由真实面平均矩及受约束局部多项式再现确定，局部辅助系数不成为全球状态。

RZ采用真实r权重的streamfunction表达及轴正则补空间，2:1使用原parent/subface真测度。必须证明面均值再现、局部rank、normal连续、相内正则与D相容；本提案不把XY示意公式当成已完成RZ/AMR设计或测试。

几何推进需要C1延拓：同一source从同一准备系数提供 `R_ext`，在Γ与物理normal及其切向导数一致，单相局部与物理场一致；相内延拓/混合jet按声明的多项式约束唯一确定。原characteristic积分和Jacobian公式保持。**这一点明确改变原“所有查询均为一个全域C1物理速度”的假设**；物理分相场与几何C1延拓是同一owner的具名派生视图，消费者不得混用。不得以第二套隐式source或失败fallback掩盖不相容。

## 3. phase-cut局部动能Hodge

原 `M_diag` 保留为真实dual材料质量及历史动能诊断。改变的是运动方程惯性度量与物理动能，不改ρ、相体积或native均值的定义。

按唯一几何标记与Γ相交的momentum组件dual以及固定重构halo。不同同组件dual仍分割物理测度；非标记部分沿用原对角项，标记部分用同源分相basis积分：

`H_rho = sum_unmarked M_f*e_f*e_f^T + sum_marked_duals integral rho*R_d(q)^T*R_d(q) dV`。

`K_phys=1/2*U^T*H_rho*U`。常向量要恢复真实总材料质量；正权求积和均值可识别性必须给出free DOF正定性及有效条件数。这里不是按频率调质量、密度平滑或全局拟合。H局部稀疏，准备阶段静态消元，保留全球U/p数量；生产saddle使用H，不把H退回M作为fallback。非cut单相光滑区域必须与原基线保持规定精度和效率。

## 4. 压力、载荷、黏性与材料能量一起闭合

压力保持 `C=-D_I^T`、原壁面/开边界及压力nullspace。力仍由完整同源P2/FirstJet曲率和实际 `J=-sigma*kappa` 装配，完整常应力 `J_f(J0)=J0*I` 保持；Γ法向basis必须来自上述同一物理/几何normal。divergence-conforming场的积分散度恒等式与α、C、界面normal之间要用真实测度验证，不能再靠独立几何或平滑力消差。

动量时间离散使用 `H_mid*(U_new-U_old)/dt`；几何/密度使H随时间变化，必须把 `1/2*U^T*dot(H)*U` 的贡献与同源材料通量和真实连续性余项一起装配/诊断。对流的skew/KEP配对按H扩展，材料F仍唯一，共享面反号/AMR来源约束保持；R和uR显式保留，不强制R=0，也不偷偷恢复独立库存。

黏性项从同一分相物理梯度和真实μ建立非负双线性形式；界面/trace项如需要必须随真实μ消失，exact μ=0得到零算子，不加黏度地板、界面人工阻尼或默认upwind。量纲、常场零耗散与原受控耗散比较保持。上述是必要联合改造范围，不能只接H而沿用未经配对验证的旧能量项。

## 5. 文件预算与执行边界

请求同一H13的联合例外：**至多24个production文件、10份已有直接相关短契约**；新增局部detail文件计入24。生产仅限Hydrodynamics及其只读几何查询接线：source/SpatialAdvector、MomentumGeometry/SinglePhaseSaddle、NormalLoad、Phase/EnergyViscosity、MaterialTransport/DualMaterial、KepConvection及必要Hydro集成。不得更改Interface accepted字段或characteristic/Jacobian算法、电场、字典语义与无关模块。大规模研究仍在benchmarks和D:证据目录。

实际开工前在现有TODO登记精确文件清单；若24文件不足，先交付具体依赖理由，不能边实现边自动超预算。正常实现/修错在本批准方案内自主推进，不为每个局部选择重复询问。下一步只做同一H13的联合候选构造与原验收，不拆r11/source/Hodge独立任务。

## 6. 新算子验收与原物理门槛

- T0：标准构建无新warning，architecture/repository_layout、唯一owner/来源及native未知量数量。
- T1：phase内P2/P3制造场、normal连续/切向跳跃、准确nativeface均值、shared/master/2:1/RZ轴真测度、相跨越连续性；H的对称正定/常量总质量/真实分相动能独立积分；C负转置、完整常应力矩阵与gauge shift；实际FirstJet方向功、真实H变化与F/R/uR/对流/黏性配对。不能仅验证scalar平面参考。
- T2：本文36精确平面case采用原误差测度，要求偏移和高密度normal/响应至少达到原1.8阶；原72几何/1728静滴全范围及全部六组回升关闭；原r11/H08a相关短契约及28原窗口/空间/时间/R范围、H10b原84无黏场景、H10/H11对流耗散范围回归。
- T3：原完整实际毛细波初值、周期长度、三档空间/时间、频率误差与>=1.8门槛、strict μ=0/低μ/high-density/AMR/RZ、原效率要求全部保持。H13全范围通过以后才允许路线图下一任务。

**能量定义变更公开登记**：新主物理动能为K_phys；必须独立用相分割体积分证实其连续物理一致性与收敛，并验收物理动能加原表面能的真实收支。原K_diag、旧原能量gate与所有历史FAIL继续原样输出、单列比较，不用新名称改判旧FAIL，也不隐藏旧基线。批准本提案意味着授权这一离散能量定义变更，而不是授权放宽频率、收敛、质量漂移或效率门槛。

## 7. 可比性、风险及明确授权内容

旧源码/hash/配置/场数据/对角能量全部封存；新身份单列。相几何、材料量、native面均值、pressure、应力、初值/周期/误差测度可直接比较；K_phys与K_diag不同，必须同时报告，不声称逐位等价。

尚未解决的核心风险为局部富集rank与条件数、一般曲面的C1延拓与FirstJet方向功、RZ/AMR、移动几何H变化、无黏高密度效率。它们是批准后需要完成的工作；失败保留FAIL，无case特供或自动fallback。

请求授权：在**同一H13**按以上预算联合改造phase-aware唯一source、局部物理动能Hodge及对应弱力/材料/黏性能量配对，并公开采用K_phys。现有两份小维护补充可纳入本预算。若保持原M与原唯一对角能量冻结，则本提案不能实施，现有反例不能按本路线关闭。

授权依据：[开发总纲](../00_LAMP_Development_Charter.md)§2.1规定“不得由 Codex 在实现过程中自行改变方法”；[原H13授权](H13_Native_Flux_Port_Method_Change.md)§3明确“原常应力完整矩阵、原C负转置、材料M和原native动能保留”。因此本次需要的是具体数学方法决定，现有自动续跑不等于批准。

## 8. 2026-10-09 / 1030 source可识别性审计

授权后求积核的八组Hodge契约和完整momentum-geometry契约已通过，但实施phase-enriched source前做了一个必要的代数可识别性审计。方格内直界面 `x=3/8`，液/气速度分别取 `(0,5/8)` 与 `(0,-3/8)`。所有四个native normal face means都为0，和静止场相同；两相散度0、interface normal连续且为0，但相间切向跳跃为1。该P0场属于T1要求的P2/P3再现空间，可对应周期y方向的自由无黏涡片。精确观测矩阵及复算在[审计证据](../../runs/checks/H13/20261009_081136/phase_source_identifiability_audit.md)。

该精确反例只证明一个cut格的四面均值不足以识别其内部切向模式。103518接续检查发现实际prepared source使用邻格观测：四个XY/RZ、uniform/AMR组合均有19或23个轴向观测、局部矩阵rank1而实际支持矩阵rank2，常场相系数复现误差至多2.22e-16。因此“整个网格U相同”的推论不成立，不能据单格例否定已有邻域重构。当前C1 source不能表示跳跃仍是实现事实，但不等于观测不可能。独立切向模态的修订已由用户明确授权，不再等待方法答复。

原单格有理数结果保留；过宽的全source阻断结论撤回。未解析内部模式与可由邻域识别的相多项式分别验收，不降低T1。H13继续FAIL，H14/H15不启动。

## 9. 已授权独立切向模态的实现与能量配对

组合速度为 `x=(U,q)`，U仍是原native面均值，q是同一source拥有并随时间推进的切向模态。物理场 `R_phys(x)=R_base(U)+sum(q_m*T_m)`；每个T由同一Q4的分相streamfunction产生，界面法向迹为0、native面法向均值为0，故压力块为 `C_ext=(C,0)^T`，D不新增第二份。几何C1延拓取共享base normal；零normal的切向补场不修改界面normal及其沿界面的导数。

第一步实现native cut-cell内部模态：XY `T=curl(phi_+*B)/g`，B在cell四边按三阶零消失（保证velocity/first-jet在同相cell边界为0），g由真实界面gradient正尺度确定。RZ使用 `Psi=r^2*phi_+*B/(g*h_r)`、`T=(Psi_z/r,-Psi_r/r)`，r=0格的radial B采用偶多项式，轴迹用解析极限。该内部子空间用于未观测的独立内部运动；它不能替代一般切向跳跃的共享trace空间，后者及phase-P2/P3再现仍为最终T1必需，不能用bubble-only结果收缩目标。

真实动能为 `K=1/2*x^T H_rho*x`，H含U-U、U-q、q-q完整分相Gram块。q初值必须显式提供，不能由U=0推断q=0。固定几何的方程为 `H*(x_new-x_old)/dt + K_mu*x_mid + C_ext*p = f_mid`；所有free native与q行共同求解。纯法向牵引对T的功严格为0，由curl零normal恒等式证明；体力/切向牵引采用同一T弱测量。黏性从同一phase first-jet构成非负Gram，有限mu的跳跃trace一致性仍需完成，mu=0精确零。

移动几何必须输出 `H`变化、端点时间余项及唯一材料F/真实R/uR的完整组合能量配对；q的cut变化/AMR转移必须以物理场投影并记录能量余项，不能丢弃旧模态或只改名旧K_diag。仅固定几何的代数演化通过不能宣称完整H13。现有24生产/10契约范围内优先实现上述owner、Hodge和saddle块；若确切接线依赖超出该范围，先提交具体必要文件与责任理由，不以正常实现选择为由反复询问。
