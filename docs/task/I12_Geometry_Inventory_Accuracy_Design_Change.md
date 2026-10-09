# I12 / 几何精度与独立守恒库存分离

**状态：用户已授权，2026-10-01；数值实现分项推进，I12 / Gate I 仍 FAIL。** 用户明确要求：“如果严格每个cell体积守恒和曲率精度很难同时满足，那就 放弃每个 cut cell 的 Hermite 几何必须机器精度地匹配守恒体积”。依据[总纲 §2.1](../00_LAMP_Development_Charter.md)，本文记录这一数学变更，不重复请求批准。既有程序路线、状态与 owner 优先，不改冻结 PDF。

## 变更原因与证据范围

[r12 条件性试验](I12a11r12_Corner_Cut_Conditioning_Report.md) 的独立解析积分表明：小相靠近角点时，固定边界 first jets 的 Q4 center bubble 体积灵敏度按 cut 线性尺度五次方退化。同一 `128*epsilon*capacity` 目标差，在最弱样例导致约 `0.00369` 的 center 改动和约 `0.945` 的 center Hessian 改动，而边界 C1 数据不变。这证明局部病态，不能据此声称所有保守高阶方法均不可能。

[r11 实际差额分解](I12a11r11_Weak_Cut_Closure_Gap_Audit_Report.md)、[r9 曲率阶段审计](I12a11r9_Translation_Curvature_Stage_Audit_Report.md) 和 [r10 pullback 审计](I12a11r10_Translation_Pullback_Audit_Report.md) 表明，极小几何/库存差会参与下一步重构并影响曲率；同时重构本身也有误差，不能把全部问题归因于最后一次 bubble。[r13 时间控制](I12a11r13_Translation_Time_Step_Control_Report.md) 中 dt 减半使同终点曲率误差增大约 44 倍，继续减小 dt 又在 t=.025 后拒绝，未证明简单缩步可解此问题。

因此停止把更紧的逐 cut 几何体积闭合当成曲率精度的代理指标。r14 精确标量细化与 [r15 受限库存转移](I12a11r15_Conservative_Refinement_Inventory_Report.md) 保留为已验证能力；后续不再以加密恢复机器精度匹配为先决条件。

## 改变的离散性质

对容量 `C_K` 定义 `e_K=V_K(phi_h)-M_K`。`M_K` 始终为唯一独立守恒相体积，仍由唯一共享区间通量更新：

`M_K^{n+1}=M_K^n-sum_f(s_Kf*F_l,f)`。

内部面贡献成对抵消；coarse flux 是实际 fine subflux 之和；实际存储库存严格满足 `0<=M_K<=C_K`。不允许库存 clipping、解析重置或用重新积分的 `V_K` 覆盖 M。

取消 I-R2 中把一维 bubble 方程强制求到库存 roundoff 预算的要求。bubble 变为局部、受几何精度约束的体积校正：预测几何已满足独立几何误差预算时保持不变；需要修正时，仅在既有拓扑/位移证书和额外导数扰动预算的交集中减小几何—库存差。不能为了更小体积残差越过导数预算。两种预算无法同时满足时明确拒绝并报告，不能伪报接受。

唯一几何、first-jet 直接输运、派生 mixed、master/hanging 和 B/J coherence 保持。法向、曲率、cut 与相面测量仍来自该几何。`e_K` 是同一状态的派生诊断，不是第二套界面或额外自由输运变量。局部库存守恒残差与 `e_K` 必须分别输出。

## 几何误差与校正预算

几何误差使用具名 `abs(e_K)/C_K` 门槛，独立于守恒 roundoff 和求积估计。第一项实现沿用现有用户请求的 `relativeVolumeTolerance`（默认 `1e-10`）作为几何门槛，停止暗中压到 `16/64*epsilon`；这个默认值已有来源，但它本身不保证曲率收敛。原 I12 接口位置、梯度、曲率及整体体积误差门槛不变。

对当前 cardinal bubble `B=256*xi^2*(1-xi)^2*eta^2*(1-eta)^2`，修正 `delta_c` 有全单元上界：

- `abs(delta_phi)<=abs(delta_c)`；
- `abs(delta_phi_x)<=4*abs(delta_c)/h_x`，y 类似；
- `abs(delta_phi_xx)<=32*abs(delta_c)/h_x^2`，yy 类似；
- `abs(delta_phi_xy)<=16*abs(delta_c)/(h_x*h_y)`。

这些是保守多项式界，来自一维 `16*t^2*(1-t)^2` 的导数；不是经验曲率滤波。可由调用者逐单元提供具名物理 Hessian 分量预算，换算成 center 修正区间，并与原 Core 证书取交集。预算必须显式给出并报告；不从解析真解倒推，不随失败临时放宽。该界只控制新增 scalar 扰动；曲率还含梯度分母及 RZ 轴向项，必须维持既有正规性检查并通过实际解析收敛验收，不能把 Hessian 界冒充完整曲率误差证明。后续实际预算选取须通过固定配置的 h/dt 矩阵验证。

第一项算子仅放宽 cut cell 的几何匹配；纯相、未解析成分、拓扑拒绝不自动放宽。完整运动中可能出现几何纯相与微量库存不一致，其可接受误差和跨步传播必须在库存协调任务中明确处理，不能在此文中默许。

## 库存有界性是接线前提

历史上旧 M 与几何液量约 `5e-12` 的相对差曾在单元排空时留下真实负库存；现有 `HermiteGalsInventoryAssembly` 只允许原 roundoff 内的共享通量修复。因此单独降低闭合强度再直接切换完整驱动是不完整的实现。

保持现有 raw 几何时空通量、边界历史与身份检查。后续必须建立有来源、守恒、有界且误差受控的共享通量处理：raw/final 通量及两者之差分别可查，边界与总通量约束不被隐藏改写；内部校正仍只存一次并用于库存和 reflux。不能仅把现有 roundoff 常数放大到几何容限，也不能因失败静默切换 donor/upwind。若需要新的离散限制公式，先补充其可行性、精度和 T1 说明，再实施；本文不提前宣称具体尚未实现的公式已获验证。

## 实施顺序与验收

1. r16：本设计变更、规范例外与入口更新；代码保持不变。
2. r17：局部几何精度预算闭合算子。共用原 I10 几何/拓扑/求根内核，显式区分新契约和归档对照；无隐式 fallback。T1 覆盖弱 cut 保留、真正校正、有限导数预算、无可行交集拒绝、XY/RZ 与 AMR、M/identity 不变。既有旧路径保持显式 A/B 可比。尚不切换完整驱动。
3. 库存协调：基于实际 `M-V` 缺陷，完成共享区间通量的有界性与精度分析/实现，必要时分拆算子任务。测试排空/充满、正反相、内部抵消、边界、counterflow、AMR reflux、失败隔离；原守恒门槛不变。
4. 接线与验收：半步/全步使用同一新契约；发布分别检查库存与几何误差，记录校正限制与 raw/final 通量。重跑原 polynomial、非多项式 h/dt、rotation/vortex/AMR crossing/曲率矩阵。只有全部满足原 gate 才宣告 I12 通过。

每项先 TODO，遵守 <=5 production + <=2 test-related 的预算；算子 PASS 不替代完整推进 PASS。不因本次授权提前进入 Hydro 或改变 Electric/耦合方法。

## 可比性与诊断

冻结文件、原方法 PDF、历史可执行及归档保留。A/B 使用相同初态、mesh、dt、速度、终止时间与原误差门槛，变化的闭合契约和精度预算显式列出。至少记录：库存局部/全局 ledger、几何 `e_K` 的有符号和/L1/Linf、修正前后 center、导数扰动界、受限单元、求积估计、raw/final 通量差、reflux 与耗时。预算、诊断与物理守恒验收不混用；Report 不参与数值分支。
