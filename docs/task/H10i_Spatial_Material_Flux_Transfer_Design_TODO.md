# H10i / H08–H09 同源空间材料通量重设计

状态：**FAIL（同源几何查询契约不完整）/ engineering PASS**，2026-10-04；本轮已结束。[报告](H10i_Spatial_Material_Flux_Transfer_Report.md)记录64布局的常规来源核验和4个实际重合面反例；后者使真实dual质量与候选相段查询均被拒绝。下一唯一任务[H10i1 重合相迹/裁剪契约设计](H10i1_Coincident_Phase_Trace_Clip_Contract_TODO.md)。已核验 H10h 的115 payload/28当前交付、1,042源文件/346输入；本轮0 production/0 test/1 research，仅做被动观察。H10未关闭，H11/H08a等待。

目标：给唯一材料通量 owner 定义相分割感知的 primal/dual 空间求值，消除旧体积分配把 cut pressure 材料散度强行输入纯气体 dual 区域的问题。只改变这一空间材料转移事实；不在同一任务实施新时间积分、完整流体时步或动量算子调参。

允许修改：本设计/TODO/报告、当前耦合设计/API 状态说明及现有 `cases/benchmarks/cch_convection/driver.cpp` 的被动证据观察；0 production / 0 test / 1 research 文件。实现前单独写设计，冻结后续实现任务的 typed API、文件清单和不超过 5 production / 2 test 的预算；本任务不修改生产材料定义。临时分析归 `runs/checks/H10i/<时间戳>`，关闭证据归 `runs/baselines`。

禁止修改：Pure FirstJet/Hermite–GALS primary、几何发布、库存/bubble/逐 cell 体积拟合、材料质量修复、第二套 momentum `rho*u` 通量、原 D_I/C、冻结 PDF、原 App/字典/研究门槛。无 AMR 特例、算例名称分支、解析数据进入求解、增加全局高阶状态或隐式 fallback。旧材料实现继续冻结，不能在它上面静默改变字段语义。

## 仅四项材料耦合硬约束

1. **唯一共享材料通量。** 一个 owner 从同一实际几何和 advector 定义 `F_rho = rho_l Q_l + rho_g Q_g`；primal、dual 和 momentum 消费该对象在各真实面的求值。不同几何面允许有不同积分值；同一几何面只有一份值。不用旧 primal 标量 F 的面积分配代替相内空间积分，也不独立插值另一份 momentum 密度通量。
2. **内部面反号共享 / AMR 子面守恒。** 每条真实共享段只计算一次，两邻域取相反外法向号；粗面值由实际子面求和。这是面通量的代数抵消，不能推论非守恒 GALS 的下一层几何质量必须满足 FV 更新。
3. **来源一致。** 相区域、物理 XY/RZ 测度、实际 native advector、geometry/layout identity、完整 clock/候选和 AMR/regrid 来源相同。真实质量仍直接测量唯一 FirstJet 几何；不以通量反推质量或修正界面。
4. **KEP/skew 结构。** 动量使用该 owner 的实际 `q_dual` 和 `R_dual`，保持材料动能恒等式、显式边界功和既有 `R/uR` 一致性。常速度和能量核验须使用修改后的实际输入，不能继续套用旧映射缺陷。

压力伴随、真实几何质量、非负黏性等其他模块现行约束保持；本授权只解除材料空间转移的附加过约束。

## 不再强制的关系与必须报告的差异

几何 quarter 体积权重、tensor 分配、`q_dual = T q_primal`、几何质量等于 `T m_primal`、`R_dual = T R_primal + G/dt` 及逐 cell `R=0` 均不作为新材料求解的等式约束。它们作为历史实现/诊断比较保留；精度、R 的收敛和系数尺度仍需实际验收，不能用解除等式来免除 PDE 精度要求。

固定 mesh 上保留几何映射 T 仅供诊断，定义

\[
E_{sp}=q_{dual}-Tq_{primal},\qquad
G=\Delta m_{dual}-T\Delta m_{primal},
\]

\[
R_{dual}=TR_{primal}+G/\Delta t+E_{sp}
=\Delta m_{dual}/\Delta t+q_{dual}.
\]

新硬检查是上述**含 E_sp 的定义自洽**，不是 `E_sp=0`。不同速度分量的实际 dual 质量和缺陷分别核算；regrid 后重建实际来源。半离散对流继续使用 `N=diag(q_dual)/2+K`、`K^T=-K`；常速度一致性和材料动能预算从真实 F/R 推导，不依赖旧 quarter 交换关系。

## 实现前必须解决的设计问题

- 同一材料 owner 如何取得真实 primal 子段与 dual 内部段的液/气区间及正物理求积，如何保证交集和 AMR 子面的共享归属；禁止用体积分数按面积分摊 cut face。
- 实际 native advector 的局部 compatible 重构如何保持原 D_I 的运动学意义、常数/affine 再现和轴正则；是否需要一般 one-ring patch。旧 H07 面常值材料率及旧 primal F 不额外冻结为新空间求值必须复制的硬数值，不另立材料通量 owner。
- 在同一源上先定义相通量，再推导 dual q/R 和 final convection 的空间矩目标。不能先任意修正目标再寻找输入，也不能把纯气体行的大 q/M 当作合理 transport 尺度。
- 明确哪些矩是运动学再现/精度条件，哪些只作为残差诊断；不同时强制旧 quarter 交换、所有 cut 行精确 affine 与非守恒材料 KE 的不相容条件。非均匀材料须给出实际二阶余项及质量误差范数，不能把有限网格下界比值当成普遍不可能证明。
- 定义来源拒绝、所有权/生命周期、所有 native/轴/端口行、重建成本及后续短契约；全局只保留原 native 速度与压力 DOF。

验收：T0 当前 protected 源码/输入/可执行程序/PDF 哈希和架构/目录通过，确认只做设计与被动观察。T1 完整推导四硬约束及上式缺陷分解，独立核验实际共享段/子面/XY–RZ 物理测度和来源，不靠旧 commuting identity 证明正确。设计必须列出后续实现的固定数字门槛：原恒等式 1e−12 gross、常速度 1e−10、内部 affine 先精确/二阶通过，再原 H10b 空间阶 ≥1.7、最细 L2≤0.05、Linf 下降及全部原 KE/时间/donor gate。

材料验证计划必须包含 H10h 的实际纯气体反例（XY 2:1、rho=1e6、n16、dof103：三层液体体积和质量变化均为零，旧 q=R=12499.9875）、n4/8/16/32 四网格组合、rho=1/1e3/1e6、移动/crossing/regrid。纯相常速度闭合误差按**该相的实际通量尺度**核验，不以邻近高密度材料的 gross 掩盖泄漏；混相 R、材料空间矩和导数行尺度的收敛在实现前定目标，不能要求严格 R=0。全局共享抵消不能代替局部相准确性。

失败规则：设计不闭合、缺少同源空间数据或与四硬约束冲突均 FAIL 并拆出明确接口任务；只有真正外部依赖才 BLOCKED。不得以提高 H10h 的固定预算系数、替换范数或豁免原研究 gate 收口。前置不通过不运行原完整 H10b，也不进入 H11/H08a。

交付：单一 owner 的数学/API 设计、旧→新性质对照、真实反例及独立诊断证据、固定的后续实现文件清单/预算/数值门槛。原 H08/H09 PASS 是原范围的历史验收，H10g/h FAIL 保持；本任务的设计 PASS 也不能宣称完整材料 PDE 或流体模块通过。
