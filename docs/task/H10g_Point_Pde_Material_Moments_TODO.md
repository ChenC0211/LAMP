# H10g / 点采样PDE材料矩与最终有界对流

状态：**FAIL，诊断关闭**，2026-10-04。真实质量范数目标在两个材料组合的n16→32空间阶不达原1.7，候选退役；最终QP未实施。48组必要下界及原材料矩来源审计通过，仅保留被动观察器。启动独立校验H10f的148个冻结payload及当前18份交付，保护源码1,042份、输入346份；见`runs/baselines/H10g_20261004`。[实现前设计](H10g_Point_Pde_Material_Moments_Design.md)、[报告](H10g_Point_Pde_Material_Moments_Report.md)。下一单任务[H10h](H10h_Final_Operator_Joint_Feasibility_TODO.md)直接联合最终pair与PDE矩可行性，H10未关闭，H11/H08a等待。

目标：最终native convection从实际点采样PDE与真实质量定义出发，联合constant/affine空间矩、material-KE/skew及有界二次/三次余项。不能继续将flux-centroid几何矩当成PDE目标，也不能仅换权重使一条功恒等式通过。均匀材料的中心通量已有H10f可行证书，应保留其正确低阶/有界行为。

允许修改：启动时独立核验[H10f报告](H10f_Compatible_Bounded_Moments_Report.md)与冻结链，先写数学设计。预算仍为CchKepConvection.h/.cpp、现有短契约、现有研究driver/README及对应文档，2 production +1 test +1 research；临时定位runs/checks/H10g/时间戳，关闭证据runs/baselines。生产只保留通过的一条路线。

禁止修改：Interface、H08/H09质量/F/R定义、旧App/输入/PDF、CMake/通用工具。无库存/bubble、人工黏度、donor fallback、问题名/AMR分支、解析状态或载荷修正，不删native边界/轴行。不提前H11/H08a；若必须一般Hodge/材料矩重设计，拆出明确数学任务，不偷偷改变KE身份。

数学约束：`a=(-N+u*q)/M`，C要求`K*1=q/2`。PDE一次目标应由同一真实质量与实际native advector生成，而不是默认几何b/M正确。显式单周期端口的横向坐标z至少需两个skew相容条件：`sum(b+z*q)=0`、`sum(z*b+q*z²/2)=0`。相容修正须同时量化实际PDE余项；H10f数据已证明仅消除几何Psi不足。唯一F/q/R/uR与原材料KE保持。

实现要求：

- 保留H10f均匀n=16原s预算的明确可行系数，不再把其迭代失败解释为交集不可行；一般中心KEP可作为最终最小改动锚点。
- 先以真实质量定义实际PDE误差范数，推导材料矩相容修正。研究候选可从最小化`sum M*(delta_a/h²)²`出发；其等式约束最小改动权重为`M*h⁴`，是二阶PDE余项范数的推导，不是继续调整失败惩罚指数。先证明Gram可识别、实际误差尺度与阶，再考虑生产。该候选目前未实施或验收。
- 若需系数尺度S包含真实M*V/h，明确它只是算子预算，不是另一材料F；预算须从PDE Taylor展开在试验前固定。不能提高H10f固定C_aff=8来重标其失败结果；新设计须给出独立数学依据并保留原裁决。
- 一般one-ring上的共享反对称pair联合满足最终矩；约束最小改动须有可行性/正确工作集及KKT证据。局部segment不足统一one-ring，不能独立拟合再负转置配对。系数、二次/三次矩与实际apply必须独立审计。
- 实际速度来源为原H07 `totalFaceRates/faceMeasures`及可识别局部functional，绑定同一geometry/clock；无解析速度输入、持久高阶速度或额外主状态。周期/native等值拓扑显式，非法关系拒绝，无旧路线fallback。

验收：T0标准Release、来源/架构/目录、原App/输入/PDF保持。T1旧124场景、1,516内部矩、真实非零R/uR/跨界/regrid、原1e−12载荷/KE与1e−10常速度保持；新增真实最终端口/材料矩/拒绝契约。T2先原内部affine反例n=8/16/32四组合，再全部native行与rho=1/1000/1e6真实材料PDE余项的h阶、行系数尺度及高阶矩；不把48组观察器exit 0当作PDE PASS。T3前置全部通过才原H10b28窗口/15,116步/84gate完整重跑，整体截断/终值阶≥1.7、最细L2≤0.05、Linf下降、原KE/时间门槛及4个donor控制保留。

失败规则：目标不对应PDE、相容性/系数/余项/身份/完整研究失败均FAIL并退役候选，不是BLOCKED。不能用几何矩或材料功PASS替代空间一致性，不提前进入H11。

交付：实现前PDE材料矩设计、唯一通过的最终operator API及独立契约（若有）、可行/不可行及全部行/材料/内部反例证据、原完整研究与明确裁决、冻结与退役档案。本轮仅诊断/目标候选，不增加生产API；原完整H10仍FAIL。
