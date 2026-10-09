# LAMP 开发总纲

**版本：** 2026-09-30 / V2 收口版  
**状态：** 现行开发总纲（normative）  
**适用范围：** Interface、Hydrodynamics、Electricity 及其 EHD 耦合。Plasma/PIC 暂不纳入本轮改造。

**2026-10-02 用户授权变更：** 当前 Interface 主线转为纯 Hermite–GALS（非守恒），见[设计变更](task/I12_Pure_Hermite_Gals_Design_Change.md)。与独立 M、体积闭合、共享守恒相通量有关的下述旧要求仅适用于历史守恒比较；新路线保留 first-jet、同源几何、曲率/AMR 准确性，体积漂移必须报告。Hydro/电场要求不因此自动放宽。

## 1. 项目目标

**2026-10-03 用户授权：** Interface/Hydro 耦合采用[完全非守恒材料输运](task/H07_Pure_Gals_Hydro_Coupling_Design.md)。材料质量、密度、黏度直接由当前唯一 GALS 几何分相积分取得，不另设守恒材料库存或反向修正界面。质量漂移、材料连续性缺陷及其动量/动能项必须显式诊断并验收；下文守恒材料输运要求仅供历史比较。压力相容、黏性非负、电荷/Gauss 约束保留。

**2026-10-04 用户授权：** H08/H09 重设计仅将唯一共享材料通量、内部面反号共享/AMR 子面守恒、来源一致、KEP/skew 结构作为材料耦合硬约束。旧 quarter 分配及 pressure→dual 精确交换改为诊断，不再约束新材料求解；R 不要求严格为零，仍须核验精度、收敛和稳定性。Pure FirstJet 独立，无库存/bubble；详见[现行耦合设计](task/H07_Pure_Gals_Hydro_Coupling_Design.md)和[H10i](task/H10i_Spatial_Material_Flux_Transfer_Design_TODO.md)。其他模块硬约束保留，历史失败不改判。

**2026-10-05 用户授权：** r11 为 common regular source/characteristic 基础设施最后关口；除非出现可复现且直接阻断 H08a 的数学反例，不再拆 source/Jacobian/低秩精修任务，不新增硬约束、特殊 case 或放宽门槛。通过后立即回到 [H08a 同源时空相通量积分](task/H08a_SpaceTime_Material_Flux_Upgrade_TODO.md)，以 R 成为可收敛的小余项为目标；保持 Pure FirstJet 独立、无 M_K/bubble、无 R=0 或 quarter commuting 强制。

LAMP 当前阶段只围绕一个目标收敛：**可靠、准确、高效地模拟电喷雾液体界面的强变形、射流、颈缩与破碎，并能够在低黏乃至无黏条件下稳定运行。**

本轮开发采用三份冻结方法作为数学主线：

1. [Conservative AMR Hermite–GALS](methods/QUANTUM_Conservative_AMR_Hermite_GALS_Method.pdf)：保留 first-jet/AMR/同源几何基础；相体积守恒部分按 2026-10-02 授权转纯路线；
2. [CCH–EHD](methods/QUANTUM_CCH-EHD_Method.pdf)：流体、质量—动量一致输运、压力—速度相容离散与电—毛细载荷；
3. [SI-HHO-Q](methods/SI_HHO_Q_QUANTUM_Electricity.pdf)：sharp-interface HHO 电势、守恒电荷与 Maxwell 应力。

本轮**不做**：Plasma/PIC、3D、hp-adaptivity、额外湍流模型、离散能量 shape derivative Maxwell force、非线性介电模型，以及任何与上述主目标无关的功能扩展。

## 2. 文档权威关系

**最新BF-02裁决（2026-10-09）：PASS（仅本任务）**。[TODO](task/BF-02_TODO.md) / [报告](task/BF-02_Report.md)：同一C的恒曲率实际RHS、真实静态自由面/gauge/壁反力完成，200次saddle、2/2新契约+10/10回归通过；3 production/2测试及完整cell/face数据归档。BF-01已封存PASS，下一次仅执行BF-03，本轮BF-03/H14/H15未启动。当前真实源码/输入新A/B基线、旧H13 FAIL/缺失/timeout及全部物理gate保持，sharp材料/唯一F/真实R-uR/C-D和旧消费者保持。下列旧BF-00入口为历史检查点。

**2026-10-09 当前BF-CSF路线入口：** 用户更新目标，按 [BF-CSF任务书](roadmap/LAMP_BF_CSF_Codex_Tasks.md) 和 [BF-CSF V1](roadmap/LAMP_BF_CSF_Method_V1.md) 逐任务实施。当前仅 [BF-00](task/BF-00_TODO.md)：0 production变更，确认方法与冻结基线；原payload丢失状态及待审查恢复方案见 [报告](task/BF-00_Report.md)。旧H13独立切向模态联合改造保留为独立未采用候选，不继续其开发或混入BF-CSF。新方法仅在同一Pure FirstJet几何上派生只读color和新的毛细RHS；native U/p、C/D、sharp材料/唯一F/真实R-uR、对流/黏性保持。原H13 FAIL和历史授权保留为各自旧范围；BF-00未完成前不进入BF-01，H14/H15不启动。

**2026-10-09 用户修订授权（同一H13独立切向模态）：** 用户明确授权加入有初值和演化方程的独立切向模态，并将其动能、弱力、黏性和材料能量一并耦合。组合速度状态为同一owner的 `x=(U,q)`；此明确授权覆盖旧“所有辅助系数由U静态消元、全球仅U/p”的相关限制。q必须显式初始化和推进，完整H含U-q交叉项，压力扩展为 `(C,0)^T`；唯一材料F、真实R/uR、严格mu=0及全部原物理gate保持。详细方法和执行范围见[方法修订](task/H13_Phase_Cut_Kinetic_Hodge_Method_Change.md)§9与[现有TODO](task/H13_Weak_Traction_Acceptance_TODO.md)。内部bubble子空间不能替代共享trace/P2-P3验收；仍仅H13，FAIL，不启动H14/H15。

**2026-10-09 用户授权（同一H13联合动能改造）：** 用户明确回复“授权”，批准[分相source与局部物理动能Hodge方案](task/H13_Phase_Cut_Kinetic_Hodge_Method_Change.md)。范围至多24生产文件、10份已有直接短契约；允许改变旧对角惯性/唯一native动能及全域C1物理场假设，采用同一owner的分相物理场与C1几何延拓、K_phys，并配对弱力/材料/黏性。原M_diag作为真实材料质量及旧能量诊断保留，原C/D、accepted FirstJet状态与characteristic/Jacobian公式不变；F唯一、R/uR显式、严格mu=0与全部频率/收敛/效率门槛保持。两项已验证KEP算术及reference维护纳入该预算。H13仍FAIL，授权不代表验收，不允许H14/H15或无关改造；下述2026-10-08旧冻结范围由本次明确范围替代。

**2026-10-08 用户授权（H13联合改造）：** 用户明确同意[共同native面通量端口与P2弱牵引提案](task/H13_Native_Flux_Port_Method_Change.md)，允许在同一H13内修改此前冻结的common source master顶点normal值/切向导数准备，并与同源完整P2弱牵引联合验证。预算不超过四个production文件、两份已有短契约；原M/C、characteristic/Jacobian积分、F/R/uR、KEP、FirstJet状态与全部原物理门槛保持。该授权不重开独立r11精修任务，不改判历史结果，不代表候选通过或允许H14。

### 2.1 数学规范

数学方法以“三份冻结方法 + [数值硬约束](01_LAMP_Numerical_Invariants.md)”为准。任何实现如果需要偏离，必须先形成单独的设计变更说明，解释：

- 为什么原方法无法实现或验证；
- 改动改变了哪条离散性质；
- 新的算子级验收是什么；
- 旧结果是否仍可比较。

不得由 Codex 在实现过程中自行改变方法。

### 2.2 工程规范

工程实现以本总纲、[架构规范](02_LAMP_Architecture_Rules.md)、[字典规范](03_LAMP_Dictionary_Rules.md)和[输出诊断规范](04_LAMP_Diagnostics_Output_Rules.md)为准。

### 2.3 执行顺序

实际开发以[分模块任务清单](roadmap/05_LAMP_Module_Task_Roadmap.md)为唯一任务入口。旧 `TODO`、阶段报告和性能试验记录保留用于追溯，但**不再自动产生新任务**。

## 3. 六项总体验收目标

| 目标 | 要求 |
| --- | --- |
| 鲁棒性 | 高密度比、强表面张力、2:1 AMR、低黏/无黏都必须是正常工作状态，不靠“最小黏度”、全局 upwind 或隐式 fallback 存活。 |
| 准确性 | 先通过算子级恒等式和 MMS/解析收敛，再讨论完整电喷雾。曲率、质量、压力梯度、电通量不得由互不一致的几何/插值链路产生。 |
| 守恒性 | Interface 与材料输运接受非守恒；质量漂移、材料连续性缺陷及动量/动能一致性项显式核算并收敛验证。动量使用唯一同源材料通量，电荷/Gauss closure 保留守恒要求。 |
| 效率 | 全局未知量只保留方法真正需要的自由度。高阶能力优先局部重构、局部 Hodge、静态消元和 prepared operator；禁止为了局部高阶精度引入持续膨胀的全局高阶状态。 |
| 可维护性 | 一种物理量一个 owner、一种 production 方法一条主线；无隐式 fallback、无重复 stencil、无数值内核直接读字典或写结果文件。 |
| 可观测性 | 每个接受步都能从控制台和 `output/series` 判断守恒、残差、时间步、迭代数、能量与耗时；每次运行能还原实际配置和程序版本。 |

## 4. 当前方法收口决定

### 4.1 Interface

当前生产候选目标为 **纯 AMR Hermite–GALS（非守恒）**：

- 直接输运 first jet：`phi + grad(phi)`；
- `phi_xy` 是紧凑 Hermite closure，不是独立二阶动力学状态；
- 不再携带独立守恒 cell phase volume；体积从唯一几何测量，作为误差诊断；
- Q4 内部值直接采用 characteristic 预测，不做守恒体积闭合；
- AMR 保留 exact restriction / master constraints；纯路线不推进相通量或执行相库存 reflux。

旧独立界面路线已按 G02 退役并保存快照；现有 Hydro/EHD 必需适配移入 Compatibility，禁止新依赖或作为失败后备。

### 4.2 Hydrodynamics

生产目标固定为 **native staggered CCH**：

- 速度在 native face，压力在 cell；
- 拓扑散度 `D_I` 唯一；压力载荷由 `C = -D_I^T` 生成；
- AMR 精度进入局部 compatible Hodge/reconstruction，而不是增加第二套 pressure-gradient stencil；
- 材料来自同一非守恒 GALS 几何；动量使用唯一材料通量并显式处理材料连续性缺陷，无独立守恒质量库存；
- 对流采用 KEP/中心主项 + 局部受控耗散；
- 黏性采用能量型 `K_mu = E^T H_mu E`，`mu=0` 时严格消失；
- 时间推进采用 midpoint predictor–corrector。

当前大型 `CompatibleMoment/cchMomentMidpoint` 路径保留为实验 baseline，**停止继续扩大全局高阶速度状态**。

### 4.3 Electricity

Electricity V2 直接新建 **SI-HHO-Q / sharpHHO2** 主线，不从 `sharpHermiteP3` 逐步打补丁演化：

- `P2` cell/face data；
- `P3` local potential reconstruction；
- cell DOF 局部静态消元；
- 全局只保留 skeleton/interface trace；
- equilibrated flux 同时服务 Gauss、电流和 Maxwell trace；
- surface charge 使用 BDF2–Schur。

`sharpHermiteP3` 仅作 A/B 回归。

## 5. 开发纪律

1. **一次只推进一个模块。** 顺序为 Interface → Hydro → Interface/Hydro → Electricity → EHD。
2. **一次任务只改变一个数学事实。** 不在同一任务同时做方法重构、性能优化、字典迁移和多算例调参。
3. 一次任务原则上修改 `<= 5` 个 production 文件和 `<= 2` 个测试文件；超出时先拆任务。
4. 每个任务结束只有 `PASS / FAIL / BLOCKED` 三种状态。`BLOCKED` 时说明缺失接口，不准顺手加 fallback。
5. 新方法未通过 T1/T2 验证前，不用于 Taylor cone 或电喷雾调参。
6. 无调用者的失败路线按用户授权退役并归档；仍有现行消费者的兼容代码单独隔离、冻结依赖，不能成为隐式后备。
7. 性能优化不得改变离散方程、验收容限或实际物理参数；若必须改变，视为新的数值任务。
8. 任何“跑通一次”的结果都不是方法完成证据。必须同时满足对应 invariant、基准和收敛验收。

## 6. 验证分层

| 层级 | 内容 | 触发频率 |
| --- | --- | --- |
| T0 | 编译、静态接口、配置解析、结构检查 | 每个任务 |
| T1 | operator/unit：代数恒等式、局部重构、守恒闭合 | 每个任务 |
| T2 | module benchmark：单模块短算例 | 一个小阶段结束 |
| T3 | convergence + AMR：网格/时间收敛、粗细接口 | 一个模块 milestone |
| T4 | coupled production：静滴、Taylor cone、jet breakup | 大阶段结束 |

T0/T1 未通过时，禁止用 T4 结果“反向证明”实现正确。

## 7. Definition of Done

一个模块只有在以下条件同时满足时，才可声明进入 production：

- 冻结方法已完整落地，没有隐藏的旧 stencil/fallback；
- 对应的数值硬约束全部通过；
- 关键解析/MMS/AMR 验证通过，并记录误差和收敛阶；
- 字典只暴露受支持选项，启动时输出实际生效配置；
- 控制台与结果文件包含规定诊断；
- 性能没有出现结构性自由度/几何复杂度失控；
- 对旧 baseline 的 A/B 差异有说明；
- 文档和测试与代码同一次任务更新。

完成上述条件后，才允许进入下一模块。

2026-10-02 工程约束补充：[目录与文档](05_LAMP_Repository_Rules.md)、[代码与测试](06_LAMP_Code_And_Test_Rules.md)、[使用与调试](07_LAMP_Run_And_Debug_Rules.md)。程序名统一为 LAMP（Levelsetbased Adaptive Multiphysics Platform）；历史证据保留原名称。
