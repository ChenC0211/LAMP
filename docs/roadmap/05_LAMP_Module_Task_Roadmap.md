# LAMP 分模块任务清单

**版本：** 2026-10-09（同一H13接续）
**用途：** 后续给 Codex 的唯一任务入口。每次只发一个 task ID；完成验收后才能进入下一项。

**最新BF-02裁决（2026-10-09）：PASS（仅恒曲率load任务）**。[TODO](../task/BF-02_TODO.md) / [报告](../task/BF-02_Report.md)：新两组契约86.7/9.6秒、200次实际saddle及10/10回归，完整cell/face与来源hash归档。下一次仅执行BF-03；本轮BF-03/H14/H15未启动。BF-01已封存PASS，旧H13 FAIL/缺失/180秒timeout与全部后续gate保持。下面BF-00 BLOCKED为历史检查点。

**2026-10-09 方法入口变更：** 用户已更新当前目标，按 [LAMP_BF_CSF_Codex_Tasks](LAMP_BF_CSF_Codex_Tasks.md) 和 [LAMP_BF_CSF_Method_V1](LAMP_BF_CSF_Method_V1.md) 逐BF-ID实施。当前唯一任务为 **BF-00**，[TODO](../task/BF-00_TODO.md) / [报告](../task/BF-00_Report.md)。原long-wave/H12历史payload丢失，用户已确认；按BF-00原规则暂BLOCKED，当前源码/输入及现存证据冻结和回归继续完成，恢复方案待决定。BF-01尚未启动；旧H13保持FAIL、切向模态候选及以下授权/检查点保留为旧路线记录，不继续其方法开发。H14/H15不启动。

**当前授权（2026-10-09）：** 用户已批准同一H13的[分相source/局部物理动能Hodge联合方案](../task/H13_Phase_Cut_Kinetic_Hodge_Method_Change.md)，并明确授权独立切向模态的初值、演化及动能/弱力/黏性/材料能量耦合。[现有TODO](../task/H13_Weak_Traction_Acceptance_TODO.md)登记24生产/10已有契约上限及全部原验收；本范围内继续实施，无待授权阻断。候选仍FAIL，H14/H15未启动。以下待授权标签为历史检查点。

**最新103518（2026-10-09）：** 仍仅H13，**FAIL / 已授权候选实施中**。实际邻域观测rank复核纠正了1030单格反例的过宽推论。当前内部切向模态已显式初始化并接入完整[U,q] Hodge、压力saddle、体力弱功、bulk黏性和唯一phase/material atom通量；移动H的交叉能量、真实R及端点时间余项显式输出。新增专项及完整momentum-geometry契约通过，完整契约101.74秒。内部bubble只是第一子空间；共享trace/P2-P3再现、有限mu界面一致性、完整KEP/uR动量推进及模态转移尚未闭合，原T2/T3未通过。详见[执行TODO](../task/H13_Weak_Traction_Acceptance_TODO.md)，H14/H15不启动。

**授权前045122检查点（2026-10-09）：** [动能条件下界](../task/H13_Phase_Cut_Kinetic_Metric_Certificate.md)完成36条件响应+36真值投影；准确normal/伴随力和原M的假设族存在一阶反馈下界。之后用户已明确授权[联合方法方案](../task/H13_Phase_Cut_Kinetic_Hodge_Method_Change.md)；此历史检查点不改写原失败，不启动H14/H15。其余较早标签属于历史检查点。

**最新035647（2026-10-09）：** H13仍为唯一任务，FAIL。[精确平面直接反例](../task/H13_Exact_Plane_Cut_Response_Counterexample.md)已独立核验：一般cut的normal插值约一阶，rho=1e6/n64 middle cut的回复模态比=-10.4379；n8约1秒最小复现。当前normal面矩/弱牵引联合候选不能采用，H14不启动。先推导共同source一般内部normal与原真实M/C/应力相容的修正；不得自动扩大现行授权。189次原分项解、36组MMS及原静滴/周期失败证据全部保留。[报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-09--035647原载荷分项与精确平面cut反例)区分本轮与以下旧检查点。

**2026-10-09接续：** 仍只做H13。[执行TODO](../task/H13_Weak_Traction_Acceptance_TODO.md)及[003903报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-09--003903原系数与完整静滴周期审计)：隔离KEP准备精修与新source reference同步通过原完整契约，永久采用待第五production/第三已有测试两个明确预算例外；当前生产候选完整静滴72/1728完成，六组寄生速度末档回升，H13仍FAIL。n8新目录完整256步无黏周期自然完成、独立核算一致，单网格频率/能量改善不代表三档空间/时间/AMR/低黏与材料长回归通过；H14不启动。此前86步停止及以下检查点保持历史适用范围。

**当前 Hydro 状态（2026-10-07）：** [H13 同源Hermite毛细载荷](../task/H13_Hermite_Capillary_Traction_Report.md) **FAIL，仍为唯一当前任务**。原静滴七组伪流例外、无黏波能量/频率未收口；所有隔离载荷参考均未进入生产；Gram族及等密度1.557阶的简单加性候选均未采用；三个实际观察确认当前常权应力平均的信息损失；局部应力偏差及已知压力lifting均未采用；后者的等密度负刚度、0.716阶与能量升高已由原native/source分解复现；二相传输81组及九次FirstJet方向更新已完成；共享一圈应力矩候选末档频率阶-0.435、原能量峰值缺陷随加密增至0.452，独立矩阵指数/功核算复现，已拒绝；18组光滑传输及原实际算子复现确认native压力/加速度近二阶、最终source界面运动约一阶，固定二次应力profile预测后余量二阶以上；局部normal最小代价消差虽可达，但邻段扰动1.506倍，逐段叠加留下全部误差，已排除；完整P2压力矩差=tangent+jump/h已由原全矩阵核验，实际可达jump端口相对尺度约0.335且不随加密消失；仅改弱牵引虽频率二阶但能量缺陷约5.56；假设面通量矩匹配运动/牵引的线性参考频率1.987阶、原能量舍入闭合；具体共同source与P2弱牵引联合改造已获用户明确授权；真实source候选契约通过，完整弱牵引及全物理门槛待完成，候选未验收；一般曲面/AMR/实际动态及联合精度仍待解决。定位发现[直接阻断H08a的实际Jacobian反例](../task/H13_H08a_Transparent_Picard_Jacobian_Counterexample.md)，已按用户限定例外在一个现有文件中先代数消去C1内部边界相消项；25/25受影响检查通过，新版本完整H08a的28窗口及原空间/时间/R历史门槛通过；H13单周期复测仍失败。r11其他来源与原门槛保持，不另拆精修任务。历史H08a/H10/H11/H12规定范围结果保留；这一处限定正确性修复已收口。Pure FirstJet独立，无M_K/bubble、R=0或quarter commuting强制；唯一共享F、真实R/uR、KEP/skew保持。H14未启动。

| ID | 当前任务 | 验收 |
| --- | --- | --- |
| I12p1 | [纯状态/发布/原驱动](../task/I12p1_Pure_Hermite_Gals_Publication_TODO.md) | **PASS**：[报告](../task/I12p1_Pure_Hermite_Gals_Publication_Report.md)，12 原完整算例 / 1,036 步、28 CTest 通过；无独立 M、无体积闭合、唯一几何 |
| I12p2a | [纯路线收敛与步内 crossing](../task/I12p2a_FirstJet_Convergence_TODO.md) | **PASS**：[报告](../task/I12p2a_FirstJet_Convergence_Report.md)，空间/曲率/时间原阶门槛及 48 crossing 通过 |
| I12p2b | [纯路线动态 AMR](../task/I12p2b_FirstJet_AMR_Transfer_TODO.md) | **PASS**：[报告](../task/I12p2b_FirstJet_AMR_Transfer_Report.md)，16 exact/8 projection、12 粗化收敛组、12 动态完整组（1,036 步/36 regrid）、35 CTest 通过 |
| I12p3 | [首次纯路线 Gate I 审计](../task/I12p3_FirstJet_Gate_Audit_TODO.md) | **FAIL（历史）**：[报告](../task/I12p3_FirstJet_Gate_Audit_Report.md)，原主要缺口已由 p3a–p3d 关闭；当前裁决见 p3r2 |
| I12p3a | [纯 first-jet 时间步 producer](../task/I12p3a_FirstJet_Step_Producer_TODO.md) | **PASS**：[报告](../task/I12p3a_FirstJet_Step_Producer_Report.md)，16场景/80拒绝、36CTest、固定/动态24完整组（2,072步），原完整输出保持 |
| I12p3b | [纯路线独立 App/字典入口](../task/I12p3b_FirstJet_App_Entry_TODO.md) | **PASS**：[报告](../task/I12p3b_FirstJet_App_Entry_Report.md)，37 CTest、24 实际完整算例/2,072步/36 regrid；与原模块数值逐字相同 |
| I12p3c | [纯路线正式输出诊断](../task/I12p3c_FirstJet_Output_TODO.md) | **PASS**：[报告](../task/I12p3c_FirstJet_Output_Report.md)，37 CTest、18 频率对照、ParaView 实读、24 完整组/116帧；数值及几何identity保持 |
| I12p3d | [纯路线效率与规模验证](../task/I12p3d_FirstJet_Performance_Scale_TODO.md) | **PASS**：[报告](../task/I12p3d_FirstJet_Performance_Scale_Report.md)，37 CTest、24 完整组/2,072步、58 次规模运行；数值和 VTK 保持，随后复核 Gate I |
| I12p3r1 | [前次纯路线 Gate I 复核](../task/I12p3r1_FirstJet_Gate_Reaudit_TODO.md) | **FAIL**：[报告](../task/I12p3r1_FirstJet_Gate_Reaudit_Report.md)，数值/主要交付 PASS；两项诊断已由p3e/p3f关闭，当前裁决见p3r2 |
| I12p3e | [字典范围错误上下文](../task/I12p3e_FirstJet_Input_Diagnostics_TODO.md) | **PASS**：[报告](../task/I12p3e_FirstJet_Input_Diagnostics_Report.md)，15 精确拒绝/6 合法边界、原 App/输出契约；判据和合法数值保持 |
| I12p3f | [同源曲率时序诊断](../task/I12p3f_FirstJet_Curvature_Diagnostics_TODO.md) | **PASS**：[报告](../task/I12p3f_FirstJet_Curvature_Diagnostics_Report.md)，30实际组/150行、90频率实验统计、96VTK保持；原数值/输入/输出契约通过 |
| I12p3r2 | [当前最终 Gate I 复核](../task/I12p3r2_FirstJet_Gate_Final_TODO.md) | **PASS / Gate I PASS**：[报告](../task/I12p3r2_FirstJet_Gate_Final_Report.md)，11archive/11,438payload、当前身份及完整证据链通过；下一项H01 |

下述 I07–I10 的守恒专项与旧 I12 失败保留历史记录，不作为新纯路线必须恢复的前提。H07/H-R3 在新路线接入前须单独设计质量—动量耦合，不能沿用已失效的守恒相通量假设。

## 0. 任务执行格式

每个 task 必须包含：

```text
任务 ID / 标题
目标：只解决一个数学或架构事实
允许修改：文件/目录范围
禁止修改：明确不碰的模块
数学约束：必须保持的关系
实现要求：最小必要改动
验收：可执行的 T0/T1/T2 条件
失败规则：何时返回 BLOCKED
交付：修改文件 + 测试结果 + 未解决问题
```

状态只有：`PASS / FAIL / BLOCKED`。

默认变更预算：`<=5` 个 production 文件 + `<=2` 个测试文件。若超出，先拆 task。

---

## S. Phase 0：冻结与立法

本次文档包已完成 `S02–S07`；在开始数值代码前还需完成 `S01` 与 `S08`。

| ID | 任务 | 验收 |
| --- | --- | --- |
| S01 | 冻结当前源码/可执行/关键 case 为 baseline | 记录 archive/hash、build 类型、代表性 case 与现有性能；不改数值代码 |
| S02 | 建立开发总纲与方法权威关系 | 本文件组可从 `docs/README.md` 找到；旧 TODO 不再是自动任务源 |
| S03 | 建立架构规范 | owner/依赖/fallback 规则明确 |
| S04 | 建立数值硬约束 | Interface/Hydro/Electricity/Coupling invariants 明确 |
| S05 | 建立输出诊断规范 | accepted/rejected、series/meta、关键残差定义明确 |
| S06 | 建立验证矩阵与 task gate | T0–T4 和 I/H/E/C 顺序明确 |
| S07 | 历史文档降级 | 旧长 TODO 顶部标注 historical/reference，避免 Codex 继续从中自行扩展任务 |
| S08 | 冻结实验路径 | `parametricGals`、大型 CompatibleMoment、`sharpHermiteP3` 标成 baseline；不再加 V2 新功能 |

**Gate S：** S01–S08 完成后才进入 I01。

---

## I. Interface V2：纯 AMR Hermite–GALS（当前）；守恒原任务保留历史

| ID | 任务 | 唯一验收点 |
| --- | --- | --- |
| I01 | 定义 V2 state/API | first-jet skeleton、conserved cell phase volume、geometry snapshot ownership 清楚；无数值行为变化 |
| I02 | tracer 返回 `B_h + J_Bh` | Jacobian 与同一离散 characteristic map 的 numerical differentiation 一致 |
| I03 | vertex `phi + grad(phi)` first-jet transport | affine/analytic velocity test |
| I04 | edge-midpoint first-jet transport | edge value/normal derivative transport test |
| I05 | uniform `phi_xy` compact closure | target polynomial reproduction + convergence |
| I06 | AMR master/hanging skeleton closure | coarse/fine C1 trace 一致；hanging 不独立推进 |
| I07 | conserved cell phase volume state | analytic geometry initialization/volume integration |
| I08 | unique phase-face flux | internal face cancellation + bounded flux constraints |
| I09 | 2:1 reflux | coarse flux = sum(fine subflux) |
| I10 | interior bubble scalar closure | monotonicity、safeguarded root、topology/amplitude guard |
| I11 | authoritative Hermite geometry publication | volume/normal/curvature/cut geometry 同源；无第二份 production geometry |
| I12 | standalone benchmark | static、translation、rotation、vortex、AMR crossing、curvature convergence 全通过 |

**Gate I（当前）：** 保留 I01–I06/I11 的有效 first-jet、同源几何与 AMR 约束；I07–I10 的独立库存/守恒闭合由顶部授权设计变更替代。当前[I12p3r2](../task/I12p3r2_FirstJet_Gate_Final_Report.md)为 **PASS**，显式 `hermiteGals` / FirstJet 已成为规定速度独立 Interface 候选，允许按顺序开始 H01。该gate不代表Hydro/EHD或低黏耦合通过。

**I12 当前授权调整（2026-10-01）：** 按用户要求取消逐 cut cell 几何体积与独立库存的机器精度硬匹配，见 [设计变更](../task/I12_Geometry_Inventory_Accuracy_Design_Change.md) 与 [r16 TODO](../task/I12a11r16_Geometry_Inventory_Separation_TODO.md)。顺序为几何精度预算闭合 → 共享通量/库存有界性协调 → 原 I12 完整验收；停止以加密恢复机器精度闭合作为先决路线。守恒、严格界限、唯一几何及原准确性 gate 不变，I12 仍 FAIL。

---

## H. Hydro V2：native staggered CCH

| ID | 任务 | 唯一验收点 |
| --- | --- | --- |
| H01 | [native face/cell DOF layout](../task/H01_Native_Face_Cell_Dof_Layout_TODO.md) | **PASS**：[报告](../task/H01_Native_Face_Cell_Dof_Layout_Report.md)，23 场景 / 29 拒绝、7/7 CTest；Nu+Np 与 native mesh 规模线性增长，无全局 C1 高阶速度状态 |
| H02 | [唯一 integral divergence `D_I`](../task/H02b_Integral_Divergence_TODO.md) | **PASS**：[H02b 报告](../task/H02b_Integral_Divergence_Report.md)，75 场景/2,400 组字段/1,308 粗面求和/32 拒绝，最大局部缺陷1.88e-15；H02a 编译边界和 H02b 数学均通过 |
| H03 | [`C=-D_I^T` 自动生成](../task/H03_Negative_Transpose_Pressure_Load_TODO.md) | **PASS**：[报告](../task/H03_Negative_Transpose_Pressure_Load_Report.md)，75 场景/2,400 对 u,p，XY/RZ/2:1 最大伴随缺陷7.39e-16，32拒绝、25零轴模式，1/1新契约及8/8回归 |
| H04 | [single-phase saddle solve](../task/H04_Single_Phase_Saddle_TODO.md) | **PASS**：[报告](../task/H04_Single_Phase_Saddle_Report.md)，封闭域逐连通域真体积压力均值、39 algebra/21拒绝；三档 XY 连续 pressure/velocity MMS 末档阶2.026/2.003，1/1新契约及9/9回归；projection 不作物理算法 |
| H05 | [RZ true measure](../task/H05_Rz_True_Measure_Axis_Regularity_TODO.md) | **PASS**：[报告](../task/H05_Rz_True_Measure_Axis_Regularity_Report.md)，72测度场景/216解析场/9轴极限/216奇偶反射/16拒绝；RZ散度MMS末档阶2.004，1/1新契约与10/10回归；strain/hoop黏性归H06 |
| H06 | [energy viscosity `K_mu`](../task/H06_Energy_Viscosity_TODO.md) | **PASS**：[报告](../task/H06_Energy_Viscosity_Report.md)，实际native P2 strain、21场景/55解析场/6谱/105黏度缩放/25拒绝；XY/RZ uniform/2:1 strain阶2.059–2.356、耗散阶1.811–1.882；零黏度精确消失、1/1新契约和11/11回归 |
| H07 | 同源相分割/材料；拆为 H07a 与 H07b | **PASS**：[H07a](../task/H07a_Same_Geometry_Phase_Rate_Report.md)同源 cell/face 测度、瞬时率；[H07b](../task/H07b_Phase_Quadrature_Viscosity_Report.md)相内正权求积/H06黏性，22场景/75,632点/30拒绝，XY/RZ uniform/2:1耗散阶1.828–2.053；无中心判相或材料库存 |
| H08 | [唯一材料通量与非守恒缺陷](../task/H08_Material_Flux_Continuity_Defect_TODO.md) | **PASS**：[报告](../task/H08_Material_Flux_Continuity_Defect_Report.md)，三层几何质量/中点材料积分/实际非零缺陷，309场景/28拒绝；光滑时间阶1.997–2.000，1/1新契约+14/14回归；无库存或第二套 `rho_f u_f` |
| H09 | 同源 dual 材料与动量一致性；拆为H09a/H09b | **PASS**：[H09a真实dual几何/质量](../task/H09a_Momentum_Dual_Geometry_Report.md)与[H09b共享流量/非零缺陷/uR](../task/H09b_Dual_Material_Consistency_Report.md)均通过；H09b 273场景/24实际平移/20拒绝，常速度误差3.58e−12、光滑时间阶2.000；1/1新契约+17/17回归 |
| H10 | KEP/skew最终native空间对流 | **本次空间/KE gate PASS**：H10o闭合入口及v2完整84项通过；原一般trace与切界面范围见报告，下一H08a |
| H10c | [最终native联合矩](../task/H10c_Native_Convection_Moments_TODO.md) | **FAIL**：原affine dof911误差7.42e−14、1,516独立矩及19/19回归通过；完整84门槛中10空间+4原donor控制失败，uniform终值阶改善但整体空间不达标 |
| H10d | [整体空间一致性收口](../task/H10d_Final_Operator_Consistency_Closure_TODO.md) | **FAIL**：[报告](../task/H10d_Final_Operator_Consistency_Report.md)，15,344行定位，闭合全P2不相容；惩罚候选系数增长/最细准备失败并已退役。生产恢复H10c，原反例与4/4受影响检查通过；前置失败，未重跑完整研究 |
| H10e | [有界最终算子一致性](../task/H10e_Bounded_Operator_Consistency_Report.md) | **FAIL**：联合求解未收敛并退役，生产保持H10c；48组真实材料精确affine/skew缺口及加密趋势已记录 |
| H10f | [相容有界空间矩](../task/H10f_Compatible_Bounded_Moments_Report.md) | **FAIL**：均匀原有界交集可行；真实材料几何矩/PDE目标差异明确，率加权修正4组超界并退役，无新生产实现 |
| H10g | [点采样PDE材料矩](../task/H10g_Point_Pde_Material_Moments_Report.md) | **FAIL**：M*h^4目标两个材料空间阶失败并退役；48组必要下界/原材料来源核验通过，生产保持H10c |
| H10h | [最终算子联合可行性](../task/H10h_Final_Operator_Joint_Feasibility_Report.md) | **FAIL**：固定区间四组不可行；48 组/46,680 行独立证书通过，定位纯气体 dual 的输入 q/R 尺度；production 保持 |
| H10i | [同源空间材料通量设计](../task/H10i_Spatial_Material_Flux_Transfer_Report.md) | **FAIL / engineering PASS**：64常规布局来源核验通过；4个合法平面publication拒绝真实dual质量/重合零迹；0 production/0 test，先补source契约 |
| H10i1 | [重合相迹/clip契约](../task/H10i1_Coincident_Phase_Trace_Clip_Contract_Report.md) | **设计PASS / engineering PASS**：clip内部迹与material因果迹明确；16实际场景source仍FAIL，0 production/0 test |
| H10i1a | [同源native Q4 trace port](../task/H10i1a_Native_Quartic_Trace_Port_Report.md) | **PASS / engineering PASS**：4 production/1既有测试；1,408真实段/40零迹/80侧区间通过，mass/clip/publication未迁移 |
| H10i1b | [闭合裁剪矩](../task/H10i1b_Closed_Clip_Moments_Report.md) | **PASS / engineering PASS**：5 production/2既有测试；真实dual内部8组mass/16盒成功；16/16回归，native publication未迁移 |
| H10i1c | [正规native零边拓扑](../task/H10i1c_Native_Zero_Edge_Topology_Report.md) | **PASS / engineering PASS**：5 production/2既有测试；Core 16/16，native8组/40支撑，16/16回归，旧publication仍拒绝 |
| H10i1d | [只读半面trace cache](../task/H10i1d_Readonly_HalfFace_Trace_Cache_Report.md) | **PASS / engineering PASS**：2 production/1既有测试；16实际来源/1728半面，17/17回归，typed细化/普通视图/事务通过 |
| H10i1e | [FirstJet派生测量/实际发布](../task/H10i1e_FirstJet_Trace_Publication_Measure_Report.md) | **PASS / engineering PASS**：5 production/2既有测试；实际16组/352CV/native8组/40支撑，两方向账本、17/17回归通过；primary/Q4/transport/clock独立 |
| H10i1f | [native trace正式消费链](../task/H10i1f_Native_Trace_Phase_Consumers_Report.md) | **PASS / engineering PASS**：4 production/2既有tests；实际16组H07/真实dual mass，864 face/864 mass/32clip、17/17回归；没有新材料F |
| H10i2 | [独立来源/空间设计复核](../task/H10i2_Spatial_Material_Source_Reaudit_Report.md) | **设计/来源 PASS，engineering PASS**：0 production/0 tests，16实际source/1472overlap、39CSV/7回归；明确kinematics复用接口及实施预算 |
| H10j0 | [实际compatible advector](../task/H10j0_Spatial_Advector_Reconstruction_Report.md) | **PASS**：3文件/1既有test；64布局/128更新、共享/真实D_I/lift/复用通过，7/7回归；无phase/material F |
| H10j | [同源空间phase-rate owner](../task/H10j_Spatial_Phase_Rate_Owner_Report.md) | **PASS**：2 production/1 test；64布局/128factory、独立causal相率与21段gas witness、8/8回归；旧material consumer未迁移 |
| H10k | [唯一材料owner/source迁移](../task/H10k_Material_Spatial_Owner_Migration_Report.md) | **PASS**：2 production/2 tests；384准备/三密度、单份F/actual stamp、8/8回归；dual仍旧quarter |
| H10l | [同一F的dual迁移/quarter退役](../task/H10l_Same_F_Dual_Material_Migration_Report.md) | **FAIL / 已恢复H10k**：5 production/2 tests候选；3 basis/4 closed失败，576准备的独立身份核验通过 |
| H10l1 | [actual advector closure与不可穿透来源](../task/H10l1_Spatial_Advector_Closure_And_Boundary_Report.md) | **PASS**：2 production/2 tests；208 basis、32 boundary布局、64材料布局及8回归 |
| H10l2 | [same-F dual迁移重试](../task/H10l2_Same_F_Dual_Material_Retry_Report.md) | **PASS**：5 production/2 tests；退役quarter/legacy，576准备/208 basis/四closed/8回归 |
| H10m | [同一F下最终native对流复核](../task/H10m_Final_Native_Same_F_Convection_Report.md) | **FAIL**：0 production/0 test；原affine通过，空间问题仍在，独立donor四gate失败 |
| H10n | [donor耗散率验收定义](../task/H10n_Donor_Dissipation_Gate_Contract_Report.md) | **PASS**：0 production/0 test；12原窗口不变，独立D身份/四新gate通过，四旧FAIL单列 |
| H10o | [有界最终native空间对流](../task/H10o_Bounded_Final_Native_Convection_TODO.md) | **PASS**：2 production/1 test；原内部反例和完整v2通过，成本及边界范围见报告 |
| H08a | [同源时空相通量积分](../task/H08a_SpaceTime_Material_Flux_Report.md) | **PASS（规定范围）**：24 短检查、28 时空窗口及原 H10b 84 gate；唯一 F/真实 R/uR/KEP，小余项收敛；Pure FirstJet 独立；完整 App/全耦合二阶/效率未验收 |
| H08a1 | [同一实际速度分支jet](../task/H08a1_Actual_Advector_Cell_Jet_TODO.md) | **PASS**：2 production/1 test；144组合/69,804点及轴/真实regrid，原4,432行输出保持；局部导数不代表跨界特征线 |
| H08a2 | [实际分片速度特征线来源](../task/H08a2_Piecewise_Characteristic_Source_TODO.md) | **FAIL**：真实跳跃/288差分通过，24时间族中7个RZ Jacobian失败；候选退役，H08a1代码/8项回归保持 |
| H08a2r1 | [轨迹/J联合时间精度](../task/H08a2r1_Partition_Consistent_Characteristic_Report.md) | **PASS**：原24族、288差分、真实跳跃及9回归通过；原4,727行保持，尚未接FirstJet |
| H08a3 | [实际运动到FirstJet消费桥接](../task/H08a3_Resolved_FirstJet_Motion_Bridge_Report.md) | **FAIL（历史）**：原19,120 DOF审计保持；必要消费接入并入H08a，不再作为下一独立先决任务 |
| H08a2r2 | [原单侧选支精度先决](../task/H08a2r2_Skeleton_Motion_Selection_Report.md) | **FAIL**：32族/56,992实际DOF，uniform原单侧梯度一阶下界；AMR剪切前提部分失败。审计已退役，1379源码/输入恢复、4/4回归通过 |
| H08a2r3 | [共用正则advector可行性](../task/H08a2r3_Common_Regular_Advector_Report.md) | **PASS（未部署）**：144结构组、32族96窗口、sin六窗口；C1/面矩及原精度门槛通过，1379源码恢复、4/4回归保持 |
| H08a2r4 | [生产正则重构准备](../task/H08a2r4_Regular_Advector_Preparation_Report.md) | **PASS（准备层）**：144结构组/12+6 source布局、9/9回归；原4727+3364行保持，生产B未切换 |
| H08a2r5 | [材料/dual参考积分次数](../task/H08a2r5_Material_Reference_Polynomial_Degree_Report.md) | **PASS**：720新检查/672高精度复算、原1154行保持；5/5回归，production 0/test 2 |
| H08a2r6 | [benchmark多项式schema](../task/H08a2r6_Benchmark_Polynomial_Schema_Report.md) | **PASS**：34表/716228行保持（计时单列），18条三/五项回读、3/3架构目录；production 0/benchmark 1 |
| H08a2r7 | [完整法向多项式表示](../task/H08a2r7_Complete_Normal_Polynomial_Report.md) | **PASS**：1632高精度复算、9/9回归；原8813行/716228导出行保持，实际B仍旧；production 4/test 2 |
| H08a2r8 | [正则场共同解析内核](../task/H08a2r8_Regular_Advector_Evaluation_Report.md) | **PASS**：144结构组、C1/面矩/Gauss、完整atom及32拒绝；9/9回归、原9849行保持；production 2/test 1 |
| H08a2r9 | [空间来源消费者契约](../task/H08a2r9_Spatial_Source_Consumer_Contracts_Resume_Report.md) | **PASS（恢复）**：原64/144/208基覆盖保持；48布局/432输入两来源物理积分与完整次数导数、三负控及9/9回归通过；production 0/test 2 |
| H08a2r10 | [缺秩准备保留可辨识低阶族](../task/H08a2r10_Identifiable_Lower_Space_Preparation_Report.md) | **PASS**：六反例关闭，56/56原门槛、448独立输入、144结构组及9/9回归；满秩权重保持；production 2/test 1 |
| H08a2r11 | [首次同一正则来源部署](../task/H08a2r11_Common_Regular_Source_Deployment_Report.md) | **FAIL（历史）**：Q4系数预算/常数精确J及研究旧patch依赖未闭合，原候选退役；当前裁决见第二次恢复 |
| H08a2r12 | [正则来源浮点消费者契约](../task/H08a2r12_Regular_Source_Arithmetic_Consumer_Report.md) | **PASS**：48布局/1152输入、单位/Decimal复算、严格事件与数值J及负控；隔离候选phase/characteristic通过，production 0/tests 2 |
| H08a2r13 | [benchmark来源观察器](../task/H08a2r13_Benchmark_Source_Observer_Report.md) | **PASS**：旧patch导出退役，原37表/733364行保持，完整双端/五点参考、14负控及9/9回归；production 0 |
| H08a2r11恢复 | [同步部署审计](../task/H08a2r11_Common_Regular_Source_Deployment_Resume_Report.md) | **FAIL**：4原窄域affine/常数、4闭域零功消费者及非线性J参考未收口；6/9候选，原运动7组通过；候选退役 |
| H08a2r14 | [缺秩常数与可辨识子空间](../task/H08a2r14_Identifiable_Constant_Subspace_Report.md) | **PASS**：原4 dual反例关闭，64布局536输入/144结构/9回归；121842条权重逐位保持，仅216条缺秩改变，未部署 |
| H08a2r15 | [闭域真实材料动能参考](../task/H08a2r15_Closed_Material_Energy_Reference_Report.md) | **PASS**：原四closed关闭，124场景/6144 DOF与三负控、9/9回归；10773原行保持，production 0 |
| H08a2r16 | [非线性运动可信参考](../task/H08a2r16_Resolved_Nonlinear_Motion_Reference_Report.md) | **FAIL**：参考已解析，原两J超1e−11；24时间族/288FD/2304分区/54负控及9/9回归保持，production/test 0 |
| H08a2r17 | [跨分片point/J精度修复](../task/H08a2r17_Nonlinear_Crossing_Jacobian_Accuracy_Report.md) | **PASS（隔离算法）**：原两J反例关闭，原完整覆盖保持；该冻结算法已由r11第二次恢复部署 |
| H08a2r11第二次恢复 | [共同正则来源部署完成](../task/H08a2r11_Common_Regular_Source_Deployment_Resume2_Report.md) | **PASS**：4 production/2 test，同源场/运动、9回归、原内部affine/native及28窗口84gate全部通过；基础设施收口，立即回H08a时空相积分 |
| H11 | [局部受控耗散](../task/H11_Controlled_Dissipation_Report.md) | **PASS（规定范围）**：19 短检查、36 窗口/140 gate、原 H10b 84 gate；同一 F 上局部非负 D、实际材料 KE 和光滑二阶；非完整 Euler/NS |
| H12 | [同源常量法向应力](../task/H12_Constant_Normal_Stress_Report.md) | **PASS（静平衡范围）**：13检查，36几何/432块/864实际解及独立复算；无黏/高密度比、AMR及原残差门槛保持 |
| H13 | [同源Hermite毛细载荷](../task/H13_Weak_Traction_Acceptance_TODO.md) | **FAIL / 冻结A/B基线**：旧切向模态候选及未完成项保留；当前方法入口转BF-CSF，见本页方法变更注记。[原报告](../task/H13_Hermite_Capillary_Traction_Report.md)，不进入H14 |
| BF-00 | [方法切换与冻结基线](../task/BF-00_TODO.md) | **PASS（仅立法/冻结任务）**：用户明确解除旧档案缺失前置；当前真实源码/输入基线冻结，旧FAIL/缺失/timeout和gate保持。[报告](../task/BF-00_Report.md) |
| BF-01 | [Regularized color](../task/BF-01_TODO.md) | **PASS（仅本任务）**：[报告](../task/BF-01_Report.md)，123.714秒新契约及13/13回归、660行真实cell数据；sharp材料与旧H13保持 |
| BF-02 | [同一C的恒曲率CSF载荷](../task/BF-02_TODO.md) | **PASS（仅本任务）**：[报告](../task/BF-02_Report.md)，200次实际saddle、2/2新契约+10/10回归，完整free-row/wall reaction与来源审计；下一次BF-03 |
| H14 | midpoint predictor–corrector | variable `dt` 二阶时间；所有材料/力属于同一 midpoint candidate |
| H15 | AMR Hodge + production boundary set | AMR gradient MMS、free-slip/no-penetration/pressure traction/open backflow tests |

**Gate H：** H02/H03 不过不做 convection；H08/H09 不过不做高密度比；H12 不过不做真实 capillary production。

**当前 Hydro 状态（2026-10-07）：** [H13 同源Hermite毛细载荷](../task/H13_Hermite_Capillary_Traction_Report.md) **FAIL，仍为唯一当前任务**。原静滴七组伪流例外、无黏波能量/频率未收口；所有隔离载荷参考均未进入生产；Gram族及等密度1.557阶的简单加性候选均未采用；三个实际观察确认当前常权应力平均的信息损失；局部应力偏差及已知压力lifting均未采用；后者的等密度负刚度、0.716阶与能量升高已由原native/source分解复现；二相传输81组及九次FirstJet方向更新已完成；共享一圈应力矩候选末档频率阶-0.435、原能量峰值缺陷随加密增至0.452，独立矩阵指数/功核算复现，已拒绝；18组光滑传输及原实际算子复现确认native压力/加速度近二阶、最终source界面运动约一阶，固定二次应力profile预测后余量二阶以上；局部normal最小代价消差虽可达，但邻段扰动1.506倍，逐段叠加留下全部误差，已排除；完整P2压力矩差=tangent+jump/h已由原全矩阵核验，实际可达jump端口相对尺度约0.335且不随加密消失；仅改弱牵引虽频率二阶但能量缺陷约5.56；假设面通量矩匹配运动/牵引的线性参考频率1.987阶、原能量舍入闭合；具体共同source与P2弱牵引联合改造已获用户明确授权；真实source候选契约通过，完整弱牵引及全物理门槛待完成，候选未验收；一般曲面/AMR/实际动态及联合精度仍待解决。定位发现[直接阻断H08a的实际Jacobian反例](../task/H13_H08a_Transparent_Picard_Jacobian_Counterexample.md)，已按用户限定例外在一个现有文件中先代数消去C1内部边界相消项；25/25受影响检查通过，新版本完整H08a的28窗口及原空间/时间/R历史门槛通过；H13单周期复测仍失败。r11其他来源与原门槛保持，不另拆精修任务。历史H08a/H10/H11/H12规定范围结果保留；这一处限定正确性修复已收口。Pure FirstJet独立，无M_K/bubble、R=0或quarter commuting强制；唯一共享F、真实R/uR、KEP/skew保持。H14未启动。

**H07 拆分及用户授权（2026-10-03）：** 用户选择“完全接受非守恒材料输运”，见[设计](../task/H07_Pure_Gals_Hydro_Coupling_Design.md)。不恢复 Interface M，也不新增 Hydro 守恒材料库存；几何材料质量和非零连续性缺陷显式保留，动量/动能处理其一致性项。先执行 [H07a](../task/H07a_Same_Geometry_Phase_Rate_TODO.md) 的原生相测度/瞬时率，再执行 [H07b](../task/H07b_Phase_Quadrature_Viscosity_TODO.md)真实分相体积求积/H06 接入；两项均已验收，H07 PASS。H08/H09 按新路线验收，原质量机器精度闭合仅归历史比较；高密度比与无黏 gate 不提前豁免。

**H02 拆分（2026-10-03）：** 建立独立库/architecture 分类并新增散度测试会超过默认两份测试文件预算，因此先做 [H02a TODO](../task/H02a_Cch_Build_Boundary_TODO.md)，已 [PASS](../task/H02a_Cch_Build_Boundary_Report.md)：真实 CCH 库不依赖旧 Hydro/Compatibility，7/7 CTest、133 架构检查/31 拒绝，原 H01 23 场景/29 拒绝保持。其后 [H02b](../task/H02b_Integral_Divergence_Report.md)在此库实现唯一 D_I，75 场景的 random/解析积分、selected 边界和真实 2:1 局部闭合均通过；两项共同关闭 H02。随后 [H03](../task/H03_Negative_Transpose_Pressure_Load_Report.md)直接负转置生成 C 并通过伴随恒等式，[H04](../task/H04_Single_Phase_Saddle_Report.md)完成封闭域单体参考求解与均匀 XY 连续 MMS，[H05](../task/H05_Rz_True_Measure_Axis_Regularity_Report.md)验证同源真实物理散度及轴规则；[H06](../task/H06_Energy_Viscosity_Report.md)已从实际native strain组装能量型黏性并验证无黏极限；H07a/H07b、H08及H09a/H09b均已通过，当前下一项H10。AMR compatible Hodge、一般速度transfer与完整production边界仍归H15。

---

## IH. Interface/Hydro 无电联合验收

| ID | 任务 | 验收 |
| --- | --- | --- |
| IH01 | rigid/translation coupling | Interface phase flux 与 Hydro mass flux 同源 |
| IH02 | static drop | analytic curvature 与 Hermite curvature 两阶段分别通过 |
| IH03 | inviscid/low-viscosity oscillation | `mu -> 0` 连续，不靠算法切换 |
| IH04 | high-density-ratio moving interface | `rho_l/rho_g=1e3` 常规；`1e6` 压力测试 |
| IH05 | AMR-crossing interface | mass/pressure/capillary 无系统 coarse–fine source |
| IH06 | neutral jet/neck breakup | 真实 necking/topology 前后的材料漂移、动量/动能缺陷核算与稳定性达到 production gate |

**Gate IH：** 通过后才开始 Electricity V2 与最终 EHD 接线。

---

## E. Electricity V2：SI-HHO-Q

| ID | 任务 | 唯一验收点 |
| --- | --- | --- |
| E01 | `P2/P3` basis + quadrature | straight planar/RZ basis/mass exactness |
| E02 | weighted local potential reconstruction | random `P3` reproduction 到舍入误差 |
| E03 | stabilization | symmetric、nonnegative、constant mode exact |
| E04 | cell static condensation | local Schur symmetry；stable recovery；不显式 inverse |
| E05 | single-phase global skeleton solve | symmetry defect `<1e-12` + MMS/analytic convergence |
| E06 | Hermite phase element/subface geometry | curved volume/face quadrature convergence；先 `q_s=0` |
| E07 | shared interface trace + charge moments | `D_n+ - D_n- = q_s` 的 `P2` moments 满足 |
| E08 | equilibrated flux recovery | local/internal/global Gauss closure |
| E09 | direct panel current | current 由同一 flux residual 获得；停止调用旧 remap current path |
| E10 | BDF2–Schur | fixed geometry，`dt*sigma/epsilon = 0.1,1,10,100`，稳定性与时间阶 |
| E11 | 2:1 AMR + small-cut agglomeration | interface translation / small-cut sweep，condition 与误差可控 |
| E12 | Maxwell trace | `D_n`, `E_t`, `tau_n`, `tau_t` 独立收敛 |

**Gate E：** V1–V12 类 verification 全部完成后才接 Taylor cone。

---

## C. EHD / 电喷雾联合验收

| ID | 任务 | 验收 |
| --- | --- | --- |
| C01 | fixed geometry electrostatic coupling | geometry/material identity 与 electric result 对齐 |
| C02 | manufactured Maxwell traction | traction deposition 独立收敛 |
| C03 | constant electro-capillary jump | `Cp + F_Gamma` 达线性/舍入误差 |
| C04 | static conducting drop | `sigma*kappa + T_E,nn` 与 parasitic current |
| C05 | midpoint EHD corrector | geometry/charge/field/traction 同一 midpoint candidate |
| C06 | Taylor-cone steady state | 分别记录 `R_p,R_Gamma,R_mu,R_c`，不只看最终速度 |
| C07 | low-viscosity EHD sequence | `mu_l -> 0` 连续稳定 |
| C08 | electrospray jet/neck breakup | 相体积/质量/电荷/能量/性能完整记录 |
| C09 | AMR production run | 长时间、动态 AMR、restart/output reproducibility |

---

## P. 最后才做性能工程

只有对应模块已通过数值 gate 后，才允许做性能 task：

- prepared operator / factorization reuse；
- cache lifetime；
- local quadrature/vectorization；
- parallel scheduling；
- sparse solver/preconditioner；
- memory reduction。

性能 task 必须在相同 case、相同网格、相同 dt、相同 tolerance 下 A/B；不得靠降低精度换速度。

---

## 下一步建议

**当前 Hydro 状态（2026-10-07）：** [H13 同源Hermite毛细载荷](../task/H13_Hermite_Capillary_Traction_Report.md) **FAIL，仍为唯一当前任务**。原静滴七组伪流例外、无黏波能量/频率未收口；所有隔离载荷参考均未进入生产；Gram族及等密度1.557阶的简单加性候选均未采用；三个实际观察确认当前常权应力平均的信息损失；局部应力偏差及已知压力lifting均未采用；后者的等密度负刚度、0.716阶与能量升高已由原native/source分解复现；二相传输81组及九次FirstJet方向更新已完成；共享一圈应力矩候选末档频率阶-0.435、原能量峰值缺陷随加密增至0.452，独立矩阵指数/功核算复现，已拒绝；18组光滑传输及原实际算子复现确认native压力/加速度近二阶、最终source界面运动约一阶，固定二次应力profile预测后余量二阶以上；局部normal最小代价消差虽可达，但邻段扰动1.506倍，逐段叠加留下全部误差，已排除；完整P2压力矩差=tangent+jump/h已由原全矩阵核验，实际可达jump端口相对尺度约0.335且不随加密消失；仅改弱牵引虽频率二阶但能量缺陷约5.56；假设面通量矩匹配运动/牵引的线性参考频率1.987阶、原能量舍入闭合；具体共同source与P2弱牵引联合改造已获用户明确授权；真实source候选契约通过，完整弱牵引及全物理门槛待完成，候选未验收；一般曲面/AMR/实际动态及联合精度仍待解决。定位发现[直接阻断H08a的实际Jacobian反例](../task/H13_H08a_Transparent_Picard_Jacobian_Counterexample.md)，已按用户限定例外在一个现有文件中先代数消去C1内部边界相消项；25/25受影响检查通过，新版本完整H08a的28窗口及原空间/时间/R历史门槛通过；H13单周期复测仍失败。r11其他来源与原门槛保持，不另拆精修任务。历史H08a/H10/H11/H12规定范围结果保留；这一处限定正确性修复已收口。Pure FirstJet独立，无M_K/bubble、R=0或quarter commuting强制；唯一共享F、真实R/uR、KEP/skew保持。H14未启动。
