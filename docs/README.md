# LAMP 开发文档

**当前 BF-02 PASS（2026-10-09）：** 同一C上的恒曲率实际CSF RHS已完成；[TODO](task/BF-02_TODO.md) / [报告](task/BF-02_Report.md) / [API](api/BF-02_Cch_Balanced_Csf_Load.md)。新两组契约86.7/9.6秒、200次实际saddle及10/10回归通过，完整9,528条cell/22,002条face数据归档。BF-01已封存PASS；下一次仅执行BF-03，本轮BF-03/H14/H15未启动。原H13 FAIL/缺失/timeout与完整动力学gate保持，以下旧入口为历史检查点。

**2026-10-09 当前任务：BF-00 / 方法切换确认与冻结基线。** 用户已将目标入口更新为 [BF-CSF任务书](roadmap/LAMP_BF_CSF_Codex_Tasks.md) 与 [BF-CSF V1方法](roadmap/LAMP_BF_CSF_Method_V1.md)。每次只执行一个BF-ID；当前 [BF-00 TODO](task/BF-00_TODO.md) / [报告](task/BF-00_Report.md)。本任务0数值变更，冻结当前native H12/H13源码/输入及现存证据。原long-wave/H12历史payload丢失，用户已确认；按BF-00原规则暂BLOCKED，恢复方案待明确决定。旧H13 FAIL、未采用的切向模态候选和旧报告保留；以下当前/最新标签为该旧路线的历史检查点，H14/H15不启动。

**当前授权081136（2026-10-09）：** 用户明确“授权”同一H13的[联合分相source/物理动能Hodge方案](task/H13_Phase_Cut_Kinetic_Hodge_Method_Change.md)，精确范围及执行格式在[现有TODO](task/H13_Weak_Traction_Acceptance_TODO.md)。开始授权范围内实施，H13仍FAIL；以下待授权/blocked说明属于历史，H14/H15未启动。

**最新045122（2026-10-09）：H13仍为唯一任务，FAIL。** [同一TODO](task/H13_Weak_Traction_Acceptance_TODO.md)、[独立动能证据](task/H13_Phase_Cut_Kinetic_Metric_Certificate.md)及[待授权联合方法](task/H13_Phase_Cut_Kinetic_Hodge_Method_Change.md)：72次原public saddle核验，准确normal仍有一阶反馈；半格原对角动能缺失导致条件Rayleigh下界。phase-aware唯一source与局部物理动能Hodge涉及冻结方法，尚未实施。生产/永久测试0修改，旧失败保留，H14未启动；下列“最新”段落为历史检查点。

**最新接续035647（2026-10-09）：H13仍FAIL。** [执行TODO](task/H13_Weak_Traction_Acceptance_TODO.md)及[直接反例](task/H13_Exact_Plane_Cut_Response_Counterexample.md)：189次原载荷分项解逐位复现原场；36组精确平面制造解、真实alpha/M及独立PDE核验确认一般cut的normal插值约一阶、高密度回复方向错误。n8最小反例约1秒复现，当前联合候选不能采用，H14未启动；须先推导一般cut的共同source/真实M/C/应力相容修正。下列003903数据与两个待批预算补充保留各自范围。

**最新接续（2026-10-09）：H13仍为唯一当前任务，FAIL。** [执行TODO](task/H13_Weak_Traction_Acceptance_TODO.md)与[003903报告](task/H13_Hermite_Capillary_Traction_Report.md#2026-10-09--003903原系数与完整静滴周期审计)：隔离KEP两行准备精修、独立新source reference同步均通过原完整契约，永久采用待两个明确预算例外。当前生产候选完整静滴72几何/1728解自然完成，独立原场/source载荷核算一致，六组寄生速度末档回升仍FAIL；新目录n8完整256步无黏周期完成，频率误差2.25%、原峰值能量缺陷1.2151e-4。三档频率/时间/AMR/低黏及材料长回归未验收，H14不启动。此前86步停止及旧检查点保留历史结论，当前API以最新说明为准。

**当前 Hydro 状态（2026-10-07）：** [H13 同源Hermite毛细载荷](task/H13_Hermite_Capillary_Traction_Report.md) **FAIL，仍为唯一当前任务**。原静滴七组伪流例外、无黏波能量/频率未收口；所有隔离载荷参考均未进入生产；Gram族及等密度1.557阶的简单加性候选均未采用；三个实际观察确认当前常权应力平均的信息损失；局部应力偏差及已知压力lifting均未采用；后者的等密度负刚度、0.716阶与能量升高已由原native/source分解复现；二相传输81组及九次FirstJet方向更新已完成；共享一圈应力矩候选末档频率阶-0.435、原能量峰值缺陷随加密增至0.452，独立矩阵指数/功核算复现，已拒绝；18组光滑传输及原实际算子复现确认native压力/加速度近二阶、最终source界面运动约一阶，固定二次应力profile预测后余量二阶以上；局部normal最小代价消差虽可达，但邻段扰动1.506倍，逐段叠加留下全部误差，已排除；完整P2压力矩差=tangent+jump/h已由原全矩阵核验，实际可达jump端口相对尺度约0.335且不随加密消失；仅改弱牵引虽频率二阶但能量缺陷约5.56；假设面通量矩匹配运动/牵引的线性参考频率1.987阶、原能量舍入闭合；具体共同source与P2弱牵引联合改造已获用户明确授权；真实source候选契约通过，完整弱牵引及全物理门槛待完成，候选未验收；一般曲面/AMR/实际动态及联合精度仍待解决。定位发现[直接阻断H08a的实际Jacobian反例](task/H13_H08a_Transparent_Picard_Jacobian_Counterexample.md)，已按用户限定例外在一个现有文件中先代数消去C1内部边界相消项；25/25受影响检查通过，新版本完整H08a的28窗口及原空间/时间/R历史门槛通过；H13单周期复测仍失败。r11其他来源与原门槛保持，不另拆精修任务。历史H08a/H10/H11/H12规定范围结果保留；这一处限定正确性修复已收口。Pure FirstJet独立，无M_K/bubble、R=0或quarter commuting强制；唯一共享F、真实R/uR、KEP/skew保持。H14未启动。

根目录仅保留现行关键约束，阅读顺序：

1. [00 开发总纲](00_LAMP_Development_Charter.md)
2. [01 数值硬约束](01_LAMP_Numerical_Invariants.md)
3. [02 架构](02_LAMP_Architecture_Rules.md)
4. [03 字典](03_LAMP_Dictionary_Rules.md)
5. [04 诊断输出](04_LAMP_Diagnostics_Output_Rules.md)
6. [05 目录与文档](05_LAMP_Repository_Rules.md)
7. [06 代码与测试](06_LAMP_Code_And_Test_Rules.md)
8. [07 使用与调试](07_LAMP_Run_And_Debug_Rules.md)

| 分区 | 内容与入口 |
| --- | --- |
| [roadmap](roadmap/05_LAMP_Module_Task_Roadmap.md) | 当前开发顺序；历史 TODO 不自动产生新任务 |
| [guides](guides/Run_And_Debug.md) | 运行、构建、调试和结果查看 |
| `api` | 接口设计；当前 FirstJet 从 [五文件字典 schema](api/G04_Interface_Dictionary_Schema.md)进入；Hydro V2 从 [H01 原生布局](api/H01_Cch_Native_Dof_Layout.md)、[H02a 编译入口](api/H02a_Cch_Build_Boundary.md)、[H02b 积分散度](api/H02b_Cch_Integral_Divergence.md)、[H03 压力载荷](api/H03_Cch_Pressure_Load.md)、[H04 封闭域鞍点参考求解](api/H04_Cch_Single_Phase_Saddle.md)与 [H05 真测度/轴规则](api/H05_Cch_Physical_Divergence.md)及[H06 能量黏性](api/H06_Cch_Energy_Viscosity.md)进入 |
| `task` | 任务书、报告、设计变更；[H06 PASS](task/H06_Energy_Viscosity_Report.md)、[H05 PASS](task/H05_Rz_True_Measure_Axis_Regularity_Report.md)、[H04 PASS](task/H04_Single_Phase_Saddle_Report.md)、[H03 PASS](task/H03_Negative_Transpose_Pressure_Load_Report.md)、[H02b / H02 PASS](task/H02b_Integral_Divergence_Report.md)、[H02a PASS](task/H02a_Cch_Build_Boundary_Report.md)、[H01 PASS](task/H01_Native_Face_Cell_Dof_Layout_Report.md)、[长涡旋基线](task/Interface_Long_Vortex_Verification_Report.md)、[G04 字典职责固化](task/G04_Interface_Dictionary_Ownership_Report.md)、[G03 通用输入报告](task/G03_Interface_Generic_Input_Report.md)、[G02 整理任务](task/G02_LAMP_Engineering_Conventions_TODO.md) |
| `methods` | 三份冻结数学方法 PDF；保持原文件字节和原名称 |
| `assets` | 插图 |
| `archive`、`history_todo`、`manuals` | 旧文档、备份与历史手册；不覆盖现行约束 |

[Gate I 报告](task/I12p3r2_FirstJet_Gate_Final_Report.md)记录当前界面的实际验收范围。程序更名及文件退役不修改历史报告的数字、源码哈希或结论。
