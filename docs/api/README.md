# 接口索引与状态

**当前 Hydro 状态（2026-10-07）：** [H13 同源Hermite毛细载荷](../task/H13_Hermite_Capillary_Traction_Report.md) **FAIL，仍为唯一当前任务**。原静滴七组伪流例外、无黏波能量/频率未收口；所有隔离载荷参考均未进入生产；Gram族及等密度1.557阶的简单加性候选均未采用；三个实际观察确认当前常权应力平均的信息损失；局部应力偏差及已知压力lifting均未采用；后者的等密度负刚度、0.716阶与能量升高已由原native/source分解复现；下一步只诊断局部二相压力传输的最终native/FirstJet耦合；一般曲面/AMR/实际动态及联合精度仍待解决。定位发现[直接阻断H08a的实际Jacobian反例](../task/H13_H08a_Transparent_Picard_Jacobian_Counterexample.md)，已按用户限定例外在一个现有文件中先代数消去C1内部边界相消项；25/25受影响检查通过，新版本完整H08a的28窗口及原空间/时间/R历史门槛通过；H13单周期复测仍失败。r11其他来源与原门槛保持，不另拆精修任务。历史H08a/H10/H11/H12规定范围结果保留；这一处限定正确性修复已收口。Pure FirstJet独立，无M_K/bubble、R=0或quarter commuting强制；唯一共享F、真实R/uR、KEP/skew保持。H14未启动。

当前 Interface 入口为纯 Hermite–GALS / FirstJet。阅读从 [五文件生产 schema](G04_Interface_Dictionary_Schema.md)、[纯路线](I12p1_Pure_Hermite_Gals_API.md)、[时间步](I12p3a_FirstJet_Step_API.md)、[AMR 事务](I12p2b_FirstJet_AMR_Transfer_API.md)、[输出](I12p3c_FirstJet_Output_API.md)开始。旧 p3b 命名 benchmark 和 G03 单文件布局已退役；基础 first-jet、trace、mixed/master closure 对应 I02–I06。

I07–I12a 系列中的相库存、reflux 和守恒闭合描述历史 GALS 研究接口；当前 App 不启用这些输运路线。旧 InterfaceState/参数曲线耦合只在 Compatibility 中维持已有 Hydro/EHD 消费者。API 文件的存在不表示相关方法已成为生产入口，也不自动产生开发任务。

Hydro V2 从 [H01 原生 DOF 布局](H01_Cch_Native_Dof_Layout.md)、[H02a 编译入口](H02a_Cch_Build_Boundary.md)、[H02b 唯一积分散度](H02b_Cch_Integral_Divergence.md)、[H03 负转置压力载荷](H03_Cch_Pressure_Load.md)、[H04 封闭域鞍点参考求解](H04_Cch_Single_Phase_Saddle.md)和 [H05 真测度/轴规则](H05_Cch_Physical_Divergence.md)及[H06 能量型黏性](H06_Cch_Energy_Viscosity.md)进入。真实 CCH 库不链接旧 Hydro/Compatibility；同源散度/压力功伴随及封闭域 reference solve 已 PASS，均匀 XY 连续 MMS 验证 pressure/velocity 二阶空间精度。H07b真实分相求积/黏性、H08中点材料流量/非零缺陷与H09实际dual材料/常速度一致性已验收。H08a真正时空材料积分、H10最终KEP与H11受控耗散已在规定范围通过；[H12常量法向载荷](H12_Cch_Constant_Normal_Stress.md)完成实际静平衡。[H13毛细候选API](H13_Cch_Hermite_Capillary_Traction.md)已实现但完整验收FAIL；尚无完整CCH流体App或全耦合时间推进。

现行约束在 docs 根目录，下一项开发以 [路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md)为准，实际通过范围见 [Gate I](../task/I12p3r2_FirstJet_Gate_Final_Report.md)。
