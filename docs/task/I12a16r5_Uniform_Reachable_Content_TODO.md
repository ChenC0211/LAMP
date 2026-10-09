# I12a16r5 / 可达来源具有同一含液量时的代数积分

目标：只解决完整可达来源上含液量严格为常量时，混合内部／边界来源的相通量可直接求值这一事实。先读 [设计](I12_Uniform_Reachable_Content_Design.md)。这是对r4下一项范围的明确修正：不要求先求出对常量乘积无影响的完整B/J，不声称证明任意轨迹终止。

允许修改：PolynomialVelocitySlab.h/.cpp、HermiteGalsTemporalPhaseFlux.h/.cpp（4production）；bounded_spacetime_flux_contract.cpp及至多一个现有Core contract（<=2tests）。禁止修改actual tracer、速度evaluate、库存装配/reflux、I10/I11、benchmark案例/阈值、冻结和其他模块。

数学约束：独立M及几何精度预算保持分离。对同一owned polynomial source、原物理domain及old→arrival历史，严格包围有限速与原舍入操作可能到达的来源；旧native Q4闭覆盖采用同一snapshot的wl/wg，所有可能边界采用同一显式history。只有联合含液量上下界完全相等才提取常因子L±=cT±。不按epsilon合并不同权重，不从几何gas推断M=0，不凭点采样判相。

实现要求：独立Core可达来源对象不含B/J，不声称事件唯一或原tracer一定终止；Interface在原完整回溯认证前使用显式常量规则。保留原immediate-boundary分支；新分支必须验证全域速度/CFL、FP、时间、闭leaf覆盖、source identity及完整边界历史。RZ新能力先限严格零径向速度。非恒定来源继续原算法，不能按小E丢弃pending。

验收：T0原39targets、26CTest、81freeze/3PDF/r4archive保持，无新增warning；未新增能力的原输出保持。T1 XY/RZ轴向、uniform/AMR、gas/liquid、混合命中/内部、时间反向与角点常量，独立积分参考；实际成功轨迹必须落在可达包围中；库存权重差异、隐藏界面、缺失历史、CFL/FP/axis/opaque/source等拒绝或回退。T2原四个driver smoke与full实际运行，保留原30/400步、预算与全部native faces；真实记录新失败阶段或完整通过，不将局部认证冒充I12通过。

失败规则：若新包围遗漏actual来源、非常量被合并、旧结果改变或guard被绕过则FAIL并修复。算法尚未完成不标外部BLOCKED；只有真实外部依赖阻断才报告。I12完整gate未通过前保持FAIL，不进入Hydro。

交付：设计、代码、全部运行/比较/审计/归档、report及未解决问题。

执行结果：**PASS（本项）**。24积分场景、800成功actual traces、64原map拒绝、30极值来源点、13拒绝检查通过；旧integral完整前缀保持。39targets及26CTest全部通过，24项旧完整输出不变；81freeze/3PDF/r4archive保持。四个原smoke第一步均通过，原full在第二步库存装配真实越界并拒绝，I12仍FAIL。见 [报告](I12a16r5_Uniform_Reachable_Content_Report.md)。
