# H13 Hermite 曲率毛细研究（当前候选FAIL）

最新035647：[精确平面MMS直接反例](../../../docs/task/H13_Exact_Plane_Cut_Response_Counterexample.md)确认一般cut的normal插值约一阶和高密度法向回复方向错误；当前联合候选不能采用。n8最小反例约1秒复现，程序/fields归该任务证据目录，未改本入口。原静滴/波/两项待批修正结果保留，H14未启动。

2026-10-09接续：当前入口已绑定同一共同source的完整局部弱牵引候选，API见[现行说明](../../../docs/api/H13_Cch_Hermite_Capillary_Traction.md)。完整静滴72/1728自然完成，独立原场/72几何source载荷核算一致，但六组末档寄生速度回升，仍FAIL。n8新目录无黏256步完整周期自然完成，频率误差2.25%、原能量峰值缺陷1.2151e-4；三档/AMR/低黏尚未验收。隔离KEP准备修正与新source reference同步通过原完整契约，永久采用待两个独立预算例外。详见[003903报告](../../../docs/task/H13_Hermite_Capillary_Traction_Report.md#2026-10-09--003903原系数与完整静滴周期审计)，不能声称完整毛细能力。此前86步停止及以下日期条目保留各自旧范围，不迁为新候选验收。

```powershell
.\tools\build.cmd -Preset local-benchmarks -Target LAMP_cch_capillary_benchmark
.\tools\run_case.cmd benchmarks/cch_capillary
.\tools\run_case.cmd benchmarks/cch_capillary/static
```

当前主配置运行实际毛细波候选，XY uniform/AMR、n8/16/32，各推进一个解析参考周期；`static` 运行72个二次/非多项式圆球几何与1728次解析/Hermite载荷求解。退出0仅表示数据完成，尚不代表[H13完整任务](../../../docs/task/H13_Hermite_Capillary_Traction_TODO.md)通过；完整扫描、空间/时间与低黏序列仍须独立验收。

2026-10-06结果见[H13报告](../../../docs/task/H13_Hermite_Capillary_Traction_Report.md)：静态全部数据/装配审计通过，但7组伪流末档回升；实际无黏波完整n8周期的峰值相对扰动能量缺陷45.2792，短窗n8/16/32及减半dt尚未建立能量收敛。剩余原长扫描已按失败候选主动停止，保留程序、完成窗和AMR前缀；不得把此入口当作已验证毛细生产能力。

研究用毛细子系统的对称 kick–drift–kick：半步真实native鞍点力/黏性CN，整步由中间实际速度驱动已有H08a同源时空材料积分及FirstJet输运，再作终点半步力/黏性更新。rho_l=rho_g=1，保留实际材料缺陷和R的动能量级作为小振幅线性化省略项诊断；不把此研究积分器冒称H14完整NS中点corrector。

侧壁采用显式phi法向对称数据，边界phi及切向一阶量从原始几何沿实际运动回溯，仅法向导数按声明边界取零；没有解析未来波幅、频率或形状重置。法向速度使用既有impermeable接口。参考频率只用于周期长度与事后比较，不进入更新。

原始历史包含实际几何模态、表面长度、动能、毛细/黏性功及真实材料缺陷。低黏、实际曲率波动、边界误差与能量收敛都需据数据判断，不能用瞬时回复方向代替完整动态验证。

当前版本另保存每步surface、native faces/cells、压力、FirstJet primary/closure与瞬时功交换，均为被动观察。研究CLI的 `probe [n [steps [periods]]]` 仅供独立runs/checks目录中的失败定位，不改变默认完整研究。只有可复现且直接阻断H08a的数学反例允许限定修复；本轮已确认一处透明内部片Jacobian相消错误并作等价代数修正，见[反例说明](../../../docs/task/H13_H08a_Transparent_Picard_Jacobian_Counterexample.md)。其他来源保持，不以毛细失败自动产生基础设施任务。

`work [n]` 为只读载荷范围研究：用原公开运动查询逐基向量导出B、物理求积及严格native无散循环，不更新生产状态。2026-10-06后续隔离功配对原型仅位于 `runs/checks/H13/20261006_070108`，虽改善能量，频率未达原目标，未采用；标准wave/static入口仍使用原候选，完整结论见H13报告。
