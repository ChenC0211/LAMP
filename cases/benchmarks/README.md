# 完整研究与 benchmark

这里保存完整时间窗、网格/时间收敛研究及专属驱动；小契约检查在 [tests](../../tests/README.md)。驱动只编排生产数值算子和独立参考，不向生产模块提供算法。

```powershell
.\tools\build.cmd -Preset local-benchmarks
.\tools\run_case.cmd benchmarks/hermite_gals
.\tools\run_case.cmd benchmarks/hermite_gals_regrid
.\tools\run_case.cmd benchmarks/hermite_gals_space
.\tools\run_case.cmd benchmarks/cch_convection
.\tools\run_case.cmd benchmarks/cch_constant_normal_stress
```

其他名称用 `run_case.cmd -List` 查询。研究配置在各目录的 `benchmark.json`；每次运行的来源、日志和状态仍归 `runs/<研究名称>/<时间戳>`。完整研究不属于默认 CTest。

[CCH光滑无黏对流](cch_convection/README.md)的规定闭合范围已通过H10o/H08a；[局部受控耗散](../../docs/task/H11_Controlled_Dissipation_Report.md)已通过H11研究，四项旧donor失败继续保留。完整CCH流体App仍在开发。

[H12常量法向应力](cch_constant_normal_stress/README.md)验证36个固定几何/864次真实鞍点静平衡，覆盖XY/RZ、AMR、密度比10^6和严格零黏，已按原门槛通过。解析常应力研究不代表真实曲率毛细流或时间推进已通过；复现命令只改研究名称。

[H13实际Hermite毛细载荷](cch_capillary/README.md)已有完整静滴扫描和真实动态候选，**当前FAIL**：伪流收敛例外及无黏波能量问题未闭合。入口为 `benchmarks/cch_capillary/static` 和 `benchmarks/cch_capillary`，退出0仅表示数据完成；完整H13未验收，不进入H14。
