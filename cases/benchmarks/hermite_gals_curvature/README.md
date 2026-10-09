# benchmarks/hermite_gals_curvature

完整 FirstJet 研究；驱动只调用生产 GALS 算子与独立参考。完整时间窗、收敛研究分别验收。

```powershell
.\tools\build.cmd -Preset local-benchmarks
.\tools\run_case.cmd benchmarks/hermite_gals_curvature
```

研究模式在 benchmark.json；记录与 CSV 行写入 runs/benchmarks__hermite_gals_curvature/<时间戳_标识>/run.log。程序和输入来源在 run.json；本驱动不输出 ParaView 场，查看场请运行对应 hermite_gals 配置。

统一操作见 [运行与调试](../../../docs/guides/Run_And_Debug.md)；配置与能力状态见 [算例清单](../../catalog.json) 。
