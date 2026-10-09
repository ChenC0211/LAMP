# H12 解析常量法向应力研究

标准构建与运行：

```powershell
.\tools\build.cmd -Preset local-benchmarks -Target LAMP_cch_constant_normal_stress_benchmark
.\tools\run_case.cmd benchmarks/cch_constant_normal_stress
```

结果自动进入 `runs/benchmarks__cch_constant_normal_stress/<timestamp>/`。完整输入矩阵及验收见 [H12 TODO](../../../docs/task/H12_Constant_Normal_Stress_TODO.md)。固定几何的真实 `A=M/0.125+K_mu` 鞍点解，共 36 几何、432 个矩阵、864 个正/负常应力求解；每个参数组合独立裁决。原始网格、相体积、质量、黏性矩阵、压力、速度、载荷、壁面反力及残差写入 `output/series`。

载荷实际调用 `CchInterfaceNormalLoad`；使用同一 Pure FirstJet alpha 与压力算子 C，RHS 为负载荷。参考压力仅供事后核对，不参与求解。密度比覆盖 1 至 10^6、黏度包含严格零；RZ 包含轴和真实粗细子面。dt 仅定义静态质量块，不表示时间推进。

本研究不验证曲率驱动的毛细流、动态界面、真实电应力或完整 Euler/NS；一般时变法向力、切向力与全耦合分别由后续任务验收。`probe` 仅供 runs/checks 中短定位，不替代完整 `full`。
