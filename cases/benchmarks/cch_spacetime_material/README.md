# 同源时空材料通量研究

H08a 验证驱动。2026-10-06 的 28 窗口研究及配套契约通过规定的积分、一致性与小余项收敛范围，见[报告](../../../docs/task/H08a_SpaceTime_Material_Flux_Report.md)。完整流体 App、全耦合二阶及效率未验收。调用生产 Pure FirstJet/CCH 时空积分，以平移圆/轴上球的解析 swept-box 体积独立核对通量。初始标量使用非多项式 `exp(r²−a²)−1`；参考只用于输出误差及给定的物理外边界，不改内部几何、F 或 R。

```powershell
.\tools\build.cmd -Preset local-benchmarks -Target LAMP_cch_spacetime_material_benchmark
.\tools\run_case.cmd benchmarks/cch_spacetime_material
```

固定物理窗口 T=0.04，XY/RZ、uniform/2:1，空间 n=4/8/16，时间 n=8 的 1/2/4/8 步。当前研究 owner 的密度比为 1e6；三密度共享恒等式另由短契约检查。时间连续窗口消费上一实际 end publication，不能用解析状态重置。记录通量、几何误差、primal/dual 积分缺陷及缺陷率；不以积分缺陷随 dt 变小替代率的收敛。

`full` 是正式全数据；`space`/`time` 选择一个研究矩阵；`probe` 是开发者 n4 两布局定位，应在新的 `runs/checks` 中运行。成功退出仅表示数据生成完成，不表示 H08a 通过；验收须另核对预设收敛门槛和原 H10b 全窗口。
