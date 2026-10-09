# H08a2r6 / benchmark 法向多项式导出 schema 报告

2026-10-05，**PASS（导出格式）**。三处 benchmark CSV 的系数列名和值统一按实际容器生成，当前三项格式保持，五项高次和尾部字段已独立回读验证。production 0、benchmark 1；生产 B/atom/J/F 未改。

[任务书](H08a2r6_Benchmark_Polynomial_Schema_TODO.md)、[实施前设计](H08a2r6_Benchmark_Polynomial_Schema_Design.md)。本项仅关闭表示升级前的输出列错位风险，不代表 H08a 时空输运完成；H11 未启动。

## 修改与验证

`cases/benchmarks/cch_convection/driver.cpp` 增加共同 header/data writer。header 直接使用 `CchNormalPolynomial` 的 begin/end coefficient 类型长度；零高次项也完整写出，无固定三项假设。三处是 `spatial_zero_rate.csv`、`advector_atoms.csv` 和 `boundary_atoms.csv`。现有模式、配置、gate 及物理算法未改。

标准 Release 分别构建并捕获修改前和最终程序，完整执行原三个模式。独立 CSV parser 检查所有 header 唯一、行宽匹配、文件集一致及全部物理/来源字段：

| 模式 | CSV 数 | 数据行 | 比较结果 |
| --- | ---: | ---: | --- |
| trace-contract | 25 | 20,460 | 全部文件逐字相同 |
| advector | 4 | 342,244 | 除原 prepare/evaluate 计时两列，全部字段相同 |
| boundary | 5 | 353,524 | 全部文件逐字相同 |

共 **34 份 CSV / 716,228 行**。advector 覆盖原 64 布局/128 revisions，boundary 覆盖原 32 布局/64 revisions；没有缩短模式。两计时列单列为有限非负耗时，不当作物理变化。

另以实际共同 writer 生成三/五项、单端/双端及三种尾部字段共 **6 份合成表 / 18 行**。独立 parser 从真实旧 header 推导预期新列名，按名称逐项核对非零 c3/c4/e3/e4 和所有元数据/尾部值。审计源码与最终源码只相差归档中的合成入口，helper 完全相同。临时入口已移除，最终程序重新标准构建并用于上述全部复跑；没有新增永久审计命令。

架构、架构自检、目录检查 **3/3 通过**。1,381 项源码/输入清单只有预算内 driver 和配套 README 变化；生产、测试、工具与构建文件未改。六份现行 phase/momentum/material/dual/KEP/characteristic 短契约 exe 与 H08a2r5 哈希相同，本项不重复执行其不受影响的数值全场景。

## 剩余事项与证据

下一唯一任务 [H08a2r7](H08a2r7_Complete_Normal_Polynomial_TODO.md)升级完整 `CchNormalPolynomial`/phase rate 表示及 phase/momentum 参考。实际 B/atom/J 与 characteristic 证书需后续同步部署；来源切换后仍必须重跑原 native、内部 affine、H10b-v2 全部 28 窗口/84 项。本格式任务未运行完整 H10b，也不改变历史 FAIL。

四材料硬约束保持；不增加 R=0、逐 cell 库存或 bubble 要求。当前材料仍为空间相率配合时间步，尚未完成同源时空相积分。

检查：`runs/checks/H08a2r6/20261005_143319_179`；冻结：`runs/baselines/H08a2r6_20261005`。`schema_final.json` 保存逐表统计；`verify_schema.py` 可复算独立列/值和全部 before/after 比较；三个 runtime 目录绑定修改前、可退役审计、最终源码与 exe；`verification.json` 绑定源码、构建和未改六份 exe。临时脚本仅留证据归档。
