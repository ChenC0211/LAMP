# 通用工具

| 入口 | 用途 |
| --- | --- |
| `build.cmd` / `build.ps1` | CMake presets 配置、构建；`-Test` 执行短检查 |
| `run_case.cmd` / `run_case.ps1` | 按算例名称运行、预检、准备调试；`-List` 列举配置和研究 |
| `postprocess/analyze_hydro_standard.ps1` | 已保存 Hydro 基准的只读分析；`-Run` 指向 runs 时间戳目录 |
| `postprocess/analyze_interface_reference.ps1` | 已保存 Interface 的独立解析参考检查与图；验证输入/程序来源，读取 reference/verification.json |
| `docs/` | 历史方法/手册构建工具；用途见子目录说明 |

普通使用无需浏览研究历史脚本。已退役的专属运行、绘图与分析脚本保存在 G02 源码快照，不再是入口。新增工具按 [规范](../docs/07_LAMP_Run_And_Debug_Rules.md)归类并登记。

运行器实时显示程序的 stdout/stderr，并同时刷新 UTF-8 `run.log`；预检输出同样实时保存到 `check.log`。程序负责按 `Log.level` 控制 summary/detail/debug，运行器不延迟或吞掉输出。成功/失败退出码、预检与准备状态仍记录在 `run.json`。

含 `retired.json` 的历史算例不进入 `-List`；直接选择时返回退役理由，不启动程序。原输入和历史运行结果保留。
