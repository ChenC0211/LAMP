# LAMP 使用、构建与调试约束

状态：现行工程约束，2026-10-02。具体操作见 [运行与调试指南](guides/Run_And_Debug.md)。

用户入口固定为：

```powershell
.\tools\build.cmd
.\tools\run_case.cmd -List
.\tools\run_case.cmd hermite_gals_vortex_xy_dynamic_amr
```

更换算例只改末尾名称。算例配置只读，修改模板后在新的时间戳目录运行；结果固定在 `runs/<算例>/<时间戳_标识>`。启动器保存输入/可执行程序与哈希、线程设置、预检/物理运行退出码和耗时。失败保留日志与已接受前缀，不声称完成。

独立 Interface 的新工况使用 [五文件 schema](api/G04_Interface_Dictionary_Schema.md)，按职责修改 run/geometry/mesh/fields/interface 字典；禁止为运行算例往生产入口添加问题名称分支。解析误差 gate 在独立 benchmark 驱动中执行，生产输出 COMPLETED 只表示配置的步骤完成。

标准构建使用 CMake presets。`build.cmd -Test` 构建并执行默认短检查；`local-debug` 构建带符号程序到 `bin/relwithdebinfo`；`local-benchmarks` 按需构建 `cases/benchmarks` 驱动。机器依赖路径只放 `CMakeUserPresets.json`。

调试先用启动器 `-PrepareOnly -Program bin/relwithdebinfo/LAMP.exe`，再用 F5 加载该运行目录。禁止通过编辑生产代码写死 case、关闭接受条件、提高黏度或跳过断言来“调通”。临时计时/输出只能由被动观察器产生，不能进入数值分支；正式保留时归 Diagnostics/IO，并有明确字段与消费方式。

新脚本必须有用途、输入/输出、失败退出码和索引。脚本只负责构建、复制、运行或后处理；不得藏入另一套求解器。专属研究参数放对应 `cases/benchmarks/<研究>`；不得再为每个算例复制一份 `run_*.ps1`。

ParaView 打开运行目录内 `output/fields/solution.pvd`。后处理只读取已保存的结果及其来源，不重新构造一套“更好看”的界面作为解。
