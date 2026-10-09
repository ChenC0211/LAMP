# LAMP 运行与调试

命令在仓库根目录执行。Windows 需 Visual Studio 2022 C++ 工具、CMake 和 Eigen/SuperLU，本机路径在 `CMakeUserPresets.json`。`.cmd` 包装器使用进程级 PowerShell 执行策略，无需改全局策略。

新电脑先按 [GitHub 与 VS Code 开发流程](Git_Workflow.md) 安装依赖，并将 [本机 preset 示例](CMakeUserPresets.example.json) 复制到根目录的 `CMakeUserPresets.json`，填写实际 D: 路径。该本机文件不提交 Git。

```powershell
.\tools\build.cmd
.\tools\run_case.cmd -List
.\tools\run_case.cmd hermite_gals_rotation_xy
.\tools\run_case.cmd hermite_gals_vortex_xy_dynamic_amr
```

程序为 `bin/LAMP.exe`。标准命令是 `cmake --preset local`、`cmake --build --preset local --target LAMP`；包装器会自动查找 VS 附带的 CMake。其他机器调整本机 preset 路径。

六类 GALS 基准：`hermite_gals_static_xy`、`hermite_gals_static_rz`、`hermite_gals_translation_xy`、`hermite_gals_translation_rz`、`hermite_gals_rotation_xy`、`hermite_gals_vortex_xy`。无后缀为固定均匀网格，`_amr` 为固定 AMR，`_dynamic` 为从均匀网格开始的动态切换，`_dynamic_amr` 为从 AMR 开始的动态切换。

每次目录为 `runs/<名称>/<YYYYMMDD_HHMMSS_mmm_标识>`，嵌套名称用双下划线。包含输入副本、`runtime` 程序/DLL 快照、`run.json` 状态与哈希、`check.log`、`run.log` 和 `output`。`latest.txt` 指向最新启动；以 `run.json` 判断是否完成。失败保留已接受结果和真实退出码。

```powershell
.\tools\run_case.cmd hermite_gals_translation_xy -CheckOnly
.\tools\run_case.cmd hermite_gals_translation_xy -PrepareOnly
.\tools\run_case.cmd hydro_static_drop/xy
```

预检检查配置且不推进时间；准备模式只复制输入/程序。默认 2 线程，可加 `-Threads 4`；BLAS 默认 1 线程。独立 Interface 的 dt/endTime 在 `system/run.dict`，网格/重网格在 `system/mesh.dict`；步数由 dt/endTime 推导，模块方法/边界在 `system/interface.dict`。修改后重新运行，不覆盖旧目录。

```powershell
.\tools\build.cmd -Test
.\tools\build.cmd -Preset local-benchmarks
.\tools\run_case.cmd benchmarks/hermite_gals_smoke
.\tools\run_case.cmd benchmarks/hermite_gals
.\tools\run_case.cmd benchmarks/hermite_gals_regrid
.\tools\run_case.cmd benchmarks/hermite_gals_space
```

完整研究配置在 `cases/benchmarks/<名称>/benchmark.json`，默认 CTest 不执行它们。`space`、`curvature`、`time`、`coarsening` 用于独立收敛研究；完整窗口验收与收敛通过分别记录。历史 Hydro 配置属于冻结/实验耦合，不表示新 GALS–Hydro 已验收。

```powershell
.\tools\build.cmd -Preset local-debug
.\tools\run_case.cmd hermite_gals_rotation_xy -PrepareOnly -Program bin/relwithdebinfo/LAMP.exe
```

VS Code 也可选择准备调试任务。复制终端打印的运行目录，F5 选择 LAMP 调试配置并粘贴。断点可设在 `InterfaceApp.cpp`、FirstJet step 或主入口。Release 与调试程序分开输出，使用相同数值路线。

ParaView 打开 `output/fields/solution.pvd`，界面在 `output/fields/interface`，误差历程在 `output/series`，来源与帧清单在 `output/meta`。短检查与临时调试结果在 `runs/checks`。
