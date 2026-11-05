# LAMP

- [App 工作流组织规范](docs/api/App_Workflow_Contract.md)
**Levelsetbased Adaptive Multiphysics Platform**：基于水平集的自适应多物理平台。

```powershell
.\tools\build.cmd
.\tools\run_case.cmd -List
.\tools\run_case.cmd hydro_static_drop
```

构建使用同一缓存目录增量编译，复用未受源码或头文件改动影响的目标文件；默认只构建主程序，不自动跑测试。固定使用同一 preset、Release 配置和工具链，日常无需 clean 或 Rebuild。修改公共头文件后，其依赖文件仍需重新编译；切换编译参数或文件并行模式也可能触发一次重编译。

已有成功配置的构建目录时，日常可用 `.\tools\build.cmd -SkipConfigure` 跳过显式配置，直接增量编译；首次构建或更换工具链/配置选项时仍用不带该参数的命令。

换算例只改末尾名称。配置在 `cases`；每次结果在 `runs/<算例>/<时间戳_标识>`；ParaView 打开该目录的 `output/fields/solution.pvd`。

已保存检查点的算例，中断后可交互续跑：

```powershell
.\tools\run_case.cmd vof_a_liquid_gas_drop_deformation_amr -Resume
```

终端会列出该算例的旧运行，输入序号或粘贴旧运行/检查点路径即可；续跑结果写入新目录。完整参数见[检查点与中断续跑说明](docs/guides/Run_And_Debug.md#vof-检查点与中断续跑)。

在另一台 Windows 10（1903 起）/11 x64 电脑运行，使用 `dist` 下的便携 ZIP，完整解压到可写目录后双击 `Start.cmd`，选择算例和“新运行 / 续跑 / 检查字典”。包内包含程序、所需 DLL 和算例输入，无需安装 VS、CMake、Python 或数值库。保留整个目录结构；配置在 `cases`，结果在 `runs`。制作与使用见[便携发布](docs/guides/Run_And_Debug.md#便携发布windows-x64)。

编译好程序后，运行以下命令，按提示选择要携带的算例；脚本直接打包现有 `bin/LAMP.exe`，预检后生成新的便携 ZIP：

```powershell
.\tools\release_portable.cmd
# 也可直接指定一个或多个算例：
.\tools\release_portable.cmd -Cases "vof_a_liquid_gas_drop_deformation_amr,vof_b1_charged_drop"
```

只收录选中的算例及其输入资源、运行程序、必要 DLL 和许可文件，不带构建目录、历史结果或未选算例。运行库缓存位于 `.portable/runtime`，后续自动复用。打包脚本不调用编译器或 CMake，无需添加跳过编译的参数。

- [使用与调试](docs/guides/Run_And_Debug.md)
- [GitHub 与 VS Code 开发流程](docs/guides/Git_Workflow.md)
- [算例索引与能力状态](cases/README.md)
- [短检查](tests/README.md) · [完整 benchmark](cases/benchmarks/README.md)
- [开发约束](docs/README.md) · [当前唯一主线](docs/roadmap/LAMP_Two_Phase_Solver_Mainline.md)
- [G02 整理任务](docs/task/G02_LAMP_Engineering_Conventions_TODO.md)

项目名称统一为 LAMP；历史验收与冻结方法保留原名称。当前本地工作区为 `D:\Chen\LAMP`；构建、运行和临时输出保存在 D:。
