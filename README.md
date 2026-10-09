# LAMP

**Levelsetbased Adaptive Multiphysics Platform**：基于水平集的自适应多物理平台。

当前已验收的界面入口是纯 Hermite–GALS / FirstJet，支持 XY/RZ、固定网格和动态 AMR 的规定速度验证。完整低黏/无黏流体、电场与射流破碎耦合仍按路线图开发。

```powershell
.\tools\build.cmd
.\tools\run_case.cmd -List
.\tools\run_case.cmd hermite_gals_vortex_xy_dynamic_amr
```

换算例只改末尾名称。配置在 `cases`；每次结果在 `runs/<算例>/<时间戳_标识>`；ParaView 打开该目录的 `output/fields/solution.pvd`。

- [使用与调试](docs/guides/Run_And_Debug.md)
- [GitHub 与 VS Code 开发流程](docs/guides/Git_Workflow.md)
- [算例索引与能力状态](cases/README.md)
- [短检查](tests/README.md) · [完整 benchmark](cases/benchmarks/README.md)
- [开发约束](docs/README.md) · [现行任务路线图](docs/roadmap/05_LAMP_Module_Task_Roadmap.md)
- [G02 整理任务](docs/task/G02_LAMP_Engineering_Conventions_TODO.md)

项目名称已由 QUANTUM 改为 LAMP；历史验收与冻结方法保留原名称。当前本地工作区为 `D:\Chen\LAMP`；构建、运行和临时输出保存在 D:。

