# hermite_gals_rotation_xy_dynamic

纯 Hermite–GALS / FirstJet 规定速度算例：通用表达式输入，基础 8×8，初始 uniform，三次动态网格切换；dt=0.01，30 步，终时 0.3 code units。

在项目根目录运行：

```powershell
.\tools\run_case.cmd hermite_gals_rotation_xy_dynamic
```

配置在 [system](system)，可复制整个算例并改名以创建新工况。结果在 `runs/hermite_gals_rotation_xy_dynamic/<时间戳_标识>/output`，ParaView 打开 `fields/solution.pvd`。预检加 `-CheckOnly`，只准备输入/程序快照加 `-PrepareOnly`。

[通用输入与算例说明](../interface_hermite_gals/README.md) / [当前验收范围](../../docs/task/I12p3r2_FirstJet_Gate_Final_Report.md)。

[G04 schema](../../docs/api/G04_Interface_Dictionary_Schema.md): run / geometry / mesh / fields / interface.
