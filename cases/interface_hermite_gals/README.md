# 纯 Hermite–GALS 水平集配置

当前 `method hermiteGals` / FirstJet 独立规定速度入口；六类 × 固定/动态 × 初始 uniform/amr 共 24 组完整通过记录见 [Gate I 报告](../../docs/task/I12p3r2_FirstJet_Gate_Final_Report.md)。不求解流体或电场，几何体积漂移作诊断。本目录是索引，实际配置是 cases 下的 24 个独立目录。

## 通用运行命令

```powershell
.\tools\run_case.cmd hermite_gals_vortex_xy
.\tools\run_case.cmd hermite_gals_translation_xy
.\tools\run_case.cmd hermite_gals_vortex_xy_dynamic_amr
.\tools\run_case.cmd -List
```

结果统一在 `runs/<名称>/<时间戳_标识>/output`；控制台打印 ParaView 的 `fields/solution.pvd` 路径。预检加 `-CheckOnly`，准备调试快照加 `-PrepareOnly`。完整规则见 [固定运行指南](../../docs/guides/Run_And_Debug.md)。

## 六类基础配置

| 基础名称 | 配置的初始界面/速度 | dt | 步数 | 终时 |
| --- | --- | ---: | ---: | ---: |
| `hermite_gals_static_xy` | 圆 / 零速度 | .005 | 4 | .02 |
| `hermite_gals_static_rz` | 球 / 零速度 | .005 | 4 | .02 |
| `hermite_gals_translation_xy` | 圆 / constant (.2 0) | .01 | 40 | .4 |
| `hermite_gals_translation_rz` | 球 / constant (0 .2) | .01 | 40 | .4 |
| `hermite_gals_rotation_xy` | 偏心圆 / polynomial 刚体旋转 | .01 | 30 | .3 |
| `hermite_gals_vortex_xy` | 偏心圆 / 时变 polynomial 涡旋 | .001 | 400 | .4 |

每个基础名称均有四种配置：基础名称为固定均匀网格，追加 `_amr` 为固定局部加密，追加 `_dynamic` 为初始 uniform 动态切换，追加 `_dynamic_amr` 为初始 amr 动态切换。例如 `hermite_gals_translation_rz_dynamic_amr`。基础分辨率均为 8×8；动态切换在总步数的 1/4、1/2、3/4 处向下取整，细化 exact、粗化 FirstJetProjection。

模板域为 [0,1]²，code units；XY 为面积，RZ 为完整轴对称体积。输入按 [五文件 schema](../../docs/api/G04_Interface_Dictionary_Schema.md)固化，不使用 problem 名称。字段含义见 [输出 API](../../docs/api/I12p3c_FirstJet_Output_API.md)。初始/末态始终输出，其余由 run.dict 的 Output.Fields.writeEverySteps 决定。更改分辨率须同步选择满足实际半单元 CFL 的 dt，输入预检不能替代实际运行检查。

## 调整工况

编辑或复制对应配置的 system：run.dict 修改 Time.dt/endTime 和输出，geometry.dict 修改域/坐标/Initial，mesh.dict 修改基础网格及加密，fields.dict 修改 Constants/Velocity，interface.dict 只改 method/Boundary。总步数由 dt/endTime 推导，不另设 steps；不要添加生产问题枚举。程序失败保留真实最后接受态和原因，不能通过缩短窗口声称完整回逆通过。

生产输出的解析误差为 NA/null；独立 [benchmark](../benchmarks/README.md)负责解析误差与精度 gate。旧 `run_hermite_gals.cmd` 已退役，使用上面的固定入口。原归档批量结果在 `runs/interface_hermite_gals`，当前新结果按具体配置名称保存。
