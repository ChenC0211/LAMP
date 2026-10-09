# G03 通用界面输入与 benchmark 分离报告

完成日期：2026-10-03。生产 Interface 不再按 `problem` 名称选择硬编码问题；现有 24 个模板已迁移为 [通用输入](../api/G03_Interface_Transport_Input.md)，目录名称、时间窗和运行命令保留。

输入通过 `Domain`、`coordinates`、`Initial.phi`、`Velocity`、`Boundary`、`Constants` 和显式 `Refinement/Regrid` 描述。Initial 的空间梯度与 Boundary 的固定时间空间梯度由同一 Core 表达式链式求导；constant 和 Q3 polynomial 速度使用既有不可变 endpoint provider。新增椭圆或支持范围内的速度工况只修改字典，主推进仍调用 `prepareHermiteGalsFirstJetStep`，网格切换仍使用现有 transfer 事务。Core/Geometry 与 Modules/Interface 算法未修改。

已退役三个生产文件：`InterfaceFirstJetBenchmark.cpp/.h` 和 `InterfaceBenchmarkInput.cpp`。原文件、24 个旧字典、旧检查等共 36 项及 SHA-256 清单保存在 [G03 快照](../../runs/baselines/G03_20261003/source_manifest.json)。迁移脚本归档在同一基线的 engineering_scripts 中，不加入 tools 或生产树。当前 App 的配置/数据 contract、输入 reader、推进和 writer 分别归 InterfaceTransport.h、IO/InterfaceTransportInput.cpp、InterfaceTransport.cpp 和 IO/InterfaceFirstJetOutput.cpp。

独立 `cases/benchmarks/hermite_gals/driver.cpp` 与测试参考继续提供解析误差和精度 gate；生产没有专用解析解或静态/平移精度断言。没有解析参考时，CSV 的 analytic_volume/absolute_volume_error 为 NA，accuracy.json 的参考及解析误差为 null，真实体积漂移、曲率和性能保留；生产末态标识为 COMPLETED，不能据此声称 benchmark 精度 gate 通过。effective_config 与 configuration_trace 输出通用输入、实际值和默认系数；导出的 Interface 块已验证能重新读取。

| 验证 | 结果与证据 |
| --- | --- |
| Release 构建、全套 CTest | 70 项注册，68 项有效检查全部通过，0 失败；[构建/检查日志](../../build/g03_full_build_ctest.log) |
| 原有缺失 fixture | legacy_nodal_state_contract、pic_checkpoint_contract 继续 Disabled，不计作通过 |
| App 短契约 | 12 组，48 接受步、36 次 regrid；41 个无效输入、7 个精确范围上下文、2 个多余输入、1 个实际 CFL 拒绝；[验证记录](../../runs/checks/interface_first_jet_app_contract/20261003_101428_c9a5f2ba/verification.json) |
| 表达式契约 | 既有 s/x/y/r/z 及新增时间变量，独立值/偏导、非法语法和奇点检查通过 |
| 通用新工况 | 非单位域椭圆：constant、等价 polynomial、非单位 origin/scale 的 affine shear 三组，独立界面方程/面积检查及 effective config 回读通过 |
| 输出契约 | 18 组输出频率，48 帧、96 个 VTK 文件、6264 个界面样本；数值频率独立、SHA-256、拒绝前缀、现有输出保护、IO 失败、Windows 原子替换重试通过；[记录](../../runs/checks/interface_first_jet_app_contract/20261003_101428_c9a5f2ba/output_contract/verification.json) |
| 全部模板预检 | 24/24；[记录](../../runs/checks/G03/preflight.json) |
| 完整独立 benchmark | 12/12，通过完整窗口；12 行最终误差结果与 G02 冻结运行逐字节相同；[比较](../../runs/checks/G03/benchmark_comparison.json)、[实际运行](../../runs/benchmarks__hermite_gals/20261003_101542_962_c82664fc/run.log) |
| 动态 AMR 涡旋 | 实际 400 步、3 次 regrid；401 行接受态的 14 个数值字段逐项相同，11 个场文件逐字节相同；[比较](../../runs/checks/G03/migration_comparison.json)、[结果](../../runs/hermite_gals_vortex_xy_dynamic_amr/20261003_101336_450_36dd3788/output/fields/solution.pvd) |
| 调试版本 | RelWithDebInfo 重建；RZ 初始 AMR 的 4 步、3 次 regrid 实际通过；[运行日志](../../runs/hermite_gals_static_rz_dynamic_amr/20261003_101805_964_g03_debug_ce6eb014/run.log) |
| 工程约束 | architecture、自测与 repository_layout 通过；目录检查拒绝生产 Interface 中命名 benchmark 定义及 cases 中旧 schema |

动态涡旋峰值相对体积漂移仍为 `2.1375898866666891e-5`。此次将未提供的解析误差显式置为 NA/null，不放宽任何数值算子的接受条件或 benchmark 精度门槛。首次新增双分量常速检查触发既有“角点同时跨两壁未解决”拒绝，该限制保持并已写入 API；随后新增椭圆平移/剪切分别在现有支持范围内验证。

没有重跑全部 24 个模板的完整窗口，已完成 24 组预检、短契约、独立完整 12 组和复杂动态 AMR 400 步比较。空间/时间/曲率/粗化的全部收敛扫描不在此次输入迁移中重复执行。第三方 SuperLU/BLAS 库缺 PDB 的原有调试链接警告保留；LAMP 调试构建及实际推进成功。

范围仍为独立规定速度输运。任意非 polynomial 空间速度输入、周期边界、一般动态 AMR 策略、原生 Hydro/EHD 耦合与低黏/无黏射流破碎物理验收属于后续对应任务，不通过增加专用 problem 分支接入。

```powershell
.\tools\build.cmd
.\tools\run_case.cmd hermite_gals_rotation_xy
.\tools\build.cmd -Preset local-benchmarks
.\tools\run_case.cmd benchmarks/hermite_gals
```
