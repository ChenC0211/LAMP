# G04 Interface 字典职责固化报告

完成日期：2026-10-03。独立 Interface 已从 G03 单文件输入迁移为五文件输入，遵循原字典规范的核心思想：按数据职责分文件、每个值只有一个声明位置、严格拒绝无效输入、记录真实生效配置与来源。当前规范见 [字典约束](../03_LAMP_Dictionary_Rules.md)和 [Interface schema](../api/G04_Interface_Dictionary_Schema.md)。

| 文件 | 当前唯一职责 |
| --- | --- |
| system/run.dict | Run.app、Time.units/dt/endTime、Output.Fields；总步数由公共 fixedStepCount 一次推导 |
| system/geometry.dict | Geometry.coordinates、Domain、Initial.phi |
| system/mesh.dict | Mesh.resolution、initialMesh、Refinement、Regrid |
| system/fields.dict | Fields.Velocity 及共享输入表达式参数 Constants |
| system/interface.dict | Interface.method 与 Boundary |

24 个完整模板、六组最小 fixture、输入与输出契约已同步迁移。没有重复 key 的兼容别名，旧单文件配置、problem/Benchmark 和错放字段直接拒绝。独立 Interface 不消费的材料、solver 或其他物理模块文件也不能作为占位输入。Fields.Constants 只声明输入表达式共享参数，供 Initial、Boundary 和 Velocity 使用，不承载网格、主时钟、材料或求解器控制。

App reader 按文件解析后汇总 typed config，数值内核仍只接收已验证数据。沿用 Core 表达式、constant/Q3 polynomial endpoint provider、Hermite–GALS/FirstJet 推进及现有 AMR transfer 事务。时间输入改为 dt/endTime，复用公共固定时步计数规则；未改变 dt、实际步数、重网格时刻或数值接受条件。

五份输入都记录 SHA-256；configuration_trace.csv 的 file/line 对应真实文件，默认值与 consumer 保持统一追踪。跨文件范围错误定位到实际 geometry/mesh/run 字段。effective_config.dict 按文件汇总真实生效值，output/meta/effective 下同时导出可直接回读的五份同名字典，由同一个 writer 生成。

薄开发约束 [03](../03_LAMP_Dictionary_Rules.md)重新明确全项目文件职责与迁移规则；完整字段、范围、默认与可运行示例放在 G04 API。索引、运行说明与模板说明已更新；G03 及旧命名 benchmark API 标记为历史布局，不作为当前输入示例。

| 验证 | 结果与证据 |
| --- | --- |
| Release 构建 | 成功；[日志](../../build/g04_build.log) |
| 受影响 CTest | architecture、自测、repository_layout、runtime_registry、interface_first_jet_app_contract 共 5/5 通过，0 失败；[日志](../../build/g04_validation.log) |
| App 短契约 | 12 组、48 接受步、36 次 regrid；49 个无效输入、7 个跨文件范围上下文、4 个缺失模块输入文件、1 个多余文件和 1 个实际 CFL 拒绝；[记录](../../runs/checks/interface_first_jet_app_contract/20261003_104955_c0260515/verification.json) |
| 输出与通用输入 | 18 组输出频率、3 组非单位域椭圆，独立验证 48 帧、96 个 VTK 文件、6264 个界面样本；五输入独立 SHA-256、五份生效字典回读、数值频率独立、拒绝前缀、现有输出保护、IO 失败及原子替换重试通过；[记录](../../runs/checks/interface_first_jet_app_contract/20261003_104955_c0260515/output_contract/verification.json) |
| 全部模板预检 | 24/24；[记录](../../runs/checks/G04/preflight.json) |
| 动态 AMR 涡旋完整窗口 | 400 接受步、3 次 regrid；与 G03 运行相比，401 行接受态的 14 个数值字段完全相同、11 个场文件逐字节相同；[比较](../../runs/checks/G04/migration_comparison.json)、[运行](../../runs/hermite_gals_vortex_xy_dynamic_amr/20261003_105542_410_a1d6db2f/run.log) |
| 调试版本 | RelWithDebInfo 构建成功，RZ 初始 AMR 的 4 步、3 次 regrid 实际通过；[构建](../../build/g04_debug_build.log)、[运行](../../runs/hermite_gals_static_rz_dynamic_amr/20261003_105659_978_g04_debug_358cf347/run.log) |
| 文档示例 | 从 G04 API 的五段字典原文生成输入，预检及 30 步旋转实际完成；[记录](../../runs/checks/G04/doc_example_verification.json)、[运行日志](../../runs/checks/G04/doc_example_run.log) |

迁移前 63 项输入、reader、writer、检查及相关文档保存在 [G04 快照清单](../../runs/baselines/G04_20261003/source_manifest.json)。迁移脚本归档到同一基线的 engineering_scripts，不进入 tools 或生产目录；活动 cases 仅保存输入模板。

本次验证针对输入布局与 IO 修改，未重复运行全套 CTest、全部 24 个模板的完整窗口或全部空间/时间/曲率收敛扫描。第三方 SuperLU/BLAS 原有缺 PDB 调试链接警告保留；LAMP 调试构建和实际推进成功。解析参考与精度 gate 仍属于独立 benchmark，生产 COMPLETED 仅表示配置窗口完成，解析误差继续为 NA/null。

此次固化范围为独立规定速度的 Interface 入口。Hydro/EHD 冻结 reader 没有全局替换；后续按同一职责规则逐项迁移并验证。原生耦合、任意外部速度 reader、一般动态 AMR 策略和低黏/无黏射流破碎物理验收仍属于后续开发任务，既有角点同时跨两壁的轨迹限制保留。

运行命令保持不变：

```powershell
.\tools\build.cmd
.\tools\run_case.cmd hermite_gals_rotation_xy
.\tools\run_case.cmd hermite_gals_vortex_xy_dynamic_amr
```
