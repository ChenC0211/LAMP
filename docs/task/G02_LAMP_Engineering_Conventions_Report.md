# G02 / LAMP 工程规范与清理报告

完成：2026-10-03。任务依据是用户对测试分类、旧界面路线退役、目录约束和 LAMP 更名的明确要求；承接 G01。

## 当前入口与目录

程序为 **LAMP — Levelsetbased Adaptive Multiphysics Platform**，构建目标 `LAMP`，发布程序 `bin/LAMP.exe`，带符号程序 `bin/relwithdebinfo/LAMP.exe`。工作区路径仍为 `C:\Chen\QUANTUM`；现存 `Quantum` C++ 名字空间、`quantum_*` 库目标及历史格式标识保留源码/格式兼容，冻结证据不改名改值。

```powershell
.\tools\build.cmd
.\tools\run_case.cmd hermite_gals_vortex_xy_dynamic_amr
.\tools\build.cmd -Test
.\tools\build.cmd -Preset local-benchmarks
.\tools\run_case.cmd benchmarks/hermite_gals
```

普通运行、完整研究、预检和调试准备使用同一启动器，选择名称即可。`cases/catalog.json` 列出 51 个入口及能力状态；运行目录为 `runs/<名称>/<时间戳_标识>`。直接对 cases 模板执行物理求解会以 2 退出，预检仍允许读取模板。

docs 根目录现有 README 和 00–07 八份关键约束。新增三份薄约束明确目录/文档、生产/兼容/测试、运行/调试的归属；API、任务、路线图、方法 PDF、使用指南各归子目录。单次调试和归档不进入生产树。

tests 根目录仅有 CMakeLists、README、布局检查；实现归 contracts 的模块子目录，辅助代码归 support。完整 FirstJet 时间窗和收敛驱动归 cases/benchmarks，默认不进入 CTest。Hysing 原始参考归对应 cases；原始数据未修改。22 份旧算例详细说明归 docs/history_todo/cases，当前 README 统一提供有效命令与能力状态。

tools 根目录仅保留 build、run_case 的 cmd/ps1 和索引；后处理、文档构建归各自子目录。旧专属运行/绘图脚本退役，文档工具生成物归 runs/checks。

## 代码退役与兼容边界

- 删除 37 个生产源码/头文件，包括旧独立 Interface App 的 levelSet、保守几何和参数曲线分发/输出/专属驱动，以及 InterfaceComponentMotion、InterfaceConservativeGals 和三种旧 InterfaceQuartic 实现。
- 当前 Interface App 仅保留四个实现文件：入口、严格输入、FirstJet 驱动、同源输出。模块根目录保留 19 个实现、19 个公开头；内部辅助归 Detail。
- Hydro/EHD 仍需的状态、边界、载荷和输运适配归 Compatibility，独立 owner 为 lamp_interface_compatibility。没有伪造尚未实现的新 GALS–Hydro 耦合来替代它。
- 通用几何测量的旧 InterfaceState 适配拆到兼容层；原生积分算术仍只有一份。当前 GALS 库不能链接兼容层，不能新增兼容消费者；76 条既有 include 边逐条冻结，其中 71 条 Interface、5 条旧 Electricity。
- 旧大型 harness 拆成短 FirstJet 契约与完整研究驱动。43 个旧测试/脚本文件退役；无配置的旧 App 对比注册移除。历史 GALS 守恒研究 API 仍有局部算子消费者，但没有生产 App 选择入口。

整理前快照在 `runs/baselines/G02_20261002/workspace_source_before.zip`。migration、test_migration、case_retirement、case_documentation_migration 和 runtime_migration 清单记录原路径与哈希。旧 QUANTUM 二进制/PDB 归 runtime_legacy。迁移时校验原字节；随后当前 README、入口和品牌的有意修订不冒充字节保持。冻结方法 PDF 和历史结果不改写。

`source_audit.json` 对整理前 918 个生产文件审计：873 个除 include 路径/行尾外完全一致，37 个明确退役，8 个有意修改。8 项限于入口收敛、品牌、输入声明和测量适配拆分；HermiteGals 数值实现未改动。

## 验证结果

| 验证 | 结果与证据 |
| --- | --- |
| Release 全目标构建 | PASS；build/g02_release_build_final.log；最终品牌重建另记 g02_release_brand_final.log |
| RelWithDebInfo 构建 | PASS；build/g02_debug_build.log、g02_debug_brand_final.log；静态 RZ 动态 AMR 实跑 4 步、3 次网格事务通过 |
| 架构与布局 | PASS；13 个 owner、370 个生产实现、884 个文件；120 条策略自检、23 个端到端错误变体均按预期处理 |
| 默认 CTest | 注册 69 项，67 项实际执行均已通过；2 项缺原始 fixture 为 Disabled，不计通过 |
| CTest 复核范围 | g02_ctest.log 中 66 项通过，1 项旧方法接受断言失配；按已退役行为改为明确拒绝后，g02_ctest_app_final.log 单项通过。没有重复运行不受影响的 66 项或放宽数值门槛 |
| App 契约 | 12 个短算例、48 步、36 次 regrid；24 个非法输入、15 个上下文诊断、6 个范围预检；1 个物理拒绝保留接受状态；3 个退役方法明确拒绝 |
| 通用入口实跑 | 静态 XY、静态 RZ 动态 AMR、XY 平移、动态 AMR 涡旋通过；预检与 PrepareOnly 不推进时间 |
| 数值对照 | 动态 AMR 涡旋 400 步、3 次 regrid；accuracy.json 与 G01 运行逐字一致 |
| benchmark 驱动 | 短动态运行 12/12、完整固定网格时间窗 12/12 通过；完整研究接受 1036 步，退出码 0 |
| 工具与配置 | JSON/PowerShell 语法通过；Hydro 后处理 helper 编译通过；模板直接物理求解拒绝且未生成 cases/output |

动态 AMR 涡旋对照：位置误差 6.1108866740340337e-6、梯度 2.4299317790353199e-4、曲率 8.3748682055109214e-3、绝对体积误差 3.8134720931753963e-7、峰值相对体积漂移 2.1375898866666891e-5。

主要运行位置：

- `runs/hermite_gals_vortex_xy_dynamic_amr/20261003_001819_865_7909cdea`（完整场输出与数值对照）。
- `runs/benchmarks__hermite_gals/20261003_002007_219_a1a04fcc`（完整研究，CSV 行在 run.log）。
- `runs/benchmarks__hermite_gals_smoke/20261003_001946_333_ffe20956`（短研究）。
- `runs/hermite_gals_static_rz_dynamic_amr/20261003_003936_375_debug_69567b1f`（带符号实跑）。

两个缺失的原始 fixture 检查为 legacy_nodal_state_contract 和 pic_checkpoint_contract；恢复相应原文件后才启用。既有 Hydro 索引转换警告和第三方 SuperLU/BLAS 缺 PDB 警告仍存在，不影响本次构建。未实际启动 IDE 断点会话；已验证带符号构建、真实运行和 F5 准备目录。

本次没有执行 space/time/curvature/coarsening 全部收敛研究，也没有完成新的 Hydro/EHD 物理 gate。工程收口不表示低黏/无黏、射流或破碎的完整平台目标已实现；下一项仍按路线图推进 H01。

临时工程脚本已在交付时归档至 runs/baselines/G02_20261002/engineering_scripts；不是后续用户入口。
