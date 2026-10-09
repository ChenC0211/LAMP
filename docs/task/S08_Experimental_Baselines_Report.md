# S08 / 实验基线冻结交付报告

**结论：PASS；Gate S 通过，可在下一任务开始 I01。** 执行日期：2026-10-01。任务定义见 [S08 TODO](S08_Experimental_Baselines_TODO.md)。

## 冻结边界

冻结的是三条旧方法的开发身份：仅保留显式 A/B、旧 case 回归和比较用途，不再扩展 V2 功能。方法仍可按原配置运行；本任务没有禁用程序、改变默认值、重命名字段或添加 runtime 分支。

| 实验 baseline | 已记录的专有实现文件 | V2 目标与边界 |
| --- | --- | --- |
| `parametricGals` | 20 | `conservativeHermiteGals` 必须另建 first-jet skeleton 与独立 cell phase volume，不延伸全局参数曲线/H² 状态 |
| `CompatibleMoment/cchMomentMidpoint` | 38 | native staggered CCH 使用 face/cell 状态与局部 compatible Hodge，不继续扩大全局 C1/Q8 速度状态 |
| `sharpHermiteP3` | 11 | `sharpHHO2` 按 SI-HHO-Q 独立建设，不把旧 P3 FV 修正逐步包装成 HHO |

可核对的文件清单、冻结前后 SHA-256、方法选择字段及代码标记点见 [S08_Frozen_Baselines.json](S08_Frozen_Baselines.json)。69 个专有文件涵盖 Interface 的 `InterfaceParametric*` 实现/应用适配、Hydro 的 `HydroCompatibleMoment*` 与 `HydroMomentCase`、Electricity 的 sharp/fitted/旧电荷路径及专有 Detail 头文件。

另列 12 个 reader、dispatch 与共享参数曲线几何文件作为接线/身份参考；这不把整个 Core/Solvers 或 reader 的其它分支冻结。后续 V2 必须继续按唯一 owner 和独立状态规则实施，不能从本次清单推导新的全局参数曲线生产路径。

本清单与代码注释是开发约束及审计依据，不是运行时锁或自动拒绝所有后续代码修改的机制。新 V2 路径不能隐式回退到旧路径，也不能共享旧方法的可写状态。

## 实际修改

以下五个 production 文件仅增加 `S08_FROZEN_BASELINE_BEGIN/END` 注释块：

1. [InterfaceParametricProjection.h](../../include/Modules/Interface/Compatibility/InterfaceParametricProjection.h)：标记全局参数曲线/H² 状态。
2. [InterfaceParametricTransport.h](../../include/Modules/Interface/Compatibility/InterfaceParametricTransport.h)：标记参数曲线输运/通量路径。
3. [HydroCompatibleMomentFlow.h](../../include/Modules/Hydrodynamics/HydroCompatibleMomentFlow.h)：标记旧速度历史 owner。
4. [HydroCompatibleMomentSpline.h](../../include/Modules/Hydrodynamics/HydroCompatibleMomentSpline.h)：标记全局 C1/Q8 速度表示。
5. [SharpDielectricElectrostatics.cpp](../../src/Modules/Electricity/SharpDielectricElectrostatics.cpp)：标记普通和 fitted sharp P3 求解实现。

没有修改接口声明、状态布局、算子、配置默认值、reader、dispatch、容限或输入。所有注释移除后，原内容与任务前快照一致（仅忽略 CRLF/LF）。未添加测试文件。

## 验证与证据

证据目录：[runs/baselines/S08_20261001](../../runs/baselines/S08_20261001/)。任务前五个文件和原始运行程序保存在 `before/`。

| 层级 | 检查 | 结果 |
| --- | --- | --- |
| T0 | Release 主程序及三个旧方法 contract 目标构建 | PASS：退出 0，411.805 s |
| T1 | 注释以外的代码不变；其它源码/输入/测试/工具/构建配置不变；冻结清单身份匹配 | PASS：5 个注释文件、1,086 个其它文件；文件清单未变，69 个专有实现及 12 个引用文件哈希匹配 |
| T2 | 旧参数曲线状态、旧 Hydro 时间验证、旧 Electricity 生命周期/解析场检查 | PASS：3/3，均退出 0 |

`source_before.csv`、`comment_only_verification.csv`、`static_verification.json` 记录 T1 细节。构建命令、退出码和耗时记录在 `build.log` / `build_result.json`。短验证的选择在执行前已经写入 TODO，不以结果临时改选更容易通过的检查。

| 原有验证命令 | 墙钟 | 观察 |
| --- | --- | --- |
| `quantum_interface_state_contract --parametric-state-control` | 0.425 s | 参数曲线唯一表示、native identity、物理时钟、旧状态及失败隔离通过 |
| `quantum_hydro_module_contract --moment-time` | 95.910 s | 两组时间收敛阶分别为 2.000142、2.075674，满足原先 `1.8 < order < 2.2`；接受/拒绝及 variable-dt 检查通过 |
| `quantum_electricity_evolution_contract` | 3.264 s | 旧电荷输运策略、电学生命周期、解析场及发布对象寿命检查通过 |

测试执行文件保存在 `runtime/`，各自哈希、参数、线程与退出码见 `checks.json`；完整日志和 Hydro 时间历史在 `checks/`。固定 OMP 2 线程、BLAS/MKL 1 线程，详细设置见 `verification_environment.json`。电学日志保留 OpenMP team-size 提示，程序退出 0。以上单次墙钟不是性能优化的 A/B 证据。

构建时 `HydroCompatibleMomentSplineTransport.cpp` 出现 Eigen Index 到 int/const int 的 C4244 提示。该数值文件本轮哈希未变，已另存旧 `batched_momentum_compile.log` 中的同类警告为 `historical_warning_build.log`。本任务不扩展范围修复该既有数值文件，也不声称构建零警告。

S01 已知的 `hydro_moment_xy` 首步失败仍保留；本次短验证不能将它改记为通过，也不能证明任何 V2 数值 gate 完成。

## Gate S 逐项核验

S02–S06 是既有规范包的只读核验，本轮没有另开这些任务。详见 `gate_S_prior_audit.csv`。

| 条目 | 当前证据 | 结论 |
| --- | --- | --- |
| S01 | 基线源码/二进制/输入归档、构建类型及测量报告；本轮重验两份归档哈希 | PASS：已知失败如实保留 |
| S02 | docs/README 可直接访问六份规范与三份冻结方法；总纲 §2 指定权威与 05 唯一任务入口 | PASS |
| S03 | 架构规范 §1–4、§6 明确依赖、owner、只读视图与禁止隐式 fallback | PASS |
| S04 | 数值规范包含 I-R1..5、H-R1..7、E-R1..7 与跨模块时层/接受态要求 | PASS：规范已建立，不代表实现通过 |
| S05 | 输出规范 Rev.1 包含 accepted/rejected、VTK XML/PVD、series/meta、残差和 failure 语义 | PASS：规范已建立，旧输出尚未迁移 |
| S06 | 总纲 T0–T4 验证层级、05 的任务格式、S/I/H/IH/E/C 顺序及各 gate | PASS |
| S07 | 本轮开工前重验 23 份交付文档与证据哈希；14 份历史 Markdown 的顶部标记 | PASS |
| S08 | 五处代码身份、冻结清单、不变性、构建及旧路径回归 | PASS：5 个纯注释改动，构建和 3/3 检查通过 |

三个冻结方法 PDF 的实际内容哈希与 S01 相同。最终逐项表见 `gate_S_final_audit.csv`。**Gate S 已满足**，其含义只是允许开始 I01 的 state/API 任务；不允许跳过 I02–I12，亦不证明低黏/无黏电喷雾主目标已实现。

## 交付与范围

production 文件 5，测试文件 0；新增本 TODO、报告、JSON 清单，更新文档索引。`build/`、`bin/` 仅由正常构建更新，前序归档与报告保持原样。`evidence.zip` 保存本次运行证据与交付副本，哈希及逐文件身份见 `delivery_manifest.json`。

本轮只做 S08。下一项为 I01，须另写 TODO；本轮未创建 V2 API 或修改算法。
