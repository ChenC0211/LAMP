# I01 / Interface V2 state/API 交付报告

**结论：PASS。** 执行日期：2026-10-01。本轮只完成 I01；下一任务为 I02，Gate I 和完整项目目标均未完成。

任务定义见 [I01 TODO](I01_Interface_V2_State_API_TODO.md)，类型及使用边界见 [API 说明](../api/I01_Interface_V2_State_API.md)。开工前已核对 S08 的 9 个交付文件及 evidence.zip 哈希，Gate S 前置条件满足。

## 实际交付

| 文件 | 变化 |
| --- | --- |
| [HermiteGalsState.h](../../include/Modules/Interface/HermiteGalsState.h) | 新建 vertex/edge first-jet、独立 cell inventory、派生 closure、时间/身份和 candidate/只读 state/geometry view 类型 |
| [HermiteGalsState.cpp](../../src/Modules/Interface/HermiteGalsState.cpp) | 新建数据复制、唯一 snapshot 所有权、身份失效、结构校验及失败隔离 |
| [interface_hermite_gals_state_contract.cpp](../../tests/interface_hermite_gals_state_contract.cpp) | 新建独立生命周期及输入契约测试，使用解析 fixture 和真实 Core snapshot |
| [tests/CMakeLists.txt](../../tests/CMakeLists.txt) | 注册一个新 contract 目标 |

production 文件 2、测试文件 2，符合任务预算。另交付 TODO、本报告、API 说明和文档索引。未修改 Core、旧 InterfaceState、Hydro/Electricity、Apps、reader 或算例输入，没有将 V2 方法接入运行分支。

## 已实现的架构事实

- 主量只有 master vertex 的 `phi/dx/dy`、master edge 中点的 `phi/dn` 和独立 `M_K`；`phi_xy`、bubble 单独存放。XY/RZ 测度明确标识，不 clipping、不从几何覆盖库存。
- 候选复制输入，完成态只读；geometry view 保有同一 snapshot 和状态的寿命，不能通过外部可写 shared_ptr 改写 authority。
- 主量更新分配新身份并清除旧派生结果。时间相同的两轮候选仍有不同身份；过期结果和重复绑定被拒绝。
- `finish()` 仅消费完成的 Interface 容器。App/Coupling 持有并决定替换 accepted handle，不在 Interface 内推进主时钟或耦合提交。
- 相体积负值/非有限值、非法索引/时间、错误 snapshot 家族、mesh 或活动域/原生角点不匹配有明确拒绝；失败不会污染旧接受态。

## 验证结果

证据目录：[runs/baselines/I01_20261001](../../runs/baselines/I01_20261001/)。

| 层级 | 实际检查 | 结果 |
| --- | --- | --- |
| T0 | Release 构建主程序、新 contract、原 Interface contract、architecture 目标 | PASS，最终构建退出 0，9.334 s |
| T0 | `quantum_help`、`architecture` | PASS，0.08 s / 0.73 s |
| T1 | `interface_hermite_gals_state_contract` | PASS，0.06 s；静态只读/类型检查，输入复制、寿命、时间/身份、候选拒绝、失效及绑定失败原子性 |
| T2 | `quantum_interface_state_contract --parametric-state-control` | PASS，退出 0，0.357 s；原有状态生命周期短回归 |
| 范围 | 原有文件哈希及 S08 冻结实现核验 | 1,538 个原有源码/输入/测试/工具等文件未变；唯一既有代码树改动为 tests/CMakeLists.txt；另增 3 个源码/测试文件；69 个冻结专有实现全部匹配 |

最终 CTest 3/3 加 legacy 独立检查 1/1。没有运行或声称通过新的输运 benchmark、收敛测试和耦合生产算例。没有新增 compiler warning；构建保留既有 SuperLU `/LTCG` 自动重链接信息，这不是新 warning。主程序行为未接入新容器，不能用本次可执行文件重新链接后的字节变化推断数值方法变更。

首次 contract 运行未通过：Core 的 moved-from snapshot 使用了与新模块不同的标准异常子类，测试 helper 最初只捕获 `invalid_argument`。修正为捕获标准异常并继续严格核对消息后通过，未放宽生产验收或隐藏错误。日志 `ctest_initial.log` 保留该次失败。复核另修正了 bounds 检查：直接比较原生角点，避免 `center +/- size/2` 的舍入差；新增 `xmin=0.1` fixture 后重建并完成最终检查。

构建脚本直接调用受本机 PowerShell execution policy 限制，因此实际采用脚本对应的 CMake configure/build 命令；最终命令、环境及退出结果见 `verification.json`、`configure.log`、`build_final.log` 和各结果文件。测试执行文件和 DLL 保存在 `runtime/`，由 `runtime_hashes.json` 记录身份。未修改机器 execution policy。

冻结方法 PDF 第 4–5 页公式已渲染核对，并结合配套 TeX 和现行数值规范确认数据拆分。现有 MiKTeX 缺少 Adobe-GB1 CMap，图片中中文缺失，不能把其文本抽取作为完整中文正文证据；该工具限制及日志如实保留。没有修改或重导出权威 PDF。

## 未完成项与边界

`bindGeometry()` 是接收外部闭合结果和 snapshot 的所有权接口，**没有验证 polynomial/volume 的数值一致性**。同源生成、master 完整性/AMR C1、相体积上界、compact closure、flux/reflux、bubble topology guard 等仍由 I02–I11 实现和验证。当前完成态表示容器已装配，不表示数值 candidate 已被接受。API 说明列明生产者前置条件，测试没有把任意同尺寸几何当作数值通过证据。

未注册 `conservativeHermiteGals`；未发布 `alpha/Q_l`，未接 Hydro 新主线。S01 的 `hydro_moment_xy` 首步失败仍是已知问题，本轮短回归没有将其改记为通过。低黏/无黏强变形、射流和破碎目标仍需按路线图推进。

归档 `evidence.zip` 保存构建/测试/校验与交付副本；`delivery_manifest.json` 记录归档及交付 SHA-256。Git 元数据不可用，未伪造 commit。本轮不开始 I02。
