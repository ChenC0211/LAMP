# I01 / 定义 Interface V2 state/API

任务依据：[05 分模块任务清单](../roadmap/05_LAMP_Module_Task_Roadmap.md) I01。执行日期：2026-10-01。前置条件为 S08 报告中的 Gate S；开工时复核交付文件与归档哈希。

## 目标

只解决一个架构事实：V2 的主状态由 transported first-jet skeleton 与独立守恒相体积组成，派生闭合量及唯一只读 Hermite snapshot 由 Interface 持有，候选态不能修改接受态。现有运行路径无数值行为变化。

## 允许修改

- `include/Modules/Interface/HermiteGalsState.h`、`src/Modules/Interface/HermiteGalsState.cpp`：类型、所有权、时间层、身份及输入结构校验。
- `tests/interface_hermite_gals_state_contract.cpp`、`tests/CMakeLists.txt`：独立契约验收及注册。
- 本 TODO、`I01_Interface_V2_State_API_Report.md`、`docs/I01_Interface_V2_State_API.md`、`docs/README.md` 的任务索引。
- 证据存入 `runs/baselines/I01_20261001/`；允许正常构建更新 `build/`、`bin/`。

预算：2 个 production 文件、2 个测试文件。

## 禁止修改

- S08 冻结实现、旧 InterfaceState、Core/Solvers、Hydro/Electricity/Plasma、Apps、reader、默认方法和 cases。
- 本轮不实现 I02–I12 的 tracer、输运、compact closure、AMR restriction、FV/reflux、bubble 求根、几何生成及 benchmark。
- 不把 V2 方法注册为可运行选项，不加入旧算法 fallback，不放宽现有容限，不修改前序证据。

## 数学约束

- master vertex 主量为 `phi, phi_x, phi_y`；master edge 中点主量为 `phi, partial_n phi`，导数为物理导数。hanging 数据由 I06 生成，不能作为独立推进未知量。
- `M_K` 是独立负相体积，不从 snapshot 重算覆盖，不 clipping。具体 XY/RZ 测度显式标识；体积初始化、上界及几何闭合验收属于 I07/I08/I10。
- `phi_xy` 和 `c_K` 独立于主状态，分别为 compact closure 和 scalar-volume closure 的结果。
- snapshot 继续使用 Core 的 native enriched-quartic 表示，只向消费者提供携带 geometry identity/time layer 的只读句柄。
- App/Coupling 决定候选是否接受；Interface 不推进主时钟、不私自提交耦合状态。

## 实现要求

- [x] 复核 Gate S、方法 PDF 中式 (11)–(24) 及六份规范；记录开工前源文件身份。
- [x] 新建独立类型，明确 native 索引、master edge 正向法向、稀疏活动几何与 cell inventory 的归属。
- [x] 候选构造复制输入；只读完成态和几何视图保有对象寿命；不接受外部可写 shared_ptr 别名。
- [x] 主量替换产生新身份并清除旧派生结果；过期/异源身份、未绑定候选、重复绑定、非法时间或 native mesh 不匹配明确报错。
- [x] 绑定接口只接收外部已算出的闭合量与 snapshot，不生成几何，不声称验证它们的数值一致性。I11 必须建立和验证同源发布；本轮用标明用途的解析 fixture 验证所有权。
- [x] 完成 Release 构建、契约测试、架构检查和既有 Interface 状态短回归，记录实际范围。

## 验收

- **T0：** 构建 `QUANTUM`、`quantum_interface_hermite_gals_state_contract`、`quantum_interface_state_contract`、`quantum_architecture_check`；新头文件可作为测试首个 include；无新增 warning；architecture 与 quantum_help 通过。
- **T1：** 主量与闭合量分离的编译期检查；输入/候选/完成态/保留几何句柄的寿命与隔离；M 独立保存；不同 candidate 即使时间相同也有不同身份；更新主量使派生量失效；错误身份、时间、索引、尺寸、非有限值、负体积及错误 snapshot 家族被拒绝。保留主量原值，不用 clipping 修复。
- **T2：** 运行原有 `quantum_interface_state_contract --parametric-state-control`，验证 legacy 短流程未改变；核对旧源码及输入哈希。没有新增输运算法，本轮不声称完成 V2 数值 benchmark。

## 失败规则

- 若必须修改旧算法、跨模块数值接线或超出预算才能定义 API，记 `BLOCKED` 并说明接口缺口，不加 fallback。
- 本次引入的编译/所有权/回归错误记 `FAIL`，在本任务范围内修复；数值闭合没有实现不能伪装为通过。
- S01 已知 Hydro 首步失败继续保留，本任务不修复也不重新标为通过。

## 交付

类型与生命周期实现、独立契约测试、API 边界说明、TODO/报告、构建/测试/文件身份证据。结束状态仅使用 `PASS / FAIL / BLOCKED`；下一任务是 I02，不在本轮开始。

**最终结论：PASS。** 2 个 production 文件、2 个测试文件；Release 构建、最终 CTest 3/3 和 legacy 检查 1/1 通过。数值算法、几何同源发布及 Gate I 未完成。详见 [I01 报告](I01_Interface_V2_State_API_Report.md)。
