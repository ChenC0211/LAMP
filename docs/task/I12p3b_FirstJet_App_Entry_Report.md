# I12p3b / 纯路线正式 App 与字典入口报告

**PASS；完整 Gate I 仍待 I12p3c/p3d。** [任务](I12p3b_FirstJet_App_Entry_TODO.md)已让主程序通过严格字典运行六类 FirstJet 问题和动态网格，数值内核保持不变。可运行字典和 CLI 见 [API](../api/I12p3b_FirstJet_App_API.md)。

修改 5 个 production 文件、2 个 test 相关文件：新 `InterfaceFirstJetBenchmark.h/.cpp`，既有 `InterfaceGeometryBenchmark.h/.cpp` 与 `InterfaceBenchmarkInput.cpp`；新 `interface_first_jet_app_contract.cmake` 和 tests/CMakeLists。没有修改 Core/Modules、数值参数、旧 benchmark 或冻结输入。新入口显式使用 `hermiteGals`，复用 p3a/p2b；原方法默认值与旧 reader 分支保留。实际接受态的网格、publication、时钟由 App 成对提交。

| 验证 | 结果 |
| --- | --- |
| 39 targets Release 构建 | PASS，23.134 s；主程序与先前 App 专项测试的 SHA 相同；无新增 compiler warning |
| 原 36 CTest | 全通过，368.692 s；35 个数值测试完整输出逐字相同 |
| 新 App CTest | PASS，21.841 s；12 短算例、48 接受步、36 regrid、24 输入拒绝、1 运行 CFL 拒绝、3 旧入口预检 |
| 实际主程序完整矩阵 | 24 组、2,072 接受步、36 regrid；463.425 s（包含每组预检和运行） |
| 全矩阵与模块档案对照 | 所有比较的数值 token **逐字相同**，未引入容差 |
| 历史 crossing probe | 仍退出 1，完整输出相同；不冒充纯路线通过证据 |
| S08 / PDF | 69 个冻结实现、11 个未改接线参考和 3 PDF 哈希相同；1 个 reader 参考按授权新增分支，旧行为单独核验 |

实际主程序矩阵为 fixed/dynamic × 六问题 × uniform/AMR，resolution=8，Q2 初值：静态 dt=.005/4步、平移 .01/40步、旋转 .01/30步、涡旋 .001/400步。dynamic 在四分之一附近的原三步序号切换，保持原步长、终时和参数。逐个核对每个接受步的 step/time/V0/V/dV/V，最终位置/梯度/曲率/解析体积误差/峰值漂移，以及每次重网格的 step/time/kind/目标网格/细化粗化数量/域归一化体积变化。固定参考为 p1，动态参考为 p2b。原逐 cell 体积误差仍由未改的模块回归检查，本 App 最终摘要没有另复制该计算。

App 短测的 AMR 初值使用非多项式，覆盖六类问题的解析边界与真实内部输运。运行 CFL 负例在首步实际 tracer 拒绝，accepted_step/time 保持 0，没有 `[Step]`；`--check-case` 全程不建网格或输出。输入负例保留未知/重复/no-op 配置、类型和数字错误、非法重网格步序、额外 Run 字段或文件等。`writeFields=true` 目前明确失败；stdout 有几何体积与漂移，inventory/flux ledger 明确为 NA。

架构统计由 386 TU/914 文件/3133 includes 变为 **387/916/3143**；12 targets、5 条既有 legacy edges 不变。39 个运行文件中仅 QUANTUM.exe 改变，38 个检查可执行文件哈希不变。p3a 归档哈希保持。

审计纠正：初版 TODO/audit 沿用“81 freeze”简写，把 S08 的 12 个接线参考误当成不可修改实现，导致脚本第一次报错。S08 [报告](S08_Experimental_Baselines_Report.md)原文明确“这不把整个 Core/Solvers 或 reader 的其它分支冻结”。本次 reader 原本就在允许的 5 文件内；审计现区分 69 冻结实现与 12 参考，从 S01 原档案取出 reader，验证其哈希及旧函数主体逐字保持（移除新增纯分支并还原共用验证调用后）；共用 Run/文件验证去除追踪附加语句后与原逻辑一致。没有改变历史清单，也没有绕过冻结方法限制。

开发中失败如实保留：首次构建误用 publication report 成员名，改为既有 `maxAbsCurvature` 后通过；首次 App CTest 的未知方法负例混入新模式专有 dt/steps，严格 reader 更早拒绝未知字段，未到该测试期望的方法枚举错误；已用旧兼容字段隔离该负例，所有输入仍严格拒绝。初次及修正后测试目录、日志和退出码均入档，未删除失败运行。

机器证据在 [本任务目录](../../runs/baselines/I12p3b_20261002/)：`verification.json`、`source_audit.json`、`regression_comparison.json`、`full_comparison.json`；`full_cases.json` 指向 24 组独立输入/日志。`package_evidence.ps1` 归档全部交付、嵌套算例证据、两次 App CTest 目录及运行文件，逐条复核 live/archive SHA；最终标识在 `archive_verification.json`。

这是独立规定速度 App，不等于电喷雾或低黏 Hydro 验收。下一项按单任务顺序实施 I12p3c 的同源只读输出，随后 I12p3d 效率与规模证据，再复核 Gate I。
