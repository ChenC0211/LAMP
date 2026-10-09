# I03 / master vertex first-jet 输运交付报告

**结论：PASS。** 执行日期：2026-10-01。本轮只完成 I03，下一任务为 I04。Gate I 与低黏/无黏电喷雾完整目标仍未完成。

范围与预设验收见 [TODO](I03_Vertex_First_Jet_Transport_TODO.md)，公式、来源与失败语义见 [API 说明](../api/I03_Hermite_Gals_Vertex_API.md)。开工前核验了 I02 的 8 个交付文件及归档 SHA-256，保留了开工 TODO 和 1,545 个既有文件的哈希。

## 实现

| 文件 | 变化 |
| --- | --- |
| [HermiteGalsVertexTransport.h](../../include/Modules/Interface/HermiteGalsVertexTransport.h) | 显式 master 子集输入、只读带来源 vertex result |
| [HermiteGalsVertexTransport.cpp](../../src/Modules/Interface/HermiteGalsVertexTransport.cpp) | native mesh/source/时间校验、I02 回溯、单一 snapshot 采样及转置拉回、整批失败隔离 |
| [interface_hermite_gals_vertex_contract.cpp](../../tests/interface_hermite_gals_vertex_contract.cpp) | 真实 Q4 源、解析与数值微分、XY/RZ 全顶点局部收敛和负面契约测试 |
| [tests/CMakeLists.txt](../../tests/CMakeLists.txt) | 注册 I03 contract |

计算直接使用 `phi_new=phi_old(B_h)`、`grad_new=J_Bh^T grad_old(B_h)`。phi 和 gradient 来自同一旧 snapshot 的同点查询；J 来自 I02。固定点仍计算非平凡切线；结果只含显式目标的 vertex first jets，不冒充新的完整状态。

旧 primary、edge、相体积、phi_xy、bubble 系数、snapshot 身份及时间保持不变。非法/hanging/缺失 target、错误来源/时间/mesh/domain/axis/尺度、physical boundary hit、离开 active geometry、非有限拉回及中途 callback 失败均经拒绝测试。真实 AMR fixture 只用于验证 hanging node 不得独立更新；一侧 corner-sector fixture 明确拒绝，不代表 AMR/sector 输运完成。

## 数值结果

原始证据：[runs/baselines/I03_20261001](../../runs/baselines/I03_20261001/)。主测试源是严格正的 tensor Q4 scalar，因此真实负相库存为零；这是输运算子验收，不是界面动力学 benchmark。

| 检查 | 覆盖 | 最大绝对误差 |
| --- | --- | --- |
| 解析平移/非对称剪切 | 49 个内部 master vertices，phi | 5.55112e-16 |
| 解析平移/非对称剪切 | 两个 gradient 分量 | 2.85882e-15 |
| 实际复合 scalar 中央差分，epsilon=1e-5 | 49 个内部顶点、两个分量 | 2.35965e-11 |
| 静止和零步 | 全部 81 顶点，与旧 snapshot 采样比较 | 逐位相等 |
| 静止采样与原 primary 比较 | phi / gradient | 0 / 2.77556e-17 |

affine 原预设容限为 `2e-12*(1+|reference|)`，独立中央差分容限为 `2e-8`；均保持不变。完整数值见 `affine.csv`、`scalar_difference.csv` 和 `stationary.csv`。

两个非线性解析速度均覆盖全网格 81 个顶点，包含位置固定但梯度改变的边界点；RZ 轴线保持偶标量的零径向导数。`dt=0.2,0.1,0.05,0.025`，以下是最细一级；完整 8 行见 `local_convergence.csv`。

| 速度 | phi 误差 | gradient 误差 | phi 局部阶 | gradient 局部阶 |
| --- | --- | --- | --- | --- |
| XY logistic | 3.75415e-8 | 8.07436e-7 | 2.993119 | 3.006337 |
| RZ cubic radial + axial logistic | 1.17793e-7 | 4.32528e-6 | 2.973183 | 3.012723 |

最后两级全部满足开工时指定的 `[2.8,3.2]`。这是 RK2 单步局部误差的三阶表现，不能替代多步完整 Interface 的时间精度、守恒或曲率验收。

## 首次失败与定位

首次 CTest 为 4/5 通过，I03 在静止断言失败。初版测试将“无输运漂移”写成与原始 primary first jet 逐位相等。诊断显示节点 3 的 dy 为 `-0.22937744140625002`，primary 为 `-0.22937744140625005`；旧 snapshot 的直接采样已经是前一个值。

Core 的 `QuarticCell::approximateJetAtChart` 使用归一化缓存进行浮点 metric 计算，因此存在 `2.77556e-17` 的舍入差。最终测试要求结果与同一 snapshot 采样逐位相同，并用 `8*machine_epsilon*(1+|primary|)` 检查缓存与 primary 的舍入量。生产公式未改，未绕过旧几何，未放宽解析误差或收敛阶容限。

`ctest_initial_failure*.log/json`、`stationary_diagnostic.log` 和 `TODO_before_implementation.md` 保留首次失败与原始验收；`vertex_refined_oracle.log` 和最终 CTest 日志保存定位后的结果。

## 构建、回归与范围

| 检查 | 结果 |
| --- | --- |
| Release：主程序、I01/I02/I03、旧 Interface contract、architecture | PASS；初次 11.696 s，最终增量 5.956 s；无新增 compiler warning |
| CTest：help / architecture / I01 / I02 / I03 | 5/5 PASS，0.01 / 0.75 / 0.01 / 0.01 / 0.10 s |
| 原 `quantum_interface_state_contract --parametric-state-control` | PASS，退出 0，0.383 s；精确计时见 legacy_result.json |
| 既有代码树/输入/工具 | 1,544 文件哈希不变；仅既有 tests/CMakeLists.txt 改动，另增 3 个源码/测试文件 |
| S08 冻结实现 / 方法 PDF | 69/69 与 3/3 哈希一致 |

预算为 2 个 production、2 个测试文件。Core、I01/I02、旧数值实现、Apps、reader、默认配置及 case 输入均未修改；SuperLU 既有 `/LTCG` 自动重链接信息保留。Git 元数据不可用，没有伪造 commit。

`verification.json` 记录命令、工具链、线程及数值摘要，`ctest_details.log` 保留测试原文；`runtime/`、`runtime_hashes.json` 保存执行文件身份。交付副本及 SHA-256 见 `delivery/`、`delivery_files.json`、`evidence_inventory.csv`、`evidence.zip` 和 `delivery_manifest.json`。

## 未完成

I04 edge first jet、I05 closure、I06 AMR、I07–I10 守恒闭合、I11 同源几何发布及 I12 完整 benchmark 均未开始。`conservativeHermiteGals` 未接入 reader/生产 App，也未接入 Hydro 新主线。

物理边界的 boundary-history tangent 与 corner-sector 一侧输运仍需专门实现；本阶段遇到这些来源明确失败。速度 provider 的同源 value/derivative、不可变性与生命周期仍是调用前置条件。本轮未构造 Hydro 插值器，也未进行大网格性能验收。

S01 已知的 `hydro_moment_xy` 首步失败继续保留；本次短回归不能将其改记为通过。后续按路线图执行 I04。
