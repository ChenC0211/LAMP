# I04 / master edge-midpoint first-jet 输运交付报告

**结论：PASS。** 执行日期：2026-10-01。本轮只做 I04，下一任务为 I05；Gate I 及低黏/无黏电喷雾完整目标仍未完成。

依据与预设验收见 [TODO](I04_Edge_First_Jet_Transport_TODO.md)，公式、目标身份和 normal 约定见 [API 说明](../api/I04_Hermite_Gals_Edge_API.md)。前一轮归类为 progress；开工逐个核验 I03 的 8 个交付文件和归档哈希，并保留 1,548 个既有文件的 source baseline。

## 交付

| 文件 | 内容 |
| --- | --- |
| [HermiteGalsEdgeTransport.h](../../include/Modules/Interface/HermiteGalsEdgeTransport.h) | 显式 native edge 选择、只读带来源 result |
| [HermiteGalsEdgeTransport.cpp](../../src/Modules/Interface/HermiteGalsEdgeTransport.cpp) | native master side 验证、midpoint、完整 gradient 拉回后的 normal 投影、批量失败隔离 |
| [HermiteGalsFirstJet.h](../../src/Modules/Interface/HermiteGalsFirstJet.h) | vertex/edge 共用的私有来源校验和单点 first-jet 计算 |
| [HermiteGalsVertexTransport.cpp](../../src/Modules/Interface/HermiteGalsVertexTransport.cpp) | 改用共享代码，公有 I03 API 和原测试不变 |
| [edge contract](../../tests/interface_hermite_gals_edge_contract.cpp) / [测试注册](../../tests/CMakeLists.txt) | 独立算子、native 身份与拒绝路径验收 |

预算为 4 个 production 文件、2 个测试文件。Core、I01/I02、冻结实现、Hydro/Electricity、Apps、reader、字典和 case 输入未修改。

edge value 与完整梯度来自同一旧 snapshot 的同一点，使用 I02 的 `J_Bh^T` 拉回后才提取 normal；水平边取 +y，竖直边取 +x。不会把旧 edge dn 当作完整 gradient。source 的 vertex/edge/M/closure/identity/time/snapshot 均保持不变，结果仅代表显式选定 edge 的阶段更新。

## 算子数值

证据目录：[runs/baselines/I04_20261001](../../runs/baselines/I04_20261001/)。精度 fixture 是真实 native Q4 scalar，严格正，因此负相库存为零。这是输运算子验证，不是界面动力学 benchmark。

| 检查 | 覆盖 | 最大绝对误差 |
| --- | --- | --- |
| 平移、x/y 两方向 affine shear 的 phi | 56 条水平 + 56 条竖直内部边 | 5.55112e-16 |
| 同上，normal derivative | 同上 | 8.93730e-15 |
| 对实际复合 scalar 沿 normal 中央差分，epsilon=1e-5 | 全部 112 条内部边 | 2.09633e-11 |
| 静止/零步与旧 snapshot 采样比较 | 全部 144 条边 | 逐位相等 |
| 静止 snapshot 与原 edge DOF 比较 | phi / dn | 0 / 5.55112e-17 |

所有原定容限保持不变：affine 为 `2e-12*(1+|reference|)`，中央差分为 `2e-8`，静止缓存舍入为 `8*epsilon*(1+|DOF|)`。固定 midpoint 的非单位 J 测试覆盖两个边方向；Midpoint 时间层与同物理时刻 End 层的计算一致。

XY logistic 和 RZ 奇对称 cubic radial 解析速度覆盖全部 144 条边，包含物理边界及 8 条 RZ axis 边。`dt=0.2,0.1,0.05,0.025` 的完整 8 行保存在 `local_convergence.csv`，最细一级如下：

| 速度 | phi 误差 | dn 误差 | phi 局部阶 | dn 局部阶 |
| --- | --- | --- | --- | --- |
| XY logistic | 3.62402e-8 | 7.90967e-7 | 2.993087 | 3.006694 |
| RZ cubic radial + axial logistic | 1.34483e-7 | 4.22623e-6 | 2.994178 | 3.013038 |

两例最后两级均满足 `[2.8,3.2]`。这是 RK2 单步局部误差阶，不能当作完整 Interface 的多步三阶精度。axis 径向 dn 保持严格零。

## native 身份与隔离

测试真实 2:1 mesh 上的 44 条 master sides，其中 4 条以 hanging vertex 为端点的垂直方向细边仍被正确接受；8 条位于 coarse side 上的 fine subedges 被拒绝，即使测试的 primary 列表故意将其列为独立 edge。此检查证明 I04 不独立推进受约束子边，不代表 I06 AMR closure 已完成。

其它拒绝检查覆盖空/重复/逆序/反向/未知 endpoints、非 native 长边、缺失 source edge、错误 geometry identity/time/domain/axis/尺度、mesh 变化、boundary-history tangent 缺失、离开旧活动几何、corner-sector 和 callback 中途失败。重复目标在任何 velocity query 之前拒绝；中途失败发生在第一条边已完成之后，source 仍保持不变。共享拉回的非有限结果拒绝也由未修改的 I03 overflow contract 回归覆盖。

不可表示 midpoint 的网格在 Core 的 `DyadicGridCoordinates` 构造时已报错，测试核对该明确上游拒绝。I04 内部的 midpoint 防御分支没有通过伪造无效 mesh/source 来强行触发，不声称该分支被直接覆盖。

## 构建与回归

| 检查 | 结果 |
| --- | --- |
| Release：QUANTUM、I01–I04、旧 Interface contract、architecture | PASS，全目标重建 9.205 s；无新增 compiler warning |
| CTest：help / architecture / I01 / I02 / I03 / I04 | 6/6 PASS；0.01 / 0.73 / 0.01 / 0.01 / 0.05 / 0.12 s |
| 原 `quantum_interface_state_contract --parametric-state-control` | PASS，退出 0，0.383 s |
| I03 原始 affine/FD/stationary/order 数值 | 11 行 17 位精度输出与 I03 归档逐字一致 |
| 既有源码/输入/工具哈希 | 1,546 不变；仅 I03 cpp、tests/CMakeLists.txt 改动，另增 4 个源码/测试文件 |
| S08 冻结文件 / 方法 PDF | 69/69 与 3/3 哈希一致 |

首次编译在新测试的多变量 `auto` 声明处失败，修正为独立声明。首次 CTest 为 5/6：所有 I04 数值和 native side 测试通过，随后测试试图构造的极窄网格被 Core 提前拒绝。最终测试明确核对这个实际拒绝位置，没有改变生产计算、Core 或数值容限。日志保存在 `build_initial_failure*` 与 `ctest_initial_failure*`，开工 TODO 副本也保留。

构建/测试命令、工具链、线程及摘要见 `verification.json`；原始测试输出见 `ctest_details.log`，I03 对照见 `I03_numeric_comparison.json`。`runtime/`、`runtime_hashes.json` 保留实际执行文件；`delivery/`、`delivery_files.json`、`evidence_inventory.csv`、`evidence.zip`、`delivery_manifest.json` 保存交付副本与 SHA-256。Git 元数据不可用，没有伪造 commit。

## 未完成

I05 mixed closure、I06 AMR constraint/transfer、I07–I10 守恒闭合、I11 同源几何发布和 I12 完整 benchmark 均未在本轮开始。V2 仍未进入 reader/生产 App 或 Hydro 新主线。

物理边界历史 tangent、corner-sector 一侧数据以及完整 active-band 管理尚未实现；当前按明确契约拒绝。速度插值 value/derivative 的同源性和不可变生命周期由 provider 保证。未声称此测试证明曲率、守恒、破碎、低黏稳定性或大网格性能。

S01 已知 `hydro_moment_xy` 首步失败继续保留，短回归没有将它重新分类。下一项按路线图执行 I05。
