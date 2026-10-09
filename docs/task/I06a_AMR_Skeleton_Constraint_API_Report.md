# I06a / AMR skeleton 精确约束接口交付报告

**结论：PASS（仅 I06a）。** 执行日期：2026-10-01。I06 整体、Gate I 和完整电喷雾目标均未完成，下一项 I06b。

任务边界、拆分理由和预设阈值见 [TODO](I06a_AMR_Skeleton_Constraint_API_TODO.md)，接口语义见 [API 说明](../api/I06a_AMR_Skeleton_Constraint_API.md)。完整 I06 需要至少 6 个 production 文件，依架构规范 §8 先拆接口与数值接线；本轮没有修改路线图顺序。

前序接口调研明确了可复用的 exact constraint owner，归为 progress；本轮重新核验 I05 的 10 个交付文件和归档，全部哈希一致，保存 1,556 文件的开工清单及实施前 TODO，再执行代码修改。

## 修改范围

| 文件 | 变化 |
| --- | --- |
| [ActiveQuarticHermiteSkeleton.h](../../include/Core/Geometry/ActiveQuarticHermiteSkeleton.h) | 中立 skeleton input、只读 owner、native/master 元数据与 node/edge trace 查询 |
| [ActiveQuarticHermiteSkeleton.cpp](../../src/Core/Geometry/ActiveQuarticHermiteSkeleton.cpp) | 复用现有约束引擎，拥有输入/mesh identity，包装失败上下文 |
| [QuarticLeafData](../../include/Core/Geometry/Detail/LevelSetHermiteQuarticLeafData.h) | 增加 skeleton-only 构造出口和只读查询；共享原有 exact trace/linear combination |
| [CMake](../../CMakeLists.txt) | 新 cpp 加入已有 strict FP / no-LTO geometry object target |
| [skeleton contract](../../tests/interface_hermite_gals_skeleton_contract.cpp) / [注册](../../tests/CMakeLists.txt) | 公有 API、解析/任意数据/ownership/拒绝与旧 geometry 对照 |

预算：4 production + 2 test。没有修改 I01–I05、mesh、Modules、Apps、字典、case、三份方法 PDF 或 S08 冻结专有实现。

输入中只有 regular vertex 的四个 Hermite 参数和 maximal side 的两个 midpoint 参数进入 Core 依赖图；hanging/fine data 完全由原有 exact master constraints 生成。Core 中参与插值的 mixed 参数仍是 Interface 派生量，未新增独立输运状态。构造不包含任何 center datum 或 cell geometry。

## Q4 与 coarse/fine 验证

证据目录：[runs/baselines/I06a_20261001](../../runs/baselines/I06a_20261001/)。最终测试原文为 `skeleton_final.log`，表格为 `q4_reproduction.csv`、`c1_constraints.csv`。

Q4 在平移、非二进制坐标、各向异性 XY/RZ 域上验证：uniform、单格细化、相邻区域细化、0–3 层 2:1 嵌套细化。单格细化测试 XY 全部 25 个 monomial + 1 个线性组合，RZ 全部 15 个 even monomial + 1 个组合；其余拓扑测试混合 Q4。

| 检查 | 覆盖 | 最大绝对误差 |
| --- | --- | --- |
| 解析 Q4 node / midpoint / trace jet | 8 组 XY/RZ × 拓扑，20,242 次 node/trace 查询，另核对所有 midpoint | 1.42109e-13 |
| shifted Q4 coarse/fine trace | 水平/竖直、三级细化 | 1.77636e-15 |
| axis parity | 978 次轴线 node/trace 查询，含细化轴线 | dx/dxy 严格零 |
| 任意 master DOF coarse/fine | 6 组 XY/RZ × 3 拓扑，每条边 9 点 | 0 |
| 独立 power-basis quartic 参考 | 同上，含端点为 hanging 的 perpendicular master edges | 0 |
| 原有 native Q4 snapshot 的 cell-side jet | 每个 cell 的每条边 3 点，任意 cell-center 输入 | 3.41061e-13 |

全部使用预设 `3e-11*(1+|reference|)` 容限，没有放宽。任意数据测试用固定整数随机序列和 dyadic 网格；零误差只针对该组实际测量，不声称所有浮点网格位级相同。shifted 网格的非零 defect 来自最终 double 查询，内部约束仍使用同一 exact homogeneous graph。

三级 XY：52 leaves、24 hanging vertices、48 constrained edges、396 个局部独立插值参数；三级 RZ：52 leaves、18 hanging vertices、36 constrained edges、432 参数。计数逐网格验证 `4*N_master_vertex+2*N_master_edge`，没有 hanging 或 cell center 进入独立计数。

hanging/fine 输入置 NaN、完全省略 fine midpoint、不修改 owner 的并发查询（8 个任务各 20 次）、修改/析构源 input 与 mesh、owner copy/move、uniform 与 coarse-support 完整的 AMR active subset 均通过。缺失 coarse support、非法 identity/mask/数值、axis parity 与非法查询均明确拒绝。

## 回归与工程验证

| 检查 | 结果 |
| --- | --- |
| Release 主程序、I01–I05、新 contract、旧 Interface、architecture 构建 | PASS；首次全目标 27.402 s，无新增 compiler warning |
| CTest help、旧 Interface state、I01–I05、新 skeleton | 8/8 PASS，7.94 s |
| architecture | 1/1 PASS，0.76 s |
| 增补 AMR active subset 后最终 skeleton CTest | PASS，0.79 s |
| 旧 `--parametric-state-control` / Core GALS control | 退出均为 0；Core GALS 完整输出逐字不变 |
| 原有 native geometry/AMR contract 前后 | 全部检查通过；去除 CTest 线程限制引出的两行 OpenMP 提示后，完整输出逐字相同 |
| I03/I04 输出 | 23 行数值记录与 I05 归档逐字相同 |
| 实际编译选项 | 新 cpp `/fp:strict /GL-` 已核对并保存 |
| 源码/输入/工具边界 | 1,553 既有文件不变；仅根 CMake、Core detail header、tests CMake 三个既有文件改变，另新增三个文件 |
| 冻结文件 / 方法 PDF | 69/69、3/3 哈希一致 |

首次测试中，多层 refine box 未命中 root cell center；第二次构造出非 2:1 面，Core 正确拒绝。最终改为逐层内缩的 refinement boxes，并断言实际最高层为 3。另修正测试期望中的 `nonfinite` 为上游真实消息 `non-finite`。这些都是测试 fixture/消息检查修正，生产数学、异常规则与阈值不变；三次失败日志保留。开始验收前还对 shifted endpoint 查询使用真实端点，避免测试算术把端点舍入到边外。

`source_audit.json`、`compile_contract.json`、各结果 JSON、runtime 哈希、交付副本及归档记录验证来源。Git 元数据不可用；没有伪造 commit。

## 未完成与下一步

下一项 I06b：在同一 primary revision 上实现 AMR 非均匀 compact mixed closure，接入 I03/I04 只推进 master 的结果，并完成 I06 最终验收。尤其要解决局部细化短 skeleton 的合法支持点及 closure/constraint 先后关系；不能用缺少 Q4 支持的低阶 stencil 或任意 mixed seed 替代。

本轮提供的是已给定 master mixed 参数后的精确约束接口，不证明 AMR mixed 收敛、完整时间推进、曲率、phase-volume conservation 或低黏稳定性。I07–I12 和后续 Hydro/Electricity/EHD gate 均保留。S01 的 `hydro_moment_xy` 首步已知失败未被此次回归重新分类。
