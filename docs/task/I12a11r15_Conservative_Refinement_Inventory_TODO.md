# I12a11r15 / 精确细化的守恒库存转移

## 目标与修改范围

只解决一个事实：给 r14 的同一精确细化标量分配子叶独立库存，使每个父叶的原始 `M_P` 保留为守恒目标，子叶严格有界，且不把库存重置成重新测量的几何体积。

允许新增 2 个 production 文件：`include/Modules/Interface/HermiteGalsRefinementInventory.h`、`src/Modules/Interface/HermiteGalsRefinementInventory.cpp`；2 个 test-related 文件：`tests/interface_refinement_inventory_contract.cpp`、`tests/CMakeLists.txt`。另含本 TODO/报告、API/README 和 `runs/baselines/I12a11r15_20261001/`。Interface 的已有源文件、Core、I06/I10/I11、flux/B/J、全部冻结文件与数值工作预算禁止修改。

## 数学约束

输入旧 mesh、同一只读 `HermiteGalsState`、目标 mesh 与原 I07 测量配置。要求完整 native Q4 支持和全域库存；同一时间层、XY/RZ 度量不变。r14 工厂提供唯一精确目标 scalar，复用现有 `InterfaceCellMeasure` 求积与纯相证书；按真实叶边界建立唯一父子归属，不反演显示坐标。

旧库存必须严格位于容量界限内，且与同一旧 scalar 的 canonical 几何液量相差不超过原 `64*epsilon` 相对容量预算。cut cell 使用现有局部 Green 液量，纯相使用原严格相位证书的 0/容量。目标几何与父库存不相容时明确拒绝，不改变旧来源。

未变叶直接复制原 `M_P`。细化叶以子几何液量为分配参考，但原 `M_P` 始终是目标。对子叶定义严格容量界限与原 64 epsilon 几何偏差预算的交集；仅在同一个父叶的这些允许增减量之间按剩余预算比例分配总差 `M_P - sum(V_C)`。分配前检查可行性，分配后检查实际存储值、父子总量及几何偏差。不得事后 clipping，不得跨父叶借库存，不得忽略未分配差额；无法表示/不相容则抛出带父叶/子叶上下文的异常。

父子求和误差只允许原共享浮点 roundoff budget，必须报告原目标、几何和、实际子库存和与残差。T1 另要求常规夹具父子缺陷在 4 epsilon 父容量内。保持现有数值容限、迭代/求积工作上限，不为失败场景扩大预算。

## 架构及实现要求

新增只能由转移函数生成的只读结果，拥有旧来源身份/时间、旧新 mesh identity、精确目标 Snapshot、子库存、几何测量、父子映射与审计报告。禁止返回可写别名；copy/move 与来源/网格销毁后的读访问明确受检。结果不是 `HermiteGalsState` 或 `HermiteGalsPublication`，不通过 `bindGeometry/finish` 伪造 accepted regrid。

复用原 r14 限制、I07 度量和同源几何积分；不新增多项式重构、mixed closure、bubble、相区间/Green stencil 或平移解析参考。新源文件由现有 Interface target glob 正常编译。后续 I11 发布和局部重试另立任务。

## 验收与失败规则

T0：Release 无 warning；新库存契约与相邻 r14、I07/I10/I11、时空库存/reflux 契约通过；原输出不变。81 冻结/参考文件、3 PDF 及 r14 archive 不变。采用 r14 已核验的测试墙钟设置，数值工作预算不变。

T1：XY/RZ、正反相、uniform/AMR/多层、未变叶和纯相；实际父子总量与容量界限；原 M 含可接受正负 roundoff 偏差时仍保留原目标；重复调用一致、只读及来源生命周期。独立检查测量来自结果的同一 scalar，原 M/identity/time/primary 不变。覆盖错误 mesh/domain、coarsening、部分覆盖、过界库存、超原几何预算、无可行分配、moved-from 等明确拒绝。

PASS 只代表库存转移算子。任何数值/表示/来源断言不满足则 FAIL 并保留证据；只有具体缺失依赖才 BLOCKED。父任务 I12 仍须真实推进和全部原矩阵通过，不进入 Hydro。

## 交付

源码、算子场景/拒绝证据、实际测试日志与耗时、父子守恒及同源测量审计、冻结与 runtime 哈希、API、报告和可复核 archive。
