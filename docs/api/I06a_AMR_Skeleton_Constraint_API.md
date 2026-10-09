# I06a / AMR master–hanging skeleton 约束 API

状态：I06a 已验收，2026-10-01。后续 [I06b](I06b_Hermite_Gals_AMR_Closure_API.md) 已完成非均匀 mixed closure、nodeSensitivity 扩展与 first-jet transport 接线，父任务 I06 已验收。本文记录 I06a 的基础约束接口。接口：[ActiveQuarticHermiteSkeleton.h](../../include/Core/Geometry/ActiveQuarticHermiteSkeleton.h)，实现：[ActiveQuarticHermiteSkeleton.cpp](../../src/Core/Geometry/ActiveQuarticHermiteSkeleton.cpp)。

## 职责与输入

```cpp
Quantum::Geometry::ActiveQuarticHermiteSkeleton skeleton(mesh, input, radialAxisSymmetry);
```

这是中立 Core 的不可变代数约束 owner，没有 Interface 时间、相体积、速度、bubble、cell centers 或 geometry snapshot。Interface 后续调用它时仍须绑定 I01 primary identity/time；Core 的 `validateMesh(mesh)` 只验证 native geometry/connectivity 身份。

`ActiveQuarticHermiteSkeletonInput` 包含：

- `activeCells`：完整 native cell-index 长度的 0/1 mask。
- `nodes`：完整 native vertex-index 长度的物理 `(phi,phi_x,phi_y,phi_xy)`。只消费 active regular vertices；hanging 和 unused payload 可为 NaN。
- `edgeMidpoints`：实际 active leaf side，端点索引按正 x（水平）或正 y（竖直）排列；值为 exact parameter midpoint 的 phi 和正坐标方向 normal derivative。水平 normal 为 +y，竖直为 +x。

每个 maximal active side 必须提供唯一 midpoint 输入。fine subedge 可不提供；提供时其数值也被忽略，但 side identity 仍须合法且唯一。active hanging vertex 的实际 coarse owner 必须在支持集中，缺失时明确失败，不将其提升为独立节点。

`phi_xy` 在此是调用方已经恢复的派生量。Core 将它视为局部 Hermite 参数参与约束，**不把它变成独立时间推进状态**。此基础接口自身不执行 AMR mixed closure；不能用任意零 mixed 值代替 I06b 的实际闭合解，也不能单凭此接口宣称 V2 已能进行 AMR 时间推进。

`radialAxisSymmetry=true` 要求域下界 r=0；axis master 的 dx/dxy 和 axis 竖直 master edge.dn 必须严格零。违反时拒绝，绝不自动修正。无轴的非负半径域可作为普通 Cartesian skeleton，由 Interface 的 measure 语义另行约束。

## 所有权与查询

构造深拷贝输入并执行现有 `QuarticLeafData` 约束；它在生成全部 node/edge 约束后停止，未注册任何 cell-center datum，也未构造 cell polynomial、sector、cut geometry 或 snapshot。mesh 和 input 析构后，已有 owner 的查询仍有效；owner 可以廉价复制且支持并发只读查询。

| 方法 | 结果 |
| --- | --- |
| `vertices()` | native index 递增的 used vertices；`masterSide={-1,-1}` 表示 independent regular vertex，否则为实际 coarse/master side |
| `edges()` | lexicographic side keys；`parentSide={-1,-1}` 表示 independent midpoint，否则为 coarse side |
| `independentDofCount()` | `4*N_regular_vertex+2*N_master_edge`，无 hanging/fine/cell-center DOF |
| `nodeJet(v)` | physical `(phi,dx,dy,dxy)` |
| `edgeMidpoint(edge)` | `(phi,+normal derivative)`，使用实际端点的精确 parameter midpoint |
| `edgeJet(edge,s)` | physical `(phi,dx,dy,dxy)`，s 为该边上的物理 x 或 y 坐标，含端点 |
| `validateMesh(mesh)` | 校验调用方当前 mesh 仍与构造时一致；不是 regrid 操作 |

DOF count 是局部插值参数数量，不能当作 Interface 独立输运未知量数。node/edge 查询返回最终转为 double 的数值；**不能将这些舍入结果当成新的独立 hanging data 重新导入另一份 production geometry**。

## 约束来源与精确性

构造继续使用原有 `resolveVertex/resolveMidpoint/traceJet`：

\[
\mathcal D_{hang}=C_{AMR}\mathcal D_{master}.
\]

同一 coarse edge 的 value trace 和 normal-derivative trace 各由五个 quartic Hermite 条件决定。fine node/edge data 在 exact homogeneous dependency graph 中生成，fine trace 从这份精确参数图求值；直到公开查询返回时才转为 double。

本任务把原有 trace coefficient 访问和 homogeneous linear combination 抽为共享内部函数。旧 geometry 构造和新只读查询调用同一实现；依赖登记顺序、basis、乘法/累加顺序不变。查询不向 dependency graph 添加数据，不修改 owner，也没有复制另一套 polynomial stencil。

新 cpp 单独加入已有 `quantum_hermite_geometry` strict FP / no-LTO object target。MSVC 的实际编译命令已验证含 `/fp:strict /GL-`；不是仅在文档中声明编译语义。

## 失败与验证边界

非法 active mask、缺 coarse support、非 2:1 native layout、缺 master midpoint、重复/反向/非 side key、非有限独立输入、不合格 axis 数据，以及域外/非有限查询均失败。错误包装包含 `[ActiveQuarticHermiteSkeleton]`、构造/查询阶段，node/edge 查询还包含对应索引。

无默认未绑定实例；移动后的源对象查询失败。候选构造失败不发布部分 owner，输入、mesh 及已有 owner 保持不变。

[报告](../task/I06a_AMR_Skeleton_Constraint_API_Report.md) 给出 shifted XY/RZ Q4、0–3 层拓扑、任意 master 数据 C1、旧 snapshot 边界对照、支持完整的 AMR active subset 与拒绝路径证据。这里的 C1 是 value/gradient trace 连续及其 tangential derivative 约束，不声称全 Hessian 为 C2 连续。

本接口没有新增 regrid/全 cell refinement/coarsening、mixed closure、守恒通量、bubble 或 geometry publication。mixed closure 见 I06b；I07–I12 仍需按路线图完成。
