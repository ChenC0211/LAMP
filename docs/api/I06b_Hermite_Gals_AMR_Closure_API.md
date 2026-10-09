# I06b / AMR mixed closure 与 skeleton API

状态：I06b 与父任务 I06 已验收，2026-10-01。接口：[HermiteGalsAmrClosure.h](../../include/Modules/Interface/HermiteGalsAmrClosure.h)，实现：[HermiteGalsAmrClosure.cpp](../../src/Modules/Interface/HermiteGalsAmrClosure.cpp)。仍未完成 I07–I12，不是 production Interface 发布入口。

## 输入与结果

```cpp
HermiteGalsAmrClosureResult closeHermiteGalsMixedAmr(
    const QuadtreeMesh& mesh, const HermiteGalsCandidate& candidate,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes = {});
```

输入必须是 I01 的未绑定 geometry 的 primary revision，mesh 身份一致。`geometryCells` 指定非空 native active support；primary.vertices 恰好包含其中所有 regular masters，primary.edges 恰好包含所有 maximal active master sides。hanging vertex 不得进入 primary，fine subedge 不能冒充 master。

支持完整 native 2:1 mesh，以及 coarse 约束与五点 stencil 支持完整的 active subset。拓扑由 I06a Core owner 校验；缺 coarse source、缺 master、额外未使用数据、非法 side、非 2:1/native layout 均失败。不为缺失支持补低阶公式、跨 inactive gap 取点或另建 scalar。

结果拥有：

- `vertexMixedDerivatives()`：物理 phi_xy，与 candidate.primary().vertices 同序。
- `skeleton()`：完成 mixed closure 后的只读 exact Core skeleton；hanging/fine 数据仅由这个 master trace 限制生成。
- `primaryIdentity()` / `timeLayer()`：来源 candidate 的 revision/time；之后 replacePrimary 不会重新绑定旧结果。
- `report()`：master/hanging vertex、master/constrained edge、连通直线及 compatibility matrix nonzero 数量；两个残差和 directional disagreement。

不修改 candidate、first jet、phase volume、mesh 或既有 snapshot，不生成 center/bubble、cut geometry 或新的 geometry snapshot。I11 publisher 仍须核验结果 identity 并处理 I10 的 bubble。结果的 Core skeleton 是代数边界算子，不能用 rounded hanging 查询值创建第二个独立 geometry owner。

## 非均匀五点支持

水平与竖直 master sides 分别按物理 line coordinate 和正向端点排序。只有前一段终点与后一段起点为同一 native vertex 才连成同一个 component；gap 两侧不共享 stencil。

每个 component 形成 `endpoint, master midpoint, endpoint, ...` 的交错序列：

- x 方向取节点 phi_y 与水平 master.dn；
- y 方向取节点 phi_x 与竖直 master.dn；
- master endpoint 直接读本次 first jet；perpendicular hanging endpoint 读取同一个 Core master trace 的切向导数；
- coarse side 内部的 collinear hanging 与 fine midpoint 不增加样本或未知量。

每个 regular target 使用连续五个样本；起点索引 `clamp(target-2,0,count-5)`，物理边界/短组件端部使用同一四次插值的单侧求导。每个 component 至少五点；即使其它 disconnected component 有足够点，也不能借用。

求导复用 I05 的 [CompactQuarticDerivative](../../include/Core/Geometry/CompactQuarticDerivative.h)：真实坐标、平移缩放、减去目标值、degree-4 interpolant 求导。两个方向固定平均，disagreement 不选择权重。普通 uniform 输入与 I05 的结果在本次测试中位级一致。

## 约束依赖与闭合顺序

Core 新增 `nodeSensitivity(vertex,component)`，返回 node jet 对源 master vertex/edge 参数的浮点链式系数。它复用已有 cell-sensitivity 的 reverse dependency traversal，未复制 restriction 或 transpose。component 为 phi/dx/dy/dxy = 0..3，edge index 对应 Core input.edgeMidpoints 的顺序。

记 master mixed 向量为 q，约束后的 first-jet 样本为 `g(q)=a+Cq`，五点求导并固定平均为 L，则原闭合关系是：

\[
q=L(a+Cq),\qquad (I-LC)q=La.
\]

实现只在装配 affine 常量项时使用 q=0；它不是闭合结果、不对外发布。所有非零 mixed 依赖均由同一 Core DAG 提取，不按大小丢弃。暂存稀疏系统通过 Eigen SparseLU 求解，随后使用最终 q 构造唯一返回的 exact skeleton，再直接代回原五点公式核查一致性。

这个系统的未知量只是本来就需要恢复的派生 mixed 向量，不增加独立 primary、时间历史或持久化全局矩阵，也不修改 phi/gradient 做 KKT correction。当前测试的合法 2:1 组件中，所选 hanging **切向**样本与 mixed 无关，系统实际退化为单位阵。测试同时测到完整 hanging jet 的 40 项非零 mixed 依赖；不能据此宣称闭合矩阵的非对角情形已被覆盖。

## RZ 与残差

AxisymmetricVolume 的非负半径由 I01 校验；域含 r=0 时 Core 要求 axis vertex.dx 和 axis edge.dn 严格零。x 向反射导数及最终 axis mixed 严格取零，y 向 axis 数据也为零。非零 axis 输入被拒绝，不作修正。

`algebraicResidual` 为：

\[
\max_i |(Aq-b)_i|/\max(1,|b_i|,|q_i|).
\]

`compatibilityResidual` 为最终 skeleton 直接求得的 `q-(Dx+Dy)/2`，逐节点除以 `max(1,|q|,|Dx|,|Dy|)`。两者分别校验 `<=1e-12`；前者不能替代后者。矩阵、右端、解、残差以及 final query 都检查有限性；因子分解/求解失败明确抛错。

`maxRelativeDirectionalDisagreement` 使用 I05 相同定义 `|Dx-Dy|/max(1,|Dx|,|Dy|)`，只作诊断。三个量先在 kernel 中计算和校验，再复制到只读 Report；Report 不驱动算法分支。

## 精度与边界

I12a10r2 增加显式外边界 symmetry，沿原 RZ 反射约束路线。调用方将与 I07/I10/I11 相同的声明传给 I06；默认空列表保持原行为。每个 plane 必须是有限的真实外域截面、有 active master trace，不能重复或重复声明自动 RZ 轴。对应 primary vertex 法向导数和 master edge.dn 必须严格为零，极小非零也拒绝，不调整 primary。

法向反射方向的紧凑导数取精确零，对应派生 q 的兼容系统行施加 q=0；其余行仍原五点固定平均。求解后只在这些已装配的约束行写入精确零，再构造同一 Core exact hanging skeleton。实际边界 node 的 normal/mixed 再次检查严格零，最终五点代回残差仍按原 1e-12 检查。没有按数值大小推断 symmetry、没有改变未声明边界、没有对 first jets 做投影。

Report 新增 `symmetryConstrainedVertices`，统计由显式 plane 新增的不同 master 约束节点；已由自动 RZ 轴约束的节点不重复计数。该计数不驱动算法。此功能的完整推进证据单独见 [I12a10r2 TODO](../task/I12a10r2_Declared_Symmetry_Closure_TODO.md)。

Q4 再现、真实 I03/I04 master-only 输运后的独立 Vandermonde 代回、C1、active subsets 与隔离检查见 [报告](../task/I06b_AMR_Mixed_Closure_Report.md)。AMR 收敛预设至少约三阶：quartic tangential trace 的 O(h^4) 求导误差再沿 skeleton 求导时可产生 O(h^3) 项。因此保留 uniform 四阶，同时分别报告 AMR master/hanging 误差，不能直接承诺所有 AMR 点四阶。

本轮两个 MMS 的最细 master 阶分别为 3.2966、3.5251，满足预设 `[2.7,4.3]`；这不是完整界面、曲率或输运时间阶的验收。构造测试样本时的库存不作为 I07 几何体积初始化证据。

失败消息包含 `[Interface/HermiteGalsAmrClosure primary=... vertex=...]` 和阶段；整批成功才返回。corner-sector/contact、physical inflow history、动态 regrid/coarsening 的守恒事务、I07–I12 与 production App 接线仍在后续范围。
