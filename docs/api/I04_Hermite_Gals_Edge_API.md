# I04 / Hermite–GALS master edge-midpoint API

状态：I04 已验收，2026-10-01。依赖 [I01 state](I01_Interface_V2_State_API.md) 与 [I02 tracer](I02_Hermite_Gals_Characteristic_API.md)，与 [I03 vertex](I03_Hermite_Gals_Vertex_API.md) 共用一阶输运执行链。完整闭合、几何发布和 Gate I 尚未完成。

公有接口：[HermiteGalsEdgeTransport.h](../../include/Modules/Interface/HermiteGalsEdgeTransport.h)。实现：[HermiteGalsEdgeTransport.cpp](../../src/Modules/Interface/HermiteGalsEdgeTransport.cpp)。私有单点公式与来源校验：[HermiteGalsFirstJet.h](../src/Modules/Interface/HermiteGalsFirstJet.h)。

## 公式与 normal

```cpp
HermiteGalsEdgeTransportResult transportHermiteGalsEdges(
    const QuadtreeMesh& mesh,
    const HermiteGalsState& source,
    const std::vector<std::array<int,2>>& masterEdges,
    const LevelSetBoundaryTrace::Domain& domain,
    const HermiteGalsVelocityView& previousVelocity,
    const HermiteGalsVelocityView& arrivalVelocity);
```

对边的 midpoint `x_m`，按冻结方法 §4.2 式 (22)–(24)：

\[
\phi_m^{new}=\phi^{old}(B_h(x_m)),\qquad
g_m^{new}=J_{B_h}(x_m)^T\nabla\phi^{old}(B_h(x_m)),\qquad
(\partial_n\phi)_m^{new}=n_e\cdot g_m^{new}.
\]

phi 和完整旧 gradient 来自 source 持有的同一个 native Q4 snapshot 的一次同点查询；B/J 来自 I02 的同一离散 midpoint map。完整 gradient 先拉回，再取 normal 分量。不能只输运旧 dn，因为剪切会把旧切向导数贡献到新 dn。

| 边 | 端点方向 | 固定 normal | 保存的分量 |
| --- | --- | --- | --- |
| 水平 | +x | (0,1) | g_new.y |
| 竖直 | +y | (1,0) | g_new.x |

这些 normal 是 skeleton 的坐标方向，既不是 cell outward normal，也不是液体界面 normal。RZ 中 x 表示 r，轴线竖直边的 dn 为径向导数。固定 midpoint 也必须使用真实 J；位置不动时 dn 仍可能变化。

## 目标身份与 native topology

masterEdges 由原生 vertex indices 组成，不接受 vertex label、face ID 或调用方传入的坐标。目标非空、严格按 endpoint pair 的字典序排列、无重复；每条边的两个端点必须沿正坐标方向，且 pair 已存在于 source.primary().edges。

I01 只保证 edge 列表的结构合法，因此本入口还独立核验 mesh 的真实 leaf sides。所有 leaf sides 按方向、所在直线、起点和终点分组排序；严格包含于更长 native side 的 fine subedge 不独立推进。此检查遍历全部 native leaves，不因 active band 的选择而漏掉 coarse 约束。

以 hanging vertex 为端点不自动意味着该 edge 被约束：垂直于 coarse side 的细边仍可以有独立 midpoint DOF。测试同时覆盖合法的这种边和必须拒绝的 coarse-side 子边。这里只识别目标身份；AMR restriction、hanging trace 生成、coarsening/regrid 和完整 C1 closure 属于 I06。

midpoint 由实际端点作 `a+(b-a)/2`，固定坐标保持原值，变化坐标必须严格位于两个端点之间。内核不把无法表示的 midpoint 吸附到端点。当前 Core mesh constructor 会提前拒绝测试中的 binary64 坐标坍缩网格；I04 的检查仍作为入口防御保留。

## 来源与结果

source 必须仍绑定相同 mesh。previous velocity 的 geometry identity 以及完整 step/begin/end/position 与 source 一致；arrival layer.begin 等于 source 物理时刻，实际到达时刻可以是显式 Midpoint 或 End。零间隔只允许 Initial。

domain bounds 必须与 mesh 和 snapshot 精确一致，axis 语义与 snapshot 一致；cellSize 为有限正值且不超过各方向最小 native leaf size。速度 provider 遵守 I02 的不可变生命周期及同源 value/derivative 前置条件，不得在 callback 中修改 mesh 或 source。

返回只读 `HermiteGalsEdgeTransportResult`，没有默认空构造：

- `edges()`：本次选定的 `{vertices, phi, dn}`，顺序与请求一致；
- `sourceIdentity()`、`sourceTimeLayer()`；
- `arrivalVelocityGeometryIdentity()`、`arrivalTimeLayer()`。

arrival identity 表示所用速度的耦合几何，可与 source 相同或不同；它不是新的 Interface state identity。只返回选定 edge stage 数据；没有自动更新未选边，也没有发布新几何、修改相体积或推进 accepted 主时钟。

## 失败隔离与支持边界

整批目标的来源、方向、native master 身份及 midpoint 先验证，随后才进入 velocity callback。结果在局部 vector 中累积；全部成功并再次检查 mesh 身份后才返回。异常含 source/arrival identities、edge endpoints 和阶段。

source 的 vertices、edges、M、phi_xy、bubble、snapshot、identity/time 均不修改。失败不会发布部分结果。缺少物理边界历史 tangent、departure 离开旧活动几何、非有限一阶结果或 callback 失败时明确退出；corner-sector 的一侧几何也明确拒绝。没有反射、外推、裁剪、Hessian 输运或 legacy fallback。

静止/零步结果逐位保留旧 snapshot 的采样；已有 Core normalized metric cache 相对原 edge DOF 可有舍入差。I04 的实测最大 dn 差为 `5.55112e-17`，独立用 `8*epsilon*(1+|DOF|)` 检查。此处理与 I03 使用同一来源。

本阶段每次调用重新建立 native side 身份表，排序复杂度 `O(N_leaf log N_leaf)`，存储随 mesh 规模线性增长；未新增全局高阶状态。prepared topology 与并行优化留给数值 gate 之后，不在本次引入缓存生命周期。

## 验证范围

[edge contract](../tests/interface_hermite_gals_edge_contract.cpp) 使用真实 Q4 snapshot，在 56 条水平和 56 条竖直内部边上检查平移、两个方向 affine shear 及实际复合 scalar 的 normal 数值微分；在全部 144 条边上检查 XY/RZ 单步局部收敛与轴线。主精度 fixture 为严格正 Q4，负相库存真实为零。

这些证据只证明 edge first-jet 算子；不证明多步界面、曲率、体积闭合、破碎或电喷雾稳定性。原始数值、前序回归与未完成事项见 [I04 报告](../task/I04_Edge_First_Jet_Transport_Report.md)。下一项为 I05 compact mixed closure。
