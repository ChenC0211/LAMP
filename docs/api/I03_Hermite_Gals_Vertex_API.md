# I03 / Hermite–GALS master vertex first-jet API

状态：I03 已验收，2026-10-01。只实现选定 master vertices 的一阶输运；完整 skeleton、守恒闭合、几何发布与 Gate I 尚未完成。

接口：[HermiteGalsVertexTransport.h](../../include/Modules/Interface/HermiteGalsVertexTransport.h)。实现：[HermiteGalsVertexTransport.cpp](../../src/Modules/Interface/HermiteGalsVertexTransport.cpp)。前置类型见 [I01 state](I01_Interface_V2_State_API.md)、[I02 characteristic](I02_Hermite_Gals_Characteristic_API.md)。

## 单一源与公式

```cpp
HermiteGalsVertexTransportResult transportHermiteGalsVertices(
    const QuadtreeMesh& mesh,
    const HermiteGalsState& source,
    const std::vector<int>& masterVertices,
    const LevelSetBoundaryTrace::Domain& domain,
    const HermiteGalsVelocityView& previousVelocity,
    const HermiteGalsVelocityView& arrivalVelocity);
```

对 mesh 中的真实 arrival 坐标 `x_v`，调用 I02 得到同一离散 midpoint 回溯的 `B_h(x_v)` 与 `J_Bh(x_v)`。随后从 `source.geometry().snapshot()` 一次查询 `scalarJetNative(B_h)`，计算：

\[
\phi_v^{new}=\phi^{old}(B_h(x_v)),\qquad
\nabla\phi_v^{new}=J_{B_h}(x_v)^T\nabla\phi^{old}(B_h(x_v)).
\]

`J` 的行是 departure 分量、列是 arrival 坐标；因此 `dx_new=J.xx*dx_old+J.yx*dy_old`，`dy_new=J.xy*dx_old+J.yy*dy_old`。梯度使用物理坐标，无归一化、有限差分替代或另一份 scalar 重建。

固定点也走完整映射与拉回计算：`B_h(x_v)=x_v` 不意味着 `J_Bh=I`。静止或零步时结果逐位保留旧 snapshot 的采样；Core 的浮点 metric cache 与 primary 导数可以相差舍入量，内核不绕过 snapshot 去取另一来源。

## 输入与来源

| 输入 | 约束 |
| --- | --- |
| source | I01 已发布只读状态，仍绑定相同 native mesh；整个调用期间保持有效 |
| masterVertices | 非空、严格递增、无重复的 native vertex indices；每个属于 source 的 master 列表且 `TVertType::non` |
| previousVelocity | geometry identity 与 source 相同；完整 step/begin/end/position 与 source time layer 相同 |
| arrivalVelocity | 显式 time layer，其 begin 等于 source 物理时刻；前进间隔必须可表示，零间隔只允许 Initial |
| domain | 矩形端点与 mesh/snapshot 精确一致，axis 语义与 snapshot 一致；正且有限的 cellSize 不超过最小 native leaf 尺度 |

masterVertices 是明确的更新子集。调用方负责选出本阶段所需 master 集合；部分结果不意味着未选中的节点已经推进。hanging nodes 不可独立更新。全网格 AMR master/constraint 闭合仍属于 I06。

arrival velocity 的 identity 表示它所依赖的耦合几何，可与 source 相同或不同；它不是本函数生成的新 Interface identity。速度 callback 必须遵守 I02 的不可变、成对 value/derivative 插值前置条件，不得通过捕获别名修改 mesh 或来源状态。

I02 的半单元采样位移、子步控制、boundary events、activeSegment 和 RZ parity 规则继续适用。minimum leaf scale 在本次调用开始计算，保守约束所有选定顶点的 tracer；这里尚无局部 CFL 优化或并行 callback 查询。

## 结果与失败隔离

返回值拥有独立的 vertex vector，通过 const reference 读取，不能默认构造。它带有：

- 实际更新的 `{vertex, phi, dx, dy}`；
- `sourceIdentity()`、`sourceTimeLayer()`；
- `arrivalVelocityGeometryIdentity()`、`arrivalTimeLayer()`。

所有目标先完成合法性检查，所有数值结果暂存于局部 vector；只有整批成功且最终 mesh 身份复核通过才返回。异常消息包含 `source`、`arrival`、`vertex` 及失败阶段。source 的 primary、edge、相体积、closure、身份、时间和 snapshot 均不修改；没有候选状态、snapshot 或 accepted step 的提前发布。

以下情况明确失败：非法目标或来源、旧几何已不匹配 mesh、物理边界命中却无 boundary-history tangent、departure 位于旧活动几何之外、非有限采样/拉回结果、速度 callback 失败。含 corner-sector 的几何也拒绝，因为单个共享 vertex first jet 不足以表达已有的一侧数据。此处没有反射、外推、裁剪或 legacy fallback。

## 验证与边界

[contract](../tests/interface_hermite_gals_vertex_contract.cpp) 使用真实 Core native Q4 snapshot：解析平移/剪切覆盖 49 个内部顶点；复合 scalar 的中央差分独立检查 gradient；XY logistic 与 RZ 奇对称 cubic radial 解析速度覆盖全部 81 个顶点，包含边界固定点与轴线。

主精度测试使用严格正 Q4 scalar，`M_K=0` 是其真实负相库存，检验的是 first-jet 算子。它不检验界面破碎、守恒通量、曲率或完整多步误差。单步误差随 dt 减半趋向三阶，对应 RK2 的局部截断阶；完整 Interface 时间收敛仍需 I12。

此次不输运 `phi_xy` 或 Hessian，不推进 edge/mass/bubble，不构建新的 geometry。I01 目前只做结构绑定，primary 与 polynomial 的生产一致性还需 I05–I11 完成。定量结果、失败路径和未解决问题见 [I03 报告](../task/I03_Vertex_First_Jet_Transport_Report.md)。
