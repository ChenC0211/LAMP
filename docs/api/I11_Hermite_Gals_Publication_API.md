# I11 / 权威 Hermite 几何发布 API

**状态：PASS，2026-10-01。** 本接口完成模块级发布，不执行 coupled accepted-step commit。Gate I 仍须 I12。

## 输入与结果

入口位于 `Modules/Interface/HermiteGalsPublication.h`：

```cpp
auto mixed = closeHermiteGalsMixedAmr(mesh, candidate);
auto bubble = closeHermiteGalsBubble(mesh, candidate, mixed, prediction);
auto publication = publishHermiteGalsGeometry(mesh, std::move(candidate), bubble);
```

prediction 属于调用方提供的当前 primary revision。bubble 必须为 I10 返回的不可伪造结果；publisher 核对 mesh、identity、未绑定状态。主量是 master vertex/edge first jet 与 M；mixed/c 仍为派生量。I10/I11 使用 zero-center cardinal extension，故在本基准下 c 等于 center scalar。不能把其它 extension 下的系数直接代入。

成功后 candidate 被消耗，publication 保持其原 identity、time layer、measure 和所有 primary 值。`std::move` 只绑定右值引用；函数内全部检查通过前不移动或绑定调用方候选。失败抛出带 geometry identity、stage 的异常，候选仍可检查、修正后生成新 revision 或以同一输入重试。失败不修改旧 accepted handle。调用方须在所有模块验收后自行提交 accepted state。

I01 `bindGeometry/finish` 保留为结构性低层端口与测试 fixture；它们不是 I11 数值验收证明，也不能构造 `HermiteGalsPublication`。后续生产接线须使用受检 publisher。

## 一个不可变 scalar owner

publisher 在私有同 revision candidate 副本中构造唯一 native Q4 snapshot；副本不 mint 新 identity。所有验收直接作用于最终返回的 snapshot，不先验收一份再重建另一份。发布结果的复制共享 Data；移动后的源句柄拒绝访问。

| 只读入口 | 含义与所有权 |
| --- | --- |
| `state()` | I01 不可变状态；M 保留原始值 |
| `volumes()` | I07 在同一 snapshot 上的几何测量；不是更新后的 M |
| `domain()` | 与候选 active cells/XY 或 RZ 度量对应的 native leaf domain |
| `topology()` | 同源完整 `NativeComponents`，允许多个分量和纯相空集合 |
| `facePartitions()` | 实际 vertex-midpoint 半面的 retained cuts；source 指针与 state 相同 |
| `sample(component, graph, fraction)` | 同一 native chart 上的 raw normal/curvature，附 identity/time |
| `quadrature(component, graph)` | 缓存同源 GL8 节点、`ds` 和 XY `ds` / RZ `2*pi*r*ds` 权重 |
| `validateMesh(mesh)` | 检查状态 mesh identity 与缓存 domain 所有权 |
| `report()` | 事后统计；不作为数值分支输入 |

MappedFacePartitions 保留 domain 地址，因此只在 Data 的最终存储位置构造 cuts。Data 不再移动，用户 copy/move 只改变 shared handle。结果拥有 snapshot 与所需派生数据；网格和输入变量销毁后仍可读取 surface samples。对外部 mesh 的消费仍须显式校验 identity。

half-face 的 Range 保留 `fromEnd/offset/length/negative` 语义，Root 保留 native cell/chart/enclosure。不要将它们转换为两个全局显示端点再相减。I08 消费 `publication.state()` 时也保留相同 snapshot 指针；本任务验证了 published midpoint 的零速度库存恒等性。

## 发布门槛与曲率符号

- 首先使用 I07 的既有 ledger、capacity、moment estimate 验收；逐单元再次检查 M 严格有界及 `abs(rawLiquid-M)/capacity <= relativeInventoryTolerance`，默认 `1e-10`。M 不裁剪、不回写。
- 完整 native topology 与所有实际半面 cuts 必须成功构造，且 source 与唯一 snapshot 相同。
- 每个 graph 两端和所有 quadrature nodes 检查 interface residual、unit normal、曲率有限性。残差门槛为 `relativeInterfaceResidualTolerance * min(cell width)`，默认 `1e-9`；normal 范数缺陷默认 `1e-12`。公开 `sample()` 也执行这些检查。
- `Frame.normal` 从负相指向正相，`Frame.curvature` 是同一 scalar gradient/Hessian 给出的 `+div(n)` 平面/子午曲率，未平滑或拟合。液滴内负外正时，XY 圆为 `+1/R`，RZ 球总曲率为 `+2/R`；反号时两者同时反号。
- RZ 非轴处 `curvature = Frame.curvature + n_r/r`；r 从保留的 endpoint-relative chart 计算，不反演舍入后的全局点。真正 axis endpoint 由 topology endpoint kind 与精确 chart 同时识别，要求 `r=0,n_r=0`，采用偶对称极限 `2*Frame.curvature`，没有近轴阈值。
- quadrature 权重有限且严格为正；纯相无 interface graph、sample 或伪曲率。报告汇总 component/graph/node/half-face/axis 数和已测最大误差。

配置值均须有限且在 `(0,1)`；体积配置中的 boundary contacts/symmetry planes 沿用 I07 语义。GL8/subdivision discrepancy 仍是误差估计，不是严格误差上界。发布检查不证明每个连续点上的全域数值误差界。

## 验证边界

见 [I11 报告](../task/I11_Authoritative_Geometry_Publication_Report.md) 与 `runs/baselines/I11_20261001/`。该任务没有新增 App reader、主时钟、Hydro/Electricity 接线或 VTK 输出流。

I12r5 支持 native Q4 面内孤立重根的相区间；I12r6 支持同一标量的正则孤立接触完整连接。实际边与内部 chart 边保留全部不同根，根据原标量相位及根包围内非零法向导数证明单元内曲线关联，负相片只跨正长度面区间合并。面上的无重因子仅用于零集身份，不能替代原梯度/Hessian；transverse integration sections 仍要求原标量简单根。普通负相变分在需要重根认证的面查询上仍拒绝。

新增解析验证涵盖 XY/RZ、uniform/AMR、正反相、角点、近相切两侧及二/四重非二进制接触位置，使用 I07/I10/I11 实际发布链。发布的固定 GL8 节点没有自适应弧长精度保证；需要受控弧长/面积时，可用同一 `NativeGraphChain` 的既有 `LevelSetInterfacePanels` material adapter。见 [I12r6 报告](../task/I12r6_Native_Isolated_Contact_Report.md)。

正长度零迹、真正奇异零点及未解析 closed subcell cycle 仍拒绝；I10 signed-partial/位移证书是充分条件。任意步内跨面通量的时间精度、完整 GALS 步进、regrid、实际拓扑变化与低黏电喷雾稳定性须后续任务验证。体积/拓扑/切面派生缓存仍分别构造，尚未实施跨模块缓存复用的性能优化。
