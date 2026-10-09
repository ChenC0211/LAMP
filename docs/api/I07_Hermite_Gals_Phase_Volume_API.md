# I07：相体积测量与守恒库存初始化

入口为 `Modules/Interface/HermiteGalsPhaseVolume.h`，namespace `Quantum::Interface`。I07 不实现时间推进；最终 production gate 仍见路线图 I08–I12。

## 数据流

```cpp
const auto measured = integrateHermiteGalsPhaseVolumes(mesh, initialTemplate, config);
auto candidate = measured.initializeCandidate(mesh);
// candidate.primary().phaseVolumes 是新 revision 的 M。
// candidate.hasGeometry() == false；I11 仍需验证并发布对应几何。
```

`initialTemplate` 是 I01 提供的、尚未经过 App 接受的模块 State：它拥有一份初始 native Q4 snapshot。初始模板的 M 不参与积分；测试中以零占位，不把它冒充满足几何体积约束的 accepted solution。这样先测量初始几何，再形成库存，不需要倒置为“先求解未知初始 M 的 bubble”。最终 primary 与 snapshot 的数值一致性及 production 初始化接线由 I11 验证。

积分函数只读，可以在 Initial / Midpoint / End 用作几何诊断；结果绑定原 State 的 identity/time，并保留同一不可变几何的所有权。`initializeCandidate` 只允许 `step=0, position=Initial, begin=end=0`。它拷贝 first-jet 主量并替换初始 M，返回新 identity、没有 geometry/closure 的 candidate，不替换 App 的 accepted handle。旧模板及结果不变。重复调用只会产生不同的未发布候选，不会推进或修改旧状态；restart、regrid 与 advancing layer 明确拒绝。

## 积分与边界

所有库存 cell 必须恰好由 `primary.geometryCells` 覆盖。函数使用相同 immutable snapshot 的完整 Core 零集发现、纯相证书和实际 native half-face 分段。`discoverVolumes` 不要求旧 InterfaceState。I12a11r8 将 native Q4 cut 体积接入与 PhaseRegion 共用的局部 Green 积分：原点为单元左下角，液/气直接积分后采用较小相及其补量，保留原 retained cuts、求积及细化限制。共享测量一致性和回归已验证，但运动曲率收敛仍失败，见 [报告](../task/I12a11r8_Shared_Local_Volume_Measure_Report.md)。

全局 fixed/free moments 仍以共同零原点计算，保留原共享面反对称账本；其浮点和可与局部体积不同，不用于覆盖局部结果。legacy `evaluateVolumes` 与 sensitivity 路径保持原行为。局部积分额外检查双 Green/容量相对缺陷及 moment estimate 不超过原区域默认 `1e-10`，cut 仍严格有界。

XY measure 为 `dx dy`；RZ 为 `2*pi*r dr dz`。实际 Cartesian leaf capacity 来自 `NativeLeafDomain`，包括 2:1 边的实际细分。没有第二套 scalar、PL 面积或 cell-center 判相；没有材料/速度依赖。

返回 `cells()` 为同 native 排序的 cell、capacity、原始 Green integral 和证书（-1 液体、+1 气体、0 cut）。`phaseVolumes()` 是这份几何**可用于新初始化的**体积值，不是原 State 已输运的 conserved inventory：

- 有完整纯相证书的 cell 使用严格 0 / capacity，并检验原始积分与证书一致。
- cut cell 直接使用原始积分，必须严格处于 `[0,capacity]`；越界拒绝，不 clipping。
- `rawLiquid` 始终保留，包括纯相积分的舍入差。

`config.boundaryContacts` 与 `symmetryPlanes` 传给现有 Core 完整拓扑发现。未声明的外边界相交、不满足 symmetry trace、退化/未解析的零集均拒绝，不删组件，也不切到旧方法。Core 当前要求每个 component 至少有三段 retained graphs；一段图的微小角部 cap 明确拒绝。已验证八段图、总面积 `2^-30` 的薄膜可保留。当前径向 snapshot 要求完整域从 `r=0` 开始，未实现环形截断域；离轴 torus 可以由轴上偶对称 Q4 描述并已验证。

## 容限与诊断

两个具名相对容限默认均为 `1e-10`，必须有限且在 `(0,1)`：

| 字段 | 归一化 / 含义 |
| --- | --- |
| `relativeLedgerTolerance` | 单元 boundary capacity / 纯相积分缺陷除以本 cell capacity；全局 fixed-boundary 组装闭合除以总 capacity |
| `relativeMomentEstimateTolerance` | 所有 free graph GL8 与细分估计差之和除以总 capacity |

这些是独立验收 gate。改变它们不调整 Core 求积规则；估计差不是严格误差上界，不能代替解析验证。cut 严格有界检查不使用这两个容限放行越界。

`report()` 包含液/气/cut cell 数、component/graph/root incidence 数、capacity 总量、raw 与初始化总量、误差估计、容量和纯相缺陷、边界组装闭合。额外保留原 primary 的 `totalSourceInventory`、有符号 `geometricInventoryDifference = sum(raw)-sum(M_source)` 及 `maxCellRelativeInventoryDifference = max |raw-M_source|/capacity`；这些诊断不修改库存或决定 bubble 修正。

## 验证边界

解析 circle/product-of-circles 在当前 Q4 fixture 中可再现，mixed 来自 I06，而 center samples 由测试解析函数供应；这不等于 I10 bubble closure 或 I11 publisher。非多项式水平集的体积误差包含表示误差，单独做 N=8/16/32 收敛。粗细网格逐 cell 积分求和一致只证明几何 measure 的可加性，不证明 I09 flux reflux 或动态 conservative remap。

见 [TODO](../task/I07_Conserved_Phase_Volume_TODO.md)、[报告](../task/I07_Conserved_Phase_Volume_Report.md) 与 `runs/baselines/I07_20261001/` 原始证据。
