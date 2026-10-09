# I12a14 / 特征相与库存含量覆盖 API

消费 [I12a13 characteristic enclosure](I12a13_RK2_Characteristic_Enclosure_API.md)，遵循[跨叶覆盖设计](../task/I12_Characteristic_Phase_Coverage_Design.md)。结果描述同一旧 Q4 标量在已认证 departure box 中的严格相及有界液体含量，不更新 M 或 geometry，不生成通量。

两个入口均返回不可默认构造的 `HermiteGalsCharacteristicPhaseCoverage`：

```cpp
coverHermiteGalsCharacteristicPhase(mesh, enclosure);
coverHermiteGalsCharacteristicPhase(mesh, enclosure, inventoryRegionSource);
```

第一种入口对应几何指标函数，负相含量 1，正相含量 0。第二种显式采用 r18 库存区域的 wl/wg；其 old state identity 与 snapshot 指针必须都与 enclosure 相同。接口不会从相同采样值推断来源相同，也不会重新计算这些权重。

| 查询 | 含义 |
| --- | --- |
| `characteristic()` | 同一 I12a13 包围、interval、old geometry、velocity 与 characteristic owner |
| `cells()` | 按 native cell 排序的所有 closed departure-box 交集，含零宽度 seam/point |
| `phase()` | 所有叶严格同号时为 Liquid / Gas，否则 Unresolved |
| `minimumLiquidContent()` / `maximumLiquidContent()` | 所有叶含量范围的整体 hull，位于 [0,1] |
| `inventorySource()` | 保留的只读库存 source；几何入口为 null |
| `validateMesh(mesh)` | 沿用精确 native 几何/拓扑身份校验 |

每个 `HermiteGalsCharacteristicPhaseCell` 含 `cell`、旧叶与 departure box 的闭交集 `lower/upper`、`phase` 及含量上下界。Liquid 的范围为 [wl,wl]，Gas 为 [wg,wg]，Unresolved 为 [wg,wl]。几何入口使用 wl=1,wg=0。纯气相的库存含量可以大于零；跨叶几何同相时，整体含量也未必是一个常数。

`HermiteGalsCertifiedPhase::Unresolved` 表示尚未证明整个 box 为严格单一相。它可能包含真薄层/相切，也可能仅因保守 box、Bernstein 区间过宽而未确定；不能当成零通量、气相、已命中或已证明无事件。若之后细分，必须保留同一完整来源。

并集覆盖以 x 端点划分闭条带，再检查每条带的 y 闭区间并集。不累加面积、不取 donor centroid、不忽略零测度的闭接触，也不依赖容差。每叶严格相使用 Core 原生 `nativeQuarticBoxSign`。结果对象只读且保留 owner，复制共享；源的局部变量销毁后仍可查询，移动后的对象拒绝访问。错误的 mesh、来源、活跃标量支持、浮点环境或覆盖空隙明确拒绝。

此 API 不处理 I12a13 尚未认证的 boundary/substep 区域，不认证连续 ODE，不执行时空 quadrature，也不证明一般库存通量更新有界性。完整 I12 gate 尚未通过。
