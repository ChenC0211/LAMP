# I12a11r18 / 有界库存区域积分 API

## 当前适用范围

2026-10-02：用户已授权转[纯 Hermite–GALS 非守恒主线](../task/I12_Pure_Hermite_Gals_Design_Change.md)。本 API 保留给显式历史守恒库存路径，纯 FirstJet 推进不消费此库存区域权重，不以其认证完成作为接受条件。


`HermiteGalsInventoryRegionSource(mesh, regions, config)` 在原 `HermiteGalsPhaseRegionSource` 上准备库存权重。新 owner 共享原不可变几何/state，并持有各 native leaf 的容量、几何液量、原 M、两相权重及恢复残差；没有新几何、库存更新或接受步发布。

输入要求完整的 native Q4 几何与库存支持，原 M 严格位于原 native 容量界限内。`config.relativeGeometryInventoryTolerance` 默认 `1e-10`，限制 `abs(V_geometry-M)/C_native`；它不改变守恒 roundoff。原/几何容量差与整格库存恢复缺陷均须不超过原 `64*epsilon*C_native`。详细公式与范围见[设计说明](../task/I12_Bounded_Inventory_Transport_Design.md)。

| 入口 | 语义 |
| --- | --- |
| `source()` / `regions()` | 同一不可变 state 与原 prepared 几何区域算子 |
| `cells()` | native 顺序的容量、V/G、原库存、液/气几何区权重、整格表示量 |
| `report()` | 实际几何容限、原/表示库存、几何差 L1、最大容量/几何/表示缺陷与调整叶计数 |
| `measure(mesh,lower,upper)` | 对同一矩形返回原几何测量、库存 liquid/gas、有符号与 L1 修正 |
| `validateMesh(mesh)` | 原严格 mesh identity 检查；等价网格允许，拓扑/坐标改变拒绝 |

`HermiteGalsPhaseRegionMeasure::leaves` 新增每个实际相交旧叶的 cell、裁剪边界、容量和 canonical 液/气量。它们来自本次同一裁剪/求积，使用既有 smaller-phase/complement 规则；原区域标量结果保持不变。库存查询消费这些贡献，不重复积分或重建几何。

准备权重时，构造器先完整校验 source/mesh，再通过仅 friend 可见的私有入口批量测量旧叶，避免 N 次全网格扫描。公开查询仍逐次校验网格，没有公开 unchecked API。结果拥有其来源；copy 共享只读 Data，moved-from 访问拒绝。查询按值返回的数据不构成 prepared source 的可写别名。

库存权重严格处于 `[0,1]`，子区域的库存 liquid/gas 严格有界。整格恢复缺陷为单独的 roundoff 诊断，不能当成把 M 重置为测量值的许可。全为几何恒等权重时直接保留原区域标量结果。

本算子可以表达几何纯相中的微量库存差，但不改变 r17 的 pure 几何闭合门槛；现有完整驱动未切换。后续仍需将它接入唯一时空 face owner、库存装配与 reflux，并通过原 I12 验收。

参见 [TODO](../task/I12a11r18_Bounded_Inventory_Region_TODO.md) 与[报告](../task/I12a11r18_Bounded_Inventory_Region_Report.md)。
