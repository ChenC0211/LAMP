# H10i1b / 同源 native Q4 闭合裁剪矩

显式 Core 入口：`LevelSetMappedInterface::integrateNativeQuarticRegionMoments(source,nativeCellIndex,graphs,faceCuts,regions,axisymmetric,greenOrigin)`。[任务](../task/H10i1b_Closed_Clip_Moments_TODO.md)、[实施前设计](../task/H10i1b_Closed_Clip_Moments_Design.md)。Interface 正式消费入口是 `HermiteGalsPhaseRegionSource::measure`，Hydro 的真实 `CchMomentumMass` 沿既有 region owner 消费，不增加第二材料来源。

范围为一个实际 active native Q4 leaf 内的完整、无正面积重叠/间隙的轴对齐矩形分区。每个 Region 的四个顶点必须严格凸且循环排列，每条完整边都有唯一同源 NativeQuarticTrace；原输入方向可反转。相同 regionId 的矩形组成一个区域，内部 seam 不算 fixed 边界，内部界面仍计 free。不同 regionId 的重合相边界分别取各区域内部相迹；经原完整 zero-set normal-partial 证书证明唯一正规 graph 后仅计 fixed。没有 epsilon 距离判定、owner 采样投票、第二 Q4 或 rounded endpoint 重造。

调用者必须提供同源完整正规 graph inventory。普通曲线继续使用原 retained root ordinal、region enclosure checks、bounded refinement 和 controlled GL。无修正 Jacobian；图方向矩受控积分，另一方向由共同 shifted Green primitive 在同一 retained endpoint 上分部积分恢复。XY potentials 是 `(x-x0)` / `(y-y0)`；RZ 是 `pi*(r-r0)*(r+r0)` / `2*pi*r*(z-z0)`。误差是 floating quadrature estimate，不是严格区间或守恒库存证明。

结果 `regions` 按实际 regionId 返回 capacity、两相 axial/radial closed moments、liquid fixed/free contributions、errorEstimate 和 freeArcs。整体返回 exterior fixed liquid、内部不同区域的 coincident liquid 项、coincidentArcs、sourceRootIncidences 和 localRootRefinements。localRootRefinements 计实际 canonical endpoint 与 transverse root 的 refinement 工作；外部共享 cut 细化仍由 consumer 单独累加。重合项是几何账本：

`sum(region.fixedLiquid) = exteriorFixedLiquid + coincidentLiquid`。

它不能当作材料 F 的反号共享规则，也不改变旧 BoundaryMoments::assemblyClosure 的含义。同一区域人工 seam 的 fixed 和 coincident 项都排除；拼成一个区域后其相界面由 free 项唯一恢复。

`RootResolutionError` 的 segmentIndex 指 `faceCuts`，rootIndex 指该 typed cut 的 `events()` 中 value event。consumer 在 owning immutable typed query 上细化、重启所有相关矩，失败不修改原 source/cuts。来源不同、moved-from、缺边/缺 inside side、坏分区、未认证重合/临界 graph 或原有数值工作预算耗尽均拒绝。不存在裁负值、half 边界相、容差扩大或质量修复。

该接口不支持斜 polygon、不迁移 native 零边 topology publication 或 legacy sensitivity，也不构造材料通量/时空积分。Pure FirstJet 的 primary/Q4/transport/clock、库存为空的路线、四项材料硬约束和真实 R/uR/KEP 保持。H10 最终 PDE 与原 H10b 收口另行验收。
