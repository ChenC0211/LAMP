# I12a2 / 平面跨面时间积分 API

`HermiteGalsTemporalVelocitySource::constant(old, arrivalId, arrivalTime, value)` 持有常速度值并生成实际 previous/arrival first-jet queries；Jacobian 为零。复制共享原 owner。既有任意 callback 构造器保留，但 `constantValue()` 返回空，不会靠采样自动获得常速度资格。

`integrateHermiteGalsPlanarFace(mesh, interval, nativeFace)` 消费 I12a1 区间并返回 `HermiteGalsIntegratedFaceTransfer`。输入旧 snapshot 是唯一标量来源；不创建时间节点几何，不更新库存，不接 driver。方向来自实际 native edge，两个邻居消费同一正 x/正 y 面量。

Core 的 `Snapshot::integrateNativeAffineTranslatedFace(a,b,velocity,referenceTime,begin,end)` 只作几何积分。它在现有 strict TU 内：

1. 以 exact expansion 计算从旧时刻至 end 的整段扫掠盒，并以 outward enclosure 查找全部相交原生叶子。
2. 从每个叶子的 exact homogeneous Q4 DOFs 逐项证明全部 25 个约束符合 affine polynomial，且跨叶子物理多项式完全相同。rounded diagnostic DOFs 不参与证明。含非仿射扰动的近似平面仍拒绝。
3. 将 face 参数 s 与区间参数 tau 映射为单位正方形。原标量在此域为 affine；用精确符号、边交点的有理表达裁剪负相半平面，完整处理所有 crossing events。
4. 以多边形解析面积及一阶 s 矩积分 XY/RZ 测度。RZ 权重为 2*pi*r；接近纯液的比例用同一 exact rational 的气相补集求值，避免近一舍入越界，不作库存裁剪。

Core 返回总时空测度和液相比例；Interface 乘以常法向速度的绝对值，分别写入对应正/负 gross volume。signed volume、平均 rate 和受检区间组合沿用 I12a1。

能力边界：仅保留域内部、active native Q4、全扫掠盒同一平面；RZ 仅常轴向平移。未知 `activeSegment` callback、外部边界历史、非平面或正时空测度恒零 trace 明确拒绝。只在单个时刻整面过零的事件为零测集，正常积分。这里的 exact 指平面证明、事件谓词和有理矩；返回的浮点量及 pi 仍有表示舍入。

这不是一般运动界面的时空求积器，也不证明库存/bubble 闭合或 I12 已通过。测试与证据见 [TODO](../task/I12a2_Planar_Crossing_Integral_TODO.md) 和对应报告。
