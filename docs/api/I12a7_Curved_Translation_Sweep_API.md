# I12a7 / curved Q4 常速度扫掠积分 API

`integrateHermiteGalsSweptFace(mesh, interval, nativeFace, regions, config)` 使用同一区间的 owned constant velocity、characteristic owner 和旧 accepted scalar。`regions` 是 I12a6 的 prepared source，必须同时具有相同的 primary identity 和 Snapshot 地址。

支持无切向分量的常速度正常穿面平移：XY 的 x/y 两方向、RZ 的轴向。零法向速度的结果严格为零。非零切向速度、opaque velocity/coverage callback、RZ 径向运动明确拒绝。原 I12a2 精确平面入口不变；两个入口互不 fallback。

Core 的 `Snapshot::nativeNormalTranslationSweep` 从精确 expansions 计算旧坐标 `x − u*(time−referenceTime)`，给出内/外 double 矩形，几何 totalVolume 与 coordinateVolumeBand。内矩形包含于精确扫掠，精确扫掠包含于外矩形。暂要求原 native Q4 域完全 active，当前面与外扫掠矩形在原域内且有严格正宽度的可分辨内矩形。出域需要未来的 BC 历史接口，本入口不外推。

I12a6 用原 scalar 的完整曲线分解分别测量这两个矩形；包围重合时只测一次。积分相量使用两次直接测量的平均，并用较小相稳定表达互补量。这里“完整”指保留旧多项式的全部曲线/裁剪根，不能靠几个时间采样全为气相就跳过事件。

返回 `HermiteGalsSweptFaceResult`：`transfer` 是原 I12a1 四个 gross integrated volume，保留 interval/face authority；`report` 包含两次直接相测量、矩形查询/arc/根细化计数，以及分开的 coordinateVolumeBand 和 momentErrorEstimate。法向速度不变号，只有对应一个方向的 gross total/liquid 非零。平均 rate 仍由 transfer/该 interval 的真实 dt 派生。

坐标包围带只界定几何位置不确定性；两次液量是浮点求积，不能把 innerLiquid/outerLiquid 叫作严格液量下上界。momentErrorEstimate 是两次 GL8 估计的平均，也不是严格误差界。默认两个相对门槛均 1e-10，超过门槛或 gross 越界直接失败；这些门槛不修改 I12a3 的库存舍入门槛。

此入口不更新 M、不发布 midpoint/end、不建立新时间节点几何、不更改原 RK2/first jet，也不提供一般变化速度或边界历史。区间切分沿用同一旧源，子区间仍有真实 begin/end。

任务与验收状态见 [TODO](../task/I12a7_Curved_Translation_Sweep_TODO.md) 和 [报告](../task/I12a7_Curved_Translation_Sweep_Report.md)。
