# I09 / Hermite–GALS 2:1 reflux API

状态：I09 PASS，2026-10-01。对应冻结方法式 (34)、I-R3。Gate I 尚未完成。

## 调用与权威结果

```cpp
#include "Modules/Interface/HermiteGalsReflux.h"

const auto flux = Quantum::Interface::advanceHermiteGalsPhaseVolumes(
    mesh, previous, midpoint, velocity);
const auto reflux = Quantum::Interface::assembleHermiteGalsReflux(mesh, flux);
const auto& updatedM = reflux.source().phaseVolumes();
const auto& nativeFlux = reflux.source().faces();
const auto& sides = reflux.sides();
```

`source()` 保留原 `HermiteGalsPhaseFluxResult` 的同一不可变 allocation；faces、updated M、previous identity、midpoint geometry identity、时间区间和 snapshot 都不复制成另一份 authority。销毁调用方的 source handle 后仍有效。新 `FluxResult::validateMesh(mesh)` 转发已有 midpoint state 的精确 mesh identity 校验；组装前必须通过。`RefluxResult::validateMesh` 使用同一入口。

本方法共用一个 midpoint interval；I08 已将每个实际子面的同一 `dt*Q_l` 反号写入 coarse/fine cells。I09 显式提供完整粗面的细面汇总，没有独立粗面近似通量可再修正，也不得再次给 M 加 reflux correction。未提供局部时间步/subcycling flux register。

## 完整侧边语义

`sides()` 按源库存的升序 cell index，每 cell 按 `West, South, East, North` 排列；cell 为 native index，不是 label。

| 字段 | 定义 |
| --- | --- |
| `subfaceCount` | 1：conforming/fine side；2：2:1 coarse side |
| `subfaces[0..count)` | `source().faces()` 中的索引，按正切坐标升序；不是 packed native face ID |
| `outwardSign` | west/south 为 -1；east/north 为 +1 |
| `totalPositiveRate/totalNegativeRate` | 同一方向子面的最终 `Q+ / Q-` 之和 |
| `liquidPositiveRate/liquidNegativeRate` | 同一方向子面的最终 `Q_l+ / Q_l-` 之和；包含 I08 微修复 |
| `totalRate/liquidRate` | 子面已缓存的 signed rate 之和，仍以全局 +x/+y 定向 |
| `liquidVolume` | 子面已缓存 `dt*Q_l` 的直接和，不重新计算 `dt*sum(Q_l)` |
| `integrationErrorEstimate` | 子面原估计的直接和；不涵盖时间、标量重构或根位置误差 |

例如外向液量为 `side.outwardSign * side.liquidVolume`。外向 gross 出流在 sign=+1 时用 positive，在 sign=-1 时用 negative。signed 总和与两个 gross 总和之差可能因求和次序有末位差异；有限体积账本应使用已有 `liquidVolume`。

所有汇总只读取 I08 最终面条目，不读 `rawLiquid*`。细侧引用一个实际子面；粗侧引用同样的两个实际子面。因此 `Q_coarse = Q_fine,0 + Q_fine,1` 按相同顺序直接成立。这里是积分通量相加，不是 normal velocity 或相体积分数的面积平均。RZ 的 `2*pi*r` 度量已在 I08 积分中进入结果，I09 不重复乘 metric。

## 拓扑、诊断与失败

复用 `FaceSegmentGraph` 的 cell side/实际 leaf edge 关联，不另做坐标近邻判定。组装核对每侧 1/2 个子面、面 owner/方向、完整 incidence、split-side 邻居层级为 coarse+1，以及选定域内 fine side 的反向关联。选定域外的邻居仍按实际子面处理，只有域内完整两细侧参与 internal residual。当前数值测试使用完整计算域。

结构 gate 通过后才生成 `HermiteGalsRefluxReport`。`coarseSides` 统计所有两子面侧边；`interiorCoarseSides` 统计两个细邻居均在选定域的粗面；`coarseFinePairs` 统计选定域内粗细配对。三项 `max*Residual` 分别为最大绝对

```text
abs(coarse.outward + fine0.outward + fine1.outward)
```

对应总速率、液相速率和液相输运体积，未归一化；没有隐藏 tolerance 或通过 Report 改变通量的分支。

改变网格、moved-from handle、缺失/重复/不匹配拓扑或非有限汇总均拒绝；源结果与原状态不变。内核不读字典、不写文件、不接受可写 face 输入或积分回调。存储规模为每 cell 四侧和每侧最多两条索引，不增加新的物理未知量。

## 验证边界

见 [I09 报告](../task/I09_Coarse_Fine_Reflux_Report.md)。XY/RZ、uniform/两类 AMR、四方向、counterflow、修复后的粗细通量、纯相/零速度及 step-end-aligned 多步 AMR 穿越通过。

多步穿越测试使用解析 midpoint geometry，并把整面穿越时刻放在步端，验证库存通过 actual subfaces 正确跨网格传递。它不证明任意穿越时刻的时间精度或 I12 完整 GALS 输运。早期 off-step 穿越探针误差约 `3.9e-5`，暴露出中点 Heaviside 时间求积并不精确；该限制仍需 I12 的时间收敛/有界性验收。Core 的半面内部精确相切/重根、正长度零迹拒绝限制仍在；I10 bubble、I11 同源几何发布及最终低黏电喷雾目标未完成。
