# I05 / 均匀网格 Hermite mixed closure API

状态：I05 已验收，2026-10-01。公有接口：[HermiteGalsMixedClosure.h](../../include/Modules/Interface/HermiteGalsMixedClosure.h)，实现：[HermiteGalsMixedClosure.cpp](../../src/Modules/Interface/HermiteGalsMixedClosure.cpp)。此阶段从 primary first jet 恢复派生 `phi_xy`，不推进 Hessian、相体积或几何。

## 来源与输出

```cpp
HermiteGalsMixedClosureResult closeHermiteGalsMixedUniform(
    const QuadtreeMesh& mesh,
    const HermiteGalsCandidate& candidate);
```

candidate 是 [I01](I01_Interface_V2_State_API.md) 的未发布 primary revision。调用前必须满足：

- native mesh 身份仍与 candidate 一致；候选尚未绑定 geometry；
- mesh 为完整、非周期、同层 Cartesian uniform leaves，每个方向至少两个单元；支持各向异性、平移域和整体同层 refinement；
- primary 包含全部 native vertices 和真实 native leaf edges；拒绝缺失、假长边、AMR/hanging；
- vertex first jet 与 edge normal derivatives 属于 candidate 的同一时间层。水平 edge.dn 为 phi_y，竖直 edge.dn 为 phi_x。

通常由调用方将 [I03](I03_Hermite_Gals_Vertex_API.md)、[I04](I04_Hermite_Gals_Edge_API.md) 同阶段结果装入新的 I01 candidate，然后调用本函数。闭合不接受分离的旧 snapshot/closure，也不自行合并来源不明的数组。

结果拥有不可变 mixed vector，`vertexMixedDerivatives()` 与 candidate.primary().vertices 同序，值为物理坐标下的 phi_xy。`primaryIdentity()` 和 `timeLayer()` 标记它所依据的 revision。candidate 之后 replacePrimary 会产生新 identity；旧结果不会自动重新绑定，后续 publisher 必须核对 identity。

candidate、first jet、edge values、相体积和几何状态完全不修改。输出不是完整 `HermiteGalsClosure`：cell bubble coefficients 仍需 I10 决定；同源 geometry publication 属于 I11。

## 离散公式

冻结方法 §5.2 式 (27)–(28) 指定：

\[
D_x\phi_y=\frac{\phi_y(-h_x)-8\phi_y(-h_x/2)+8\phi_y(h_x/2)-\phi_y(h_x)}{6h_x},
\quad
\phi_{xy}=\tfrac12(D_x\phi_y+D_y\phi_x).
\]

x 线上读取 vertex.dy 和水平 edge.dn；y 线上读取 vertex.dx 和竖直 edge.dn。不差分 phi，不再构造 gradient，不使用速度二阶导数。

实现使用 [CompactQuarticDerivative.h](../../include/Core/Geometry/CompactQuarticDerivative.h) 的五点四次插值导数。内部点取 `[-h,-h/2,0,h/2,h]`，中心权重为零，等价于上述公式。物理下边界取 `[0,h/2,h,3h/2,2h]`，权重为 `[-25,48,-36,16,-3]/(6h)`；上边界用镜像 `[3,-16,36,-48,25]/(6h)`。这是同一四次插值的单侧求导，不读取域外值，也不指定物理 inflow history。

节点使用真实 native 坐标，midpoint 由实际端点产生；权重在以目标点为原点、支持区跨度为尺度的坐标中计算，因此不会用两个较大绝对坐标直接构造高次幂。求和前减去目标值，常数导数严格为零。uniform 坐标因 binary64 表示存在舍入时，仍按实际点求导，保持 Q4 再现。

该纯数值函数从已有 Core GALS predictor 原样提取。旧 predictor 通过原有 require 回调保留错误上下文；V2 只调用这个求导算子，不调用旧 predictor，也不与它共享可写状态。

## RZ 轴线

AxisymmetricVolume 要求非负半径。域含 r=0 时，所有轴线 vertex.dx 与竖直 edge.dn 必须严格零；非零值明确失败，不作修正。

偶标量对应偶 phi_z；其轴线中心反射导数 `D_r phi_z` 严格为零。轴线上 phi_r 的 vertex/edge 数据也严格为零，因此 `D_z phi_r=0`，最终 phi_rz 精确为零。非轴线顶点按相同四阶求导闭合。该 parity 是 RZ 正则性条件，不是额外可调混合权重。

## 诊断、失败与范围

`report()` 返回顶点数和：

\[
\max_v\frac{|D_x\phi_y-D_y\phi_x|}{\max(1,|D_x\phi_y|,|D_y\phi_x|)}.
\]

该 disagreement 只做诊断；始终采用 1/2 与 1/2，不因诊断值改变权重、接受条件或算法。诊断采用先缩放再相减，避免有限的大数在差值中溢出。

不完整或非均匀输入、已发布候选、mesh 身份改变、轴线 first jet 不合格、奇异支持点及非有限数值结果明确失败。错误包含 primary identity、阶段和 vertex index；整批成功才返回 mixed vector，不提前写入任何 state。I01 已经拒绝负半径候选，I05 同时保留非负半径入口防御。

存储量随 skeleton 规模线性增长；native edge 映射使用有序表，当前总工作量为 `O(N log N)`，没有全局高阶未知量或额外 scalar。缓存/并行优化不在本任务范围。

本入口要求完整 uniform skeleton；active-band 局部缺失、AMR master/hanging closure、contact/sector 一侧数据、守恒体积/bubble 及新 geometry 均由后续任务处理。[验收报告](../task/I05_Uniform_Mixed_Closure_Report.md) 给出了 Q4、四阶空间收敛、输运接入及旧实现不变证据，不能据此声明完整 Interface/曲率/破碎已通过。
