# H10l1 / 同一实际速度场的均值闭合与不可穿透来源

2026-10-04；对应[设计](../task/H10l1_Spatial_Advector_Closure_And_Boundary_Design.md)、[报告](../task/H10l1_Spatial_Advector_Closure_And_Boundary_Report.md)。原 [H10j0](H10j0_Cch_Spatial_Advector.md) 的 geometry/actual owner、atoms/views、完整 clock 与 provenance 保持，新增接口仍由 `CchSpatialAdvector.h/.cpp` 和原 CCH target 持有。

## 物理均值与实际积分

在每个原 native face 的局部体积切向坐标上，令 `s=t-t0`、`L=t1-t0`，使用 `B=b+k*(s-L/2)`。XY 坐标不变；RZ 使用 `X=r²/2`、`Y=z`，径向 `B=r*u_r`。b 来自原物理面均值，k 来自原 one-ring weights。原 support/rank/selected tangent 缺口不变。

同一 face 的前缀积分为 `w*(b+k*(w-L)/2)`。完整面 `w=L` 的零均值切向贡献精确为零；cell 两侧分别在真实 face 上积分并反号求差，仍构成原 lift。原 D_I 独立核验实际 field divergence，误差门槛为 `1e-12*sum(abs(A*u))`，gross 为零时要求精确零，无经验 floor。实际 cellVelocity、atom polynomial 和 phase/material consumers 均来自此同一场，没有第二套散度或后置通量修补。

## 显式静态边界

```cpp
auto law = CchAdvectorBoundary::impermeable(mesh, momentumGeometry, boundaryDofs);
auto actual = reconstruction.evaluate(mesh, nativeValues, geometryIdentity, fullClock, law);
```

`boundaryDofs` 必须由调用方明确给出，是原 selected domain 上真实 native boundary DOFs；包括明确选择的 selected-domain 边界及 RZ axis。重复、越界、内部面、错误 mesh/regrid 或另一 momentumGeometry instance 均拒绝。每次 evaluate 要求指定面的输入 normal mean **严格为零**，即使非零值为 denormal 也拒绝。

law 只绑定静态 geometry 和排序后的 DOFs，副本保留同一不可变身份；它不持有速度、phase、clock 或接受状态。`momentumGeometry()`、`velocityDofs()`、`sameIdentity()` 与 `validateFor()` 用于只读来源核验。actual owner 保留 law 的共享副本，每次实际求值仍获得独立 advector identity 和原完整 clock。

`actual.boundary()` 返回只读 law 指针，未绑定时为 null；指针生命周期不超过 actual owner。原四参数 evaluate 继续表示未指定该 normal law。只有明确指定的 face 将 normal slope 固定为零，配合输入零均值使整个 normal polynomial 精确为零。mean==0 本身不代表墙面，其他内部或一般边界的变号 trace/counterflow 保留。

phase 与唯一材料 F 直接消费这些实际 polynomials；边界 atoms/dual subsegments 为零是上游同一场的结果，不在 H08/H09 另压 F。此 law 不提供外部 phase 数据，也不豁免原 selected/exterior causal 流入拒绝。

## 范围

本项验证 homogeneous impermeable normal 条件与 closure，不代表一般 H14/H15 边界/接受事务。[H10l2](../task/H10l2_Same_F_Dual_Material_Retry_Report.md)已单独完成same-F dual迁移并删除旧quarter；closed fixture显式消费此law。材料硬约束仅唯一共享 F、内部反号/实际 AMR 子面和、来源一致、KEP/skew；R、E_sp、G 可以非零，Pure FirstJet 不受库存或 bubble 约束。
