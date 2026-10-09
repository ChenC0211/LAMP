# H07b / 同源正权分相求积与 CCH 黏性

状态：2026-10-03，**H07b / H07 PASS**；[任务书](../task/H07b_Phase_Quadrature_Viscosity_TODO.md)、[报告](../task/H07b_Phase_Quadrature_Viscosity_Report.md)。采用[完全非守恒材料路线](../task/H07_Pure_Gals_Hydro_Coupling_Design.md)，[H08 材料流量/非零缺陷](H08_Cch_Material_Transport.md)已 PASS，下一项 H09：同源dual材料与动量一致性。本项是几何/黏性算子接线，不是完整两相 NS 或时间步。

公开头：[HermiteGalsPhaseQuadrature.h](../../include/Modules/Interface/HermiteGalsPhaseQuadrature.h) 属于正式 `quantum_interface`；[CchPhaseViscosity.h](../../include/Modules/Hydrodynamics/CchPhaseViscosity.h) 属于 `lamp_hydrodynamics_cch`。依赖只读地从 CCH 指向正式 Interface；没有旧 Hydro/Compatibility。

```cpp
auto geometry = CchPhaseGeometry::build(mesh, metrics, publication);
auto quadrature = HermiteGalsPhaseQuadrature::build(mesh, publication);
auto model = CchPhaseViscosity::build(mesh, geometry, quadrature, muLiquid, muGas);
auto force = model.apply(mesh, nativeValues, geometryIdentity, timeLayer);
auto power = model.dissipation(mesh, nativeValues, geometryIdentity, timeLayer);
```

## 几何求积

唯一 immutable FirstJet publication 保留原生 Q4、完整 geometry identity、`{step,begin,end,position}` 及实际 mesh。求积共享同一 snapshot 寿命，不复制 scalar、不重拟合或更新材料库存。`points()` 每项包含 native cell、物理位置、严格正的 `physicalWeight` 和相标签；`phi<0` 为液相。`cells()` 是未改动权重累加的液气容积和误差诊断，不是新的材料状态。

每条 x/r 积分线从同一 Q4 的正式 MappedSegment 取得全部液气区间。短区间按独立 offset/length 嵌入，内层 y/z 用 GL4，外层自适应 GL8；RZ 权重含 `2*pi*r`。每个物理点必须通过原生 Q4 的严格相证书，未分辨相、不可表示内部点或非正/非有限权重立即拒绝。纯相用完整相证书；不能用中心或顶点判相。

拓扑图端点、中点与已括住的 x 投影转折位置只提供数值分段，不能替代原 scalar 或承担判相。控制原单元缩放坐标的各相 P4 矩；非轴 RZ 另控制乘 `(r_c/r)`、`(r_c/r)^2` 的矩，覆盖 H06 局部 P2 hoop 耗散。轴邻域沿用 H06 奇偶空间。

GL8 与两个子区间 GL8 的差是误差**估计**，不是严格上界。初始分段按宽度分配局部预算；接受叶段的估计还必须累计通过整单元检查。独立 Green 相容积检查阻止粗、细求积共同漏掉小相区；两侧正容积必须各自有正权节点。不得通过缩放权重、裁剪或拟合界面使检查通过。

默认 `relativeMomentTolerance=1e-9`、`relativeVolumeTolerance=1e-8`，均相对真实 `V_c`；合法范围为有限 `(0,1e-6]`。默认深度24、每单元8192次 slice 查询、整体200万节点；深度合法0–40，其余上限必须正。重试最多5轮，所有 slice 查询计入同一工作限额，超限或最终不满足精度明确拒绝。`report()` 返回实际点数、查询数、最大相容积差和累计矩估计。相容积差是同一静态几何的数值求积一致性，不是材料守恒要求。

## 黏性消费

`CchPhaseViscosity` 检查相同 geometry/full time/native mesh，允许 CCH 明确选中 publication 域的子集，且每个选中单元必须有覆盖。每点按真实相取 `muLiquid` 或 `muGas`，复用 H06 的同一 `CchStrainOperator` 和 `CchEnergyViscosity`，构造 `H_mu=2*mu_phase*dV`、`K_mu=E^T H_mu E`。不按 alpha 平均全单元点黏度；没有第二应变实现、第二散度或黏度 floor。

材料值必须有限且非负。正权与非负黏度保证矩阵对称半正定；两相黏度为零时矩阵、载荷和耗散精确零。即使无黏仍检查来源、layout、mesh 和严格零径向轴值；近轴正则性来自 H06 奇偶重构。

prepared 对象及其副本共享不可变 data，保留几何、求积和 metrics 的寿命。坐标/拓扑改变后旧对象失效，须显式 rebuild。`apply/dissipation/validateFor` 检查来源；只读 `viscosity()` 提供 H06 矩阵以组装 H04 的速度块，消费该低层矩阵前调用者仍须核验当前来源。完整时步的候选/接受绑定归后续任务。

## 复现与范围

标准构建 `tools/build.cmd -Target quantum_cch_phase_viscosity_contract`，然后 `ctest --preset local -R "^cch_phase_viscosity_contract$" --output-on-failure`（CMake 在 PATH 时）。XY/RZ uniform/2:1、纯相、平面、圆/球、极小液气相、选中半域及非多项式几何均有检查。

两相耗散 MMS 使用实际 native 面值，参考由独立解析相区/高阶积分计算；末档阶1.828–2.053。它验证黏性算子收敛，不验证完整 PDE。[H08](H08_Cch_Material_Transport.md)已实现唯一中点材料流量和非零连续性缺陷；H09质量—动量一致性、完整边界/Hodge、时间推进和低黏EHD尚待实现。
