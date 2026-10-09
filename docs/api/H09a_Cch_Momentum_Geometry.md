# H09a / 原生动量控制体与同源几何材料质量

状态：2026-10-03，H09a PASS，验收见[任务书](../task/H09a_Momentum_Dual_Geometry_TODO.md)与[报告](../task/H09a_Momentum_Dual_Geometry_Report.md)。沿用[完全非守恒材料路线](../task/H07_Pure_Gals_Hydro_Coupling_Design.md)。本接口只准备当前层几何/质量；[H09b](H09b_Cch_Dual_Material.md)已实现流量转移、真实非零dual缺陷及uR项并通过[验收](../task/H09b_Dual_Material_Consistency_Report.md)，两项共同关闭H09；下一项H10。

公开头：[CchMomentumGeometry.h](../../include/Modules/Hydrodynamics/CchMomentumGeometry.h)，唯一实现 owner 为 `lamp_hydrodynamics_cch`。依赖正式 Interface 的只读区域积分，不链接旧 Hydro/Compatibility。

```cpp
auto dual = CchMomentumGeometry::build(mesh, physicalDivergence);
auto mass = CchMomentumMass::build(mesh, dual,
    actualPhaseGeometry, rhoLiquid, rhoGas);
mass.validateFor(mesh, actualGeometryIdentity, actualTimeLayer);
const auto& Mu = mass.masses();
const auto& intersections = dual.volumes();
const auto& commonSides = dual.segments();
```

## 实际控制体

`volumes()` 按现有 velocity DOF 保存矩形范围、真实物理体积、一或两个实际 native cell 的交集，以及四边的 segment 索引。每个 cell 的某速度方向由该方向两侧 faces 各拥有半矩形；2:1 双子面沿切向分开。normal 深度取实际 cell 中心，不取压力算子中的虚拟 fine 中心。selected 边界只保留内侧半体积；孔洞/不连通域同样使用实际选择。

`overlaps` 的 `pressure` 是同一 layout 的压力 DOF，`cellFraction=V_fc/V_c` 是纯几何转移权重，不能解释为该交集的 phase fraction。RZ 使用真实 `pi*(r_hi+r_lo)*dr*dz`，因此径向两半通常不同。各方向交集精确分区同一 cell；控制体本身没有新的全局自由度。

`segments()` 每条共同边只存一次，以正坐标法向定向；`negativeVelocity/positiveVelocity` 为同速度方向的 owner，外域为−1。每个 dual 四边的端点精确覆盖。RZ 轴边面积为0，其动量半控制体仍有正体积；轴速度正则规则由 H05/后续方程消费，几何对象本身不设置速度。

每段具有唯一几何来源：`primalVelocity` 指向实际 primal face，或 `ownerPressure` 指向实际 cell 内部。它们为后续唯一材料流量转移提供来源，**本接口还没有生成 dual 材料流量**。共享边段不能另插值一套 rho*u。

`metrics()` 保留原 H05 metrics 和同一 D_I；没有第二 primal 散度。`sameIdentity()` 比较同一不可变几何 data；复制共享只读数据，重新构建产生新的几何 identity。

## 同源分相质量

`CchMomentumMass` 保留同一 `CchPhaseGeometry` publication 与 dual geometry。构建要求相同 layout identity、相同保留 D_I、实际 mesh 和 XY/RZ 测度。FirstJet 不含 M。

正式 `HermiteGalsPhaseRegionSource` 从 publication 的原 Q4 snapshot 裁剪完整 dual 矩形，返回各实际叶交集的分相体积。开放端点的 contacts/symmetry 声明从 retained curve metadata 去重读取；重新发现的是同一 scalar 的派生拓扑，不复制 snapshot、不拟合替代几何、不初始化库存。

`phaseVolumes()` 按 velocity DOF 分别保留液气体积；`overlapPhaseVolumes()` 的内层顺序对应几何 overlaps。纯相精确处理；小相独立测量后由容量求补相。`masses()` 为有限正 double，可由后续动量块消费：

\[
m_{geo,f}=\rho_l\,V_{l,f}+\rho_g\,V_{g,f}.
\]

两相密度有限正数，无 density floor/真空。质量由真实 cut 区域积分，通常不等于 `sum(V_fc/V_c * m_cell,c)`。保存的几何权重不能拿来平均密度并冒称 dual 材料积分。

`report()` 给出每速度方向的液气总体积/质量、最大逐 cell 相分区误差、区域 Green ledger 缺陷、相对 moment error estimate 与 cut 交集计数。每方向 dual 相积分应覆盖同层 cell 相积分；这里的误差门槛为1e−10，是空间积分自洽检查，不要求时间连续性缺陷为零，也不校正质量或体积。moment estimate 不是严格误差上界。

## 来源、拒绝与能力范围

`geometryIdentity()/timeLayer()` 原样返回 publication 来源；`validateFor` 比较当前实际来源、完整 `{step,begin,end,position}` 和 mesh。相同时刻不等于相同候选区间。错误 layout/D_I、密度、来源/时钟、变更 mesh 和正材料乘积下溢/溢出、不可表示的正动量质量明确拒绝；已持有对象不被失败构造修改。

对象不借 mesh 指针；等效实际 mesh 和原 mesh 析构后均可验证。实际 refine 后旧对象拒绝使用，重建新 publication/layout/D_I/dual/mass；不会把旧数组重新贴上新标签。本项未实现旧→新速度 transfer、兼容 AMR Hodge 或接受事务。

区域查询目前逐次执行正式 mesh 验证；本项没有新增 unchecked 访问或修改 Interface 性能接口。大规模批量裁剪的成本仍须在后续完整算例测量；算子空间阶不是完整流体 PDE 阶、时间阶或高密度比稳定性证明。

复现：`tools/build.cmd -Target quantum_cch_momentum_geometry_contract`，随后 `ctest --preset local -R "^cch_momentum_geometry_contract$" --output-on-failure`（CMake 在 PATH 时）。
