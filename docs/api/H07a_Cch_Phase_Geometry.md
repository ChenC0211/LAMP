# H07a / CCH 同源相测度与瞬时率

状态：2026-10-03，验收见[报告](../task/H07a_Same_Geometry_Phase_Rate_Report.md)、[任务书](../task/H07a_Same_Geometry_Phase_Rate_TODO.md)。材料输运采用用户授权的[完全非守恒路线](../task/H07_Pure_Gals_Hydro_Coupling_Design.md)。[H07b 分相体积求积/黏性](H07b_Cch_Phase_Quadrature_Viscosity.md)已 PASS，H07 关闭；下一项 H09。

公开头：[CchPhaseGeometry.h](../../include/Modules/Hydrodynamics/CchPhaseGeometry.h)，实现归 `lamp_hydrodynamics_cch`，正式只读依赖 `quantum_interface`。不链接旧 Hydro/Compatibility，不复制 geometry snapshot。

```cpp
auto phase = CchPhaseGeometry::build(mesh, metrics, publication);
auto reconstruction = CchAdvectorReconstruction::prepare(mesh, momentumGeometry);
auto actual = reconstruction.evaluate(mesh, nativeValues,
    phase.geometryIdentity(), phase.timeLayer());
auto rate = CchPhaseRate::prepareSpatial(mesh, phase, actual);
```

`metrics` 为 H05 保留的同一 layout/D_I/真物理测度；`publication` 必须为无 M 的 FirstJet。CCH 可选择 publication 支持域的子集；缺少选中单元的几何必须拒绝。XY/RZ 模式、真实单元容积必须一致。复制对象共享不可变数据，保留 publication/metrics 寿命，不借调用者 mesh 指针。

`cellVolumes()` 按 pressure DOF 返回液气几何容积；cut 使用原始同源 Green 测量，pure 使用 publication 已认证的精确容量/零值。raw pure Green 值仍由 publication 作为舍入误差 ledger 保留，不能误当 pure 容积 authority。`liquidFractions()` 为 V_l/V_c；它不成为独立库存。

`faceMeasures()` 按 velocity DOF 返回液气面积，直接积分已保留的两条真实 vertex–midpoint trace。`phi<0` 为液相；每个保留区间分别按独立 offset/length 积分，RZ 面含 `2*pi*r`。液气面积分别累计，绝不由两次舍入端点之差构造短区间或用全相减法抹掉小气相。零面积轴面两相面积精确零。闭合 gate 为 `128*epsilon(double)` 相对真实面面积，失败拒绝，不剪裁面积。

`CchPhaseRate::faceRates()` 为液气瞬时体积流率，正号沿 native 正坐标法向；当前由[同一actual B的空间相率工厂](H10j_Cch_Spatial_Phase_Rate.md)生成，ordinary/重合trace消费实际normal polynomial与causal侧。原face-constant `evaluate` 在H10l2退役，`isSpatial()`恒为true；调用方包含`CchSpatialAdvector.h`，momentumGeometry使用同一metrics。`totalFaceRates()` 是 A_f*u_f。`cellOutwardRates()` 只从保留 D_I 的非零项符号生成，没有第二份散度；相邻单元以相反符号消费同一面率，2:1 聚合使用真实 fine subface 求和。轴上径向速度必须严格零。

来源校验比较完整 identity 及 `{step,begin,end,position}`，只比较 time() 不足以绑定同一候选。调用者提供速度实际来源标签，不能把旧速度重新贴成当前标签；此 API 无法识别调用者伪造的 provenance。`validateFor` 也核验实际 mesh；等价内容 mesh 可用，任何坐标/拓扑变更后旧对象失效，须显式 transfer/rebuild。不同尺寸一致的 layout 也不能混用。

时间单位：XY 体积率为面积/时间（单位深度），RZ 为体积/时间。该对象是某一瞬时的 rate，没有 dt、时空积分或 accepted commit。共享面抵消不证明下一层 GALS 几何质量守恒；本项不生成 F_rho、材料质量更新、dual correction 或物理时步。[H08](H08_Cch_Material_Transport.md)已核算非零cell材料连续性缺陷，H09后续处理dual材料/动量一致性。

复现短检查：标准 `tools/build.cmd -Target quantum_cch_phase_geometry_contract`，然后 `ctest --preset local -R "^cch_phase_geometry_contract$" --output-on-failure`（CMake 已在 PATH 时）。完整运行仍使用统一 case 启动器；目前没有 CCH 生产流体 App。
