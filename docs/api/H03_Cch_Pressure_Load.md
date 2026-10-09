# H03 / CCH 负转置压力载荷

状态：2026-10-03，**H03 PASS**；[任务书](../task/H03_Negative_Transpose_Pressure_Load_TODO.md)、[验收报告](../task/H03_Negative_Transpose_Pressure_Load_Report.md)。本接口实现冻结 CCH §6.2 的积分压力载荷；[H04](H04_Cch_Single_Phase_Saddle.md)已消费它并完成封闭域单相参考求解。[H05](H05_Cch_Physical_Divergence.md)同源物理散度/轴规则已通过，[H06 能量型黏性](H06_Cch_Energy_Viscosity.md)已 PASS，[H07b 分相求积/黏性](H07b_Cch_Phase_Quadrature_Viscosity.md)已 PASS，下一项 H09；尚无可运行 CCH 流体 App。

## 入口和来源

公开头为 `include/Modules/Hydrodynamics/CchPressureLoad.h`，类型位于 `Quantum::Hydrodynamics`，消费者链接 [lamp_hydrodynamics_cch](H02a_Cch_Build_Boundary.md)。

```cpp
const auto layout = CchDofLayout::build(mesh, domain);
const auto divergence = CchIntegralDivergence::build(mesh, layout);
const auto pressureLoad = CchPressureLoad::build(divergence);
const CchNativeValues values(layout, nativeNormalVelocity, cellPressure);
const auto load = pressureLoad.apply(mesh, values);
const auto& C = pressureLoad.matrix();
const auto& D = pressureLoad.divergence().matrix();
```

build 只接受已准备的权威 [CchIntegralDivergence](H02b_Cch_Integral_Divergence.md)，保留同一不可变 handle，以 `-divergence.matrix().transpose()` 直接生成 compressed C。它不查询实际面、不独立装配 pressure gradient。稀疏 C 是派生代数矩阵；唯一 D_I 和布局仍由 retained divergence 持有。

| 接口 | 契约 |
| --- | --- |
| build(divergence) | 从同一 D_I 生成 C，不接收其他 gradient/geometry 数据 |
| divergence() | 返回 retained 权威散度的 const 引用，供后续 solver 消费同源 D_I/C |
| layout() | 返回上述散度的同一 CchDofLayout |
| matrix() | 返回 Nu×Np 的 const sparse matrix，nnz 与 D_I 相同 |
| apply(mesh, values) | 验证 actual mesh 与 values 的 layout identity，使用 C 乘 pressure，返回 Nu 个有限积分载荷 |

复制共享只读 C 和来源。没有借用 mesh/domain 或原 divergence 对象的生命周期；原对象析构后仍可读取矩阵，apply 时提供内容完全等价的实际 mesh。matrix() 的低层代数消费不检查字段来源，携带状态计算使用 apply。

## 压力功与边界语义

\[
C=-D_I^T,\qquad u^TCp=-(D_Iu)^Tp.
\]

每项系数精确为同一 D_I 的负转置，沿用 H01 的正坐标方向。内部面 `(Cp)_f=A_f(p_positive-p_negative)`；2:1 接口沿用实际 fine subfaces，不插值 virtual child 或修改 C。RZ full 2πr 测度直接继承，原零轴列变为 C 的精确零行，原 velocity DOF 保留。

本接口的 D_I/C 均为边界消元之前的完整算子。selected 边界（包括未选邻居）只有一项负转置贡献；常压只在内部面严格抵消，不是这个完整矩阵的全局零模。[H04](H04_Cch_Single_Phase_Saddle.md)在 retained 完整算子之上消元封闭域给定 normal velocity，并按压力连通域施加真体积零均值。RZ 物理正则和完整外部 pressure traction 等仍归后续任务，不能提前删掉边界载荷来制造零模。

返回值是 pressure×area 的积分载荷：XY 为单位厚度载荷，RZ 为全旋转体载荷。它不乘 dt、不除密度/dual volume，不是压力加速度；该转换及粗细接口高阶精度由后续 compatible Hodge 处理。

## 拒绝与验证范围

不同 layout identity、stale mesh、XY/RZ 或 selected 域来源不符均由同源 values/layout 检查拒绝。CchNativeValues 继续拒绝 NaN/Inf；有限 pressure 导致 C·p 溢出时抛出 `[CchPressureLoad] nonfinite integrated pressure load`，保留不可变算子和调用者值。velocity 改变不影响结果。

75 场景、2,400 组速度/压力的最大 η_adj 为 **7.38771e-16**，低于 H-R2 的 1e-12 和本验证的 1e-13；独立端点/有向图面载荷缺陷最大 2.25179e-15。25 个轴速度模式具有精确零散度/压力功，32 个拒绝和来源/寿命检查通过。此验收不证明 pressure/velocity MMS、Hodge、物理 AMR transfer 或低黏流体稳定性，Gate H 继续待验收。
