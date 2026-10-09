# H05 / 真测度物理散度与轴边界规则

状态：2026-10-03，**H05 PASS**；[任务书](../task/H05_Rz_True_Measure_Axis_Regularity_TODO.md)、[验收报告](../task/H05_Rz_True_Measure_Axis_Regularity_Report.md)。[H06 能量型黏性](H06_Cch_Energy_Viscosity.md)已 PASS，[H07b 分相求积/黏性](H07b_Cch_Phase_Quadrature_Viscosity.md)已 PASS，下一项 H09。当前无 CCH App/完整流体时间推进；AMR pressure Hodge 和完整 production 边界仍归 H15。

## 入口与来源

公开头 `include/Modules/Hydrodynamics/CchPhysicalDivergence.h`；类型位于 `Quantum::Hydrodynamics`，消费者链接 [lamp_hydrodynamics_cch](H02a_Cch_Build_Boundary.md)。

```cpp
const auto physical = CchPhysicalDivergence::build(mesh, integralDivergence);
const auto cellAverageDivergence = physical.apply(mesh, nativeValues);
```

build 验证实际 mesh，保留同一 [H02 D_I](H02b_Cch_Integral_Divergence.md)/layout，按原生 velocity 顺序准备 Core face areas 与轴 DOF。cellVolumes()直接引用 retained domain 的真实 V_c；faceMeasures()/axisVelocityDofs()/divergence()/layout() 均只读。复制共享 prepared 数据，不借用 mesh/domain；重新建 mesh/material 来源不能只按尺寸复用。

apply 拒绝不同 layout 和非零径向轴面速度，再调用 retained D_I.apply 验证实际 mesh 并计算唯一积分通量，最后逐 cell 除 V_c。结果是 cell-average physical divergence，单位为速度/长度；不是 integrated net outflow、时间相通量或动量载荷。

RZ 采用 full 2π 旋转体测度：V_c=π(r_+²-r_-²)Δz，径向 A_f=2πr_fΔz，轴向 A_f=π(r_+²-r_-²)。XY 是单位厚度。没有额外归一化稀疏矩阵、虚拟粗面、经验1/r补项或另一个散度 stencil。

**压力载荷继续为 C=-D_I^T。不能把 volume-normalized divergence 的负转置当作新的 C。** [H04 鞍点](H04_Cch_Single_Phase_Saddle.md)继续消费原积分算子，源对象及矩阵没有改变。

## 轴正则与派生反射值

径向轴面保留原 native DOF，A_f精确零，但 apply 要求其 u_r精确零。任何非零值（包括极小非零）都拒绝，不改成零；积分 D_I 的零轴列不能作为正则证明。

```cpp
using Axis = CchDofLayout::NormalAxis;
const double radialGhost = physical.reflectAxisNormalVelocity(Axis::First, radialAtPositiveR);
const double axialGhost = physical.reflectAxisNormalVelocity(Axis::Second, axialAtPositiveR);
const double pressureGhost = physical.reflectAxisPressure(pressureAtPositiveR);
```

输入是 caller 在正半径的有限样本，返回相应负半径伙伴：径向速度奇，轴向速度/压力偶。轴向对称中心导数因而满足∂r u_z=0；radial odd反射与轴面u_r=0一致。这些是边界派生值，不存 ghost primary state，不构造另一套全局速度场。

只有实际选中 RZ 轴面的域可调用反射；XY、离轴域和未选中轴面的域拒绝。非法 component、非有限样本、归一化溢出也明确拒绝。任意失败不修改来源或已返回的数值数组；旧 mesh 的几何/AMR 改变后 apply/build 拒绝，独立等价 mesh 可继续使用。

本接口实施轴边界奇偶条件，不据此证明任意内部场光滑。轴上 u_r/r 的 strain 正则极限和黏性已由[H06](H06_Cch_Energy_Viscosity.md)实现；没有径向epsilon、人工黏度、速度clipping或fallback。

## 已测范围

72组XY/RZ轴上/离轴、n=4/8、depth=0/1/2、全域/左半域/不连通域/带孔域；216组独立解析通量，9组半径1/1e-3/1e-6的有限轴极限，216组独立正负半径反射，16项拒绝通过。解析体积平均散度的最大归一化误差3.33e-15。

均匀 RZ 三档非多项式 divergence MMS 末档阶2.00362，n=16的物理体积加权相对L2约0.2325%。本项不是RZ压力/速度PDE收敛、AMR pressure Hodge或完整流体稳定性验收。
