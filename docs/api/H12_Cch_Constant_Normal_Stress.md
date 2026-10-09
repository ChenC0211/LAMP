# H12 压力相容常量法向应力

`CchInterfaceNormalLoad` 属于 `lamp_hydrodynamics_cch`，生产头位于 `include/Modules/Hydrodynamics`，实现位于对应 `src`。这是只读载荷 owner，不是毛细时间步或接受状态。

```cpp
const auto pressure = CchPressureLoad::build(phase.metrics().divergence());
const auto normal = CchInterfaceNormalLoad::constant(mesh, pressure, phase, J0);
normal.validateFor(mesh, pressure, phase.geometryIdentity(), phase.timeLayer());
auto rhs = normal.balanceLoads();
for (double& value : rhs) value = -value; // A*u + C*p = rhs
```

## 唯一来源与符号

工厂保留同一 `CchPressureLoad` 和 `CchPhaseGeometry`，要求相同 layout identity 及同一 D_I owner。alpha 取保留 publication 的派生 `liquidFractions()`；仅构造临时 `J0*alpha` 并调用现有 `pressure.apply`，不重组梯度、法向或几何。返回 `F=C(J0 alpha)`，它是 `Cp+F=0` 中的平衡项，传给现有鞍点 RHS 时应减去 F。

J0 为给定的有限法向总应力，不读取曲率、电场或材料库存；几何、材料通量、R/uR、KEP 与 Interface 接受状态不更新。非有限 J0、溢出载荷、错误 layout/D_I、mesh/source/full clock 明确拒绝。验证还要求相同 C owner；按值复制共享不可变存储，重新装配 C 不视为原 owner。失败不改变已返回对象。

只读入口为 `pressureLoad()`、`geometry()`、`layout()`、`constantJump()`、`balanceLoads()`。完整 `{step,begin,end,position}` 沿用 phase 验证；regrid 后需实际 FirstJet transfer，再重建 metrics、phase、C 和载荷。不可重贴来源或仅按相同尺寸复用。Hydro 可选含洞、断开或离轴子域；标量支撑范围仍服从既有 FirstJet 契约。

## 静平衡与边界

全原生行在 `p=-J0 alpha` 时平衡。封闭域求解使用逐连通域真体积均值 gauge，解析参考为 `p=J0 mean_V(alpha)-J0 alpha`；自由面平衡，固定边界由壁面反力承担。全原生 C 在外边界一般不消灭常数，因此不能要求 gauge 解的全行反力为零。RZ 轴载荷和径向速度为零；固定边界反力功为零。

实际速度块由 `CchMomentumMass` 与 `CchPhaseViscosity` 生成 `A=M/dt+K_mu`；零黏时黏性矩阵严格消失。解析压力只用于事后检查，不替换 `solve`。dt 只定义静态质量块，不声明时间积分精度。

## 原系统的近零残差校正

H12 暴露原一次 SparseLU 解的局部连续性不足：力平衡已接近舍入误差，但部分单元通量尺度小至 1e-32，普通增广点积丢失所需残差。`CchSinglePhaseSaddle` 现在最多对同一增广矩阵、同一 LU 因子作 8 次残差校正，没有第二个物理算子或解析静滴分支。

`b-Kx` 使用乘积 FMA 误差及补偿求和。若原代数和动量门槛均已通过，校正的动量 RHS 取零，使 `A du+C dp=0` 保持已有力平衡；已通过的压力均值同样取零校正 RHS。每次仍解完整增广系统，并重检原代数、动量、逐单元连续性和真体积压力均值。容差与零尺度规则不变，超过次数或非有限值仍抛错。首次解通过则不校正；report 的 `linearRefinements` 仅记录结果。

实际验收与原失败见 [H12 报告](../task/H12_Constant_Normal_Stress_Report.md)，完整研究见 [benchmark](../../cases/benchmarks/cch_constant_normal_stress/README.md)。该接口尚不代表真实曲率毛细力、一般时变应力、切向力或完整流体 App 通过。
