# H06 / native CCH strain 与能量黏性

状态：2026-10-03，**H06 PASS**；[任务书](../task/H06_Energy_Viscosity_TODO.md)、[报告](../task/H06_Energy_Viscosity_Report.md)。[H07b 分相求积/黏性](H07b_Cch_Phase_Quadrature_Viscosity.md)已 PASS，下一项 H09；完整流体边界归 H15。本项不代表 CCH App、完整 NS 或低黏耦合已通过。

## 使用与所有权

```cpp
#include "Modules/Hydrodynamics/CchEnergyViscosity.h"
using namespace Quantum::Hydrodynamics;
auto points = CchStrainOperator::volumeQuadrature(mesh, physicalDivergence);
auto strain = CchStrainOperator::build(mesh, physicalDivergence, points);
auto viscous = CchEnergyViscosity::build(strain, std::vector<double>(points.size(), mu));
auto integralViscousLoad = viscous.apply(mesh, nativeValues); // Momentum uses its NEGATIVE.
double dissipation = viscous.dissipation(mesh, nativeValues);
```

prepared 对象共享不可变 data，复制不重建矩阵，不借 mesh 寿命。strain 保留 [H05](H05_Cch_Physical_Divergence.md) 同一 metrics/layout/D_I，矩阵只有 Nu 列；H06 不重新装配散度/压力载荷，也没有全局高阶速度未知量。当前原生速度仍只有 Nu，局部 P2 系数在准备时消元进 sparse E。

`CchStrainPoint` 包含 native cell DOF、物理坐标和物理 dV 权重；build 复制点列。内置 GL2×GL2 是单相体求积：XY 每点 V_c/4，RZ 每点 (V_c/4)(r_q/r_c)，真实 V_c 来自 retained Core domain。体积、r/z 一阶矩均有独立验证。零权重点可用于边界极限查询；不作为体积分贡献。

H06 检查 native 来源、点位于单元、非负有限权重及材料值，**不验证 Hermite geometry/time/material identity、切分覆盖或相求积精度**。[H07a](H07a_Cch_Phase_Geometry.md)已提供同源 cell/face 测度和瞬时率，[H07b](H07b_Cch_Phase_Quadrature_Viscosity.md)已用唯一几何 owner 构造/核验真实 phase split 体积样本及其身份，并以typed CCH对象接入本算子；本H06低层port自身仍不判相。材料采用[完全非守恒路线](../task/H07_Pure_Gals_Hydro_Coupling_Design.md)，不以独立库存约束几何。

## 实际应变

每个求积点的行顺序：

| 测度 | 行 |
| --- | --- |
| XY | ∂x u_x，∂y u_y，(∂y u_x+∂x u_y)/sqrt(2) |
| RZ | ∂r u_r，∂z u_z，u_r/r，(∂z u_r+∂r u_z)/sqrt(2) |

支撑通过同一 selected native face 邻接走两层；不从空间距离跨孔洞/不连通域借值。缩放 P2 重构使用 QR；非轴 patch 的原生 normal faces 作精确插值约束，SVD 只生成局部约束 nullspace，不显式逆全局矩阵。规则网格 normal strain 在单元中心等于原两面的中心差分；其它导数及 AMR patch 通过同一局部多项式得到，非多项式二阶由实际求积验证。

触轴支撑复用 H05 的奇偶反射，基底原点移到 r=0，并显式采用径向奇 P2={r,rz}、轴向偶 P2={1,z,r²,z²}。在 r=0 令 hoop=∂r u_r，剪切为正则零极限；无需计算 u_r/0，也不靠大数抵消满足轴条件。所有 apply/dissipation 都要求 native axis u_r 精确为零。

支撑不足 P2 rank 时明确拒绝，不能降阶或自动扩大到全域；新几何/AMR/layout 要重建 prepared operator。独立等价 mesh 可通过同一 Core identity 检验。

## 能量关系与无黏极限

黏度数组按同一 point 顺序复制；每点各 strain 行权重都为 2μ_q w_q。K_mu 只由 E^T diag(weights) E 形成，权重非负，因此 u^T K_mu u≥0。`dissipation` 从同一 E 的实际应变计算非负加权平方，使用 long double 累加并检查返回有限；它不修改物理状态。

全部有效权重为零时 K_mu 是压缩的空稀疏矩阵，apply 和 dissipation 精确为零。dissipation 在验证来源/轴值后直接返回零，不先计算可能溢出的 E u。正黏度没有 floor；已测缩放到 1e-300。非法材料、权重溢出、非有限组装/应用或耗散均明确拒绝。

没有添加 no-slip tangential penalty、open traction 或稳定化。H15 必须在同一能量双线性型里完成物理边界，且罚项随 μ 消失；此处 μ=0 是算子性质，不能据此证明不可压 Euler 时间推进稳定。

## 验证范围

21 个实际网格/selected 域场景、55 个多项式/刚体场、6 个小矩阵谱、105 个黏度缩放、12 个 MMS 点、3 个轴半径、6 个 H04 单体消费和25项拒绝均通过。三档均匀/2:1 XY/RZ 的末档 strain 阶2.059–2.356、积分耗散阶1.811–1.882。6个已测谱仅保留XY三种刚体/约束后RZ轴向平移零模；未把它推广为任意复杂域的通用谱证明。
