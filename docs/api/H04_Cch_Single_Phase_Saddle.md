# H04 / 单相封闭域 CCH 鞍点参考求解

状态：2026-10-03，**H04 PASS**；[任务书](../task/H04_Single_Phase_Saddle_TODO.md)、[报告](../task/H04_Single_Phase_Saddle_Report.md)。已有模块级封闭域 reference solve；[H05](H05_Cch_Physical_Divergence.md)已验证 RZ 物理散度/轴规则，[H06 能量型黏性](H06_Cch_Energy_Viscosity.md)已 PASS，[H07b 分相求积/黏性](H07b_Cch_Phase_Quadrature_Viscosity.md)已 PASS；后续H08a/H11与[H12常应力](H12_Cch_Constant_Normal_Stress.md)已在规定范围验收，当前下一项为H13。AMR Hodge、完整 production 边界集和流体 App/时间推进尚未实现。

## 输入和入口

公开头 `include/Modules/Hydrodynamics/CchSinglePhaseSaddle.h`，类型位于 `Quantum::Hydrodynamics`；链接 [lamp_hydrodynamics_cch](H02a_Cch_Build_Boundary.md)。

```cpp
const auto velocityBlock = CchSaddleVelocityBlock::build(layout, integralAu);
const auto saddle = CchSinglePhaseSaddle::buildClosed(pressureLoad, velocityBlock);
const CchSaddleRhs rhs(layout, integratedMomentum, prescribedBoundaryVelocity);
const auto result = saddle.solve(mesh, rhs);
const auto& values = result.values();
const auto& report = result.report();
```

CchSaddleVelocityBlock 复制有限 Nu×Nu 的 A_u，绑定同一 layout，检查相对 symmetry≤1e-12、LLT 正定和有限 factor。A_u 是后续物理 owner 提供的积分速度块；本接口不推断密度/dt、不生成 Hodge/黏性，也不替换系数。修改材料或 A_u 时重新构造 block 和 prepared saddle，不以相同尺寸复用旧 factor。

CchSaddleRhs 复制两个长度 Nu 的有限数组：integralMomentum 为 r_u，prescribedNormalVelocity 只给边界 DOF；free face 的该 payload 必须精确为零。RHS 绑定 layout identity，solve 验证实际 mesh。结果返回只读 CchNativeValues 和只读 CchSaddleReport，没有主时钟、commit 或另一个物理 state owner。

## 封闭域消元与压力均值

buildClosed 显式选用所有边界 normal velocity 已知的封闭压力问题。以同源 C 的两项行识别内部 free face，一项行识别 selected 边界，零行识别 RZ 轴并固定零速度；不独立构造 gradient 或 geometry。所选域内的洞/未选邻居同样成为给定速度边界。给定边界可非零，但每个压力连通域的净流量必须兼容。

令 f/b 为 free/prescribed velocity，G 的每列为一个连通域内的真实 cell volumes，则一次求解

\[
\begin{bmatrix}A_{ff}&C_f&0\\C_f^T&0&G\\0&G^T&0\end{bmatrix}
\begin{bmatrix}u_f\\p\\\lambda\end{bmatrix}
=\begin{bmatrix}r_f-A_{fb}u_b\\-C_b^Tu_b\\0\end{bmatrix}.
\]

G^T p=0 是逐连通域 ΣV_c p_c=0，无任意 cell gauge pin。C_f 与 C_f^T 直接来自 retained [H03](H03_Cch_Pressure_Load.md)；该来源仍保留同一 H02 D_I。λ 是临时压力均值约束 multiplier，不是物理源或 primary state。返回恢复后的 Nu+Np 原生数据；检查实际 D_Iu，不能用 λ 吸收不兼容边界流量。

固定速度处的动量方程由边界反力承担，本接口只验收 free momentum。pressure/traction anchor、free-slip 切向条件、open backflow 等完整边界归 H15；不能将 buildClosed 自动解释为这些模型。

## Prepared 数据与线性后端

只读端口包括 pressureLoad()/velocityBlock()/layout()、增广 matrix()、freeVelocityDofs()、pressureComponents()/pressureGaugeCount()。复制共享 prepared 数据和 factor，原 mesh/domain/input 对象析构后可提供独立等价 mesh 继续求解。

增广矩阵 compressed / column-major，采用 [Eigen SparseLU](https://libeigen.gitlab.io/eigen/docs-3.4/classEigen_1_1SparseLU.html) 作冻结 CCH §11.1 允许的 direct reference 后端，分解在 buildClosed 完成，后续 RHS 使用同一 factor。没有显式 inverse、投影物理步骤、迭代 fallback 或隐式压力修正。这是开发参考求解，不是大规模 MINRES/AMG 性能交付。

## 残差与拒绝

| config / report | 定义与默认 gate |
| --- | --- |
| algebraicRelativeTolerance / algebraicRelativeResidual | 增广完整 ||Kx-b||₂/||b||₂≤1e-11 |
| momentumRelativeTolerance / freeMomentumDefect | free momentum 的 L2 缺陷 / 各行 |A_u u|+|Cp|+|r_u| 的 L2 尺度≤1e-11 |
| continuityRelativeTolerance / maxCellContinuityDefect | max_c |(D_Iu)_c|/Σ_f |A_fu_f|≤1e-11 |
| pressureMeanRelativeTolerance / maxPressureMeanDefect | 每个连通域 |ΣV_c p_c|/Σ|V_c p_c|≤1e-11 |
| boundaryFluxRelativeTolerance / maxBoundaryFluxDefect | 每个连通域 |ΣsA_fu_b|/Σ|A_fu_b|≤1e-12 |

零尺度要求缺陷精确零。各 tolerance 有限且处于 (0,1)；代数残差与物理 closure 分开校验。report 另记录 freeVelocityCount、pressureGaugeCount、maxAbsGaugeMultiplier 和 linearRefinements；它们不决定数值分支。

2026-10-06 [H12](H12_Cch_Constant_Normal_Stress.md)补充同一增广矩阵/同一 LU 因子的有界线性残差校正（最多8次）：补偿计算 b-Kx，已满足门槛的力平衡和压力均值通过零校正 RHS 保持。原矩阵、gauge、全部残差门槛及零尺度规则不变，首次直接解通过则不校正；失败仍拒绝，不返回半成品。

错 layout/stale mesh、非法矩阵/配置、free prescribed payload、非零轴速度明确拒绝；component flux 不兼容、非有限解/尺度或任一 residual 未通过抛出明确异常，不返回半成品 result。完整零 RHS 精确零，失败后的合法调用及已返回结果保持。

39 个独立 algebra 场和 21 个拒绝通过；均匀 XY 的三档连续 reaction/incompressibility MMS 实测速度阶2.003、pressure阶2.026。H04 本身对RZ/AMR只提供代数/真测度gauge测试；[H05](H05_Cch_Physical_Divergence.md)单独验证物理散度/轴规则，不将归一化算子反转置为压力载荷。RZ pressure/velocity、AMR Hodge和完整边界仍待对应验收。
