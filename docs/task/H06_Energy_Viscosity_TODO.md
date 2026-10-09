# H06 / native CCH 能量型黏性算子

**状态：PASS，2026-10-03；验收见[报告](H06_Energy_Viscosity_Report.md)。** 前置 H01–H05 PASS。依据[路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md) H06、冻结 CCH §8/§12.1 与 H-R5。下一项 H07。

## 目标

只解决一个数学事实：从实际 native 面速度构造 XY/RZ 应变 E_h，在非负物理求积下形成 K_mu=E_h^T H_mu E_h，使黏性耗散非负且 mu=0 精确消失。不是由调用者提供任意 E 的代数包装。

## 允许修改

- `include/Modules/Hydrodynamics/CchEnergyViscosity.h`、`src/Modules/Hydrodynamics/CchEnergyViscosity.cpp`：prepared 局部应变和黏性，归独立 CCH owner。
- root `CMakeLists.txt`：显式 source 注册。
- `tests/contracts/hydrodynamics/cch_energy_viscosity_contract.cpp`、`tests/CMakeLists.txt`：短契约与注册。
- 本 TODO、API、报告、现行路线图/架构/索引；证据 `runs/baselines/H06_20261003`，临时工作 `runs/checks/H06/<时间戳>`。

预算：3 production/build + 2 tests。

## 禁止修改

H01–H05 实现/测试、Core/Solvers、旧 Hydro/Compatibility、Interface/Electricity/App、cases、冻结 PDF 和历史报告。禁止全局高阶速度 DOF、人工最小黏度、顶点平均判相、黏性 upwind、轴 epsilon、clipping、静默 fallback 或第二套 D_I/C。

H07 接入唯一 Hermite 相分割/材料身份，H15 完整 no-slip/free-slip/traction/open 边界与 AMR pressure Hodge；本项不得宣称这些已实现。先提供带 native 来源的实际物理求积端口，不假造相几何；未接 Hermite 的输入仅称 single-phase/operator 验证。

## 数学约束

- XY E=(∂x u_x,∂y u_y,(∂y u_x+∂x u_y)/sqrt(2))；RZ E=(∂r u_r,∂z u_z,u_r/r,(∂z u_r+∂r u_z)/sqrt(2))。
- P2={1,x,y,x²,xy,y²} 局部重构只作用于准备阶段的 native 面中心值；规则内部为对称中心模板、二阶一致，2:1 patch 为缩放局部多项式再现。无新 primary coefficient state。
- 支撑从同一 selected-domain 面邻接局部两层得到，不跨不连通域/孔洞凭距离借值；不足 P2 rank 明确拒绝，不降阶。非轴邻域约束重构精确保留本单元的原生 normal face 值；轴邻域显式限制到径向奇/轴向偶的正则多项式子空间，避免靠舍入抵消满足奇偶。
- 轴邻域复用 H05 奇偶反射；axis native u_r=0 精确校验。r=0 的 hoop 用 d_r u_r 正则极限，绝不除零。
- 求积点、物理 dV 权重来自调用者并复制绑定同一 layout；内置单相 GL2×GL2 使用真 V_c（RZ 权重含 2πr）。H_mu 行权重为 2 mu_q w_q≥0；mu 为非负有限数，不设 floor。cut phase quadrature 必须由后续唯一 Hermite owner 供应。
- K_mu 仅由同一 E 和 H_mu 形成，动量侧使用 -K_mu u。零黏度无非零矩阵项、输出精确零；任意 μ 缩放连续线性。

## 实现要求

- [x] TODO 先于源码；核验 H05 交付、冻结源码/输入/程序/PDF身份。
- [x] 不可变 prepared strain 与 viscosity，公开头自足；同一来源、复制/寿命/错来源拒绝，有限性与轴校验。
- [x] XY/RZ 轴上/离轴、均匀/2:1 实际面值 P2 应变再现；轴正则及 rigid translation/XY rotation 零耗散。
- [x] 对称、PSD、独立逐点耗散等式；6 个已测谱排除非刚体零能量模式（在明确轴约束后）。
- [x] μ=0、稀疏零权、1至1e-12及1e-300线性缩放；无人工残留。实际 K_mu 进入 H04 同源单体参考求解，原连续性/压力均值通过。
- [x] 非多项式 XY/RZ 应变及物理黏性耗散三档 n=4/8/16，均匀和实际2:1网格，误差递减、末档阶≥1.8；独立解析参考不复用 E。
- [x] 标准 build/CTest、架构与目录/前置/App回归、来源和交付审计；报告真实通过范围，关闭 H06 后下一项 H07。

## 验收

**T0：** local Release 无新增 warning，architecture/self_test/repository_layout/lamp_help 通过；CCH 仍只依赖 Core/Eigen，cpp 只有一个 CMake owner。

**T1：** P2 应变归一化误差≤1e-11，刚体耗散≤1e-11；对称相对缺陷≤1e-12、u^T K u 与独立 2mu D:D 求积归一化差≤1e-11、最小特征值≥-1e-11*max(1,谱半径)，且没有额外数值零模。零黏度矩阵 nonzeros=0、apply/energy 精确零；缩放误差≤1e-11（含极小正黏度）。错布局/mesh/轴值、非法点/权重/mu/长度、rank不足、构造或应用溢出明确拒绝，失败不改变已有对象。

**T2：** XY u=sin(pi x)cos(pi y), v=cos(pi x)sin(pi y)；RZ u_r=r sin(pi z), u_z=cos(pi r)sin(pi z)。准确导数/hoop独立计算，在实际求积点使用物理 dV 加权 L2；单相 mu=1 的积分耗散另用独立高阶求积参照。4/8/16均匀及2:1序列均须递减且末档阶≥1.8。完整黏性 NS/PDE 边界与时间精度不是本项算子级验收，不能把该项 PASS 当作 Gate H。

## 失败规则

应变一致性、轴极限、PSD/额外零模、无黏或收敛失败即 FAIL，在允许范围修复，不调 gate。若需改前置方法、相 owner、完整边界或超过预算，先拆缺失 contract 并保持 H06 未关闭；实际外部依赖无法解决时才 BLOCKED。局部输入不支持须显式拒绝，不以只测通用矩阵恒等式替代 strain 验证。

## 交付

生产 strain/K_mu、短契约、API/TODO/报告、原始数值/阶/回归/哈希证据；列明 Hermite 相切分、完整边界、CCH App、NS 时间推进、低黏/EHD 的未完成范围。总目标保持 active。
