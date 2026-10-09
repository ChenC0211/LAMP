# H04 / 单相封闭域 native CCH 鞍点参考求解

**状态：PASS，2026-10-03；验收见[报告](H04_Single_Phase_Saddle_Report.md)。** 下一单任务 H05。

依据：[现行路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md) H04、H-R1/H-R2、冻结 CCH §6.4/§11.1。H01/H02/H03 已 PASS。本项只解决一个数学事实：同源 D_I/C 和 SPD 速度块进入一次单体鞍点求解，压力与速度同时满足动量和不可压缩约束。

## 目标

实现封闭域 prepared reference saddle solver，显式消元给定 normal velocity，按每个连通域采用真实 cell-volume 加权零均值压力；以独立解析单相 reaction/incompressibility MMS 验证压力与速度空间精度。保留接口供后续材料/Hodge/黏性提供 A_u，不新造全局高阶速度状态。

## 允许修改

- `include/Modules/Hydrodynamics/CchSinglePhaseSaddle.h`、`src/Modules/Hydrodynamics/CchSinglePhaseSaddle.cpp`：同布局的 SPD velocity block / RHS、prepared 组装与 reference solve、只读 result/report。
- root `CMakeLists.txt`：加入 CCH 唯一 source list。
- `tests/contracts/hydrodynamics/cch_single_phase_saddle_contract.cpp`、`tests/CMakeLists.txt`：短 algebra/MMS/拒绝契约与注册。
- TODO/API/报告、当前索引与路线图；证据归 `runs/baselines/H04_20261003`，临时检查归 `runs/checks/H04/<时间戳>`。

预算：3 production/build + 2 tests。

## 禁止修改

H01/H02/H03、Core/Solvers、旧 Hydro/Compatibility、Interface/Electricity/Plasma/App、case/字典、冻结方法、现有容限。不得装配另一套 D_I/C，不引入任意 cell gauge pin、压力稳定化、projection 物理算法、人工黏度或 fallback。H05 RZ 精度/正则、H06 黏性、H07–H14 物理耦合/输运/时间层及 H15 AMR Hodge/完整 production 边界不在本项实现。

## 数学约束

\[
\begin{bmatrix}A_u&C\\C^T&0\end{bmatrix}
\begin{bmatrix}u\\p\end{bmatrix}
=\begin{bmatrix}r_u\\0\end{bmatrix},\qquad C^T=-D_I.
\]

- A_u 是 caller 提供的同布局有限 SPD 积分速度块，封装复制并验证，不推断材料、dt 或 Hodge；不显式 inverse。
- buildClosed 显式选择封闭域：边界面固定 normal velocity，RZ 零轴行必须为零。内部面保留唯一 native velocity；known velocity 消元同时进入 free momentum RHS 和 continuity RHS。
- 每个 C_f 连通压力域施加 Σ V_c p_c=0，使用 retained Core 真测度；零均值通过约束 multiplier 增广，不 pin cell，不作为放宽 continuity 的 source。
- complete matrix 保留 C_f 和 C_f^T 的同源关系。使用可靠稀疏 direct reference solve，不返回 projection 中间状态。
- 求解后必须核对实际原生 free momentum、全部 cell D_Iu、每个域的压力均值和增广完整残差；代数残差与物理 closure 分开具名。
- 错来源/不兼容边界净流量/非法 SPD/非有限结果/物理残差失败明确拒绝，prepared source 和调用者输入保持。

## 实现要求

- [x] 开工冻结源码/测试/程序/PDF/case 身份，核对 H03 交付。
- [x] 同源 velocity block/RHS 与独立 residual config；prepared 只读 CCH 来源、边界编号、pressure components、增广稀疏矩阵和 reference factorization。
- [x] 正确 normal velocity 消元、component-wise volume mean 与边界流量 compatibility；RZ 零轴固定，不假造面积。
- [x] 只返回有限、通过全部残差的不可变 native values/report；零输入、重复 RHS、错误后恢复及独立等价 mesh 生命周期通过。
- [x] 独立实际边图构造 divergence-free algebra 场，验证 XY/RZ、AMR/不连通域、非对角 SPD 和非零兼容 prescribed velocity；这不声称 AMR/RZ PDE 收敛。
- [x] 均匀 XY n=4/8/16 的独立连续 reaction + pressure MMS，通过 pressure/velocity 全域物理加权 L2 收敛；不由离散矩阵反造 MMS RHS。
- [x] 标准构建、前置/架构/App 定向回归、来源/文档/交付审计完成，裁决 H04，下一项 H05。

## 验收

**T0：** local Release 无新增 compiler warning，公开头首 include；architecture/self_test、repository_layout、lamp_help PASS，新来源唯一 CCH owner，无旧 Hydro/Compatibility 依赖。

**T1：** A_u 尺寸/有限性/相对 symmetry≤1e-12、LLT SPD 检查；同源 C 与其转置、边界消元、连通域均值、只读/寿命正确。默认 frozen config：algebraic / free momentum / cell continuity / pressure mean 分别≤1e-11；component boundary net flow / gross prescribed flow≤1e-12。零尺度精确零，不以增广 multiplier 掩盖原 D_Iu。独立 algebra solution 速度/pressure 归一化误差≤1e-9、零边界压力零功/速度能量关系≤1e-10。

**T2：** 独立连续方程 a u+∇p=f、div u=0，a=2.5；ψ=sin(πx)sin(πy)，p=cos(πx)cos(2πy)。准确 native face averages、cell-average pressure、均匀 XY 真 dual volume 和解析 face-average forcing（经同一 uniform dual-volume Hodge 投影），n=4/8/16；全部档位 pressure/velocity 加权 L2 递减，末档各阶≥1.8。不把短三档扫描称作论文完整 benchmark 或时间收敛。H01/H02/H03、旧 gradient load、FirstJet App 保持；旧源码、程序、case/PDF 哈希不变。

## 失败规则

SPD、鞍点、边界 compatibility、全 cell 连续性、压力均值、MMS 阶或来源/寿命失败即 FAIL，范围内修复，不调容限或 MMS 参数。若需更改前置/Core/Hodge 或超过预算，拆出缺失 contract 后记 BLOCKED；不得以完成 algebra recovery 代替连续 MMS。

## 交付

prepared 单相封闭域 reference solver、独立短契约、API/TODO/报告、标准构建/回归/MMS/哈希证据。H04 不代表 AMR gradient MMS、RZ 物理精度、open traction/backflow、完整无黏流体推进或 Gate H；这些仍须按路线图验收。
