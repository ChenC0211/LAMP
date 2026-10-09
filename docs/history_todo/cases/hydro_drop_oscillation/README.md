此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_drop_oscillation/README.md。

运行入口：tools\run_case.cmd hydro_drop_oscillation 。结果：runs/hydro_drop_oscillation/<时间戳_标识>/。程序：LAMP。

# B3：无黏小振幅液滴振荡

**状态（2026-09-29）：V7 接受 3 步至 t=0.003 后，按用户要求优先处理耗时，已主动停止，完整振荡未验收。** 第三步用 12 轮达到原方程残差 6.73648e−11，液量漂移为 0；这仅为原 1000 步中的前 3 步。停止发生在第四步尚未接受时，进程实际退出 −1，属于性能调度中止，不是数值发散或到达终时；见[停止记录](output/20260929_170341_415_V7_normal_shape_44e95f/performance_priority_stop.json)。V1/V2 首步自由 H² 超限；V3–V6 为诊断/方法比较主动停止，均接受 0 步。V7 使用原冻结程序，不包含后续材料组装与历史初值改进。

[V7](output/20260929_170341_415_V7_normal_shape_44e95f/run.json)按法向形变构造毛细迭代主部，保持同一加速尺度与原完整方程验收。仍为原 64×64、密度比 1000、零黏度、dt=0.001、T=1，五份字典哈希一致，OpenMP 4 线程、BLAS 单线程；真实状态以运行记录和接受态历程为准。方法、相关控制及限制见[流体 §81](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#normal-shape-capillary-iteration)。非线性候选改善不能代替物理频率、阻尼与能量验收。

**耗时处理（2026-09-29）：** [流体 §85](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#shared-complete-flux-cost)直接重放首步真实候选；两组半步＋整步界面准备耗时减少约 26%/30%，几何逐字节一致，最终六步局部耦合报告与前版一致。本轮截获运行 `20260929_191834_319_Perf_candidate_capture_56034d` 在取得非零候选后主动停止，退出 −1、接受 0 步，不与 V7 接受态拼接。原字典不变，完整振荡和原规模整步效率仍待验收。

**当前性能诊断（流体 §90）：** `output/20260929_215857_147_Perf_assembly8_a40c2d` 保持原五份字典和 8 线程，已接受至 t=0.002；第 2 步用 9 轮迭代、10 次完整候选评价，主要计时阶段累计 685.72 秒。第 3 步取得性能问题证据后主动停止，退出 −1，属于诊断中止，原 T=1 未完成。根据用户要求先分析原因，不继续零散改法或新增长跑。详见[时间步耗时诊断报告](../../docs/QUANTUM_TimeLoop_Performance_Diagnosis_20260929.md)。

## 目的与物理模型

二维不可压缩两相 NS：ρ内=1、ρ外=0.001，μ内=μ外=0，σ=1，无重力、电场或供排。参考 [Basilisk oscillation.c](https://basilisk.fr/src/test/oscillation.c) 的 D=0.2、5% 二阶角向扰动、T=1 问题，评价毛细回复、频率及数值阻尼，不加入人工黏度。

## 初始与边界条件

初速为零，负相是滴内，初始曲线 r(θ)=R[1+0.05cos(2θ)]，R=0.1。字典采用 θ=2πs 的坐标公式生成唯一 Hermite 曲线，使用实际数值曲率；不以椭圆替换。解析初始面积 πR²(1+0.05²/2) 与实际几何面积的差计为初始化误差。

参考使用 [0,0.5]² 象限；本例为其完整镜像域 [−0.5,0.5]²，中心 (0,0)，四边不可穿透、零切向牵引。外壁距滴心仍为 0.5，没有压力出口。解保持镜像对称时与象限问题对应；本程序不额外钳制内部对称面，验收需报告对称性误差，不能把域宽 0.5 的中心整滴称作同一个边界条件。

闭域压力采用零均值规范，初始 pressure_available=0；后续压力属于中点，速度与几何属于步末。

## 离散与工况

完整 C1 速度、共同中点、同源相/质量/动量通量、中心对流和数值曲率牵引。64×64 均匀网格、h=1/64、D/h=12.8，64 个初始参数弧段，动态 AMR 关闭。dt=0.001、T=1，共 1000 步，每 10 步写场，能量与分相历史逐步保存。

空间比较用 128×128、128 个初始弧段；代表档另用 dt/2，不改变零黏度和密度比。几何初始误差、形状误差及 FV 库存漂移分别统计。

## 参考和评价

无限域线性理论为 ω²=(n³−n)σ/[(ρ内+ρ外)R³]，n=2；形状周期为 2π/ω，动能主频为 2ω。有限振幅和有限域的偏差不一概归为离散误差。

2026-09-29 已接离线诊断，直接读取每个保存时刻的 Hermite 多项式。在圆心 (0,0) 周围验证规则、单值的极坐标表示，再以 1024 个等角射线求半径 r(θ)，计算 a₂=∫r cos(2θ)dθ/(πR)、b₂=∫r sin(2θ)dθ/(πR)。与隔点 512 角结果对照，单列采样差。失去该表示时明确标记模态不可用，不排序轮廓点来伪造极坐标界面。已有 shape_moment 是面积矩不对称度，不作为 a₂ 的替代。

上下/左右镜像误差为相应极径差的均方根除以 R。a₂ 同方向过零间隔提供独立形状周期，过零时间按相邻输出帧线性插值，另列区间宽度及相对于线性理论的相位滞后。中心漂移不会被重置圆心消去，亦会反映在此形状/镜像观察中。

动能沿用参考的零初相模型 K=A exp(−bt)[1−cos(2ωt)]，在全部接受步时窗上拟合；输出频率、衰减率、相对 RMS 残差及搜索是否触边。A 对给定 b/ω 做最小二乘，频率搜索范围为 [0.5,1.5]ω₀，b∈[−2/T,10/T]，允许能量增长而非强制衰减；触边或残差过大不能把频率视为可信结论。模型与形状过零独立交叉判断，不只看一条拟合曲线。

本程序 kinetic_energy 为全域 1/2∫ρ|u|²dA；参考脚本输出象限 ∫ρ|u|²dA。因此严格镜像对称时，用本程序动能的 1/2 比较参考曲线，避免混淆域倍数与动能定义。原始能量误差直接取 K+σL−E₀，同时除以 E₀ 和初始变形能 E₀−2σ√(πA₀) 列出，后者避免误差被基础圆周表面能掩盖；离散方程余额另列，不从物理能量中扣掉质量时间项后宣称守恒。

在运行前固定以下**项目判据**，不作为文献统一容差：最细档频率相对线性理论误差 <5%，动能拟合相对 RMS <5% 且未触及搜索边界；每个线性形状周期的拟合振幅衰减 <5%，全时窗原始能量变化绝对值不超过初始变形能的 10%，保存帧镜像相对 L2 <10⁻³。保留空间加密和代表 dt/2 对照，有限振幅/有限外壁误差独立说明；单次摘要不自动宣布模块通过。形变时相对于圆的曲率差含真实物理变化，不能套用静滴的 1% 曲率门槛。

## 运行与结果

输入：[公式](../../../../cases/hydro_drop_oscillation/system/geometry.dict)、[物性与边界](../../../../cases/hydro_drop_oscillation/system/hydro.dict)、[界面](../../../../cases/hydro_drop_oscillation/system/interface.dict)、[网格](../../../../cases/hydro_drop_oscillation/system/mesh.dict)、[运行](../../../../cases/hydro_drop_oscillation/system/run.dict)。

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/run_case.ps1 -Case hydro_drop_oscillation -Threads 8
```

结果在 output/<时间戳>/，其 output/series/ 包含 hydro_moment_energy.csv、hydro_moment_phase_integrals.csv 和 hydro_moment_interface.csv；output/fields/ 保存物理场与数值曲率。[本次预检](output/20260928_230936_712_analytic_input_744364/run.log)没有推进物理时间，不能说明无黏稳定性或完整振荡已通过。

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/analyze_hydro_standard.ps1 -Run runs/hydro_drop_oscillation/<时间戳> -Benchmark oscillation
```

仅分析已终止的实际运行；失败或不足完整 T=1 时加 `-AcceptedPrefix`。输出位于该次运行的 `output/analysis/<时间戳_类型>/`：history.csv、frames.csv、geometry.csv、可用时的 kinetic_fit.csv/shape_crossings.csv，以及带程序、输入、数据、分析工具身份的 summary.json。前处理核对生效工况、初态误差、接受时钟及步中压力时层，缺少新诊断输出的旧结果明确拒绝。当前只完成主程序和后处理代码编译，尚无本版计算数据或拟合结果。
