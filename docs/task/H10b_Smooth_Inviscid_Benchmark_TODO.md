# H10b / 光滑无黏对流完整周期研究

状态：**FAIL**，2026-10-03，见[报告](H10b_Smooth_Inviscid_Benchmark_Report.md)。28完整窗口/15,116步、独立时间研究和19/19回归已完成，空间一致性失败；H10未关闭，下一单任务H10c设计与修正。前置[H10a](H10a_Compatible_Center_Convection_Report.md)算子恒等式PASS保持；依照[路线图§0](../roadmap/05_LAMP_Module_Task_Roadmap.md)和[H10设计](H10_Compatible_Center_Convection_Design.md)，本轮只执行H10b。TODO先于研究源码，材料保持用户授权的完全非守恒路线。

**当前重跑（2026-10-05）：** [H10o](H10o_Bounded_Final_Native_Convection_Report.md)的实际有界闭合算子完成原28窗口，显式H10b-dissipation-v2的84项通过；原四累计donor损失阶仍FAIL，原H10b历史裁决不改。下一单任务H08a。

## 目标

只核验一个数学事实：H10a实际相容对流载荷在真实XY/RZ uniform/2:1网格上，是否准确输运光滑周期向量并保持无黏动能预算；功恒等式不能代替PDE收敛。

## 允许修改

- 新增 `cases/benchmarks/cch_convection/driver.cpp`、`benchmark.json`、`README.md`，修改 `cases/benchmarks/CMakeLists.txt` 注册及研究索引。
- 本TODO、研究报告、当前路线图/文档索引及必要的现行架构/H07耦合/H10设计状态。
- 证据 `runs/baselines/H10b_20261003`；临时核查位于start.json登记的 `runs/checks/H10b/<timestamp>`；通用启动器产生 `runs/benchmarks__cch_convection/<timestamp>`。

预算2 research/build源码文件、0测试文件；不改任何生产数值实现，不新增临时生产开关。

## 禁止修改

H01–H10a源码/旧tests、Interface/Core/Solvers/旧Hydro/Compatibility/App、既有算例输入、工具、冻结PDF和历史TODO/报告。不得注入解析解修正状态、人工黏度、滤波、投影或新材料库存。H11/H14/H15不在本项。

## 数学约束与能力边界

- 研究方程为固定散度为零的规定平移速度v下的无黏向量输运：`partial_tau w + v dot grad(w)=0`。XY：v=(1,0)，w=(cos(2π(x−tau)), sin(2π(x−tau)))；RZ：v=(0,1)，w=(r sin(2π(z−tau)), 0.4+cos(2π(z−tau)))，轴径向分量为0。完整周期T=1，mu=0、无界面力。解析参考独立于生产重构/散射。
- 三层真实纯FirstJet发布、实际相分割、H08材料通量、H09真实dual质量/共享流量和H10a算子一次构成**自治半离散算子实验**。所有apply始终验证这一不可变prepared owner的原clock/identity；研究时间tau不是生产candidate时钟，不能伪装实际Interface/NS接受事务。
- 规定advector v与被输运探针w明确区分；本研究检验线性对流算子，不能称完整非线性Euler/NS。H14/H15仍负责实际流体时间推进、压力/速度耦合、一般Hodge及生产边界。
- 周期端口由外部相反边界的同分量、同切向区间唯一配对，正负流量检查后用两owner的同一中心trace；零流量边界显式给值，RZ轴径向trace0。保持所有native自由度，包括边界半dual；不隐藏边界误差，不把研究trace当成生产周期mesh。
- RZ轴径向native值为H05规定的本质条件0，研究积分器只在自由分量推进，轴反力不做功；全部轴DOF仍计入误差，原始未约束轴载荷另列诊断。该几何轴约束不是解析解修正或生产H15边界实现。
- `dw/dtau=momentumLoads/m_mid` 只调用实际生产载荷，保留其中的uR；实际质量变化、R/G、边界功均诊断。研究RK4仅编排同一算子的四次调用，不成为生产时间算法。
- donor对照只属于研究driver，同一流量/质量/几何及周期配对；不链接为生产fallback。

## 实现要求

- [x] 冻结核验H10a交付/证据，记录源码、旧输入、PDF/App哈希和任务开始位置。
- [x] 单一研究driver，经标准CMake/通用run_case调用，输出完整时序、误差/阶、动能和原始门槛结果。
- [x] 四种测度/网格组合，n=8/16/32，中央[.25,.75]^2真实一级2:1；每组跑完整T=1、RK4步数32n，同窗donor对照。
- [x] 独立时间研究XY uniform n=8，20/40/80步对照640步的同算子终值；时间误差和解析空间误差分别报告。
- [x] 如周期输运暴露空间错误，增加只读affine载荷定位：解析外边界trace仅用于瞬时参考，不进入时间积分；逐native自由度输出实际/解析导数及coarse–fine/exterior标记，区分研究周期端口与内部装配问题。
- [x] Release构建，架构/目录/受影响H10a及前置契约回归，源/输入/证据审计归档。

## 验收 T0/T1/T2/T3

**T0：** 标准Release构建无新增warning/error；原生产/App/PDF/旧输入不变；研究只链接正式CCH及其已允许依赖，默认CTest无长研究。

**T1：** 真实周期配对唯一、对置面积/流量归一化误差≤1e−12；每次实际apply的boundary KE功/材料功归一化≤1e−12，质量不变和R/G归一化≤1e−12；全部native误差参与统计。

**T2：** 四组均完成全部时间窗；最细网格归一化weighted L2终值误差≤0.05，Linf随加密下降；最细中心相对动能漂移绝对值≤1e−6，并≤同组donor能量损失的1%；donor损失>0且末档loss的h阶在[0.7,1.4]，验证可检出的O(h)系统耗散对照。

**T3：** 解析瞬时导数的weighted L2截断误差和完整周期L2误差，四组末档空间阶均≥1.7（时间步随h缩小使RK4误差高阶）；独立时间误差两档阶均≥3.5。不能用时间自收敛或全域功平衡豁免AMR/RZ空间精度。

## 失败规则

任何固定门槛失败即H10b FAIL / H10未关闭，保留所有失败行/时序；不放宽门槛、删去native边界、缩短时窗或加入黏度。若暴露H10a/前置装配缺陷，冻结反例并先写下一单修正任务，本项不得越界修改生产数值代码。仅外部依赖确实无可推进时BLOCKED；数值失败不是goal BLOCKED。

## 交付

可复现研究入口、完整原始结果、明确PASS/FAIL报告及未解决问题；单任务结束后再进入明确的下一任务，不能以研究入口完成代替H10验收。
