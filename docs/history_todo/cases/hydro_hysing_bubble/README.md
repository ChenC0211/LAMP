此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_hysing_bubble/README.md。

# B4：Hysing 二维上升气泡

**状态（2026-09-29）：TC1 原工况已进入集中调试，TC2 尚未启动，均未验收。** TC1 初始圆的相切测度问题修正后已完成初始化，首个耦合步仍因体积约束细化失败，接受步为 0；具体版本与原始记录见 [tc1](../../../../cases/hydro_hysing_bubble/tc1/README.md)。物理矩形域、封闭自由滑移、按相重力和分相积分输出已接代码，TC1/TC2 原工况输入及六份作者参考曲线已准备。用于评价变密度两相流的耦合精度和强形变。参考原文件与比较约定见 [reference](../../tests/reference/hysing_bubble/README.md)。

## 模型、初始与边界条件

不可压缩两相 Navier–Stokes 方程，含黏性、重力和表面张力，无电场。域 [0,1]×[0,2]，初始气泡圆心 (0.5,0.5)、半径 0.25，初始速度零；重力 (0,−0.98)，终时 T=3。上下壁无滑移，左右壁不可穿透且自由滑移，无压力出口，采用封闭域压力规范。

| 工况 | 外部液体密度 | 内部气体密度 | 外部动力黏度 | 内部动力黏度 | 表面张力 |
| --- | ---: | ---: | ---: | ---: | ---: |
| TC1 | 1000 | 100 | 10 | 1 | 24.5 |
| TC2 | 1000 | 1 | 10 | 0.1 | 1.96 |

来源：[Hysing 等（2009），Quantitative benchmark computations of two-dimensional bubble dynamics](https://www.karlin.mff.cuni.cz/~hron/NMMO403/1934_ftp.pdf)，DOI: 10.1002/fld.1934。字典须将曲线内部对应到气相，不能照搬液滴的内外物性。ρg 按真实相材料积分，不使用恒密度体力替代。

## 观测量与判据

气泡区域 Ωb 的面积为 A、周长为 P：质心 yc=∫Ωb y dA/A，平均上升速度 vc=∫Ωb vy dA/A，圆度 c=2√(πA)/P。记录完整历程、面积守恒、分辨率和能量收支，明确重力势能与边界功符号。

TC1 原文参考：最小圆度 0.9012±0.0001，最大平均上升速度 0.2419±0.0002，终时质心 1.081±0.001。这些是文献报告范围，不自动等于任意粗网格的通过容差。作者公开的三组最高层级曲线已原样保存并固定哈希；其中采样极值与上述概括范围有小差异，分别保留，详见[参考数据说明](../../tests/reference/hysing_bubble/README.md)。

正式输出 `hydro_moment_phase_integrals.csv` 每步有负相、正相两行；本例取 `phase=negative` 的 `centroid_y`、`mean_u_y`、`circularity`、`volume_geometry`，并同时报告 `volume_fv`、`geometry_minus_fv` 和 `volume_fv_drift`。这里 negative 是气泡，不等于字段旧命名中的“液体”。能量历史新增实际几何材料的重力势能及 K+σP+Eg；机械能变化仍须结合黏性、外力、网格和事件收支解释，不直接当成守恒残差。

TC2 的晚期拓扑在原文各方法间存在差异，不设唯一破碎轮廓或统一断裂时刻。比较形状、积分量、质量与网格敏感性，区分已解析演化和声明的终端事件模型；与 TC1 定量精度分别给出结论。

## 离散、运行与输出

采用同一参数 Hermite 界面、共享相通量、完整速度和共同中点，开跑前固定网格、时间步、AMR/事件和初始压力。不增黏、不增加出口换取通过。原则上两档网格，确需判断收敛阶时增加一档。

两组变体分别为 [TC1](../../../../cases/hydro_hysing_bubble/tc1/README.md) 和 [TC2](../../../../cases/hydro_hysing_bubble/tc2/README.md)，各有独立 `system/` 与 README。均保留原始物理域和完整 T=3 终时，提交的基础网格为 16×32、Δt=0.002、初始几何 64 段；分辨率是否足够由冻结后的加密对照判断。通用运行器选择 `-Case hydro_hysing_bubble/tc1` 或 `tc2`；结果分别保存在该变体的 `output/<时间戳_标识>/`，包含实际输入、程序和完整日志。执行[集中验收计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)，目前不启动物理计算。

左右壁使用 `hydro.dict` 的 `boundaryPatches`：`type freeSlip`、`side west/east`、`start 0`、`end 2`；法向速度固定为零，`value` 可省略。上下壁不列 patch，沿用无滑移。封闭壁面不配置 `portPhases`、`incomingVelocity` 或压力出口。重力设置为 `accelerationX 0; accelerationY -0.98;`；`densityNegative/viscosityNegative` 对应内部气泡，`densityPositive/viscosityPositive` 对应外部液体。当前代码和输入尚未经完整运行，不能据此关闭浮升精度或强形变/拓扑验收。

## 参考比较工具

一次运行结束后，在项目根目录分析其实际时间戳目录：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/analyze_hydro_standard.ps1 -Benchmark hysingTC1 -Run cases/hydro_hysing_bubble/tc1/output/<本次运行目录>
# TC2 使用 -Benchmark hysingTC2 及 tc2 的实际运行目录。
```

工具核对程序/输入快照、工况、边界、接受时钟及三组参考文件哈希；预检和仍在运行的任务不能用于比较。失败轨迹须显式加 `-AcceptedPrefix`，只分析实际接受区间，摘要继续标明失败及是否到达 T=3。结果写入该次运行的 `output/analysis/<时间戳>_hysingTC1/` 或 TC2 对应目录。

`history.csv` 同列几何面积、FV 库存及两者之差、圆度、质心、平均速度和能量/事件记录；三份原参考分别生成对齐曲线，`reference_errors.csv` 列出共同区间的时间加权 RMS、相对 L2、采样最大误差、采样极值及发生时间、最后比较时刻的值。`summary.json` 保存实际比较范围、完成状态和来源哈希。具体范数与不可用面积的处理见[参考说明](../../tests/reference/hysing_bubble/README.md)。比较工具已编译，尚无本项目物理结果；单次误差表不替代网格/时间敏感性及模块验收。
