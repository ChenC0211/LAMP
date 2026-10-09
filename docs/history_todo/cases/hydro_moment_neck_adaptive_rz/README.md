此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_moment_neck_adaptive_rz/README.md。

运行入口：tools\run_case.cmd hydro_moment_neck_adaptive_rz 。结果：runs/hydro_moment_neck_adaptive_rz/<时间戳_标识>/。程序：LAMP。

# RZ 局部 AMR 颈缩

<!-- case-definition-20260928 -->

**算例类别：** 无电破碎与网格耦合算例。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

联合曲率/间隙分辨率与显式颈部模型，考察局部网格支持下的真实毛细颈缩。

## 计算域、初始条件与边界条件

轴线正则、封闭外壁。现有长程结果仍存在终端区域外的曲率分辨率拒绝；不能把拒绝或未分裂标为通过。

```text
Domain { min (0, 0); max (1, 1); dimension 2; }
```

```text
ParametricInterface
{
    name liquid;
    type neckedSpheroid;
    center (0, 0.5); radii (0.35, 0.3);
    neckRatio 0.005; segments 32;
}
```

初始速度类型为 `zero`；具体初场幅值与初始压力设置见下表。

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/hydro.dict) |
| 材料模型 | `twoPhase` | [hydro.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/hydro.dict) |
| 负相密度 ρ− | `5` | [hydro.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/hydro.dict) |
| 正相密度 ρ+ | `1` | [hydro.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/hydro.dict) |
| 负相动力黏度 μ− | `0.02` | [hydro.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/hydro.dict) |
| 正相动力黏度 μ+ | `0.01` | [hydro.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/hydro.dict) |
| 表面张力系数 σ | `0.1` | [hydro.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/hydro.dict) |
| 时间格式 | `cchMomentMidpoint` | [hydro.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/hydro.dict) |
| 对流格式 | `central` | [hydro.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-10` | [hydro.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/hydro.dict) |
| 背景网格划分 | `(8, 8)` | [mesh.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/mesh.dict) |
| 启用网格自适应 | `true` | [mesh.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/mesh.dict) |
| 动态加密最高层级 | `3` | [mesh.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/mesh.dict) |
| 动态单元预算 | `4096` | [mesh.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/mesh.dict) |
| 时间步长 Δt | `0.00025` | [run.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/run.dict) |
| 计算终时 T | `0.02` | [run.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/run.dict) |
| 场输出间隔 | `10` | [run.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/run.dict) |
| 界面输运方法 | `parametricGals` | [interface.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/interface.dict) |
| 重新初始化间隔 | `0` | [interface.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/interface.dict) |

网格分辨率设置：

```text
Adaptation
{
    maximumLevel 3;
    cellBudget 4096;
    maxCurvatureCells 0.5;
    minimumGapCells 2;
}
```

完整输入文件：[geometry.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/hydro.dict)、[interface.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/interface.dict)、[mesh.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/mesh.dict)、[run.dict](../../../../cases/hydro_moment_neck_adaptive_rz/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_moment_neck_adaptive_rz
```

每次运行创建 `runs/hydro_moment_neck_adaptive_rz/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

界面场文件 `*_interface_step*.dat` 按同一 Hermite 曲线的曲率特征点及每段八等分点输出，保留材料节点两侧的曲率。`curvature_meridional` 为经向曲率，`curvature_azimuthal` 为环向曲率，`curvature_one_sided` 为二者之和；网格分辨率采用两主曲率绝对值的最大值。输出节点数随曲线改变，读取时应使用文件中的连通关系。轮廓 CSV 的固定材料参数采样适合看形状，不能用于确定窄峰的最大曲率。[方法说明 §2.45](../../docs/QUANTUM_Conservative_AMR_Hermite_GALS_TODO.md#surface-extrema)属于代码更新，本次未运行该算例，不改变旧失败轨迹的结论。

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# 局部 AMR 下的 RZ 颈缩与有限尺度分裂

原冻结轨迹只接受 6 步，第 7 步因界面近似空间容量停止，尚未分裂。§40 已将同一 H2/库存重构改为按局部支撑存储，并完成独立容量与首步对照；原字典无需改动。新的 80 步运行位于 `runs/hydro_moment/20260928_moment_material_support_v4_ordinary`，完整结果仍待收齐。方法、资源上限与精度限制见[流体 §40](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#40-按局部支撑存储同一个材料近似空间2026-09-28)。[实际初始网格图](../../cases/module_history/output/20260928_184817_build_archive/local_terminal_resolution_initial.png)仅表示 t=0。

本例从静止、连通的 Hermite 界面开始，用 Compatible-Moment 两相流体的共同中点推进计算毛细运动。初始颈半径为 0.00175；8×8 背景网格允许三级局部加密，最细间距为 0.015625，事件的必要半径条件为 0.0015625。因此初态不会仅凭半径触发事件。预先声明 dt=0.00025、T=0.02，共 80 步；实际是否分裂取决于几何、库存、动量和能量的共同可行性。

初始界面先经过分辨率审查与必要加密，然后在所得原生网格上设置初始速度。各处曲率、相对界面间距和速度重构所需的邻域提出局部层级要求；一个细颈不再直接指定整条界面采用最高层级。两相材料积分、压力—速度负伴随关系、共享质量通量和黏性应变能形式均沿用同一流体方法。

有限网格不能同时分辨无限增大的颈部曲率和数学上的零半径。启用 `Moment.neckRadiusCells` 后，程序明确记录终端区域：在本例中，以最细间距和 `maxCurvatureCells` 定义径向边界 R*=0.03125，允许的内部细颈及具有实际分裂来源的轴端进入这一范围，区域及必要流体支撑仍加密到最细级别。原有自然轴端、整个欠分辨小液滴和无关液滴的接近不自动获得豁免。终端曲率仍来自同一 Hermite 曲线，不能解读为已达到曲率精度要求。

在项目根目录的开发者 PowerShell 中运行：

```powershell
cmake --build build/windows --config Release --target QUANTUM --parallel 1 -- /nr:false
$caseRun='runs/hydro_moment/my_adaptive_neck'
if(Test-Path -LiteralPath $caseRun){throw '请使用新的运行目录'}
New-Item -ItemType Directory -Path $caseRun
Copy-Item -LiteralPath cases/hydro_moment_neck_adaptive_rz/system -Destination $caseRun -Recurse
$env:OMP_NUM_THREADS='4'
$env:OMP_DYNAMIC='FALSE'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
.\bin\LAMP.exe --case $caseRun --check-case
if($LASTEXITCODE -ne 0){throw '算例预检失败'}
.\bin\LAMP.exe --case $caseRun
```

主要输出位于 `output/series`：

- `hydro_moment_resolution.csv`：每次初始、物理步开始及事件后的分辨率审查，记录网格数、饱和情况、终端范围与曲率/间距请求。
- `hydro_moment_terminal_regions.csv`：终端区的材料参数范围、轴向范围和事件来源。这里的 `terminal_geometric_volume`（在前一文件）是整个声明区域的几何体积，不是独立液量，也不是封口实际改变的两格模型体积。
- `hydro_moment_events.csv`：实际接受的分裂与两格模型能量、体积；文件不存在表示未接受分裂。
- `hydro_moment_components.csv`、`hydro_moment_energy.csv`：各分量独立库存、完整动能/表面能账以及网格转移账。
- `hydro_moment_interface.csv`：实际轮廓。`curvature_model` 只标识由封口产生的材料来源，与更广的终端分辨率区域不是同一标签。

本例是声明有限尺度下的耦合演化验证。完整完成时窗仍不能替代曲率精度、空间收敛或自然断裂时间的独立验收。当前审查发生在已接受物理步的开始和事件之后，不是一步内连续误差控制。

</details>
