此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_moment_neck_thinning_rz/README.md。

运行入口：tools\run_case.cmd hydro_moment_neck_thinning_rz 。结果：runs/hydro_moment_neck_thinning_rz/<时间戳_标识>/。程序：LAMP。

# RZ 自由颈缩与断裂

<!-- case-definition-20260928 -->

**算例类别：** 无电破碎主线算例。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

从连通颈部出发求解实际 Navier–Stokes 毛细演化，目标为颈缩、分裂和事件后液滴运动。

## 计算域、初始条件与边界条件

封闭域和轴线正则；不能人为指定颈缩速度或强制断裂。模型事件尺度、表面能变化与物理时间推进分开记录。

```text
Domain { min (0, 0); max (1, 1); dimension 2; }
```

```text
ParametricInterface
{
    name liquid;
    type neckedSpheroid;
    center (0, 0.5); radii (0.35, 0.3);
    neckRatio 0.022857142857142857; segments 32;
}
```

初始速度类型为 `zero`；具体初场幅值与初始压力设置见下表。

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/hydro.dict) |
| 材料模型 | `twoPhase` | [hydro.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/hydro.dict) |
| 负相密度 ρ− | `5` | [hydro.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/hydro.dict) |
| 正相密度 ρ+ | `1` | [hydro.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/hydro.dict) |
| 负相动力黏度 μ− | `0.02` | [hydro.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/hydro.dict) |
| 正相动力黏度 μ+ | `0.01` | [hydro.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/hydro.dict) |
| 表面张力系数 σ | `0.1` | [hydro.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/hydro.dict) |
| 时间格式 | `cchMomentMidpoint` | [hydro.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/hydro.dict) |
| 对流格式 | `central` | [hydro.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-10` | [hydro.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/hydro.dict) |
| 背景网格划分 | `(16, 16)` | [mesh.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/mesh.dict) |
| 启用网格自适应 | `false` | [mesh.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/mesh.dict) |
| 最大单元数 | `4096` | [mesh.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/mesh.dict) |
| 时间步长 Δt | `0.00025` | [run.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/run.dict) |
| 计算终时 T | `0.02` | [run.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/run.dict) |
| 场输出间隔 | `10` | [run.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/run.dict) |
| 界面输运方法 | `parametricGals` | [interface.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/interface.dict) |
| 重新初始化间隔 | `0` | [interface.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/interface.dict) |

完整输入文件：[geometry.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/hydro.dict)、[interface.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/interface.dict)、[mesh.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/mesh.dict)、[run.dict](../../../../cases/hydro_moment_neck_thinning_rz/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_moment_neck_thinning_rz
```

每次运行创建 `runs/hydro_moment_neck_thinning_rz/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# 从连通颈部开始的 RZ 毛细演化

**最新结果（流体 §36–37）：** 本模板的 N16 运行完整完成 80 步，一次模型分裂后推进 23 步两滴运动。相同程序的时间减半组完成 160 步，终态动能改变约 0.161%；全域 N32 对照也完成 80 步，但本时窗仍未分裂。[尺度对照图](../../cases/module_history/output/20260928_184817_build_archive/parametric_neck_scale.png)。本模板用于有限分辨率功能和尺度敏感性展示，不能作为已收敛的断裂时间或曲率基准。

本例从静止流体和连通的 `neckedSpheroid` 开始。颈半径 0.008，均匀网格 16×16，事件尺度为 `0.1 h = 0.00625`，所以初态不允许分裂。密度 5/1、黏度 0.02/0.01、表面张力 0.1；普通两相 NS 的共同中点推进决定后续运动，没有规定颈缩速度。

预先声明 dt=0.00025、T=0.02，共 80 步。此输入用于检验从正则连通状态进入局部事件可行域的链路；是否发生断裂以实际输出为准，不能用半径阈值推定断裂，也不代表自然断裂时刻已空间收敛。

在项目根目录的开发者 PowerShell 中：

```powershell
cmake --build build/windows --config Release --target QUANTUM --parallel 1 -- /nr:false
$caseRun='runs/hydro_moment/my_neck_thinning'
if(Test-Path -LiteralPath $caseRun){throw '请使用新的运行目录'}
New-Item -ItemType Directory -Path $caseRun
Copy-Item -LiteralPath cases/hydro_moment_neck_thinning_rz/system -Destination $caseRun -Recurse
$env:OMP_NUM_THREADS='4'
$env:OMP_DYNAMIC='FALSE'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
.\bin\LAMP.exe --case $caseRun --check-case
if($LASTEXITCODE -ne 0){throw '算例预检失败'}
.\bin\LAMP.exe --case $caseRun
```

`output/series/hydro_moment_necks.csv` 记录每个物理步开始时的实际径向极小值与事件评估：`scale_pending` 表示半径或所选截面仍太宽，`cap_pending` 表示当前独立库存与 C2 接缝条件不能产生正则、分离的帽面，`momentum_pending`/`energy_pending` 表示流体接续约束尚不满足，`ready` 才持有可接受的共同候选。未执行的评估项由 `geometry_assessed`、`caps_assessed`、`fluid_assessed` 标明，对应数值 0 不能当作物理测量。

有限、正则而尚未就绪的颈部继续原方程；无效几何、非有限值或求解故障仍报错。事件不处理普通输运失败。成功分裂另写 `hydro_moment_events.csv`，其模型体积和能量与普通物理步分账；分量库存、真实轮廓及毛细/动能账分别保存在相应 CSV 中。

当前事件支持为两个相邻轴边格中的单分支颈部。§36 起使用直接参数四次 Hermite 帽面，保留侧曲率、轴端正则性与独立体积一起定形；不需要修改本例字典。新增接缝子午曲率、局部单位径向切向、端点参数速度尺度 L 和长度方程判别式，随颈部评估输出。一般时空事件定位、断裂时间和曲率的尺度验收仍开放，详见[流体待办](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md)。

2026-09-28 的实际运行 `runs/hydro_moment/20260928_neck_thinning` 已完整完成 80 步至 T=0.02，预检和运行退出 0。终态轮廓采样的最小颈半径约 0.00499476，比初始缩小约 37.6%；液量最大漂移 6.94e−18，方程最大残差 9.81e−11。在第 57 步后第一次低于半径尺度，此后 23 次中央颈部评估均为 `cap_pending`。**分裂次数为 0**，本例不能作为完整断裂验收。

[实际颈部轮廓与半径历史](../../cases/module_history/output/20260928_184817_build_archive/moment_neck_thinning.png)显示真实连通形状及未就绪区间。所记录帽面轴端是候选诊断，未就绪帽面从未成为实际界面。下一步先分析真实接缝与帽面表示的适用性，详见[流体 §35](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#35-事件的共同可行域与从连通初态开始的毛细演化2026-09-28)；不通过扩大阈值强行分裂。

上述是 §35 旧冻结版本的记录。新完整结果位于 `runs/hydro_moment/20260928_neck_parametric_sided`，程序为 `build/parametric_neck_cap_sided.exe`，所有字典文件与旧运行相同。前 57 步动能和液量历史逐项相同，区别从事件开始：模型区域约占父液量 0.458%，初始间隙 0.01531、终态 0.02772。事件模型耗散约占本时窗表面能释放的 15%，不能仅凭小模型体积认为精度已经足够。下一阶段的尺度核对见[流体 §36](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#36-直接参数化的守恒-hermite-封口2026-09-28)；§35 的旧记录和提前停止的开发版本均保留，不与本次完整结果混用。


§37 两组输入和完整结果分别保存在 `runs/hydro_moment/20260928_neck_parametric_space` 与 `runs/hydro_moment/20260928_neck_parametric_time`。空间组只将背景划分改成 32×32；时间组只减半时间步并保持相同物理输出间隔。本模板字典未改变。当前参数路径的 AMR 曲率判据与终端事件尺度尚待统一，因此不能直接开启自动加密就认为破碎接续已完整；下一实现包见[流体 §37.4](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#374-下一主包统一局部分辨率与终端事件模型)。

</details>
