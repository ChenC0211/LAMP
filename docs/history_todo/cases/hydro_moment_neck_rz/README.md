此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_moment_neck_rz/README.md。

运行入口：tools\run_case.cmd hydro_moment_neck_rz 。结果：runs/hydro_moment_neck_rz/<时间戳_标识>/。程序：LAMP。

# RZ 近临界连通颈部

<!-- case-definition-20260928 -->

**算例类别：** 断颈事件开发算例。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

轴对称连通液体的局部断颈模型，考察父子分量、库存与材料动量的接续。

## 计算域、初始条件与边界条件

零初速、封闭外壁与轴线正则。完整速度场的断颈接续已在[流体 §55](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#native-field-transfer)接入并编译通过，尚未用本算例进行数值验收；旧低阶事件结果仅作历史对照。

```text
Domain { min (0, 0); max (1, 1); dimension 2; }
```

```text
ParametricInterface
{
    name liquid;
    type neckedSpheroid;
    center (0, 0.5); radii (0.35, 0.3);
    neckRatio 0.005714285714285714; segments 32;
}
```

初始速度类型为 `vortex`；具体初场幅值与初始压力设置见下表。

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_moment_neck_rz/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_moment_neck_rz/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/hydro_moment_neck_rz/system/hydro.dict) |
| 材料模型 | `twoPhase` | [hydro.dict](../../../../cases/hydro_moment_neck_rz/system/hydro.dict) |
| 负相密度 ρ− | `5` | [hydro.dict](../../../../cases/hydro_moment_neck_rz/system/hydro.dict) |
| 正相密度 ρ+ | `1` | [hydro.dict](../../../../cases/hydro_moment_neck_rz/system/hydro.dict) |
| 负相动力黏度 μ− | `0.02` | [hydro.dict](../../../../cases/hydro_moment_neck_rz/system/hydro.dict) |
| 正相动力黏度 μ+ | `0.01` | [hydro.dict](../../../../cases/hydro_moment_neck_rz/system/hydro.dict) |
| 表面张力系数 σ | `0.1` | [hydro.dict](../../../../cases/hydro_moment_neck_rz/system/hydro.dict) |
| 时间格式 | `cchMomentMidpoint` | [hydro.dict](../../../../cases/hydro_moment_neck_rz/system/hydro.dict) |
| 对流格式 | `central` | [hydro.dict](../../../../cases/hydro_moment_neck_rz/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-10` | [hydro.dict](../../../../cases/hydro_moment_neck_rz/system/hydro.dict) |
| 初速幅值 | `0.03` | [hydro.dict](../../../../cases/hydro_moment_neck_rz/system/hydro.dict) |
| 背景网格划分 | `(8, 8)` | [mesh.dict](../../../../cases/hydro_moment_neck_rz/system/mesh.dict) |
| 启用网格自适应 | `false` | [mesh.dict](../../../../cases/hydro_moment_neck_rz/system/mesh.dict) |
| 最大单元数 | `4096` | [mesh.dict](../../../../cases/hydro_moment_neck_rz/system/mesh.dict) |
| 时间步长 Δt | `0.00025` | [run.dict](../../../../cases/hydro_moment_neck_rz/system/run.dict) |
| 计算终时 T | `0.005` | [run.dict](../../../../cases/hydro_moment_neck_rz/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_moment_neck_rz/system/run.dict) |
| 场输出间隔 | `5` | [run.dict](../../../../cases/hydro_moment_neck_rz/system/run.dict) |
| 界面输运方法 | `parametricGals` | [interface.dict](../../../../cases/hydro_moment_neck_rz/system/interface.dict) |
| 重新初始化间隔 | `0` | [interface.dict](../../../../cases/hydro_moment_neck_rz/system/interface.dict) |

完整输入文件：[geometry.dict](../../../../cases/hydro_moment_neck_rz/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_moment_neck_rz/system/hydro.dict)、[interface.dict](../../../../cases/hydro_moment_neck_rz/system/interface.dict)、[mesh.dict](../../../../cases/hydro_moment_neck_rz/system/mesh.dict)、[run.dict](../../../../cases/hydro_moment_neck_rz/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_moment_neck_rz
```

每次运行创建 `runs/hydro_moment_neck_rz/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# RZ 局部颈部模型与实际两相流体接续

初态是一条已经很细的连通液体，颈半径约 0.002。`neckedSpheroid` 的子午线为 `r=R sin(πs)[1−(1−neckRatio)sin²(πs)]`、`z=z₀−H cos(πs)`，随后只保留统一 Hermite 几何。`Moment.neckRadiusCells 0.1` 显式开启两轴边格的欠解析颈部模型；默认省略时关闭。

本例让检测器在接受态寻找实际径向极小值，若满足声明的模型尺度，在下一物理步开始前共同接续界面、逐格库存和流体动量。初态已经满足小半径条件，因此预期在 T=0 分裂；不是由本例预测得到自然断裂时刻。随后由同一两相 Navier–Stokes、中点、质量通量和数值毛细力推进 20 步至 T=0.005。液/气密度 5/1，黏度 0.02/0.01，表面张力 0.1；网格 8×8 为功能运行范围，不是精度设置。

从开发者 PowerShell 在项目根目录运行：

```powershell
cmake --build build/windows --config Release --target QUANTUM --parallel 1 -- /nr:false
$caseRun='runs/hydro_moment/my_neck_rz'
if(Test-Path -LiteralPath $caseRun){throw '请使用新的运行目录'}
New-Item -ItemType Directory -Path $caseRun
Copy-Item -LiteralPath cases/hydro_moment_neck_rz/system -Destination $caseRun -Recurse
$env:OMP_NUM_THREADS='4'
$env:OMP_DYNAMIC='FALSE'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
.\bin\LAMP.exe --case $caseRun --check-case
if($LASTEXITCODE -ne 0){throw '算例预检失败'}
.\bin\LAMP.exe --case $caseRun
```

事件单独写入 `hydro_moment_events.csv`：父子编号、尺度、模型体积、动量/散度残差和事件能量。`hydro_moment_energy.csv` 的事件动能变化、表面能变化和模型耗散单列，普通步能量账不把它们误算成黏性或毛细功。各滴的 `volume_change` 相对其初始化或出生时的库存；轮廓/自由面场的 `curvature_model=1` 标记来自事件帽面的材料区间，随输运和后续事件祖先关系保留。

当前支持显式声明的两格轴向单分支颈部。§36 起统一使用直接参数四次 Hermite 帽面，避免逆径向坐标，并匹配实际保留侧的曲率；此前冻结程序的数值表属于其各自历史版本。§35 起，有限正则颈部未达到截面尺度、守恒帽面分离条件或流体接续条件时，记录 pending 状态并继续同一 NS 方程；无效几何与数值求解故障仍报错。事件的观察与接受共用一个候选构造，不补量、不重置速度。各步评估写入 `hydro_moment_necks.csv`；从初始未触发状态开始的输入见 [hydro_moment_neck_thinning_rz](../../../../cases/hydro_moment_neck_thinning_rz/README.md)。

2026-09-28 的实际运行位于 `runs/hydro_moment/20260928_neck_rz`，程序冻结为 `build/moment_neck_owned.exe`，预检及运行退出 0。在 T=0 发生一次事件，随后完成全部 20 步；各子滴液量漂移 ≤3.47e−18，最大方程残差 4.28e−12，事件模型耗散 2.784e−4。模型区域约占父液量 4.19%，初始断口约 0.09646，必须进一步检验空间尺度依赖。21 行能量、41 行分量、1 行事件和六帧轮廓完整保存；[实际界面图](../../cases/module_history/output/20260928_184817_build_archive/moment_neck.png)与[方法及范围](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#33-颈部断开后材料动量与能量一起接续2026-09-28)可供核对。

随后同版完成 16×16 的 20 步、8×8 半时间步的 40 步及只细化颈部的 70 叶 20 步，均到达 T=0.005。时间步影响很小；局部与均匀细化虽给出相同初始帽面，后续动能增量仍差约 29%。当前结果尚未空间收敛。独立帽面代数计算还显示同一初态取 h=1/32 时无法得到两个分离端点，尽管半径触发条件通过；这不是 N32 的流体运行结果。详见[尺度分析与后续实施](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#34-断口的尺度影响与四次封口的可行域2026-09-28)和[实际对比图](../../cases/module_history/output/20260928_184817_build_archive/moment_neck_scale.png)。运行变体由 `build/run_moment_neck_scale.ps1 -Study space` 生成空间对照；将 `space` 换为 `time` 或 `local` 生成另两项；已有运行目录不会覆盖。

</details>
