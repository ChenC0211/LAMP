此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_moment_supply_rz/README.md。

运行入口：tools\run_case.cmd hydro_moment_supply_rz 。结果：runs/hydro_moment_supply_rz/<时间戳_标识>/。程序：LAMP。

# RZ 固定湿润口供液液帽

<!-- case-definition-20260928 -->

**算例类别：** 两相供液耦合示例。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

参数 Hermite 液帽与变物性流体共同演化，液体库存由同一入口相通量增加。

## 计算域、初始条件与边界条件

南侧 r∈[0,0.5] 的外法向速度 −0.1，北侧压力 0.2；入口液相和出口外部相由 interface.dict 指定；接触圆位置固定。

```text
Domain { min (0, 0); max (1, 1); dimension 2; }
```

```text
ParametricInterface
{
    name liquid;
    type supplyCap;
    radius 0.5;
    height 0.3;
    segments 8;
}
```

初始速度类型为 `zero`；具体初场幅值与初始压力设置见下表。

显式边界输入：

```text
boundaryPatches
    {
        supply { side south; type normalVelocity; start 0; end 0.5; value -0.1; }
        outlet { side north; type pressureOutlet; start 0; end 1; value 0.2; incomingVelocity (0, 0); }
    }
```

显式边界输入：

```text
portPhases { supply liquid; outlet gas; }
```

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_moment_supply_rz/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_moment_supply_rz/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/hydro_moment_supply_rz/system/hydro.dict) |
| 材料模型 | `twoPhase` | [hydro.dict](../../../../cases/hydro_moment_supply_rz/system/hydro.dict) |
| 负相密度 ρ− | `5` | [hydro.dict](../../../../cases/hydro_moment_supply_rz/system/hydro.dict) |
| 正相密度 ρ+ | `1` | [hydro.dict](../../../../cases/hydro_moment_supply_rz/system/hydro.dict) |
| 负相动力黏度 μ− | `0.02` | [hydro.dict](../../../../cases/hydro_moment_supply_rz/system/hydro.dict) |
| 正相动力黏度 μ+ | `0.01` | [hydro.dict](../../../../cases/hydro_moment_supply_rz/system/hydro.dict) |
| 表面张力系数 σ | `0.1` | [hydro.dict](../../../../cases/hydro_moment_supply_rz/system/hydro.dict) |
| 时间格式 | `cchMomentMidpoint` | [hydro.dict](../../../../cases/hydro_moment_supply_rz/system/hydro.dict) |
| 对流格式 | `central` | [hydro.dict](../../../../cases/hydro_moment_supply_rz/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-10` | [hydro.dict](../../../../cases/hydro_moment_supply_rz/system/hydro.dict) |
| 背景网格划分 | `(4, 4)` | [mesh.dict](../../../../cases/hydro_moment_supply_rz/system/mesh.dict) |
| 启用网格自适应 | `false` | [mesh.dict](../../../../cases/hydro_moment_supply_rz/system/mesh.dict) |
| 最大单元数 | `4096` | [mesh.dict](../../../../cases/hydro_moment_supply_rz/system/mesh.dict) |
| 时间步长 Δt | `0.005` | [run.dict](../../../../cases/hydro_moment_supply_rz/system/run.dict) |
| 计算终时 T | `0.1` | [run.dict](../../../../cases/hydro_moment_supply_rz/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_moment_supply_rz/system/run.dict) |
| 场输出间隔 | `5` | [run.dict](../../../../cases/hydro_moment_supply_rz/system/run.dict) |
| 界面输运方法 | `parametricGals` | [interface.dict](../../../../cases/hydro_moment_supply_rz/system/interface.dict) |
| 重新初始化间隔 | `0` | [interface.dict](../../../../cases/hydro_moment_supply_rz/system/interface.dict) |

完整输入文件：[geometry.dict](../../../../cases/hydro_moment_supply_rz/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_moment_supply_rz/system/hydro.dict)、[interface.dict](../../../../cases/hydro_moment_supply_rz/system/interface.dict)、[mesh.dict](../../../../cases/hydro_moment_supply_rz/system/mesh.dict)、[run.dict](../../../../cases/hydro_moment_supply_rz/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_moment_supply_rz
```

每次运行创建 `runs/hydro_moment_supply_rz/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# Compatible-Moment 轴对称两相供液帽

本例已接通新矩方法的主动两相供液，并完成原定 20 步。原生单位圆柱域上，液体从底部 r∈[0,0.5] 进入，帽面固定在入口边缘并在轴上规则闭合。初始自由面 z=0.3[1−(r/0.5)²]，液/气密度 5/1、黏度 0.02/0.01、表面张力 0.1；北侧给定完整压力牵引 0.2。其余正面积边界无滑移。该粗网格例用于供液、两相惯性与毛细反馈的接续，不是喷嘴精度或破碎验收。

`geometry.dict` 的 `supplyCap` 只构造初态，此后使用唯一 Hermite 几何。入口圆盘属于液相体积的闭合边界，不是自由面，不计毛细面积。`interface.dict` 沿用 `portPhases` 声明外部纯相；入口为液体、出口环境为气体。出口 `incomingVelocity (0, 0)` 明确指定反流时的外部动量，不同时给定出口流速。流入/流出方向由实际共同 C1 迹决定。

`initialVelocity zero` 是零时间约束初始化的参考速度。真实初态必须满足非零入口与不可压条件，初始化使用真实两相材料质量范数；初帧没有已经求解的物理压力。

输入使用原中点方程，dt=0.005、20 步至 T=0.1。供液率为 π×0.5²×0.1，液量增加 6.67%。自由面离域、移动接触线、动态开口网格转移和破碎尚未完成；当前固定 4×4 网格不用于精度验收。`Restart { model none; }` 保持现有输入，附着曲线暂不支持诊断相检查点 V1。

在项目根目录运行，选用尚不存在的运行目录：

```powershell
cmake --build build/windows --config Release --target QUANTUM --parallel 1 -- /nr:false
$caseRun='runs/hydro_moment/my_supply_rz'
if (Test-Path -LiteralPath $caseRun) { throw '请换一个尚不存在的运行目录' }
New-Item -ItemType Directory -Path $caseRun
Copy-Item -LiteralPath cases/hydro_moment_supply_rz/system -Destination $caseRun -Recurse
$env:OMP_NUM_THREADS='4'
$env:OMP_DYNAMIC='FALSE'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
.\bin\LAMP.exe --case $caseRun --check-case
if ($LASTEXITCODE -ne 0) { throw '算例预检失败' }
.\bin\LAMP.exe --case $caseRun
```

若当前终端没有 `cmake`，本机程序位于 `C:/Program Files/Microsoft Visual Studio/2022/Community/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe`，用 `& '完整路径'` 替代命令名即可。复用当前已配置的 `build/windows`，不需重新安装依赖。

`output/series/hydro_moment_energy.csv` 每步记录液量、两相质量、动能、表面能及实际边界/材料功；`hydro_moment_interface.csv` 与 `output/fields` 按输入每五步输出轮廓/场。`liquid_volume_drift` 是相对初态的体积增量，应减去累计供液量后判断守恒。能量不能只看是否下降：入口持续做功，液帽的表面能和动能可以增长。边界动能通量与入口差项已经包含在输运功中，不要重复计入总账。

2026-09-28 的实际运行在 `runs/hydro_moment/20260928_supply_rz`，预检/运行退出码均为 0。最大液量收支差 5.56e−17、方程残差 9.36e−11；动量离散动能账差 5.47e−13，包含实际表面能变化后的总账差 7.69e−8。后一差项仍需后续精度收口，不能用前者替代。[实际结果图](../../cases/module_history/output/20260928_184817_build_archive/moment_supply.png)与[方法、限制及下一步](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#31-同源两相动量毛细边界功与普通供液算例2026-09-28)。

</details>
