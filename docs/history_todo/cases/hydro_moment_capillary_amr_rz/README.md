此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_moment_capillary_amr_rz/README.md。

运行入口：tools\run_case.cmd hydro_moment_capillary_amr_rz 。结果：runs/hydro_moment_capillary_amr_rz/<时间戳_标识>/。程序：LAMP。

# RZ 粗细网格液滴毛细运动

<!-- case-definition-20260928 -->

**算例类别：** 静态 AMR 耦合示例。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

使用预设 RefineBoxes 的轴对称两相液滴，考察界面经过粗细网格支撑时的相通量与受力。

## 计算域、初始条件与边界条件

模板的 adaptiveMesh=false，属于预先局部加密，不是运行时自动 AMR；轴线正则，其他外壁无滑移。

```text
Domain { min (0, 0); max (1, 1); dimension 2; }
```

```text
ParametricInterface
{
    name liquid;
    type spheroid;
    center (0, 0.5);
    radii (0.22, 0.17);
    segments 16;
}
```

初始速度类型为 `zero`；具体初场幅值与初始压力设置见下表。

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/hydro.dict) |
| 材料模型 | `twoPhase` | [hydro.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/hydro.dict) |
| 负相密度 ρ− | `1` | [hydro.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/hydro.dict) |
| 正相密度 ρ+ | `0.1` | [hydro.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/hydro.dict) |
| 负相动力黏度 μ− | `0.01` | [hydro.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/hydro.dict) |
| 正相动力黏度 μ+ | `0.001` | [hydro.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/hydro.dict) |
| 表面张力系数 σ | `1` | [hydro.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/hydro.dict) |
| 时间格式 | `cchMomentMidpoint` | [hydro.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/hydro.dict) |
| 对流格式 | `central` | [hydro.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-10` | [hydro.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/hydro.dict) |
| 背景网格划分 | `(4, 4)` | [mesh.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/mesh.dict) |
| 启用网格自适应 | `false` | [mesh.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/mesh.dict) |
| 最大单元数 | `4096` | [mesh.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/mesh.dict) |
| 时间步长 Δt | `0.002` | [run.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/run.dict) |
| 计算终时 T | `0.24` | [run.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/run.dict) |
| 场输出间隔 | `5` | [run.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/run.dict) |
| 界面输运方法 | `parametricGals` | [interface.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/interface.dict) |
| 重新初始化间隔 | `0` | [interface.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/interface.dict) |

网格分辨率设置：

```text
RefineBoxes
{
    upper_axis { min (0, 0.5); max (0.5, 1); level 1; }
}
```

完整输入文件：[geometry.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/hydro.dict)、[interface.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/interface.dict)、[mesh.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/mesh.dict)、[run.dict](../../../../cases/hydro_moment_capillary_amr_rz/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_moment_capillary_amr_rz
```

每次运行创建 `runs/hydro_moment_capillary_amr_rz/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# 固定粗细网格上的轴对称毛细运动

单位 RZ 域，4×4 背景的左上四分之一区域细化一级，共 28 个实际叶单元；网格在本算例中保持固定。初始椭球中心轴高 0.5，径向/轴向半径 0.22/0.17，跨过 z=0.5 的粗细接缝。密度 1/0.1、黏度 0.01/0.001、表面张力 1，从静止开始；dt=0.002，默认请求 120 步至 T=0.24。沿用 [均匀 RZ 算例](../../../../cases/hydro_moment_capillary_rz/README.md)的物性、边界和输出定义。

同一 Compatible-Moment 中点推进器消费 Hermite 界面的实际相通量。粗细网格共同速度由连续的完整 C2 多项式空间作保矩最小动能重建，界面力经同一映射的转置返回；独立液量直接由共享子面通量推进。`RefineBoxes` 指定初始固定加密，与运行中的 `Adaptation` 不同；本例保持固定网格；动态细化现见 [自适应 RZ 输入](../../../../cases/hydro_moment_adaptive_rz/README.md)，自动粗化、接触壁面、开口及拓扑变化仍未接通。

在项目根目录使用尚不存在的运行目录：

```powershell
cmake --build build/windows --config Release --target QUANTUM --parallel 1 -- /nr:false
$caseRun = 'runs/hydro_moment/my_capillary_amr_rz'
New-Item -ItemType Directory -Path $caseRun
Copy-Item -LiteralPath cases/hydro_moment_capillary_amr_rz/system -Destination $caseRun -Recurse
$env:OMP_NUM_THREADS = '4'
$env:OMP_DYNAMIC = 'FALSE'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
.\bin\LAMP.exe --case $caseRun --check-case
.\bin\LAMP.exe --case $caseRun
```

`hydro_moment_energy.csv` 记录每个接受步的液量、动能、表面积及能量收支；`hydro_moment_interface.csv` 和场文件记录实际子午线。单元速度输出是 P1 显示投影，界面运动使用完整共同速度。六步模块控制、输入检查和完整 120 步运行是不同范围的证据；实际状态见运行的接受历史、退出码及[流体待办](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md)。本算例本身不证明网格收敛、动态 AMR 或断裂完成。

2026-09-28：同一正式两相模块的 XY/RZ 固定 AMR 六步控制均完成，液量漂移分别不超过 1.39e−17/6.94e−18，逐格几何—库存分数差约 1.001e−12。普通运行副本 "runs/hydro_moment/20260928_native_c1_amr_rz" 已通过冻结 App 的输入预检；本轮没有另开默认 120 步长程。方法与[实际图](../../cases/module_history/output/20260928_184817_build_archive/moment_native_c1_coupled.png)见流体待办 §26。

</details>
