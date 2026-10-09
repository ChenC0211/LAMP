此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_moment_capillary_rz/README.md。

运行入口：tools\run_case.cmd hydro_moment_capillary_rz 。结果：runs/hydro_moment_capillary_rz/<时间戳_标识>/。程序：LAMP。

# RZ 椭球液滴毛细运动

<!-- case-definition-20260928 -->

**算例类别：** 轴对称毛细耦合示例。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

轴对称椭球液滴由静止开始的毛细运动，考察含方位曲率的界面受力。

## 计算域、初始条件与边界条件

轴线正则，其余外边界无滑移；曲率应包含子午曲率和方位曲率，不能直接采用平面圆滴压差。

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
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_moment_capillary_rz/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_moment_capillary_rz/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/hydro_moment_capillary_rz/system/hydro.dict) |
| 材料模型 | `twoPhase` | [hydro.dict](../../../../cases/hydro_moment_capillary_rz/system/hydro.dict) |
| 负相密度 ρ− | `1` | [hydro.dict](../../../../cases/hydro_moment_capillary_rz/system/hydro.dict) |
| 正相密度 ρ+ | `0.1` | [hydro.dict](../../../../cases/hydro_moment_capillary_rz/system/hydro.dict) |
| 负相动力黏度 μ− | `0.01` | [hydro.dict](../../../../cases/hydro_moment_capillary_rz/system/hydro.dict) |
| 正相动力黏度 μ+ | `0.001` | [hydro.dict](../../../../cases/hydro_moment_capillary_rz/system/hydro.dict) |
| 表面张力系数 σ | `1` | [hydro.dict](../../../../cases/hydro_moment_capillary_rz/system/hydro.dict) |
| 时间格式 | `cchMomentMidpoint` | [hydro.dict](../../../../cases/hydro_moment_capillary_rz/system/hydro.dict) |
| 对流格式 | `central` | [hydro.dict](../../../../cases/hydro_moment_capillary_rz/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-10` | [hydro.dict](../../../../cases/hydro_moment_capillary_rz/system/hydro.dict) |
| 背景网格划分 | `(8, 8)` | [mesh.dict](../../../../cases/hydro_moment_capillary_rz/system/mesh.dict) |
| 启用网格自适应 | `false` | [mesh.dict](../../../../cases/hydro_moment_capillary_rz/system/mesh.dict) |
| 最大单元数 | `4096` | [mesh.dict](../../../../cases/hydro_moment_capillary_rz/system/mesh.dict) |
| 时间步长 Δt | `0.002` | [run.dict](../../../../cases/hydro_moment_capillary_rz/system/run.dict) |
| 计算终时 T | `0.24` | [run.dict](../../../../cases/hydro_moment_capillary_rz/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_moment_capillary_rz/system/run.dict) |
| 场输出间隔 | `5` | [run.dict](../../../../cases/hydro_moment_capillary_rz/system/run.dict) |
| 界面输运方法 | `parametricGals` | [interface.dict](../../../../cases/hydro_moment_capillary_rz/system/interface.dict) |
| 重新初始化间隔 | `0` | [interface.dict](../../../../cases/hydro_moment_capillary_rz/system/interface.dict) |

完整输入文件：[geometry.dict](../../../../cases/hydro_moment_capillary_rz/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_moment_capillary_rz/system/hydro.dict)、[interface.dict](../../../../cases/hydro_moment_capillary_rz/system/interface.dict)、[mesh.dict](../../../../cases/hydro_moment_capillary_rz/system/mesh.dict)、[run.dict](../../../../cases/hydro_moment_capillary_rz/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_moment_capillary_rz
```

每次运行创建 `runs/hydro_moment_capillary_rz/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# 新矩轴对称毛细运动

单位 RZ 域内的椭球从静止开始，由表面张力驱动。径向/轴向半径为 0.22/0.17，中心轴高为 0.5；密度 1/0.1，黏度 0.01/0.001，表面张力 1。8×8 原生流体网格、16 段初始 Hermite 子午线，dt=0.002，请求 120 步至 T=0.24。使用 code 单位，不含电场。

这使用同一 Compatible-Moment 中点推进器。RZ 惯性由保矩局部最小动能重建定义，GALS 的共同 C1 速度及施力转置由同一个约束最小动能问题决定。质量与动量共用相通量；轴线只封闭液体体积，不承受表面张力。当前完整运动是否成功，以实际运行的退出码、接受步数及历史为准；六步模块控制不替代本例验收。

在项目根目录编译，并选一个尚不存在的运行目录：

```powershell
$cmakeExe = 'C:/Program Files/Microsoft Visual Studio/2022/Community/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe'
& $cmakeExe --build build/windows --config Release --target QUANTUM --parallel 1 -- /nr:false
$caseRun = 'runs/hydro_moment/my_capillary_rz'
New-Item -ItemType Directory -Path $caseRun
Copy-Item -LiteralPath cases/hydro_moment_capillary_rz/system -Destination $caseRun -Recurse
$env:OMP_NUM_THREADS = '4'
$env:OMP_DYNAMIC = 'FALSE'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
.\bin\LAMP.exe --case $caseRun --check-case
.\bin\LAMP.exe --case $caseRun
```

`geometry.dict` 的 `type spheroid` 只用于 RZ；`center.x` 必须为零，两个 `radii` 分别为径向和轴向半轴。轴端必须位于上下外壁之间，最大半径小于域半径；不能把 XY 的闭合 `ellipse` 轮廓直接移入这个算例。半轴相等时为数值曲率球滴初态，不预置解析压力跳跃。

`hydro.dict` 的 `cchMomentMidpoint`、`central` 与 `parametric/sharp`，以及 `interface.dict` 的 `parametricGals` 共同选择这一条路线。本例使用固定、完整、封闭的单位均匀方格，无旋流；同一入口的固定粗细叶用法见[AMR 算例](../../../../cases/hydro_moment_capillary_amr_rz/README.md)。动态 AMR、出入口、接触壁面及拓扑变化仍不支持。依赖复用本机 SuperLU 和现有 OpenMP。

输出沿用 `output/series/hydro_moment_energy.csv` 和 `hydro_moment_interface.csv`，分别记录每个接受步的物理收支与写场时刻的真实子午线。`liquid_volume` 使用 2πr 体积，`surface_energy` 使用旋转表面积；界面场中的曲率是包含环向项的总曲率，规则轴端采用解析极限。`shape_moment` 对比液体的一个横向二阶矩 `∫r²/2 dV` 与围绕液体轴向重心的 `∫(z−zc)² dV`，球为零，径向扁球为正。

场文件的单元速度是 P1 显示投影；GALS 使用完整 C1 重建。原生面文件保留两个拥有的矩，其两项代表式不是完整高阶面迹。动量功收支 `energy_balance` 与包含实际表面积变化的 `total_energy_balance` 分开记录，两者不能混为一个精度结论。方法、失败证据和后续只在 [CCH 流体待办](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md)维护。

首次正式运行位于 `runs/hydro_moment/20260927_rz_riesz`，使用冻结 `build/moment_rz_riesz_final.exe`。2026-09-28 归档快照完成 19/120 步至 t=0.038，进程仍在运行；[阶段图](../../cases/module_history/output/20260928_184817_build_archive/moment_rz_riesz_run.png)的最近轮廓为第 15 步。这些固定快照不随运行更新，最终状态需读取运行目录中的实际历史和 `run.exit`，不能以阶段图代替全部 120 步。

最终更新：同一冻结程序现已完成原定 120 步至 T=0.24，`run.exit=0`。[完整结果图](../../cases/module_history/output/20260928_184817_build_archive/moment_rz_capillary_completed.png)包含第 120 步轮廓及全部接受历史；液量最大漂移 2.082e−17、逐格分数差 ≤1.002e−12、含表面积的累计能量账最大差 9.071e−7。上述 19 步和后来的阶段图保留其历史身份；此完整运行仍不构成网格/时间收敛或破碎验收。

</details>
