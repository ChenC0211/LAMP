此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_moment_rz/README.md。

运行入口：tools\run_case.cmd hydro_moment_rz 。结果：runs/hydro_moment_rz/<时间戳_标识>/。程序：LAMP。

# RZ 相容矩单相涡流

<!-- case-definition-20260928 -->

**算例类别：** 流体功能与精度入口。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

轴对称不可压缩恒物性流动，包含径向、轴向及环向应变耗散，采用局部静态加密。

## 计算域、初始条件与边界条件

r=0 为轴线，其他外边界为封闭无滑移边界；初始涡流幅值为 1。不是一般三维流动。

```text
Domain { min (0, 0); max (1, 1); dimension 2; }
```

初始速度类型为 `vortex`；具体初场幅值与初始压力设置见下表。

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_moment_rz/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_moment_rz/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/hydro_moment_rz/system/hydro.dict) |
| 材料模型 | `constant` | [hydro.dict](../../../../cases/hydro_moment_rz/system/hydro.dict) |
| 密度 ρ | `1.3` | [hydro.dict](../../../../cases/hydro_moment_rz/system/hydro.dict) |
| 动力黏度 μ | `0.02` | [hydro.dict](../../../../cases/hydro_moment_rz/system/hydro.dict) |
| 时间格式 | `cchMomentMidpoint` | [hydro.dict](../../../../cases/hydro_moment_rz/system/hydro.dict) |
| 对流格式 | `central` | [hydro.dict](../../../../cases/hydro_moment_rz/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-11` | [hydro.dict](../../../../cases/hydro_moment_rz/system/hydro.dict) |
| 初速幅值 | `1` | [hydro.dict](../../../../cases/hydro_moment_rz/system/hydro.dict) |
| 背景网格划分 | `(8, 10)` | [mesh.dict](../../../../cases/hydro_moment_rz/system/mesh.dict) |
| 启用网格自适应 | `false` | [mesh.dict](../../../../cases/hydro_moment_rz/system/mesh.dict) |
| 最大单元数 | `2048` | [mesh.dict](../../../../cases/hydro_moment_rz/system/mesh.dict) |
| 时间步长 Δt | `0.01` | [run.dict](../../../../cases/hydro_moment_rz/system/run.dict) |
| 计算终时 T | `0.4` | [run.dict](../../../../cases/hydro_moment_rz/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_moment_rz/system/run.dict) |
| 场输出间隔 | `10` | [run.dict](../../../../cases/hydro_moment_rz/system/run.dict) |

网格分辨率设置：

```text
RefineBoxes
{
    center { min (0.25, 0.25); max (0.75, 0.75); level 1; }
    core { min (0.4375, 0.375); max (0.625, 0.5625); level 2; }
}
```

完整输入文件：[geometry.dict](../../../../cases/hydro_moment_rz/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_moment_rz/system/hydro.dict)、[mesh.dict](../../../../cases/hydro_moment_rz/system/mesh.dict)、[run.dict](../../../../cases/hydro_moment_rz/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_moment_rz
```

每次运行创建 `runs/hydro_moment_rz/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# 相容矩空间的轴对称涡流

本例使用 `cchMomentMidpoint`，直接调用共同的流体模块求解无旋流的轴对称 Navier–Stokes 方程。半径、轴向坐标分别为 x、y，区域 `[0,1]×[0,1]` 包含轴线。8×10 根网格经两层静态局部细化后为 188 个叶单元；密度 1.3、黏度 0.02，dt=0.01，40 步至 T=0.4。轴上径向速度为零，其余外壁无滑移。没有运动界面或电场。

在项目根目录的 PowerShell 中编译并运行：

```powershell
$cmakeExe = 'C:/Program Files/Microsoft Visual Studio/2022/Community/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe'
& $cmakeExe --build build/windows --config Release --target QUANTUM --parallel 1 -- /nr:false
New-Item -ItemType Directory -Path runs/hydro_moment/my_rz_run
Copy-Item -LiteralPath cases/hydro_moment_rz/system -Destination runs/hydro_moment/my_rz_run -Recurse
$env:OMP_NUM_THREADS = '4'
$env:OMP_DYNAMIC = 'FALSE'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
.\bin\LAMP.exe --case runs/hydro_moment/my_rz_run --check-case
.\bin\LAMP.exe --case runs/hydro_moment/my_rz_run
```

无需改字典或安装依赖。请为下一次计算使用新的运行目录。现有 OpenMP 设置不代表局部装配和 SuperLU 分解已全部并行。

`hydro.dict` 的关键设置是 `coordinateSystem axisymmetricRZ`、`timeScheme cchMomentMidpoint`。`Moment.convection central` 是当前默认且为本例配置；`traceResidual` 仅显式启用非负面速度差耗散。二者使用同一状态和时间推进，旧运行结果保留原耗散设置。`Moment.maximumIterations` 是非线性迭代上限，`linearResidualTolerance` 同时约束实际非线性残差。`initialVelocity zero` 给出静止初态；`vortex` 由下述流函数生成，要求区域从 r=0 开始。

令 R 为区域最大半径，H 为高度，L=min(R,H)，X=r/R，Y=(z−z_min)/H，A=`velocityAmplitude`。涡流初值为

\[
\Psi=8ALR\,X^2(1-X^2)^2Y^2(1-Y)^2,\qquad
u_r=\Psi_z/r,\quad u_z=-\Psi_r/r.
\]

程序直接计算轴上有限的多项式速度；A 是尺度系数，不是最大流速。实际面法向矩和六个独立单元量由同一速度积分，质量、压力、黏性和输运共用这些量。选取静止初态和轴向 `accelerationY` 时，均匀重力由同一压力势平衡。

输出约定与 [XY 算例](../../../../cases/hydro_moment_xy/README.md)相同：

- `output/series/hydro_moment_energy.csv` 每个接受态一行，包含真实轴对称质量、离散动能、累计耗散、功和能量差。
- `hydro_stepNNNN.dat` 保存原生叶单元上的一次速度与压力。速度属于步终点，压力属于中点；初帧压力不可用。`Coordinates` 标明 x=r、y=z，速度列 `u_x/u_y` 分别表示径向/轴向分量。
- `hydro_faces_stepNNNN.dat` 保存实际法向面多项式和体积流率。轴面面积、径向速度和流率严格为零，未引入轴面未知量。

默认在 0/10/20/30/40 步输出场。显示用的一次速度投影不等同于体内的正则输运重构；不能将显示梯度直接用作 GALS 的光滑回溯速度。

当前入口限于封闭、恒定单相物性、静态至多 2:1 网格。变质量两相、运动界面、动态 AMR、供液/出流和带电耦合仍待接通；本例不是射流断裂算例。实际验证结果统一维护在[流体待办 §16](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md)。

2026-09-27 的[实际运行](../../../../runs/hydro_moment/20260927_rz_stable_coordinates/run.log)完成 40 步，动能由 0.19350333949 降至 0.10453096029，下降 45.98%；最大相对累计能量差 1.553e−13。见[初末速度和能量图](../../cases/module_history/output/20260928_184817_build_archive/moment_rz_case_evolution.png)。另一静水控制在轴向加速度 −0.2 下完成相同时段，速度严格为零、压力最大误差 4.164e−17。

</details>
