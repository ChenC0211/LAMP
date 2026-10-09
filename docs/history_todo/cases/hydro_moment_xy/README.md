此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_moment_xy/README.md。

运行入口：tools\run_case.cmd hydro_moment_xy 。结果：runs/hydro_moment_xy/<时间戳_标识>/。程序：LAMP。

# 二维相容矩单相涡流

<!-- case-definition-20260928 -->

**算例类别：** 流体功能与精度入口。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

平面不可压缩恒物性 Navier–Stokes 流动，中心对流与能量形式黏性；局部静态加密用于观察粗细接口行为。

## 计算域、初始条件与边界条件

封闭无滑移边界；初始涡流幅值为 1。存在 RefineBoxes，不能把 adaptiveMesh=false 理解为全域均匀网格。

```text
Domain { min (0, 0); max (1, 1); dimension 2; }
```

初始速度类型为 `vortex`；具体初场幅值与初始压力设置见下表。

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_moment_xy/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_moment_xy/system/run.dict) |
| 坐标系 | `planar2D` | [hydro.dict](../../../../cases/hydro_moment_xy/system/hydro.dict) |
| 材料模型 | `constant` | [hydro.dict](../../../../cases/hydro_moment_xy/system/hydro.dict) |
| 密度 ρ | `1.3` | [hydro.dict](../../../../cases/hydro_moment_xy/system/hydro.dict) |
| 动力黏度 μ | `0.02` | [hydro.dict](../../../../cases/hydro_moment_xy/system/hydro.dict) |
| 时间格式 | `cchMomentMidpoint` | [hydro.dict](../../../../cases/hydro_moment_xy/system/hydro.dict) |
| 对流格式 | `central` | [hydro.dict](../../../../cases/hydro_moment_xy/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-11` | [hydro.dict](../../../../cases/hydro_moment_xy/system/hydro.dict) |
| 初速幅值 | `1` | [hydro.dict](../../../../cases/hydro_moment_xy/system/hydro.dict) |
| 背景网格划分 | `(8, 10)` | [mesh.dict](../../../../cases/hydro_moment_xy/system/mesh.dict) |
| 启用网格自适应 | `false` | [mesh.dict](../../../../cases/hydro_moment_xy/system/mesh.dict) |
| 最大单元数 | `2048` | [mesh.dict](../../../../cases/hydro_moment_xy/system/mesh.dict) |
| 时间步长 Δt | `0.01` | [run.dict](../../../../cases/hydro_moment_xy/system/run.dict) |
| 计算终时 T | `0.4` | [run.dict](../../../../cases/hydro_moment_xy/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_moment_xy/system/run.dict) |
| 场输出间隔 | `10` | [run.dict](../../../../cases/hydro_moment_xy/system/run.dict) |

网格分辨率设置：

```text
RefineBoxes
{
    center { min (0.25, 0.25); max (0.75, 0.75); level 1; }
    core { min (0.4375, 0.375); max (0.625, 0.5625); level 2; }
}
```

完整输入文件：[geometry.dict](../../../../cases/hydro_moment_xy/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_moment_xy/system/hydro.dict)、[mesh.dict](../../../../cases/hydro_moment_xy/system/mesh.dict)、[run.dict](../../../../cases/hydro_moment_xy/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_moment_xy
```

每次运行创建 `runs/hydro_moment_xy/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# 相容矩空间的 XY 涡流

这个算例直接调用 `HydroCompatibleMomentFlow` 的非线性中点推进。8×10 根网格经两层静态局部细化后为 188 个叶单元；恒密度 1.3、黏度 0.02，闭域无滑移，dt=0.01，40 步至 T=0.4。无界面、无电场。当前默认流体路径没有切换，本例显式选择 `cchMomentMidpoint`。

在项目根目录用 PowerShell 编译：

```powershell
$cmakeExe = 'C:/Program Files/Microsoft Visual Studio/2022/Community/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe'
& $cmakeExe --build build/windows --config Release --target QUANTUM --parallel 1 -- /nr:false
```

复制到新的运行目录后运行，无需改动已有字典：

```powershell
New-Item -ItemType Directory -Path runs/hydro_moment/my_run
Copy-Item -LiteralPath cases/hydro_moment_xy/system -Destination runs/hydro_moment/my_run -Recurse
$env:OMP_NUM_THREADS = '4'
$env:OMP_DYNAMIC = 'FALSE'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
.\bin\LAMP.exe --case runs/hydro_moment/my_run --check-case
.\bin\LAMP.exe --case runs/hydro_moment/my_run
```

已复用本机 SuperLU，无需安装新库。线程设置沿用项目配置，不表示此矩方法的所有局部装配和 SuperLU 分解都已经并行。

`hydro.dict` 中的关键设置：

| 键 | 含义 |
| --- | --- |
| `timeScheme cchMomentMidpoint` | 新面/单元矩状态；与旧 `cchMidpoint` 的面平均状态不同 |
| `initialVelocity zero/vortex` | 静止或由解析流函数积分的无散度初值；不接受旧 `smooth` |
| `velocityAmplitude` | 下述流函数的速度尺度系数，并非最大流速 |
| `Moment.convection central/traceResidual` | 默认 `central` 中心通量；`traceResidual` 仅显式启用非负面速度差耗散。本例现使用中心格式，旧运行结果保留原耗散设置 |
| `Moment.maximumIterations` | 非线性中点迭代上限，默认 24，可取 1–100 |
| `linearResidualTolerance` | 同时约束线性及实际非线性方程残差，本例 1e−11 |
| `accelerationX/Y` | 可选恒定加速度，默认 0；恒密度下通过压力势表示 |

域宽、高为 W、H，归一坐标为 X、Y，L=min(W,H)。涡流初值来自 `psi=32*velocityAmplitude*L*X²(1−X)²*Y²(1−Y)²`，速度取 `(∂y psi,−∂x psi)`；所有面与单元矩在模块内积分。均匀加速度 a 对应零均值压力势 `p_a=rho*a·(x−domainCenter)`，求解剩余压力后输出总物理压力，不由微小假速度平衡体力。

`mesh.dict` 的 `RefineBoxes` 只在初始化时应用：按单元中心是否位于矩形内、逐级细化至 `level`。允许最多 32 个矩形、level 为 1–8，矩形必须在域内；实际网格仍须满足求解器的至多 2:1 分裂边条件。`maxCellNum` 必须足够。`adaptiveMesh false` 表示运行中不改变网格，不等于没有初始细化。这个入口尚未实现动态 AMR 转移。

输出读取：

- `output/series/hydro_moment_energy.csv`：每个接受态一行，含初态；保存质量、真实 H 内积动能、累计面差/黏性耗散、体力/压力/中心对流功、能量差、迭代数和残差。
- `output/fields/hydro_stepNNNN.dat`：原生叶单元上的速度、压力及其一次梯度。速度属于步终点，压力属于实际中点；`PressureTime` 明确时刻。初帧 `PressureAvailable=0`，压力列零只是占位，不代表已求得初始物理压力。
- `output/fields/hydro_faces_stepNNNN.dat`：实际法向均值、一次系数、切向中心/尺度、面积和体积流率。法向为竖直面的 +x、水平面的 +y，`u_n(s)=mean+linear*(s−center)/scale`。它们不是由显示速度插值而来。

默认在第 0/10/20/30/40 步输出上述两类场。迭代失败保持旧状态；写文件失败可能发生于新物理态已经接受之后，此时程序停止并报告 accepted/output 两个步数，不重算该步。当前没有新矩重启文件。

2026-09-27 实际运行见 [run.log](../../../../runs/hydro_moment/20260927_xy_refined_final/run.log)：40 步完成，动能 0.04021644202 → 0.02087010804，下降 48.11%；累计黏性耗散 0.01924809877、面差耗散 0.00009823520，最大相对能量收支差 5.437e−15。独立完整 NS 制造流动的空间精度见[流体待办 §13](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md)。本例是流体模块接入证据，不是喷雾破碎验收。

本例为封闭 XY、恒定单相物性及静态细化。同一新矩推进器的轴对称入口现见 [hydro_moment_rz](../../../../cases/hydro_moment_rz/README.md)，使用独立 RZ 流动证据。新矩两相、开放边界、动态 AMR、GALS 光滑速度回溯和带电耦合仍需完成。

</details>
