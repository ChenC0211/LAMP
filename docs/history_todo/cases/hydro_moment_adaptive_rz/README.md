此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_moment_adaptive_rz/README.md。

运行入口：tools\run_case.cmd hydro_moment_adaptive_rz 。结果：runs/hydro_moment_adaptive_rz/<时间戳_标识>/。程序：LAMP。

# RZ 自适应网格上的毛细液滴

<!-- case-definition-20260928 -->

**算例类别：** 动态 AMR 功能示例。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

以实际接受界面的曲率和相向界面间距驱动网格分辨率，考察相体积库存与完整流体状态的联合转移。当前输入显式启用有滞回和速度误差约束的粗化；这是功能示例，三步时窗不保证发生粗化，也不作为网格自适应精度验收。

## 计算域、初始条件与边界条件

轴线正则、外壁无滑移；必要细化超过最高层级或单元预算时停止。完整速度转移与自动粗化见[流体 §57](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#automatic-coarsening)。粗化候选若不满足动量、不可压、能量或速度误差约束，则保留当前网格继续推进。2026-09-28 的本次更新只接通代码和输入，未运行物理时窗，数值验收仍待集中完成。

```text
Domain { min (0, 0); max (1, 1); dimension 2; }
```

```text
ParametricInterface
{
    name liquid;
    type spheroid;
    center (0, 0.5);
    radii (0.40, 0.36);
    segments 16;
}
```

初始速度类型为 `zero`；具体初场幅值与初始压力设置见下表。

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_moment_adaptive_rz/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_moment_adaptive_rz/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/hydro_moment_adaptive_rz/system/hydro.dict) |
| 材料模型 | `twoPhase` | [hydro.dict](../../../../cases/hydro_moment_adaptive_rz/system/hydro.dict) |
| 负相密度 ρ− | `1` | [hydro.dict](../../../../cases/hydro_moment_adaptive_rz/system/hydro.dict) |
| 正相密度 ρ+ | `0.1` | [hydro.dict](../../../../cases/hydro_moment_adaptive_rz/system/hydro.dict) |
| 负相动力黏度 μ− | `0.01` | [hydro.dict](../../../../cases/hydro_moment_adaptive_rz/system/hydro.dict) |
| 正相动力黏度 μ+ | `0.001` | [hydro.dict](../../../../cases/hydro_moment_adaptive_rz/system/hydro.dict) |
| 表面张力系数 σ | `1` | [hydro.dict](../../../../cases/hydro_moment_adaptive_rz/system/hydro.dict) |
| 时间格式 | `cchMomentMidpoint` | [hydro.dict](../../../../cases/hydro_moment_adaptive_rz/system/hydro.dict) |
| 对流格式 | `central` | [hydro.dict](../../../../cases/hydro_moment_adaptive_rz/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-10` | [hydro.dict](../../../../cases/hydro_moment_adaptive_rz/system/hydro.dict) |
| 背景网格划分 | `(4, 4)` | [mesh.dict](../../../../cases/hydro_moment_adaptive_rz/system/mesh.dict) |
| 启用网格自适应 | `true` | [mesh.dict](../../../../cases/hydro_moment_adaptive_rz/system/mesh.dict) |
| 动态加密最高层级 | `1` | [mesh.dict](../../../../cases/hydro_moment_adaptive_rz/system/mesh.dict) |
| 动态单元预算 | `128` | [mesh.dict](../../../../cases/hydro_moment_adaptive_rz/system/mesh.dict) |
| 时间步长 Δt | `0.002` | [run.dict](../../../../cases/hydro_moment_adaptive_rz/system/run.dict) |
| 计算终时 T | `0.006` | [run.dict](../../../../cases/hydro_moment_adaptive_rz/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_moment_adaptive_rz/system/run.dict) |
| 场输出间隔 | `1` | [run.dict](../../../../cases/hydro_moment_adaptive_rz/system/run.dict) |
| 界面输运方法 | `parametricGals` | [interface.dict](../../../../cases/hydro_moment_adaptive_rz/system/interface.dict) |
| 重新初始化间隔 | `0` | [interface.dict](../../../../cases/hydro_moment_adaptive_rz/system/interface.dict) |

网格分辨率设置：

```text
Adaptation
{
    maximumLevel 1;
    cellBudget 128;
    maxCurvatureCells 0.5;
    minimumGapCells 2;
    Coarsening
    {
        hysteresis 0.5;
        relativeVelocityError 0.001;
    }
}
```

粗化每次最多合并一代完整兄弟单元。`hysteresis=0.5` 表示合并后的曲率尺度、间距指标须留出一倍余量；保护带按父单元范围计算，包含已有预留带。`relativeVelocityError=0.001` 是质量加权速度范数的相对误差预算，分别约束当前端点、上一端点及中点，不是动能损失率或曲率误差。移除整个 `Coarsening` 块即可恢复只细化的策略。

完整输入文件：[geometry.dict](../../../../cases/hydro_moment_adaptive_rz/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_moment_adaptive_rz/system/hydro.dict)、[interface.dict](../../../../cases/hydro_moment_adaptive_rz/system/interface.dict)、[mesh.dict](../../../../cases/hydro_moment_adaptive_rz/system/mesh.dict)、[run.dict](../../../../cases/hydro_moment_adaptive_rz/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_moment_adaptive_rz
```

每次运行创建 `runs/hydro_moment_adaptive_rz/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

`output/series/hydro_moment_coarsening.csv` 记录候选及逐时层结果：`ready` 表示候选满足转移条件；`accuracy_pending`、`momentum_pending`、`energy_pending` 表示本次不粗化。`hydro_moment_adaptation.csv` 只记录实际接受的网格变更，`operation` 区分细化与粗化，并给出最大相对速度误差。未出现粗化候选时，不产生粗化记录文件；不能据此声称粗化已经验证。

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# Hermite 界面与 Compatible-Moment 流体的联合动态加密

这是一个三步接续算例：单位 RZ 域，4×4 根网格，初始椭球半轴 0.40/0.36，轴向中心 0.5；密度 1/0.1、黏度 0.01/0.001、表面张力 1，从静止开始。dt=0.002，终时 0.006。当前入口在 t=0 先审查初始几何并完成 16→64 格加密，然后在最终初始网格上设置速度，三个物理步均从这一网格开始。以前“先推进一个根网格步”的结果只保留为历史证据。

从项目根目录运行，选择尚不存在的运行目录：

```powershell
cmake --build build/windows --config Release --target QUANTUM --parallel 1 -- /nr:false
$caseRun = 'runs/hydro_moment/my_adaptive_rz'
New-Item -ItemType Directory -Path $caseRun
Copy-Item -LiteralPath cases/hydro_moment_adaptive_rz/system -Destination $caseRun -Recurse
$env:OMP_NUM_THREADS='4'
$env:OMP_DYNAMIC='FALSE'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
.\bin\LAMP.exe --case $caseRun --check-case
.\bin\LAMP.exe --case $caseRun
```

`mesh.dict` 的 `adaptiveMesh true` 与 `Adaptation` 配套，容量只由 `cellBudget` 指定；不要再填写 `QuadtreeMesh.maxCellNum`。`maxCurvatureCells` 约束最大主曲率与局部网格宽度的乘积，`minimumGapCells` 约束沿采样法线找到的相向界面间距。这里的 0.5 和 2 是小型接续演示所用的粗分辨率，不能作为精确曲率或破碎计算的推荐精度。加密带、2:1 分配及资源限制仍由现有 Interface/Core 策略处理，超出最高层级或预算会明确失败，不默默放宽判据。

一次网格变化不移动 Hermite 曲线，也不改变其曲率/表面能。独立液量守恒转移；已有流场的动态加密同时转移所需速度时层，满足动量、不可压与动能不增加约束。初始加密发生在流场建立之前，无需转移粗网格流体历史。物理推进仍使用同源质量—动量通量、共同中点、负伴随压力和应变能黏性。当前普通自适应策略只细化，联合转移模块另有实际粗化控制；自动粗化和开放边界下的联合 AMR 仍待实现。

输出位于运行目录的 `output`：

- `series/hydro_moment_energy.csv`：每个物理接受步的液量、能量和方程残差；`cumulative_mesh_energy_change`、`leaf_cells` 为累计网格动能变化与叶数。能量收支单独扣除网格变化，不将其算作物理黏性。
- `series/hydro_moment_resolution.csv`：t=0 及各物理步开始时的几何/网格审查；初始 16→64 格变化记在这里。
- `series/hydro_moment_adaptation.csv`：每次实际接受的零时间网格变化，记录前后叶数、液量差、动能变化、可行最小动能及动量/散度残差。这里的 `step,time` 是转移发生时的已接受物理时刻。
- `series/hydro_moment_interface.csv` 与 `fields`：实际接受曲线、原生叶上的流场和面通量矩。流场在下一实际求解后输出，压力时间仍是该步中点。

如果网格候选不可行，准备过程保留旧状态并停止本次运行；如果候选已接受而随后物理步失败，日志会保留该次零时间网格变化，物理时钟不会前进。输入预检不执行加密判据或物理步，不等于运行成功。实际结果与方法限制统一见[流体待办 §27](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md)。

</details>
