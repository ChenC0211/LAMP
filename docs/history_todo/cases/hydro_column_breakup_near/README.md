此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_column_breakup_near/README.md。

运行入口：tools\run_case.cmd hydro_column_breakup_near 。结果：runs/hydro_column_breakup_near/<时间戳_标识>/。程序：LAMP。

# 近断裂轴对称液柱

<!-- case-definition-20260928 -->

**算例类别：** 旧保守路线开发对照。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

无黏两相毛细破碎；近断裂初态缩短颈缩等待时间，研究显式断颈及子滴运动。

## 计算域、初始条件与边界条件

轴线正则，轴向两端为反射截面；模板为半波长，旧专用脚本的 HalfWaves=2 才生成整波长。不得与当前完整速度路线的验收混用。

```text
Domain { min (0, 0); max (4, 5.2359877559829888); dimension 2; }
```

```text
AxisymmetricColumn { levelSet liquid; radius 1; relativeAmplitude 0.954754144451598; }
```

初始速度类型为 `zero`；具体初场幅值与初始压力设置见下表。

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_column_breakup_near/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_column_breakup_near/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/hydro_column_breakup_near/system/hydro.dict) |
| 材料模型 | `twoPhase` | [hydro.dict](../../../../cases/hydro_column_breakup_near/system/hydro.dict) |
| 负相密度 ρ− | `1` | [hydro.dict](../../../../cases/hydro_column_breakup_near/system/hydro.dict) |
| 正相密度 ρ+ | `0.0001` | [hydro.dict](../../../../cases/hydro_column_breakup_near/system/hydro.dict) |
| 负相动力黏度 μ− | `0` | [hydro.dict](../../../../cases/hydro_column_breakup_near/system/hydro.dict) |
| 正相动力黏度 μ+ | `0` | [hydro.dict](../../../../cases/hydro_column_breakup_near/system/hydro.dict) |
| 表面张力系数 σ | `1` | [hydro.dict](../../../../cases/hydro_column_breakup_near/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-8` | [hydro.dict](../../../../cases/hydro_column_breakup_near/system/hydro.dict) |
| 初始压力跳跃 | `1` | [hydro.dict](../../../../cases/hydro_column_breakup_near/system/hydro.dict) |
| 背景网格划分 | `(64, 84)` | [mesh.dict](../../../../cases/hydro_column_breakup_near/system/mesh.dict) |
| 启用网格自适应 | `false` | [mesh.dict](../../../../cases/hydro_column_breakup_near/system/mesh.dict) |
| 最大单元数 | `12000` | [mesh.dict](../../../../cases/hydro_column_breakup_near/system/mesh.dict) |
| 时间步长 Δt | `0.0005` | [run.dict](../../../../cases/hydro_column_breakup_near/system/run.dict) |
| 计算终时 T | `0.2` | [run.dict](../../../../cases/hydro_column_breakup_near/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_column_breakup_near/system/run.dict) |
| 场输出间隔 | `10` | [run.dict](../../../../cases/hydro_column_breakup_near/system/run.dict) |
| 界面输运方法 | `conservative` | [interface.dict](../../../../cases/hydro_column_breakup_near/system/interface.dict) |
| 重新初始化间隔 | `0` | [interface.dict](../../../../cases/hydro_column_breakup_near/system/interface.dict) |

完整输入文件：[geometry.dict](../../../../cases/hydro_column_breakup_near/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_column_breakup_near/system/hydro.dict)、[interface.dict](../../../../cases/hydro_column_breakup_near/system/interface.dict)、[mesh.dict](../../../../cases/hydro_column_breakup_near/system/mesh.dict)、[run.dict](../../../../cases/hydro_column_breakup_near/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_column_breakup_near
```

每次运行创建 `runs/hydro_column_breakup_near/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# 近断裂液柱：主线开发算例

最新状态（2026-09-26）：近断裂初态和 OpenMP 12 / OpenBLAS 1 保持，无需新增字典设置。分裂后的局部库存约束同步已接入并通过集中包；新无电主线在运行，完整后程未确认。旧无电轨迹在分裂后继续 194 步；带电轨迹终止于第 168 步、t=0.0823046875，分裂候选的轴上气隙仍缺全域电场闭合。局部库存同步尚不解决带电事件电场问题，本轮没有重复启动带电长跑。以下历史记录保留，当前状态见[进度说明](../../docs/QUANTUM_Hydro_EHD_Progress.md)。

当前内部分裂主线使用**整波长**：下面命令自动生成 64×168 网格和中央细颈，无需手改字典。初始颈半径约 0.0375R（0.6h），实际收缩跨过 0.4h 才触发 1→2 数值事件，初速仍为零。默认 OpenMP 12 / OpenBLAS 1；不带 `-HalfWaves 2` 时保留原半波长截面脱离算例。

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/run_hydro_column_breakup.ps1 -HalfWaves 2 -Threads 12 -Mainline
```

已接通随相通量传递的子滴库存、分子滴毛细力和下一步 NS；`interface_components.csv` 输出材料液量/残差，保存场包含 `C_liquid_<id>`。[上一条整波长运行](../../../../runs/hydro_column_breakup/20260926_094340_interior_bounded_material_58fe5f/breakup_analysis.json)在 t=0.0815 分裂后继续 77 步，最终于 t=0.1159921875 后因端部积分失败退出（rc=3）。当前源码已接入端部积分坐标分段和同一步曲率共享；解析测度/电荷及集中耦合检查通过，[一条修正后的主运行](../../../../runs/hydro_column_breakup/20260926_103327_interior_chart_coords_26929c/run.log)已跨过旧积分失败时段，最终在第 333 步、t=0.15721875 停止，分裂后继续 170 步；第 334 步折返尖端竞争曲率投影导致 rc=3，尚未达到终时。上述启动命令即可使用当前程序，无需新增字典项。 带电整波长开发路径现可在命令中加入 `-Charged`；脚本生成所需电学字典，默认初始面电荷密度 0.001（代码单位），并输出 `conservative_electric_components.csv`。首条带电整波长运行在第 141 步后因细颈 P3 窗口缺少液相样本退出，未到达事件。当前已接入细颈调和电场及集中检查，并启动一条新带电主线，实际跨越仍待验证。当前无电程序已接入事件区多面元积分受力，集中检查通过；一条新的完整主线运行中，尚未形成终态验收。外侧轴向边界仍为反射截面，完整终时、网格/时间验收及喷嘴自由下游运动尚未完成。

当前加速版保留生产积分误差要求，把局部最大深度与 4096 叶工作量上限分开，并记录 `max_spatial_depth`、`max_spatial_leaves`。使用下述原命令即可，不需修改字典。解析输运、并行一致性与 Hydro 接续检查通过；[实际主复算](../../../../runs/hydro_column_breakup/20260926_033558_near_local_spatial_be9057/breakup_analysis.json)已越过原第 360 步积分缩步点，保持 dt=0.00025。该轨迹最终接受至 t=0.18775，第 519 步因折返尖端毛细投影歧义退出；步长减半也在 t=0.18490039 遇同类问题。终时 0.2 尚未完成，正在处理共同的后段受力问题；不据局部改善声明全程成功或全程加速比。

后续启动默认启用 `StepRetry.recoverAfter 8`：连续 8 步无需重试后，候选步长增加 25%，上限仍为初始 dt；所有物理门限和失败回滚保留。`-RecoverAfter 0` 保持原来的只减不增策略；与现有冻结主档作时间/空间比较时应显式选择这一相同设置。策略加入前的冻结轨迹均未启用恢复，不受后续模板变化影响。

2026-09-26 后段改进：旧截面前的 2h 气相带由独立库存确认后，事件区恢复库存驱动几何重构，避免永久保留预测标量。无需新增字典。[无电复算](../../../../runs/hydro_column_breakup/20260926_045406_near_inventory_handoff_05bd01/breakup_analysis.json)使用 `-RecoverAfter 0`，在第 203 步交接，最后接受第 290 步、t=0.112625；第 291 步因远场极小库存尾迹的曲率投影歧义退出，rc=3。它证明交接进入了实际 NS，尚未完成终时或解决全部后段受力问题；结果已保留，不原样重跑。

无电荷、轴对称半波长液柱，供界面断裂及断裂后动力学接续开发使用。采用 64×84 均匀网格，R/h=16，从初始化起就具有所需细网格；本算例不启用 AMR。

初始扰动幅度 δ=0.954754144451598，解析颈半径为 0.0375R。2026-09-26 轴邻体积矩重构修正后的原生颈半径约 0.03749986113R（0.6000h）；旧冻结版本为 0.04190R。原 δ=0.85 算例约为 0.12856R。初始速度为零；细颈进入已有终端曲率模型区域，仍须由实际收缩步跨过 **0.4h=0.025R** 才能接受断裂事件。该初值用于加快事件接续开发，不替代增长率、网格收敛或连续体断裂时刻验证，也不是旧轨迹的重启动。

在项目根目录运行：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/run_hydro_column_breakup.ps1 -Threads 12 -Mainline
```

脚本冻结程序、字典和输入哈希到新的 `runs/hydro_column_breakup/` 子目录，默认 dt=0.0005、终时 0.2，保留原缩步机制。使用 `-EndTime 0.4` 可在启动前指定更长演化时间。OpenMP 使用 12 线程，OpenBLAS 使用 1 线程；空闲 OpenMP 工作者休眠。默认 `-ThreadBinding false` 让系统调度各线程，避免多个进程的串行主线程被默认 `close` 策略集中到同一核；专用机器需要绑定时可显式选择 `close` 或 `spread`。`-Threads 2 -EndTime 0.003` 可运行短程性能对照。

集中验收所需的步长减半与一档空间加密可通过同一脚本设置：`-RadiusCells 16 -TimeStep 0.00025`、`-RadiusCells 32 -TimeStep 0.00025`，分别生成 64×84 和 128×168 网格。解析初态、物性及其他门限保持相同，断裂截断尺度仍为 0.4h，随加密缩小；两档初始化后的实际几何也需纳入比较。完整基准已跑通；首次步长减半对照在 t=0.19896875 后因封口曲率选支歧义退出，修复后的对照推进至 t=0.199572265625，又因普通曲率精度门限拒绝，仍未完成。原生积分更新后的旧初态步长对照最后接受 697 步至 t=0.143271484375，第 698 步因转折封口的毛细力投影歧义退出；旧 R32 已在 t=0.122625 脱离，随后接受至 t=0.1332265625，第 1250 步因封口跨图段投影歧义退出。新轴邻重构版在第 162 步、t=0.081 脱离，事件后接受 182 步至 t=0.14496875，第 345 步普通曲率支撑方向检查拒绝。支撑转弯纳入已有 2h 短模板选择后，集中检查通过，[复算](../../../../runs/hydro_column_breakup/20260926_012903_near_support_turn_0cdef8/run.log)在第 162 步、t=0.081 脱离，最后接受 377 步、t=0.1532734375，事件后 215 步；越过原第 345 步阻塞，但第 378 步另一查询的局部曲率一侧弧支撑不足，rc=3。已修正终端区域选择把继承的事件模型可用性误作普通曲率解析性的缺陷：选择前清除旧标记，以当前普通曲率重定连续模型区。受控肩部的 6 格缺失普通支撑在重复选择时均被保留，库存与几何不变，实际受力及 AMR 接续检查通过。[新复算](../../../../runs/hydro_column_breakup/20260926_020927_near_provider_scope_eae334/run.log)已在第 162 步、t=0.081 脱离，最后接受 377 步、t=0.1532734375，事件后 215 步；第 378 步仍因普通曲率支撑不足退出，rc=3，未跨过原阻塞；模型区范围和误差继续单列。旧轨迹保留，不把不同初始化/程序版本混作网格收敛，也不把断裂时刻变化当作纯性能加速。

结果在运行目录的 `output/fields/`；逐步质量、事件和流动记录在 `output/series/`，分阶段计时在 `run.log`。`build/interface_event_mainline_run.txt` 指向当前主线运行。实时分析可用：

```powershell
$run = Get-Content build/interface_event_mainline_run.txt
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/analyze_hydro_column_breakup.ps1 -Run $run -InProgress
```

2026-09-25：逐面输运积分、独立供体拟合及毛细力曲率投影均已并行化。单/多线程通量、供体形状、力场和力积分对照通过，保留原选支与曲率门限、串行力汇总顺序及失败回滚。实际同初值前六步物理记录和初始/第一步几何与旧程序一致；第 2–6 步平均计算约 17.09→8.23 秒，约加速 2.08 倍。这是共享本机短程观测，后续耗时随界面变化。较慢的旧对照前 165 步记录及 71 个保存文件完全一致后已退役，结果保留。

早期完整并行版本在第 218 步、t=0.109 接受脱离，继续 14 步后因远场极小库存尾迹的曲率投影歧义退出。当前改为检查两投影点之间的实际连接弧，同初值复算已正常结束至 **t=0.2、404 步**，脱离后继续 186 步；末段一次库存舍入越界缩步后恢复，无曲率拒绝。库存相对漂移 ≤2.23e−15，几何—库存相对 L1 ≤1.18e−4，尖端较事件时回缩约 0.51618R，共耗时约 26.4 分钟。见 [唯一待办 §8.1.3](../../docs/QUANTUM_Modular_Consolidation_Plan.md)与[完整终态图](../../output/figures/20260925_190556_near_subarc_9e1501_step0404/interface_shape.png)。这是一次运行完成，空间/时间验收仍开放。

线程绑定调整的同程序六步对照完成，三份物理历程及七个初始/首步保存场完全一致；共享本机短程约 2.73 秒/步。原绑定样本约 5.94 秒/步，但辅助任务优先级不同，不能据此声明独立的绑定策略加速比。三个旧辅助轨迹保留进度，仅放宽现有主线程 CPU 绑定；新启动的本算例默认不绑定。

2026-09-26：普通区域与终端模型的衔接改为固定边界曲率节点，区域选择同时检查实际节点；减少逐压力点拟合，保留原精度门限。集中受力、电荷及 AMR 接续检查通过。[当前主复算](../../../../runs/hydro_column_breakup/20260926_024515_near_boundary_knots_cf71fa/run.log)已越过旧 t≈0.15327 的曲率支撑阻塞：第 400 步、t=0.153375 的原压力面使用明确标记的终端模型曲率，已接受其受力。当前局部时空积分深度上限仍引发缩步，完整终时 0.2 尚未达到。已启动同方法 R16、初始 dt=0.00025 的单条时间对照，尚未启动新 R32。

</details>
