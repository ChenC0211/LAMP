> **状态：historical/reference（历史资料，仅供追溯）。**
> 本文原始正文保留为历史记录，其中“当前”“下一步”、TODO、旧规则与通过记录均指写作当时；不再自动产生新任务，也不能证明 V2 gate 通过。
> 现行任务唯一入口：[分模块任务清单](../roadmap/05_LAMP_Module_Task_Roadmap.md)。文档权威见[开发总纲](../00_LAMP_Development_Charter.md)，完整入口见[文档索引](../README.md)。
> 旧字典/职责说明可供查询既有实现；V2 新增或修改仍遵守现行六份规范。历史正文中的链接保留原样，目标缺失时不能作为本机已验证证据。
> 本目录中的旧设计、TeX/PDF 和索引同属历史参考；总纲指定的三份冻结方法 PDF 仍以 `docs/` 下的权威文件为准，不因本目录的历史标记而降级。

# Hydro / EHD 标准算例使用指南

更新：2026-09-27。本页说明当前迁移工作区的用法；输入能力以[字典参考](QUANTUM_Dictionary_Reference.md)为准，模块化状态以[主计划](QUANTUM_Modular_Consolidation_Plan.md)为准。原机算例及长程结果没有完整迁入。

## 1. 选择算例

可运行输入统一见 [cases 索引](../cases/README.md)。普通 Hydro 可从 `hydro_interface_smoke` 开始，EHD 可用 `ehd_diffuse_smoke`；喷雾有 `spray_explicit_smoke`、`spray_implicit_smoke` 和 `spray_dielectric_smoke`。独立 Interface、Electricity、Core 网格以及显隐 PIC 也在该索引中列出。

这些是本机重建的小算例，不能证明原喷嘴的长期演化或网格收敛。固定状态喷雾验证复用同一份介电输入，按其 [README](../cases/spray_dielectric_smoke/README.md)切换组合，不复制另一套求解器。

### 新 CCH 中点流体的直接运行

从尚未触发事件的连通状态开始，使用 [hydro_moment_neck_thinning_rz](../cases/hydro_moment_neck_thinning_rz/README.md)。它从静止流体、颈半径 0.008 开始；N16 上半径尺度为 0.00625，且实际事件还须满足独立库存、正则封口、动量和能量条件。每步 pending/ready 原因写入 `hydro_moment_necks.csv`。该 README 给出完整运行命令。§36 的直接参数帽面版本已完成全部 80 步，含一次模型分裂及 23 步滴后运动；物理初态和字典未改，数值精度仍须尺度核对。方法与实际结果见流体待办 §36。

近断裂的新输入见 [hydro_moment_neck_rz](../cases/hydro_moment_neck_rz/README.md)：`neckedSpheroid` 初态与 `Moment.neckRadiusCells` 显式开启局部颈部模型，已自动 1→2 后完成 20 步。事件的动量、库存与模型能量同步接续；事件 CSV 单列模型尺度和耗散，轮廓的红色帽面来源见[实际图](../cases/module_history/output/20260928_184817_build_archive/moment_neck.png)。它在 T=0 已满足触发条件，不是自然断裂时间的预测。该 README 提供完整编译/运行命令。

该例的同版空间/时间对照现已完成，空间敏感性仍明显，不能作为网格独立的液滴运动结果。半径阈值也不能单独保证四次帽面可行；方法限制和下一实施项见[流体 §34](QUANTUM_CCH_Hydrodynamics_TODO.md#34-断口的尺度影响与四次封口的可行域2026-09-28)。

多液滴的新输入见 [hydro_moment_two_drops_rz](../cases/hydro_moment_two_drops_rz/README.md)：两个液滴共用原生流场、分别守恒，已完成 20 步；每步独立液量写入 `hydro_moment_components.csv`。该说明也链接“附着液帽＋独立液滴”共同供液的 10 步运行。二者均为给定多分量初态，尚不包含颈部断裂事件。

新矩空间显式选择 `timeScheme cchMomentMidpoint`。单相例见 [hydro_moment_xy](../cases/hydro_moment_xy/README.md) 与 [hydro_moment_rz](../cases/hydro_moment_rz/README.md)，均直接使用模块的非线性推进器。[hydro_moment_capillary](../cases/hydro_moment_capillary/README.md) 接均匀闭域 XY 的真实两相惯性、共享质量通量、唯一 Hermite 界面和中点相容毛细力；冻结运行完成 200 步至 T=0.4，[RZ 参数椭球](../cases/hydro_moment_capillary_rz/README.md)完成原 120 步。闭域已接联合动态细化，见[对应输入](../cases/hydro_moment_adaptive_rz/README.md)；[附着两相 RZ 供液](../cases/hydro_moment_supply_rz/README.md)已接真实材料、方向质量矩和毛细边界功，完成 20 步。各例都有编译/运行命令、时层及输出说明。这些运行不构成一般收敛或破碎验收；自由面离域、开放域动态网格和拓扑仍待完成。下述原面平均入口保留其历史范围，勿互用验证结论。

原面平均流体的参考入口为 [hydro_midpoint_rz](../cases/hydro_midpoint_rz/README.md)：将 `system` 复制到新运行目录，依次执行 `bin/QUANTUM.exe --case <目录> --check-case` 和 `bin/QUANTUM.exe --case <目录>`。该历史入口为 16² 轴对称涡流、40 步至 0.4，显式选择 `timeScheme cchMidpoint`。当前新方法的独立流体入口使用上面的 `hydro_moment_xy/rz`，其后续两相与破碎保持同一 Compatible-Moment 路线。

看 `output/series/hydro_energy.csv` 中的动能和累计耗散；速度帧属于步终点，压力属于步中点，实际时刻同时写入 CSV 和场元数据。支持范围目前为封闭均匀 XY/RZ 的恒定单相物性；参数两相、开口边界、AMR 和带电仍另行推进。设置及输出含义见[字典参考 §15.0](QUANTUM_Dictionary_Reference.md#150-cch-中点单相流动2026-09-27)。下方旧算例/长程研究保留各自版本语境，其当时状态不代表当前方法完成状态。

### 新守恒 Interface 的直接运行

构建后可复用原平移驱动，不必手工新增 case：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/run_interface_geometry.ps1 -Method conservative -Problems translationRZ -Resolution 32 -EndTime 0.003 -Fields
```

脚本生成新的运行目录和 `method conservative` 字典。库存/几何误差、曲率质量和耗时分别写入 CSV，初末态另有 C 场及几何。`completed` 仅表示指定区间运行完毕，数值验收仍未完成。N=32 在此 RZ 问题对应 `64×128` 网格；XY 平移对应 `128×128`。参数及输出定义见[字典参考 6.3](QUANTUM_Dictionary_Reference.md#63-保守界面入口2026-09-21能力仍在收敛)。此前 t=0.03 的固定网格涡旋记录属于历史版本；当前自适应轨迹与真实误差以[实施计划](QUANTUM_Modular_Consolidation_Plan.md#812-a-包保守界面基线)为准。生产喷嘴仍用原字典路径。

2026-09-24：保守涡旋已增加可选原生自适应入口。下面是实际完成过的短程接续配置，未代表完整回逆验收：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/run_interface_advection.ps1 -Method conservative -Resolution 32 -AdaptiveMaxLevel 2 -CellBudget 12000 -EndTime 0.03125 -Fields
```

初态采用显式解析三次 jet，时间步按配置最细层级固定。每步从已接受状态决定是否细化；达到层级/单元上限或迁移预算不足时，保留旧状态并明确停止。输入、V2 原生几何存档和回归库存汇总指标的边界见[字典参考 6.4](QUANTUM_Dictionary_Reference.md#64-保守涡旋的原生自适应入口2026-09-24)。

2026-09-25 启动的完整回逆研究使用下面配置，**目前仍在运行，尚未验收**：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/run_interface_advection.ps1 -Method conservative -Resolution 64 -AdaptiveMaxLevel 5 -CellBudget 250000 -TimeStepBasis currentLeaves -ReserveBandCells 4 -EndTime 2 -Fields
```

它按实际最细叶选步长，仅在必需网格发生变化时增加四个目标单元宽度的预留，减少支撑带的小幅连续调整。省略 `ReserveBandCells` 保持原行为；预算不足时先撤掉可选预留，必需网格仍无法容纳才停止。完整周期可能耗时较长，每批冻结程序以便持续观察同一版本，状态与来源见实施计划。

### 正式静滴 B4

需要失败步缩步重试时，可在保守 Hydro 或 `ConservativeColumn` 带电液柱的 `run.dict/Time` 中显式添加 `StepRetry { minDt 0.000125; maxRetries 3; }`，并保证 `minDt <= dt`。失败候选恢复后重算整个区间，默认沿用较小 dt；可加 `recoverAfter 8;`，在连续 8 个无需重试的接受步后，把下一步候选 dt 增加 25%，最多恢复到初始 dt，所有物理门限保留。末步截短到 endTime；省略整个块仍保持固定时间日程。适用范围、真实输出时刻及限制见[字典参考第 15 节](QUANTUM_Dictionary_Reference.md#15-普通-hydroehd-分文件输入2026-09-19)。现有正式 B4 脚本仍生成原固定 dt 研究，不自动把旧网格趋势研究改成变步长研究。

2026-09-24 修正后的脚本按原设计生成输入：`R=σ=ρl=1`、`Oh=0.1`，非轴边界距滴心至少 `4R`，零初速与 `p=σ κ C` 初压。`Resolution` 表示 **8R 轴向长度上的单元数**，只取 64/128/256，对应 `R/h=8/16/32`；XY 网格为 N×N，RZ 为 (N/2)×N。旧脚本 `-B4 -Resolution 32/64/128` 的近边界、密度比 10 结果保留为内部对照，不属于这组正式定义。

```powershell
# 粗档 RZ，默认计算至少五个毛细时间，保存稀疏空间帧及每步收支
powershell -NoProfile -ExecutionPolicy Bypass -File tools/run_conservative_static_drop.ps1 -B4 -Coordinate RZ -Resolution 64 -WriteInterval 100
# 大密度比输入准备：只生成字典、冻结程序和研究定义，不启动求解
powershell -NoProfile -ExecutionPolicy Bypass -File tools/run_conservative_static_drop.ps1 -B4 -Coordinate both -Resolution 128 -DensityRatio 1000 -ViscosityRatio 0.01 -PrepareOnly
```

先做密度比 1，再做 1000；气液黏度比 `ViscosityRatio` 的 1 与 0.01 分开记录。每次调用建立新目录，不覆盖旧输出。默认步长按 `h^1.5` 缩放；显式 `-Steps 3` 仅用于短时诊断，不能作为完整 B4。`-Detail` 可观察未耦合保守 Hydro 各阶段耗时，普通运行保留 summary。

`study-input.json` 记录物性、尺度、时间和输入/程序哈希。完成后的 `results.json` 汇总全历程峰值寄生速度、库存收支、排除 2h 界面带后的两相压差，以及末态形状和动能。形状来自真实接受界面的端点/中点采样，采样 Linf 不是连续曲线的严格上界；动能采用单元物性、物理体积和四角节点速度平方的数值积分。`completed` 与数值验收分开，完整时长、最细档压差 `<1%`、无量纲峰值速度 `<1e-3` 和网格下降趋势仍需由实际结果证明。

## 2. 编译和运行

在项目根目录执行：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/build.ps1
```

脚本使用本机已有 SuperLU/BLAS；路径在 `CMakeUserPresets.json`。不需要安装依赖。默认只做增量构建；整包完成后运行相关检查，阶段验收才加 `-Test` 运行全部已注册 CTest。迁移来的 `.venv` 不用于构建或检查。

复制所选 `system/` 到新运行目录；下面命令使用现有的 schema 2 导电小算例：

```powershell
$caseName = 'spray_explicit_smoke'
$caseRun = Join-Path (Get-Location) 'runs/nozzle_001'
if (Test-Path -LiteralPath $caseRun) { throw '请换一个尚不存在的运行目录' }
New-Item -ItemType Directory -Path $caseRun | Out-Null
Copy-Item -LiteralPath "./cases/$caseName/system" -Destination $caseRun -Recurse
$env:OMP_NUM_THREADS = '2'
$env:OMP_THREAD_LIMIT = '2'
$env:OMP_DYNAMIC = 'false'
$env:OPENBLAS_NUM_THREADS = '1'
./bin/QUANTUM.exe --case $caseRun --check-case
if ($LASTEXITCODE -ne 0) { throw '输入检查失败' }
./bin/QUANTUM.exe --case $caseRun *> (Join-Path $caseRun 'run.log')
$runExit = $LASTEXITCODE
Set-Content -LiteralPath (Join-Path $caseRun 'exit-code.txt') -Value $runExit
if ($runExit -ne 0) { throw "运行失败，退出码 $runExit；保留日志排查" }
```

`--case` 指向包含 `system/` 的目录，`--check-case` 只检查输入，不能替代实际运行。普通 Hydro/EHD、正式 Interface/Electricity、网格验证和 schema 2 喷雾支持预检查；旧喷雾和 PIC 不能照搬这一步。线程设置只影响当前进程环境，PIC 自身的线程设置仍按原输入生效。

### 电喷雾启动与覆盖

使用新目录可保留已有结果。schema 2 喷雾在已有输出时需要交互选择，或显式 `--start new` 清空该算例的整个 `output/` 后重算；`system/`、`cache/` 不清理。日志放运行目录根部。`--check-case` 不询问、不清理，`--start resume` 目前明确报不支持。独立 Interface/Electricity 的新入口要求新输出目录；PIC 的检查点规则另见字典参考，不混用喷雾启动选项。

## 3. 时间设置和输出

schema 2 喷雾及普通 Hydro/EHD 的分文件格式在 `run.dict` 中设置 `Time/Output/Log`；旧普通 Hydro/EHD 扁平输入保留 `hydro.dict` 中的运行设置。两种格式不混填；其他 App 按字典参考各自章节选取。

| 位置 | 内容 |
| --- | --- |
| `output/fields` | 空间场、界面、力及按需的原生明细；喷雾默认每帧为 flow/interface/forces DAT |
| `output/series` | 随步记录的体积、电荷、流量等收支及 AMR 事件；喷雾采用 DAT，普通 Hydro/EHD 与 PIC 保留各自 CSV |
| `output/meta` | 配置来源、单位、字段说明和已写完帧的文件清单；不是独立的物理成功证明 |
| `cache` | 可恢复或重建的数据；显示用 DAT/CSV 不能代替检查点 |

默认 `physical` 场输出省去重复电场网格和原生面明细；需要它们时选 `diagnostic`，不通过日志级别控制。关闭空间场不关闭必要收支检查。不默认同时写 DAT/PLT 或 SI/代码单位两份副本。PIC 验证和移动边界的默认输出同样按上述用途分类，显式自定义路径保留。

### 逐步收支用 XY Line 查看

在 Tecplot 中单独加载需要的 `output/series/*.dat`，选 XY Line，X 为 `time_s`，Y 按变量名选取，例如 `max_speed_mps`、`relative_charge_closure` 或 `liquid_volume_m3`。不同量纲分开查看；按变量名匹配，不按列位置合并。

`*_stepNNNNNN_boundaries.dat` 是单帧边界分组表，X 为 `boundary_id`，不能将边界编号当时间。AMR 表只有实际事件；仅有表头表示没有事件。`predictor_relative_balance_defined=0` 表示比值未定义，对应显示占位零不是零误差。

原 `tools/post/electrospray_history.py` 未迁入；本机直接读取原生 DAT/CSV，不把缺失的后处理命令列为必需步骤。配置快照、单张表的末行及帧清单均不能单独代替实际退出码和守恒结果。

## 4. 在 Tecplot 中查看连续帧

按类型加载同一系列的全部 `*_flow.dat`，使用文件中的 `SOLUTIONTIME` 和 `STRANDID`。喷雾坐标为 `r_m,z_m`，压力为 `p_Pa`；普通域按文件变量名及单位使用，不能擅自标成 SI。interface 和 forces 可分别加载后叠加，需要合并时按变量名匹配。

界面折线及节点插值用于显示；不要把线性零等值线当精确自由界面。界面法向/曲率是原面元中点值；交错速度面上的力密度分量也不是同点牵引矢量。普通 EHD 的电力时刻由 `forceTime` 标明，按实际时间对齐。本包未验证 Tecplot GUI 动画。

## 5. 单位和物理限制

代码时间乘所选输入的 T0 才是秒；不要照抄另一算例的换算常数。RZ 积分已含完整回转测度，不重复乘 2πr；电荷 C 与电流 A 不混用，电极导电电流不等于发射电流。支持的几何、AMR 和接触边界以实际 App 能力为准。

小算例跑完不代表形成成熟射流。接触点启动、网格敏感性、长期喷射及完整发射闭合仍有未完成项，见[物理进度与限制](QUANTUM_Hydro_EHD_Progress.md)。

<details>
<summary>迁移前使用记录（2026-09-11 至 12，仅供追溯）</summary>

以下保留当时的数值、失败和命令；其中的原机算例、脚本和结果可能未迁入，不是当前运行指南。

## 6. 输出与控制台改进的运行记录（2026-09-11）

目录及控制台新版五个入口均实际完成，退出码0；29个DAT均通过本机Tecplot 2017 R3的preplot读取转换。没有打开GUI验收动画。以下目录中的PLT直接与DAT同放fields，转换日志放meta，不再建立格式专用目录。

| 入口 | 完成步数 / 输出时间 | 可直接查看的场 |
|---|---|---|
| 单相供液 | 20步；代码时间 .1/.2，2帧 | DAT / PLT（未迁入；原路径 `../runs/history/organization/output_layout_console_20260911/hydrodynamics_inlet_outlet_planar/output/fields/`） |
| 毛细静滴 | 20步；0/.01/.02，3帧 | DAT / PLT（未迁入；原路径 `../runs/history/organization/output_layout_console_20260911/hydrodynamics_capillary_planar/output/fields/`） |
| 平面运动液面 | 100步；0/.02/.04/.06/.08/.1，6帧 | DAT / PLT（未迁入；原路径 `../runs/history/organization/output_layout_console_20260911/hydrodynamics_levelset_planar/output/fields/`） |
| sharp介电滴 | 20步；流场0/.01/.02，电力0/.009/.019，6帧 | DAT / PLT（未迁入；原路径 `../runs/history/organization/output_layout_console_20260911/ehd_console_fixed/output/fields/`） |
| 原喷嘴运动介电短程 | 20步；step0/1/10/20各3类，12帧 | DAT / PLT（未迁入；原路径 `../runs/history/organization/output_layout_console_20260911/electrospray_moving_dielectric_short_rz/output/fields/`） |

本轮沿用SuperLU，目录/日志构建成功（未迁入；原路径 `../runs/history/organization/output_layout_console_20260911/build.log`）。首次五例的存档程序为 QUANTUM_output_layout.exe（未迁入；原路径 `../runs/history/organization/output_layout_console_20260911/QUANTUM_output_layout.exe`），SHA256：`901AA66F6BBB55145669C868E94CA6584C78BE9C8B32622E5C6BD6FD235336D6`。

实跑发现EHD结束提示误读内部临时标志，显示−1步；实际20步计算和数据正常，错误日志保留。修正后在上表ehd_console_fixed目录用原输入重跑20步、退出0，30份数据文件逐字节相同，6个DAT再次通过preplot。本批结束时bin/QUANTUM.exe及 QUANTUM_console_fixed.exe（未迁入；原路径 `../runs/history/organization/output_layout_console_20260911/QUANTUM_console_fixed.exe`） SHA256为 `408AA3699A58688CBA2BAA0E688881CCC9E600E583EDD1C71A55B8ED8CB61F43`；修正构建日志（未迁入；原路径 `../runs/history/organization/output_layout_console_20260911/build_console_fix.log`）。所有运行均已结束。

目录迁移前后：喷嘴66份原CSV及12份DAT、四个平面入口的44份原数据文件均逐字节相同。现有分析入口已支持新/历史路径和POINT/BLOCK；14份新流场通过原读场检查，静滴通过原压力响应检查，见读场日志（未迁入；原路径 `../runs/history/organization/output_layout_console_20260911/field_reader.log`）。**本轮只改输出组织与日志，没有修复物理问题；24.73%的体积缺陷仍在。**

下一步回到原喷嘴局部体积守恒，之后接导电液面电荷；不以扩大无关测试替代主线。

## 7. 公共坐标修复后的原喷嘴记录（2026-09-12）

同一七份输入、同一 SuperLU 后端和 20 步设置再次运行完成，退出 0。实际喷唇与钉扎点严格一致；初始及第 1 步的局部体积诊断均无负体积或超填。程序仍为 `movingDielectric`、q=0，尚未推进导电电荷。

- 本次 Tecplot DAT（未迁入；原路径 `../runs/history/electrospray/electrospray_shared_grid_20260912/moving_probe/output/fields/`）：step 0/1/10/20 的 flow / interface / forces，目录和字段沿用上文约定；本批未另做 preplot 转换或 GUI 动画检查。
- 逐步诊断（未迁入；原路径 `../runs/history/electrospray/electrospray_shared_grid_20260912/moving_probe/output/series/`）与运行日志（未迁入；原路径 `../runs/history/electrospray/electrospray_shared_grid_20260912/moving_probe/run.log`）：终步体积缺陷 `−2.1858828601374292e−27 m³`，相对累计净液体输运量 **24.726432%**，仍未通过守恒验收。
- 本次运行程序为 `runs/history/electrospray/electrospray_shared_grid_20260912/QUANTUM_probe.exe`，仅在初始和第 1 步增加只读局部体积输出。标准 App 中探针已移除，构建及七项直接相关组件检查通过；未另跑移除探针后的同一短程。

本批有意修复原生网格坐标的末位，因此不承诺与第 6 节流场逐字节一致，也不能把很小的误差变化解释为守恒问题已解决。旧输出保留。不一致旧 QMC2 缓存现在会明确拒绝，需从 case 重建，不能直接沿用旧几何状态；本批未删除任何用户缓存。详细修改和限制见[主计划第 78 节](QUANTUM_Electrohydrodynamics_Implementation_Plan.md)。

## 8. 当前开发版：局部体积修正的原样复跑（2026-09-12）

启用原喷嘴局部体积修正的最新两阶段迭代原样复跑（未迁入；原路径 `../runs/history/electrospray/electrospray_volume_correction_20260912/moving_feasibility_phase/run.log`）已完成 **20 步、退出 0**。七份输入字典与永久短程 case、第 7 节旧运行均相同。每步 3–4 次修正、零回退，沿用 SuperLU；第 6、7 节未修正旧结果不覆盖。

本次可查看 step 0/1/10/20 的 12 份 DAT / PLT（未迁入；原路径 `../runs/history/electrospray/electrospray_volume_correction_20260912/moving_feasibility_phase/output/fields/`），全部经过 preplot 转换，日志在 `meta`。四份逐步 CSV 各 20 行，帧清单与终态一致；4 个流场均为 44001 节点 / 43520 活动单元，全部 12 个 DAT 字段有限。尚未打开 GUI 检查动画。前一停止顺序版本第 2 步失败及其 6 个首步文件保留，不混作完整记录。

沿用原全局体积诊断，终步缺陷从 `−2.1858828601374292e−27` 改善为 `+1.1692559034358223e−32 m³`，占累计净液体输运量的比例从 **24.726432%** 降至 **0.000132265%**；本次全过程该百分比最大 **0.0516337%**（第 1 步）。局部 L1 相对收支最大 `8.59315e−5`，最大单元容量相对残差 `1.33321e−13`，满足本批预定工作准则。末步液体在供液端仍净流出，这是实际短时响应，不是已经建立供液射流；压力边界和原 dt 没有为得到喷流而修改。

新版本增加 `output/series/*_local_volume_diagnostics.csv` 与 `output/fields/*_volume_balance_cells.csv`，分别记录修正收敛过程和逐单元 SI 收支；仍保留原全局体积与 Tecplot 输出。逐轮控制台会说明候选计算、回退和残差。“接受修正候选”不是“完成物理时间步”；失败时可能留下部分候选场，应一起核对退出码、完整逐步 CSV 和帧清单。修正细节、实际构建和失败证据见本批记录（未迁入；原路径 `validation/electrospray/QUANTUM_Electrospray_Volume_Correction_20260912.md`）。

## 9. 原喷嘴导电运动液面：20 步实跑（2026-09-12）

新入口 electrospray_conducting_short_rz（未迁入；原路径 `../cases/electrospray_conducting_short_rz/README.md`） 的七份输入与本次正常退出的运行（未迁入；原路径 `../runs/history/electrospray/electrospray_boundary_chart_20260912/original_conducting_twenty_steps/run.log`）逐项相同。使用原物性、45.2 Pa 供液、300 V 喷嘴、0 V 顶部以及 SuperLU，显式选择 `coupled / sharpHermiteP3`。从项目根目录运行：

```powershell
# 先按第 2 节将 electrospray_conducting_short_rz/system 复制到新的 $caseRun。
./bin/QUANTUM.exe --case $caseRun
```

本次已完成 20 步，运行循环约 619 s。fields 中的 16 份 DAT / PLT（未迁入；原路径 `../runs/history/electrospray/electrospray_boundary_chart_20260912/original_conducting_twenty_steps/output/fields/`）包括 step 0/1/10/20 的 flow、interface、electric、forces，各 preplot 转换退出 0；原生定向面电流 CSV 保留。四类场按变量名合并，不按列位置混合。

导电输出新增 `output/series/*_conducting_charge_diagnostics.csv`，逐步记录实际总面电荷、旧/新电流、守恒残差及重映射。四帧场积分与全部 20 步电荷历史的复核报告（未迁入；原路径 `../runs/history/electrospray/electrospray_boundary_chart_20260912/history_review_twenty_steps.json`）通过：终步 Q 约 `3.36715e−16 C`、边界流入电流 `18.3464 µA`；最大相对电荷收支 `3.23294e−13`，原门 `1.6e−9` 不变。体积缺陷相对累计净液体输运量终步 `0.000152247%`、全过程最大 `0.0366802%`。这是短时启动，供液口仍回吐，不是形成射流的结论。

当前帧电流与载荷供下一步使用，不把它们与当前压力视作同时层静力平衡。界面文件的面元电荷在显示线段上重复，须按 `panel_id` 去重；载荷文件另给出投影点电荷和两侧重构电场。源单元及边界相对坐标说明、修复范围和程序哈希见本批记录（未迁入；原路径 `validation/electrospray/QUANTUM_Electrospray_Boundary_Chart_20260912.md`）。

</details>
