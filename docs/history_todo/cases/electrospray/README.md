此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/electrospray/README.md。

运行入口：tools\run_case.cmd electrospray 。结果：runs/electrospray/<时间戳_标识>/。程序：LAMP。

# 实际喷嘴电喷雾

<!-- case-definition-20260928 -->

**算例类别：** 工程应用模板。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

轴对称不可压缩两相流、电准静态场与漏电介质表面电荷输运；用于喷嘴供液、液面变形及最终带电脱落。

## 计算域、初始条件与边界条件

压力供液 45.2 Pa；电极 300 V，收集极 0 V；喷口唇缘接触线钉扎，固壁无滑移；环境气体允许回流。新界面/流体主线完成前，带电计算暂缓。

```text
Domain
{
    min (0, 0) um
    max (5, 16) um
}
```

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `Electrospray` | [run.dict](../../../../cases/electrospray/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/electrospray/system/hydro.dict) |
| 背景网格划分 | `(20, 40)` | [mesh.dict](../../../../cases/electrospray/system/mesh.dict) |
| 启用网格自适应 | `true` | [mesh.dict](../../../../cases/electrospray/system/mesh.dict) |
| 最大单元数 | `200000` | [mesh.dict](../../../../cases/electrospray/system/mesh.dict) |
| 时间步长 Δt | `5e-6` | [run.dict](../../../../cases/electrospray/system/run.dict) |
| 计算终时 T | `20.0` | [run.dict](../../../../cases/electrospray/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/electrospray/system/run.dict) |
| 场输出间隔 | `10` | [run.dict](../../../../cases/electrospray/system/run.dict) |
| 界面输运方法 | `activeHermiteMaterial` | [interface.dict](../../../../cases/electrospray/system/interface.dict) |

两相物性与电学材料常数：

```text
Materials
{
    Phases
    {
        Liquid
        {
            phaseId 2
            density 1240 kg/m^3
            dynamicViscosity 0.038 Pa*s
            relativePermittivity 10
            conductivity 1.4 S/m
        }
        Gas
        {
            phaseId 3
            density 1.184 kg/m^3
            dynamicViscosity 1.849e-5 Pa*s
            relativePermittivity 1
            conductivity 0 S/m
        }
    }
    Interfaces
    {
        LiquidGas
        {
            phases (2 3)
            surfaceTension 0.0452 N/m
        }
    }
}
```

完整输入文件：[electricity.dict](../../../../cases/electrospray/system/electricity.dict)、[geometry.dict](../../../../cases/electrospray/system/geometry.dict)、[hydro.dict](../../../../cases/electrospray/system/hydro.dict)、[interface.dict](../../../../cases/electrospray/system/interface.dict)、[materials.dict](../../../../cases/electrospray/system/materials.dict)、[mesh.dict](../../../../cases/electrospray/system/mesh.dict)、[normalization.dict](../../../../cases/electrospray/system/normalization.dict)、[run.dict](../../../../cases/electrospray/system/run.dict)、[solvers.dict](../../../../cases/electrospray/system/solvers.dict)、[units.dict](../../../../cases/electrospray/system/units.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd electrospray
```

每次运行创建 `runs/electrospray/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

关注电势、电场、表面电荷、Maxwell 应力、速度压力、液量及电荷收支。确认数值表面电荷、极化电荷和总电荷定义；短程正常退出不等同于形成了可信的电喷雾射流。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# 原喷嘴 Electrospray：导电液面与动态局部网格

本目录对应 [example.jpg](../../docs/assets/electrospray/example.jpg) 的喷嘴，不是旧 App 的硬编码圆液滴。按现行[字典规范](../../docs/03_QUANTUM_Dictionary_Rules.md)组织；2026-09-14 接入运行时 AMR，物性、边界数据、几何、代码时间和场输出间隔保持原设置。

**默认终时对应 4,000,000 步，不是开发短程，也不是已经完成物理验收的射流算例。** 2026-09-15 已接通本喷嘴需要的公共边界粗细网格动量通量及压力耦合：液体最低 3 级，界面等重点区域最高 4 级，删除侧面和顶部开放边界的专门加密带。原目标输入完成动态三步后应用到本目录，同配置另完成十步接续，实际退出均为 0；没有执行完整长时间计算。

## 物理设置



| 项目       | 当前输入                                                                                     |
| -------- | ---------------------------------------------------------------------------------------- |
| 几何       | 无旋流 RZ；r=0–5 µm、z=0–16 µm；喷唇 (1,3) µm，初始液面轴端 (0,4) µm                                    |
| 液体 / 气体  | 相 2：ρ=1240 kg/m³、μ=0.038 Pa・s、εr=10、σ=1.4 S/m；相 3：ρ=1.184 kg/m³、μ=1.849e−5 Pa・s、εr=1、σ=0 |
| 表面张力     | 0.0452 N/m                                                                               |
| 压力边界     | 边界 5：45.2 Pa 供液；边界 8、9：0 Pa 环境气体开放边界，允许回流；不指定供液速度或流量                                     |
| 固壁 / 接触线 | 无滑移；钉扎在几何点 6 的喷唇处，接触角可变化                                                                 |
| 电极 / 电荷  | 导体相 1：300 V；顶部边界 9：0 V；初始自由面电荷为零，随后推进导电电荷                                                |

几何中的 `vacuum` 是保留的历史相名称；本算例实际使用上表的气体物性，不是真空流动模型。初态速度、压力为零；供液压力在 NS 求解中施加。

原输入把默认导电推进与 `cellHarmonic` 初态电学选择混在一起。此次**显式改为导电路径要求的&#x20;**`sharpHermiteP3`；保留旧 coupled 的显式电荷时间推进 `explicitConductionMaterialTransport`，不在字典整理时切换隐式算法。因此这不是电离散完全不变的文字改名。

需要旧工程欠松弛电荷 — 电场迭代时，使用已有的[隐式导电短程](../electrospray_conducting_implicit_short_rz/README.md)。本模板仍保留原显式电荷模式和时间步；动态网格调度与电荷隐式迭代是独立设置，不因接 AMR 而同时切换。

## 去哪里修改



| 文件                                                                              | 唯一职责                                       |
| ------------------------------------------------------------------------------- | ------------------------------------------ |
| [run.dict](../../../../cases/electrospray/system/run.dict)                                                     | App、Time、Output、Log、Restart 和耦合更新频率        |
| [materials.dict](../../../../cases/electrospray/system/materials.dict)                                         | 两相密度、黏度、介电常数、电导率及表面张力，显式 SI 单位             |
| [hydro.dict](../../../../cases/electrospray/system/hydro.dict)                                                 | RZ 两相 NS、压力边界、界面力带宽及流体收支判据；不再内联界面块         |
| [interface.dict](../../../../cases/electrospray/system/interface.dict)                                         | 液面名称、移动边界、材料输运与钉扎接触线位置；由 Interface/IO 直接读取 |
| [electricity.dict](../../../../cases/electrospray/system/electricity.dict)                                     | 电模型、电极、初始面电荷、电荷推进及电学收支判据                   |
| [solvers.dict](../../../../cases/electrospray/system/solvers.dict)                                             | 五个实际线性系统的 SuperLU 后端及各自残差容限                |
| [geometry.dict](../../../../cases/electrospray/system/geometry.dict)、[mesh.dict](../../../../cases/electrospray/system/mesh.dict)             | 公共几何／编号；初始局部网格、运动边界加密带及运行时 AMR 调度          |
| [units.dict](../../../../cases/electrospray/system/units.dict)、[normalization.dict](../../../../cases/electrospray/system/normalization.dict) | 参考尺度及完整 RZ 回转体积测度，沿用原格式                    |

界面配置已按 schema 2 分文件格式拆到独立 `interface.dict`（液面 levelSet `"1"`、移动边界 4、`activeHermiteMaterial` 输运、点 6 钉扎），`hydro.dict` 只保留流体设置；`interface.dict` 与 `hydro.dict` 内联块两处不得同时存在，否则输入被拒绝。2026-09-20 同输入三步短程对照确认：拆分前后 12 个 Tecplot 场与 9 张逐步收支表逐字节一致，唯一差异为 AMR 表随运行变化的 `wall_seconds` 计时列；`meta/electrospray_resolved.json` 中五个界面字段来源由 `hydro.dict:Hydrodynamics.Interface.*` 变为 `interface.dict:Interface.*`。

不再使用旧 `ElectrosprayPhysical`、run 内旧时间字段或无实际消费者的 `PoissonCellCentered/useCG/useIC`。五个线性残差容限和原流体 / 电学收支判据仍为 1e−10，电荷收支门限仍为 1.6e−9，彼此独立；没有放宽验收标准。

网格保留原 20×40 背景划分，最高提高到四级，采用以下分工；运行证据与适用范围见[分辨率与边界离散记录](../../docs/validation/electrospray/QUANTUM_Electrospray_Mesh_Policy_20260914.md)。



| 区域 / 判据        | 分辨率要求                                                           |
| -------------- | --------------------------------------------------------------- |
| 远离界面的低速气体      | 最低 1 级，`Δr×Δz=0.125×0.2 µm`                                     |
| 当前液体内部         | 最低 3 级，`0.03125×0.05 µm`，随当前液面判相；两个方向均比原最低二级减半                  |
| 运动液面（含 tip）    | 边界 4，`C_Gamma=5`，最高 4 级，`0.015625×0.025 µm`                     |
| 无滑移壁面、喷口肩部、供液口 | 边界 2/3 的距离带及具名局部加密区，4 级                                         |
| 轴线             | 保留约 8 个最细单元宽的同层加密带                                              |
| 侧面 / 顶部压力开口    | 不设专门加密带；与轴线 / 壁面相交处仍受相应局部规则及 2:1 过渡约束                           |
| 高速 / 强速度梯度区    | `Velocity > 0.3 m/s` 或 `VelocityGradient > 1e6 1/s` 请求加密，最高 4 级 |

速度和速度梯度是实际已接入的指标，不是只写在注释里的开关；判据使用两相活动流场。`hysteresis=0.25` 表示降到门限的 75% 以下才允许粗化，几何最低分辨率及液面 / 边界保护仍优先。它们是本算例的初始分辨率策略，不是经网格无关性验证得到的普适最优门限。参考旧 `passive_fed` 的液体最低层级、界面距离带、壁面 / 唇缘与速度加密思想，不照搬其算例归一化常数。

`mesh.dict/RuntimeAdaptation` 当前为 `enabled true / everySteps 1 / maxSweeps 2`。每个接受步后、下一步开始前更新运动边界 4 的加密带，复用原四叉树细粗化，保持物理时刻不变；不额外推进电荷或 NS。当前液面单元和必要模板受保护，只在已确认的单相区改层级。**四级目标的三步实跑退出 0：第 1 步后总单元数 28841 → 23567，活动流体单元 21905 → 18554；AMR 的液体体积和面电荷变化均为零，后续完整耦合接续。** 共享节点支撑已在复制前纳入目标层级判定；无需调整时不生成候选网格，需要调整时直接复制原生数据、不重建背景网格。原完整 Q4 多项式和守恒检查不变。运行期再加密能力来自同一原 `adapt()`，本次短程实际发生的是粗化，不宣称长距离尖端运动或破裂已经验收。

此前有开放边界加密带的 18611 → 15773 和更早的 28358 → 25028 结果继续保留，不与当前网格混作同输入对照。关闭运行时调整时，将 RuntimeAdaptation 和 PhysicsRefinement 两块都改为只有 `enabled false`；保留几何 / 相区规则即可使用固定局部网格。频率和限制见[字典参考 5.6](../../docs/QUANTUM_Dictionary_Reference.md#56-运行时网格自适应meshdict)。本次只改正式 `mesh.dict`；物性、边界条件、时间步、终时、输出间隔及 `electrospray` 前缀不变。

## 时间、输出与运行

`Time.units code`，`dt=5e−6`，`endTime=20`，对应物理终时约 **3.31262 µs**；代码时间须乘 `T0=1.6563100047570824e−7 s`。用户已将场输出改为 `writeControl stepInterval / writeInterval 10`，即每 10 步，另保留初态、首步和末步。本次启动逻辑修改保留该设置。

输出采用 `profile physical`：Tecplot DAT 的 flow/interface/forces，位置仍为 `output/fields`；逐步收支放 `output/series`，生效配置和帧清单放 `output/meta`。`run.dict/Output/prefix` 为 `electrospray`，例如 `electrospray_conducting_step000001_flow.dat`，不再把整个过程命名成 `electrospray_init_*`；既有历史文件不改名。不生成 `tecplot/` 目录或重复 DAT/PLT。控制台 `Log.everySteps=1`、`level summary`；详细诊断按需开启，不靠关闭写盘跳过守恒检查。

逐步收支直接写 `output/series/*_diagnostics.dat`，Tecplot 中选择 XY Line，横轴选 `time_s`，纵轴可选 `max_speed_mps`、`liquid_volume_m3`、`total_charge_after_C` 或相应残差。各表分开加载，或按变量名称合并；不按变量位置混合。每行写完即可读取，运行中刷新数据不会自动改变已有图。默认不重复写 CSV；旧结果不转换或覆盖。`meta` 中的帧清单仍为 CSV。完整变量与未定义比值说明见[XY Line 指南](../../docs/QUANTUM_Hydro_EHD_Case_Guide.md#逐步收支用-xy-line-查看)。

每次实际 AMR 仅向 `output/series/<prefix>_amr.dat` 增加一行网格、体积／电荷／流量及转移误差记录；没有 AMR 事件时只有表头，不把它画成零误差曲线。后续场帧直接包含新网格，不在同一时刻多写一份空间场。不要按节点编号跨帧匹配。`*_stepNNNNNN_boundaries.dat` 则是单帧的边界分组表，横轴选 `boundary_id` 查看 `outward_J_A`；它不是逐步历程。

在 QUANTUM 根目录可以直接做只读检查：



```
.\bin\LAMP.exe --case .\cases\electrospray --check-case
```

首次求解只复制 `system/` 到新的运行目录，不在模板里累积结果；需要保留旧结果时也使用新目录：



```
\$caseRun = 'D:\Chen\QUANTUM\runs\electrospray\_001'

if (Test-Path -LiteralPath \$caseRun) { throw '请换一个尚不存在的运行目录' }

New-Item -ItemType Directory -Path \$caseRun | Out-Null

Copy-Item -LiteralPath '.\cases\electrospray\system' -Destination \$caseRun -Recurse

\# 先审核运行副本 system/run.dict 的时间范围；短程只修改副本。

.\bin\LAMP.exe --case \$caseRun --check-case

if (\$LASTEXITCODE -ne 0) { throw '配置检查失败，请先修正输入' }
```

确认运行副本的步数和输出间隔后，再单独执行：



```
.\bin\LAMP.exe --case \$caseRun 2>&1 | Tee-Object -FilePath (Join-Path \$caseRun 'run.log')

\$LASTEXITCODE
```

新版程序在交互终端会先询问：`1 从初态重新计算 / 2 接上次计算 / 0 取消`。选 `1` 清空当前算例的整个 output/ 再启动；旧结果不能由程序恢复，system/、cache/ 和算例其他文件保留。选 `2` 暂时明确提示续算尚未支持，旧结果不动；回车、输入结束或 `0` 取消启动。

脚本重复运行需显式加 `--start new` 确认覆盖，或加 `--start ask` 强制询问；未确认的非交互运行遇到旧输出会在初始化前停止。只读检查不询问也不清理。日志写在根目录 run.log，不写进待清理的 output/。完整规则见[字典参考](../../docs/QUANTUM_Dictionary_Reference.md#启动选择与覆盖旧结果)。

快速复核可直接选择已有[20 步显式导电算例](../electrospray_conducting_short_rz_v2/README.md)或[8 步隐式导电算例](../electrospray_conducting_implicit_short_rz/README.md)。构建、依赖、OpenMP 环境和 Tecplot 加载方式统一见[算例指南](../../docs/QUANTUM_Hydro_EHD_Case_Guide.md)，不在此维护另一套规范。

## 状态与来源

无开放边界加密带、界面四级的[十步 Tecplot 场](../../runs/checks/electrospray_amr_preflight_20260915/requested_ten/output/fields/)、[物理 / 网格收支](../../runs/checks/electrospray_amr_preflight_20260915/requested_ten/output/series/)、[日志](../../runs/checks/electrospray_amr_preflight_20260915/requested_ten.log)保存在隔离运行目录；已完成 10/10 步、退出 0，沿用每 10 步输出，保存 0/1/10 帧。物理终时约 8.282 ps，本次同输入约 13.5 s / 步（此前 16.3 s / 步），尚未达到加速目标或成熟喷射状态。flow 文件已有 `mesh_level` 和 `phase_id`，无需重复导出网格；物理输出与此前十步逐字节相同，可复用已有的 [Tecplot 时间曲线](../../runs/checks/electrospray_boundary_flux_20260915/requested_ten/post/series/)可辅助查看充电和流量过程。原三步失败输入回放、公共边界修复和同输入回归见[网格策略记录](../../docs/validation/electrospray/QUANTUM_Electrospray_Mesh_Policy_20260914.md)。下述为历史来源，不与当前新网格混作同输入对照。

字典迁移前输入保存在[迁移记录目录](../../runs/maintenance/electrospray_dictionary_schema2_20260914/)；AMR 修改前九份输入和 README 保存在[本批输入快照](../../runs/checks/electrospray_amr_20260914/app_runtime/production_before/)。[legacy](legacy/) 的原始几何 / 网格和示意图保留。

动态三步的[完整 Tecplot 场](../../runs/checks/electrospray_amr_20260914/app_runtime/dynamic_nozzle_fixed/output/fields/)、[逐步与 AMR 收支](../../runs/checks/electrospray_amr_20260914/app_runtime/dynamic_nozzle_fixed/output/series/)和[日志](../../runs/checks/electrospray_amr_20260914/app_runtime/dynamic_nozzle_fixed.log)可直接查看。该副本仅缩短终时并每步写场；原式的独立后处理核对了 0/1/2/3 帧及逐步记录。固定网格的接口拆分对照另有 23 份数值文件逐字节一致；失败复现、修复、输入范围和程序身份见[动态 AMR 第 8 节](../../docs/validation/electrospray/QUANTUM_Electrospray_Dynamic_AMR_20260914.md#8-第四批电喷雾时间循环的联合准备与提交)。

导电耦合的短程完成不等于整个长终时能可靠运行。曲率与局部体积修正的兼容性、运行速度、长期稳定性及网格 / 时间收敛仍需继续处理；本次不修改这些数值算法。开发进度与历史物理证据统一见[进度页](../../docs/QUANTUM_Hydro_EHD_Progress.md)和[主计划](../../docs/QUANTUM_Electrohydrodynamics_Implementation_Plan.md)。

</details>
