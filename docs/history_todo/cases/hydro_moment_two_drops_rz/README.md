此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_moment_two_drops_rz/README.md。

运行入口：tools\run_case.cmd hydro_moment_two_drops_rz 。结果：runs/hydro_moment_two_drops_rz/<时间戳_标识>/。程序：LAMP。

# RZ 双液滴毛细运动

<!-- case-definition-20260928 -->

**算例类别：** 多分量功能示例。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

两个独立椭球液滴共用流场和相材料积分，分别记录分量身份与液体库存。

## 计算域、初始条件与边界条件

封闭域，r=0 为轴线。初始即为两个分量；本例不能作为动态分裂或合并已经实现的证据。

```text
Domain { min (0, 0); max (1, 1); dimension 2; }
```

```text
ParametricInterface
{
    name liquid;
    type components;
    components
    {
        lower { type spheroid; center (0, 0.42); radii (0.18, 0.12); segments 16; }
        upper { type spheroid; center (0, 0.67); radii (0.12, 0.09); segments 16; }
    }
}
```

初始速度类型为 `zero`；具体初场幅值与初始压力设置见下表。

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_moment_two_drops_rz/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_moment_two_drops_rz/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/hydro_moment_two_drops_rz/system/hydro.dict) |
| 材料模型 | `twoPhase` | [hydro.dict](../../../../cases/hydro_moment_two_drops_rz/system/hydro.dict) |
| 负相密度 ρ− | `5` | [hydro.dict](../../../../cases/hydro_moment_two_drops_rz/system/hydro.dict) |
| 正相密度 ρ+ | `1` | [hydro.dict](../../../../cases/hydro_moment_two_drops_rz/system/hydro.dict) |
| 负相动力黏度 μ− | `0.02` | [hydro.dict](../../../../cases/hydro_moment_two_drops_rz/system/hydro.dict) |
| 正相动力黏度 μ+ | `0.01` | [hydro.dict](../../../../cases/hydro_moment_two_drops_rz/system/hydro.dict) |
| 表面张力系数 σ | `0.1` | [hydro.dict](../../../../cases/hydro_moment_two_drops_rz/system/hydro.dict) |
| 时间格式 | `cchMomentMidpoint` | [hydro.dict](../../../../cases/hydro_moment_two_drops_rz/system/hydro.dict) |
| 对流格式 | `central` | [hydro.dict](../../../../cases/hydro_moment_two_drops_rz/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-10` | [hydro.dict](../../../../cases/hydro_moment_two_drops_rz/system/hydro.dict) |
| 背景网格划分 | `(4, 4)` | [mesh.dict](../../../../cases/hydro_moment_two_drops_rz/system/mesh.dict) |
| 启用网格自适应 | `false` | [mesh.dict](../../../../cases/hydro_moment_two_drops_rz/system/mesh.dict) |
| 最大单元数 | `4096` | [mesh.dict](../../../../cases/hydro_moment_two_drops_rz/system/mesh.dict) |
| 时间步长 Δt | `0.002` | [run.dict](../../../../cases/hydro_moment_two_drops_rz/system/run.dict) |
| 计算终时 T | `0.04` | [run.dict](../../../../cases/hydro_moment_two_drops_rz/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_moment_two_drops_rz/system/run.dict) |
| 场输出间隔 | `5` | [run.dict](../../../../cases/hydro_moment_two_drops_rz/system/run.dict) |
| 界面输运方法 | `parametricGals` | [interface.dict](../../../../cases/hydro_moment_two_drops_rz/system/interface.dict) |
| 重新初始化间隔 | `0` | [interface.dict](../../../../cases/hydro_moment_two_drops_rz/system/interface.dict) |

完整输入文件：[geometry.dict](../../../../cases/hydro_moment_two_drops_rz/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_moment_two_drops_rz/system/hydro.dict)、[interface.dict](../../../../cases/hydro_moment_two_drops_rz/system/interface.dict)、[mesh.dict](../../../../cases/hydro_moment_two_drops_rz/system/mesh.dict)、[run.dict](../../../../cases/hydro_moment_two_drops_rz/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_moment_two_drops_rz
```

每次运行创建 `runs/hydro_moment_two_drops_rz/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# 两个轴对称液滴的共同流体演化

本例为新矩方法的多分量入口，已完成原定 20 步。两个不同大小的椭球液滴共用一套原生 RZ 流体速度与压力；每个液滴保留独立守恒库存，材料密度/黏性由两者的实际液相区域共同积分。表面张力 0.1，液/气密度 5/1、黏度 0.02/0.01，初速为零。原生 4×4 网格上有单元同时包含两个液体分量。

`ParametricInterface.type components` 的子块按输入顺序取得稳定分量编号。每个子块使用既有 `ellipse`、`spheroid` 或 `supplyCap` 初态定义；多分量必须互不接触且互不包含。所有分量仍属于同一种液体，`interfaceLevelSet liquid` 是整个液相的名字。不同液体材料、嵌套气泡、合并和颈部断裂事件不由这个初态入口自动实现。

`hydro_moment_components.csv` 每步给出各分量的独立液量、相对初态变化和几何体积；轮廓 CSV 与 Tecplot 自由面增加 `component`。单元液相分数是各分量的总和，不可把同一个单元误认作只能含一个液滴。流体能量仍在原 `hydro_moment_energy.csv` 中记账。

从根目录构建并复制到新运行目录：

```powershell
cmake --build build/windows --config Release --target QUANTUM --parallel 1 -- /nr:false
$caseRun='runs/hydro_moment/my_two_drops_rz'
if (Test-Path -LiteralPath $caseRun) { throw '请换一个尚不存在的运行目录' }
New-Item -ItemType Directory -Path $caseRun
Copy-Item -LiteralPath cases/hydro_moment_two_drops_rz/system -Destination $caseRun -Recurse
$env:OMP_NUM_THREADS='4'
$env:OMP_DYNAMIC='FALSE'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
.\bin\LAMP.exe --case $caseRun --check-case
if ($LASTEXITCODE -ne 0) { throw '算例预检失败' }
.\bin\LAMP.exe --case $caseRun
```

输入使用 dt=0.002、20 步至 T=0.04；粗网格仅用于多分量演化和守恒接续，不作为曲率精度设置。2026-09-28 的实际运行在 `runs/hydro_moment/20260928_two_drops_rz`，预检/运行退出 0。两滴液量最大漂移分别为 3.47e−18 和 8.68e−19；轴向/径向尺度比分别从 0.6667/0.7500 增至 0.6761/0.7736，逐渐回圆。最大方程残差 8.43e−11，包含实际表面能变化的总能量账差约 7.80e−9。

[实际形状与能量图](../../cases/module_history/output/20260928_184817_build_archive/moment_components.png) · [方法、证据及下一步](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#32-多液滴的独立库存与共同流体演化2026-09-28)。初始化两滴不代表已完成颈部断裂；分裂时的子滴库存和身份分配仍是下一主项。

另有“附着液帽＋独立液滴”供液变体，复用 `hydro_moment_supply_rz`，只改几何初态与终时；实际字典见 [geometry.dict](../../../../runs/hydro_moment/20260928_components_supply_rz/system/geometry.dict)，生成过程见 [脚本](../../build/run_moment_components_supply.ps1)。它已完成 10 步：附着液体增长 3.3333%，与入口累计量的差不超过 2.17e−17，独立滴库存保持原值。该原始运行目录已存在，重跑时需在脚本副本中另选未使用的目录。

</details>
