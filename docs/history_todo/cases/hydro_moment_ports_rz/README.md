此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_moment_ports_rz/README.md。

运行入口：tools\run_case.cmd hydro_moment_ports_rz 。结果：runs/hydro_moment_ports_rz/<时间戳_标识>/。程序：LAMP。

# RZ 单相供排流动

<!-- case-definition-20260928 -->

**算例类别：** 开放边界功能示例。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

恒物性轴对称流动，检查规定入口通量、压力出口及边界动量/功收支。

## 计算域、初始条件与边界条件

南侧外法向速度 −0.1，北侧压力 0.2；r=0 为轴线，外径向边界为无滑移壁面。初始化会满足入口与不可压缩约束。

```text
Domain { min (0, 0); max (1, 1); dimension 2; }
```

初始速度类型为 `zero`；具体初场幅值与初始压力设置见下表。

显式边界输入：

```text
boundaryPatches
    {
        supply { side south; type normalVelocity; start 0; end 1; value -0.1; }
        outlet { side north; type pressureOutlet; start 0; end 1; value 0.2; }
    }
```

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_moment_ports_rz/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_moment_ports_rz/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/hydro_moment_ports_rz/system/hydro.dict) |
| 材料模型 | `constant` | [hydro.dict](../../../../cases/hydro_moment_ports_rz/system/hydro.dict) |
| 密度 ρ | `1` | [hydro.dict](../../../../cases/hydro_moment_ports_rz/system/hydro.dict) |
| 动力黏度 μ | `0.02` | [hydro.dict](../../../../cases/hydro_moment_ports_rz/system/hydro.dict) |
| 时间格式 | `cchMomentMidpoint` | [hydro.dict](../../../../cases/hydro_moment_ports_rz/system/hydro.dict) |
| 对流格式 | `central` | [hydro.dict](../../../../cases/hydro_moment_ports_rz/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-10` | [hydro.dict](../../../../cases/hydro_moment_ports_rz/system/hydro.dict) |
| 背景网格划分 | `(8, 8)` | [mesh.dict](../../../../cases/hydro_moment_ports_rz/system/mesh.dict) |
| 启用网格自适应 | `false` | [mesh.dict](../../../../cases/hydro_moment_ports_rz/system/mesh.dict) |
| 最大单元数 | `2048` | [mesh.dict](../../../../cases/hydro_moment_ports_rz/system/mesh.dict) |
| 时间步长 Δt | `0.01` | [run.dict](../../../../cases/hydro_moment_ports_rz/system/run.dict) |
| 计算终时 T | `0.2` | [run.dict](../../../../cases/hydro_moment_ports_rz/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_moment_ports_rz/system/run.dict) |
| 场输出间隔 | `5` | [run.dict](../../../../cases/hydro_moment_ports_rz/system/run.dict) |

完整输入文件：[geometry.dict](../../../../cases/hydro_moment_ports_rz/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_moment_ports_rz/system/hydro.dict)、[mesh.dict](../../../../cases/hydro_moment_ports_rz/system/mesh.dict)、[run.dict](../../../../cases/hydro_moment_ports_rz/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_moment_ports_rz
```

每次运行创建 `runs/hydro_moment_ports_rz/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# Compatible-Moment 轴对称供液/出流

单位 RZ 圆柱域，8×8 原生叶，密度 1、黏度 0.02。下端给定向内法向速度 0.1，上端给定完整法向应力 `(-p I + 2 μ D)n = -0.2 n`；轴上规则，外圆柱壁无滑移。时间步 0.01，20 步至 0.2。它使用新矩方法原有的散度、负伴随压力、中心动量、应变能黏性和中点推进。

`normalVelocity.value` 使用外法向约定，负数表示供液；这个入口只规定法向速度，切向应力为零。进入域的动量以该法向速度、零切向速度作为给定数据。`pressureOutlet.value` 是完整应力边界，既不钉住邻近单元压力，也不同时规定出口流速。当前出口要求实际向外流动；反流需要明确的外部动量/相数据，程序会拒绝未定义的反流，不暗中复制内部值。

初始 `zero` 是约束初始化的参考场。非零供液下，它被投影到满足所有法向边界与 `D u=0` 的最小物理质量范数状态；程序打印该修正量，时间仍为零，初帧没有已求解的压力。能量账以这份真正接受的初始速度为起点。后续没有压力投影替代物理联立求解。

从根目录运行，选用尚不存在的目录：

```powershell
cmake --build build/windows --config Release --target QUANTUM --parallel 1 -- /nr:false
$caseRun='runs/hydro_moment/my_ports_rz'
New-Item -ItemType Directory -Path $caseRun
Copy-Item -LiteralPath cases/hydro_moment_ports_rz/system -Destination $caseRun -Recurse
$env:OMP_NUM_THREADS='4'
$env:OMP_DYNAMIC='FALSE'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
.\bin\LAMP.exe --case $caseRun --check-case
.\bin\LAMP.exe --case $caseRun
```

`output/series/hydro_moment_energy.csv` 除原流场收支外，还输出累计出口牵引功、入口约束反力功、净流出动能、入口动量差耗散和流入/流出体积率。入口动量差项来自给定外部速度的边界通量，不是内部黏度。应检查进出流量相等，以及动能变化、黏性与边界能量交换共同闭合；只看动能是否下降不适用于持续供液。

这是单相流体边界入口。另有[附着两相供液算例](../../../../cases/hydro_moment_supply_rz/README.md)接入真实材料和毛细反馈；自由面离域、开放域动态网格及破碎仍待完成。本例的结果见[流体待办 §28](../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md)，当前共同进度见该文档顶部。

2026-09-28：冻结 `build/moment_ports_final.exe` 在 `runs/hydro_moment/20260928_ports_rz_final` 完成原定 20 步，预检/运行退出码均为 0。进出体积率均为 0.314159265359，累计能量账最大绝对差 1.47e−11。[实际流场与能量图](../../cases/module_history/output/20260928_184817_build_archive/moment_ports_rz.png)展示发展过程；它不替代空间/时间收敛或充分发展管流精度验证。

</details>
