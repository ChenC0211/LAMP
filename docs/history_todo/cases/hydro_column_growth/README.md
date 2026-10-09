此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_column_growth/README.md。

运行入口：tools\run_case.cmd hydro_column_growth 。结果：runs/hydro_column_growth/<时间戳_标识>/。程序：LAMP。

# 轴对称液柱小扰动增长

<!-- case-definition-20260928 -->

**算例类别：** 旧路线线性阶段对照。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

无黏毛细不稳定性，初始相对扰动幅值 0.001；计算早期增长率，达到幅值 0.03 可按字典提前终止。

## 计算域、初始条件与边界条件

轴线正则，轴向反射截面。线性阶段的增长率不代表非线性断裂时刻已收敛。

```text
Domain { min (0, 0); max (4, 5.235987755982989); dimension 2; }
```

```text
AxisymmetricColumn { levelSet liquid; radius 1; relativeAmplitude 0.001; }
```

初始速度类型为 `zero`；具体初场幅值与初始压力设置见下表。

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_column_growth/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_column_growth/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/hydro_column_growth/system/hydro.dict) |
| 材料模型 | `twoPhase` | [hydro.dict](../../../../cases/hydro_column_growth/system/hydro.dict) |
| 负相密度 ρ− | `1` | [hydro.dict](../../../../cases/hydro_column_growth/system/hydro.dict) |
| 正相密度 ρ+ | `0.0001` | [hydro.dict](../../../../cases/hydro_column_growth/system/hydro.dict) |
| 负相动力黏度 μ− | `0` | [hydro.dict](../../../../cases/hydro_column_growth/system/hydro.dict) |
| 正相动力黏度 μ+ | `0` | [hydro.dict](../../../../cases/hydro_column_growth/system/hydro.dict) |
| 表面张力系数 σ | `1` | [hydro.dict](../../../../cases/hydro_column_growth/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-8` | [hydro.dict](../../../../cases/hydro_column_growth/system/hydro.dict) |
| 初始压力跳跃 | `1` | [hydro.dict](../../../../cases/hydro_column_growth/system/hydro.dict) |
| 背景网格划分 | `(32, 42)` | [mesh.dict](../../../../cases/hydro_column_growth/system/mesh.dict) |
| 启用网格自适应 | `false` | [mesh.dict](../../../../cases/hydro_column_growth/system/mesh.dict) |
| 最大单元数 | `3000` | [mesh.dict](../../../../cases/hydro_column_growth/system/mesh.dict) |
| 时间步长 Δt | `0.005` | [run.dict](../../../../cases/hydro_column_growth/system/run.dict) |
| 计算终时 T | `20` | [run.dict](../../../../cases/hydro_column_growth/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_column_growth/system/run.dict) |
| 场输出间隔 | `20` | [run.dict](../../../../cases/hydro_column_growth/system/run.dict) |
| 界面输运方法 | `conservative` | [interface.dict](../../../../cases/hydro_column_growth/system/interface.dict) |
| 重新初始化间隔 | `0` | [interface.dict](../../../../cases/hydro_column_growth/system/interface.dict) |

完整输入文件：[geometry.dict](../../../../cases/hydro_column_growth/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_column_growth/system/hydro.dict)、[interface.dict](../../../../cases/hydro_column_growth/system/interface.dict)、[mesh.dict](../../../../cases/hydro_column_growth/system/mesh.dict)、[run.dict](../../../../cases/hydro_column_growth/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_column_growth
```

每次运行创建 `runs/hydro_column_growth/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# 半波长液柱增长入口

这是 B5a 的对称域前置算例，使用真实 Interface、守恒质量—动量、毛细力和 Hydro 联立求解。完整轴向周期、断裂事件及子滴接续仍未完成。

`geometry.dict` 的 `AxisymmetricColumn` 定义

\[
R(z)=R_0\frac{1+\delta\cos[\pi(z-z_{min})/L_z]}{\sqrt{1+\delta^2/2}}.
\]

初始化直接生成同一原生 Hermite 标量 `phi=r²−R(z)²`，不创建辅助闭合轮廓。南北端是可径向滑动的反射截面；轴端及外径采用现有零法向速度/自由滑移边界。所有量为代码单位。`radius` 是体积等效半径，`relativeAmplitude` 是公式中的 δ。

本例 `R0=rho_l=sigma=1`、`rho_g=1e-4`、两相黏度为零，`kR0=0.6`、外径 4、初始速度零，初压为 `C*sigma/R0`。网格 32×42，径向 `R0/h=8`。C 按实际离散几何初始化；不以解析体积覆盖其微小插值差。

从项目根目录启动独立运行并冻结当前程序、DLL 和输入：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/run_hydro_column_growth.ps1
```

默认粗档 `dt=0.005`、终时上限 20、振幅上限 0.03。省略 `-TimeStep` 时，脚本按实际最小单元尺寸的 `h^(3/2)` 缩放该毛细步长，再向下调整为终时的整数分步；`-RadiusCells 16` 自动得到约 `dt=0.001768`。这是统一的固定毛细尺度，实际输运 CFL 仍由求解器检查；显式 `-TimeStep` 覆盖此日程，供时间减半等研究使用。实际值、策略和原毛细步长均写入 `study-input.json`。`-PrepareOnly` 仅生成输入。脚本还支持波数、外径、初始幅值和明确记录的 Oh 变体；这些选项不会自动组成完整验收矩阵。也可直接运行 `bin/LAMP.exe --case cases/hydro_column_growth`，但推荐上面的独立目录方式。

`output/series/column_mode.csv` 通过同一原生面元的 GL8 弧长积分及 `dz=n_r ds` 给出平均半径和基频余弦振幅。`single_radius_graph=1` 才表示此时仍可作单值半径分析。`interface_conservative.csv` 和 `hydro_diagnostics.csv` 分别保留体积收支及流体求解诊断。

`run.dict` 的停止配置为：

```text
Time {
    units code; dt 0.005; endTime 20;
    Stop { model columnAmplitude; relativeAmplitude 0.03; }
}
```

阈值是 `|A|/R0`，须严格大于初始幅值且小于 1。Interface 从已发布几何计算模态，Hydro App 在整个物理步联合接受后停止；不插值事件时刻，允许一个时间步的阈值越过。若 `Fields.enabled` 和 `writeFinal` 均为 true，该步强制写末帧，普通输出间隔不影响停止。`output/series/termination.json` 记录实际步号、时刻、振幅和原因 `columnAmplitude` 或 `endTime`。写盘失败仍报告失败，不能伪装成成功终止。

脚本的 `-AmplitudeLimit` 可改阈值；`-NoAmplitudeStop` 显式恢复仅按终时停止。此能力目前只接通保守 Hydro 的 `AxisymmetricColumn`，EHD 不接受这个选项。20260925_025941 的旧冻结长程没有此停止机制，仍按原终时运行，不能从新模板推断其能力。

零初速度应比较 `A(t)/A(0)≈cosh(omega*t)`；稳定波数比较余弦振荡。只拟合从初态开始连续的 `|A|/R0<0.03` 单值半径窗口：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/analyze_hydro_column_growth.ps1 -Run runs/hydro_column_growth/<运行目录>
# 读取正在运行或未完成计算的已写完整行，必须明确标成局部快照：
powershell -NoProfile -ExecutionPolicy Bypass -File tools/analyze_hydro_column_growth.ps1 -Run runs/hydro_column_growth/<运行目录> -InProgress
```

分析输出 `growth_analysis.json`，核对冻结程序、三个物理历程的共同接受时钟，并记录拟合窗口、增长率、库存/几何差和动量残差。短窗口可辨识性与单网格 5% 指标只是诊断；三档网格、时间减半、外域敏感性及完整周期规格仍需完成，`numerical_acceptance` 保持 false。

若求解器已退出而外层结果汇总中断，可从实际退出码、冻结输入和原输出重建 `result.json`，不启动求解器；此选项不能与运行参数混用：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/run_hydro_column_growth.ps1 -SummarizeRun runs/hydro_column_growth/<运行目录>
```

汇总统一以求解器的 binary64 类型严格比较 JSON/CSV 时刻，避免 Windows PowerShell 将 JSON 解析为 Decimal 后误报不同；没有扩大时钟容差。缺少求解器实际退出码的运行不能使用此选项。

</details>
