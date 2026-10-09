此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_conservative_rz/README.md。

运行入口：tools\run_case.cmd hydro_conservative_rz 。结果：runs/hydro_conservative_rz/<时间戳_标识>/。程序：LAMP。

# 高密度比 RZ 开放域液滴

<!-- case-definition-20260928 -->

**算例类别：** 旧保守输运功能示例。以下设置根据当前 system 字典整理；本次目录整理未运行物理时窗。

## 计算目的与物理模型

不可压缩变物性两相流及表面张力；密度比 1000，研究气相供排下的液滴输运和守恒。

## 计算域、初始条件与边界条件

南侧法向速度入口、北侧压力出口；供排口外部相均为气相。负的外法向速度表示流入。

```text
Domain { min (0, 0); max (1, 1); dimension 2; }
```

初始速度类型为 `zero`；具体初场幅值与初始压力设置见下表。

显式边界输入：

```text
boundaryPatches
    {
        supply { side south; type normalVelocity; start 0; end 1; value -0.02; }
        outlet { side north; type pressureOutlet; start 0; end 1; value 0; }
    }
```

显式边界输入：

```text
portPhases { supply gas; outlet gas; }
```

## 工况、网格与数值设置

数值按字典原单位列出。`Time.units code` 表示代码时间；存在 `units.dict` 时，物理时间按其 T0 换算。未显式声明的参数不在本文猜测为默认值。

| 项目 | 当前字典值 | 来源 |
| --- | --- | --- |
| 运行模块 | `SingleModule` | [run.dict](../../../../cases/hydro_conservative_rz/system/run.dict) |
| 子模块 | `(Hydrodynamics)` | [run.dict](../../../../cases/hydro_conservative_rz/system/run.dict) |
| 坐标系 | `axisymmetricRZ` | [hydro.dict](../../../../cases/hydro_conservative_rz/system/hydro.dict) |
| 材料模型 | `twoPhase` | [hydro.dict](../../../../cases/hydro_conservative_rz/system/hydro.dict) |
| 负相密度 ρ− | `1000` | [hydro.dict](../../../../cases/hydro_conservative_rz/system/hydro.dict) |
| 正相密度 ρ+ | `1` | [hydro.dict](../../../../cases/hydro_conservative_rz/system/hydro.dict) |
| 负相动力黏度 μ− | `0.01` | [hydro.dict](../../../../cases/hydro_conservative_rz/system/hydro.dict) |
| 正相动力黏度 μ+ | `0.001` | [hydro.dict](../../../../cases/hydro_conservative_rz/system/hydro.dict) |
| 表面张力系数 σ | `0.01` | [hydro.dict](../../../../cases/hydro_conservative_rz/system/hydro.dict) |
| 线性方程相对残差容限 | `1e-8` | [hydro.dict](../../../../cases/hydro_conservative_rz/system/hydro.dict) |
| 背景网格划分 | `(16, 16)` | [mesh.dict](../../../../cases/hydro_conservative_rz/system/mesh.dict) |
| 启用网格自适应 | `false` | [mesh.dict](../../../../cases/hydro_conservative_rz/system/mesh.dict) |
| 最大单元数 | `1024` | [mesh.dict](../../../../cases/hydro_conservative_rz/system/mesh.dict) |
| 时间步长 Δt | `0.001` | [run.dict](../../../../cases/hydro_conservative_rz/system/run.dict) |
| 计算终时 T | `0.02` | [run.dict](../../../../cases/hydro_conservative_rz/system/run.dict) |
| 时间单位 | `code` | [run.dict](../../../../cases/hydro_conservative_rz/system/run.dict) |
| 场输出间隔 | `10` | [run.dict](../../../../cases/hydro_conservative_rz/system/run.dict) |
| 界面输运方法 | `conservative` | [interface.dict](../../../../cases/hydro_conservative_rz/system/interface.dict) |
| 重新初始化间隔 | `0` | [interface.dict](../../../../cases/hydro_conservative_rz/system/interface.dict) |

完整输入文件：[geometry.dict](../../../../cases/hydro_conservative_rz/system/geometry.dict)、[hydro.dict](../../../../cases/hydro_conservative_rz/system/hydro.dict)、[interface.dict](../../../../cases/hydro_conservative_rz/system/interface.dict)、[mesh.dict](../../../../cases/hydro_conservative_rz/system/mesh.dict)、[run.dict](../../../../cases/hydro_conservative_rz/system/run.dict)。

## 运行方法与结果目录

在项目根目录运行：

```powershell
.\tools\run_case.cmd hydro_conservative_rz
```

每次运行创建 `runs/hydro_conservative_rz/<年月日_时分秒_毫秒_标识>/`，保存 `system/` 输入副本、`runtime/` 程序副本、`run.json` 和 `run.log`。求解器原生结果位于该运行目录的 `output/` 下，场文件在 `fields/`，时间历程在 `series/`，元数据在 `meta/`；部分独立基准另在原生输出根目录写 `metrics.csv/summary.json`。保留这层运行快照结构，使现有求解器和后处理接口一致。

仅准备输入副本使用 `-PrepareOnly`；支持输入预检的应用可使用 `-CheckOnly`。修改工况应创建新的运行副本，不能覆盖先前结果。线程设置不等于该模块每个环节都已并行。

## 输出量与结果判据

关注速度、压力、离散散度和动能历程；两相算例另记录界面轮廓、液体体积、曲率、表面能及边界质量通量。相容矩方法的速度属于步末、压力属于中点，比较参考解时须匹配时层。断颈算例记录颈半径、父子分量、事件时间和事件模型引起的能量变化。

程序退出码、输入有效性、全时窗完成和数值验收是不同结论。本 README 不因已有模板或历史通过记录而宣称当前版本通过；集中验收遵循[主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md#concentrated-benchmarks)。网格、时间误差和参考范围须同时交代，失败输出仍保留。

<details>
<summary>原方法说明与历史运行记录（旧命令及路径仅用于追溯）</summary>

# 保守 RZ 球滴与实际 NS 耦合

完整根网格上的球滴处于轴向供流中，密度比 1000，包含非零表面张力。默认 16×16、20 步到 `t=0.02`。这是可用的主程序接入示例，不是射流破碎、长期寄生流或网格收敛验收。

在项目根目录编译和运行：

```powershell
cmake --build --preset local --target QUANTUM --parallel 1 -- /nr:false
.\bin\LAMP.exe --case .\cases\hydro_conservative_rz --check-case
.\bin\LAMP.exe --case .\cases\hydro_conservative_rz
```

`cmake` 未在 PATH 时使用本机 Visual Studio 的 CMake。此 Hydro 入口不使用 Electrospray 的 `--start` 参数；保留多次结果时，将 `system` 复制到新的运行目录，再把该目录传给 `--case`。

- `interface.dict`：`transport conservative`、`advectInterface true`、`reinitializeEvery 0`；`portPhases` 显式指定每个入口/出口的全步纯相，`liquid` 表示负相、`gas` 表示正相。初始及运动后边界层必须保持该相。
- `hydro.dict`：`interfaceForceGeometry conservative` 与 `interfaceForceLayer inventory` 配套。密度/黏度来自 C，毛细力使用共享界面曲率与 `sigma*kappa*G(C)`。`interfaceThickness` 在此模式中是曲率投影搜索半宽，不是材料平滑宽度。
- `run.dict`：唯一时间步和输出计划。每步冻结当前求得的完整 NS 速度推进界面，再从推进后的 C/界面组装毛细力，用同一限幅相通量与旧/新对偶质量解 NS；黏度取步起点，毛细力取步终点候选，当前为一阶分区推进。不需要新增字典设置。

`output/fields/hydro_step*.dat` 包含同一接受时刻的 `phi、ur、uz、C、rho、mu、p、divU、C_Gamma`。`C_Gamma` 是共享界面的几何体积分数，不覆盖守恒 C。`output/series/interface_conservative.csv` 记录库存、几何差异、质量闭合、材料时间和载荷时间；`hydro_diagnostics.csv` 保留求解器、毛细力与边界收支诊断。此模式的 `negative_volume` 为 C 库存，`interface_max_change` 是相邻步标量最大变化，不是界面位移。

写接受步场时，同时保存 `hydro_capillary_source_stepNNNN.csv`：这是该步实际使用的毛细载荷，列出原生面、投影界面点、C 梯度、总/子午/环向曲率、载荷与模板标记。当前保守路径先推进界面再组装载荷，因此 `force_time=accepted_time`；此前保存的旧时间离散结果仍保留早一步的源时间，应按文件实际时间读取。初始化尚无载荷，不生成该文件。配置 `writeFirst true` 可取得第一步实际施力来源。

轴端模板额外输出 `implicit_axis`：值 1 表示三列环带体积的偶隐式图面；`even_axis=1` 且 `implicit_axis=0` 表示保留的偶多项式模板。两者均只估计同一发布界面的曲率，不改写守恒 C。

普通高度曲率现使用五列精确四次体积矩及原轴端/同弧后备。旧七列过定曲率和试验性的隐式弧段方案因平界面高频恢复符号错误已移除；无需修改字典。离散算子证据与完整耦合稳定性分别记录在[唯一待办](../../docs/QUANTUM_Modular_Consolidation_Plan.md#levelset-implementation)，不能把单算子检查或短时运行完成视为破碎或长期稳定验收。

当前接入范围是完整 XY/RZ 根网格、显式零切应力边界或无 patches 的零法向闭箱、闭合液滴（含轴线端点）和纯相端口。XY 使用 `coordinateSystem planar2D` 及相应平面几何；闭箱不配置 patches 或 portPhases，初始速度自动作零时间相容投影。非轴外壁需与界面保持纯相层隔离；分段入口法向速度在公共顶点不连续时拒绝。AMR、壁面接触事件、拓扑破碎及电场/电荷接续尚未接入此路径。

闭箱保守算例可在 `hydro.dict` 指定 `initialPressureJump`（代码压力单位，负相减正相，可带符号，默认零）。初始化使用 `p_i=initialPressureJump*C_i`，与实际守恒体积分数一致，初始速度相容投影保留该物理压力。它只规定初始条件，不修改后续曲率/毛细力；不适用于本目录原有的开放端口配置，带 `boundaryPatches` 或非保守输运时显式拒绝此键。

`tools/run_conservative_static_drop.ps1` 可在独立目录生成并运行 XY/RZ 闭箱静滴；`-InitialPressure Laplace` 指定圆的 `sigma/R` 或球的 `2*sigma/R`，默认 `Zero` 保持历史零压力启动及旧程序兼容。结果记录初始压力、`R/h` 与实际毛细时间；当前固定半径 0.23，32² 的 `R/h=7.36`，并非正式 B4 的 `R/h=8/16/32`。已保存结果仍有寄生流，不能作为长期稳定验收，见[最近的图像和说明](../../output/figures/conservative_static_drop_20260923_224343/README.md)。

</details>
