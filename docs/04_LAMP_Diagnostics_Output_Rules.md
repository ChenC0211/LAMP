# LAMP 输出与诊断规范

**2026-10-02 授权变更：** [纯 Hermite–GALS](task/I12_Pure_Hermite_Gals_Design_Change.md) 没有独立 M 或相库存 ledger。该模式必须报告同源几何体积、相对初值漂移及解析体积误差；守恒库存/通量残差输出“不适用”，不能填零冒充守恒。位置/梯度/曲率/AMR 误差和 accepted/rejected 诊断保持，Hydro/电荷诊断未被取消。I12p3c 已验证独立纯路线的 VTK XML/PVD、accepted/rejected、series/meta、配置追踪及 App 阶段计时，见[输出 API](api/I12p3c_FirstJet_Output_API.md)。Hydro/Electricity 尚无的物理量为 NA；[I12p3d](task/I12p3d_FirstJet_Performance_Scale_Report.md) 已验证核内细分计时和规模研究；计时包含关系与 regrid 源/目标规模见[性能 API](api/I12p3d_FirstJet_Performance_API.md)。[I12p3f](task/I12p3f_FirstJet_Curvature_Diagnostics_Report.md)已补齐全部接受态曲率 RMS/unavailable count，物理求积权重及 NA 语义见[API](api/I12p3f_FirstJet_Curvature_Diagnostics_API.md)，原最大曲率与 VTP 同源点值保持。

**版本：** 2026-09-30 Rev.1  
**目标：** 用户无需翻 debug 日志，就能判断计算进行到哪里、为什么接受/拒绝、数值守恒是否正常、性能瓶颈在哪里；空间场与界面结果以 ParaView 可直接读取的 VTK XML 系列格式输出，便于直接查看时间演化、AMR 网格、界面几何和多物理场分布。

## 1. 基本规则

**2026-10-03 材料耦合：** 按用户授权的[非守恒材料路线](task/H07_Pure_Gals_Hydro_Coupling_Design.md)，后续 Hydro 输出须区分几何材料质量漂移、逐单元/dual 材料连续性缺陷、由该缺陷产生的动量一致性项及动能项。它们不要求为零，不能填零声称守恒；尚未实现的字段为 NA。瞬时 phase rate 与时间积分流量的单位/时间层分开记录，电荷守恒诊断保持。

1. 数值 kernel 不写文件；App/IO 在结果返回后写文件。
2. 只有 accepted step 可以进入正常 frame manifest。
3. rejected candidate 必须记录原因，但不得伪装成物理时间已经推进。
4. 输出频率不能改变计算路径或验收结果。
5. performance timer 只做诊断，不参与数值 acceptance。
6. 所有关键残差必须注明归一化方式；禁止只输出一个没有定义的 `residual`。
7. 空间场与界面帧统一采用 ParaView 原生支持的 VTK XML 数据格式；CSV/JSON 继续承担标量时序、诊断和元数据职责。
8. 可视化输出是只读派生结果，不得反向参与时间推进、守恒更新或 acceptance 判据。

## 2. 目录

每个唯一运行根目录下固定：

```text
output/
    fields/
        volume/
        interface/
        rejected/
    series/
    meta/
```

目标文件：

```text
output/fields/solution.pvd              # ParaView 时间序列入口，推荐用户直接打开
output/fields/volume/volume_000000.vtu  # 体网格、AMR 网格及体场
output/fields/interface/interface_000000.vtp
                                        # Hermite 零水平集/界面及界面场
output/fields/rejected/                  # 可选：拒绝步快照，不进入 solution.pvd
output/series/conservation.csv
output/series/solver.csv
output/series/energy.csv
output/series/performance.csv
output/series/events.csv
output/meta/frames.csv
output/meta/effective_config.dict
output/meta/configuration_trace.csv
output/meta/manifest.json
output/meta/build_info.txt
output/meta/failure.json                 # 仅失败时
```

ParaView 输出格式统一如下：

- **体网格/体场：** 默认使用 VTK XML Unstructured Grid，扩展名 `.vtu`。quadtree/AMR 叶单元以真实叶单元拓扑写出，不为可视化重新规则化网格。
- **界面：** 使用 VTK XML PolyData，扩展名 `.vtp`。二维/轴对称计算中，Hermite 零水平集以线段/高阶采样后的 polyline 表示。
- **时间序列：** 使用 `.pvd` collection 文件记录 accepted frames 与物理时间；ParaView 打开 `solution.pvd` 即可浏览时间序列。
- **并行输出：** 若以后采用 MPI，每步分别允许使用 `.pvtu` / `.pvtp` 作为并行索引文件，各 rank 写对应 piece 文件；顶层仍由 `solution.pvd` 管理时间序列。
- **不推荐作为 V2 默认格式：** legacy `.vtk`。它可作为临时调试格式，但 production 输出统一使用 VTK XML。

`solution.pvd` 只登记 **accepted step**，其 `timestep` 必须等于真实物理时间。volume 与 interface 可作为同一时间点的不同 group/part 写入同一 collection。`frames.csv` 继续保留，作为人可读和脚本可读的 accepted-frame 清单，不由 `.pvd` 取代。

暂未迁移的模块可继续写现有等价文件，但新 V2 任务不得再创造新的平行目录体系。

## 3. 控制台：只显示运行主状态

每个 accepted step 默认一条 summary：

```text
[Step 00120] t=1.200000e-3 dt=1.0e-5 CFL=0.31 NL=2 LIN=14 \
  div=2.1e-12 dV/V=3.8e-14 dM/M=4.5e-14 umax=1.72 \
  kappaMax=8.13e1 wall=0.84s
```

Electricity 接入后追加简短项：

```text
  Esym=2.0e-14 Gauss=3.1e-12 dQ/Q=5.0e-13
```

拒绝步：

```text
[REJECT] step=120 t=1.19e-3 reason=phase_flux_boundedness \
  cell=3812 residual=2.7e-8 dt=1.0e-5->5.0e-6
```

详细 assembly、patch、quadrature、Krylov history 只在 `Log.level=detail/debug` 时输出，并优先写日志而不是刷屏。

## 4. ParaView 场输出约定

### 4.1 Volume `.vtu`

每个输出帧至少应包含网格拓扑和能够判断主要物理状态的基础场。字段按实际模块可用性逐步增加，但命名一旦进入 production 不应随意改变。

建议基础字段：

```text
CellData:
    level                   # AMR level
    cell_id
    phase_volume_fraction   # alpha
    density
    viscosity
    pressure
    divergence

PointData 或派生 CellData:
    velocity
```

Hydro 的 authoritative velocity 仍然是原生 staggered face DOF。为了 ParaView 显示而生成的 point/cell `velocity` 仅是**可视化重构量**，不得反向替代数值求解中的 face velocity，也不得用于守恒或验收计算。

Electricity 接入后可增加：

```text
potential
electric_field
electric_displacement
```

如果 cut cell 内两相存在不同的电势 polynomial，不允许为了输出而把二者无说明地平均成一个“权威电势”。可以输出清楚命名的 diagnostic view，或在后续需要时使用 phase-resolved block。

### 4.2 Interface `.vtp`

界面输出至少保存几何位置，并按模块进度增加：

```text
PointData / CellData:
    normal
    curvature
    phase_flux
    surface_charge
    normal_D_minus
    normal_D_plus
    tangential_E
    maxwell_traction
    J_Gamma
```

其中 `curvature`、`normal` 必须来自当前 accepted `HermiteGeometrySnapshot`；不得为了画图另外拟合一条独立界面后再计算曲率。

### 4.3 帧命名与时间一致性

统一采用零填充帧号，例如：

```text
volume_000120.vtu
interface_000120.vtp
```

帧号对应 accepted-frame index；物理时间由 `.pvd` 的 `timestep` 和 `frames.csv` 明确记录，不允许从文件序号推断物理时间。

生产输出推荐使用 VTK XML appended binary；若实现阶段暂时使用 ASCII XML 便于调试，必须保持完全相同的数据语义与字段命名。压缩属于 IO 选项，不得改变数值结果。

### 4.4 rejected step

若启用失败/拒绝步场快照：

```text
output/fields/rejected/volume_rejected_step000120_attempt02.vtu
output/fields/rejected/interface_rejected_step000120_attempt02.vtp
```

这些文件仅用于诊断：

- 不进入 `solution.pvd`；
- 不增加 accepted-frame index；
- 必须能从文件名或 sidecar metadata 看出 step、attempt 和 reject reason。

## 5. `conservation.csv`

每个 accepted step 至少包含：

```text
step,time,dt,
phase_volume,phase_volume_boundary_flux,phase_volume_closure,
mass,mass_boundary_flux,mass_closure,
divergence_norm,momentum_balance_norm
```

EHD 后增加：

```text
surface_charge,boundary_current,charge_closure,global_gauss_residual
```

原则：**原始量 + closure 同时保存。** 不只写“修正后误差”。

## 6. `solver.csv`

至少记录：

```text
step,time,module,stage,status,
nonlinear_iterations,linear_iterations,
linear_residual,physical_residual,
rejected_reason
```

不同系统分开记录，例如 `hydro_saddle`、`electric_skeleton`、`charge_schur`。禁止把三个系统的迭代数加成一个无法解释的数字。

## 7. `energy.csv`

Hydro 开始 production benchmark 后至少记录：

```text
step,time,kinetic_energy,viscous_dissipation,
surface_energy,external_power,numerical_dissipation_estimate
```

在尚未建立严格全离散能量恒等式的阶段，字段必须写成“diagnostic/estimate”，不得把诊断量包装成已证明的守恒律。

Electricity 增加：

```text
electrostatic_skeleton_energy,reconstructed_field_energy,stabilization_energy_difference
```

## 8. `performance.csv`

按 accepted step 记录 wall time，并使用稳定 stage 名：

```text
interface_trace
interface_flux
interface_closure
hydro_assembly
hydro_linear_solve
hydro_corrector
electric_geometry
electric_assembly
electric_solve
charge_schur
io
step_total
```

同时记录结构规模：

```text
leaf_cells,native_faces,cut_cells,interface_panels,
hydro_dofs,electric_skeleton_dofs,nnz,agglomerated_elements
```

这样可以区分“算法变慢”与“表示规模失控”。

## 9. 模块硬诊断

### Interface

每步或规定频率记录：

- phase-volume closure；
- geometric volume 与 conserved inventory mismatch；
- max/mean bubble correction；
- cut-cell/interface element count；
- coarse–fine reflux residual；
- curvature max/RMS 与 unavailable sample count。

### Hydro

- local/global continuity；
- `C=-D^T` adjoint defect（startup/mesh-change 时必查）；
- dual-mass closure；
- kinetic energy 与 controlled dissipation；
- `u^T K_mu u`；
- linear/nonlinear residual；
- AMR/Hodge conditioning summary。

### Electricity

- symmetry defect；
- min local factor pivot / max `cond(Kcc)`；
- agglomeration count；
- local/internal/interface/global Gauss residual；
- panel current closure；
- charge-Schur residual；
- skeleton vs reconstructed field energy。

### Coupling

- geometry identity/time layer；
- `J_Gamma` range；
- `Cp + F_Gamma` static-balance diagnostic；
- midpoint corrector change in geometry / mass flux / traction。

## 10. `manifest.json`

至少包含：

```text
run_id
case_path
start_time
program_version / git_commit (if available)
build_type
compiler
selected_methods
mesh_identity
initial_geometry_identity
input_file_hashes
accepted_steps
final_time
exit_status
```

没有 Git 仓库信息时写 `unavailable`，不得伪造 commit。

## 11. failure 语义

发生数值失败时：

- 退出码非零；
- accepted state 不前移；
- `failure.json` 记录 module/stage/step/time/reason/diagnostics；
- 如果启用失败场快照，统一写入 `fields/rejected/` 下的 `.vtu/.vtp`，且**不进入 `solution.pvd` 与 `frames.csv` accepted list**；
- 不能因为输出失败而把数值失败改成成功，也不能因为数值已接受而吞掉 IO failure。

## 12. 生产输出的完成标准

一个新模块只有在：

- 控制台 summary 可读；
- `solution.pvd` 可由 ParaView 直接打开，并能正确浏览 accepted time series；
- volume `.vtu` 与 interface `.vtp` 的基础字段、时间和 geometry identity 一致；
- `conservation.csv`、`solver.csv` 与 `performance.csv` 有真实数据；
- effective config 和 build/meta 可追溯；
- rejected step 可定位原因；

之后，才允许进入长时间电喷雾运行。
