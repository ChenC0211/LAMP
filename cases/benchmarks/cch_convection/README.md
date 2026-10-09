# CCH 光滑无黏对流研究

本研究调用生产H07–H10算子，检验规定常速度下的周期向量输运；不是完整Euler/NS求解器。材料来自三层真实纯FirstJet几何，没有守恒库存、质量修复、物理黏度或界面力。一次prepared owner定义自治半离散算子，研究RK4时间不冒充生产候选clock；实际CCH时间推进及生产边界仍归H14/H15。

```powershell
.\tools\build.cmd -Preset local-benchmarks -Target LAMP_cch_convection_benchmark
.\tools\run_case.cmd benchmarks/cch_convection
```

配置只有 `benchmark.json`，driver只属于本研究目录。每次输出保存至 `runs/benchmarks__cch_convection/<时间戳>`；默认CTest不运行长研究。

H11 受控耗散研究使用同一程序的显式配置，2026-10-06 已通过规定的算子、光滑收敛和无黏研究验收，见[报告](../../../docs/task/H11_Controlled_Dissipation_Report.md)：

```powershell
.\tools\run_case.cmd benchmarks/cch_convection/controlled -Threads 1
```

该配置完成 XY/RZ、uniform/2:1、n8/16/32 的36窗口：12个光滑 controlled 收敛窗口（T=1、32n步）和24个局部欠分辨波包中心/controlled对照（T=1、128n步）。调用生产 `buildControlled` 和原共享材料 owner，无黏度或隐藏 donor。输出 `controlled_summary/trajectory/gates.csv` 及 `controlled_edges/nodes.csv`，保存实际同阶段D积分、R/材料功、边界功、触发局部性和轴反力；节点/连接快照支持独立重组。`controlled_layout/overlaps/material.csv` 保存原始拓扑、真实质量、端口和材料率；`controlled_steps.csv` 保存全部步的四阶段功与步前后能量，以复算同一 RK4 积分收支。`layout_elapsed_seconds` 包含该布局的准备及已结束窗口累计时间，不是单个窗口独占计时。固定轴按原本质条件保持零速度，约束反力的功为零。完整范围与预设门槛见 [H11任务书](../../../docs/task/H11_Controlled_Dissipation_TODO.md)；这是规定advector研究，不能当作完整Euler/NS。开发者 `controlled-probe` 仅生成四布局n8的12窗口，归新的 runs/checks，不代表完整研究通过。原默认 H10b 28窗口/84 gate及旧四FAIL保持。

H10n起验收版本为 `H10b-dissipation-v2`，见[实施前设计变更](../../../docs/task/H10n_Donor_Dissipation_Gate_Design_Change.md)。四项 donor 对照改为初始离散耗散率 `D/E(0)` 的网格阶，仍用 `[0.7,1.4]`；原累计损失阶在 `legacy_gates.csv` 保留原公式和失败。`gate_contract.csv` 明确四组映射及其余80项full gate保持；新版本不表示历史84项通过。原轨迹、时间窗、中心空间/动能及全部边界门槛保持。

`donor_power.csv` 记录初值/终值的独立共享边耗散、原rhs能量功、真实材料项与外边界功；`donor_energy_segments.csv` 与 `donor_energy_nodes.csv` 保存重组输入。segment kind=0/1/2分别为内部、计一次的显式周期连接、未配对外段。固定轴保留并检查零功。`donor-audit` 完成12个原donor窗口，按新四gate及独立身份返回0/1（异常2）；同时输出旧四gate，不把该模式通过当作中心PDE通过。

研究固定运行XY/RZ、uniform/中央真实2:1、n=8/16/32，每组完整周期T=1，中心生产载荷和显式donor对照使用同一质量/流量；全部native自由度计入误差，RZ轴径向值按本质条件固定0、轴未约束载荷另列。研究显式声明相反外段的同分量/同区间和平移关系，生产H10o闭合入口验证并处理周期等值类；保留全部原half-dual。这是算子端口契约，不是完整生产边界/时间推进App。独立时间研究使用n=8的20/40/80步与640步参考，时间自收敛不能抵消空间误差。

结果包括：

- `run.log`：SPACE、SPACE_ORDER、TIME、GATE和总裁决。
- `output/series/space.csv`：终值/瞬时截断误差、动能漂移、donor损失、边界/材料功、R/G和轴诊断。
- `output/series/trajectory.csv`：每个完整窗口的定期误差/动能时序。
- `output/series/gates.csv`：预先定义的数值门槛，失败仍保留全部完整研究结果。

`probe`为开发者n=8 XY AMR完整窗口定位，`audit`保留原n=8/16逐native瞬时affine反例矩阵，`audit32`补充最细四组，均只能在 `runs/checks` 中调用。audit的解析外边界只用于瞬时参考，不进入时间积分；输出 `affine.csv` 区分内部coarse–fine和外边界位置，并在日志记录最终算子的矩残差/共享pairs/迭代数。内部误差超过1e−10明确失败，外边界误差仍如实保存。

`rows`与`rows32`是被动最终算子审计，分别运行原n=8/16八组合和n=32四组合，只能在`runs/checks`中调用。独立稀疏矩阵包含原K、共享J、实际显式研究端口及uR，与真实apply核对≤1e−12；内部/外部/固定轴全部行的误差平方贡献须加回整体。`rows.csv`保存实际误差、真实质量/h、无量纲行系数、wrapped一次矩、二次矩/h和三次矩/h²；日志给出各组统计。诊断参考不参与系数、载荷、时间积分或门槛，退出0只表示矩阵与观察器核验通过，不代表空间一致性通过。

`compatibility`与`compatibility32`为被动材料空间矩诊断，分别n=4/8/16与n=32，四网格组合、rho=1/1000/1e6，共48组。真实圆形Pure FirstJet按dt=0.1/n推进，消费原三层几何/质量与唯一F/R。`material_compatibility.csv`记录横向admissible affine样本的材料KE必要条件Phi、几何一次目标的skew必要条件Psi、导数质量RMS下界及R/G；native求和与独立segment表达式核验≤1e−12。模式只能在`runs/checks`调用，退出0表示诊断身份有效，没有构造或更新convection系数、状态或通量，不代表PDE通过；非零条件缺口也不能直接推论二阶不可能。完整参数和限制见[H10e报告](../../../docs/task/H10e_Bounded_Operator_Consistency_Report.md)。

H10f扩展同一观察模式：`material_affine.csv`保存全部native的真实质量/h、几何一次矩及其对应点导数、同一材料率的q/gross、从原体积率恢复的actual advector和`M*V`目标。summary追加几何目标的RMS/Linf与两个PDE/skew必要条件缺口。参考横向导数为0；固定轴原始值保存，RMS使用原轴约束。它是目标审计，**不是新convection实际apply或完整轨迹**。失败的相容投影代码已删除，只在H10f归档；参见[报告](../../../docs/task/H10f_Compatible_Bounded_Moments_Report.md)。

`certificate/certificate32`在`runs/checks`构造同一F的中心KEP可行证书，分别n=8/16及n=32四组合。`center_certificate_pairs.csv`给出共享skew系数，`center_certificate_rows.csv`给出quotient采样/真实质量/原固定C/P1目标、P2区间和l1检验。六个uniform组合可行，AMR中心seed不满足原矩预算仍明确输出0；模式exit 0只表示证书审计完整，不将系数用于full/RK4或宣称最小改动/完整PDE通过。原half-dual/native DOF未删除。

`bounds/bounds32`在`runs/checks`消费同样48组真实材料，追加两条横向affine/skew条件的必要RMS和二阶缩放RMS下界，不返回修正target、pair或状态。`material_moment_bounds.csv`给出原缺口、真实质量Gram、独立2×N SVD核验，以及原pressure→dual映射和primal面F scatter两种独立来源表达式；`material_moment_sources.csv`保留全部pressure的真实映射坐标矩和primal q。固定轴仍参与条件，误差限定为0。exit 0只表示证书核验通过；有限网格下界比值不能证明任意算子不可能二阶。H10g的M*h^4最小目标修正仅在冻结数据上评估，两个原空间阶门槛FAIL，公式和逐行witness只归档，不进入本driver的数值路径。

`intervals/intervals32` 在 `runs/checks` 观察同样 48 组实际材料，不构造最终算子。`material_interval_summary.csv`、`material_interval_normals.csv` 给出 H10h 固定一次残差区间的全部二维支持法向和最终共享 pair 消去证书；`material_transport_scale.csv` 保存 46,680 行的真实质量/三层相体积、q/R、actual advector 及输入强制导数行尺度 `h*abs(q)/M`。原 affine/compatibility/bounds/sources 输出保持。exit0 只表示观察器与身份检查完整；四组 XY 2:1、rho=1e3/1e6、n8/n32 区间不可行，候选 **FAIL**。同源纯气体 dual 的大 q/R 定位旧 quarter/散度转移的相准确性不足，不能靠共享 skew 修正移除固定对角项，也不能据此断言全部二阶方法不可能。见[H10h 报告](../../../docs/task/H10h_Final_Operator_Joint_Feasibility_Report.md)。

当前验收见[H10o](../../../docs/task/H10o_Bounded_Final_Native_Convection_Report.md)：原28窗口/15116步完成，H10b-dissipation-v2的84项通过，四旧donor gate仍FAIL。下述H10d–n说明为对应阶段的观察器和历史证据，不代表当前任务入口；下一唯一H08a。observer退出0本身仍不代表PDE通过。

`segments/segments32`是H10i被动来源观察，使用原16网格组合与4类selected布局，记录真实FirstJet的native/dual段液气物理矩、GL3参考、one-ring几何functional rank/缺口及原pure-gas witness三层边界。仅在`runs/checks`调用；没有新速度/通量系数、材料更新或convection apply。64布局共157,810查询行，38,770 native/119,040 dual，56条极粗selected切向functional不可识别须如实报告。exit0只表示来源观察核验完整。

`coincident`复现n4 XY/RZ uniform/2:1实际平面FirstJet，平面offset分别0.375/0.125，位于pressure cell内部dual面。记录合法publication、真实dual质量构建状态和全部704段的查询原因；H10i原冻结时4个dual质量均失败，20个正长度零迹拒绝。观察器exit0表示已完整记录该能力缺口，**不表示平面材料可用**。没有为零迹选择相分配或扰动几何。对应`coincident_source.csv`与H10i来源报告。

`trace-contract`是被动source能力库存，仅在`runs/checks`调用。n4 XY/RZ uniform/2:1×dual内部/native面×正负phi共16场景。[H10i1b](../../../docs/task/H10i1b_Closed_Clip_Moments_Report.md)后，8个dual位置的真实 publication/phase/mass 成功，16个clip盒与432个dual质量按独立平面参考核验；H10i1b冻结时8个native位置仍在旧geometric inventory出版阶段拒绝，不能计作mass成功。`trace_contract.csv`、`trace_clip.csv`、`trace_mass.csv`保存实际阶段、两方向矩、质量与独立参考。退出0只表示观察完整，不代表最终PDE或原H10b通过。

`native_trace_port.csv`/`native_trace_ranges.csv`保留正式Core只读typed查询、normal阶数/native来源及两侧相区间；H10i1a–e冻结的实际1,408段/40零迹/80侧范围通过，与H10i1a冻结结果逐字节相同。零段的point-only endpoint incidence单独核验，普通段与legacy ranges逐字段相同。H10i1b只迁移闭合geometry clip和真实质量消费，没有新材料F或接受状态；H10i1c完成Core拓扑，完整native publication仍待H10i1d–f。

H10i1c 在原 `trace-contract` 16 场景记录完成后追加显式 Core topology 观察，保留原 identity 次序及五份旧 trace/clip/mass CSV。`trace_topology.csv` / `native_side_support.csv` 标明独立 FirstJet **seed** 来源；16/16 Core 成功，native 面8组/40 side supports 的 length、XY/RZ physical measure、Green geometry moments按独立平面参考≤1e−12实际尺度核验，dual内部8组与旧Core graph/endpoint/frame/moments逐字段相同。这里没有接受 publication、材料F或convection更新，旧native publication仍明确拒绝。16/16 CTest、64普通布局22CSV和内部affine保持；下一任务H10i1d/e/f逐项接通真实消费链。[完整结果](../../../docs/task/H10i1c_Native_Zero_Edge_Topology_Report.md)。

H10i1d 在上述同一 seed/场景追加 `half_trace_cache.csv` 和 `half_trace_cache_summary.csv`：实际16组/1,728半面/80零半面/576 value events，普通和两侧物理相测度按独立平面参考≤1e−12实际面积核验，最大误差0。生产typed cache owns原Core查询，ordinary view不复制来源，zero legacy view拒绝；不是材料causal rate。原7份trace/clip/mass/Core topology表、64普通布局22CSV、内部affine与82行旧flux数值输出保持；17/17 CTest、normal/value/重复细化与失败事务、OMP2/496半面通过。完整FirstJet publication/H07仍待H10i1e/f；没有新F/R或convection更新。[H10i1d结果](../../../docs/task/H10i1d_Readonly_HalfFace_Trace_Cache_Report.md)。

H10i1e 验收时实际 16/16 publication 成功。`publication_trace_cells.csv` / `publication_trace_summary.csv` 记录最终返回的 publication identity/step/time、相同 snapshot/cache authority、实际 352 个 CV 的派生 phase volume 与 raw assembled axial/radial ledger、独立平面参考；最大误差8.23e−16，native8组有40 side graphs/80 zero halves。原8份clip/mass/port/Core/cache表、22普通来源CSV及affine保持。`trace_contract.csv` 的 native8组如今明确在 H07 `phase` 阶段拒绝，尚未完成 RegionSource/native dual mass，不能计为来源全链或 PDE PASS。当时下一任务H10i1f；原H10b失败历史、full/RK4/donor/gate保持。[报告](../../../docs/task/H10i1e_FirstJet_Trace_Publication_Measure_Report.md)。

H10i1f 当前 `trace-contract` 的原16场景全部完成publication/H07/真实dual mass。新增 `spatial_face_trace.csv` 保存864 actual face及typed kind/selected availability/两侧物理迹，`native_region_source.csv` 保存16原Snapshot来源和独立两方向体积；32clip盒、864质量、2,816 typed段通过。六份原publication/Core/cache表字节相同；四份扩展表的原ordinary8组原始行相同。空间两侧平均只用于诊断；旧rate对重合面非零actual velocity明确拒绝，零速度率为零，尚无新causal材料owner。原22普通来源CSV、affine与82行flux保持，17/17回归通过。下一H10i2单独复核，再H10j；不能以geometry消费PASS关闭H10或重写原H10b失败历史。[H10i1f报告](../../../docs/task/H10i1f_Native_Trace_Phase_Consumers_Report.md)。

H10i2在原trace-contract中只追加source_chain_faces/duals/overlaps三份坐标表：16原实际场景、864 face/864 true dual/1472 overlap由独立坐标公式重算，两方向region/clip、原64布局157810查询及三层pure-gas witness核验通过。39原CSV、内部affine与旧数值路径完整保持。独立公式试验不进入本driver或production；新kinematics/causal phase/material F尚未实施。现H10j0先准备同一actual advector与可复用geometry，H10j再消费该owner的phase积分，不把对称空间迹当材料率。[H10i2结果](../../../docs/task/H10i2_Spatial_Material_Source_Reaudit_Report.md)。


H10j0 的只读实际速度 owner 研究使用同一程序的 `advector` 模式（`LAMP_cch_convection_benchmark.exe advector`），在新的运行目录执行。它保存64布局/128实际更新的face、atom、cell与复用/耗时表；固定几何准备只做一次。输出不进入相率、材料F或原convection算法，也不计作H10b通过。原full/probe等命令与gate保持；接口见[H10j0 API](../../../docs/api/H10j0_Cch_Spatial_Advector.md)。

H10j的独立只读 `spatial` 模式在新的 `runs/checks/H10j/<timestamp>` 运行，保存64布局/128实际factory的spatial_phase_atoms/views/summary及三层pure-gas dof103的21段spatial_gas_witness。它消费正式GALS midpoint和同一advector，不生成材料F、不替代full。原 `trace-contract` 增加spatial_zero_rate，观察原16来源的四种actual velocity（正/负/零/变号）；既有其他23份CSV字节保持。phase/API验收不代表当前材料consumer已迁移，实际接口见[H10j API](../../../docs/api/H10j_Cch_Spatial_Phase_Rate.md)。

`material` 模式在新的 `runs/checks/H10k/<timestamp>` 消费正式spatial phase/actual advector与唯一CchMaterialFlux，输出material_spatial_atoms/views/cells/summary及三层gas witness。64布局×平面/移动圆球×三密度，共384材料准备；独立后处理重算相交物理积分、三层质量、实际incidence/R。几何规定运动与Hydro MMS不同，这只验证F/诊断身份，不能说明R小余项或耦合PDE通过。新增模式不改full/RK4/donor/gate；原20份source CSV字节保持，affine4份中3份保持、affine.csv数值最大变化3.290e−13（行/字段/门槛不变）。研究helper已消费实际B，不再拟合第二份native advector；此为H10k冻结时状态；H10l2后生产dual直接共享同一F。

H10l失败候选的`dual`观察模式、576准备及真实R/E_sp/G表只保存在冻结归档，该次driver恢复H10k；H10l2后重新接入经当前源码/独立数据验证的同名观察器。初版研究runtime与最终失败candidate2测试runtime分别绑定；见[失败报告](../../../docs/task/H10l_Same_F_Dual_Material_Migration_Report.md)。当前标准material/full等入口保持，不把observer exit0视为任务或PDE通过。

H10l1新增两个只读模式，只能在新的`runs/checks/H10l1/<timestamp>`运行：`closure`输出n4 XY/RZ uniform/2:1全部208非轴basis的实际faces/cells/status，记录原DI与实际场；`boundary`输出n4/8/16/32×XY/RZ×uniform/2:1×full/left-selected共32布局的绑定/未绑定64个actual revisions，保存boundary atoms/faces/cells/segments/summary。显式law选择真实boundary DOFs，未绑定对照不自动视作墙面；独立后处理核验物理积分/DI/cellVelocity/共享。原closed流函数没有内部零均值变号见证，该项由既有material短契约的selective-law平面fixture验证（含AMR真实共同侧求和）。两模式不改full/RK4/donor或gate，退出0不代表PDE通过；见[H10l1报告](../../../docs/task/H10l1_Spatial_Advector_Closure_And_Boundary_Report.md)。

H10l2当前`dual`模式在新`runs/checks/H10l2/<timestamp>`输出dual_cells/segments/primal_cells/primal_faces/overlaps/summary与gas_witness。64布局×三族×三密度共576准备；三族分别为静止平面/MMS、移动圆球/MMS及同向平移coupled制造族。独立脚本重算三层质量、同源相积分/F、真实incidence/R/E_sp/G及uR/材料动能，强制几何族不证明R小余项。gas dof103保存三层×三密度63段。原full/RK4/donor/gates不变，退出0仍不等于完整PDE通过；见[H10l2报告](../../../docs/task/H10l2_Same_F_Dual_Material_Retry_Report.md)。

H10m冻结的`donor-audit`曾按旧四gate返回1，细节见[报告](../../../docs/task/H10m_Final_Native_Same_F_Convection_Report.md)。H10n保持原`donor_audit.csv`、`trajectory.csv`和`donor_fourier.csv`，另加上述独立耗散身份与显式版本；只能在新的`runs/checks/<任务>/<timestamp>`调用，不覆盖历史输出。新旧退出语义有明确版本区分，原数值轨迹不变。


H10o当前结果见[报告](../../../docs/task/H10o_Bounded_Final_Native_Convection_Report.md)：显式闭合算子已完成原28窗口，H10b-dissipation-v2为84/84 PASS；旧四gate仍单列FAIL。full/probe/rows使用同一生产闭合入口，原F、rhs、RK4、donor及全部门槛定义保持。audit/audit32明确使用非周期解析trace重载，新闭合入口另有允许affine/独立矩短契约，不能把非周期x样本强塞进周期alias条件。

rows/rows32的独立矩阵包含周期类按真实总质量分配的最终载荷，BOUNDED_MOMENTS另列实际迭代/精化、物理矩残差、系数界和支撑。certificate/certificate32现在导出完整允许one-ring边集（包括零seed系数），行表增加同源M*B/h点PDE目标；只作独立可行性核验，不把外部优化解送入生产。最细准备成本仍较高，普通不透明trace及任意sharp-cut场景的全native有界PDE尚不由本研究证明。下一任务H08a，不启动H11。

H08a2r6 将 `trace-contract/spatial_zero_rate.csv`、`advector/advector_atoms.csv`、`boundary/boundary_atoms.csv` 的 c/e 列统一为按实际 coefficient 容器完整输出。当前三项 header/物理数据保持；未来五项表示会包含 c0…c4（双端表另含 e0…e4），下游应按列名读取尾部物理量。高项为零也不省略，不改变上述模式或验收语义；见[格式验证报告](../../../docs/task/H08a2r6_Benchmark_Polynomial_Schema_Report.md)。

H08a2r13 当前观察器读取公共物理来源，`advector_faces.csv` 已删除旧 rank/deficit/support；`closure_faces.csv` 已删除旧 mean_x/mean_y/scale/slope。当前多项式完整输出 c0…c4/e0…e4；旧 H10j0/H10l1 的三系数/patch 参考只用于历史基线。`advector_summary.csv` 的准备报告计数仍是当前 owner 的原样诊断，不能据此筛选精度验收。

`closure` 保留原208非轴单位面通量基输入及 cells/status，新增四张表：`closure_atoms` 保存完整双端物理多项式与积分；`closure_views` 保存真实 primal/dual→atom 映射及邻接；`closure_sources` 保存身份、完整时钟和规模；`closure_samples` 每atom保存五个实际物理速度值（GL3节点再加两个独立内点）。primal 邻接和 closure_cells 的 p 是 native cell ID，dual 邻接是 velocity DOF，samples 的 pressure 是 pressure DOF；这些索引不能混用。

独立参考入口为 [verify_source.py](verify_source.py)，仅需 Python 标准库。对已完成的观察运行执行，例如：

```powershell
python cases/benchmarks/cch_convection/verify_source.py closure runs/checks/<任务>/<时间戳>/closure --output runs/checks/<任务>/<时间戳>/closure_reference.json
```

首参数可换成 `advector` 或 `boundary`，末尾运行目录须包含对应 `output/series`。成功输出 JSON，失败非零退出；只读 CSV，不修改数值结果。它以完整次数物理 GL3、双端图表、五点物理查询、原 native 输入与真实子面覆盖验收，原56/64适用解析布局由原fixture固定指定。`boundary_cells` 的旧体积坐标割线仅为历史诊断，不用作新场点散度。零均值反号及极小相片仍由原 phase/material 短契约验证，不能由 closed 流函数观察器代替。该参考不证明时空材料耦合或完整 PDE；`full` 原28窗口/84gate及标准运行命令保持。详见 [r13 设计](../../../docs/task/H08a2r13_Benchmark_Source_Observer_Design.md)。
