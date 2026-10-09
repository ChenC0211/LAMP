> **状态：historical/reference（历史资料，仅供追溯）。**
> 本文原始正文保留为历史记录，其中“当前”“下一步”、TODO、旧规则与通过记录均指写作当时；不再自动产生新任务，也不能证明 V2 gate 通过。
> 现行任务唯一入口：[分模块任务清单](../roadmap/05_LAMP_Module_Task_Roadmap.md)。文档权威见[开发总纲](../00_LAMP_Development_Charter.md)，完整入口见[文档索引](../README.md)。
> 旧字典/职责说明可供查询既有实现；V2 新增或修改仍遵守现行六份规范。历史正文中的链接保留原样，目标缺失时不能作为本机已验证证据。
> 本目录中的旧设计、TeX/PDF 和索引同属历史参考；总纲指定的三份冻结方法 PDF 仍以 `docs/` 下的权威文件为准，不因本目录的历史标记而降级。

# QUANTUM（电）流体模块恢复、通用化与电喷雾解耦实施计划

**最新入口（2026-09-12）：**日常跟进看[简明进度页](QUANTUM_Hydro_EHD_Progress.md)，运行与 Tecplot 看图看[算例指南](QUANTUM_Hydro_EHD_Case_Guide.md)。输出目录与标准入口见第72节，原样 20 步体积修正见第79节，导电液面接线见第80节；下文前部保留历史审计，不代表当前状态。

> 状态：实施基线（2026-09-08）；首批见第 12 节；材料输运见第 13 节；表面张力见第 14 节；修复前静滴验证见第 15 节；圆弧距离/压力恢复修复与原样复测见第 16 节；通用 EHD App 首批接入见第 17 节（2026-09-09）  
> 范围：水平集几何、四叉树 AMR、不可压流体、界面静电/表面电荷、电流体耦合与电喷雾应用  
> 依据：当前 `QUANTUM` 源码、现有 hydrodynamics/electrospray cases、编译与最小运行探针，以及旧工程 `D:\Chen\IGMPP-ILET` 中原始实现  
> 与既有文档的关系：本计划继承 `QUANTUM_Modular_Multiphysics_RefactorPlan.tex` 的分层方向，但以当前源码和实跑结果补齐实施顺序、接口契约和验收门；如二者冲突，以本计划的“先恢复、后抽象”顺序为准。

## 1. 结论先行

本节至第 9 节保留 **2026-09-08 改造前审计基线**，其中“当前”、代码行号和探针结果不代表修改后的状态。第 10 节定义首批范围，第 12 节记录首批交付，第 13 节记录材料输运接线，第 14 节记录表面张力与界面力接线，第 15 节记录修复前静滴物理验证，第 16 节记录圆弧距离/压力恢复修复和同参复测；未列为交付的后续阶段仍未实施。

当前问题不是“流体算法已经丢失”，而是一次未完成的框架迁移把入口、配置、模板实例化、共享状态和算例模型割裂了：

- 流体 XY/RZ 的主体离散仍在，关闭 OpenMP 的诊断构建中均完成了一个最小时间步；最大速度散度分别约为 `1.47e-17` 和 `3.48e-05`。这只证明求解主体尚可执行，不代表边界、自由界面或电耦合正确。
- 正常 Release 构建虽然成功，但工程内没有生产路径调用 `Hydrodynamics::solve/solveRZ`，模板函数没有被编译，属于“假绿”。显式实例化并使用生产 `/openmp:llvm` 后，XY/RZ 均因对成员变量做 OpenMP reduction 而编译失败。
- hydrodynamics 和 electrospray 原案例都会在进入模块前被 `RunController` 的 `system/` 白名单拒绝。即使临时只保留白名单文件，Hydrodynamics 又必须读取被白名单禁止的 `hydro.dict`，形成不可满足的运行契约。
- 当前 `ElectrosprayApp` 只初始化几何、材料和电场并写一次输出；`finalTime` 被显式弃用，Hydrodynamics/Electricity 的实际 `solve` 均未进入时间循环。
- 原始工程的耦合时间循环、Hydrodynamics 与 Electricity 主体实现仍可作为行为参考，但其中含有大量 `levelSets[1]`、固定电位、EMI/BF4、场蒸发、轴对称假设和硬编码几何，不能原样回填为“通用模块”。

因此采用以下总策略：

1. **恢复优先**：先让现有求解函数在正式编译配置中真实实例化、可调用、可失败、可测试。
2. **行为钉住**：用小而明确的单物理回归测试保护旧算法，再修改架构；旧缓存只能做参考，不能代替验证。
3. **共享状态单一所有者**：网格、几何、拓扑派生量、时间和字段迁移由 Runtime/App 创建一次，模块不得各建一套。
4. **物理与应用分离**：Hydrodynamics 和 Electricity 只接受材料、边界、界面条件与源项；电喷雾的离子种类、蒸发律、发射与特定时间步准则留在 Electrospray App/Model。
5. **分阶段替换**：第一阶段保留旧 staggered Hydro 与 IIM Electricity 内核并加适配层；通过基准后，才考虑合并 XY/RZ 或改为新的 FV/EB 离散。

## 2. 已验证的当前状态

### 2.1 构建和运行探针

| 检查 | 实测结果 | 含义 |
|---|---|---|
| `cmake --build build --config Release --parallel 8` | 成功 | 只说明当前已引用的代码可构建，不覆盖流体推进模板 |
| 现有 CTest | 共 6 项，无流体/电流体测试；抽跑 2 项 Poisson/几何测试通过 | 公共 Poisson/几何已有部分保护，Fluid/EHD 为测试盲区 |
| 原 hydrodynamics case | 返回码 2，白名单拒绝 `app.dict/example.jpg/gas.dict/hydro.dict/normalization.dict/solvers.dict/units.dict` | 尚未进入 Hydrodynamics |
| 原 electrospray case | 返回码 2，白名单拒绝 `electricity.dict/electrospray.dict/example.jpg/normalization.dict/solvers.dict/units.dict` | 尚未进入 ElectrosprayApp |
| 临时白名单化 hydrodynamics case | 随后因找不到 `./system/hydro.dict` 失败 | 白名单与模块必需配置相互矛盾 |
| 临时补新格式 `Domain` 的 electrospray case | 返回码 0，只初始化并输出 | 不是耦合求解成功 |
| 显式实例化 Hydro，生产 OpenMP 配置 | XY 在 `Hydrodynamics_impl.h:106`、RZ 在 `HydrodynamicsRZ_impl.h:80` 报 C3028 | 成员 `maxDivURes_`/`uMagMax_` 不能直接做 MSVC OpenMP reduction |
| 显式实例化 Electricity::solve | 编译成功，有 `size_t -> int` 警告 | 只能说明模板可生成，未验证物理结果 |
| 关闭 OpenMP 的最小一步 Hydro | XY/RZ 均返回 0 | 主体可作为恢复起点；测试中 `noSlipNodes=0`，未覆盖真实边界、移动界面、电力或守恒 |
| 生产 OpenMP 的一步 Electricity | 返回 0；表面电荷残差约 254 次降到 `9.73e-7`，法向导数残差仍约 `7.67e-2` | 证明动态入口可执行，同时暴露需要定义/验证的界面残差 |
| 无 OpenMP 的 Electricity→Hydro 单向一步 | 返回 0，Hydro 最大散度约 `1.22e-17` | 仅为零初速、无粒子、无 no-slip 的单向 smoke，不证明耦合正确 |

Hydro smoke 使用 `[0,1] x [-1,1]` 的 20x40 uniform 网格：LS0 在域内无切割，半圆 LS1 产生 29 个界面节点；`applyDirectForce` 虽被调用，但因为 no-slip 集为空，其残差直接为 0。所有临时探针和临时 case 均位于 `runs/history/hydro/template_interface_probe_20260908/`，日志在其 `logs/` 子目录，没有修改正式源码或发布 cases。

### 2.2 运行时和字典契约断裂

1. `src/Apps/Runtime/RunController.cpp:52-94` 的白名单只允许：
   `run.dict`、`pic.dict`、`geometry.dict`、`mcc.dict`、`mesh.dict`。
2. `HydrodynamicsModule.cpp:24-29` 又硬性要求 `./system/hydro.dict`；`ElectricityModule` 同样读取 `electricity.dict`。这与白名单直接冲突，也与 `QUANTUM_Developer_Guide_ArchitectureRules_v2.tex` 中允许 hydro/electricity dict 的规则冲突。
3. RunController 仍会回退读取 `solvers.dict`、`units.dict`、`normalization.dict`、`gas.dict`，但这些文件也会先被白名单拒绝；同时它把多个共享配置优先从 `pic.dict` 提取到 cache，导致 Runtime 的公共配置入口被 PIC 特化。
4. 旧 fluid cases 的 `geometry.dict` 没有当前解析器要求的 `Domain { min ... max ... }`；域范围仍写在旧 `mesh.dict/BackgroundMesh` 中。
5. 旧 `Boundary` 使用 `domainBoundary {4,...}` 和 `line 0 { ... }` 形式；当前 `DictionaryGeometry` 的实体解析只接收 `boundary <id> { ... }`。所以旧边界不会进入 `boundaries/boundaryById`，后续无法可靠绑定壁面、出口、电极和自由界面。旧 `Phase` 的 `levelSet + operation` 写法也缺少当前 v3 的 `constraints`，实跑中三个 phase 都会被跳过。
6. `cases/electrospray/system/run.dict` 已含 `Electrospray { ... }`，而同目录 `electrospray.dict` 重复同一组驱动参数；后者当前未被 App 使用。`app.dict` 也已被新规则废弃。

### 2.3 App/Module 还没有共享仿真上下文

- `IApp` 与 `IModule` 都只有 `run()`；注释中提到未来的 `SimContext`，但尚未实现。
- `SingleModuleApp` 只是依次构造并执行模块。HydrodynamicsModule 和 ElectricityModule 会分别重读 mesh/geometry 并重建对象，因此 `(Hydrodynamics Electricity)` 并不表示两个模块耦合或共享同一网格。
- ElectricityModule 为获得材料属性临时构造 Hydrodynamics，并提供了一个 `initOnExisting` 式的局部绕行。这是对象所有权未定的症状，不应继续扩展。
- `ElectrosprayApp.cpp:189-192` 明确停用了旧时间循环；`223-231` 又用圆和水平线覆盖字典水平集，`250` 一带固定选择 level set 1。当前 App 是初始化 demo，不是电喷雾求解器。

### 2.4 Hydrodynamics 的接口和实现问题

Hydrodynamics 的主要数值流程仍在：外力扩散、非一致网格辅助系数、动量方程、中间速度、压力投影、最终速度、界面退移速度和直接力边界。但至少存在以下恢复阻塞：

1. **模板未覆盖**：生产目标没有调用 `solve/solveRZ`，因此上述 OpenMP 编译错误长期隐藏。
2. **变量步长接口回归**：旧接口为 `solve(mesh, dt, preDt, ...)`，对流项系数使用 `dt/preDt`；当前接口移除了 `preDt` 并固定成 `1.5/0.5`。而旧 driver 每步又会重新计算 `dt`，两者数学上不一致。恢复时必须选择“严格固定步长”或恢复变量步长 AB2/BDF 系数，不能静默混用。
3. **索引类型错误**：XY 约 `1314-1325`、RZ 约 `1113-1123` 等处把 `supSonCellIndex` 声明成 `double` 后用于数组索引，关闭 OpenMP 时会产生大量 C4244。
4. **材料硬编码**：`updateProperties` 直接读取 `Core/Base/Base.h` 中的液/气密度、黏度和 Reynolds 数。旧版还把两个特定 level set 的结果相加再截断；新版只接收一个 level set。这两种方式都不是通用的 phase-to-material 映射。
5. **魔法序号**：`applyEvapVel`、`applyDirectForce` 等路径固定使用 level set 1。
6. **无滑移链路失效**：`SharpInterface::initialize` 中 `setIgnore/createNoSlipNodes/buildNoSlipStencil` 均被注释，`noSlipNodeNum` 保持 0；Hydro 的直接力路径因而通常是空操作。恢复后其迭代上限又被注释，存在不收敛时死循环风险。
7. **边界条件特化**：压力方程只在固定方位钉 `p=0`，其余边界分支有大量注释；不能表达按 boundary ID/name 配置的速度/压力条件。
8. **失败语义错误**：SparseLU 失败多处调用 `exit(0)`，会把求解失败报告为成功并绕过上层清理/诊断。
9. **多界面外力覆盖**：`calculateExternalForce` 在每个 moving level set 上重新清零总外力，多移动界面时前一个界面的贡献会被后一个覆盖；应在 step 开始清零一次、随后累加 provider 的贡献。
10. **RZ 行为漂移**：旧实现每步做 4 次压力修正，当前只做 1 次；压力参考/边界位置也发生改变。这是数值算法变化，不是单纯 API 改名，必须用回归与解析基准决定最终行为。
11. **维护风险**：XY 实现约 6600 行，后半有大段注释掉的重复实现；RZ 又复制约 3700 行。第一轮只做必要修复，不在没有回归测试时同时重写离散。

### 2.5 Electricity 混合了通用电学与电喷雾模型

1. 当前不只是 `Electricity::initialize` 输出了制造解误差；默认矩阵/RHS 本身已被 `k=3`、两侧 `beta=1` 和 `cos(r^2)cos(z)`/`cos(r^2)sin(z)` 的 manufactured Helmholtz 问题侵入。旧物理实现中的两侧介电系数、空间电荷、表面电荷 jump 和电极边界项已被替换或删除。MMS 必须整体迁移到 verification source/BC policy，再恢复物理 operator 为默认。
2. `Electricity::solve` 同时负责：
   - 两介质电势与电场；
   - 界面表面电荷迭代和电流；
   - Maxwell 法/切向应力；
   - EMI/BF4 场蒸发；
   - 粒子生成、推进、收集统计。
3. 主路径反复硬编码 `levelSets[1]`、`nodeIndexToLabel[1]` 和轴对称 `1/r` 项，无法复用于任意命名界面或平面几何。
4. `Core/Base/Base.h` 中存在真空/液体介电率、298 K、EMI/BF4 质量、液体密度/黏度/表面张力、1 μm 特征尺度、±8 无量纲电位、导电率和场蒸发参数。这些是 case/model 配置，不是 Core 常量。
5. `emitParticle` 当前又把计算得到的发射数置为 0，实际发射被强制关闭。旧 species fraction、蒸发电流时间平均、collector/extractor 累计和若干粒子边界行为也未等价迁移，不能把“函数被调用”当作粒子发射已实现。
6. 旧 Poisson RHS 使用的 `interfaceSpaceRhoe` 及相关诊断状态在当前头文件中消失，可能切断体空间电荷到界面 jump 的反馈；恢复物理 operator 时必须逐项建立状态映射表。
7. 现有 `PoissonFVSolverCellCentered`、`PolarizationPoissonFVSolverCellCentered`、`FaceSegmentGraph` 和 `CutCellMetrics` 是重要的公共资产；但当前 PIC Poisson 主要处理单活动区域和固体 EB。液体—真空内部界面需要两侧系数与电势/通量跳跃，未经制造解验证前不能直接替换旧 IIM。

补充运行探针中，`Electricity::solve` 在生产 OpenMP 下可以运行：29 个圆界面节点、零速度、空粒子种类时，表面电荷迭代在 254 次后降到约 `9.73e-7`，但法向导数残差停留在约 `7.67e-2`。继续单向执行 Electricity→Hydro 也能返回 0，Hydro 最大散度约 `1.22e-17`；这仍是无 no-slip、无真实 phase/BC、无时间循环的 smoke，反而说明电学 jump 残差必须进入正式验收。

### 2.6 Mesh/Geometry 中仍有物理反向依赖

- `QuadtreeMesh::adapt` 的公共签名硬编码 `magGradU`，无显式 target level 时还依据速度梯度阈值细化；Core/Mesh 因而知道“流体速度梯度”。目标接口应只接收 refinement marks/target levels，物理指标由 `QuadtreeMeshRefiner` 的 provider 提供。
- AMR 字段迁移依赖调用者手工传入 vertex/h-edge/v-edge/cell 四组指针。Electrospray 只登记少数字段，新增字段很容易在 adapt 后漏迁移或被重置。
- `MovingEBController` 虽位于 Core/Geometry，却依赖 Particle、Normalization 和 `Solvers/Poisson/CutCellMetrics`，内部还包含离子撞壁、产额/侵蚀和 CSV 诊断。只能复用其中纯几何运动部分；粒子撞击和侵蚀律应回到 Plasma/Coupling 或具体 App Model。
- `SharpInterface` 可以按任意 level-set 数量保存界面场，这是可保留的能力；但“移动界面节点”和“固定无滑移节点”目前用不同且不完整的全局接口，应改为按语义 boundary/interface handle 构造的视图。
- SharpInterface 的曲率/表面张力路径含固定的轴对称 `-n_r/r` 修正和全局无量纲系数；平面与 RZ 必须通过 GeometryMetric/CoordinateSystem 显式选择，不能让公共界面类永远假定 RZ。
- `InterfaceNode` 和 `Vert` 还直接保存 IIM 系数、蒸发速度等 solver/model sidecar 状态；长期应由 solver-owned stencil/cache 按 label/revision 保存，几何实体只保留几何与拓扑信息。

ExplicitPIC 已实际走通 `QuadtreeMesh -> ImmersedGeometry -> QuadtreeMeshRefiner/provider -> FaceSegmentGraph -> CutCellMetrics -> cell-centered Poisson`，这是共享层最可信的参考路径。但 FaceSegmentGraph 表达的是 cell-face 通量，不能直接覆盖旧 Hydro 的 edge-centered momentum 双控制体；若保留交错布局，需要在 Solvers 层新增 `StaggeredDualGraph/DualControlVolumeStencil`。同理，固体 EB 的单侧 cut metrics 与两侧都参与求解的内部介质界面不同，后者需要 `InternalInterfaceMetrics/JumpCondition`。

## 3. 从旧工程恢复的行为基线

旧 `D:\Chen\IGMPP-ILET\src\Main.cpp` 的电喷雾流程可作为第一版耦合顺序的参考，而不是目标架构模板。旧输出和约 232 MiB/份的 checkpoints 表明确实至少运行到 loop 619；但发射、collector/extractor 和推力记录始终为零，因此只能用作链路/场趋势参考，不能视作完整物理真值。

可保留为 legacy trend 的三个观测点如下（不设为解析验收值）：

| loop | time | max div(U) | max speed | volume |
|---:|---:|---:|---:|---:|
| 1 | 0.000386838 | 4.33814e-05 | 9.05641e-05 | 1.18108 |
| 262 | 0.100047 | 5.01392e-02 | 5.7368e-02 | 1.18306 |
| 524 | 0.200088 | 4.92662e-02 | 6.06998e-02 | 1.18591 |

loop 1 的时间从初值 `dt=5e-6` 跳到约 `3.87e-4`，也直接暴露了下述旧时间记账错误。

初始化顺序：

1. 创建 mesh、ImmersedGeometry、Hydrodynamics、SharpInterface、Electricity 和粒子种类；
2. 注册部分 AMR 字段，循环 `initializeLevelSets -> adapt`；
3. `markPhases -> reinitializeLevelSets -> checkCutInfo`；
4. 更新流体材料属性；
5. 初始化 SharpInterface；
6. 初始化 Electricity；
7. 先做一次 `hydrodynamics.solveRZ(dt, preDt, ...)` 作为流场启动。

每个时间步的旧顺序：

1. 将 staggered 速度插值到顶点；
2. 用速度推进 level sets；
3. AMR adapt 并迁移登记字段；
4. 重新插值速度，重初始化 level sets，刷新 cut info；
5. 更新材料属性；
6. 更新 SharpInterface 的节点、历史量和 stencil；
7. Electricity 求解表面电荷、电场、Maxwell 应力及粒子；
8. Hydrodynamics RZ 求解；
9. 由流速、毛细、介电松弛和粒子发射限制计算下一步 `dt`；
10. 写诊断、场和 cache。

旧时间记账本身有 bug：先用 `dt = compute_dt(...)` 覆盖当前步长，再执行 `t += dt; preDt = dt`，导致时间推进使用“下一步 dt”，且 `preDt` 总等于 `dt`。新 `TimeController` 必须原子地区分 `dtUsed` 与 `dtNext`：先 `t += dtUsed`、保存 `previousDt = dtUsed`，最后提交 `dtNext`。

恢复版先保持这个**显式分裂顺序**，确保可比较；后续再增加可选的 E↔flow 子迭代。旧流程的以下内容不能进入通用层：交互式重启询问、固定 level set 1、EMI/BF4、特定电极/收集器统计、硬编码输出探针和场蒸发时间步约束。旧 `physicalField.dict` 在源码中没有读取者，不能视为已有 BC 系统；应把旧硬编码行为先固化成测试，再接入新的 BC registry。

## 4. 目标架构和所有权

```text
quantum_main
  -> RunController / CaseLoader
       -> SimulationContext（唯一创建 mesh/geometry/units/time/IO）
            |- QuadtreeMesh + ImmersedGeometry
            |- Geometry/Boundary/Phase registry
            |- FieldRegistry + topologyRevision + geometryRevision
            |- GeometryViewCache（按 active-region 选择器缓存）
            |    `- SharpInterface views / FaceSegmentGraph / CutCellMetrics
            |- OutputManager + CheckpointManager
            |
            +-> HydrodynamicsModule
            |     |- FlowState
            |     |- IncompressibleFlowSolver (Planar2D | AxisymmetricRZ)
            |     |- MaterialModel / FlowBC / body & interface force ports
            |
            +-> ElectricityModule
            |     |- ElectrostaticState
            |     |- InterfaceElectrostaticSolver
            |     |- SurfaceChargeTransport（可选、与 Poisson 分开）
            |     `- MaxwellTractionProvider
            |
            `-> ElectrosprayApp（只负责编排）
                  |- CouplingSchedule / TimeStepController
                  |- IonEvaporationModel + SpeciesConfig
                  |- ParticleEmission/Collection（可选）
                  `- electrospray diagnostics
```

对应 CMake 目标最终拆为 `quantum_core -> quantum_solvers -> quantum_{hydrodynamics,electricity,plasma} -> quantum_runtime -> QUANTUM`。用链接依赖阻止 Core 反向引用 Solver/Module，也让各 physics library 的模板在自己的显式实例化目标中真正编译。

### 4.1 `SimulationContext` 的最小职责

建议新增 Runtime 层对象（名字可取 `SimulationContext` 或 `SimulationDomain`）：

```cpp
struct SimulationContext {
    CasePaths paths;
    std::shared_ptr<const UnitsContext> units;
    std::shared_ptr<const NormalizationContext> normalization;

    QuadtreeMesh mesh;
    ImmersedGeometry<QuadtreeMesh> geometry;
    GeometryRegistry geometryRegistry;   // name/id -> phase/boundary/interface handle
    FieldRegistry fields;

    std::uint64_t topologyRevision = 0;
    std::uint64_t geometryRevision = 0;
    SimulationClock clock;               // step, t, dt, previousDt
    OutputManager output;
    CheckpointManager checkpoint;

    GeometryViewCache derived;           // 按区域选择器缓存 interface/graph/metrics
};
```

约束：

- App/Runtime 创建一次；Module 构造函数接收引用或窄接口，禁止自己读取 mesh/geometry dict。
- `topologyRevision` 每次网格拓扑变化递增，`geometryRevision` 在移动界面/边界改变时递增。所有矩阵、stencil、cell-label 映射和粒子搜索缓存都记录依赖的 revision，不匹配就必须重建。
- 不假定全系统只有一个 active-cell graph：PIC 活动域、流体相域和全域静电可能不同。GeometryViewCache 以 region selector + 两个 revision 为 key 管理 FaceSegmentGraph/metrics。
- Core 不引用 Hydrodynamics/Electricity/PIC 名字；物理模块通过 provider 提交细化指标、源项和边界策略。
- `RunConfig` 可保留共享 solver/units 配置，但不再充当隐式 PIC 配置抽取器。

### 4.2 Module 生命周期

不要求一次改完所有现有 PIC 模块。先为 Fluid/Electricity 引入新生命周期，并用 adapter 保留旧 `run()`：

```cpp
struct StepInfo {
    int step;
    double time;
    double dt;
    double previousDt;
    std::uint64_t topologyRevision;
    std::uint64_t geometryRevision;
};

class IPhysicsModule {
public:
    virtual Status configure(const ModuleConfig&, const SimulationContext&) = 0;
    virtual void registerFields(FieldRegistry&) = 0;
    virtual Status initialize(SimulationContext&) = 0;
    virtual StepReport advance(SimulationContext&, const StepInfo&) = 0;
    virtual Status onMeshChanged(SimulationContext&, std::uint64_t oldTopologyRevision) = 0;
    virtual Status onGeometryChanged(SimulationContext&, std::uint64_t oldRevision) = 0;
    virtual Status save(CheckpointWriter&) const = 0;
    virtual Status load(CheckpointReader&) = 0;
};
```

`Status/StepReport` 至少携带阶段名、线性/非线性残差、迭代数、守恒量、是否有限、失败原因；底层不得 `exit()`。

### 4.3 语义绑定而不是魔法序号

字典可继续保存整数 ID，但算法入口使用启动时解析出的 handle：

- `PhaseHandle liquid` / `PhaseHandle gas`；
- `BoundaryHandle emitterWall` / `outlet` / `axis`；
- `InterfaceHandle freeSurface`；
- handle 中可缓存 level-set index、正负侧 phase、boundary IDs 和解析时的 geometry revision。

所有名字在初始化时一次性校验：不存在、重名、正负相不一致、RZ 轴不在 `r=0`、电学区域无介电率等均立即报配置错误。生产算法中禁止直接出现 `levelSets[1]` 及“东侧就是压力出口”等位置假设。

### 4.4 AMR 与字段迁移协议

用 `FieldRegistry` 取代每个 App 手工维护的四组字段指针。每个字段注册：

- 所有者、名字和位置（vertex/cell/H-edge/V-edge/interface-node/particle）；
- 角色（primary/conserved/history/derived/scratch）；
- refine/coarsen 策略（插值、保守限制、重算或丢弃）；
- checkpoint 策略与单位；
- 当前 topology/geometry revision。

至少提供 `Interpolate`（连续标量）、`ConservativeCellAverage`（质量/体电荷）、`FaceFluxConservative`（交错速度/通量）、`LevelSetTransfer`（当前及历史 LS）、`GeometryProject`（表面电荷等界面节点场）、`Recompute`（E/梯度/曲率/metrics）和 `Discard`（临时工作区）策略，不能用同一种插值处理所有字段。

统一拓扑变更事务：

1. 各模块提交 refinement indicator，Mesh 层只合并为 target levels/marks；
2. 冻结并校验 registry，保存旧拓扑及 primary/history 字段；
3. adapt；
4. 按策略迁移守恒量、速度、压力初值、level-set 历史和表面电荷；
5. 重新初始化 signed distance，刷新梯度、cut info 和 phase 标记；
6. 重建 interface nodes/stencils、FaceSegmentGraph、CutCellMetrics；
7. 递增相应 revision，调用所有模块 `onMeshChanged/onGeometryChanged` 重建矩阵/cache；
8. 做有限性、质量/电荷变化和悬空 label 检查。

若界面移动但网格未 adapt，则跳过 2--4，只增加 geometry revision，并执行 5--8 中与几何有关的部分。

### 4.5 配置文件所有权

目标发布 case 建议保留：

```text
system/
  run.dict          # App、时间、耦合、输出、重启
  mesh.dict         # 纯网格与 AMR 限制
  geometry.dict     # Domain、线/水平集、phase、boundary 实体与名称
  units.dict        # 跨模块单位与归一化（建议合并旧 normalization.dict）
  solvers.dict      # 跨模块线性求解器默认值
  hydro.dict        # 仅流体材料、BC、界面和数值选项
  electricity.dict  # 仅电学区域、BC、jump、表面电荷选项
  pic.dict/mcc.dict # 仅启用相应模块时存在
```

这需要同步修订 ArchitectureRules v2。为了不破坏已调通的 PIC：迁移期仍接受 `pic.dict` 内嵌的 Units/Normalization/Solvers，但若共享文件与内嵌块同时出现，应报“重复来源”而不是静默按优先级覆盖。`app.dict`、重复的 `electrospray.dict` 和 `system/example.jpg` 从发布 case 移除；图片可放 `docs/assets/`。

当前 `DictionarySolvers.cpp` 主要只解析 `PoissonCellCentered`。建议让 `solvers.dict` 提供有名字的线性求解 profile（如 `hydroMomentum`、`hydroPressure`、`interfaceElectrostatic`），Hydro/Electricity dict 只引用 profile 并保留少量物理算法选项。这样线性 backend、容差、预条件和失败策略可以共用，但材料/BC 不会倒灌到 Solvers 配置。

`hydro.dict` 通过 phase/boundary/interface 名字表达物理，例如：

```text
Hydrodynamics
{
    coordinateSystem axisymmetricRZ
    timeIntegration variableStepAB2Projection

    phases
    {
        liquid { geometryPhase liquid; density 1240 [kg/m3]; viscosity 0.038 [Pa.s]; }
        gas    { geometryPhase vacuum; density 1.184 [kg/m3]; viscosity 1.849e-5 [Pa.s]; }
    }

    boundaries
    {
        emitterWall { geometryBoundary emitterWall; velocity noSlip; pressure zeroGradient; }
        outlet      { geometryBoundary outlet; velocity zeroGradient; pressure fixed 0 [Pa]; }
        axis         { geometryBoundary axis; velocity axisSymmetry; pressure zeroGradient; }
    }

    interfaces
    {
        freeSurface
        {
            geometryBoundary liquidVacuum;
            surfaceTension 0.0452 [N/m];
            advect true;
        }
    }
}
```

`electricity.dict` 只描述通用电学：区域介电率/电导率、按 boundary name 的 Dirichlet/Neumann 条件、界面 jump 与是否启用表面电荷输运。离子质量、活化能、发射律、收集器和电喷雾专属诊断放入 `run.dict/Electrospray` 的 model/species 块。

## 5. 通用化后的物理边界

### 5.1 HydrodynamicsModule 应拥有

- `FlowState`：`p, rho, mu, ux/uy` 及明确命名的历史/预测字段；
- 平面 2D 与轴对称 RZ 的不可压推进策略；
- phase 驱动的材料场构建；
- boundary handle 驱动的速度/压力条件；
- 体力和界面 traction 输入端口；
- 速度插值、散度、质量和线性求解诊断；
- 网格变化后的字段迁移与离散缓存重建。

它不应知道“电喷雾”“EMI/BF4”“阳极/阴极”，也不应直接读取 Electricity 对象。

### 5.2 ElectricityModule 应拥有

- `ElectrostaticState`：potential、E、可选体电荷；
- 多区域 `epsilon/conductivity`；
- domain/embedded/interface 边界条件；
- 两侧电势和法向通量 jump；
- 可选、独立的 `SurfaceChargeTransport`；
- 输出通用 `MaxwellTraction`，供任意流体 App 连接。

制造解残差属于 tests；粒子发射/推进属于独立模型或 PIC/Lagrangian 服务，不属于静电 Poisson 求解器。

### 5.3 ElectrosprayApp 应拥有

- 显式分裂或子迭代的耦合顺序；
- 电喷雾专属材料组合与电极角色检查；
- 场蒸发模型、species、粒子 emission/collection；
- 组合式时间步限制（对流、毛细、电荷松弛、发射等）；
- 推力、质量流率、电流等应用诊断；
- checkpoint/output 调度。

几何形状仍完全来自 `geometry.dict`；App 不写坐标公式覆盖 level set。

## 6. 分阶段实施与退出条件

各阶段都应形成可独立审查、可回退的小提交。后续阶段不得以跳过前一阶段的验收测试为代价并行“向前冲”。

### P0：建立真实编译基线和失败语义

工作项：

1. 新增编译探针，显式实例化 QuadtreeMesh 上的 Hydro XY/RZ 全部生产入口以及 Electricity 的 initialize/solve/cache 入口；加入 Release + OpenMP CI/CTest。
2. 用局部变量 reduction 后回写成员，修复 C3028；把所有 `supSonCellIndex` 改为整数并开启转换告警审计。
3. 把 SparseLU/配置/拓扑错误从 `exit(0)` 改成失败 `Status/StepReport`；所有迭代器增加 `maxIterations`、非有限检查和未收敛返回。
4. 保存本次无 OpenMP 一步探针输出，作为“能执行”的诊断基线，不作为物理 golden result。
5. 当前 QUANTUM 目录没有可用的 `.git` 历史；正式改算法前应建立可追踪的源码基线，并记录旧 `IGMPP-ILET` 的来源版本/校验值。

退出条件：

- 正式 `/openmp:llvm` 配置确实生成所有 Fluid/Electricity 求解函数；
- 编译探针无 C3028、索引窄化警告；
- 线性求解失败返回非零且 CTest 能观察到，不再出现 `exit(0)`；
- 不改变现有 PIC 回归结果。

### P1：修通 case、字典和最小单模块入口

工作项：

1. 修正白名单/配置装载：允许启用模块声明的 dict 和共享 dict；unknown file 给出精确诊断。
2. 统一共享配置唯一来源，保留 PIC 兼容读取并检测重复定义。
3. 为旧几何提供一次性迁移工具或明确报错，把 domain 移到 `geometry.dict/Domain`，把旧 `line N` boundary 转为带 ID/name 的 `boundary N` 实体。
4. 清理两个 case 中废弃/重复文件；当前 electrospray mesh 域为 `(0,-1)-(1,1)`，geometry 点却延伸到 `(5,16)`，是验证域与旧几何的混合物，不能直接当作验证 oracle。
5. 分建尺度一致的 `hydrodynamics_smoke`、`electrostatics_interface_mms` 和 `electrospray_legacy_regression`；前两者验证通用模块，后者只对照旧链路/字段趋势。Hydrodynamics 单模块实际执行可配置的 1/N 步，而非只写 `rho/mu`。

退出条件：

- 两个新 case 均能从 `QUANTUM --case <path>` 进入预期 App/Module；
- `--check-case`（或等价测试）能报告所有 phase/boundary/interface 绑定结果；
- 缺少 Domain、重复配置来源、未知 boundary name 均在求解前失败；
- Hydrodynamics smoke 的日志明确包含动量解和压力投影，不以“生成输出文件”代替成功。

### P2：引入共享 SimulationContext、FieldRegistry 和拓扑事务

工作项：

1. Runtime 创建并持有 mesh/geometry/units/clock/output；Hydro/Electricity 的构造与模块入口改为注入 context。
2. SingleModuleApp 对多个模块只构建一份 context；短期仍可顺序运行，但必须共享同一 mesh 地址和 topology/geometry revisions。
3. 建立 GeometryRegistry，将名字解析为 typed handles；去除新代码中的魔法 level-set/boundary 序号。
4. 建立 FieldRegistry；迁移 level-set 全部历史、Hydro 全部 primary/history、Electricity potential/charge 及 interface fields。
5. `QuadtreeMesh::adapt` 去除 `magGradU` 概念，改为 target-level/mark 输入；Hydro/PIC 各自实现 refinement provider。
6. 拆出 MovingEBController 的纯几何运动接口；保留现有 PIC 行为的 adapter，避免在本阶段改动其物理。

退出条件：

- 联合加载 Hydro+Electricity 时只构造一次 mesh/geometry；
- 强制 refine/coarsen 后所有已注册字段 revision 一致，无悬空 label；
- 常量场迁移保持常量，守恒场的总量误差在定义阈值内；漏注册 primary/history 字段会硬失败；
- 现有 PIC cases 仍能通过兼容 adapter 运行。

### P3：恢复通用单相/两相 Hydrodynamics

工作项：

1. 先保持旧 staggered 离散，拆出 `FlowState`、`MaterialModel`、`FlowBoundaryConditions` 和 `FlowForcing`。
2. 为 edge-centered momentum 建立 `StaggeredDualGraph/DualControlVolumeStencil` 适配层；压力/散度可复用 FaceSegmentGraph，不能假设两种控制体完全等价。
3. 恢复 `previousDt` 或将配置限制为固定步长；用时间积分单元测试验证系数。
4. 用 phase 映射生成材料场；solid/electrode mask 与 liquid/gas 混合规则分开，禁止两个 level set 相加截断的 case 逻辑。
5. 将四周硬编码速度/压力分支改为 boundary handle 条件；RZ 增加 axis 条件和 `r -> 0` 专门测试。
6. 修复并按 boundary view 重建 no-slip nodes/stencil；直接力迭代设置容差、上限、松弛、残差报告。必要时先提供强制边界/ghost 方案作对照。
7. 分离平面/RZ 的 metric policy，先减少新增重复代码；旧数值主体在回归稳定前不强行合并。
8. 删除或移入历史参考文件的大段注释重复实现。

退出条件：

- 静止场保持静止且无 NaN/Inf；
- 平面 manufactured/通道/腔体之一通过网格收敛测试，投影后散度随容差/网格合理收敛；
- 四叉树粗细交界的速度/压力测试通过；
- RZ 轴线上有限、零法向穿透，最小基准达到记录的误差阶；
- 无滑移测试确认节点数非零、残差收敛；未收敛返回失败而非死循环。

### P4：恢复移动自由界面、表面张力和 AMR

工作项：

1. 将 level-set advection 作为 Geometry 服务，由 Flow velocity provider 驱动；界面退移速度通过可选 provider 注入。
2. 明确 `advect -> adapt/remap -> reinitialize -> cut/phase -> interface stencil` 的唯一生命周期。
3. Surface tension 只依赖界面几何和材料参数，输出通用 traction/body-force view。
4. 对界面质量、体积、signed-distance 质量、曲率和寄生流建立诊断。

退出条件：

- 平移圆/形变界面测试中体积误差受控，AMR 前后无拓扑悬空；
- 静态液滴满足 Laplace 压差，寄生速度随网格细化下降；
- 关闭退移 provider 时无隐藏 level set 1 行为；开启常量退移时速度/体积变化符合解析关系。

### P5：拆分并验证通用界面静电

工作项：

1. 先建立“旧物理 operator 与当前 MMS operator”逐项映射测试，把解析源/边界/jump 完整移到 verification policy，避免只删除 residual 输出却仍在求解 Helmholtz MMS。
2. 把 Electricity 拆成 `InterfaceElectrostaticSolver`、`SurfaceChargeTransport`、`MaxwellTractionProvider`；粒子代码移出。
3. 用 region/interface/boundary handles 替换所有 `[1]`；平面与 RZ metric 明确选择。
4. 第一版可把旧 RZ IIM 包成 `LegacyIIMInterfacePoisson` 以恢复物理行为，但对外只暴露通用输入/输出。
5. 在 Solvers 层实现或扩展真正的 two-sided variable-coefficient interface Poisson：支持 `epsilon^- / epsilon^+`、potential jump、normal-flux jump 和表面电荷。
6. 评估复用 FaceSegmentGraph/CutCellMetrics 与 cell-centered Poisson 的装配/线性求解部分；不得把固体 EB 的单侧条件冒充内部介质界面。

退出条件：

- 平面与 RZ 两介质 manufactured solution 给出 potential/E/jump 残差及网格收敛；
- 常介电率退化到普通 Poisson，并与现有 solver 结果一致；
- 界面方向翻转只改变约定中的 jump 符号，不改变物理解；
- normal-derivative/jump residual 随网格或求解迭代收敛，不能重现探针中长期停在约 `7.67e-2` 而仍报告成功的行为；
- `initialize` 不再含硬编码解析解，`Electricity` 不再创建/推进粒子。

### P6：通用电流体耦合

工作项：

1. 定义稳定的数据接口：Electricity 输出界面 Maxwell traction/体电力，Hydro 只消费 `FlowForcing`；Flow 输出速度，SurfaceChargeTransport 和 Geometry 消费速度 view。
2. ElectrosprayApp 先实现旧显式分裂顺序；配置可选耦合子迭代、松弛和最大次数。
3. 统一时间层：所有模块使用同一 `StepInfo(dt, previousDt)`，明确电应力用于哪个时间层。
4. 建立质量、体积、表面/体/粒子电荷和电流平衡表；每步输出机器可读诊断。

退出条件：

- 关闭电力时结果退化到 Hydro-only；相同介电率且零电荷时界面电应力跳跃为零；
- 静态介质界面/液滴的 Maxwell traction 与解析或文献基准吻合；
- 10--100 步 coupled smoke 全字段有限，投影、界面电荷和耦合残差均有界；
- 调换无关 boundary/level-set 的输入顺序不改变结果。

### P7：电喷雾模型回归 App 层

工作项：

1. 将 EMI/BF4、温度、活化能、导电率、特征尺度、电位和发射参数迁入 electrospray 配置；Core/Base 只保留数学/普适物理常数（若确需）。
2. 独立实现 `IonEvaporationModel`，输入局部场、温度和界面状态，输出连续电流/质量源；以单元测试验证数量级和正负离子分支。
3. 独立实现 charge-to-macro-particle 转换，移除 `particleNum = 0`；定义随机抽样、宏粒子权重、守恒回写和可复现实验种子。
4. 复用已有 Lagrangian/PIC 粒子推进能力，而不是在 Electricity 中维护第二套 pusher。
5. 把旧 `compute_dt` 拆成可组合 constraints；无粒子时发射限制返回 infinity，每个限制在日志中单独报告。
6. 从 `geometry.dict` 读取真实喷口/电极/自由界面；删除 App 中的圆和 `y+7.98`。

退出条件：

- 禁用 evaporation/emission 时 EHD case 仍可运行；
- 启用连续蒸发但不生成粒子时，表面电荷/质量源守恒；
- 启用粒子后宏粒子电荷与从界面扣除的电荷在容差内相等；
- 电喷雾 smoke 完成多个输出周期，电流/质量流率/推力的单位和符号经过检查。

### P8：重启、发布回归与遗留清理

工作项：

1. Checkpoint 增加 magic、schema version、case/config hash、坐标系、topology/geometry revisions、字段清单和单位；停止裸写整个对象内存布局。当前 Electricity cache 已有 `ELC2/version 2`，旧 cache 没有 header，二者应由显式 legacy importer 转换，主 reader 不猜格式。
2. 建立从固定 checkpoint 继续 N 步与连续运行的对比测试。
3. 将旧行为需要保留的 adapter 标注 deprecation；只有新回归全部覆盖后才删除 legacy IIM/旧接口。
4. 更新 ArchitectureRules、用户 case 指南和字典 schema；删除重复 case 文件和失效注释。

退出条件：

- 不兼容 cache 明确拒绝并说明原因；兼容 cache 的续算结果在规定容差内等同连续运行；
- Release/OpenMP、Debug 和无 OpenMP 至少各有一个 fluid/EHD CI job；
- 发布 cases 通过 schema lint，`system/` 无图片、脚本和废弃 dict；
- 下述测试矩阵全部进入 CTest/CI。

## 7. 最小验证矩阵

| 层级 | Case/Test | 必查量 |
|---|---|---|
| 编译 | `hydro_template_instantiation` XY/RZ + OpenMP | 所有生产模板生成；零关键告警 |
| Schema | valid/invalid fluid dictionaries | Domain、名字绑定、重复来源、单位、坐标系 |
| Mesh | constant/linear field refine-coarsen | 精确性、守恒、revision、无悬空 label |
| Hydro | zero field / manufactured / channel or cavity | 速度、压力、投影散度、线性残差、收敛阶 |
| Hydro-AMR | 粗细界面压力投影 | flux balance、速度连续性、散度 |
| RZ | axis regularity / manufactured | `r=0` 有限性、轴条件、带 metric 的守恒 |
| Free surface | static circle/drop | Laplace 压差、曲率、寄生流、体积 |
| Geometry motion | translating/deforming level set | 形状/体积误差、重初始化质量 |
| Electrostatics | planar/RZ dielectric interface MMS | potential、E、jump residual、收敛阶 |
| Surface charge | zero/source/advection manufactured | 表面电荷守恒和时间收敛 |
| EHD | static Maxwell-stress benchmark | traction 符号/大小、零场退化 |
| Coupled smoke | 10--100 步、含一次 AMR | finite fields、质量/电荷平衡、残差有界 |
| Electrospray | evaporation off/on、particles off/on | 模型退化关系、发射电荷/质量守恒 |
| Restart | continuous vs checkpoint restart | 字段、历史时间层、粒子 RNG、diagnostics |

验收阈值不要凭当前 smoke 数字拍定。先用制造解/解析解做网格序列，记录观测阶和离散误差，再在测试中固定“绝对/相对误差 + 收敛阶 + 守恒”三类阈值。

## 8. 文件级改造地图

| 区域 | 第一批修改 | 后续归宿 |
|---|---|---|
| `src/Apps/Runtime/RunController.cpp` | 修白名单、消除 PIC 专属共享配置入口、case 预检 | CaseLoader + schema/manifest 驱动 |
| `IApp.h` / `IModule.h` / `RunConfig.h` | 注入共享 context 的 adapter | 明确 App 编排与 Physics lifecycle |
| `SingleModuleApp.cpp` | 一次构建 context，多模块共享 | 单模块调试器而非伪耦合器 |
| `ElectrosprayApp.cpp` | 恢复受测时间循环，删除硬编码几何 | 只保留 coupling/model/output 调度 |
| `DictionaryGeometry.cpp` | 旧 schema 明确迁移/报错，完整实体校验 | GeometryRegistry typed handles |
| `QuadtreeMesh.h/_impl.h` | adapt 接受 marks；禁止 `magGradU` 参数 | 通用拓扑操作 + FieldRegistry transaction |
| `QuadtreeMeshRefiner` | 接收模块 provider | 多物理 refinement policy 合并器 |
| `SharpInterface.h/_impl.h` | 恢复并测试 no-slip，按 handle 建 view | 通用界面几何/stencil 服务 |
| `MovingEBController.h` | 提取纯几何运动核心 | impact/erosion 下沉到相应 Module/App |
| `Modules/Plasma/ExplicitPIC/Geometry/LevelSetQuery.h` | 提取无 PIC 依赖的采样部分 | `Core/Geometry/LevelSetSampler` |
| `Solvers/IO/KVBlockParser.h` | 统一 Hydro/Electricity 的重复行解析与错误定位 | 下沉 `Core/IO`，所有字典共用 schema 基础 |
| `src/Solvers/IO/DictionarySolvers.cpp` | 增加命名 linear-solver profiles 与严格校验 | 多物理共用 backend/容差/失败策略 |
| `Hydrodynamics*.h/*_impl.h` | 模板编译修复、失败返回、`preDt`、去魔法索引 | FlowState + solver/material/BC/forcing |
| `HydrodynamicsModule.cpp` | 实际单步/多步入口，使用共享 context | 可独立验证的 physics module |
| `Electricity*.h/*_impl.h` | MMS/粒子拆出，去 `[1]`、配置化区域和 BC | electrostatic/surface-charge/traction 三层 |
| `ElectricityModule.cpp` | 移除临时 Hydro 和私有 mesh/geometry | 可独立验证的 electrostatics module |
| `Core/Base/Base.h` | 标记并迁移 case 常数 | 仅保留真正通用常数/基础设置 |
| `Solvers/Poisson/*` | 复用装配基础并补 two-sided interface solver | 单一通用 Poisson 家族 |
| PIC 私有 `PicCutCellMetrics/EdgeCutLengthHelper` | 与 Solvers 版本做行为比对 | 删除无消费者的重复实现，避免漂移 |
| `cases/hydrodynamics`, `cases/electrospray` | 迁移 schema、去重复/非字典文件 | 发布级回归 case，非硬编码输入壳 |
| CMake/CTest | 显式模板、schema、fluid/EHD tests；拆 physics library targets | 依赖方向与 CI 质量门 |

## 9. 实施决策与禁止事项

已确定的方向：

- 旧工程是**行为参考**，不是可整包复制的架构。
- 第一版恢复保留旧数值内核；架构接口先稳定，离散替换后置。
- Runtime/App 拥有共享网格几何；Module 拥有自己的状态和物理配置；Case 只提供数据。
- 几何实体用名字/ID 解析成 handle；不允许算法依赖输入排列顺序。
- App 决定耦合与时间推进；Hydro/Electricity 不互相构造或读取对方字典。
- 所有拓扑/几何相关 cache 都受 revision 管理；所有 primary/history 字段都经 registry 迁移。

明确禁止：

- 直接把旧 `Main.cpp` 粘回当前 main；
- 只扩大白名单后就宣布 fluid case 已恢复；
- 用“程序返回 0/写出了 dat 文件”替代数值验收；
- 在通用 Electricity 中保留 EMI/BF4、场蒸发或粒子 pusher；
- 未做两介质 jump MMS 就用当前 PIC Poisson 替代 IIM；
- AMR 后靠调用者记忆手工补字段列表；
- 求解失败 `exit(0)`，或无上限的非线性/直接力循环；
- 在没有测试保护时同时合并 XY/RZ、改离散、改单位和改耦合顺序。

## 10. 建议立即执行的首个开发批次

实施记录：本批次已接入正式代码和 CTest，具体支持边界、命令及实测结果见第 12 节；这不等于 P0/P1 全部完成，也不代表完整 EHD 已恢复。

首批只做“让失败可见、让最小 Hydro 真跑”，控制范围如下：

1. 把临时模板探针正式化并加入 CMake/CTest；
2. 修复 Hydro OpenMP reduction、整数索引、solver failure return；
3. 修正 Runtime 对 hydro/electricity/shared dict 的契约，并增加 case schema 预检；
4. 新建干净的 planar 与 RZ Hydro smoke cases，不沿用电喷雾常数和陈旧几何；
5. HydrodynamicsModule 读取 `steps/dt/coordinateSystem`，真实调用 1--3 步 solve；
6. 加入 zero-flow、projection divergence、RZ-axis 三个最小测试；
7. 保持 Electrospray 时间循环关闭，直到 P0/P1 的 Hydro 单物理门通过。

完成这个批次后，项目状态会从“主构建假绿、case 不可达”变成“流体入口真实编译、错误可观察、单物理可回归”。随后再进入共享 SimulationContext/AMR 字段注册；这比直接恢复电喷雾大循环更容易定位数值和接口问题。

## 11. 完成定义（最终）

“（电）流体模块已实现”必须同时满足：

1. Hydro 与 Electricity 可独立运行，也可在同一 `SimulationContext` 中组合；
2. 平面/RZ、uniform/quadtree、固定/移动界面均有自动回归；
3. 材料、BC、界面和源项由 phase/boundary/interface 语义配置，不含魔法序号；
4. EHD 通过解析/制造/守恒基准，不只通过 smoke；
5. Electrospray 的几何、材料、species、发射律和电位全部来自 case 配置；
6. AMR 和 restart 保存所有当前/历史状态，拓扑/几何 cache 不会跨 revision 误用；
7. 任何线性/非线性失败都返回可诊断的非零状态；
8. 现有 Plasma 回归不退化，且公共 Geometry/Mesh 不反向依赖特定物理。

达到以上条件后，水平集几何与四叉树网格才真正成为 Hydro、Electricity、Plasma 可共同依赖的基础，而 Electrospray 只是这些通用能力之上的一个可替换应用工作流。

## 12. 首批交付：让失败可见、让最小 Hydro 真跑（2026-09-09）

### 12.1 已实现的范围

- 正式目标 `quantum_ehd_template_compile` 实例化 Hydro XY/RZ 的新旧推进入口及缓存入口、Electricity 的初始化/推进/缓存入口，使用与生产目标相同的编译选项（本机 MSVC `/openmp:llvm`）。Electricity 仅为编译覆盖，不是电流体物理验收。
- 修复 Hydro 成员变量 OpenMP reduction 和浮点辅助单元索引；新 `HydroSolveReport` 返回阶段、原因、动量/压力相对线性残差、投影前后散度及最大速度。检查矩阵、右端项、分解、解及真实残差，失败不再 `exit(0)`。旧四参数入口保留，失败抛异常；直接力循环增加上限及非有限值检查，但本批案例不启用直接力。
- 新增不需要 `SharpInterface` / `ImmersedGeometry` 的固定步长单相入口，以及显式常密度/动力黏度初始化。独立路径不调用电喷雾材料、蒸发或表面力逻辑；成功步之后更改 `dt` 会被拒绝。旧可变步长驱动尚未修复，本批不重开该驱动。
- `HydrodynamicsModule` 真正执行 1--3 步动量方程和压力投影，读取 `steps/dt/coordinateSystem/density/viscosity`；可选零初速或平滑非零初速。每步写 CSV 诊断，可选写顶点 `.dat` 场。
- 最小闭域压力改为齐次 Neumann 加 `p[0]=0` 参考，匹配四周零法向速度。新诊断从最终速度重新计算**所有单元**的离散通量散度，包括压力参考单元、边界和 RZ 轴邻近单元。旧几何耦合入口保持原压力边界选择。
- Runtime 白名单放行 Hydro、Electricity 及已有读取逻辑所支持的共享字典。`Solvers/Units/Normalization` 同时出现在独立字典和 PIC 内嵌块时明确报冲突，PIC-only 原有内嵌路径保留；气体/反应等其他共享配置的来源重构不在本批。
- `--check-case` 仅支持 `SingleModule` 中唯一的 Hydro/Hydrodynamics：验证字典、Domain、网格尺寸/容量和支持范围，不创建网格、不推进、不创建 `output/` 或 `cache/`。正常和预检共用模块列表解析，支持 inline / 跨行列表；执行前先检查全部模块名。

### 12.2 首批契约与明确限制

发布案例为 [`hydrodynamics_smoke_planar`](../cases/hydrodynamics_smoke_planar/README.md) 和 [`hydrodynamics_smoke_rz`](../cases/hydrodynamics_smoke_rz/README.md)，各仅含 `run.dict`、`hydro.dict`、`geometry.dict`、`mesh.dict` 四个输入。

- **只接受数字代码单位**。不支持 SI 换算，也不使用共享 `Solvers/Units/Normalization` 参数；独立文件或 PIC 同名内嵌块均明确拒绝，避免配置被忽略或输出虚假的 SI 报告。通用共享配置仍留待后续接入。
- 几何仅为 `Domain { min (...) max (...) dimension 2 }`；RZ 使用 `x=r, y=z`，要求径向下界为 0。不接受 Phase、Boundary、LevelSet、Line、运动或旧几何 schema。
- 复用 `QuadtreeMesh` 的**固定均匀网格**；`adaptiveMesh false` 显式声明，不进行初始化细化。`divide` 每维为 `[4,256]` 内整数，`maxCellNum` 至少容纳 `(Nx+1)*(Ny+1)` 个顶点且不超过 1000000，`maxInitialIterations` 为 0。坐标与网格分量支持逗号或空格分隔。
- 支持闭域零法向速度；切向黏性边界仍沿用旧离散。**不声称已恢复通用 no-slip、入口/出口或收敛阶正确性**。
- 单相、固定材料、固定 `dt`、1--3 步；无 AMR、自由界面、几何运动、电力或 EHD 多模块耦合。Hydro 必须独立运行，Electrospray 时间循环继续关闭。
- 求解失败时返回非零，CLI 立即中止。**失败步不保证回滚**；库调用者须重建对象或恢复全部主场和历史场后重试。尚未到达的诊断值为 `NaN`，CSV 的失败行表示尝试推进的时间，而非已接受时间步。
- `linearResidualTolerance` 是相对线性残差门限（右端项为零时用绝对残差），不是离散误差或散度收敛标准。`max_speed` 汇总单元插值速度模与交错面速度分量的最大值；`axis_max_radial_velocity` 仅对 RZ 有轴条件含义。

### 12.3 构建、运行与回归命令

从仓库根目录执行；当前 CMake 将主程序输出到 `bin/QUANTUM.exe`：

```powershell
cmake --build build --config Release --parallel 8
ctest --test-dir build -C Release -R "^hydro_" --output-on-failure
.\bin\QUANTUM.exe --case cases/hydrodynamics_smoke_planar --check-case
.\bin\QUANTUM.exe --case cases/hydrodynamics_smoke_planar
.\bin\QUANTUM.exe --case cases/hydrodynamics_smoke_rz
```

正常运行将写入所选 case 的 `output/hydro_diagnostics.csv` 和 `hydro_step0001.dat` 等文件；重复运行会覆盖同名结果。Hydro CTest 自动复制输入到 `runs/history/checks/tests/hydro_cases/<test>/run_<id>/`，不改动源 case，每次保存 `test_stdout.log`、`test_stderr.log` 及输出。

### 12.4 实测与验收记录

两个正式 smoke 均使用 12x12 网格、非零初速、`rho=1`、`mu=0.01`、`dt=0.001`，完成 3 步至 `t=0.003`：

| 指标 | planar2D | axisymmetricRZ |
|---|---:|---:|
| 第 1 步投影前最大散度 | `1.6901e-1` | `2.2926e-1` |
| 第 1 步投影后最大散度 | `3.70e-15` | `5.64e-14` |
| 第 3 步投影后最大散度 | `6.25e-17` | `4.77e-16` |
| 三步最大动量相对残差 | `2.70e-16` 以下 | `3.08e-16` 以下 |
| 三步最大压力相对残差 | `1.18e-14` 以下 | `2.29e-14` 以下 |
| 第 3 步最大速度 | `3.0645e-2` | `2.8430e-2` |
| 轴径向速度 | 不适用 | `0` |

非零速度仍被保留，避免“清零整个速度场也通过散度测试”。独立数值回归从返回的交错速度重算散度，要求投影后最大散度不超过 `1e-8`，首步至少下降六个数量级，零流保持至 `1e-13`；这些是首批数值健全性门限，不是网格收敛阶证据。

将 `linearResidualTolerance` 故意设为 `1e-30` 时，真实动量相对残差约 `2.69e-16`，程序以返回码 **3** 退出，CSV 保留 `success=0, stage=momentum-solve` 和超限原因；非法配置/不支持范围返回码 **2**。这区分了正常完成、求解失败和输入失败。

最终验证结果：

- `cmake --build build --config Release --parallel 8` 成功。使用 MSVC/LLVM OpenMP 与 Eigen SparseLU；当前环境未找到外部 SuperLU。遗留 Geometry 缓存/Electricity 模板仍有 `size_t -> int`（C4267）警告，本批未把工程宣称为零警告，也未完成全量 Debug/无 OpenMP 构建矩阵。
- `ctest --test-dir build -C Release -R "^hydro_" --output-on-failure`：**25/25 通过**。覆盖模板编译、零流、非零投影、RZ 轴、非法步长/材料、真实求解残差失败、两个三步案例、预检无写入、inline/跨行模块列表、缺失/重复/不支持的几何配置及共享来源冲突/不支持配置。
- 现有 4 个 Poisson/几何案例与 ImplicitPIC Turner periodic-Y 的 20 步 smoke，使用未改动输入的隔离副本运行，均返回 0。cut-circle MMS 电势 `L2≈2.68e-6`，polarization MMS `L2≈1.67e-7`；其中三个旧测试为 `RegressionHarness mode record`，这里只记作运行与指标记录通过，**不是固定历史 baseline 对比通过**。
- 隔离旧案例产物保存在 `runs/history/hydro/hydro_first_batch_baseline_65c251a4/`；Turner 的兄弟目录碰撞数据按原相对路径复制，未覆盖源案例 output/cache。原第六项 `regression_implicitpic_eb_intercept_smoke` 引用的 `cases/implicitpic_eb_intercept_smoke` 在当前工作区不存在，因此未运行；没有将完整既有 CTest 套件标为全绿，也没有为通过测试伪造该案例。

首批之外的下一步仍为共享 SimulationContext、字段生命周期与独立边界/材料契约，不在本次继续实施。

## 13. 快速接线批：流体—命名水平集—材料闭环（2026-09-09）

### 13.1 优先级调整

用户明确要求“先把流体模块快速接起来跑起来，后面再复合验证”。因此暂不把 Poiseuille / MMS / 三档网格收敛验收作为继续接线的前置条件，也不展开完整 SimulationContext、FieldRegistry 或 AMR 重构。本批交付目标为可配置的持续多步运行，以及现有流体与水平集几何之间的实际数据交换；保留编译、有限值、线性残差、非零失败状态和现有 smoke 回归。

### 13.2 接通的运行链路

1. `geometry.dict` 的命名 LevelSet 经严格的受支持 schema 适配，生成公共 `DictionaryGeometry`，与 Hydro 共用同一个 `QuadtreeMesh` 构造 `ImmersedGeometry`；初始化使用公共水平集算法，不在 App/Module 内写死圆或电喷雾坐标公式。
2. `materialModel twoPhase` 通过 `interfaceLevelSet` 名称绑定一个材料界面，使用配置的负侧/正侧密度、动力黏度和过渡厚度构建材料场。复用现有 tanh 平滑与顶点到单元插值，不读取 `nonDim_*`、Reynolds 或电喷雾材料常数。
3. 使用这些材料场执行现有 XY/RZ 动量与压力投影。`accelerationX/Y` 为代码单位体加速度，装配的体力密度为 `rho * acceleration`；默认均为 0。
4. `advectInterface true` 时，把新速度及前一步速度转为顶点场，调用公共 `ImmersedGeometry::updateLevelSets`，再刷新导数、cut 信息和输出材料。固定拓扑下历史存储无需 resize，派生数据刷新不会覆盖两层水平集历史。
5. CSV 每步记录流体与材料/几何状态；`.dat` 可按 `writeEvery` 输出，并强制保留末步，`writeInitialFields` 可额外输出初始状态。`logEvery` 独立控制进度日志。原单相三步案例输入和已有用户结果不改动。

`steps` 已解除 3 步上限，接受小于 `INT_MAX` 的正整数，仍要求 `dt>0`、终止时间有限且全过程固定步长。原三步 smoke 继续保留作为兼容回归，并不是运行能力上限。

### 13.3 新输入契约及限制

`materialModel` 默认为 `constant`，沿用 `density/viscosity` 与 Domain-only 几何；选择 `twoPhase` 时改用以下参数，不同时接受常材料键：

```text
materialModel twoPhase
interfaceLevelSet "fluidInterface"
densityNegative 2
densityPositive 1
viscosityNegative 0.02
viscosityPositive 0.01
interfaceThickness 0.02
advectInterface true
reinitializeEvery 0
```

此处 `twoPhase` 表示一个水平集分隔的**两种流体材料**，不是固体掩膜，也不是已经通过完整两相流物理验证。所有物理量仍为代码单位。材料平滑使用 `H=0.5*(1+tanh(phi/(2*interfaceThickness)))`，负侧权重 `1-H`，正侧权重 `H`。

- 快速几何子集为 `Domain → Point → Line → LevelSet`，只支持一个命名、闭合、位于矩形域内部的材料轮廓。点/线采用现代正整数 1-based ID；绑定使用 LevelSet 名称，不要求 dummy LS0 或固定 LS1。
- 支持直线和正向圆弧（`radian` 单位为 pi）；本批不接受反向/零编号线引用、样条、Phase/Boundary、固体或接触线。这样不会把尚未接通的固壁参数静默当成流体配置。
- `advectInterface true` 要求所选 LevelSet 声明 `moving true`；关闭该开关时材料界面保持固定。`reinitializeEvery 0` 不重初始化；正整数表示在相应步调用已有精确重距离能力，只能与界面输运一起启用。
- 继续使用固定均匀网格和闭域零法向速度。RZ 中 `x=r`，域从 `r=0` 开始；离轴闭合截面绕轴旋转对应环形材料区域，不应解读为轴心球形液滴基准。
- **本批不启用表面张力、SharpInterface 固壁直接力、蒸发、电力或 EHD 耦合。** 旧 SharpInterface 的默认 RZ 曲率、Weber 数常数和被注释的 no-slip 构造均未被冒充为已恢复功能。
- 未实现保守两相输运、体积修正、移动拓扑或 AMR 字段迁移。体积变化是诊断，不代表已证明守恒；压力/速度与材料更新采用显式分裂时间顺序，不作完整多物理收敛声明。

### 13.4 使用与检查

新案例：`cases/hydrodynamics_levelset_planar`、`cases/hydrodynamics_levelset_rz`。从项目根目录运行：

```powershell
cmake --build build --config Release --parallel 8
.\bin\QUANTUM.exe --case cases/hydrodynamics_levelset_planar --check-case
.\bin\QUANTUM.exe --case cases/hydrodynamics_levelset_planar
.\bin\QUANTUM.exe --case cases/hydrodynamics_levelset_rz
ctest --test-dir build -C Release -R "^hydro_" --output-on-failure
```

运行会写入所选案例的 `output/`；需要保留已有结果时，先复制 `system/` 到新案例目录。CTest 始终使用 `runs/history/checks/tests/hydro_cases/` 中的新副本。原 `cases/hydrodynamics` 的旧电喷雾相关输入尚未整体迁移，本批应使用这里的新案例入口。

### 13.5 本批实跑结果

- Release/生产 OpenMP 构建成功；`ctest --test-dir build -C Release -R "^hydro_" --output-on-failure` **33/33 通过**（约 37 秒），保留原 25 项，并补入多步、输出节拍、材料/加速度消费和水平集接通检查。没有开展新的解析解、网格收敛或 EHD 复合验证。
- planar/RZ 两个新材料界面案例均为 32x32、100 步、`dt=0.001`，至 `t=0.1`。均输出初始及 20/40/60/80/100 步共 6 帧，材料非均匀、水平集实际变化；第 100 步仍有 44 个 cut cells。

| 第 100 步指标 | planar2D | axisymmetricRZ |
|---|---:|---:|
| 动量相对残差 | `3.45e-16` | `3.52e-16` |
| 压力相对残差 | `1.39e-14` | `1.23e-13` |
| 投影后最大散度 | `2.22e-16` | `3.04e-15` |
| 最大速度 | `3.511e-2` | `3.265e-2` |
| `max(abs(phi(t)-phi(0)))` | `3.054e-3` | `2.867e-3` |
| 单元密度范围 | 约 `[1.000004,1.977486]` | 约 `[1.000004,1.977433]` |
| 轴径向速度 | 不适用 | `0` |

`interface_max_change` 是**水平集标量场相对初始场的最大差**，不是精确的几何界面位移；`negative_volume` 是按平滑相分数积分的面积（平面）/带 `2*pi*r` 度量的体积（RZ），不是精确切割体积。材料诊断对应输运后的材料场，流体残差对应该步用输运前材料完成的显式分裂求解。失败行的几何诊断可能是最后一次成功收集的状态，应先检查 `success/stage/reason`。

- 两个常材料案例额外完成 100 步；关闭写场时无 `.dat`。7 步输出检查仅产生第 0/3/6/7 步场，确认末步不会因输出间隔而丢失。
- `--check-case` 在新 planar 案例通过且未创建 output/cache；未知界面名与零厚度预检返回非零。几何初始化、输运、材料更新及输出异常可进入 CSV 失败记录；运行时若负相达到外域边界则明确停止，因为尚无接触线模型。
- 另在隔离输入上各跑 6 步：`advectInterface false` 时水平集差严格为 0；`reinitializeEvery 2` 能完成已有精确重距离调用，字段与残差有限。后者平滑负相体积从约 `0.11818` 变为 `0.11700`，说明**重距离仍有体积漂移，不能宣称体积守恒**；发布案例默认关闭重距离。
- 现有 ImplicitPIC Turner periodic-Y 20 步在独立副本中返回 0；未修改 PIC 物理或公共几何内核。这里仅记录本批 smoke 兼容性检查，不宣称完整 PIC 基准比较通过。

本轮正式几何案例日志/场位于 `runs/history/checks/tests/hydro_cases/hydro_levelset_planar/run_4ff29a8985ee/` 和 `hydro_levelset_rz/run_2c9a3ea614a5/`；手动分支检查与 PIC 副本位于 `runs/history/hydro/hydro_runtime_checks_a5d026b2/`。用户原先在两个三步 smoke 案例里的 output/cache 保留。

后续先沿此运行链路补表面张力、真实壁面和更广泛几何，再安排物理/耦合验证；是否开启这些项需单独确定范围，本批不继续扩展。

## 14. 表面张力与界面力接线（2026-09-09）

### 14.1 本批范围与接口

按用户“接表面张力与界面力”的要求，沿第 13 节固定均匀网格、单一域内闭合材料界面的链路继续接线，不扩展 AMR、电求解、蒸发、固壁或接触线。

- 公共 `SharpInterface` 新增受检的单界面几何重建与节点牵引扩散入口 `rebuildHydroInterface` / `spreadHydroTraction`，不依赖 Hydro 材料配置或电喷雾常数。旧 `initialize/update/IBM_surfaceForceSpreading` 保留原行为，未把旧电喷雾驱动重新打开。
- 原曲率实现默认包含 RZ 环向项，旧扩散缺少节点半径与单元半径的比值，而且邻接/stencil 存在未经检查的 `-1` 索引。新入口独立检查固定网格、节点连接、导数、曲率、扩散支持与有限值，不直接调用旧的未受检邻接构造。
- `HydrodynamicsModule` 负责界面物理模型：把配置的表面张力、附加法向/切向牵引转换为节点上的力面密度。Core 只提供几何与插值/扩散。
- Hydro 增加显式单元体力密度重载 `solve/solveRZ(mesh, dt, options, forceX, forceY)`，装配 `externalF = interfaceForce + rho * acceleration`。入口检查尺寸、有限值及输入与 `externalFx/Fy` 输出的别名冲突；不启用旧 `surfaceForces/evaporation/directForce` 标志，不读取 `number_We_inv/nonDim_g`。
- 每步顺序为：`phi^n → 材料/界面节点/曲率 → 节点牵引 → 单元体力 → 动量/压力投影 → 可选 phi^(n+1) 输运 → 派生量与材料刷新`。静态材料界面和移动材料界面使用相同受力入口；下一步从当前水平集重新构建，避免沿用旧节点。

### 14.2 配置、符号与度量

在原 `materialModel twoPhase` 配置内新增以下可选项（默认均为 0）：

```text
surfaceTension 0.01
interfaceNormalTraction 0
interfaceTangentialTraction 0
```

`surfaceTension` 是非负常量表面张力系数；另两项是可正可负的、空间均匀的附加界面牵引，不是体加速度，也不是电场或 Maxwell 应力的自动计算。三项必须有限，且即使显式填 0，也只能在 `twoPhase` 中声明。数值仍为代码单位：表面张力乘以曲率得到力面密度，扩散后得到每体积力（XY 按单位出平面厚度）。没有 SI/Weber 数自动换算，也未接入变表面张力或 Marangoni 模型。

法向约定 `n = grad(phi)/|grad(phi)|`，从负相指向正相；切向 `t=(-n_y,n_x)`。曲率及牵引为：

```text
XY: kappa = -div(n)
RZ: kappa = kappa_in_plane - n_r/r
traction = (surfaceTension*kappa + interfaceNormalTraction)*n
           + interfaceTangentialTraction*t
```

因此负相为圆内部时，XY 曲率为负，正表面张力指向圆内；正的附加法向牵引指向正相。RZ 的 `x/y` 分别是 `r/z`；当前 CLI 不接受接轴/接壁轮廓，RZ 示例是离轴环状材料区域，不是球形液滴。

扩散复用 Peskin 核并归一化每个节点的有效单元权重 `sum(w)=1`。XY 节点测度为 `ds`、单元测度为 `cellArea`；RZ 分别为 `2*pi*rNode*ds` 和 `2*pi*rCell*cellArea`。每个贡献为 `traction * nodeMeasure * w / cellMeasure`，保证该离散规则下节点输入积分与单元输出积分一致。这只是**力扩散积分一致性**，不是两相质量、能量或完整离散力平衡的证明。RZ 的径向积分分量也不能当成三维 Cartesian 全局净合力。

### 14.3 运行与诊断

新增 [`hydrodynamics_capillary_planar`](../cases/hydrodynamics_capillary_planar/README.md) 和 [`hydrodynamics_capillary_rz`](../cases/hydrodynamics_capillary_rz/README.md)：固定 32x32、零初速、`surfaceTension=0.01`、`dt=0.001`、20 步，默认 `advectInterface false`，用于先确认曲率—力—压力响应。旧两个 level-set 输运案例保持零界面力输入，不改已有物理设置。

```powershell
cmake --build build --config Release --parallel 8
.\bin\QUANTUM.exe --case cases/hydrodynamics_capillary_planar --check-case
.\bin\QUANTUM.exe --case cases/hydrodynamics_capillary_planar
.\bin\QUANTUM.exe --case cases/hydrodynamics_capillary_rz
ctest --test-dir build -C Release -R "^hydro_" --output-on-failure
```

正常运行会覆盖所选 case 中的同名输出；CTest 自动使用隔离副本。`--check-case` 仅做输入和支持范围预检、不构造节点，不能替代实际运行时的可解析拓扑/曲率检查。

CSV 保留原 21 列，追加 12 列：`surface_tension`、`interface_node_count`、`interface_surface_measure`、`max_abs_curvature`、`force_time`、`interface_force_x/y`、`spread_force_x/y`、`max_interface_force_density`、`interface_normal_traction`、`interface_tangential_traction`。力积分只计界面贡献，不含 `rho*acceleration`。

`force_time=(step-1)*dt` 表示**实际送入本步动量求解的旧时间层力**，而场文件中的流速/压力与输运后水平集属于输出时层；不将旧时层力混写到该 `.dat` 场。全零系数不构建力节点，新增力诊断为 0。几何、牵引、扩散失败分别记录 `interface-geometry`、`interface-traction`、`interface-force-spreading`；失败行的未完成力诊断可能为 0 或部分结果，必须先看 `success/stage/reason`。

本批仍是显式受力、固定 `dt`；未加入毛细时间步自动限制。更大表面张力、更小网格或更大材料比需要重新选择步长。静态圆滴 Laplace 压差幅值、寄生流衰减、网格收敛、形变/振荡频率、Marangoni/EHD 耦合仍留待后续复合验证，不以“跑完”替代这些验收。

### 14.4 构建与实跑记录

- Release 全构建成功，使用生产 MSVC/LLVM OpenMP 和 Eigen SparseLU；仍有已有 Geometry/Electricity 的 C4267 警告，外部 SuperLU 不可用。
- `ctest --test-dir build -C Release -R "^hydro_" --output-on-failure`：**44/44 通过**（22.33 秒）。原 33 项保留，新增 2 项库接口检查和 9 项案例检查：力输入的尺寸/非有限/别名/溢出拒绝、体力累加、XY/RZ 曲率及扩散积分、两个毛细力案例、移动受力、独立法向/切向牵引、显式零系数兼容及三种非法配置。
- `--check-case` 对新 planar 毛细力案例通过，源 case 未产生 output/cache。正式实跑与 CTest 都使用隔离副本，未覆盖用户已有 smoke 结果。

零初速、固定界面两个案例均完成 20 步至 `t=0.02`，有初始及第 10/20 步三帧；第 20 步记录如下：

| 指标 | planar2D | axisymmetricRZ |
|---|---:|---:|
| 界面节点数 | `32` | `32` |
| 界面长度/旋转面积 | `1.129315` | `3.547849` |
| 最大绝对曲率 | `5.630028` | `7.051045` |
| 最大界面体力密度 | `0.946465` | `1.175162` |
| 压力相对线性残差 | `2.03e-15` | `1.98e-15` |
| 投影后最大散度 | `6.98e-17` | `1.01e-15` |
| 最大速度 | `8.16e-4` | `1.25e-3` |
| `max(abs(phi(t)-phi(0)))` | `0` | `0` |

两种坐标的节点与单元积分力差均在 `1e-16` 量级以内。平面输出中，以距圆心 `<0.1` 和 `>0.35` 的顶点采样，平均压力分别约 `0.055527` 和 `-1.09e-5`，确认正表面张力的圆内增压方向；自动测试只检查正负方向及非零响应，不设 Laplace 幅值误差门限。上述静态界面仍有非零流速，**不宣称寄生流消除或静态力平衡已通过**。

另在两套隔离输入上启用 `advectInterface true`、平滑初速幅值 `0.05`，保留 `surfaceTension=0.01`，各运行 100 步至 `t=0.1`：

- XY：水平集场值最大变化 `3.05371e-3`，最终最大曲率 `5.66635`、体力密度最大值 `0.943382`，投影后散度 `2.22e-16`。
- RZ：水平集场值最大变化 `2.69394e-3`，最终最大曲率 `7.10364`、体力密度最大值 `1.183374`，投影后散度 `1.83e-14`，轴径向速度为 0。

力诊断随几何实际改变，不是重复使用初始曲率。这里仍只是固定网格上的接线/健全性检查，没有完成形变解析基准或严格体积守恒验证。

现有 ImplicitPIC Turner periodic-Y 20 步在新隔离副本返回 0；没有改变其输入物理或碰撞数据。这是 smoke 兼容性记录，不是全量 Plasma 历史 baseline 比较。手动毛细力、100 步移动受力与 PIC 日志/场保存在 `runs/history/hydro/hydro_force_runtime_1d8c5b50/`；自动回归副本在 `runs/history/checks/tests/hydro_cases/`。

### 14.5 后续复核时优先处理的遗留风险

本次接口审查还发现两项原求解主体的风险，未借本轮表面力接线悄悄改动离散：

1. 当前 RZ 径向黏性装配疑似缺少环向应力项。`HydrodynamicsRZ_impl.h` 中活动矩阵/RHS 与同文件保留的旧注释实现不同，后者显式含 `-2*mu*u_r/r^2` 贡献；后续应以 RZ 制造解/黏性算例核对完整算子，不能仅凭本批小散度认定 RZ 动量物理正确。
2. 当前 AB2 对流历史使用本步密度乘历史速度（如 `_qNorthPre = _rhoNorth * _vNorthPre * AreaNS`），不是独立的历史密度。移动变密度问题的时间层一致性、动量/质量守恒必须专项复核。

显式单元力入口也仍要求调用方保证力场按**同一网格的 cell index 顺序**排列；当前仅检查长度、数值和输出别名，不提供跨网格重排或 revision/AMR 自动迁移。本模块在同一个 mesh 上初始化力场，不跨网格传递。完整 FieldRegistry/AMR 生命周期契约仍未实施。

## 15. 静态压差、寄生流和收敛验证（2026-09-09）

本轮按用户要求只增加可复现验证工具、独立几何探针和报告，未修改生产几何/流动离散。详见 [静滴验证完整报告](validation/hydro/QUANTUM_Capillary_Verification_20260909.md)，其中固定测量口径、工程门限、原始结果路径及重分析命令均已落文档。

共完成十个正式流动实验：原材料与等材料高黏各 N32/64/128 的空间序列、N64 两个额外时间步、N32 一个黏性时间长跑和 N64 开启输运对照；另完成 16 项 XY/RZ 解析距离/字典距离几何测量。结果不是全绿：

- 原材料固定 dt=0.0025、T=1.28 的压差误差依次为 **0.043515%、0.064373%、0.292553%**。均满足 1% 幅值筛查，却随加密增大，压差收敛未通过。
- 寄生流 Umax 为 `2.12249e-3 → 5.13102e-4 → 3.59557e-4`，表观阶从 2.048 降到 0.513；细网格 Ca 满足 1e-3，但不能认定稳定二阶。N32 延长至 T=12.96 后 Umax 仍为 `2.11505e-3`，确认持续寄生流。
- 等材料高黏三档的压差误差为 2.48922%、1.79453%、1.35126%，共同 dt 的细网格仍超过 1%。N64 将 dt 从 0.002 减至 0.001、0.0005，误差变为 **6.59507% → 1.79453% → 0.301941%**；速度已呈平台，压差仍明显依赖 dt，支持压力分裂偏差解释。
- 几何探针中解析 phi 的曲率误差随网格下降，实际字典圆弧路径却从 0.536% 增至 2.057%。已定位到圆弧先折线化、再取到弦段距离的初始化瓶颈。积分力一致性约机器精度不代表局部力平衡。
- 高黏自由界面对照的分片线性面积相对漂移为 `5.23e-7`，只通过单网格小位移筛查，不宣称完整质量守恒或自由界面收敛。RZ 只做离轴环面几何测量，没有静态球滴验收。

Release 构建成功，既有 Hydro 回归 **44/44**、分析器自测 **24/24** 通过；正式实跑均正常完成。这些是运行/回归通过，不替代物理失败项。本轮未重跑 Plasma 基准，源案例已有输出保留。

下一批建议先修**公共圆弧距离/曲率初始化**，再明确**物理压力与投影势的恢复契约**，随后复测局部受力平衡、寄生流和自由界面；AMR 和电耦合扩展建议安排在这些基础误差得到控制之后。以上修复本轮未实施。

## 16. 公共圆弧距离、压力恢复修复与原样复测（2026-09-09）

第 15 节的前两项修复现已实现。完整接口、兼容性边界、测试和数值见[修复与复测报告](validation/hydro/QUANTUM_Capillary_Fix_Revalidation_20260909.md)。

- 公共初始化对有限圆弧使用解析距离与一致的解析射线符号；保留 innerPoint、样条采样、闭合和近零保护，修复反向圆弧采样漏乘 pi。R=0.18 的字典 phi 与解析 phi 最大差降至 `5.55e-17`，N32/64/128 曲率相对 L2 误差恢复为 `0.283657% → 0.0777031% → 0.0274988%`。
- 通用固定 dt 闭域 XY/RZ 使用增量投影：面预测减去历史压力梯度，本步私有投影势只用于速度校正，物理 p 再累加。共用面梯度模板；不把常黏度旋转修正套到任意两相物性。
- 旧电喷雾入口保留非增量及原东侧压力边界；禁止同实例投影模式切换，旧 Hydro 缓存保留 legacy-only。原 `y≈0` 电喷雾强制行限制为旧入口，不再污染通用平移域。
- 十组实验的 **40 个输入文件 SHA256 与前轮完全相同**，分析器/工程门限不变，全部在新目录完成。原材料压差误差变为 **0.109114%、0.0237513%、0.0119114%**，恢复单调下降；等材料为 **0.110144%、0.0238316%、0.0118893%**，三档均通过 1% 幅值门。
- 等材料 N64 的 dt=0.002/0.001/0.0005 压差误差由 **6.59507%/1.79453%/0.301941%** 降至 **0.0236909%/0.0238316%/0.0238316%**。时间解差进入打印精度底的项不计算表观阶；不宣称动态二阶压力。
- 原材料 N128 Umax 从 `3.59557e-4` 降至 `2.74746e-4`，但 Umax 空间阶仍为 2.226、0.704；N32 长时间仍有 `2.08623e-3` 的持续流动。**寄生流未消除、稳定二阶未证明**。自由界面近似面积相对漂移 `4.87e-7`，仍仅为单网格小位移筛查。

Release 全构建、Hydro CTest **50/50**、Python 分析器自测 **24/24** 通过。新压力回归包含变量 rho/mu、已知/未知压力、gauge、跨内部 y=0、别名输入、坏历史、模式切换与旧缓存边界。增量投影仍有启动/时间分裂误差，旧几何缓存不会自动获得新解析初始化；未实施 Hydro AMR 或 RZ 球滴支持。

下一批优先检查界面力到单元/面和压力梯度的**局部离散平衡**，以持续寄生流及稳定收敛阶为目标；不再把本轮已修的圆弧距离或旧压力滤波作为尚未消除寄生流的唯一原因。

## 17. 通用电流体 App 首批接入（2026-09-09）

根据用户“可以的话开始调电流体 app”的后续安排，在保留第 16 节流体未验收项的前提下开始最小耦合。新增 `Electrohydrodynamics`（别名 `EHD`）App，详细模型、接口、输出时间层和实跑结果见[首批理想介电耦合报告](validation/ehd/QUANTUM_EHD_App_First_Coupling_20260909.md)。

- 单个 Hydro 驱动拥有 mesh/geometry/flow/clock，电学通过 `HydroCoupling` 借用同一状态并提供受检单元体力，不用多个独立 Module.run 假装耦合。
- 首批仅 planar2D、固定均匀网格、单命名闭合界面、无自由电荷的理想介电两相。电场通过公共变系数 FV Poisson 真正求解，介电率随公共水平集刷新；电体力与现有毛细力相加，推进速度/压力和可选的界面运动。
- `run.dict` 只选择 App；几何、流体、电材料与电压各归对应字典。严格 `--check-case` 无写入；失败非零退出，电场载荷和步末流体场分别标注时间。
- 旧 `ElectrosprayApp` / `Electricity` 的初始化演示与 RZ 制造解没有当作物理求解复用；保留旧行为并明确提示。新 App 不含电导率、自由表面电荷、发射或粒子，不宣称旧导电电喷雾已恢复。

本批不消除第 16 节的寄生流和局部离散平衡风险，也不解除 AMR、RZ 及重启限制。后续需先给理想介电分支补电场/形变与网格—厚度收敛，再分批实现表面电荷守恒、导电耦合和喷雾边界/模型。

本批接线现已实跑完成：Release 全构建、新增 EHD **21/21**、既有 Hydro **50/50**、毛细分析器 **24/24** 通过。发布案例固定界面 20 步至 T=0.02，另一个隔离副本启用界面运动、将电压差改为 0.5 后运行 100 步至 T=0.1，均返回 0。后者 LS 最大场值变化 `1.62486e-4`，电场/电力随几何刷新，末步最大流速 `2.68736e-3`、投影后散度 `1.38778e-17`。这证明最小理想介电耦合链路可运行，不证明导电电喷雾、锐界面电应力或形变收敛已验收。

### 17.1 SuperLU 后端更正与复测

已找到本机 `C:/vcpkg-master/installed/x64-windows` 中的 SuperLU 7.0.0；此前是跨 vcpkg 前缀导致依赖发现失败，不是环境未安装。CMake 现统一传播并链接 SuperLU/OpenBLAS/METIS，部署对应运行库，缺 SuperLU 则配置失败。Hydro 动量/压力不再硬编码 SparseLU；新 EHD 电势强制 SuperLU，失败不回退。实际成功日志分别报告三个后端，CSV 和原数值门限保持。

Release 全构建、EHD **21/21**、Hydro **50/50**、分析器 **24/24** 再次通过；固定 20 步、运动 100 步在 10 个字典哈希完全相同的新副本中完成，全部活动后端日志为 SuperLU。旧 Plasma/MMS 调用者的默认后端策略保留，未全局移除其他算法；隔离的既有 polarization MMS 兼容性 smoke 返回 0。更正依据、构建命令与验证边界见[首批报告第 7 节](validation/ehd/QUANTUM_EHD_App_First_Coupling_20260909.md#7-superlu-启用与原样复测2026-09-09)。

## 18. 理想介电液滴首轮定量验证（2026-09-09）

按用户“先定量验证，再扩展导电电喷雾、RZ、AMR”的顺序，本轮只补验证工具/独立参考，不改生产电学、力或输运离散。完成 13 组固定圆滴电场/电力测量、8 组短时运动/零电压及 dt 对照，详细方法、失败项、门限和复现命令见[定量报告](validation/ehd/QUANTUM_EHD_Dielectric_Validation_20260909.md)。

- 同方盒边界独立圆谐波参考通过外边界、圆上 V/Dn/Et 跳跃及模式自收敛检查，不用无限域或三维球滴公式代替实际二维边界。
- 原 smoke 的 F2 误差约 31.2%，首轮 N128/t0.01 仍为 11.2%（5% 门未过）；固定 t 加密趋于平台。联合细化至 N512/t0.0025，F0/F2 误差降至 0.768%/2.971%，通过幅值筛查，保留粗网格失败产物。
- 两档运动网格均有近似 V² 的沿场形变；N64 时间步减半后 D 差约 0.009%。完成面积、质心和零电压对照诊断，不宣称稳态形变或空间收敛已验收，零电压寄生流仍存在。
- Release 全构建、EHD 23/23、Hydro 50/50 通过；新验证探针及所有活动电/流求解均要求 SuperLU。源案例参数未变，未扩展导电、RZ、AMR。

下一步优先：验证并改善电力/毛细力到面速度和压力梯度的局部平衡，明确 t/h 的可承受误差，补第三档运动网格及较长时间的自由边界形变参考；不能把本轮单个细化点通过工程幅值门当作整个 EHD 模型完成物理验收。

## 19. 界面力与压力的相容面离散（2026-09-09）

按用户后续要求，本批只修正界面力局部离散并保留原验证口径复测。方法、输入哈希、完整指标及未通过项统一记录于[界面平衡报告](validation/ehd/QUANTUM_EHD_Interface_Balance_20260909.md)，第 18 节结果继续作为旧基线保留。

- 新增受检 XY/RZ 面力端口，与历史压力使用相同 dual volume；保留旧 cell 入口，不扩大为 AMR 支持。
- 毛细力改为 `sigma*kappa_face*G(Hcap)`，Hcap 是公共 φ 的紧支撑三角裁剪几何面积分数，曲率仍取实际重建节点。未拟合圆或调整 σ，未把材料 tanh 厚度改作固定毛细厚度。
- 电力由实际 Poisson 位移通量构造面载荷，平层恒 Dn/Et 时满足离散 Maxwell 链式关系；cell E 从 D/ε 重建。电力与 σ 均只加一次，真实非梯度形变载荷保留。
- 原 CSV 列/参考解/掩膜/门限保留；旧 node/Peskin 力列现在明确作为参考诊断，另记实际面载荷和 `F(t^n)-G p^(n+1)`，后者不包含完整惯性/黏性平衡判据。
- Release 构建、Hydro 52/52、EHD 24/24 及原毛细分析器 24/24 通过。13 组原静态电输入、8 组原运动输入已复测；所有活动电/流求解均使用 SuperLU。10 组原毛细输入的完整结果见报告 §4.3。

结论边界：平层/恒曲率制造测试的机器精度平衡，不等于真实曲率、零初压、任意 t 的液滴都已静止。原固定 t 的电 F2 偏差平台仍在，N512/t.0025 新 F2 误差约 3.05%（原 5% 门通过）；高黏细网格虽远场压差通过，界面局部压力恢复仍较慢。下一步应先处理/验证这一压力启动与时间分裂问题，再补薄层解析度、运动网格/时间及较长时间形变对照，不急于扩展导电电喷雾、RZ 电算子或 AMR。

## 20. 高黏零初压恢复与耦合校正（2026-09-09）

已继续完成第 19 节的压力启动优先项，详见[压力恢复与原样复测报告](validation/ehd/QUANTUM_EHD_Pressure_Recovery_20260909.md)。旧报告和输出保留为基线。

- 固定均匀闭域在原 predictor/Poisson 之后，使用同一已装配动量矩阵、压力梯度和 FV 通量作速度—压力稀疏校正；不改初压、dt/T、物性、力或历史时间层。旧非增量与非均匀路径不启用该校正。
- momentum、pressure、coupled 分别缓存实际 SuperLU LU，仅完整 CSC 逐字节一致时复用；矩阵变化或求解失败明确失效，不回退 SparseLU。
- 保留原 38 列及残差口径，追加 8 列。新动量量是已装配增量形式残差，不能冒充独立总载荷重建；原 `F−Gp` 也不含完整惯性/黏性平衡。
- Release 全构建及 Hydro 54/54、EHD 24/24 通过。原 18 组运动/毛细、80 个字典、10,384 步原样复测成功；原 13 组静态电输出逐字节一致。参考、掩膜、门限和物理参数不变。
- 高黏 N128/T=.256 的末态 `max(F−Gp)` 从 `.195338` 降至 `.000576599`，Umax 仍为 `2.61840e−7`。N64 三档 dt 的压力差已低于原打印精度，但不宣称时间二阶或寄生流消失。

下一步应集中在充分解析的材料薄层、真实曲率误差、第三档运动网格和较长时间形变参考。原固定厚度的电力幅值偏差、部分稳态指标及稳定空间阶仍未解决；继续暂缓导电电喷雾、RZ 电算子、AMR。

## 21. 材料薄层解析度与真实曲率的独立验证（2026-09-09）

按用户最新要求，本批只补验证，不改生产离散。完成 11 次固定厚度电场、四档网格×实际/恒曲率的 2048 步 NS 对照，以及 12 组子网格平移几何检查，详见[薄层与曲率报告](validation/ehd/QUANTUM_EHD_ThinLayer_Curvature_20260909.md)。原参考、掩膜、门限和 SuperLU 要求不变，旧结果保留。

- 固定 t=.01，N512 的 F2 锐界面误差仍为 10.6647%，256→512 增量仅为参考值的 .0950%，局部力限制后 L2 差 .8885%；支持厚度相关偏差占主导，但不把有限网格差当成严格误差界。
- t=.005/.0025 的 N512 F2 误差为 5.7103%/3.0496%，局部力跨网格差却为 3.30%/10.42%；后者通过原 5% 积分门仍不等于薄层充分解析。
- 同 Hcap、零初始 p/u/history 的恒曲率测试支路全程 Umax≤1.12e−17；真实曲率 N32/64/128 与旧生产末态 Umax 完全一致，新增 N256 为 5.95e−8。当前圆滴的剩余寄生流应优先归因于曲率空间变化，而非再修改已修正的圆弧距离或恒曲率压力平衡。
- κ RMS 的中心组表观阶为 1.868/1.499/2.066；三个对齐在最后一级分别为 2.066/1.972/1.930。均值的小偏差可能来自抵消，不能代替局部误差，也不宣称任意界面稳定二阶。
- 两个新增测量探针构建成功，原 Hydro 54/54、EHD 24/24 回归通过。本批未重新重放全部旧运动组，也未将新单元压力对照替换旧固定掩膜压差指标。

下一小批先评估通用曲率重建/面采样候选，并用非圆光滑界面防止圆特化；电学继续分开固定厚度解析度与锐界面模型偏差。两者按新报告协议复测后，再补运动网格和长时间形变；导电、RZ 电算子、AMR 仍暂缓。

## 22. 持续推进目标与代码结构约束（2026-09-09）

用户已要求在曲率问题处理后持续推进，直到流体—电流体—电喷雾在约定支持范围内可靠正确运行，并明确要求代码卫生和合理结构。第 21 节的“暂缓”是实施顺序，不再需要逐批等待是否继续的确认；仍不得跳过物理模型选择、验收或失败可见性。

按以下关口推进，未通过的能力保持明确受限，不以 smoke、制造解或某个算例专用修补替代通用物理验收：

1. **公共几何与 Hydro**：完成真实曲率误差分解和通用候选，保持现有节点/拓扑、力压离散及基线；验证圆/平移/非圆几何、静滴压差与寄生流、运动界面及原时间/空间口径。候选失败必须保留记录，不能只展示选中的方案。
2. **理想介电 EHD**：区分真实有限厚度模型与理想锐界面目标，验证界面电势/位移通量跳跃、电载荷及与压力的相容性；再验证短时和较长时间形变、空间/时间敏感度和守恒。保持独立参考，不用三维球滴公式替代二维有限方盒。
3. **导电与电荷输运**：明确电荷守恒、体/表面电荷及导电电流的模型与时间层，建立独立松弛/输运/界面电流基准；在没有守恒验证前不宣称导电电喷雾恢复。
4. **RZ 与电喷雾**：补齐轴线、球面曲率、电算子及应力度量，再处理喷嘴/入口、接触界面和流出边界。实际电喷雾的材料、供液/电边界和电荷注入/发射模型必须有明确依据；缺失且影响结果时向用户询问，不从旧算例硬编码推断。
5. **AMR 与综合可靠性**：实现并验证粗细网格间通量、物性、界面及电荷转移，确认重网格前后守恒和解的收敛；补完整工作流的失败行为、输出时间层及重启合同。最终以支持能力矩阵、可复现基准、原输入回归和收敛证据审计完成状态。

代码结构约束：通用几何/网格与离散基础放 Core；Hydro、电学、电荷及界面物理各归所属模块；App 只编排共享状态与时间推进；材料/几何/电极/供液等算例参数留在字典。测试解析解、圆/椭圆特殊参数和候选筛选不进入生产算法。避免并行保留重复生产求解路径、无用配置开关或悄悄回退；每批检查接口职责、资源生命周期、失败契约和遗留调用者影响。SuperLU、原参考/掩膜/门限及既有输入/结果归档继续保持。

## 23. 真实曲率重建接入与原样复测（2026-09-10 审计）

已完成[真实曲率修正报告](validation/hydro/QUANTUM_Curvature_Reconstruction_20260909.md)。生产只新增公共局部标量重建头，并在受检平面界面路径选择张量积双四次求导；保留旧投影、闭环、Hcap、面延拓、压力和电学。首轮退化的总三次/总四次候选保留测量，没有按圆滴半径或网格切换方法。

- 三批各 24 组圆/非多项式椭圆、2,957 节点几何逐项一致，选中方法在全部网格/移位的曲率 RMS 改善；原 12 组平移曲率及法向误差均降低。
- 原四档公共圆弧 actual/oracle 共 2,048 步复测，真实曲率 N32/64/128/256 末态寄生速度降低 2.82/11.33/27.32/55.44 倍；恒曲率舍入量级平衡不变。没有解析初压或恒曲率生产替换。
- 原 18 组、80 输入、10,384 步均完成，三个活动流体求解和电学使用 SuperLU；保留原压力/形变/时间口径。两组 N128 平均压差误差增大、全域压力 L2 略增，以及未过的稳态指标均记录，不能称全部指标都改善。
- RZ 两组共 120 步、13 份输出逐字节一致；此为兼容性证据，不是新增 RZ 曲率或电算子精度验收。曲率阶段 79/79 CTest、Hydro Python 113/113、EHD Python 82/82 通过。

本关口支持保留平面曲率修正，不意味着整体流体/电流体已经无剩余误差。按第 22 节持续推进，下一批已启动[锐界面电学重建合同](validation/ehd/QUANTUM_EHD_Sharp_Interface_Trace_20260910.md)：先验证一般斜界面 V/Dn/Et 与 Maxwell 牵引，再接公共守恒电算子、原圆滴和运动验证。实验候选先在 tests 验证，未把仅无厚度的一维谐均值替换当成二维锐界面模型完成。

## 24. 锐界面电学局部 trace 候选（2026-09-10）

已完成测试层两相共同 VΓ/Et/Dn 重建和 **972/972** 个多角度、亚格位置、物性比、纯法/切向及混合场的解析 patch；20 项标量不变性和 32 项非法输入自测通过，加入后原 Hydro/EHD 合计 **81/81 CTest** 通过。详见[局部合同、原始误差与限制](validation/ehd/QUANTUM_EHD_Sharp_Interface_Trace_20260910.md#5-局部候选的实跑结果)。

局部候选与独立 Maxwell 应力/分相面通量一致，距离谐均值基线在 920/3888 个面样本上有超出局部舍入界的偏差。此批仅新增 tests 与 CMake 目标，250 个生产 C++ 文件和生产 exe 未改变，没有把局部 QR 宣称为全局 SuperLU 电场求解，也没有切换生产模型。

下一具体实施项是“保留 TPFA 骨架 + 共同 trace 的切向通量修正”，先复用同局部矩阵，再核验共享内面反号、常数零模、锚定矩阵秩/棋盘格及完整离散残差；随后接公共几何、圆滴独立参考与 EHD 原样运动复测。全局稳定性、曲面与运动收敛仍未证明，原 diffuse 结果作为模型对照保留，目标未完成。

## 25. 公共守恒通量端口与全局平面验证（2026-09-10）

已完成[全局通量报告](validation/ehd/QUANTUM_EHD_Global_Interface_Flux_20260910.md)。界面专用候选仍在 tests；公共 Poisson 只新增通用积分面通量端口，复用同一 SuperLU 求解主体，严格校验双向守恒、非法输入和真实完整残差，新端口禁缓存冻结及替代求解后端。App 与生产 diffuse 电模型未切换。

- 原 972 组局部输入的 3,888 个修正面通量通过；旧默认 trace 输出逐字节不变。
- 216 组全局平面问题、基线/候选共 432 次实际 SuperLU 求解通过；候选最大电势误差 `6.11e-16`、物理面通量误差 `3.61e-15`。72 个 N8 锚定矩阵通过奇异值/棋盘格检查；这不是曲界面收敛或一般 SPD 证明。
- 矩形 8×6 回归发现公共 `BackgroundMesh` 竖边长度误用 hx，已源头修成 hy，未放宽校验。边端点/面矩/Gauss、仿射电场及 99 次失败无写回检查通过。
- 整项目 Release 构建通过，Hydro/EHD 84/84 CTest 通过。独立 before/after 原样复制的四个旧 Poisson 算例 15 份输出一致；PIC 20 步前后完成，23/25 份输出一致，仅运行清单/性能计时文件不同。线程限额导致的前后相同 PIC 失败也保留，没有把未存在的 EB smoke 算例计入通过。

紧接着先处理公共缓存读回/旧 parent 恢复仍可能保留矩形错边长的问题，正常缓存与既有结果必须保留；然后推进真实 level-set 曲界面电场及 Maxwell—Hcap—压力的相容性验证。原 10,384 步长程历史组本批未再重放；导电、电喷雾、RZ 电学和 AMR 仍未验收，持续目标未完成。

## 26. 公共缓存、静态加密初始化与曲面电场测量（2026-09-10）

详见[曲面通量与缓存报告](validation/ehd/QUANTUM_EHD_Curved_Flux_Cache_20260910.md)。保持生产模型和旧物理口径不变，新增 36 组二维圆柱独立解析测量，108 次实际 SuperLU 求解；同时关闭公共缓存的错误边长恢复和 AMR 存储/旧索引初始化缺口。

- 公共 QMC2 仍为 version 2，独立暂存、校验后再提交。旧矩形错误边长及必要拓扑失配明确拒绝，原网格与输入文件不变，不静默修补历史几何。
- 实跑先发现 AMR 新实体水平集数组被赋为空；修复细化与粗化回收槽。进一步的阶段标记暴露静态加密缺失旧索引准备，造成旧插值历史可能访问 `parent=-1`；已抽取 `prepareTransferIndices()`，静态加密与动态 adapt 共用，避免复制初始化逻辑。
- 同样例复测：9 组正常缓存往返、27 组读取拒绝、22 组写前保护通过；包括两层父子链、粗化和再次细化，保留位级快照及完整 QMC2 文件字节比较。所有初始失败/停滞日志仍保留。尚不保证任意缓存的分类字段已经可直接供 FD 使用，也不宣称完整物理重启或 AMR 守恒转移完成。
- 曲面参考自身 2,287 项独立检查通过，原 216 组平面输出逐字节不变。圆柱 affine 候选细网格明显改善，但 ε=10 的 N16/N32 电势误差仍可能比 TPFA 大 1.17–2.78 倍；精确点 Et 也未消除曲面及边界误差，不能按小线性残差宣布精度通过。
- 最终整项目 Release 构建及 Core/Hydro/EHD 86/86 回归通过；原局部 trace/面通量、全局平面及本批曲面重放均逐字节一致。旧四个 Poisson 算例的 15 份输出不变；PIC 20 步前后及各一次重复运行成功，但其近零 `dQrel_stage` 有 `2.09e-14` 的首轮跨 exe 差异，旧 exe 自身重跑也变化 `2.10e-14`，如实保留而不计为逐字节通过。未重跑原 10,384 步全组。

下一具体实现为测试层曲率感知分相 P2 与完整 `Q(Pfit)−B_TP(Pfit)` 面修正，先保持原 36 组和退化记录，再补非圆独立参考及真实 level-set 几何误差。之后仍须接 Maxwell—Hcap—压力平衡并重放旧液滴/运动口径，满足门限后才能切换生产锐界面模型；导电、电荷输运、电喷雾、RZ 与 AMR 的后续关口不变。

## 27. 曲率感知电学 P2 与单侧 Maxwell 载荷（2026-09-10）

已完成[二次电学重构报告](validation/ehd/QUANTUM_EHD_Quadratic_Trace_20260910.md)。候选、圆解析解和诊断各自隔离在 tests；公共求解器仍用同一通用面通量端口及实际 SuperLU，没有切换生产模型或在 App 中添加算例分支。

- 局部分相二次重构的 168 组独立多项式合同、67 个非法输入反例通过；48 组局部圆测量明确没有物理通过门，保留两个粗网格面通量退化和部分电势退化。
- 原 216 组平面、36 组圆问题及旧局部 trace/面通量输出原样重放、逐字节相同。P2 平面组通过；圆问题新增第四种算法，精确样本和实际求解误差分列，不调整原参考、区域或残差门。
- N128 新 P2 界面求解通量密度最大误差：ε比2 为 `1.20577e-4`，ε比10 为 `6.54538e-3`；比原 affine 最差 shift 汇总约改善 19.59/9.79 倍。但 ε比10 粗网格仍有不如 TPFA 的情况，bulk 和外边界误差并未全面改善。
- 同组固定 128 角度、完整 5×5 窗口测单侧场/Maxwell 应力。固定同一 global P2 求解电势，P2 后处理的最细应力 RMS 比 P1 降低约 13.65/12.80 倍；但固定 P2 后处理时，ε比2 原 affine 全局解反而比新 P2 的应力 RMS 更小。按构造为零的切向跳跃不计独立物理通过，不能仅以电势改善宣称界面载荷验收。
- 整项目 Release 构建及 Core/Hydro/EHD **88/88** 回归通过；252 个生产 C++/头文件及主程序与上一批完全相同，最终输入与 exe 哈希稳定。未重跑旧 10,384 步组或 Poisson/PIC 兼容性组，不重复计入本批成绩。

下一批先做公共水平集几何的显式测试适配：将几何与独立参考拆开，保留解析几何、只替换法向/曲率、完整 LS 几何三层对照；补任意线段分相、交点/相别/长度/位置/法向/曲率误差，避免用旧 EB 聚合量或 Hydro 投影节点冒充同一精确界面。之后仍按第 22 节补非圆、Maxwell—毛细力—压力与运动验收；导电、电喷雾、RZ 电学和 AMR 不跳关。

## 28. 公共水平集接入与 PL 几何误差归因（2026-09-10）

已完成[水平集几何报告](validation/ehd/QUANTUM_EHD_LevelSet_Geometry_20260910.md)。只在 tests 增加独立几何快照、公共圆弧初始化 fixture 和合同回归；生产 252 个 C++/头文件及主程序不变，App 不新增诊断分支。

- 原 36 组参数扩为八算法，每次 288 个实际 SuperLU 解，两次全组 JSON 完全一致；旧平面、旧四算法圆输出也逐字节一致。边界、解析真值、3h 掩膜、门限不变。
- 36 组公共圆弧初始化与解析节点距离全为零差，未触发零值保护；只换同位置 Q4 法向/曲率时，电场与原解析几何非常接近。
- 完整 PL 位置/相别/长度配 Q4 导数明显退化。N128 的 PL 位置 RMS 为 `3.37241e-5`；重建点曲率 RMS `1.07341e-3`，大于真界面点的 `1.74162e-4`。face/link 的切割并集份额 RMS 为 `.00536721/.00898891`，末级约一阶；保留粗细网格反增与相别/切割计数差异，不以全域 RMS 稀释局部误差。
- 即便相对原 affine 的求解界面通量 24/24 改善，同一解析几何 P2 后处理的应力 RMS 仍 24/24 更差；不能将通量改善等同电载荷通过。当前 PL+Q4 是诊断混合几何，不进入生产。
- 整项目 Release 构建及 Core/Hydro/EHD **89/89** 通过。新增独立几何合同 9,162 断言、56 个拒绝输入通过；保留首次 LNK1104 与无进程可见后的成功重链日志，不推测瞬时文件访问原因。未重新重放旧 10,384 步组。

下一小批先拆出公共原始标量导数合同，再测试共享节点双三次 Hermite 连续标量，让值、交点、法向和曲率最终来自同一插值。先证明 C1 和多项式再现，明确其非 C2 与根/切触/最近点限制；不得以单次 Newton 或 PL 回退冒充完整高阶几何。之后补非圆电学、完整 LS 载荷与原液滴/运动验收，再推进第 22 节的后续能力。

## 29. 公共原始导数与连续标量候选（2026-09-10）

已完成[公共标量导数与 Hermite 报告](validation/ehd/QUANTUM_EHD_Hermite_Scalar_20260910.md)。生产只拆分公共重构的原始导数与法向/曲率计算，常量/驻点可求导，原几何退化检查不变；没有复制 QR、增加 App 分支或切换物理模型。

- 独立标量合同 224 组多项式、376 项拒绝通过；测试层 Hermite 的 64 组 Q3、410 个双侧共享边、25 项拒绝通过，明确 C1 非 C2。
- 24 组真实圆/旋转椭圆点测量显示标量细网格 RMS 约 `3e-9`，但 Hermite 曲率 RMS 24/24 大于同点 Q4，N128 约大 3.5–4.3 倍，细网格法向也退化。保留无物理通过门的测量，不按连续性优势宣称几何/载荷验收。
- 原 24 组曲率和 12 组平移在旧/新 exe 下逐字节相同；旧 PL 合同、平面、P2 圆、八算法 LS 圆均与上批逐字节相同。SuperLU 和旧参考/掩膜/门限不变。
- 整项目 Release 与 Core/Hydro/EHD **91/91** 通过；保留首次 LNK1104 与无同名进程、目标不存在的诊断，重链成功但瞬时文件访问原因未确定。未重新重跑原长程 10,384 步组。

下一步实现同一连续标量上的线段全部根/相区与重建点几何，再原样评估全局电学与完整 LS 载荷；不能用真点标量误差代替实际零点/分相误差。根据完整对照决定双三次是否足够或需更高阶/C2，随后依第 22 节推进，不跳过原液滴/运动与守恒验收。

## 30. Hermite 几何接入、实测改善与未通过合同（2026-09-10）

见[Hermite 几何与载荷报告](validation/ehd/QUANTUM_EHD_Hermite_Geometry_20260910.md)。生产代码和主程序不变；新增标量限制、根隔离、几何组合均在 tests，圆 fixture 与共享误差统计分离。这里完成了电学测量，但**几何合同尚未完成**，不切换生产方法。

- 原 36 组扩为 12 算法，两次各 432 个实际 SuperLU 解，完整 JSON 逐字节相同。旧平面、P2 圆、八算法 LS 圆及标量/PL 合同原样重放一致，原参考、门限与掩膜不改。
- N128 的实际界面位置 RMS 降到 `2.97290e-9`，cut-union face/link 份额 RMS 为 `5.05541e-7/7.83248e-7`；较 PL 明显改善。但同点曲率仍差于 Q4，细网格法向也未改善，部分通量及 Maxwell 指标保留退化。
- 修复根工具对两个/六个充分分开的简单根的过保守误拒。其后 35 个根测试、576 个 Q3 限制及 72 个局部投影完成；116 个分相检查为 108 成功、1 个独立证明的表示限制、7 个仍失败，JSON `status:fail`、exit 2，未扩大例外名单。
- 已证实一个合法 Hermite 端点严格为零，转换后 double 幂系数的精确和却为 `-3/2^57`。另有条带内部相区冲突，不能都归因成端点或几何欠解析。下一步直接修共享节点权威标量与保留误差的限制/全根表示，不继续堆物理层或按算例改系数。
- 整项目 Release 构建成功；Core/Hydro/EHD **91/92**，旧 91 项通过，新几何测试失败并保留。272 个源码/CMake 和 11 个 exe 在最终复测前后哈希稳定，252 个生产文件及主 exe 与前批相同。历史长程组和 Poisson/PIC 本批未重跑。

先关闭上述分相失败，再补非圆电学参考及完整载荷归因，随后回到 Maxwell—毛细力—压力和原液滴/运动验收。导电、电荷输运、电喷雾、RZ 和 AMR 后续关口不变；持续目标仍未完成。

## 31. 权威 Hermite 限制与原分相失败关闭（2026-09-10）

见[精确限制与分相修复报告](validation/ehd/QUANTUM_EHD_Authoritative_Hermite_20260910.md)。以原共享节点 DOF 为唯一权威标量，限制保留浮点展开尾项；同一根算法返回完整包络，相别取自无根间隙。统一点/线段坐标转换，仅合并同相格线细分，不删除异相短层。

- 原 116 项分相与新增 6 项近格线检查全通过，移除上一批的表示冲突例外；根位置、相别、长度和拒绝门限不改。
- 新独立整数回归核对 784 项算术、23 组根、27 组限制、384 个单元界、6,144 个 Bernstein 控制量及 2,688 个点值界；335,417 断言通过。严格区间快筛只在已证明单号且满足闭单元范围时生效，近零仍精确回退。
- 整项目 Release 与 Core/Hydro/EHD **93/93** 通过。原 36 组/12 算法两遍各 432 次 SuperLU 解，完整 JSON 一致；旧八方法算子结果、旧平面/P2/LS、标量/PL 回归保持一致。没有 SparseLU 回退。
- 修复后的圆问题主要量只发生舍入量级变化，不宣称真实曲率或 Maxwell 误差随分相合同变绿而解决。旧退化、无物理精度门的测量及有限拓扑/局部投影支持范围保留。
- 本批仅测试候选/CMake 修改，生产 252 个 C++/头文件与主程序不变；274 个输入和 12 个 exe 在最终复测前后哈希稳定。保留两次 LNK1104、主动停止的无快筛性能尝试及其未完成状态。历史长程与 Poisson/PIC 组未计入本批实跑。

下一具体批次为非圆独立电学参考及“几何 × 电势来源”载荷归因；参考/fixture/通用候选分离，非 SDF 椭圆隐式函数不得替代物理距离掩膜。随后按第 22 节继续局部力平衡、运动液滴、导电/电荷输运、电喷雾、RZ 与 AMR，不提前切换生产模型。

## 32. 非圆独立参考、共享诊断与载荷对照（2026-09-10）

已完成[旋转椭圆参考及测量报告](validation/ehd/QUANTUM_EHD_Ellipse_Reference_20260910.md)。参考、形状 fixture、共享 trace/geometry 诊断内核分层，候选不接收椭圆参数；生产 252 个 C++/头文件与主程序不变，仍实际 SuperLU，不切换生产电模型。

- 独立参考 16,927 项断言、7/7 段通过，两次 JSON 一致；检查界面条件、势差分场、独立分相积分、Gauss 与最近距离。中轴非唯一法向、一般不可分离近切触及不可表示参数仍明确拒绝，不伪装支持。
- 36 组椭圆各 108 次实际 SuperLU、27 次矩阵检查，三轮完整输出逐字节相同。非 SDF 节点函数与真实欧氏 3h 掩膜严格分开，旧圆/平面/LS 及标量几何合同重放字节不变。
- N128 重建点曲率 RMS 最大 `1.42139e-3`；位置末级约四阶不等于切割面份额同阶，face 份额 Linf 在 N64→128 几乎停滞。精确点值的 P2 Maxwell 误差仍明显，几何路径差异不足以解释全部载荷误差。
- 24 组非均匀介质中，Hermite-P2 相对解析 P2 的界面通量 11 组变差；相对 TPFA 仍有 4 组退化，最坏 2.08949 倍。保持物理精度门 null，不按残差小宣称 EHD 通过。
- 全量 CTest 100/101，唯一失败是既有 PIC EB smoke 缺算例目录；原 Core/Hydro/EHD 范围最终 95/95。保留 CLI 重链 LNK1104 和无可见同名进程/目标不存在的诊断，重试构建成功。十项 CLI 拒绝通过，历史 10,384 步长程组未重跑。

下一小批直接复用 Hydro staggered 面力入口：固定既有 Hcap，测试 `tau_face*G(H)` 的恒牵引压力平衡和真实变牵引首步响应；不把变牵引写成 `G(tau*H)`，不把椭圆强行当静滴。以独立原动量方程及 `F→F+Gq` 动态不变性验接线，再接 P2 载荷、毛细力与原运动/厚度口径。导电、电荷输运、电喷雾、RZ 电学和守恒 AMR 后续关口不变，持续目标未完成。

## 33. Maxwell 面牵引与真实 Hydro 接线合同（2026-09-10）

已完成[面牵引—Hydro 报告](validation/ehd/QUANTUM_EHD_Maxwell_Hydro_20260910.md)。新增形状无关的测试面力 adapter，复用生产压力面度量；旧独立动量差分抽成共享测试头，旧回归前后输出逐字节不变。生产模块与主程序不变，不新写电学或 NS 求解器。

- 两遍各 100 次真实 Hydro 单步，momentum/pressure/coupled 三处实际 SuperLU；5/5 section、249,513 本测试断言通过，两遍 JSON 一致。电场为解析参考，本批不计电学 FV 求解。
- 16 组零初始 p 的恒牵引，最终最大速度 `4.21482e-18`、压力误差 `1.11022e-16`；另16组预平衡 p 保持平衡。原门限不改。
- 32 对圆/椭圆变牵引保留非零速度，加入 Gq 后最大速度差 `2.05998e-17`、压力移位误差 `3.43475e-16`，最大独立首步原动量相对残差 `2.64988e-15`。这不是曲面动力响应解析精度或时间收敛验收。
- 4 个零载荷控制、8 组电场/规范/相标号变换、矩形 face 单位检查及20个无部分写回反例通过。Release 构建和 Core/Hydro/EHD 96/96 通过，CLI五项拒绝通过；历史长程组未重跑，既有缺目录 PIC 测试未声称修复。

下一批直接闭合真实 FV 电势→P2 单侧 Maxwell→同一 H 面力→Hydro 的误差对照，复用现有装配与局部采样，不复制另一套solver。当前仍是解析tau+PL H，尚未完成同源高阶几何、毛细联合载荷或原液滴运动/厚度验收；后续导电、电荷输运、电喷雾、RZ 与守恒 AMR 目标保持不变。

## 34. 实际 FV—P2 Maxwell—Hydro 链路及误差传递（2026-09-10）

已完成[实际电势—载荷—Hydro 报告](validation/ehd/QUANTUM_EHD_FV_Maxwell_Hydro_20260910.md)。原全局装配/求解和两种局部采样政策抽成共享测试实现；形状参数留在 fixture，生产模块/App 不变。不新写 solver，不混用 clipped faceCenter 与 strict xΓ 窗口。

- 两遍各14组、28次实际FV、12次N16矩阵检查、98次零初态Hydro，三处Hydro后端均SuperLU，完整JSON一致；最大独立原动量相对残差 `1.07665e-15`、散度 `1.43115e-16`，原门不改。
- 同一PL H上比较oracle及解析/Hermite局部几何×精确点值/VA/VH共七路；保存实际数组，21对先差再范数。独立重算5,488项统计最大差 `4.44089e-16`，28次FV的旧共同诊断全部相同。
- 完整H_VH在N64的速度差占同网格oracle响应2.77%–4.59%；解析几何+精确电势点值仍为2.60%–4.50%。局部P2截断误差已传入流场，不能因残差达标宣布物理验收，也不能靠范数相减分摊误差。
- 旧完整圆432解/椭圆108解/平面及100步解析载荷回归重放字节不变；Release与Core/Hydro/EHD97/97通过。保留一次LNK1104、只读诊断和成功重试；CLI七项拒绝通过。生产252个文件/主exe不变，旧缺目录PIC及历史长程组未声称解决或重跑。

下一步在局部重建层处理精确点值对照中仍明显的误差，保留七路来源、固定窗口和原门限作原样复测；随后合入毛细联合载荷并恢复原静滴/运动/厚度与空间/时间口径。固定同网格oracle不是连续NS真值，同源高阶面积、连续界面力积分、导电/电荷输运、电喷雾、RZ和守恒AMR仍有未完成关口。

## 35. 局部 P3 transmission trace 与真实 FV 误差复测（2026-09-10）

已完成[局部 P3 报告](validation/ehd/QUANTUM_EHD_Cubic_Trace_20260910.md)。新增形状无关七状态局部重构及同一 Hermite cell 的三阶导/κ_s；全局 FV 保持 P2，生产/App 不变。默认七路字节不变，可选 `--cubic` 复用同一25点增加六路，不复制求解器。

- 独立局部回归152,783条断言、96行无改善门测量；真实耦合两轮各14组/28 FV/182 Hydro、12次矩阵核验，全部实际 SuperLU，JSON一致。Core/Hydro/EHD 99/99通过。
- N64完整 C_H_VH 的速度差相对同网格oracle从原2.77%–4.59%降至0.20%–1.07%；三条主路径108个F/u/p比较全部改善，36条网格序列下降，但不是连续NS验收。
- 弱矩15/108项退化，圆ε−=10 N64净力差模增至原3.94倍；真实FV增量已与局部误差竞争。保留全部反例、参考/掩膜/门，未凭总L2掩盖抵消效应。
- 独立审计17,836项统计、1,456,000个有限数组值；可选模式旧828,032个数组元素逐值相同。旧圆/椭圆/默认七路与几何回归字节重放不变。保留首次测试未初始化读取的真实失败及最小修复、LNK1104与成功重试、12项CLI拒绝。

下一步固定局部P3/H/Hydro，隔离全局cut-face的高阶通量修正；P3位移沿每段为二次式，需二点Gauss而非旧中点积分，仍保留TPFA主项和旧bulk/边界/矩阵门。之后联合毛细载荷并回到原静滴、运动、厚度和时空收敛。导电、电荷输运、电喷雾、RZ与守恒AMR仍未完成，持续目标保持活动。

## 36. P3 切面通量与真实 FV—Hydro 原样复测（2026-09-10）

已完成[全局切面修正报告](validation/ehd/QUANTUM_EHD_Cubic_Global_Flux_20260910.md)。抽取P2/P3共用有序面段校验；P3每段两点Gauss，在七维基底内Q−TP后单次投影。全局仍共享原装配/公共SuperLU，仅contrast cut-face升级，bulk/边界/参考/3h掩膜/门限不变；生产App未切换。

- 独立局部通量576主例、20协变、48平面P2退化、2中点反例、30拒绝合同，两轮325,435条断言通过，JSON一致。真值是独立解析线积分，不是同Gauss自证。
- 可选17路两轮各14例/56FV/24dense/238Hydro/1904配对，所有实际后端SuperLU；完整JSON一致，Core/Hydro/EHD101/101通过。旧七路/十三路/完整圆与椭圆及局部P3字节重放不变，旧P2合同通过。
- N64圆/椭圆高介电比界面通量密度误差降至旧值8.64%/11.81%；圆净力差7.54e-5→1.20e-5。N64完整新路径速度差占同网格oracle约0.18%–0.99%。但全网格F/u/p有56/144项退化、弱矩60/144项退化，N16椭圆ε−=10压力误差约增21%，不得宣称全面精度通过。
- 独立审计29,988项统计及1,904,000个有限原始值；旧1,500,032个数组元素逐值一致。保留一次LNK1104及成功重试、14项CLI拒绝；290源码/18exe最终复测前后不变，生产252文件/主exe不变。

下一里程碑直接推进同一φ/H/face frame上的Maxwell—毛细联合载荷与实际App接线：先用常曲率毛细纯梯度核速度不变/压力移位，再换真实曲率；复用现有HydroCoupling时间循环，提取公共不可变几何上下文而非生产依赖测试。须同时补公共修正通量的原两电极Dirichlet+侧壁零Neumann支持，保留旧扩散后端为显式对照，不能把全边解析Dirichlet当原App验收。之后原静滴、电极V/E/F0/F2/净力、厚度和运动时空矩阵原样复测；导电/电荷/电喷雾/RZ/AMR目标不变。

## 37. 公共修正通量的静态混合边界（2026-09-10）

已完成[混合边界报告](validation/ehd/QUANTUM_EHD_Mixed_Boundary_20260910.md)。仅扩展公共校验和声明：允许静态D/严格零N，至少一侧D；N面每个修正系数及常量必须为零。原装配cpp、all-D/nullptr、实际SuperLU和失败事务性不变，不推断介质相别，不作纯N定规。

- 固定62成功、96拒绝共158调用，两轮18,360断言通过，JSON一致。独立完整矩阵残差≤3.22e-15，边界积分通量误差≤1.62e-14，原1e-11门不变；分层真值仅针对严格cell-midface对齐。
- 原7/13/17路、完整圆/椭圆JSON逐字节不变；旧99项通量拒绝及stdout/stderr相同。Release重试成功，Core/Hydro/EHD102/102，保留首次两个LNK1104及7项CLI拒绝证据。
- 本批是真实生产端口修改，但不是sharp App、电极力矩/薄层或联合毛细物理验收。下一步公共几何快照与受控浮点实现边界，再接同源Maxwell—毛细联合载荷及现有时间循环。原物理参考/掩膜/门限、运动/厚度/时空收敛与后续导电/电喷雾/RZ/AMR关口不变。

## 38. 共享 Hermite 几何与窄公共快照（2026-09-10）

已完成[公共几何迁移报告](validation/ehd/QUANTUM_EHD_Shared_Hermite_Geometry_20260910.md)。四引擎集中Core/Geometry/Detail，旧测试仅alias；公开move-only深快照，生产不包含测试和shape参考。4头逆还原与冻结源码逐字节一致，未改运算/门限。内部链接引擎与strict实现对象隔离，其他Core/App浮点选项不变。

- 实际CL/符号检查通过；两轮57 public captures、160预期拒绝、11,357断言及全部字段位比较通过，JSON一致。明确保留Q4共享近似DOF、局部投影/C1、当前x64及禁FP traps的边界。
- 首次整项目Release构建成功，Core/Hydro/EHD103/103。原7/13/17路、完整圆/椭圆、标量/几何/精确限制/混合边界JSON逐字节不变，300sources/20exe最终复测哈希稳定。
- 尚未自动切换生产物理模型。下一批公共快照直接进入126次真实Hydro的Maxwell—毛细联合载荷矩阵，再接原混合电极App与运动/厚度/时空物理验收；后续导电/电喷雾/RZ/AMR目标不变。

## 39. 同源 Maxwell—毛细联合面力（2026-09-10）

已完成[联合载荷报告](validation/ehd/QUANTUM_EHD_Joint_Load_20260910.md)。同一公共快照、原PL H、每活动面一次P3 QR与两个真实电势评价；`(tau-sigma*kappa)G(H)` 共9路，公共κ=+div(n)。共享测试载荷/响应辅助代码，旧campaign不改，不复制solver。

- 两轮各14组/28FV/12dense/126Hydro，实际SuperLU；2119fit/4238evaluation，JSON相同。常毛细静态速度≤8.69e-19、压力移位误差≤5.08e-15，原门保持。
- 467,461个新旧E-only原始数值逐binary64相同；独立重算全部联合数组关系及126组staggered弱矩通过。旧7/13/17、圆/椭圆等JSON原哈希不变，首次Release和Core/Hydro/EHD104/104通过；304源码/21exe最终复测哈希稳定。
- 真实圆曲率仍产生非零首步毛细速度，N16/32/64为7.44e-6/1.73e-6/3.27e-7；真实椭圆运动不是静滴寄生流控制。不得把离散接线通过当作实际电极、稳态或连续物理验收。

第40批已按[生产接入合同](../runs/history/ehd/ehd_sharp_app_20260910/experiment_contract.json)开始：明确选择sharp电学及Hydro Hermite几何，纯Hydro与zero-V保留同毛细路径，共用原HydroCoupling时钟和源快照；公共数值内核不依赖tests/形状参考。保持原finite-box参考、掩膜/门限与默认diffuse对照，先实际App固定界面，再原运动/厚度/时空验证；导电/电喷雾/RZ/AMR持续目标未完成。

## 40. 实际 sharp App 接入及原物理门复测（2026-09-10）

已完成[生产接入报告](validation/ehd/QUANTUM_EHD_Sharp_App_20260910.md)。Core不可变共享几何、Electricity唯一数值内核、Hydro唯一时钟及App独立诊断I/O分层；新算例显式选择sharp，旧diffuse默认不变。公共P1/P3五个头机械逆还原与stage39源码哈希相同，不生产依赖测试，不复制solver。

- 两轮生产合同各12成功/12失败事务性/6浮点拒绝、119,226断言；实际SuperLU。原App check/smoke/zeroV/equalε/runtimeFailure通过，原零电压逐场相等不改。
- Release与Core/Hydro/EHD110/110通过；11个旧完整数值JSON哈希不变，337输入/22exe重放前后稳定。历史PIC与长期运动验证未声称通过。
- 原有限方盒参考、掩膜和门限用于实际App N32/64/128。V/bulkE/F0/F2通过；净力比2.15e-5→2.46e-6→2.87e-7仍未达1e-10，切向泄漏13.48%/12.74%/12.84%仍未达5%。三个分析明确exit2，物理验收未完成。
- 原生面力与单元邻面平均的净力一致，净力偏差来自tau反演非对称；同G(H)常系数诊断仍有约13%切向指标，不可误当局部Maxwell跳跃切向应力。这不改变原验收，也不允许通过改输出、减净力或换参考掩盖问题。

下一窄修先隔离全局face-center窗口方向偏置及电势/牵引反演误差，再处理充分解析载荷层及局部离散平衡。修复后原样复测，再推进运动、时空/守恒、导电电荷输运、电喷雾、RZ与AMR；持续目标未完成。

## 41. 全局物理面窗口方向偏置（2026-09-10）

已完成[窄修及复测报告](validation/ehd/QUANTUM_EHD_Symmetric_Sampling_20260910.md)。sharp全局切面显式选择法向6×切向5的对称实际单元窗口，旧GlobalClipped/TraceStrict及局部25点不变；共享采样实现，无第二装配/solver。

- 独立方形/矩形网格836面/30拒绝、177,964断言通过；两轮真实App的12份电场/面力/φ/metadata文件逐字节相同。原算例、参考、掩膜、门不变。
- N32/64/128净力比降至2.95e-6/9.44e-8/7.71e-10，改善7.27/26.0/372倍，仍未过1e-10；切向指标仍约13%。原其余四门通过，电势总误差略增，不能声称全面改善。
- Release与Core/Hydro/EHD111/111通过；原生产合同和11个旧完整JSON哈希不变，338输入/23exe重放稳定。

下一定位真实全局切面frame/κ_s及最终修正权重的反射缺陷；使用现有装配的可选只读诊断，失败事务性和关闭诊断时原计算保持，不复制装配、不引入shape分支。之后继续载荷层局部平衡及原运动/薄层/守恒验收，整体持续目标未完成。

## 42–43. 曲率导数修正及实际运动推进（2026-09-10）

已完成[曲率导数修正与运动记录](validation/ehd/QUANTUM_EHD_Curvature_Derivative_Motion_20260910.md)。实际装配定位后直接修公共三阶导数重构：节点五点差分加双线性插值，电学全局P3和共享载荷同源；原phi/相别/投影/κ、压力、PL力层及canonical查询保留，不作圆/轴特判。

- 原样N32/64/128净力比降至5.51e-14/2.85e-14/1.72e-13，全部通过原1e-10；切向仍约13%，原分析仍exit2，不放宽门限。
- 真实sharp App已开启界面运动：N32、dt0.001、T0.2，V0和V0.1各200步全部成功，逐步SuperLU。终态D为0和2.03e-4；原轮廓面积漂移为-1.32e-6和-6.07e-6。实际动量平衡残差小，但零场仍有约2.06e-5寄生速度，不能替代连续物理验收。
- 本轮仅构建主程序与两个相关目标、运行两个定向CTest（均通过），完整回归留在里程碑，不再为每个局部改动扩张诊断体系。

接下来用实际运动时空对照推进，并保留PL力层切向误差待办；导电/电荷输运/电喷雾/RZ/AMR尚未实现验收，持续目标不变。

## 44. 充分解析载荷层与两网格真实运动（2026-09-11）

已完成[平滑载荷层实现与结果](validation/ehd/QUANTUM_EHD_Resolved_Load_Layer_20260911.md)。显式regularized选择使用同公共phi的紧支撑C2指标、同压力G；原PL面积份额、电场装配和默认算例保留，宽带触壁须严格证书、不截尾。

- 原N32/64/128切向指标从约13%降至7.834%/2.221%/0.586%，观测阶1.818/1.923；N64/N128原全部门通过，N32仍失败不排除。N128半宽从.04降至.02仍全部通过，F0/F2相对变化约3.22e-6/4.91e-6。
- regularized真实N32/N64、V0/.1各200步，全部SuperLU成功；有场终态D=2.0203e-4/2.0682e-4，网格变化2.319%；面积漂移约-6.05e-6/-2.66e-6。零场寄生速度随加密下降约5.43倍，但不为零。
- PL基线半dt对照完成，D变化0.01039%；不冒充新力层时间收敛。两项C++定向回归及28项原Python测试通过，未扩新目标或完整矩阵。

下一主线接有限电荷弛豫/界面输运，并保留regularized时间/进一步空间及独立动态基准待办。旧导电电势实际上是RZ制造解、发射被禁用，不能靠启用旧循环交付；需在受检公共设施上实现真实电荷守恒闭合，材料/发射参数与算例分离。整体目标未完成。

## 45. 带电泊松、完整 Maxwell 牵引与守恒电荷更新（2026-09-11）

本批直接推进导电所需的生产代码，不再扩展理想介电扫描。regularized N64半dt对照已结束，各电压400步全部SuperLU；终态形变在原输出精度下相同，面积漂移变化约6.04e-9，详见[原结果补记](validation/ehd/QUANTUM_EHD_Resolved_Load_Layer_20260911.md)。这不是正式时间收敛阶，也不替代剩余独立动态基准。

已实现：

- `SurfaceChargeTransport`保存闭合、有序界面各段总电荷Q，以迎风共享面通量更新 `Q_new=Q+dt*(I+F_in-F_out)`；I约定为负相流入减正相流出。伸缩通过新段长恢复密度，绝不全局缩放/减净电荷。显式对流CFL、有限数值与下溢检查、守恒残差及失败不改状态。它不另建网格/水平集/时钟；段身份、几何和传导电流仍由上层提供。
- 现有P3增加给定 `q=[Dn]` 及沿界面一、二阶导数的仿射跳跃；同一全局FV/SuperLU求解加入单元内界面总电荷 `Q/area` RHS和仿射面通量。独立高斯检查从无电荷的sumD=0扩展为sumD=Q，默认nullptr原路径保留。`SurfaceChargeSource`提供同一电荷场的单元积分和局部jet；二者几何/积分一致性由将来的电荷场适配层保证，不假装此接口已经构造了电荷面板。
- 使用两侧实际重构电场的完整Maxwell张量牵引，验证法向恒等式和 `tau_t=q*Et`。Core的可选`includeTangentialLoads`提供同一H、同一压力G的四面横向插值，电力为 `tau_n*GH+tau_t*R(GH)`。它是交错梯度插值近似，不把离散GH方向冒称为精确几何法向。默认关闭，原面集合/算术不变；边界零延拓需要严格H平台证书。Hydro毛细力允许新面主梯度为零，不向切向力混入毛细项。

验证仅覆盖本次实现：主程序与相关目标Release构建成功；[三项定向CTest](../runs/history/ehd/ehd_surface_charge_20260911/targeted_ctest.log)全部通过，总运行约5.5秒。原sharp生产12成功/12失败/6浮点环境拒绝、regularized原覆盖保留；新增加一个实际带电制造场及5项输入拒绝，局部P3增加16个带电补丁/8项非法输入。没有重跑完整矩阵。

带电制造场在两相均为`V=.1y`、介电常数不同，所需 `q=-.1*(epsilon_positive-epsilon_negative)*n_y`；每格Q独立由公共上下边界切分及高斯定理积分，不输入解析电势给求解器。实际求解电势最大误差约6.31439e-16，切向牵引最大约0.00499977；两侧Dn、Et、独立法向/切向Maxwell真值、横向G和最终面力通过。[首轮实际结果](../runs/history/ehd/ehd_surface_charge_20260911/sharp_charge_first.json)保留；随后仅补独立法向牵引断言，最终三项CTest仍通过。输运测试中的指数弛豫使用外加电流，仅证明该显式更新的一阶时间行为，不能当成自洽充电验收。

**当前边界：给定表面电荷→全局电场→完整界面电力已接通；电荷输运内核已实现，但尚未形成导电App充电闭环。** 下一交付是从同一公共几何构造有序电荷积分段，将一个守恒电荷场同时供给Q和局部q jet，再接有限电导的法向电流反馈与固定界面充电算例（先关闭发射）。不得用稳态电导方程取代有限弛豫，外相零电导也须保留电场。随后接运动/对流、电喷雾发射、RZ和AMR；真实发射模型或参数需要用户选择时再明确报告。整体持续目标未完成。

## 46. 固定界面有限充电 App 闭环（2026-09-11）

**已跑通实际充电反馈，不再只是给定电荷的静电求解。** 新算例见[cases/ehd_leaky_charging_planar](../cases/ehd_leaky_charging_planar/README.md)：负相sigma=1、正相sigma=0，两相epsilon分别2/1，保留外相电场；从q=0开始，实际传导电流更新Q并反馈下一步电场，完整Maxwell力进入原Hydro。参数仍全部在算例字典中。

生产分层保持：

- Core新增`LevelSetInterfacePanels`，使用同一公共Hermite快照的共享边切分构造单闭环，每切单元一段，按t定向；GL8在实际零曲线的单调图上积分，不用解析圆或折线替代。结果只保留几何值，不跨步借用快照。目前明确拒绝碰壁、顶点/沿边零、多交点、多分支、多闭环及未认证隐藏拓扑。
- Electricity的`SurfaceChargeField`在固定段上保存Q；三邻段的二次平均值重构同时给出单元积分和局部q/导数。局部弧坐标用保留求积数据的Hermite插值，二次电荷重构是分片近似，不冒称全局C2。类禁止复制/移动以保护source回调生命周期；旧source保持旧时间层。
- 传导电流来自同一物理FV面：混相P3分相E通量复用原QR，再分配实际D与P3 D的差，确保两相epsilon加权重组回实际D。`J=sigma_negative*Eflux_negative+sigma_positive*Eflux_positive`，每格I=−sumJ，与电极电流闭合；bulk电流不移到界面，只允许并记录求解舍入量。每次成功报告单独返回只读电流，失败没有旧电流缓存可误用。
- App使用Hydro唯一时钟。t^n电势/电力/电流先计算，Hydro步成功后由通用`completeStep`回调提交Q^(n+1)；不提前写成功电荷状态。首版显式更新，检查dt*max(sigma/epsilon)<=.25，明确固定界面、`chargeConvection false`及无重初始化。流体速度仍求解，但本分项算例关闭其电荷对流反馈，不能称为完整运动导电液滴。

实际运行：N32、44个界面积分段、T=.2；原dt=.01共20步与半dt=.005共40步均exit0，全部60步Hydro成功、电学/耦合SuperLU、诊断有限。证据分别在[dt=.01](../runs/history/ehd/ehd_charging_20260911/final_dt001/output/surface_charge_diagnostics.csv)及[dt=.005](../runs/history/ehd/ehd_charging_20260911/half_dt0005/output/surface_charge_diagnostics.csv)。

| 终态或全过程指标 | dt=.01 | dt=.005 |
|---|---:|---:|
| 最大面板平均电荷密度 | .013385044610742493 | .013375050455977610 |
| 总净电荷 | −4.62954e-17 | −4.78133e-17 |
| 最大单步电荷—电极电流平衡误差 | 3.81842e-18 | 2.08709e-18 |
| 最大bulk绝对电流（未分配） | 2.91434e-15 | 2.82413e-15 |

两档终态电荷密度的弧长加权相对L2差为7.5276943e-4（0.07528%）。这证明本短时步长敏感性较小，不是正式收敛阶、稳态或独立连续充电基准验收。外相sigma=0时电场仍非零，dt=.01末个力时刻maxE约.135056。

必要回归：主程序Release构建成功；[两个定向CTest](../runs/history/ehd/ehd_charging_20260911/final_targeted_ctest.log)通过（约6.6秒），包含原sharp覆盖、独立线性制造场逐格传导电流与非零电极电流符号、快照销毁后使用及重构q的GL8积分等于每段Q。未扩展参数矩阵。另[故意压低流体残差门](../runs/history/ehd/ehd_charging_20260911/fail_momentum/run.log)使第1步在momentum-solve失败、exit3；电荷CSV只有表头，未提交新电荷状态。原失败输出保留。

诊断修正：上述两档数值运行的早期sharp metadata仍有旧`perfectDielectric`硬编码标签（实际输入、运行日志及电荷文件均为leakyDielectric）；原文件保留、不据该错误标签识别物理模型。现已改为实际model，并在带电face输出追加正侧Dn/q/tau_t/横向GH，旧无电荷列不变。[修正后的两步实际输出](../runs/history/ehd/ehd_charging_20260911/metadata_check/output/ehd_charge_electric_sharp_step0002_metadata.json)确认leakyDielectric；324个载荷面最大Dn跳跃误差1.39e-17、tau_t=qEt误差0，头两步电荷峰值与原运行相同。这是I/O修正检查，不额外宣称一次物理精度验收。

下一主线：先接固定界面上的真实速度电荷对流，再实现运动后积分段变化的守恒重映射，使电荷与运动界面真正双向耦合；补独立充电参考及必要空间/时间证据，但不以大规模扫描阻塞实现。导电喷雾发射、RZ和AMR仍未完成，整体持续目标保持。

## 47. 真实速度对流与运动导电液滴闭环（2026-09-11）

**已接通并实跑运动导电液滴，包含跨网格线后界面分段数量变化。** 本轮以生产实现为主，没有扩展理想介电参数矩阵。可直接运行[固定界面对流算例](../cases/ehd_leaky_advection_planar/README.md)及[运动导电液滴算例](../cases/ehd_leaky_moving_planar/README.md)。

实现分层：Core新增独立的交错速度快照采样；公共`ImmersedGeometry::getNextCoord`反解原水平集数值回溯，复用其算术，严格路径拒绝裁剪/无效查格/超半格位移/不收敛，旧回溯默认行为保留。Electricity在已有`SurfaceChargeTransport`中增加直接积分Q工厂和周期弧长重叠重映射；`SurfaceChargeField`复用自身弧长映射，检查单调、局部投影、方向和一周绕数，候选新场完整成功后才替换旧场。App只连接真实前后速度、实际新水平集和唯一Hydro时钟，没有几何特判、第二个流体求解器或独立时钟。

固定界面用实际native-face旧/新速度平均值的端点u·t驱动迎风通量。运动路径则使用与水平集一致的顶点速度和回溯逆映射：先将t^n电流加入旧段总电荷，再随材料端点移动并投到实际新界面，按新弧长上的重叠长度分配总Q；**不再叠加固定界面的u·t通量**。因此运动输出的Euler通量/CFL为0，但`remapped=1`及材料位移非零。分段总Q不经density往返，也不全局减净电荷/缩放；新密度使用原GL8几何长度。供体为分片常数近似，不能当作高阶重映射。

实际App结果（均N32、dt=.01、20步至T=.2，实际SuperLU）：

| 运行 | 界面/速度初值幅度 | 终态最大平均电荷密度 | 最大单步电荷—电极平衡误差 |
|---|---|---:|---:|
| [关闭对流](../runs/history/ehd/ehd_charge_advection_20260911/off/run.log) | 固定 / .1 | .0133850446107425 | 3.82e-18 |
| [开启对流](../runs/history/ehd/ehd_charge_advection_20260911/on/run.log) | 固定 / .1 | .0133884054294053 | 3.82e-18 |
| [运动](../runs/history/ehd/ehd_charge_advection_20260911/moving/run.log) | 移动 / .1 | .0134022598827729 | 6.51e-18 |
| [运动并跨网格线](../runs/history/ehd/ehd_charge_advection_20260911/moving_crossing/run.log) | 移动 / .2 | .0134188006807043 | 5.39e-18 |

开/关对流仅一项配置不同，终态电荷密度弧长加权相对L2差为.00744105（0.7441%）；最大实际CFL=.0876666，绝对电荷通量和最大约.00270791。两档运动不是精度收敛对照，而是集成运行：幅度.2的第12步面元数44→46，该步重映射电荷残差5.42e-20；全程最大重映射残差1.63e-19、最大材料位移.00123053、最大新界面投影修正9.01e-5（h=.03125）。幅度.1的最大重映射残差2.71e-19。投影修正被显式记录，并不等于物理轨迹误差已经验收。

必要验证：[最终transport定向CTest](../runs/history/ehd/ehd_charge_advection_20260911/final_transport_ctest.log)通过，约1.6秒，保留原充电/守恒/GL8/失败无写回检查，新增实际面速度采样、非零切向输运及超CFL拒绝，两例signed-Q变段数周期重叠及失败不变，以及一个old/new仿射速度的独立数值轨迹逆映射。首轮链接发现Core遗漏Field源目录，已补构建收录并成功重新构建，不是运行时换求解后端。

当前界限：这80步证明闭环、分段转移和离散电荷守恒已运行，不代表导电液滴物理精度或全部流体/电喷雾验收完成。smooth初速度并非无散场；本例是接线/分段变化验证，不能替代零初速物理基准。尚需独立充电/运动参考与必要时空证据，重初始化/接壁/多分支/破裂拓扑仍明确拒绝；发射、RZ和AMR也未完成。下一主线在此完整App上做紧凑的导电物理验收，然后推进RZ及电喷雾入口/电荷输出边界；真实发射模型/材料参数需要用户选择时明确报告。保持旧参考、掩膜、门限和失败记录，不再以扩展测试框架代替功能推进。

## 48. 物理电荷弛豫与 RZ 共享载荷接通（2026-09-11）

本轮完成一项物理时间基准，并继续落地RZ公共载荷代码；没有扩展理想介电扫描。详见[电荷弛豫报告](validation/ehd/QUANTUM_EHD_Charge_Relaxation_20260911.md)及[永久算例](../cases/ehd_charge_relaxation_planar/README.md)。

两相sigma/epsilon均为1，初始均匀q0=.02，固定界面且关闭对流。由q=[Dn]、J=D，独立连续解为q=q0·exp(−t)，电极净向外传导电流为Q；外加电压保持实际求解。真实App以dt=.05/.025/.0125分别完成20/40/80步至T=1，密度连续误差2.553423%/1.263186%/.628276%，观测阶1.01537/1.00760。全140步电极电流相对Q_before的最大绝对差3.81e-15，单步电荷平衡最大误差1.41e-16。17个已输出面元帧共748行与独立Euler递推最大偏差1.23e-14，连续参考与离散参考没有混用。此基准支持同弛豫率的物理电荷更新和时间收敛，不证明一般电导率比的空间精度或运动电喷雾。

RZ生产推进：Core的`InterfaceLoadGeometry::Options`显式选择axisymmetric，仅开放离轴内部regularized界面。保留frame子午线曲率与ks、原GH和原二维PL分数；新增`surfaceCurvature=kappa_m+n_r/r_interface`与`dualMeasure=2*pi*r_face*dualArea`。Hydro同一毛细力入口使用三维总曲率及物理积分测度，力密度不再乘r，继续走原RZ NS和唯一水平集推进。算例参数与几何仍在字典，见[hydrodynamics_hermite_rz](../cases/hydrodynamics_hermite_rz/README.md)。原默认planar的新增字段只是原值赋值，其梯度、frame和曲率导数不变。

[RZ实际运动运行](../runs/history/hydro/hydro_hermite_rz_20260911/final_moving/run.log)完成N32、dt=.001、20步至T=.02，动量/压力/耦合修正均SuperLU。最大散度4.23e-16、最大完整动量平衡相对残差2.29e-16、轴上径向速度0；终态maxSpeed=.00110291、界面最大变化1.06278e-5，原体积口径相对漂移−2.36e-6。环面总曲率不恒定，本例不能按静滴速度为零验收。首个运行目录`moving`的App选择拼写错误导致exit2，配置已改为既有SingleModule/Hydrodynamics入口，原失败记录保留；未覆盖后续`final_moving`。

**新增RZ曲率精度门仍失败，不计为验证通过。** [完整定向回归记录](../runs/history/hydro/hydro_hermite_rz_20260911/curvature_detail_ctest.log)保留旧12成功/12失败/6环境、regularized3/5/1984、带电1/5全部覆盖；新增227个RZ面完成共享几何逐位一致、径向度量、毛细面力/独立积分以及两项范围拒绝检查。独立解析环面总曲率的最大相对误差为2.24595%、最大绝对误差.0566871，未达到预设1%。最差相对误差点约(r,z)=(.322903,.532199)：子午线κ=5.61197（圆真值5.55556）、周向项−3.04703，总曲率2.56494而参考2.50860；证据显示主要误差来自现有粗网格子午线曲率，在环面内侧相减后放大。门限、N32、参考未放宽或替换，CTest仍报告失败。

RZ电学仍未接通，新增库层检查也拒绝把RZ载荷交给平面电场，失败不写回。下一实现项是公共FV的独立径向face/cell度量和满足`V_rr+V_zz+V_r/r=0`的局部P3基底，再接真实RZ理想介电App；不能只给平面通量乘r。有限电荷后续还需区分弧长ds与表面积2*pi*r*ds，接轴界面的非周期面板/正则极限另行实现。N32曲率精度失败继续保留待处理；发射、RZ电流体、AMR及完整物理验收均未完成，持续目标保持。

## 49. RZ 电流体 App 实际运动接通（2026-09-11）

**RZ 电场—Maxwell 面力—NS—运动界面已接通，真实完成 20 步。** 详见[RZ App 报告](validation/ehd/QUANTUM_EHD_RZ_App_20260911.md)及[永久环面算例](../cases/ehd_dielectric_rz/README.md)。公共 FV 的物理面/体度量、局部 RZ P3 方程与逐点 r 加权面通量均已实现，不是给平面解事后乘半径；共享几何、原 Hydro 时钟、实际 SuperLU 及默认 planar 路径保持。坐标制只取 Hydro 配置，轴电极错误在预检拒绝。

N32、dt=.001 至 T=.02，全 20 步电学/流体成功且诊断有限；最大电通量相对残差 1.61e-15、电极通量不平衡 3.56e-14、投影后散度 3.31e-16、完整动量平衡残差 2.18e-16，轴上径向速度0。终态 maxSpeed=.00102786448、界面最大变化9.92051e-6。环面非静态球滴，不声称零寄生流或物理精度已验收。

仅运行三个已有定向目标（约8秒）：公共FV非零RHS制造解、局部RZ二/三次调和多项式通过；sharp新增三个RZ实际求解通过，等介电电势误差1.12e-15、物理电极通量误差2.42e-14、独立原生面力积分差1.30e-18。**原N32曲率2.24595%仍超过1%，sharp总目标继续失败，门限与失败未改。** 未扩展参数扫描或完整回归矩阵。

下一实现为RZ有限电荷/运动输运的面积测度，再处理接轴面板、入口/出口及电喷雾边界；现阶段RZ电荷、接轴/接壁界面、发射与AMR仍未完成。曲率误差和独立时空/动态基准继续作为明确待办，不以本次短时接线运行代替完整目标。

## 50. RZ 导电电荷与运动面积重映射接通（2026-09-11）

**RZ 有限导电—完整电力—运动界面—电荷守恒闭环已实跑。** 详见[RZ 电荷报告](validation/ehd/QUANTUM_EHD_RZ_Charge_20260911.md)、[固定弛豫算例](../cases/ehd_charge_relaxation_rz/README.md)和[运动导电环面](../cases/ehd_leaky_moving_rz/README.md)。Core 保留原 ds 几何并追加同一零曲线的 GL8 物理面积；电荷 Q/平均值矩/运动重映射采用 `2πr ds`，q 的导数仍对弧长。固定端点通量为 `2πr u_t`，运动不叠加第二次 Euler 对流。Sharp 使用 Q/物理体积与物理分相电流，源/几何/电场坐标制必须一致。共用原网格、水平集、NS、时钟及 SuperLU，不另造喷雾求解路径。

两个真实 App 各完成20步且数值诊断有限。RZ 同弛豫率算例 dt=.05 至T=1，实际q=.007169718448187798，对连续指数解误差2.553423%，对独立Euler递推的面元密度误差至多1.70e-14；两种误差明确分开。运动算例dt=.01至T=.2、初始非零总Q=.071061176640578289，终态Q=.071061176640578746；第14步44→46面元，最大单步电荷平衡误差4.16e-17、最大重映射残差1.39e-17。短时守恒不代替运动/空间物理验收。

仅运行两个已有定向目标：transport RZ面积/矩/端点通量/伸缩重映射与失败检查通过；sharp新RZ带电制造场电势误差9.37e-16，逐格传导电流与完整Maxwell检查通过。**原N32曲率2.24595%仍未过1%，sharp总目标继续失败。** 首轮测试lambda捕获编译错误已修正，原记录保留；未扩展完整参数矩阵。

下一实现转向接轴开放界面链与轴正则极限，再接喷嘴入口/出口和发射条件；当前离轴闭环面积输运将直接复用。接壁/破裂/多分支、重初始化、发射、AMR和完整曲率/时空/动态物理验收仍未完成，持续目标保持。

## 51. 固定接轴导电球滴 App 接通（2026-09-11）

**真实接轴球滴已完成电场—完整 Maxwell 面力—NS—有限电荷弛豫的20步运行。** 见[交付与复测记录](validation/ehd/QUANTUM_EHD_Axis_Sphere_20260911.md)及[永久算例](../cases/ehd_charge_relaxation_sphere_rz/README.md)。本批围绕实际App补齐公共Hermite径向偶性、真实零面积轴电通量、两端接轴开放面元及极点q_s=0重构、输入镜像轮廓与Hydro轴边界接口。共用原网格、GL8面积、FV、压力兼容面力和SuperLU，没有新建喷雾专用求解器。

N32、R=.18、dt=.05至T=1，22个面元；初/终Q=.008143021279328281/.002919158494512298，保存66行面元密度对独立Euler递推最大差4.80e-15，而连续指数解误差仍2.553423%。最大电荷平衡误差8.56e-16、散度3.45e-15、完整动量相对残差1.75e-15；界面固定、流体仍运动，不当成静滴或动态精度验收。

接轴拓扑不定义旧闭环节点/Peskin参考，因此只在新轴路径将旧八项参考诊断留空并标记不可用，活跃面力和完整动量诊断保持；原离轴算例原样复跑20步，其hydro/electric/charge全部原列逐项一致。仅两个已有定向目标复测，新增接轴覆盖通过；**原N32环面曲率2.24595%>1%仍是唯一失败，门限未改。** 旧弥散negative_volume不作为解析球体积证据。

下一主线为接轴水平集/材料端点一致运动与开放链面积重映射，再接喷嘴供液/出口；接轴运动目前明确预检拒绝。曲率精度、动态及时空基准、发射、破裂、重初始化及AMR仍未完成，整体目标不变。

## 52. 接轴运动带电球滴接通；下一批转向供液/出口（2026-09-11）

**运动球滴完整耦合已实跑20步，不再只做固定形状电荷更新。** 见[本批交付](validation/ehd/QUANTUM_EHD_Axis_Moving_20260911.md)与[永久算例](../cases/ehd_leaky_moving_sphere_rz/README.md)。公共回溯/材料轨迹补轴边界模式，开链追踪N+1端点并按2πr ds守恒重映射；首轮实跑发现旧面→顶点无约束外推破坏轴ur=0，已在Core共用速度接口施加径向奇/轴向偶重构，Hydro与EHD共同调用，默认旧路径不变。

原N32、dt=.01至T=.2真实完成：两极r严格0且z移动，step20面元22→24，总Q最大漂移8.67e-18，最大重映射残差1.73e-18，散度7.04e-14，完整动量相对残差1.38e-15；电场/流体SuperLU。原离轴运动算例原样复跑20步，全部原有流体/电场/电荷数值列及保存面元值逐项一致。只运行两个已有定向目标；transport通过，**sharp唯一失败仍为N32环面曲率2.24595%>1%**。短时运行及守恒不替代物理精度验收，首轮失败保留。

下一批直接补通用速度入口与压力/应力出口，接入原动量、压力与coupled连续性装配，先交付单相实际供液/排出，再接喷口界面入流与电荷端通量。现新App依然封闭域；不能靠放开几何预检或复用旧legacy入口宣称电喷雾已通。喷嘴/发射、破裂、重初始化、AMR及曲率/动态/时空物理验收仍未完成。后续按实际功能里程碑推进，不扩展无关测试矩阵。

## 53. 单相供液/压力牵引出口实跑；补全公共 RZ 周向黏性项（2026-09-11）

**平面与RZ单相供液—排出各完成20步，已突破新App封闭域限制。** 见[本批实现与独立证据](validation/hydro/QUANTUM_Hydro_Open_Boundary_20260911.md)及[平面](../cases/hydrodynamics_inlet_outlet_planar/README.md)/[RZ](../cases/hydrodynamics_inlet_outlet_rz/README.md)永久算例。静态patch由字典给出；入口指定向外法向速度，出口给完整压力牵引(-pI+2μD)n=-p_out*n，切向牵引0。出口真实半控制体、边界连续性、压力负载和角点动量通量接入原CN/AB2—projection—coupled，不固定出口速度，不删有出口时的cell0连续性或减p0。两相入口仍明确拒绝。

N12、dt=.01至T=.2、入流速度.05、外部压力.3，实际进出流量平面.05/RZ .05π；最大全域体积流量误差7.81e-18/6.25e-17，全部SuperLU。现有coupled-pressure目标的原闭域/动态检查保留，并通过固定规模通流、完整非线性斜向角点通流、两幅值独立线性RZ应变/黏性牵引检查；未扩大网格时间矩阵。轴入口预检exit2、缺失回流状态运行exit3，失败可见。

**重要修正：旧受支持RZ径向动量漏了-2μur/r²，已对所有checked固定均匀RZ补齐，不能只给新patch路径补丁。** 原运动接轴球滴原样复跑20步，终态maxSpeed从.05591247480变为.03397537998，φ变化从.01614532794变为.01223564439，面元24变为22；Q最大漂移3.47e-18。旧stage49–52的RZ运动/黏性数值需按该缺项重评，不以小线性残差或历史短时成功充当物理正确性证书，历史日志保留。

下一步接喷口相别、水平集入流/唇口条件和非轴开链电荷端通量，继而接壁电学；真实喷口参数及边界模型需在算例接线前明确。内部固体喷嘴、无滑移/接触角、真实出口回流、发射、破裂、重初始化和AMR仍未完成。原N32曲率2.24595%>1%未修改，本批没有复跑整个sharp目标；RZ径向出口及完整空间/时间/动态精度验收也未完成。持续目标保持，不把本次供液跑通视为完整电喷雾完成。

## 54. 喷口非轴端点与电荷边界输运接入公共几何/Field（2026-09-11）

本批是喷口接口实现，不是新一轮液滴参数扫描，也**尚不是运动喷口 App**。原电喷雾输入没有可沿用的动态唇口/接触线规律，已询问首个可配置模型是否采用“固定唇口＋给定供液速度，接触角随流动变化”；未把初始几何角度暗设为永久接触角。

- Core 的 `LevelSetInterfacePanels(..., boundaryContacts)` 支持显式声明的非轴边界端点与轴端连接。精确网格顶点使用共享节点身份，不移动点或修改 φ；只接触零角点的相邻单元须通过共享 Hermite DOF 的精确 Bernstein 同相证书。完整网格边求根单独保留零端点并区分整边重合，不再用相分段推断完整零点清单。未声明接壁、沿边重合、分支及不明确拓扑仍拒绝；空声明保持原路径。声明只标识几何，不施加钉扎物理约束。
- `SurfaceChargeTransport` 新显式边界状态重载允许非零端面积通量：入流必须给外部 q，出流取端面元 q，按真实相对面积通量上风更新。边界电荷输入/输出/净输入与法向导电电流分别报告，守恒式为 `ΔQ = dt*(ΣI_normal + I_boundary)`，尺度包含两端绝对输运量，CFL 包含端点损失；失败不提交状态。旧重载仍要求开链两端精确零通量。
- 能力已接入实际 `SurfaceChargeField`，不是孤立 kernel。非轴唇口使用普通单侧二次电荷重构，只有 Axis 端使用偶四次及零面积极限；面元平均值、局部 q 及导数来自同一重构。即使显式给了外部 q，轴端非零通量仍拒绝。声明链的材料重映射、自动壁面速度采样尚未实现，入口明确失败，不能借旧“两材料极点”的覆盖约束冒充喷口运动。

主程序和既有 transport 目标[最终构建](../runs/history/electrospray/nozzle_endpoint_interfaces_20260911/build_field.log)成功；只运行这一个[定向回归](../runs/history/electrospray/nozzle_endpoint_interfaces_20260911/transport_field.log)，10.85秒通过。新增固定 N32 球冠几何使用测试节点值 `φ=r²+z²+.25z−.25`，唇口为精确节点(.5,0)，得到27个面元；同一零曲线 GL8 面积对独立球冠面积的相对误差2.22e-16。在实际 Field 上给唇口相对切向速度.01、入流q=.25、dt=.001、法向电流0，新增总Q=7.853981633974483e-6，与 `dt*2πr*u_t*q` 相符，逐面元重构积分一致。另检查正负电荷、反向端流、缺入流状态、非零轴流和失败不改状态。**这是保留几何上的边界输运与重构检查，不是距离初始化、运动喷口、电场接壁或真实曲率精度验收。** 首次 kernel/几何构建和测试也保存在同目录，未覆盖。

原[接轴运动球滴](../runs/history/electrospray/nozzle_endpoint_interfaces_20260911/rz_original/run.log)以原5个 system 文件、原20步重新完成，电学/动量/压力/coupled 均为 SuperLU。对照第53节 hoop 修正后的基线，4份诊断文件及3份保存的电荷面元文件全部逐字节相同；没有用新端点功能改写旧运动结果。

下一实现仍是**无电喷口供液时界面真实生长**：接通公共水平集入口回溯/边界值及唇口精确零值保护、接壁界面力与压力边界收支，然后接声明端点的运动电荷通量和接壁电学。当前 App 的两相入口限制没有放开；上面测试直接给节点 φ，不声称旧初始化的零值扰动已修复。材料/供液/润湿或钉扎条件留在算例配置，等待必要物理选择；本批未更改或重跑原 N32 曲率门，最近失败仍为2.24595%>1%，AMR及完整物理/时空验收仍未完成。

## 55. 水平集入流与 NS 壁面速度接线（2026-09-11）

本批继续补无电喷口的真实入口链路，没有扩展液滴参数扫描。**公共水平集边界与速度转换已实现；两相喷嘴 App 尚未贯通，不能把下列回归当成喷嘴供液验收。**

- Core 的 `LevelSetBoundaryConditions` 显式声明入流 φ(point,time) 与固定零值顶点。初始化和后续检查保留声明节点的精确值；回溯实际到达入口时使用到达时刻的边界值，不裁回域内，也不采样域外速度。内部完整步保留原 RK2 算术；只有出域中点才缩短子步，剩余时间继续推进。入口时间、固定拓扑及所有 moving LS 的候选值先校验，再提交 current/pre/prePre 和诊断，提交复用已有历史缓冲。`refreshLevelSetDerivatives` 只更新导数，保留已有时间层；新模式明确拒绝不支持的旧 cut-length、redistance/reinit 和 cache 路径。该约束不等于接触角模型，也不把唇口当材料端点。
- Core 的 `StaggeredVelocitySampler` 增加四条显式切向壁面迹，填补交错网格缺失的半格条带；物理边界值由调用方给出，Core 不选择滑移模型。Hydro 新 `HydroBoundaryVelocity` 从实际 native 速度重构这些迹，与原动量装配共用零切向牵引公式；角点使用同侧单侧差，RZ 轴保留 ur=0 和原轴向偶延拓。Hydro 输出/界面推进与 EHD 运动电荷回调使用同一转换入口，无 patch 时仍走原路径。没有新增网格、时钟或喷雾专用求解器。

[主程序与既有 transport 目标构建](../runs/history/hydro/levelset_inflow_20260911/build_final.log)成功；只运行一个[定向 CTest](../runs/history/hydro/levelset_inflow_20260911/transport.log)，11.28秒通过。新增固定 N32 检查包括实际公共圆弧初始化、入口精确零值、两步仿射 φ/三层历史、非零到达时间0.0032668（独立二次方程参考）、失败不写回及合法角点映射被更改后的拒绝；另以独立非恒定零剪切仿射速度核对四壁、角点、旧/新时间层及原 native 插值算术。它们是运动学和离散接线检查，不是曲率/接触线物理基准。

实际 App 原样复跑：

- [平面供液](../runs/history/hydro/levelset_inflow_20260911/planar_supply/run.log)与[RZ供液](../runs/history/hydro/levelset_inflow_20260911/rz_supply/run.log)各20步、T=.2、SuperLU成功，新壁面转换用于实际场输出；各自输入与第53节最终基线一致，两份 hydro_diagnostics 逐字节一致。
- [原接轴运动带电球滴](../runs/history/hydro/levelset_inflow_20260911/rz_original/run.log)20步成功、实际电学/流体SuperLU；与第54节相比5份输入、4份诊断及3份保存的电荷面元文件全部逐字节一致。

下一交付仍按实际能力推进：接壁面毛细力与压力/边界收支，再把上述 LS 入口声明接入两相 Hydro 的字典和时间循环，交付**无电喷嘴持续供液与界面增长**；之后在同一算例开启电场和运动端点电荷输运。现有两相入口/接壁电学限制未移除，固定唇口等物理选择尚待确认。原 N32 曲率2.24595%>1%的失败、AMR、破裂、重初始化和完整时空/物理精度验收继续保留，本批不另开全矩阵复测。

## 56. 接壁毛细力、端点 KKT 投影与独立压力恢复（2026-09-11）

本批将共享 Hermite 毛细载荷接到显式给定法向速度的边界，复用原局部压力梯度 G、RZ 物理测度和 Hydro 求解器，没有新增喷嘴专用求解器。`InterfaceLoadGeometry::Options` 的 `normalVelocityBoundaries` 与 `boundaryContacts` 必须共同显式给出；接触拓扑由真实零曲线开链检查，触壁载荷带只能落在声明的完整受约束 native 边界面区间。未声明边界仍要求原闭单元带外证书，压力出口不因此获准接触，带切向载荷的新路径仍拒绝；默认空声明保留原路径。

共享 `Snapshot::boundaryProjectedFrame` 区分内部最近点与声明端点。端点使用向可行曲线方向的单侧 KKT 条件，不强迫普通切向 stationarity 为零；输出分别记录 `boundaryEndpointProjection`、`endpointKktResidual` 和原始 stationarity。没有把投影点裁到墙上、用解析圆替换实际 Hermite 曲线，或把几何端点声明当成机械接触线定律。

[最终主程序构建](../runs/history/hydro/capillary_boundary_20260911/build_main_final.log)成功；既有 [coupled-pressure 回归](../runs/history/hydro/capillary_boundary_20260911/coupled_pressure_roundoff_fix.log)的 `all` 模式通过，CTest 总计4.69秒，原测试与门限保留。新增检查只用固定 N32 球冠 `φ=r²+z²+.25z−.25`、唇口(.5,0)、`σ=.03125`，south 法向速度0、north 外压.875，初始压力与速度均为0；载荷半宽1/32、1/8各实际求解2步。参考独立由解析 φ、展开的 C2 quintic H 和 `κ=2/√.265625` 构成 `p*=p_out+σκ(1−H)`，不由生产 indicator/曲率、矩阵或求解结果反造；比较绝对压力，不减 gauge。正则过渡层是离散载荷参考，只有纯相压力跳跃对应 Young–Laplace 真值。

- 半宽1/32：最大速度5.22016e-16，绝对压力误差2.97540e-14，原 local-G 误差4.86722e-13，独立散度6.59799e-16；没有端点投影。
- 半宽1/8：实际1个载荷面使用端点投影，KKT残差0而原始 stationarity 为0.0075792382823855441；独立解析单侧切向不等式成立。最大速度2.15269e-16，绝对压力误差2.12053e-14，local-G误差2.19491e-13，独立散度3.74453e-16。两档均实际使用 predictor/projection/coupled SuperLU，最大完整动量相对残差1.36213e-15。

首轮静态平衡曾被出口回流门误判：[保留的失败诊断](../runs/history/hydro/capillary_boundary_20260911/open_failure_diagnostic.log)中真实最终速度仅4.22012e-16，出口最小向外速度−3.74945e-17。原门把舍入尺度也乘在已经相消的最终速度上；修复仅为回流判定加入机器 epsilon 乘未相消压力冲量尺度，求解后再计入实际 predictor 尺度。没有放宽线性求解残差门、添加有量纲常数下限或裁剪速度，原失败日志保留。[真实缺失回流状态反例](../runs/history/hydro/capillary_boundary_20260911/reject_backflow/run.log)仍在第1步 `coupled-balance-check` 退出3：出口向外速度−0.049999999999999996、旧速度0，判定容差5.000331624505635e-10，明确要求 incoming velocity/material data。

用同一最终 exe 对[原接轴运动带电球滴](../runs/history/hydro/capillary_boundary_20260911/rz_original_final/run.log)完成20步、T=.2并退出0；与第55节 `runs/history/hydro/levelset_inflow_20260911/rz_original` 对照，5份 system 输入、4份诊断CSV和3份保存的电荷面元CSV全部逐字节一致。两个回放均写入此前不存在的独立 run.log，没有覆盖历史目录。

**本节交付是接壁共享几何载荷与自由流体动量行的压力平衡，不是完整机械壁反力、接触线 law 或真实喷嘴供液验收。** 两相 Module 的入口字典/时间循环接线属于后续工作，不能计入上述实跑。原 N32 曲率2.24595%>1%的失败继续保留，本批未重跑或改动该门；接壁电学、运动端点电荷通量、AMR、破裂、重初始化及完整时空/物理精度仍待完成。

## 57. 真实两相 Hydro 供液与运动界面接通；体积预算尚不合格（2026-09-11）

**生产 SingleModule/Hydrodynamics 已完成带毛细力的两相 RZ 持续供液20步，不再只在独立测试驱动中接线。** 永久输入为 [hydrodynamics_prescribed_inlet_rz](../cases/hydrodynamics_prescribed_inlet_rz/README.md)，使用原公共网格、圆弧距离初始化、NS/压力恢复、水平集与 SuperLU。它是显式给定运动学条件的理想供液例，不是用户真实喷嘴参数或完整接触线模型的替代。

实现：字典新增 `interfaceBoundary { model prescribedScalar; fixedContacts; inflow; }`，各负向法向速度入口必须显式给定 `phi=offset+slope*(s-origin)`；只有声明接触侧允许初始辅助闭轮廓延伸出域，未声明墙面、轴镜像、真实零根/单链、完整带覆盖等检查保留。`HydroInterfaceBoundary` 适配器负责公共标量边界/载荷声明、实际几何体积和边界相通量；Module 在初始化前注册条件，后续走真实速度—RK2—导数刷新—材料更新时序，避开不支持该拓扑的旧 cut-length/Peskin 参考。接壁电学和运动端点电荷尚未接通，新模式对 coupled App 明确拒绝。

按同一份 N32、dt=.01、20步、半球 R=.5、rho=10/1、mu=.05/.01、sigma=.01、入口速度.05 的输入，连续处理了两处真实断点，未挪唇口、改输入或换成二次测试 φ：

- 初始化在 `r=.5,z∈[0,1/32]` 检出唇口零点与 `z=1.3130800568474189e−5` 的第二个真实 Hermite 根。声明开链现保留同边全部根，按共享节点 DOF 的精确导数符号判定角点属于哪一单元，单元仍只许一个规则图段；同边小段使用原截线求根/GL8，不用错误的弦中点跨段投影。原不带声明的路径及分支/不明确拓扑拒绝保持。
- 第一遍已完成 NS，但从墙点(.53125,0)回溯左移约1 ulp 后产生8.88e−17假法向入流。原因是入口和墙面速度经 native→vertex→双线性跨唇口混合。Hydro 增加持有 native/壁面快照的分段速度回调：首个法向网格层由本段真实法向约束连接到首内层，出口及内域保留原采样。Core 新可选旧/新速度回调只替换插值，原 RK2、时间、标量和历史提交不变；不扩大标量入口、裁轨迹或抹掉真实入口速度。求解后场输出也使用该回调。

[最终构建](../runs/history/hydro/hydro_prescribed_inlet_20260911/build_final.log)成功。[实际供液日志](../runs/history/hydro/hydro_prescribed_inlet_20260911/supply_patch_velocity/run.log)记录20步至T=.2、exit0；最大散度7.29590e−16、完整自由动量/约束行残差3.24983e−16、全域体积流量不平衡1.73472e−17，轴 ur=0。终态 maxSpeed=.057780964276221702、phi最大变化.00896345650756911。几何液体体积从.26179940779932309增至.26894150702634473，界面真实增长。两相净质量流量不再套单相稳态零通量门，原单相门和不可压体积通量门保留。

**体积预算仍失败，不能标成可靠物理验收。** 南侧精确入口量 `Q=.05*pi*.5²`，T=.2累计输入.007853981633974482；同一实际零曲线 GL8 散度定理几何体积增量比它少.000711882406952845，即供液量的9.06397%。四个新诊断列如实报告初值、当前几何体积、累计负相入流和预算缺陷；原 `negative_volume` 弥散口径未改，没有全局 φ 补偿。源码显示首步旧速度只有入口非零而体内未投影，旧/新平均含非散度零启动场；`-Q*dt/2≈-1.96350e−4`只是该因素的量级估算，不是已测误差分解。后续每步仍亏约2.5–2.7e−5，现 level0 的半拉格朗日双线性 φ 点值插值不保体积，连续重构速度也未获局部保体积证明。下一批先修相容 U0/历史，再处理局部界面输运预算，原输入复测后再开启接壁电耦合。

仅使用已有目标：本批 [coupled-pressure 与 transport 定向回归](../runs/history/hydro/hydro_prescribed_inlet_20260911/targeted_regression.log)通过；分段速度接线后的[最终 transport 回归](../runs/history/hydro/hydro_prescribed_inlet_20260911/transport_patch_velocity.log)14.03秒通过，包含真实公共圆弧初始化双根/小段和两步回调输运。新增半球体积对独立 `2*pi*R³/3` 相对差7.6395e−8只检几何，不替代原真实曲率门或动态守恒。首轮链接时 MSBuild 未重载 GLOB 新源、初始化双根失败、回溯失败的日志均保留在同一 build 批次，无参数扫描扩张。

[原接轴运动带电球滴回放](../runs/history/hydro/hydro_prescribed_inlet_20260911/rz_original/run.log)也用最终主程序完成20步、T=.2、exit0；5份输入及7份原CSV与第56节最终基线全部逐字节一致。本批运行均已结束，未覆盖历史证据。

原 N32 环面曲率2.24595%>1%未修改或重跑，仍待解决；机械壁反力/接触线物理、真实喷嘴、电场接壁/电荷移动边界、重初始化、AMR、破裂和完整时空/物理精度验收均未完成，持续目标保持。

## 58. 优先恢复原 electrospray 喷嘴输入与真实三相初始化（2026-09-11）

用户要求停止扩展验证支线，优先跑通 `cases/electrospray`，并明确选择 `geometry.dict/example.jpg` 的喷嘴，而非旧 App 硬编码的 R=.5 圆滴。**本节实际交付是原喷嘴几何/三相/初始细化恢复，尚未接受物理时间步，不能报为电喷雾跑通。** 供液 U0 与标量插值的分步对照暂停；其已有源码改动保留但未完成独立实跑，不计入本节验收。

原样[运行](../runs/history/electrospray/electrospray_original_20260911/baseline/run.log)首先被 `system/electrospray.dict` 与 `example.jpg` 白名单拒绝；重复配置已归档到 case 的 `legacy/`，图片移到 case 根目录，未删除资料。[再次运行](../runs/history/electrospray/electrospray_original_20260911/layout_fixed/run.log)确认缺失 v3 Domain。旧 mesh 的 `[0,1]×[-1,1]`、App 的硬编码圆与图中 `[0,5]×[0,16]` 喷嘴不一致；旧电学又是 manufactured Helmholtz 且使用 SparseLU，不能恢复调用后就称物理电喷雾。

现将原喷嘴9点12线、3相、固定电极壁和移动弯月面迁移到当前 v3 schema；等价共线 spline 改为 straight，原文本归档。保留20×40背景、全域level3、200000容量与3次初始细化预算，通过公共 `QuadtreeMeshRefiner` 提供 target levels。App使用真实公共距离初始化及相标记，不另建喷雾网格或假定整条LS1都是自由界面。输出只含实际两条LS和phase ID，不写伪电场/速度。

直接修复两项公共初始化断点：显式 `LevelSet.signClosure` 使用有向几何线封闭符号路径，辅助线不参与距离/cut/force；声明端点须连续闭合，域墙符号取域内单侧迹，坐标及距离不移动。`markPhases` 原先在每个中间条件就写相别，会将液帽覆盖成真空；现先检查引用/运算和有限性，只以各相最终复合值生成新标记，完整计算后提交，空phase保留旧行为。

[最终构建](../runs/history/electrospray/electrospray_original_20260911/build_geometry_phase.log)exit0；直接在正式 `cases/electrospray` 目录[实跑](../runs/history/electrospray/electrospray_original_20260911/nozzle_geometry_run.log)约3.2秒，得到51200叶单元、51681顶点，电极/液体/真空顶点数7766/2481/41434，并写入 [electrospray_init.dat](../runs/history/manual_cases/electrospray/output/electrospray_init.dat)。随后 `physical-coupling` 明确exit3：**物理步数0，电场/NS/导电运动耦合仍未接通**。本轮未开参数矩阵、未把旧exe测试算作新回归；公共改动对其它生产路径的回归与原曲率门仍需保留后续验收。

首轮输出的全51681节点核对又定位右壁(5,.05)、(5,.35)两点错相：竖直线段的加权横坐标计算将常量5舍入1 ulp，与域内迹查询重合。公共 `segmentRayCrossingParity` 现对竖线直接保留其恒定横坐标，未扩大容差/移动网格/修改距离。再次[构建](../runs/history/electrospray/electrospray_original_20260911/build_vertical_identity.log)exit0，[原样实跑](../runs/history/electrospray/electrospray_original_20260911/nozzle_final_run.log)约3.46秒，电极/液体/真空修正为7768/2481/41432；仍51200单元、物理步0、physical-coupling exit3。原错误输出完整归档 `nozzle_geometry_before_vertical_fix_output`，最终输出回到正式case目录。

最终场原口径全节点核对：两LS非各自零节点的错符号均0、非真实材料界面错相0、非finite0；参考phi0=z−3，phi1在z≤3为r−1，其余为sqrt(r²+(z−3)²)−1。距离误差8.89e−16/5.00e−5，后者符合DAT六位有效数字舍入。LS零节点161/62，并集222（真正材料界面190，另32在ignored横线）；仍为原±1e−8保护值，不算精确零值/曲率/动态验收。此次只核对已运行App的实际场，没有追加算例矩阵。

下一主线直接在该喷嘴接真实电极、固体流体约束、导电界面与时间循环。已询问用户是否沿用旧IGMPP物性、300/0V、45.2Pa压力供液；当前图没有旧工程的内部引出极，已另询问顶部z=16接地或补回引出极。当前units的1m/1s/1V是移植占位，不能静默解释成微喷嘴物理尺度。原硬编码粒子发射、场蒸发和接触线规律不擅自带入公共代码。详情及运行方法见 [case README](../cases/electrospray/README.md)。

后续只读接口核实：当前公共 Poisson 的可复用能力分属互斥入口，不能把参数确认等同于已接通。`FaceSegmentGraph::build(activePredicate)` 与 `CutCellMetrics` 的复合LS/边界ID映射已支持选择未知区域；`PolarizationPoissonFVSolverCellCentered` 的旧EB入口可接分相系数和多电极Dirichlet，但明确拒绝RZ。严格RZ修正通量入口则拒绝cutMetrics/embeddedBC、删去单元及有父子身份的细化网格；本case虽全域level3后叶单元等大，也不满足原level0契约。下一电学实现必须在共用算子中真正接通RZ固体/电极边界与真实细化拓扑，不能通过伪改Nx/Ny、抹掉parent/level或回退到制造解绕过。此次没有启动新物理运行，物性/供液/接地选择仍待用户明确。

## 59. 原喷嘴物理定标与真实初始电场接通（2026-09-11）

用户随后明确确认：沿用旧工程物性、45.2 Pa 压力供液，喷嘴300 V，顶部z=16为0 V收集极。**这一选择已落入正式case，不再等待参数确认。当前主程序已解出实际喷嘴初始电场，但NS/导电时间循环仍未接通，接受物理步数0，退出码仍为3。** 本批不另开参数矩阵，不将t=0电场或旧制造解作为电喷雾成功。

### 59.1 输入与物理尺度

新增纯输入值对象 `Apps/Runtime/IO/ElectrosprayPhysicalConfig`，与网格、场和求解器解耦。`system/electricity.dict` 使用独立必填 `ElectrosprayPhysical` 块，拒绝未知键/重复/缺项/非有限值。液体为phase2：rho=1240 kg/m³、mu=.038 Pa·s、epsr=10、conductivity=1.4 S/m；gas为phase3：rho=1.184、mu=1.849e-5、epsr=1、conductivity=0；sigma=.0452 N/m，初始自由面电荷0。固体phase1为300 V，top边界9为0 V；活动域其余5/8/10/11/12显式电学零法向通量。压力供液5=45.2 Pa，出口9=0 Pa；只将boundary4视为运动弯月面，整条LS1还包含固定固壁3。

物性与尺度来源为旧 `D:/Chen/IGMPP-ILET/include/Base/Base.h:17–28`，气体参数在62/65行；旧 `HydrodynamicsRZ_impl.h:2280` 给定无量纲压力.001，用p0=gamma/L0=45200 Pa得到45.2 Pa。旧2277行的x<1筛选已注释，不照搬整个底边强设cell压力。旧半径尺度为L0=1e-6 m，Phi0=50 V，u0=sqrt(.0452/(1240*1e-6))=6.037517114114527 m/s，T0=1.6563100047570824e-7 s。R=1 µm半球的纯毛细压差90.4 kPa，远大于45.2 Pa，不宣称供液本身能平衡初始弯月面。

几何坐标/圆心/innerPoint与refineBox均显式`um`，App实日志仍为代码域(0,0)--(5,16)、全域level3。旧radius解析器仅stod，保留1代码长度并明确注释；不把1解释为1 m。RZ完整回转量采用Normalization dimension3，不另乘1 m slab。`run.dict`时间数值未改，仍是代码时间，dt=5e-6对应约8.28e-13 s；还没有实际时间推进。原5份输入已备份至 `runs/history/electrospray/electrospray_physical_20260911/config_before/`。

### 59.2 接公共算子，而非重建喷嘴数值内核

`ElectrosprayElectrostatics`仅做App绑定：用实际水平集的cell中心插值和各相最终复合表达式分类，构造真实active `FaceSegmentGraph`；不平均phase ID，不改mesh.Nx/Ny、level、parent或label。固体未知被排除，其实际邻接面仍保留，绑定到固定几何边界2/3的300 V。每个面检整段覆盖，并核对声明固体边界的总长度，缺失、曲/斜固壁或非贴面对齐均明确报错，不用阶梯固体冒充原几何。

共享 `PolarizationPoissonFVSolverCellCentered::solveFittedBoundary` 增加显式active-leaf贴面边界端口，复用原TPFA装配、RZ面积/体积、SuperLU、完整A*x-b残差及连通分量锚定。私有validator按真实active predicate重建公共graph校验身份、覆盖与反向系数；固体Dirichlet使用owner到实际面的距离，不用到固体cell中心的整格距离。输出预初始化；失败不改原phi，成功保留inactive值与标签/存储。不加gauge/safety pin，不回退CG/SparseLU。旧fluxCorrections/sharp端口的限制与验收门保留，未借新端口放宽它们。

液气介电系数当前明确为`cellHarmonic`的距离加权harmonic epsilon；初始rhs/q=0。它是物理初始介电问题的cell材料近似，**不是sharp-P3、导电松弛或真实Maxwell载荷验收**。电导率已读取但没有假装推进电荷。积分电位移输出单位C，不是电流A。

### 59.3 原case实跑证据与下一主线

[最终主程序构建](../runs/history/electrospray/electrospray_physical_20260911/build_initial_electric_final.log)成功；随后直接运行正式 `cases/electrospray`，约5.66秒：[真实初始电场日志](../runs/history/electrospray/electrospray_physical_20260911/nozzle_initial_electric_run.log)。实际51200叶单元、51681顶点，三相cell数7680/2422/41098，43520个电势未知，SuperLU=1、完整相对残差1.3957440758304249e-15。电势cell范围[.57928697449778199,300] V；0 V施加在顶部实际面，不要求其相邻cell中心为0。

电极气侧/液侧向活动域外的积分电位移分别−1.4708285968159063e-14/−1.4054606759980231e-15 C，顶部1.6113746644246534e-14 C，其余域边界0；相对边界闭合误差2.7754985322853464e-12。最大cell净通量7.8370557594587605e-28 C，离散电能2.4170619966336129e-12 J。面有效法向电场最大5.1981444028788611e7 V/m，仅为TPFA面压降/距离，不当作自由界面两侧sharp迹。

权威SI产物为case `output/electrospray_init_electric_cells.csv` 与 `electrospray_init_electric_faces.csv`，几何DAT仍保留代码坐标/LS/phase口径。上批几何输出完整归档 `geometry_only_output_before_SI/`。首轮MSBuild GLOB重载漏新cpp的链接失败和新局部epsilon0与旧全局常量重名的编译失败分别保留 `build_initial_electric.log` / `build_initial_electric_reload.log`；实际运行使用修复后的新exe。

只对本次修改的共享Poisson补做既有[flux contract回归](../runs/history/electrospray/electrospray_physical_20260911/poisson_contract.log)：目标重新构建，0.83秒通过，未改测试/门限、未跑其它参数矩阵。这与实际case通量核对都不是完整物理收敛验证。

最终将units的rho0从未指定时的默认1显式设为canonical `epsilon0*Phi0/L0²=442.70939064 C/m³`，与Normalization保持一致，避免后续电荷换算混尺度。本次q=0电势与显式eps0诊断不受影响；正式case[再次实跑](../runs/history/electrospray/electrospray_physical_20260911/nozzle_initial_electric_final_run.log)约4.43秒，实际日志两处rho0一致，仍physical-coupling exit3。三份场文件与前次逐字节相同：几何DAT SHA256 `F3788869B01103705F1E162B1D5F7D1DB87FD39351A4F67314FBD021A760DFCE`；cell电场CSV `722670DDA25D47BB09D4BF6551E455CAA743ECC57BB2B23ABE2672F2B911D39E`；face CSV `6E311F28C2D04C7F6B3D76EA327CD2DDBF8DD369CD699C750D91F6C1C07B41B3`。最终system完整快照`inputs_final/`，前次电场移存`electric_before_charge_scale_output/`。

全量实际产物独立重算：174080个有向面中86560对内部面逐对精确反消；非有限0；导体cell全部exact300 V；边界2/3/9用`Vowner−E_normal*distance`逐面回推300/300/0 V，误差exact0，HN与轴通量全0。全量面重算离散能量与`0.5*300V*电极流出电位移`相对差4.15801e-12。这里的cell材料近似、自由界面迹和未推进范围保持，不把线性/预算闭合解释为空间精度已收敛。

下一步仍直推原喷嘴：现有checked Hydro仅接受完整level0矩形流体域，内部固体面/active行映射缺失；pressureOutlet又会拒绝真正入流，不能硬套45.2 Pa供液。接线顺序与最小实际验收写入 本计划第 60 节的贴面 Hydro 实现记录（原临时切片计划已退役）：统一真实活跃域描述，补固定壁无穿透/切向无滑移及明确压力入口物性与对流donor，贯穿现有CN/AB2、projection、coupled与诊断，再绑App。不得伪改mesh或另造NS求解器。原曲率门、供液体积预算、运动电荷、时步/空间收敛、AMR等旧未完成项保持，不因本次初始电场成功而标为通过。

## 60. 原喷嘴贴面 NS 压力供液首次完成短程（2026-09-11）

**原喷嘴上的压力驱动 NS 已真实完成 20 步；尚非完整电喷雾。** 用户确认外侧 boundary8（r=5,z=3..16）按0 Pa环境开放，顶部9为0 Pa出口，固体无滑移。正式 `electricity.dict` 已补 `ambientPressure_SI 8 0`、`solidWallModel noSlip` 和显式 `incomingMomentumModel interiorContinuationZeroShear`。底部5保留45.2 Pa压力供液，喷嘴300 V/顶部0 V及旧物性不变。入口/环境口各自声明液体/气体 incoming 物性；正常压力牵引不等于把边界cell压力设成常数，也不指定流量。顶部没有 incoming 状态，真实回流仍报错。

### 60.1 共享代码接线与范围

新增 `HydroFittedDomain` 作为真实活跃cell/native面及边界描述，由已有初始电学的实际graph和几何边界ID绑定；验证固定同层叶拓扑、实际邻接/面积/半控制体、图覆盖及网格身份，保留level/parent/label，不伪造Nx/Ny。固体单元不作流体未知，固壁/轴法向速度精确消元；每个活动cell仍保留连续性方程，压力口提供物理锚定，无任意gauge或pin。

共享 `Hydrodynamics` 新显式 `solveFittedRZ` 入口，保持原CN/AB2 predictor—增量投影—coupled残差修正；不是另一个BE/Stokes求解器。远壁规则RZ行抽出 `HydroNativeRZRow_impl.h`，新旧入口复用原算术；`HydroFittedMomentum` 负责实际半面/半dual全牵引、法向黏性、交叉剪切、no-slip切向通量及独立hoop，`HydroFittedProjection` 只负责同一D/C映射的代数装配。压力入口采用内侧延拓动量迹与零剪切切向迹，面质量通量仍使用求得的边界速度，不施加 ub=ui。没有读取固体伪速度/物性。

三个线性阶段仍统一 checked SuperLU；完整动量残差和逐cell连续性按输入 `linearResidualTolerance=1e-10` 检查后才提交物理场/历史。边界拓扑变更、混用旧历史、可变dt、未定义的压力half/full阶梯剪切连接等明确拒绝。当前不是一般粗细交界AMR端口，固定同层level3并不等于动态AMR通过。

App仅负责字典边界、SI归一化、材料场及循环。rho0=1240、U0=L0/T0，p0=rho0*U0²=45200 Pa，因此45.2 Pa→代码0.001；RZ实际面积已有2πr，流量只乘U0*L0²。开发模式 `physicsMode frozenMaterialHydro` 显式冻结初始rho/mu，不推进水平集、不施加毛细/电力、不更新电荷；默认 `coupled` 不会自动降级，仍在未实现的完整耦合处失败。

### 60.2 实跑及最小回放证据

[构建](../runs/history/electrospray/electrospray_hydro_20260911/build_main_retry.log)成功；首轮已有MSBuild状态文件访问拒绝保留在 `build_main.log`，只读检查确认其它活跃构建位于独立论文目录后以批准权限重试，未终止用户构建、删除状态文件或覆盖失败日志。

[真实喷嘴短程](../runs/history/electrospray/electrospray_hydro_20260911/frozen_material_run/run.log)退出0。它是正式case输入的独立副本：只在副本run.dict选择开发模式、finalTime=1e-4、writeDt=5e-5，原dt=5e-6与其余输入不变，正式case默认模式/终止时间不改。实际51,200个level3叶单元，43,520个活动压力未知、87,012个自由速度（84,654条规则动量行、2,358条边界分片行）。完整物理时长仅1.6563100047570826e-11 s，不能代表稳态、喷射发展或长时稳定。

20步均三阶段SuperLU=1。[逐步诊断](../runs/history/electrospray/electrospray_hydro_20260911/frozen_material_run/output/electrospray_init_frozen_hydro_diagnostics.csv)的最大momentum/projection/coupled线性相对残差分别3.71550e-16、7.86714e-15、3.10240e-15；最大完整动量残差9.69693e-16，最大代码散度3.76213e-22。第20步maxSpeed=4.1101118017602977e-7 m/s；向外Q（m³/s）为供口−5.6627214656092695e-19、外侧5.6297884289064513e-19、顶部3.293303670281255e-21，净量−3.0359639535481072e-34。20步最大 `abs(netQ)/sum(abs(portQ))=6.24749e-16`。这些是不可压体积通量收支，不是移动界面的液体质量/体积验收。

保存step0/10/20的实际SI cell/native-face场；初始几何DAT及两份电场CSV与第59节最终产物SHA256全部一致。三帧全量核对：43520活动/7680固体、全部level3；15360个solid-only native、188固壁法向native及320轴径向native均exact0，7份Hydro CSV无非有限值；边界面位置/类型/覆盖均正确。CSV未单列固壁切向trace，不能只凭零法向输出宣称全部no-slip精度通过。逐面独立RZ流量核对与更多数值见[产物核对记录](../runs/history/electrospray/electrospray_hydro_20260911/artifact_audit.md)。原正式case输出和旧失败产物均未覆盖。

只补做既有[接轴运动带电球滴回放](../runs/history/electrospray/electrospray_hydro_20260911/rz_original_replay/run.log)，同一新exe完成20步至T=.2，exit0；5份system输入、4份诊断CSV及3份电荷面元CSV与第57节 `runs/history/hydro/hydro_prescribed_inlet_20260911/rz_original` 全部逐字节一致。未改参考/门限或扩展参数矩阵。

下一主线仍在本喷嘴：接固定壁3与自由弯月面4的分离运动、同一贴面动量/压力映射上的毛细和电载荷，再接有限电荷与导电推进。不得整条LS1随动、用有效TPFA面电场冒充sharp两侧迹，或把冻结材料模式当完整喷雾。原N32曲率2.24595%>1%、第57节供液几何体积亏损9.06397%、接触线、时空收敛、动态AMR与破裂等未完成项不因此清零，持续目标保持。

## 61. 唇缘钉扎、真实叶布局与自由液面初值分离（2026-09-11）

用户确认：接触线先钉扎在原喷口唇缘，接触角可变，不沿固壁爬移。`electricity.dict` 现显式为 `contactLineModel pinned`、`pinnedContactPoint 6`，仅引用几何ID。未给压力-only旧快照强加接触线要求；完整界面入口则要求该物理选择完整存在。

代码审查确认不能继续以同一C1零线表示固定wall3+自由cap4：非零梯度给唇口唯一切线，会把壁的切线错误地固定给自由液面；regularized H也会沿固定壁生成假界面带。仅过滤固体面或projection候选不足以修复。本批先分离初始源，不用删除guard冒充运动/载荷已连接。

- 公共 `UniformLeafGrid` 提取原Hydro实际leaf lattice校验及label映射；Hydro删除重复装配前校验，保留RZ测度/D/C和原精确identity。曲率Reconstruction与公共/内部Hermite Snapshot构造链新增显式 `CartesianGridScope::UniformLeaves`，默认RootOnly；读取真实leaf尺寸和祖先身份，不伪改mesh.Nx/Ny。Q4、共享DOF缩放及零点/曲率算术保持。
- 公共 `ImmersedGeometry::initializeUnsignedDistanceToLines` 复用从原初始化机械抽出的同一straight/arc/spline原语，只对选定物理线求距离，不增加closure、近零位移、相标记或历史提交。
- `ElectrosprayInterface` 仅作App绑定：明确原复合LS各线fixed/moving所有权，独立freeLevelSet只取cap4；活动节点由实际fitted fluid cells定义；固壁处符号来自液/气单侧trace，inactive payload禁止送全域Snapshot且CSV写空值。原复合LS和三相初始电场均不改。

[构建](../runs/history/electrospray/electrospray_interface_20260911/build_canonical_pin.log)成功；[原输入副本](../runs/history/electrospray/electrospray_interface_20260911/nozzle_interface_final/run.log)约5.97秒，实际160×320、全level3，44,001活动/7,680 inactive节点，free/fixed/pin分别1/1/1，并写出独立 [free-interface CSV](../runs/history/electrospray/electrospray_interface_20260911/nozzle_interface_final/output/electrospray_init_free_interface_initial.csv)。7份system输入与正式case相同，原几何DAT/两份电场CSV与第59/60节逐字节一致。**此后仍physical-coupling exit3，未新增物理时间步，不能称运动界面或电喷雾已通。**

首轮pin失败保留 `nozzle_coupled/run.log`；精确诊断 `nozzle_pin_diagnostic/run.log` 定位Point6=(1,3)而实际网格节点=(1,3.0000000000000004)。现要求声明坐标逐位等于canonical leaf格点，再按已验证lattice identity查真实label1319；无最近点容差、无节点/曲线位移。另保留原浮点距离==0判断，pin定义在声明的canonical DOF上，输出两套坐标。必须区分：hypot(1,2^-51)会舍入成1，故computed distance=0不是actual点的解析零根证书（精确圆方程残差2^-102）；本批不作该错误宣称。

仅追加一个真实4×4→8×8适配契约到原scalar目标，默认拒绝、显式叶布局/点值和原mesh身份不变均通过；[结果](../runs/history/electrospray/electrospray_interface_20260911/scalar_contract_adapt.json)exit0，原64 polynomial cases/410 edge samples及原门保留，未运行 `--measure` 扫描。初版fixture走 `applyRefineBoxes` 未刷新boundary metadata而被正确拒绝，失败保留 `scalar_contract.log`；改用与App一致的真实 `adapt`（含defineTypes），没有删检查。`applyRefineBoxes` 的metadata刷新缺口不在本批修复，保留待办。

公共改动后的[压力供液回放](../runs/history/electrospray/electrospray_interface_20260911/pressure_replay/run.log)与[原RZ回放](../runs/history/electrospray/electrospray_interface_20260911/rz_original_replay/run.log)均20步exit0；分别7份输入/9份CSV及5份输入/7份CSV逐字节等于第60节。它们使用 `build_main_final.log` 后的exe；随后仅App接触点绑定/输出/报错文字改变，未再改公共数值源。首轮缺少显式QuadtreeMesh头的编译失败保留 `build_main.log`，未覆盖任何历史产物。

全量[实际产物核对](../runs/history/electrospray/electrospray_interface_20260911/artifact_audit.md)：活动44,001节点均有有限phi，inactive7,680节点全部空phi；wall3除pin的60点全负、wall2除pin的128点全正，不再沿固定壁取零距离。独立四分之一圆弧距离最大差5.08220e-21 m，属浮点/单位回转量级，不代表端点曲率验收。[最终诊断文字增量构建](../runs/history/electrospray/electrospray_interface_20260911/build_final_diagnostic_text.log)exit0，未为纯诊断文字重复物理回放。

下一处真实接线：活动流体域的自由cap单侧Hermite/实际零链、pin与axis端点投影，再将其载荷接入现有fitted NS。独立开弧distance在端点Voronoi区域非C1，active mask本身不能消除距离锥；不可直接用跨solid或跨固定壁的居中Q4，也不可将initialization success替代曲率、毛细平衡、运动体积、Maxwell及导电推进验收。第60节列出的全部未完成门保持。

## 62. 原喷嘴唇缘单侧距离导数接通；面元/载荷复用真实叶布局（2026-09-11）

本批继续原喷嘴，不改物性、几何、网格、dt或终止时间。**已从原始自由液面距离算出唇缘单侧导数；仍未把整条活动域 Hermite 零链和毛细载荷接入 NS，完整模式物理步数仍为0。** 上一批的初始化分离是实际进展，本批没有以新的测试矩阵代替喷嘴接线。

公共 `LevelSetCurvature::Reconstruction` 新增显式 native-index 二值 sampleMask。只拷贝入选节点值，允许被排除的 payload 为 NaN；默认空mask保持原算术。若原5×5模板不全入选，选择包含查询点的最近完整入选5×5，仍复用同一QR与25项Q4，不外插、不补零、不按计算出的曲率选模板。它只规定导数取样支持，不声称生成一个活动域C1标量。

公共初始距离原语可同时返回每条有向曲线两端指向曲线内部的单位切向。圆弧/直线复用原曲线，spline沿原采样折线；退化端点报错。App仅以初始自由cap切向确定端点内侧半平面，并与实际active节点取交集，再调用公共重构；不在App塞入圆半径或解析曲率，也不把初始切向固定成演化接触角。

[主程序构建](../runs/history/electrospray/electrospray_contact_jet_20260911/build_main.log)成功；[原输入副本实跑](../runs/history/electrospray/electrospray_contact_jet_20260911/nozzle_initial/run.log)约5.92秒。唇缘(1,3)的42,021个入选节点位于cap侧；导数实际只使用选中的25点。新 [pinned_initial_jets.csv](../runs/history/electrospray/electrospray_contact_jet_20260911/nozzle_initial/output/electrospray_init_pinned_initial_jets.csv) 给出法向 `(0.99999999999013989,-4.4407727235295914e-6)`、子午曲率 `1000353.1348116726 m^-1`、RZ曲率 `2000353.1348018125 m^-1`。用初始半球独立参考 `2/R=2e6 m^-1`，该点相对误差约0.0176567%。这只是单侧Q4距离导数诊断，不是整条曲率收敛或共享bicubic端点曲率验收；Q4的rr/zz不能作为bicubic独立节点DOF硬塞进去。

原几何DAT、两份电场CSV和上批free-interface CSV四份产物全部逐字节一致，原7份system输入不变。日志明确 `physical-coupling` / App rc=3；直接PowerShell重定向stderr的外层工具显示exit1，不将其解释为其他物理失败或成功。尚无毛细/电荷/运动新步，正式case的历史输出未覆盖。

`LevelSetInterfacePanels` 构造尾参scope与 `InterfaceLoadGeometry::Options::gridScope` 默认RootOnly，显式UniformLeaves贯通共享布局、Snapshot、Panels与validateSource；删除重复root布局构造，保留native面遍历顺序、原精度门、RZ测度、端点/拓扑/完整带/投影守卫。**它们仍不支持内部固壁的active-fluid载荷，不能只因接受level3就宣称原喷嘴载荷接通。**

仅在原 scalar 回归追加一个单侧多项式契约和同一4→8真实叶网格上的Panels/LoadGeometry接线检查；[最终结果](../runs/history/electrospray/electrospray_contact_jet_20260911/scalar_contract_final.json)exit0，64 polynomial cases、410 edge samples、32拒绝、145163断言，原误差门不变，未运行measure。首个新增布局圆心(.521625,.533875)、R=.23在N8产生同边双根，被原Panels守卫拒绝，保留 `scalar_contract.log`；将这份新布局fixture圆心改为(.5,.5)使其符合原单交点契约，未改生产守卫、旧fixture或喷嘴输入。

共享几何修改后的既有[原RZ带电液滴回放](../runs/history/electrospray/electrospray_contact_jet_20260911/rz_original_replay/run.log)20步、T=.2、exit0，SuperLU；5份输入和7份CSV与第61节逐字节一致。本批不再重复压力-only20步，它的上一批结果仍保留，不冒充本批重跑。

下一实际实现切片是同一Hermite内核的active patch：实际cell mask与每个节点唯一的 `{phi,phi_r,phi_z,phi_rz}`，只读active四角，不以inactive零多项式占位；pin只固定phi=0，导数仍随液面更新。精确边/顶点统一选择实际active incident cell，seed跳过inactive、Newton不得穿固体、segment不能丢弃inactive区间后称完整；原隐藏5×5三阶重构也须明确支持域。然后以实际active边界incidence替代外矩形端点假设，将自由cap的同一cell曲率、GL8与KKT载荷送入现有fitted D/C动量映射。第60/61节未完成的全局曲率、守恒、运动、Maxwell/导电、收敛、AMR等门均保持。

## 63. 原喷嘴活动域 Hermite 与完整自由液面面元接通（2026-09-11）

**原喷嘴已实际形成唇缘→轴线的完整51面元自由液面；毛细/电荷/运动时间循环仍未接通，物理步数仍0。** 本批交付是同一个公共Hermite标量上的active零链，不是只输出Q4曲率或将固体填零后继续计算。

公共 `ActiveHermitePatch` 是native-index的active-cell二值mask与每个节点唯一的 `{value,dx,dy,dxy}`，Snapshot深拷贝，仅验证/读取活动cell四角，缩放一次供相邻cell共用。实际网格身份仍由UniformLeafGrid验证。inactive cell/unused node的直接查询、空活动域、非finite已用DOF、非零轴奇导数均明确拒绝。原默认全域构造/算术保留；restricted域的额外5×5重构三阶导尚无闭合，明确拒绝，不读取缺失节点。常规同cell三阶导仍可用。

公共Geometry构造/查询贯通patch：精确grid边/顶点可选实际active incident cell；不按近零位置跨入活动区。种子只来自active，原物理边根仍只生成一份；Newton trial不在inactive采样，segment的正长度片必须完整位于活动cell，不能略去固体段。Panels对restricted路径使用实际active/inactive native边界识别声明端点，支持原内部唇缘；inactive无零多项式、无假相别。原root默认闭链路径及其guard保持，active路径仍需显式开链声明、全部根/单cell规则图段/节点度数/整条链验收。

### 63.1 实跑定位并修复初始端点相容性

[原始Q4节点实跑](../runs/history/electrospray/electrospray_active_geometry_20260911/nozzle_raw/run.log)约4.89秒、exit3，断在active Hermite候选单元的Cartesian-partial符号检查；实际DOF完整保留。独立核对定位唇缘下左SW `(31,59)`：两个方向的Bernstein偏导范围都跨0，NW/NE的径向偏导下界却均>.998。原pin的 `phi_z=-4.440772723573383e-6` 确实令壁3上紧邻pin下方出现错误正相，实际Hermite边多项式另有向下 `1.15720132877e-13 m` 的零根；不是单纯放宽检查就能接受的场。

初始声明自由曲线本来就满足端点切向 `dphi/ds=0`。现仅t=0将pin梯度投影到该输入曲线法向，保留Q4法向分量和mixed导数，不填解析曲率，也不把初始角当演化接触角条件。raw Q4诊断完全保留。raw/compatible的51681×6项DOF逐字符串核对仅pin1319的phi_z一项由上述负值变为0；所有原phi、其它导数、activity与身份不变。

修正后SW仍不能通过原偏导充分条件，但已有精确 `cellSignExceptVertices` 证明其16个Bernstein控制为14负2零、每边均有负控制、其它三角严格负，故闭cell只可能在NE pin为零。restricted buildSeeds仅在这个更强证书成立时允许纯零角接触；否则仍原偏导门。`cellSigns`仍0、全部边和零根照常收集，Panels仍核角点incidence，不删根、不把零角冒充严格闭cell同号。原全域路径和门不改。

### 63.2 原输入运行结果与保留范围

[最终构建](../runs/history/electrospray/electrospray_active_geometry_20260911/build_contact_compatible.log)成功；[原输入副本运行](../runs/history/electrospray/electrospray_active_geometry_20260911/nozzle_compatible/run.log)约7.69秒，实际43520活动cell、44001活动节点，Hermite 54候选cell/52种子，完整51面元，pin/axis端点保留。随后仍在 `physical-coupling` App rc=3，未接受新物理步。正式case的7份system与该副本逐字节相同，历史正式输出未覆盖。

新 [free_panels_initial.csv](../runs/history/electrospray/electrospray_active_geometry_20260911/nozzle_compatible/output/electrospray_init_free_panels_initial.csv) 的51编号与cell index均唯一、459个数值均finite、面积全正。独立求和RZ面积 `6.283185323952918e-12 m²`，相对原半球 `2*pi*R²` 差 `+2.6695588e-9`。仅在51个面元中点，圆距离最大误差 `6.03603e-14 m`，RZ曲率范围 `[1999657.88034,2000164.56667] m^-1`、相对 `2/R` 最大误差约0.017106%。面积使用GL8，但此曲率统计不是全GL8节点、端点或收敛验收；原N32环面曲率门没有改变。

原几何DAT、两电场CSV、free-interface初值和raw pinned-Q4 CSV共5份与第62节逐字节相同。已有[scalar契约](../runs/history/electrospray/electrospray_active_geometry_20260911/scalar_contract.json)exit0，64 polynomial/410 edge samples保持，54拒绝、145266断言，新增1个8×8 L形active patch契约，包含固体NaN未读、精确壁选格、1ulp域外拒绝和两端active但中间穿solid的segment拒绝；原门限不改，无新目标/measure扫描。[既有原RZ带电滴回放](../runs/history/electrospray/electrospray_active_geometry_20260911/rz_original_replay/run.log)20步、T=.2、exit0、SuperLU，5份输入及7份CSV逐字节等于第62节。

下一步直接同一喷嘴的fitted毛细载荷/NS：公共载荷仍须支持active图，内部端点KKT不能调用仍仅认外矩形的旧boundaryProjectedFrame。端点曲率还应明确来自cap incident cell的同一bicubic，不能使用另一侧canonical Hessian或独立Q4数值替代。运动时须重建一次共享DOF并只固定pin的phi，不能冻结本批t=0切向；还需闭合端点力/壁反力及液体体积预算。导电、电荷、Maxwell、原曲率与体积缺陷、时空收敛、AMR等未完成项保持。完整产物记录见[本批审计](../runs/history/electrospray/electrospray_active_geometry_20260911/artifact_audit.md)。

## 64. 原喷嘴毛细力—贴面 NS 接通，顶部气体回流显式配置（2026-09-11）

**实际原喷嘴完成20个带毛细力的NS步，exit0；完整电喷雾尚未完成。** 本批仅使用一个原喷嘴短程：保持实际level3叶网格、物性、45.2 Pa供口、300 V喷嘴及0 V顶部；数值dt仍5e-6，开发副本终止时间1e-4、输出间隔5e-5。冻结液面/物性并省略电力/电荷输运的范围由显式 `frozenInterfaceCapillary` 模式声明，不能冒充 `coupled` 喷雾轨迹。实际[运行目录](../runs/history/electrospray/electrospray_capillary_20260911/frozen_capillary_open_top/)和[诊断](../runs/history/electrospray/electrospray_capillary_20260911/frozen_capillary_open_top/output/electrospray_init_frozen_capillary_diagnostics.csv)已保存。

### 64.1 公共几何与载荷接线

公共Snapshot新增 `jetInCell` 和 restricted-only `patchBoundaryProjectedFrame`。端点记录由已验证开链首尾提供实际cap incident cell及内向切线符号；验证真实活动边界、权威零点与单元入射方向，曲率使用该cell的同一bicubic Hessian，不换成canonical另一侧、解析值或独立Q4。stationary/KKT与原外矩形端点接口共用选胜算术，原入口契约保留；KKT残差与端点不为零的stationarity残差分别记录。

`HydroFittedCapillary` 位于流体模块，不含算例物性/SI输出/时间循环。它验证snapshot activity与实际fitted压力域一致，使用同一公共Panels识别开链；H采用原公共compact C2正相指示函数。自由内部面的载荷为

`G_f(H)=[C H]_f/V_f`，`C=-D^T`，`f_sigma,f=-sigma*kappa_RZ,f*G_f(H)`。

这里C与V均取实际fitted NS；原动量RHS乘回同一V，故积分载荷恰为 `-sigma*kappa*C(H)`，没有另乘半径/密度，也没有另造矩形压力梯度。底部/外侧/顶部压力口必须整owner cell严格离开载荷带；固壁/轴约束行不加自由载荷。原 `HydroBalancedForce` 抽出最后一步通用面力算术供两种适配器共用。App只负责 `sigma_code=sigma_SI/(p0*L0)=1`、数值半宽及SI输出；两个冻结模式共用原NS时间循环，没有第二套Stokes/NS求解器。

新运行数值项 `capillaryHalfWidthCells=3`，半宽为 `3*max(dr,dz)=.15`代码长度（.15 µm），全宽覆盖6个轴向cell。默认值有明示，副本也显式填写；它不是物理界面厚度、薄层外推或旧厚度验证的通过判据。正式 `run.dict` 与原geometry/mesh保持不变。

### 64.2 原样实跑揭示的两处阻塞

首次[原输入初始化](../runs/history/electrospray/electrospray_capillary_20260911/nozzle_initial/run.log)在毛细载荷几何查询exit3，定位增强后的[日志](../runs/history/electrospray/electrospray_capillary_20260911/nozzle_query_context/run.log)给出native面62120、query=(.828125,3.3999999999999999)。真实点处于活动j67，但除以hy=.05后舍入成68，而canonical z68=3.4000000000000004；旧restricted选格遗漏实际包含点的相邻格。现只在restricted分支检查邻格的真实闭区间包含关系与activity，不移动点、不加容差、不接纳solid点；原unrestricted路径不改。新增最小非二进制间距fixture同时验证下侧1ULP活动点接受、上侧1ULP非活动点拒绝。

修后[首轮毛细NS](../runs/history/electrospray/electrospray_capillary_20260911/frozen_capillary_run/run.log)已完成683面载荷与首步求解，但顶部9回流而无incoming state，原保护明确拒绝，accepted steps=0、诊断CSV只有表头。用户随后确认顶部允许0 Pa环境气体回流，电学仍0 V。正式 `electricity.dict` 新增可选 `outletIncomingPhaseId 3`，绑定已有reservoir与气相密度/黏度/原incoming momentum闭合，不删除回流保护。旧历史输入未改，未声明该项时仍为仅出流PressureOutlet。代码还让拒绝日志包含实际阶段/native面/回流速度及舍入阈值，不把失败候选提交为accepted状态。

### 64.3 实际产物与短程结果

载荷共683个唯一native面，29个唇缘端点投影；H为43520活动cell和7680空inactive payload，317个过渡cell。所有载荷数值finite、dualVolume正、两侧cell均active；SI `f=-sigma*kappa*G(H)` 与 `integrated=f*V` 的独立相对误差均<4.4e-16。实际C(H)与独立A(Hhigh-Hlow)的最大绝对差4.25007e-17代码单位。RZ曲率范围 `[1999657.73766,2000654.24162] m^-1`；不插入2/R。29个端点κm=1000574.0975105042 m^-1，可从其同cell9189的共享节点DOF独立复算，区别于原Q4 κm=1000353.1348116726；端点KKT/phi残差为0，stationarity最大1.5e-7 m未被伪置零。

最大面力密度5.60668377118e11 N/m³；径向/轴向自由内部面载荷积分分别 `-4.59058352224e-7` / `-2.97551096111e-7 N`。径向数值是RZ径向自由度的工作量对应积分，不是完整三维净径向矢量力；这两个和均不含固壁反力或接触线支撑，不能单独作为界面总力验收。

用户确认回流后的[20步运行](../runs/history/electrospray/electrospray_capillary_20260911/frozen_capillary_open_top/run.log)exit0，三个阶段全SuperLU，完整动量相对残差最大1.76445390333e-14，最大代码散度1.51038065471e-15，均保留原1e-10线性/局部连续性检查。末步速度8.21437280429e-4 m/s；供口/外侧/顶部向外Q为 `+1.13193961068e-15` / `-1.12536564959e-15` / `-6.57396109448e-18 m³/s`，净Q=-1.89197038279e-31 m³/s。供口此时实际流出，不应称为供液喷射；初始半球90.4 kPa毛细压差与45.2 Pa供口不平衡，且尚无电力。物理时间仅1.65631e-11 s，不能代替原finalTime20（4百万步）或稳态/喷射时间尺度验收。

原循环诊断的20步总边界Q相对闭合误差最大1.05949179169e-12（step1，分母为三类口Q绝对和）；仅末步保存native面补偿求和得到1.15734631409e-17，两种求和口径分别保留，不以较小末步值替代原20步指标。step0/10/20全部15868个受约束native速度均exact0，包括188固壁normal、320轴normal及15360 solid-only；不等于额外的切向wall trace或壁反力证明。

默认完整模式[实跑](../runs/history/electrospray/electrospray_capillary_20260911/full_mode_open_top/run.log)完成初始电场、活动cap与毛细载荷后仍在未实现的电力/运动/导电循环exit3，accepted物理步0，不自动调用冻结模式。正式历史输出未覆盖。原7份初始几何/电场/自由液面/DOF/panel产物与第63节逐字节相同；顶端回流前后新增两份初始载荷/H也相同。完整[机读审计](../runs/history/electrospray/electrospray_capillary_20260911/run_audit.json)给出逐项数值与SHA。

现有scalar目标[契约结果](../runs/history/electrospray/electrospray_capillary_20260911/scalar_contract_final.json)exit0，64 polynomial/410 edge samples及原门不改，61拒绝、145306断言；新增的端点和非二进制选格契约均通过，无新target或measure扫描。共同几何/力修改后的一轮[原RZ回放](../runs/history/electrospray/electrospray_capillary_20260911/rz_original_replay/run.log)20步exit0，5份system及7份CSV与第63节逐字节相同；此回放使用active-lookup版本，随后仅喷嘴输入/绑定和fitted失败诊断改变，未再改RZ数值路径。

### 64.4 下一主线与未通过项

继续同一喷嘴接Maxwell电载荷、共享活动域液面运动/材料更新及导电电荷推进，不新增宽参数测试矩阵。本批的phi、初始接触切向和材料冻结仅限显式开发模式；运动时pin只固定phi，不固定接触角。需要量化真实液体体积与电荷收支、端点通量及壁/接触线反力，再谈喷射轨迹。原N32环面曲率2.24595%>1%门、RZ供液体积相对输入9.06397%缺陷、薄层厚度依赖、时空精度、AMR及重初始化/破碎均未由本批解决；不因短程exit0而降门或判整个目标完成。完整证据见[本批审计](../runs/history/electrospray/electrospray_capillary_20260911/artifact_audit.md)。

## 65. 原喷嘴毛细力＋无自由电荷电极化力接通同一 NS（2026-09-11）

**同一原喷嘴已完成20个带毛细力和电极化力的NS步；完整运动/导电电喷雾仍未完成。** 沿用第64节用户确认的顶部0 Pa环境气体回流、0 V收集极，以及45.2 Pa供口、300 V喷嘴和旧工程物性；本批没有改正式7份输入，也没有新增测试矩阵。实际[短程目录](../runs/history/electrospray/electrospray_electric_force_20260911/frozen_ehd_run/)显式选择 `frozenInterfaceEhd`，仅冻结初始液面、材料和电场，不能替代默认 `coupled`。

### 65.1 数值模块与 App 分层

公共电学模块新增 `FittedDielectricPolarization`，消费初始求解器的实际cell电势、harmonic面介电系数、活动图及边界条件。逐个验证真实native邻接、owner/neighbor距离、双向面位移、RZ面积和fitted流体域身份；固体epsilon占位值不参与，真实零面积轴必须显式齐次Neumann，不能以面积乘零隐藏错误Dirichlet。

力使用已有介电模块抽出的同一 `DielectricPolarizationForce::density` 算术：

`f_e,f = -0.5 * [Et_f^2 + Dn_f^2/(epsilon_low*epsilon_high)] * [C(epsilon)]_f/V_f`。

上式为一致电学代码单位，交给Hydro前另乘 `epsilon0*(Phi0/L0)^2/p0=0.48972277725663715`。`Dn`来自实际求解的harmonic FV通量；cell E以实际面电位移/RZ面积重构后除本cell介电常数，`Et`为相邻cell切向分量平均。C和V直接取现有贴面NS；只在active-active自由速度面施加，排除压力口、固壁及轴约束行。均匀介电常数必须给出精确零极化梯度，不把端口上单侧C项当物理极化力，也不另加一遍Maxwell牵引。

这是 **q=0 cellHarmonic极化近似**，不是由台阶电势推得的sharp双侧界面迹，也未证明一般二维/RZ离散电能变分一致。当前52个介电跃迁面与毛细半宽0.15 µm、683面的数值带不同；虽共享压力映射，仍不能宣称两力厚度一致、局部联合静力平衡或原薄层问题已解决。初始无自由电荷守卫仍保留；液体1.4 S/m电导率尚未进入时间推进。

App新增 `ElectrosprayForces` 只负责SI归一化、输出和两力相加；`frozenInterfaceEhd` 与已有压力/毛细模式复用同一CN/AB2 predictor—增量投影—coupled残差修正循环，没有新增NS求解器。旧介电力算术抽取后的既有[局部力契约](../runs/history/electrospray/electrospray_electric_force_20260911/local_force_balance.log) `PASS all`，包括原层/曲线非梯度力及失败不提交字段的检查；没有改变旧参考或门限。

### 65.2 原喷嘴实跑

[最终构建](../runs/history/electrospray/electrospray_electric_force_20260911/build_final.log)成功，当前主程序与本批冻结exe的SHA256相同，均为 `6D2A9823001C0ECA711B51A60726FCB97589ECBCFB8EE398D99B93072FFBED40`。原输入短程副本只选新开发模式并令 `finalTime=1e-4`、`writeDt=5e-5`、显式 `capillaryHalfWidthCells=3`；原dt5e-6、真实level3叶布局和所有物理输入不变。

初始电载荷52面，最大分量 `6.3453911601051111e11 N/m³`；径向/轴向自由面载荷积分 `+3.6155877244285639e-8` / `+2.5036698328391978e-8 N`，不包含固壁反力或接触线支撑，径向也不等于三维净矢量力。实际 `C(epsilon)` 与独立 `A*(epsilon_high-epsilon_low)` 最大绝对差 `4.4408920985006262e-16`（代码单位）；实际Gauss相对残差 `1.3844065069984089e-15`、边界电位移闭合 `2.7753598199027292e-12`。原初始电场输出及其既有诊断不改，新适配器计算乘法顺序略有不同，不能以新残差替换旧电场精度口径。

[20步日志](../runs/history/electrospray/electrospray_electric_force_20260911/frozen_ehd_run/run.log)exit0，三阶段全SuperLU；[逐步诊断](../runs/history/electrospray/electrospray_electric_force_20260911/frozen_ehd_run/output/electrospray_init_frozen_ehd_diagnostics.csv)末步时间 `1.6563100047570826e-11 s`、最大速度 `0.0020935451895587945 m/s`、完整动量相对残差 `1.6884059909753345e-16`、最大代码散度 `1.3252311172764652e-18`。末步供口/外侧/顶部向外Q分别 `+1.0616859852335748e-15` / `-1.0557625965983341e-15` / `-5.923388635240493e-18 m³/s`，净Q `+2.8900083448299527e-31 m³/s`。**供口仍回吐，不能据此声称液体供入或喷射。** 本短程也不是原finalTime20的4百万步、长时稳定性或运动液体体积验收。

全20步而非仅末步：momentum/projection/coupled线性相对残差最大分别 `2.2459906432268753e-15` / `1.0920978683592757e-15` / `2.1419406196191728e-14`；完整动量最大 `1.0218509697881676e-14`、最大代码散度 `8.0161376861435471e-16`。沿用三类压力口绝对Q之和归一化，`abs(netQ)/sum(abs(portQ))` 最大 `1.412651469590172e-12`（step1），不以末步较小数替换。全部诊断有限。

已有[毛细模式原样回放](../runs/history/electrospray/electrospray_electric_force_20260911/capillary_replay/run.log)同样20步exit0，末步速度仍 `0.00082143728042927827 m/s`。默认完整模式使用正式7份输入的[独立副本](../runs/history/electrospray/electrospray_electric_force_20260911/full_mode/run.log)，完成初始电场、活动cap和两力装配后在 `physical-coupling` 返回3、accepted步0；不自动调用冻结模式。未覆盖正式历史输出。本批未重跑RZ，不把第64节旧回放冒充新exe结果。全量力/速度/输入SHA核对见[产物审计](../runs/history/electrospray/electrospray_electric_force_20260911/artifact_audit.md)。

原样回放7份输入/15份旧CSV逐字节等于第64阶段；新EHD短程的9份旧初始产物也不变。step0/10/20每帧的188固壁法向、320轴径向及15360 solid-only native速度均exact0，六份场CSV全部数值finite，不额外宣称壁切向trace或壁反力验收。独立重算52电力面SI公式最大相对差 `6.02e-16`，683面total=cap+electric最大相对差 `3.36e-15`；全活动cell E由实际有向电位移重构的最大相对差 `3.42e-15`。这验证实际输出/单位/映射一致性，不是新的物理收敛结果。

### 65.3 下一实现切片：同一喷嘴的真实液面运动

接下来不是增加冻结力参数扫描，而是把已求得的速度接到独立freeLevelSet，再更新液气材料/电场并接导电电荷。复用已有两时刻速度回调/RK2和同一个active Hermite标量；不能平移固定wall3、从固体伪零值插值、把初始pin切向当演化角度，或只删除RootOnly守卫。

需要同时补齐实际UniformLeaves活动域速度采样与固壁/轴/压力口物理迹；由旧快照反追踪更新活动节点，重建唯一共享 `{phi,phi_r,phi_z,phi_rz}`，先完整验收候选Panels/端点再提交状态。pin固定phi与位置，角度随真实液面更新。现有端口incoming相ID不自动构成全局phi入流闭合，不能随意填充域外phi；运动电荷也必须使用同一轨迹和端点/面积收支，不能继续调用仅支持root且拒绝pin的旧路径。

真实液体体积、电荷收支与端点通量随这个主线实现直接核对；原N32环面曲率门、9.06397%供液体积缺陷、厚度依赖、壁/接触线反力、时空收敛、动态AMR及重初始化/破碎问题仍保留。冻结q=0短程的线性收敛不解除这些门，也不判持续目标完成。

## 66. 原喷嘴运动接线已编译，首步定位到尖角接触表示限制（2026-09-11）

**冻结界面 EHD 仍能完成20步；运动模式尚未接受第1步，不能称为电喷雾跑通。** 本批只用原喷嘴同一短程输入排除实际接线故障，没有新增参数矩阵。顶部继续采用用户确认的0 Pa环境气体回流、0 V收集极；六份物理/网格配置与正式算例不变，开发副本仅选择 `movingDielectric`、显式写默认同值半宽3、缩短终止/输出时间，dt仍5e-6。完整证据保存在 [electrospray_motion_20260911](../runs/history/electrospray/electrospray_motion_20260911/) 及[审计](../runs/history/electrospray/electrospray_motion_20260911/initial_and_failed_attempts_audit.md)。

### 66.1 接线与代码结构

新增公共 `ActiveLevelSetTransport`，从保留的活动域共享DOF构造旧Snapshot，复用原两时刻RK2回溯。实际节点更新采用 `old_phi + phi_old(foot)-phi_old(actual_start)`，导数采用旧jet加 `Q4(delta_phi)`，不在零位移时重拟合初始相容jet。pin只固定phi值，导数允许变化；轴奇导数仍为精确零。新 `HydroFittedVelocitySampler` 留在流体模块，采样实际叶单元/边/节点及固壁、轴、压力口速度迹，不从solid占位值插值，不把实际坐标改成规则乘积格。

`ElectrosprayInterfaceMotion` 负责算例的pin→axis链、固壁、相别及端口验收；App复用既有NS循环，按候选几何→材料/电场→毛细/电力→体积通量的顺序检查后才提交。新 `HydroFittedPhaseMeasure` 负责公共RZ几何体积与纯相压力口通量，App只做SI转换。完整液体体积使用自由面GL8加实际湿壁/域边闭合，不把自由帽体积当全部液体。

回溯的数值速度延拓在 `|phi|<=0.35` 全速、`|phi|>=0.55` 为零，中间复用已有稳定C2函数；这些是代码phi幅值宽度，不是距离或物理界面厚度。旧/新零链的重构支撑整cell必须落在全速带内，压力口整owner cell必须离开截止带及额外halo；没有填充域外phi、移动固定wall3或给固体造零标量。电学仍明确q=0，不推进配置中的1.4 S/m电导率。

### 66.2 同输入四次实跑及修复

| 副本 | 首步NS后的结果 | 处理 |
| --- | --- | --- |
| `moving_run` | 实际节点不满足严格Cartesian乘积坐标假设 | 速度采样保留实际四边形/边中点，规范网格仅用于拓扑搜索；不移动查询点 |
| `moving_actual_mesh` | 极短回溯段的相位分段无可表示内部点 | 新公共 `segmentInsideActiveDomain` 用精确参数区间检查整段活动域覆盖；不调用phi/根/中点，原相位分段保护不改 |
| `moving_transport_domain` | C2直接多项式在端点附近舍入略大于1 | 复用已有尾部稳定的 `regularizedPositiveIndicator`；不clamp、不放宽范围检查 |
| `moving_stable_extension` | 新共享jet在唇缘左下湿液格产生额外零线 | **尚未修复的尖角表示问题**；拒绝候选，不提交运动步 |

四次日志与可执行文件均保留；四次运行的运动及NS两类接受步诊断均只有表头。最终程序SHA256为 `4882AE7F6AFCAD479A3AD7F34F9D339D9092542AE557136DEE976AC116E97963`。原scalar契约通过，64 polynomial/410 edge samples、原门不改，63拒绝、145317断言；新增检查仅覆盖实际活动域短段/穿固体语义。最终程序的[冻结EHD原样回放](../runs/history/electrospray/electrospray_motion_20260911/frozen_ehd_replay/run.log)20步exit0，末步速度仍 `0.0020935451895587945 m/s`，原完整动量与连续性指标不变。没有再跑RZ或薄层矩阵。

冻结回放的7份输入及第65节全部18份CSV逐字节一致；20步三个求解阶段均为SuperLU，完整动量残差最大 `1.0218509697881676e-14`、代码散度最大 `8.016137686143547e-16`，没有改原1e-10门。这仅证明原冻结路径未被本批接线破坏。

运动模式初始化从共享Hermite的实际cell中心判相，恰有cell9194由旧顶点平均法的气相改为液相；中心 `(0.984375,3.175)` 的原圆参考半径平方约0.99961914，确在初始液内。旧7初始产物保持不变，真实新模式step0电场/载荷另存，不能把旧冻结电场冒充运动初态。初始闭合液体体积为 `1.1519173071546478e-17 m³`。由于候选尚未提交，后续材料/电场更新、运动体积/通量误差与动态力闭合都**没有**通过实跑验收。

### 66.3 当前主线：尖角两侧的单侧延拓，而不是更多参数测试

最终候选实际改变1237个phi节点，最大增量 `3.6266288583330919e-11` 代码长度、最大回溯/网格比 `1.1644178690874548e-9`。pin的 `(phi,phi_r,phi_z,phi_rz)` 为 `(0,0.99999999996925437,-7.2553641335589707e-11,2.1719998055625101e-5)`；左下活动cell `(31,59)` 不能取得原偏导或角点单相证书。

这不只是充分条件过严：令唇缘局部坐标 `X=r-1, Z=z-3`，共享非退化C1标量满足 `phi=alpha*X+beta*Z+O(|(X,Z)|²)`。当前alpha>0、beta<0，零线 `X=-(beta/alpha)*Z+O(Z²)` 必在Z<0时进入X<0的湿液格，内壁下方也会出现局部正相。beta来自本次cap侧 `Q4(delta_phi).dy`，不是重新引入初始raw切向；现有日志未保存全部候选增量/NS后速度，不能断言其负号属于真实连续运动还是重构误差。但即使改进Q4，只要允许真实beta<0，一个跨尖角各侧的共同C1 jet仍有同一表示限制。

下一实现切片明确为：在公共几何中区分终止自由面的单侧jet与湿固壁侧的单相延拓，显式声明尖角各侧的连续性范围；自由面保留演化角度，壁侧不产生假气相/额外物理界面。材料、H、体积与零链必须使用一致的这一表示，不能只让Panels挑选想要的一条分支。**不将phi_z强置0、不删根、不跳过湿液格，也不把构造失败当已接受运动。** 这处完成并原样跑过当前喷嘴后，才继续导电电荷和非零q电力接线。原曲率、体积、厚度、时空收敛、壁反力、AMR及破碎未通过项全部保留。

## 67. 原喷嘴运动步：单侧接触延拓与局部坐标（2026-09-11）

本批继续同一 `movingDielectric` 短程，物性、几何、网格、压力/电压、dt及20步终止条件不变，不增加参数矩阵。顶部已按确认配置为允许气体回流的0 Pa环境口，电学仍0 V。**最终原喷嘴运动介电短程接受20步、exit0；完整导电电喷雾与体积守恒验收仍未完成。**

公共 `ActiveHermitePatch::CornerSector` 显式区分湿壁侧延拓与自由面jet：只在三活动格/一非活动格的声明尖角允许一个单侧覆盖，值及沿公共边的多项式保持C0，湿壁格须通过“除声明零角外严格单相”的精确证书。自由面本身的切向导数继续演化；这不是固定接触角或接触线反力模型。App只绑定实际壁面/唇缘及两侧采样支持，通用延拓和验证仍在公共模块。

运动后的自由面确有小于全局坐标ULP的圆弧，不能删掉。网格边交点保留根参数，局部Frame保留cell和局部坐标，Panels的小同边圆弧在局部坐标中构造及GL8积分；全局point仅是舍入后的显示坐标。显式端点相对求根策略用于这一路径，原默认求根/残差门不变。[既有标量契约及同fixture小圆弧检查](../runs/history/electrospray/electrospray_contact_sector_20260911/scalar_contract_panels.log)通过：64 polynomial、410 edge、76拒绝、145380断言，原门限不变；不是整体曲率或EHD验收。

### 67.1 同一首步暴露并修复的两处下游问题

- [带具体区间的失败](../runs/history/electrospray/electrospray_contact_sector_20260911/moving_partition_context/run.log)：实际边端点 `(1,3.0000000000000004)→(1.03125,3)` 的湿区小于全局坐标的可采样间隔。新公共 `partitionSegmentParametric` 复用原限制多项式、根包络及判相，但返回实际子段锚点和根参数；反向保留同一坐标图、不做会抹掉小量的 `1-t`。体积边界项直接积分 `abs(dr*(t1-t0)*(2*r0+dr*(t0+t1)))`，不先把根转成全局坐标再相减。旧 `partitionSegment` 的物理内部点保护不变；网格切分本身不可表示的情形仍拒绝。
- 极小运动造成约0.94%的轴向毛细载荷跳变：水平面9321、9340误将湿壁格内的近零残差点作为自由面投影，取到约 `-7.575e7 m^-1` 的湿壁延拓曲率，而不是约 `+2.000574e6 m^-1` 的cap端点曲率。投影接受现在服从已有精确单相证书：若该格只在声明角点有零值，内部近零点不能作为零面接受；仍由原端点KKT和cap侧jet处理真实端点，不扩大容差或移动点。

批次目录为 [electrospray_contact_sector_20260911](../runs/history/electrospray/electrospray_contact_sector_20260911/)。首次 `moving_nozzle` 在局部根精度契约失败后人为终止，不能记为完成运行；`moving_local_roots` 与 `moving_partition_context` 均在首候选体积检查返回3、完整接受步数0。它们已完成候选材料/电场/力更新，但step1候选载荷文件不是已提交场。原冻结EHD20步基线及历史失败均保留。

### 67.2 最终20步结果与未通过项

[修复后的实际运行](../runs/history/electrospray/electrospray_contact_sector_20260911/moving_parametric_projection/run.log)完成20步并正常退出。最终程序SHA256为 `4DF045B39B715BC84AD66E81F375C0C494A636001F26EF1A13F89C12A5A4DCE0`，保存为批次下 `QUANTUM_parametric_projection.exe`。本批四个运动尝试的7份字典逐字节一致；除开发模式/短程时长所在run.dict外，6份物理输入均与正式算例逐字节一致，正式run.dict未改。

[NS诊断](../runs/history/electrospray/electrospray_contact_sector_20260911/moving_parametric_projection/output/electrospray_init_moving_dielectric_diagnostics.csv)与[运动诊断](../runs/history/electrospray/electrospray_contact_sector_20260911/moving_parametric_projection/output/electrospray_init_interface_motion_diagnostics.csv)各20行，全部数值有限；每步三阶段NS与电场均为SuperLU。完整动量相对残差最大 `8.010856348872897e-15`，代码散度最大 `1.1535468447915813e-15`，末步最大速度 `0.0020940575129193497 m/s`。物理时长仅 `1.6563100047570826e-11 s`，是接线短程，不是成熟喷射轨迹。

两处错误投影均恢复为cap端点，首步RZ曲率 `2000574.113879055 m^-1`；总轴向毛细载荷由 `-2.97551096110843e-7` 变为 `-2.97551096230613e-7 N`，相对变化仅 `4.03e-10`，原约0.94%跳变消失。首步小面元保留8个有限局部GL8节点，弧长 `1.450240223721055e-16 m`，RZ面积 `9.11212806556497e-22 m²`；不能因全局半径显示为pin值便认为该段不存在。step1/10/20均保存完整节点DOF、sector、面元、局部GL8、材料/电场/速度/压力。末步pin仍phi=0，cap侧phi_z为 `-5.920844786677655e-9`，湿壁sector的phi_z为0，表明固定的是位置而非自由面角度。

**体积收支未通过，不能用“小绝对误差”掩盖。** 初始液体体积 `1.1519173071546478e-17 m³`，末体积 `1.1519173060520303e-17 m³`；时间积分液体输入为 `-8.840267908937839e-27 m³`（实际净流出），体积缺陷 `-2.185907512040547e-27 m³`，占该净流出量 **24.7267%**，全程最大比值 **25.6295%**。末步供口向外流量 `+1.0594455440385104e-15 m³/s`，侧面/顶部气体进入，仍不是供液形成喷流。NS离散连续性通过不等于当前插值速度/水平集几何体积守恒通过。

材料判相与电场确实每步重算，但本短程尚无cell跨相，20步 `changed_material_cells=0`；q始终0，1.4 S/m电导率尚未参与推进。[最终现有契约](../runs/history/electrospray/electrospray_contact_sector_20260911/scalar_contract_final.log)通过：原64 polynomial/410 edge及门限不变，77拒绝、145392断言；新增项复用同一L形fixture检查参数湿区、反向小量与湿壁假投影。最终程序的[冻结EHD回放](../runs/history/electrospray/electrospray_contact_sector_20260911/frozen_ehd_replay/run.log)也完成20步exit0，7份字典及全部18份CSV与第65节逐字节一致。

### 67.3 下一主线

不再追加薄层/网格测试矩阵。下一实现分为两个明确工作包：一是让运动采样与同一native面通量的RZ散度相容，修正当前体积收支，而不是全局移动phi抹平账目；二是给实际fitted电学接导电电流、守恒自由电荷状态和非零q电力。应优先评估复用现有cellHarmonic/fitted FV/SuperLU的体积电荷表征，明确与sharp面电荷模型的区别，不能把q=0的现有极化项当完整带电力，也不能重置随材料变化的电荷。旧sharp电荷链尚有UniformLeaves、活动域、局部Frame及钉扎端点接口缺口，不能直接解除限制后调用。

本批只完成原喷嘴的运动介电接线及所述表示修复。原真实曲率、体积守恒、界面厚度/力平衡、时空收敛、接触线/壁反力、导电、电荷、AMR和破碎未通过项全部保留。

## 68. Native RZ通量相容速度已接入，运动体积缺陷仍未消除（2026-09-11）

本批在 `HydroFittedVelocitySampler` 增加显式 `RZFluxCompatible` 重构，保留真实活动叶单元的四个子四边形、原生面平均速度和物理RZ测度。加权Piola映射保留每面的积分通量，局部流函数修正给出共享法向迹、整段无滑移壁和pin零速度；轴列采用 `psi=Phi(xi²,eta)` 的正则延拓。实际面/体积测度还须与原生FV测度在坐标舍入范围内一致。App仅在运动模式的两时刻速度采样调用处选择新重构，仍复用原NS、公共水平集输运、共享Hermite几何与材料/电力路径，没有修改物性、NS离散或时间输入。**面通量相容不等于点值水平集输运已经几何体积守恒。**

[构建日志](../runs/history/electrospray/electrospray_flux_motion_20260911/build.log)与[唯一新增契约](../runs/history/electrospray/electrospray_flux_motion_20260911/flux_contract.log)已保存。目标 `quantum_hydro_fitted_flux_velocity_regression` 使用同一N8网格的完整RZ圆柱和L形活动域：解析流函数端点差生成native Q，逐真实半边GL8积分后核对整面通量；另检查整壁/pin、共享法向、轴正则性、非零离散散度和固体侧拒绝。结果为5319断言、4拒绝，最大面通量相对差 `4.5933553274841072e-16`，独立有限差分RZ散度最大误差 `5.6501185674573406e-11`，未放宽 `1e-8` 门限。它是采样器契约，不是新的NS或液体体积验收。

### 68.1 同输入实际20步

[moving_flux运行](../runs/history/electrospray/electrospray_flux_motion_20260911/moving_flux/run.log)完成20步、exit0。保存程序 `QUANTUM_flux_velocity.exe` 的SHA256为 `6C03627B957A1852F03D72D364ACB432841C7A3F7D9189B512221C6BA2179CE6`。其 `system/` 的 electricity、geometry、mesh、normalization、run、solvers、units七份字典均与第67节 `moving_parametric_projection/system/` 逐字节一致，仍是 `dt=5e-6`、`finalTime=1e-4`、半宽3的同一短程，没有调整压力、电压、材料或网格。

[NS诊断](../runs/history/electrospray/electrospray_flux_motion_20260911/moving_flux/output/electrospray_init_moving_dielectric_diagnostics.csv)和[运动诊断](../runs/history/electrospray/electrospray_flux_motion_20260911/moving_flux/output/electrospray_init_interface_motion_diagnostics.csv)各20行，全部数值有限。三阶段NS各20次、逐步电场20次均为SuperLU。完整动量相对残差最大 `8.010856348872897e-15`、代码散度最大 `1.1535468447915813e-15`，与第67节相同；三端口 `abs(netQ)/sum(abs(portQ))` 最大 `3.856806893791906e-13`。末步速度 `0.0020940575129016624 m/s`，第67节为 `0.0020940575129193497 m/s`。输出存在细小变化：59份同名CSV中38份逐字节相同、21份不同，不能将本次运动运行称为原样回放。

### 68.2 体积与物理边界仍未通过

末液体体积仍为 `1.1519173060520303e-17 m³`，累计液体输入仍为 `-8.840267908937839e-27 m³`（净流出），原定义体积缺陷仍为 `-2.1859075120405466e-27 m³`。相对于累计净流出量的缺陷为 **24.7267111648%**，全程最大 **25.6294971098%**，均未改善。末步供口仍向外 `+1.0594455440385104e-15 m³/s`，外侧/顶部气体进入；不能称为供液喷射，也不能用小绝对体积量或NS连续性通过替代此失败。下一步仍应处理同一水平集输运与几何体积预算，不通过移动全局phi抵消账目。

本次物理时长仅 `1.6563100047570826e-11 s`。材料/电场每步重算，但20步跨相cell数仍为0；q始终0，液体1.4 S/m电导率和自由电荷尚未推进，因此不是完整导电电喷雾跑通。**本批没有新增冻结EHD回放**，第67节的旧回放不计作本exe的新验证。所有既有曲率、体积、厚度/力平衡、收敛、壁/接触线反力、导电电荷、AMR和破碎未通过项继续保留。

## 69. 实际原生面的局部相通量账本已接入，尚未施加几何守恒投影（2026-09-11）

本批新增 `HydroFittedPhaseMeasure::phaseFluxBudget`：沿同一速度采样器的实际 vertex→midpoint→vertex 半边分相，复用解析 `integratedNormalFlux`，保留嵌套局部参数的极小区间；每个负相通量只计算一次，以正负号加入两个活动owner。这里的 `positiveFluxByNativeFace` 指朝+r/+z定向的**负相**通量，不是正相体积。新诊断使用候选新界面与当前速度的瞬时通量，单位m³/s；没有将其作为几何修正量写回。公共 `scalarJetAtCellLocal` 直取同cell的原始value、gradient和Hessian，允许零梯度，为后续共享DOF约束提供接口。

### 69.1 两次失败保留与显式数值契约

三次运行的7份dict全部相同，也与第68节 `moving_flux/system/` 逐字节相同。失败均发生于首步NS之后、完整候选提交之前，接受步数0、rc3；没有将其输出算作已接受步。

| 保存副本及日志 | 结果 | 保存exe SHA256 |
| --- | --- | --- |
| [moving_budget](../runs/history/electrospray/electrospray_local_budget_20260911/moving_budget/run.log) | face9254，`(1,3.0250000000000008)`，极短Cartesian子段无可表示内部点 | `ABAAB007006B3B5BFB1CB1E9E4052690422CF3726E64343B499AB6839B324910` |
| [moving_budget_box](../runs/history/electrospray/electrospray_local_budget_20260911/moving_budget_box/run.log) | face10526，`(.96875,3.2250000000000005)`，含界面根的半边还有ULP网格尾片；整半边符号证书不足 | `0DEA64B16940EF08D7DBCC4DF2CCDFB266FB9CC524A7DDA2AB9F3EFC89F00D11` |
| [moving_budget_certified](../runs/history/electrospray/electrospray_local_budget_20260911/moving_budget_certified/run.log) | 20步、exit0 | `DAE741EC3DB605D6361410E867754C4C259F35A29C3031170A2D5F84A51E2411` |

三个exe分别保存在本批目录的 `QUANTUM_phase_budget.exe`、`QUANTUM_phase_budget_box.exe`、`QUANTUM_phase_budget_certified.exe`。[最终构建](../runs/history/electrospray/electrospray_local_budget_20260911/build_certified_pieces.log)已完成。新增 `certifiedSegmentSign` 以精确活动cell参数区间覆盖和向外舍入的Bernstein包围子盒证明整段严格相别；不能证明时返回0。显式 `CertifyUnresolvedGridPieces` 仅在网格子片两anchor仍不同、整个闭片被证明严格单相时保留原anchor及完整参数(0,1)。含根/触零而无法证明的尾片仍拒绝；默认partition仍拒绝原不可表示中点，原collapsed-gridcrossing门不变，没有删根或平均相别。

[定点契约](../runs/history/electrospray/electrospray_local_budget_20260911/flux_contract_certified.log)通过10182断言、8拒绝，仍只用原两个N8域和原1e-8门：实际co-normal半边GL8、部分/反向/可加/anchored极小通量、线性相区账本、raw jet及显式尾片完整保留。[原scalar契约](../runs/history/electrospray/electrospray_local_budget_20260911/scalar_contract_certified.log)仍为64 polynomial、410 edge samples、145392断言、77拒绝，原归一化最大误差 `1.1510792319313623e-12`，原测量与门限未改。没有新目标或参数矩阵。

### 69.2 最终20步产物核对

第68节全部59份同名CSV与最终运行**59/59逐字节一致**，包括原NS、运动、材料/电场、界面和力输出；仅新增7份相通量CSV。NS与运动诊断仍各20行且全finite，NS三阶段各20次、电场20次SuperLU。原完整动量残差最大 `8.010856348872897e-15`、代码散度最大 `1.1535468447915813e-15`，未因账本接线改变。

[相通量逐步诊断](../runs/history/electrospray/electrospray_local_budget_20260911/moving_budget_certified/output/electrospray_init_phase_flux_diagnostics.csv)20行全finite，每步87520活动面、52个进入根分段的半边。`abs(assemblyClosure)/sumAbsFaceQ` 最大 `6.4659612520122125e-18`，改用 `abs(boundaryQ)` 分母为 `5.138253048239798e-16`；绝对闭合最大 `4.9094168786560375e-31 m³/s`。账本边界与原纯相压力口通量的差，按原门 `abs(diff)/(abs(ledgerQ)+abs(pressureQ))` 归一化最大 `1.0915878165108841e-16`，绝对差最大 `1.9721522630525295e-31 m³/s`。

step1/10/20各有43520行cell ledger、87520行face ledger，其中86560共享面；ID均唯一、数值全finite。从保存的面值按owner±独立重组各cell，并与同帧cell表及逐步诊断核对（以下m³/s，由原CSV的SI数值重算）：

| step | cell表总和 | 面表外边界定向总和 | 逐cell重组最大相对差（分母为该cell面绝对通量和） |
| --- | ---: | ---: | ---: |
| 1 | `5.427140749765718e-17` | `5.4271407497657235e-17` | `2.1029685841267766e-16` |
| 10 | `5.35725345555645e-16` | `5.357253455556448e-16` | `1.605311030569109e-16` |
| 20 | `1.0594455440385106e-15` | `1.0594455440385104e-15` | `1.874877772093479e-16` |

三帧cell总和与对应diagnostic差最大 `9.860761315262648e-32 m³/s`，边界总和差最大 `3.944304526105059e-31 m³/s`；这是逐项SI输出再求和的舍入差，不要求与先求和再缩放逐位一致。此核对只证明账本的共享面/符号/输出一致性。

### 69.3 体积验收仍失败，下一步是同源局部几何量

末步体积缺陷原样为 `-2.1859075120405466e-27 m³`，占累计净液体流出量 **24.7267111648%**，全程最大 **25.6294971098%**。当前仍未完成实际四子四边形控制体的界面重裁、逐cell几何体积或非线性Hermite守恒投影；不能把canonical Cartesian的panel归属无声当成实际控制体。下一切片复用Panels原GL8及同源cell-local子弧求实际控制体体积，再由raw gradient构造共享DOF体积变化约束，不全局shift phi抹账。

本批无新冻结EHD回放。物理时长仍仅 `1.6563100047570826e-11 s`，q=0、1.4 S/m电导率未推进，供口仍回吐；不是导电喷雾验收。旧曲率、体积、厚度/力平衡、收敛、壁/接触线反力、导电电荷、AMR及破碎未通过项全部保留。

## 70. 同源子弧与实际控制体固定边界积分接口（2026-09-11）

本批完成三项公共接口，尚未把几何守恒约束写回水平集：`AxisymmetricBoundaryIntegral::axialVolume` 直接积分保留anchor/参数的直边片，使用真实右侧co-normal，保留微斜竖边及sub-global-ULP片的非零贡献；`HydroFittedPhaseMeasure::boundaryMoments` 汇总实际vertex→midpoint→vertex八边控制体容量及其湿润固定边界项；Panels新增显式共享不可变Snapshot的保留查询模式，以原cell-local单调图提供GL8子弧、正的ds/RZ权重及同cell原始梯度模。当前仅既有定点测试消费这些新入口，喷嘴App尚未调用新boundaryMoments或保留源构造。**固定边界项并非完整负相体积**，还缺自由曲线在实际控制体中的重裁积分。

[配置](../runs/history/electrospray/electrospray_cell_geometry_20260911/configure.log)与[构建](../runs/history/electrospray/electrospray_cell_geometry_20260911/build.log)成功，SuperLU为required target；保存副本 `runs/history/electrospray/electrospray_cell_geometry_20260911/QUANTUM_cell_geometry.exe` 与主程序SHA256均为 `75998217DC868BF3370E0A58CD9FAAB05AE554FE6B25571C8C31C707CD62425A`。[既有flux/geometry契约](../runs/history/electrospray/electrospray_cell_geometry_20260911/flux_geometry_contract.log)通过10681断言、12拒绝，仍只有原两个N8域、原1e-8门，没有新增目标或参数矩阵。新增检查包括独立多项式部分/反向/可加积分、极小参数片、微斜竖边、闭合环柱体积不依赖轴向参考点；同一L域的平面 `phi=z-.6875` 以解析自由平面项补全各cell湿体积，合计湿体积 `5*pi/16`、容量 `5*pi/8`。这些固定平面结果不能推广成所有曲线的overlay验收。

在该固定平面上，保留源构造与旧Panels构造的原端点、测度及GL8位置/法向/曲率/权重逐值相同。新子弧检查使用梯度模为3的同源标量，未假定距离场梯度为1；全弧/部分/反向和sub-global-ULP正长度均保留。借用源的旧构造不开放新查询；局部图本身无法容纳GL8内部点的非零片仍拒绝。GL8子弧积分是原数值求积，不是精确可加原函数。[原scalar契约](../runs/history/electrospray/electrospray_cell_geometry_20260911/scalar_contract.log)仍为145392断言、77拒绝、64 polynomial、410 edge samples，原归一化最大误差 `1.1510792319313623e-12`，原门未改。

**本批没有新喷嘴20步运行或冻结EHD回放。** 第69节最新实跑的体积缺陷 `-2.1859075120405466e-27 m³`、24.73%失败仍保留，不能从接口契约宣称改善；q=0及电导率未推进的边界不变。下一步仍是实际控制体自由曲线重裁和共享DOF体积约束。当前global/cell-local物理embedding并非严格互逆，GL8还有独立求积误差，两者都不能以ULP小量为由忽略；讨论中的控制体外壳界尚未实现，也不作为验收证据。不全局shift phi抹账，不借本批扩大物理结论。

## 71. Tecplot 输出、标准算例和主线进度入口（2026-09-11）

按用户要求，本轮先交付可查看的计算过程，不扩大物理测试矩阵。日常跟进改看[简明进度页](QUANTUM_Hydro_EHD_Progress.md)；运行方法、Tecplot字段和实际结果路径集中在[算例指南](QUANTUM_Hydro_EHD_Case_Guide.md)。参考工程路径已确认为 `D:\Chen\passive_fed`，借鉴其“场云图—液面曲线—过程诊断”的组织，不移植旧算例常量或未接通的电荷/电流结果。

输出职责分开：`Core/IO/TecplotWriter` 只负责 FE/BLOCK 格式，`TecplotMeshOutput` 只适配原网格和已有节点数组；`Modules/Hydrodynamics/IO/HydroFittedOutput` 整理实际活动流体单元和速度重构；`Apps/Runtime/IO/ElectrosprayTecplotOutput` 组合SI流场、电学、界面、面载荷和时间序列。未改 PIC 的旧网格输出函数。普通Hydro/EHD保持原文件名、变量名和节点插值算术，DAT文本改为17位DOUBLE/BLOCK及时间标签；因此不能再用旧FEPOINT文本解析器或逐字节DAT比较来推断物理变化。

新增 [cases/electrospray_moving_dielectric_short_rz](../cases/electrospray_moving_dielectric_short_rz/README.md)，完整复制第69节20步输入，正式 `cases/electrospray` 的完整模式未改。喷嘴输出 `output/tecplot/*_{flow,interface,forces}.dat` 及已完整写帧索引；保留全部原CSV。界面为原面元端点及GL8点的显示折线，曲率/法向取原面元中点值；毛细力显示列为SI总力减极化电力。文件明确这些来源、单位、时间以及未推进导电电荷的限制。

本轮使用SuperLU重新构建，成功日志为 [configure_retry.log](../runs/history/organization/tecplot_workflow_20260911/configure_retry.log)、[build_verified.log](../runs/history/organization/tecplot_workflow_20260911/build_verified.log)。程序归档 `runs/history/organization/tecplot_workflow_20260911/QUANTUM_tecplot.exe`，与实际运行的 `bin/QUANTUM.exe` SHA256均为 `86DDB22D63BD5DB1766FE08181302592F12F3767F5B3675A9E8BBC1F6D60218C`。本批此前的编译名称冲突、源文件注册遗漏、数组适配错误日志均保留，没有覆盖为成功日志。

五个原参数算例实际完成：单相进出口20步、毛细静滴20步、运动界面100步、sharp介电滴20步、原喷嘴20步，退出码均0。全部 **29份新时间序列DAT** 用本机 Tecplot 360 EX 2017 R3 的 `preplot` 实际读取并转成非空PLT，全部退出0；这是数据读取检查，没有声称打开GUI完成视觉验收。普通EHD流场时间0/.01/.02，电力时间0/.009/.019，保留真实受力时刻差异。首次单相启动被PowerShell把OpenMP警告当异常中断，原目录保留；实际完成结果在 `_retry` 目录，不把中断计作求解成功。

[喷嘴本轮日志](../runs/history/organization/tecplot_workflow_20260911/electrospray_moving_dielectric_short_rz/run.log)完成20步，step0/1/10/20各输出三类DAT。与第69节基线相比，**7/7份dict与66/66份原CSV逐字节一致**。末步速度仍 `0.0020940575129016624 m/s`，体积缺陷仍 `−2.1859075120405466e−27 m³`，占累计净液体流出量 **24.73%**；本轮输出工作没有修正这一失败项。时长仍 `1.6563100047570826e−11 s`，不是完整导电喷射过程。

### 71.1 暂存的局部体积修正工作，不作为本轮物理成果

用户转向输出与跟进之前，已新增同源线段根来源查询、实际控制体内自由子弧归属/体积组装、16个Hermite自由度的体积灵敏度，以及独立的SuperLU线性体积修正组件。前两类由 `HydroFittedInterfaceGeometry` / `ActiveHermiteVolumeSensitivity` 提供，后者为 `HydroFittedVolumeCorrection`；App只改为共享保留不可变几何源，**没有调用线性体积修正，更没有将其写回液面**。新线性组件已随本轮主程序编译，但尚未单独验证或接入非线性接受流程。

此前保存的[既有两域定点检查](../runs/history/electrospray/electrospray_cv_partition_20260911/geometry_contract.log)通过11101断言、12拒绝，[原标量检查](../runs/history/electrospray/electrospray_cv_partition_20260911/scalar_contract.log)仍145392断言、77拒绝，原门限不变。增加的解析平面体积/灵敏度核对不能替代实际曲线、非线性体积守恒或喷嘴验收；浮点几何嵌入与GL8求积限制仍在。接下来回到同一原喷嘴，完成局部体积修正及实际更新后检查，再接导电液面电荷，不扩展无关算例矩阵。

## 72. 沿用Plasma输出目录，改善控制台跟进（2026-09-11）

按用户追加要求，Hydro/EHD/喷嘴不再新增 `tecplot/` 格式目录。公共 `Core/IO/OutputLayout.h` 提供按需建目录的统一入口：`fields/`为空间DAT/CSV/TSV，`series/`为逐步残差、流量、体积和电荷诊断，`meta/`为清单、JSON说明及sharp装配明细。没有改Plasma代码或搬动旧结果。喷嘴帧清单的文件引用相对于 `output/` 根目录；原文件名、字段、数值及单位不变。

控制台默认使用中文物理摘要。Hydro/EHD沿用原logEvery并在首个成功步汇总实际求解器；喷嘴增加 `logEvery`（默认1）与 `verboseDiagnostics`（默认false），第一步/最后一步总会报告，逐步CSV不受日志频率影响。摘要包括真实物理时间、最大流速、连续性误差、动量残差、供液/侧面/顶部流量、液体体积缺陷及墙钟耗时。供液以流入为正、侧面和顶部以流出为正；缺陷百分比相对累计净液体输运量，不是液滴总体积损失比例。失败仍报告实际步骤/阶段/底层位置与原因，不以完成文案代替CSV验收。

输出与日志版 [构建](../runs/history/organization/output_layout_console_20260911/build.log)成功，SuperLU未变。原五个标准输入实际完成20/20/100/20/20步、退出0，存档 `QUANTUM_output_layout.exe` SHA256=`901AA66F6BBB55145669C868E94CA6584C78BE9C8B32622E5C6BD6FD235336D6`。四平面入口44份原输出、喷嘴66份原CSV及12份DAT按文件名核对全部逐字节相同；29份DAT全部通过本机Tecplot2017R3 preplot。PLT直接放fields，转换日志放meta，当前目录没有tecplot子目录。各案例结果路径见[使用指南](QUANTUM_Hydro_EHD_Case_Guide.md)。

实跑发现EHD结束文案曾将会被消费清空的 `completedStep_` 临时标志误当累计步数，打印−1；原日志保留。这是新增日志问题，不是物理步数失败。删除错误getter，结束文案改为时间循环完成，实际步数仍由Hydro结束行给出。[修正构建](../runs/history/organization/output_layout_console_20260911/build_console_fix.log)后仅在 [ehd_console_fixed](../runs/history/organization/output_layout_console_20260911/ehd_console_fixed/run.log) 用原输入重跑20步、退出0，30份数据逐字节相同，6DAT再次通过preplot。当前bin与存档 `QUANTUM_console_fixed.exe` SHA256=`408AA3699A58688CBA2BAA0E688881CCC9E600E583EDD1C71A55B8ED8CB61F43`。

同步修正活跃分析入口的路径和读场兼容性：共享Python路径选择/POINT-BLOCK解析器、已有C++读场器及CMake入口保留原变量、采样掩膜和物理门限；新分类目录一旦出现，不能用旧根文件补缺。历史只有meta/ls_diag时仍按旧布局。10个Python文件语法及29个BLOCK/1个历史POINT只读检查通过；已有 `quantum_hydro_regression` [构建](../runs/history/organization/output_layout_console_20260911/build_field_reader.log)与14帧 `verify_fields`、原静滴 `verify_pressure inside_positive` [检查](../runs/history/organization/output_layout_console_20260911/field_reader.log)通过。没有新物理参数矩阵。

### 72.1 主线的实际下一阻断

在用户追加输出要求前，曾将 `HydroFittedInterfaceGeometry::evaluate` 接到原喷嘴初始化。[真实运行](../runs/history/electrospray/electrospray_local_volume_20260911/initial_geometry/run.log)在boundary moment face60857、中心(0.984375,3)遇到 `PolynomialSegmentRootsCandidate` 近多重/临界根符号无法判定，0步、退出3。程序保存为 `runs/history/electrospray/electrospray_local_volume_20260911/QUANTUM_initial_geometry.exe`，SHA256=`1A823EEDE6560C4D3467230BC0D1DC5C517ECB157ACCCBCADB28823447F13C08`；[构建日志](../runs/history/electrospray/electrospray_local_volume_20260911/build_initial.log)保留。临时初始化探针已从标准入口移除，失败没有隐藏或计作物理通过；没有创建新的非线性投影组件，也没有改变原NS/输运算法。

下一步先修这处真实唇缘半面的几何分段，再接每个控制体的体积约束：新旧几何和新旧速度分别形成液体通量，满足 `V_new-V_old+dt/2*(F_old+F_new)=0`，每次修正后重新检查全部单元与几何边界。当前仍有24.73%体积缺陷，仍q=0，导电液面尚未接入原喷嘴运动链；本轮输出和日志改进不解除这些未通过项。

## 73. 真实喷嘴近端分相参数化修复与原样短程复测（2026-09-11）

### 73.1 不是重根：保留末端极短相段

再次在原喷嘴初始化调用 `HydroFittedInterfaceGeometry::evaluate`，取得[实际限制多项式诊断](../runs/history/electrospray/electrospray_volume_roots_20260911/initial_root_context/run.log)。face60857 的后半边从 `(0.984375,3)` 到 `(1,3.0000000000000004)`，标量 cell `(31,60)` 的局部端点为 `(.5,0)` 与 `(1,7.1054273576010019e-15)`。展开系数保留了正的末端值，约 `P(1)=6.31451030642266e-32`；全段导数严格正，根的末端参数距离约 `1-t=4.04128659611050e-30`。原失败来自近 t=1 的相消和参数分辨率，不是真实重根，不能将其吸附为端点或删除极短相段。

公共参数分相接口现在按精确端点值的绝对值选择近端作为参数原点，交换原 anchors 后重新形成同源精确限制多项式。返回顺序仍沿调用方向；保留参数方向，不计算 `1-root`。根隔离器、临界根拒绝和原精度门不变，新增检查禁止将非零端点误认作根。因隔离器只接纳已证明的内部单根，参数片的相号可由精确端点符号逐根翻转，并检查两端符号与根数一致，避免在极小参数上进行不必要的高阶精确乘积而下溢。两端均为零时保留原有间隙采样和失败处理。

交点来源查询另允许宽策略在同一活动闭单元中，对没有物理可表示中点、但实际和局部 anchors 均不塌缩的网格片完成局部多项式隔离，保留两侧和仅触角点的全部记录。严格策略仍拒绝无中点；参数分相对触零网格片的原严格同相证书要求没有放宽。物理坐标、网格身份和钉扎位置没有修改。

### 73.2 已通过的范围与仍然失败的位置

修复后，[实际初始化诊断](../runs/history/electrospray/electrospray_volume_roots_20260911/initial_phase/run.log)已完成所有 `boundaryMoments`，进入自由界面面元的实际控制体归属阶段；该版本随后在无物理中点的交点查询处失败。完善该查询后，[下一份实际诊断](../runs/history/electrospray/electrospray_volume_roots_20260911/initial_inventory/run.log)进一步定位到 `panel=0, sourceCell=9189, face=60857, half=1`：交点来源查询还在原方向隔离近末端根。两份诊断都为 **0 个接受步、退出 3**，不是局部体积或守恒通过。

后续不能只删除 `localPoint` 塌缩检查。当前链路依次以 double 局部点生成切点键、GL8 子弧点和控制体归属位置；位移若已丢失，之后再转 ExactExpansion 也不能恢复。另有更早的映射问题：`gridPoint(physical)-cell` 与 `origin+h*local` 不是精确互逆，此处 z 的两种嵌入相差约 `8.88e-17`，比近根位移大得多。**下一步先统一实际半边到标量局部坐标的映射，再让保留精度的根/切点/求积点贯穿实际控制体积分。** 在此之前不声称完成真实控制体的几何相交计算，也不启用线性或非线性体积修正。

局部体积临时探针保存在对应诊断程序中，未留在标准运行入口；原失败记录完整保留。归档程序：`QUANTUM_phase.exe` SHA256=`D4A6F300EFB6B42CFC2800A33787C1DA7FE0E1998DEE65B0BA142176CB12DF5E`；`QUANTUM_inventory.exe` SHA256=`042933810EB0D40FC8EEE70F2D07AC796C51FFCA110A590952518575689FE9D1`，均位于 `runs/history/electrospray/electrospray_volume_roots_20260911/`。

### 73.3 原门限下的最终验证

[最终构建](../runs/history/electrospray/electrospray_volume_roots_20260911/build_final.log)成功；当前 `bin/QUANTUM.exe` 与归档 `QUANTUM_final.exe` SHA256=`546EE6D44B11875282354654B4E611BFE36CAF9053FC78FCBD362AB1DF3D8BD3`。三项直接相关的[现有回归](../runs/history/electrospray/electrospray_volume_roots_20260911/contracts_final.log)全部通过，[详细输出](../runs/history/electrospray/electrospray_volume_roots_20260911/contracts_final_details.log)保留：

- `hydro_fitted_flux_velocity`：11101 断言、12 拒绝，仍为原两个 N8 域和原门限。
- `ehd_hermite_scalar_contract`：145431 断言、81 拒绝；原 64 个多项式、410 条边采样不变，原最大归一化误差仍为 `1.1510792319313623e-12`。新增同一 N8 网格的精确 `phi=x-1+y^2` 反例，以独立稳定二次公式核对近末端两相与反向参数；另核对无中点小边的四个零端点来源及严格策略拒绝，不修改原门。
- `ehd_hermite_geometry_contract`：31702 断言，109 个非法输入和 15 个不支持情况继续拒绝；原 122 个分相检查通过，最大归一化误差 `1.5660758744363164e-12`。

[原喷嘴标准短程](../runs/history/electrospray/electrospray_volume_roots_20260911/moving_20/run.log)按原七份字典完成 **20 步、退出 0**；7/7 输入逐字节一致。与第72节喷嘴结果相比：13/13 DAT（12份时间帧和1份初始几何）及59/66 CSV逐字节一致；7份相通量CSV存在参数重定向后的浮点末位差异：

| 时刻 | 单元相流量最大绝对差（m³/s） | 面相流量最大绝对差（m³/s） |
|---|---:|---:|
| step1 | 1.5407439555097887e-33 | 1.5407439555097887e-33 |
| step10 | 1.2325951644078309e-32 | 1.2325951644078309e-32 |
| step20 | 2.5422275265911513e-32 | 2.4651903288156619e-32 |

逐步相通量统计最大绝对差 `3.944304526105059e-31 m³/s`，原共享面闭合/压力口检查仍通过。完整流动诊断和界面运动诊断逐字节相同：末步最大流速仍 `0.0020940575129016624 m/s`，体积缺陷仍 `-2.1859075120405466e-27 m³`，占累计净液体流出量 **24.726711%**。这次没有修正该守恒缺陷，也没有接入导电液面；总物理时间仍约 `1.65631e-11 s`。

没有新曲率/厚度/收敛矩阵；DAT 与已通过 preplot 的旧文件相同，本轮未重复进行 Tecplot 转换或 GUI 检查。构建过程中一次旧诊断进程占用 bin 导致链接失败，失败日志 `build_chart.log` 与当次 `initial_chart/` 均保留，不计作修复后的验证；重新链接成功后各版本均归档再运行。最终构建第一次权限审核超时未启动，重试后成功，不视为数值失败。

### 73.4 结构整理的用户确认顺序

用户已确认：先调通原喷嘴算例，再实施[代码与字典整理方案](QUANTUM_Code_And_Dictionary_Organization_Plan.md)，固化正式规范、仓库级 `AGENTS.md` 和可自动执行的依赖/配置检查。当前不提前展开目录与字典重构，也不把介电短程成功当作这一前置任务全部完成。

## 74. 局部守恒几何的一致映射：先完成清分母截线端口（2026-09-11）

### 74.1 本轮修改的准确范围

在公共标量内部实现 [RationalSegmentCoordinate 与 segmentPolynomialCleared](../include/Core/Geometry/Detail/LevelSetHermiteScalar.h)：保留局部坐标的线性分子与严格正的常分母，直接生成 `Q(t)=Dx^3*Dy^3*P(Nx(t)/Dx,Ny(t)/Dy)` 的精确展开系数。`P` 仍由原来一次缩放的 16 个单元 Hermite 自由度定义，包括已声明的 corner sector；不重拟合、不先计算舍入后的局部商、不吸附交点。正比例因子使根与相号不变，但 **Q 的数值和导数不能直接用作物理 phi 或梯度**。

这只是后续一致映射的代数前置，**没有改旧 point query、旧分段、Panels 求积或 App，也没有接入体积修正**。旧 `segmentPolynomialExact` 的实现未动。非法分母、越界局部端点、非活动单元以及不能表示的精确运算尺度仍应明确拒绝。内部代码保留在 `Core/Geometry/Detail`；未向 App 暴露精确展开类型，未开始目录/字典迁移。

### 74.2 已明确的接线约束，尚未实现

名义 `h` 与逐格舍入的坐标边界不能同时视为精确无缝映射。例如原 z 方向 `h=.05`，j=59、60 的相邻名义单元分别有约 `1.8041e-16` 的重叠、`2.6368e-16` 的间隙。它们是浮点表示差别，不是允许修改网格或钉扎位置的理由。

后续局部守恒几何采用显式、同源的共享边界映射：`x=Xi+(Xi+1-Xi)*u`、`z=Zj+(Zj+1-Zj)*v`，其中 X/Z 为现有 canonical 坐标，宽度保留精确差。真实 native 顶点、边中心、控制体和速度自由度不动。该映射保留同一局部 P 及共享边上的连续值，但与旧名义 h 的物理嵌入不是逐位等价；跨格物理法向导数也不能再宣称严格 C1。求积的长度、RZ 权重、单元内梯度与体积灵敏度必须采用同一映射的 Jacobian，不能直接照搬旧 physical jet。

必须同时接好以下消费关系，不能只替换其中一项：

1. 每个旧状态/候选状态生成一份只读的实际半面分段缓存，绑定同一个不可变 Snapshot 和实际网格身份；保留相段参数、根的来源与包围区间。
2. 固定边界体积积分 `boundaryMoments` 与相流量 `phaseFluxBudget` 共用这些湿润相段；流量仍用原 `integratedNormalFlux` 和真实面测度，不额外乘几何映射的 Jacobian。
3. 自由界面弧的切点、GL8 求积、实际控制体归属及灵敏度也来自同源映射；不能先压回 double 局部点再假称恢复精度。
4. 在同一映射下组装 `V_new-V_old+dt/2*(F_old+F_new)`，每次候选修正后重新检查全部活动单元。旧总体积/边界流量及其验收 CSV 继续作为原口径的独立诊断。

若只换自由弧体积而保留旧固定边界分相，初始体积就不是同一闭边界的积分；若只换体积而保留旧相流量，局部修正会对应两条不同的液面。共享面正负装配闭合、压力口收支通过都不能排除这两类不一致。

### 74.3 本轮构建与有界核对

[主程序与两项原几何目标的构建](../runs/history/electrospray/electrospray_cleared_restriction_20260911/build.log)、[原 exact-restriction 目标的构建](../runs/history/electrospray/electrospray_cleared_restriction_20260911/build_exact.log)均退出 0，沿用 SuperLU 和原浮点编译属性。当前 `bin/QUANTUM.exe` 与新归档 `runs/history/electrospray/electrospray_cleared_restriction_20260911/QUANTUM_cleared.exe` SHA256 均为 `4D8D5FCDD9ED9F9AA2695DD7C2EEB3AAF0A1F226ED918A38ADE669BC2C9EC92C`。

- [两项原检查](../runs/history/electrospray/electrospray_cleared_restriction_20260911/legacy_contracts.log)通过：标量仍为 145431 断言、81 拒绝，几何仍为 31702 断言、109 非法输入和 15 不支持情况拒绝；原误差、122 项分相检查与第 73 节一致。
- [精确限制检查](../runs/history/electrospray/electrospray_cleared_restriction_20260911/exact_contract.log)通过：原分类统计仍为 784 个算术核对、23 个根问题、27 个限制问题、26 个非法输入与 2 个不支持情况拒绝。新增使用独立任意长度整数 dyadic 运算的 4 个小型契约例、79 次精确核对、1 个近端根问题和 8 个拒绝；没有新增网格或时间步矩阵。总计 339980 断言；原解析 Q3 最大绝对误差仍为 `2.2204460492503131e-16`。

新增核对包含：7 个不同 dyadic 参数上的齐次 Hermite 因式公式、反向恒等、分母为 1 时与旧接口逐系数等价、实际 canonical 宽度与名义 `.05` 的区别、近端二次根包围区间的独立符号验证，以及同一 N8 唇缘 sector 的单元局部自由度来源。错误分母、越界、非活动单元、不可表示运算尾和非最近舍入环境均保持明确失败。

没有为尚未被 App 消费的端口重复运行原 20 步或生成新 Tecplot 帧。最近一次真实喷嘴结果仍是第 73 节归档：24.726711% 的累计净流出量相对体积缺陷尚未修复，仍为 q=0 介电短程。本轮不作为局部守恒或导电电喷雾验收；代码/字典整理继续遵守“原算例调通后实施”的用户顺序。

## 75. 同映射半面分段接入体积边界项与相流量（2026-09-12）

### 75.1 接线范围

- [LevelSetMappedSegment](../include/Core/Geometry/LevelSetMappedSegment.h) 保留同一个不可变 Snapshot，显式采用第 74 节的共享边界映射。网格切点为正分母的精确有理参数；在其局部参数上进行原有近端根隔离，并保留根包围区间、单元来源及仅触角点的来源。切点排序、相段覆盖与单元关联不先形成舍入后的全局点。共享边上不同单元的非零限制多项式须精确正比例，才共用根。
- [内部实现](../include/Core/Geometry/Detail/LevelSetMappedSegment.h) 留在 Core/Geometry/Detail，公开包装仍编译在原严格浮点几何 TU 内，避免内部展开类型跨编译单元。给积分消费者的结果是分开保留的近端偏移和正长度；它们是末端浮点近似，不是区间误差证书。内部仍保存有理位置，供后续自由弧切分使用。
- [MappedFacePartitions](../include/Modules/Hydrodynamics/HydroFittedPhaseMeasure.h) 为实际活动 native 面的两个半面生成一份只读缓存，检查源、网格和固定流体域身份。新的固定边界体积积分与相流量重载共用这一缓存，不各自重算液面。真实面中心、顶点、速度重构和正负面装配均不改变。
- [轴对称体积边界项](../src/Core/Geometry/AxisymmetricBoundaryIntegral.cpp) 与 [面流量积分](../src/Modules/Hydrodynamics/HydroFittedVelocity.cpp) 消费同一负相区间，不先拼出两个舍入端点再相减。`fromEnd` 只是积分局部坐标反向，不反转物理余法向；通量仍含原 native RZ 测度，不额外乘 Jacobian。面流量系数的原表达式提取为共享函数，旧积分顺序保留。

这实现了同映射的**固定边界体积项与相流量**，还不是完整单元液体体积。自由弧的切分、求积和灵敏度尚未切换到该映射，不能混用旧自由弧后启用守恒修正。没有迁移目录、启用新字典语法或接入导电电荷。

### 75.2 实际喷嘴中发现并处理的尺度问题

第一次[原输入诊断](../runs/history/electrospray/electrospray_mapped_faces_20260912/moving_probe/run.log)完成初始全部活动半面分段与固定边界积分，但在第 1 步 NS 后的候选几何中失败：face9243 half1 从 `(1,2.9750000000000001)` 到 `(1,3.0000000000000004)`，精确乘积尾项超出 binary64 指数范围。退出 3，0 个接受步；该失败记录保留。

局部坐标的分子、分母随后同时按 2 的同一整数次幂精确缩放，使正分母处于单位量级，再清分母。这样不改变坐标比值，避免人为的清分母比例因子取三次幂时造成下溢；没有删除相段、改变根精度判据或允许丢弃运算尾。输入本身超出可表示范围时仍失败。

[缩放后的原样运行](../runs/history/electrospray/electrospray_mapped_faces_20260912/moving_scaled/run.log)完成初始和第 1 步候选的上述诊断，并继续按原物理路径完成 20 步、退出 0。临时诊断只读，不将候选几何校正或新体积定义写回；各数字仍是代码单位的原始装配统计：

| 同映射诊断 | 初始化 | 第 1 步候选 |
|---|---:|---:|
| 根的单元来源记录数（不等于独立交点个数） | 101 | 98 |
| 固定边界体积项装配差 | 3.3306690738754696e-14 | 3.5665914666690488e-14 |
| 相流量装配差 | 未在初始化计算 | 5.082197683525802e-21 |
| 液体边界向外流量 | 未在初始化计算 | 8.989027521061825e-06 |

这些数值不是完整单元守恒通过的判据。原最终体积缺陷仍为 `-2.1859075120405466e-27 m³`，占累计净液体流出量 **24.726711%**；仍为 q=0，物理终时约 `1.65631e-11 s`。

诊断程序与对应 App 源文件单独归档：未缩放 `QUANTUM_probe.exe` SHA256=`09A41DE8472CD116CEEE77FBDB4278020EC5E115DC10E08E371BC57A2EAE6859`；缩放版 `QUANTUM_scaled_probe.exe` SHA256=`615CD134FC425B281813FAB91F5FDC87BE9A7C3231FBDA705DE73D212742D6FA`，均在 `runs/history/electrospray/electrospray_mapped_faces_20260912/`。临时初始/第 1 步探针随后从标准 App 移除，不让未完成的局部守恒流程变成默认运行分支。

### 75.3 原样对照与下一项工作

[逐字节对照](../runs/history/electrospray/electrospray_mapped_faces_20260912/comparison.json)：七份输入字典未变；相对第 73 节原喷嘴记录，15 份 DAT 与 67 份 CSV 全部同名且 SHA256 一致。其中 DAT 为 fields 的 13 份及 meta/ls_diag 的 2 份；CSV 为 fields/series 的 66 份及 meta 帧清单的 1 份。原输出布局不变，没有 tecplot 子目录。没有将控制台或程序文件的哈希当作流场等价证明。

[直接相关检查](../runs/history/electrospray/electrospray_mapped_faces_20260912/contracts_scaled.log)通过：原两个 N8 域的 `hydro_fitted_flux_velocity` 为 19712 断言、17 拒绝，仍使用原 `1e-8` 门；逐半面、面、单元及全域核对新增同源积分与独立解析平面，检查反向余法向、错误反向半面通量及近端小湿段。标量检查为 145474 断言、81 拒绝，原 64 多项式/410 边采样及最大归一化误差不变；新增同一 N8 上的近端两相、无中点四来源和非二进有理网格切点核对，没有新参数矩阵。

**下一步：把这些保留参数的交点用于自由界面弧在真实控制体内的切分、GL8 求积与灵敏度，完成同映射的完整单元体积，再启用候选液面的体积修正。** 在此之前继续保留原 24.73% 失败项，不接入导电液面来掩盖这一问题；结构/字典整理按用户确认顺序留待原算例调通后实施。

### 75.4 最终标准入口与构建

移除临时 App 探针后的[最终构建](../runs/history/electrospray/electrospray_mapped_faces_20260912/build_final.log)退出 0；该批次结束时 `bin/QUANTUM.exe` 与归档 `runs/history/electrospray/electrospray_mapped_faces_20260912/QUANTUM_final.exe` SHA256 均为 `D5E1449F649EC93DFC28EB25FE9BB8CF9AED59BCF2F85DE72E36468CE10BD6A1`。SuperLU 与严格浮点属性保持。

[最终四项检查](../runs/history/electrospray/electrospray_mapped_faces_20260912/contracts_final.log)全部通过：上述 Hydro/标量检查，以及原 `ehd_hermite_geometry_contract` 和 `ehd_exact_restriction_contract`；后两项的统计和误差与第 74 节一致。20 步的实际运行证据属于上面归档的缩放诊断程序；移除纯诊断调用后没有再重复同一 20 步，也未重复转换逐字节相同的 Tecplot 数据。

首轮构建的自动权限审核曾超时，编译未启动；按工具提示重试后成功，不计作数值失败。此后的构建均完成，未遗留运行中的诊断进程。

## 76. 同源自由液面积分接线与切点定位精度（2026-09-12）

### 76.1 本次实现的范围

新增公共 [LevelSetMappedInterface](../include/Core/Geometry/LevelSetMappedInterface.h)，实现仍位于原严格浮点的 `LevelSetHermiteGeometry.cpp` 编译单元。输入为同源标量单元、图线端点身份、保留参数的实际边界切分及凸控制体片；输出自由弧的控制体归属、轴对称体积矩和 16 个原物理节点自由度的体积灵敏度。Core 不认识喷嘴、Hydro 配置或文件输出。

`PanelGraph` 新增端点的 W/S/E/N 边身份及完整边根库存编号：直接从原拓扑库存保留，不按接近的坐标重新匹配，不改变旧面元及其求积值。Hydro 的 [mapped evaluate 重载](../src/Modules/Interface/InterfaceCellMeasure.cpp)将固定边界矩、自由界面矩和节点灵敏度装配到真实原生控制体索引；同一控制体的四片子四边形共用 owner，不把内部片缝当成控制体边界。旧 nominal-h 接口保留，实际控制体构造由两个路径共用。

本次仍是**局部体积修正的几何准备**。没有调用 `HydroFittedVolumeCorrection::solve` 改写候选液面，未推进导电液面电荷，也没有实施目录/字典迁移。

### 76.2 原喷嘴暴露的下溢及修复

四次定位运行分别保存到 `runs/history/electrospray/electrospray_mapped_interface_20260912/`；每份运行的七份字典均与 `cases/electrospray_moving_dielectric_short_rz/system` 的当前原输入 SHA256 相同。仍为 45.2 Pa 压力供液、喷嘴 300 V、顶部接地及原物性，均使用 SuperLU。

1. [首个实跑](../runs/history/electrospray/electrospray_mapped_interface_20260912/moving_probe/run.log)在初始 `panel=0, sourceCell=9189` 报 exact product tail underflow。
2. [增加阶段定位后的复现](../runs/history/electrospray/electrospray_mapped_interface_20260912/moving_context/run.log)明确失败在第一个 GL8 点的基函数/jet 求值，而不是此前的求交或控制体归属判断。局部 r 坐标的分子为 `{−6.2233267924194347e−34, 0.5}`、分母 `0.5`：完整展开张量基函数会生成超出 binary64 指数范围的尾项，但实际求积系数并不为零。
3. 修复仅让**本来就是浮点 GL8 的终端求值**使用因子形式：分别从保留比值求出 `q`、`1−q` 及可能抵消的线性因子；绝不用 `1−rounded(q)`。切点、排序、横向隔根和归属仍保留原 Ratio。横向简单根的导数符号由完整隔根库存的相符号变化确定，浮点导数必须与该证书一致，并通过 `γ128·Σ|term|` 的舍入抵消筛查；该筛查使用 double epsilon，不是严格积分包络，也不是改动原物理/守恒门限。
4. 同时对局部线段两端分母相同的情形做恒等化简，去掉多余的 `D/D`。这是额外的等价运算整理，不将它误报为本例 jet 下溢的直接原因。

### 76.3 原样复跑后的当前失败点

[修复后实跑](../runs/history/electrospray/electrospray_mapped_interface_20260912/moving_factored/run.log)越过前述下溢，在初始 `panel=18, sourceCell=11713` 的极短子弧遇到 GL8 点改变实际控制体归属。[进一步检查切点包络](../runs/history/electrospray/electrospray_mapped_interface_20260912/moving_resolution/run.log)确认：第一个求积点位于左切点自身的定位包络内，因而尚不能确定它在该边界的哪一侧。

该弧的局部 r 坐标端点代表约为 `0.6000002297422129` 和 `0.600000229742216`，长度约 `3.11e−15`。第一个 GL8 点距左代表仅约 `6.17e−17`，而左根包络的上界距该代表约 `4.44e−16`。**代表点可区分，不等于交点已足够精确到可以切出并积分这条小弧。** 新检查明确拒绝落在切点包络内的求积点，没有删除小弧、归并到邻居、移动网格或放松旧阈值。

四次定位运行都退出 **3**，接受步数均为 **0**；本次没有完成初始完整单元体积，也未进入第 1 步候选液面的体积检查。不能将此前第 75 节的 20 步记录写成本次新体积流程通过，更没有新的体积守恒通过结果。此前 **24.73%（相对累计净液体流出量）** 的缺陷及 q=0 限制仍保留。

### 76.4 有界回归、标准入口与下一步

[最终构建](../runs/history/electrospray/electrospray_mapped_interface_20260912/build_final.log)退出 0，[最终四项直接相关检查](../runs/history/electrospray/electrospray_mapped_interface_20260912/contracts_final.log)全部通过。原两个 N8 域、`1e−8` 门及三项原 Hermite 检查保留；没有重新开展网格/曲率/厚度矩阵。

- Hydro 检查现为 **20113 断言、17 拒绝**。新增解析平面检查 359 项：逐单元容量/液体体积、固体 NaN、原节点编号、独立 GL4 对全部 16 列灵敏度，以及常量 phi 扰动的解析体积变化。
- 在同一 N8 网格上新增抛物线 `r=1−z²` 的薄弧检查 42 项。分界 `z=2^-30` 处横向根距右单元边界仅 `2^-57` 个局部单位；薄弧体积矩 `3.383677780714189e−27`、value 列之和 `−5.8516723170686385e−9` 符合解析积分，西侧节点的微小系数非零，反向边查询不改变结果。这不是喷嘴或 NS 守恒验收。

临时 App 调用已移除；`ElectrosprayHydrodynamics.cpp` 与本次修改前保存副本 SHA256 同为 `991A676436BF50DA54107204C15B00BAE87867F8F6D2FE5C0E082325FC8ED413`。保留的公共积分组件和回归没有变成默认的未完成守恒步骤。本批结束时 `bin/QUANTUM.exe` 与[归档标准程序](../runs/history/electrospray/electrospray_mapped_interface_20260912/QUANTUM_final.exe) SHA256 同为 `AEB2361D0BF6F18569A6218A5F341D4C5E1FA82CF37B397752F69B59E810A7DC`。未再次运行移除诊断后的同一 20 步，也未宣称新 Tecplot 流场对照已通过。

诊断程序、临时 App 源码和各阶段 Core 头文件均保留在上述 D: 目录。诊断程序依次为 `QUANTUM_probe.exe`、`QUANTUM_context.exe`、`QUANTUM_factored_probe.exe`、`QUANTUM_resolution_probe.exe`；其 SHA256 分别为 `44FD968ACD33F69A50D1036E0DDC177693FDF09E4469E82F0BCF92AC6E7461EC`、`D2979891A5F134FD40A5D4A4C0D7F5D732FAE284457FDA15C6CCFBD0C550338B`、`0120BAF1D312DCB74B6311331EB058DE17BAFAE613037F93BF926F772143C6B3`、`13844116FDD596330B87876044590FA97A39C4010A3CFE904821C973F3F23E50`。所有本次运行/构建进程已结束，旧结果未覆盖；输出布局仍无 `tecplot/`。

**下一步直接处理共用切点的定位精度。** 在保留参数的局部隔根区间内继续细化，而不是把根转成全局 double 后再求差；将细化结果同时用于固定边界矩、相流量和自由弧切分。不能只给自由弧换一套交点而保留面通量的旧交点。先复现并解决本节短弧，再检查原喷嘴初始体积和第 1 步；随后接完整单元收支和非线性体积修正。导电电荷接线及代码/字典规范固化仍按已确认顺序推进，不以本节局部检查替代完整算例通过。

## 77. 按需交点细化：初始逐单元体积算通，第 1 步边界归属仍失败（2026-09-12）

### 77.1 实现与职责

- `Core/Geometry/Detail/LevelSetMappedCoordinates.h` 提取共用保留坐标运算；`LevelSetMappedRootRefinement.h` 在原简单根包络内重建同一多项式的局部限制并隔根。新包络不得扩大且必须严格缩小，端点根检查非零导数；不转成全局 double、不丢精确尾项、不改隔根门限。
- `LevelSetMappedSegment::refined(rootIndex)` 返回新的不可变查询，同点的闭单元关联与相邻相区间同步更新；正反输入索引保持对应。原查询在成功或失败后仍有效。
- `LevelSetMappedInterface` 只在分段/积分判定根包络不足时请求细化。规范边根以同源多项式的唯一边序号识别；同身份多来源取非空包络交集，不以几何容差合并。没有面来源的图端点可局部细化。不同参数化的浮点代表仍受各自包络约束，不宣称精确水密。
- Hydro 的 `MappedFacePartitions::refineRoot` 成功后才替换半面查询；`evaluateRefined` 丢弃旧结果并重新计算完整边界矩、自由弧和系数。后续相流量必须使用更新后的缓存。256 次面根和 16 次局部端点的限制是工作量上限，耗尽即失败，不是放宽精度。整个准备过程后段失败时，缓存可能保留前面成功的细化，但不返回部分体积，也不接受物理时间步。

### 77.2 原喷嘴同参实跑

输入仍为 `electrospray_moving_dielectric_short_rz` 原七份字典，SuperLU、2 线程、原 20 步时间设置；没有改物性、供液、电极、回流、钉扎或原验收口径。

1. [统一细化尝试](../runs/history/electrospray/electrospray_root_refinement_20260912/moving_nested/run.log)在喷唇极近端点根引发精确尾项下溢。此策略已撤除；当时还有一项旧测试因公开根所属 chart 改变而失败，默认构造恢复后旧检查重新通过。未以忽略下溢或修改旧门限处理。
2. [按需细化运行](../runs/history/electrospray/electrospray_root_refinement_20260912/moving_adaptive/run.log)完成初始全部单元体积：101 个根的单元关联、3 次面根细化、61 段子弧；体积 `11.519173071546486 code³`，与旧全局积分差 `7.1054273576010019e−15 code³`；负体积和超填单元均为 0。此前第 19 个初始面元的失败已越过。
3. 同次运行在第 1 步 NS 后的候选液面 `panel=0, sourceCell=9264` 失败：极短弧的 GL8 点没有实际活动控制体归属。[增强定位复跑](../runs/history/electrospray/electrospray_root_refinement_20260912/moving_ownership/run.log)保存局部/实际保留坐标、横向根包络及所有候选四边形；算术和物理输入未改。两次均退出 3、接受步数 0，未完成第 1 步完整单元体积。

后者位于喷唇右侧，子弧轴向跨度约等于规范 `z=3` 与实际顶点 `z=3.0000000000000004` 的差。将日志17位组件恢复为原 binary64，再用独立 BigInt 有理运算检查：实际半壁从 `(1,3+2^-51)` 到 `(1.015625,3)`，故 `z_wall(r)=3+2^-51−2^-45(r−1)`；整个横向根包络均有 `4.3527e−16 < z_wall−z < 4.3528e−16`。两端叉积严格为负（约 `−6.80112e−18`）；包络宽仅约 `1.0762e−42`，不能改变壁内外关系。其他控制体由径向/轴向范围排除。因此已确认是**规范自由液面与实际固壁的坐标嵌入不一致**，不是横向根精度问题。下一步须修正同一实际喷唇上的几何约束，不能继续盲目细化或删除短弧。**局部体积修正尚未启用，原 24.73% 缺陷、q=0 及导电液面未接入的限制不变。**

### 77.3 回归与标准程序

[最终构建](../runs/history/electrospray/electrospray_root_refinement_20260912/build_final.log)退出 0，[四项最终回归](../runs/history/electrospray/electrospray_root_refinement_20260912/contracts_final.log)全部通过。Hydro 仍 20113 断言/17 拒绝；标量 145555 断言/81 拒绝（新增单根8项和共享边双根73项）；几何仍31702；精确限制340338（新增358项断言，包括4次细化、24项独立整数有理比较、4项拒绝）。原两个 N8 域、所有旧门限和旧用例未改，没有新物理矩阵。

临时 App 探针已归档并移除，标准 `ElectrosprayHydrodynamics.cpp` 恢复为 SHA256 `991A676436BF50DA54107204C15B00BAE87867F8F6D2FE5C0E082325FC8ED413`。本批结束时，[标准程序](../runs/history/electrospray/electrospray_root_refinement_20260912/QUANTUM_final.exe)与 `bin/QUANTUM.exe` SHA256 同为 `24848C6F89C338A46792797AAB4DF82E490A0AF05125D43F222054E08127570F`。诊断程序 `QUANTUM_nested_probe.exe`、`QUANTUM_adaptive_probe.exe`、`QUANTUM_ownership_probe.exe` 独立保留，未覆盖旧结果；本批进程均已结束。未再次重跑标准程序的20步，也未声称生成新的通过验收的 Tecplot 动画。目录/字典规范仍在原喷嘴调通后按确认方案实施。

## 78. 修复公共网格坐标，越过喷唇归属失败（2026-09-12）

### 78.1 复用范围与实际修改

遵守用户确认的“复用旧求解器，恢复耦合链，修复数值问题并通用化”，本批没有换 NS、电场、电荷方程或时间分裂。复用现有 `BackgroundMesh` / `QuadtreeMesh`、编号连接、细化/粗化和场量转移流程，修复其不同位置使用不同浮点加减顺序生成同一坐标的问题。这是**有意改变原生几何末位的数值修复**，不是仅整理接口，也不是把旧结果原样保住。

新增 `Core/Mesh/DyadicGridCoordinates` 小型公共坐标函数：由已有根/父子拓扑恢复整数地址 `(i,j,L)`，中心使用下一层的奇数地址；所有构造点与细分中的 21 类新点共用严格浮点的坐标路径，域端点按输入原值保留。不增加网格实体状态、另一套网格或新字典。精确二进制缩放后同址坐标相同；无法精确缩放的次正规步长、地址超限或相邻点坍缩明确拒绝。名义网格宽度并不承诺等于两个舍入端点相减。

QMC2 **读取**在已有拓扑检查后，追加全部存活单元、角点和边的严格坐标检查，包括非叶父单元、父边和悬挂节点。缓存格式不变；不静默改动旧坐标、cut 数据或历史场，不一致旧缓存要求从 case 重建。写前仍为原 storage 检查，本批不声称其也执行了新增坐标 gate；独立 BackgroundMesh 的 BMC2 读取尚未迁移此检查。相同细化决策下编号/连接流程不变，但不承诺任意依赖几何的后续 AMR 决策逐位相同。

喷嘴接触点绑定追加实际顶点必须与声明点完全一致的检查；不再只把自由液面钉到另一套近似坐标，也不单独挪动唇缘节点。

### 78.2 原喷嘴复跑结果

七份输入字典运行前后均与 `cases/electrospray_moving_dielectric_short_rz/system` 逐字节一致。仍为旧物性、45.2 Pa、300 V / 顶部 0 V、原回流与钉扎边界，SuperLU、原 dt 与 20 步终时，未改变验收口径。[实跑日志](../runs/history/electrospray/electrospray_shared_grid_20260912/moving_probe/run.log)退出 **0**、接受 **20 步**，末步物理时间 `1.6563100047570826e−11 s`。

只在初始化和第 1 步候选液面增加只读局部体积诊断，不修改候选状态、不调用体积修正：

| 项目 | 初始液面 | 第 1 步候选液面 |
|---|---:|---:|
| 同源根的单元关联数 | 118 | 115 |
| 按需面根细化次数 | 0 | 0 |
| 自由积分子弧数 | 51 | 52 |
| 液体体积 / code³ | 11.519173071546488 | 11.519173071518946 |
| 与旧全局积分之差 / code³ | 8.8817841970012523e−15 | 1.7763568394002505e−15 |
| 负体积 / 超填单元数 | 0 / 0 | 0 / 0 |

实际喷唇与钉扎点均为 `(1,3)`，此前首步极短弧无活动归属的失败消失。新增局部体积只检查上表两个状态，不把 20 步完成说成所有后续局部体积均已验证。

终步最大速度 `0.0020940575138516048 m/s`，体积缺陷 `−2.1858828601374292e−27 m³`，累计净液体输入 `−8.8402679089376678e−27 m³`，缺陷占累计净输运量 **24.726432305603126%**。旧基线约 **24.726711%** 保留，末位变化不构成守恒改善。**局部体积修正仍未启用，导电液面电荷仍未推进，不能据此宣称完整电喷雾通过。** [Tecplot DAT](../runs/history/electrospray/electrospray_shared_grid_20260912/moving_probe/output/fields/)保留 step 0/1/10/20 的流场、界面、载荷，仍为 `output/fields` / `series` / `meta`；本批未另做 preplot 或 GUI 动画验收。

### 78.3 构建、直接相关检查与归档

[标准 App / Poisson 构建](../runs/history/electrospray/electrospray_shared_grid_20260912/build_final_app.log)及[五个其余目标构建](../runs/history/electrospray/electrospray_shared_grid_20260912/build_final_tests.log)均退出 0。[最终七项组件检查](../runs/history/electrospray/electrospray_shared_grid_20260912/contracts_final.log)全部通过：网格缓存、Hydro 面通量/速度、两个 Poisson 通量/混合边界、三个原 Hermite 检查；原门限与旧用例未改，没有新的物理收敛矩阵。

缓存检查现为 13 次完整往返、30 次事务式读拒绝、22 次原写前拒绝；新增 `h=0.4` 根网格→L3→粗化→再细化四阶段独立整数坐标对照，总计 398810 条坐标断言。三个单 ULP 扰动均位于原宽松几何容差内，但被新严格读取 gate 拒绝；失败保留接收对象及输入文件。极小步长另验证精确缩放接受、丢位缩放明确拒绝。Hydro 和三项原 Hermite 计数仍分别为 20113 / 145555 / 31702 / 340338。

因修改公共网格，隔离复制并运行三个既有 cut-circle App，均退出 0；原 `RegressionHarness mode record` 未改，**只算运行完成，不算基线精度通过**。PIC 的 EB smoke 缺原输入未运行；已有 Turner 短程在初始化报 `CouplingWeightsCache::build: inconsistent local offsets in parallel build`、退出 2。修改前归档程序在同输入、同两线程环境也以同错退出 2，确认该失败已存在；真正原因未诊断，不能把 PIC 算作通过。详情及原样对照见[记录](../runs/history/electrospray/electrospray_shared_grid_20260912/regression_cases/RESULTS.md)。没有转去修改 Plasma。

只读诊断 App 已归档后移除，标准 App 源码 SHA256 仍为 `991A676436BF50DA54107204C15B00BAE87867F8F6D2FE5C0E082325FC8ED413`。实跑[诊断程序](../runs/history/electrospray/electrospray_shared_grid_20260912/QUANTUM_probe.exe) SHA256 为 `917DE5BD942597EF33BEB22F6364B49073016BB8E0ADA24E666D8A8E37D53AFD`。后续补入极小步长拒绝检查并移除探针，最终标准程序与[归档程序](../runs/history/electrospray/electrospray_shared_grid_20260912/QUANTUM_final.exe) SHA256 为 `348C45CD301060A028F932CCCEE02A4910EAE634A3A169AE907F859626B9B023`；没有再次重跑最终标准程序的同一 20 步。最终涉及源文件存于本批 `source_final/`，各构建/运行进程均已结束，旧算例结果和缓存未覆盖。

下一步沿现有面通量和局部体积组件完成候选液面的守恒修正，然后适配已有电场/表面电荷求解链到真实喷嘴边界。必须区分已有算法复用、耦合适配与必要数值修复；结构/字典规范仍在原算例调通后实施。

## 79. 原喷嘴逐单元体积修正接线（2026-09-12，原样 20 步已完成）

本批把既有实际面通量、同源液体体积与灵敏度、SuperLU 线性步接入运动液面的候选状态。NS、电场、物性、45.2 Pa 供液、300 V / 顶部接地、回流、唇缘钉扎及原 dt / 20 步输入均不变。旧工程未找到可直接恢复的逐单元体积修正，故新外层修正明确归为数值补充，不冒称只是恢复旧接口。详见[本批完整记录](validation/electrospray/QUANTUM_Electrospray_Volume_Correction_20260912.md)。

每轮用同一候选几何重算实际控制体液体体积和相通量，检查 `R_K=(V_new−V_old)+dt/2*(Q_old+Q_new)`；通过时才交给后续材料、电场与受力更新。新增工作迭代准则为局部 L1 残差 / 逐单元输运量不超过 `1e−3`，最大 |R_K| / 单元容量不超过 `1e−10`。两者不是原物理验收门；原全局体积、净流量诊断及约 24.73% 的失败基线全部保留。不裁剪体积、不删除小行或小界面段。

只修节点值的前两次实跑均在第 1 步失败：喷唇新湿单元对允许修正方向极不敏感，线性值修正达到 3.05 code。追加接触点沿湿壁切向的一阶导数后，位置仍固定、角度仍可变、湿壁单侧延拓不改；[第三次实跑](../runs/history/electrospray/electrospray_volume_correction_20260912/moving_contact_gradient/run.log)最大线性值修正降至 `9.16e−11 code`，但局部相对收支仅从 0.429 降至 0.371，未收敛，退出 3、接受 0 步。末次回退在轴端上方单元 12800 的微弧物理坐标分段处拒绝，随后复用已有 cell-local 曲线接口贯通所有显式 sector 面元，保留原根库存与几何检查。

最终内部迭代先降低完整 L1 残差；达到原定工作准则后保持合格，并严格降低梯形目标的越界总量，直至逐单元目标全部有界。只改变内部进展判断；不裁剪实际或目标体积、不放宽最终门限、不增大 12 轮/8 次回退上限。原每轮相通量重算保留；雅可比仍只含体积灵敏度，未声称完整牛顿法。多次中间失败及内部策略修订全部保留在本批记录。

[两阶段版本](../runs/history/electrospray/electrospray_volume_correction_20260912/moving_feasibility_phase/run.log)七份字典哈希原样，实际完成 **20 步、退出 0**，同一终步原全局缺陷从 `−2.1858828601374292e−27` 降至 `+1.1692559034358223e−32 m³`，占累计净液体输运量从 **24.726432%** 降至 **0.000132265%**；本次全过程最大百分比 **0.0516337%**。20 步共 64 次修正、零回退；局部相对收支最大 `8.59315e−5`，最大单元容量相对残差 `1.33321e−13`。这是本次原短程守恒改善，不替代更长时间或原物理收敛验证。

最终[四项直接相关组件检查](../runs/history/electrospray/electrospray_volume_correction_20260912/contracts_stopping_order.log)全部通过，Hydro **20651 断言 / 30 拒绝**，三个原 Hermite 检查仍为 145555 / 31702 / 340338；它们不包含非线性 Step，该部分用上述原喷嘴实际运行验证。新局部收支写入 `output/series` 与 `output/fields`；模块只发布计算数据，App IO 负责逐轮刷新。step 0/1/10/20 共 [12 份 DAT / PLT](../runs/history/electrospray/electrospray_volume_correction_20260912/moving_feasibility_phase/output/fields/)均通过 preplot，完整输出核对见 [output_review_final.log](../runs/history/electrospray/electrospray_volume_correction_20260912/output_review_final.log)。

下一步恢复导电耦合链：先适配已有 `SurfaceChargeField` 的显式叶网格/活动域/局部液面来源，再接固定物质端点、相同运动轨迹、电流源与完整 Maxwell 载荷；继续复用现有积分电荷输运和重映射。导电液面电荷仍未推进，尚不增加导电通过声明；目录/字典整理仍按用户确认顺序在原算例调通后执行。

## 80. 导电液面的共同几何与局部电荷重映射（2026-09-12，接线中）

新增能力在原 `SurfaceChargeField` / `SurfaceChargeTransport` 内，未另建电荷求解器：显式共享 Snapshot / UniformLeaves / 活动域，局部 Frame 来源核对，以及以“面元编号＋面元内面积”表达重叠区间的 open-chain 重映射。沿用原积分 Q 更新、三点重构和守恒门限，旧接口和拒绝保留。其几何和局部面积组件已通过原电荷测试目标；它们尚未作为一条运动电荷链接入 App。

正式原喷嘴实际接线暴露了非 sector 内部 GL8 点与端点的两类坐标转换问题，均按原样初始化证据修复：内部点使用同单元原点逆映并严格保留物理点；端点保留原局部边根参数，不强求新旧显示坐标逐位相同。sector 来源、原 Panel/积分和 Core 根库存检查不变。最终电荷组件检查退出 0，新增非二进制 h 的定向例检出 71 个旧回映失配并验证新接口，未增加物理矩阵。

[最终正式喷嘴初始化](../runs/history/electrospray/electrospray_charge_geometry_20260912/original_nozzle_charge_retained_endpoint/run.log)已完成 51 个电荷面元初始化，逐面元身份和回转面积与同一流体液面一致；七份正式输入哈希不变。之后在 `physical-coupling` 明确退出 **3**，接受 **0 步**，不自动降级为 q=0 短程。运动物质端点重映射、带电电场/导电电流与完整 Maxwell 载荷仍待接通；未放宽旧 bulk 电流残差门，不把非界面电流任意搬到液面。

详细范围、构建/组件/原样运行记录、程序哈希及后续接口见[导电液面接线记录](validation/electrospray/QUANTUM_Electrospray_Conducting_Interface_20260912.md)。正式代码/字典规范仍按用户要求，在原算例完整接通后实施。

### 80.1 共用物质轨迹、局部运动电荷与电场组合接口

后续同批推进已完成 `ActiveLevelSetTransport` 原回溯的持有式物质逆映射、公共局部曲线投影、固定物质接触端/轴端的 Field 运动重映射，以及同一 Poisson 求解器的拟合边界／积分通量修正组合入口。均复用已有算法与 SuperLU，不重建网格、另写 NS 或新增电荷方程。内部面积舍入为端点时明确拒绝，尚未将任意极小剩余子面积列为已支持范围。

原电荷目标包含独立非恒定带符号电荷重映射、原数值轨迹的解析逆、微面元及失败事务，最终 **18.47 s、通过**；原 Poisson 目标新增同网格 L 形活动域的非零源/仿射通量独立检查，**0.70 s、通过**。原参考、门限、检查和拒绝计数保留，新增拒绝另列。它们是组合接口验证，不是导电液面物理验收。

共用轨迹版本在 [moving_trace_first](../runs/history/electrospray/electrospray_material_charge_motion_20260912/moving_trace_first/run.log)原样完成 **20 步、退出 0**，七份字典不变；全部 **86 个 DAT/CSV** 与第 79 节体积修正版逐字节一致，12 个 Tecplot 场已转换为 PLT。后补面积身份拒绝和电场组合接口的最终程序另完成原正式喷嘴初始化（SuperLU、51 个共同液面电荷面元），随后仍在未接完的 coupled 检查处退出 **3**，接受 **0 个导电时间步**；未把两个程序的运行证据混用。

下一步不展开其他物理矩阵：适配已有 sharp/P3 电荷跳跃和分相电流到真实喷嘴活动域、固壁电极及局部曲线，让正式 App 消费上述运动电荷与电场入口，再接完整电应力。当前 q=0 cellHarmonic 短程基线保持，P3 数值选择与接口恢复须明确区分。本批各命令均已结束；详细日志、边界、输出和哈希见同一[接线记录第 6–7 节](validation/electrospray/QUANTUM_Electrospray_Conducting_Interface_20260912.md)。

### 80.2 原喷嘴 P3 带电入口与分相导电电流

新 `FittedSharpDielectricElectrostatics` 入口已经复用同一 Sharp TU 的 P3 跳跃、通量装配和分相电流，使用实际活动叶网格、真实壁面电极及共享局部曲线。App 负责 SI 归一化，不复制电荷或电场算法。取样排除固体，曲率导数明确使用同单元原始 Hermite D3；原完整域重构 D3 和 q=0 cellHarmonic 路径保留，不混称同一离散。

[原喷嘴运行副本](../runs/history/electrospray/electrospray_fitted_sharp_20260912/original_nozzle_p3_first/run.log)已实际求解 43520 个活动单元、74 个修正面，SuperLU 相对残差 `1.3988232658696913e−15`；同次分相电流给出总边界流入约 `2.22411e−5 A`。唯一输入改动为明确选择 `sharpHermiteP3`，其余六份字典哈希一致，物性/电压/压力/网格/时间设置均不改。原 case 输入未覆盖。

初始 51 个共同液面面元仍为 q=0，但电导率已参与 J 的求解。程序在完整 Maxwell 载荷及运动电荷尚未接入循环处退出 **3**，接受 **0 步**，因此下一段是实际时间耦合，不再重复把“电场入口适配”列为未开始。旧 RZ 环面曲率误差 `2.24595% > 1%` 在本轮再次出现，整条首轮结果 JSON 与历史日志一致，保留未通过状态，不调门。具体日志、代码作用域与新增解析验证见[接线记录第 8 节](validation/electrospray/QUANTUM_Electrospray_Conducting_Interface_20260912.md)。

### 80.3 原喷嘴完整 Maxwell 与毛细载荷

现已完成活动域载荷专用局部投影、真实壁面标量梯度与完整横向载荷集合，复用原 P3/Maxwell 及毛细公式，使用同一 Frame 和实际压力梯度。壁面显式采用 ActiveHermiteWallGradient，不把无滑移当作标量梯度零，不宣称近壁离散平衡已经通过。旧 q=0 接口和材料投影门不变。

原喷嘴按上一轮七份字典原样初始化，实际完成 **794 个载荷面、111 个横向新增面、8 个非零固壁来源和 38 个端点投影**；法向电应力峰值约 **12.32 kPa**。初始 q=0，切向应力为零；独立非零带符号电荷的平面/RZ 新检查已通过原门。全部旧电学报告属性（除总断言数）保持，旧 N32 环面曲率 2.24595% 超过 1% 仍导致整体回归失败。

新完整受力和电势/电荷初态均已输出并转换 Tecplot PLT，文件、构建、原样运行和输入/源码哈希见[本轮记录](validation/electrospray/QUANTUM_Electrospray_Maxwell_Loads_20260912.md)。仍在 `physical-coupling` 退出 **3**、接受 **0 个导电步**，不能称新载荷已被 NS 消费。下一段保留一个实际运动循环，拆出纯材料更新后串联 NS → 液面/体积修正 → 原电流与唯一物质轨迹的电荷重映射 → 新 P3 电场与完整受力，按全部耦合检查登记接受步；不另写第二份长驱动。

### 80.4 实际导电时间循环与首步近喷唇投影

已按 80.3 的顺序接入原驱动：纯材料更新、持有式导电时层、电流只积分一次、原物质轨迹重映射、新 P3 电场/电流和完整载荷，以及具名导电 Tecplot 帧和电荷收支。没有第二份 NS 或另一套电荷方程；初态借用作用域与提交后的场引用已核对。当前正式目录/字典迁移仍未开始。

原喷嘴实际运行请求两步；dt、物性、边界、网格、参考和门限不变。第 1 步 NS、液面体积修正、电荷候选和新 P3 电场已执行，但在完整载荷的速度面 60840、坐标 (0.828125,3) 投影失败，实际退出 **3、接受 0 个导电步**。增加诊断后按同七份字典再次运行，复现原八种子零收敛、喷唇端点单侧条件不合格；不能以钉死接触角、放宽门限或吸附端点替代修复。

旧 q=0 分支另按原输入截为单步，实际退出 **0**：41 个完整文件哈希一致、5 个截短序列对应行一致。导电初态四类场的单位/积分关系完成核对并转换 PLT，但未验任何已接受的导电跨步收支。具体日志、源码职责、失败数值、归档程序与下一步公共近边界坐标/投影修复见[本轮记录](validation/electrospray/QUANTUM_Electrospray_Conducting_Time_Loop_20260912.md)。原真实曲率、薄层/局部受力平衡、网格/时间收敛和 AMR 验收均未被本批替代。

### 80.5 近喷唇坐标表示修复与已接受的导电推进

公共边界相对坐标保留 `1−微小偏移` 的内部点身份，同一 Hermite 标量继续供曲率、电场重构、电荷取样和界面力使用。复用原截面根及驻点搜索，原端点条件、残差门和物质运动限制均保留；人工分半截面求根的中点重根明确拒绝。不是另一套几何或 NS 求解器，也不宣称全局最近点完备证明。

原样两步输入已完成 **2 步、退出 0**；两步电荷相对收支残差分别约 `2.70e−13`、`3.23e−13`，均小于原 `1.6e−9` 门。两步完整场、真实定向面电流及跨步电荷检查通过，12 个 Tecplot DAT 已转换 PLT。旧 q=0 单步回归的完整场/首行诊断保持原结果。后补中点拒绝的公共几何检查通过 145755 个断言、101 项拒绝检查，原误差不变。

同一原喷嘴的 20 步延长已完成、退出 0，物性、边界、网格和 dt 不变，仅延长终时及恢复 0/1/10/20 场输出频率。每个接受步电荷、流体及体积历史齐全并核对，四帧场通过导出恒等式检查、16 份 PLT 转换成功。终步 Q 约 `3.36715e−16 C`、边界流入电流 `18.3464 µA`；体积缺陷相对累计净输运量终步 `0.000152247%`。新增 `cases/electrospray_conducting_short_rz`，七份输入与该次运行逐项相同。证据与范围见[本批记录](validation/electrospray/QUANTUM_Electrospray_Boundary_Chart_20260912.md)。

原喷嘴短程耦合链已接通，可以本批为冻结基线进入已确认组织方案的首批整理；正式目录/字典迁移目前尚未开始。不以短程接线完成代替喷射、收敛、薄层/局部平衡或 AMR 验收，后续时间演化仍需继续。
