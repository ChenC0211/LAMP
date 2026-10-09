> **状态：historical/reference（历史资料，仅供追溯）。**
> 本文原始正文保留为历史记录，其中“当前”“下一步”、TODO、旧规则与通过记录均指写作当时；不再自动产生新任务，也不能证明 V2 gate 通过。
> 现行任务唯一入口：[分模块任务清单](../roadmap/05_LAMP_Module_Task_Roadmap.md)。文档权威见[开发总纲](../00_LAMP_Development_Charter.md)，完整入口见[文档索引](../README.md)。
> 旧字典/职责说明可供查询既有实现；V2 新增或修改仍遵守现行六份规范。历史正文中的链接保留原样，目标缺失时不能作为本机已验证证据。
> 本目录中的旧设计、TeX/PDF 和索引同属历史参考；总纲指定的三份冻结方法 PDF 仍以 `docs/` 下的权威文件为准，不因本目录的历史标记而降级。

# QUANTUM 代码组织规范

更新：2026-09-19，沿用 A2 + B1/B2 职责约定。**本文件第 1–4 节是目录、依赖与新增代码约束的规范入口；实际输入语法统一见字典参考。** 当前完成状态只在[模块化主计划](QUANTUM_Modular_Consolidation_Plan.md)第 8 节维护，第 5 节保留各批当时的实施记录，不把其中旧待办视为当前任务。历史缺失文件明确标为未迁入，不作为本机验证证据。[组织方案](QUANTUM_Code_And_Dictionary_Organization_Plan.md)及旧架构 TeX/PDF 保留为历史资料，不重新定义现行规则。

2026-09-29 补充：新算例与数值结果目录遵循下表及第 4 节的时间戳规则。全项目模块、字典、支路和输出的集中整理安排在流体与界面调试完成后，唯一待办为[主计划阶段 P](QUANTUM_Modular_Consolidation_Plan.md#project-usability-cleanup)；该阶段的输出设计尚未实现，不把规划当作现有能力。

## 1. 先复用，再修复与通用化

参数界面的当前边界补充（§2.30）：`InterfaceParametricReconstruction` 在同一个连续材料空间内以能量归一化坐标、直接体积约束分解构造唯一曲线；材料节点细化不等于流体网格 AMR。`HydroSplineVelocity` 在 Hydro 中从累积原生面流量构造 C1 相容速度及同源通量；`Hydrodynamics::splineVelocity` 只从实际 Current/Previous 层发布只读自有快照。`HydroInterfaceTransport` 借用这份快照，Core 沿同一 C1 场推进 Gauss 位置及导数，不在样条结点重启。旧 C0 分片入口保留为显式对照，不能作为 C1 的失败后备。当前支持完整均匀 XY 零法向/零剪切闭域；一般时变、RZ、开口和混合 AMR 仍待统一扩展，禁止在不支持区域拼接旧 C0 场。Interface 仍拥有唯一几何和独立库存，App IO 拥有输出；两档完整冻结往返及真实 NS 两层发布已通过，不能据此宣称毛细耦合或拓扑完成。详见[唯一待办](QUANTUM_Conservative_AMR_Hermite_GALS_TODO.md)。

§2.31 时间扩展：Hydro 的 `timeLinearParametricFlow` 借用布局相同的实际两层 C1 快照；Core 的 `advanceTimeDependent` 微分同一非自治 Gauss 阶段方程，包含阶段物理时刻。Interface 的 `TimeLinearFlux` 保留两个端点通量，原冻结入口显式传入两个相同端点，复用一份时空通量实现；`prepareTimeLinear` 使用真实步长推进既有状态。轨迹二阶、反向通量和短 NS 更新通过，原生顶点停滞处的完整时变交线分区仍未完成；不能按组件通过切换默认或称毛细耦合完成。

复用旧求解器，恢复耦合链，修复有证据的数值问题并通用化；不另建平行的 NS、电场或电荷输运求解路线。共享一套权威网格、水平集及界面身份。算例的喷嘴尺寸、相 ID、物性、压力和电势不得写入公共算法。

纯结构批次不得同时改离散、时间分裂、参考值、掩膜或误差门。不同自由度/边界离散的 EHD 与喷嘴路径不可因目录相同就视为等价。物理失败保持可见；短程运行成功不代替守恒、独立基准和网格/时间收敛验收。

## 2. 目录职责

§2.37 参数毛细力：Core 的 `sampleSurface` 消费唯一参数弧段，返回位置、局部导数、法向与曲率；它不管理物性、力或另一份几何。`HydroSplineVelocityAdjoint` 是原自然样条/累积面流量构造的转置，与前向重构共用布局和自然边界定义。`buildParametricCapillaryLoad` 从同一 InterfaceState 的独立库存及几何形成常曲率压力模态和变曲率界面积分；末端转为原生面力密度，壁面约束行保留反力。来源由 `HydroForceSource` 的参数几何句柄及真实时间绑定，并对照实际质量步的旧/新候选；App 不重标旧力或生成旧标量后备。当前只支持均匀 XY 闭域，不把空间虚功或起始回圆控制称为全离散能量、RZ/AMR 或完整 CCH 验收。

§2.36 参数质量交接：Interface 在同一个 Green 时空积分中计算净量与正向部分，并发布实际正/负累计相通量；`InterfaceState` 持有该记录及旧几何身份、网格绑定和物理时钟。Hydro 的 `HydroConservativeTransport` 新重载直接消费此记录，复用原质量映射；`HydroParametricMaterial` 从同一独立库存生成物性，不借旧标量几何或再次限制通量。公共保留多项式的求值/求根/积分仍由 Core 拥有，Interface 不包含私有展开算术。不以四步无毛细质量交接称为完整 CCH。

§2.35 无电 CCH：相容 `D/C` 仍由 `HydroConservativePressure` 持有，`buildCompatibleScalarGradientLoad` 借用同一算子将应力乘相分数梯度转为真实对偶体积上的力密度；不构建另一套界面。`compatibleNativePressure` 是现有 Hydro 原生保守联立流程的显式装配选择，当前只支持均匀根网格完整 XY/RZ 闭域，接受历史中不能切换。规范压力与速度边界由同一联立装配消元，最终物理压力取体积加权零均值。旧供体输运/黏性/力时层保留，参数界面真实通量及受力来源的正式交接、完整 CCH 中点/AMR 尚未完成；不得因解析曲率平衡通过就切换默认或声称动态毛细验收。

`include/` 与 `src/` 的应用目录对应，只有头文件或只有源文件的条目按实际需要保留，不创建空壳文件。

| 目录 | 职责 |
| --- | --- |
| `Core` | 公共几何、网格、场、单位、通用 IO；不依赖物理模块和 App |
| `Solvers` | 通用算子、线性代数及求解器配置；数值内核不自行读取 case；IO 适配器只解析显式传入的字典，不依赖物理模块和 App |
| `Modules/Hydrodynamics`、`Modules/Electricity` | 方程、离散、边界与状态；接收明确的参数及共享几何/场，不读取 App 配置 |
| `Modules/Interface` | 自由液面状态、相体积/相通量及灵敏度、完整运动与体积修正、材料映射；复用 Core 几何，不管理喷嘴字典编号；当前已接入的几何能力范围见第 5 节 |
| `Modules/*/IO` | 模块字典适配与输出；不承担 App 时钟或耦合步骤 |
| `Apps/Runtime` | App/Module 接口、注册、调度及通用运行支持；不再放具体 `*App` 实现 |
| `Apps/RunCase.cpp` | 应用装配入口：注册具体 App、绑定已选中 case 的输入与启动策略，再交给通用注册器执行；不拥有物理时间循环 |
| `Core/IO/DictionaryUnits`、`Modules/Plasma/IO/DictionaryNormalization` | 前者读取旧通用单位语法，后者解释参考物种与归一化输入、解析物理参考尺度；几何闭合来源和报告时机仍由 App 绑定 |
| `Apps/SingleModule` | 独立模块选择与注册入口；不承担模块数值求解 |
| `Apps/Interface` | 规定输入下的独立界面运行入口；Setup 绑定公共固定域，IO 读输入和写场，主循环只调用同一 Interface 运动服务 |
| `Apps/Hydrodynamics` | 完整域 Hydro 独立运行适配（可选闭域保守 AMR）；`HydroCase` 为 Hydro/EHD 共用的无主循环算例状态，`Coupling` 为 App 耦合阶段契约；严格几何/网格解析归 Hydro/IO，App 只传路径并消费记录，见[域输入归属](validation/organization/QUANTUM_Modular_M6_Hydro_Domain_Input_20260919.md) |
| `Apps/Electrohydrodynamics` | EHD 入口及自身主循环；`Coupling` 组织电场/电荷时层，`IO` 组合跨模块诊断 |
| `Apps/Electrospray` | 喷嘴工作流入口及唯一时间循环 `ElectrosprayTimeLoop` |
| `Apps/Electrospray/Setup` | 把 case 的相、电极、固定壁、压力口、接触线与模块接口绑定 |
| `Apps/Electrospray/Coupling` | 应用单位适配、模块间场交换、候选导电时层及界面运动接线；不另造求解器 |
| `Apps/Electrospray/IO` | 喷嘴场输出、体积/电荷收支和控制台报告 |
| `Apps/Electrospray/Compatibility` | 现有旧字典 reader、输入结构及明确兼容行为；不是新物性模块 |
| `Apps/Verification/Pic` | PIC 验证工作流、验证专用配置及 reader、回归基线处理；Plasma 的 PicConfig 不保存 App 测试开关 |
| `cases/<name>/system` | 可运行输入；不存放另一套求解算法或算例专用运行程序 |
| `cases/<case>/output/<时间戳_标识>/` | 新计算的实际输入/程序快照、运行日志，以及 `output/fields、series、meta`；不覆盖算例模板或其他运行 |
| `runs/` | 保留已有历史结果及引用；新数值结果采用上述 case 输出目录，不批量迁移旧数据 |
| `build` | 编译配置、缓存、库、可执行程序及构建日志；不保存算例物理场、研究结果或文档编译中间文件 |
| `docs` | 顶层为现行规范、指南和主计划；阶段报告按主题放 `docs/validation/`，不混作当前待办 |
| `tests` | 验证脚本、判据与结构检查；数值输出写入独立运行目录 |
| `tools/post` | 从实际运行输出提取可视化或统计；不推进方程、不改参考/门限，不覆盖原始诊断，明确数据来源与完成状态 |

不新增含混的 `Utils/Common/Manager` 汇总目录。公共能力按物理/数学职责归属；跨模块接线留 App，不为方便而让 Electricity 依赖喷嘴配置或让 Core 认识某个液滴。

喷雾的运行组合由 `Setup/ElectrosprayWorkflow::problem` 明确表示：导电、运动介电、固定材料压力驱动、固定界面毛细和固定界面介电；它们组合同一模块服务，不建立另一套 NS。旧 `physicsMode` 名称只由 Compatibility 适配器选择组合，时间循环和导电预检查不再解析字符串。版本化 reader 从互斥的 `Coupling` / `Verification` 选择并持有组合，实际 App 与生效报告消费同一结果，不再由报告重新推测。IO 可输出原历史标签。固定验证显式限制为 perfectDielectric/cellHarmonic，不静默更换模型或开放 AMR/续算；参考[验证入口记录](validation/organization/QUANTUM_Modular_M7_Fixed_Verification_20260919.md)。旧消费者及兼容退出条件仍见[运行组合记录](validation/organization/QUANTUM_Modular_M7_Spray_Workflow_20260919.md)。

公共编译选项由中立的 `quantum_compile_options` 接口目标提供，模块和独立测试不得读取主程序的编译属性。此目标只承载编译设置和外部编译器运行时，不拥有生产源码、项目包含路径或物理/应用库依赖；数值后端仍由 `quantum_linalg` 声明。本机路径留在 `CMakeUserPresets.json`。

独立电学不应为了计算无人消费的流体物性创建 `Hydrodynamics`。旧 `SingleModule/Electricity` 兼容入口已移除该临时对象和无调用者的 `initOnExisting` 重复初始化接口；旧参数仍解析以保持输入兼容。该入口仍使用旧顶点电学及 `LegacyNodalState`/毛细初始化，不把这次清理描述为旧模型全部独立化或替换成正式电学。见[清理记录](validation/organization/QUANTUM_Modular_M4_Legacy_Electric_Initialization_20260919.md)。

## 3. 接口、所有权与 IO

拟合域的毛细与外部载荷由 `HydroFittedForce::combine` 合成，完成后私有持有原生面场及生成来源。App 负责电学模型、共同几何/时层和单位尺度绑定，Hydro 负责逐项转换、合成及原生行检查。诊断通过标量视图读取，数值步同步借用 `HydroForceDensity`；借用期间须保持结果存活且不得修改场。复制结果句柄共享已完成的数据，不复制整场；替换工作空间不改变保留的旧结果。冻结来源的生成时刻与当前求值时刻仍按既有 Hydro 规则区分，AMR 需在目标域重新装配。

完整域合成载荷由 `HydroInterfaceForce::Assembly` 私有持有；App 不直接累加毛细力和耦合力。该对象同步借用当前网格，仅活于一次载荷/求解阶段：先准备原界面载荷，再按原顺序接收一次耦合单元结果和一次面结果，最后发布借用的只读 `HydroForceDensity`。返回 false 的面提供者只能交付零场；合成失败即失效，不回滚重试。发布后的源身份仍由既有 Hydro 来源捕获绑定，几何检查、唯一电场求解和物理步接受仍由 App 编排；数值模块不持有提供者或 IO。旧 `build` 仅为文件内核，禁止恢复 App 可写合成场入口。范围见[载荷合成记录](validation/organization/QUANTUM_Modular_M2_Load_Assembly_20260919.md)。

完整域 zero/smooth 初始速度由 Hydro 的 `initializeCompleteVelocity` 生成；App 只把已校验的设置绑定到具名 profile、幅值和 XY/RZ 域，不保存第二份初始场或历史。该入口复用 `initializeState` 的一次发布，保持零物理压力、原面采样和非零散度，不求解或推进时间；边界相容投影仍为显式后续阶段。外部给定原生场及拟合域继续使用既有 `initializeState`。范围和对照见[初始速度记录](validation/organization/QUANTUM_Modular_M2_Velocity_Initialization_20260919.md)。

- `initialize` 只做初始化；每步用 `advance/assemble/evaluate` 等实际语义。`Config` 是输入，`State` 是状态，`Report` 是诊断，`Writer` 或 IO 函数负责输出。
- 数值计算不打开输出文件，不由 `writeFields` 决定计算路径。计算结果先返回，再由驱动明确调用 writer；只允许必要的显式只读诊断回调，并说明失败时刻。
- 每个头文件自行提供所需声明，不借用无关头文件的传递声明；源文件首先包含自身接口。减少不必要的完整类型依赖，不靠包含顺序让编译成功。
- 一个运行过程只有一个主时钟。回调借用的界面快照、面元和电流必须标明时间层及有效期；不能在对象移动或提交后留下悬空引用。
- 对称截面的轴/坐标、原生根身份和精确零法向导数证明归 Core 几何，随 Interface 拓扑保留；不得把它转换成固定接触点或 RZ 轴端。截面仍参与域内体积闭合，不生成自由表面端盖或镜像体积。Interface 以实际域内列体积拟合偶函数，共享标量节点强制零法向/混合导数；物理边界必须显式绑定为不可穿透对称面，RZ 物理截面仅限轴向。材料与电荷重映射核对新旧截面身份，不能据此允许拓扑事件。Hydro 的 `AxisymmetricColumn` 仅负责初值与反射声明，后续复用同一保守 Interface/载荷/NS 服务，不保留解析几何所有者。
- 保守 Hydro 的可选 `StepRetry` 时钟由 `HydroCase` 持有接受时刻，Runtime 只读取参数；其他 App 必须拒绝该扩展。`HydroStepTransaction` 拥有同网格的字段/历史/材料恢复，App 保留旧 Interface，所有物理准备成功后才共同接受。失败恢复后整步重新装配，不能复用被拒绝区间的相通量、载荷或线性分解。保守一阶动量显式选择变步长策略，普通 CN/AB2 保留固定步长；输出失败不回滚物理接受。闭域保守 AMR 由 App 暂存目标网格、Hydro 和双时层 Interface，在目标下一物理步成功后一次交换；拒绝时旧拥有者保留。Interface 统一解析 Adaptation 和规划两层支持并集，Hydro 提供双时层守恒迁移，App 不复制其数值算法。拓扑事件和电荷恢复另须实际闭合，不能沿用本能力声明。
- 普通电学通过 `DielectricSolveResult/DielectricEvaluation` 发布只读且持有自身数据的结果，不公开可写 Field。旧场明确绑定实际源电荷；电荷完成后的新快照不能被当作旧场的同迭代电荷。IO 使用结果保留的标量，历史成功结果不随下一步更新或服务销毁而改变；网格生命周期仍由调用方负责。拟合带电 P3 使用模块构造的 `FittedChargeState`，绑定电荷、电势、电流、两层历史和源时层；App 只替换完整句柄，源查询显式持有时层。cellHarmonic 使用 `FittedCellHarmonicState` 持有固定图/边界、实际材料、只读电势及来源时刻；模块区分未解种子与已解结果，承担求解、界限及面电通量，App 绑定编号和单位。联合接受仍待处理，见[普通域记录](validation/organization/QUANTUM_Modular_M4_Electric_Publication_20260916.md)、[拟合域记录](validation/organization/QUANTUM_Modular_M4_Fitted_State_20260916.md)和[cellHarmonic 记录](validation/organization/QUANTUM_Modular_M4_Cell_Harmonic_20260916.md)。
- 只有所有物理检查通过的候选才能标作接受步。场文件写入失败也必须传播；单个 DAT/CSV 存在不能证明一个时间步完成。保留完整帧清单与进程退出码。
- 喷嘴显式电荷用旧电流推进；显式选择的新隐式模式用同一步内收敛电流和两层材料历史。两者的新载荷都供下一次 NS；纯目录整理不得改变时间式，也不能把不同时层场量当作静力平衡。欠松弛时间内核、电场拟合缓存、电荷/材料历史和电极收支归 Electricity；命名 App 绑定单位、边界、材料轨迹与联合候选接受，IO 不计算时间推进。
- 保留 SuperLU 实际调用；缺失时明确失败，不回退 Eigen SparseLU。保留几何严格浮点、相关精确浮点和禁用 LTO 的现有编译属性。
- 并行计算须明确只读共享数据和独立写入对象；不得让异常直接穿越 OpenMP 并行区，汇合后明确传播失败，不发布不完整候选。确定性收支保留约定的累加顺序；复用缓存必须核验所属几何、拓扑和实际系数，不能只看稀疏结构就冻结变化的物理矩阵。性能计时不参与数值接受条件。
- 显式固定算子的右端项入口必须拥有已经完整验证的数值矩阵、原生自由度映射及装配数据，不借用可能失效的网格数组；它求解的是该固定问题，不自动代表当前可变网格。几何、系数、边界或身份改变须从完整入口重新准备。变化的右端项仍做有限性、边界通量、完整残差及失败失效检查；不能仅凭场标签相同就认定网格没有变化。
- 线性系统残差与动量、连续性、电通量、电荷等物理收支检查必须使用独立的具名参数。旧统一容限只允许在 Compatibility 读入时展开；新增调用不再据 solver 精度隐式改变物理验收门限。仅某条路径消费的显式选项，在其他路径必须拒绝，不能静默忽略。
- 普通 Hydro/EHD 的 `HydroCompleteConservation` 分别持有完整动量、初始连续性、出口可接受性及边界流量容限；正式 reader 与直接模块调用采用独立默认值。旧扁平 reader 才将原统一容限展开，App 只传递已解析值。拟合域拒绝 `completeConservation`；完整域仍拒绝 `fittedNumerics`。原判据中的范数、16 倍流量系数及舍入项保持，不能借参数拆分放宽守恒或改变原生非均匀投影策略。具体语法与迁移见[字典参考第 15 节](QUANTUM_Dictionary_Reference.md)。
- 普通电学的边界位移通量、锐界面重构/Gauss 收支及电荷闭合由 `DielectricConfig` 的独立参数持有；正式 Electricity 和分文件 EHD 共用物理 reader，旧 EHD reader 才展开原线性容限及其 16 倍电荷门限。三种电荷推进路径与电极收支使用同一个显式 chargeClosure，生效报告写实际值；拟合电学拒绝普通域专有选项，继续使用自身已发布的检查参数。范数、舍入项、实际求解和失败停止语义保持，语法与默认值迁移见[字典参考第 7.1 节](QUANTUM_Dictionary_Reference.md#71-普通电学的独立守恒参数2026-09-19)。

原生非均匀完整域的完整动量残差必须对最终速度与本步实际装配评价，不能用预测或 Poisson 的线性残差代替。当前原投影策略以 `coupledCorrectionUsed=false` 标明没有联合校正，完整残差是诊断，不新增接受门，也不把成功步解释为完整 CN 闭合；有黏性分裂缺陷单独公开。XY/RZ 压力加载沿用虚拟子单元原算子和对偶体积，均匀/拟合策略保持各自门限。见[残差记录](validation/organization/QUANTUM_Modular_M2_Adaptive_Momentum_Residual_20260919.md)。

2026-09-25：`HydroConservativeTransport` 持有实际半对偶质量、已接受相通量的相容映射及未投影的零时间动量转移，复用 `MomentumDuals`。动量迁移必须消费真实目标 C 谱系，不拥有网格修改、压力/黏性算子、边界投影或 App 时钟；返回速度不是 Hydro 接受态。每个单元内的小图系统只是几何通量提升，非负质量转移也不增加速度自由度或另一套 NS。闭域保守混合叶已使用专属物理装配策略接入原 Hydro 联立服务；开放边界仍限原根网格，不能因账本可构建就放开旧非均匀算子。

同日接续：联合投影归 `Hydrodynamics::projectConservativeMomentumRemap`，消费实际接受源的边界/材料/时钟和真实目标 C 谱系，借用原 `HydroBoundaryProjection`、受检 SuperLU 及临时可复用分解；两列约束 SVD 只做秩与相容响应，不是另一个 NS 求解器。成功/失败均恢复源报告和原压力缓存，乘子不得写成物理压力。`prepareConservativeMeshTransfer` 进一步要求两个实际 Interface 时间层及其目标迁移谱系，以各层 C 分别迁移动量，创建独立 Hydro 拥有者并重建目标物性；当前低阶压力转移保留父单元平均和积分。App 仍须把目标网格、Interface 双历史及 Hydro 统一发布，模块不代行该交换或启动时间循环。

## 4. 输入、输出与验证

2026-09-25：`HydroConservativePressure` 归 Hydro，持有不可变原生域身份、实际 `MomentumDuals` 及共享压力迹装配；RZ 几何压力源也由该服务负责。`HydroConservativeStress` 复用同一对偶分区及现有原生速度迹拟合，提供共享黏性和 RZ 环向项；接受与受检联立求解仍归 Hydrodynamics。`HydroUniformLoad` 的 native 入口复用同一压力矩阵及原 Interface 曲率查询，不能独立重建另一条施力界面。实际混合 NS 报告压力功、黏性功和动能变化率；非零压力功必须保留为缺陷，不能以谱检查或短程通过声称精确能量守恒。既有均匀入口保持。

正式 App 在自己的输入边界声明并检查文件集合，主程序路由不能成为唯一保护层。Interface、Electricity、MeshVerification 和 schema 2 Electrospray 直接调用也要拒绝无消费者的文件及子目录；旧入口仍使用原兼容白名单，普通 Hydro/EHD 的分文件检查仍归 HydroCaseInput。新增文件必须先有所属 reader 和消费者，不能只向总入口名单追加名称。

schema 2 喷雾的 `interface.dict` 由 Interface reader 独立解析，Hydro reader 只接受流体块并返回 `FittedFlowCaseConfig`。原 `hydro.dict/Interface` 仅在 App 输入边界兼容拆分，直接交给同一模块 reader，不恢复组合物理配置或新旧来源覆盖。独立文件与内联块互斥且至少存在一份；来源文件、字段路径和未消费标记随实际读取记录。标准输入和迁移方法见[字典参考第 5.3 节](QUANTUM_Dictionary_Reference.md#53-材料和流体)。

普通 Hydro/EHD 的分文件输入由 App IO 组合：`run.dict` 的 Time/Output/Log 直接生成 `RunControls`，`interface.dict` 由 Interface 的节点 reader 解析，`hydro.dict` 的物性/边界由 Hydro reader 解析。新入口不经旧 `DictionaryHydrodynamics` 对象中转；旧扁平输入只在兼容分支适配，内联界面语法也复用 Interface reader。两种来源不可混填，未消费的文件和选项明确拒绝。两个 App 共享只读 `HydroCaseSchedule` 的输出选择；EHD 保留首步求初始电场及原场时间层。能力和示例见[字典参考第 15 节](QUANTUM_Dictionary_Reference.md#15-普通-hydroehd-分文件输入2026-09-19)。

固定步长运行的 `endTime/dt` 校验归共享 `RunControls::fixedStepCount`；App 为自己的循环确定总步数，Interface/Electricity 输出器接收只读终止步号，不反算另一份步数，也不以浮点 `time >= endTime` 决定末帧。保持原整数倍容差及各 App 能力上限，普通 Hydro/EHD 的旧整数步数不经时间相除恢复。循环按已完成数 `0..N-1` 选取下一步 `1..N`，结束时不递增已达 `INT_MAX` 的步号。物理时间仍由原乘法生成，模块接受与输出完成的区别保持，见[固定步数记录](validation/organization/QUANTUM_Modular_M7_Fixed_Step_Schedule_20260919.md)。

2026-09-25：液柱 Fourier 模态由 Interface 的 `InterfaceColumnMode` 对同一发布几何求积，结果绑定来源；输出器仅消费它，不再拥有数值公式。Runtime 解析 `Time.Stop`，Hydro App 在模块联合接受后消费模态并决定停止/强制末帧；物理接受和输出完成仍是两个时钟。停止是接受步上的观测，不是插值得到的精确事件时间，也不开放拓扑变化。具体范围见字典参考第 15.1 节。

上述四类运行的场输出选择统一调用 `RunControls::writesFieldsAtStep(step, finalStep)`；IO 不另存下一次输出时刻，强制首末帧不移动定期间隔。时间间隔从零计，跨过间隔边界时选择对应接受步；小于或等于 dt 的间隔每步最多一帧，不能因区间编号溢出漏帧。物理接受仍由原模块/耦合流程决定，选择写场不代表接受计算步；普通 EHD 的首步初始电场例外留在其 App 适配器，场时间层不变。公共运行控制的实际语义见[字典参考第 5.2 节](QUANTUM_Dictionary_Reference.md#52-一个主时钟输出和日志)。

拟合液面选择和接触设置由 Interface 的 `FittedInterfaceCaseConfig` 拥有，既有 `DictionaryInterface` 解析；不得重新放回 `FittedFlowCaseConfig`。`hydro.dict` 的组合 reader 只确定块位置与读取顺序，返回分属 Hydro/Interface 的值。App 绑定几何编号、相身份和跨模块条件；Interface 数值服务仍消费解析后的坐标/域，不持有字典编号。力带宽度和流动守恒门仍归 Hydro，文件语法与配置归属分别描述。正式输入直接组合模块记录，旧输入只单向适配；见[拟合界面输入记录](validation/organization/QUANTUM_Modular_M7_Fitted_Interface_Input_20260919.md)。

正式独立 Electricity 的物性、固定域、标量名称、材料宽度和网格能力由 Electricity/IO 的 `DielectricDomainInput` / `readDielectricDomainInput` 负责。App 组合该记录与 RunControls，负责文件集合和运行控制，不重复展开物理设置。Core 几何 reader 的显式路径入口、网格 reader 的双路径入口不得修改工作目录；网格与几何必须消费同一来源，旧默认入口只作转发；文件或解析失败用异常交给调用方报告，不终止宿主进程。解析仍直接写入记录，失败后的记录应丢弃，不承诺回滚。普通 EHD 在预检绑定所需界面厚度，电学输出只接收已绑定坐标信息，不保存或传递整份 Hydro 配置。范围见[电学域输入记录](validation/organization/QUANTUM_Modular_M7_Electric_Domain_Input_20260919.md)。

CLI 仍为 `QUANTUM.exe --case <case>`。当前字典按[字典参考](QUANTUM_Dictionary_Reference.md)及[算例指南](QUANTUM_Hydro_EHD_Case_Guide.md)使用。schema 2 支持两条运动喷嘴路径，以及 perfectDielectric/cellHarmonic 的三种固定状态 Verification 组合：独立 materials/hydro/electricity/solvers 与通用 Time/Output/Log 已接通；其他 App 按各自正式参考选择输入，历史格式保留兼容入口。新字典必须同时落地解析、校验、实际消费者及参考文档，不能只增加示例或开关。

新增喷嘴算例优先使用 schema 2；新通用能力按其拥有者落在 Core/IO、Modules/*/IO、Solvers/IO、Runtime/IO 或命名 App 接线层。不得把新物性和所有边界继续塞进 Compatibility；兼容层只用于旧输入或明确标注的迁移适配。新物理字段必须有量纲检查、具名消费者、来源记录及直接检查。应用支持范围不因更换版本号自动扩大；未支持的模型、求解器和混合新旧来源必须失败。

输出固定按用途组织：`output/fields`、`output/series`、`output/meta`，不建立 `tecplot/`。喷嘴 Tecplot 场采用具名 SI 单位；保留原生面/单元诊断、物理时间及来源。控制台用中文 CFD 表述，报告当前阶段、步数、物理时间、流速、残差及体积/电荷收支；详细信息与摘要频率分开。电喷雾控制台不输出计算耗时统计，详细日志级别也不恢复这些统计；保留模块内部计时及既有文件计时列，不为此增加重复文件或字典分支。

2026-09-19 默认路径统一：PIC 验证的空间采样、切割单元/面及网格文件归 `fields`，回归报告和自力汇总归 `meta`；MovingEB、MotionSchedule 和显式 PIC AMR 报告归 `series`。旧 Electricity 初始化的默认场为 `fields/electricity_init.dat`。显式指定的旧输出路径仍按请求使用，不复制一份旧目录结果；恢复相同运行目录时若需追加迁移前的几何诊断，应显式保留旧路径。低层旧 writer 的相对前缀接口仍由调用者传入；新 App 必须按用途绑定。Windows 未处理异常的应急日志保留 `output/quantum_crash_unhandled.log`，不为崩溃处理增加目录创建或正常物理输出副本。

默认文件输出应服务物理分析，不同时重复保存相同电学量的多套网格、所有中间初始化状态或全部原生数组。schema 2 使用独立的 `Output.Fields.profile physical/diagnostic`：前者默认精简，后者按需启用原生面/单元明细；小体积守恒汇总与失败信息保留，检查本身不由写盘开关决定。不要用 `Log.level` 控制文件冗余，不默认同时生成 DAT/PLT 或 SI/代码量两份副本。历史基线不删除；有意改变的输出清单单列，保留物理数值按原输入核对。

电喷雾 `series` 的数值表统一采用 Tecplot XY Line 可读的 ordered POINT DAT，不默认重复写 CSV；`meta` 帧清单中的字符串路径仍用 CSV，`fields` 的 diagnostic 原生明细不随之改格式。复用 Core/IO/TecplotWriter 的 `History` 序列化，App IO 决定列、单位及何时记录，不让格式 writer 承担时间推进。每行写完关闭，完整行可在运行中读取；没有事件的表头不代表有一个零值事件。时间历程以 `time_s` 为横轴；单帧边界分组表以 `boundary_id` 为横轴，不伪装成逐步电流历程。字符串常量存 AUXDATA；原未定义比值须显式有效性标志，不静默替换成零误差。旧 CSV 基线保留，后处理明确支持两种输入并拒绝同名双份来源。

验证节奏（2026-09-19，按用户要求调整）：以完整职责工作包为单位推进，把该职责的实现、调用者迁移和旧接口清理一起完成后统一验收。工作包开始时明确结束条件；不能把每个函数抽取或字段迁移都作为独立验收批次，也不无限推迟集成验证。

工作包内按需要增量编译受影响目标；只有行为变化、具体风险或失败线索才运行少量直接相关检查。纯文档修改只检查修改内容和链接；纯移动、改名及简单委托不新增镜像实现的测试。数值公式/装配顺序、守恒与时间层、AMR 接受、检查点恢复或数据格式发生变化时，及时运行能定位该风险的最小检查，其余对照仍合并到工作包结束。

工作包结束做一次全目标增量构建、受影响模块的既有 CTest 和少量代表性短程对照；优先选一至两个能覆盖变更的输入，有不同实际执行路径时再补充。跨模块公共契约/构建依赖变化或阶段里程碑运行完整 CTest，不在每个中间修改重复全套。通过后仅因新的相关修改、失败或未解决疑点扩大或重跑。完整物理参数矩阵、长时间运行和收敛研究只在触及相应方法/精度或专门物理验收时安排。

纯结构对照沿用原输入和原数值判据；有意变化的路径/日志单列，不放宽门限或替换失败参考。每个工作包只保留一套必要的前后基线与结果，已有匹配基线直接复用，不逐小步复制整套源码、程序、输出或另写验证脚本。复用既有测试；只为新行为、已定位缺陷或确有价值的覆盖缺口补测试。结项在主计划简记范围、已跑/未跑检查和剩余问题，需要详细证据时整包共用一份报告。历史报告中的高频验证次数不是后续执行要求。本节调整开发验证频率，不取消程序运行时的守恒、状态有效性和失败停止检查。

迁移后现有检查见 [tests/README](../tests/README.md)；需要完整回归时执行（工作包内用 `-R` 选择相关检查）：

```powershell
ctest --preset local
```

原机的架构扫描包含 119 个拒绝/允许检查，脚本未迁入，不能宣称已恢复这份覆盖率。当前 [architecture 检查器](../tests/architecture/check.cpp) 与 [CMake 构建图导出](../tests/architecture/check.cmake) 已接入 CTest，不新增解释器依赖：扫描 `src/include` 的源码布局与所有字面量 include，检查实际目标链接、包含目录和每个生产编译单元的唯一归属。注释及原始字符串不算依赖；非活动条件分支仍检查；相对路径和 `Runtime/` 别名按真实文件解析，宏形式 include 和无法解析的项目头明确拒绝。外部系统头由编译器验证。本门不推断函数体的物理职责、字典消费者、时间层或 IO 副作用。

方向约束为 Core→Core，Solvers→Core/Solvers，Interface→Core/Interface，Electricity 与 Plasma→自身/Core/Solvers，Hydrodynamics→自身/Core/Solvers/Interface。Interface 构建目标另允许既有 `quantum_linalg`。App 负责组合；具体 App 注册与输入绑定现归 `runCase`，通用 `AppRegistry::run` 不认识具体 App，六项旧注册例外已删除并改为拒绝检查。Module→Runtime/App 和 Runtime→具体 App 的例外均为零。旧 Electricity 仍保留三个现存 Legacy 文件、五条精确跨域/入口 include，以及显式兼容库；禁止新消费者、整目录豁免和已无消费者的过期例外。见[分派记录](validation/organization/QUANTUM_Modular_M6_App_Dispatch_20260918.md)及[原架构门记录](validation/organization/QUANTUM_Modular_M7_Architecture_20260918.md)。

两个喷嘴 Setup 文件可包含 `DictionaryGeometry.h`，仅用来绑定已加载的几何定义，不允许调用 `readDict()`。模块 IO 可以读取公共字典，数值文件不可以；这些函数级规则尚未由本次精简门覆盖。现存 PIC 大驱动内的装载、输出与数值职责仍待后续拆分。

物种质量/电荷、等离子频率与热速度参考的输入解释使用 Plasma 的 `resolvePlasmaReferenceScales`；App 不再复制这些计算。保留旧 manual/物理策略回退、派生尺度重置及 Core 最终构造规则，不因归属变化调整有效单位。Core 通用单位 reader 不依赖 Plasma；schema 2 的严格 manual reader 保持独立输入契约。实际 SI/代码单位、独立运行和装载后脱离字典目录的对照见[归一化输入记录](validation/organization/QUANTUM_Modular_M6_Normalization_Inputs_20260918.md)。

Core 独立网格验证由 `Apps/Verification/Mesh/MeshVerificationApp` 持有输入能力过滤、解析探针及运行/输出顺序，数值适配、场转移、面图和 QMC2 读写继续消费 Core 原实现。主程序只注册并把目录准备交给该 App；独立验证程序仅链接 Core，不借物理模块或 Solvers 建网格。三个已有字典的实际 reader 与支持表必须同步，未知/未消费设置拒绝。四阶段只是验证序号，没有物理时钟；网格缓存读回不能表述为完整场/历史续算。见[网格算例](../cases/mesh_refinement/README.md)。

验证 App 的 `Regression*` / `SelfForceTest3` 设置与基线文件处理已退出 Plasma/Runtime。`DictionaryPIC` 保留物理 reader 及原词法规则，提供同步、不保留的附加记录消费入口；验证 App 只解释自己的记录，不复制物理配置或括号解析器。Runtime 字面量 include 限于 Runtime/Core/Solvers，不允许引用任何物理模块。`RunConfig` 仍前置声明并传递既有服务拥有者，此门不等于完整类型解耦。范围及验证见[验证输入记录](validation/organization/QUANTUM_Modular_M6_Verification_Inputs_20260918.md)。

原 `architecture_layout` 及 `QUANTUM_TEST_OUTPUT_ROOT` 属于未迁入的测试布局。现有 `architecture` / `architecture_self_test` 检查源码层级、含混工具目录、字面量 include、实际 CMake 依赖和每个生产编译单元的唯一归属；失效兼容例外也会失败。它们不替代数学、默认输出清单或用户自定义路径的检查；这些由相关运行契约及整包代表算例核对，不声称恢复了原历史覆盖率。新的 App 算例通过 `tools/run_case.ps1` 创建 `cases/<case>/output/<时间戳_标识>/`，复制实际需要的输入后运行；物理探针、图和后处理结果随该次运行保存。模块级集中检查的研究数据归入 `cases/module_history/output/<时间戳_标识>/`。build 仅保留编译配置、产物和日志；仍采用旧结果路径的驱动在所属工作包内迁移。旧兼容输入和独立参考若仍被回归消费，不因名称相似便删除。

本机新证据随对应算例或上述模块检查目录集中保存，已有 `runs/checks/` 及工作包链接保留；不为目录收口批量复制或搬动已链接的结果。旧机的 `runs/history/<主题>/<原批次名>/` 是历史归档约定，未迁入的目录不能当作本机现有数据。归档中的源码/程序仅为来源快照，不作为第二套开发入口；复用每次运行已保存的输入和程序身份，不为每个小改动重复归档完整源码/程序。清理须记录精确目标并同步活跃引用，来源不明的结果先保留，不能伪造旧日志。入口见[文档索引](README.md)、[算例索引](../cases/README.md)和[测试说明](../tests/README.md)。

启动与覆盖（2026-09-14）：新版 Electrospray 的启动选择归 `Apps/Electrospray/IO/ElectrosprayRunSession`，运行层只传递 CLI 意图。先完成配置/几何绑定校验，再询问或执行显式启动选择；在任何结果 writer 和网格初始化前处理旧输出。同目录重算必须由交互选择或 `--start new` 确认，清理范围仅为该算例普通 `output/` 目录，不动 system/cache。整个清理树先核验，拒绝符号链接和 Windows 重解析点；运行锁覆盖整个 App 生命周期。取消、未支持的续算、非法输入和只读预检查不修改旧输出。不得取消 writer 的同次运行重复写入保护来替代启动管理。

真正断点续算单独接入：当前菜单和 CLI 的 resume 只明确报“不支持”，不追加旧 CSV、不把 Tecplot 显示场当作完整状态。后续须保存并恢复已接受步的时间、网格/液面自由度、NS 状态、电荷及两层历史、守恒累计量和输出进度，并做连续运行与中断续算的同输入对照。

## 5. 实施记录（各批当时范围，不是当前待办）

2026-09-27 整体重构：InterfaceParametricReconstruction 在 Interface 中拥有一个全局连续 Hermite 修正空间、H2 目标和非线性逐格液量约束；旧固定端点体积投影已退出正式 prepare，单弧预测仅作构件。预测/修正合成后只发布一条参数曲线。时空定位按局部材料—时间片细化，轨迹通过 coordinateResolution 传递已有原生坐标关联尺度；它不是 ODE 误差认证。拒绝输入档案归 App IO，不在数值模块中写文件。完整结果与限制见[待办 §2.24](QUANTUM_Conservative_AMR_Hermite_GALS_TODO.md)，下文保留各包历史范围。

2026-09-27 规则返穿通量：`InterfaceParametricFaceTimeFlux` 归 Interface，拥有临时时空定位图、原映射交点和一份有向原生面通量；接受态仍只有唯一参数几何及独立库存。Core `LevelSetGaussTrace` 保留积分位移，供相对格面残差使用，不在几何投影后修液量。两种局部交线坐标是同一积分的图，不是失败后切换的求解链。N16 原标准涡旋全回转已完成；加密网格、当前通量的实际 Hydro 整段更新及 G1 验收仍开放。见[待办 §2.23](QUANTUM_Conservative_AMR_Hermite_GALS_TODO.md)。以下为各包当时的结果。

2026-09-27 参数界面实际入口：`InterfaceState` 增加明确的参数表示，私有拥有几何—库存候选、原生网格/域身份和完成时钟；该表示不分配旧节点场或标量 Snapshot。`InterfaceParametricMotion` 负责网格/面索引绑定、时层连续性及候选准备，旧标量消费者明确拒绝这种状态。Core `LevelSetGaussTrace` 拥有从既有验证中提取的同源 Gauss 轨迹/导数；`InterfacePrescribedData` 提供同源解析速度与精确可分离时间换算，App 只编排时钟、调用和 IO。`Apps/Interface` 显式选择 `parametricGals`，场输出直接消费唯一几何。主/独立入口在相同第 6 步失败且输出相同：结构接线完成，格面返穿的数值范围尚未补齐，不作为完整涡旋、NS 或 G1 验收。见[待办 §2.22](QUANTUM_Conservative_AMR_Hermite_GALS_TODO.md)。以下保留前包当时的范围。

2026-09-27 参数界面组件：Core `ParametricHermiteCurve` 拥有共用多项式、实际相测度及同源距离/法向/曲率查询；Interface `InterfaceParametricProjection` 持有不可变几何—库存候选，`InterfaceParametricTransport` 在旧态只读条件下准备 first-jet 源、共享原生面库存及同一联合投影。初始化建账与后续 FV 更新分开，禁止投影覆盖库存；普通 C1 节点不发布唯一二阶框架。当前推进器限均匀平面单位域、冻结同源流映射/面通量及已解析横向穿越。正式组件已由验证驱动消费，但字典 App、网格/时层身份和旧 InterfaceState/标量 Snapshot 消费链尚待接入；不得以第二份持久节点标量场掩盖表示变化。详见 [GALS 待办 §2.21](QUANTUM_Conservative_AMR_Hermite_GALS_TODO.md)。

2026-09-26 GALS 预测器：`Core/Geometry/LevelSetBoundaryTrace` 的值与切线入口共用同一离散轨迹；新 `ActiveQuarticHermiteGals` 只拥有 first-jet 与紧凑 skeleton 闭合，输出不可变几何候选及原回溯。它不推进有限体积库存、不提交 InterfaceState、不决定拓扑事件，也不读取字典。当前仅完整均匀 XY/RZ，物理流入导数、接触和混合层级另行接入；不得把该候选直接标为保守物理步。复用既有 jet 检查目标，状态与后续接线见[GALS 待办](QUANTUM_Conservative_AMR_Hermite_GALS_TODO.md)。

同一 G1 中心闭合：Core `Snapshot::nativeQuarticCenterInterval` 从原 Q4 精确系数拥有零集邻域的偏导/实际 bubble 位移允许范围；`Modules/Interface/InterfaceQuarticBubble` 拥有给定库存的独立中心求根、残差与不可达分类。原共享 skeleton 和库存不可变，不把体积误差改写成物理通量或材料速度，也不接受 Interface 时钟。当前每轮试探批量重建全域几何是待消除的开销，不作为最终局部复用实现。

同一 G1 输运候选：`InterfaceQuarticFlux` 消费旧 Q4、实际分时刻回溯及冻结法向通量原函数；`InterfaceConservativePhase` 继续拥有独立有界库存；`InterfaceConservativeGals` 统一候选顺序、身份和失败隔离。`Apps/Interface/Verification/InterfaceGalsBenchmark` 仅供应规定流动、调用和输出，不能在 App 中补改中心或库存。此入口当前只用于完整均匀闭域；Hydro、AMR 和材料事件接线尚未完成，未改变旧默认路径。

G1 回溯精度修订：`Core/Geometry/LevelSetClosedRk4Trace` 拥有闭域 RK4 位置与同阶段切线，GALS 的预测和全部分时刻通量统一选择它。旧 RK2 只作显式对照，不在新路径失败后回退。速度仍只提供同源值/Jacobian，没有额外 HessU、单独 gradU 拟合或库存补偿器。方法选择留在 Core/Interface 配置中，App 默认明确输出 `closedRk4`。

G1 精度与时间分段：实际半面端点的旧 Q4 特征零点及相通量积分区间由 `InterfaceQuarticFlux` 负责，App 不拼补相通量。旧分段策略只允许显式诊断选择，生产默认无失败回退。`InterfaceConservativeGals::Report` 可选择持有只读预测快照和未限幅库存诊断，但独立 FV 接受库存仍是唯一守恒状态。解析参考曲线及其二阶角向变分只在 `Apps/Interface/Verification/InterfaceGalsAccuracy`，不把参考 Hessian、形状或曲率送回物理更新。法向匹配是带残差/采样检查的诊断，未解析样本必须显式记录。

G1 几何复用：不可变 Core Q4 leaf 以 `call_once` 保留精确 Bernstein 系数的向外包围；张量限制采用精确展开的分方向收缩，不增加第二份近似标量。Interface 相通量可在声明的闭域/半格 CFL 契约下排除严格纯相的整个回溯支撑，仍消费原法向通量原函数及方向账本。中心闭合产生的最终标量快照经私有候选入口复用，拓扑发现、材料身份和接受时钟仍正常执行；全域试探积分尚待局部化。分项计时归候选报告，App 只输出。相同短输入的几何与原指标一致，详情只在统一 G1 证据中维护。

纯几何拆分批：Core `InterfaceGeometry` 拥有同源标量、原 PL 比例与精确来源校验；旧 `InterfaceLoadGeometry` 已整体迁为 Hydro 的 `HydroInterfaceLoadGeometry`，不保留 Core 别名。普通电学服务只消费纯几何和可选物理牵引查询，独立 sharp 入口不创建机械加载带，旧带宽参数有明确弃用提示与未消费记录。编译归属、独立运行、来源/失败检查和 118+13 个数值文件对照见[几何职责记录](validation/organization/QUANTUM_Modular_M4_Interface_Geometry_20260916.md)。后续普通域与拟合带电域只读结果、电荷两层历史及其网格转移、cellHarmonic 求解状态已接入；Hydro/Interface 候选、联合接受及完整 AMR 仍按主计划推进。

普通完整域电力批：`DielectricElectrostatics` 仅发布实际 P3 牵引或弥散极化系数、材料与来源；`HydroUniformLoad` 拥有法向/切向加载、标量梯度、边界自由行和力积分。`SharpDielectricDiagnostics` 已从 Electricity/IO 迁到 EHD App/IO，读取两个模块结果。24 项检查及 118+13 个数值文件对照见[普通域加载记录](validation/organization/QUANTUM_Modular_M4_Uniform_Load_20260916.md)。该批遗留的 Core 综合几何依赖已由后续纯几何批消除；场结果的可变旧视图尚待收紧。

拟合域 P3 牵引批：Core `InterfaceTraceGeometry` 只保存同源局部查询，Electricity 返回实际场与法向/切向牵引，Hydro 持有完整主/横向梯度载荷集合并加载。`HydroFittedInterfaceLoadGeometry::faces()` 为只读，加载必须使用其同一采样集合及接收域网格身份；不能删面、重排或以另一份几何冒充同源输入。App 保留单位和结果组合；普通完整域旧电力已由后续批移入 Hydro，旧 Core 机械载荷字段也已由后续纯几何批迁入 Hydro。精确范围、22 项检查及输出对照见[Maxwell 牵引记录](validation/organization/QUANTUM_Modular_M4_Maxwell_Traction_20260916.md)。

cellHarmonic 极化批：Electricity 从实际 FV 场与 Core 域生成物理系数；`HydroFittedScalarGradientLoad` 持有压力相容梯度、自由行、对偶体积与加载运算。接收方校验源网格身份、选中单元和 XY/RZ 度量；App 只连接结果及单位。电学源码全部进入独立库，不再有 FittedDielectricPolarization 的主程序例外。两侧独立解析/积分检查、11 个标准输入及详细输出对照见[极化加载记录](validation/organization/QUANTUM_Modular_M4_Polarization_Load_20260916.md)。拟合域锐界面牵引/加载由后续批接续，普通完整域及完整接受来源仍待收口。

正式独立电学入口使用 `Apps/Electricity/ElectricityApp` 和生产 `DielectricEvolution`，固定界面运行不创建 Hydro，不加载 Plasma 配置；同一 App 已通过 Electricity/Solvers/Core 独立链接和实跑。旧 SingleModule/Electricity 的运行器、reader/config 归 `Apps/Electricity/Compatibility`，保持原模型语义，旧 SparseLU 不被当作正式服务后端。电学物理输入归 `Modules/Electricity/IO`，数值服务头不再声明读文件接口。`FittedDielectricPolarization` 的 Hydro 加载依赖已由后续极化批拆分；详见[独立入口记录](validation/organization/QUANTUM_Modular_M4_Electricity_App_20260916.md)。下文先前阶段描述按本条接续。

拟合域电学子步：`FittedChargeEvolution` 接收归一化材料、电图/边界、最终几何及材料图，执行原显式或 BE/二阶历史电荷—P3 场子步及同迭代电极检查。App 组合一份电学状态，只保留 SI/编号适配、耦合与共享载荷组装。迭代 observer 只观察已完成场求解，不提供数值电流；候选失败不得改变旧状态。独立链接、回滚及新喷嘴对照见[拟合电学记录](validation/organization/QUANTUM_Modular_M4_Fitted_Charge_20260916.md)。AMR 历史转移、完整状态来源与旧运行入口仍待后续收紧。

普通 EHD 电学服务：电场/电荷所有权、待消费电流与时层检查、显式充电/对流/重映射和电极收支归 `DielectricEvolution`；App 只安排原耦合顺序并输出模块结果。场量在完成电荷子步后仍属于旧受力时刻，不得冒充新时间同迭代解；数值失败后模块拒绝取场及继续推进。`quantum_electricity` 独立链接 Solvers/Core；仍依赖 Hydro 的旧 `ElectricityModule` 和 `FittedDielectricPolarization` 显式保留为待迁移消费者。实际 XY/RZ/规定边界及独立电学证据见[本批记录](validation/organization/QUANTUM_Modular_M4_Dielectric_Evolution_20260916.md)。

迁移后 M3/M6 首批：普通 Hydro/EHD 的节点液面、历史、输运及材料映射由 `InterfaceNodalState` 持有，App 直接借用同一边界报告；不得恢复可变副本。`quantum_hydrodynamics` 独立编译并通过不含 App/Electricity 的实际求解检查，已有浮点/SuperLU 契约随源码保持。原测试源码未迁移，目前检查范围以 [tests/README](../tests/README.md)和[本批证据](validation/organization/QUANTUM_Migration_M3_Hydro_20260916.md)为准；历史条目不代表本机复现。

M6 PIC 进程输入：65 项旧隐式环境覆盖统一由 Plasma IO 捕获为 `ImplicitPicOverrides`，随 `PicCaseInput` 持有；数值驱动、极化求解和轨道分段不得重新读取这些环境变量。程序化默认输入不隐式套用进程覆盖，需要时由调用者显式装载；缺失值仍遵循原字典/默认策略，路径仍相对于执行目录。保留旧布尔/数值解析差异，禁止用模块内逐步查询恢复第二套输入来源。双入口与输入快照见[环境输入记录](validation/organization/QUANTUM_Modular_M6_Pic_Environment_Inputs_20260918.md)。`ignoreHostHints` 必须对初始和每次跨单元定位都忽略提示编号，交给现有 Core 位置查询，再按原策略选祖先单元；不得沿用仅因编号有效就未经位置核对的宿主。真实粗细网格、嵌入边界和续算见[轨道定位记录](validation/organization/QUANTUM_Modular_M6_Pic_Segmenter_Lookup_20260918.md)。

M6 Core 网格系数缓存：`DiscCoeff` 默认初始化全部有效成员；QMC2/BMC2 共用 `MeshCacheIO` 的成员序列化，保留原生记录长度及偏移，将填充字节写为零。不得在 Plasma 或按零计数模板额外清除 RHS，旧成员位模式原样恢复，不能用新默认值覆盖历史有效常量。旧原生 ABI 兼容不等于跨平台格式；冻结小样本、真实新旧恢复及确定性范围见[系数缓存记录](validation/organization/QUANTUM_Modular_M6_Mesh_Cache_Coefficients_20260918.md)。

M6 Core 空闲标签缓存：QMC2 继续保存四类全容量数组，但只读取其有效空闲前缀，未使用尾部写零。前缀长度由实体占用确定，包含高水位以下的孔洞，并保留非叶父记录；不得以活动叶数或高水位代替占用统计，不重排或修复有效前缀，不修改实时分配表。旧 reader 保留原校验和提交流程，输入文件不被改写。真实细化/粗化、双向旧缓存及重复快照证据见[标签缓存记录](validation/organization/QUANTUM_Modular_M6_Mesh_Cache_Label_Tails_20260918.md)。

M6 Plasma 边界：`PlasmaRunInputs` 只包含两个 PIC 驱动实际消费的八项物理服务输入，App 适配必须保持原只读拥有者及输出开关，不重读或改写物理设置。模块不得继承 Runtime 接口、保存 `RunConfig` 或依赖 App 包含路径；`quantum_plasma` 只通过明确库依赖连接，禁止恢复主程序对 Module 源文件的通配编译兜底。宏步/子循环、随机数、碰撞及恢复语义属于 Plasma，拆分时按实际输入验收；本批明确范围见[Plasma 边界记录](validation/organization/QUANTUM_Modular_M6_Plasma_Boundary_20260916.md)。

M6 碰撞输入后续：四类物理 reader 与优先级选择、依赖校验、报告写出归 `Modules/Plasma/IO`，Runtime 只保存完成装载的 `PlasmaCollisionInputs`，不得恢复四类对象的独立解释或部分发布。通用单位/归一化仍复用原服务；保持既有 mcc/pic/独立字典逐块优先级。实际 MCC、AMR、续算及错误码修复见[输入记录](validation/organization/QUANTUM_Modular_M6_Plasma_Inputs_20260916.md)，本条接续下文旧阶段的“Runtime Plasma 装载待迁移”说明。

M6 碰撞运行服务：`PicCollisions` 持有物种数据库、MCC 算子及其输入拥有者；禁止每步重新初始化或无拥有者地转交其借用引用。驱动保留显式整步/隐式半步调用顺序和已有粒子生成合并规则，服务统一单位换算，显式初始化诊断归 IO。实际生命周期、三种 RNG 模式、电离与 perParticle 弹性续算范围见[服务记录](validation/organization/QUANTUM_Modular_M6_Pic_Collisions_20260916.md)。

M6 碰撞阶段：`PicCollisions::Stage` 按运行持有原始/待提交生成缓冲、批次及累计统计，只发布 const 视图；初始化算子、物种索引和默认权重须覆盖阶段寿命。显式在每物种推送后碰撞/过滤，全部推送后才提交；隐式每半步遍历全部投射物，之后稀疏化、核算和提交。未知物种警告、宿主修复和无效半步的旧差异保留。累计清空与批次发布仍由调用方在原时机执行；原 apply 接口不自行提交，`advanceImplicitHalf` 明确执行隐式半步后的稀疏化/记账/提交，不重置 RNG 或提前累计批次。借入实时阶段视图，不保存宏观步内引用；禁用或无效 dt 保留原批次和调用者统计。见[阶段记录](validation/organization/QUANTUM_Modular_M6_Pic_Collision_Stage_20260918.md)。 完整半步接线见[隐式新增粒子阶段记录](validation/organization/QUANTUM_Modular_M6_Pic_Implicit_Spawn_Stage_20260919.md)。

M6 碰撞报告：显式 MCC/离子/中性粒子文件由 `IO/PicCollisionDiagnostics` 持有，能量汇总归 `Diagnostics/PicStepDiagnostics`；驱动保留采样阶段和时间换算。通道种类必须查询初始化后的启用通道并复用 `PostCollisionModel` 的分类，不得按报告中的字面名称另立判据。统计和粒子只读借用，不修改 RNG 或推进碰撞。原二维电子代理与三维中性粒子能量、恢复截断及输出阶段语义见[报告记录](validation/organization/QUANTUM_Modular_M6_Pic_Collision_Reports_20260918.md)。

M6 壁面流失统计：`PicBoundaryLosses` 按宏步拥有原事件统计，只发布 const 结果；调用方保持物种编号、默认权重、准备后的速度尺度及事件阶段。几何 ID 解析、SEE、侵蚀、电极/能谱观察者和 IO 保留原调用位置，统计服务不反向调用它们。显隐模式的非有限值/周期处理和分类门明确保留；隐式域电荷仍按原物种顺序归约，不能改用显式逐事件和。验证及尺度接线修正见[壁面流失记录](validation/organization/QUANTUM_Modular_M6_Pic_Boundary_Losses_20260918.md)。

边界通量/能量公式归现有 Diagnostics 的值类型；IO 只序列化这些值。隐式侧面积、电荷和电流由 `PicBoundaryLosses` 在冲击记录后同步归约，返回不含借用的结果，保留 SI 面积和时间回退。显式累计窗口、瞬时输出单位及隐式单步报告不得混用；整数计数、零面积能量 NaN 与空行策略保持。原事件归约及实际连续/续算验证见[边界通量记录](validation/organization/QUANTUM_Modular_M6_Pic_Boundary_Rates_20260919.md)。

PIC 活粒子统计归 Diagnostics：`PicSpeciesCount::sample` 读取真实粒子容器，正权重覆盖默认值，逐粒子累加，不改成缓存 `number` 或额外有限性过滤；总宏粒子计数保持整数。显式报告窗口只在新区间读取首快照。显式体内能谱保留按恢复后物种名称索引查默认权重，隐式汇总保留按物种位置取值，不能为合并入口改变这一区别。观察回调同步消费、失败原样传播，不保存粒子/回调借用，不改变 RNG 或物理推进。旧字段/单位/采样时机及实际对照见[粒子统计记录](validation/organization/QUANTUM_Modular_M6_Pic_Particle_Samples_20260919.md)。

M6 粒子归属：粒子状态与缓存、合并/分裂归 `Modules/Plasma/Particles`，粒子壁面冲击/侵蚀执行归 `Modules/Plasma/Boundary`，禁止 Core 反向包含 Plasma。旧顶点电学粒子混合实现归 `Modules/Electricity/Legacy`，仅由 `quantum_legacy_electricity` 显式链接 Plasma；正式电学独立库不带该依赖。粒子缓存及验证限制见[粒子归属记录](validation/organization/QUANTUM_Modular_M6_Particle_Ownership_20260916.md)。MovingEB schema/reader 后续已归 Plasma，Core 用不含物理类型的扩展块保留输入，App 通用开关检查不链接 Plasma；控制器持有不可变配置，具体兼容与原失败见[壁面配置记录](validation/organization/QUANTUM_Modular_M6_MovingEB_Inputs_20260917.md)。MotionSchedule schema/reader 后续也已归 Plasma；`LegacyGridMotion` 只接收并拥有明确线段绑定，不能在计算服务内查找场景边界编号。驱动保留原网格重建/更新频率，原连续/续算差异不得当作正确性通过，见[栅极运动记录](validation/organization/QUANTUM_Modular_M6_Grid_Motion_20260917.md)。

M6 分阶段运动状态：`LegacyGridMotionEvolution` 拥有初始参考和最后接受的运动，驱动不得把已移动检查点轮廓重新当作初始参考。恢复按已完成步数及原频率还原生效状态，提前核验并恢复内点；不相容输入不得删除已有输出。该移动路径必须求解当前极化矩阵，不得仅凭稀疏结构冻结旧系数。实际固定/细化恢复和数值对照见[恢复记录](validation/organization/QUANTUM_Modular_M6_Grid_Restore_20260917.md)。

M6 检查点 IO：`PicCheckpointStore` 统一 PIC 目录选择、元信息、RNG 与保留规则，接收明确组件清单；隐式七组场历史由 `ImplicitPicFieldCheckpoint` 序列化。禁止两条驱动重新复制文件格式实现；显式输出预检与两种自动选择要求的旧差别按实际能力保留。组件拥有者及读写顺序不转交 Core 通用目录服务，具体兼容和未闭合的事务/完整性契约见[检查点记录](validation/organization/QUANTUM_Modular_M6_Pic_Checkpoint_20260917.md)。

M6 Core/Lagrangian 后续：固定均匀网格的闭合轮廓归 Interface 的 `UniformInterfaceContour`，几何快照不混存牵引力；节点牵引积分/分配归 Hydro 的 `HydroNodalTraction`。App 不保存旧 `SharpInterface` 来调用该路径，Core 不保留物理方法转发。旧快照保留原节点，复用前检查完整网格身份；全部载荷校验后才累加输出。活跃毛细面载荷与节点参考积分的区别保持。闭合轮廓范围见[节点牵引记录](validation/organization/QUANTUM_Modular_M6_Nodal_Traction_20260916.md)。后续旧毛细力、IBM 插值/分配、固定速度力已归 Hydro 的 `LegacyNodalForcing`，旧 Maxwell 叠加归 Electricity 的 `LegacyNodalMaxwell`；Core `initializeGeometry` / `updateGeometry` 不隐式计算物理载荷，旧 App 显式组合。旧物理核的检查及限制见[旧节点核记录](validation/organization/QUANTUM_Modular_M6_Legacy_Nodal_20260916.md)。后续几何与状态已分离：Core 的 `LegacyNodalGeometry` 不拥有速度/力数组，Hydro 的 `LegacyNodalState` 组合几何并拥有速度历史、受力和旧 SIC2/混合输出；旧电学兼容库显式链接 Hydro，正式电学不带此依赖。禁止把物理状态或方法转发放回 Core。状态拆分范围见[节点状态记录](validation/organization/QUANTUM_Modular_M6_Nodal_State_20260917.md)。后续蒸发速度及三层历史也已归 Hydro，旧电学显式返回并交付速度，Core `InterfaceNode` 不再保存该字段；通用 IIM 系数/模板仍保留。原重建清零、固定保留及 SIC2 原生格式必须按冻结旧缓存核验，不以物理状态归属完成宣称可移植序列化或完整旧瞬态链通过，具体见[蒸发状态记录](validation/organization/QUANTUM_Modular_M6_Nodal_Evaporation_20260917.md)。

M2 规则 RZ 行：完整域和拟合域的度量/物性准备共用 `Detail/HydroNativeRZRowPreparation.h`，不得恢复两份对偶体积、面密度或黏度取值实现。调用者先完成原拓扑/活动邻域/边界筛选，再采样物性；拟合域保留活动流体顶点校验，纯矩阵装配不增加对流密度读取。原生与拟合的合法模板差异不以规则行合并为由删除，实际范围和细化输入限制见[规则行记录](validation/organization/QUANTUM_Modular_M2_RZ_Row_Inputs_20260917.md)。

M2 固定步接口约束（2026-09-19）：完整 XY/RZ 与拟合 RZ 的正式求解只通过 `advance/prepareStep`，先显式初始化并通过 `updateMaterials` 提供物性来源；禁止恢复只带 dt 的旧 options 重载或无来源记录的物性写入口。拟合网格转移必须捕获已接受且非停止/非待接受的状态，完成步号/时间与接受记录相同，必需物性记录按原冻结或时变规则核对；转移只保留该历史，不增加物理步。旧节点非增量耦合与缓存仍是显式隔离的兼容能力，不据其存在放行拟合 AMR 的来源检查。原装配/投影/转移算式不因接口退出改变，范围及对照见[固定步接口记录](validation/organization/QUANTUM_Modular_M2_Managed_Entry_Retirement_20260919.md)。

M2 拟合投影生命周期：`HydroFittedProjection::System` 持有一次准备的只读压力方程及私有逆质量；时间步和 AMR 共用相同校正接口，重复历史修正只更新右端。域只在当前同步调用内借用且必须保持地址、算子和定义不变，禁止跨步或跨 AMR 保存；不得通过校正参数换入其他域。AMR 相容路径不装配，修正次数、物理压力和历史时刻保持。

M2 耦合方程生命周期：完整域与拟合域共用 `HydroCoupledBalance::Assembly/Operator`；本步 A/C/D 与压力参考策略确定后一次完成块装配，平衡修正只重新计算右端。算子仅暴露 const 矩阵，稀疏策略保留本步 C/D，禁止通过修正请求更换几何算子。A/C/D 的拥有者须覆盖整个同步求解且保持不变，算子不能跨步或跨 AMR 保存；不为此复制物理场。原生加载积累顺序与稀疏表达式的浮点边界、物理门限保持。直接方程契约及完整输出对照见[耦合算子记录](validation/organization/QUANTUM_Modular_M2_Coupled_Operator_20260916.md)，不代表其余动量行/投影已合并。

M2 拟合喷嘴接受约束：App 在界面/电学/载荷来源时层核对、既有通量检查和材料更新后，通过不抛异常的句柄/记录交换提交电学、载荷、材料及体积/通量记录，随后发布新时间诊断。禁止恢复可能复制场数组的提交赋值或在写盘后才更新接受记录；交换后须刷新后续载荷借用。错误日志区分物理接受与输出完成，失败立即停止；这不是 Hydro 回滚、重试或完整 AMR 联合保证。候选曲率/迭代输出必须保留未接受含义。实际检查见[拟合域提交记录](validation/organization/QUANTUM_Modular_M2_Fitted_Commit_20260916.md)。

M2 普通 App 接受约束：`HydroCaseCoupling::completeStepWithVelocity` 只完成物理子步，`publishCompletedStep` 只发布已完成的拥有型结果；不得恢复在电荷完成前写新流场的顺序。App 在 Hydro/Interface/材料/电荷及检查后接受，再写新时间结果。`HydroCase::progress()` 的 `acceptedStep/outputStep` 分别记录物理接受和输出完成，沿用同一时钟；任何计算或 IO 失败都停止 Case，接受步号不提供旧状态回滚。电荷发布核对原步号/受力时刻/dt 且只消费一次，不能保留速度/材料轨迹借用或重算电流。失败报告只用已采集诊断，不查询停止的模块；诊断写入失败不能掩盖原始错误。正常与故障的实际验证见[接受与发布记录](validation/organization/QUANTUM_Modular_M2_Coupled_Acceptance_20260916.md)。该约束尚未代替拟合域/AMR 的完整联合接受。

M3 普通节点状态约束：公共 `FieldValues` 在 const 场上只返回 const 元素/指针，不能通过公开存储取得可写数组；可变拥有者仍使用原分配、填值、复制和别名检查。`InterfaceNodalState` 接管原几何/历史，不为只读接口复制整场。几何/边界报告借用截至推进或销毁，数值失败后拒绝新的几何/边界查询，已有借用必须结束；材料轨迹继续逐次检查有效期。预检查失败不消耗状态，实际数值失败停止服务；这不提供全状态回滚。独立检查、全量构建和实际对照见[普通界面只读记录](validation/organization/QUANTUM_Modular_M3_Nodal_Const_20260916.md)。普通 EHD 接受及发布由上述后续批接续，单个场文件仍不能认定整组输出完成。

M3 普通节点材料运动约束（2026-09-19）：请求材料映射时，`InterfaceNodalState` 同时持有最终标量的 Hermite 查询快照；`motionFor` 在原网格/拥有者/索引/步时/拓扑核对后，将消费者保留的节点标量与既有上一层历史逐字节比较，再借出配套快照和原轨迹。App 不得另配目标快照或重建速度轨迹；EHD 使用实际电场结果保留的标量作来源，Electricity 仍只依赖通用 Core 输入。提取只扫描一次来源；逐点查询继续检查有效期。借用必须在映射、网格及状态均存活且未推进的同步范围内消费；推进/数值失败/析构使回调失效，预检查拒绝保留有效期。该快照是查询表示，不是另一份可推进标量；不提供速度绝对时层证书或回滚。验证见[节点映射记录](validation/organization/QUANTUM_Modular_M3_Nodal_Material_Binding_20260919.md)。

M3 普通节点输入约束（2026-09-19）：运动选择、重初始化周期及规定标量的接触/入口数据归 `InterfaceNodalConfig`，不得重新放进 `HydrodynamicsConfig`。值类型只依赖 Core 坐标和标准容器，不依赖 Hydro、Runtime 或 App；原仿射标量由它计算。旧综合 reader 可保留原解析与跨配置兼容检查，但分别返回各模块记录；Hydro 边界适配只绑定命名标量入口与显式流动边界，机械加载只接触点和所需边界数据，IO 只接收列/取场所需的选择。App 继续将输入绑定为现有 Interface 数值选项，不增加一份可推进状态。输入语法仍为旧 `hydro.dict`，实际限制和验证见[界面输入记录](validation/organization/QUANTUM_Modular_M3_Nodal_Input_20260919.md)。

M3 状态/初始化约束：`InterfaceState` 私有持有显示值、约束、接触初始导数、重构及工作缓存；App 只拥有可移动句柄和几何编号，不得恢复直接填场、写缓存或替换来源/面元的入口。`InterfaceInitialization` 拥有初始接触模板与几何发布，App 绑定输入身份并调用已有 Core 距离服务。`InterfaceMeshTransfer` 拥有 Q4 几何/面元转移，Electricity 只消费匹配的 Core 面元并检查电荷/历史。Prepared 不是可用几何，Ready 不是联合接受；失败候选清除几何且不得继续推进。只读数据同步借用，保存不可变几何须按值持有共享句柄。完整接受来源和联合提交仍待迁移；25 项检查、实际三步及原生拒绝范围见[状态与转移记录](validation/organization/QUANTUM_Modular_M3_Interface_State_20260916.md)。

M3 独立界面入口：`Apps/Interface` 直接调用生产 `InterfaceMotion::preparePinnedToAxis`；不得反向调用电喷雾初始化或另写输运/修正。原固定域单元判相、直线完整面覆盖与几何边界匹配归 `Core/Geometry/FittedDomainGeometry`，不搬入电学、流场边界条件。`Modules/Interface/IO` 负责 `interface.dict`，可以依赖公共字典设施；Interface 数值文件仍禁止依赖 IO、App、Hydro 或 Electricity。`InterfacePrescribedData` 只提供解析初值/速度/同源法向通量，不求解 NS。当前独立初值适配只支持代码单位、固定均匀 RZ 和单个湿壁钉扎端到轴线液面；未知/无消费者配置拒绝，不据此限制原喷嘴原生双四次及 AMR 能力。新入口只接受新输出目录，原 App 的输入和启动政策不变。接口、字典、原机独立链接与实跑见独立界面算例记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M3_Interface_Case_20260916.md`）；本机重建输入的三步证据见上述状态与转移记录。

M3 相数据归属（2026-09-19）：`InterfaceState` 私有持有 `InterfacePhaseState`，体积测度与相通量由 `phase()` 只读借用。Ready 只表示几何可用，`hasPhase()` 才表示相数据已准备；初始化失败不发布部分账本，成功后拒绝重复初始化。运动必须消费同一旧状态的相数据，最终测度/通量归新候选；几何 AMR 不复制旧索引账本，需以目标 Hydro 通量重新初始化。App 提交交换完整状态，保留旧状态到本步诊断结束，不新增整场复制。回调的实际速度时层仍须由提供者保证，不把所有权封装称为完整时间来源证书；联合接受与 SI 全局诊断仍归 App。原数学模板和数值门限保持，见[相状态记录](validation/organization/QUANTUM_Modular_M3_Phase_State_20260919.md)。

M3 材料映射消费（2026-09-19）：`InterfaceMotionResult` 私有保存原轨迹、旧/最终几何共享身份和实际步时；使用 `materialMapFor(source,target,time)` 验证后提取同一回调，禁止从可编辑的诊断报告重建来源。导电接线用旧电学电势来源、当前最终液面和两端时间核对，模块仍不反向依赖 Electricity。提取后的回调保留原速度捕获，数据所有权/借用期限继续遵守原提供者契约，不重算输运或把材料点再次投影到修正液面；任意外部速度回调的语义和完整联合接受仍须各自保证。几何身份与原轨迹验证见[映射来源记录](validation/organization/QUANTUM_Modular_M3_Material_Map_Binding_20260919.md)，完整时间参数以以下时间历史条目为准。

M3 时间历史（2026-09-19）：`InterfaceState::progress()` 只返回 Ready 几何的本地完成步号/时刻；初始为零、几何 AMR 原样传递，失败或移动后不可查询。`InterfaceStepTime` 由 App 唯一时钟提供步号、begin/end 和原 dt，禁止用 end-begin 改写求解 dt 或另累加时钟；运动入口核对连续性后才查询速度/准备候选。`materialMapFor` 的第三参现为完整 `InterfaceStepTime`，替代前批单独 dt。喷雾适配器从匹配的 pending Hydro 直接构造当前/历史速度，AMR 对 Interface/电学步时作相同来源检查，独立场输出拒绝错标时间。记录表示模块完成，不表示联合接受或可恢复历史；详见[时间历史记录](validation/organization/QUANTUM_Modular_M3_Time_History_20260919.md)。

M3b 完整运动规范：`InterfaceMotion::preparePinnedToAxis` 拥有原 RZ 钉扎端到轴线液面的支撑带、输运准备、边界约束、修正自由度、最终状态与材料轨迹；App 仅绑定原 Hydro 速度/通量及喷嘴编号。禁止在 App 恢复这些算法或进行第二次材料输运。速度回调副本由材料映射保留，须拥有其数据或保证生命周期；只读并行查询由提供者明确承诺。域/通量同步借用，结果不保留借用。新状态只在终结阶段填充，异常后丢弃，返回不代表耦合接受。等层三次与原生双四次能力保留，不扩称为任意闭合/多接触/XY 运动；规定速度独立检查、原输入对照和剩余模块入口见完整运动记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M3b_Interface_Motion_20260916.md`）。下列早期状态/通量条目的后续职责以本条和较新记录为准。

M3 运动支撑：内部 `InterfaceMotionSupport` 集中支撑/压力口单元集合、带宽选择及逐阶段全闭单元检查；几何和布局仅同步借用。输运及体积修正沿用旧标量选定的带宽，压力口固定自由度消费同一集合，不另建 App 副本。原生选择保留已有效的初始/翻倍结果，仅在翻倍冲突时以同一整单元证书二分收缩；不得试跑输运、减小原支撑/压力口余量或跌破初始带宽。后续阶段不重新选带宽以掩盖液面接近压力口。实现、真实冲突拒绝及静态细化短程范围见[带宽记录](validation/interface/QUANTUM_Interface_Band_Selection_20260917.md)，不等于动态 AMR 或空间收敛验收。

M3c 通量与修正规范：相通量及其灵敏度、线性和完整非线性体积修正归 `InterfacePhaseFlux/InterfaceVolumeCorrection`。通过 `InterfaceTransportDomain/InterfaceNormalFluxView` 借用公共域、界面约束和冻结通量，不接收 Hydro 具体类型，不复制网格或另建速度算法；Hydro 适配只调用已有重构。四个通量查询须属于同一个原生面原函数，所有者在同步调用内存活且不变，返回结果不得保留这些借用。原压力口平均相流量诊断仍在 Hydro。Interface 库现依赖 Core 与公共线性代数，保持 SuperLU 和原精确浮点属性；禁止恢复旧 Hydro 相通量/体积修正文件或转发层。完整运动、跨调用材料轨迹及旧数组时间来源尚未统一，不能以本批独立修正检查替代完整输运验收；当前证据见相通量与修正记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M3c_Interface_Transport_20260916.md`）。下列 M3a/M3b 条目是各批历史边界，库依赖以本条为准。

M3b/M3c 相体积边界：原生半面切分、根细化、相体积与灵敏度归 `InterfacePhaseGeometry/InterfaceCellMeasure`，只接收 Core 的 `NativeLeafDomain` 与同源标量/面元；不得重新引入 Hydro 域或 NS 边界依赖，也不保留旧同名转发实现。原生面视图借用现有网格，不另存面数组；原逐面算式和精确浮点属性随实现保留。实际对照及尚未迁移的速度通量/完整修正流程见界面相体积记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M3b_Interface_Measure_20260916.md`）。

M3a 状态首批：自由液面数据/掩膜/钉扎坐标/几何/面元/模板缓存归 `Modules/Interface/InterfaceState`；喷嘴组合该对象，只另存几何编号和初始接触诊断。不得恢复 App 内第二套标量类型或独立几何/面元指针；实际构造点通过 `retainGeometry` 核验来源对。保留三次/双四次表示及原构造、数值检查和联合接受时序；显示字段不是权威标量。`quantum_interface` 为依赖 Core 的真实库，独立链接不带 App/Hydro/Electricity；后续输运服务须继续满足规定输入独立运行，不以反向依赖完成接线。可写准备数据及完整运动/体积约束尚未迁完，范围和实际对照见界面状态记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M3a_Interface_State_20260916.md`）。

M5 AMR 支撑：Interface 依据已发布标量的全闭单元证书提供液面保护，Hydro 提供物理边界两侧实际单元保护，Core 执行面邻居扩展；App 合并后调用原 Q4 拓扑保护。原标量验收掩膜、无变化先判、历史转移和联合提交顺序保持；不得在 App 恢复物理保护计算或把拓扑保护替代标量验收。实际动态粗化及隐式二阶历史续算的本机范围见[支撑记录](validation/organization/QUANTUM_Modular_M5_Adaptation_Support_20260917.md)；动态细化及完整联合拒绝已由下述联合转移批接续。

M5 每轮拓扑请求由 Core 的 `QuadtreeMeshRefiner::buildMappedTargetLevels` 组装：细化继承源叶指标，粗化保留所有覆盖旧叶的最大值，调用原层级策略后锁定已做拓扑闭合的保护掩码，并沿用原无变化判断。源网格、指标和保护在同步调用期间保持不变；相区分类由调用方按目标中心提供，Core 不认识流体相或单位。App 只绑定这些输入，Hydro 仍计算速度/梯度指标。请求不克隆网格或转移物理历史，只有实际变化才准备候选；输出是工作缓冲，拒绝后不得使用部分结果。

M5 联合来源：AMR 在无变化先判之前核对已接受的 Hydro/电学步号和时刻，以及 Interface、电势、电荷面元和载荷保留的共同标量。接受观察只在物理接受和输出完成后同步借用，不提供数值输入或保存字段副本。喷嘴 App 库由主程序与检查共用；五阶段拒绝、原体积门拒绝、有效候选丢弃及后续真实提交的边界见[联合转移记录](validation/organization/QUANTUM_Modular_M5_Joint_Transfer_20260917.md)。本记录接续上一段的动态细化与联合拒绝待办；已提交后的 IO 失败仍停止，不承诺回滚。

M5 候选提交（2026-09-19）：`ElectrosprayMeshCandidate` 私有拥有目标网格、域、电学绑定、Interface、Hydro、导电/载荷状态及材料；外部只能查询报告和 const 目标状态，禁止逐项提取或替换。准备源必须存活且保持不变到提交；提交入口再次检查拥有者、原生网格身份、Hydro/Interface 历史和电学/力来源，然后按原顺序整体交换，检查完成后的交换无求解、回调或分配。候选只能提交一次，提交后目标查询拒绝，旧拥有者保留到原诊断区间结束；报告仍可读。无变化先判不创建候选记录或复制历史。该检查不复制所有源数组，不替代同步源不变约定；不承诺并发安全、任意 App 数组改写检测或提交后的 IO 回滚。验证范围见[候选提交记录](validation/organization/QUANTUM_Modular_M5_Candidate_Commit_20260919.md)。

M2b 连续性约束：完整域速度/散度诊断、拟合域 NS/AMR 的净流量和尺度、初始投影的原全域判定归 `HydroContinuityBalance`。这是体积流量平衡，不在 App 另算一套残差；不同阶段的原流量不等式、相对残差或全域范数判据不能以“统一”为由互换。速度每次修正后重新评价，已有结果只用于该速度候选，不作为跨步缓存；保留原生半面累加顺序、粗细行、XY/RZ 测度、并行归约及失败传播。实际对照和剩余装配职责见连续性记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M2b_Continuity_20260916.md`）。

M2b 载荷来源约束：原生载荷须携带装配时生成的 `HydroForceSource`，记录域/共享网格身份、源时刻、固定或时变属性以及可用的不可变液面。消费端不得通过重新标时或临时生成来源冒充已接受载荷；AMR 的新载荷绑定新域/新液面，保留原物理时刻。来源预检查不修改物理场，复用既有精确网格校验，不增加全场散列、网格复制或第二套时钟。原力密度单位、时层与受力公式保持。可变标量、公开数组及剩余装配职责的限制和实际复测状态见载荷来源记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M2b_Load_Source_20260916.md`）。

M2c 运行约束：完整均匀域 Hydro 与普通 EHD 各自的命名 App 拥有主循环，共用 App 层不含主循环的 `HydroCase`；不得恢复 `HydrodynamicsModule::runCoupled`。现有电学回调协议归 `Apps/Hydrodynamics/Coupling`，不放数值模块；生产数值入口仍是公共 `Hydrodynamics::advance`。原受力及流量/压力平衡计算归 `HydroInterfaceForce/HydroFlowDiagnostics`，只接收所需参数，不依赖 case/Runtime 或写文件。原场与诊断序列化归 Hydro/IO，App 按原时层调度；NS 成功不代表随后界面、电荷和写盘阶段已成功，失败实例不得继续推进。均匀域能力、原综合字典、界面接线与公开取场的剩余迁移见运行职责记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M2c_Runtime_20260916.md`），不得把本批当作 M2–M7 全部完成。

普通 Hydro/EHD 运行设置约束（2026-09-19）：`HydrodynamicsConfig` 不持有步数、主时钟或输出设置；旧 reader 的 `RunInput` 只是输入适配记录。App 一次生成只读 `HydroCaseSchedule`，Hydro 与 EHD 共享同一对象，复用 `RunControls` 并保留旧格式的精确整数步数，不由 `endTime/dt` 反算。模块 IO 接收明确的报告时间、步长和文件前缀，不保存运行配置或派生第二套时钟。普通域暂仍读取旧 `hydro.dict`，通用控制类型的存在不代表支持全部选项；实际输入/节奏与验证见[运行设置记录](validation/organization/QUANTUM_Modular_M7_Hydro_Run_Controls_20260919.md)。

M2 状态查询约束：`Hydrodynamics` 私有持有物理场、历史及工作数组。App 和 IO 使用 `HydroStateView`，不得公开可写指针或经 const Field 改写 pointee；该视图只在下一次模块变更/销毁前同步借用，不自动保持旧时层。需要保留的速度由 `vertexVelocity/boundaryAdvectionVelocity/averagedNativeVelocity` 或拟合域采样器返回拥有型结果；节点/场输出使用显式独立副本。Current/Previous 来自真实模块历史，predictor 不能当作更早物理层。求解失败停止服务，已有拥有型结果保持；本批未引入全状态回滚，见[验证与有效期](validation/organization/QUANTUM_Modular_M2_Hydro_Fields_20260916.md)。

M2b 初始场/物性约束：生产 Hydro/EHD 与喷嘴须通过 `initializeState` 和带域/来源时间的 `updateMaterials` 准备状态，再调用 `advance`；不得在运行适配中直接发布初值/历史或退回旧物性入口。冻结物性可跨步复用，时变物性须匹配步起点；AMR 保留这些记录并拒绝尚未更新到完成子步的物性。物性计算失败停止实例，不能换入口或经 AMR 绕过。原初始速度相容投影、物性公式、NS、字典和门限保持。旧公开字段及运行器债务见初始化记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M2b_Initialization_20260916.md`）；旧接口存在不授权恢复生产直接写场。

M2b 联合修正规范：完整域与拟合域共用 `HydroCoupledBalance` 的块装配和完整动量残差；原压力/散度行策略负责系数与载荷，不再在域驱动中复制联合装配或范数。压力基准必须显式指定；开放压力条件不删除连续性行，粗细面不强制假设 `C=-D^T`。保留原逐面运算顺序及 Eigen 稀疏表达式求值，不以中间向量替换为纯接口整理；模块原验收门限和 SuperLU 不变。该服务不拥有物理状态、主时钟、字典或求解缓存；构建、对照及未完成职责见装配记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M2b_Balance_20260915.md`）。

M2b 几何约束：选中的完整叶域、原生面几何、XY/RZ 测度和面图/网格身份校验归 `Core/Mesh/NativeLeafDomain`，复用原网格与图服务，不创建另一张网格或第二份面数组。`HydroFittedDomain` 组合该对象；压力口、无滑移/回流闭合、速度行消去、虚拟子单元对偶和方程矩阵仍归 Hydro，不搬入 Core。逐面查询不隐式扫描全网格；网格可能改变后，批量复用前必须核验精确身份，尺寸检查不能替代身份。不得在 Hydro 重新实现已迁移的几何构造，也不得让 Core 推断物理边界条件。检查、原生矩阵对照及喷嘴复测见几何职责记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M2b_Geometry_20260915.md`）。

已实施：具体 App 离开 `Runtime/Apps`；EHD 耦合类与普通诊断 writer 提取；喷嘴毛细/介电/导电力装配与文件/控制台输出拆分；原电学专用 writer 归 Electricity/IO；旧喷嘴输入集中在 Compatibility；CMake 路径与依赖检查同步。

以下 A2 职责拆分已实施，并完成 SuperLU 构建与结构检查；两条原喷嘴均完成 20 步、退出 0，导电 76 个、旧介电 86 个完整数值文件分别与 A1 逐字节一致。证据及未完成项见 A2 记录（未迁入；原路径 `validation/organization/QUANTUM_App_Organization_A2_20260912.md`）。

| 原存量位置 | 当前归属 |
| --- | --- |
| `ElectrosprayTimeLoop.cpp` 的区域、材料绑定 | `Setup/ElectrosprayHydroDomain` 与 `Setup/ElectrosprayMaterials`；边界、材料和历史输出共用原机械归一化值对象 |
| 同文件的场/相通量输出和三个 CSV 历史流 | `IO/ElectrosprayHydroDiagnostics`；驱动按原先位置显式调用，仍由驱动接受状态 |
| `Setup/ElectrosprayElectrostatics.cpp` 的输出 | `IO/ElectrosprayElectricDiagnostics`；取面值仍共用 Setup 中原拟合边界距离接口，不再求解电场 |
| `Setup/ElectrosprayInterface.cpp` 的初始输出 | `IO/ElectrosprayInitialInterfaceDiagnostics`；通过明确的只读观察接口，在 Snapshot 构造前写实际节点，在初始化提交前写面元 |
| `Coupling/ElectrosprayInterfaceMotion.cpp` 的运动输出 | `IO/ElectrosprayInterfaceMotionDiagnostics`；输出与计算调用同一个链验证函数；体积修正进度 writer 由驱动传入 |

初始化观察接口是为保留失败诊断设置的同步借用，不是第二个时钟或状态发布接口。不得在回调中修改几何、保留借用数据或提前标记接受步。候选 CSV 仍可能先于耦合接受写出；接受历史、Tecplot 帧清单与实际退出码共同判断运行状态。初始几何绑定/网格准备仍由 App 调度已有 Core 能力，不因本次职责拆分重写公共几何算法。

随后按组织方案 B 批实现版本化字典、新旧等价配置与求解器实际消费；Runtime 的存量 Plasma 数据装载属于 D 批，不随本次改写碰撞或粒子算法。新增代码应直接遵守前四节，不复制上述存量混合形式。

B1 已新增数值参数的独立调用接口：`HydroFittedNumerics` 和 `DielectricFluxChecks` 归各自模块；喷嘴 Setup 仅组合它们与跨模块检查参数；旧 reader 展开原值；IO 独立输出来源与消费者。构建、直接接线检查及三个 20 步入口完整对照均已完成，证据见 B1 记录（未迁入；原路径 `validation/organization/QUANTUM_Dictionary_B1_Numerical_Controls_20260912.md`）。这一步不启用新字典语法，旧共享 `PoissonCellCentered` 选项仍明确标为喷嘴未使用，不能把它称为 B 批全部完成。

B2 已接入版本化喷嘴的严格读取、公共 SI 材料表、模块各自配置、五个真实线性系统、输出/日志消费者、生效配置记录和只读几何绑定检查。旧 Hydro lexer 通过别名保留原行为，不全量替换旧 geometry/PIC 解析器；原兼容喷嘴仍可运行。具体实现、归档 EXE、原短程复测与限制见 B2 记录（未迁入；原路径 `validation/organization/QUANTUM_Dictionary_B2_Versioned_Input_20260912.md`）。

新旧配置等价检查和 CLI 预检查分别为 `electrospray_versioned_input`、`electrospray_versioned_cli`；实际数值参数调用继续用 `electrospray_numerical_controls`。结构变更仍需 `architecture_dependencies`，八项例外没有扩大。schema 迁移比较须核对生效输入和原生数值文件，不能放宽原同输入比较器。Runtime 的历史选择探测和 Plasma 装载、App 的只读兼容绑定对象尚未全部退役，后续按实际触及范围处理，不因此另拆一轮阻塞原喷嘴推进。

M1a（2026-09-15）进一步收紧输入：新版 reader 保留模块拥有的类型，几何身份检查、配置报告和 Hydro 压力/回流材料绑定直接消费这些类型。新增代码不得再让这些入口经过 `ElectrosprayPhysicalConfig`；剩余旧签名的临时转换只放 `Compatibility/ElectrosprayInputAdapters`，不在各调用处另做一套映射。七个迁移边界已加入依赖检查，检查自例由 16 增至 22，原八项 Runtime 例外不变。实现、实际对照状态与仍待退出的消费者以 M1a 记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M1a_20260915.md`）为准；这不表示 M1b 或模块独立运行已经完成。

M1b（2026-09-15）完成正式喷嘴的配置消费收口：Setup/受力/场输出只接收实际所需的模块配置、材料或尺度；时间循环、电荷和 AMR 接线直接组合模块输入，不再接受 `ElectrosprayConfig/PhysicalConfig`。主时钟及输出只由 `RunControls` 提供，`ElectrosprayWorkflow` 仅表达 App 模式、初始几何选择和旧日志行为的显式兼容信息；不得向其中重新添加物性、力带宽或求解容限。旧输入单向适配到同一调用链，旧结构仅供旧 reader 与来源报告使用，不恢复反向构造函数。线性系统值类型与字典 reader 分离；34 个迁移边界受依赖检查约束。实际原输入对照及仍未完成的模块服务归属见 M1b 记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M1b_20260915.md`）。

M2a（2026-09-15）新增公共 `Hydrodynamics::advance`：完整域 XY/RZ 与拟合域 RZ 的生产驱动均使用明确的单步输入；旧重载只作现有数值回归/旧几何适配，不恢复生产分支。载荷标明位置、代码单位下的力密度及来源时刻；一个 Hydro 实例核验步号/时间连续性，AMR 同时转移此身份，不自行推进主时钟。数值失败须报告物理/历史场是否可能改变，并停止该实例，不能把成功 Hydro 子步当成完整耦合接受。拟合域材料赋值及顶点闭合归 Hydro，App 只做分类/尺度绑定；精确网格身份归 Core。优先使用不暴露可写指针的 `HydroStateView`；既有公开字段与运行器迁移债务、原输入对照见 M2a 记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M2a_20260915.md`）。

M2a 补充：同一 Hydro 对象的公共 `advance` 与旧 `solve*`/cache 入口不得混用；约束在进入内核时固定，AMR 必须保留，数值失败不能通过换入口重试。请求预检查失败不占用入口。喷嘴时间循环、AMR 与流场 writer 已改用 `state`、具名当前/历史速度的 `fittedVelocity` 和 `fittedPhaseOutwardFlux`；新增消费者不得直接读取公开场数组或另组一套重构。借用视图与拥有型速度快照的有效期分开说明；原生数组已由后续[私有化批](validation/organization/QUANTUM_Modular_M2_Hydro_Fields_20260916.md)封闭，普通域消费者也已迁移。实际实现、原输入复测及下一段 M2b 范围见入口与状态记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M2a_State_20260915.md`）。

M4/Hydro 原生加载后续：App 的毛细装配统一调用 Hydro 服务，来源在成功后绑定，失败工作区必须丢弃；原生入口复用现有加载几何，均匀入口保持原算式。Electricity 校验真实分裂面关联，Hydro 的粗细面标量一致性检查消费共享虚拟子单元模板，实际力仍来自原压力矩阵。常量插值舍入不以两端相同就强制为零，检查扩展须单列；Interface 显示值不得用旧节点数组代替。静态细化介电、解析参考和原数值对照范围见[加载记录](validation/organization/QUANTUM_Modular_M4_Native_Capillary_20260917.md)，不代表动态 AMR 或一般二维粗细电学收敛。

M2 Hydro 接受生命周期：耦合 App 调用 `prepareStep`，完成 Interface/电学/材料及收支检查后调用候选的 `accept(mesh)`，之后才发布新时间输出。未接受候选必须阻止下一步和 AMR 捕获；放弃候选停止模块，不能用子步成功报告冒充接受或恢复旧状态。只读查询在待接受期读取实际候选层，拥有型采样保留原生命周期；销毁拥有者后候选不得再接受。`advance` 保留独立调用者的同步接受行为。数值内核和主时钟不归候选对象，不增加整场副本，当前范围见[候选记录](validation/organization/QUANTUM_Modular_M2_Step_Candidate_20260917.md)。

M2 动量装配生命周期：矩阵条目和右端向量由本次 `HydroMomentumAssembly` 持有，原生和拟合行策略、边界及载荷显式写入；完成后预测求解与联合校正消费同一份矩阵和右端向量，不恢复求解器的跨步装配成员。结束过程保留调用方稀疏存储方式、条目合并和运算顺序，不借结构整理更换合法离散公式。装配对象不拥有物理状态或接受权；验证范围见[装配记录](validation/organization/QUANTUM_Modular_M2_Momentum_Assembly_20260917.md)。

M6 PIC 输入边界：App 或独立调用者先装载 `PicCaseInput`，两种 PIC 驱动按值消费准备好的配置记录；不在物理运行内部读取系统字典或从输出目录推断输入是否存在。诊断来源、几何文件存在标记与实际读取结果分别保留，检查点要求消费准备输入。原显隐缺省/异常策略由共同 IO 服务保留；不恢复隐式 reader 空包装或第二套语法。原数值与恢复对照见[PIC 输入记录](validation/organization/QUANTUM_Modular_M6_Pic_Case_Inputs_20260917.md)。

M6 PIC 单步边界：App 在 `withSession` 的同步有效期内驱动 `PicStepSession`；模块保留原初始化、整步和状态历史，不能重新开一层宏观循环。步号来自恢复范围，成功后仅推进一次，失败或异常后停止；驱动提前返回不能声称完成。会话不可复制/移动，不得逃逸或当作长期快照，原步内输出/检查点顺序仍保持，拆分剩余 IO 时不得改物理格式。检查范围见[单步会话记录](validation/organization/QUANTUM_Modular_M6_Pic_Step_Session_20260917.md)。

PIC 会话使用恢复后的有效 `dt` 和步号生成 `PicStepTime`，显式、隐式整步及后续完成阶段借用同一代码单位区间；起止时间沿用整数步号乘 dt，不累计浮点时间，不改变原 SI 换算顺序或粒子子步。`advance` 返回可保留的值类型 `PicStepResult`，区分 MacroStep、Completion、Accepted；异常向外传播，同时在 `lastStep()` 记录阶段和 `threw`，不能以零状态码单独判断成功。完成阶段失败时接受步号/时刻不前进，但物理字段和已写输出可能改变；MacroStep 内仍有原交错 IO，不宣称这两个阶段已分离全部物理与输出职责。非法步数/时间区间在调用会话驱动前拒绝；已越过结束步的恢复仍为空运行。

M6 PIC 检查点编排：App 按准备输入的开关、频率和末步请求 `advance(saveCheckpoint)`，数值驱动不得复制自动调度策略。原物理步成功后才执行完成回调，完成失败也停止会话，步号不得提前增加。`PicCheckpointWriter` 统一原公共写入顺序，借用实时 RNG，物理组件继续使用各自 writer；禁止增加隐式状态快照或重复格式实现。元信息/RNG 及组件写入失败必须传播；未完成目录不能选择，也不计入保留数量。组件的打开、写入、刷新和关闭都成功后才移除标记，随后清理旧检查点；几何要求与实际构造条件一致。会话停止但不回滚物理步或同一步目录，不宣称原子恢复/掉电持久性，见[检查点完成记录](validation/organization/QUANTUM_Modular_M6_Pic_Checkpoint_Completion_20260917.md)。

M6 PIC 检查点接线：`PicCheckpointStore::forRun` 持有准备输入导出的路径/组件要求，保留显式水平集名字门及隐式场历史要求；构造不访问文件。显式 `prepareExplicitRestart` 在输出预检阶段选择一次目录，返回私有持有目录和输出恢复步号的 `PreparedExplicitRestart`；稍后完整元信息读取和 dt 覆盖必须使用该对象，不重新选择 latest。宽松输出解析、碰撞预检、完整解析和 dt 生效仍在原阶段；准备后出现新检查点不能改变恢复来源，已选目录缺失或出现未完成标记时停止，不回退。这只固定同一工作目录下的来源路径，不锁定缓存内容或回滚已处理的输出。`PicCheckpointWriter::LiveState` 在写出时同步借用实时运行域、壁面、物种、碰撞及可选隐式历史，按原顺序调用组件 writer，不保存这些引用或在早期冻结碰撞指针。不得为统一接口改变主 RNG、组件格式或拒绝恢复时的输入缓存行为，原阶段依据见[接线记录](validation/organization/QUANTUM_Modular_M6_Pic_Checkpoint_Wiring_20260918.md)。

M6 PIC 恢复读取：`PicCheckpointReader` 拥有一次选定的目录与元信息，在原阶段恢复实时网格、RNG 和粒子，不持有物理状态副本，不重新选择更新目录。缺失网格的显隐差别、元信息阈值、RNG 警告和后续物种失败时的部分读取语义保持；失败调用者必须停止，不能声称回滚。几何、壁面和场历史仍按各自初始化依赖恢复，不为合并函数改变顺序，见[恢复读取记录](validation/organization/QUANTUM_Modular_M6_Pic_Checkpoint_Reader_20260917.md)。

M6 PIC 粒子初始化：`initializePicSpecies` 填充按配置身份新建的空物种，直接使用调用者的 RNG 和网格；尺度换算、采样顺序、正态分布的逐物种生命周期及恢复时不抽样的语义保持。显隐壁面排除策略由各调用者明确提供，不能为合并采样改变容差、把剔除改成重采样或新增随机流。物种分类和诊断保持原顺序，见[粒子初始化记录](validation/organization/QUANTUM_Modular_M6_Pic_Particle_Initialization_20260917.md)。

M6 PIC 场文件：`PicFieldOutput` 的底层 writer 只序列化已准备的坐标/变量/掩膜，不插值、不做单位换算或保存物理副本。显式阶段入口同步借用当前输入，调用已有帧/剖面服务准备数据后写出；不另立数值实现。三种原顶点布局、路径、排序容差及连接规则分别保持，不能为共用 writer 改变原显示值或字段顺序；队列与 AMR 刷新仍遵循调用方原时序。直接 IO 见[场输出记录](validation/organization/QUANTUM_Modular_M6_Pic_Field_Output_20260917.md)，显式阶段见[显式输出记录](validation/organization/QUANTUM_Modular_M6_Pic_Explicit_Output_Stage_20260918.md)。

M6 PIC 待写队列：`PicFieldOutputQueue` 拥有准备字段、掩膜、坐标及隐式输出配置，借用当前网格和不可变单位上下文；物理输入只在同步准备帧时借用，不进入待写任务。调用者必须在网格变更前和正常结束时刷新，AMR 重建后在原位置更新坐标；错误网格或未刷新任务先拒绝，不能让旧帧引用新网格。立即路径不入队，批量路径移动字段，完整写出后才清空；部分失败保留原已写文件及待写帧并由运行方停止，析构不得补写或重试。原队列与失败契约见[输出队列记录](validation/organization/QUANTUM_Modular_M6_Pic_Output_Queue_20260917.md)，配置/坐标及完整隐式准备阶段见[隐式输出记录](validation/organization/QUANTUM_Modular_M6_Pic_Implicit_Output_Stage_20260918.md)。

M6 PIC 输出重建：`PicOutputReconstruction` 共用原掩膜插值并处理已准备显示场的周期顶点配对，不修改物理单元场或模拟边界条件。输入采用当前网格索引；仅大小匹配的掩膜生效，稀疏邻域回退、限幅、最终填充值及 Y 后 X 的周期次序保持。不能把显示用重建替换粒子取场或物理求解，见[场重建记录](validation/organization/QUANTUM_Modular_M6_Pic_Output_Reconstruction_20260917.md)。

M6 PIC 粒子轨迹：`PicParticleTrajectories` 拥有运行内编号计数、轨迹批次和样本；配置归 `PicParticleTrajectoryConfig`，旧 `OutputConfig` 名称仅作类型别名。每次借用粒子和模拟 RNG，选样顺序及调用阶段不得因诊断拆分而改变。启动清理和正常结束写出显式调用，析构不补写；运行编号和样本不写入现有检查点，保留原续算重新建轨迹的语义。CSV/轨迹标记共用 `TimeSeriesIO::sanitizeCsvToken`，见[轨迹记录](validation/organization/QUANTUM_Modular_M6_Pic_Trajectories_20260917.md)。

M6 PIC 物种与边界 CSV：`PicSpeciesSeries` 私有拥有列标记及 primary/detail/budget、边界计数/能量/通量六种文件流，调用者不借用列名或流。诊断数组按调用借用，收支和边界记录保留原整数/浮点类型，按物种及边界顺序逐项写入；计算方仍负责原残差、时间窗口、通量及单位换算，IO 不重算守恒量。续算步、零通量筛选和失败后的调度由驱动按原规则提供，不能因流归属变化合并显隐的刷新、门控或浮点格式差异。见[物种 CSV 记录](validation/organization/QUANTUM_Modular_M6_Pic_Species_Series_20260917.md)及[边界 CSV 记录](validation/organization/QUANTUM_Modular_M6_Pic_Boundary_Series_20260917.md)。

M6 PIC 逐步汇总：`PicStepSeries` 拥有总览、过程与几何耦合 CSV，接收准备好的标量和同步借用名称，不复制诊断状态或重算能量/几何统计。显隐布局、计数类型、缓冲、续算步及原外层门控分别保留；不能为了合并格式引入字段或隐式刷新，见[逐步汇总记录](validation/organization/QUANTUM_Modular_M6_Pic_Step_Series_20260917.md)。

M6 PIC 逐步诊断：`PicStepDiagnostics` 计算过程主导项和运动几何统计，借用原物种/反应累加量与只读边界控制器，不复制物理状态、不更改 RNG 或写盘。显式累计窗口与隐式单步由调用方提供，不为共用诊断重排通道、累加或活动区查询。选定名称由结果拥有，IO 同步借用；原 alpha 阈值、零面积回退和并列规则保持，见[逐步诊断记录](validation/organization/QUANTUM_Modular_M6_Pic_Step_Diagnostics_20260917.md)。

M2 压力方程生命周期：`HydroPressureSystem` 只拥有一次完成装配的矩阵/右端。完整域的 `HydroBoundaryProjection` 同时持有本次面连接、面积、逆质量和固定速度；方程、校正及通量核算消费同一只读结果，不重新查询活动几何。XY/RZ 度量和密度插值由原适配准备，稀疏条目/求和顺序保持。初始投影采用尺度 1，保留已给定的固定面值并在全部检查通过后发布；乘子不写入物理压力、不占用压力分解缓存或时间历史。时间步采用原 dt，强制固定面值，按原顺序写投影势和速度；不承诺失败回滚。结果仅限当前同步调用，不能跨步/AMR 复用；只增加局部面方程记录，不复制网格或物理场。拟合/原生自适应投影保留各自压力基准和缩放，不据共享存储替换合法离散。见[边界投影记录](validation/organization/QUANTUM_Modular_M2_Boundary_Projection_20260918.md)。

M2 单元压力行：完整域 XY/RZ 共用 `assembleCellPressure` 的邻域、插值、基准和条目构建；坐标策略只提供原面系数与通量表达式，不另存网格字段。保留 XY 倒数纵横比和 RZ 东西面分别相乘的浮点顺序，不能为表面统一改写算式；原生自适应面装配及拟合/边界投影按各自合法策略保留，见[压力行记录](validation/organization/QUANTUM_Modular_M2_Pressure_Rows_20260917.md)。

M2 原生速度校正：`correctNativeVelocity` 共用压力面校正、速度到顶点插值及单元梯度循环。坐标分支只保留原 XY 边界筛选和 RZ 面积/体积散度；不得改变压力面梯度、密度除法、历史发布或旧缓存梯度幅值的计算顺序。受管理状态仍不写未版本化旧缓存，见[速度校正记录](validation/organization/QUANTUM_Modular_M2_Velocity_Correction_20260917.md)。

M2 原生动量矩阵：`buildNativeMomentumMatrix` 共用两种坐标的邻域、物性插值与条目构建。`HydroNativeMomentumDetail` 仅提供逐行局部度量及应力系数计算，保留 XY 倒数比例、RZ 先乘面积后除长度及径向黏性两项之和的运算次序；不得预先合并交叉应力比例。RZ 规则行、行掩膜和递归叶面策略显式保留，不能将坐标合并当成 XY 任意深度细化的等价证明，见[动量矩阵记录](validation/organization/QUANTUM_Modular_M2_Momentum_Matrix_20260917.md)。

M2 原生动量右端：`buildNativeMomentumRHS` 共用 CN/AB2 的邻域、当前/历史取值、迎风和受力累加；原生面面积由矩阵与右端共同调用，保持虚子单元半径。XY 黏性通量先作差再乘面长，RZ 径向两项分别乘面积；不能将这两个表达式代数统一。边界值、行掩膜清零、规则行调用和借用力来源按原坐标/域策略保留，不改变步状态或接受顺序，见[动量右端记录](validation/organization/QUANTUM_Modular_M2_Momentum_RHS_20260917.md)。

M2 原生旧路径退出：XY/RZ 与拟合内部行直接调用 `buildNativeMomentumMatrix/buildNativeMomentumRHS` 的原坐标特化，不恢复四个旧转发函数或已注释掉的装配实现。公共 XY/RZ 求解入口统一放在 `Hydrodynamics_impl.h`；原生 RZ 规则行和叶面遍历归 `Detail/HydroNativeRZRowPreparation.h`，可独立包含。完整域失败即停与拟合域候选发布的差异仍由各自策略负责，不能据此声称全状态回滚，见[退出记录](validation/organization/QUANTUM_Modular_M2_Native_Retirement_20260917.md)。

M6 PIC 粒子场准备：`PicParticleFields` 同步借用粒子、沉积模板及归一化参数，只写调用方提供的显示场，不创建 RNG、缓存或修改物理状态。CSR 沉积保留粒子/条目顺序，数密度保留切割体积和 SI 换算；电子能量保留显式有限正权重、隐式回退后只查正值及宿主单元取样策略，不能以共用为由替换。模板必须对应当前粒子顺序和网格，输出零初始化且不与输入别名；选择、命名、缓存生命周期与调度仍由驱动负责，见[粒子场准备记录](validation/organization/QUANTUM_Modular_M6_Pic_Particle_Fields_20260917.md)。

M6 PIC 输出准备：`PicOutputReconstruction` 统一合法四叉树上的活动掩膜和标量转换循环，切割活动区优先于面图，alpha 仅作额外过滤。调用者提供单位规则并保留显式“顶点插值后换算/截断”和隐式“单元换算/截断、显示填补后再插值”的顺序；不能代数交换这些步骤。隐式零势区域标记及各显示标量的填补保留顶点邻接重复计数、原地遍历和 20 次上限，只写准备好的场，不修正物理电势；掩膜与字段必须对应当前网格，见[输出准备记录](validation/organization/QUANTUM_Modular_M6_Pic_Output_Preparation_20260917.md)。

M6 PIC 字段规则：`PicFieldSelection` 只持有输出元数据，统一别名/默认/命名和隐式稳定排序，不接字典 reader 或物理场。显式保留空列表默认、始终输出物种密度和固定字段顺序，隐式保留显式选择与用户列表顺序；不得为共用而新增不支持的字段。`PicFieldKeys` 每个物种输出组重置基础名计数，保留两种符号拼写和原后缀碰撞语义。重排只移动名称及对应字段，数量不匹配仍留给原 writer 校验。实际输入规则归字典参考第 8 节，验证见[字段选择记录](validation/organization/QUANTUM_Modular_M6_Pic_Field_Selection_20260917.md)。

M6 PIC 输出帧：`PicFieldFrame` 拥有准备后的名称、顶点/单元值和有效掩膜，组装期间借用当前网格与只读选择；待写任务只持有准备结果，不保留物理场或粒子引用。共同粒子字段组织保留显式密度/能量分别构造缓存、隐式单次缓存复用及快速模式的 Nstar/宿主单元策略；密度回退到 CSR 不得顺便改变快速能量来源。显式插值后换算/截断，隐式换算/填补后插值或直接写单元，不能为共用而交换顺序。隐式字段/单位门控、缓存策略与完整帧准备归输出队列；运行方通过共同输出调度选择正的完成步号后调用，队列不再持有另一份场频率。驱动保留物理时层、AMR 和结束刷新。显式阶段归 `PicFieldOutput::writeExplicit`，保留 Profile1D 先写、字段固定排列与立即提交；驱动仍在初始半步后、原频率及末步调用，平均剖面留在求解阶段。帧内规则见[输出帧记录](validation/organization/QUANTUM_Modular_M6_Pic_Field_Frame_20260917.md)，完整隐式阶段见[隐式输出记录](validation/organization/QUANTUM_Modular_M6_Pic_Implicit_Output_Stage_20260918.md)，显式阶段见[显式输出记录](validation/organization/QUANTUM_Modular_M6_Pic_Explicit_Output_Stage_20260918.md)。

M6 PIC 电学/几何输出：电势/电荷单位换算、电场分量拆分、水平集单元采样及隐式极化系数代理归现有输出帧；只同步借用当前物理数据，不保留几何/状态引用。`muFace` 保留逐有向面段对 K/L 的累计，内部面两向条目都参与，不能改成面积加权或只取 owner 范围。可选掩膜诊断归 `PicFieldOutput`，保留原列名、原始坐标/电势和默认精度；文件打开失败仍静默，目录异常仍传播。此诊断的 x_m/y_m 旧列名不意味着自动转换，正式 packed fields 仍由输出帧换算为 SI。见[验证记录](validation/organization/QUANTUM_Modular_M6_Pic_Geometry_Fields_20260917.md)。

M6 PIC 剖面输出：`PicProfiles` 拥有瞬时文件和平均窗口；同步 `Sample` 只借用调用时的场、单位和物种，不跨步保留引用。显式平均消费求解缓存，瞬时输出重建诊断缓存；隐式两者在原求解阶段消费缓存和 EcellNew，禁止为接口一致改变采样时层。恢复准备瞬时 CSV 时按原规则删除边界步及其后行，平均窗口从空开始；未就绪时保持旧平均快照文件。字段输出频率统一传入，缓存数量不符仍跳过平均累积和写出。诊断耦合的位置/宿主提取现归 `PicCurrentWeights`，按次复制选项并替换临时 hostCellLabels，不改变调用者选项或求解缓存。原缓存调试策略和剖面数学核保持，见[剖面输出记录](validation/organization/QUANTUM_Modular_M6_Pic_Profiles_20260917.md)。

M6 PIC 边界诊断：底层边界测量只借用当前场/网格并返回代码单位标量及自有调试文本；显式任意配对、隐式切向对齐配对分别声明，不能为共用服务改变切割过滤、面段顺序或采样时层。周期 CSV 归既有 `PicStepSeries`，电极/电功的 SI 采样及历史归下述诊断对象，驱动保留原门控及阶段。见[边界诊断记录](validation/organization/QUANTUM_Modular_M6_Pic_Boundary_Diagnostics_20260917.md)。

M6 PIC 选定电极电流：`PicElectrodeCurrent` 拥有电极身份、同步场样本和差分历史，只发布只读标量结果；不持有网格/粒子引用，不从波形 CSV 或旧检查点推断历史。边界电压按调用者传入的 SI 时刻求值；无单位时显式保留 E0 乘法（包括零），隐式保留原场值。写出开关及频率归 `PicStepSeries`，不得让文件可用性改变其他诊断中的电流可用性。四侧功积分和波形出流口径不得擅自合并，见[电流波形记录](validation/organization/QUANTUM_Modular_M6_Pic_Current_Waveform_20260917.md)及[采样记录](validation/organization/QUANTUM_Modular_M6_Pic_Boundary_Sampling_20260918.md)。

M6 PIC 四侧/嵌入电功归 `PicBoundaryWork`，仅持有 SI 样本、四侧历史和标量结果；法向场、四侧电压及嵌入电压在原阶段分别同步采样，服务不保存借用引用。显式按原 SI 时刻求值，隐式读取已更新的边界存储值；不合并两条路径的门控或求值时层。保留四侧顺序、嵌入 ID/默认选择、原电荷映射遍历及零/非有限电荷过滤。域边界电功与选定电极波形的注入和零 dt 历史规则分别保持，见[边界电功记录](validation/organization/QUANTUM_Modular_M6_Pic_Boundary_Work_20260917.md)及[采样记录](validation/organization/QUANTUM_Modular_M6_Pic_Boundary_Sampling_20260918.md)。

M6 PIC 序列窗口归 `PicDiagnosticWindow`，私有拥有显式窗口的源汇、边界、能量及起始状态，读取者只借用 const 结果；数据长度和物种下标遵守初始化后的固定物种契约。`PicSpeciesBudget` 的累计/残差计算归诊断，IO 只序列化；隐式单步统计不得隐式改为多步窗口。见[统计窗口记录](validation/organization/QUANTUM_Modular_M6_Pic_Diagnostic_Window_20260917.md)。

M6 PIC 能量收支：`PicPowerWindow` 私有拥有显式报告间隔的七项传递累计、场/总能量历史及时间；`PicStepEnergyHistory` 拥有隐式逐步场/动能历史。显式按原场输出采样且报告后清空传递项，隐式按每个原诊断步及当前 dt 更新；有限性门、求和和电荷分类规则不得合并改写。IO 的 `PicStepSeries` 只序列化已准备值，保留三套列布局、格式、恢复裁剪及原刷新点。见[能量收支记录](validation/organization/QUANTUM_Modular_M6_Pic_Energy_Balance_20260917.md)。

M6 PIC 粒子/电荷收支：`PicStepParticleBalance`（原 PicStepChargeBalance）拥有隐式步首物种数、三类保留生成量/能量与阶段/物理电荷闭合，只发布 const 结果。两个 MCC 半步直接逐粒子累加同一物种和，不改浮点分组；MCC 只计数量，SEE 在原提交前核算并在提交后合入源电荷小计，普通源在原提交后核算。总注入能量为源加 SEE，保留原三维 SI 动能函数和权重回退；未知物种跳过。这里统计原阶段的保留缓冲，后续提交修复/拒绝、固体剔除及粒子管理不改定义。物种末态仍在原报表门内采样；新步重新构造，重启不覆盖历史诊断。原 OpenMP 归约、吸收、电荷归一化和 CSV 布局/精度不变。显式电流纯汇总仍只在原输出门通过时消费，IO 不计算物理收支。见[隐式粒子收支记录](validation/organization/QUANTUM_Modular_M6_Pic_Implicit_Particle_Balance_20260918.md)。

M6 PIC 边界初始化归 `Boundary/PicBoundarySetup`：共同场映射消费已准备的配置/几何，仅返回拥有型结果；不读文件、不保留输入引用。显隐共享线归属、域边位置/名称推断、波形比较和冲突检查；模式只选择原警告/错误文本。显式额外周期配对和域 SEE 在原启动位置完成，隐式 SEE 保留恢复后的原调用位置与模板选择。运行期时变 BC 和墙面物理继续消费原配置，不能把设置阶段映射当成物理推进。见[边界初始化记录](validation/organization/QUANTUM_Modular_M6_Pic_Boundary_Setup_20260917.md)。

M6 PIC 嵌入几何归 `PicEmbeddedGeometry`：选择、组合场、CutCellMetrics、物理 alpha 和边界标签由一个不可复制的模块对象私有拥有，只发布 const 视图。单水平集视图借用现有 ImmersedGeometry，组合视图由服务拥有；几何或网格变化后须按原顺序重建，视图不得跨重建/重选/销毁使用。原活动单元余量、alpha 优先筛除、移动线坐标变换分别按模式保留，不能以共用为由改变物理区域。原始几何与 Refiner 由 `PicDomain` 持有；运动/恢复时机仍由原驱动编排，不宣称完整几何生命周期已合并。见[嵌入几何记录](validation/organization/QUANTUM_Modular_M6_Pic_Embedded_Geometry_20260917.md)。

M6 PIC 的网格、细化配置副本、细化器、原始/派生几何和面图统一归不可复制/移动的 `PicDomain`（替代 `PicGeometry`）。成员顺序保证细化器/派生视图先于原始几何和网格销毁；场、粒子及运动服务只在所属会话内借用。网格对象在恢复/AMR 后地址不变，叶数组、几何派生量和求解缓存仍按原阶段转移、重建及失效。构造不建面图；启动细化后和原 AMR 阶段显式重建，不承诺自动追踪任意可变借用。几何仅初始化一次，失败停止会话；空几何、显隐恢复差异和规定运动交接保持。细化器借用的几何字典仍由调用者保证寿命。见[运行域记录](validation/organization/QUANTUM_Modular_M6_Pic_Domain_20260918.md)。

PIC AMR 的几何刷新策略、适配前标量诊断和选择性重初始化连续次数归 `PicAmrGeometry`，驱动只发布实际适配上下文并读取 const 报告。必须在 mesh.adapt 前 capture，在原 MovingEB 网格通知之后 refresh，随后按原顺序重建 face graph、cut metrics、物理场并使求解缓存失效。诊断只复制原有按水平集的小型统计，不新增网格/场副本；报告不得跨下一次 capture/refresh 使用。连续次数仍为运行内历史，检查点不新增该字段。原始几何的寿命由 `PicDomain` 覆盖；完整 AMR 联合提交仍由原消费者协调。见[AMR 几何记录](validation/organization/QUANTUM_Modular_M6_Pic_Amr_Geometry_20260917.md)。

规定几何运动由 `PicScheduledGeometry` 持有计划、当前观测状态和同步事件观察函数；几何重建与求解缓存失效完成后接受运动。网格、原始几何、派生几何、面图、细化器和求解器采用有明确生命周期的借用绑定，服务存活期间对象身份须稳定。未更新步返回已接受几何对应的观测状态并置 changed=false；强制刷新即使 changed=false 仍重建。恢复启动使用已恢复状态，不重复移动。失败沿用中止语义，不承诺部分几何回滚。trace 过滤/文件/控制台及几何指标 CSV 归 `PicStepSeries`，观察函数借用的报告对象必须覆盖运动服务生命周期，不能推进或修改几何；报告打开及同步事件位置保持。见[规定运动记录](validation/organization/QUANTUM_Modular_M6_Pic_Scheduled_Geometry_20260917.md)及[运动诊断记录](validation/organization/QUANTUM_Modular_M6_Pic_Geometry_Reports_20260918.md)。

PIC 的 MovingEB 控制器、按步撞壁查询上下文和最近更新标记归 `PicMovingGeometry`，驱动仅通过 const 控制器视图提供派生几何/诊断读取。beginStep 清理撞壁上下文及控制器计数，保留上次 updated 供原时刻诊断使用；advance 发布本次结果。隐式粒子剔除须在 advance 与 refreshAfterUpdate 之间。精确移动后刷新 AMR 参考线，平流后清除陈旧参考线；失败中止，不新增回滚。残差诊断必须使用本次求解实际采用的边界表。见[MovingEB 记录](validation/organization/QUANTUM_Modular_M6_Pic_Moving_Geometry_20260917.md)。

PIC detail 流及初始物种列布局归每次运行的 `PicStepSeries`，不保留进程级静态文件缓存；`PicDiagnostics` 只计算结果。恢复追加前应验证可观察的表头/单位布局，拒绝不兼容表头时不得先截断该文件。历史中已丢失的物种值不能凭当前状态补造。见[detail CSV 记录](validation/organization/QUANTUM_Modular_M6_Pic_Detail_Series_20260917.md)。

MCC 的应用计数归实际算子，不能由 PIC 步号推导；检查点 IO 只读写其版本化状态，模式/种子验证成功后才恢复。重写同一步目录先使旧碰撞组件失效，写完前保留未完成标记。通用直接求解缓存须同时核对完整稀疏结构及系数，重新分解失败不能保留有效因子标记。见[碰撞续算记录](validation/organization/QUANTUM_Modular_M6_Pic_Collision_Restart_20260917.md)。

M6 PIC 物理场：`PicFieldState` 绑定会话网格，私有拥有五个公共场；`ImplicitPicFieldHistory` 借用该拥有者并持有下一步/半步电场，寿命不能越过场与网格。两条 PIC 数值驱动均只读电势、电荷量、密度、RHS 和电场，不再拥有场状态/电场历史的友元权限。写入通过沉积、求场和历史模块操作完成；AMR 后的数组借用必须重新取得。半步使用当前与下一步的原平均，接受只复制下一步值，不能交换存储或覆盖半步。快速 AMR 电势必须携带当前网格大小和标签；无效输入在修改前拒绝。缓存读取仍遵循原 IPF1 的头拒绝/部分读取语义，读取失败后的重建由原驱动执行，不声称事务回滚。见[场状态记录](validation/organization/QUANTUM_Modular_M6_Pic_Field_State_20260917.md)。

M6 隐式 PIC 极化阶段：`ImplicitPicPolarization` 拥有 Q*/N*/rho*、mu/kappa 及沉积工作区，仅发布 const 查询；粒子、轨道缓存、归一化和面图按调用借用。沿用完整单元体积，不能在组织调整中改用切割体积；粒子权重回退、OpenMP 归约分组、域边界控制的均值修正和面调和规则保持。驱动保留边界判断、诊断/计时、Poisson 求解及 AMR 的分段重置时刻；沉积缓存的粒子/网格时层仍由构建方保证，不把大小预检扩称拓扑身份验证。见[极化阶段记录](validation/organization/QUANTUM_Modular_M6_Pic_Polarization_20260917.md)。

M6 隐式 PIC 轨道阶段：`ImplicitPicOrbit` 私有拥有取场、预测端点/宿主/分段及平均权重工作区；查询 CSR 为 const，下一次构建、重置或网格变化后必须更新借用。粒子/几何按调用借用，预测不推进粒子或 RNG。分段预算、积分模式、并行开关和鞘层单位换算现由轨道服务持有配置快照；每次构建接收原调用者的切割/拦截门并绑定当前几何，不持有跨重建的几何指针。patch 强制 fused、分段下限及隐式原始水平集编号保留；轨道权重仍启用原统计，不能套用当前位置物理权重的关闭策略。可选 AMR 清理、初步 CSR 验证与故障停止位置保留。OpenMP 工作块不得假定请求线程全部形成：平均权重与普通耦合权重都按全部逻辑分块执行 parallel for，固定区间及拼接次序，避免嵌套串行化时漏粒子。见[轨道阶段记录](validation/organization/QUANTUM_Modular_M6_Pic_Orbit_20260917.md)及[配置/普通权重分块修正](validation/organization/QUANTUM_Modular_M6_Pic_Orbit_Configuration_20260918.md)。

M6 PIC 当前位置权重：`PicCurrentWeights` 绑定原网格并私有拥有物种 CSR 与复用输入数组，物种按原顺序同步构建；调用者显式安排重建/清空。物种内粒子顺序、宿主标签及显隐几何选项不得混用。诊断使用独立缓存，不能覆盖求解缓存；const 借用不是跨粒子移动/AMR 的快照，服务不自动判断拓扑版本，也不承诺构建异常回滚。见[当前位置权重记录](validation/organization/QUANTUM_Modular_M6_Pic_Current_Weights_20260918.md)。

M6 PIC 耦合配置：同一服务内的 `Policy` 拥有配置快照和一项最小叶宽度；隐式轨道只借用 const 宽度，AMR 在原阶段原位刷新，不自动探测拓扑变化。物理/诊断选项每次从当前几何生成，不能保存跨重建的几何指针。显式保留 metrics/尺寸门及物理调试，隐式保留调用者的 `useCutMetrics`、原始水平集编号和物理统计关闭；诊断保留自己的复合编号、尺寸检查及调试默认值。不得为统一接口抹平这些差异。错误网格的 policy 在更换缓存前拒绝，数值构造内核保持原样，见[耦合配置记录](validation/organization/QUANTUM_Modular_M6_Pic_Coupling_Policy_20260918.md)。

M6 PIC 细化指标：Core 只定义 `MeshRefinementIndicator.h` 中的数据记录；PIC 的 `PicRefinementIndicators` 拥有本会话值，调用细化器时同步提供。禁止恢复全局注册表或跨会话查找。流体 `CurrentFields` 与 PIC 指标不能在同一次调用混合；缺失/错误长度指标仍保留既有警告及忽略规则。提供方负责原网格/时层与重建时机，不能把长度检查宣称为拓扑身份验证。见[细化指标记录](validation/organization/QUANTUM_Modular_M6_Pic_Refinement_Inputs_20260918.md)。

指标更新的外层布尔返回值只表示是否完成本次重算；内部空间平滑通过输出数组返回结果，保持显式 `void`。提取时不得把外层早退替换扩散到嵌套 lambda；本次已恢复原平滑体并核对解析结果，见[平滑修复记录](validation/organization/QUANTUM_Modular_M6_Pic_Amr_Smoothing_Return_20260918.md)。

M6 PIC 细化状态：显式和隐式策略分别归 `ExplicitPicRefinement` / `ImplicitPicRefinement`，对象绑定稳定网格且不可复制/移动，只发布 const 结果。显式每次调用更新 EMA、按原间隔构造指标，适配时转移 EMA；隐式保留原相对变化门控、数量限制、持续性及几何请求语义，字段重置不清空门控历史。样本与配置仅在调用期间借用；查询不得当作跨网格快照。实际适配继续调用 Core，完成标记仍在原几何/场重建阶段之后清理；不新增检查点格式或回滚保证。见[细化状态记录](validation/organization/QUANTUM_Modular_M6_Pic_Refinement_State_20260918.md)。

M6 PIC 场边界：`PicFieldBoundaryState` 私有持有域边、嵌入默认和按 ID 边界；构造只复制初始记录，更新只能由拥有者完成。隐式引用/指针跨更新保持地址，读取到本次更新值；不得把引用当作旧时刻快照。显式每次求场保留原局部快照寿命。只更新启用的 Dirichlet 波形，保留 Neumann/Robin/Periodic、禁用和缺失边界的原策略；电功选择中启用的非 Dirichlet ID 不退回默认电极。输入映射继续归 `PicBoundarySetup`，字典单位解析和求解器边界记录不变。见[场边界记录](validation/organization/QUANTUM_Modular_M6_Pic_Field_Boundaries_20260918.md)。

M6 PIC 当前电荷沉积：`ExplicitPicChargeDeposition` 绑定场拥有者及原网格，私有持有电子数、工作区和最近电荷统计；网格/场寿命覆盖它，配置、物种、CSR、归一化及几何只在调用内借用。显式采用原有效体积、alpha 优先/中心判定回退和固体清零；隐式初始化/重建经 `PicFieldState::depositFullCellCharge` 保持全单元体积及叶单元判断。复用 `DepositQK`，保留粒子权重回退、累加次序和单个缓存粒子数不符时的原跳过行为。缺少物种缓存或字段大小未随网格重置先于写入拒绝；大小检查不代表拓扑身份验证。RHS 复制、电子数字段重置仍在原节点，既有隐式轨道平均沉积归属不变。见[沉积记录](validation/organization/QUANTUM_Modular_M6_Pic_Charge_Deposition_20260918.md)。

M6 PIC 普通求场：`PicPoissonFieldSolve` 不可复制或移动，私有持有普通 Poisson 求解器、选项快照和显式 setup 状态；绑定场及网格必须覆盖其寿命，面图/切割指标/边界只在调用内借用。显式 setup、缓存求解、失败回退和报告合并保留；几何 `clearCache` 不重置 setup 状态，AMR `resetAfterMeshChange` 同时重置，两者不可合并。电场重建的配置策略与隐式普通重建默认策略保持区别；写下一步历史前核对场拥有者身份，半步形成及接受仍单独执行。隐式极化求解及周期对称处理现由 `ImplicitPicFieldSolve` 写入，完整初始化/步内编排仍需收紧。见[求场记录](validation/organization/QUANTUM_Modular_M6_Pic_Poisson_Stage_20260918.md)。

M6 PIC 隐式求场：`ImplicitPicFieldSolve` 绑定当前场，私有持有极化求解器、选项快照和 auto 模式迭代抑制状态。场/网格寿命覆盖服务，借入的极化沉积结果、面图、切割指标和边界必须对应当前网格；不复制数值状态。原模式/环境覆盖规则在构造时求值，定期几何运动禁用矩阵冻结的规则保留。求解失败保持原异常；auto 模式在迭代回退后抑制下一次迭代，显式 cg/ic 模式不锁存抑制。几何 `clearDirectCache` 只使求解缓存失效，AMR `resetAfterMeshChange` 还解除抑制，两者不可合并。求解计时仅覆盖原 solver 调用；报告、周期洁净单元配对投影、电场重建、半步形成按原顺序执行。移除两条驱动对私有场的访问，不扩称任意不匹配网格输入的校验或失败事务回滚。见[隐式求场记录](validation/organization/QUANTUM_Modular_M6_Pic_Implicit_Field_Solve_20260918.md)。

M6 PIC 二次电子发射：`PicSecondaryEmission` 私有持有配置、换算尺度、边界选择、恢复后目标物种索引/默认权重/质量、正态缓存、每步撞壁/显式发射报告和跨步调试预算；初始化仅一次，不留配置容器借用。撞壁阶段同步借用当前网格/几何、物种、RNG 和待提交缓冲，不自行稀疏化或提交。显式先统计选中离子撞壁，数量抽样后再检查发射质量；隐式沿用已有权重/能量，在启用/边界/物种/质量门通过后才计撞壁。原产额表尾端差异、整数期望数抽样、正态缓存寿命和失败落点的随机消耗保持。beginStep 清报告但不清跨步调试预算，隐式不重置正态缓存。显式域边来源由发射服务逐粒子持有：0–3 对应四边，内嵌发射为 -1；`thinExplicit` 通过原 `PicSpawnManagement` 稳定筛选粒子及来源，外部不得改排/裁剪缓冲。`accountRetained` 重算保留权重、能量矩及域边电荷，空缓冲归零而保留撞壁统计；域边电功读取这一只读报告。来源仅为步内状态，不改检查点。隐式稀疏化后的逐物种源收支仍在原调用阶段。见[保留发射修复](validation/organization/QUANTUM_Modular_M6_Pic_Retained_Emission_20260918.md)。见[发射阶段记录](validation/organization/QUANTUM_Modular_M6_Pic_Emission_Stage_20260918.md)。

M6 PIC 边界粒子源：`PicParticleSources` 持有准备好的线段、法向、源配置、原域界快照及显式待提交粒子的边界 ID；不保留 reader/原始配置引用，不拥有 RNG。空缓冲传入 emit 开始新批次，非空缓冲可追加但必须保持来源对应；通过 `thinExplicit` 使用原 `PicSpawnManagement` 稳定筛选粒子与标签。accountRetained 借入当前物种/权重/名称索引/速度尺度和域边 ID，返回只读使用的能量与四域边电荷报告；未知物种跳过，未映射/内嵌来源只计能量，重复域边 ID 按左/右/下/上首次匹配。拒绝落点或全部丢弃贡献零，调用者不得自行改排/裁剪缓冲。原几何准备、显隐法向与放置差异、数量舍入/截断、每次发射的高斯缓存和随机顺序保持；隐式保留原收支阶段，不跟踪这些显式来源标签。筛选、提交、固体剔除仍在原调用位置，来源标签不入检查点。见[粒子源收支记录](validation/organization/QUANTUM_Modular_M6_Pic_Source_Accounting_20260918.md)。

M6 PIC 新增粒子管理：`PicSpawnManagement` 拥有模式和配置快照；借用的物种名映射与默认权重容器必须覆盖其寿命，实时物种/网格/RNG/缓冲仅在调用内借用。保留 soft/hard cap 与最小概率的原规则、每候选一次均匀抽样、稳定顺序及 1/p 权重补偿，不统一显式一次全局概率与隐式分阶段概率。显式提交采用原目标计数/预留/追加并标记触及物种，不补权重或修宿主；隐式修复宿主和缺失权重、跳过域外候选、非空提交更新所有物种数量，仅 MCC 保留缺物种警告与全种群标签检查。提交不消费诊断缓冲，不承担完整性回滚或后续剔除。显式报表继续写共同概率，隐式按两次 MCC 输入数量加权；非有限旧权重的统计运算次序保持。见[新增粒子管理记录](validation/organization/QUANTUM_Modular_M6_Pic_Spawn_Management_20260918.md)。 `acceptImplicit` 明确区分 Collision/Secondary/Source：碰撞统计累加，其余覆盖；碰撞/SEE 先记账再提交，源先提交再记账，SEE 提交后才并入 sourceC。空 SEE 整段跳过；所有空保留缓冲的 commit 都不修正元数据。异常保留已完成的统计、粒子与收支，不能为了统一流程挪动这些顺序。 `acceptExplicit` 只在全部物种推送和源发射后调用：先根据全部候选计算共同概率，按 MCC/SEE/源顺序稀疏化及维护来源标签，核算 SEE/源后按 MCC/源/SEE 顺序提交。结果拥有注入能量、源四侧电荷及触及掩膜，借入引用不跨调用；调用方在原调试消息后筛除触及物种。碰撞阶段只桥接私有缓冲，不提前累计或清空。两种接受入口不能互换，见[显式新增粒子阶段记录](validation/organization/QUANTUM_Modular_M6_Pic_Explicit_Spawn_Stage_20260919.md)。

M6 PIC 合并/拆分：`ParticleManager::Step` 借用当次配置与网格，拥有保护掩膜，不跨网格/几何更新复用。显式仅从切割单元标记直接邻居；隐式从有限的 alpha < 1−1e−12 单元出发，保留面遍历顺序和原地传播；掩膜同时抑制合并与拆分。共用统计保留默认权重回退、求和/排序顺序和显隐诊断门控。`PicStepSeries` 持有管理流、21 列格式和恢复步包含式裁剪；调用者仍决定写出/刷新时机及新增粒子保留概率。见[合并拆分记录](validation/organization/QUANTUM_Modular_M6_Pic_Resampling_20260918.md)。

M6 PIC 全物种管理：`ParticleManager::Step::applyAll` 顺序消费调用者的实时物种，归约计数并在隐式模式累加带符号损失，再即时通知只读观察者；不缓存物种/报告引用，不在此做整步调度。空物种和缺省权重回退保留；短损失数组独立做范围判断，显式不写这项差额。观察者异常后不回滚当前及此前粒子变更，也不继续后续物种。驱动保留执行/报告频率、计时和显隐不同的固体筛除门。见[全物种管理记录](validation/organization/QUANTUM_Modular_M6_Pic_Population_Stage_20260918.md)。

M6 PIC 边界粒子：`PicSecondaryEmission` 从当前几何/字典解析法向，沿用选中水平集、边界覆盖和切割指标回退；发射先抽能量/方向，再查宿主，失败仍消耗原随机序列。落点使用调用者提供的域框、边距、偏移与尺度，不捕获可失效的网格视图。`ParticleManager` 的固体筛除保留顺序、可选物种掩膜和逐删除权重累加；显式不刷新宿主，隐式刷新但不顺带丢弃域外幸存粒子。截止统计保持筛除前后求和再相减，不改为逐删除累加；各项均同步借用输入，收支/提交和日志留在调用阶段。见[边界粒子记录](validation/organization/QUANTUM_Modular_M6_Pic_Particle_Boundaries_20260918.md)。

M6 PIC 物种状态：`PicSpeciesState` 不可复制/移动，拥有实时物种、主 RNG、配置默认权重/分类和恢复后名称索引；数值/IO 接口仍同步借用可变物种及引擎。配置分类保留正、负及近中性三条独立规则，恢复不重贴配置来源。名称索引保留重复名称首次插入，在恢复后建立；若外部更改运行物种名称/次序，需要显式重建且重新检查下游借用，不宣称可变借用能阻止身份更改。初始化仅一次，观察回调位于原采样后、追加前；失败停止会话，不承诺回滚。见[物种状态记录](validation/organization/QUANTUM_Modular_M6_Pic_Species_State_20260918.md)。

M6 PIC 性能输出：计时采样留在原积分阶段，`PicStepDiagnostics` 拥有计时记录及回退/排名汇总，`PicStepSeries` 拥有性能与回退 CSV；输出可用性不得改变物理推进。显式性能时间标为 `time_code`，隐式为 `time_s`；旧显式表头/短初始行仅按已知布局迁移，检查整份表的结构后再改写/裁剪。原显隐频率和刷新策略保持，见[性能输出记录](validation/organization/QUANTUM_Modular_M6_Pic_Performance_Series_20260918.md)。

M6 隐式 PIC 控制台：现有 `IO/PicStepSeries` 中的 `PicImplicitConsole` 持有本次运行的输出配置、上次采样诊断、告警防抖/冷却与 20 步耗时窗口。驱动按原阶段观察诊断并明确传入已完成步号；报告没有推进时钟，恢复后历史从空开始。当前电势/三维电子能量/碰撞分类由 `Diagnostics/PicStepDiagnostics` 的只读归约返回自有值；完整步输入只同步借用。保留原诊断快照与当前粒子数/状态的采样差别，以及隐式域电荷独立归约结果；不能借用未按隐式规则累计的显式逐事件电荷表。MCC 控制台保留原通道名字分类、未知投射物归中性组的兼容口径，不影响模型派发。控制台频率不得门控物理求解或 CSV，详见[控制台报告记录](validation/organization/QUANTUM_Modular_M6_Pic_Implicit_Console_20260918.md)。

M6 显式 PIC 控制台：`PicExplicitConsole` 在同一 IO 中持有配置并同步借用当前物理输入和原共享诊断，报告不重新采样、保存物理引用或改变步进。Diagnostics 的 `PicConsoleValues` 共用显隐电势/三维能量和旧通道归类；显式投射物按配置名字集合，隐式按当前首个同名物种电荷，不能为共用而抹去恢复后的分类差别。显式粒子数量、沉积电荷和四侧电流由专属只读归约提供，未知/周期损失、单位和面积门保持。当前状态与共享诊断允许原来的采样时层差别，非频率末步不强制打印，控制台开关不能影响 CSV 或物理。范围见[显式控制台记录](validation/organization/QUANTUM_Modular_M6_Pic_Explicit_Console_20260918.md)。

M6 PIC 输出元数据：`PicOutputCadence` 从准备配置捕获文件输出频率和结束步号，显隐驱动及 `PicOutputManifest::write` 消费同一只读对象，不再另传可互相矛盾的清单频率数组。场输出、带末步的历程、纯周期诊断、数量历程和显式详细表刷新分别具名；显式末尾补场仍在原结束阶段，隐式场仍不强制补末帧，原非正频率兼容行为保持。它不管理物理步骤、专用控制台/粒子管理频率或逐流打开阶段。`PicOutputManifest` 同步借用输入和该对象，在原启动点生成运行清单/列说明；文件角色、频率键和列声明归 IO，不在两个积分器中各保留一份序列化代码。它不持有物理状态、不更改输出门控，清单声明与实际启用输出保持区别；现有格式/顺序由对照约束。见[输出清单记录](validation/organization/QUANTUM_Modular_M6_Pic_Output_Manifest_20260918.md)。

M6 PIC 输出目录：同一 IO 文件中的 `PicOutputDirectory` 持有根目录及派生路径，构造不进行文件操作。驱动在原阶段分别请求创建、字段历史准备和旧报告清理；轨迹重置、恢复提示及报告打开不能为合并接口而换序。新运行清空整个 fields 子树，正恢复阈值只递归裁剪尾部步号不小于阈值的文件；显式 nextStep 与隐式 nextStep + 1 由驱动原样提供。清理只覆盖对应模式的旧名称，非递归 remove 保留非空目录；创建异常仍传播，不承诺回滚。见[输出目录记录](validation/organization/QUANTUM_Modular_M6_Pic_Output_Directory_20260918.md)。

M6 PIC 物理诊断：`ElectrodeCurrentGroups` 只拥有分组和每步统计，`Plasma/IO/PicPhysicsDiagnostics` 持有电极/反应/能谱输出流及配置转换。流必须实际打开且状态正常才允许采样/阻止延迟打开；嵌入电功需要的边界 ID 由物理启用条件保证，不能借用诊断或控制台频率。原采样阶段保持，见[物理诊断记录](validation/organization/QUANTUM_Modular_M6_Pic_Physics_Diagnostics_20260918.md)。 后续 `PicDiagnosticsOutput` 统一拥有分项配置快照、频率、实际流和每步门控，同时持有每步都更新的独立控制台桶；不保存反应数据库/碰撞算子指针或粒子引用。显式在 MCC 初始化后尝试并在原步内阶段重试反应打开，隐式仅在原启动阶段尝试，包括无 MCC 的空反应输出；不得相互套用。调用者保持原采样/发布阶段，最终步不强制命中诊断频率，见[诊断生命周期记录](validation/organization/QUANTUM_Modular_M6_Pic_Physics_Diagnostics_Stage_20260919.md)。

M6 通用块读取：`Core/IO/DictionaryBlockText` 持有文件读取与词法层级发现，调用者指定 scope；Solvers 的 KVBlockParser 仅保留键值解释。碰撞 loader 与 RunCase 共享来源只匹配顶层，诊断/扩展中的同名块不得提升为物理输入。别名及文件优先级、缓存名和物理解析归调用者；UTF-8 BOM 不参与首个块名。见[块作用域记录](validation/organization/QUANTUM_Modular_M6_Dictionary_Block_Scope_20260918.md)。

M6 隐式 PIC 诊断序列：physics/amount/amount_species/eb_band_rms 四类流归现有 `IO/PicStepSeries`；阶段粒子数量与嵌入带电场归约归 `Diagnostics/PicStepDiagnostics`，只借用当前场/面图并返回自有标量。驱动保留原采样、就绪联动、频率/末步及刷新顺序，不重新序列化表头和行。嵌入带的缺失/非有限电场仍计入原所选单元数分母；数量直接累计原阶段场，不补乘体积，不借结构迁移改变统计定义。真实空带、复合几何、两条 AMR 与旧检查点范围见[诊断序列记录](validation/organization/QUANTUM_Modular_M6_Pic_Implicit_Diagnostic_Series_20260918.md)。

M6 隐式 PIC 场子步：生产驱动通过现有 `ImplicitPicFieldSolve::advance` 调用轨道沉积至 next/half 电场的完整场阶段，不恢复逐项数值编排。普通场求值器/历史必须绑定同一场所有者，极化状态必须绑定同一网格；当前尺寸核对在沉积前完成。粒子与缓存只同步借用，域四侧中性判定、阶段收支、周期投影及计时顺序保持。该调用更新电势和派生层，不接受 current 电场、不推进主时钟；求解失败保留已完成沉积/诊断，不承诺事务回滚。独立对照、真实 AMR、回退及旧检查点范围见[场子步记录](validation/organization/QUANTUM_Modular_M6_Pic_Implicit_Field_Stage_20260918.md)。

M6 PIC 粒子推送：生产驱动通过 `PicSpeciesState` 的共同匹配插值/推送入口消费原场与缓存，推进策略及冲击缓冲归该状态对象并只发布 const 结果。推进策略一次配置、自有值与水平集选择；域只在配置时借用，显隐嵌入资格/选择/符号规则分别保持。初始半步、永久推进和隐式轨道预测消费同一快照，不能在每步重新传入或替换；未配置拒绝发生在清空冲击或关闭半步资格前。显式保留逐物种推送与壁面/SEE/MCC 交错，隐式保留先全部推送后处理；不能用统一批量调用改变显式时序。新步清空含空物种的旧冲击，结果不跨后续清空/同物种推送保存；失败不承诺回滚。`EmbeddedImpact` 四参数构造必须明确绑定 crossing，不再依赖可变成员顺序；原生产穿越门保持关闭。实际吸收、并行、AMR、续算和旧标记错误范围见[推送记录](validation/organization/QUANTUM_Modular_M6_Pic_Particle_Push_20260918.md)，策略所有权与新旧对照见[推进策略记录](validation/organization/QUANTUM_Modular_M6_Pic_Push_Policy_20260918.md)。

M6 显式 PIC 场阶段：初始化与普通步通过 `PicPoissonFieldSolve::advanceExplicit` 执行当前位置权重、切割体积沉积/RHS、SI 时刻边界、普通求势及当前 E。场/沉积必须同一所有者，缓存/显式耦合策略必须同一网格；所有场及电子数量存储在 AMR 后先由原阶段重置，调用前核对。面图/几何/宿主的拓扑时层仍由调用方保证，不扩称任意输入身份验证。波形/电功采样、AMR 指标、诊断和粒子半步保持原后续顺序；evalE 计时为求值与原采样两段之和。无额外场副本、数值回调或主时钟；失败保留已完成阶段，不承诺回滚。范围见[显式场阶段记录](validation/organization/QUANTUM_Modular_M6_Pic_Explicit_Field_Stage_20260918.md)。

M6 PIC 粒子初始化/重定位：`PicSpeciesState` 拥有初始匹配插值、反向半步及一次性资格；检查点准备、已执行半步或已开始永久推送的状态禁止再初始化。缓存/粒子数与场尺寸预检在修改前，开始数值操作后的失败终止会话，不承诺回滚。原半步吸收/压缩与不漂移规则保持，不加入完整推送的周期回绕。AMR 宿主/背景宿主/物种数量更新也归粒子状态，显式配置开关与隐式固定查找顺序不能混同；关闭提示时只对负标签执行回退的旧规则保持。几何/缓存时层及后续权重重建仍由调用方负责。并行、真实细化/粗化与续算范围见[粒子初始化记录](validation/organization/QUANTUM_Modular_M6_Pic_Particle_Initialization_20260918.md)。

M6 隐式 PIC 普通场重建/恢复：`PicPoissonFieldSolve` 统一启动、回退与 AMR 的当前位置权重、完整体积沉积/RHS 和普通场重建；可选历史只在求 E 成功后播种，启动不提前创建/播种。映射快路径单独求 E 并播种，不沉积或求势；有效 IPF1 恢复原样保留七个字段及不同时间层，不重建权重。恢复只先检查历史归属，备用重建的所有者/策略/尺寸在实际回退时检查，允许有效 reader 重建旧尺寸存储。原自动缺必需文件拒绝、显式步号读失败回退及拒绝前碰撞输入缓存保持；格式/选择完整性不在本次接线中另改。拓扑时层由调用方保证，失败不承诺回滚。范围见[隐式重建记录](validation/organization/QUANTUM_Modular_M6_Pic_Implicit_Field_Rebuild_20260918.md)。

M6 PIC 诊断快照：`PicDiagnostics::sampleExplicit/sampleImplicit` 共用有效体积、物种人口/密度及三维电子归约，显式保留电势极值，隐式保留变系数高斯输入和原四子槽体积判定。名字优先、权重回退、alpha 裁剪/缺项及串行算式保持；不以共享代码抹去已有判定差异。`summarizeImplicitPicParticles` 在原位置执行电荷→可选能谱观察→动能累计；观察者只同步消费准备值，不修改借用粒子或保留引用。收支完成仍在粒子汇总后、场诊断前；IO 只序列化，不扫描物理状态。驱动保留原采样阶段、频率、上次诊断和后续 CSV 调用。显式范围见[诊断采样记录](validation/organization/QUANTUM_Modular_M6_Pic_Explicit_Diagnostic_Sample_20260918.md)，共享归约及隐式观察见[隐式诊断记录](validation/organization/QUANTUM_Modular_M6_Pic_Implicit_Diagnostics_20260918.md)。

M6 显式 AMR 状态协议归 `ExplicitPicRefinement`：网格/几何转移、壁面通知、派生视图、场重置、宿主重定位及缓存/耦合更新作为同一模块操作。到期后先核对私有网格/场归属和模式；friend 仅服务该核对，不恢复驱动对私有状态的访问。服务持有报告，调用方保留整步阶段/控制台输出；失败不承诺回滚。原首子节点叶判定保持。范围见[显式 AMR 记录](validation/organization/QUANTUM_Modular_M6_Pic_Explicit_Amr_Stage_20260918.md)。

M6 显式 AMR 配置：已有输入 IO 的 `ExplicitPicAmrSetup` 私有持有 AMR 和细化器输入，只读视图随该对象存活，不保存网格/粒子引用。构造只解释准备输入；网格初始化后执行层级启用判断，物种恢复和源准备后再执行电子/指标判断，各阶段只允许一次。报告打开、后续关闭、调试消息和温度求值保留原阶段，不能在配置读取时提前关闭报告或选择未恢复物种。保留配置块存在/启用差别、无配置块的细化器回退、原指标名单和温度的 NaN/Inf 行为；实际指标及拓扑操作仍归物理服务。范围见[显式 AMR 配置记录](validation/organization/QUANTUM_Modular_M6_Pic_Explicit_Amr_Setup_20260918.md)。

M6 隐式 AMR 配置：已有输入 IO 的 `resolveImplicitPicAmr` 只解释已准备的 `PicCaseInput`，产生自有指标/几何选项及报告来源；不得在此读取当前环境、字典或持有运行态引用。实际几何细化门和电子编号由拥有者提供，报告仍在驱动原阶段写出。保持 mesh/PIC 优先级、类型化覆盖、不同默认/下限与原警告规则；冲突去空白和运行模式不去空白不能混同。旧 exact-refresh 选择仅影响兼容报告，物理精确刷新仍依据当前移动几何状态，不把迁移变成新物理开关。见[AMR 配置记录](validation/organization/QUANTUM_Modular_M6_Pic_Amr_Setup_20260918.md)。

M6 隐式 AMR 状态协议归 `ImplicitPicRefinement`：重建/映射分别保留，输出坐标通过只读网格的同步输出回调在场重建后刷新，再完成 pending 和原可选瞬态清理。到期操作先核对活状态/视图归属；失败不回滚或发布完成报告。原 AMR 诊断之后的求解器缓存/重试重置仍由外层在原位置调用。不得把输出回调用于额外物理推进。范围见[隐式 AMR 记录](validation/organization/QUANTUM_Modular_M6_Pic_Implicit_Amr_Stage_20260918.md)。

M6 PIC AMR 报告归现有 `PicStepSeries`：驱动只选择原开关、频率和采样阶段，IO 持有文件流、表头、重启裁剪及序列化。显式 defaultfloat 与隐式 fixed16 格式分别保留；level-set 索引仍从 1 开始，掩码只读取当前几何/报告，不保留跨 AMR 的引用。控制台仅格式化已准备数据，不扫描网格或推进物理。指标关闭、失败流、原刷新与求解器失效时机均保持。范围见[AMR 报告记录](validation/organization/QUANTUM_Modular_M6_Pic_Amr_Reporting_20260918.md)。

M6 PIC 运行域启动归 `PicDomain`：先 `prepareMesh`，在原位置映射场边界，再 `initializeForRun`；拥有者统一几何恢复、细化、静态规则卸载及派生面图/指标准备，返回已解析的运行开关。网格准备和几何初始化各执行一次，开始后的失败停止会话，不回滚或重试。显式缺网格可继续、恢复后可启动细化；隐式缺网格失败，恢复/已加载快照跳过启动细化，检查点优先于快照。原快照读写异常处理与显隐 phase/composite 顺序保留；MovingEB 仍由其拥有者在原后续位置初始化。refiner 借用的几何输入须保持有效。范围见[运行域启动记录](validation/organization/QUANTUM_Modular_M6_Pic_Domain_Startup_20260918.md)。

M6 PIC 物种启动归 `PicSpeciesState::initializeForRun`：采样或读取 RNG/粒子后才执行嵌入剔除，成功返回首个循环步；域按次借用，模式只读。部分恢复失败保留已读入状态，不提前剔除、回滚或重试。运行中的剔除共用同一入口，显式保留存活宿主、隐式刷新宿主；掩码和损失数组按原规则处理。名称索引、轨迹编号、半步修正、场和主时钟仍在原阶段，禁止因组合入口提前执行。范围见[物种启动记录](validation/organization/QUANTUM_Modular_M6_Pic_Species_Startup_20260918.md)。

M6 隐式切割单元诊断：`PicCutcellStateView` 在 Diagnostics 中按当前网格/场同步采样，IO 只选择行并序列化；原法向场内核、单位来源和边界 map/fallback 优先级保持。视图与行中的边界指针不跨输入有效期、AMR 或状态更新保留，不复制整套物理场。驱动保留原场子步后的门控和 current/next 字段组合。目录创建、法向场计算、打开/截断、逐行输出的顺序及原异常/警告规则不能因迁移改变。范围见[切割诊断记录](validation/organization/QUANTUM_Modular_M6_Pic_Cutcell_Dump_20260918.md)。

M6 碰撞初始化/恢复（2026-09-19）：`PicCollisions::create` 可同步借用已选定的 `PicCheckpointReader`，先完成原物种/截面/算子初始化，再按 `saveRng` 恢复计数，最后返回拥有者；不保存 reader 引用。创建失败或禁用时不读检查点，恢复失败沿用 IO 异常且不返回部分对象。`PicCollisions` 不再公开 RNG 重置，reader 的对应读取只供工厂调用；只读 RNG 查询、原显式整步/隐式两半步、报告与输出预检时机保持。低层算子、检查点格式和旧 global/perParticle 缺计数组件规则不改；counterBased 仍要求该组件。范围见[初始化记录](validation/organization/QUANTUM_Modular_M6_Collision_Initialization_20260919.md)。

