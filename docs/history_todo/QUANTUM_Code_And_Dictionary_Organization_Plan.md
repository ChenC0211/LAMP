> **状态：historical/reference（历史资料，仅供追溯）。**
> 本文原始正文保留为历史记录，其中“当前”“下一步”、TODO、旧规则与通过记录均指写作当时；不再自动产生新任务，也不能证明 V2 gate 通过。
> 现行任务唯一入口：[分模块任务清单](../roadmap/05_LAMP_Module_Task_Roadmap.md)。文档权威见[开发总纲](../00_LAMP_Development_Charter.md)，完整入口见[文档索引](../README.md)。
> 旧字典/职责说明可供查询既有实现；V2 新增或修改仍遵守现行六份规范。历史正文中的链接保留原样，目标缺失时不能作为本机已验证证据。
> 本目录中的旧设计、TeX/PDF 和索引同属历史参考；总纲指定的三份冻结方法 PDF 仍以 `docs/` 下的权威文件为准，不因本目录的历史标记而降级。

# QUANTUM 代码结构与字典规范整理方案

2026-09-28 用户补充：界面与流体集中调试完成后，再检查全项目模块化、字典规范、代码卫生与多余支路，并整理控制台及结果输出。当前唯一排期见[主计划阶段 P](QUANTUM_Modular_Consolidation_Plan.md#project-usability-cleanup)，本页不另开一套执行清单。

后续执行衔接（2026-09-15）：用户已要求启动全项目模块化收紧，新的分批安排见[模块化收紧计划](QUANTUM_Modular_Consolidation_Plan.md)。本页保留 A1/A2、B1/B2 等已实施整理的历史与证据；不将下文旧排期和新计划并行维护成两套当前任务。已生效规则仍以正式规范为准。

日期：2026-09-11；实施更新：2026-09-12。A1 阶段状态：**原喷嘴导电短程已接通，App 目录、EHD 适配和界面力 IO 拆分已实施；三组原输入对照全部通过，EHD 28 个、导电喷嘴 76 个、旧介电喷嘴 86 个完整数值文件分别与整理前逐字节一致。** 详见 [A1 记录](validation/organization/QUANTUM_App_Organization_A1_20260912.md)。最新输入能力见下方 B2 更新；正文中的拟议设计不替代正式参考。

**A2 更新：** 剩余流体区域/材料绑定、初始/运动界面诊断及时间循环 CSV 历史已拆分，构建和结构检查通过；两条原喷嘴均完成 20 步、退出 0，导电 76 个、旧介电 86 个完整数值文件分别与 A1 逐字节一致，见 [A2 记录](validation/organization/QUANTUM_App_Organization_A2_20260912.md)。

**B1 更新：** 线性系统与守恒门限的独立接线、生效数值控制报告已完成；导电 76 文件、介电 86 文件、EHD 28 文件的原输入 20 步对照全部一致，原 RZ 曲率失败保留。见 [B1 记录](validation/organization/QUANTUM_Dictionary_B1_Numerical_Controls_20260912.md)。

**B2 更新：** 两条运动喷嘴已接入 schema 2：独立材料表、模块物理字典、五个实际线性系统、Time/Output/Log、只读几何绑定和生效配置记录。新旧生效输入逐项一致，实际 CLI/数值接线检查通过；原 20 步复测及 Tecplot 结果状态见 [B2 记录](validation/organization/QUANTUM_Dictionary_B2_Versioned_Input_20260912.md)。实际语法以[字典参考第 5 节](QUANTUM_Dictionary_Reference.md#5-schema-2材料模块和运行控制分开)为准。该支持不推广到其他 App，也不等于 D 批或全部物理验收完成。

本计划保留整体目标与尚未实施的设计；已生效部分以[正式代码组织规范](QUANTUM_Code_Organization_Standard.md)为准，并由仓库 `AGENTS.md` 引用及 CTest 检查约束。实际迁移和验收状态见规范第 5 节；不把 A1/A2 称为 A/B/C/D 全部完成。数值主线见[当前进度](QUANTUM_Hydro_EHD_Progress.md)，现有用法仍以[算例指南](QUANTUM_Hydro_EHD_Case_Guide.md)为准。

## 1. 建议结论

保留 `Core / Solvers / Modules / Apps / cases` 的基本分层，不另造一套框架。整理重点不是把文件夹改好看，而是让四个问题有唯一答案：

- **运行哪个过程？** 看 `Apps/<应用名>/` 和 `run.dict`。
- **方程、离散和边界条件在哪里？** 看相应物理模块；跨模块的调用顺序属于 App。
- **某个物理量在哪里配置？** 相物性统一一处，流体和电学边界各归自己的字典。
- **一次时间步怎样推进、何时接受结果？** 看一段清楚的驱动流程，不必追踪散落的初始化函数和输出副作用。

用户确认的顺序保持：**先接通原喷嘴算例，再分批整理并固化规范。** 导电 20 步已经完整运行并保存收支/场基线，结构整理现已开始；它不替代长时喷射物理验收。先拆职责、再迁配置，结构批次不改数值算法或算例参数。

## 2. 当前问题与代码依据

本节保留提出方案时的问题清单，部分已由 A1/A2 解决；当前路径、职责及尚未完成项以正式规范和各批实际记录为准，不把旧目录描述当作新代码模板。

| 当前情况 | 实际影响 | 处理方向 |
|---|---|---|
| `Apps/Runtime/Apps` 同时放 App 入口、喷嘴界面运动、毛细力包装、电场组装和开发模式 | 文件名都像 App，难以判断哪些能共用、哪些只是接线 | 应用按目录分组；计算实现和应用适配分别归位 |
| `RunController.cpp` 约 1061 行，识别 Hydro/EHD、解析 PIC/MCC 相关共享数据、处理单位和报告 | 加 App 需要继续修改总控，普通流体也被迫接触等离子体配置逻辑 | 总控只负责选择与调度；物理配置由对应装载器负责 |
| `ElectrohydrodynamicsApp.cpp` 内部有较完整的 `DielectricCoupling`；喷嘴在另一套驱动中组织流动、电场、界面和输出 | 两条路径的生命周期、单位与场量交换关系不直观 | 先显式列出差异，再提取相同的耦合适配；不能直接认定两套离散等价 |
| 喷嘴 `electricity.dict/ElectrosprayPhysical` 同时包含密度、黏度、供液压力、固壁、接触线、电势和残差容限；该短程没有 `hydro.dict` | 修改流体条件要去找电学字典，配置和模块名称不对应 | 物性、流体、电学、线性求解设置分别归属 |
| 同一个 `electricity.dict` 在旧独立 Electricity、通用 EHD、喷嘴入口中有三种读取方式 | 不能从文件名判断支持哪些键；电学物理解析还混在求解实现中 | 新版本使用同一电学配置定义，旧格式分别适配；解析实现归 Electricity/IO |
| `initializeElectrosprayCapillarity`、`initializeElectrosprayEhdForces` 每步也被调用，内部同时计算和写文件 | `initialize` 名称误导；计算函数受输出开关影响，候选状态也可能先写文件 | 拆成无 IO 的计算函数与明确的输出函数 |
| Hydro 的时间/输出在 `hydro.dict`，喷嘴的时间/输出在 `run.dict`；残差容限又分散在不同配置 | 换一个 App 就要重新学习“在哪里改步长、输出、求解器” | 新规范统一运行控制；旧写法通过兼容层映射 |
| Hydro/EHD 当前拒绝共享 `units/normalization/solvers.dict`，喷嘴则走另一条装载路径 | 字典同名不等于模块已经消费设置；当前不能仅移动参数就宣布通用 | 迁移字典必须同步接通实际消费者，并输出生效配置 |
| 喷嘴加载了 `solvers.dict`，但实际电场与流体仍直接指定 SuperLU，读取 `ElectrosprayPhysical.linearResidualTolerance` | 用户修改 solver 文件可能不改变实际求解；同一容限还被用于相通量收支检查 | 为每项参数标明消费者；分别配置线性求解和守恒检查，迁移时保留原值 |
| `DictionaryElectrospray` 与严格 Hydro/物性解析器行为不同 | 未知键、重复键、数值尾随字符、无效布尔值等检查不统一 | 公共语法读取 + 各自类型校验，错误带文件、行号、字段路径 |

主要依据：

- [当前应用装配入口](../src/Apps/RunCase.cpp)、[其头文件契约](../include/Apps/RunCase.h)、[通用注册/执行](../src/Apps/Runtime/AppRegistry.cpp)。原 RunController 已退出；原白名单、共享装载策略与时序迁入应用入口，见[分派记录](validation/organization/QUANTUM_Modular_M6_App_Dispatch_20260918.md)。
- [喷嘴入口](../src/Apps/Electrospray/ElectrosprayApp.cpp)、[喷嘴时间推进](../src/Apps/Electrospray/ElectrosprayTimeLoop.cpp)、[提取后的 EHD 适配](../src/Apps/Electrohydrodynamics/Coupling/DielectricCoupling.cpp)。
- [毛细力适配](../src/Apps/Electrospray/Coupling/ElectrosprayCapillarity.cpp)、[电力合成](../src/Apps/Electrospray/Coupling/ElectrosprayForces.cpp)、[已分离的输出](../src/Apps/Electrospray/IO/ElectrosprayForceDiagnostics.cpp)。
- [喷嘴物理配置](../cases/electrospray_moving_dielectric_short_rz/system/electricity.dict)、[旧喷嘴运行配置读取](../src/Apps/Electrospray/Compatibility/DictionaryElectrospray.cpp)。

这不是否定已有分层：`HydroFitted*`、公共 Hermite 几何、Electricity 的电荷输运以及新建的分层输出已有可复用基础，应在其上整理。

## 3. 目标代码结构

`include/` 与 `src/` 保持对应。下图表示目标职责，不要求一次创建所有目录或空壳类；未涉及的 Plasma 内部结构保持原样。

```text
include/ 与 src/
  Apps/
    quantum_main.cpp                 # 仅 src：命令行、选择 case、返回退出码
    Runtime/                         # App 接口、注册表、case 路径、通用运行控制
      IO/                            # run 配置与兼容输入装载
    SingleModule/                    # 独立运行入口；模块运行适配
    Electrohydrodynamics/
      ElectrohydrodynamicsApp.*      # 通用 EHD 工作流
      Coupling/                      # 电荷/电场/界面力与流动之间的适配
      IO/                            # 跨模块诊断和输出组合
    Electrospray/
      ElectrosprayApp.*              # 电喷雾工作流，不实现 NS 或 Poisson 装配
      Setup/                         # 将 case 的相、边界、接触线绑定到共用模块
      IO/                            # 应用配置、喷雾过程输出
      Compatibility/                 # 迁移期旧输入与旧开发模式的入口映射
    Verification/
      Pic/                           # 现有验证 App；解析解/误差判据不混入生产入口
  Modules/
    Hydrodynamics/
      Flow/                          # NS、压力恢复、流动状态和求解报告
      Boundary/                      # 开放边界、固壁、轴线等流体条件
      Interface/                     # 界面输运、毛细力、相体积与体积修正
      IO/                            # hydro 字典、流体输出字段组装
      Detail/                        # 不供外部直接使用的实现
    Electricity/
      Electrostatics/                # 电势、电场、介电界面离散
      SurfaceCharge/                 # 界面电荷、导电源项与电荷输运
      Forces/                        # 电应力/极化力计算
      IO/                            # electricity 字典和电学输出
      Detail/
    Plasma/                          # 本轮不迁移 PIC/MCC 数值实现与输入
  Core/
    Geometry/  Mesh/  Field/          # 保留现有公共基础
    Materials/                       # 拟新增：相属性表/相间属性表，不含求解流程
    IO/                              # 通用字典语法、格式写入、OutputLayout
    Units/                           # 单位和参考尺度，不含喷嘴物性常数
  Solvers/                           # 通用数值算子、线性代数及求解器配置
```

### 3.1 App 的边界

- `SingleModuleApp`：独立运行某一模块。现有 `modules (...)` 保留兼容；模块顺序运行不等于共享网格的物理耦合。
- `ElectrohydrodynamicsApp`：可复用的流体—电场—界面电荷耦合工作流。平面/RZ、介电/导电是受支持能力组合，不是为每个算例建一个 App。
- `ElectrosprayApp`：组织供液、喷嘴/电极、接触线以及未来发射或粒子交换；使用相同的流体、电学能力。喷嘴尺寸、压力、电压、物性只在 case 中。
- `Verification`：解析参考、制造解和误差判据。标准物理算例仍通过正常 Hydro/EHD/App 运行，不复制一套求解器。

各 App 不调用另一个 App 的 `run()` 来实现耦合，而是共用模块和必要的耦合适配。暂时不强行合并矩形域 Hydro 与实际喷嘴活动流体域的离散：网格与自由度、力的位置、边界处理均须保持原含义。

### 3.2 依赖与所有权

1. `Core` 不依赖物理模块或 Apps；`Solvers` 不读取 case 文件、不认识喷嘴和电喷雾。
2. 物理数值对象接收类型化参数、网格/几何视图和场，不读取 `run.dict`，不依赖 `ElectrosprayConfig`。
3. 模块之间通过场和有单位约定的数据交换；跨 Hydro/Electricity 的具体接线在 App 层。纯电学算法不为使用方便而整体放进 Hydro。
4. 现有 `HydrodynamicsModule` 等 `IModule` 入口先保留，逐步拆为薄运行适配和可复用的单步求解能力；不要求本轮同时改动全部 PIC 入口。
5. 每个运行过程只有一个主时钟和一套权威网格/几何。界面快照、面元和通量数据说明所属时间层与有效期；不能让一个模块私自重建另一套界面或推进第二个时钟。
6. App 只把检查通过的候选结果标成“接受步”。必要的失败诊断可以输出，但不能进入接受帧清单；IO 失败与物理求解失败分别报告。

现有 `RunConfig` 同时装运行选项和已加载的气体/反应数据库。逐步区分“运行配置”和“已加载上下文”，保留兼容接口；不要把所有模块状态都搬进一个新的全局 Context。单位尺度也只生成一份明确的值对象，不能在毛细力、电力、输出中分别推导出不同版本。

### 3.3 具体文件怎么拆

| 现有文件/内容 | 建议归属与处理 |
|---|---|
| `Runtime/Apps/*App.*` | 分别放入相应应用目录；入口只保留配置、初始化、推进和错误处理的主干 |
| `ElectrohydrodynamicsApp.cpp` 内 `DielectricCoupling` | 提取到 `Apps/Electrohydrodynamics/Coupling/`，保持调用次序和时间层 |
| `ElectrosprayHydrodynamics.cpp` | case 边界绑定留 `Electrospray/Setup`；公共流动求解仍在 Hydro；时间循环留 App；CSV 字段组装移到 IO |
| `ElectrosprayInterface.cpp` / `ElectrosprayInterfaceMotion.cpp` | case 边界/唇缘绑定留 Setup；与喷嘴无关的输运和几何检查提取到 Hydro/Interface；纯几何算法留 Core/Geometry |
| `ElectrosprayElectrostatics.cpp` | 相/电极 ID 绑定留 Setup；可复用的材料映射和电场求解移 Electricity；写场移 IO |
| `ElectrosprayCapillarity.cpp` / `ElectrosprayForces.cpp` | 计算函数无文件副作用；毛细计算归 Hydro，电力计算归 Electricity，总力组合归耦合适配 |
| `Runtime/IO/ElectrosprayPhysicalConfig.*` | 旧格式 reader 留 Compatibility；新配置不再用一个结构包揽所有模块物理量 |
| `Runtime/IO/SharpDielectricDiagnostics.*` / `SurfaceChargeDiagnostics.*` | 仅依赖电学数据的部分归 Electricity/IO；确实跨流体、电荷的部分留 EHD/IO |
| `Runtime/IO/DictionaryGas/CrossSections/Reactions/Collisions` | 后续随 Plasma 配置职责整理；先从总控物理逻辑隔离，不顺便重写已跑通的碰撞输入 |

命名使用 `initialize / advance / assemble / evaluate / validate / write` 的实际语义；`initialize` 不能兼作每步更新。`Config` 只放输入，`State` 只放演化状态，`Report` 放诊断，`Writer` 只负责输出。不新增无明确职责的 `Utils/Common/Manager` 大杂烩，也不按行数机械拆文件。

## 4. 字典的目标职责

保持 `cases/<case>/system/*.dict` 和 `QUANTUM.exe --case <case>`。不恢复 `app.dict`，也不为每个子功能增加一个字典。

**建议仅新增一份共享 `materials.dict`。** 这是对白名单的有意扩展，需要装载器、校验器、文档一起落地；它是共用材料数据，不是新的物理模块。其余优先整理现有字典的块结构。

| 字典 | 唯一职责 | 不应再放的内容 |
|---|---|---|
| `geometry.dict` | 现有几何、相、水平集、边界实体与编号 | 压力、电势、黏度；本轮不改旧几何格式 |
| `mesh.dict` | 网格层级、初始细化、AMR 控制 | App 再定义一套网格细化参数 |
| `run.dict` | App 选择、时间、重启、耦合调度、输出、日志 | 密度、黏度、电导率、毛细系数、线性求解器参数 |
| `materials.dict`（拟新增） | 按几何相 ID 绑定材料：密度、动力黏度、相对介电常数、电导率；按两相组合绑定表面张力系数 | 边界条件、时间步、力的离散厚度 |
| `hydro.dict` | 流体模型、初始流场、速度/压力边界、界面输运、毛细模型及其离散设置、接触线条件 | 电压、电荷模型；新格式不再重复材料常数及运行输出设置 |
| `electricity.dict` | 电学模型、电势边界、电荷初值/输运模型、电界面离散 | 供液压力、无滑移、液体黏度、接触线钉扎 |
| `solvers.dict` | 按实际方程/线性系统配置后端、线性容限、迭代控制 | 物理边界、界面厚度、以线性容限替代守恒验收 |
| `units.dict` / `normalization.dict` | 保留现有单位和参考尺度入口；分清量纲换算与几何积分约定 | 用粒子数密度代替流体质量密度；把 RZ 的三维体积归一化误称为三维网格 |
| `pic.dict` / `mcc.dict` | 现有 Plasma 规范与兼容路径 | 本轮不要求迁入新 materials 或改写已运行的 Plasma 算例 |

按能力加载所需文件：单相 Hydro 不要求 `electricity.dict`；EHD 不加载反应截面或 MCC。不是每个 case 都要补齐全表。

材料表可以同时记录流体和电学属性，但模块只消费所需字段。比如理想介电计算保留材料电导率作为资料，不代表已推进导电；生效配置必须列出“本模型未使用”，不能与请求了却未实现的电荷输运混为一谈。

### 4.1 块结构与统一语法

- 新规范使用明确的配置版本，例如 `Run.schemaVersion 2`；该版本仅用于本次运行/物理配置迁移，不重解释旧 geometry/mesh/PIC 的格式版本。未写版本按现有兼容路径读取。
- `run.dict`：`Run`、`Time`、`Output`、`Log`、可选 `Restart`，以及被选 App 的同名工作流块。
- `hydro.dict`：`Hydrodynamics` 下分 `Model / Initial / BoundaryConditions / Interface / Numerics`。
- `electricity.dict`：`Electricity` 下分 `Model / Initial / BoundaryConditions / SurfaceCharge / Numerics`。
- 块名统一大驼峰，键名统一小驼峰；新文档只展示一个规范写法。历史别名在兼容 reader 处理，不扩散到求解代码。
- `model` 用经过校验的枚举；只支持明确列出的模型组合，不依靠多个松散布尔量组合出“看起来可运行”的分支。
- 边界条目按现有 geometry 边界 ID 绑定；相、水平集、接触点保留原身份和 1-based 约定，进入内部后集中转换。先用条目标签/注释提高可读性，不修改几何编号来配合目录整理。
- 物理参数新格式优先 `数值 + 单位`，如 `pressure 45.2 Pa`；旧 `_SI` 和代码单位数值通过兼容 reader 明确转换。单位支持必须实际实现量纲检查后才发布，不只修改示例。
- 标量键不允许重复；边界等集合允许多条，但同一对象的条件不能重复或冲突。未知键、未知块、缺少必填值、数字尾随垃圾、无效布尔值均报错。
- 错误信息至少包含 `文件:行号 + 块/字段路径 + 原值 + 期望`。解析语法归通用 reader，数值范围归各模块，跨字典 ID/模型组合归 App 的配置检查。

### 4.2 时间、输出、日志的统一示意

以下是**拟议语法，不是当前可直接运行的输入**。为复现原喷嘴，示意中仍保留原来的代码时间值，不改成更长的物理过程。

```text
Run
{
    schemaVersion 2
    app Electrospray
}
Time
{
    units code
    dt 5e-6
    endTime 1e-4
}
Output
{
    prefix electrospray
    Fields
    {
        enabled true
        format tecplot
        writeControl timeInterval
        writeInterval 5e-5
        writeInitial true
        writeFirst true
        writeFinal true
    }
}
Log
{
    everySteps 1
    level summary
}
Electrospray
{
    Coupling
    {
        electricUpdate everyStep
    }
}
```

`timeInterval` 的单位继承 `Time.units`；采用按步输出时写 `writeControl stepInterval`，间隔必须是整数。`dt/endTime` 是新格式唯一主时间来源；旧 `steps` 由兼容层换算并校验，不允许两处同时控制主时钟。新 `endTime` 不整除固定 `dt` 的行为必须明确，首批保持当前拒绝规则，不偷偷截短末步。

输出位置固定遵守 Plasma 的用途分类：`output/fields`、`output/series`、`output/meta`，**没有 `tecplot/`**。`format` 只决定文件格式。空间场可降频，逐步收支诊断和失败信息不能因控制台降频而消失。

控制台显示选中的 App、**实际启用**的物理模型、单位、网格/界面身份、求解器、步数与物理时间，以及每步连续性、动量、体积/电荷收支。详细装配信息只在详细日志中出现；电学与流场若不是同一时间层，必须分别注明。

### 4.3 原喷嘴参数的归属示例

| 原参数 | 新归属 | 必须保持的物理含义 |
|---|---|---|
| 液体 `1240 kg/m³`、`0.038 Pa·s`、`εr=10`、`1.4 S/m` | materials：几何相 2 | 不因改配置改变物性 |
| 气体 `1.184 kg/m³`、`1.849e-5 Pa·s`、`εr=1`、`0 S/m` | materials：几何相 3 | 回流材料仍为环境气体 |
| 表面张力 `0.0452 N/m` | materials：相 2–3 的相间属性 | Hydro 只配置启用模型与数值处理 |
| 边界 5 的 `45.2 Pa` | hydro / BoundaryConditions | 压力供液，不改成指定速度/流量 |
| 边界 8、9 的 `0 Pa` | hydro / BoundaryConditions | 两处均允许气体回流，入流材料和动量闭合显式保留 |
| 固壁无滑移、点 6 钉扎 | hydro / BoundaryConditions、Interface | 固定接触位置，不能顺便规定接触角 |
| 导体相 1 为 `300 V`，边界 9 为 `0 V`，其余电通量条件 | electricity / BoundaryConditions | 流体出口与电学接地各自独立，完整迁移所有电边界 |
| `capillaryHalfWidthCells 3` | hydro / Numerics 的界面力设置 | 数值力带半宽，不当作物理界面厚度 |
| `linearResidualTolerance` 的线性求解用途 | solvers 的对应线性系统条目 | 分清电势、动量、压力与耦合系统；原值不能用统一默认值覆盖 |
| 同一旧容限的相通量守恒检查用途 | hydro / Numerics / Conservation | 与线性容限分开命名，仍保留原值和原分母定义 |

目前 `ElectrosprayPhysical.model=leakyDielectric` 与 `physicsMode=movingDielectric` 并存，实际运动链为 `q=0` 介电。**新配置不能继续靠这种二义性表达模型。** 新理想介电算例应明确选择理想介电模型；导电液面必须明确请求电荷演化，未接通就失败。兼容模式保留旧计算含义、原始请求和“实际未推进导电”的报告，不能自动把正式 `coupled` 算例改为介电短程。

## 5. 配置读取和代码可读性的落地规则

### 5.1 只解析一次，配置来源可追踪

装载顺序：选择 case/App → 读取所需字典 → 类型/单位检查 → 相与边界绑定检查 → 生成只读生效配置 → 初始化求解对象。计算函数不重新打开字典，也不再次按字符串判断 App 模式。

公共语法 reader 放 Core/IO；`DictionaryHydrodynamics`、`DictionaryElectricity` 等仍由各自模块拥有，不新建一个全知的物理配置大类。现有 `KVBlockParser`、`HydroInput::Reader` 先保留兼容接口，再提取共同的词法和错误定位；不要在第一批替换全部 geometry/PIC 解析器。

新旧配置同时提供同一项时**报冲突而不规定静默优先级**。迁移前后比较类型化配置和换算后的数值；不能仅比较字典文字。

普通运行在 `output/meta/` 记录实际生效值、单位/参考尺度、来源字段、默认值与未使用材料属性；`--check-case` 只读检查且不写这些文件。当前预检查只支持独立 Hydro/EHD，喷嘴支持必须真正实现后再写入使用指南，不能靠放开白名单假称支持。

### 5.2 求解器设置不能只“读到了”

`solvers.dict` 目标按实际系统配置，例如电势、动量、压力与拟合域耦合系统；条目名称须与代码中的真实调用对应。Hydro/EHD 的共享 solver 配置接入前，继续如实拒绝不支持的设置。

SuperLU 作为本项目要求的后端必须真正使用；缺少依赖时明确报错，不静默回退 Eigen SparseLU。直接求解器不接受看似生效的 CG/IC 组合。几何判据、守恒误差阈值与线性残差是三类不同设置，不能因为整理成一份字典而混为一个 tolerance。特别是目前喷嘴物理配置中的同一个容限有多种用途，兼容层应分别填入原值，而不是只迁移其中的线性求解用途。

#### 5.2.1 B 批接线前核对（2026-09-12，尚未启用新输入）

A2 等待原样对照期间核对了真实调用。`ElectrosprayPhysical.linearResidualTolerance` 不能简单搬成一个“线性求解精度”后删掉其他用途：

| 当前用途 | 实际位置/约束 | 迁移要求 |
| --- | --- | --- |
| 初始、q=0 更新及 P3 电势求解 | `Setup/ElectrosprayElectrostatics`；P3 模块还据此检查离散高斯收支，并有原机器精度下限 | 分清线性系统与电通量检查，旧输入分别映射原值；不删掉派生下限 |
| 拟合域 NS、局部体积修正的线性子问题 | `ElectrosprayTimeLoop` 的 `HydroSolveOptions`、`ElectrosprayInterfaceMotion` 的 `correctionOptions.linearTolerance` | 按真实系统接入 solver 配置；体积修正的非线性停止准则不是同一个参数 |
| 相通量装配、压力口通量一致性 | `ElectrosprayTimeLoop` 的两项收支检查 | 归流体守恒设置，保持原分母及原值 |
| q=0 极化载荷的电场/通量一致性 | `buildFittedDielectricPolarization` 的 `relativeTolerance` | 归电学数值检查，不冒充重新求解 Poisson 的容限 |
| 运动电荷、旧电流及电极收支 | `ElectrosprayConducting` 使用 `16 * linearResidualTolerance` | 保留该原派生值及原绝对电荷/电流归一化，不变成另一个默认门限 |

同样，给 `DielectricConfig` 成员赋值不等于该条求解/重构分支消费了它，生效配置应逐项沿调用核实。当前 `UnitParser` 返回 SI 值和单位 token，未给出字段量纲类型；新语法要复用换算并补足字段量纲约束，不能把“能识别单位”当作“物性单位正确”。旧 `HydroInput::Reader` 是严格 token reader，旧喷嘴物性另有严格的扁平行 reader；优先提取已有词法/定位能力并保留兼容适配，不再复制第三套宽松 parser。`RunController` 的现有白名单尚无 `materials.dict`，喷嘴 `--check-case` 也尚未支持；两者必须与新装载/检查路径一起实现，不能只改白名单或示例。

### 5.3 让主流程读起来像数值流程

B1 接线实施记录见 [数值参数分离](validation/organization/QUANTUM_Dictionary_B1_Numerical_Controls_20260912.md)。进一步沿 NS 调用核对后，已将动量预测、压力修正、耦合修正的线性门限，以及完整动量平衡、局部连续性、出口可接受性分开；旧值和分母保持。新版字典尚未启用，后续配置迁移应接入这些真实独立成员，不再通过一个旧统一容限控制所有检查。

驱动主干按“初始化 → 本步计算 → 检查 → 接受状态 → 输出”组织，具体时间分裂和场量时刻保持现状。不要为视觉整齐改变电荷、界面、材料或载荷的更新顺序。

报告使用具名字段：阶段、时间、单位、自由度位置、残差、接受状态。与求解相关的函数注释优先解释离散约定、所有权和适用范围；过时的“仅平面”“SparseLU”“只负责选择 App”等注释随对应实现修正。

目录移动同步修改 CMake 源列表、include 路径和必要的测试引用。当前有依赖 `/fp:strict`、`/fp:precise` 和禁用 LTO 的几何实现，必须保留编译属性；不能因为迁移路径而改变浮点行为。目录整理不等于马上重建整个链接库体系。

## 6. 分批实施与停止条件

| 批次 | 工作范围 | 验收/停止条件 |
|---|---|---|
| A：先理清喷嘴主线 | App 按目录分组；提取 EHD 耦合适配；喷嘴计算与写文件分离；修正 initialize/update 命名；新增简短目录职责说明 | 可顺着 App 看清一步流程；原输入不变，原喷嘴短程与一个受影响 EHD 短程结果不变；不扩大物理测试矩阵 |
| B：规范新增物理配置 | 实现版本化输入、materials 与 hydro/electricity 分工、Time/Output/Log；接通实际使用的 solvers 配置；保留旧格式 reader | 新旧输入生成等价内部配置；原数值结果不变；错误配置在求解前清楚失败；CLI、输出目录不变 |
| C：随物理主线提取公共能力 | 在修局部体积守恒、接导电液面时，将通用输运/电荷/受力实现归模块，保留 case 绑定；明确状态接受与时间层 | 原喷嘴先推进，不等待全项目迁移；原体积/电荷验收继续执行，不因结构调整掩盖失败 |
| D：其余清理 | 处理剩余 Runtime 的 Plasma 配置装载和验证 App；更新旧架构说明，必要时逐步细分模块子目录 | 单独安排，不阻塞前面两项物理任务，也不强制批量迁移历史算例 |

批次 A 内部先做无 IO 计算与显式 writer 的拆分，再提取适配和移动应用目录；不要反过来先大批搬文件。Hydro/Electricity 的全部子目录细分属于目标结构，按实际触及的文件逐步迁移，不作为 A 的完成条件。

每批拆分“纯结构变化”和“算法变化”：不在一次改动里同时搬文件、改字典、换离散、换验收口径。纯结构批次优先比较原生数值文件是否一致；有意改变的路径/元数据另列，不拿日志哈希代替数值验证。

验证控制在直接相关范围：少量解析错误检查、受影响目标编译、原喷嘴既有短程与对应 EHD 入口；触及公共组件时才加已有相关回归。无需为目录整理重跑全部曲率、厚度、网格收敛矩阵。

## 7. 算例和文档也要分清职责

- 保留现有 `cases/<name>/system` 路径，不一次搬动历史算例；用一个算例索引区分可运行标准算例、开发短程、尚未接完的完整模型。
- README 固定说明：物理模型、已实现/未实现部分、单位、运行命令、输出字段与观察重点、已知失败。介电短程不能在名称或说明上冒充完整导电电喷雾。
- 脚本和误差判据放 `tests/` 或已有工具目录，运行产物放新运行目录；case 不掺入另一套算法。
- 生效规范集中维护：架构说明、字典参考、算例指南、当前进度。历史日期报告作为证据归档，不继续充当多个互相冲突的规范入口。
- 现有[架构规则 v2](QUANTUM_Developer_Guide_ArchitectureRules_v2.tex)的分层原则继续沿用；其中旧 `apps/` 布局和字典白名单与当前实现已有差异。本方案获准实施后更新规范并标注替代关系，不在提案阶段宣称新格式已经生效。

## 8. 用户确认的实施顺序与规范固化

1. 先完成当前原喷嘴算例接线和所需验证；保留物性、边界、模型支持范围与验收口径，不能以介电 smoke 替代完整电喷雾通过。
2. 然后实施 A、B，按实际涉及的共用能力落实 C，再处理非阻塞的 D。分批迁移、原样对照，不一次搬动整个仓库。
3. 形成一份正式架构规范和一份正式字典参考，记录版本、支持状态、兼容关系及迁移方法；本计划保留为实施记录，旧规范明确标注被替代的条款。
4. 在仓库级 `AGENTS.md` 固定开发约束并链接正式规范：目录职责、依赖方向、物性与 case 分离、配置唯一来源、计算与 IO 分离、SuperLU、失败传播、输出目录和验收要求。避免另抄一套会漂移的完整规范。
5. 将能自动执行的约束纳入现有构建/测试入口：禁止新增反向层依赖和数值模块依赖 App 配置；新 schema 严格检查未知键、重复来源、单位及不支持的设置；求解器生效配置与实际调用一致；保留相关回归。存量例外要逐项登记，不能用整目录豁免掩盖新增违反。
6. 后续新增 App、模块或字典设置必须同时补职责说明、配置参考和针对性检查。不能只增加开关而不实现消费者，也不能靠注释声明已经遵守。

本节保留全套固化要求。A1 已生成正式代码组织规范、现行兼容字典参考、仓库 `AGENTS.md` 和 `architecture_dependencies` 检查；A2 已实施上述职责拆分并扩展结构检查，两条原喷嘴完整数值对照也已通过。B 批新语法、对应参考与消费者验证尚未完成，不能把现有八个运行适配例外扩大为整目录放行。

### 8.1 恢复开发的约束：复用旧求解器（2026-09-12 用户确认，即日起执行）

**复用旧求解器，恢复耦合链，修复数值问题并通用化；不是大半物理算法都重新发明。** 这是当前恢复工作的执行要求，不等待目录迁移后才生效。

- 以 `IGMPP-ILET` 的已有 NS、电场/表面电荷和耦合时序为行为基线；`passive_fed` 可作补充实现与输出对照。旧代码存在不等于该模型已在当前 App 接通，也不等于旧输出是解析真值。
- 每批明确区分“沿用的旧实现”“接口/配置适配”“有证据的数值修复”“确实缺失而新增的能力”，优先接入现有模块，避免另建平行求解器或第三套同功能电荷算法。
- 算法改变必须说明旧实现的具体问题、修改范围、时间层/守恒量是否改变及直接相关验证。不能把更换离散或分裂时序说成单纯恢复接口；旧专用常数、制造解、禁用分支和已确认错误也不照搬。
- 原喷嘴是当前主线。先解决已定位的真实边界/几何和耦合问题，再完成导电液面；不因通用化重新扩展一轮与当前阻塞无关的验证矩阵。已约定的原参考、掩膜、门限和失败项继续保留。
- A/B 结构整理仍在原算例调通后实施；正式架构规范和仓库 `AGENTS.md` 应引用本条复用原则，不以重构为由另起物理求解路线。
