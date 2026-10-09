> **状态：historical/reference（历史资料，仅供追溯）。**
> 本文原始正文保留为历史记录，其中“当前”“下一步”、TODO、旧规则与通过记录均指写作当时；不再自动产生新任务，也不能证明 V2 gate 通过。
> 现行任务唯一入口：[分模块任务清单](../roadmap/05_LAMP_Module_Task_Roadmap.md)。文档权威见[开发总纲](../00_LAMP_Development_Charter.md)，完整入口见[文档索引](../README.md)。
> 旧字典/职责说明可供查询既有实现；V2 新增或修改仍遵守现行六份规范。历史正文中的链接保留原样，目标缺失时不能作为本机已验证证据。
> 本目录中的旧设计、TeX/PDF 和索引同属历史参考；总纲指定的三份冻结方法 PDF 仍以 `docs/` 下的权威文件为准，不因本目录的历史标记而降级。

# QUANTUM 字典参考

2026-09-28：参数界面新增两种解析初值。`type analytic` 在 s∈[0,1] 上填写 XY 的 x/y 或 RZ 的 r/z 公式；`type implicit` 填写 F、seed 及可选 traceTolerance，提取所选规则零集。二者均支持 constants 和 segments，并汇入同一 Hermite 初始化，不新增运行期界面表示。完整语法、分量/边界限制和可运行示例见[解析几何输入](QUANTUM_Analytic_Geometry_Input.md)。当前接通 cchMomentMidpoint 两相 Hydro，编译及输入预检通过，物理精度待集中验收。

2026-09-19：schema 2 喷雾的界面设置使用独立 `interface.dict`，见第 5.3 节；原 `hydro.dict/Interface` 继续作为兼容输入，不能与新文件混填。模型、编号及数值含义保持，不新增版本号。

2026-09-19：普通 Hydro/EHD 新增分文件输入，运行设置归 `run.dict`、节点界面设置归 `interface.dict`，见第 15 节。旧扁平输入继续兼容；不使用喷嘴的 schema 2。

2026-09-19 默认输出目录收口：PIC 验证的 CSV 空间采样及网格 DAT 归 `output/fields`，回归报告和自力汇总归 `output/meta`；MovingEB 三项默认诊断归 `output/series`，MotionSchedule 默认为 `output/series/geometry_step.csv`，显式 PIC AMR 报告为 `output/series/amr_report.csv`。旧 SingleModule/Electricity 默认写 `output/fields/electricity_init.dat`。显式配置的输出路径保持原义，旧结果不迁移、不生成兼容副本；需要续接旧几何诊断文件时显式填写原路径。仅改变默认位置，字段、数值和采样时刻不变。

2026-09-17 补充 PIC 原格式的字段和剖面输出，见第 8–9 节；没有升级 PIC 输入版本。

2026-09-16 新增独立 Interface 的规定输入入口，实际支持范围见第 6 节；其他 App 不因此改用 `interface.dict`。

版本：B2，2026-09-12。**schema 2 已接入两条运动喷嘴路径：导电液面和理想介电液面。** 第 1–4 节保留无版本号的兼容输入，第 5 节是新输入规范；两者不能混填。普通 Hydro/EHD 的旧格式及第 15 节分文件格式、Plasma 原格式均不支持把版本号直接改为 2。B2 短程复核状态见实施记录（未迁入；原路径 `validation/organization/QUANTUM_Dictionary_B2_Versioned_Input_20260912.md`）。

## 1. 兼容格式：入口与文件归属

PIC 的 `output/series/perf_step.csv`：显式时间列为 `time_code`（代码时间），隐式为 `time_s`（秒）；计时列均为实际墙钟秒。显式初始行与普通步均为 20 列；恢复读取已知旧表头及 17 列初始行时自动兼容迁移，旧数值不作单位换算。本项不增加字典设置，见[性能输出记录](validation/organization/QUANTUM_Modular_M6_Pic_Performance_Series_20260918.md)。

`Run.app PicVerification`（别名 `PicVerificationApp`）继续从 `pic.dict` 读取 `RegressionHarness`、三种 unit-box/电场检查、三种 cut-circle 检查、`RegressionEmbeddedConductor` 和 `SelfForceTest3`。这些配置由验证 App 的 `DictionaryPicVerification` 持有，PIC 时间推进的 `PicConfig` 不再包含验证开关；普通 PIC 的物理设置仍由原 `DictionaryPIC` 解析。旧块名、PIC 内别名和物理默认值保持，未增加输入版本；默认输出位置按本页 2026-09-19 迁移说明分类。词法复用同一次扫描，回归基线仍由验证 App 记录/比较。注意旧验证 reader 的单位失败回退仍保留；本批不扩充单位表，例如 `pC` 未支持，明确库仑量应写 `1e-12 C`。见[验证输入记录](validation/organization/QUANTUM_Modular_M6_Verification_Inputs_20260918.md)。

运行适配说明（2026-09-16）：普通 Hydro 的 `SingleModule` 选择及 `Hydrodynamics/Hydro` 注册名、普通 EHD 的 App 名和 `hydro.dict/electricity.dict` 语法不变。主循环已迁到各自 App，模块字典 reader 保留在原 IO 目录；不新增版本号或运行参数。职责与对照见运行迁移记录（未迁入；原路径 `validation/organization/QUANTUM_Modular_M2c_Runtime_20260916.md`）。

普通 Hydro/EHD 的旧兼容输入仍在 `hydro.dict` 中读取 `advectInterface`、`reinitializeEvery` 和 `interfaceBoundary`；reader 单独返回 Interface 配置，运行设置由 App 单向适配。新分文件入口见第 15 节，两套来源不可混填；独立规定速度 Interface App 的格式仍见第 6 节。早期内存归属迁移见[界面输入记录](validation/organization/QUANTUM_Modular_M3_Nodal_Input_20260919.md)。

算例仍为 `cases/<name>/system/*.dict`，运行仍为 `QUANTUM.exe --case <case>`。新导电短程用 electrospray_conducting_short_rz（未迁入；原路径 `../cases/electrospray_conducting_short_rz/README.md`），旧介电短程另用 `electrospray_moving_dielectric_short_rz`。不要直接启动原正式喷嘴的长终时来代替开发短程。

| 文件 | 兼容模式实际用途 | 与新版的区别 |
| --- | --- | --- |
| `geometry.dict` | 几何、相、水平集、边界 ID；由公共几何读取器消费 | 不在本次改变格式/编号 |
| `mesh.dict` | 初始网格及细化规则；由公共网格读取器消费 | 一般动态跨级 AMR 仍未在喷嘴验收 |
| `run.dict` | `Run.app` 选择 App；旧 `Electrospray` 块控制时间、模式和输出 | 新版改用 Time/Output/Log 等块 |
| `electricity.dict` | 旧 `ElectrosprayPhysical` 同时存物性、流体/电学边界及原数值容限 | 新版分归 materials/hydro/electricity/solvers |
| `units.dict` | 单位模式和参考尺度；喷嘴实际使用 SI 转换 | 不把参考粒子数密度当作流体质量密度 |
| `normalization.dict` | 归一化体积测度；本 RZ case 用 dimension 3 | 不是三维 Cartesian 网格 |
| `solvers.dict` | 旧共享配置由运行层加载 | 喷嘴 fitted/P3/NS 当前直接指定 SuperLU，并使用物理配置中的原容限；并非所有此文件开关都被这些调用消费 |

本表只描述兼容模式。迁移到第 5 节的新格式时，必须同时迁移运行、材料、流体、电学和线性系统设置，不能只改文件名；几何、网格及其编号保持原样。

## 2. 喷嘴运行控制：run.dict

```text
Run
{
    app Electrospray
}
Electrospray
{
    physicsMode coupled
    outputName electrospray_init
    writeFields true
    dt 5e-6
    finalTime 1e-4
    writeDt 5e-5
    restartLoop -1
    initialAdaptIterations -1
    initializeLevelSets true
}
```

| 字段 | 默认值 / 当前含义 |
| --- | --- |
| `physicsMode` | 默认 `coupled`；实际选择下表的推进范围，不只是标签 |
| `dt`, `finalTime`, `writeDt` | 默认分别 `5e-6`, `20`, `0.1`；均为代码时间，**默认终时不是 20 步** |
| `capillaryHalfWidthCells` | 默认 3；数值界面力带半宽除以 `max(dr,dz)`，不是物理界面厚度 |
| `outputName`, `writeFields` | 默认 `electrospray_init`, true；输出按 fields/series/meta 分工 |
| `logEvery`, `verboseDiagnostics` | 默认 1, false；摘要间隔与详细报告，不降低逐步收支 CSV 频率 |
| `restartLoop` | 默认 −1；当前喷嘴拒绝非负重启请求，不把目录缓存当作可用续算 |
| `initialAdaptIterations` | 默认 −1，使用 mesh 既有初始迭代设置；旧别名 `initialIterations` 保留 |
| `initializeLevelSets` | 默认 true；当前喷嘴要求从算例几何初始化，false 会拒绝 |

| 模式 | 实际推进 | 电学选择 |
| --- | --- | --- |
| `frozenMaterialHydro` | 固定初始材料的压力驱动 NS，无毛细/电力推进 | `cellHarmonic` 初态绑定 |
| `frozenInterfaceCapillary` | 固定液面的毛细力—NS | `cellHarmonic` 初态绑定 |
| `frozenInterfaceEhd` | 固定液面的毛细力＋q=0 极化力—NS | `cellHarmonic` |
| `movingDielectric` | NS、运动液面、物性和 q=0 介电场；不推进导电电荷 | `cellHarmonic` |
| `coupled` | NS、运动液面、物性、旧电流驱动的液面电荷及新 P3 电场/Maxwell 力 | 必须显式 `sharpHermiteP3`，SI、`leakyDielectric` |

2026-09-19：上表五个名称在兼容输入边界转换为明确的运行组合，原拼写、默认 `coupled` 和输出标签保持；无效名称在适配阶段拒绝，早于物理输入和求解。运行循环不再以字符串选择组合。默认首步输出仍只对 `coupled`、`movingDielectric` 开启；旧配置对象的显式 `writeFirst` 覆盖仍生效，但旧 `DictionaryElectrospray` 并未读取该键，不应把它写作旧字典的新能力。

第 5 节 `schemaVersion 2` 的 `Coupling` 仍按电学模型选择：`leakyDielectric` 为运动导电，`perfectDielectric` 为运动介电。固定材料/固定界面的三个组合现有独立 `Verification` 定义，见第 5.2 节；它显式要求 `perfectDielectric/cellHarmonic`，不会将导电模型静默替换为介电模型。旧入口继续兼容，消费者清单及退出条件见[运行组合记录](validation/organization/QUANTUM_Modular_M7_Spray_Workflow_20260919.md)，新旧小算例的实际对照见[验证入口记录](validation/organization/QUANTUM_Modular_M7_Fixed_Verification_20260919.md)。

读入位置是 [Compatibility/DictionaryElectrospray](../src/Apps/Electrospray/Compatibility/DictionaryElectrospray.cpp)，模型/支持范围由 [App](../src/Apps/Electrospray/ElectrosprayApp.cpp)及[导电预检查](../src/Apps/Electrospray/Coupling/ElectrosprayConducting.cpp)核验。旧运行 reader 对未知键、重复项、部分布尔值和数字尾随文本的处理尚未统一为严格语法；不要依靠拼写错误触发默认值。严格检查已接入第 5 节的新格式，不对全部历史 reader 作追溯性改写。

## 3. 喷嘴物理输入：electricity.dict

旧块名为 `ElectrosprayPhysical`。相、边界、几何点 ID 都为 **1-based**；带 `_SI` 的量使用名称所示 SI 单位，无后缀的相对介电常数和误差容限为无量纲量。不要把这个块直接用于独立 Electricity/EHD 的另一种 reader。

| 字段组 | 含义与当前原喷嘴值 |
| --- | --- |
| `model`, `interfaceDiscretization` | 导电短程为 `leakyDielectric`, `sharpHermiteP3`；介电兼容模式的实际推进范围以上节为准 |
| `fluidInterfaceLevelSet`, `movingBoundary` | 选水平集及真正运动的自由边界；本例为 `"1"`, 4，不移动固定湿壁 3 |
| `liquidPhase`, `liquidRelativePermittivity`, `liquidConductivity_SI`, `liquidDensity_SI`, `liquidDynamicViscosity_SI` | 相 2；10、1.4 S/m、1240 kg/m³、0.038 Pa·s |
| 相应 `gas*` 五项 | 相 3；1、0 S/m、1.184 kg/m³、1.849e−5 Pa·s |
| `conductorPhase`, `conductorPotential_SI` | 导体相 1、300 V |
| `surfaceTension_SI`, `initialSurfaceCharge_SI` | 0.0452 N/m、0 C/m²；当前喷嘴初态要求零自由面电荷，之后导电模式实际推进 |
| `domainBoundaryPotential_SI <ID> <V>` | 本例边界 9，0 V 收集极 |
| `domainBoundaryHomogeneousNeumann <ID>` | 其余活动流体域电通量边界逐项声明；本例 5、8、10、11、12 |
| `supplyPressure_SI <ID> <Pa>` | 边界 5、45.2 Pa；指定压力，不指定入口速度或流量 |
| `outletPressure_SI <ID> <Pa>` | 边界 9、0 Pa；与电学接地是独立条件 |
| `outletIncomingPhaseId` | 本例 3，允许环境气体回流；缺省 −1 的旧 pressureOutlet 拒绝回流 |
| `ambientPressure_SI <ID> <Pa>` | 边界 8、0 Pa 环境气体开放边界，可有多条 |
| `solidWallModel`, `incomingMomentumModel` | `noSlip`, `interiorContinuationZeroShear`；不等同于给定开口速度 |
| `contactLineModel`, `pinnedContactPoint <ID>` | `pinned`、几何点 6；接触位置固定，接触角不被规定 |
| `linearResidualTolerance` | 本例 1e−10；B1 兼容 reader 在读入时展开为独立线性/守恒参数，电荷收支仍取原 16 倍。当前文件仍只有这个旧输入字段，不能在旧块中填写新内部成员名 |

此块 reader 会拒绝未知键、重复/冲突数据、缺少必填值和不合法数值；几何 ID 是否存在及模型组合能否运行由绑定层核验。[实际解析器](../src/Apps/Electrospray/Compatibility/ElectrosprayPhysicalConfig.cpp)仍在 Compatibility，尚未迁成公共材料配置。

## 4. 单位、求解器及其他入口

2026-09-18 职责更新：旧通用 `units.dict` reader 现位于 Core/IO，旧 `normalization.dict` reader 和参考物种/等离子频率/热尺度解析现位于 Plasma/IO；仍由应用入口选择独立文件或 PIC 嵌入块、绑定几何闭合并输出原报告。输入版本、默认值和回退规则不变，见[归一化输入记录](validation/organization/QUANTUM_Modular_M6_Normalization_Inputs_20260918.md)。schema 2 的严格 manual 输入不因此开放物理归一化策略。

旧 SI `manual` 保留显式派生尺度；`plasmaFrequency`（含原别名）按参考数密度、质量和电荷确定时间尺度；`thermal` 按参考温度确定电势及时间尺度。电子/质子名称及正的显式电荷、质量优先级保持；未知物种警告后回退电子，无效参数沿用原警告/回退，CodeUnits 下物理策略警告后忽略。非 manual 策略仍重置派生尺度，由 UnitsContext 完成原构造；其中 UnitsContext 的旧电荷尺度回退和 NormalizationContext 的派生电荷尺度是不同字段，不能在迁移中合并或暗改。具体数值与失败范围见上述记录。

当前喷嘴 `L0=1e−6 m`、`T0=1.6563100047570824e−7 s`、`Phi0=50 V`；代码时间乘 T0 才是秒。`rho0` 在该单位字典中是参考**电荷密度**，不是液体质量密度；流体密度来自 `liquidDensity_SI/gasDensity_SI`。RZ 体积含完整回转测度，不重复乘 2πr。

兼容模式的 `PoissonCellCentered` 旧块中即使保留 `useCG/useIC`，也不能宣称 fitted/P3 调用了这些算法。本项目实际要求 SuperLU；修改旧 solvers 字典不等于改变了喷嘴所有线性系统。新格式已按第 5 节分别接入实际线性系统和守恒门限。

B1 新增 `output/meta/<outputName>_numerics.json`：记录 14 项数值参数的有效值、来源、倍率、默认状态、是否由选中模式消费及调用说明；控制台同步提示旧共享 solver 选项在喷嘴中未使用。该文件只是求解前的**数值控制**记录，不是完整物性/边界生效配置，也不是求解成功证明。内部独立参数和实际验证状态见 B1 记录（未迁入；原路径 `validation/organization/QUANTUM_Dictionary_B1_Numerical_Controls_20260912.md`）。

独立 Hydro、普通 EHD 的当前文件/时间控制见[算例指南](QUANTUM_Hydro_EHD_Case_Guide.md)。**旧格式喷嘴**仍不支持只读 `--check-case`，退出 2；**schema 2 喷嘴**已实现该检查。Plasma 的现行输入不在本次迁移范围，不强制批量补 materials 或修改 PIC/MCC 算例。

## 5. schema 2：材料、模块和运行控制分开

本机可运行输入见[导电喷嘴](../cases/spray_implicit_smoke/system/)和[理想介电喷嘴](../cases/spray_dielectric_smoke/system/)，均为迁移后重建的三步小算例。原机长程输入和物理验收状态另见主计划。

### 5.1 支持范围及语法

`Run { schemaVersion 2; app Electrospray; }` 显式选择新入口。未写版本号保留原兼容路径；不支持未知版本或版本 2 的其他 App/别名。分文件格式使用以下 10 份字典；保留内联 Interface 的旧 schema 2 输入仍使用原 9 份。不能掺入 PIC/MCC 或旧 `ElectrosprayPhysical` 块。

文件集合由 Electrospray 的版本化 reader 检查；通过主程序或直接调用该 reader，均拒绝额外文件和子目录。各物理块继续由下列所属 reader 解析，旧九文件入口的 Interface 来源互斥规则保持。

| 文件 | 内容 / 实际读取者 |
| --- | --- |
| `run.dict` | Run、Time、Output、Log、可选 Restart，以及 Electrospray 中互斥的 Coupling 或 Verification；Runtime 读取通用控制，App 绑定工作流 |
| `materials.dict` | 相物性和相间表面张力；Core/Materials 值对象、Core/IO 读取器 |
| `hydro.dict` | 流体模型、初态、压力边界、载荷带宽及守恒门限；Hydrodynamics/IO |
| `interface.dict` | 液面名称、移动边界、输运与钉扎位置；Interface/IO，具体格式见第 5.3 节 |
| `electricity.dict` | 电模型、电极、初始面电荷、电荷推进及电通量门限；Electricity/IO |
| `solvers.dict` | 五个具名线性系统；Solvers/IO 读取，App 传到已有真实调用 |
| `units.dict`、`normalization.dict` | 显式 SI 参考尺度及 manual/RZ 体积闭合；复用 UnitsContext/NormalizationContext |
| `geometry.dict`、`mesh.dict` | **仍用原公共格式和解析器**；不重解释弧半径、编号或网格设置 |

新块名用大驼峰、字段名用小驼峰，区分大小写。标量一行一个，或用分号分隔；块的左花括号可以放下一行。支持 `#`、`//` 注释及带引号的名称；新 ID 列表写 `(2 3)`，不加逗号。未知项、重复项、缺必填项、错误布尔值、数字尾随垃圾、错误量纲均拒绝。布尔值只接受 `true/false`。新解析错误给出文件、行号、字段路径、原值和期望；跨几何绑定错误另明确对象 ID/条件，旧几何解析器的诊断格式不在本批重写。

当前能力为无旋流 RZ、固定拟合固壁上的一条开放自由液面，端点在轴上或已声明的固定接触点。同层叶网格及第 5.6 节所列原生 Q4 粗细网格/运行时 AMR 已接通。**本喷嘴所需的开放边界粗细动量闭合已接通，不再要求侧面/顶部整条同层加密带。** 公共 Hydro 按实际网格构造共有通量和压力相容的动量行，无新增字典求解分支；原固体隔离、压力模板、2:1 与几何身份检查仍执行，不能据此宣称任意壁面角区均已支持。仅从字典初态启动；不支持重启、额外运动固壁或周期边界。`--check-case` 验证配置和几何身份，不构造网格、不求解、不写 output/cache；网格精确对齐、曲率质量和时间推进稳定性仍在实际运行中检查。

### 5.2 一个主时钟、输出和日志

```text
Run { schemaVersion 2; app Electrospray; }
Time { units code; dt 5e-6; endTime 1e-4; }
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
Log { everySteps 1; level summary; }
Electrospray { Coupling { electricUpdate everyStep; } }
```

上述通用控制字段均必填。`Electrospray` 必须且只能包含一个 `Coupling` 或下述 `Verification` 块。`Time.units code` 接受裸代码时间；选 `SI` 时 `dt/endTime` 必须写 `数值 + 时间单位`。`endTime` 必须是固定 `dt` 的整数倍，沿用原 `1e-12 * endTime` 检查，不截短末步。`timeInterval` 的量纲继承 `Time.units`；`stepInterval` 的间隔必须是正整数。

总步数上限仍为 `INT_MAX`；独立 Electricity 为初始求场另留一个编号，上限为 `INT_MAX-1`。`writeFinal true` 按 App 确定的最后一个接受步写出，前提是 `Fields.enabled true`；允许范围内的结束时间舍入不会漏掉末帧。例如 `dt 0.01; endTime 0.030000000000000002;` 仍推进三步并在第 3 步写末帧。数值时间使用原来的 `step * dt`，不会为了匹配输入文字补短步或移动物理状态。修复与对照见[固定步数记录](validation/organization/QUANTUM_Modular_M7_Fixed_Step_Schedule_20260919.md)。

固定状态验证使用以下定义替换整个 `Electrospray` 块，其余材料、单位、几何、电学、流体和求解器字典保持明确配置：

```text
Electrospray { Verification { problem fixedInterfaceDielectric; } }
```

| `Verification.problem` | 实际求解 | 保持固定的状态 |
| --- | --- | --- |
| `pressureDriven` | 压力驱动 NS，不加载毛细或电力 | 初始网格、材料与电场 |
| `fixedInterfaceCapillary` | 初始毛细载荷下的 NS | 初始界面、材料与毛细载荷 |
| `fixedInterfaceDielectric` | 初始毛细＋q=0 介电极化载荷下的 NS | 初始界面、材料、电场及载荷 |

三者仍执行原初始电场求解，要求 `perfectDielectric`、`cellHarmonic`、零初始面电荷及 `SurfaceCharge.model none`；不推进界面或电荷，不支持运行时 AMR 和检查点续算。它们消费同一 NS 时间循环，不能解释为成熟喷射或界面平衡验证。`Coupling` 与 `Verification` 同时存在、两者均缺失、未知/重复 problem 或混入旧 `physicsMode` 均拒绝。

正式输入的 `writeFirst` 始终显式填写，不按验证组合隐式改变；若要复现旧冻结入口，设为 false，并按原 `writeDt` 配置时间输出间隔。生效报告记录 `Electrospray.Verification.problem` 的来源，标出不参与该固定组合计算的运动/体积修正等设置；输出文件继续使用历史标签以便前后比较。参考[复用现有喷嘴输入的三种验证](../cases/spray_dielectric_smoke/README.md)。

三个写场标志分别控制初态、强制首步和强制末步；`false` 不会取消另一个定期输出规则已经命中的帧。`enabled false` 只关闭空间场，逐步收支和失败信息继续保存。目录仍为 `output/fields`、`output/series`、`output/meta`；`tecplot` 是格式，不是目录。

2026-09-19 输出调度统一：Electrospray、独立 Interface/Electricity 和普通 Hydro/EHD 共用按接受步号选择输出的规则。`timeInterval` 从零开始计间隔，每次跨过一个或多个间隔边界时写该接受步，最多一帧；间隔不大于 dt 时每步写一帧，不补算中间物理状态。边界比较沿用电学路径的 `1e-12` 个间隔舍入容差，强制首末帧不改变间隔原点。例如 dt=0.01、间隔=0.030000000000000002 时，第 3 步即按定期间隔写出，不依赖 writeFinal。极小但合法的间隔也不会因时间/间隔除法溢出而漏帧。普通 EHD 首步初始电场的额外输出选择见第 15 节；旧喷嘴的极近浮点边界输出也按此公共规则收齐，数值步进不变。

`Output.prefix` 命名整个计算过程，不表示初态。正式 Electrospray 使用 `electrospray`，例如 `electrospray_conducting_step000001_flow.dat`；旧兼容默认值和历史文件名不追溯修改。

`Output.Fields` 另支持可选 `profile physical;`（默认）或 `profile diagnostic;`。physical 保留 flow/interface/forces；导电电势、积分电荷及导电电荷累积率（A）已经在 flow，界面两侧电场在 forces，不另写重复 electric 网格、原生面/单元和初始化中间态明细。diagnostic 按同一写场时刻额外保存这些明细。逐步物理收支、求解失败和必要边界分组电流仍保留；该选择独立于 `Log.level`，不改变求解算法或检查门限。旧无版本 reader 的原输出行为暂不变，推荐使用 schema 2。实际批次核对见输出精简记录（未迁入；原路径 `validation/electrospray/QUANTUM_Electrospray_Performance_And_Output_20260912.md`）。

physical 导电帧清单的 `electric_file` 指向同帧 flow（该文件确实包含电学量），`face_flux_file` 为空表示未请求原生面明细；不能据此推断计算没有电通量。Tecplot 可直接读取 DAT，默认不另生成一份 PLT 副本。

电喷雾的逐步收支和 AMR 表统一写入 `output/series/*.dat`，采用 Tecplot XY Line 的 ordered POINT 格式，不再默认同时写 CSV；不需要新增字典设置，旧兼容 Electrospray 的这些表也使用相同 writer。`Output.Fields` 仍只控制空间场，不控制逐步收支。`*_stepNNNNNN_boundaries.dat` 保持原保存帧频率，是单帧边界编号分组表；不增加电流采样或积分。`meta/*_series.csv` 等字符串文件清单仍保留 CSV，历史输出不改写。读取说明见[算例指南](QUANTUM_Hydro_EHD_Case_Guide.md#逐步收支用-xy-line-查看)。

`Log.everySteps` 控制流体接受步摘要，首末步和失败始终报告。`summary` 保留阶段/接受状态/收支；`detail` 另启用逐轮体积修正进度和详细场诊断。两种级别均不向控制台输出计算耗时统计，物理时间照常显示；AMR 表中的 `wall_seconds` 列保留。日志降频不改变计算或逐步收支 DAT。可选 `Restart { model none; }`，省略也为 none；其他重启模型拒绝。初始几何细化次数只来自 `mesh.dict`，不在新 run 重复设置；液面始终按几何初始化。

#### 启动选择与覆盖旧结果

新版 Electrospray 在交互终端启动时询问：`1` 从初态重新计算、`2` 接上次计算、`0` 取消；回车或输入结束按取消处理。选 `1` 会先清空**当前算例 output/ 内的全部旧结果**，包括旧的后期帧和其他前缀，不仅覆盖同名文件；删除后不能由程序恢复。system/、cache/ 及算例其他文件不动。需要保留结果时先使用新的 runs 目录。日志放算例根目录的 run.log，不要把正在写入的日志放到待清理的 output/ 内。

脚本可用 `QUANTUM.exe --case <运行目录> --start new` 显式确认重算覆盖；`--start ask` 强制询问，也可从标准输入读取选项。不传启动参数的非交互运行只在没有旧输出时自动开始；已有输出则在初始化前报错，不默认删除。`--check-case` 始终只读，即使同时传入 `--start new` 也不清理。其他 App 和旧兼容喷嘴不消费这些新启动参数。

`2` 或 `--start resume` 当前明确报“断点续算尚未支持”并保留旧结果；`Restart { model none; }` 的语法和含义不变。DAT/CSV 不是完整续算状态，真正续算另批实现。新版本在同一算例上使用 `.electrospray.lock` 进程锁；Windows 正常退出或进程结束自动释放并移除锁文件，POSIX 保留零字节锁文件但释放系统锁。不要手动删除正在使用的锁来并发重算同一目录。

### 5.3 材料和流体

`Materials/Phases/<标签>` 每相必填 `phaseId`、`density`、`dynamicViscosity`、`relativePermittivity`、`conductivity`；当前必须恰为 Hydro 指定的两个流体相。密度、黏度、相对介电常数为正，电导率非负。导体固相在电极绑定中指定，不冒充第三种流体物性。

`Materials/Interfaces/<标签>` 必填 `phases (2 3)` 和非负 `surfaceTension 0.0452 N/m`；当前恰为所选液气相对，同一相对不能反向重复。

| Hydrodynamics 块 | 必填设置 |
| --- | --- |
| Model | `formulation incompressibleTwoPhase`、`coordinateSystem axisymmetricRZ`、`liquidPhase`、`gasPhase` |
| Initial | `model rest`；初始速度和压力为零，指定压力由 NS 施加 |
| BoundaryConditions/Supply | `boundaryId`、`pressure 45.2 Pa`；入口回流材料为液体，不指定流量 |
| BoundaryConditions/Outlet | `boundaryId`、`pressure 0 Pa`、`incomingPhase 3` 或 `none`；none 保留拒绝未配置回流的出口 |
| BoundaryConditions/Ambient/各标签 | `boundaryId`、`pressure 0 Pa`；可零项或多项，入流材料为气体 |
| BoundaryConditions | `solidWallModel noSlip`、`incomingMomentumModel interiorContinuationZeroShear` |
| Numerics/InterfaceForce | `halfWidthCells 3`；数值力带半宽，不是物理界面厚度 |
| Numerics/Conservation | `momentumBalance`、`localContinuity`、`outletAdmissibility`、`phaseFluxBalance`；本例均 `1e-10` |

`interface.dict` 只包含一个 `Interface` 块，由 Interface/IO 直接读取。以下值对应本机三个喷雾小算例，实际编号须匹配各自 `geometry.dict`：

```text
Interface
{
    levelSet "liquid";
    movingBoundary 4;
    transport activeHermiteMaterial;
    ContactLine { model pinned; points (3); }
}
```

这五项均必填：液面名称、正的移动边界 ID、既有输运模型、pinned 接触模型及无重复的正 ID 列表。钉扎位置不指定接触角，未开放其他输运/接触模型。App 绑定实际几何，AMR 与生效配置报告使用同一 `FittedInterfaceCaseConfig`；`Numerics/InterfaceForce` 的力带宽度仍归 Hydro。

迁移时将原 `Hydrodynamics.Interface` 块原样移入新文件，其他键不变。没有 `interface.dict` 时仍读取旧内联块；两处同时存在、两处都缺失或新文件含未知块均拒绝。新报告记录 `interface.dict:Interface.*`，兼容输入保留 `hydro.dict:Hydrodynamics.Interface.*` 来源。此格式适用于 schema 2 的 Coupling 和三种 Verification；后者仍将未消费的 transport 标为未使用。独立规定速度 Interface App 和普通 Hydro/EHD 各自使用第 6、15 节的格式，不按文件名混用模型。

所有误差容限是 `(0,1)` 内的无量纲数。动量、局部连续性、相通量、出口条件各用独立参数，保持原分母和机器精度项。所有压力口已给入流状态时，`outletAdmissibility` 在配置报告中明确标为本模式未用。

### 5.4 电学和线性系统

| 选择 | 理想介电运动液面 | 导电运动液面 |
| --- | --- | --- |
| Electricity/Model/model | `perfectDielectric` | `leakyDielectric` |
| Model/interfaceDiscretization | `cellHarmonic` | `sharpHermiteP3` |
| SurfaceCharge/model | `none` | `explicitConductionMaterialTransport` 或 `implicitConductionMaterialTransport` |
| SurfaceCharge/closureTolerance | **不可填写** | 必填，本例 `1.6e-9`；不再由线性容限派生 |
| Numerics | `potentialBounds`、`polarizationConsistency` | `potentialBounds`、`reconstructedPoisson`、`gaussClosure` |

Numerics 本例均 `1e-10`。不支持交叉混搭；初态 `Initial { surfaceChargeDensity 0 C/m^2; }` 当前只接受零自由面电荷。显式模式仍用旧电流推进；隐式模式在同一新液面上做电荷—电场欠松弛迭代，首步后向 Euler、后续旧工程两层时间式。两者都沿共用材料轨迹搬运电荷，新电场/载荷供下一 NS 步。

仅隐式模式必须填写 `SurfaceCharge/Iteration`；其他模式填写该块会报错：

```text
SurfaceCharge
{
    model implicitConductionMaterialTransport
    closureTolerance 1.6e-9
    Iteration
    {
        underRelaxation 0.1
        maxIterations 500
        densityChangeTolerance 4.4270939064e-10 C/m^2
        relativeEquationTolerance 1e-10
    }
}
```

`underRelaxation` 在 `(0,1]`，`maxIterations` 是正整数，`densityChangeTolerance` 为正 SI 面电荷密度，`relativeEquationTolerance` 在 `(0,1)`。密度容限在命名 App 中按 `epsilon0*Phi0/L0` 转为代码量，不沿用无量纲默认值冒充 SI。相邻变化量和完整离散方程必须同时通过；电极收支独立使用 `closureTolerance`。设置及消费者写入配置来源；计算中不输出内迭代临时场。原 `dt*max(sigma/epsilon)<=0.25` 限制暂保留，不能把该模式视为大步长精度已经验收。标准入口见隐式导电短程（未迁入；原路径 `../cases/electrospray_conducting_implicit_short_rz/README.md`）。

`BoundaryConditions/Conductor` 填 `phaseId 1; potential 300 V;`。`BoundaryConditions/Domain/各标签` 用 `boundaryId` 和 `type` 绑定：`dirichlet` 必须提供 `potential 0 V`；`homogeneousNeumann` 不允许附带 potential。每个活动流体域电边界必须显式声明，不能与压力边界设置混为一项。

`solvers.dict` 格式为 `Solvers { <系统名> { backend SuperLU; relativeTolerance 1e-10; } }`，恰有五项：

| 系统名 | 实际消费者 |
| --- | --- |
| ElectricPotential | 初始、q=0 更新或 P3 电势线性求解 |
| Momentum | 拟合域 NS 动量预测线性求解 |
| Pressure | 拟合域 NS 压力修正线性求解 |
| CoupledVelocityPressure | 拟合域 NS 速度—压力耦合修正 |
| InterfaceVolumeCorrection | 原逐单元体积修正中的线性子问题 |

只接受 SuperLU，不接受 SparseLU、CG、IC 或无消费者的迭代次数字段；未发生某个条件分支时不代表其门限改变了其他系统。体积修正的非线性停止准则仍保持原值，没有偷换成线性残差。

### 5.5 单位、生效记录及迁移步骤

物理字段需显式单位。常用支持：长度 `m/cm/mm/um/nm/km`，时间 `s/ms/us/ns/min/h/hr`，电势 `V/mV/kV`，压力 `Pa/kPa/MPa/bar/atm/Torr/mTorr`，密度 `kg/m^3`、`kg/m3`、`g/cm^3`，黏度 `Pa*s`、`N*s/m^2`，电导率 `S/m`，表面张力 `N/m`，面电荷 `C/m^2`、`C/m2`。单位随字段量纲检查，不能用 V 代替 Pa。

当前新入口使用 `Units.mode SI`，必填 `L0/T0/Phi0/rho0/n0/writeSIOutputs`，其中 writeSIOutputs 必须 true；派生速度/电场尺度复用公共 UnitsContext。`rho0` 支持 `C/m^3`、`C/cm^3`，是保留的体电荷参考；`n0` 支持 `1/m^3`、`1/cm^3`，是保留的粒子参考。喷雾计算不使用这两项，来源报告均标为 `used=false` 并说明原因；仍按现有格式读取和校验，改变它们不会调整喷雾电荷尺度。面电荷密度实际使用 `epsilon0*Phi0/L0`，完整 RZ 电荷使用 `epsilon0*Phi0*L0`；公共 NormalizationContext 的派生体电荷尺度为 `epsilon0*Phi0/L0^2`。`Normalization { dimension 3; strategy manual; }` 表示完整 RZ 回转体积，不是三维笛卡尔网格。其他归一化策略、额外手动派生尺度本入口拒绝。

正常运行写 `output/meta/<prefix>_resolved.json`：schema、原请求电模型、实际模型、转换后的生效值、参考尺度，以及新输入每项的来源文件/行号/原值/单位/默认状态/消费者。理想介电模式未消费材料电导率；该事实与未用的体电荷/粒子参考、出口检查单列。该文件是**求解前配置快照**，不是时间步成功证明；目录或写入失败会传到 App 退出码。旧格式另保留 B1 的 `_numerics.json`，新格式不伪造旧统一容限来源。

迁移时先复制到新 case：保留 geometry/mesh/units/normalization 四文件，按上述归属拆分其余设置；不同时保留旧时间、旧物性块或旧 PoissonCellCentered。旧统一 `1e-10` 要逐项填入原消费者，电荷门限明确填 `1.6e-9`。先 `--check-case`，再实跑；核对转换后配置和原生结果。本机复用 `interface_state_contract`、`electrospray_mesh_contract`、各 `spray_*` 与 `architecture` 检查；原机未迁入的测试名称不代表本机覆盖。字典文字和来源位置允许改变，生效参数与物理输出仍使用原比较判据。

### 5.6 运行时网格自适应：mesh.dict

`RuntimeAdaptation` 是公共网格字典的可选顶层块，省略时保持固定拓扑。当前实际消费者为新版导电 Electrospray 的原生 Q4 粗细网格路径；不是给任意 App 打开的通用开关。接线与实跑状态见动态 AMR 第 8 节（未迁入；原路径 `validation/electrospray/QUANTUM_Electrospray_Dynamic_AMR_20260914.md#8-第四批电喷雾时间循环的联合准备与提交`）。

```text
RuntimeAdaptation
{
    enabled true
    everySteps 1
    maxSweeps 2
}
```

- `enabled` 必填，只接受 `true/false`。关闭时块内只填写 `enabled false`，不保留无效控制项。
- `everySteps` 为正整数，默认 1：每接受这么多物理步、且还有下一步要推进时，检查一次移动加密带。不是时间步长，也不是场输出间隔。
- `maxSweeps` 为 1–32 的整数，默认 2：一次检查最多调用原四叉树细粗化多少轮；每轮仍可在原 `adapt()` 内处理多个层级。它不是物理迭代次数，也不保证必须改变网格。
- 块及字段不可重复；未知键、非法数值、缺括号和嵌套块明确拒绝。当前采用多行块、每行一个键值的写法。

必须同时启用 `adaptiveMesh true`、`GeometryRefinement/enabled true`，并在原 `refineBoundary` 块中恰好绑定一条运动液面规则，例如本喷嘴的 `boundaryId 4 C_Gamma 5.0 maxLevel 4`（整条规则写在一行）。固定固壁、轴线、供液入口的局部细网格区域继续用原 `refineRegion`；开放气体边界不要求专门加密带。`C_Gamma` 按原网格对角尺度定义加密带，不是毛细力数值厚度。运行时用当前液面替换该边界的初始折线，不重新初始化水平集。

当前要求 `leakyDielectric + sharpHermiteP3`、明确回流物性的压力储库、初始化后仍有粗细叶单元的 Q4 表示。`--check-case` 只检查配置与几何绑定，Q4／网格质量须在初始化和实际候选中检查；同层 Q3 初态不能只加一个开关就变成 Q4。没有显式接入运行时调度的其他工作流会拒绝启用该块，不能静默忽略。

2026-09-14 增加当前相区与流场判据，示例见正式 `cases/electrospray/system/mesh.dict`：

```text
# 写在 GeometryRefinement 内；一条相规则一行，相 ID 从 1 开始。
refinePhase
{
    phase 2 minLevel 3 maxLevel 4
    phase 3 minLevel 1 maxLevel 4
}

# 独立顶层块；不再填写另一套 adaptEverySteps 等调度。
PhysicsRefinement
{
    enabled true
    indicators (Velocity VelocityGradient)
    velocityThreshold 0.3 m/s
    velocityGradientThreshold 1e6 1/s
    hysteresis 0.25
    caps
    {
        Velocity maxLevel 4
        VelocityGradient maxLevel 4
    }
}
```

`refinePhase` 的相 ID、最低/最高层级必须有效且不可重复。初始化由公共几何标相；运行时由**当前权威液面**判定活动流体相，固体仍由固定几何绑定。不读取旧初始顶点相标签来指示已移动的液体。当前含界面单元及保护层仍优先保持不变。

Electrospray 的流场指标仅支持 `Velocity`、`VelocityGradient`，每项必须同时有对应的正 SI 门限及层级上限；可单独启用一项，不可留下另一项无消费者的门限/cap。速度支持 `m/s、cm/s、mm/s、km/s`，速度梯度支持 `1/s、s^-1`，裸数或错误量纲拒绝。`hysteresis` 必填且位于 `(0,1)`。

指标由 `Modules/Hydrodynamics/HydroFittedRefinement` 从已接受的原生速度重构读取：速度取中心及四个内部采样点的模最大值；梯度为四分之一单元偏移的内部割线估计，RZ 包含 `u_r/r`，不是另建 NS 算子或误差估计器。两项各自除以门限：比值 `r>1` 时按原 `ceil(log2(r))` 层级请求细化并受各自 cap 限制；`1-hysteresis <= r <= 1` 保持；更小时允许每轮粗化一级。流场判据取较细要求，再应用原几何规则：距离带的 `maxLevel` 限制带内上限，固定区域提供最低层级，相规则限制最低/最高层级，最后还受 2:1 邻接及保护层约束。正式喷嘴的各项最高层级统一为 4，避免互相压制；提高某项上限时须同时检查相区和边界规则。初态速度为零，只按几何/相区建网格；首个接受步后才有流场指标，不做全域最大值归一化。

每轮候选网格必须重新绑定指标编号；细化子单元继承已接受母单元的指标，粗化单元取覆盖旧叶单元的最大值，避免中心采样漏掉局部高速区。这只控制分辨率，物理速度仍走原模块转移。SI 门限、cap、消费者和首次生效时刻进入原 `_resolved.json`。关闭流场指标时使用只有 `enabled false` 的块；若关闭 RuntimeAdaptation，也必须关闭 PhysicsRefinement。Plasma 原有物理判据和调度未改。

本批保护所有含液面的单元、边界单元以及必要差分模板邻域，只在已确认的单相区改变层级；当前液面单元不拆分／合并，不支持破裂与接触线迁移。网格变化后的速度两层历史、电荷历史、压力、材料、电势和界面力联合准备，通过检查才提交。零时间连续性修正复用 `Pressure` 的 SuperLU 设置，不累加到物理压力，不额外推进 NS 或电荷。

共享节点支撑在复制前由 Core 的 `topologyProtection` 纳入几何保留掩膜，再调用原 `mayAdaptToTargets` 判定是否需要生成候选。没有可执行调整时，不复制网格或物理历史；需要调整时直接复制原生网格数据，不重新生成背景网格。原标量保护掩膜、完整多项式一致性、体积和电荷验收仍在实际转移中检查；失败明确返回，不再先生成一个已知缺少支撑的候选后重试。没有新增输入开关。AMR CSV 的 `protected_cells` 记录包含共享支撑的几何保留单元数，不是实际细化/粗化数量；现统一计入，不再仅在重试时计入。实跑与诊断差异见复制前判定记录（未迁入；原路径 `validation/electrospray/QUANTUM_Electrospray_Iteration_Cost_20260913.md#11-网格直接复制与共享节点保护前置2026-09-15`）。

实际网格事件仅增加 `output/series/<prefix>_amr.dat` 的一行：步数、SI 时刻、单元数量、体积／电荷／边界液体通量变化、两层散度和速度修正、动量／压力积分变化及耗时。原体积和电荷累计量不重置；场输出频率仍只由 `run.dict/Output/Fields` 决定，后续 Tecplot 帧直接包含新网格。配置来源和调度值记入原 `_resolved.json`，不另建重复空间场目录。

## 6. 独立 Interface 规定输入算例

消费者为 `Apps/Interface/InterfaceApp`。标准输入见 [interface_prescribed_cap_rz](../cases/interface_prescribed_cap_rz/README.md)。该入口不创建流场或电场，不支持直接迁入喷嘴的 SI、导电设置；保守供液的局部网格细化见第 6.5 节，其语法与第 5.3 节喷雾不同。

`system/` 只接受 `run.dict、geometry.dict、mesh.dict、interface.dict、solvers.dict` 五个文件。`run.dict` 使用 `Run { app Interface; }`（不写 `schemaVersion`），其 Time/Output/Log/Restart 复用现有 `RunControls` 读取与消费者；Time 只接受 `units code`，Restart 只接受 `model none`。Fields 只接受 Tecplot/physical，可使用 stepInterval 或 timeInterval，关闭场输出不关闭守恒历程。

`pinnedToAxis` 方法需要这五个文件；`conservative` 方法只需要前四个，必须省去未使用的 `solvers.dict`。额外文件和子目录由 Interface 输入 reader 拒绝，主程序和独立程序遵守同一规则。`--check-case` 的这类失败返回 2，不创建 output/cache；不能通过混入 hydro/electricity/units 等字典启用本入口没有的能力。

当前本机标准输入是迁移后重建的三步检查：`dt=0.01`、终止时间 `0.03`；以下 Point ID 和数值与该输入一致，不代表原机十步算例已复现。

界面字典由 `Modules/Interface/IO/DictionaryInterface` 严格解析；未知键、重复键、非法量值及重复边界组均失败。编号列表用空格分隔：

```text
Interface
{
    units code
    excludedPhase 1
    Initial { model parabolicCap; pinnedPoint 3; height 0.125; }
    Transport
    {
        model pinnedToAxis
        halfWidthCells 6
        Velocity { model axialColumn; peakSpeed 0.05; }
    }
    Boundaries
    {
        stationaryWall (3 4)
        symmetry (8)
        singlePhasePort (5 6 7)
    }
}
```

| 设置 | 数学含义、单位及约束 | 实际消费者 |
| --- | --- | --- |
| `excludedPhase` | `geometry.dict` 中被排除的固定固体相，正的 1-based ID | App Setup 调用 Core 判相和完整叶域构造 |
| `Initial.pinnedPoint` | 固定湿壁唇缘的 Point ID；实际网格上必须有该活动顶点 | Setup 绑定坐标，Interface 初始化及运动约束 |
| `Initial.height` | 正的代码长度；初始 `phi=z-z_pin-h(1-r²/R²)`，`R=r_pin` | `InterfacePrescribedData::initializeParabolicCap` |
| `Transport.halfWidthCells` | 正的支撑带宽度；该初态不是距离函数，实际模板仍由模块严格检查 | 同一 `InterfaceMotion::preparePinnedToAxis` |
| `Velocity.peakSpeed` | 有限代码速度；`u_z=U(1-r²/R²)^5`，`r>=R` 为零；`u_r=0` | 规定速度与同源完整 `2*pi*r` 原生面通量；不是 NS 解 |
| 三类 `Boundaries` | 非空且互斥的几何边界 ID；固定固体只能为 stationaryWall，symmetry 必须在 r=0 | `InterfaceTransportDomain` 壁面/对称/单相端口约束 |

几何沿用 Core 原有**多行**语法，只定义固定固体与有效域；不重复提供运动液面。初始适配器要求无周期、无 MovingEB、无运动水平集，采用均匀 RZ 叶网格。默认 `mesh.dict` 使用 `BackgroundMesh { divide (...) }` 和 `adaptiveMesh false`；保守方法可按第 6.5 节显式启用后续细化。不接受 QuadtreeMesh/RefineBoxes/GeometryRefinement/PhysicsRefinement/RuntimeAdaptation 等无消费者的块；公共固定网格默认容量仍由原读取器提供。

`solvers.dict` 只配置实际系统，沿用喷嘴相同的名称与 backend 约束：

```text
Solvers
{
    InterfaceVolumeCorrection { backend SuperLU; relativeTolerance 1e-10; }
}
```

线性容限实际传入同一体积修正服务，不改变原局部/全局物理收支门限。输入来源与消费者写到 `output/meta/configuration.txt`；完整场帧清单为 `meta/frames.csv`，不另写重复物理 CSV。`--check-case` 不验证尚未构建的离散网格；实际面覆盖、接触点对齐、液面链和支撑带在初始化/输运时检查，失败明确返回且不标作接受步。

独立钉扎 RZ 算例现提供前后时刻的规定法向通量，使用步内积分和对应守恒修正导数，不新增字典键。`Output.Fields` 选择的时刻还输出 `output/series/interface_cap_accuracy.csv`，比较输运前、输运后、修正后的解析曲率与形状；`geometry_time` 与 `accepted_time` 分别表示几何实际时间和接受步末时间。三档 B3 研究入口及限制见[算例说明](../cases/interface_prescribed_cap_rz/README.md)。这些设置不自动改变 Electrospray 的旧时刻通量能力。

### 6.1 B0 独立几何与制造拓扑（2026-09-19）

`Run { app Interface; }` 还可显式选择以下模式，与上面的规定速度输入互斥：

```text
Interface
{
    Benchmark { problem cassiniRZ; resolution 16; writeFields true; }
}
```

该模式只接受 `system/run.dict` 和 `system/interface.dict`。run 中只接受 `Run.app`，不填 `Time` 或物理求解器；制造采样时刻由问题定义。三个 Benchmark 字段都必填，未知字段、额外字典和超出 `[2,256]` 的整数分辨率会拒绝。`resolution=N` 决定网格，`writeFields` 控制界面 Tecplot 输出；指标、分量和事件始终记录。问题列表、各网格定义、返回码及运行脚本见 [B0 样例](../cases/interface_geometry_b0/README.md)。

输入由 Apps/Interface/IO 读取，解析形状与参考积分属于验证驱动；重构、连通域、状态所有权与制造父子关系由 Core/Interface 正式服务完成。它不接受或消费电喷雾材料/电荷模型，不能作为物理破碎开关。

同一 `Benchmark` 输入还接受 P2 的 `componentStepXY`、`stationaryPairXY`、`translationPairXY`、`translationXY`、`translationRZ`，定义和实测边界见[多分量运动检查](../cases/interface_components_step/README.md)。前两项验证小步运动/静止保持，后三项在 N=16 已通过三步平移；它们的时间和速度由问题定义，不能用作任意物理时间控制器。

### 6.2 B1/B2 规定速度输运（2026-09-20）

`Benchmark.problem` 可取 `zalesakXY`（B1 开槽圆盘，周期 1）、`vortexXY`（B2 单涡，周期 2）或 `vortexStressXY`（周期 8）。`resolution` 表示单位方域的每向格数，当前接受 `[32,256]`；正式研究取 64、128、256。必填字段仍为 `problem`、`resolution`、`writeFields`。仅这些问题额外接受 `cfl`（默认 0.25，范围 `(0,0.25]`）和 `endTime`（默认完整周期，必须大于零且不超过周期）。CFL 按规定速度的模长上界限制时间步；完整周期步数取 4 的倍数，包含四分之一周期输出。缩短终止时间的结果明确标为 `partial`；原 B0/短程平移输入不接受这两个字段。

该输入仍只有 `run.dict` 与 `interface.dict`，不消费 `Time`、物性或电学字典。B1 的旋转速度及解析气相流入、B2 的静止外壁均由标准问题定义；不增加喷雾字典开关。物理定义、运行命令、诊断局限与实施状态见 [B1](../cases/interface_zalesak_b1/README.md)、[B2](../cases/interface_vortex_b2/README.md)。

### 6.3 保守界面入口（2026-09-21，能力仍在收敛）

`Benchmark` 新增可选 `method`，默认 `levelSet` 保留上面原路线；显式填写 `conservative` 时调用同一 Interface 库的守恒库存、体积一致曲面/PLIC 通量、jet 预测及公共几何同步，不再调用旧几何体积修正。示例：

```text
Interface {
    Benchmark {
        problem translationRZ;
        method conservative;
        resolution 32;
        writeFields true;
        cfl 0.25;
        endTime 0.003;
    }
}
```

当前平移问题为 `translationXY`、`translationPairXY`、`translationRZ`、`stationaryPairXY`，保留其原几何/速度定义；`resolution≥16`，`endTime` 默认 `0.03`、范围 `(0,0.25]`，`cfl` 默认 `0.25`、范围 `(0,0.25]`。步数按速度模长上界和实际最小网格宽度确定。**这些问题的 N 不是每向格数**：XY 为 `4N×4N`、RZ 为 `2N×4N`，摘要另写实际单元数和网格宽度。原 `levelSet` 平移仍使用其既有固定时间控制。

`vortexXY`／`vortexStressXY` 沿用单位方域及原周期，当前通过同一空间解析通量形状与时间幅值的分离表示接入时空积分。三角函数在原解析零点附近评价，保留精确的壁面/中线零值，幅值反号处单独切分时间；未将三角速度伪装成多项式。任意随时间改变形状的非多项式回调仍没有该能力，不能仅凭采样猜测零点。**完整周期尚未验收**：原 N=32 首步直边扫掠失败已越过，固定高度模板及形状约束后，N=32、`dt=0.0075` 已运行四步至 `t=0.03`，局部图面支撑改造后末步几何/库存相对 L1 为 `0.00174513`、全步曲率有效覆盖 100%。N=64 也完成八步至 `t=0.03`，末态误差 `0.00127701`、全步覆盖 100%；两档都仍超 `1e-3` 目标，没有完整空间/时间收敛结论。此前 C 梯度法向试验已撤回，当前结果来自固定图面支撑改造。此前第四步几何失败作为历史证据保留。`zalesakXY`、制造拓扑和接触边界问题仍不接受 `method conservative`，不静默回退旧方法。

该路线输出：

2026-09-24 默认物理曲率更新无需修改字典：有效高度/同弧拟合的节点曲率经固定叶尺度 h 的弧长三阶正则化后用于供体、AMR 和保守毛细力。支持同尺度闭合 XY/RZ、双轴端 RZ、声明的物理端及单轴—物理端分支。物理端保留原始端点曲率，有限区间不强制零斜率；轴端才使用偶反射。缺失/跨尺度支撑或未声明端点保留原查询并报告能力状态。原 C、零线、法向及接触条件不因正则化改写；这不扩展接触端 AMR 能力。`curvature.csv` 的 curvature 为实际选择值，method/height/fit 仍说明原始提供者，不能据此推断曲率未经正则化。`metrics.csv` 追加 `regularized_donor_cells`；Hydro 的 `hydro_capillary_source_step*.csv` 追加 `regularized_curvature,secondary_curvature_weight,regularization_status`。末项依次为 0 未请求、1 已应用、2 不支持几何、3 曲率不完整、4 跨尺度、5 未解析转角。长期数值验收仍以实施计划为准。

- `metrics.csv`：库存与几何体积、`CΓ−C` 有符号/L1/最大单元误差、C 范围、边界收支与库存闭合、散度/CFL/限幅、曲面供体数、未归属库存、曲率覆盖及计算/输出耗时。另记录多图面加权节点数/最大图面差异、单调形状约束限制的节点数/最小导数缩放、角点数据非单调的混合单元数；这些几何诊断初始未同步时留空。另有低阶基线修复轮数、修复面数及最大单方向通量修复量（体积单位），初态为零。几何未分辨单元数及其少数相体积另由 `unresolved_phase_cells`、`unresolved_plic_volume_defect` 记录；这部分库存未删除。`blended_local_graph_nodes` 统计含局部图面的加权节点，`max_local_graph_rms`、`max_local_graph_residual`、`max_local_graph_chart_spread` 分别保留长度加权的线段矩方程 RMS、PLIC 中点到拟合图面的最大偏差和图面间标量差（按网格宽度归一化）；RMS 与最大中点偏差不是同一种残差，后者不能被平均矩拟合掩盖。曲率 L2 仅在平移解析曲率可用时填写，按有效面元物理测度加权，必须连同覆盖率阅读。
- `components.csv`：同步后的逐液体分量库存/几何差异；初始时刻尚无同步报告，不伪造分量误差行。
- `curvature.csv`：每步图段中点和轴端的实际组合曲率、方法（0 不可用、1 高度函数、2 同弧拟合）及原始质量枚举；不可用曲率留空。质量枚举见 `InterfaceHeightCurvature.h` 与 `InterfaceCurvature.h`。
- `writeFields=true`：初态/末态的 `inventory_*.csv`（原生单元 C、CΓ、体积差）和显示用 `geometry_*.dat`。后者明确命名 `raw_levelset_curvature`，不能冒充组合曲率。
- `summary.json`：成功完成请求区间为 `completed`，**`numerical_acceptance=false`**；返回 0 只说明请求步数完成。失败为 `failed` 或 `needs_refinement`，返回 3/4，保留最后已接受步/时间及 `failure.txt`。阶段指标可能超过精度门限，不以运行完成自动验收。

复用 `tools/run_interface_geometry.ps1`／`tools/run_interface_advection.ps1` 的 `-Method conservative`；旧默认不变。入口只接受 `run.dict` 和 `interface.dict`，尚未迁移生产喷嘴、NS/Maxwell 闭环、固定接触与实际流入。

### 6.4 保守涡旋的原生自适应入口（2026-09-24）

本节描述 `method conservative` 的 `vortexXY`、`vortexStressXY` 所用 `Adaptation` 块；省略时仍使用原固定均匀网格路径。参数 Hermite 方法的共享字段与粗化扩展见 §6.6，不能套用本节旧标量迁移预算和输出格式。示例：

```text
Interface {
    Benchmark {
        problem vortexXY; method conservative;
        resolution 32; writeFields true; cfl 0.25; endTime 0.03125;
        Adaptation { maximumLevel 2; cellBudget 12000; }
    }
}
```

| 字段 | 默认/范围 | 含义 |
| --- | --- | --- |
| `maximumLevel` | 必填，1–8 | 相对初始 N×N 背景网格的最大 dyadic 层级 |
| `cellBudget` | 必填，N² 以上、至多 2000000 | 候选实际叶单元数量上限；包含父节点/边/顶点的分配容量另计 |
| `timeStepBasis` | `configuredFinest` 或 `currentLeaves`，默认前者 | 固定按配置最细网格选步长，或按当前实际最细叶选步长；仅此自适应验证入口消费 |
| `maxCurvatureCells` | 0.125，(0,0.5] | 原生图段端点/中点处的最大主曲率尺度：XY 为 `abs(kappa_meridional)*h`，RZ 为 `max(abs(kappa_meridional),abs(kappa_azimuthal))*h`。不是曲率误差上界；不可用时仍请求质量加密，并保留原标量估计 |
| `minimumGapCells` | 8，≥2 | 已报告相反分支采样距离所需的网格尺度；不是全局最小间隙保证 |
| `maximumGeometryFractionDifference` | 0.02，(0,1] | 已同步状态逐单元 `abs(C_Gamma-C)` 的加密阈值 |
| `reserveBandCells` | 0，[0,64] | 必需网格确需改变时，额外预留的目标层级单元宽度；0 保持原策略 |
| `remapScalarTolerance` | 1e-12，≥0 | 零时间迁移绝对标量差上界预算 |
| `remapRelativeGeometryTolerance` | 1e-3，(0,1] | 迁移后几何—库存 L1 的预算，相对初始液体库存 |

初态显式从原解析有符号距离生成三次 nodal jets，并由该几何初始化 C；摘要标记 `initial_scalar=analytic_cubic_jets`。未截断已有 Q4 状态，固定网格默认仍使用原解析 Q4 初值，因此比较两条路径时应同时报告初始误差。默认时间步仍按配置最细层级和速度上界确定，摘要为 `time_step_basis=configured_finest_level`。

显式设置 `timeStepBasis currentLeaves` 后，按照当前最细叶和同一全时段速度上界选步长。原最细网格时钟保留为整数刻度，粗网格一步可跨多个刻度；真实加密后重新缩短该步，绝不超过加密规划已留的一步位移缓冲。完整周期的四分之一时刻及终点精确落点，不用速度正负抵消放大时间步。摘要为 `time_step_basis=current_leaves`，此时 `dt=null`；`base_dt` 是基础刻度而非固定物理步长，实际范围见 `minimum_accepted_dt/maximum_accepted_dt`。`metrics.csv` 追加每步 `dt,finest_leaf_level,clock_tick,planned_maximum_travel`。运行脚本可传 `-TimeStepBasis currentLeaves`。此选项改变时间离散和加密时刻，需要自己的实际精度证据；不自动继承固定步长轨迹的数值验收。每步仍复查已接受状态，保持共同层级界面支撑；目前只细化、不粗化。

候选使用原生 `adapt`、保守状态迁移及同一生产整步，成功后才同时替换网格/域/状态。达到层级/叶预算上限返回 `needs_refinement`；迁移预算失败返回 `failed`；旧状态、时钟和网格保留。物理候选失败不会自动无限加密或缩步。记录 `adaptation.csv` 的 requested/prepared/accepted/saturated/rejected 阶段及迁移误差；`required_level=maximumLevel+1` 仅是超上限需求的下界。

2026-09-25：`reserveBandCells` 将必需支撑与预留分开。每步仍按原曲率/间隙/几何差、模板宽度、下一步位移和 2:1 平衡求必需网格；只有该网格需要改变时才扩大预留带，不因预留带自身移动再次触发加密。预留超出 `cellBudget` 时明确退回必需计划；必需计划仍超预算则按原规则停止。它不会跳过每步检查、缩减支撑、推迟必要的升阶或触发粗化。`adaptation.csv` 追加 `required_leaf_count,planned_leaf_count,reserve_applied,reserve_budget_limited`；摘要记录 `reserve_band_cells`。脚本参数为 `-ReserveBandCells 4`，省略仍为 0；预留会增大单次网格和迁移成本，实际收益须由完整轨迹比较。

2026-09-25：RZ 同时检查环向主曲率，避免漏掉子午曲率接近零的细液柱。使用与毛细载荷相同的有效曲率来源；真正 `r=0` 轴端使用偶对称有限极限，不能用 `1/r` 判据令球滴端点无限加密，也不使用可能在颈部相互抵消的总曲率。`max_resolution_kappa_h` / `max_curvature_cells` 现记录最大主曲率尺度；两种自适应入口均追加 `max_azimuthal_kappa_h,azimuthal_requests`。改变影响 RZ 分辨率规划，不自动赋予 Hydro 保守 AMR 动量或拓扑事件能力，仍只是端点/中点采样指标。

2026-09-24 曲率判据与联合迁移更新：有效高度/局部拟合提供者与物理载荷共用曲率；原标量 Hessian 尖峰另存，不能覆盖有效提供者而单独触发加密。`adaptation.csv` 增加 `max_resolution_kappa_h,max_raw_kappa_h,raw_only_curvature_cells`；无网格变化但存在原标量单独超限时记 `raw_curvature_disagreement`。原 0.125 阈值、模板/近邻分支请求和层级上限不变。

共同曲线实际用于混合父单元细化时，当前独立 Interface 入口联合生成新子叶 C 和零时间同步几何。每个原父库存、纯父单元、时钟、表面/液体 ID 和祖先记录保持；新子叶相位可随约束曲线改变，随后复用物理步的库存驱动重构、纯相门及拓扑身份检查发布单一 Γh。双向端点/中点投影要求位移不超过当地 h/4、法向同向；这是有限采样检查，不是连续 Hausdorff 证书。没有共同曲线应用时保留原标量迁移，不因纯支撑带扩展额外同步。

新增 `geometry_synchronized,max_sync_scalar_change,max_sync_interface_shift` 分别记录是否同步、节点标量变化和采样几何位移；标量变化不能当作距离。`remapScalarTolerance` 继续约束同步**之前**的多项式网格转移，`remapRelativeGeometryTolerance` 约束最终 Γh 与 C 的 L1；`remap_shared_curves.csv` 的 `changed_pure_children` 记录共同曲线相对旧标量迁移改变了多少纯子叶参考。未增加字典开关。该迁移尚不提供电荷、动量或表面材料转移，未接触接触端/湿角 AMR；第 187 步保存状态的曲率仍约 5.34%，不声明长程精度验收。

`inventory_*.csv` 追加实际 `level,xmin,ymin,xmax,ymax`。有场输出时每次接受网格变化也写帧；`scalar_*.state` 使用 `QUANTUM_CURVATURE_CONTROL_NATIVE_CUBIC_V2`，下一行指向相对路径的 QMC2 `mesh_cache.bin`，随后保留原生节点 jet（悬挂数据为 NaN）、单元容量和 C。同一不可变网格的多帧共用首次保存的缓存，移动文件时须保留该相对依赖。它可作几何重放，**不是包含全部时间/材料历史的重启文件**；V1 均匀读取器必须明确拒绝 V2。现有 `quantum_interface_state_contract --cubic-mesh-transfer-replay <file>` 已支持 V2 的原生网格/约束身份重放。

同步后保存帧另写 `reconstruction_<step>.csv`，包括节点坐标、来源、最近 PLIC 单元、标量、混合/纯相限幅前后的导数、被拒候选值及标量修正量；这是原同步报告的诊断输出，不修改数值状态。首次出现缺失曲率时额外保留其真实上一接受帧及当前帧，方便定位阶段差异，无需逐步保存全部网格。

启用场输出且实际细化混合父单元时，另写 `remap_curvature.csv` 的 source/target 两阶段，以及 `fields/remap_inventory_<step>.csv`。这里的 step 是**待执行物理步号**，time 是迁移时刻（上一接受步的时间），target 是零时间迁移候选，并非已接受物理帧。曲率行保留查询点、原生测度、方法和实际高度模板支撑；库存保留 C/CΓ 与真实叶边界。用 `tools/analyze_interface_vortex_curvature.ps1 -Run <目录> -Step <待执行步号> -Stage source` 或 `-Stage target` 在正确物理时刻评价；不传 Stage 时仍评价已接受物理步。`tools/analyze_interface_remap_columns.ps1` 接受相同参数，限 XY 五列有效模板，独立复现列矩并作几何体积/三列对照；这些对照不替换生产库存或曲率。

回归库存保留按初始背景单元汇总的 `initial_grid_aggregate_inventory_return_l1`；自适应路径原 `relative_inventory_return_l1` 留 null，不能把汇总误差当成细叶 L1。另以解析圆直接积分每个当前叶矩形，summary 提供 `initial_native_disk_inventory_l1`、`return_native_disk_inventory_l1` 和 `return_native_disk_roundoff_indicator`，均按解析圆面积归一化；`native_inventory_reference` 明确参考为 `analytic_disk_current_leaves`。只有初始/完整返回时刻适用，短程末态的返回字段留 null。`return_metrics.csv` 追加 `native_disk_inventory_l1`、未归一化的 `native_disk_reference_volume`、归一化的 `native_disk_roundoff_indicator` 和 `max_native_disk_fraction_error`；舍入指示量不等同严格误差界。旧离散初值比较字段保持原义。

中途曲率/库存精度须用独立输运参考。`tools/analyze_interface_vortex_reference.ps1 -Run <算例目录> -Step <已保存步数>` 按实际叶边界裁剪参考多边形；`tools/analyze_interface_vortex_curvature.ps1` 使用相同参数评价保存的原生曲率及覆盖率。两者从已接受 metrics/summary 读取时刻和周期，不运行求解器、不改原结果。几何回归对称差仍从当前真实零面计算；自适应全周期验收尚待完成。`completed` 只说明请求区间结束，`numerical_acceptance` 仍为 false。

2026-09-24 细化分配更新无需新增字典项：对混合父单元，目标原生几何指导子叶分配，但每个父单元的原库存始终是等式约束。几何纯子叶在舍入闭合阶段也保持其相位；没有合适的混合接收者时报告预算内残差，当前支撑无法容纳父库存时回到原来的 C 匹配供体。`adaptation.csv` 追加 `geometric_donors`、`geometric_fallbacks`、`geometric_correction_L1`，分别报告本次迁移使用几何分配/回退的父单元数，以及几何参考为满足父库存所需的绝对体积调整量。后者不是总库存漂移；原迁移预算与饱和停止规则保持。闭合修复后同一两层轨迹完成 16 步且曲率覆盖完整，末态独立曲率误差仍为 1.291%，没有完成完整周期验收，详见实施计划最新记录。

复用驱动：`tools/run_interface_advection.ps1 -Method conservative -Resolution 32 -AdaptiveMaxLevel 2 -CellBudget 12000 -EndTime 0.03125 -Fields`。这是一项入口/接续检查；实际结果与精度限制见[唯一待办](QUANTUM_Modular_Consolidation_Plan.md#levelset-implementation)。

### 6.5 保守供液的局部细化（2026-09-24）

本节用于正式 `InterfaceApp` 的固定湿角、规定流场供液。`interface.dict` 选择 `Transport.model conservative`，`system` 只保留 run/geometry/mesh/interface 四份字典。初始网格仍均匀，随后由当前曲率、库存—几何差和支撑带共同请求局部细化。`mesh.dict` 示例：

```text
adaptiveMesh true
BackgroundMesh { divide (32, 64) }
Adaptation
{
    maximumLevel 2;
    cellBudget 16000;
    maxCurvatureCells 0.04;
}
```

`maximumLevel` 必填且为 1–8；`cellBudget` 必填、大于初始 Nx×Ny、至多 2000000，计入固体支撑叶。`maxCurvatureCells` 默认 0.125，范围 `(0,0.5]`。其余可选项为 `minimumGapCells`（默认 8，≥2）、`maximumGeometryFractionDifference`（默认 0.02，范围 `(0,1]`）、`remapScalarTolerance`（默认 `1e-12`，≥0）和 `remapRelativeGeometryTolerance`（默认 `1e-3`，范围 `(0,1]`）。前三个分辨率指标含义同第 6.4 节；几何迁移预算相对此次迁移前的液体库存。

同一 `Adaptation` 块也接受 `reserveBandCells`（默认 0，[0,64]），与第 6.4 节共享必需网格触发及预算回退语义。供液记录 `interface_adaptation.csv` 同样追加四个预留计划字段；只有后续物理步成功，预留网格才随域、通量和状态一起接受。

湿角纯扇区允许重新延拓标量，前提是原独立库存严格纯相且精确符号检查证明全部子叶相位不变；`remapScalarTolerance` 约束其余单元。原始全域标量差仍完整记录，不能把纯延拓差当作零面位移。固定接触坐标、真实父库存与时钟保持。每次候选网格按原几何重新绑定固定壁、轴及端口，端口相位继承初始声明的纯相；同一几何端口内初始相位不一致会明确拒绝。

本入口继续使用 `run.dict/Time.dt` 的固定步长，**不接受** benchmark 的 `timeStepBasis`。需要按允许的最细网格设置 dt；例如本例 `maximumLevel=2` 的完整供液段使用 `dt=0.005; endTime=0.3`。物理 CFL 和原精度门保持，失败不自动缩步或重试。省略适配时设置 `adaptiveMesh false` 并省去整个 `Adaptation` 块；旧 `pinnedToAxis` 不支持适配。

`output/series/interface_adaptation.csv` 区分 `prepared`（迁移完成、尚未接受物理步）、`accepted`（物理步通过）和 `rejected`。记录源/接受时刻、叶数、所需层级、曲率尺度、迁移误差和纯延拓单元数。只有迁移及后续物理步完成后，才同时替换网格、域、通量和状态；预算、层级或物理步失败返回 3，仍保留上一接受状态和输出。尚不支持粗化、移动接触线或拓扑事件；此能力不等于耦合 NS/电场喷嘴的 AMR 已接通。

### 6.6 参数 Hermite–GALS 开发入口（2026-09-27）

在规定速度的 `Interface.Benchmark` 中显式设置 `method parametricGals`。它使用唯一参数四次几何、独立 FV 液量和联合 H2 投影，与旧 `conservativeGals` 选择分开；不会改变 Hydro 默认。只读取 `run.dict` 和 `interface.dict`，其他文件、开槽圆盘及 RZ benchmark 输入不支持。`run.dict` 只含 `Run { app Interface; }`；该入口由标准问题定义网格和速度，不读取普通 Hydro 几何字典。

`problem` 为 `vortexXY`（周期 2）、`vortexStressXY`（周期 8）或 `reversedXY`（中心方域的强拉伸回逆问题，周期 15）。前两者沿用 [0,1]² 的原闭壁速度；`reversedXY` 使用 [−0.5,0.5]² 及零法向、非零切向壁面速度，仅支持 `parametricGals`，不能以旧 T=2 定义代替。工况、解析速度和参考见 [B1 README](../cases/interface_reversible_vortex/README.md)。

`resolution` 为每向根网格数，范围 16–256，非二进制分辨率未验收；`interfaceSegments` 可选，默认等于 resolution，显式值为 8–512，表示初始材料分片数，不是运行期分片上限。`cfl` 默认 0.25，范围 (0,0.25]；`endTime` 默认完整周期，允许 (0,周期]；`writeFields` 是必填布尔值。没有新的体积容差、失败回退或曲率滤波开关。

**2026-09-28 动态网格扩展（代码与预检，未数值验收）：** 三种参数涡旋均接受 `Adaptation`。`maximumLevel`、`cellBudget`、`timeStepBasis`、`maxCurvatureCells`、`minimumGapCells`、`maximumGeometryFractionDifference`、`reserveBandCells` 的默认/范围沿用 §6.4；几何与曲率观察取当前参数曲线，不能把旧标量的观察方式当作此方法的实现。参数方法拒绝 `remapScalarTolerance` 和 `remapRelativeGeometryTolerance`，网格转移沿用物理重构的同一 H2 与库存约束。可选 `Coarsening { hysteresis 0.5; }` 的 hysteresis 范围 (0,1)，省略则不粗化；此块只允许参数方法。只在无必要细化时尝试一代粗化，并重新核对目标支撑。示例见 [B1 AMR 字典](../cases/interface_reversible_vortex/amr/system/interface.dict)。

算例及完整命令见 [interface_parametric_vortex](../cases/interface_parametric_vortex/README.md)和 [B1](../cases/interface_reversible_vortex/README.md)。输出为 `output/fields/interface_NNNNNN.dat`、`cells_NNNNNN.dat`、`output/series/interface.dat` 和 `output/meta/summary.json`；返回形状/曲率只在完整周期结束后计算，之前为 null，缩短终时标为 partial。B1 终时曲率与法向距离用弧长加权，材料对应点位移另列；旧涡旋的曲率权重保持材料参数，摘要明确区分。场图在初始、四分之一周期、请求终时以及失败后的最后接受态写出，不是可续算状态。AMR 另写 `output/series/amr.dat`，只登记已随物理步共同提交的迁移；单元液量按各叶真实面积归一化。

**数值证据与实现范围：** 旧 T=2 的 N16/N32 完整回转在 §2.26 通过，返回曲率相对误差约 0.105%/0.043%；后续证据以[界面待办](QUANTUM_Conservative_AMR_Hermite_GALS_TODO.md)为准，不能移作新 B1 的通过结论。新 B1 均匀/AMR 入口只完成编译和字典预检，未运行物理时步。物理 prepare 拒绝时附 `output/meta/rejected_input.dat` 供诊断，不增加续算选项；物理矩形域或未提交 AMR 候选使用 V2，包含实际叶矩形和迁移标记，旧 V1 回放器不能直接读取。G1–G4 整体验收仍开放。

## 7. 独立 Electricity App（2026-09-16）

2026-09-19 错误边界更新：公共几何/网格 reader 原九处进程退出已改为异常，保留实际输入路径；调用方统一报告。正式 Electricity 的这类输入错误现返回 2，`--check-case` 不创建 output/cache。合法输入、旧宽松解析及弃用提示保持，见[错误边界记录](validation/organization/QUANTUM_Modular_M7_Dictionary_Failure_20260919.md)。

2026-09-19 输入归属更新，语法不变：Electricity/IO 的 `readDielectricDomainInput` 读取 `electricity.dict` 的物理块和 Domain，并通过 Core 的显式路径读取同源 geometry/mesh。App 只组合该模块记录与 RunControls；`dt` 作为原显式充电限制的预检参数传入，不由电学 reader 保存运行时钟。Sharp 的旧 `interfaceHalfWidth` 仍只作弃用审计，Diffuse 仍实际消费材料宽度；原先的能力限制和拒绝顺序保持，见[验证记录](validation/organization/QUANTUM_Modular_M7_Electric_Domain_Input_20260919.md)。

`Run { app Electricity; }` 调用正式 `DielectricEvolution` / P3 / SuperLU 服务，不构造 Hydro，不加载 Plasma 字典。标准输入为 [electricity_uniform](../cases/electricity_uniform/system/) 与 [electricity_charging_rz](../cases/electricity_charging_rz/system/)。历史 `Run { app SingleModule; }` 配合 `SingleModule { modules (Electricity); }` 仍走 Apps/Electricity/Compatibility 中的旧顶点初始化；没有将其原输入默默解释成新物理模型。

2026-09-19：旧初始化入口也已移除不被电场消费的临时 Hydro 物性计算。旧 `computeHydroProperties`（别名 `computeHydro`、`hydroProperties`）和 `propertyLevelSetIndex`（别名 `levelSetIndex`、`levelSet`）继续按原规则解析，但不再触发被丢弃的计算；默认值及非法整数拒绝保持。`writeFields` 和 `outputName` 仍控制原顶点场输出。此兼容入口使用原 SparseLU 初始化及旧模型，不能作为正式电学服务的物理验收，见[清理及对照记录](validation/organization/QUANTUM_Modular_M4_Legacy_Electric_Initialization_20260919.md)。

本入口只接受 `run/geometry/mesh/electricity.dict` 四份文件，全部代码单位。run 的 Time/Output/Log 使用第 5.2 节的公共运行控制，App 名改成 Electricity，Time.units 必须 code，不填 schemaVersion 或 Electrospray 块。只支持 profile physical、Restart.model none；支持 stepInterval/timeInterval、首末帧和关闭场输出，逐步收支始终保留。现有非空 output 目录会拒绝重算，使用新的 runs 目录；不消费 --start。

electricity.dict 包含两个块。物理 Electricity 块复用普通 EHD 的严格 reader、默认值与数值服务，例如：

```text
Electricity
{
    model leakyDielectric
    interfaceDiscretization sharpHermiteP3
    permittivityNegative 2
    permittivityPositive 1
    electrodeAxis y
    lowerPotential 0
    upperPotential 0.01
    linearResidualTolerance 1e-9
    conductivityNegative 0.1
    conductivityPositive 0.01
    initialSurfaceChargeDensity 0.001
    chargeConvection false
}
Domain
{
    coordinateSystem axisymmetricRZ
    interfaceLevelSet liquid
    interfaceHalfWidth 0.08
    motion fixed
}
```

Domain 由 App 绑定器消费：coordinateSystem 选 cartesianXY/axisymmetricRZ，interfaceLevelSet 必须选中 geometry.dict 唯一的静止水平集，interfaceHalfWidth 为正的正则化半宽（代码长度），motion 当前只支持 fixed。这里没有另一套材料或几何状态。物理块选 perfectDielectric 时不填电导率/初始电荷/chargeConvection；XY 支持原 diffuseTanh 或 sharpHermiteP3，RZ/导电要求 sharpHermiteP3。固定界面没有速度提供者，chargeConvection 必须 false。电极轴选 x/y；RZ 域始于 r=0 时仅允许轴向 y 电极，其余边界齐次 Neumann。数值求解始终 SuperLU；当前直接消费 linearResidualTolerance，不接受无消费者的 solvers.dict。

geometry/mesh 沿用 Core 多行格式：完整矩形域、单一固定水平集、无 Phase/Boundary/周期/固体运动，几何 dimension 为二维坐标网格的 2；RZ 度量由 coordinateSystem 显式指定。mesh 只允许 adaptiveMesh false 与 BackgroundMesh/divide，Nx/Ny >= 8；不接受 AMR/细化设置。正面积外边界必须满足原 sharp 几何分离和模板检查，r=0 的平滑轴接触由原对称路径处理。预检查不建网格，实际几何质量仍在运行时验证。

每个导电步先消费旧场真实电流更新电荷，再求当前电荷的电场，成功后发布同一时刻的电势/电荷。保留原显式 CFL <= 0.25 和电极闭合门限，不支持该独立入口的隐式/材料运动策略。固定网格场输出包含原单元电势、epsilon、Ex/Ey、节点水平集；interface_panel_mean_charge 为界面单元内的面元平均面电荷密度，其余单元为零，不是体电荷。series/electric_electric.dat（默认样例 prefix=electric）记录接受步、真实 SuperLU 残差、总电荷、电流来源时刻与电极收支。meta 保留完整输入及生效数值选项。

验证范围见[独立入口记录](validation/organization/QUANTUM_Modular_M4_Electricity_App_20260916.md)：解析均匀介质场、显式总电荷衰减、同一 App 独立链接及旧耦合输出对照。该结果不宣称旧顶点初始化、任意多界面、SI 材料表、规定运动或 AMR 已迁移。

### 7.1 普通电学的独立守恒参数（2026-09-19）

正式 Electricity App 和第 15 节分文件 EHD 共用 `electricity.dict/Electricity/Conservation`。`linearResidualTolerance` 只控制线性求解及其代数残差；省略下列物理参数时采用独立默认值，不随线性精度变化：

| 键 | 默认值 | 实际消费者及支持范围 |
| --- | --- | --- |
| `boundaryFlux` | `1e-9` | diffuseTanh 的全域边界位移通量收支；不用于 sharpHermiteP3 |
| `reconstructedPoisson` | `1e-9` | sharpHermiteP3 重构通量的逐单元方程残差 |
| `gaussClosure` | `1e-9` | sharpHermiteP3 的全域 Gauss 收支，含已给定面电荷 |
| `chargeClosure` | `1.6e-8` | leakyDielectric 的未分配体电流及电荷增量—电极电流闭合；该值直接使用，不再乘 16 |

四项均须有限且位于 `(0,1)`；未知、重复、非法值及没有对应模型/离散消费者的显式设置拒绝。导电锐界面示例为 `Conservation { reconstructedPoisson 1e-9; gaussClosure 1e-9; chargeClosure 1.6e-8; }`，diffuseTanh 示例只填写 `Conservation { boundaryFlux 1e-9; }`。原 Gauss/边界通量的网格规模舍入项、面电荷输运的机器精度检查及显式 CFL 限制保持。

旧扁平 EHD 输入不接受 Conservation，兼容 reader 将原 `linearResidualTolerance` 展开为实际使用的通量容限，并将 16 倍值交给电荷收支消费者。正式输入若以前填写了不同于 `1e-9` 的线性容限，迁移时需显式填写原通量容限及原 16 倍电荷容限（仅填写适用项），才能保留原物理验收门。直接模块调用也使用独立默认值；拟合喷嘴继续消费第 5 节既有 Numerics/SurfaceCharge 参数，拒绝普通域的 boundaryFlux/chargeClosure 选项。生效配置仍写入原 `output/meta/electricity_effective.txt`，列出实际通量和电荷容限，不增加报告文件。

## 8. PIC 场输出选择（2026-09-17）

ExplicitPIC / ImplicitPIC 沿用 `system/pic.dict` 的原格式。`DictionaryPIC` 读取字段列表，由 [PicFieldSelection](../include/Modules/Plasma/IO/PicFieldSelection.h) 解释；不是喷嘴 schema 2 的 `run.dict/Output.Fields`。例如：

```text
Output
{
    Fields
    {
        every 1
        list (phi rho n E Te_e electron_mean_energy)
    }
}
```

字段名不区分大小写；列表决定准备哪些已有诊断场，不改变求解模型或单位归一化。

| 选择名 | 内容与范围 |
| --- | --- |
| `phi`、`rho` | 电势、体电荷密度；两种 PIC 均支持 |
| `E` / `efield` / `electricfield` / `electric`，或 `Ex`、`Ey` | 两个电场分量，或单独选分量 |
| `n` | 各运行物种的 SI 数密度；显式无论列表是否含 n 都写出 |
| `Te_e` / `te` / `electron_temperature` / `electron_temp` / `temperature` | 电子温度，列名 `Te_e_eV` |
| `electron_mean_energy` / `electron_mean_energy_ev` / `emean` / `mean_energy` | 电子平均动能，列名 `electron_mean_energy_eV` |
| `psi` / `levelset` / `level_set` / `ls` / `level-set` | 原水平集显示场；隐式采用组合水平集 |
| `psi_composite` 及 `levelset_composite` / `level_set_composite` / `ls_composite` / `composite_levelset` / `composite_level_set` | 仅隐式，组合水平集 |
| `psi_all` 及 `levelset_all` / `level_set_all` / `ls_all` / `all_levelsets` / `all_level_sets` | 仅隐式，逐个输出已绑定的嵌入水平集 |
| `muCell`、`muFace` | 仅隐式，现有单元及面代理极化系数诊断 |
| `w_mean`（`w` / `wmean` / `weight` / `weights`） | 仅显式，单元平均粒子权重 |
| `w_sum`（`wsum` / `weight_sum` / `weightsum` / `wphys`）、`w_max`（`wmax` / `weight_max`） | 仅显式，单元权重总和、最大值 |
| `n_macro`（`macro` / `macrocount` / `particle` / `particlecount` / `particle_count` / `n_particles` / `np`） | 仅显式，单元宏粒子数 |

显式空列表使用 `phi rho psi E n Te_e n_macro`，保持原固定输出顺序；隐式空列表不写场。隐式按列表中最早出现的对应名称/别名排序，同一物种组或电场分量保持原稳定次序。未识别名称不启用字段，显式仍输出物种密度；隐式只有未知名称时不生成场文件。逐步收支仍由其独立配置控制。

物种与水平集名称转为小写字段键并合并/裁掉多余下划线。显式将 `+ - *` 分别写为 `plus minus star`，隐式将这些符号按下划线处理，例如 `Ar+` / `Ar-` 的密度列分别为显式 `n_arplus_SI_1pm3` / `n_arminus_SI_1pm3`、隐式 `n_ar_SI_1pm3` / `n_ar_2_SI_1pm3`。同一基础名重复时沿用 `_2` 等后缀；计数在每个物种字段组重新开始，仍不提供与独立带后缀基础名的全局去重。

字段是否实际产生仍取决于原场景条件，例如电子能量需要存在有粒子的电子物种，逐水平集输出需要已绑定的嵌入几何。本节不扩展这些能力；实际对照见[字段选择记录](validation/organization/QUANTUM_Modular_M6_Pic_Field_Selection_20260917.md)。

## 9. PIC 瞬时与平均剖面（2026-09-17）

两种 PIC 沿用 `system/pic.dict` 的 `Output` 格式，由 `DictionaryPIC` 解析，运行中的文件与窗口归 [PicProfiles](../include/Modules/Plasma/IO/PicProfiles.h)。示例：

```text
Output
{
    Fields
    {
        every 2
        list (phi rho n E)
    }
    Profile1D
    {
        enabled true
        nbinsX 4
    }
    Profile1DAvg
    {
        enabled true
        avgWindowSteps 2
    }
}
```

- 瞬时剖面输出到 `output/fields/profile1d.csv`，每个输出步、每个 x 分箱一行；默认关闭，`nbinsX` 的 PIC 字典默认值为 256。
- 平均剖面输出到 `output/fields/profiles_avg.csv`，窗口就绪后覆盖单份快照。默认关闭、窗口 0；只有启用且窗口大于 0 才累积。统计采用滑动窗口的均值和总体标准差。
- 两者的写出频率均来自 `Fields.every`，驱动按最小 1 归一化；平均在每次有效采样时累积，写出还需满足频率和窗口就绪。剖面块内的 `every` / `path` 不作为当前 reader 的有效独立设置，文件名按上述固定。
- `Profile1D` 和 `Profile1DAvg` 也可嵌在 `Output/Fields` 下；`ProfilesAvg` / `profiles_avg` 是平均块的别名，`enable` 是 `enabled` 的别名。`nbins` 是 `nbinsX` 的别名；`window` / `windowSteps` 是 `avgWindowSteps` 的别名；读取这些数值时，如尚未明确指定 enabled，会自动启用对应功能。
- 保留原非正分箱数行为：瞬时输出至少 1 箱，平均采用 128 箱回退；窗口小于等于 0 禁用平均。
- 恢复时瞬时 CSV 按传入的输出恢复边界删除该步及之后行，保留更早数据与表头。平均窗口重新从空开始，未就绪期间原平均快照保持；窗口没有随检查点恢复。

显式瞬时剖面位于场输出窗口，使用重建的诊断缓存；平均使用求解阶段缓存。隐式两类剖面都在原求解阶段读取缓存及 EcellNew。采样点及原分箱、密度/能量计算规则保持；这一配置不改变 PIC 时间积分。[验证记录](validation/organization/QUANTUM_Modular_M6_Pic_Profiles_20260917.md)覆盖每两步输出、两步平均窗口和旧检查点恢复。

## 10. PIC 周期诊断与电流波形（2026-09-17）

```text
Output
{
    Series
    {
        every 2
    }
    PeriodicDiagnostics
    {
        enabled true
    }
    CurrentWaveform
    {
        enabled true
        side left
    }
}
```

- 两项默认关闭，写出频率均取 `Output/Series.every`（最小 1）；固定文件分别为 `output/series/periodic_step.csv` 和 `output/series/waveform_step.csv`。各功能块内的 `every` / `path` 不作为当前 reader 的独立设置。
- `PeriodicDiagnostics` 也可嵌在 `Output/Console`；`console` / `printConsole` 控制显式摘要，尚未显式设置 enabled 时会自动启用诊断。隐式配对详情由环境变量 `QUANTUM_IMPLICIT_PERIODIC_DIAG_DEBUG` 开启。
- 周期诊断开关本身不设置物理边界。对应两侧 Field BC 必须为 Periodic，且通过边界 ID 映射到 geometry.dict 域边界；`PeriodicPairs` 必须定义这两个 ID 的配对及平移。只在 PIC Domain 写 Periodic 不足以启用完整周期闭合，应检查 CSV 的 periodicX/periodicY 标志。
- `CurrentWaveform` / `Waveform` 也可嵌在 `Output/Series`；`side` / `electrode` 指定 left/right/bottom/top/auto。auto 优先时间驱动的 Dirichlet 侧，否则首个 Dirichlet 侧；默认 auto，未明确设置 enabled 时指定 side 会自动启用。
- 周期文件保留六列 step/time_s/phiJumpX_V/phiJumpY_V/periodicX/periodicY；波形保留 step/t_si_s/V_drive_V/J_cond_Am2/J_disp_Am2/J_total_Am2。开启诊断不更改原求解/电流差分时层；续算按传入输出恢复边界裁剪旧 CSV。

完整可运行输入及前后对照见[边界诊断记录](validation/organization/QUANTUM_Modular_M6_Pic_Boundary_Diagnostics_20260917.md)。本节记录当前 reader 和消费者行为，未引入新语法。

波形历史规则：电极选择在初始化时固定；位移电流每个原统计步根据相邻电场样本计算，`Series.every` 只控制写出。当前检查点未保存这份诊断历史，也不从 CSV 重建，因此恢复后的第一个采样步位移电流为零。关闭波形文件不等价于取消其他诊断中的选定侧电流；无有效侧才标记为不可用。波形传导项沿用墙面物种出流，四侧驱动功的注入修正口径另行保持，见[验证记录](validation/organization/QUANTUM_Modular_M6_Pic_Current_Waveform_20260917.md)。

边界电功沿用 `power_step.csv`：显式需要 `system/units.dict` 的 `Units { mode CodeUnits writeSI true }`（或既有 SI 配置）及已建立的单位/归一化上下文；隐式按原上下文条件启用。域侧显式扣除源/二次电子注入，隐式使用粒子出流，均含位移电流及嵌入 Dirichlet 边界功。零/非正 dt 不推进这份四侧历史，恢复首步无前一电场参考；它与波形历史分别管理。未增加字典项，见[验证记录](validation/organization/QUANTUM_Modular_M6_Pic_Boundary_Work_20260917.md)。

PIC Series 统计口径：显式物种预算、边界/过程/几何诊断累计到下一次 `Series.every` 输出（含末步不足窗口），报告后下一步清空；隐式沿用采样步本身的数据。`Fields.every` 独立控制场输出，不能用其频率推断预算窗口。恢复从空窗口重新累计，原检查点未保存部分窗口，CSV 裁剪边界仍各自保持。没有新增字典项，见[验证记录](validation/organization/QUANTUM_Modular_M6_Pic_Diagnostic_Window_20260917.md)。

## 11. PIC 检查点失败规则（2026-09-17）

2026-09-18 网格系数缓存修复没有新增字典项：新建 RHS 初始化为零，QMC2/BMC2 写出时将结构体填充字节置零；原版本、布局及旧成员值保持。旧缓存中零计数条目的非零 RHS 也会保留，不要求删除缓存重新生成。双向读取和真实续算范围见[系数缓存记录](validation/organization/QUANTUM_Modular_M6_Mesh_Cache_Coefficients_20260918.md)。

`Checkpoint` 原键保持：enabled 默认 false、every 默认 100、keep 默认 2、path 默认 `../cache`、saveRng 默认 true、restart 默认 false、restartStep 默认 -1。启用且 every 大于 0 时按周期和末步请求保存；keep 小于等于 0 保留全部。相对 path 以运行目录为基准，目录名为 `pic.<下一步号>`；负 restartStep 自动选择最新可用项，非负值指定目录。

新程序在写入前建立 `.pic_checkpoint_incomplete`，全部必需文件成功写完并关闭后移除。带标记目录不能自动或显式选择；自动选择还要求必需组件是非空普通文件。keep 仅计算可用目录，失败目录保留用于诊断。写入失败结束运行，显式返回 2，隐式沿已有完成回调返回 -3；已执行物理步、普通输出及同一步覆盖不回滚，其他旧检查点不因这次失败清理。

几何要求按实际构造条件判断，仅 Domain 的 geometry.dict 不要求 iGeometry 缓存；MovingEB 开启且保存恢复状态时，默认 `moving_eb_cache.bin` 和相对 restartPath 都基于选定检查点目录，绝对路径按原设置使用。成功格式不变，没有新增输入键，见[验证记录](validation/organization/QUANTUM_Modular_M6_Pic_Checkpoint_Completion_20260917.md)。旧格式显式选择/宽松 reader 规则保留；非空检测不等于任意历史损坏验证，也不提供掉电持久性或旧程序的标记识别能力。

PIC 功率表 `output/series/power_step.csv` 沿用既有语法。显式按原场输出频率报告，源汇能量累计到报告后才清零；隐式历史每个诊断步更新，Series 频率只决定写行，能量变化率始终使用当前物理 dt。显式采用相邻报告时间差；恢复后两种历史均为空，不从保留 CSV 重建。显式 CodeUnits 列名和驱动功置零规则保留。上述两种时间口径与 Series 预算窗口各自独立，没有新增输入键，见[能量收支记录](validation/organization/QUANTUM_Modular_M6_Pic_Energy_Balance_20260917.md)。

`Output/CurrentBalance` 的启用语法和 `Series.every` 沿用。显式只写整除步，默认浮点格式；隐式额外写最后一步，电荷列固定 16 位小数。微小电荷可能舍入为零，旧归一化和恢复裁剪边界保持；本批未改变输出精度或新增配置。见[电荷收支记录](validation/organization/QUANTUM_Modular_M6_Pic_Charge_Balance_20260917.md)。

PIC `BC/Field` 的边界 ID 映射已由显隐共同服务处理，语法/默认/容差不变。域边优先由线端点位置判断，无法分类才用 left/right/bottom/top 名称；跨两个域边的同一实体仍拒绝。显式周期 wrap 仍同时要求 Field 两侧 Periodic 与匹配几何配对，隐式沿用原 Field 判定。`BC/Particle` 的全局与局部 SEE 参数选择仍保留两条路径的原差异；无新字典版本。见[边界初始化记录](validation/organization/QUANTUM_Modular_M6_Pic_Boundary_Setup_20260917.md)。

PIC 的嵌入边界/水平集选择和内部组合缓存改由模块对象管理，字典语法、默认符号和退化处理未变；显式无有效水平集时仍关闭 EB，隐式保持原 ready/启用区分。既有 alpha 采样数 6、显隐包含容差和移动壁面策略均保持，无新增配置。见[嵌入几何记录](validation/organization/QUANTUM_Modular_M6_Pic_Embedded_Geometry_20260917.md)。

PIC 的几何启动由模块共享实现，GeometryRefinement、maxInitialIterations、MotionSchedule、初始化快照环境开关及检查点语法不变。启动细化仍屏蔽物理指标，只在目标层级有变化时调用 adapt，重建原几何后再次判断收敛；静态几何的运行期细化关闭规则保持。无新增配置。见[几何启动记录](validation/organization/QUANTUM_Modular_M6_Pic_Geometry_Startup_20260917.md)。

PIC GeometryRefresh 及既有 AMR 几何环境开关的解析、优先级、数值回退与告警保持；已解析值映射到 `PicAmrGeometry` 的类型化策略。full/adaptive/skip、移动壁面强制刷新、活跃/混合选择及连续次数上限仍按原语义执行，无新增字典键。整数环境读取的既有最小值 1 规则未改。见[AMR 几何记录](validation/organization/QUANTUM_Modular_M6_Pic_Amr_Geometry_20260917.md)。

MotionSchedule 与 `QUANTUM_IMPLICIT_MOVING_TRACE_STEP/ALL/PATH` 的解析、频率、阶段切换及追加写入语义保持；仅将运行状态和几何更新流程交给 `PicScheduledGeometry`。无新增字典键、检查点字段或配置优先级。见[规定运动记录](validation/organization/QUANTUM_Modular_M6_Pic_Scheduled_Geometry_20260917.md)。

MovingEB 原字典、单位、更新/输出频率和检查点格式保持，配置交给 `PicMovingGeometry` 的控制器持有。restart=false 仍跳过缓存读写；启用恢复但缓存无效仍中止。单水平集侵蚀的切割标签生命周期已修复，显式时变边界残差改用求解器的当前边界值，无新增输入键。见[MovingEB 记录](validation/organization/QUANTUM_Modular_M6_Pic_Moving_Geometry_20260917.md)。

显式 PIC `overview_detail_step.csv` 续算现在补齐原表头中的物种密度/宏粒子数/物理数列，输出频率、SI/代码单位和精度不变，无新增输入键。同一运行保持路径和单位模式；恢复首行的生成表头必须与已有表头相同，不兼容时明确报错并保留该文件。已有历史行按原恢复步截断规则保留，不追补过去已漏写的值。见[detail CSV 记录](validation/organization/QUANTUM_Modular_M6_Pic_Detail_Series_20260917.md)。

启用 MCC 且 `Checkpoint.saveRng true` 时，检查点增加 `pic_collisions_state.txt`，记录随机数模式、种子及实际碰撞调用次数；没有新增字典键。`counterBased` 恢复必须有该组件，旧检查点不能仅靠步号补算，缺失/损坏/模式或种子不匹配会拒绝恢复；显式预检在裁剪输出历史前执行。旧 `global`/`perParticle` 检查点允许缺少该文件。`saveRng false` 同时跳过这项读写，不能据此声称确定性续算。极化求解的 `freezeDirectMatrixUntilPatternChange` 原键保留，实际分解仅在矩阵结构和系数均不变时复用；系数变化会重新分解。见[验证记录](validation/organization/QUANTUM_Modular_M6_Pic_Collision_Restart_20260917.md)。

## 12. PIC 电极、反应与能谱诊断（2026-09-18）

以下为 `system/pic.dict` 已有语法，未增加输入项：

```text
Output
{
    Diagnostics
    {
        every 2
        baseDir physics_diagnostics
        ElectrodeCurrent
        {
            enabled true
            includeDomain true
            perBoundary true
            group electrode (1 100 101)
        }
        ReactionDiagnostics
        {
            enabled true
            summary true
            perChannel true
            perProjectile true
            spawnedSpecies true
        }
        EnergyHistogram
        {
            enabled true
            bins 80
            maxEnergyEV 200
            sample electron_bulk electron bulk
            sample wall_all all wall
        }
    }
}
```

- 这三项默认关闭；功能块内的正 `every` 优先，否则取 Diagnostics 的正 every，再回退 Series 输出频率。按整除步采样/写出，不累加被跳过的诊断步，也不强制补写末步。示例文件位于 `output/series/physics_diagnostics/` 下的 `electrodes`、`reactions`、`energy` 子目录；baseDir 默认 `diagnostics`。
- 电极正边界 ID 重复时计入首个配置分组，未分组边界仍可逐边界输出；includeDomain 控制是否计入四侧域边界。原分组表 boundary_count 按配置列表长度计数，不据去重后的边界数重写。
- 非零反应统计需要已建立的 MCC 配置，诊断开关本身不启用碰撞。显式现在从第一个满足频率的步写行，已有后续行口径保持；恢复仍按原输出恢复边界裁剪，不追补过去漏写的记录。`ReactionDiagnostics` 和 `Reactions` 两个既有诊断块名均可使用；碰撞输入发现现在只选择顶层，嵌套报表不再导致缺少截面的误报。
- 能谱 sample 依次为名称、物种选择、位置。位置 bulk 为体内样本，wall 包含 domain/embedded；物种可用实际名称或 all/electron/ion/neutral，后三种按电荷符号分类。没有 sample 时每个配置物种生成一个 bulk 样本；达到/超过能量上限的样本计入末箱，每次采样前清零。
- 输出流只有实际打开且状态正常才参与采样。开启、关闭或改变这些可选报表的频率不再改变嵌入边界电功；相应修正会改变旧程序遗漏电功的 power/process 行，场和粒子推进保持。

新旧实际 MCC、每步/隔步、电功开关及旧检查点的验证范围见[物理诊断记录](validation/organization/QUANTUM_Modular_M6_Pic_Physics_Diagnostics_20260918.md)。

## 13. PIC 内嵌输入的块层级（2026-09-18）

物理碰撞输入逐块按以下优先级选择：mcc.dict 的顶层块、pic.dict 的顶层块、对应的独立 system 字典。Gas、CrossSections、Reactions、Collisions 可分别来自不同来源；既有别名和同一文件内的首个匹配规则保持。

`Output/Diagnostics/Reactions` 是报表设置，不是物理 Reactions。扩展块中嵌套的 Gas/CrossSections/Collisions 同样不参与碰撞装配。仅启用反应报表无需补建一套碰撞输入。

Solvers、Units、Normalization 内嵌于 pic.dict 时也必须是顶层块；同时存在顶层块与同名独立字典仍报重复来源，嵌套扩展中的同名项不会触发该检查。独立字典和内嵌块的实际物理参数解释保持，文件开头的 UTF-8 BOM 可继续使用。

块名采用现有未加引号的关键字语法，名字后只能有空白/行注释再接 `{`。提取器不会把标量值、列表成员或引号字符串中的同名文字提升为块。这里没有新增物理键，也未将所有旧 reader 改为统一严格解析。

[验证记录](validation/organization/QUANTUM_Modular_M6_Dictionary_Block_Scope_20260918.md)包含旧误报、真实重复来源、两种诊断别名、混合物理输入、旧失败缓存和 BOM 的实际运行对照。

## 隐式 PIC 环境覆盖的读取时机（2026-09-18）

旧 `QUANTUM_IMPLICIT_*` 选项以及 `QUANTUM_SEGMENTER_TRACE/NOHINT` 的本批 65 项名称、优先级、解析及默认值保持，没有新增字典键。它们现在于 `readPicCaseInput` 准备输入时一次读取，运行中修改环境不再改变该次隐式 PIC；下一次运行重新读取。直接使用模块库时默认不读取进程环境，可由调用者提供 `ImplicitPicOverrides`，或显式调用 IO 的 `readImplicitPicOverrides()`。未提供的覆盖项继续采用原字典/默认策略，相对输出/快照路径仍按执行目录解释。

TRACE/NOHINT 仍按变量是否存在启用，`0` 和 `false` 也算启用；要关闭应取消该环境变量。NOHINT 让初始及后续宿主单元都按位置查询，不使用粒子或上一段的提示编号；原跨单元分段超限问题已修复，真实分段预算限制仍保留，见[轨道定位记录](validation/organization/QUANTUM_Modular_M6_Pic_Segmenter_Lookup_20260918.md)。选项清单及输入捕获规则见[环境输入记录](validation/organization/QUANTUM_Modular_M6_Pic_Environment_Inputs_20260918.md)。


## 14. Core 网格独立验证（2026-09-19）

`Run { app MeshVerification; }` 选择网格验证 App，无版本号；不加载 Plasma、单位归一化或物理求解设置。主程序与 `quantum_mesh_standalone` 消费同一 App，后者只链接 Core。标准输入和运行命令见 [mesh_refinement](../cases/mesh_refinement/README.md)。

| 文件/块 | 实际支持与消费者 |
| --- | --- |
| `run.dict/Run` | 仅 `app MeshVerification`；不接受 Time、Restart、Output、schemaVersion 或旧模块配置 |
| `geometry.dict/Domain` | 仅 `min (x,y)`、`max (x,y)`，有限代码坐标及正面积；域由既有 DictionaryMesh reader 读取 |
| `mesh.dict/BackgroundMesh` | 仅 `divide (Nx,Ny)`，两个整数各在 2–64 范围 |
| `mesh.dict/adaptiveMesh` | 必须显式 `true` |
| `mesh.dict/QuadtreeMesh` | 必须显式 `maxCellNum`、`maxInitialIterations 0`，验证 App 自行安排适配阶段 |
| `mesh.dict/GeometryRefinement` | `enabled true`，一个 `refineRegion/entry`；必填 name、shape box、xmin/xmax/ymin/ymax、maxLevel；区域具有正面积并位于域内，级别 1–3 |

沿用原 mesh reader 的逐行块/记录语法；实际读取值须与严格能力过滤一致。最细完整网格 `Nx*Ny*4^maxLevel <= 4096`，容量须满足 `2*(Nx*2^maxLevel+1)*(Ny*2^maxLevel+1) <= maxCellNum <= 65536`。这里只提供一个有界验证流程，不声称支持所有 Core mesh reader 选项。多区域、边界/相细化、物理指标、运动、几何对象、其他字典文件和未列出的字段均拒绝，不能静默忽略。

固定执行初始、区域细化、QMC2 网格读回后粗化、再次细化四阶段；调用现有 refiner/adapt/field-transfer/face-graph/cache 服务。`--check-case` 无网格及文件写入；实际运行只接受不存在的 output/cache，保留旧结果。输出按 fields/series/meta 分工，不将阶段号当物理时间。缓存只含网格，解析探针场继续使用内存副本；不支持物理检查点续算。解析误差门与结果清单见算例指南，失败保留已写出的证据。

## 15. 普通 Hydro/EHD 分文件输入（2026-09-19）

适用于固定均匀网格的 `SingleModule { modules (Hydrodynamics); }` 和 `Run { app Electrohydrodynamics; }`。保持原 App 选择、代码单位、物性/边界和求解器能力；不要添加 `schemaVersion 2`。这是普通节点界面格式，不是第 6 节规定速度 Interface App 的格式。

在 `run.dict` 中出现 Time、Output、Log 或 Restart 时，整套采用分文件输入：

| 文件 | 拥有者与实际用途 |
| --- | --- |
| `run.dict` | App/Runtime 的 Time、Output、Log；主时钟和共享输出选择 |
| `hydro.dict` | Hydro 物性、初速度、压力/速度边界、界面载荷、线性残差容限及独立 Conservation；Hydrodynamics 块 |
| `interface.dict` | Interface 的 `advectInterface`、`reinitializeEvery`、`interfaceBoundary`；twoPhase 必须提供，constant 不接受 |
| `geometry.dict`、`mesh.dict` | 保持原域、水平集和均匀网格格式 |
| `electricity.dict` | 仅 EHD 需要，保持原电学物性、离散与电荷设置 |

`run.dict` 示例；独立 Hydro 将首行替换为 `Run { app SingleModule; }`，再加 `SingleModule { modules (Hydrodynamics); }`：

```text
Run { app Electrohydrodynamics; }
Time { units code; dt 0.001; endTime 0.003; }
Output
{
    prefix hydro;
    Fields
    {
        enabled true; format tecplot; profile physical;
        writeControl stepInterval; writeInterval 3;
        writeInitial true; writeFirst false; writeFinal true;
    }
}
Log { everySteps 1; level summary; }
Restart { model none; }
```

Time/Output/Log 必填；Restart 可省略，省略等价于 none。`endTime` 须为整数步，步数为 1..INT_MAX-1；`Time.units` 仅 code。输出支持 `stepInterval`（正整数步）和 `timeInterval`（正代码时间），复用公共运行设置。只支持 `format tecplot`、`profile physical` 和 `Restart.model none`。日志默认 `summary`；仅独立保守 Hydro（`Interface.transport conservative`、未耦合电学）另支持 `Log.level detail`，输出每步材料、界面、毛细力、质量通量、NS、接受及写场耗时，界面再分供体、通量、库存、预测、同步、材料映射。计时不参与数值决策；普通旧方法和 EHD 仍只支持 summary。`enabled false` 只关闭场文件，逐步收支仍保留。首帧、第一步和末步的强制输出分别由三个布尔项控制；常规间隔命中的帧仍会写出。EHD 初始电场在第 1 步求值，`writeInitial` 也会选中其第 1 步电学输出；电场/载荷和完成流场仍保留各自原时间层，不改写为同一时刻。

上述整数步及 `step*dt` 规则是省略 `Time.StepRetry` 时的默认路径。保守 Hydro 及已接入联合接受的 `ConservativeColumn` 带电液柱可显式使用（闭域 Hydro AMR 的候选接受见 15.2）：

```text
Time
{
    units code;
    dt 0.001;
    endTime 0.0033;
    StepRetry { minDt 0.000125; maxRetries 3; }
}
```

该可选块必须同时提供正 `minDt <= dt` 和整数 `maxRetries`（1–30）；`endTime/minDt` 必须小于 `INT_MAX-1`。`Interface.transport` 须为 `conservative`，电学耦合须提供联合候选接受协议。普通节点 Hydro、普通漫界面 EHD、Electrospray 及独立 Interface/Electricity 不支持此块。旧格式仍只支持原固定步长。

启用后 `dt` 是初始最大步长；物理候选失败且已恢复原接受态时，将本步 dt 减半，重新准备界面、相/质量通量、毛细力和 NS。每个接受步至多减半 `maxRetries` 次，减半后低于 `minDt` 则停止。默认成功后沿用实际接受的 dt。可选 `recoverAfter N;`（整数 0–1000000，默认 0 关闭）：连续 N 个无需重试的完整接受步后，下一步候选 dt 提高至原来的 1.25 倍，最多恢复到 `Time.dt`；发生重试时稳定步计数清零，物理接受后的输出失败仍停止。新候选照常执行原 CFL、精度、守恒及回滚门，不是误差放宽。末步截短到精确 `endTime`，允许这一个剩余区间短于 `minDt`，但不能再把失败的末步缩到该下限以下。空间欠分辨、拓扑闭合等错误可能不能靠减步解决，预算耗尽后明确停止。近断裂无电/带电模板现取 `recoverAfter 8;`；启动脚本的 `-RecoverAfter 0` 关闭恢复并省略该键，便于保留旧冻结程序和时间对照的原策略。

物理步号仅在共同接受后增加。`hydro_diagnostics.csv` 写实际时刻和实际 dt，Interface 历程与场文件绑定相同接受时层。`timeInterval` 按相邻接受时刻是否跨过从零起算的输出间隔选帧，`stepInterval` 按接受步号；末帧仍受 `writeFinal` 控制。被放弃的尝试只写 `[Hydrodynamics retry]` 日志，不冒充一个接受步；最终未能恢复推进时仍保留失败诊断。物理接受后的输出错误停止运行，不重复该步。

保守动量显式选择 `HydroTimePolicy::ConservativeOneStep`：每步使用 `M^{n+1}, M^n` 及该区间唯一积分质量通量，隐式迎风对流替代旧 AB2，黏性梯形项只依赖当前区间。完整耦合压力修正使用该步实际矩阵。它仍是冻结输运速度的总体一阶方案，不能当作变步长二阶 AB2；没有更改普通/拟合 Hydro 的固定步长保护。闭域 Hydro AMR 使用未发布的独立候选（15.2）；已接入的终端截面脱离及 `ConservativeColumn` 电荷使用各自联合接受路径，电流体共同 AMR 尚未接通。

保守 Hydro 在相同写场时刻输出 `<prefix>_interface_geometry_stepNNNN.csv`：记录已接受原生曲面的表面 ID、图段及端点/中点位置、法向和分配的物理权重。位置直接来自唯一已接受 `Γh`，不由显示用节点 phi 作线性等值线替代。每段端点/中点采用 1/6、2/3、1/6 的采样权，三项之和为原段物理测度；用于形状诊断的离散求积，不能据此声称连续最大误差上界。

detail 的界面计时行另含 `constant_phase_chords`（整个供体带严格同一纯相的弦数）、`primitive_requests` 和 `primitive_evaluations`（区间原函数请求数和局部复用后的实际计算数）。这些仅为成本观测，不改变库存/几何门限，也不随物理状态保存；次数减少不直接等价于同幅度的总耗时下降。

最小 `interface.dict` 为 `Interface { advectInterface true; reinitializeEvery 0; }`。`advectInterface` 必填；重初始化周期默认 0，非零要求实际运动。`interfaceBoundary` 可省略；提供时沿用 prescribedScalar 的 model/fixedContacts/inflow 块，命名入口仍必须绑定真实负法向速度流入边界。接触入口仅 Hydro 支持，EHD 不因此获得壁面电荷输运能力。

分文件 `hydro.dict/Hydrodynamics` 中，`linearResidualTolerance` 只控制线性求解，默认 `1e-8`，允许范围仍为 `(0, 0.01]`。可选 `Conservation` 块单独控制物理检查：

| 键 | 默认值 | 实际消费者及原判据 |
| --- | --- | --- |
| `momentumBalance` | `1e-8` | 联合修正后的完整动量相对残差；不代替预测/压力线性残差 |
| `initialContinuity` | `1e-8` | 初始边界相容投影后的全单元体积流率范数，保留原舍入项；要求 boundaryPatches |
| `outletAdmissibility` | `1e-8` | 初始化及时间步的出口回流检查，保留原速度/压力冲量尺度；要求 pressureOutlet |
| `boundaryFluxBalance` | `1e-8` | 每步全边界体积流量及常物性质量流量闭合，保留原 `16 × 容限 × 流量尺度`；要求 boundaryPatches |

四项均须有限且在 `(0,1)` 内。省略项使用各自默认值，不随 `linearResidualTolerance` 改变。闭域示例为 `Conservation { momentumBalance 1e-8; }`；边界算例可按上表补全其他项。未知/重复键、非法值和没有对应边界的显式设置拒绝。普通域与喷雾拟合域的 Conservation 格式及范数不同，不能混填。分文件 EHD 的电通量/电荷门限单独归 `electricity.dict`，见第 7.1 节，不借用 Hydro 容限。

旧扁平 Hydro/EHD 输入继续由兼容 reader 将 `linearResidualTolerance` 展开为原来的四项同值，保留原数值判据；旧格式不接受新 Conservation 块。迁移旧输入且原容限不是 `1e-8` 时，应在新 Conservation 中显式填写原值（只填写实际使用的边界项），才能保持原物理门限。直接模块调用的 `completeConservation` 也使用独立默认值；拟合域明确拒绝该选项，继续使用 fittedNumerics。原生非均匀完整域仍只报告原投影的完整动量缺陷，不因新增参数暗中增加联合修正或收敛门。

迁移时从 `hydro.dict` 移除 steps/dt/outputName/writeFields/writeInitialFields/writeEvery/logEvery，并将三个界面键移到 Interface 块。时间、输出、界面键在新旧位置重复时拒绝，不覆盖或忽略。没有上述运行块时继续按旧扁平 hydro.dict 读取，此时不接受 interface.dict。完整示例复用 [hydro_contact_smoke](../cases/hydro_contact_smoke/system/run.dict) 和 [ehd_diffuse_smoke](../cases/ehd_diffuse_smoke/system/run.dict)，无需新增标准算例。

### 15.0 CCH 中点单相流动（2026-09-27）

普通 `SingleModule/Hydrodynamics` 可显式选择 `timeScheme cchMidpoint;`，使用正式中点推进器。当前 case 范围为 `materialModel constant`、封闭均匀 XY/RZ、无旋流及固定步长；拒绝两相、EHD、AMR、`boundaryPatches` 和 `Time.StepRetry`。省略该键取 `cnAb2`，保留原运行路径。已有模块级参数两相中点入口不因此成为普通 case 的支持能力。

```text
Hydrodynamics {
    coordinateSystem axisymmetricRZ;
    timeScheme cchMidpoint;
    materialModel constant; density 1; viscosity 0.02;
    initialVelocity vortex; velocityAmplitude 0.2;
    linearResidualTolerance 1e-10;
    Midpoint {
        maximumIterations 30;
        relativeVelocityTolerance 1e-12;
        absoluteVelocityTolerance 1e-13;
    }
}
```

`Midpoint` 可省略，默认最多 16 次迭代、相对速度容限 `1e-9`、绝对速度容限 `1e-11`。次数限为 1–100，容限须正且有限；该块只能配合 `cchMidpoint`。这些控制检查同一中点速度的非线性一致性，不替代线性残差或物理守恒检查。迭代未收敛时整步拒绝并保留旧物理态；接受后的文件写出失败则停止，已接受的新物理态不回滚、不重算。

`initialVelocity` 支持 `zero`、原 `smooth` 和新增 `vortex`。`vortex` 从流函数差生成真实面平均速度，闭域初速满足离散散度约束；XY 为双正弦模态，RZ 为轴上正则的多项式/正弦涡流。RZ 初值及用法见 [hydro_midpoint_rz](../cases/hydro_midpoint_rz/README.md)。`smooth` 在中点路径启动时作相容投影。该单相输入不需要 `interface.dict`。

中点路径新增 `output/series/hydro_energy.csv`，包含实际质量、动能、累计黏性耗散、累计外力功及绝对收支差 `E+D−W−E0`，并记录中点迭代次数与速度缺陷。质量/动能使用原生物理体积和对偶面速度，不使用插值节点显示场。字段 `time` 为完整步时刻，`pressure_time` 为中点时刻。Tecplot 场同样带 `PressureTime`、`PressureTimeRole` 和 `VelocityTimeRole`；初帧压力仅为初始条件，不作为动态中点压力。原 `hydro_diagnostics.csv` 和场输出节奏保持。

### 15.0a 新矩空间的单相中点入口（2026-09-27）

**物理域内自由表面（2026-09-28，离域接续的部分代码）：** 无新增字典项，尚未开放普通两相 NS 的界面穿出口。正式表面能与 `interface_measure` 只积分原 Hermite 曲线位于实际 Domain 内的参数区间；虚拟相闭合边不计自由表面积。曲面被域边界截断时圆度/球形度为 NaN。常值压力分解的外边界项使用真实相区域从域内侧取得的湿润区间；弱毛细力保留截断端点项。[流体 §59](QUANTUM_CCH_Hydrodynamics_TODO.md#physical-surface-restriction)定义这一测度。

**出口运动与库存约束接续（2026-09-28）：** 完整 C1 场已从原 pressureOutlet 区间提供流函数法向二阶延拓，仅供材料轨迹使用；壁面及指定速度入口不允许越过。体积约束导数已限于物理域内，RZ 虚拟轴闭合可以随极点延伸到域外；固定供液口仍使用原边界定义。非固定闭合的出口相通量已接运动世界面。表面能穿出口输运与分量退出尚待接通，当前耦合求解器仍对中点和步末保留统一的离域拒绝；不能仅添加 pressureOutlet 就作为完整离域算例运行。没有新增用户参数或通过放宽字典绕过此限制。方法见[流体 §60](QUANTUM_CCH_Hydrodynamics_TODO.md#outlet-kinematics-inventory)。

**空界面与退场记录（2026-09-28，§61 接续）：** 模块已支持零活动界面的背景相状态及显式吸收性退场操作，普通耦合出口尚未启用，仍无新增字典参数。分量必须在物理域内既无自由曲面、又无负相区域才可退场；固定供液分量和微小域内液滴不会因此被删除。`hydro_moment_components.csv` 追加 `status/retired_at`：活动行为 `active/nan`；退场行为 `retired/接受步末时刻`，几何体积和弧数为零，`liquid_volume` 仍保留该分量的有符号 FV 数值余量。这些余量计入总库存及误差预算。没有活动界面时不产生该时刻的 `_interface*.dat` 线区，流场和相积分输出照常保留。方法、吸收性边界的含义与尚未接通的表面能出口项见[流体 §61](QUANTUM_CCH_Hydrodynamics_TODO.md#empty-phase-retirement)。

**表面能输运输出（2026-09-28，§62 接续）：** 无新增字典项。`hydro_moment_energy.csv` 追加 `cumulative_surface_advection`（向外为正的表面能输运）、`cumulative_surface_endpoint_work`（同轨迹的弱毛细端点功）、`cumulative_surface_flux_error`、`cumulative_surface_flux_budget`、`cumulative_surface_critical_estimate` 及本步 `surface_work_defect`。这些量已乘实际 σ；XY 为单位厚度的能量，RZ 为旋转后的物理表面能。`total_energy_balance` 增加前两项之差；`energy_balance` 仍是原动能方程余额。误差列只描述该求积的数值估计，不包含 H2 重构、时间推进和全部几何误差，也不构成任意 Flow 回调的严格误差证书。普通离域仍受保护，待实际出口 patch 准入与共同退场提交完成后开放。公式和临界交点处理见[流体 §62](QUANTUM_CCH_Hydrodynamics_TODO.md#surface-outlet-transport)。

**分相积分与能量观测（2026-09-28，完整参数两相路径，代码已接、待集中验收）：** 无新增字典项。`output/series/hydro_moment_phase_integrals.csv` 每个接受端点写两行，`phase=negative/positive` 分别对应字典负/正相，不能一概解读为液体/气体。字段含义如下：

| 字段 | 定义及范围 |
| --- | --- |
| `time`、`present` | 所属端点物理时间；几何体积大于零则 present=1。初态也有记录。 |
| `volume_fv`、`volume_geometry` | 独立有限体积库存与同一几何材料求积的体积；XY 为面积，RZ 为完整旋转体体积。 |
| `geometry_minus_fv`、`volume_fv_drift` | 几何体积减库存、库存减初值。开放域库存变化包含真实供排，不能直接当成守恒误差。 |
| XY `centroid_x/y`、RZ `mean_radius,centroid_z` | 位置的相体积平均；RZ 平均半径不是三维笛卡尔质心。 |
| `integral_u_*`、`mean_u_*` | 相区域内完整速度的体积分及除以几何体积后的平均值；不是界面节点平均、质心差分或最大速度。 |
| `central_xx,central_xy,central_yy` 或 `central_rr,central_rz,central_zz` | 物理坐标中心二阶矩的积分，未除以体积；RZ 使用 2πr 测度。 |
| `interface_measure` | 所有分量位于物理域内的自由界面长度（XY）或旋转表面积（RZ）；不包含闭合用壁面。两相行共享该值。 |
| `circularity`、`sphericity` | 负相封闭体的 2√(πA)/P（XY）或 π^(1/3)(6V)^(2/3)/S（RZ）；另一维度、正相、有固定接触端点或曲面被域边界截断时为 NaN。多分量给出整体指标。 |

相不存在时位置均值、速度均值和中心二阶矩为 NaN，不给零值制造有效观测。`hydro_moment_energy.csv` 中既有 `shape_moment` 改为物理中心矩之差除以和；RZ 横向矩为 ∫r²dV/2、轴向矩为 ∫(z−zc)²dV。旧输出使用旧定义，不能直接拼接新旧历程。

能量历史追加 `potential_energy`、`initial_potential_energy`、`mechanical_energy`、`mechanical_energy_change`。Eg=−∫ρg·x dV，使用同一几何材料，物理坐标原点为零势参考；机械能为 K+σS+Eg，change 为相对初态变化。它不等于计入黏性、供排、网格和事件后的一整步收支残差，也不替代原 `cumulative_external_work`；单相路径当前四列为 NaN。字典与压力中点时层保持原定义。详见[流体 §58](QUANTUM_CCH_Hydrodynamics_TODO.md#phase-observables)及 [Hysing 输出对照](../cases/hydro_hysing_bubble/reference/README.md)。

**参数曲率与分片（2026-09-28，代码更新，待集中数值验收）：** 无新增字典项。细化/粗化的主曲率指标读取同一四次曲线的数值驻点及叶面交点；特征线源与毛细积分共用相应几何分片。`maxCurvatureCells` 和 `minimumGapCells` 原物理含义不变，后者仍为法线间距观察。`*_interface_step*.dat` 新增 `curvature_meridional`、`curvature_azimuthal`，既有 `curvature_one_sided` 在 RZ 为二者之和；每段输出节点数可能改变。详见[界面 §2.45](QUANTUM_Conservative_AMR_Hermite_GALS_TODO.md#surface-extrema)。

**物理矩形域（2026-09-28，代码已接、未数值验收）：** 参数界面的网格、相容量、面通量、初始几何及毛细牵引现在共用实际 `Domain`。支持平面物理矩形和含轴 RZ 矩形，根网格两方向分割数可以不同；中心和半轴使用物理坐标/长度，不把域归一化为正方形。供液帽以实际南边界 z_min 为底，`radius/height` 为物理长度。均匀及混合原生叶消费相同的真实单元/面段；既有单位域输入无需修改。新输入见 [Hysing TC1](../cases/hydro_hysing_bubble/tc1/README.md)、[TC2](../cases/hydro_hysing_bubble/tc2/README.md)。积分观测量及作者参考曲线见本节更新；代码编译不代表矩形域或浮升数值验收。

**标准观测输出（2026-09-29，代码编译，未物理验收）：** 没有新增字典开关。普通 moment 场输出沿用原 Fields 频率，追加 `speed_sampled_max_5/9`（完整速度在每叶 5×5/9×9 节点的采样最大值，不是连续上界）、`p_volume_mean`（P1 压力按真实 XY/RZ 单元测度的均值）。原格心 `speed/p` 保留。两相每帧另写 `<prefix>_arcs_stepNNNN.csv`，保存接受曲线的局部 [0,1] 幂基系数、材料区间与端点种类；不是续算状态。逐步分相历程追加同一自由表面的平均曲率、曲率平方积分、求积估计与裁剪标记；两相 meta/hydro_moment_case.json 记录生效域、物性、边界及压力/几何时层。B2/B3 离线分析读取这些实际量，见各自 README；不对旧输出猜补完整速度或几何。

**封闭自由滑移与两相重力（2026-09-28，代码已接、待集中数值验收）：** `cchMomentMidpoint` 新增 `boundaryPatches` 类型 `freeSlip`，表示不可穿透、零切向牵引；例如 `leftWall { side west; type freeSlip; start 0; end 1; }`。`value` 可省略，填写时必须为零；不接受 `incomingVelocity` 或对应的 `portPhases`。未列出的边界仍无滑移，RZ 轴沿用轴正则条件；当前一般 RZ 边界迹仍要求计算域含轴。仅滑移壁不需要压力出口，压力采用封闭域体积均值规范；实际通流 patch 仍需压力出口配对。两相 `parametricGals` 现接受 `accelerationX/Y`，通过同一中点的分相材料积分组装 ρg，并计入 `cumulative_external_work`，不使用均匀密度代替。无新增算例结果；Hysing 的 1×2 物理矩形域及两组输入见本节物理矩形域更新。

**常物性供液/出流更新（流体 §28、31）：** `cchMomentMidpoint` 的 `materialModel constant` 接受既有 `boundaryPatches`，支持固定原生 XY/RZ（含 2:1 混合叶）。例如 `supply { side south; type normalVelocity; start 0; end 1; value -0.1; }` 和 `outlet { side north; type pressureOutlet; start 0; end 1; value 0.2; }`。前者给定外法向速度、零切向应力，后者给定完整应力 `(-p I+2μD)n=-0.2 n`；它不是出口单元压力赋值。未列出边界仍无滑移，RZ 轴不能配置开口。可在开口中指定 `incomingVelocity (x, y);` 作为实际流入时的外部动量数据，仅新矩路径接受此键。给定速度入口省略时取指定法向速度及零切向分量；压力出口省略时实际反流会明确拒绝。给出口设置该项不同时约束其流速。省略边界块保留原闭域行为。

非零供液时，声明的 `initialVelocity zero/vortex` 作为初态参考，程序在质量范数中投影到满足给定法向速度和不可压约束的初态；打印修正量，时钟仍为零，初始化乘子不输出为压力。两相初始化直接使用真实 cut-material 质量度量。后续仍直接解同一个中点鞍点系统。普通[常物性 RZ 供液](../cases/hydro_moment_ports_rz/README.md)及[附着两相 RZ 供液](../cases/hydro_moment_supply_rz/README.md)均已完成各自 20 步；两相可用范围见下段，不能将固定湿润口推广为自由面离域或任意接触线模型。

能量历史追加 `cumulative_boundary_traction_work`、`cumulative_boundary_reaction_work`、`cumulative_boundary_kinetic_flux`、`cumulative_inflow_dissipation`、`inflow_volume_rate`、`outflow_volume_rate`、`net_boundary_volume_rate`。前两项是应力边界与给定速度约束的功，第三项以向外为正，第四项来自入口动量差。原 `cumulative_central_work` 包含边界输运贡献；计算总账时不可再重复加减第三、第四项。初帧尚无实际中点，流量诊断列为 0，初始速度本身已满足边界。

**动态网格（流体 §27、52、57）：** [输入与运行说明](../cases/hydro_moment_adaptive_rz/README.md)。运动两相的 `mesh.dict` 使用 `BackgroundMesh`、`QuadtreeMesh { maxInitialIterations 0; }`、`adaptiveMesh true` 和 `Adaptation { maximumLevel 1; cellBudget 128; maxCurvatureCells 0.5; minimumGapCells 2; }`；`cellBudget` 替代 `maxCellNum`。这些分辨率数值只用于粗网格接续示例。初态先满足几何分辨率再初始化速度，后续每步开始前审查接受几何。曲率、相向间距和库存差额来自唯一 Hermite 曲线，保留原加密带/预算与 2:1 分配；`reserveBandCells` 默认 0。省略下面的 `Coarsening` 块保持只细化。旧标量的 `remapScalarTolerance`/`remapRelativeGeometryTolerance` 不改变本路径的 Hermite 库存预算，也不用于动量或能量约束。转移保留曲线、物理时钟与独立库存，当前/上一端点/中点均使用各自完整速度和材料；必要细化不可行或资源不足时明确失败。

**可选粗化（2026-09-28，代码接通，待集中验收）：** 在 `Adaptation` 内加入 `Coarsening { hysteresis 0.5; relativeVelocityError 0.001; }`。两项可省略，默认值如示例，允许范围均为 `(0,1)`；只有完整的 `cchMomentMidpoint + twoPhase + parametricGals` 消费者接受。`hysteresis` 规定合并后的曲率/间距裕量，父单元范围内的几何支撑、预留带、终端最细区和边界分段受保护；每次最多合并一代完整兄弟单元。`relativeVelocityError` 约束三个物理时层各自的全域质量加权速度相对误差，包含浮点消去余量，不是动能损失率或局部空间误差。动量、能量或误差不满足则保留当前网格；数值异常仍报错。`hydro_moment_coarsening.csv` 记录候选状态及已评估时层，`hydro_moment_adaptation.csv` 的 `operation` 和 `maximum_relative_velocity_error` 只记录接受变更。候选拒绝不回退已完成的必要细化。详见[方法 §57](QUANTUM_CCH_Hydrodynamics_TODO.md#automatic-coarsening)。

参数界面的普通 AMR 入口现在允许原已支持的壁面和供排 patch，保留分段端点及原弱法向约束，尚未数值验收。单相常物性动态网格、自由界面离域和一般接触线模型仍未交付；已接的 RZ 断颈模型不等于任意拓扑变化均已支持。物理矩形网格检查 `hx≈hy`，不要求 `Nx=Ny`；例如 1×2 域的 16×32 根网格满足要求。

`hydro_moment_adaptation.csv` 单列实际接受的网格变化和端点动能变化；`hydro_moment_energy.csv` 追加 `cumulative_mesh_energy_change`、`leaf_cells`，能量账扣除累计网格动能变化，不将其混入黏性。零时间变化的历史压力是投影值；普通场输出在下一实际中点求解后发生，仍使用真实 `pressure_time`。

`timeScheme cchMomentMidpoint;` 显式选择新的面一次矩、单元速度矩、线性压力和弱黏性状态，直接调用模块的非线性推进器。接受分离字典、封闭 `planar2D/axisymmetricRZ`、`materialModel constant`、`initialVelocity zero/vortex`；新增两相范围见下段。常物性可用本节所述 `boundaryPatches`；拒绝旧 `smooth`、EHD、`Time.StepRetry`、旧 `Midpoint` 与 `Conservation` 块。RZ 的涡流初值要求 Domain 从 r=0 开始。没有修改省略 `timeScheme` 时的旧默认。

**运动两相入口（验收中）：** `materialModel twoPhase` 配合 `Interface.transport parametricGals`、`advectInterface true`、`reinitializeEvery 0`，在物理矩形 XY 或含轴无旋 RZ 使用，支持均匀网格及由 `RefineBoxes` 指定的固定 2:1 混合叶。闭域可用 `adaptiveMesh true` + `Adaptation` 启用联合动态细化；开口暂不支持动态网格转移。XY 的 `geometry.dict` 使用 `ParametricInterface { name liquid; type ellipse; center (0.5 0.5); radii (0.23 0.17); segments 32; }`，其名字与 `interfaceLevelSet` 一致；没有额外标量水平集。闭合 RZ 使用 `coordinateSystem axisymmetricRZ` 和 `type spheroid; center (0 0.5); radii (0.22 0.17); segments 16;`，径向中心必须为零，半轴为径向/轴向物理长度，子午线两端落在轴上且远离上下壁。Hydro 指定 `interfaceForceGeometry parametric`、`interfaceForceLayer sharp`、中心动量格式和所需两相密度/黏度、`surfaceTension`，不填写 `interfaceThickness`。当前不接受任意接触、重距离或预设牵引；均匀加速度见本节的按相重力更新，固定湿润供液口见下面更新。闭域完整输入与输出见 [XY 算例](../cases/hydro_moment_capillary/README.md)、[RZ 算例](../cases/hydro_moment_capillary_rz/README.md)及[固定 AMR 的 RZ 算例](../cases/hydro_moment_capillary_amr_rz/README.md)。固定 AMR 实际反馈见流体待办 §26，联合动态网格转移与普通三步接续见 §27。XY 冻结运行已完成原 200 步，RZ 冻结运行已完成原 120 步。两者不构成空间/时间收敛或破碎验收。

**附着两相供液更新（流体 §31；物理域扩展见界面 §2.44）：** 含轴 RZ 域可用 `ParametricInterface { name liquid; type supplyCap; radius 0.5; height 0.3; segments 8; }` 初始化 z=z_min+height×[1−(r/radius)²]。`radius`、`height` 须为正且小于实际径向/轴向域长度，`segments` 为 8–512；不填写 `center/radii`。自由面入口边缘固定，另端在轴上规则闭合；底部须有 `start 0; end radius; value <0;` 的 `normalVelocity` 液体入口，边缘外侧紧邻未列出的无滑移壁面。比如 `supply { side south; type normalVelocity; start 0; end 0.5; value -0.1; }`，并在 `Interface` 中填写 `portPhases { supply liquid; outlet gas; }`，覆盖每个声明开口。出口可用 `outlet { side north; type pressureOutlet; start 0; end 1; value 0.2; incomingVelocity (0, 0); }`。实际相/动量进出由共同 C1 速度迹、几何湿润区域和同源方向质量矩决定，纯相声明不能覆盖不相容的进入几何。入口圆盘只参与相测度及边界功，不计自由面面积。边界覆盖须对齐实际原生面；预检会检查。

该输入已接两相惯性、黏性、毛细反馈和开放域中点功，完整字典见 [hydro_moment_supply_rz](../cases/hydro_moment_supply_rz/README.md)。`liquid_volume_drift` 是液量相对初态的变化，在供液问题中须减累计入口液量才可判断守恒。`energy_balance` 为动量离散的动能账，`total_energy_balance` 还包含实际自由面能量变化，两者不能混同为精度误差。自由面离域、移动接触线、多连通域破碎和开口动态 AMR 仍开放；附着端点不支持诊断相检查点 V1，当前输入用 `Restart { model none; }`。

**多个液体分量更新（流体 §32）：** `ParametricInterface` 可改为下列形式。每个子块沿用同一 XY `ellipse`、RZ `spheroid` 或附着 `supplyCap` 初态；原平坦的单分量格式继续有效。

```text
ParametricInterface {
    name liquid;
    type components;
    components {
        lower { type spheroid; center (0, 0.42); radii (0.18, 0.12); segments 16; }
        upper { type spheroid; center (0, 0.67); radii (0.12, 0.09); segments 16; }
    }
}
```

外层仅 `name/type/components`；每个子块包含 `type/center/radii/segments` 或 `type/radius/height/segments`，子块名称不得重复。编号按初始输入顺序从 1 开始，之后随材料保留。它们共用 Hydro 的同一种液体密度/黏度；不是多种液体材料。每个分量独立推进库存，流体消费其总和；多滴可同时切到同一单元。各曲线须是互不接触、不互相包含的正则有效边界，程序会拒绝无法分离的几何，不自动合并。附着子块仍须满足本节供液口与外部相设置。

普通 [hydro_moment_two_drops_rz](../cases/hydro_moment_two_drops_rz/README.md)已完成 20 步。新增 `output/series/hydro_moment_components.csv` 每步输出 `component/liquid_volume/volume_change/geometry_volume/arcs`；原轮廓 CSV 和自由面场增加 `component`。能量历史仍是共同流体与全部自由面的总账。初始化多个分量不触发断裂事件；显式两格 1→2 见下段，嵌套气泡、合并及诊断相检查点 V1 的多分量表示尚未完成。

**显式颈部模型更新（流体 §33）：** RZ 参数两相可设置 `Hydrodynamics.Moment { neckRadiusCells 0.1; }`；默认 0/省略为关闭，正值须 ≤0.5。按接受态 Hermite 的实际径向极小值及相邻轴边格的共用截面检测，在下一物理步前尝试局部事件。§35 将该半径条件与实际候选可行性统一：有限正则颈部尚未满足截面尺度、守恒 C2 封口、动量或能量条件时，记录 pending 原因并继续同一物理方程；只有 ready 候选可发布。无效几何、非有限量或求解故障仍报错，不作为普通输运失败的后备。`hydro_moment_necks.csv` 的评估标志区分未执行项与真实数值；事件 CSV 只记录已接受分裂。 §36 已将封口统一改为直接参数四次 Hermite 构造，匹配保留侧几何 C2、轴端正则性与独立体积，参数根取零曲率连续分支；没有新增用户调参或第二种运行封口。评估输出的 `lower/upper_tangent_r` 为指向各自帽面的单位径向切向，`lower/upper_curvature` 为保留侧有向子午曲率，`lower/upper_cap_length` 为端点参数速度尺度 L（不是整段弧长），`lower/upper_cap_discriminant` 为长度方程的判别式。`ParametricInterface { name liquid; type neckedSpheroid; center (0, 0.5); radii (0.35, 0.3); neckRatio 0.005714285714285714; segments 32; }` 可初始化近断裂液体；`neckRatio` 在 (0,2/3)，其余中心/半轴/分段规则同 RZ spheroid，也可用作 `components` 子块。

普通 [hydro_moment_neck_rz](../cases/hydro_moment_neck_rz/README.md)已自动分裂并完成 20 步。事件单独记录到 `hydro_moment_events.csv`；能量 CSV 增加 `cumulative_event_kinetic_change`、`cumulative_event_surface_change`、`cumulative_event_dissipation` 和 `event_count`。子滴 `volume_change` 从出生库存起算。轮廓和自由面场的 `curvature_model=1` 表示材料来自事件帽面，不代表物理曲率已验收。该 case 初态已欠解析、事件在 T=0，不预测自然断裂时刻。

```text
Hydrodynamics {
    coordinateSystem planar2D;
    timeScheme cchMomentMidpoint;
    materialModel constant; density 1.3; viscosity 0.02;
    initialVelocity vortex; velocityAmplitude 1;
    linearResidualTolerance 1e-11;
    Moment { maximumIterations 24; convection central; }
}
```

`Moment` 可省略；默认迭代上限 24、中心通量 `central`，上限取 1–100。单相、固定物性 Stokes 与运动两相均使用完整 C1 速度；旧低阶时间推进的 `traceResidual` 已退出，选择时会明确报错，现有标准 case 已使用 central。本入口使用 `linearResidualTolerance` 同时检查线性方程与重新计算的实际非线性方程，不使用旧 Midpoint 的速度变化容限。均匀 `accelerationX/Y` 统一在完整物理测试空间积分密度乘加速度，两相使用真实中点分相密度；输出总物理压力，不再走单相专用压力势接线。固定非均匀密度仅支持 Stokes，不冒充有质量输运的变密度 NS。这些接入仍待集中数值验收，静水平衡不能由编译推定。

可在 `adaptiveMesh false` 的 `mesh.dict` 中加入 `RefineBoxes { name { min (x0,y0); max (x1,y1); level n; } }`，只在初始时细化。最多 32 框、level 为 1–8、坐标必须在 Domain 内；足够的 `maxCellNum` 和实际至多 2:1 分裂边仍为运行条件。这不是动态网格转移。旧流体路径不接受这个块。

详细初值公式、输出列与编译/运行命令见 [hydro_moment_xy](../cases/hydro_moment_xy/README.md) 和 [hydro_moment_rz](../cases/hydro_moment_rz/README.md)。`hydro_moment_energy.csv` 使用实际非对角 H 动能和对应物理体积；RZ 面积、体积含旋转权重。单元场提供中心值及一次梯度，面场提供法向均值和一次矩；轴面径向速度与流率为零。终点速度与中点压力分别标记时刻，初始压力 `PressureAvailable=0`。RZ 的 `u_x/u_y` 表示径向/轴向分量。显示投影不替代 GALS 光滑速度；新矩两相例另外输出实际接受界面、两种相分数、表面能及变质量功。当前仍无新矩重启。

### 15.1 保守 RZ 半波长液柱初值（2026-09-25）

`geometry.dict` 可用下面的初值块代替 Point/Line/LevelSet 辅助轮廓：

```text
Domain { min (0, 0); max (4, 5.235987755982989); dimension 2; }
AxisymmetricColumn { levelSet liquid; radius 1; relativeAmplitude 0.001; }
```

三项都必填；`levelSet` 与 `hydro.dict` 的 `interfaceLevelSet` 同名。半径为 `R0*(1+delta*cos(pi*(z-zmin)/Lz))/sqrt(1+delta^2/2)`，`radius=R0>0`、`|relativeAmplitude|<1`，最大半径须小于外径。坐标制须为 RZ，Interface 须开启 conservative 输运、关闭重初始化；毛细配置为 conservative/inventory，不允许混入固定接触点、portPhases、流入/流出 patch 或另一套轮廓。

南北端声明为反射截面，端点半径随真实流体运动。此输入没有增加周期配对或拓扑断裂能力。物性、主时钟和输出仍归原字典；新诊断 `column_mode.csv` 记录原生界面基频幅值。用法与剩余数值验收见[液柱算例说明](../cases/hydro_column_growth/README.md)。

`run.dict/Time` 可选振幅停止，仅用于上述保守 Hydro 初值：

```text
Time {
    units code; dt 0.005; endTime 20;
    Stop { model columnAmplitude; relativeAmplitude 0.03; }
}
```

`model` 和 `relativeAmplitude` 必填；唯一模型为 `columnAmplitude`，阈值须为有限数、位于 `(0,1)` 且严格高于实际初始 `|A|/R0`。App 在首次达到阈值的整个物理步联合接受后停止，`endTime` 仍是上限；不插值阈值时刻。初态/后续几何必须可作为跨两声明截面的单值半径图分析。`Fields.enabled && writeFinal` 时该步强制写末帧；关场输出仍按同一物理轨迹停止。`series/termination.json` 记录原因、实际时钟和振幅，输出失败不生成成功终止记录。省略 Stop 保持原终时行为；普通非保守 Hydro、其他初值和 EHD 拒绝此设置。

### 15.2 保守 Hydro 运行时加密（2026-09-25）

完整 XY/RZ 闭域的 conservative Hydro 可在 `mesh.dict` 显式启用运行时加密，复用独立 Interface 的 `Adaptation` 语法：

```text
adaptiveMesh true;
BackgroundMesh { divide (32, 32); }
Adaptation {
    maximumLevel 3;
    cellBudget 12000;
    maxCurvatureCells 0.1;
    reserveBandCells 2;
}
```

`maximumLevel` 和 `cellBudget` 必填，分别限于 1–8 和大于根网格单元数、至多 250000。`cellBudget` 是候选叶单元上限，内部存储由它分配；不要同时设置 `QuadtreeMesh.maxCellNum`。初始化仍为声明的根网格，`maxInitialIterations` 如提供必须为 0。第一物理步必须能在根网格上解析，拥有真实双时层后才开始加密。

其余可选键：`maxCurvatureCells` 默认 0.125，`minimumGapCells` 默认 8，`maximumGeometryFractionDifference` 默认 0.02，`reserveBandCells` 默认 0，`remapScalarTolerance` 默认 1e-12，`remapRelativeGeometryTolerance` 默认 1e-3。前两项关注主曲率尺度和相向界面间距，不能解释为曲率误差或连续最小间距保证。几何迁移 L1 的绝对预算为相应源时层液相库存乘以 `remapRelativeGeometryTolerance`；不改变守恒门。

App 每个接受步后检查当前与上一时层的界面，取共同界面层级与支持区域的并集；必要的 2:1 平衡及可选 reserve 都计入同一单元预算。reserve 仅随必要加密扩展，超预算时放弃额外 reserve；必要区域超预算或层级饱和则明确拒绝。当前只加密、不粗化。此策略检查已接受几何，未将节点最大速度冒充拟合速度场的连续运动上界；原输运支持和物理 CFL 检查仍有效。固定步长需按预期细网格选择；可选 StepRetry 不能修复空间欠分辨或资源不足。

迁移分别保留 C 当前/上一层的库存和谱系，两层速度用各自时层质量做守恒无散投影；物理压力保留父单元平均与体积积分，物性从目标 C 重建。随后在候选网格上执行下一真实 Interface—NS 步，全部物理检查成功后才交换网格、Hydro、Interface 三个拥有者。准备或物理失败保留旧接受态；输出失败保留已接受物理态并停止，不重复该步。

加密事件写入 `output/series/hydro_adaptation.csv`，包含源时刻、接受/拒绝状态、单元数、两层动量闭合/散度/动能变化、压力积分差、几何库存 L1 和下一物理步散度。普通场文件写实际原生叶单元，Interface 历程跨加密连续追加；未启用 AMR 的算例不新增该事件文件。

半波长液柱保留两端反射截面和偶函数模板，目标每个截面面段必须仍是物理对称边界；移动反射端点不改为固定接触点。现有闭曲线拟合不适用于这种端点，迁移保留原生标量/几何参考及其独立误差检查。宽压力行的毛细尾迹判定检查所有非零系数支撑的实际 C，不能只看面两侧；相分辨阈值、几何纯相、曲率带宽及分支拒绝保持原定义。

目前仅接受 `transport conservative; advectInterface true;` 的闭域完整流体网格；普通节点输运、流入/流出 patch、固定壁面接触和电学耦合不因此获得 AMR 能力。混合网格的压力功、黏性精度与全方法收敛仍须独立验收，运行时接线不等于破碎验证完成。


### 15.3 无电液柱终端颈缩模型（2026-09-25，开发中）

仅用于 conservative Hydro 的 RZ 半波长反射液柱；默认不启用。当前开发主线先使用固定均匀网格，实际喷嘴与电学耦合尚未接通，不能据此宣称已完成电喷雾破碎。

在 interface.dict 中显式指定：

    Interface {
        transport conservative;
        advectInterface true;
        reinitializeEvery 0;
        neckBreakup {
            maximumLevel 0;
            radiusCells 0.4;
        }
    }

两个键均必填。maximumLevel 为 0–24 的声明最细层级，0 表示根网格；若同时配置 Adaptation，必须与其 maximumLevel 一致。radiusCells 必须处于 (0,1)，是上反射端颈半径相对局部径向单元宽度的数值封口阈值，不是曲率细化门限。初始颈半径必须大于该阈值；候选必须在实际物理步内继续缩颈且已达声明最细层级，才允许准备封口。它与 Time.Stop 的小振幅停止模式互斥。

初始化先用普通 C 同步建立共同标量尺度。最细网格的末端薄颈以显式局部模型继续受力：连续模型带延伸至首个半径至少为 2h、曲率可用的径向截面，使用当前共享原生标量的曲率及 RZ 轴极限，方法号为 4（EventModel）。该 2h 条件用于模型带覆盖，不能触发封口。光滑区保持原曲率提供者；弧长正则化在模型带处分段。同一原生径向图上的近距离投影竞争可按最近点继续，并记录候选曲率差。对启用终端模型的状态，当库存尾迹不超过现有相分辨阈值、整个压力支撑均为同一几何纯相时，还允许跨连续单值 r(z) 图段选择最近点；这可跨过普通/模型曲率提供者的边界，但各候选须有自己的可用曲率，不在两者之间插值或扩大原生模型带。折叠、反向图段、不可用候选和不同表面仍拒绝，可分辨界面不因此获得跨图延拓。尾迹分类由 C 和完整压力支撑决定，与投影点是否位于搜索带内无关。

封口逐格匹配原独立 C，保留远区标量并生成截面脱离身份；事件时间记录本次物理步区间，不插值物理奇异时刻。事件带两轴单元之外必须为几何纯气，允许其低于现有相分辨阈值的库存尾迹原值保留，并记入几何—库存差及未分配库存；不删除尾迹，也不把它当作已解析液滴。断裂前即使采用终端曲率，几何仍进行 C 同步；只有实际截面脱离后，模型区才保留物理 jet 预测，避免仅凭相同高度库存重新接回原截面。有截面脱离记录且仍保持原生连通性、轴端和最细层的端部，可在后续运动中局部转弯，曲率保留其当前符号；封口前仍要求径向液体图。该放行用于模型选区和 AMR 终端走廊，跨折返的多分支力投影竞争仍按原规则拒绝。C 始终由共享面通量更新，模型区的几何—库存差另记。目标物性取同一 C，当前步毛细力取候选 Γh，经真实保守 NS 成功后与界面共同接受；失败沿用整步事务恢复。已有保守动量为一步质量通量格式，黏性只使用当前区间，不在此处虚构 AB2 或电荷历史。

只有启用该模型时，interface_conservative.csv 才增加模型区单元数、库存占比、几何 L1 占总库存比、事件数、实际模型力面数、最大力及投影候选力差。场输出的 capillary_source CSV 追加同图候选数和曲率差。投影力差是局部敏感性观测，不是完整误差界；模型受力也不代表已解析奇异曲率。普通算例原 CSV 列保持。

最新构建还写 output/series/interface_events.jsonl：每个接受事件一行，含物理步、时间区间、事件类别、液体父子 ID、脱离截面、当前表面及其父 ID、各液体的已分配 C 库存/几何体积，以及总库存和未分配尾迹。时间区间不是插值出的连续奇异时刻；体积属于当前计算域，半波长反射域的 1→1 截面脱离需按镜像关系解释。事件记录仅在联合接受后的输出中追加，跨 AMR 保留写入位置。

保守 Hydro/EHD 还写 `output/series/interface_components.csv`，每个接受时刻按实际液体 ID 分列已分配 C 库存、几何体积和触边标记；最后一行的保留 ID=0 只记录未分配库存，并标识该时刻记录完整。各液体库存加未分配项才是全域总库存。触边包含轴线/反射截面，不能用来判断喷嘴脱离；半域镜像换算只适用于声明的对称算例。该输出不受场输出频率控制，AMR 继承同一输出流，不增加字典项。旧冻结程序没有此文件时，分析脚本明确报告缺失，不补造逐滴体积。

当 Fields.enabled 为 true 时，新构建额外保存事件前的原接受态、事件步及下一接受步的场文件，沿用原字段和文件名，不重复追加时间序列。前帧输出失败直接停止；失败候选不写入事件日志。旧冻结副本只按原间隔保存场；首次事件运行已在第 199 步、t=0.4975 的真实截面脱离中触发输出，第 198/199/200 步原生界面与事件时钟、身份记录已核对一致。完整后续运动与数值验收仍开放。

该主包用 `tools/analyze_hydro_column_breakup.ps1 -Run <运行目录>` 汇总冻结输入、共同接受时钟、库存/几何差、事件区占比、截面端点及事件后 NS 历程，写入运行目录的 `breakup_analysis.json`。运行中加 `-InProgress`，脚本只读取完整行；终止失败保留实际退出码，不能当作完成。半域镜像后的总 C 含未分配尾迹，不能直接冒充各滴库存；端点回缩也不等于整滴下游平移。该汇总供两档网格及一个时间代表比较，不自动宣布数值验收通过。

与 Adaptation 同时启用时，终端径向连通区域仍请求细化到声明最细层；该层无法满足的曲率尺度另记为终端模型请求，不冒充已解析曲率。区域以最细径向尺度及 maxCurvatureCells 界定，和实际受力模型带、封口阈值分别记录；普通区域曲率、拟合质量、相向间距与几何—库存拒绝保留。hydro_adaptation.csv 增加 terminal_curvature_requests 和 maximum_terminal_curvature_cells；没有换网格但存在该请求的成功物理步记 terminal_model。

普通曲率仍优先采用有效高度函数，其次采用同分支局部拟合。2026-09-25 增加有界模板选择：3.5h 模板仅在几何拟合失配时尝试 2h，双侧各至少 h 支撑和原质量门限保持；其他拒绝原因不经此策略放行。无需新增字典项，也不扩大事件曲率范围。

事件后的 AMR 可迁移当前/上一时层，但模型单元必须保持相同物理范围、层级及独立 C，只改变周围支持区。迁移保留原生标量和事件谱系，不用整条曲线重拟合封口；仍检查原标量误差与几何—库存预算。已有制造封口后三步 NS、160→208 格双时层迁移及下一步 NS 控制。模型单元的继续细分/合并、零时间重同步、表面材料和电荷转移尚未实现，相关调用拒绝；真实断裂后的长程与网格/步长验收仍开放。[此前 AMR 开发输入，已按用户要求停止](../runs/hydro_column_breakup/20260925_162212_r16_amr_refined_projection/study-input.json)使用 32×42 根网格、Adaptation.maximumLevel=1、cellBudget=6000、reserveBandCells=2，neckBreakup.maximumLevel 同为 1；最细 R/h=16，dt=0.0008838834765、终时 1。radiusCells=0.4 对应最细物理封口尺度 0.025R；这与旧 R8 轨迹分别记账，网格/时间精度尚未验收。AMR 迁移速度现用同一质量投影做有界残差校正，使下游保守体积通量可用；仅该投影可做最多三次线性残差精化，原线性容差及动量/能量与 C 门限保留，无新增字典项；[本包证据](../cases/module_history/output/20260928_184817_build_archive/interface_event_coupling_evidence.json)区分原失败运行、AMR 扩展及尾迹修复的程序身份和检查范围。

**近断裂主算例（2026-09-25）：** 现使用 [cases/hydro_column_breakup_near](../cases/hydro_column_breakup_near/README.md)。`geometry.dict` 将 `relativeAmplitude` 改为 `0.954754144451598`，`mesh.dict` 使用 `divide (64, 84)`、`adaptiveMesh false`，`neckBreakup.maximumLevel=0`，`radiusCells=0.4`；细网格从初始化生效，实际径向触发尺度为 0.025R。`run.dict` 为 dt=0.0005、endTime=0.2，保留原 StepRetry，`Log.level=detail` 输出分阶段计时。零初速，初始化后原生颈半径约 0.04190R；解析输入半径 0.0375R 与保守重建后的半径分开记录。终端曲率模型从初始细颈区域参与计算，接受事件仍要求真实收缩物理步。该初值用于缩短开发过程，不是增长率/连续体断裂时刻验收输入。

运行 `powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/run_hydro_column_breakup.ps1 -Threads 12 -Mainline`。脚本冻结程序和字典到新运行目录，以环境变量配置 OpenMP 12 线程、OpenBLAS 1 线程、空闲工作线程休眠；不新增物理字典参数。`-EndTime` 在启动前指定终时，`-Threads` 控制计算线程。Interface 通量并行仅由显式承诺查询并发安全的调用者启用，任意用户回调仍默认串行；当前保守 Hydro 路径已经接通。

