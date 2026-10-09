# LAMP 字典约束

现行版本：2026-10-03。目标是让用户能明确知道算什么、哪些配置生效、由谁消费。当前独立 Interface 的冻结输入见 [G04 schema](api/G04_Interface_Dictionary_Schema.md)。G03 的通用数据思想保留，单文件布局退役；历史配置与报告不作为当前可运行示例。

## 1. 固定职责

按数据的职责分文件，不能因为某个 App 最先使用它，就塞进该 App 的模块字典。

| 文件 | 唯一职责 | 独立 Interface 当前使用 |
| --- | --- | --- |
| `run.dict` | App 选择、主时钟、输出、日志、restart 和耦合调度 | Run.app、Time、Output.Fields |
| `geometry.dict` | 计算域、坐标/测度、相和边界 identity、初始几何 | Geometry.coordinates、Domain、Initial.phi |
| `mesh.dict` | 基础网格、初始加密、AMR 层级、网格变更策略 | Mesh.resolution、initialMesh、Refinement、Regrid |
| `fields.dict` | 独立给定场、外部场输入及其表达式参数 | Fields.Velocity、Constants |
| `interface.dict` | 界面模型/方法、界面数值选项、该模块的边界条件 | Interface.method、Boundary |
| `materials.dict` | density、viscosity、permittivity、conductivity、surface tension 等共享物性 | 未使用，不提供占位文件 |
| `hydro.dict` / `electricity.dict` | 对应物理模块的模型、初值、边界和数值选项 | 未启用，不读取 |
| `solvers.dict` | 线性/非线性系统、solver、preconditioner 及求解容限 | 当前 Interface 无此配置 |
| `units.dict` / `normalization.dict` | 参考尺度和单位转换 | 当前独立 Interface 仅 code units；Time.units code |

Run 块本身只选 App；Time、Output 等是 run.dict 内的同级块。每个 App 只实现实际支持的子集，文件职责保持一致，其他字段明确报错。不能接受一个字段却不消费它，也不能加入尚未实现的占位模型。

本次落实到独立 Interface。尚未迁移的 Hydro/EHD 等模块保留其冻结 reader，不因更新本约束就声称已符合新布局；后续对应开发任务按此职责归位，不做无验证的全仓替换。

## 2. 字段与 owner

- 每个输入值只在一个文件声明、解析一次；由 App 汇总已验证配置，数值内核不读取字典或路径。
- 时间与输出不能重复出现在 interface/hydro/electricity/mesh 字典中；网格与域也不能复制到模块字典。
- 初始 phi 描述初始几何，归 geometry.dict；界面边界 phi 归 interface.dict。给定速度归 fields.dict，耦合时实际速度由 Hydro 发布，Interface 消费只读 view。
- Fields.Constants 仅保存这些输入表达式共享的参数，供初始 phi、边界 phi 和给定速度解析引用，不是文本宏或文件 include。共享材料、solver 容限、网格/时钟控制不能挪入 Constants；未知、保留名称或未使用参数报错。
- 不用 problem/case 名称在生产代码选择初值、速度、解析解或验收；专用参考与 benchmark gate 归 tests/cases。

## 3. 严格读取与版本迁移

未知字段/块、重复 key、非法类型、垃圾尾缀、未使用文件或选项均报错。跨文件组合约束也必须验证，例如 endTime/dt 的整数步数、Regrid 必须发生在末步之前、Refinement 必须在 Domain 内。

默认值由 reader 的 ConfigurationTrace.defaultValue 记录，不能依赖 case 名、目录、隐藏环境或网格尺寸。派生量不能假装成默认输入；当前 steps 由 Time.dt/endTime 按公共时步计数规则一次推导，不另设 Time.steps。

旧 schema 只允许通过明确且有限期的 reader 迁移处理；不能同时接受多份等价字段、静默覆盖或失败后切换方法。G04 不兼容单文件 Interface 配置：完整模板与 fixture 已迁移，旧错放字段直接拒绝。新增/修改 schema 必须同步更新 reader、模板、测试、生效配置输出、API 和迁移说明。

## 4. 单位、物性与容限

区分物理参数、数值方法、求解器容限和物理验收容限。使用具体名称，如 linear.relativeTolerance、nonlinear.relativeTolerance、conservation.phaseVolumeTolerance、verification.adjointTolerance，不共享含义模糊的 tolerance。验收参考不参与生产修正或容差调节。

共享物性由 materials.dict 以 phase ID / phase pair 绑定，模块只消费所需属性。模型允许的精确零值，例如 viscosity、conductivity、surface tension，不能被 reader 意外禁止。code/SI 的支持范围及转换 owner 必须写明，不能把未转换的数字称为 SI。

## 5. 来源、输出与错误

统一使用 output/meta/effective_config.dict 和 configuration_trace.csv，不建立第二套追踪。生效配置写实际值与默认值，保留所属文件；独立 Interface 另导出 effective 下同名的五份字典用于回读，仍由同一个 writer/配置 owner 生成。所有消费文件都记录来源哈希。

configuration_trace.csv 至少包含 field、file、line、raw_value、effective_value、units、consumer、defaulted、used、unused_reason。错误必须给出真实 file:line、block.field、原值和期望类型/范围/consumer。跨文件校验不能把 geometry/mesh/run 错误统一定位到 interface.dict。

当前字段的“能解析”“能完成推进”“精度通过”“物理验收通过”分别报告。未实现示例标记为目标；历史 [字典参考](history_todo/QUANTUM_Dictionary_Reference.md)和 I12p3b–e/G03 文件用于追溯，当前输入以本约束及对应现行 API 为准。
