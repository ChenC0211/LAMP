> **状态：historical/reference（历史资料，仅供追溯）。**
> 本文原始正文保留为历史记录，其中“当前”“下一步”、TODO、旧规则与通过记录均指写作当时；不再自动产生新任务，也不能证明 V2 gate 通过。
> 现行任务唯一入口：[分模块任务清单](../roadmap/05_LAMP_Module_Task_Roadmap.md)。文档权威见[开发总纲](../00_LAMP_Development_Charter.md)，完整入口见[文档索引](../README.md)。
> 旧字典/职责说明可供查询既有实现；V2 新增或修改仍遵守现行六份规范。历史正文中的链接保留原样，目标缺失时不能作为本机已验证证据。
> 本目录中的旧设计、TeX/PDF 和索引同属历史参考；总纲指定的三份冻结方法 PDF 仍以 `docs/` 下的权威文件为准，不因本目录的历史标记而降级。

# 文档入口

当前模块化状态与共同排期在[主计划](QUANTUM_Modular_Consolidation_Plan.md)第 8 节维护；界面与流体分别执行 [Conservative AMR Hermite–GALS 待办](QUANTUM_Conservative_AMR_Hermite_GALS_TODO.md)和 [CCH 流体待办](QUANTUM_CCH_Hydrodynamics_TODO.md)，分别验收后推进整体破碎与带电接入。电喷雾破碎仍是最终目标。历史报告中的待办、检查数量和旧电脑路径不作为当前执行要求。构建和运行先看[项目说明](../README.md)，完整工作包收尾才集中验证。

面向使用者的三册 **LaTeX + PDF** 见[说明书索引](manuals/README.md)：[整体架构与数学模型](manuals/QUANTUM_Architecture_Models.pdf)、[运行与结果阅读](manuals/QUANTUM_User_Guide.pdf)、[字典配置指南](manuals/QUANTUM_Dictionary_Guide.pdf)。按 2026-09-19 当前实现编写，包含真实算例、参数解释和模型适用范围。

| 要找什么 | 现行入口 |
| --- | --- |
| 目录、依赖、接口及新增代码规则 | [代码组织规范](QUANTUM_Code_Organization_Standard.md)第 1–4 节 |
| 实际支持的输入、默认值和迁移 | [字典参考](QUANTUM_Dictionary_Reference.md) |
| 模块化完成条件、证据与剩余工作 | [模块化收紧计划](QUANTUM_Modular_Consolidation_Plan.md) |
| 水平集新方法的实施待办、依赖与验收 | [Conservative AMR Hermite–GALS 待办](QUANTUM_Conservative_AMR_Hermite_GALS_TODO.md)，[主计划排期入口](QUANTUM_Modular_Consolidation_Plan.md#levelset-implementation) |
| CCH 无电流体方法、独立验收与接入 | [CCH 流体待办](QUANTUM_CCH_Hydrodynamics_TODO.md) |
| 选输入、运行及查看输出 | [本机算例索引](../cases/README.md)、[使用指南](QUANTUM_Hydro_EHD_Case_Guide.md) |
| 尚未解决的数值问题 | [物理进度与限制](QUANTUM_Hydro_EHD_Progress.md) |
| 已注册检查及其范围 | [测试说明](../tests/README.md) |
| 历史计算与阶段证据 | [运行目录](../runs/)、[分类报告](validation/)；以主计划中的具体记录链接为准 |

输入语法以字典参考为准，代码组织以组织规范为准；指南只说明用法，不另立一套规则。新增能力需要同时落实模块消费者、支持范围、对应检查及这两份规范中的相关条款。结构检查的实际覆盖及保留兼容边见主计划第 8 节。

数学推导和设计资料：

- [Conservative AMR Hermite–GALS 设计 PDF](QUANTUM_Conservative_AMR_Hermite_GALS_Method.pdf)、[LaTeX](QUANTUM_Conservative_AMR_Hermite_GALS_Method.tex)：2026-09-26 当前方法改造依据；实现与验收状态见[待办](QUANTUM_Conservative_AMR_Hermite_GALS_TODO.md)，设计能力尚不代表已实现。
- [旧水平集方法设计 PDF](QUANTUM_LevelSet_Method_Design.pdf)、[LaTeX](QUANTUM_LevelSet_Method_Design.tex)：2026-09-21 的 CLSVOF 型守恒、分支曲率、拓扑与电荷接续方案，保留作旧路线追溯。
- [流体数值方法 PDF](QUANTUM_Hydrodynamics_Numerical_Methods_Revised.pdf)、[LaTeX](QUANTUM_Hydrodynamics_Numerical_Methods.tex)、[源码身份](QUANTUM_Hydrodynamics_Numerical_Methods.sources.json)。这是标注日期的数学/源码快照，不是当前模块化待办。
- [NS 时间投影笔记](QUANTUM_NS_Projection_Time_Discretization.pdf)、[LaTeX](QUANTUM_NS_Projection_Time_Discretization.tex)。
- [模块蓝图](QUANTUM_Modular_Model_Blueprint.pdf)、[LaTeX](QUANTUM_Modular_Model_Blueprint.tex)和 [Popinet 分裂设计](QUANTUM_Popinet_EHD_Split_Scheme.pdf)属于设计资料；拟议语法和模型不能直接当作已实现输入。
- [早期组织计划](QUANTUM_Code_And_Dictionary_Organization_Plan.md)和 [EHD 实施计划](QUANTUM_Electrohydrodynamics_Implementation_Plan.md)用于追溯原安排，当前执行入口仍是模块化主计划。

[主计划 PDF](QUANTUM_Modular_Consolidation_Plan.pdf)及[LaTeX](QUANTUM_Modular_Consolidation_Plan.tex)是 2026-09-16 的 v1.8 阅读快照，实时内容维护在 Markdown；旧 v1.3/v1.4 同样只作历史资料。原导出脚本 `tools/docs/export_modular_consolidation_plan.py` 未迁入本机，不能按旧文档直接运行。新的使用者说明书有独立的[构建脚本](../tools/docs/build_manuals.ps1)，不覆盖这些历史计划快照。

历史链接若目标未迁入，保留原路径并明确标注，不能作为本机已验证证据；不按旧记录重建虚假的历史结果。新的计算输入和结果放 `runs/`，标准 `cases/` 只保留输入，源码快照不作为第二个开发入口。
