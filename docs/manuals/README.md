# QUANTUM 使用者说明书

**状态：historical/reference（2026-09-19 使用快照）。** 三册说明书用于查询当时的模型、字典和运行方式，不维护 V2 开发任务、现行输入支持范围或验收状态。当前任务唯一入口是[分模块任务清单](../roadmap/05_LAMP_Module_Task_Roadmap.md)，规范与冻结方法见[文档索引](../README.md)。

| 说明书 | 内容 | PDF | LaTeX |
| --- | --- | --- | --- |
| 第一册：整体架构与数学模型 | 分层与依赖、模块状态、流体/界面/电学/PIC 方程、耦合顺序、适用范围 | [阅读](QUANTUM_Architecture_Models.pdf) | [源码](QUANTUM_Architecture_Models.tex) |
| 第二册：运行与结果阅读 | 本机编译、十五个样例、启动与恢复、输出、Tecplot、常见问题 | [阅读](QUANTUM_User_Guide.pdf) | [源码](QUANTUM_User_Guide.tex) |
| 第三册：字典配置指南 | 按 App 区分文件、字段含义/单位/组合、模型与边界、容限、完整 EHD 输入 | [阅读](QUANTUM_Dictionary_Guide.pdf) | [源码](QUANTUM_Dictionary_Guide.tex) |

运行前核对[当前算例索引](../../cases/README.md)和实际输入/reader。手册的示例和通过记录仅说明当时状态，不等于当前实现已完成验收。

旧语法与职责可查询[历史字典参考](../history_todo/QUANTUM_Dictionary_Reference.md)和[历史代码组织规范](../history_todo/QUANTUM_Code_Organization_Standard.md)；V2 新增/修改以 [03 字典规范](../03_LAMP_Dictionary_Rules.md)和 [02 架构规范](../02_LAMP_Architecture_Rules.md)为准。[旧主计划](../history_todo/QUANTUM_Modular_Consolidation_Plan.md)与[旧物理进度](../history_todo/QUANTUM_Hydro_EHD_Progress.md)仅供追溯，不产生当前任务。

历史说明提到通过 [tools/docs/build_manuals.ps1](../../tools/docs/build_manuals.ps1) 重建，三个 TeX 共享 `quantum-manual.sty`。当前目录缺少该样式文件，因此不能把旧重建命令宣称为现成可用；本任务未补造样式、重编 PDF 或重新验证手册中的算例。

原有 PDF/TeX 保持不变。它们是历史阅读快照，与总纲指定的三份冻结方法 PDF 的数学权威不同。
