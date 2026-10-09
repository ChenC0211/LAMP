# S07 / 历史文档降级交付报告

**结论：PASS（S07 文档权威与入口验收）。** 执行日期：2026-10-01。任务定义见 [S07 TODO](S07_Historical_Documents_TODO.md)。

## 目标与变更

当前任务唯一入口仍为 [05 / 分模块任务清单](../roadmap/05_LAMP_Module_Task_Roadmap.md)。本次只修正文档权威与可达入口，不改变数学方法或实现。

- `docs/history_todo/` 的 14 份 Markdown 全部增加顶部 `historical/reference` 声明，包括两个长 TODO、旧主计划、性能方案、规范/指南与旧入口。声明明确原文“当前”“下一步”仅指写作当时，不产生新任务、不证明 V2 gate 通过。
- 历史正文原始字节保留，包括其中旧链接；未把缺失旧结果包装为已验证证据。
- 新建 [docs/README.md](../README.md)，直接链接六份现行规范、三份冻结方法、历史入口与任务记录。
- 更新项目、算例和说明书 README 中的任务权威声明；旧标准算例状态注明是 2026-09-29 快照，不暗示进程现在仍在运行。
- 00 总纲修正三个冻结 PDF 的相对链接；02 架构与 03 字典规范各修正一个历史参考链接。01、04、05 不改；规范规则、公式与 gate 内容不变。

三份 `docs/` 根目录冻结方法 PDF 仍是开发总纲指定的数学主线。历史目录保存配套 TeX 不改变其权威；本任务没有重编或修改任何 PDF/TeX。

## 验证与证据

证据目录：[docs/.build/S07_20261001](../.build/S07_20261001/)。

| 层级 | 验收内容 | 结果 |
| --- | --- | --- |
| T0 | 历史顶部声明；现行入口文件与锚点可达；六份规范、三份冻结方法可访问 | PASS：14 份历史顶部声明、136 个现行/新增链接，无失效目标 |
| T1 | 14 份历史正文逐字节保留；受保护文件清单/哈希不变；规范仅链接目标变化 | PASS：14 份正文、1,154 个受保护文件；六份规范只发生五处允许的链接替换 |
| T2 | 项目入口 → 现行文档 → 路线图与旧 TODO 的权威关系复核 | PASS：现行入口指向 05，历史顶部声明明确不再产生任务 |

17 份编辑前文档快照（14 份历史 Markdown、3 份已有 README）保存于 `before/`，哈希见 `documents_before.csv`。六份规范的内容比较参照 S01 冻结归档中的 `normative_reference_S01/`，不是本次编辑前的独立快照；比较时仅允许指定的五处链接目标替换，并忽略 CRLF/LF 差异。受保护对象清单见 `protected_before.csv`；详细结果在 `verification.json`、`historical_verification.csv`、`normative_verification.csv`、`link_verification.csv`。最终文档身份见 `delivery_manifest.json`。

此任务为文档权威与入口修改，未编译程序、未新增测试、未重跑物理算例。T0/T1/T2 在本报告中分别指文档结构、不变性和入口验收，不作为任何 I/H/E 数值 gate 的证据。

## 缺失资料与边界

- 旧索引中的 `hydro_midpoint_rz`、`hydro_moment_capillary`、`mesh_refinement`、`pic_smoke` 及三个 spray smoke 目录当前缺失，已改成明确的历史索引项。
- `cases/CASE_TEMPLATE.md`、`cases/module_history` 和 `docs/manuals/quantum-manual.sty` 当前缺失；入口如实注明，不补造历史资料或声称手册可直接重建。
- 历史正文、单个 case README 和 `tests/README.md` 中的旧链接/旧记录不在本次批量修复范围。历史目录顶部声明限定该目录中的原文；项目/文档入口统一规定任务只来自 05。现行入口与新增链接单独核验。
- 本次修改了文档索引及根 README，因此 S01 的“当前所有源文件等于基线”检查后续可能报告预期的文档差异。S01 原始归档与交付证据仍应保持相同 SHA-256，不回写其历史报告。

## 交付与下一步

交付本 TODO/报告、14 份历史顶部声明、4 个现行入口（含新建文档索引）及 3 份规范的链接修正；共 23 份 Markdown。production 文件 0，测试文件 0，构建/程序/输入不变。

本轮只处理 S07。下一项为 S08，须单独定义 TODO 与验收；Gate S 尚不能据此宣告完成，I01 未开始。
