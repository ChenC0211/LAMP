# 文档构建

这些脚本用于历史数学方法和手册排版：`build_interface_design.ps1`、`build_levelset_method.ps1`、`build_manuals.ps1`。它们不是求解器构建或算例运行入口。

冻结原件在 docs/methods，历史手册在 docs/manuals。重新生成文件应进入 runs/checks 的任务目录并复核，不覆盖冻结 PDF 后沿用原哈希。当前使用说明维护 docs/guides 的 Markdown。
