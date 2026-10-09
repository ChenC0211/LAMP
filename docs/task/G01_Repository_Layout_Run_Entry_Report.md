# G01 / 目录与通用入口报告

G01 的工程工作完成，并由用户的后续要求扩展为 [G02 / LAMP 规范与清理](G02_LAMP_Engineering_Conventions_Report.md)。当前入口、规范和验证结果以 G02 报告为准。

G01 固化 24 份 Hermite–GALS 配置、通用 run_case 启动器、CMake 构建与 VS Code 调试入口；输入/程序快照、时间戳运行、真实退出码和 ParaView 场路径接通。21 个内部头从 src 迁入 include/Detail，原生产实现除 include 路径外逐字审计保持。

旧结果迁移清单见 `runs/MIGRATION_G01.md` 和 `runs/baselines/G01_20261002/result_migration.json`，历史结果与验收包保持追溯关系。G01 源码/运行快照和审计不因 G02 的退役而改写。

G01 初次完整构建补齐旧独立 App 的链接源；初次 CTest 暴露旧守恒研究失败和原始 fixture 缺失，不能称作全部通过。G02 后当前默认有效检查已全部通过，缺失两项明确 Disabled，完整研究改为显式 cases/benchmarks 入口。

工程整理按用户授权处理，不改变数学方法、容差或物理参数；原 G01 的兼容 wrapper/专属脚本已按 G02 退役，程序改为 LAMP。使用说明在 docs/guides/Run_And_Debug.md。
