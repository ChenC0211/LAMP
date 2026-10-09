# LAMP 目录与文档约束

状态：现行工程约束，2026-10-02。与旧路径约定冲突时，以本文件为准；数学规范不因此改变。

| 目录 | 允许的内容 |
| --- | --- |
| `include/{Core,Solvers,Modules,Apps}` | 生产头文件；模块内部头归 `Detail` |
| `src/{Core,Solvers,Modules,Apps}` | 对应的 C++ 实现；每个翻译单元只有一个 CMake owner |
| `tests/contracts/<模块>` | 短契约、局部算子、输入及输出格式检查 |
| `tests/support`、`tests/fixtures` | 测试专用辅助头、最小静态输入；不得成为生产内核的依赖 |
| `cases/<算例>` | 完整验证/应用算例的配置、静态输入、参考数据和说明 |
| `cases/benchmarks/<研究>` | 完整时间窗、收敛扫描、论文 benchmark；专用驱动只编排生产算子与独立参考 |
| `runs/<算例>/<时间戳_标识>` | 单次运行的输入、程序快照、日志、结果及状态；不得覆盖旧运行 |
| `runs/checks/<检查>` | 自动检查结果及临时调试；新调试工作再建时间戳子目录 |
| `runs/baselines/<任务>` | 冻结验收证据、来源/哈希、退役源码包；不作为当前入口的运行依赖 |
| `build`、`bin` | CMake 中间产物与可执行程序/DLL；不可成为算例输入 |
| `tools` | 仅通用构建/运行入口和职责明确的工具子目录 |
| `docs` 根目录 | 总纲、数值/工程约束及索引；不放单次任务、API 草案、图片、PDF |

文档分区固定为：`api`（接口）、`task`（TODO/报告/设计变更）、`roadmap`（现行路线图）、`guides`（使用）、`methods`（冻结数学方法）、`assets`（插图）、`archive`/`history_todo`（历史）。历史报告中的数字、结论、归档哈希不改写；移动后维护链接和迁移清单。

禁止在 `src/include/cases/docs` 中随手增加 `tmp`、`debug`、`output`、`old`、`backup`、`(1)` 副本或生成文件。尚未确定归属的试验先放 `runs/checks/<任务>/<时间戳>`。需要成为长期实现时，完成归属、接口、依赖、验证和文档后再进入生产树。

目录改动必须同次更新 CMake、消费者、索引和相关检查；删除实现前保存可追溯快照与清单。空的迁出目录及时移除。`repository_layout` 和 `architecture` 是必过检查。
