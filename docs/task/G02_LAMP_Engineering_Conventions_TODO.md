# G02 / LAMP 工程规范与旧路线清理 TODO

依据：用户 2026-10-02 的补充要求。本任务承接 G01 的目录、通用运行入口与头文件迁移，并将程序更名为 LAMP（Levelsetbased Adaptive Multiphysics Platform）。

## 执行范围

- [x] 固化薄约束：目录与文档、生产代码与兼容层、测试与 benchmark、运行与调试。
- [x] docs 根目录只保留关键约束与索引；API、任务、方法 PDF、使用指南各归子目录。
- [x] 小型契约测试归 tests；完整时间窗、收敛扫描与论文 benchmark 归 cases/benchmarks，默认 CTest 不跑完整研究。
- [x] 清理 tools 根目录；只保留通用入口，通用后处理、文档构建分别归子目录。
- [x] 在源码快照和逐文件清单之后，删除已放弃的独立界面路线、专属测试与脚本。保留 Hydro / 电喷雾仍调用的代码时，归独立兼容层并冻结其依赖边；不作为当前 GALS 入口的备用路线。
- [x] 项目、可执行程序、构建目标、运行信息和当前使用文档统一为 LAMP；历史证据不改写。
- [x] 自动检查目录、依赖与测试分类，避免后续开发重新混入临时文件。
- [x] Release 构建、CTest、完整动态 AMR 算例、benchmark 短运行与调试构建验证；准确记录历史缺失 fixture 和未完成的物理能力。

## 删除边界

当前生产界面主线是纯 Hermite–GALS / FirstJet。老 levelSet、PLIC、conservativeGals、parametricGals 的独立 Interface App 分发与专用研究脚本退役。仍被旧 Hydro / EHD 调用的边界、载荷和状态适配必须显式隔离；本次不伪造一个尚未实现的新 GALS–Hydro 耦合来替代它。

用户本次清理要求覆盖先前冻结文件的保留要求；原冻结验收及源码保存为历史证据，不再要求失效实现继续存在于生产目录。数值算法、容差与通过标准不因工程整理而放宽。

临时审计脚本放 runs/checks/G02，源码快照与迁移/退役清单放 runs/baselines/G02_20261002。交付记录见 G02_LAMP_Engineering_Conventions_Report.md。

完成记录：2026-10-03；见 G02_LAMP_Engineering_Conventions_Report.md。

临时工程脚本已在交付时归档至 runs/baselines/G02_20261002/engineering_scripts；不是后续用户入口。
