# LAMP GitHub 初始化记录

日期：2026-10-09。远程目标为账号 ChenC0211 下已有的私有仓库 [LAMP](https://github.com/ChenC0211/LAMP)。远程已有一个仅包含两行 README 的初始化提交 `2a3c0be`；通过合并保留该历史，README 使用本地完整项目说明。GitHub 的 size 为 0 不表示没有提交，先前的空仓库判断已更正。当前源码首次纳入 Git，快照标签为 `snapshot-20261009`，对应本地导入提交 `bb9c17e`；main 另包含远程历史合并记录。该标签保存工程现状；各物理能力仍以现行任务报告和路线图为准，首次提交不改变已有 PASS、FAIL、缺失或未验收结论。

## 纳入范围

提交 src、include、CMake 共享配置、tools、tests、cases、docs、共享 VS Code 构建/调试配置、忽略与属性规则、Git 使用说明、本机 preset 示例及 Pull Request 模板。日常 runs、build、bin、虚拟环境、本机 CMakeUserPresets 和凭据不提交。

本次未修改生产源码、数值测试、算例输入或冻结方法 PDF。工程架构检查器修复了跨盘符外部 include 路径判断，并加入 Windows 跨盘符自测。新增仓库管理文件使用 LF；其他文件保留原始字节，避免首次导入引起换行转换或历史哈希变化。

三个历史 ZIP 保留原字节：task.zip 有 5 项、task1.zip 有 5 项、task88.zip 有 2 项内容不同于现行 Markdown；并非全部重复备份。历史文档和冻结 PDF 保留现有目录及名称。

## 本机构建迁移

忽略的 CMakeUserPresets 已改用 `D:/Chen/openLib` 下的 Eigen/SuperLU/CBLAS；local、local-debug 和 local-benchmarks 分别使用新的 windows-local、windows-local-debug 和 windows-local-benchmarks 构建目录。旧缓存记录 `C:/Chen/QUANTUM`，未删除旧构建产物。新机器使用 docs/guides 的 preset 示例填写本机路径。

## 检查

初始化检查证据保存在本机 `runs/checks/github_bootstrap/`，不随日常运行目录上传。

| 检查 | 命令或范围 | 结果 |
| --- | --- | --- |
| 全新 Release 配置和 LAMP 构建 | `tools/build.cmd -Jobs 4`，新 windows-local 目录 | PASS |
| 架构检查器构建 | `cmake --build --preset local --target quantum_architecture_check --parallel 2 -- /nr:false` | PASS |
| 架构及目录 | architecture、architecture_self_test、repository_layout，独立 CTest 窗口 | 3/3 PASS，3.29 秒 |
| 应用入口 | lamp_help、interface_first_jet_app_contract，另一 CTest 窗口 | 2/2 PASS，125.07 秒 |
| 跨盘符回归 | 新外部目录正例、原 133 项策略自测及 31 项拒绝 fixture | PASS |
| 上传候选 | 构建/运行产物、本机 preset、常见凭据模式和 50 MiB 文件阈值 | 未发现误入产物、匹配凭据模式或超阈值文件 |

上述为五项不同 CTest 的两个执行窗口，未执行全部 96 项已注册测试或完整物理 benchmark。legacy_nodal_state_contract、pic_checkpoint_contract 因原始 fixture 缺失保持 Disabled，不计作通过。正式数值能力仍按各任务 gate 验收。

全新编译记录已有源码中的 C4244、C4267 类型转换警告，以及 CICMomentExactCut.cpp 的 C4828 字符编码警告；本次未改动这些生产文件。最终依赖选择为既有 Eigen 3.4.1 配置包，SuperLU/CBLAS 来自 D:/Chen/openLib，示例 EIGEN3_INCLUDE_DIR 的 3.3.9 仅为 fallback。编译警告保留在 build.log。

首次提交脚本会逐文件比对工作树原始字节对应的 Git blob SHA 与索引，发现意外换行或过滤器转换即停止；本机 candidate_manifest.json 保存候选 SHA-256 清单。

首轮 architecture 失败：`Forbidden include directory: quantum_hermite_geometry ->`。实际依赖为 C: 上的 Eigen 配置包；C++ filesystem 无法将跨盘符路径表示为项目相对路径，旧检查器把空相对路径当作内部目录。修复先规范化实际路径，再把不同盘符的目录归为外部依赖；同盘符内部层级检查保留，原所有拒绝自测继续运行，不改数值容限。

## 历史证据缺口

部分既有报告引用 `runs/baselines/...`，当前工作区没有该目录。本次不改写历史结论或伪造原始结果。旧证据的定位、恢复和独立备份仍需后续完成；源码和报告快照不等于完整历史计算数据归档。

## 分支约定

main 为唯一长期集成分支，任务使用 feat、fix、chore、docs 或 experiment 前缀，按任务创建并通过 Pull Request 合并。首次上传同时建立 `chore/reproducible-build` 作为后续工程工作的入口。使用说明见 [GitHub 与 VS Code 开发流程](../guides/Git_Workflow.md)。
