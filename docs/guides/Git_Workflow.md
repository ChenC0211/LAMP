# LAMP GitHub 与 VS Code 开发流程

本地工作区为 `D:\Chen\LAMP`。Git 管理源码、静态输入、参考数据、测试和文档；构建产物、日常计算结果及本机依赖配置留在 D:，另行备份。首次快照及其验收范围见 [初始化报告](../task/GitHub_Bootstrap_Report.md)。

远程私有仓库为 [ChenC0211/LAMP](https://github.com/ChenC0211/LAMP)，HTTPS 地址为 `https://github.com/ChenC0211/LAMP.git`。

## 在新电脑构建

安装 Git、VS Code、Visual Studio 2022 的 C++ 桌面开发工具和 CMake 3.20 或更高版本。项目使用 C++17、Eigen、SuperLU 和 BLAS；参考依赖组合为 Eigen 3.3.9 与现有 SuperLU/CBLAS Release 库。SuperLU 和 BLAS 应匹配 MSVC x64 ABI；OpenMP 构建需要相应运行库。VS Code 建议安装 Microsoft C/C++ 和 GitHub Pull Requests and Issues 扩展。

将仓库克隆到 D:，在根目录创建本机配置：

```powershell
Set-Location D:\Chen\LAMP
Copy-Item docs/guides/CMakeUserPresets.example.json CMakeUserPresets.json
```

按实际安装目录修改示例中的四个依赖路径。CMake 优先使用已安装的 Eigen3 配置包，只有未找到配置包时才使用 EIGEN3_INCLUDE_DIR；需要固定已安装包时设置 Eigen3_DIR，并检查 CMakeCache.txt 的最终依赖位置。初始化机器实际使用既有 Eigen 3.4.1 配置包，示例的 3.3.9 路径是 fallback，不表示本次构建已换用该版本。依赖路径不进入共享 CMake 配置，不上传第三方库或 `.venv`。示例为 local、local-debug 和 local-benchmarks 使用独立构建目录，避免旧 C: 工作区缓存或 benchmark 配置影响正常构建。

```powershell
.\tools\build.cmd
.\tools\build.cmd -Test
.\tools\run_case.cmd -List
.\tools\run_case.cmd hermite_gals_rotation_xy
```

包装器会寻找 VS 2022 附带的 CMake。完整研究驱动需显式执行 `tools/build.cmd -Preset local-benchmarks`，默认 CTest 不代表完整物理验收。依赖下载、临时文件和缓存使用 D:；按需将进程级 TEMP/TMP 指向 `D:\Chen\Codex\temp`。

## 分支和提交

`main` 是唯一长期集成分支。每项任务从最新 main 创建短期分支：`feat/<任务>`、`fix/<任务>`、`chore/<任务>` 或 `docs/<任务>`。未采用的数值研究用 `experiment/<任务>`，报告原始结果和失败原因；不得混入已验收实现。分支按任务划分，不按物理模块永久拆分。

```powershell
git switch main
git pull --ff-only
git switch -c fix/example-task
# 修改后逐文件检查并暂存本次任务内容
git diff
git add -p
git diff --cached
git commit -m "fix: describe the corrected behavior"
git push -u origin fix/example-task
```

`git add -p` 按改动块暂存已有文件；新增文件需明确执行 `git add 实际文件路径`。提交消息使用 `feat:`、`fix:`、`chore:`、`docs:`、`test:` 等前缀，说明实际行为。一项任务可以同时修改源码、测试、算例和说明。首次完整导入按现状保存，不补造过去的开发历史。

在 GitHub 创建目标为 main 的 Pull Request，填写改动、验证范围、已知限制和证据位置。配置/构建、architecture、repository_layout 和受影响的短契约通过后合并；算法或误差预算变化时追加对应 benchmark，禁止放宽物理断言来消除失败。推荐 Squash and merge；本地回到 main 后执行 `git pull --ff-only`。确认已合入后再清理任务分支。

## VS Code 操作

`Ctrl+Shift+G` 打开源代码管理，逐文件查看 diff 和暂存；左下角分支名称用于创建或切换分支。Commit 保存本地历史，Push 上传到 GitHub。切换任务前提交或妥善保存未提交改动。

`Ctrl+Shift+B` 执行现有 LAMP 构建任务。执行 `LAMP: prepare debug case` 后，F5 选择 LAMP 调试配置并粘贴新运行目录。详细命令见 [运行与调试](Run_And_Debug.md)。

## 数据和验收版本

`src`、`include`、`tools`、`tests`、`cases`、`docs` 和共享 `.vscode` 配置提交 Git。`build`、`bin`、日常 `runs`、算例 output、`.venv`、CMakeUserPresets 和凭据不提交。科研参考 txt/csv、冻结方法 PDF 和包含唯一历史版本的 ZIP 保留；不按扩展名笼统排除。

首次导入的 `.gitattributes` 保留现有文件字节，只为新增仓库管理文件规定 LF。任何以后统一换行的任务必须单独审核历史哈希影响。

`snapshot-20261009` 记录首次快照，不代表完整物理验收。正式验收或论文计算应在干净工作目录建立明确标签，并记录 commit/tag、依赖与构建版本、输入哈希、执行命令、结果归档路径和哈希、通过范围及限制。使用未提交改动时同时归档补丁。

GitHub 不保存被忽略的计算结果。冻结证据继续遵守 `runs/baselines/<任务>` 的项目约定，同时备份到 D: 上独立归档或选定的私有存储，并在报告中记录定位信息。历史报告的数字、结论和哈希不改写；缺失证据只记录缺失，不改判历史结果。
