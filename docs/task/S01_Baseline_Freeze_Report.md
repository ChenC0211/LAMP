# S01 / 基线冻结交付报告

**结论：PASS（基线冻结任务）。** 执行日期：2026-09-30。任务定义见 [TODO](S01_Baseline_Freeze_TODO.md)。

当前源码、原始程序、关键算例与历史证据已归档，并完成内容哈希校验、Release 构建及短验证。`hydro_moment_xy` 的现有数值失败已复现并保存；S01 的 PASS 只表示基线可追溯，不表示 Hydro 或任何 V2 数值 gate 通过。本任务未进入 S07/S08/I01。

## 1. 归档与身份

归档目录：[runs/baselines/S01_20260930](../../runs/baselines/S01_20260930/)。

| 对象 | SHA-256 |
| --- | --- |
| `baseline.zip`（29,143,503 bytes） | `b1bf9258e46db2ff95fcf48d3375cdadab408c20e5caa5f2bd3ca0ca90d5dab1` |
| `inventory.csv` | `a0b797f92cff8c162135ee3c205feb4bf2a549ad51e8cc87d2d30a1e75e4c4a7` |
| 原始 `bin/QUANTUM.exe` | `7ef81fea0bd2ae43c64ebd2922f8b071871935f81be1a433342c7b16891f4441` |
| 当前源码重建 `QUANTUM.exe` | `39da03ca01c8354bcafbc828b554ad1151dba926becfabdaac056fc3eed7d1ea` |

`baseline.zip` 在构建前生成，共 **1,687 个文件**：源码/测试/构建入口 960、原始运行文件 2、原始测试运行文件 4、构建证据 5、case 输入/说明 127、历史证据 44、规范/冻结方法文档 9、初始任务计划 1、依赖文件 501、工具 34。`inventory.csv` 逐条记录 archive 路径、原始路径、长度、时间及 SHA-256；`manifest.json` 记录范围与排除项。

原始可执行程序与 DLL 另解出至 `original_runtime/`，重建程序保存于 `rebuilt_runtime/`。原始程序是否由冻结源码完整构建无法独立证明，因此两者分别保留，不能混称为同一二进制。重建已更新工作区 `bin/QUANTUM.exe`，原始程序仍在快照中。

未复制全部历史输出、构建对象/PDB、旧 `QUANTUM_.7z` 和失效 Python 虚拟环境。后生成的验证证据存入 `evidence.zip`，其哈希及关键文件哈希见 `delivery_manifest.json`。本地归档不等于异地备份。

## 2. 构建环境与 T0

- Windows x64，Visual Studio 17 2022 generator，MSVC **19.40.33812.0**，C++17。
- `local` preset，Release；缓存记录 `/O2 /Ob2 /DNDEBUG`。原有特定源文件 `/fp:precise;/GL-` 等例外由归档的 `CMakeLists.txt` 保留。
- Eigen **3.3.9**；SuperLU **5.2.2**；配套 `CBLAS_Rel.lib`；原始 OpenMP DLL 已归档。依赖头文件/静态库归档在 `external/`，移机恢复时需在用户 preset 中显式重定向路径。
- Git commit：`unavailable`。工作区没有 `.git` 元数据，PATH 中也无 Git，不推测提交号。
- 执行 `powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/build.ps1 -Target QUANTUM`：退出 **0**，墙钟 **20.669 s**；日志 `build.log`、记录 `build_result.json`。
- 三个 contract target 构建成功，记录 `contract_build_retry.log` 与 `contract_build_retry_result.json`。首次经 `powershell.exe -File` 传逗号分隔 target 被当成单个名称，构建命令失败；已改用 CMake 独立参数，保留首次日志，未修改项目构建脚本。
- 构建含既有 SuperLU `/LTCG` 链接提示；本次没有修改 production 源码或引入编译器告警。

## 3. T1 与内容完整性

| 现有检查 | 结果 | CTest 墙钟 |
| --- | --- | --- |
| `quantum_help` | PASS | 0.179 s |
| `interface_time_flux_contract` | PASS | 0.509 s |
| `hydro_gradient_load_contract` | PASS | 0.100 s |
| `electricity_polarization_contract` | PASS | 0.120 s |

4/4 通过。完整命令与输出在 `ctest.log`、`ctest.xml`；这些检查针对现有实现，不替代未来 I/H/E 的验收。

归档 1,687 个条目内容与清单逐一匹配，1,588 个源码/输入/依赖文件与冻结时完全一致。六次运行的 20 个输入文件副本均与归档一致，见 `measurement_input_verification.csv`。未改 `include/`、`src/`、现有测试、原始 case 或容限。

## 4. T2：同输入、两种程序的短运行

运行位于 `measurements/{original,rebuilt}/<case>/`，逐项保存输入、命令、程序哈希、日志、退出码、墙钟及输出。均使用 `OMP_NUM_THREADS=2`、`OMP_THREAD_LIMIT=2`、`OMP_DYNAMIC=FALSE`、`OPENBLAS_NUM_THREADS=1`、`MKL_NUM_THREADS=1`，未缩短原始输入终时。

| 算例 / 输入来源 | 原始程序 | 重建程序 | 观察 |
| --- | --- | --- | --- |
| `cases/interface_geometry_b0` | 退出 0；5.925 s | 退出 0；5.900 s | `cassiniRZ`，resolution=16；7 regular frames、2 expected singular frames、split/merge 各 1；最大相对体积误差 `5.722619618338822e-16` |
| `cases/hydro_moment_xy` | 退出 3；0.900 s | 退出 3；0.873 s | 首个时间步失败，accepted=0、output=0；未完成原定 40 步 |
| `runs/checks/electricity_uniform/20260919_131755_7d0bf84b/main/system` | 退出 0；0.148 s | 退出 0；0.165 s | 原 `sharpHermiteP3` 输入，3 步至 `t=0.003`；残差 `5.6563416767817254e-16` |

Interface 结果明确标记 `manufactured_geometry_only`，曲率等不可用指标是 `null`；不能称为 Hermite–GALS 输运/曲率收敛验收。Electricity 的 3 步固定界面运行不能称为 SI-HHO-Q 验收。

Hydro 两种程序给出相同核心错误：

```text
moment-advance failed: parametric flux primitive is not on the bound native half charts; accepted=0, output=0
```

此项记为**现有 baseline 数值失败**。S01 不修复该数学路径、不加 fallback、不更改容限；后续应在相应路线图任务中评估其影响。

计时均为单次测量，包含启动与 IO；并非隔离性能实验。尤其不足以从两种程序的短时间差宣称加速或退化。对应源码没有改动。

## 5. 历史证据及限制

归档保留了以下原始路径中的输入、日志、series/meta 与退出记录（排除 fields/runtime/cache）：

| 历史路径 | 冻结事实 / 限制 |
| --- | --- |
| `runs/interface_parametric/20260927_044807_conforming_final16` | 旧 `parametricGals` 的 128 步至 t=2；summary 的 `g1_accepted=false`；旧几何表示与时间推进，不是 V2 Gate I |
| `runs/hydro_moment/20260927_xy_refined_final` | 当时日志完成 40 步至 t=0.4；当前输入/当前程序已不同，不能用该记录掩盖本次 Hydro 失败 |
| `runs/checks/electricity_uniform/20260919_131755_7d0bf84b/main` | 独立电学输入与历史输出；本次完整复用其输入 |
| `cases/hydro_drop_oscillation/output/20260929_195614_382_Perf_whole_step_0a024f` | 4 线程；总运行墙钟 1222.716 s；手动停止，退出 -1；只接受至 t=0.001，完整基准未完成。该总时长包含未接受候选，不能作为每步耗时 |

最后一项的程序 SHA-256 为 `12ce4764ca73e802e821cfd007477ae0920b9b99f1b1698e5d14e1c884b11eed`，与本次原始/重建程序均不同。它的最后接受态方程残差 `1.5463529340093175e-11`、液量漂移 `6.9388939039072284e-18`。停止原因与原始 `run.json` 的 failed 状态一并保留。

另归档 `docs/history_todo/QUANTUM_TimeLoop_Performance_Diagnosis_20260929.md`，仅为二级报告证据。该报告引用的 `cases/module_history/...` 原始诊断目录当前不存在；其中阶段计时不能冒充本次实测或独立核验结果。

## 6. 复核与恢复

在项目根目录执行：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/freeze_s01_baseline.ps1 -Destination C:/Chen/QUANTUM/runs/baselines/S01_20260930 -VerifyOnly
```

该命令校验 archive/inventory 哈希、全部归档内容以及当前源码/输入/依赖是否仍等于 S01。后续合法任务改动源码后，源码比较会明确失败；这不等于原归档内容损坏。工具拒绝覆盖已有基线目录。

恢复请解压到**新目录**，先核对 manifest 中的哈希；原运行入口为 `original_runtime/QUANTUM.exe`，本次重建入口为 `rebuilt_runtime/QUANTUM.exe`。用 `--case <独立运行副本>` 指定所需输入；不要直接覆盖历史输出。

## 7. 交付、预算及后续门槛

- 新增 `tools/freeze_s01_baseline.ps1`、本报告与 TODO；production 文件 **0**，测试文件 **0**。
- 生成 `runs/baselines/S01_20260930/` 的归档、清单、运行与验证证据；既有 `build/` 和 `bin/` 经正常构建更新。
- 本任务交付 S01；**Gate S 尚未完成**。路线图虽声明 S02–S07 文档包完成，当前 `docs/README.md` 缺失，旧引用路径及历史文档顶部标记仍需后续核对。本任务不替其宣告通过。
- 下一项按路线图末尾顺序核验/处理 S07，然后 S08；只有 S01–S08 的证据齐备才进入 I01。
