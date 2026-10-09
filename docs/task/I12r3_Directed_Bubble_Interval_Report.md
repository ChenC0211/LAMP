# I12r3 / 带同一证书的定向 bubble 区间交付报告

**状态：PASS（定向区间搜索及其合同）；父任务 I12：FAIL；Gate I 未通过。** 2026-10-01。

本任务按 [TODO](I12r3_Directed_Bubble_Interval_TODO.md) 完成：当原对称区间不能覆盖体积闭合所需的单侧修正时，沿目标方向重新搜索，并仍由原 signed-Cartesian-partial / 位移证明引擎认证整个区间。I12r2 找到的 4 个原始运动阻断均已越过；后续仍有闭合或库存拒绝。更长的接受时间不等于精度、收敛或 production 验收通过。

## 改动与保持的关系

- Core 增加 `QuarticCenterDirection::{Both,Decrease,Increase}`。默认 Both 保留原对称区间及 9 次倍增上限；显式单侧区间以原预测 c 为一个端点，最多 16 次倍增，仍使用原 exact Bernstein / outward subdivision 证明内核。
- I10 先检查原对称区间和端点体积的单调性。只有目标在原括区外，才按 `M(c)` 非增关系选择方向，获得新证书并重新测量两端；仍执行原 moment、体积界限、topology、monotonicity 和 safeguarded root 检查。
- 所有证书从同一个原 predictor 出发，不重定位累计位移预算。默认 `maxDisplacementCells=0.125`、体积残差 `1e-10`、求积估计 `1e-11` 均不变。`retainedDerivativeFraction` 保留原搜索种子与减半尝试语义，不将其误称为最终导数下界的固定比例。
- 原区间已可达时保持原路径。失败不修改 M、candidate identity 或发布状态；纯相不匹配、未解析几何及不可证明的区间仍拒绝。没有把 I12r2 的另一种 `changeCertificate` 合同替换进 I10。
- 报告增加原对称边界、是否定向扩区和扩区单元数；benchmark 增加 `BUBBLE_EXPANSION`。该行表示一次候选闭合成功，不能单独证明整个时间步已接受。

最初沿用 9 次单侧倍增时，AMR rotation 的 upper 停在 `0.17111038039062512`，低于 I12r2 已审计根 `0.37909467603515645`。保留 [初始范围探针](../../runs/baselines/I12r3_20261001/range_cap_probe.log)，更新 TODO 后将显式单侧工作上限改为 16；每个试探仍须通过同一证书。这增加有限搜索工作，不增加允许误差或物理位移。搜索仍是充分且有限的，并不保证找到最大安全区间。

## T0 / T1 与历史控制

Release 构建成功；MSVC `19.40.33812.0`、CMake `3.28.3-msvc11`，local preset。OMP 线程数/上限为 2，OMP_DYNAMIC=false，BLAS/MKL 线程数为 1。构建日志未发现新增 compiler warning。

15 项 CTest 全通过，CTest 总计 92.53 秒；覆盖 architecture、help、历史 state 和 I01–I11 合同。随后只在 bubble 测试增加“定向求根迭代耗尽时不污染 candidate/M”断言，重新构建并单独重跑 I10，通过，65.89 秒。生产代码与其余测试自 15 项回归后未改变。最终审计采用此最后一次 I10 输出。

新增 XY/RZ × uniform/AMR × 正/反液相共 8 个制造解配置，108 个单元实际扩区，其中 c 减小 56、增大 52。每个最终根确实越过原对称边界；独立截面积分使用解析 Q2 加 bubble、独立一维根及 adaptive Simpson，不复用 snapshot 的零根。验证边界 C1 不变、M 不变、默认 Both 等于显式 Both、证书位移预算与非法枚举拒绝。

| 几何/网格（两种液相符号） | 扩区单元总数 | 最大独立体积误差 / capacity | 最大闭合残差 / capacity | 最大位移界 / cell |
| --- | ---: | ---: | ---: | ---: |
| XY uniform | 26 | 7.975225e-11 | 7.975232e-11 | 0.112098181 |
| XY AMR | 42 | 9.680157e-11 | 9.680257e-11 | 0.112497212 |
| RZ uniform | 18 | 1.033214e-11 | 1.033129e-11 | 0.112412583 |
| RZ AMR | 22 | 1.033214e-11 | 1.033129e-11 | 0.112412583 |

14 个未增加输出的 CTest 程序输出与 I12r1 逐字一致。原 I10 输出前缀及 `SAFEGUARD,388,128,2` 不变。legacy GALS 和 parametric state 原始日志哈希分别保持 `C6D8D43979F6A5998E9C54B292EA7EE38EA627A6352E66F9CB0BDBA779AAA5D1`、`0D1289C91325D6830310BBBE623182DA30911742DBE142B75FED4D855C23975A`。

legacy conservative GALS 对照通过，实际程序输出及已有 OMP Warning #96 与 I10 一致；原始日志只有 3 行 PowerShell 调用位置/命令上下文不同。审计只排除 `At line:` 和 `+ ` 开头的 shell 上下文，保留 OMP warning 及所有程序数值；不将此旧运行警告称为编译警告或新的数值失败。

`--gross-flux-contract` 通过。`--crossing-probe` 保持实际退出码 1，其 7 行跨面反例指标与 I12r1 完全一致，没有用 WILL_FAIL 把完整 benchmark 改成通过。

## T2 / 原运动 benchmark

最终完整运行耗时 238.185 秒，退出码 1：8 个 moving case 和 2 个跨面反例仍失败。4 个 static 输出和 12 个静态曲率输出与 I12r2 完全一致。以下对比保持网格、dt、目标步数及容限；旋转目标时间为 0.3，其余运动目标时间为 0.4。

| case | 原接受时间 | 本次接受时间 | 本次失败阶段/原因 |
| --- | ---: | ---: | --- |
| translation XY uniform | 0.02 | 0.02 | end_closure：目标仍不在已认证体积区间 |
| amr_crossing XY AMR | 0 | 0 | midpoint_closure：定向搜索不能扩展所需边界 |
| translation RZ uniform | 0 | 0 | midpoint_closure：纯液库存不匹配 |
| amr_crossing RZ AMR | 0 | 0 | midpoint_closure：纯液库存不匹配 |
| rotation XY uniform | 0.01 | 0.02 | midpoint_closure：后续定向区间不足 |
| vortex XY uniform | 0.031 | 0.034 | end_closure：后续定向区间不足 |
| rotation XY AMR | 0 | 0.01 | half_flux：下一步 cell 72 库存越界超过舍入预算 |
| vortex XY AMR | 0.034 | 0.036 | end_closure：后续定向区间不足 |

最终日志记录 10 次带扩区的候选闭合。单独 `--rotation-amr-first-step` 通过，验证原 AMR rotation 首步确实解除阻断；完整运行仍在下一步失败。

| 推进更远的 case（本次最后接受态） | position 误差 | gradient 相对误差 | curvature 绝对误差 | 全局体积误差 | 逐 cell M 解析误差 / capacity |
| --- | ---: | ---: | ---: | ---: | ---: |
| rotation uniform，t=0.02 | 3.277154e-4 | 5.773681e-2 | 3.723220 | 0 | 5.717406e-6 |
| vortex uniform，t=0.034 | 2.020259e-4 | 3.298195e-2 | 0.720495 | 6.938894e-18 | NaN |
| rotation AMR，t=0.01 | 6.309709e-6 | 1.656260e-2 | **32.942469** | 6.938894e-18 | 1.350289e-5 |
| vortex AMR，t=0.036 | 1.028448e-4 | 1.493175e-2 | 1.052049 | 2.081668e-17 | NaN |

AMR rotation 首步曲率误差很大，不能因位置误差小、体积闭合与位移证书通过而宣称几何精度合格。前后终止时间不同，表中进展也不能用作同一时刻的精度改善对比。完整前后误差见 [accepted_error_comparison.csv](../../runs/baselines/I12r3_20261001/accepted_error_comparison.csv)。vortex 尚无独立逐 cell 解析体积，保留 NaN。本任务未重跑时间/网格 refinement，不声称当前 moving T3 收敛通过。

## 源码范围与证据

修改 5 个 production 文件和 2 个测试，无新增源码文件；1578 个受审计文件中 1571 个不变。

- `include/Core/Geometry/LevelSetHermiteGeometry.h`、`src/Core/Geometry/LevelSetHermiteGeometry.cpp`、`include/Core/Geometry/Detail/LevelSetHermiteQuarticLeafData.h`：显式方向及同一证书搜索。
- `include/Modules/Interface/HermiteGalsBubble.h`、`src/Modules/Interface/HermiteGalsBubble.cpp`：条件扩区、原括区保留及报告。
- `tests/interface_hermite_gals_bubble_contract.cpp`、`tests/interface_hermite_gals_benchmark.cpp`：制造解/失败隔离和原运动比较。

69 个 S08 frozen 文件、12 个 integration reference 和 3 份方法 PDF 哈希不变。上次 I12r2 归档及其交付源码已在本次开始时校验；原归档不修改。没有 Git 元数据，不虚构 commit。

证据目录：[I12r3_20261001](../../runs/baselines/I12r3_20261001/)。其中保留初始与最终 build、初始不足搜索范围、完整失败日志；[最终 I10 日志](../../runs/baselines/I12r3_20261001/ctest_bubble_final_details.log)、[完整 benchmark](../../runs/baselines/I12r3_20261001/benchmark_final.log)、[源码审计](../../runs/baselines/I12r3_20261001/source_audit.json)、[回归比较](../../runs/baselines/I12r3_20261001/regression_comparison.json)、[数值汇总](../../runs/baselines/I12r3_20261001/numerical_summary.json) 支撑以上结论。

交付源码/文档和 18 个当前 runtime 一并打包。[verification.json](../../runs/baselines/I12r3_20261001/verification.json) 区分本任务 PASS 与父任务 FAIL；归档整体 SHA256 和逐 entry 验证见 [delivery_manifest.json](../../runs/baselines/I12r3_20261001/delivery_manifest.json) / [archive_verification.json](../../runs/baselines/I12r3_20261001/archive_verification.json)。后者的 PASS 只表示字节完整性，不代表物理验收。

## 未完成项与后续边界

本轮只完成 I12r3，下一任务仍从 I12 未解决的完整输运/闭合/精度问题中选取单一事实，不进入 H01。现有统一 signed partial 是充分条件，不能把证书未找到等同于物理解不存在；后续失败需按实际候选重新审计，不能直接外推 I12r2 对旧候选的结论。

[通量时间积分设计变更](I12_Phase_Flux_Time_Integration_Design_Change.md) 仍未批准、未实现。本任务不改变冻结的时间取样公式。两个解析跨面反例仍证实该缺口；其余后续失效不能未经核验一概归因于它。精确面内相切/重根、正长度零迹、未解析 closed subcell cycle、动态 regrid、完整 V2 配置/输出、移动曲率/体积收敛与强变形 A/B 仍未完成，低黏/无黏电喷雾目标保持 active。
