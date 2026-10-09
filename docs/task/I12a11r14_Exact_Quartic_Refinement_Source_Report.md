# I12a11r14 / 精确 Q4 细化标量来源报告

**状态：PASS，2026-10-01。** 精确细化 T1 和合并审计后的 31 项回归通过。初次 bubble 的墙钟超时、同环境归档对照与最终复验分别保留。父任务 I12 / Gate I 仍为 FAIL。

## 交付的算子

新增 `Snapshot::refineNativeQuartic(targetMesh, source)`，允许 cut cell 的精确 Q4 限制。完整 active、同物理域、无 corner-sector；每个新叶必须属于一个旧父叶。不生成库存或 Interface publication，不接实际推进。

新 master DOF 保留 exact numerator/denominator，cardinal datum 直接保留已有精确表示；hanging/fine-side 继续使用原共享约束引擎。子叶通过原 Q4 限制内核共用父叶清分母因子。返回前对每叶 25 个 normalized 数据及对应 physical 共享约束作精确相等核验，未以浮点误差容限替代多项式证书。临时借用在构造结束时清除。

修改 4 个 production 文件、2 个 test-related 文件，均在 [TODO](I12a11r14_Exact_Quartic_Refinement_Source_TODO.md) 的范围内。接口边界见 [API](../api/I12a11r14_Exact_Quartic_Refinement_API.md)。旧 constructor 和已有 double transfer 保留原语义。

## T1 实测

56 个场景，XY/RZ × 正反相 × 有无非全局 center 扰动 × 7 种网格/来源组合。覆盖不变网格、均匀细化、uniform→AMR、AMR→更深 AMR、大坐标、非二进制坐标、已经精确细化来源的再次细化。每场景还重复构造并销毁旧来源及两份外部网格，再检查查询值。

共检查 4,832 个目标叶、43,488 个内部 jet 样点和 5,120 个曲线求积节点。六类输入拒绝通过：coarsening/跨父叶、域变化、partial-active、非 Q4、corner-sector、moved-from；旧来源仍可用。

| 浮点核验 | 最大误差 | 原容限 |
| --- | ---: | ---: |
| value / gradient / Hessian，除以 `1+abs(reference)` | `3.6522241312976431e-13` | `1e-10` |
| 边界 trace value / gradient，绝对差 | `1.7763568394002505e-15` | `1e-12` |
| 子叶液量之和与父叶液量，除以父容量 | `1.4135798584282298e-15` | `1e-10` |
| 同一保留位置的子午曲率消费，绝对差 | `1.7053025658242404e-13` | `1e-9` |

这些浮点检查验证消费链路；整多项式相等由构造内的 exact identity 证明。此处体积是独立几何测量，还不是守恒 `M_K` 转移。

## T0 与执行记录

首次完整构建成功，耗时 `412.5500601 s`。初次 CTest 共 31 项，30 项通过，bubble 在 `180.02 s` 超时；新增测试 `37.12 s`，56 场景及 6 拒绝全部通过。完整运行记录为 `550.6049331 s`。除了架构计数和超时截断的 bubble 输出，28 项原输出与 r8 逐字一致。

同环境归档 smoke 程序 `50.0295932 s`，当前程序 `50.31 s`，完整输出一致。归档 bubble 完整运行 `245.7850977 s`，当前 bubble `249.06328 s`，两者全部断言通过，输出与历史记录逐字一致；当前/归档耗时比约 `1.0133`。归档程序也较历史运行慢约一倍，本对照未显示约两倍的代码性能回退。

依现有长测试的设置，只把 bubble 的 CTest 墙钟保护时间从 180 s 调整到 300 s；测试源码、样本、数值容限与求解工作预算未变。初次 30 项通过结果保留，仅复验受墙钟限制影响的 bubble。最终审计明确记录两次 CTest 的来源，不把拼接结果冒充一次原始运行。

最终 bubble CTest 通过，测试计时 `243.30 s`（命令 `243.3536873 s`）；配置构建 `2.5474922 s`。31 个不同测试、32 次 CTest 执行、1 次初始墙钟超时，最终每项都有通过记录。原 29 项数值输出逐字一致；architecture 保持 12 targets / 380 translation units / 901 files，仅 include 从 3092 增至 3093。legacy crossing 的预期失败 exit=1 及原输出保持。

最终构建无 warning。81 个冻结/参考文件、3 份 PDF 和 r13 archive 哈希全部不变。运行包记录 29 个 executable 与原 OpenMP DLL；没有新增数学 stencil、容差、解析重置或库存 clipping。

执行中保留了全部失败日志：首次编译的异常消息参数类型错误；首次多层网格违反 2:1；大坐标测试错误使用显示位置反演局部坐标；非二进制细化的冗余分母乘积及重复构造中的未平衡交叉乘积触发精确算术拒绝。前两项修正编译/fixture，第三项改用原有保留 chart，后两项以共用分母、cardinal 数据复用和精确二进制幂平衡修正表示，未删去尾项或提高数值容限。

## 未完成边界

本接口仅解决精确标量来源。独立库存的父子守恒转移、I11 受检 regrid publication、局部细化重试及 I12 原平移/旋转/vortex/AMR/收敛矩阵尚需继续。不能由本项 T1 通过宣称 I12 通过，也不能进入 Hydro。

证据目录：[I12a11r14](../../runs/baselines/I12a11r14_20261001/)；[验证摘要](../../runs/baselines/I12a11r14_20261001/verification.json)、[56 场景](../../runs/baselines/I12a11r14_20261001/refinement_comparison.json)、[逐测试输出与来源](../../runs/baselines/I12a11r14_20261001/regression_outputs.json)、[回归比较](../../runs/baselines/I12a11r14_20261001/regression_comparison.json)、[源码/冻结审计](../../runs/baselines/I12a11r14_20261001/source_audit.json)、[bubble A/B](../../runs/baselines/I12a11r14_20261001/bubble_performance_comparison.json)。同目录 `evidence.zip`、`delivery_manifest.json` 和 `archive_verification.json` 保存交付文件及逐项哈希。
