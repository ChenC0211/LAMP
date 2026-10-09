# I12a16r7 / 认证通量舍入修复结果

**PASS（本任务）；原完整 I12 / Gate I 仍 FAIL。** 原三个第二步 full 证书越界已消除，四个原始 rotation/vortex uniform/AMR 均完成第二步。随后均在第三步 half 来源认证失败，完整 30/400 步没有完成。

按[设计](I12_Certified_Roundoff_Repair_Design.md)，原唯一库存核函数接收可选严格面界限。增长仅使用 upper 以下余量，减小仅使用 lower 以上余量；原 64epsilon、128 次上限、唯一共享面流量、total/boundary 固定及严格库存界保持。没有改动 M 或要求几何与 M 机器精度相等。无证书入口保持旧算式。

验证：40 个原修复会越界的缩减反例全部通过新修复，覆盖空/满、正/负、非零 lower、小体积尺度和极小流量；另有 4 个合法增长和 8 个拒绝检查。39 targets 编译通过、无 warning；26 CTest 全通过（309.69 s）。25 项旧完整输出相同，库存测试旧完整前缀相同，原 smoke 与 legacy crossing 输出相同。81 冻结文件、3 PDF 和 r6 archive 保持。生产 3 文件、测试 1 文件。

原 full（31.45 s）16 个成功 half/full 区间的 ledger/reflux residual 均为零；第二步 8 次认证选择均合法。第三步失败如下：

| 算例 | 已接受时间 | 失败面 | 原因 |
| --- | --- | --- | --- |
| rotation uniform | 0.02 | 44 | depth 40，initial boundary hit 无法排除 |
| rotation AMR | 0.02 | 155 | 4096 panels，midpoint domain unresolved |
| vortex uniform | 0.002 | 44 | depth 40，initial boundary hit 无法排除 |
| vortex AMR | 0.002 | 29 | depth 40，initial boundary hit 无法排除 |

测试调用记录保留一次 PowerShell 数组展开错误：39 targets 被合并为一个 MSBuild 名称，构建直接失败，后续测试没有执行。修正参数后重新完整编译并验证；没有将旧可执行计入证据。

用户随后授权：当前路线很难实现时转为“纯 Hermite–GALS，不守恒”。因此本项归档后停止继续扩展守恒通量认证，下一项先明确新路线的状态/几何/体积诊断契约，再实现纯 first-jet 输运与验证。此次局部 PASS 不作为新路线或完整 I12 通过的证据。

证据：[审计](../../runs/baselines/I12a16r7_20261002/verification.json)、[原驱动结果](../../runs/baselines/I12a16r7_20261002/driver_verification.json)、[源码范围](../../runs/baselines/I12a16r7_20261002/source_audit.json)、[回归比较](../../runs/baselines/I12a16r7_20261002/regression_comparison.json)。
