# I12a11 / 非多项式完整平移空间收敛报告

**FAIL：12 组中 11 组完成 t=.04；XY AMR 最细档库存拒绝，另有 XY uniform 曲率不单调。** 0 production + 1 test，遵守 [TODO](I12a11_Nonpolynomial_Translation_Study_TODO.md)。不改变原门槛；I12/Gate I 仍未通过。

## 完整推进结果

新增 opt-in `--translation-space-study`，复用 advanceTranslation。初始 profile 为既有 `exp(q)-1`，之后不解析重置；XY/RZ × uniform/AMR，n=8/16/32，固定 dt=.01、4 步。runCaseUsing 仅返回其已有的最终 Errors，失败返回空值，原行为和输出不变。

下表为实际完成 t=.04 后的最大曲率误差，以及 16→32 的梯度/曲率阶。缺失的最终数据不使用 t=.03 的 accepted error 补充。

| 配置 | n=8 曲率误差 | n=16 曲率误差 | n=32 曲率误差 | 末档梯度阶 / 曲率阶 | 结论 |
|---|---:|---:|---:|---:|---|
| XY uniform | 7.2089e-4 | 1.0019e-3 | 7.0990e-5 | 4.6191 / 3.8189 | FAIL，8→16 曲率增长 |
| XY AMR | 1.0021e-3 | 7.0990e-5 | 未完成 | 未计算 | FAIL，库存拒绝 |
| RZ uniform | 3.7214e-3 | 4.5696e-4 | 3.4170e-5 | 3.4300 / 3.7413 | PASS |
| RZ AMR | 3.7212e-3 | 4.5697e-4 | 3.4170e-5 | 3.6258 / 3.7413 | PASS |

XY uniform 的末档阶数虽足够，但 8→16 的曲率阶为 -0.4748，不满足本任务及 I12 要求的误差递减。不能删除该档或只引用末档数值宣布通过。两组 RZ 的梯度/曲率均逐档递减。11 个完成案例的实际 ledger/reflux 输出均为零；全局相体积相对解析解的误差仍含初始 Q4 近似误差，并不因此为零。完整 position/gradient/curvature/global volume/cell-M 数据保留在原日志。

XY AMR n=32 已接受 3 步，到 t=.03 后的下一半区间在 `half_integrated_flux` 拒绝：

```text
[Interface/HermiteGalsPhaseFlux geometry=150] advance: inventory violation exceeds roundoff budget cell=1472
```

accepted identity/M 隔离断言通过。t=.03 的位置、梯度、曲率误差分别为 2.2172e-10、3.3381e-8、5.9350e-6；这些是失败前状态，不作为 t=.04 收敛数据。本项 `failures=2` 表示一个推进失败和一个网格组的曲率递减失败，不是两个推进案例失败。

## 证据与后续

整个 12 组矩阵运行 1108.00 秒、退出 1。Release 构建 3.44 秒、无 warning；11 项直接回归全部通过，33.98 秒，完整输出与 I12a10r2 逐字一致。所有 production、其它 27 个可执行、81 冻结/参考文件、3 PDF 和前序 archive 哈希不变；旧 midpoint crossing probe 仍退出 1，输出相同。源码/runtime/日志及逐项审计归档于 `runs/baselines/I12a11_20261001/`。证据审计 PASS 不改变数值 FAIL。

下一独立诊断优先重放 XY AMR n=32 的库存拒绝，将原始 M 更新与同一旧 scalar 的 cell volume、平移后区域量逐项分解。已读代码中 I10 默认闭合容限为 1e-10，库存舍入预算为 64 eps；这只是待核验的尺度差异，尚不能认定为本次失败原因。曲率问题随后比较预测 scalar 与实际 bubble/publication，不能把全部失败都归因于时间积分或某一个容限。保持现有生产路径与原 gate，继续推进。
