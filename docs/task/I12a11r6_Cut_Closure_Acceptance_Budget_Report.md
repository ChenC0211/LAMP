# I12a11r6 / cut 闭合验收预算与求根启动条件报告

**FAIL：单位与 30 项回归、原四步案例和长平移通过，但空间矩阵只有 10/12 组完成，I12 仍 FAIL。** 对应 [TODO](I12a11r6_Cut_Closure_Acceptance_Budget_TODO.md)。

## 实现与单位验证

2 production + 2 tests。所有单元的验收预算为 `min(relativeVolumeTolerance,64 eps)`；预算内 cut cell 不再启动 root。只有超预算 cut cell 使用原证书、括区及 safeguarded root 解到 `min(relativeVolumeTolerance,16 eps)`，最终核验按实际是否激活 root 区分两种精度。新增 `relativeAcceptedVolumeTolerance` 报告字段，原纯相字段仍保留。M、first jets/B/J、通量和全部几何/求积/迭代限制不变，详见 [API](../api/I10_Hermite_Gals_Bubble_API.md)。

首轮 Release 构建 9.61 秒，无 warning；bubble 单位测试 75.01 秒 PASS。6 组新增 XY/RZ cut 对照确认：32 eps 默认预算内不校正，128 eps 触发 root；32 eps 配合用户 8 eps 上限也触发 root。后两种的最终残差达到实际工作精度，M/系数及 identity 的要求按情形验证。原 23 条 bubble 数值行与 r4 完全一致，Newton/bisection/endpoint 仍为 503/128/0。

## 原四步结果与精度限制

同一 XY AMR n=32、dt=.01、非多项式、四步配置，经 480.59 秒到达 t=.04。旧库存越界单元、拒绝和诊断异常均为零，ledger/reflux 为零；原 first-step pure/cut 拒绝及原 fourth-step cell 1472 库存拒绝均未重现。

首步 midpoint 曲率误差保持 8.4588e-6；end 从 r4 的 3.2321e-3 降为 6.3637e-6，对 predictor 的放大从 621.45 倍降至 1.2236 倍。独立 bubble value/gradient/Hessian 检查通过。

但最终 t=.04 的位置/梯度/曲率误差为 2.1667e-9、1.8610e-6、1.1505e-3，不能将推进完成或守恒零残差当作曲率收敛证明。此外，首步 end 的 cell 1408 原闭合体积与局部 Green 测量相差 8.7007e-18，约为库存预算的 2.51 倍；这未触发本次库存越界，但说明工作精度分配仍不构成所有同源测量误差的严格控制证明。后续需结合空间矩阵定位。

## 回归与完整空间矩阵

完整 Release 重建 31.89 秒，无 warning；30/30 CTest 通过，205.01 秒。28 项输出逐字不变。bubble 仅新增上述预算场景，原数值行不变；translation smoke 仅 XY uniform/AMR 两行几何误差改变，其余所有列与行保持。旧 midpoint crossing probe 仍退出 1，原错误不伪装为通过。

四组静止各 4 步和 XY/RZ × uniform/AMR 四组平移各 40 步全部完成，210.27 秒、退出 0。四组平移的最大 ledger 残差为 1.3878e-17，reflux 为零，几何误差完整保留于日志。

12 组空间矩阵完成运行，耗时 1230.44 秒、退出 1；10 组推进到 t=.04。XY uniform n=32 和 XY AMR n=16 在接受 t=.03 后的下一 midpoint closure 出现纯相拒绝，分别为 cell 531 和 341。XY uniform 的 n=8→16 梯度/曲率误差增长；两组 XY 因未完成全部级别，程序记为 INCOMPLETE，而非另计一次 order failure。不能把统计中的 order_failures=0 解读为收敛通过。

| 配置 | n=8 曲率误差 | n=16 曲率误差 | n=32 曲率误差 | 最后一级梯度/曲率阶 |
|---|---:|---:|---:|---:|
| XY uniform | 7.2089e-4 | 2.9548e-3 | 拒绝 | 未完成 |
| XY AMR | 2.9627e-3 | 拒绝 | 1.1505e-3 | 未完成 |
| RZ uniform | 1.0264e-2 | 9.8473e-4 | 8.1810e-5 | 5.0329 / 3.5894 |
| RZ AMR | 1.0263e-2 | 9.8472e-4 | 8.1810e-5 | 5.9058 / 3.5894 |

RZ 两组均满足本矩阵的逐级下降和阶数要求。XY AMR n=32 的独立四步诊断与矩阵最终 BENCHMARK 行逐字一致。下一项只读定位两个新拒绝的差额来源，再决定修复；不因旧失败消失而放松纯相或收敛门槛。

## 冻结及交付

源码审计确认仅 2 production + 2 tests 改动、无新增源码；81 项冻结/参考文件、3 PDF 与 r5 archive 哈希不变。架构仍为 12 targets、380 translation units、900 files、3088 includes、5 legacy edges。完整日志、逐项输出对比、源码/runtime 哈希、交付清单及可独立验证的 archive 位于 `runs/baselines/I12a11r6_20261001/`。本项数值结果 FAIL 与归档完整性 PASS 分开记录。
