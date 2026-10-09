# I12a9 / 同源常速度完整驱动报告

**FAIL（两组 XY 长程平移）；完整驱动已接入。** RZ 均匀网格与 AMR 均完成原 dt=.01 的 40 步，到达 t≈.4；XY 两组均在 t=.04 后的 full_integrated_flux 拒绝。I12 保持 FAIL，继续修复已定位的曲线裁剪归属检查。

## 实现

遵守 [TODO](I12a9_Translation_Driver_TODO.md) 的 2 production + 2 test 文件预算。新增 constantArrivalView 保留同一 immutable velocity provider，允许实际半步 candidate/time。benchmark 完整半/全区间共享旧 scalar、速度、characteristic 与显式常气相 BC history；所有面经 I12a8/I12a3/I09，再由原 I10/I11 闭合发布。发布后的 M 与最终共享通量形成的 M 逐项相同。interior 不重置解析场，B/J、bubble、容限均未修改。详见 [API](../api/I12a9_Translation_Driver_API.md)。

## 实际推进证据

`--translation-sweep-smoke`：8 组静止四步/平移首步全部通过，用时 29.98 秒；实际半步身份/时间、constant B/J、provider 生命周期及 6 类拒绝通过。

`--translation-sweep`：197.25 秒，退出 1，保留完整 8 组结果。静止四组通过，最大 first-jet/M drift 为 2.23e-16 以下。平移结果如下（误差定义沿用原 harness）：

| 算例 | accepted time | 位置误差 | 梯度误差 | 曲率误差 | 结果 |
| --- | --- | --- | --- | --- | --- |
| XY uniform | .04 | 1.64e-16 | 2.21e-15 | 1.03e-13 | 第五步 full flux 拒绝 |
| XY AMR | .04 | 9.13e-17 | 1.96e-15 | 1.36e-13 | 第五步 full flux 拒绝 |
| RZ uniform | .40000000000000019 | 1.80e-16 | 1.93e-15 | 5.69e-14 | 40 步 PASS |
| RZ AMR | .40000000000000019 | 1.80e-16 | 1.93e-15 | 6.13e-14 | 40 步 PASS |

RZ 两组最大共享库存 ledger 为 1.39e-17，reflux residual 为零，最终逐单元相对库存误差分别约 1.89e-15 和 1.02e-14。XY 拒绝后逐项验证 accepted identity/M 未变。

两组 XY 原始错误均为 `MappedInterface GL8 root-enclosure ownership: free-curve point has different overlapping region owners`，发生在 chosen=1 的宽度约 3.8e-8 的微小自由子弧。GL8 代表点归属已经通过，随后的横向根包围端点归属检查失败。它们越过了旧中点通量的原失败步，但不能因此将本项完整平移验收标 PASS。下一项检验包围端点接触裁剪边界与真实曲线跨界的区别；不得跳过小弧、放宽库存或降低 closure 门槛。

主程序及 26 个测试目标 Release 构建通过，用时 40.74 秒，无 warning；**26/26 CTest 通过**，337.56 秒。25 项旧测试完整输出逐字不变，architecture 仍为 12 targets / 380 TUs / 899 files / 3085 includes / 5 legacy edges。81 冻结/参考、3 PDF 与 I12a8 archive 哈希核对通过。旧 crossing probe 仍退出 1 且输出不变；旧完整 driver 未在本项重跑。原始证据、source/runtime 哈希和逐项归档校验保存于 `runs/baselines/I12a9_20261001/`。证据审计 PASS 不改变本项长程数值验收 FAIL。

## 未完成

XY 40 步、非恒定速度 rotation/vortex、独立跨面完整发布与最终 I12 收敛/其它 gate 仍待通过。旧默认 driver 作为明确历史对照保留，不是新入口的 fallback。
