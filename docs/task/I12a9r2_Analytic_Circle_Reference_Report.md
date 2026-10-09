# I12a9r2 / 解析参考修复与完整平移报告

**PASS。I12a9 的常速度完整平移范围现已通过，I12 整体仍 FAIL。** XY/RZ、uniform/AMR 四组平移均按原 n=8、dt=.01 接受全部 40 步，四组静止四步也通过。没有改生产推进或验收门槛。

## 实现与独立证据

符合 [TODO](I12a9r2_Analytic_Circle_Reference_TODO.md) 的 0 production + 2 test 文件预算。circleCell 的同一解析 primitive 将病态 asin(d/r) 改为 atan2(d,s)，保持原弦高、积分区域和分段；runCase 在后处理前设置 error_reference，附加诊断失败也会记录并继续其它算例，失败计数不变。

修复前新增同一 reference contract 退出 2，整圆参考误差为 **-6.4844671299191248e-11**。修复后 36 组 XY/RZ、uniform/AMR、n=4/8/16、多个平移/旋转中心全部通过：最大总面积误差 **4.16334e-17**，最大局部切分误差 **1.38778e-17**；原绝对总面积 1e-12 门槛未改。参考函数只计算误差，从未用于重新初始化数值库存或 interior first jets。

## 完整推进

`--translation-sweep` 用时 231.88 秒，退出 0，8 组结果全部 PASS。平移实际终点均为 .40000000000000019：

| 算例 | 位置误差 | 梯度误差 | 曲率误差 | 最大逐单元相对库存误差 |
| --- | --- | --- | --- | --- |
| XY uniform | 2.01e-16 | 1.69e-15 | 7.90e-14 | 4.00e-15 |
| XY AMR | 1.10e-16 | 2.42e-15 | 1.40e-13 | 1.28e-14 |
| RZ uniform | 1.80e-16 | 1.93e-15 | 5.68e-14 | 1.88e-15 |
| RZ AMR | 1.80e-16 | 1.93e-15 | 6.13e-14 | 1.02e-14 |

最大步内库存 ledger 为 1.39e-17，reflux residual 为零。静止最大 drift 为 2.23e-16 以下。所有 half/full 源身份、独立 M、原 I10/I11 与失败事务规则沿用 I12a9；本项仅消除独立参考函数的误报。

## 回归与交付

benchmark Release 构建通过，3.42 秒，无 warning；5/5 直接 CTest 通过，17.31 秒。architecture、CLI、root-enclosure contract 三项完整输出与 I12a9r1 逐字相同；translation smoke 的唯一输出变化是解析参考 cell-volume 误差列，位置/梯度/曲率、数值体积、ledger/reflux、drift 和全部其它输出不变。

原三档静态 curvature study 退出 0，6.43 秒；rotation-curvature audit 退出 0，20.73 秒，仍是旧 midpoint 路线的诊断，不代表旋转完整推进通过。旧 crossing probe 退出 1 且输出逐字不变。所有 production hash、其它 27 个可执行、81 冻结/参考、3 PDF、I12a9r1 archive 均不变；没有为测试文件改动重复执行未受影响的完整生产回归。

原始日志、scope/runtime 哈希、输出列对照与逐项 archive 校验位于 `runs/baselines/I12a9r2_20261001/`。I12a9 和 I12a9r1 原始 FAIL 归档保留，后续常速度验收以本报告为准。

## 继续推进

下一步完成任意步内平面跨面的完整发布证据，并实现非恒定速度的同源时间通量。rotation/vortex、h/dt 收敛、零迹/拓扑、动态 regrid 及最终 Gate I 审计仍须实际完成；常速度 40 步通过不代替这些任务，不进入 Hydro。
