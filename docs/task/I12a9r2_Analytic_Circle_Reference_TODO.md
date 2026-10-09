# I12a9r2 / 解析圆参考的端点消减

## 目标：一个数学事实

独立误差参考必须稳定地计算同一解析圆/球的单元交集。I12a9r1 的 XY uniform 已接受 40 步，但 circleCell 中的 asin primitive 在圆端点产生约 3.24e-11 消减误差，超过原 1e-12 总面积核对限，且 catch 再次计算同一参考使进程中止。只修复参考计算及其失败归因，不改变数值解。

## 允许与禁止修改

0 production；最多 2 测试相关文件：interface_hermite_gals_benchmark.cpp、tests/CMakeLists.txt。另含本 TODO/报告和 `runs/baselines/I12a9r2_20261001`。所有 production、B/J、时空相积分、库存、bubble、冻结路径及 PDF 禁止修改；不提高 1e-12 参考总面积检查或任何原 gate。

## 数学约束与最小实现

- 先加入独立 reference contract，在原 circleCell 上复现总面积错误并保存；再用同一半径/弦高的 atan2 替换病态 asin 原函数表达式。积分区域、切分位置与解析函数不变。
- 在 XY/RZ、uniform/AMR、多个网格尺度与平移/旋转中心上核对单元解析量总和、整圆/球、对称分区及局部可加性。测试不读取数值 phi，也不向独立 M 或 first jets 写值。
- runCase 在误差计算前设置真实 error_reference stage。失败后附加诊断即使失败也必须报告并让其它算例继续；原失败计数仍保留，不吞掉验收失败。
- 原完整 translation-sweep 40 步参数保持；所有四组都必须实际执行。原默认 driver 继续作为历史对照，不作为 fallback。

## 验收

- T1：新 `--circle-reference-contract` 修复前 FAIL、修复后 PASS；绝对总面积原 1e-12 限不变，记录最大误差。有效参数覆盖 endpoint cancellation 原反例。
- T2：原 I12a9 完整 4 组平移各 40 步及 4 组静止四步；保留全部 accepted time、位置/梯度/曲率/库存误差与后续失败。
- T0：benchmark Release 构建、新 reference contract、translation smoke、直接 Core root contract、architecture/CLI；证明所有 production hash、81 冻结/参考、3 PDF 和 I12a9r1 archive 不变。旧 crossing 应相同；原 curvature study 与 rotation-curvature audit 的受影响参考诊断实际重跑，结果按实际记录。

## 失败规则与交付

出现新的真实推进拒绝时按真实 stage/cell 保留结果并独立拆项。不能通过修改参考值去贴合数值库存，也不能忽略解析参考失败。交付最小测试修复、前后证据、日志/哈希归档及未解决问题；完整 I12 未通过前保持 FAIL。

## 执行结果

PASS：36 组独立参考、5 项直接 CTest、原 curvature/rotation 参考诊断通过；四组平移各完成 40 步、四组静止通过。所有 production 与 27 个未修改可执行 hash 不变。见 [报告](I12a9r2_Analytic_Circle_Reference_Report.md)。I12a9 常速度范围通过，I12 继续推进。
