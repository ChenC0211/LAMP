# I12r3 / 带同一 signed-partial 证书的定向 bubble 区间

## 任务 ID / 标题

I12r3：修复 bubble 对称安全区间不能覆盖实际所需单侧安全修正的问题。前置 I12r2 审计归档已核验，上一 goal 轮为 progress。父任务 I12 仍 FAIL。

## 目标

仅改变安全区间的搜索方向：在原对称区间无法括住 M 时，按 M(c) 单调性选择增大或减小 c，并用原 signed-Cartesian-partial / 位移证书覆盖预测系数至新端点的整个区间。不是放宽 volume/displacement/topology 容限，也不把不同合同的 changeCertificate 直接替换进求解器。

## 允许修改

- Core `LevelSetHermiteGeometry.h/.cpp`、`Detail/LevelSetHermiteQuarticLeafData.h`：同一 center-interval 证书引擎增加显式方向，原默认对称入口保持行为。
- Interface `HermiteGalsBubble.h/.cpp`：目标方向的 certified bracket 扩区及诊断，最多共 5 production 文件。
- `tests/interface_hermite_gals_bubble_contract.cpp`、`tests/interface_hermite_gals_benchmark.cpp`：最多 2 tests。
- TODO/报告、bubble API、docs README、`runs/baselines/I12r3_20261001/`。

## 禁止修改

时间取样/phase-flux、I02–I09、I11 发布算法、Hydro/Electricity/App/reader、三份冻结方法 PDF、S08 frozen 源。禁止裁剪 M、增大 0.125-cell 默认位移、修改已有容限、fallback/KKT 或新增 production 状态。时间积分设计变更仍待明确决策。

## 数学约束

1. 固定 candidate skeleton 和 M，仅变动 interior center；bubble/gradient 在边界恒零。
2. 单侧区间由**同一个** exact Bernstein/outward subdivision 引擎证明：整个参数区间的可能零集上，存在同一个有符号 Cartesian partial 正下界，并证明 displacement<=原预算。
3. retainedDerivativeFraction 保持原搜索语义与证书判据；默认对称搜索的尝试次数保持。单侧倍增工作上限允许单独设置为有限的 16 次，所有新区间仍须通过同一证书；无自动增大物理预算。
4. M(c) 非增；扩区前验证证书，扩区端点通过既有体积/求积/topology gate，再允许 safeguarded Newton/二分。
5. 原对称区间已可达时，原算式、Newton 路径和结果保持；原 pure mismatch、证书不足、几何未解析和不可达 M 仍拒绝。
6. 一切失败不污染 candidate/accepted state；报告区分原区间与扩区后的证书。

## 实现要求

Core 提供默认 Both 与显式 Decrease/Increase 三种方向；不复制证明内核。I10 先保留既有对称 bracket 流程，仅对不可达的 active cut cell 尝试相应单侧证书。更新完整端点 bracket 和 root 步长记录，所有接受仍经相同 volume/topology checks；不能直接接受上一任务诊断探针的结果。

执行中新增证据：`range_cap_probe.log` 证明沿用 9 次倍增时，AMR rotation 的 directed upper 停在 0.17111038039062512，而先前外部根约为 0.37909467603515645；这是工作上限不足，不能以此结束修复。仅增加显式单侧搜索的有限工作量，不改对称默认或安全判据。

## 验收

- T0：Release 主程序及相关 targets、architecture、无新增 warning；5+2 文件预算、冻结 hash。
- T1：XY/RZ × uniform/AMR 的解析 Q2+bubble 制造目标；两方向证书、确实越过原区间、独立截面积分、边界 C1、失败隔离；原 I10/I11 及 I01–I09 回归。
- T2：重跑原 I12 moving benchmark，逐项比较 accepted time、失败阶段、position/gradient/curvature/volume；验证 I12r2 暴露的失败是否解除，不将更长运行当完整 gate。
- 控制：原 symmetric Core / legacy GALS 与 parametric state 对照；原 q=1/4、3/4 midpoint 跨面反例保持失败。
- 完整 I12 若仍失败，保留原始失败和未解决项；不通过 WILL_FAIL 或缩小目标变绿。

## 失败规则

若单侧仍不能证明同一 signed partial、位移或 target bracket，明确拒绝。若需要换成不同数学证书、改变时间公式或预算，停止该方向并另写设计说明，不在本任务擅自改方法。若需超出文件预算，拆任务。

## 交付

最小代码、方向合同与诊断、解析/实际运动验证、前后误差与失败对比、源码/runtime/归档哈希、报告及尚未完成的 I12/Gate I 项。总体目标保持 active。

## 执行结果

**PASS（本任务）；父任务 I12 仍 FAIL。** 5 production + 2 tests，15 项回归通过；最终新增失败隔离断言后单独重跑 I10 通过。8 个制造解配置覆盖 108 个实际扩区单元，独立体积误差最大 9.680157e-11，位移界最大 0.112497212 cell。4 个原旋转/涡旋阻断越过，但完整 benchmark 仍有 8 个运动和 2 个跨面失败，AMR rotation 曲率误差 32.942469，不能宣称精度或 Gate I 通过。详见 [交付报告](I12r3_Directed_Bubble_Interval_Report.md)。本轮不启动下一任务。
