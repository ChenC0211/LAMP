# I12a3 / integrated face volume 共享库存装配

## 任务 ID / 标题

I12a3：将同一区间的实际面积分量装配为严格有界的守恒库存。I12a2 已 PASS 并归档；用户已明确要求持续推进直至 I12 通过。

## 目标

只解决一个事实：M 的唯一更新必须直接消费共享 integrated liquid volume，不能再从一个 midpoint snapshot 或平均 rate 重新积分。

## 允许修改

最多 4 个 production 文件：HermiteGalsTemporalPhaseFlux.h/.cpp、HermiteGalsPhaseFlux.cpp、新增私有 Detail/HermiteGalsInventoryAssembly.h。最多 2 个测试相关文件：新增 integrated inventory contract 与 tests/CMakeLists.txt。本 TODO、API/报告及 runs/baselines/I12a3_20261001 证据。

## 禁止修改

不接 driver、I09 reflux、bubble、general quadrature、first-jet、Hydro/Electricity、字典或输出。不改变已有 64*eps 库存舍入预算、128 次修复上限、CFL/散度/守恒门槛。旧瞬时算子保留显式 A/B，不成为新装配失败的 fallback。

## 数学约束

完整且有序的 actual faces 必须来自同一个 I12a1 interval；区间 begin 必须为旧 accepted 实际时刻。每个面的四个非负 gross integrated volumes 及唯一 signed liquid volume 被相邻两侧反号使用。库存只按共享面交换更新，不由几何重算，不 clipping。只允许原有预算内的内部共享 flux 舍入修复；原始与最终积分量均保留可追溯 owner，source/mesh/measure 不变。

## 实现要求

从旧 I08 提取同一个私有库存装配器，保持旧路径的运算顺序与结果；新路径以体积为单位调用，旧路径按原 dt*rate 单位调用。新增结果显式保留 interval，不借用 midpoint geometry identity。总 outflow CFL、散度、局部严格容量和整体 ledger 全部沿用既有检查。输入 integrated transfer 是代数契约，本装配器不冒称证明任意外部输入的求积/边界历史。

## 验收

- T0：Release、architecture、冻结文件及 4+2 文件预算。
- T1：新实际 I12a2 单面跨越量配合显式解析边界通量，验证 q=1/4、3/4 及相位扫描的原始 M，XY/RZ、正反方向、内部抵消、纯相、counterflow、同源区间组合；禁止把解析末态写入更新。
- T1：同源完整面/区间/mesh 校验；错误来源、遗漏/重复面、晚起始子区间、超 CFL、非零总散度、超预算库存失败，旧 authority/M 不变；可达预算内舍入通过原共享面修复。
- 回归：I12a2 的 19 项检查及旧 crossing 输出；新增装配未接旧 benchmark，因此旧失败仍显式保留。

## 失败规则

所有不兼容输入或超出原严格预算的库存缺陷明确拒绝。若需要新的数学事实或超出文件预算，先拆下一任务，不降低 gate。

## 交付

共享装配实现、API、独立测试、原始日志/逐输出比较/SHA256 归档。验收后继续 I09 的区间积分消费等必要后续任务，直至 I12 的真实验收通过。

## 执行结果

**FAIL。** 公共装配完成，19 项前序回归通过且 18 项完整输出一致；新 RZ 跨越矩阵在预算内一 ulp 超容处暴露旧单向舍入修复缺陷。独立代数检查通过不能替代此失败。详见 [报告](I12a3_Integrated_Inventory_Assembly_Report.md)。先拆出内部有界舍入修复任务，再重新完成 I12a3 验收。

**后续验收 PASS：** 独立 [I12a3r1](I12a3r1_Bidirectional_Roundoff_Repair_Report.md) 完成原预算内的内部增/减修复，148 组库存检查、20 类拒绝和 20 项 CTest 通过。初次失败归档保留，当前可继续 I09 区间消费。
