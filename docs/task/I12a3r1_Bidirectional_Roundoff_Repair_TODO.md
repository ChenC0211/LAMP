# I12a3r1 / 内部共享通量的双向舍入修复

## 任务 ID / 标题

I12a3r1：在原预算内使用内部面尚有余量的液相 gross 通量，修复 I12a3 的 RZ 一 ulp 容量缺陷。父任务失败版本已归档，用户授权持续推进至 I12 通过。

## 目标

只解决原修复“只能减小液量”的可达性限制；支持受限增加内部面的液量，并保持两侧守恒及原有 gross bounds。

## 允许修改

最多 1 个 production 文件：`include/Modules/Interface/Detail/HermiteGalsInventoryAssembly.h`。最多 1 个测试文件：`tests/interface_integrated_inventory_contract.cpp`。另含本 TODO、API/报告更新和 runs/baselines/I12a3r1_20261001 证据。

## 禁止修改

不改时间积分、Q4、boundary history、总通量、独立 M 的输入、CFL/散度/守恒门槛、bubble/reflux/driver 或冻结文件。不提高 64*eps 面修正/库存预算及 128 次修复上限。不修改边界面液量，不用 M clipping 或 global KKT，不将不收敛输入自动切到旧算子。

## 数学约束

每次进入修复前所有原始库存缺陷必须已在原预算内。仅在实际内部面增加液相正/负 gross 分量，其上界仍为该方向固定 total，且受邻居真实库存余量与相对原始面量的 64*eps*min(capacity) 预算共同限制。每个改变量被双方反号消费，重装配后再次检查。原有受限减量仍可用；未找到可行内部调整明确失败。

## 实现要求

在同一共享装配器中先检查邻接单元及方向余量，利用局部可达的增加量消除刚转纯相时的缺陷；随后保留原受限减量。用绝对变化同时审核增加/减少的面修正预算和诊断。保存原始量、最终量和修复次数；失败不提交任何 state。初次 I12a3 的 native total 与新积分 total 相同，本任务不通过更改度量计算顺序来绕过修复问题。

## 验收

- T0：Release、architecture、冻结与 1+1 文件预算。
- T1：I12a3 原失败矩阵完整通过（XY/RZ、uniform/2:1 AMR、正反速度、17 个相位及完整/组合区间）；原 q=1/4、3/4、dt=.02 的库存反例；独立 counterflow/减量/增量修复与超预算拒绝。
- T1：每个最终 face 保持 0<=liquid+/-<=total+/-，相对 raw 的每方向改变量均在原面预算；边界量与 total 不变；从最终共享面量重算 M 完全一致，严格容量和 ledger 通过。
- 回归：原 19 项检查、新 inventory contract、旧 crossing；记录任何原输出变化及其舍入修复原因，不降低标准。

## 失败规则

无相邻可行余量、不可表示的受限改变量、超出预算或达到原工作上限均明确拒绝。若仍缺一般可达性能力，记录新证据并拆独立子任务，不掩盖数值失败。

## 交付

通过或失败的完整矩阵、API/报告、源码/runtime 哈希和原始日志归档。只有 I12a3 的完整新矩阵通过后才进入 I09 区间消费。

## 执行结果

**PASS。** 1 production + 1 test 文件；148 组检查、20 类拒绝和 20 项 CTest 通过，19 项原完整输出不变。I12a3 随本修复通过，详见 [报告](I12a3r1_Bidirectional_Roundoff_Repair_Report.md)。I12 整体未通过，继续 I09 的区间积分消费。
