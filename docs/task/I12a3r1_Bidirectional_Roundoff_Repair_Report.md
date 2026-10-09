# I12a3r1 / 内部共享通量双向舍入修复报告

**PASS；I12a3 库存装配验收随本修复通过。** I12 整体仍 FAIL，新区间库存尚未接 I09/driver。按用户授权继续下一独立任务。

## 修改与约束

仅修改 1 个 production 文件 `Detail/HermiteGalsInventoryAssembly.h` 和 1 个测试文件 `interface_integrated_inventory_contract.cpp`。详见 [TODO](I12a3r1_Bidirectional_Roundoff_Repair_TODO.md)。

同一装配器在每轮确认库存缺陷仍处于原 64*eps*capacity 预算后，检查邻接内部面的剩余 gross 液量空间。超容单元可增加内部出流；负库存单元可增加内部入流。改变量同时受方向 total、邻居库存余量、相对原始面量的原 64*eps*min(capacity) 预算约束；每次改动重新装配，双方反号使用同一个面量。之后仍可使用既有受限减量。

修正预算和最大修正诊断改为绝对差，统一覆盖增/减。总通量、边界液量、旧 M、CFL/散度/守恒容限、128 次工作上限均不变；没有 M clipping、第二套库存更新或全局 KKT。无法找到可行相邻调整仍明确拒绝，不声称已求解任意网络的全局可达性。

## 验证

- 主程序及 21 个相关测试目标 Release 构建通过，24.63 秒，无 warning。
- **20/20 CTest PASS**，111.02 秒。除主动扩展的新 inventory 测试外，19 项前序完整输出与 I12a3 失败归档完全一致，包含 architecture 和旧 I08 修复数据。
- 新库存合同 **148 组检查、20 类拒绝通过**：136 组 XY/RZ、uniform/2:1 AMR、正反速度、17 相位的真实 I12a2 积分装配；原代数/counterflow/减量检查；固定边界下两个方向的内部增量；8 组 dt=.02、q=1/4 与 3/4 的 XY/RZ 正反跨越。
- 十进制反例保留 h=.25 的跨越单元及原 XY 单面 0.00375/0.00125 参考；法向保留域取 [.25,.75]，以确保原生平面在整个旧支持上精确可表示，外边界历史显式给出。它验证库存与原单面积分，不冒称是完整单位域 driver。
- 全区间及组合区间核验每个最终 gross subset bound、严格 M 容量、从共享最终面量重算的逐单元 M 完全一致、global ledger、原面修正预算，以及 boundary/total 完全不变。用例累计记录 10 次正向、11 次反向 gross 分量增加（包含完整/组合两次消费）。
- I12a3 的 RZ q=12/16、5.55e-17 超容反例通过，不再把缺陷推向满液边界。超预算和不兼容输入继续明确拒绝。
- 69 冻结文件、12 集成参考、3 PDF 身份保持。I12a3 失败 archive 原样保留。
- 旧 `--crossing-probe` 仍退出 1，完整输出不变；未接新区间 driver，不能据此声称 I12 长程 benchmark 通过。

复现：`build/windows/tests/Release/quantum_interface_integrated_inventory_contract.exe`。原始日志、构建目标、source/runtime SHA256、20 项输出比较见 `runs/baselines/I12a3r1_20261001/`；归档哈希和逐条验证见其中 `archive_verification.json`。

## 下一项

让 I09 持有新区间库存结果并汇总其最终 actual subface integrated volumes，不重新积分、不再次修正 M。一般曲面/变化速度的时空源、完整步进及 I12 其它 gate 仍须独立实施验收。
