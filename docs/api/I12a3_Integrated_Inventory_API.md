# I12a3 / 区间积分库存装配 API

`assembleHermiteGalsIntegratedInventory(mesh, interval, faces, config)` 直接消费按 native row 排列的完整 `HermiteGalsIntegratedFaceTransfer` 数组。所有面必须与同一个 interval 同源；begin 必须为旧 accepted 时刻，不能从较晚子区间直接更新旧库存。晚起始子区间可先通过 I12a1 的面量组合组成完整区间。

结果 `HermiteGalsIntegratedInventoryResult` 保留 interval、原始 `rawFaces()`、共享最终 `faces()`、新库存和 `report()`。没有 fabricated midpoint geometry identity，也不发布/提交新几何。原始输入与 accepted state 不变。邻居使用同一个 signed liquid volume 反号交换；四个 gross volumes 保持非负和液量不超过总量。

装配与旧 I08 共享私有 `Detail/HermiteGalsInventoryAssembly.h`。旧路径以原 dt*rate 单位进入，新路径以 integrated volume 单位进入；CFL、总散度、容量、整体 ledger、64*eps 内部面舍入修复预算与 128 次上限保持原状。结果检查源/时间/mesh/完整面，并不替外部提供者认证任意输入的求积或边界历史。

I12a3 初次交付为 FAIL：真实 RZ 跨面矩阵暴露旧单向舍入修复的可达性缺陷。**经 [I12a3r1](../task/I12a3r1_Bidirectional_Roundoff_Repair_Report.md) 修复，当前 I12a3 验收 PASS。** 在原预算内可增加具有余量的内部液相 gross 分量，改变量同时受邻居库存空间与该方向 total 约束；每次改动重装配。边界与 total 不变，绝对修正诊断覆盖增/减；仍保留原减量和工作上限。148 组库存检查、20 类拒绝及 20 项回归通过。此能力不代表 I09/driver 已接线，也不保证任意不可达网络都能修复。
