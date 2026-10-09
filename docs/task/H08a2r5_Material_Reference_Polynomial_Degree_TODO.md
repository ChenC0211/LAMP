# H08a2r5 / 材料与 dual 参考积分解除固定次数

状态：**PASS（参考消费者）**，见[报告](H08a2r5_Material_Reference_Polynomial_Degree_Report.md)及[实施前设计](H08a2r5_Material_Reference_Polynomial_Degree_Design.md)。720项新检查、672条独立高精度复算及5/5回归通过；原1154行保持，production 0/test 2。下一[H08a2r6](H08a2r6_Benchmark_Polynomial_Schema_TODO.md)处理benchmark系数列。输入为 [H08a2r4](H08a2r4_Regular_Advector_Preparation_Report.md) 的生产准备层；实际 B/atom/J 仍为原来源，H08a/H11 未通过或启动。

目标：在升级 production 法向多项式前，先关闭两个编译期固定三系数的材料参考积分消费者，建立独立的 quartic normal、RZ 加权 quintic 积分检查。只解决参考消费者可接受真实次数这一事实，不切换生产来源。

允许修改：production 0；test 2：`tests/contracts/hydrodynamics/cch_material_transport_contract.cpp`、`tests/contracts/hydrodynamics/cch_dual_material_contract.cpp`。审计和中间结果归本任务 `runs/checks`。若必须增加其他消费者，先拆任务，不扩大本项。

禁止修改：Interface/FirstJet、接受状态、材料 F/R/uR、KEP、生产 B/J/atom、压力/黏性/EHD、App/字典。禁止通过丢弃高次系数、隐式转回三系数或提高容差实现兼容。

数学约束：四材料硬约束保持。参考积分别于生产算法，有解析/充分阶数核验；真实 XY/RZ 测度、正负 normal 号、shifted 两端坐标和极短正区间不能退化。旧二次来源的确定性材料/dual 输出保持。

实现要求：先逐项核对 `referenceAtom` 固定数组和 work-scale 假设；按输入次数处理全部系数，用有独立来源的 quartic 值验证 RZ 测度后最高五次。保留原 short contract 入口，不增加问题枚举或生产开关。benchmark 导出列、phase/momentum 参考及真实 atom/J 迁移留在各自预算任务；它们在部署之前必须全部关闭。

验收：T0 标准 Release、架构/目录、原 material/dual 契约；两份旧确定性输出逐行保持。T1 新 quartic 的完整/部分/双端 shifted 积分、正负方向及 XY/RZ 与独立解析结果一致，极短区间非零贡献不丢失；附原二次回归。T2 新检查仍为秒级短契约，不在 tests 中引入完整 benchmark。未改生产来源，本项无需重跑完整 H10b，但不得沿用本项结果替代正式部署后的 H10o/原 affine/H10b-v2 28 窗口 84 项。

失败规则：任何旧输出改变、高次项丢失、测度/方向或短区间失败均 FAIL；失败候选退役并保留证据。不因内部实现困难使用 BLOCKED；仅真实外部依赖阻塞才可报告。

交付：两项次数独立的参考消费者、独立积分证据、回归及 source/runtime 冻结、剩余真实消费者清单。完成后再按预算升级完整多项式表示及统一 B/atom/J，然后返回 H08a3 与同源时空材料积分。一次仍只推进一项。
