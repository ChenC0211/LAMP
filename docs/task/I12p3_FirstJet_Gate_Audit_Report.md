# I12p3 / 当前纯路线 Gate I 审计

**Gate I：FAIL，内部交付缺口；无外部阻塞。** p1/p2a/p2b 的数值结果保持 PASS。审计按[总纲的完成条件](../00_LAMP_Development_Charter.md)检查，当前还不能将模块测试驱动当作 production Interface。

| 要求 | 当前证据 | 裁决 |
| --- | --- | --- |
| 唯一 scalar 服务位置/法向/曲率/cut/体积 | [p1 API](../api/I12p1_Pure_Hermite_Gals_API.md)、p1/p2b 契约中的 owner 检查 | PASS |
| 直接 first jet、同一实际 B_h/J_Bh、I06 mixed closure | [p1 报告](I12p1_Pure_Hermite_Gals_Publication_Report.md)，p2b 回归保留 trace/vertex/edge 完整输出 | PASS（已测范围） |
| 无 M、无体积拟合、NA 不冒充零守恒残差 | FirstJet 状态拒绝 M；12 原完整算例记录几何体积与漂移 | PASS |
| 静止/平移/旋转/vortex 及收敛/crossing | [p2a 报告](I12p2a_FirstJet_Convergence_Report.md)：三档空间/时间、48 crossing | PASS（原规定范围） |
| 动态 AMR 后继续实际输运 | [p2b 报告](I12p2b_FirstJet_AMR_Transfer_Report.md)：exact/projection、粗化三档、1,036 步/36 regrid | PASS |
| 正式模块时间步组合与 App 接线 | 完整纯时间步组合仍在测试文件；正式 App 无 FirstJet 入口 | FAIL |
| 字典/实际配置与 production 输出诊断 | 现有 reader 不接受纯路线；新模式未接 accepted VTK/series/meta | FAIL |
| 效率与规模增长证据 | 有完整运行 wall time、无新增全局高阶物理状态；缺当前纯步按阶段/网格规模统计 | FAIL（证据不足） |

证据中的速度是独立模块基准的规定速度，仍未接 Hydro。取消相体积守恒没有自动证明密度/动量一致、低黏稳定或电喷雾破碎准确；这些仍属于后续 H/IH/E/C 顺序。

真实缺口定位：

- [benchmark 的 `predict` 与 `finishFirstJet`](../../tests/interface_hermite_gals_benchmark.cpp) 已调用 production tracer/vertex/edge/I06/publication，但组合、中心回溯和真实外边界数据接线位于测试匿名域。App 不能复用这个时间步入口；应提取模块 producer，保持原输出和原拒绝语义，不复制一份公式。
- [正式字典 reader](../../src/Apps/Interface/IO/InterfaceBenchmarkInput.cpp) 的方法白名单仅有 `levelSet/conservative/conservativeGals/parametricGals`；[advection dispatcher](../../src/Apps/Interface/Verification/InterfaceAdvectionBenchmark.cpp) 也没有当前纯路线分支。旧名 `conservativeGals` 不等于此次 FirstJet，不能据名称宣称可用。
- p1/p2 驱动输出到测试 stdout，档案里的日志/JSON 是验证证据；它们不是规范要求的 production `output/fields/solution.pvd`、accepted-frame、series/meta 和有效字典追踪。
- p2b 动态 full 238.31 s、p1 固定网格 full 226.65 s 是原规定场景的一次 wall time，不是可比网格规模研究，也不能证明每阶段结构复杂度。暂不做数值优化或改变容限。

下一步按单任务顺序补齐：

1. **I12p3a**：提取可复用纯 first-jet 时间步 producer。只改变时间步组合归属，保持旧源/拒绝事务、真实 trace 与外边界含义，用原纯基准比较完整数值输出。
2. **I12p3b**：正式独立 App/字典入口与 typed 配置，调用新 producer；严格拒绝不支持的选项，验证接受时钟及动态网格接线。
3. **I12p3c**：接通规范的只读 VTK XML、accepted/rejected、volume drift、series/meta/effective config，验证输出频率不改变解。
4. **I12p3d**：在同一实际 App 路径上记录规模、阶段耗时、真实状态/几何数量与内存，保持误差门槛；若出现结构问题另立任务，不混做优化。
5. 回到 **I12p3** 汇总复核，全部必要项有证据后才可声明 Gate I PASS。

本审计不改 production/test/case，不重复已归档确定性计算；校验 p1/p2a/p2b archive、81 freeze、3 PDF 和当前源码范围。机器证据：[矩阵与缺口](../../runs/baselines/I12p3_20261002/gate_matrix.json)、[审计汇总](../../runs/baselines/I12p3_20261002/verification.json)。
