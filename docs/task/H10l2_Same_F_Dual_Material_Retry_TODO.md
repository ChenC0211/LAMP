# H10l2 / 同一 F 的 dual 材料迁移重试

状态：**PASS**，见[报告](H10l2_Same_F_Dual_Material_Retry_Report.md)。本轮完成5 production/2 tests；下一唯一[H10m](H10m_Final_Native_Same_F_Convection_TODO.md)，尚未实施。以下为实施前固定要求，原H10l FAIL及H10b gates保持。

目标：H09 直接消费原 `CchMaterialFlux::dualFaceRates()` 的唯一共享存储，保留真实三层 dual 质量与同源 clock/actual B/geometry；删除 quarter 延拓及最后 legacy phase/material 分支。现 KEP/uR 继续消费真实 R，不将严格材料 commuting 或 R=0 加回硬条件。

## 允许范围与禁止项

固定 **5 production / 2 tests**：`CchPhaseGeometry.h/.cpp` 仅退役 legacy evaluate，`CchMaterialTransport.cpp` 仅退役 legacy 表达，`CchDualMaterial.h/.cpp` 完成 same-F 消费及真实诊断；仅迁移现 dual/KEP 两个契约。既有 research driver/README、实施前设计/API/report/current 状态可更新；无新源码、target、tools 或变体。临时文件归 `runs/checks/H10l2/<timestamp>`，关闭时冻结 baselines。

禁止修改 Core/FirstJet、phase spatial 数值、已验收 actual advector/boundary、原 D_I/C、真实 dual 几何与密度/质量定义、KEP 生产算子、App/字典/PDF 或原 H10b gates。无库存、bubble、mass repair、AMR 特例或失败 fallback。H10l 归档候选只作参考，不能直接假定与当前来源兼容。

## 实施前数学与来源

- 对每条真实 segment 直接引用 material 原 dual view；不重新积分、乘密度、按面积分摊 primal F 或求 quarter。原 layout/D_I、momentumGeometry instance、actual advector stamp、三层 geometry/full clock 均严格核验。
- 真实 incidence 给出 `q_d=sum(s*F)`、`r_d=delta_m_d+dt*q_d`。uR/material KE 使用这同一 R。原真实 overlapFraction 的 T 仅用于诊断：`E_sp=q_d-T*q_p` 为率，`G=delta_m_d-T*delta_m_p` 为质量；独立核验 `r_d=T*r_p+G+dt*E_sp`，不要求 E_sp/G/R 为零，不伪填零。
- closed streamfunction fixture 显式绑定 H10l1 typed law，验证每条指定边界 atom/subsegment 的 normal/F 为零。不能凭零 native mean 推断不可穿透；未指定面的相反向流必须保留。全部原 basis 使用现实际 B，无新容差/豁免。
- 迁移最后两份测试的 legacy A_phase*u/quarter 专属断言为独立物理积分、真实 incidence/source/KEP 身份；原其余测试与 PDE 门槛保留。确认无调用者后删除 legacy API。原研究 full/RK4/donor/gates 不变，仅增加被动实际诊断。

## 验收 T0/T1/T2

**T0：** 标准 Release、architecture/self-test/layout 及 phase/material/momentum/dual/KEP 必要回归；1044 起始 source 成员按实际启动库存核对，5/2预算，CMake/inputs/App/PDF 不变。legacy evaluate 调用与 quarter 算法库存归零。绑定实际 compiled source、exe/lib/DLL、运行退出与 flat SHA256 manifest，原 H10l 失败证据不改。

**T1：** XY/RZ、uniform/2:1、full/selected、rho=1/1e3/1e6；zero、正负、局部变号、小相自身尺度、axis、三层跨面。独立核验原 atoms→primal/dual 同一 F、内部反号及真实 AMR 子面和；完整来源/clock/instance/copy/lifetime/regrid/nonfinite 拒绝。保留原56粗 selected tangent 缺口，不新增豁免。原208 basis 与四 closed 场景的迁移后必要回归通过。

**T2：** 独立三层 dual 质量、q/R、T/G/E_sp、上述分解及 consistencyImpulse/materialKineticChange/KEP-skew 身份，均按实际 gross/work ≤1e−12核验。gas dof103 的三层21段×三密度及 coupled 制造族必须记录真实 q/R。强制几何与速度不同的族只证明诊断身份，不宣称 R 小余项或 PDE 收敛。重跑原内部 n8/16 affine；若新 F 使最终 operator 不通过，保留具体 witness，另拆最终算子任务，不按 segment trace 补丁或放宽 gate 收口。

失败规则：必要回归、同源/唯一F/共享、独立物理积分/质量或预算失败即 FAIL 并冻结反例。需要改 KEP 或额外 production 时拆下一项；只有真实外部依赖反复阻断才 BLOCKED。

交付：实施前设计、真实 API/report、退役库存、独立 CSV、scope/runtime/frozen manifest与明确下一单任务。通过后才联合最终 native constant/affine spatial moments 与 material-KE/skew，内部反例精确/二阶通过后重跑原 H10b **28窗口/84 gates**；H11 继续等待，H08a 同源时空材料通量随后单独升级。
