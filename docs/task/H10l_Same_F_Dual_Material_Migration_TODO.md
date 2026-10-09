# H10l / 同一 F 的 dual 材料迁移与 quarter 退役

状态：**FAIL / 恢复后的 engineering PASS**，见[报告](H10l_Same_F_Dual_Material_Migration_Report.md)。候选5 production/2 tests在原basis与closed必要回归失败，已保留全部证据并恢复H10k。下一唯一[H10l1](H10l1_Spatial_Advector_Closure_And_Boundary_TODO.md)，本轮未实施。下列为实施前固定要求，未放宽原gate。

目标：H09直接消费[唯一瞬时材料F](../api/H10k_Cch_Spatial_Material_Flux.md)的dual视图，保持真实三层dual质量、完整clock/actual advector/momentumGeometry来源；删除quarter数值延拓与最后legacy相率/材料入口。原production KEP算子仍消费H09，不在本项重设计最终spatial operator。

允许修改：`CchPhaseGeometry.h/.cpp`（仅删除legacy evaluate及其专属存储/代码）；`CchMaterialTransport.cpp`（删除legacy表达，不重算F）；`CchDualMaterial.h/.cpp`（same-F消费、真实R/E_sp诊断）。固定**5 production/2 tests**：仅既有`cch_dual_material_contract.cpp`和`cch_kep_convection_contract.cpp`迁移。现research driver/README、实施前设计/API/report与当前状态可更新，原CMake/target，无新源码/tools/variant。

禁止：Core/FirstJet/Q4/phase spatial数值、advector重构、D_I/C、密度定义、原dual几何/质量定义、黏性/电场、`CchKepConvection.h/.cpp`、App/字典/PDF或原H10b gates。无库存、bubble、mass repair、AMR特例、失败fallback；不得改回quarter/精确material commuting/R=0硬条件。

## 实施前数学与来源固定

- 对原真实momentumGeometry的每条segment直接消费`material.dualFaceRates()`，不重新积分/乘密度、按面积分配primal F或求quarter。核验material原geometry instance、actual advector stamp、原layout/D_I和三层mass来源；copies保留原prepared owner。
- 以真实dual incidence计算`q_dual=sum s F`、中点积分及`r_dual=delta_m_dual+dt*q_dual`；R/uR和材料KE继续使用此同一R。mapped pressure quantities只为诊断。
- 定义`E_sp=q_dual-T*q_primal`（率单位）及`G=delta_m_dual-T*delta_m_primal`（质量单位），其中T仅为原真实overlapFraction。固定并独立核验`r_dual=T*r_primal+G+dt*E_sp`，不要求E_sp/G为零。quarter/commuting旧report字段可具名保留历史/诊断含义，禁止赋伪零冒称通过。
- 最后两份契约实际spatial factory复用原momentumGeometry；迁移其旧A_phase*u/quarter专属断言为独立physical integral、actual incidence/source和material-KE/skew断言，保留原其余测试与PDE精度门槛。无调用者后一次删除`CchPhaseRate::evaluate`和material legacy分支。现Transport公开spatial接口保持，不能隐藏增加第6生产文件。
- 研究helper已迁移原B；只允许观察新的actualF/R/E_sp。原full/RK4/donor/原gate不变。临时检查`runs/checks/H10l/<timestamp>`，冻结baselines；不把observer exit0当PDE通过。

## 验收 T0/T1/T2

**T0：** 标准Release与architecture/self-test/layout；受影响phase/material/momentum/dual/KEP必要回归。源成员/inputs/App/PDF/Cmake保持预算；最后旧evaluate及quarter算法调用库存为0。旧历史档案不改，旧失败不改判。独立绑定compiled源码、实际exe/library/DLL、退出记录和flat SHA256 manifest。

**T1：** XY/RZ、uniform/2:1、full/selected、rho1/1e3/1e6、zero正负/零/局部变号、小相、轴、真实跨面三层。唯一F/原atoms→dual共享、真实内部反号及AMR subface sums、来源/完整clock/geometryinstance/copy/lifetime/regrid/nonfinite明确拒绝，所有身份≤1e−12 gross。缺失selected/exterior流入相仍拒绝。原56粗selected缺口如实保留，不新增豁免。

**T2：** 独立三层dual质量、真实q/R、mapped量/G/E_sp及上述分解、consistencyImpulse/materialKineticChange与KEP/skew身份共同核验；新gas dof103三层21段×三密度q为正确closed sum，实际R保留不能填零。几何和advector运动一致的coupled制造族才可讨论R收敛；强制不同运动只证明诊断身份。内部原affine counterexample先重跑；若最终native受新F影响而失败，保留具体witness和最终operator后续任务，不靠trace补丁/放宽门槛修复。

失败规则：有限性、同源/单份F/actual sharing、独立真实phase/F/质量或预算/必要回归失败即FAIL并冻结反例；只外部依赖确实阻碍才BLOCKED。需改KEP/其他production时拆下一任务，不超预算，不用R=0或旧quarterfallback遮掩。

交付：真实API/report、删除库存、独立actualCSV、scope/runtime/frozen manifest及明确下一单任务。只有H10l关闭；下一任务联合设计最终native constant/affine空间矩与material-KE/skew，内部counterexample精确/二阶通过后再原H10b全部28窗口/84gate。H10/H11/H08a未自动通过。
