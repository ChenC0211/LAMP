# H10k / 唯一空间材料通量 owner 与来源迁移

**2026-10-04当前收口：** [H10l2](H10l2_Same_F_Dual_Material_Retry_Report.md)已PASS；旧quarter/legacy evaluate退役，H09直接共享H08原dual F，真实R/E_sp/G与uR/材料KE通过。现行接口见[H09 API](../api/H09b_Cch_Dual_Material.md)。下一唯一[H10m](H10m_Final_Native_Same_F_Convection_TODO.md)联合最终native空间矩与material-KE/skew，内部affine先于原完整H10b。Pure FirstJet独立；R/E_sp/G不要求零，H11/H08a继续等待。

状态：**PASS / engineering PASS**，2026-10-04，见[报告](H10k_Material_Spatial_Owner_Migration_Report.md)和[实际API](../api/H10k_Cch_Spatial_Material_Flux.md)。2 production/2既有tests、384实际准备及8/8必要回归通过；瞬时owner补充见[实施前addendum](H10k_Instantaneous_Material_Flux_Addendum.md)。下一唯一[H10l](H10l_Same_F_Dual_Material_Migration_TODO.md)；不授权提前执行H11/H08a。

目标：让现 H08 材料 owner 从同一 actual spatial phase atom 生成唯一 F_rho，并携带实际 advector identity。只解决材料来源和单份通量，不重设计 dual transfer 或 convection。

允许修改：现 `CchMaterialTransport.h/.cpp`（**2 production**）；现 `cch_phase_geometry_contract.cpp` 与 `cch_material_transport_contract.cpp`（**2 tests**）；既有 `cases/benchmarks/cch_convection/driver.cpp`/README、对应设计/API/报告与现行状态。原 CMake/target，0新库、target、tools 或变体文件。临时 `runs/checks/H10k/<timestamp>`；证据归 baselines。

禁止修改：Core/FirstJet/出版/Q4、advector 重构、CchPhaseGeometry 的数值算法、D_I/C、黏性/电场、密度定义、dual/KEP生产实现及其两份测试、App/字典/PDF、H10b gate。无库存、bubble、质量修复、AMR 特例或 fallback。

数学约束：只有唯一共享 F_rho、内部反号/真实 AMR 子面和、来源一致、KEP/skew 是材料硬约束。几何完全非守恒且独立；保留实际质量差和 R，不能强制 quarter/pressure→dual commuting/R=0。

实现要求：

- 实施前固定调用库存、actual source stamp 和一份原子材料率的所有视图/单位；报告 scope 与旧消费者的明确边界。
- `CchMaterialSources` 加 actual advector stamp；spatial 路线来自 midpointRate.advectorIdentity，不能由 caller 自造、贴给旧数据或用0充当有效实际来源。validateFor 比较三层 geometry、完整clock、原layout/D_I、same momentumGeometry instance及此 stamp；stamp 本身不是 H14 commit。
- H08 对每个 spatial phase atom 唯一计算 `rho_l*Q_l + rho_g*Q_g`，随后 primal/dual/material boundary 视图只从这份 canonical F 求和。phase 不能乘密度；H09 后续也不能再积分或另乘密度。两个相分量分别有限性/非零 representability 核验，不吸收细相；无可写缓存或备用 material owner。
- 原真实三层压力质量与 midpoint dt 积分、实际 continuity defect/balance/report 保留，报告其空间/时间误差范围，不能以 R=0 为接受条件。时间仍为原 midpoint 近似，不冒称时空高阶或接受事务。
- phase/material 两份原契约迁移到实际 spatial factory 和同一 geometry instance，用独立 physical polynomial/相交参考代替旧 face-constant 断言。研究 helper 使用其原 geometry，不另建实例；保留原 H10b full/RK4/donor/门槛与所有 native 行。
- 现 dual/KEP 两份测试仍调用旧 evaluate；本任务保留该显式历史入口及其旧 material 表达，只作为未迁移回归，不能标成新路线 PASS。新 spatial 工厂不回落 legacy。**H10l**迁移这最后两类 consumers 后一次删除 legacy evaluate 和 quarter 数值算法。禁止在本预算内暗改其 helper 或伪造新的 source stamp；移除旧入口的时序取代早期 H10j/k 草案。
- 不把旧 dual quarter 视图直接称为新同源 F。新的材料 dual rate port 只供下一 H10l 消费；当前生产 consumer 尚未全迁移。

验收 T0：标准 Release、architecture/self-test/layout 和原 phase/material/momentum/dual/KEP 必要回归；源码/输入/App/PDF库存保持预算。原 source/affine CSV变化逐项解释；不能把允许来源迁移后的数值变化当成无条件字节保持或 silently 改 gate。

验收 T1：actual spatial owner 的单份 F、primal/dual atom求和、内部反号、真实 AMR 子面和、原 geometry/full clock/actual advector stamp/source寿命及拒绝事务，≤1e−12 gross。XY/RZ、uniform/2:1、full/selected，零迹正负/零/变号、细相、轴和原未识别 tangent 的已声明范围保持。

验收 T2：正式三层移动/crossing，rho=1/1e3/1e6，独立真实几何质量、相率和材料率，actual R及全域 balance 身份核验。H10h pure gas dof103 的三层21段新 phase/F 与其 signed outward sum独立重查，不能仅用等密度或库存证明。copy/lifetime、wronggeometry/clock/layout/D_I/advector/nonfinite/regrid明确拒绝；没有外部 inflow port的zero拒绝保留。

失败规则：任何必要来源/单份F/真实相准确性/有限性/预算/回归失败即 FAIL，冻结反例并修正；不能用 legacy fallback、黏度、mass correction 或 relaxed R 接受。若发现仍需改 dual/KEP/Phase数值，拆下一任务，不突破2 production/2 tests。仅确实无法推进的外部依赖才 BLOCKED。

交付：实际 API/owner、执行与独立数值证据、全部仍未迁移 consumer、预算/compiled runtime/flat manifest。H10k后唯一下一任务 H10l：消费同一材料 F、同一 momentumGeometry，定义真实 E_sp/R 并删除 quarter；预核最多5 production/2 tests，不改最终 convection。再单独收口 native constant/affine 与 material KE，内部 affine先通过，再原H10b全部28窗口/84gate。H10/H11/H08a仍未关闭。
