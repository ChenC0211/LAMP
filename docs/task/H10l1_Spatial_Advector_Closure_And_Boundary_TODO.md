# H10l1 / 实际 advector 的零净流量精度与封闭边界来源

状态：**PASS**，见[报告](H10l1_Spatial_Advector_Closure_And_Boundary_Report.md)。仅完成2 production/2 tests范围；下一唯一[H10l2](H10l2_Same_F_Dual_Material_Retry_TODO.md)，尚未实施。下列为实施前固定要求，原gate不变；H10l FAIL保持历史结论。

目标：同一actual compatible B既保留原D_I/物理面均值，又能消费明确的不可穿透边界条件；定位零净native流量的closure反例，并从实际trace work推导和验证误差预算。不能预先认定只是1e−300标度问题。无第二速度场、F、几何或材料修复。

允许修改：现`CchSpatialAdvector.h/.cpp`，**2 production**；现`cch_phase_geometry_contract.cpp`及`cch_material_transport_contract.cpp`，**2 tests**；既有研究driver/README、实施前设计/API/report与current状态。原CMake/targets，无新源码/tools。临时`runs/checks/H10l1/<timestamp>`，关闭时冻baselines。

禁止：Core/FirstJet、phase spatial积分、材料Flux/Transport、dual/KEP生产及两tests、D_I/C、原密度/三层几何质量、App/字典/PDF及H10b gates。当前生产保持H10k，不恢复失败quarter之外的新数值分支或通过忽略failure让H10l自动PASS。

数学硬约束：actual B/trace/normal mean和原D_I来源一致，唯一atoms/真实AMR subfaces、完整clock/layout/geometry instance。材料四硬约束保持，Pure FirstJet独立，无库存/bubble。`u_f=0`只代表物理面均值零，**不能自动推断pointwise impermeability**；内部零净counterflow必须保留。

## 实施前推导与最小实现

- 从冻结的n4 XY AMR basis u3、RZ uniform/AMR basis u2出发，定位实际cell/piece、原DI、normal trace系数/物理均值及signed/gross积分。复现原reject，不先改容差。区分真实几何/系数错误与零净流量下可证的浮点work误差。
- 保持原DI为独立authority。优先稳定装配physical-mean约束及实际积分；误差标度如确需含zero-mean trace variation，必须从实际多项式L1/work与浮点舍入推导并独立核验，不能加经验floor、放宽1e−12或把field.div直接贴成expected遮蔽真实积分。全零数据仍必须精确零。
- 原streamfunction在native边界均值为0而子段F非零（XY约1.367e−3，RZ约7.230e−3），没有显式pointwise边界law。新增通用typed只读homogeneous impermeable normal-trace来源，在同一prepared B施加；明确选择真实native boundary/selected边界、核验layout/来源/有限性和实际zero native均值。不得按case名称、AMR位置或mean==0自动压斜率。
- 不可穿透指定面必须整个normal polynomial严格0；一般边界/内部mean0仍可有正负counterflow。保留普通MMS/selected56未识别tangent库存、RZ axis真实测度与所有原causal phase拒绝。
- 同一B的原速度重构/lift和atoms/primal/dual views共同消费此边界law，不能在material/H09另压F。无主clock/接受事务；H14/H15一般边界仍单独验收。

## 验收 T0/T1/T2

**T0：** 标准Release、architecture/self-test/layout、原phase/material/momentum/dual/KEP必要回归；原其余源码/inputs/App/PDF保持，明确编译/实际二进制与flat manifest。当前旧dual/KEP回归只证明H10k原范围，不称新迁移通过。

**T1：** n4 XY/RZ uniform/2:1全部非轴单面basis及线性组合，不再错误拒绝；独立physical trace integral/D_I/actual B volume divergence满足推导误差预算。零净及近相消流量、tiny/大动态范围、全零、same-source/clock/instance、wronglayout/regrid/nonfinite拒绝保留，无silent fallback。

**T2：** 原四封闭streamfunction在明确law下通过所有真实boundary atoms/subsegments零normal；原DI与实际B散度、phase/material通量及内部sharing独立核验。未绑定law的原相同数据仍显式展示zero-mean变号trace，不能误归封闭。内部/一般边界counterflow、原tiny phase own-scale与64布局/56缺口无新豁免。wrong binding/lifetime/copy、selected外部phase流入限制不变。

失败规则：任何必要身份/支持范围/精度或预算失败为FAIL并保留counterexample；需要改变phase/material/dual/KEP时先拆任务，不暗改原gate。只有真实外部依赖反复阻断才BLOCKED。

交付：实施前设计、实际API/来源、独立witness与运行报告、scope/compiled runtime/frozen manifest。只关闭H10l1；随后另起任务按原5 production/2 tests重试H10l same-F dual迁移，在closed fixture显式绑定已验收law。再联合final native spatial moments/material KE，先内部affine再原H10b28窗口/84gate；H11/H08a仍等待。
