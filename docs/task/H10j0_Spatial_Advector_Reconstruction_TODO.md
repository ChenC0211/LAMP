# H10j0 / 可复用 native compatible advector

状态：**PASS / engineering PASS**，验收见[报告](H10j0_Spatial_Advector_Reconstruction_Report.md)与[实际API](../api/H10j0_Cch_Spatial_Advector.md)。本任务只关闭运动学 owner，不计作phase/material/PDE PASS；下一唯一任务H10j。下列为实施前固定要求。

目标：只实现与密度/phi无关的实际compatible advector：几何只准备一次，速度更新只生成不可变系数与fresh来源identity；供随后H10j与最终native/Interface readonly view消费。

允许修改：新 `include/Modules/Hydrodynamics/CchSpatialAdvector.h` 和 `src/Modules/Hydrodynamics/CchSpatialAdvector.cpp`，以及现 `CMakeLists.txt` 的CCH显式源owner登记（共3，<=5）；仅既有 `tests/contracts/hydrodynamics/cch_momentum_geometry_contract.cpp`（1 test）。既有cch_convection driver/README只能补充该实际owner观察、对应设计/API/报告。0新库/target/tools。开始前验证H10i2 flat manifest/canonical/archive并保存实施前设计、source/input目录、原数值入口。临时 runs/checks/H10j0，完成 runs/baselines。

禁止修改：Core geometry/root/trace、Interface primary/Q4/publication/transport/clock、CchPhaseGeometry/旧rate、H08/H09 F/R/uR、convection/压力/散度/黏性、App/字典/PDF与原gate。没有phi query、phase-rate、density、材料库存/bubble、质量修复、AMR特例、额外global DOF、解析输入修正或fallback。本轮不实施H10j。

数学约束：保持原native scalar总量 A_f*u_f、D_I、actual subface/shared trace；局部B的体积坐标div为原(D_I*u)/V。XY affine可识别时精确，RZ限定轴正规族及even-r² axial物理平均再现；axis不除零。它是运动学精度/来源要求，不给材料加回quarter/commuting/R=0。材料仅唯一F、内部反号/真实AMR子面和、来源一致、KEP/skew。

实现要求：按现API分 geometry preparation 与 actual evaluate；同一scaled三列physical-mean SVD最小范数，无rank fallback。原n4 selected2/3的56未识别tangent保留explicit deficit；full16及n>=8 selected不得套用此豁免。真实one-ring/面均值/axis、分片transfinite和零normal lift、actual cell/atom polynomial、native同点两分量都属于该owner。canonical overlay包括kinematic断点，按实际primalVelocity或ownerPressure去重；public immutable getters不暴露setter/裸cache。每次实际evaluate复制原输入并生成不可伪造advector identity，声明geometry/full clock原样保持；不冒称H14实际commit验证。几何副本共享storage，速度更新不重做SVD/atom图。

T0：标准local-benchmarks Release；architecture/self-test/layout与原phase/momentum/dual/KEP必要短契约、Core/FirstJet关键source hash；新cpp唯一编译owner且不进入legacy lib；原App/PDF/输入保持。无新调试程序target。

T1：16full布局（XY/RZ×uniform/真实2:1×n4/8/16/32）及原selected三类，独立face physical-mean、原D_I总量与cell divergence、normal反向/真实AMR区间加回<=1e−12 gross；所有native行/axis保留。可识别XY affine和RZ正规/even-r² axial族<=1e−12真实尺度；一般RZ point/mean差异单列O(h²)，不改原point参考或标成全affine exact。非仿射actual速度检查分片lift、atom polynomial与独立局部系数/积分参考，不以本owner自己的输出作参考。tiny interval使用endpoint-relative参数，不能靠两个root全局坐标相减。

T2：同一geometry先后evaluate两组actual值，确认geometry storage/SVD/atom库存不变而输入系数/identity改变；copied owner、源句柄释放后lifetime、wrong mesh/layout/regrid/clock及nonfinite/FP/axis拒绝、失败前后原owner数据不变。相同momentum instance来源显式核验。记录准备/更新工作量与局部storage随真实面/segment线性，不借本任务性能优化改变构造。

失败规则：来源/可识别再现/总量/divergence/复用或预算失败，冻结具体反例后修正或拆单一接口；无旧方法fallback/容限放宽。极粗selected的信息缺口只允许预列的精度限制，不能把未识别行从finite/shared/来源检查中删除。只有真正外部依赖才BLOCKED。

交付：实现/API、独立短契约和现研究actual owner证据、原protected/affine/KE比较、准备复用与来源库存、支持/拒绝范围、flat manifest。只关闭H10j0；随后H10j的causal phase factory，再分别迁移H10k/H10l与最终native/PDE。原内部affine先过，再原完整H10b全部原gate；H11/H08a仍等待。
