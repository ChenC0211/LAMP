# H10i1a / 同源native Q4零迹与两侧相port

状态：**PASS（Core 只读来源接口）/ engineering PASS**，2026-10-04；[验收](H10i1a_Native_Quartic_Trace_Port_Report.md)、[实施前 API 设计](H10i1a_Native_Quartic_Trace_Port_Design.md)。H10i1 的 108 payload 已独立核验，上一轮为 progress。本轮4 production/1既有测试；实际1,408段/40零迹/80侧相区间通过，未迁移clip消费者、publication、H08/H09材料或最终convection。真实质量/裁剪/native出版仍拒绝；下一唯一任务[H10i1b](H10i1b_Closed_Clip_Moments_TODO.md)，完整来源链仍按H10i1b–f推进。

目标：保留正长度恒零restriction，认证实际两侧normal Taylor的相区间，给clip内部迹和后续共享材料因果迹一个统一来源。Legacy bool phase Range/根/敏感度仍拒绝零迹，不以half分配、采样epsilon或rounded Q4替代。

允许production四文件：`include/Core/Geometry/LevelSetMappedSegment.h`、`include/Core/Geometry/Detail/LevelSetMappedSegment.h`、`include/Core/Geometry/Detail/LevelSetMappedCellLayout.h`、`src/Core/Geometry/LevelSetHermiteGeometry.cpp`。测试仅既有`tests/contracts/interface/interface_root_enclosure_contract.cpp`（归quantum_core）。可更新既有研究driver/API/任务文档。0新工具/0新CMake target；临时结果`runs/checks/H10i1a/<时间戳>`，冻结证据baselines。超过4 production/1 test须先明确拆任务，不静默扩张。

禁止改变：Interface primary/Q4/transport/publication、库存/bubble/体积拟合或质量修复、D_I/C、H07–H10生产消费者、App/字典/PDF/原门槛、全局DOF。不得穿透Core Detail给Hydro或复制第二套精确多项式。新port不是接受状态、clock或材料通量owner。

实施前写API设计：在既有owner中定义独立typed immutable NativeQuarticTrace、显式factory及normalAxis；source owning/usesSource、range/root/chart/native side来源、availability、normal order/sign和同阶refinement。普通与零range严格分型；regular legacy view只在无零range时有效。source源active支持不等于Hydro selected支持；后者由后续consumer显式限制。

数学规则采用H10i1：value限制所有exact coefficient为0才是恒零；每侧按真实native incidence取得0–4阶同源restrictLocal，首条非恒零normal系数隔离所有实际根并在正长度开区间认证sign*side^order。derivative dyadic/refine也必须限制同一阶偏导。事件overlay使用retained rational charts；不从rounded root/全局端点差造小相。全零normal Taylor、uncovered支持、非轴对齐、FP/环境/moved-from、预算耗尽明确拒绝；缺外侧是显式unavailable，不外推。

T0：标准Release、architecture/self-test/repository_layout及受影响原root/shiftedGreen/FirstJet/phase/material/dual/KEP回归通过；App/输入/PDF保持。原普通constructor及其拒绝语义、原volume/geometry输出无数值改动。

T1：新typed实际输出核验（不是旧异常计数）：axis0/1、XY/RZ native uniform/真实2:1、内部段/原生面/hanging subface、反向与分片、外界/selected support声明；positive/negative Q4；`(x-a)*(1+y+y²)`、偶阶/三阶与normal系数过零事件、原普通横穿/相切、多重根、短相和近零非零trace。native阶数/order/sign/全零判据须exact认证，测度/多项式reference≤1e−12实际尺度；细相relative gates沿用原契约，不统一改成绝对1e−12容差。

实际覆盖H10i1已发布dual内部8组中的1,408段/40零段，ordinary ranges保持，zero ranges完整且两侧相正确；原生面不能发布时用同一正式native scalar fixture直接验新Coreport，同时保留publication拒绝记录，不能标为已修复。source寿命、wrong source/axis、moved-from/环境与失败事务检查保持；普通旧root/variation调用不能读到缺失/假根数据。

FAIL规则：来源、exact零判据、侧相/事件、反向/子面或必要旧回归任一失败即FAIL，冻结原因后修正或拆明确接口工作，不fallback、不改几何或加容差。非外部困难不BLOCKED。

交付：正式Core API/实现、同一既有短契约、实际新port输出及独立reference、所有未支持项/consumer尚未移植项、平面manifest/protected核验。下一H10i1b单独实现真实clip测量。H10i历史FAIL保持，H11/H08a等待。
