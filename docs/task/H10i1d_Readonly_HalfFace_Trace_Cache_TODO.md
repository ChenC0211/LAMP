# H10i1d / 同源只读半面 trace cache

状态：**PASS / engineering PASS**；本轮唯一任务已完成，见[报告](H10i1d_Readonly_HalfFace_Trace_Cache_Report.md)/[API](../api/H10i1d_Readonly_HalfFace_Trace_Cache.md)。已核验H10i1c的175 archived payload/61当前交付，实施前[设计](H10i1d_Readonly_HalfFace_Trace_Cache_Design.md)及目录/gate快照见`runs/checks/H10i1d/20261004_111021_237`。16实际来源/1728半面、short16布局/160零半面/64hanging及同源细化/失败事务、OMP2/496半面通过；17/17回归、原来源29CSV和内部affine保持。下一唯一任务[H10i1e](H10i1e_FirstJet_Trace_Publication_Measure_TODO.md)，H10仍未关闭。[总契约](H10i1_Coincident_Phase_Trace_Clip_Contract_Design.md)及 [Core trace API](../api/H10i1a_Native_Q4_Trace.md)。

目标：为真实 native half-face 准备唯一 owning immutable typed trace cache，复用 Core `NativeQuarticTrace`；普通 legacy phase view 保持，零迹的两侧相和 normal Taylor events 不丢失。

允许修改：仅 `include/Modules/Interface/InterfacePhaseGeometry.h`、`src/Modules/Interface/InterfacePhaseGeometry.cpp`（2 production），既有 `tests/contracts/interface/interface_hermite_gals_flux_contract.cpp`（1 test）；既有 research/docs 可更新。0新工具/target；超过预算先拆任务。开始前保存设计、来源/输入目录和原 gate，核验 H10i1c 平面 manifest。

禁止修改：Core Q4/root/regularity，Interface primary/transport/clock/volume/publication，Hydro 材料/对流/散度/压力，App/字典/PDF和验收容限。无库存/bubble/几何体积拟合/质量修复、scalar 微扰、epsilon 零判据、AMR 特例或临时调试开关。

数学约束：每条实际 half-face 的 extent、native row/half identity、同一 owning Snapshot、typed normal axis 必须一致；反向输入和实际 AMR 子面 partition 可加。事件细化继续原 immutable same-root refinement，失败保留原 cache；普通 consumer 的有限根不能混入 normal Taylor event inventory。inactive/外部侧 availability 如实保留，不外推。

实现要求：显式 `MappedFacePartitions::prepareNativeTrace` / `traceHalfFace`（最终签名先在设计核对现有 API）；cache 只保留原 owning typed 查询。`halfFace` 对 ordinary 给原 legacy view，对 positive-length zero 明确拒绝；不能用对称 half 相替代 causal 材料率。细化后重新生成依赖测度，不把原始 roots/ranges 复制成另一数值来源。

T0：标准 Release、architecture/self-test/layout及受影响旧 face/flux/FirstJet/publication/region/phase/material/dual/KEP 回归；H10i1c 普通来源、闭合 clip、实际 Core native-side 与内部 affine 保持。

T1：XY/RZ、uniform/真实 2:1/hanging、正负相、normal order/event、ordinary→zero 混合半面；同源/错源/实际 extent、selected availability、copy/move/lifetime、细化和失败事务。独立闭式物理面积/两侧液气测度、反向和 actual 子面加回 ≤1e−12 实际尺度；原严格 regularity/FP/工作预算保持。旧 publication native 八组仍拒绝，Core PASS 不能冒称完整材料链通过。

失败规则：任何来源、identity、事件、measure、反向/partition 或必要回归失败，保存具体反例后修正或拆任务；无 fallback、容限扩大或改判旧失败。

交付：API/设计/报告、既有短契约结果、实际编译源码/二进制/退出和 CSV、完整目录范围审计及平面 SHA256 manifest。随后单独登记 H10i1e FirstJet 派生出版测量；H10i1f/H10i复核后才 H10j。最终 native 内部 affine 先通过，再原 H10b；H10未关闭，H11/H08a等待。

材料硬约束仍只有唯一共享 F、内部反号/实际 AMR 子面和、来源一致、KEP/skew。quarter/pressure→dual 精确交换与 R=0 仅诊断；继续真实 R/uR，不约束或修正 Pure FirstJet 几何。
