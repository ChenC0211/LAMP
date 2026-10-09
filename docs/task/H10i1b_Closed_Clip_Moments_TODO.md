# H10i1b / 同源零迹的闭合裁剪矩

状态：**PASS / engineering PASS**；仅关闭本任务，[验收](H10i1b_Closed_Clip_Moments_Report.md)。前置冻结核验 PASS，实施前设计见[设计](H10i1b_Closed_Clip_Moments_Design.md)。[总契约](H10i1_Coincident_Phase_Trace_Clip_Contract_Design.md)、[Core port](../api/H10i1a_Native_Q4_Trace.md)。只实现真实 clipped region 的 geometry-only Green moments 与一个正式消费入口；不改 native-edge publication、材料 F、最终 convection 或时间推进。

目标：同一 native Q4 与 NativeQuarticTrace 的 region-inside 相迹共同给出闭合液/气矩，free interface 与 fixed clip 重合只计一次。首先让已成功出版的 dual 内部位置八组真实 region/momentum mass 成功，不能以单独面迹 PASS 替代闭合体积。

允许 production 五文件：`include/Core/Geometry/LevelSetMappedInterface.h`、`include/Core/Geometry/Detail/LevelSetMappedInterface.h`、`src/Core/Geometry/LevelSetHermiteGeometry.cpp`、`include/Modules/Interface/Detail/HermiteGalsRegionMeasure.h`、`src/Modules/Interface/HermiteGalsPhaseRegion.cpp`。测试仅既有 `interface_phase_region_contract.cpp` 与 `cch_momentum_geometry_contract.cpp`（2）。可更新既有 research driver/API/任务文档；0新工具/target。超过5/2先拆任务。临时证据 `runs/checks/H10i1b/<时间戳>`，冻结基线；先验 H10i1a 的 manifest，再写实施前设计及固定 gate。

禁止改变：Interface primary/Q4/transport/clock、库存/bubble/体积拟合/质量修复、legacy 敏感度、H08/H09/F/R/uR、native convection、D_I/C、App/字典/PDF、原容限。NativeQuarticTrace 只有几何来源，不得从 terminal rounded 区间自建 Q4/root/overlay；不得用采样 owner 投票或坐标 epsilon 判定 free/fixed 重合。不要求逐 cell Hermite 几何匹配任何守恒库存。

API 设计：在既有 Core owner 增加显式 `integrateNativeQuarticRegionMoments`，输入同一 snapshot 的实际 graphs、typed cuts、真实 regions、XY/RZ 与共同 Green origin；输出两相闭合 axial/radial moments、真实 arc/root refinement 计数与 error estimate，无修正 Jacobian。复用原 regular graph 隔离/controlled GL engine。zero clip 取盒内侧；原 point-only endpoint incidence 单独保留，不冒称一条零曲线；同源支撑/正规 branch 认证重合后只计 fixed 项，内部 free arc 保留。缺侧、inactive 正面积交集、未认证退化与预算不足明确拒绝，失败不污染 prepared source。

T0：标准 Release、architecture/self-test/repository layout、新 Core trace 与原 root/shifted Green/FirstJet/phase/material/dual/KEP 受影响回归保持；生产预算、原 source/input/PDF 和候选数据实际来源逐项核验。

T1：独立平面及乘正非仿射 Q4 的闭式矩形参考，XY/RZ 正负相、内部/外边界 zero clip、真实2:1/hanging subface、反向/partition 加回、axis/selected、wrong source/生命周期/失败事务。液气非负且容量闭合、物理体积及两方向 Green moments误差 ≤1e−12实际 capacity；短相/近零非零沿用原相对 gate，无 clipping负值、half 边界值或容差扩大。一般临界界面仍如实拒绝。

实际八组已出版 dual 位置必须调用正式 PhaseRegionSource 与 CchMomentumMass：16个剪裁盒和八个质量构建均成功，液/气矩按独立参考核验；不再只统计旧异常或新面迹。原生面位置八组 publication 尚未迁移，仍应记录原拒绝，不能作为本轮质量成功。再跑 H10i 64普通来源布局、真实 phase/material/dual/KEP 回归及内部 affine；工程/来源 PASS 不关闭 H10。

新几何 clip 账本须记录 `sum(fixed)=exterior_fixed+coincident_interface_moment`，不能强制内部单侧几何迹反号共享或静默重解释旧 assemblyClosure；后续 H10i1e 的全体积账本单独迁移。材料唯一共享 F、内部反号/AMR 子面和、来源一致、KEP/skew 仍是四硬约束，R 不要求为零。

FAIL：来源、归属、闭合矩、真实 mass/clip 或必要回归任一失败，冻结具体原因再修正/拆正式接口，不 fallback、不冒称全链完成。交付 API/实现、同一既有短契约、实际体积/质量输出与独立参考、限制及平面 manifest。下一 H10i1c 单独迁移正规 native 零边拓扑；H10i1d/e/f、H10j 和最终 PDE 按序，H11/H08a 保持等待。
