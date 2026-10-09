# H10i1c / 正规 native 零边派生拓扑

状态：**PASS / engineering PASS**；本轮唯一任务已完成，详见[报告](H10i1c_Native_Zero_Edge_Topology_Report.md)和[API](../api/H10i1c_Native_Zero_Edge_Topology.md)。实施前设计见[设计](H10i1c_Native_Zero_Edge_Topology_Design.md)。已前置核验 H10i1b 的219 archived payload/54当前交付，并保存本任务实施前设计、原 gate 和来源目录快照。[总契约](H10i1_Coincident_Phase_Trace_Clip_Contract_Design.md)。本任务只迁移 Core 同源正规零边拓扑及其 frame/geometry moment 读端口，不改 Interface primary/Q4/transport/clock，也不把实际 publication 全链提前算作通过。Core16/16场景、短契约16布局/80支撑/32hanging、8组失败事务及16/16受影响回归通过；普通来源22CSV与内部affine逐字节保持。下一唯一任务[H10i1d](H10i1d_Readonly_HalfFace_Trace_Cache_TODO.md)，H10仍未关闭。

目标：显式 `discoverNativeTraceComponents` 给出正规 native zero edge 的 typed side-attached graph、真实 incident owner、两侧相与端点/vertex/subedge incidence；完整 native leaf 零集与同号 normal partial 证书排除遗漏内部零分支。canonical edge arc 只由一个几何 owner 输出，不伪造无限根的 rootOrdinal。

production 上限五文件：`include/Core/Geometry/LevelSetMappedInterface.h`、`include/Core/Geometry/LevelSetInterfacePanels.h`、`include/Core/Geometry/Detail/LevelSetNativeInterfaceChain.h`、`include/Core/Geometry/Detail/LevelSetMappedInterface.h`、`src/Core/Geometry/LevelSetHermiteGeometry.cpp`。test 仅既有 `interface_root_enclosure_contract.cpp` 和 `interface_hermite_gals_publication_contract.cpp`（2）；既有 research/API/任务文档可更新，0新工具/target。超过 5/2 先拆。证据 `runs/checks/H10i1c/<时间戳>`，最终冻结基线。

禁止变更：legacy discovery/sensitivity 的恒零拒绝、scalar 与 root engine、库存/bubble/拟合/修复、材料 F/R/uR、最终 native convection、D_I/C、App/字典/PDF/原容限。不能只凭零边一阶导数忽略其他内部环/分支；不能复制 Q4、用 epsilon 或节点采样证明全叶相。inactive/selected 支持不外推；一般临界 crossing、全零叶与未认证退化仍拒绝。

T0：标准 Release、architecture/self-test/repository layout、原 root/FirstJet/publication/phase/material/dual/KEP 受影响回归、H10i1b 闭合 clip 与 64 普通来源、内部 affine；原 protected 输入/来源逐项核验。T1：XY/RZ 正负相、regular plane/positive nonaffine Q4、真实 uniform/2:1/hanging native edge、axis/边界/多 component、端点 incidence 与方向/唯一 owner、错源/生命周期/失败事务；原严格正规性及 ≤1e−12 实际尺度 gate，不增加 AMR 特例或无穷 root inventory。

首先核验 H10i 已记录的 native-face 位置八组在新显式 Core topology port 的真实能力，dual 内部八组普通图结果保持。旧 publication 尚无新 typed measurement/cache consumer，若仍拒绝应明确记录，不能以 Core 来源 PASS 冒称整个 publication/质量链修复。

FAIL：来源、完整零集、正规性、去重/端点拼接、frame/geometry moments 或必要回归任一失败，冻结具体反例再修正/拆任务，不 fallback。交付同源 typed API、既有短契约、实际八组能力与限制、独立参考、目录差异核验和平面 manifest。之后 H10i1d/e/f 单独接通 cache/派生出版/正式消费者，再复核 H10i→H10j；内部 affine 先通过再原 H10b。H10 未关闭，H11/H08a 保持等待。

材料硬约束仍只有唯一共享 F、内部面反号/实际 AMR 子面和、来源一致、KEP/skew；quarter/压力到 dual 的精确交换及 R=0 不设为硬约束。几何容量与曲线准确性按各自 gate 验收，不重建守恒库存。
