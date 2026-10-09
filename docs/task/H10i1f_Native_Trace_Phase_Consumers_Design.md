# H10i1f / 实施前设计

日期：2026-10-04。实施前 H10i1e 268 payload / 79 canonical delivery、d 189 与 c 175 archived payload 核验通过。原 source/input catalog、4 production / 2 tests / driver 的实施前快照、原 TODO、App/PDF hash 保存在 `runs/checks/H10i1f/20261004_121139_996`。本设计在数值改动前保存副本，后续不追改。

唯一事实：实际 Pure FirstJet publication 的原 native-side 几何进入 H07 空间面测度与 RegionSource/真实 dual mass。新材料 F、空间因果率、H08/H09转移与 convection 不在本轮实现。

H07 面测度显式读取 publication 的原 `traceHalfFace`。普通 half 仍用其原 `ordinary()` ranges 与原物理积分算术。zero half 使用同一 typed query 的两侧 normal phase ranges，按 Hydro 所选 domain 的真实 incident cell 过滤 availability；publication scalar 有额外活动支撑时，也不能把所选域外的 side 当作 inside-domain 相。每个共享 face 只测量一次，实际 AMR subface不合并。

新增 `CchPhaseFaceTrace`：保留 kind（Ordinary / InteriorCoincident / BoundaryCoincident）、原 selected negative/positive cell、每个可用侧的完整物理 liquid/gas trace 与 zero half 数量。内部重合诊断测度为两侧空间迹的平均；外部/selected 重合诊断仅用实际内侧。此平均是明确定义的空间诊断，不能当作 half phase 或 causal material rate。旧 `CchPhaseRate::evaluate` 若任一重合 face 的原实际 normal velocity 非零，显式拒绝；velocity精确为0时率为0，不依赖诊断相分配。没有外部边界数据或新空间rate factory被暗中创建。原 face closure 128 eps gate及旧 rate identity不放宽。

RegionSource 按原 primary transport 明确选择 discovery：FirstJet调用原 complete native trace topology，ConservativeVolume保持旧入口。保留原 Snapshot 的 owning geometry view alias，原 component/source/clock不变。native-side graph所关联的每个实际 incident leaf必须有同源完整纯内部相证书；graph留在完整 topology 中用于连通性/来源统计，不传给旧 regular region integrator。实际任意 rectangle/dual交集的纯相叶片用其原capacity/0，普通cut仍用原 local trace/root/controlled Green积分、严格非负/capacity、误差预算与coverage。没有额外界面体积项重复添加到 region内。

`nativeSideGraphs()` 为prepared region来源提供只读数量，便于验证正式消费使用的是原 trace topology；不是求解分支、material inventory或第二几何。CchMomentumMass自身算法不修改，它现有调用RegionSource即可得到真实native dual/overlap质量，必须由实际运行证明。

预算固定4 production：CchPhaseGeometry header/cpp、HermiteGalsPhaseRegion header/cpp；2既有tests：cch_phase_geometry_contract、cch_momentum_geometry_contract。既有research只增加空间face kind/侧availability、真实region/dual来源的被动证据，不修改原numerical driver。0新target/tools。必要新文件或Core修改先拆分，不能以fallback越过预算。

验收：原16 research场景实际publication/H07 geometry/dual mass全部成功；独立平面两侧face与region/真实dual/overlap质量≤1e−12实际尺度。短契约覆盖XY/RZ、uniform/真实2:1、正负相、正非仿射Q4、Initial/Midpoint/End、ordinary→side、多component、axis、selected/physical boundary、missing exterior、wrong source/time/domain、FP、lifetime/失败事务。原未认证RZ域角点/critical/allzero保持拒绝。标准17回归、原64 source布局/22 CSV/157810 queries、内部affine、clip/topology/cache、旧flux/KEP/CAP结果保持；原native8 mass由失败变成功的证据单列，不改历史冻结文件。

完成只关闭H10i1f，下一任务独立H10i来源链复核，然后H10j。仍仅唯一共享F、内部反号/真实AMR子面和、来源一致、KEP/skew四类材料硬约束；quarter/精确交换/R=0仅诊断，真实R/uR和Pure FirstJet独立性保持。原完整H10b未重新验收，H11/H08a等待最终native空间/PDE通过。
