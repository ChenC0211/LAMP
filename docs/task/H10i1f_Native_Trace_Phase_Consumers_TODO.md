# H10i1f / 原生 trace 的 H07 与 region 消费链

状态：**PASS / engineering PASS**，2026-10-04；[验收报告](H10i1f_Native_Trace_Phase_Consumers_Report.md)。原16组 publication/H07/真实dual mass全部成功，17/17回归、原普通来源/affine/flux保持；4 production/2既有tests，0新tools/targets。下一唯一任务[H10i2独立来源/设计复核](H10i2_Spatial_Material_Source_Reaudit_TODO.md)，本轮不启动H10j。实施前核验 H10i1e 的 268 frozen payload / 79 canonical delivery 与全部前置归档通过。[实施前设计](H10i1f_Native_Trace_Phase_Consumers_Design.md)保持原字节。

目标：使已实际发布的 Pure FirstJet native-side 几何进入 H07 同源面测量与 RegionSource/真实 native dual mass 消费，闭合来源链；不在本轮构造新的材料 phase rate 或 convection。

允许修改：`include/Modules/Hydrodynamics/CchPhaseGeometry.h` 及对应 `.cpp`、`include/Modules/Interface/HermiteGalsPhaseRegion.h` 及对应 `.cpp`（4 production，最多 5；任何额外文件先明确登记/必要时拆分），既有 `tests/contracts/hydrodynamics/cch_phase_geometry_contract.cpp` 与 `tests/contracts/hydrodynamics/cch_momentum_geometry_contract.cpp`（2 tests）；既有 research/docs 可更新。0 新工具/target。开始前核验 H10i1e 平面 manifest、保存实施前设计/目录/原 gate。

禁止修改：Core Q4/root/trace/topology 及其预算，FirstJet primary/phi/grad/center/Q4/transport/clock，H08/H09 F/R/uR、最终 convection、压力/散度、App/字典/PDF及原容限。无 inventory/bubble/体积拟合/质量修复、source perturbation、half 分配、AMR 特例或失败 fallback。

数学约束：每个只读测度由相同 publication Snapshot/revision/time 得到。native 面的无速度几何测度必须明确普通/重合类型；内部重合诊断可使用原两侧对称空间迹，外部使用原可用内侧并记录 kind，不能被当作 material causal rate。旧面常数 rate 不得用这种诊断伪造新的 F；未有明确因果/外部输入的非零速度重合 rate 仍拒绝。实际 wholeleaf 纯内部证书即使有 attached graph 也适用；region box 的内部体积不重复添加外部侧界面。shared/AMR subedge 与同源 clip/root 来源保持。

实现要求：先核对 H07 build 对 legacy halfFace 的调用与 RegionSource 的旧 discovery，显式迁移 FirstJet native trace。保留默认 historical Conservative/legacy scope。真实 region/dual 交集使用原完整证书或原 regular local integration，不能把 attached graph 送到尚未支持它的 sensitivity/region integrator。face kind/可用侧保持 typed，并使旧 rate 明确处理/拒绝未迁移类型。

T0：标准 Release、architecture/self-test/layout、原 FirstJet 四模式、publication/phase-region/root/half-cache/flux/H07/H08/H09/KEP 必要回归；原 64 普通来源布局、22 CSV、内部 affine、clip/mass 参考保持。冻结历史 PASS/FAIL 不改判。

T1：原 16 research 场景均实际完成 publication/H07 geometry/真实 dual mass，尤其原 8 native 面不得停在旧 phase legacy view。独立 planar XY/RZ uniform/真实 2:1、正负相 volume/mass 与真实 dual box 参考 ≤1e−12 actual capacity/material scale；三类 clock、原 Snapshot 地址、actual subedge、面 kind/两侧 availability 逐项验证。ordinary→side、多 component、轴/selected/真实边界、missing exterior、wrong source、FP、lifetime/失败事务核验。原未认证 RZ 域角点/critical/allzero 保留拒绝。

失败规则：任一来源/类型/regularity/非负容量/回归失败，冻结具体反例后修正或拆任务。必要链若需修改超范围 Core 或材料 owner，不绕过 source/gate；不扩大容限或隐藏拒绝。

交付：实施前设计/API/报告、实际 H07/region/dual mass 及独立参考、标准短契约、实际编译来源/二进制/退出/CSV、完整 source/input 范围核验与 flat SHA256 manifest。只关闭 H10i1f，下一项为独立 H10i 来源链复核，再 H10j；原 native affine 先通过再原完整 H10b，H11/H08a 等待。

材料仅四硬约束：唯一共享 F、内部反号/真实 AMR 子面和、来源一致、KEP/skew；quarter/压力到 dual 精确交换与 R=0 仅诊断，继续真实 R/uR。Pure FirstJet 几何完全独立。
