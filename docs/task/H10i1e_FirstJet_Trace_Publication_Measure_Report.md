# H10i1e / FirstJet 同源 trace 派生测量与实际发布

状态：**PASS / engineering PASS**。本轮只关闭 H10i1e。Pure FirstJet 原 native-face 8 个场景已能实际发布，同一返回 Snapshot 的 complete topology、实际 half cache、两方向 Green 账本与 phase volume 一致；未改变 primary/闭包值/clock，未构造新材料 F。原 8 个 native-face H07 消费仍在 `phase` 阶段拒绝，下一唯一任务 [H10i1f](H10i1f_Native_Trace_Phase_Consumers_TODO.md)。H10 尚未关闭，H11/H08a 未启动。

[任务书](H10i1e_FirstJet_Trace_Publication_Measure_TODO.md)、[实施前设计](H10i1e_FirstJet_Trace_Publication_Measure_Design.md)、[具体身份/边界核对](H10i1e_Trace_Publication_Scope_Note.md)、[API](../api/H10i1e_FirstJet_Trace_Publication_Measure.md)。检查根目录 `runs/checks/H10i1e/20261004_114041_022`；冻结位置 `runs/baselines/H10i1e_20261004`。标准 `local-benchmarks` Release 构建，最终运行均已终止并记录退出码。以 `checks/final2/`、`checks/runtime_final2/`、`checks/compiled_source2/`、`checks/LastTest_final2.log` 为最终证据，前一轮 `final/` 是追加实际 axial ledger 观察之前的 preflight，不能替代最终来源。

FirstJet 发布在最终 Data 内先固定 domain，再绑定原 trace topology/cache，随后用显式 `integrateHermiteGalsNativeTracePhaseVolumes` 测量。`discoverNativeTraceVolumes` 为独立只读入口；`measureNativeTraceVolumes` 复用原缓存，且在同一 scalar allocation 上重新认证公开 topology 的完整性，检查全部 phase/component/graph/side inventory，拒绝手改 mask 或遗漏 component。同步借用 alias 不复制 Snapshot。ConservativeVolume、bubble 和旧敏感度继续原入口，没有失败 fallback。

zero support 没有固定 half phase。shared fixed 面只包括原 ordinary negative ranges，一次积分后给 incident CV 反号；原 Core side moment 只计入真实可用液侧，canonical geometry owner 可在气侧。缺少液侧的外部/selected 图保留几何，free liquid moment 为零，不制造外部液体。原完整 normal zero-set 证书认证 leaf 内部相；纯内部相用实际 capacity/0，一般 cut graph 保留原 local Green/controlled/root integration。actual assembled axial/radial ledger 与 phase volume 分开保存并核验，没有体积拟合/库存/质量修复。

| 验收项 | 最终证据与结果 |
| --- | --- |
| 原 16 research 来源场景 | **16/16 actual FirstJet publication，352 个实际单元**；scalar/revision/time/cache 一致，无 M；原 native 8 组 **40 side graphs / 80 zero halves**。独立 planar phase volume 与 assembled 两方向 Green/capacity 最大相对误差 **8.224464630855156e−16**，原 1e−12 actual capacity gate 保持 |
| 实际 native publication 短契约 | **8 组** XY/RZ×uniform/真实 2:1×正负相；使用正非仿射 Q4；source/identity/clock、原 center、raw endpoint/中点与 quadrature、气侧 canonical owner、actual hanging subedge、standalone/retained 测量一致。失败保留 candidate，成功 consume；foreign scalar、手改 phase certificate 拒绝 |
| 扩展实际 publication | **16 组**通过：2 ordinary→side 拼接、2 multi-component、4 XY/RZ selected 及 physical exterior 等；**24 个 side graph 没有可用液侧**，体积未被外部延拓污染；总参考误差最大 **5.551115123125783e−17** |
| 原未认证边界 | **4 个 RZ axis/domain corner** 保留原 Core 的 undeclared-boundary 拒绝，candidate identity/geometry/clock 保持；不计作支持范围，未修改 Core。原 critical/allzero/high-order normal 拒绝及 FP 契约由原 root/cache 短契约保持 |
| 独立两方向测量短契约 | **16 layouts / 400 actual CVs / 168 zero halves / 84 side graphs**；两个 normal axis、XY/RZ、正非仿射 Q4、uniform/真实 AMR、正负相，original fixed+free 与独立 rectangle/shell 参考通过；最大相对误差约 **8.22446e−16**。完整 discovery 与 retained cache 结果一致；漏 component、legacy zero 入口仍拒绝 |
| 原 topology/clip/cache CSV | **8/8 逐字节相同**：clip、mass、typed trace port/ranges、Core topology/support、half cache/summary；仍保留原 **16 clip boxes / 432 dual mass rows**。只有 actual capability 库存 `trace_contract.csv` 的 native publication 由拒绝变成成功并在 H07 phase 明确拒绝；新增 publication 两表不修改旧参考/算法 |
| 原普通 source/material 来源 | **64 layouts / 157,810 queries，22/22 CSV 逐字节相同**；原物理/reference gate 保持。连同上项共 **30 原 CSV 相同** |
| 最终 native 内部 affine | **3,760 行逐字节相同**；内部最大误差 **1.584919338725414e−13**，原 AMR dof 911 反例误差 **7.41629e−14**，原 gate 保持 |
| 原 flux / KE / native moments | **82 行旧 flux 数值保持**；KEP SUMMARY/NATIVE_MOMENTS、CAP、typed trace、clip/mass、Core side topology 和 half cache 原数值行保持 |
| 标准 Release 必要回归 | **17/17 PASS**，最终实耗 **226.82 s**：architecture/self-test/layout、root/shifted、publication/phase-region/flux、FirstJet 四模式、H07/material/momentum/dual/KEP |

范围：**5 production / 2 既有 tests / 1 既有 research driver，0 新工具/target**。生产文件仅 InterfaceCellMeasure header/cpp、HermiteGalsPhaseVolume header/cpp、HermiteGalsPublication.cpp。测试文件为既有 publication 与 phase-region 契约；原 8 拒绝事务实际在 publication 文件，任务书初稿误列 first_jet 文件，实施前已纠正且保留 FirstJet 四模式回归。source 目录 **1,042 项，1,034 未变、8 变更**；input **346 项，344 未变，仅 driver/README 变更**，目录成员保持。Core/Q4/root/FP/regularity、primary/transport/clock、Hydro F/R/uR/KEP/压力/散度、App/字典/PDF 均保持。

最终 research driver SHA256 `91447E04362F2DE84011614263B764A9FBB7F22E2E70CFF302A7A9BE77A80560`；research executable SHA256 `F9DAEB3C5A54C1E0777B6D572226B6BDBE3D63F2F35CCB2C70DFDC55C7B82529`。四个最终 research mode 的实际 exit/source/binary identity 均核对；实际编译变更、最终测试可执行 hash、原运行与 preflight 失败日志均保留。`verify_saved_data.ps1` 从保存数据独立核验，`verify_scope.ps1` 校验目录/变更预算/被动观察/保护源；flat manifest 覆盖完整交付与证据，并重新核对 d 的 189 archived payload 和 c 的 175 archived payload。实施前设计 hash 保持，旧冻结历史未改写。

限制：当前 native 面可发布，不代表 H07 face measure/rate、RegionSource 或 native dual mass 完整消费已迁移。实际 8 native scene 在旧 H07 `halfFace` 返回 `coincident trace has no legacy phase view`，明确冻结到下一 H10i1f；未把诊断空间迹当作因果材料 F。旧长研究 H10b 未重跑，其 14 失败保持，内部 affine 通过不能代替 PDE 收口。完整流体/电喷雾与低黏/无黏联合验收仍未完成。

材料只有 **唯一共享 F、内部反号/真实 AMR 子面和、来源一致、KEP/skew** 四类硬约束。quarter/压力到 dual 精确交换与 R=0 仅诊断，真实 R/uR 保留；Pure FirstJet 几何完全独立，无守恒库存/bubble。按 H10i1f→H10i 来源链复核→H10j 与材料/最终 native 收口顺序继续；最终内部 affine 先通过，再原完整 H10b。
