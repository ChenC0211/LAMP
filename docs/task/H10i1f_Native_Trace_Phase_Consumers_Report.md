# H10i1f / 原生 trace 几何消费链验收

2026-10-04。结论：**PASS / engineering PASS**。本轮只接通已发布 Pure FirstJet native-side 几何到 H07 空间面测量、RegionSource 和真实 native dual mass；4 production / 2 既有 tests / 1 既有 research driver，0 新工具/target。[任务书](H10i1f_Native_Trace_Phase_Consumers_TODO.md)、[实施前设计](H10i1f_Native_Trace_Phase_Consumers_Design.md)、[API](../api/H10i1f_Native_Trace_Phase_Consumers.md)。下一唯一任务是 [H10i2 独立来源/设计复核](H10i2_Spatial_Material_Source_Reaudit_TODO.md)，本轮不执行 H10j。

材料继续仅四类硬约束：唯一共享 F、内部反号/真实 AMR 子面和、来源一致、KEP/skew。quarter、精确 pressure→dual commuting、严格 R=0 只作为诊断，保留真实 R/uR。Pure FirstJet primary/Q4/transport/clock 完全独立，没有守恒库存、bubble、几何拟合或质量修复。本轮没有构造新材料 F 或修改最终 convection。

## 正式消费路径

H07 按原 publication `traceHalfFace` 构造普通/内部重合/边界重合类型，保留 selected domain 的 actual incident cell 和可用侧。普通 half 的积分算术保持；内部 zero 部分提供两侧对称**空间诊断**，边界仅提供实际内侧，不伪造外部相。每个实际共享 native/AMR subface 一份测度，原 128 eps face closure 不变。

原 `CchPhaseRate::evaluate` 明确拒绝重合面任意非零 actual normal velocity，精确零速度率为零。短契约包括 ±1e−300 的非零值。没有把空间对称迹当作新因果材料率；H10j 才处理同一 advector 的 causal inflow 与外部来源。

RegionSource 为 FirstJet 使用原完整 native trace topology，historical ConservativeVolume 保持原入口。side graph 的实际 incident leaf 必须具备同源完整纯内部证书；graph 保留拓扑/来源，避免送入旧 regular integrator。真实 rectangle/dual/overlap 用 capacity/0 或原 ordinary local clip/root/Green integration，无重复外部侧体积。`CchMomentumMass` 算法文件未改，它原有调用现在直接消费正式 RegionSource。

## 实际来源与独立参考

最终 `trace-contract` 仍使用原 n4 XY/RZ、uniform/真实 2:1、dual 内部/native 面、正负 phi 的 **16 个原场景**，没有改 source、identity 次序或原数值路径。全部 publication、H07、真实 dual mass 成功；上一冻结版本的 native8组在 H07 legacy view 拒绝，本轮关闭这一消费缺口。

| 实际记录 | 数量 / 裁决 |
| --- | --- |
| 最终 publication / CV | 16 / 352；原 publication 表逐字节保持 |
| native side graphs / zero halves | 40 / 80；原 Core/cache 表逐字节保持 |
| 实际空间 face / coincident face | 864 / 40；原两侧/物理平面参考通过 |
| whole RegionSource | 16；两方向 Green 与原 Snapshot/identity/time 通过 |
| 实际真实 dual mass | 864；液/气体积和 rho=1000/1 独立平面参考通过 |
| 两方向 clip box | 32；内部侧参考与非负容量通过 |
| 实际 typed segment | 2,816；来源与原 port 认证通过 |

新保存 `spatial_face_trace.csv` 与 `native_region_source.csv` 只观察正式输出。原六份 publication/Core/cache CSV 字节保持；clip/mass/port/ranges 四份表因 native8组成功而扩展，原 ordinary8组的 header/全部原始行完全保持。旧 legacy segment 查询仍在 zero trace 上拒绝，这是原 legacy 接口能力记录，不是新工厂失败。

## 短契约与回归

既有 phase test 新增 48 个非仿射 Q4 native publication（XY/RZ、uniform/AMR、两方向、正负、Initial/Midpoint/End），8 selected 和 6 physical boundary；3,168 face、316 coincident face、632 非零 causal rate 拒绝。原可用侧、selected source 额外支撑/部分支撑、axis、ordinary→side、多 component、完整 time、FP 失败事务和 publication 句柄释放后 source lifetime 通过。

既有 momentum test 新增 48 同源 publication、288 真实 mass 构造（full/selected、rho=1/1e3/1e6）、12,024 true dual、19,872 actual overlaps、144 region box、240 side graphs。独立体积/材料质量最大相对误差 **1.26e−16**。另 20 个 exterior/missing-side、ordinary→side、双 component 场景验证 region 两方向及外部真实 dual mass，误差 **1.12e−16**；失败 region 查询后原来源保持。原未认证 RZ 轴/域角点四个拒绝及 critical/allzero regularity gate 保持。

标准 Release 构建与 **17/17 CTest** 通过（architecture/self-test/layout、FirstJet 四模式、publication/root/Green/region/flux、H07/H08/H09/dual/KEP）。原 64 普通布局、157,810 查询和 22 CSV 逐字节保持。内部 affine 3,760 行原表保持，最大误差 **1.5849193387e−13**，原 XY AMR n16 dof911 反例误差 **7.4162898045e−14**。原 82 行 flux 数值、KEP/CAP/native topology/cache/clip/mass 短结果保持。

## 范围与冻结

1,042 源文件目录成员保持，仅 4 production / 2 tests / driver 共7项改变；346输入成员保持，仅既有 driver/README改变。被动 face/region 观察段删除后，research全文与实施前相同。原 Core、primary/transport/clock、H08/H09 F/R/uR、KEP/convection、压力/散度、App/PDF保持；无新生产变体或工具。

实际最终编译源、research exe/DLL、测试 exe hash、四模式退出与 CSV、CTest 日志、实施前目录/设计和独立数据/scope核验保存至 [冻结证据](../../runs/baselines/H10i1f_20261004/README.md)。[flat manifest](../../runs/baselines/H10i1f_20261004/delivery_manifest.json)覆盖所有 delivery/evidence；[独立核验](../../runs/baselines/H10i1f_20261004/verification.json)检查 current canonical 和 e/d/c archived payload。旧历史冻结结论不改判。

本轮几何消费链 PASS 不代表 H10i 的空间设计已独立重新验收，不代表材料 causal owner 或 PDE 已实施。H10i2 先复核，随后 H10j/H10k/H10l 分别推进。原完整 H10b 的 28窗口/15,116步/14失败保持，尚未重跑；最终 native 内部 affine/空间阶先通过，再按原 H10b gate 重跑。H10 未关闭，H11/H08a 等待，长期 goal 保持 active。
