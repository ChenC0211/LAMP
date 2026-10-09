# H10i1c / 正规 native 零边拓扑验收

2026-10-04，**PASS / engineering PASS**。本轮只关闭 H10i1c；goal 保持 active/progress，H10 和原完整 H10b 未关闭，完整流体/EHD 未验收。[TODO](H10i1c_Native_Zero_Edge_Topology_TODO.md)、[实施前设计](H10i1c_Native_Zero_Edge_Topology_Design.md)、[API](../api/H10i1c_Native_Zero_Edge_Topology.md)。下一唯一任务为 [H10i1d 只读半面 trace cache](H10i1d_Readonly_HalfFace_Trace_Cache_TODO.md)，此轮未启动它。

开始前独立核验 H10i1b 的 **219 archived payload / 54 当前交付**，保存实施前设计及 1,042 source / 346 input 路径与哈希。检查目录 `runs/checks/H10i1c/20261004_103851_399`；冻结目录 `runs/baselines/H10i1c_20261004` 保存最终源码、实际研究二进制/依赖、退出记录、原始 CSV、完整回归日志、范围及数据验证脚本。平面 manifest 逐项核验 archived payload 和当前 canonical delivery；先前归档的已替代源码只核验 archive。

新增显式 `discoverNativeTraceComponents`。原默认 topology 和 legacy sensitivity 保持恒零拒绝；新路径用原 Q4/retained root 引擎认证实际 zero subedge，在每个 incident leaf 同时要求完整零集的统一有符号 normal-partial 证书与原边限制的严格一阶非零、无根。完整证书排除隐藏内部分支；不凭边上采样推断整个 leaf。

每条真实子边只由一个实际 incident cell ID 确定的 canonical geometry owner 输出。`NativeSideSupport` 只读并 owns 同一 Snapshot，保留实际 native edge/vertex、normal side、两侧 cell/phase/availability；不伪造恒零多项式的 finite root ordinal。ordinary→side 以原 actual vertex incidence 拼接，原 valence、正则性、axis/声明边界和最少 graph 规则保持。液侧负片和 geometry owner 可不同；inactive/外部 unavailable 不外推，不制造外部液体库存。

sampling/frame 使用同一 rational native cell chart；仅附着支撑的 normal coordinate 固定为精确 side 0/1，tangent 保持严格内点。原 raw scalar jet 给出法向，普通 graph 的 strict interior gate 不变。Green 几何矩取已认证直线的 exact retained endpoint primitive，XY/RZ 测度保持；旧 sensitivity 和未迁移的 region consumer 明确拒绝 attached graph，没有修正 Jacobian、材料通量或另一套 scalar。

| 验收项 | 实际结果 |
| --- | --- |
| 原 16 场景的独立 seed Core 查询 | **16/16**；8 个 native-face 场景获正式 Core typed 拓扑，共 **40 条 side supports**；原 dual 内部 8 组的 graph/endpoint/frame/geometry moments 与旧 Core 逐字段相同 |
| 独立平面 XY/RZ length/physical measure/两方向 Green 矩 | 原 **1e−12 实际尺度** gate；最大研究误差约 **4.24074e−16** |
| 既有 root 短契约扩展 | **16 布局 / 80 支撑 / 32 actual hanging 子边**；正负相、两种法向、正非仿射 Q4、uniform/真实 2:1；原 frame chart round trip、完整 normal 证书、子弧 partition、actual edge/vertex、canonical owner 通过 |
| 拓扑扩展 | 普通→side 同一链 **2 side / 2 ordinary**；两个 side components 共享一个真实 negative region；physical/selected boundary availability 正确，高阶退化/crossing/全零拒绝；wrong source/cell/public addressing、legacy sensitivity、owning lifetime 通过 |
| FirstJet 失败事务 | **8/8**；旧出版仍拒绝 zero edge，candidate identity/phi/grad/edge jets/无 M 主状态和时钟保持；Core seed 不是 accepted publication |
| 原已出版 dual 内部来源 | **8/8** publication/phase/mass；**16 clip / 432 mass**，五份 trace/clip/mass/typed range CSV 与 H10i1b 逐字节相同 |
| 原普通来源与材料诊断 | **64 layouts / 157,810 queries**；两模式 **22/22 CSV 逐字节相同**；原独立物理/reference gate 保持 |
| 最终 native 内部 affine | **3,760 行逐字节相同**；内部最大 **1.58492e−13**，原 AMR dof 911 反例 **7.41629e−14** |
| 标准 Release 受影响 CTest | **16/16 PASS**；原 KEP SUMMARY/NATIVE_MOMENTS、32 CAP、NATIVE_TRACE、COINCIDENT_CLIP/MASS 逐条相同 |

范围固定 **5 production / 2 既有 tests / 1 既有 research driver**，0新生产工具/target。五份生产文件是 public MappedInterface/InterfacePanels、两个对应 detail 头和原 strict geometry TU；测试只扩展 root enclosure 与 Hermite–GALS publication。原 research trace 观察原样保留，新 independent seed Core observer 放在原 16 场景记录之后，避免改变原 identity 次序。`trace_topology.csv`、`native_side_support.csv` 明确标记 seed 来源，不是 Hydro acceptance。

最终 source catalog **1,042：1,034 不变、8 变更**；input **346：344 不变，仅既有 driver/README 变更**，成员保持。App/PDF、字典、Hydro F/R/uR/KEP、Interface primary/Q4/transport/clock、原研究数值函数与门槛均保持。实际 final driver SHA256 `A36FF8AE2EEDCDC30F4E75BABF9771BC71B0EC7834E43621E01DF05DBC9D1D8E`，research exe `089CEBF12F4D18EBC542F06986C0B691BC7E5A0DC3405C624BE69017CC431507`，最终四模式均对应此二进制。

保留初次编译/短契约 preflight 日志：一个 legacy 拒绝检查初次插入位置错误已修正；fixture 首次使用错误字段名、退化拒绝字符串与 FirstJet config 字段名随后核对现有 API 修正；一次已有 target 名称写错后使用标准已注册 target。未扩大数值容限或改变退化接受条件。最后补充 wrong cell、完整 normal 证书、ordinary 新旧入口逐字段比较及真实 edge ID 直接复用后，重新标准构建并运行最终回归和四个研究模式。

限制：本任务仅验收 Core native zero-edge topology/frame/sampling/Green moments，projection/load/remap 未独立验收。半面 typed cache、FirstJet 派生 volume/publication 与 H07 消费仍依次归 H10i1d/e/f；旧 publication 的 native **8/8 拒绝仍完整保留**，不能据 Core PASS 宣称整个来源/质量链可用。无 material F、空间 convection factory 或时空输运迁移，原完整 H10b 没有重跑、历史 14 项失败保持。来源链闭合并复核 H10i 后才 H10j；最终 native 内部 affine 先通过，再原 H10b，H11/H08a继续等待。

材料耦合硬约束仍只有**唯一共享 F、内部反号/真实 AMR 子面和、来源一致、KEP/skew**。quarter/pressure→dual 精确交换和 R=0 只作诊断，继续实际 R/uR；Pure FirstJet 完全独立，无守恒库存、bubble、逐 cell 几何体积拟合或质量修复。
