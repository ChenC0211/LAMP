# I12a10r2 / I06 声明外边界对称性闭合报告

**PASS：完整平面跨面由 36/48 恢复为 48/48，四组 40 步平移与四组静止均保持通过。I12/Gate I 仍 FAIL。** 实现遵守 [TODO](I12a10r2_Declared_Symmetry_Closure_TODO.md) 的 2 production + 2 test 预算。

## 修复

I06 接收与 I07/I10/I11 相同的显式 SymmetryPlane 列表，默认空。声明必须对应真实 active 外边界；重复、非法方向/坐标、重复自动 RZ 轴、非零 vertex/edge 法向 first jet 都拒绝，包括 1e-30 的非零值。

沿原 RZ 轴反射路线，法向紧凑求导与 mixed 兼容系统行施加精确边界约束，随后构造原 exact hanging skeleton，并核查实际边界 normal/mixed。其它系统行、固定双向平均及两个 1e-12 残差检查保持。没有修改 first jets/M、B/J、Core restriction、bubble/topology gate，也没有根据数值大小猜测 BC。benchmark 初始化和 half/full completion 传同一声明；无声明的原始失败诊断保留且输出不变。详见 [API](../api/I06b_Hermite_Gals_AMR_Closure_API.md)。

## 验证

- 30 组非平面 Q4：XY 四外边界、平移缩放域、uniform/两种 AMR、四边同时约束、RZ 自动轴加三外边界。最大 mixed 误差 2.17e-14，AMR C1 缺陷 1.78e-15；所有声明边界的五个 normal Q4 DOF 严格零。
- 8 类非法声明/输入全部拒绝，candidate identity/primary/M 不变。
- 原 48 组 crossing publication 全部真实完成 half/full 发布；最大相对单元库存误差 1.1437146127282951e-15、穿面液量误差 3.0922059403117525e-15、first-jet/法向误差 6.1062266354383639e-16、曲率误差 5.7412822849917378e-14。原先成功的 36 组数值输出完全一致。
- 四组 XY/RZ × uniform/AMR 平移各完成 40 步至 t=.4，四组静止各 4 步。整份 full translation 输出与 I12a9r2 逐字一致；运行约 182.90 秒。
- Release 重建 28 targets，47.77 秒，无 warning；30 项指定 CTest 全通过，194.97 秒。27 项旧完整输出逐字一致；AMR contract 原有全部输出内容不变，仅追加新验证。architecture 仅增加一个 include：12 targets / 380 TUs / 899 files / 3086 includes / 5 legacy edges。

首次聚焦测试中，完整 crossing 已通过，但新增 uniform fixture 调用了只接受 AMR fine side 的 C1 测试辅助函数而失败。修正该调用条件后，uniform 仍核验 Q4 mixed 和全部边界法向 DOF，所有 AMR 继续核验 fine-side C1；未改生产实现或数值门槛。初次失败日志保留。

81 冻结/参考文件、3 PDF 和 I12a10r1 archive 哈希不变；旧 midpoint crossing probe 仍退出 1，输出不变。证据、源码、28 可执行与 1 DLL 的逐项哈希归档位于 `runs/baselines/I12a10r2_20261001/`。

## 剩余工作

本项补齐平面跨面与边界约束，不代表一般界面或非恒定速度已通过。下一项补充非多项式界面的完整常速度空间收敛，以覆盖一般 Q4 数据；然后继续非恒定速度时空源、事件覆盖、旋转/vortex 与其余 I12 验收。q=1/2 正长度零迹及其它既有未支持情况继续保持明确拒绝。不进入 Hydro。
