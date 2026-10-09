# I02 / 回溯映射与同公式 Jacobian

任务依据：[05 路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md) I02、[数值约束](../01_LAMP_Numerical_Invariants.md) I-R4、冻结 Hermite–GALS 方法 §4.1（式 18–19）。执行日期：2026-10-01。前序 I01 已 PASS，开工时复核其交付文件及归档。

## 目标

只解决一个数学事实：V2 tracer 同时返回实际离散 `B_h` 及 `J_Bh`，Jacobian 来自同一 midpoint RK2 公式的 tangent-linear differentiation，并与实际映射的独立数值微分一致。

## 允许修改

- 新建 `include/Modules/Interface/HermiteGalsTrace.h`、`src/Modules/Interface/HermiteGalsTrace.cpp`。
- 新建 `tests/interface_hermite_gals_trace_contract.cpp`，在 `tests/CMakeLists.txt` 注册该契约。
- 本 TODO、`I02_Characteristic_Map_Jacobian_Report.md`、`docs/I02_Hermite_Gals_Characteristic_API.md`、`docs/README.md` 的任务索引。
- 证据目录 `runs/baselines/I02_20261001/`，以及正常构建产物。

预算：2 个 production 文件、2 个测试文件。

## 禁止修改

- 不修改 I01 主状态或 Core 的既有回溯公式；不扩展 S08 冻结方法、Hydro、Electricity、Apps、reader 或 cases。
- 不开始 I03/I04 的 skeleton 输运、I05 的混合导数闭合或相体积/AMR 算法。
- 不用空间有限差分构造生产 Jacobian，不加入速度 Hessian、独立梯度拟合、边界裁剪/反射或旧方法 fallback。
- 不把未定义的边界命中时间导数包装为内部 Jacobian，不将 tracer 验收视为 Gate I 完成。

## 数学约束

- 固定使用通用 Core `LevelSetBoundaryTrace::traceWithJacobian` 的显式 midpoint RK2；这是冻结方法允许的 integrator，不另外复制第二套回溯 stencil。
- 一个速度采样同时返回 value 和空间 derivative；二者必须由同一不可变插值分支求值，时间插值权重相同。
- 完整子步 `x_m=x-h*u(x,t)/2`，`B=x-h*u(x_m,t-h/2)`；`J=I-h*Du_m*(I-h*Du_start/2)`。多个成功子步按实际轨迹顺序组合。
- Jacobian 行为 departure 分量，列为 arrival 坐标；它是固定插值/子步分支上的导数。
- 边界事件保留 side、命中时间和无 Jacobian 的显式状态；没有 boundary-history tangent 时不得做 first-jet pullback。
- 速度 view 自带 geometry identity 和时间层，两个端点的物理时间必须与请求的旧/新（或中点）区间一致。结果保持这些来源信息。

## 实现要求

- [x] 复核 I01、六份规范及冻结公式，记录开工前文件哈希。
- [x] 增加只读的 timed velocity view、V2 tracer 和结果接口；owned callback 的捕获对象必须由 provider 保证寿命/不可变。
- [x] 从传入物理时间推导实际 trace interval，不独立接受可能矛盾的 dt；拒绝无法分辨的中点时间层。
- [x] 复用 Core 的同公式 primal/tangent，保留 boundary/active-domain/axis/CFL 检查，错误含模块、阶段和来源身份。
- [x] 独立检查四个 Jacobian 分量和 primal 采样序列；覆盖 nonlinear、两子步、多种 dt、XY/RZ、静止/零步和边界事件。
- [x] 以解析 time-dependent rotation 的长区间组合验证 `B/J` 在均匀和交替步长序列上的二阶收敛。
- [x] 完成编译、架构、I01 回归和旧 Interface 短回归，交付误差/阶数及原始日志。

## 验收

- **T0：** Release 构建 QUANTUM、I01/I02 contract、原 Interface contract、architecture；新头文件可独立作为首个 include，无新增 warning；architecture 与 quantum_help 通过。
- **T1：** affine shear 的解析 `B/J` 误差不超过 `2e-14*(1+|reference|)`；nonlinear 离散 map 的中央差分覆盖四分量、多个位置和 dt，误差随扰动缩小二阶下降，最细扰动误差 `<2e-8`；两个成功子步的有限差分误差 `<2e-9`，primal 和 tangent 路径点值及采样顺序完全一致。
- **T1：** 零步 identity、静止点非平凡 tangent、只读来源时间/身份、四面 boundary hit 无内部 J、invalid velocity/derivative/time/domain/axis/CFL 明确拒绝，且不修改输入状态。
- **T2（tracer 专项）：** time-dependent rotation 在固定总时间、`N=8,16,32,64` 上，uniform 和 alternating dt 的位置/全 Jacobian 最后两级收敛阶各在 `[1.85,2.15]`。不能将该结果推广为完整 Interface 输运阶。
- **回归：** I01 state contract 和 `quantum_interface_state_contract --parametric-state-control` 通过；除测试注册外的原有源码/输入及 69 个 S08 专有实现哈希不变。

## 失败规则

- 如果 Core 实际公式不满足 I-R4，任务判 `FAIL`；先明确问题与所需最小修复范围，不把错误包装成通过。若必须跨模块或超预算，记 `BLOCKED` 并说明需拆分内容。
- 边界命中需要 boundary-history 导数时显式报告缺口；本轮不发明边界物理、不静默切换 RK4 或其它方法。
- 不以放宽验收容限、调整物性或改动旧 benchmark 掩盖失败。S01 已知 Hydro 首步失败保留。

## 交付

V2 tracer 接口、独立契约、TODO/报告/API 说明、原始误差和收敛阶、构建与文件身份证据。只做 I02，结束状态为 `PASS / FAIL / BLOCKED`；下一项 I03 留待下一轮。

**最终结论：PASS。** Release 构建、CTest 4/4 与旧 Interface 短回归通过；全部预定 B/J 误差及收敛阶门槛通过。production 2 文件、测试 2 文件。详见 [I02 报告](I02_Characteristic_Map_Jacobian_Report.md)。
