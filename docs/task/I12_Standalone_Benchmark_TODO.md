# I12 / standalone benchmark

## 任务 ID / 标题

I12：Conservative AMR Hermite–GALS 完整一步与 standalone benchmark。I11 archive 和 10 项交付哈希已核验；上一轮为 progress。

## 目标

检验同一 accepted state 的 first-jet 输运、共享相通量、reflux、bubble 闭合和唯一几何发布能否共同实现静止、平移、旋转、vortex、AMR crossing 与曲率收敛。不能用逐步重置解析几何/M 代替数值时间推进，也不能只用对齐步端的跨面案例宣称通过。

## 允许修改

- 新 `tests/interface_hermite_gals_benchmark.cpp` 与 `tests/CMakeLists.txt`，先建立经公开模块 API 的完整链路及解析反例。
- 本阶段不修改 production 算式；若证据要求新接口或数学修复，先记录失败、提出范围明确的 I12 子任务或设计变更，保持路线图顺序，不能扩大本任务范围修多个数学事实。
- 本 TODO、报告、`docs/README.md` 与 `runs/baselines/I12_20261001/`。预算 0 production + 2 tests。

## 禁止修改

三份冻结 PDF、S08 frozen sources、Hydro/Electricity、App/reader/输出主线；不改旧路径、不加 fallback、M clipping、全局 upwind/KKT，不调整既有验收容限以掩盖链路失败。

## 数学约束

1. 每一步的 vertex/edge master first jet 来自上一份已发布 snapshot 和同一离散 B/J。只在真实外边界提供明确的解析 boundary first jet，不在内部覆盖输运结果；hanging 数据仍由 master 约束生成。
2. M 仅经 I08 共享 face flux 更新，并由 I09 汇总实际子面；禁止更新后几何重算 M。最终 I10 闭合、I11 发布全部通过后才替换 accepted handle。
3. 半步预测采用旧 snapshot 生成私有 quarter-time first-jet geometry，用 I08 做半区间守恒更新；所得 M 与半步 skeleton 经 I10/I11 得到 full-interval midpoint。再以该 midpoint 更新全步 M，并与全步 skeleton 闭合发布。私有 quarter geometry 是阶段预测，不是 accepted 状态，也不对外供力/曲率；其暂存 M 不参与通量求积。
4. 所有阶段保留真实时间和 geometry identity；速度/face polynomial 对应同一解析场及时间，RZ 采用实际度量。预测、闭合或 publication 失败不能修改 accepted state。
5. 位移、gradient、curvature、M 的误差分别记录；守恒 closure 与相对解析误差不互相代替。固定 AMR crossing 与动态 refinement/coarsening 分别核验，缺少 regrid 证据不得声称完整 Gate I。

## 实现要求

测试 harness 只负责编排既有模块，不复制 tracer、mixed closure、FV 更新或 bubble 算式。center prediction 从同一旧 scalar 经 I02 回溯获得。初始解析 first jet/I07 库存用于建立唯一初始状态，之后不重置。全部六类 benchmark 和独立解析跨面反例保留结果；一项失败也不能通过忽略用例或降低标准将 I12 标 PASS。

## 验收

- **T0：** Release 编译、architecture、前序 I01–I11 回归；无新增 warning，冻结文件和改动预算核验。
- **T1：** 原 state 不受失败影响；flux 与 midpoint 同源、内部 cancellation/coarse–fine residual 小于既有容限；零速度 4 步 first jet/M 无系统漂移（scaled max `<=1e-11`）。
- **T2：** circle/sphere 的 static/translation、off-center rotation、可逆差分旋转 vortex 和 2:1 AMR crossing，至少记录 position/gradient/curvature/phase-volume 的最终误差、accepted time、拒绝阶段。所有非静止长程运行必须实际完成目标时间。
- **T3：** 至少三档 h 和 dt；空间 curvature/gradient 必须递减，末档阶数分别 `>=1.5`/`>=2`；独立时间细化位置误差末档阶数 `>=1.8`；M closure `<=1e-10`。舍入底限须以误差数据说明，不伪造阶数。
- **跨面反例：** 解析平面在任意步内穿越 face，检查原始 midpoint flux 与严格容量条件/全步 skeleton 可达 M；不通过时必须报告原始误差、CFL 与理论原因，不移到恰好步端以求通过。
- **Gate I 完成审计：** 六类 benchmark、AMR restriction/conservative coarsening、boundary/axis、长时强变形与旧法可比 A/B、配置/诊断接线逐项查证；缺少项保持未完成，不进入 H01。

## 失败规则

数值失败为 FAIL；缺接口或冻结方法与现有离散存在不可调和冲突则记录 BLOCKED 的具体依赖和可复现证据。若需改变冻结方法，必须先提出设计变更说明；不得自行改公式。task 状态与整个 goal 状态区分：本任务未过不代表总体目标完成，也不自动把 goal 标 blocked。

## 交付

完整链路/反例的可执行源码、每项原始通过/失败数据、前序回归、源码/runtime 哈希、TODO/报告与下一步的具体修复依据。I12 全部验收之前 Gate I 保持未完成。

## 本轮结论

**FAIL。** 完整一步驱动已执行；静止和静态曲率收敛通过，移动界面存在 phase flux/bubble 可达性失败。独立解析平面在 CFL=0.08 时复现上游库存缺口或 2% 超容；I08/I10 的拒绝保持 accepted state 不变。生产算式不改，已形成 [设计变更说明](I12_Phase_Flux_Time_Integration_Design_Change.md) 等待确认。详见 [报告](I12_Standalone_Benchmark_Report.md)。I12/Gate I 不通过，不进入 H01。
