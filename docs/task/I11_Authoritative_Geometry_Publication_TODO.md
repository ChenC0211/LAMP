# I11 / authoritative Hermite geometry publication

## 任务 ID / 标题

I11：唯一权威 Hermite 几何发布。I10 archive 与 10 项交付哈希已核验；上一轮为有效进展。

## 目标

将同一 primary revision 的 I10 closure 构造成一份不可变 native Q4 snapshot，使体积、完整零集、法向/曲率、实际半面 cut 及后续 phase flux 消费同一 geometry identity/time/source。数值检查失败不得改变候选或原 accepted handle。

## 允许修改

- 新 `include/Modules/Interface/HermiteGalsPublication.h`、`src/Modules/Interface/HermiteGalsPublication.cpp`。
- `include/Modules/Interface/HermiteGalsState.h`、`src/Modules/Interface/HermiteGalsState.cpp`：仅供 publisher 使用的私有同 revision 副本/成功消耗入口；保留 I01 structural binding 行为，更新角色/系数基准注释。
- 新 `tests/interface_hermite_gals_publication_contract.cpp`、`tests/CMakeLists.txt`。
- TODO/API/报告、README 及 `runs/baselines/I11_20261001/`。预算 4 production + 2 tests。

## 禁止修改

Core/mesh、I02–I10 数值算法、legacy geometry/curvature/transport、冻结 PDF/S08 专用源、Hydro/Electricity、App/reader/case。不得用 fitted/height curvature、另一张 PLIC/参数曲线几何或旧方法 fallback。此任务不开展 I12 全部输运 benchmark，不将现有几何不支持情况静默转换为纯相。

## 数学约束

1. 必须消费不可伪造的 I10 Result，并核对其与候选的 revision/mesh/未发布状态。使用 I10 zero-center cardinal extension，把原 primary master first jet、闭合 mixed/c 恢复为唯一原生 Q4；hanging trace 仍由 Core exact constraints 导出。
2. 在最终将返回的同一 snapshot 上复核几何库存误差 `abs(M_geom-M)/capacity<=1e-10`，M 严格有界，完整 topology 与实际 face partitions 均成功；不能以 I10 的 Report 代替重新验收，也不能通过 Report 决定数值分支。
3. 每个 frame/GL8 node 使用同源 raw gradient/Hessian 所得 normal/meridional curvature。RZ 总曲率为 `k_meridional+n_r/r`，真实 axis endpoint 用偶对称极限 `2*k_meridional`；轴身份取自 retained endpoint/chart，不能按近轴阈值猜测。
4. 零集在原 leaf chart 内采样；保留 frame provenance，禁止先舍入全局点再反推 local 坐标。cut ranges 保留现有 offset/length 语义，不以显示端点重建。
5. 发布结果的 state、体积结果、topology、cut partitions 共享同一个不可变 snapshot。新 candidate 才能有新 identity；成功消耗本候选，不能同 revision 重复发布第二 snapshot。复制/移动结果仅共享 owner。
6. 发布是 Interface module result，不是耦合 accepted-step commit；App 保持 accepted handle 直到全模块通过。本 API 失败保留 candidate 的 identity、primary、未绑定状态，旧 accepted 数据不变。

## 实现要求

私有 candidate 副本保留原 revision，不额外 mint identity；全部数值验证/内存分配完成后以不抛异常的步骤消耗原 candidate。新只读 Publication 类型只可由受检入口构造，暴露 state、I07 几何体积、完整 topology、半面 cut、同源 surface samples/正权 quadrature 和诊断。I01 bindGeometry 继续明确为结构性低层端口，不把它冒充数值发布证明。

体积复用 I07，切面复用 MappedFacePartitions，分量/采样复用 Core NativeGraphChain；这些只是同一 scalar 的派生只读数据。所有 checked frame 的 residual、unit normal、曲率有限性需通过具名 gate。内核不读字典、不写文件，不引入新主时钟。

## 验收

- **T0：** Release 主程序、I01–I11/旧 native state/architecture 构建与运行，无新增 warning；文件预算通过。
- **T1 同源与事务：** 所有 source 指针严格相同，identity/time/measure 对齐；输入复制别名、生命周期、move/copy、stale revision/mesh、重复发布、晚期体积 gate 失败与成功重试。失败不污染 candidate/旧 accepted。
- **T1 几何：** XY/RZ、uniform/2:1 AMR、正反 phase 圆/球与多分量；normal/meridional/physical curvature 对解析值，normal 误差 `<=2e-10`、curvature 误差 `<=2e-8`；RZ axis pole 同源偶极限。半面 ranges/root provenance 直接消费，同 source 恒等式检查。
- **T1 体积：** 对每 cell 原 M 不变；发布 residual `<=1e-10`，解析总量误差 `<=1e-10`（总容量尺度）。纯相允许空 topology，几何报告不产生伪接口。
- **T2 模块链路：** 已闭合非零 bubble 的发布通过；改变候选 primary 必须得到新 geometry identity。I08 使用发布的 midpoint state，输出 Q_l/M 指向同一 snapshot；零速度库存不变。此项不替代 I12 自洽时间推进/收敛。
- **回归：** I07–I10 完整输出及旧 GALS/parametric controls 不变；69 个 S08 专用源与 3 PDF 哈希不变。

## 失败规则

数值/ownership/事务条件失败为 FAIL 并在范围内修复；缺接口、需改冻结方法或超预算为 BLOCKED 并记录证据。禁止放宽容限或以几何发布 PASS 声称面内重根/相切、零迹、未解析 closed cycle、任意步内跨面输运已解决；这些仍须 I12 gate 前处理。

## 交付

源码/测试、API/TODO/报告、原始数据与早期失败、源码/runtime 哈希及逐条验证归档。I11 通过后才进入 I12；总体 goal 保持 active。

## 执行结论

**PASS，2026-10-01。** 4 production + 2 tests；15 项 CTest 全通过，69 个冻结源码/3 PDF 未变。唯一 snapshot、RZ 轴极限、非零 bubble 发布、失败候选保留与重试、published midpoint 零流通量链路均验证。首次运行发现并修复临时 domain 移动导致的 cache identity 失效；早期日志保留。详见 [报告](I11_Authoritative_Geometry_Publication_Report.md)。下一项 I12，Gate I 未完成。
