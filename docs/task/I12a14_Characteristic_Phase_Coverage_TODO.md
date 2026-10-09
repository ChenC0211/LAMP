# I12a14 / 特征回溯的跨叶相与库存含量覆盖

目标：只解决同一 old Hermite geometry 上 `H(phi^n(B_h))` 及已准备的有界库存含量 `rho^n(B_h)` 的完整 closed-box 覆盖事实，消费 I12a13 已认证的实际单步 B_h 包围，为一般 phase-flux 积分提供可靠分类和范围。

允许修改：`include/Modules/Interface/HermiteGalsTemporalPhaseFlux.h` 与 `src/Modules/Interface/HermiteGalsTemporalPhaseFlux.cpp`，共 2 production 文件。新增 `tests/interface_characteristic_phase_contract.cpp` 与 `tests/CMakeLists.txt`，共 2 tests。本文、设计/API/报告、文档入口及本项证据。

禁止修改：Core RK2、速度 slab、Q4 符号算术、原 phase/inventory region 权重与积分、M 更新、bubble、几何发布、reflux、原 benchmark 参数、冻结文件及 gate。不把未认证的相区间当成零，不由样本同号推断无事件。

数学约束：用同一已保留 old snapshot 的 native leaf bounds 与严格 Q4 Bernstein box sign。每个 closed departure box 必须完整覆盖，包含 coarse/fine 接缝、仅边/点接触及退化 box；不能用面积误差容限掩盖空隙。每叶 negative/positive 给出严格液/气相；零返回值表示未确定，不等于 phi=0 或确有界面。整体只在所有叶严格同号时认证相。采用 r18 权重时，每叶液相含量为 wl，气相为 wg，未确定为 [wg,wl]；整体取所有叶范围的 hull，保持 0<=wg<=wl<=1。纯几何分类不意味着独立 M 为 0 或 capacity。

实现要求：先写设计。复用 native leaf search 与 Q4 box sign，不重构标量、不另写 polynomial engine。对返回叶交集用纯坐标比较的矩形并集覆盖检查，禁止依赖浮点面积和采样。提供几何与显式 inventory source 两种入口，私有构造只读结果保留完整 characteristic/enclosure/interval/old geometry 和可选 inventory owner。检查 inventory 的 old state identity 与同一 snapshot 指针，沿用 native mesh identity；移动后访问和来源不符均拒绝。

验收：T0 构建、架构、冻结、来源范围与 runtime 审计；无新增 warning，旧告警必须核验归属。相关 a12/a13、phase/inventory/flux/reflux、r23 smoke/跨面旧数值输出保持。T1 XY/RZ、均匀/AMR、rotation/cubic/shear/axial/zero、跨多叶/粗细/轴/边/点；已认证相与实际回溯的权威旧相一致，库存权重逐叶/整体包含性不加容差。增加薄层在时间采样之间穿越、同号相切、空间内部隐藏相的反例，必须保持未确定，细分纯区域后才认证。测试几何纯相仍有非零 gasWeight 或小于 1 的 liquidWeight；source 不符、过期网格、moved-from、lifetime 与浮点环境检查。

失败规则：漏叶/覆盖空隙被接受、错误纯相认证、库存范围排除实际含量、来源丢失或旧行为改变，均 FAIL 并保留证据，不放宽相符号与守恒预算。未确定是有效区间结果，不声称已完成 quadrature。仅在确需外部信息而无法进展时报告 BLOCKED；当前无此阻塞。

交付：2 production + 2 test-related 文件、T0/T1 记录、设计/API/报告、入口、archive/hash。下一项在本项归档后连接 face×time 的速度/含量范围与有误差上界的时空积分；本项不同时实现通量提交，不标记 I12 或 Gate I 通过。
