# I12a11r19 / 有界库存表示接入共享时空面传递

目标：只解决一个数学事实：由 r18 的同源有界库存表示生成唯一的区间面液量，并由原库存装配与 reflux 消费；验证几何—库存差不会在测试中的排空/充满过程中变成负库存或超容量。

允许修改：`include/Modules/Interface/HermiteGalsTemporalPhaseFlux.h`、`src/Modules/Interface/HermiteGalsTemporalPhaseFlux.cpp`（2 production）；新增 `tests/interface_inventory_sweep_contract.cpp`、`tests/CMakeLists.txt`（2 test-related）；本文、Report/API、README、现行设计说明的进度链接及本项证据目录。

禁止修改：Core、r18 权重公式、原库存装配与 roundoff 常数、reflux、I10/I11、完整时间驱动、冻结 PDF/历史 archive、其它物理模块及原 I12 精度门槛。

数学约束：同一旧 M/geometry、时间区间、owned constant velocity、实际 map 与边界历史；当前能力仍为单轴恒定平移，XY 或 RZ 轴向。积分区域/坐标包围/相测量沿用原扫掠链路；用 r18 权重积分库存液量，严格 `0<=F_liquid^+/-<=F_total^+/-`。实际内部面只生产一次，原 M 更新和 coarse=sum(fine) 使用同一缓存面量。不得重置 M、clipping、扩大修复预算或退回旧几何通量。

实现要求：为库存区域 owner 增加显式 overload，与原几何 overload 共用扫掠及边界内核，不复制 stencil。报告区分库存加权/原几何、同次查询的几何液量、库存修正及 L1 诊断；Report 不控制数值分支。prescribed inflow 继续使用同一显式全相历史，不外推 ghost 几何。原纯代数 Transfer 仍不是任意求积的认证对象，实际 producer 必须由调用者显式选择。

验收：T0 无 warning，旧几何扫掠/入流/库存/reflux 与相邻契约输出不变。T1 覆盖正反速度/相位、XY/RZ、uniform/AMR、平面排空/充满、曲线库存差、零速度、区间切分、严格容量、内部抵消、原边界 ledger、coarse=sum(fine)、身份/边界/支持域拒绝及失败隔离。几何差使用原默认 1e-10；制造 4e-11 相对容量的正负差；旧几何通量在相应反例中必须仍明确拒绝，新库存通量必须通过原装配。库存与直接 departure-region 参考差不超过原 64 epsilon 容量预算，区间分割使用原扫掠 1e-10 相对总传递量门槛，reflux 应舍入量级。

失败规则：实际积分/库存/来源断言或旧回归失败即 FAIL，保留原日志；缺失具体接口才 BLOCKED。纯相闭合及完整 first-jet/geometry 发布未在本项切换，不能把面算子 PASS 当成 I12 通过。

交付：源码/API、完整 face→inventory→reflux 场景数据、对照反例、回归/耗时、冻结/hash 审计、报告和可复核 archive。
