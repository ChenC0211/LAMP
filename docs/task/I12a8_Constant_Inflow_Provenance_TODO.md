# I12a8 / 常相入流历史与区间身份

## 一个数学事实

域外入流相态必须是明确的 BC 数据源，其真实范围和 owner 必须参与 interval/face 的同源判断；同一区间不能拼接两份无关边界历史。先支持 axis-aligned owned constant velocity 下，实际整面、时间不变的纯液/纯气入流。这是完整常速度步进的边界依赖，不是一般 BC 历史插值。

## 范围

最多 4 个 production 文件：HermiteGalsTemporalPhaseFlux.h/.cpp 与 Core LevelSetHermiteGeometry.h/.cpp；最多 2 个测试相关文件：新 inflow provenance contract、tests/CMakeLists.txt。另含本 TODO/API/报告、runs/baselines/I12a8_20261001。I12a7 已 PASS 并归档。

## 合同

1. immutable constant inflow source 从一份尚未附加 BC 的 base interval 绑定；保持旧 state、velocity、characteristic 和 mesh owner/真实时间范围。显式输入完整且有序的 actual inflow face phase entries；只允许 Liquid/Gas，拒绝缺面、重复、非入流面和空洞支撑。零速度允许空条目表。不从旧标量或解析 interior 采样猜外部相态。
2. interval 新重载持有该 source，原无 BC 重载保持原语义。sameInterval、appended 和库存的 source 检查包含 BC owner；不存在/存在、独立创建的相同数据均不能冒充同源。半/全区间可共享覆盖其完整范围的同一 BC source，组合不丢失 owner。避免循环所有权：BC 的 base interval 不允许已绑定 BC。
3. 新 translation-face 入口：实际入流边界只消费显式 pure phase history；内部/出流继续 I12a7 同源扫掠，零法向量严格为零。无 BC 时入流明确拒绝，不外推、不读取解析内场补值。范围限 XY 单轴平移/RZ 轴向、完整 active 原域，无 opaque coverage。
4. 将 I12a7 已有 exact total-volume 度量提为共享 strict Core primitive，边界和内部调用同一计算路径。只统一时间/面度量，不改相积分、库存 64 eps 预算或边界液量以求容量通过。
5. 不改 first jet/B/J、bubble、I09/I11、旧平面算子或 benchmark driver。新结果明确标识 prescribed inflow；内部保留原 swept diagnostics。

## 验收

- T0：Release 主程序/相关回归、architecture、81 冻结/参考+3 PDF、前序 archive 与 4+2 预算。
- T1：XY/RZ、正反流、uniform/AMR、纯液/纯气及不同面不同 phase；四个 gross volumes、真实 dt、区间切分、BC owner 生命周期与保留、内部结果和 I12a7 逐字段相同。
- T1：无/异源 BC、不同旧 state/velocity/characteristic、覆盖时间不符、非法面/phase、非 constant/单轴速度、opaque coverage、错误/陈旧 mesh、moved-from source；组合和 I12a3 混源拒绝且状态不变。
- T1：完整常相面集合经 I12a3/I09 装配，M 只消费最终共享通量；验证严格容量、纯相/边界 ledger 与内部抵消，不把 algebra contract 代替完整 first-jet 推进。
- Core/I01–I11、I12a1–a7 受影响合同与旧 crossing probe；不宣称长程 I12 PASS。

## 后续

本项验收后接常速度完整半步/全步 driver，再据实际失败拆解余下曲率/运动不一致；变化速度与非均匀/时变边界另有明确源合同，不由本项推断支持。

## 执行结果

PASS：74 组检查、22 类拒绝、25 项 CTest；24 项旧输出逐字不变。详见 [报告](I12a8_Constant_Inflow_Provenance_Report.md)。I12 仍 FAIL，继续完整常速度驱动。
