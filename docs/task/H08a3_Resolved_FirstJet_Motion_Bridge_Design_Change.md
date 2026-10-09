# H08a3 / 消费桥接前的真实自由度可用性门槛

2026-10-05，实现前固定。[任务书](H08a3_Resolved_FirstJet_Motion_Bridge_TODO.md)。H08a2r1的80冻结payload、17现行文件及前置81payload已复核；上一goal turn为已实现且验证的进展。

现行FirstJet step真正消费主顶点、主edge中点以及每个native cell中心；挂点由master closure生成，不独立输运。其旧入口直接调用连续速度RK2。CCH来源目前仅在唯一内向初始分支上自动回溯，明确拒绝同时事件。原H08a2r1时间矩阵的32初值避开网格对齐，不能据其通过认定所有实际DOF可用。

本任务先在**不改生产算法**的条件下执行以下可用性先决检查，再决定是否可以接桥；这个检查不是新的材料约束。它检验目标消费入口所需数据是否已经存在，不能用只读内部点查询替代完整消费而宣称H08a3通过。

## 固定范围和先决门槛

预留production范围为 `HermiteGalsFirstJetStep.h/.cpp` 与 `CchCharacteristic.h/.cpp`（4文件），若确需独立中立owner再先拆任务；不得更改H08a2r1原运动/事件裁决。检查仅扩展现有 `tests/contracts/hydrodynamics/cch_characteristic_contract.cpp` 的临时显式审计模式（1 test），标准CMake原target构建，无新工具或运行脚本入tools。没有生产变更之前，审计失败就保存证据并恢复测试源码；不留下失败候选或生产开关。

**T0**：复核冻结及启动源码/输入哈希；标准Release构建原characteristic target。审计模式与当前实际二进制、源码绑定。

**T1前置**：XY/RZ、uniform/真实2:1、n4/8，枚举同一原Pure FirstJet publication的所有独立主顶点、主edge中点、cell中心。取同一实际native面均值构造的零场、轴向平移、XY对角平移/RZ径向伸缩加平移、原光滑仿射场，两个实际target时刻0.125/0.25、原maxStep=0.0125。记录每个真实DOF的自动结果/拒绝；外部/轴点单列，不能删行。

对解析轨迹全程在域内的点，必须提供可消费的唯一内部运动/J，且point/J对解析解scaled误差≤1e−4。零场应对所有内部DOF成立。若自动拒绝，逐个显式初始piece调用公开同一tracer，记录合法结果的差异及剩余拒绝；此枚举仅诊断，不以挑第一个成功分支作为生产修复。恒定平移的内部同时穿面必须单独记录，不能借助缩短时间或扰动DOF避开。

另以真实native切向跳跃且法向为零的速度进行同样审计，记录多分支运动/J及其FirstJet pullback差异；不预设任意分支都应相同，不把一般非光滑解强加为精确仿射。该反例用于判断是否能合法自动选择单侧，不能把均值或参考解注入source。phi/grad只从同一原snapshot和当前实际运动取得，用于审计，不发布候选或改变原状态。

**T2**：只有上述先决通过，才固定并实现typed桥接的源身份、完整clock、DOF/边界消费和独立解析值/梯度时间收敛门槛，随后跑原九项、FirstJet step回归及完整消费验证。若先决失败，本任务必须明确FAIL/尚未接桥，T2未执行；不能用已完成的前置审计缩小原目标改判PASS。

## 失败后行动

若失败源于原tracer的初始/同时事件语义，修改它会改变本任务明示禁止的数学事实，应先形成独立的运动语义前置任务，保留旧门槛裁决与反例。H08a3原消费目标完整保留，待前置通过后继续；不是外部BLOCKED。四项材料硬约束、独立FirstJet、真实R/uR与H11顺序保持。
