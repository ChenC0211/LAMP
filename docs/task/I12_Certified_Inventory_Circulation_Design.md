# I12 / 原认证区间内的共享库存通量选择

## 原失败证据与约束

r6首次原full诊断显示：第二步half在旧M=0或C的多个cell越界；如rotation uniform cell19，代表值更新为-1.495909855120503e-10，64epsilon预算2.220446049250313e-16，而认证净更新区间[-4.487729565361512e-10,1.495909855120504e-10]包含0。AMR满cell同样存在。所有已输出局部区间都与[0,C]相交，但这只是必要条件，不能逐cell独立选面通量。

原每面四个认证区间继续作为唯一积分误差约束。已有独立代表值只是区间内的一种选择；选择另一个区间内数值仍有相同width误差上界。不改被积分函数、速度、M、几何、time/history、gross总通量或边界通量，不增加64epsilon修复限额。本项只增加联合选择相通量代表值的算子。

## 联合可行问题

固定每个face已有total±代表值T±和物理边界已有liquid±代表值。内部face的两个方向各为一个有界弧：canonical负cell→正cell的L+，反向的L-；范围为原liquid认证区间与[0,T±]的交。两个相邻cell共享同一个变量。

增加一个纯代数reservoir节点R。每个cell有固定弧R→cell，流量为旧独立M；以及cell→R，允许流量[0,C]，表示待求新库存。物理边界两个方向固定为原通量，在R与对应cell之间连接。所有节点流量守恒恰好等价于M_new=M_old+入流-出流及0<=M_new<=C。R不表示新的物理储库或边界条件；其约束是全局库存与原边界ledger的一致性。

初始化内部流量为原代表值，库存弧为合法旧M。根据节点不平衡添加标准super-source/super-sink边，使用确定性最短增广路求可行circulation。库存弧初始化不使用clip后的新M。若所有必要流量不能送达，严格拒绝，并给出不可行状态；不宣称独立cell区间相交就是全局可行。

## 数值、归属与返回规则

所有输入均为有限binary64。残量、节点平衡、min瓶颈及是否饱和用现有Core exact expansion的精确加减/符号；无epsilon吸收、无舍入残量路径或浮点“近似可行”。Core不认识mesh/相/Interface。新TU以strict FP、无LTO独立编译；保留既有FP环境检查。BFS增广具有Edmonds–Karp有限终止界，不能通过忽略小残量结束。

Interface公开新批量入口，仅接受完整有序同一interval的真实BoundedFaceIntegral及同一旧inventory/snapshot。若原代表值已在旧装配的舍入域内，沿用原装配，避免改变已有正常结果。需要联合选择时，将精确可行流量转为原区间内的binary64代表值，交给**原**库存装配和reflux；逐项再检查最终liquid±仍处于原认证区间、原T±和边界值保持。若binary64转换或原舍入修复使证书失效则拒绝，不发布状态。

返回只读结果保留原certified inputs、选择后原inventory result及独立诊断：是否需要选择、增广数、改动面数、相对代表值改变量。误差报告仍使用原四width预算，不将选择差额计作新几何误差或舍入修复。没有第二个M或phase-face owner；消费者只使用该结果内原inventory result的共享通量。

## 验证与限制

Core小网络用独立构造可行流与割条件验证，覆盖双向边、必须沿多边传播、固定边界、不相交／全局不可行、亚正规残量、错误输入/FP。Interface验证provenance、顺序、原区间、旧可行输入不变、真实XY/RZ/AMR积分及original driver。四个smoke和原30/400步full全部实际运行；不能用算子局部PASS代替I12。若遇到原区间确实不可行，需要另项输运方法工作，不在本项扩大积分容限或clip M。
