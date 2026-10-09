# I12p2b / 纯 FirstJet 动态 AMR 传递

目标：在真实零时间网格切换后发布可继续实际 first-jet 输运的唯一 Hermite 几何，细化保留已有 exact Q4，粗化显式插值。

允许修改：新增 HermiteGalsMeshTransfer.h/.cpp；HermiteGalsPublication.h/.cpp（4 production）；interface_hermite_gals_benchmark.cpp、tests/CMakeLists.txt（2 tests）；设计/API/报告/证据。

禁止修改：Core exact restriction/几何、actual tracer、I05/I06、独立库存/通量/bubble/reflux、Hydro/Electricity、旧算例/容限、冻结资料。不得把测量体积回填 M，不引入解析重置或失败后的方法切换。

数学约束：source 为合法 FirstJet publication，完整 Q4 支持，同物理域；网格切换不推进时间或步号，新 revision 有新 identity。ExactRefinement 必须调用原 refineNativeQuartic，不能 coarsen；其 closure 的 mixed/center 是原标量限制的派生数据，不是新增独立动力学。FirstJetProjection 显式采样旧 scalar 的目标 master first jet/edge/center，再用原 I06 和纯 publication；不宣称一般分片 Q4 精确或体积守恒。失败保持 source/旧网格可用。唯一 scalar 仍服务 cut/法向/曲率/体积。

实现要求：先写[设计](I12_FirstJet_AMR_Transfer_Design.md)。两个模式共用真实目标 skeleton 的 master/hanging 归属，不复制约束算法。原 publication 提取 private 同源测量/拓扑/表面检查，旧入口行为保持；exact 模式只能由受检传递 producer 调用，不能新增公开任意 snapshot 发布旁路。粗化/投影拒绝丢失连通分量或端点类型，记录实际体积变化；该检查不冒充误差上界。

验收：T0 <=4 production/2 tests，完整39targets、原28CTest/输出保持、freeze81/PDF3/p2a archive保持，无新增warning。T1 XY/RZ、uniform↔AMR、同网格、非多项式、正反方向、来源/域/支持/模式拒绝，exact 的值/gradient/Hessian <=1e-10、trace <=1e-12、同点曲率 <=1e-9；全局Q2投影再现保留原精度。T2 实际推进中 uniform→AMR→uniform→AMR，再继续原 static/translation/rotation/vortex，保持原目标时间与参数；输出每次网格切换及接受时间、体积漂移、误差。粗化非多项式三档 n=8/16/32 梯度/曲率末阶 >=2/1.5。所有实际失败原样记录，不能删组或降低阈值。

失败规则：旧源被修改、时间/identity错误、多owner、方法隐式切换、伪造exact/coarsening守恒、后续实际推进或精度门槛未通过均FAIL；内部数值问题不标外部BLOCKED。

交付：模块 API/实现、算子与动态网格实跑证据、回归/审计/归档、报告和剩余 I12 gate 项。

实跑布置（运行前固定）：动态 full 沿用 p1 原 12 组 n=8、dt/目标步数，在 floor(N/4)、floor(N/2)、floor(3N/4) 完成步后交替 uniform/AMR，时间不变；smoke 各4步，必须发生3次网格切换。粗化收敛覆盖 XY/RZ × uniform细网格/局部AMR来源，先在细网格实际平移2步（dt=.01），投影到 n=8/16/32 粗网格，再实际平移2步，统一终止 .04。无内部解析重置。

结果：**PASS**，见[报告](I12p2b_FirstJet_AMR_Transfer_Report.md)。16 exact、8 projection、45拒绝、12粗化收敛组及12完整动态组通过；动态 full 1,036步/36次网格切换，35CTest及4直接旧契约通过。完整 Gate I 尚未通过，下一项 I12p3。
