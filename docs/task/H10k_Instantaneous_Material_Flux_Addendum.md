# H10k / 瞬时材料 F owner 的实施前补充

2026-10-04。在原H10k设计/初步迁移后，实际时钟库存确认：原CchMaterialTransport同时负责三层质量差和Midpoint积分，只接受Midpoint/End完整区间，不能在H10h的真实Initial/End层准备瞬时F。不能重贴snapshot clock，也不能让观察器另乘密度形成第二个材料通量算法。

因此在相同CchMaterialTransport.h/.cpp内定义唯一只读CchMaterialFlux，消费实际spatial CchPhaseRate和两密度，保存原phase/advector/layout/D_I/full clock，一次计算canonical材料atoms；所有primal/dual view只求和。它接受各合法瞬时时间层，但没有质量库存、dt、R或接受事务。错来源/instance/fullclock/非有限/不可表示相分量明确拒绝，旧face-constant输入不能进入此工厂。

CchMaterialTransport只是此F owner的区间消费者：新增prepare(mesh,begin,preparedMaterialFlux,end)，保留该同一Flux storage，计算原pressure三层几何质量差和dt*F/R。原prepare(mesh,begin,phaseRate,end,rhos)在spatial输入时只编排该唯一工厂，再调用同一区间装配；不重复计算密度权重。任务书暂留的legacy输入只有原显式历史表达，无失败fallback；H10l删除最后旧调用者后退役。

prepare/validateFor仍核验实际advector stamp、原momentumGeometry instance、layout/D_I/mesh与fullclock；初/中/末层瞬时F分别保留原source，区间只消费合法midpoint F。material瞬时owner与区间wrapper是同一通量的数据流，不是备用数值路线、独立global DOF或第二套几何。

文件预算仍2 production/2既有tests；新typed class属于同一material header/cpp/CMake owner，0新文件/target/tools。补充验证三真实层21段F、prepared Flux copy与区间复用，独立density加权参考。原四材料硬约束及R/uR保持，H10l/final native/H10b顺序不变。本补充先于此typed Flux实现，不回写原实施前设计。
