# I12a16r1：物理边界相历史

本项延续 [a16](I12a16_Polynomial_Inventory_Driver.md) 的来源缺口，具体范围见 [TODO](../task/I12a16r1_Boundary_Phase_Source_TODO.md)。独立库存与几何预算的既有分离保持。

新增常相四边历史绑定 unbound interval，按 west/south/east/north 存储显式 Gas/Liquid。它是应用给出的物理边界条件，不是旧几何的外推。历史的旧状态、owned polynomial velocity、原 characteristic owner 与时间覆盖必须匹配；不同历史对象即使值相同也不是同一来源。组合 transfer 继续保留同一历史。

初始命中证书使用实际速度计算的严格区间界，验证原 tracer 构造条件、正持续时间下的原半格 CFL、RZ 轴条件及唯一 outward-backtrace side。对被声明的一边，整个闭面上的法向速度必须严格入流；在其余相交边上必须排除同时 backward crossing。所有正持续时间轨迹随即从原 run 的初始边界分支返回 point、elapsed=0、side，无需 midpoint、旧 donor 或 exterior velocity。不能把这个证书应用于随后才命中的轨迹。

时间积分允许 begin=旧时间。原 tracer 在 dt=0 返回 identity；本API明确只认证 t>旧时间的立即命中。单独 t=旧时间切片对空间—时间体积积分为零测度，边界条件与旧相在该切片不一致也不会改变积分。旧 closed no-hit API 未改变，不用这一事实忽略正测度的 pending 区域、空间切向点或多壁事件。

驱动的 rotation 液滴中心距(.5,.5)恒定，支持半径上界为该距离+radius。vortex 的既有场绕(.5,.5)旋转，每条物理轨迹保持半径，因此同一上界小于 .5 时整段解析物理边界均气相。这个声明仅用于显式边界历史；内部旧几何气相依然可能具有非零 gasWeight，必须由原 inventory source 提供。完整驱动仍须逐个实际面验证原characteristic分支。
