# I12a11r18 / 同源几何上的有界库存区域积分

目标：只解决一个数学事实：几何体积与 M 不再机器精度相等时，在同一 Hermite 相区上构造有界的输运相含量，使整格积分恢复原库存至既有 roundoff 预算。属于[用户已授权的误差分离路线](I12_Geometry_Inventory_Accuracy_Design_Change.md)，具体公式先见[本项补充说明](I12_Bounded_Inventory_Transport_Design.md)。

允许修改：`HermiteGalsPhaseRegion.h/.cpp` 暴露同次区域积分的叶贡献；新增 `HermiteGalsInventoryRegion.h/.cpp`（共 4 production）；新增 `tests/interface_inventory_region_contract.cpp` 及 `tests/CMakeLists.txt`（2 test-related）。另含本文、API/Report、设计说明、README 与证据目录。

禁止修改：Core 求积/根表示、原 I07/I10/I11、现有时空通量与库存装配、完整驱动、其它物理模块、冻结 PDF/历史 archive 和原精度门槛。

数学约束：唯一几何不变，原 M 不变，无新增独立推进状态。以原 native 容量 C、同源测量 Cg/Vg 定义 a=M/C、v=Vg/Cg；若 a<=v，液/气几何区权重为 (a/v,0)，v=a=0 时约定 (1,0)；若 a>v，权重为 (1,(a-v)/(1-v))。权重必须严格在 [0,1]；几何残差受独立具名预算限制；原/几何容量差与整格表示缺陷仍不超过原 64 epsilon 预算，不将其算作可调几何容限。全单元测量和任意子区域均由原 prepared geometry owner 提供，不再重构第二份界面。

实现要求：新只读 prepared owner 持有同一来源和原区域算子、各叶权重与报告；查询复用同一次几何积分的叶贡献，返回原几何液量、有界库存液量/气量及修正诊断。所有容量、布局、mesh、来源与有限性检查显式；copy/move/lifetime 明确。原区域标量输出顺序和数值不变；不重写求积或用 clipping 修正 M。纯几何恒等权重时返回原几何积分。

准备全域权重时只校验一次完整 mesh identity，通过私有 friend 批量入口复用原区域求积内核；不公开 unchecked 查询，避免每个旧叶都再次扫描全网格的平方开销。普通查询仍逐次校验 mesh。

验收：T0 Release 无 warning；原区域、时空扫掠/库存/reflux、I07/I10/I11 及新相邻契约通过，旧输出不变。T1 覆盖 XY/RZ、uniform/AMR、正反相、正负非机器精度缺陷、纯相微量库存、任意子区域严格界限、整格恢复、分割可加性、独立分段常量参考、来源不变、生命周期及非法输入拒绝。几何匹配预算默认沿用 1e-10；代数压力夹具的大缺陷仅用明确配置，不当作生产门槛。

失败规则：预算、代数或实际积分断言失败即 FAIL 并保留证据；仅缺失具体接口时 BLOCKED。该结果不能冒充共享面通量、完整接受步或 I12 通过；后续接线另立任务。

交付：源码/API、数学说明、算子矩阵与回归日志、误差/耗时、冻结与 archive 哈希、报告和可复核 evidence 包。
