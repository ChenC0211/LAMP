# I12a11r14 / 精确 Q4 细化标量来源

## 目标及范围

只补齐一个算子事实：对同一个完整 native Q4 scalar 做纯细化，细化后每个叶的多项式严格等于原父叶的限制，支持 cut cell；不通过舍入采样另拟合新界面。现有 `ActiveQuarticHermiteTransfer::prepare` 刻意保护全部 interface cell 并拒绝细化它们，不能直接用于 r13 后的局部重试。

允许 4 production：`include/Core/Geometry/LevelSetHermiteGeometry.h`、`src/Core/Geometry/LevelSetHermiteGeometry.cpp`、`include/Core/Geometry/Detail/LevelSetHermiteQuarticLeafData.h`、`include/Core/Geometry/Detail/NativeHermiteTransferSampling.h`；2 test-related：新算子测试与 `tests/CMakeLists.txt`。另含本 TODO/报告、API 文档、README 和 `runs/baselines/I12a11r14_20261001/`。

禁止修改 Interface primary/M/flux、B/J、I06/I10/I11、冻结代码与所有容差/工作限。本项不提供 accepted regrid，不将新的 scalar 绕过 I11 接入实际推进。

## 数学约束与实现

新增 Snapshot 的显式 native Q4 纯细化工厂，输入完整旧只读 Snapshot 和目标 mesh。仅支持同一物理域、完整 active 覆盖、无 corner-sector；目标每个叶必须包含于恰好一个旧父叶，允许未细化叶，拒绝 coarsening/跨父叶重拟合、错误来源、域变化与不支持的 coverage。

复用现有 exact-chart transfer locator、homogeneous polynomial restriction、shared master/hanging constraint 引擎。新 master value/first/mixed jet、edge value/normal 和 center 直接保留原来源的 exact numerator/denominator；floating jet/DOF 仅为查询值，不作为新独立输入重新装配。源数据在构造中只读，成功结果完整拥有自己的精确数据，不借用可变外部 mesh/patch，不需要旧 Snapshot 继续存活。

构造结束逐目标叶核验 25 个 normalized Hermite 数据与旧父叶的精确限制相等（交叉相乘的 exact identity）；相等性连同 Q4 唯一性证明整多项式不变。现有 constructor 和旧 double transfer 路径保持原语义。不得重复实现另一套 restriction/stencil。

执行补充：在原限制内核外统一同一父叶的清分母因子，避免逐点分母的冗余乘积；原 cardinal datum 保留已有精确表示，交叉比较前只作精确的二进制幂平衡。所有共享约束的 physical datum 也须逐项证明相等。新增大坐标、非二进制坐标、连续两级细化和重复构造检查；这些检查不放宽原容限。继续沿用 exact expansion 的有限指数窗口，禁止丢弃不可表示的尾项。

## 验收

T0：Release 无 warning；原完整 30 项 CTest 保持通过，除新增文件导致的架构计数外原输出不变；81 冻结/参考文件、3 PDF 和 r13 archive 不变。

墙钟执行补充：初次完整 CTest 的 bubble 达到 180 s 时限。同环境原封归档程序完整运行 `245.7850977 s` 且全部断言通过，归档/当前 smoke 分别 `50.03/50.31 s`、输出完全一致。依现有 `interface_conservative_contract` 的 300 s 设置，仅将 bubble 的 CTest 墙钟保护时间改为 300 s，再复验该项；保留首次完整运行与归档对照。不得改变测试内容、求解迭代/工作上限或数值容限；其余已通过项不重复运行。T0 以初次其余 30 项和 bubble 复验的明确来源合并审计，不伪造一次全通过的日志。

T1：XY/RZ，uniform→uniform、uniform→AMR、AMR→更细 AMR、多层和未变化 mesh；全 Q4 与非全局多项式的任意兼容 C1 数据，含原 center bubble。检查 exact restriction certificate、同位置 value/gradient/Hessian、共享 trace、phase cuts/体积与曲率消费，以及源销毁后的 owner 生命周期。覆盖错误 domain、coarsening、非 Q4、partial-active、corner-sector 和 moved-from 输入的拒绝；失败不改变旧来源。

浮点几何/导数比较使用原算子级量级容差，严格多项式相等不能以浮点采样替代。若目标主量格式不能保持 exact restriction，则 FAIL 并报告，不用普通采样伪称精确。本项即使 PASS，守恒 M 细化分配、受检 I11 regrid publication、局部加密重试与 I12 数值 gate 仍须后续独立完成。

## 交付与失败规则

交付源码、执行日志/耗时、逐场景证据、完整回归差异、冻结/runtime 哈希、API 和 archive。实际数值不满足则 FAIL；只因具体缺失依赖无法执行才 BLOCKED，不提高容限。完成本算子后再推进与同一精确几何匹配的守恒库存转移，不进入 Hydro。
