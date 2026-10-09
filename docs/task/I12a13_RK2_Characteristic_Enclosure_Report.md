# I12a13 / 同源 RK2 特征时空包围报告

状态：**PASS（本项）**。I12 / Gate I 仍 **FAIL**，不进入 Hydro。依据 [TODO](I12a13_RK2_Characteristic_Enclosure_TODO.md)、[设计](I12_RK2_Characteristic_Enclosure_Design.md)与 [API](../api/I12a13_RK2_Characteristic_Enclosure_API.md)，本项完成 owned polynomial velocity source 下实际单步、无 hit RK2 分支的 B/J 区间包围，没有完成一般相通量。

## 实现与范围

旧 double tracer 与新区间求值共用 midpoint、target、Jacobian 的原代数表达，保持括号及运算顺序。中点速度按实际 arrival/old endpoint 查询平均，保留半步时间下溢为零时的 arrival 分支；不改为直接查询中间时刻的 slab 值，也不替换成连续 ODE 解。新区间计算沿用已有 strict/no-LTO Core TU 及向外舍入算术。

Core 对空间 box × arrival-time interval 给出 midpoint、departure 与 J 范围，同时检查域、原半 cell 速度限制、RZ 精确轴条件及无新 wall contact。能证明整个区域走完整单步时才返回；未知边界、缩步或 opaque active coverage 拒绝。闭时间段包含旧时刻的恒等映射；退化空间线/点与单点时间有明确行为。

Interface 的只读结果保留原 interval、唯一 old geometry、owned velocity 与 characteristic owner。没有虚构时间节点 geometry ID，没有独立第二个 tracer，没有修改 M、几何发布、bubble、相积分或 reflux。原 constant/opaque source 不自动获得 polynomial 能力。mesh 检查继续使用精确几何/拓扑数据身份。

5 个 production 文件为 `LevelSetBoundaryTrace.h`、`PolynomialVelocitySlab.h/.cpp`、`HermiteGalsTemporalPhaseFlux.h/.cpp`；2 个测试相关文件为新增 contract 与 `tests/CMakeLists.txt`。root CMake、原 benchmark 参数、冻结文件及全部库存/几何预算不变。

## T1 结果

新 `interface_characteristic_enclosure_contract` 通过。8 组速度配置包括恒定/反向旋转、恒定/反向 cubic swirl、shear、RZ strain、RZ axial translation、zero；均匀 16 叶与 AMR 28 叶各 10 个主包围区域，共 20 个，包含 4 个精确 RZ 轴线区域。

| 检查 | 结果 |
| --- | --- |
| 实际 Core/Interface 回溯 | 4,044 次，B/J 均在返回范围内，包含性不加容差 |
| 旧时刻恒等映射包含性 | 500 个点 |
| End / 真正 Midpoint view | 均通过，同一 old/velocity owner |
| 独立仿射离散 RK2 对照 | 最大误差 `1.1102230246251565e-16` |
| 主矩阵最大 L1 travel 上界 | `3.4200000000000029e-2`，小于 AMR 的 h/2=`0.0625` |
| 空间/时间子区间、退化 face/point | 包含性及范围收缩通过 |
| 拒绝检查 | 24 项，通过；不等同于 24 个不同错误类别 |

回溯计数为 20 × 202 次主矩阵/单点检查，加 3 次物理 outflow 和 1 次 subnormal 时长回溯。恒等映射 500 点单独计数，不冒充实际正时长回溯。独立参考只用于仿射离散 RK2，不声称它是连续旋转解。

边界测试同时检查实际 inflow/精确 wall-contact hit 与区间拒绝；midpoint/target 域覆盖不足、CFL、非法 RZ 轴、opaque coverage、非法 box/time/domain、rounding/FTZ 均明确拒绝。source、interval 和结果的 moved-from、复制、owner lifetime、不同 source、网格细化等契约通过，成功/失败均不改变旧 geometry/M。密集样本验证实现；覆盖的依据仍是 outward 算术、共用表达和充分分支条件。

## T0、失败记录与审计

35 个 runtime executable 目标重建成功，用时 `248.4725509 s`。21 项 CTest 全部通过，用时 `157.7019925 s`；19 项旧数值测试完整输出逐字相同，包括原 trace/vertex/edge、a12 source、共享通量/inventory/reflux、r23 smoke/跨面等。架构输出同样未变：12 targets、383 production translation units、907 files、3114 includes、5 exact legacy edges；新测试 TU 不计入该 production 图。旧 crossing probe 仍 exit 1，完整输出与 a12 相同，不把已知 FAIL 改成通过。

保留两次初始失败：首次测试编译引用了不存在的 `Snapshot.valid()`，改用实际 `nativeCellBounds()` 查询后编译通过；首次运行错误地要求独立分配的相同 mesh 对象被拒绝，修正为现有 `NativeMeshIdentity` 的精确几何/拓扑身份语义。细化/不同拓扑仍拒绝，旧规则未改。修正后的直接测试用时 `0.1967737 s`，随后完整重建与 CTest 再次通过。

完整重建触发 `ImmersedGeometry_impl.h` 的 11 条 C4267 `size_t → int` 缓存序列化告警。本项前后该文件哈希完全相同，11 个行号逐一核验；没有新增源码告警，未改 Plasma 编译选项，也未抑制告警。初始审计的“所有告警为零”条件因此失败，已按现行 02 的“无新增 compiler warning”修正并单列证据。TODO 原先的“无 warning”简写同步澄清；这不是历史版本重编译对比，也不能声称完整构建零告警。

81 个冻结文件、3 份方法 PDF 与 a12 archive 未变。源码范围精确为 5 production + 2 test-related；严格浮点/no-LTO 编译方式从生成的 vcxproj 核验。证据保留 35 个 executable 与原 DLL 的 runtime 身份，并在归档时逐项检查 archive/live payload hash。

证据目录：[I12a13_20261002](../../runs/baselines/I12a13_20261002/)。主要记录：[验收](../../runs/baselines/I12a13_20261002/verification.json)、[算子矩阵](../../runs/baselines/I12a13_20261002/operator_verification.json)、[旧输出对照](../../runs/baselines/I12a13_20261002/regression_comparison.json)、[源码范围](../../runs/baselines/I12a13_20261002/source_audit.json)、[旧告警身份](../../runs/baselines/I12a13_20261002/warning_audit.json)。目录含初始失败、构建、CTest、审计和打包脚本。

## 未解决与下一项

保守区间可能拒绝实际部分轨迹可行的 box，尤其接触轴邻域或物理边界时；本项没有 boundary-history 或缩步全过程认证。也没有建立一般 old phase 的跨叶符号覆盖、相事件完备、时空 quadrature、一般速度下的共享库存有界性或完整旋转/涡旋/时间收敛。原几何—库存分离策略保持，独立 M 不因新包围而被重置。

下一项在本归档完成后单独写 TODO，把 B_h 的 departure box 接到同一旧 Hermite 标量的 native leaf / 严格相符号覆盖；不通过采样同号推断无薄层或无事件。原 I12 完整矩阵仍是最终 gate。
