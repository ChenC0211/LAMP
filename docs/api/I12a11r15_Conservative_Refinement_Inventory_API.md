# I12a11r15 / 精确细化库存 API

入口为 `refineHermiteGalsInventory(oldMesh, sourceState, targetMesh, volumeConfig)`，返回不可伪造的只读 `HermiteGalsRefinementInventory`。它是准备结果，不是新 `HermiteGalsState`、I11 publication 或 accepted regrid。

输入必须为完整 native Q4 支持与全域库存。原状态时间层和 XY/RZ 度量不变，目标网格只能对同物理域纯细化。r14 工厂产生唯一精确目标 Snapshot，Core 的 coverage、corner-sector、coarsening 和精确表示拒绝规则保持。

原 M 必须严格有界，且与同一旧标量的几何液量差不超过现有 `beta = 64*epsilon` 相对容量预算。cut cell 使用原局部 Green 积分，纯相按原证书使用 0/容量；测量沿用 I07 的 ledger/moment 配置。这里没有可调的新库存容限。

未变叶直接复制原库存。对于细化父叶，设子几何液量为 `V_i`，容量为 `C_i`，双精度参考库存为 `v_i`，原父库存为 `M_P`。定义可行区间 `L_i=max(0,V_i-beta*C_i)`、`U_i=min(C_i,V_i+beta*C_i)`。总差 `delta=M_P-sum(v_i)` 为正时，允许增量为 `a_i=U_i-v_i`；为负时允许减量为 `a_i=v_i-L_i`。

只有 `abs(delta)<=sum(a_i)` 且所有量有效时，才赋值 `M_i=v_i+delta*a_i/sum(a_i)`。实际存储后再次检查严格容量界限、原几何偏差预算、逐父与全局守恒残差。不可行或不能表示时明确拒绝，并给出旧来源 identity、stage、parent、child 上下文。没有事后 clipping、跨父叶调量或遗漏差额；任何非法旧库存都不能借此修复。

| 只读入口 | 含义 |
| --- | --- |
| `sourceIdentity/timeLayer/source` | 原不可变来源及其时间层 |
| `geometry()` | r14 产生的同一精确目标标量 |
| `phaseVolumes()` | 子叶独立守恒库存，按 native index 排列 |
| `cells()` | 同源测量、相位证书、容量和原始 Green 液量 |
| `parentCells()` | 每个目标 native cell 唯一旧父叶 |
| `parentBalances()` | 每父容量、原 M、子几何和、实际库存和及残差 |
| `report()` | 总量、最大偏差及细化/不变叶计数 |
| `validateMeshes(old,target)` | 检查旧新 native mesh identity |

copy 共享同一不可变 Data 和 scalar；moved-from 句柄拒绝访问。结果拥有旧来源与全部目标数据，不借用外部 mesh。来源、时间、M 和 geometry 不会在失败时被修改。下游仍须实现同一 Snapshot 的 I11 受检发布，不能通过 `bindGeometry/finish` 把本结果冒充已接受状态。

详见 [TODO](../task/I12a11r15_Conservative_Refinement_Inventory_TODO.md) 与 [报告](../task/I12a11r15_Conservative_Refinement_Inventory_Report.md)。I12 总体仍为 FAIL。
