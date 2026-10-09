# I10 / Hermite–GALS interior bubble closure API

状态：I10 PASS，2026-10-01。对应冻结方法式 (35)–(37) 与 I-R2；I11 后续已实现同源发布，Gate I 仍未通过。

## 输入与系数基准

```cpp
#include "Modules/Interface/HermiteGalsBubble.h"

const auto mixed = Quantum::Interface::closeHermiteGalsMixedAmr(mesh, candidate);
const Quantum::Interface::HermiteGalsBubblePrediction prediction{
    candidate.identity(), predictedCoefficients
};
const auto bubble = Quantum::Interface::closeHermiteGalsBubble(
    mesh, candidate, mixed, prediction);
bubble.validateCandidate(mesh, candidate);
const auto& closure = bubble.closure();
```

`candidate` 必须尚未绑定 geometry，`phaseVolumes` 与 `geometryCells` 完全同序同域，每个 M 严格在 `[0,capacity]`。M 是求根目标，既不由 predictor 覆盖，也不因 pure-cell residual 而 clipping。I06 mixed closure 与 prediction 必须属于同一 primary identity；mesh 和 mixed 时间层也被校验。

`predictedCoefficients` 按 `geometryCells` 顺序，是当前候选上的临时预测，不是新增 transported state。预测的取得属于调用方的输运/初始化流程；I10 不另推进中心值。它由函数读取，结果不保留外部可写数组别名。

本入口明确选择 **zero-center cardinal extension**：Core Q4 的前四个一维基函数为 `q_i = H_i-H_i(1/2)*q4`，`q4(t)=16*t^2*(1-t)^2`，故 `q_i(1/2)=0`。令 center DOF=0 所得多项式为 `phi_boundary`，则

```text
phi = phi_boundary + c*q4(xi)*q4(eta)
c = phi(cell parametric center)   // 仅在此明确基准下
```

任意其他 boundary extension 的 center sample 不得无说明地当作 c。I11 使用该结果时必须保持这一基准。master 数据直接来自 candidate；mixed 来自同 revision 的 I06；hanging/fine trace 由同一 Core 约束引擎生成，不用其 rounded query 值重建另一套多项式。

## 体积与求根

`psi>=0` 保证增加 c 时负相集合单调缩小，XY/RZ 正测度下 `M(c)` 非增；正规零集的灵敏度是 `-integral psi/|grad(phi)| dA<=0`。理论上，固定有限 skeleton 且在开单元 `psi>0` 时，`c→+∞` 对应 M→0、`c→-∞` 对应 M→capacity；求解器只允许证书给出的有限安全区间，不用极大 c 强行逼近纯相。

`InterfaceCellMeasure::discoverBubbleVolumes` 在一个不可变 native Q4 snapshot 上进行完整分量发现及 25-DOF sensitivity。I12a11r8 将 cut 体积改为与 I12 PhaseRegion 共用的局部 Green、小相直接积分及补量路径；局部 retained cuts 来自同一 scalar，其原 256 次 root refinement 限制保持。导数仍使用既有 mapped GL8 center sensitivity 和 half-face partitions。全局 fixed/free moments 仍读该 partitions，保留原共同原点账本；局部体积与其浮点和可不同。共享测量与回归已验证，运动曲率门槛仍未通过，见 [报告](../task/I12a11r8_Shared_Local_Volume_Measure_Report.md)。

`Volumes::momentErrorEstimateByCell` 记录实际采用的逐 cell 积分估计，Q4 cut cell 现为共享局部积分的估计；全局为逐 cell 和。估计是 GL8/细分差，不是严格误差上界。导数用于建议 Newton 步；最终接受必须重新检查实际局部体积 residual。

每个需修正的 cut cell 从 predictor snapshot 调用 `nativeQuarticCenterInterval`：Core 的 exact Bernstein 数据及向外取整的子盒界覆盖整个系数区间，在可能出现零集的位置证明同一 Cartesian partial 保持符号和正下界，并限制位移。该证明排除修正途中临界零及新增内部闭合圈；它是充分条件，不是最大可达区间。

I12r3 增加可选第四参数 `Snapshot::QuarticCenterDirection::{Both,Decrease,Increase}`。默认 Both 保留原对称搜索及 9 次倍增上限；Decrease/Increase 将原预测 c 固定为右/左端点，沿单侧最多 16 次倍增。每次试探使用同一个 signed-partial/位移证明内核；增加的是有限搜索工作量，不增加允许位移、求积或守恒误差。`retainedDerivativeFraction` 保持原来的导数下界搜索起始比例及后续减半尝试语义，不是承诺最终 partial 下界恰为原梯度的该比例。

先测对称安全区间两端体积并检查与 predictor 的非增关系。若目标已在原区间内，原求解路径保持。若目标不可达，利用 `M(c)` 非增性选择 Decrease（需增加液量）或 Increase（需减少液量），从同一原 predictor 重新证明整个单侧区间。证书端点确实扩展所需方向后，再测两端体积并执行全部既有 moment/topology/monotonicity checks；仍不能括住目标即拒绝。

随后保留有符号体积括区。Newton 必须位于括区内部且步长不超过前一步的一半，否则二分；已满足 volume tolerance 的端点可以直接使用。每次评估同时检查导数非正、括区体积顺序、完整分量数/phase certificates 未变。迭代或浮点存储分辨率耗尽即拒绝，没有旧方法后备路径。单侧证书始终相对原预测 c，不通过多次重定位累积位移预算。

## 配置与结果

| Config | 默认值 / 语义 |
| --- | --- |
| `relativeVolumeTolerance` | `1e-10`；请求的每 cell 最大 `abs(M_geom-M_target)/capacity`，实际工作精度还受库存预算约束 |
| `relativeMomentEstimateTolerance` | `1e-11`；每 cell adaptive moment estimate/capacity，所有试算均检查 |
| `relativeLedgerTolerance` | `1e-10`；逐 cell 容量/纯相原始积分一致性，以及全域 boundary ledger |
| `maxDisplacementCells` | `0.125`；位移上界除以该 leaf 的 `min(hx,hy)` |
| `retainedDerivativeFraction` | `0.25`；既有 Core 区间证书的具名控制 |
| `maxIterations` | 64；批次 root 更新上限 |
| `boundaryContacts/symmetryPlanes` | 显式传给同源 topology discovery；不产生新接触角条件 |

I12a11r1 诊断发现：容限内的旧 M/几何差额仍可能是后续库存舍入预算的数百倍，在几何液量流空时形成负库存。I12a11r2 将原 I08 的 64 eps 预算保持原值移至共享 [私有定义](../../include/Modules/Interface/Detail/HermiteGalsVolumeRoundoff.h)，I10 实际相对工作容限为 `min(relativeVolumeTolerance, inventoryRoundoffBudget/4)`，默认 16 eps（约 3.55e-15）。已经启动的 cut root，其括区、求根终止和最终核验使用同一工作容限。只提高同一方程的求解精度，不修改库存/通量、位移/拓扑证书或 moment gate；仍可能因精度不可达或存储分辨率不足而明确拒绝。r2 原先将该容限同时施加给纯相，导致预算内 17 eps 尾差拒绝，见 [r2 报告](../task/I12a11r2_Closure_Working_Precision_Report.md) 与 [r3 诊断](../task/I12a11r3_Pure_Phase_Gap_Audit_Report.md)。I12a11r4 分离纯相预算为 `min(relativeVolumeTolerance, inventoryRoundoffBudget)`，默认 64 eps；cut 求根精度保持不变。r4 仍因预算内 cut cell 被过度修正而未完成原推进，见 [r4 报告](../task/I12a11r4_Pure_Phase_Roundoff_Budget_Report.md) 与 [r5 诊断](../task/I12a11r5_Working_Precision_Measure_Audit_Report.md)。I12a11r6 将相同的 `min(relativeVolumeTolerance, inventoryRoundoffBudget)` 验收上限用于所有单元；只有超预算 cut cell 才启动 root，预算内保留 predictor 系数和真实残差。启动后的 root 仍求解到上述工作精度，不放宽其证书或改变 M。验收见 [r6 TODO](../task/I12a11r6_Cut_Closure_Acceptance_Budget_TODO.md)。

纯相 geometry 由完整证书确认；纯相预算内 M 不匹配不改系数，也不改目标，保留在 residual 中。超过纯相预算则拒绝创建未标记内部相分量；原 1e-12 相对容量的尾量不再被接受，4 eps 和 32 eps 的尾量按实际残差报告。用户的更严格上限仍生效。raw Green 体积与证书体积同时记录，cut-cell 原始体积严格有界，不 clipping。

Result 不可变，提供完整 `HermiteGalsClosure`、identity/time、逐 cell report 和汇总 report。它保留 mixed skeleton 的所有权与 mesh identity，但**不返回试算 geometry**，不改 candidate、不绑定 geometry、不 commit。`validateCandidate` 在后续消费前检查 revision、mesh 及仍未发布；源对象释放后 closure 仍可读。

逐 cell 报告包含预测/最终 c、求根采用的证书区间、partial 方向/符号/下界、位移系数、initial/geometric/raw/target/min/max volume 与最终导数。I12r3 的 `directionalExpansion` 标明是否扩区，`symmetricLowerCoefficient/symmetricUpperCoefficient` 保留原对称边界；未激活的纯相/已闭合 cell 使用预测值作为两端。

汇总包含 root 更新批次、geometry evaluations、Newton/二分/端点步数、最大初始/最终/纯相相对残差、所有试算的最大求积估计、signed/L1 residual、最大系数修正和位移界；`directionallyExpandedCells` 单独计数，`relativeWorkingVolumeTolerance`、`relativeAcceptedVolumeTolerance`、`relativePureVolumeTolerance` 分别记录活动 cut root 工作精度、全部单元验收上限和纯相上限（后两者相同）。Report 在 gate 后生成，不参与数值分支。工作残差不等于对全部几何/求积误差的严格认证界。

错误消息含 primary identity、cell 和阶段；失败不返回部分 closure。内核不读字典、不写文件，无第二份生产 geometry、KKT 或额外状态更新。当前按全域批量重建临时 snapshot 评估相互独立的 cell root，没有宣称完成局部复用/性能优化。

## 验证边界

见 [I10 报告](../task/I10_Interior_Bubble_Closure_Report.md)。圆/球 XY/RZ、uniform/2:1 AMR、相反修正方向、单调扫描、独立解析截面积分、C1 边界与纯相/失败隔离均通过。

现有 Core 对未解析的 closed subcell cycle 要求物理细化；单 cell 内完整小圈也不能获得本入口所需的统一 signed-partial 证书。不能把这种充分证书缺失称作“所有这类零集都不正规”。半面内部精确相切/重根、正长度零迹等既有限制也未在本任务解决。任意步内跨面输运、真实 topology change、动态 regrid、曲率收敛及低黏电喷雾稳定性仍属后续 gate。

[I12r2 审计](../task/I12r2_Bubble_Rejection_Audit_Report.md) 已在原运动失败的 4 个单元中找到原对称区间外、满足体积容限且通过现有 Core regularity/displacement 同伦证书的根；另有 2 个纯液库存不匹配和 2 个区间外探测未解决的情况。探针只存在于显式测试模式，生产 I10 合同和求解路径未变。同伦证书与 `centerInterval` 的 signed-partial 搜索合同不同，不能据此自动互换或接受试算；后续定向区间搜索仍须独立验收，曲率准确性也须另行验证。

该段描述 I12r2 当时状态；后续 [I12r3](../task/I12r3_Directed_Bubble_Interval_Report.md) 以同一个 signed-partial 证明引擎接入上述定向搜索，没有直接换用 `changeCertificate`。小位移/体积闭合仍不提供曲率误差保证，尤其需关注小 cut cell 的大系数修正；完整运动精度和 Gate I 未因此完成。


