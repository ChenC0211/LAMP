# I12a14 / 特征回溯的跨叶相与库存含量覆盖报告

状态：**PASS（本项）**。I12 / Gate I 仍 **FAIL**，尚未进入 Hydro。依据 [TODO](I12a14_Characteristic_Phase_Coverage_TODO.md)、[设计](I12_Characteristic_Phase_Coverage_Design.md)与 [API](../api/I12a14_Characteristic_Phase_Coverage_API.md)，完成同一旧 Hermite 标量上 `H(phi^n(B_h))` 及 r18 有界库存含量的完整 closed-box 覆盖。本项没有宣告一般 phase flux 或完整 benchmark 通过。

## 实现与数学边界

消费 I12a13 的实际 RK2 departure box，使用 old snapshot 的 native leaf search 枚举所有闭交集。叶接缝、仅边/点接触与退化 box 均保留；按 x 端点划分闭条带，再核验每条带的 y 闭区间并集，确保完整覆盖。覆盖检查不使用面积误差容限，也不靠中点或采样确定 donor leaf。

每叶严格相直接来自 Core 原有 `nativeQuarticBoxSign`。所有叶严格同号时才给出整体 Liquid/Gas，否则为 Unresolved。后者不证明真有界面，也不证明无事件，不能视为气相或零通量。没有新写多项式符号引擎。

显式 inventory 入口要求同一 old identity 与同一 snapshot 指针，保留原 source owner。按每叶已有 wl/wg 给出液相 [wl,wl]、气相 [wg,wg] 或未确定 [wg,wl]，整体取范围 hull。该选择保持 0≤wg≤wl≤1；纯气相仍可有预算内的独立液体库存，跨叶纯相也可能具有不同含量权重。未重新匹配几何体积、修改 M 或生成第二份界面。

只修改 `HermiteGalsTemporalPhaseFlux.h/.cpp` 两个 production 文件，新增 contract 并修改 `tests/CMakeLists.txt` 两个 test-related 文件。Core RK2、velocity slab、符号算术、库存权重与积分、bubble、发布、reflux、root CMake、原 benchmark 参数及全部 gate 均未改。

## T1

最终 `interface_characteristic_phase_contract` 通过：

| 指标 | 结果 |
| --- | --- |
| 完整区域查询 | 131 |
| 权威旧相与库存含量样本 | 5,895；包含性容差为 0 |
| 实际正时长 RK2 回溯 | 4,788；其余样本为旧时刻恒等映射 |
| closed native leaf 交集 | 543；单查询最多 28 叶 |
| 隐藏液相反例 | 4：两种网格各 1 个时间薄层、1 个空间内部液泡 |
| 库存配置 | 12：8 个纯相尾量配置、4 个平面 cut 权重配置 |
| 拒绝检查 | 11 项，含保守 CFL 拒绝、来源、mesh、moved-from、partial support 与浮点环境 |

XY/RZ 与 uniform/AMR 四组主矩阵通过。XY 包含 zero、constant/reversing rotation、reversing cubic swirl、shear、axial translation；RZ 包含 zero、axial translation 和 reversing strain，另含精确轴线。四种网格/度量均检查跨叶、接缝线和接缝点。返回的每个闭叶集合与独立全 native leaf 扫描逐项一致，包括只有零宽交集的叶。

时间薄层在 t=0.2、0.25、0.3 三个节点均为气相，却在两节点之间有液相。完整区间正确返回 Unresolved，细分至 [0.274,0.276] 后认证液相；两侧纯时间区间认证气相。相同的同号相切例仍保留 Unresolved，并检查旧标量确有精确零点。空间内部液泡在查询 box 的四角及中心均为气相，完整 box 不被误认证；小区域细分后认证液相。

库存测试明确检查纯气相含量大于 0、纯液相含量小于 1；跨叶权重变化没有被抹平。平面 cut、纯相与未确定区域的整体和逐叶范围都包含实际 `rho^n(B_h)`。query、成功/失败路径、复制、移动及 source 原局部变量销毁后的访问均不改变 old identity、唯一 snapshot 或 M。部分 active scalar 支持在原 a13 入口拒绝，不能绕过其完整支持前提。

上述样本是实现检查。完整包含性的依据是已认证 B_h 包围、closed 叶并集、原 Q4 符号证书及既有有界权重；未以采样同号推断无事件。

## 初始失败与 T0

首轮编译成功。首次运行在 AMR 的远角 constant-rotation 查询 [0.78,0.80]×[0.72,0.74] 上，因整段时间的保守速度包围超过原半 cell 充分条件而拒绝。该查询保留为预期拒绝测试；纯气相认证使用仍在同一圆外、但包围可认证的 [0.78,0.80]×[0.49,0.51]。没有修改速度、时间范围、h、CFL 或任何已有 benchmark 参数，也未声称原远角区域已认证。

第二次运行的 partial-support 夹具错误地保留了仅属于 inactive cell 的 edge-midpoint/vertex，原 skeleton 规则将其拒绝。修正测试夹具后，得到有效的部分标量 source，并确认 a13 按原“full old scalar support”前提拒绝它。两次失败日志均保留，production 实现不需要为它们修改规则。修正后的直接测试通过，用时 `0.3750099 s`；随后增加真正 RZ 接缝、精确相切零点和 weighted plane 的分类断言，再执行完整构建与回归。

36 个 runtime executable 目标重建成功，用时 `61.3382086 s`；22 项 CTest 全部通过，20 项旧数值输出逐字相同，包括 a12/a13、原 trace/vertex/edge、phase/inventory/reflux 与 r23 smoke/跨面。架构输出仍为 12 targets、383 production translation units、907 files、3114 includes、5 exact legacy edges。新测试 TU 不计入这个 production 图。

本项构建日志无 compiler warning；a13 已单列的 legacy 告警没有在本项被修改或抑制。新 contract 的 strict/no-LTO 编译设置从生成 vcxproj 核验。81 个冻结文件、3 份方法 PDF、a13 archive 保持不变。旧 crossing probe 仍为 exit 1，输出不变，继续作为已知 FAIL 保留。

源码范围审计为 2 production + 2 test-related 文件。归档保留 36 个 executable 与原 DLL 的身份，逐条验证 archive/live payload hash。证据目录：[I12a14_20261002](../../runs/baselines/I12a14_20261002/)，包括 [验收](../../runs/baselines/I12a14_20261002/verification.json)、[算子](../../runs/baselines/I12a14_20261002/operator_verification.json)、[旧输出](../../runs/baselines/I12a14_20261002/regression_comparison.json)、[源码审计](../../runs/baselines/I12a14_20261002/source_audit.json)、初始失败、构建/CTest 日志与可复跑脚本。

## 下一项与未解决问题

下一项在本归档后单独写 TODO，将 face×time 上的同源速度范围、相/库存含量范围接入有误差上界的时空积分。Unresolved 区域需要细分或明确的误差预算，不能因粗细求积相近便跳过隐藏事件。

当前仅覆盖 a13 已证明的完整单步、无 hit 分支；边界历史、缩步过程及未认证的 CFL/轴邻域 box 仍未解决。本项不执行一般 flux commit，不证明一般库存更新有界性、旋转/涡旋的完整输运、时间阶、动态 AMR 或低黏/无黏 Hydro 稳定性。独立 M 严格守恒、几何精度预算及原 I12 全矩阵仍是后续验收要求。
