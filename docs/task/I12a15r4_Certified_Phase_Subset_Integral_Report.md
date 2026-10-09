# I12a15r4 / 认证相子区域积分报告

**本项 PASS。原 a15 uniform shear 的固定预算失败已解决，57 个主积分全部通过；I12 / Gate I 仍 FAIL。** 见 [TODO](I12a15r4_Certified_Phase_Subset_Integral_TODO.md)、[设计](I12_Certified_Phase_Subset_Integral_Design.md) 和 [API](../api/I12a15r4_Certified_Phase_Subset_Integral_API.md)。本项没有提高原误差、4096 panels 或 depth 40 预算，没有改变独立库存 M 与几何精度分离的路线。

## 实现与依据

原 a15 先认证完整 RK2 single/no-hit 分支及闭 donor 覆盖。对于原相未决 tile，新实现以同一 field 的 r2 affine departure 作用于所有已覆盖原 native Q4 叶，再由 r3 返回每个叶的 phi form。选择首个 form 的一次部分，将所有系数差与原剩余项向外合并成共同 band。实际 donor 无论是哪一叶，都被该 band 包含；没有平均/拟合第二相场，叶外延拓只用来构造对所有可能 donor 都成立的保守包络。

上包络严格负给出明确液体子区域，下包络严格正给出明确气体子区域。原 Snapshot 的精确有理运算及 unit-square 半平面裁剪提取到共享 Core detail header：交点、面积、第一切向矩均在原 Expansion/rational 算术中计算，保持旧表达式顺序。新消费者先精确组合 length、duration、径向权重，再 outward 除法及 pi 区间转换；旧 Snapshot 入口仍保留原浮点表示行为。

非恒定仿射函数的等号直线是零测集，可以在 polygon 裁剪中保留边界；恒等零函数则不能据此认证任何正测度相，明确返回空子区域。子区域计算不使用采样事件或近似根。

对正/负 gross 分别使用点态速度的非负最小界。实际线性时间重构被 begin/end 的空间范围 hull 包含，另保留原 8eps/8tiny 的原浮点混合误差包络，避免同一系数的时间依赖重复导致常速下界变成零。原 donor 的 wl/wg extrema 决定全域与子区域密度界，满足：

`L >= wg_min*T + (wl_min-wg_min)*v_min*A_liquid`

`L <= wl_max*T - (wl_max-wg_max)*v_min*A_gas`。

运算全程向外，贡献仅使用认证下界，再与旧积分界相交。将子区域下界降低到不超过 T.lower 只保持原组件界排序，不改变积分或库存值。速度反号使 v_min=0 时仍保留原 gross/refinement 路线，未将 counterflow 抵消为净均值。原四宽和 E、pending、分割选择、完整来源/mesh/库存身份检查保持。

## T1：独立子区域检查

原 r3 的 232 区域、33,498 次实际回溯与 57,470 次闭叶检查仍通过，完整原输出作为新输出前缀逐字不变。追加 **70 个独立子区域场景、6 类新拒绝检查**：

| 场景 | 检查 |
| --- | --- |
| 64 个半平面组合 | XY/RZ、水平/竖直 face、正/负 gross、精确/非零 residual、几何/库存权重、单 band/两个 donor 的共同 band；对称及偏移三角形的独立面积/第一矩参考 |
| 3 个常数 scalar | 负/零/正；恒零保持未决，不能归为纯相 |
| 可变库存范围 | wl=[0.5,0.75]、wg=[0,0.25] 的独立极值参考 |
| 不同 donor 斜率 | 两个 form 的 s/t 系数同时变化，共同包络包含参考且仍缩小积分不确定性 |
| 反向速度 | 正/负 gross 同时保留 |
| 新拒绝 | 空 forms、错误 wl/wg、越界 content、反向 residual、非有限系数、非 nearest FP 环境 |

非零 residual 的参考上下相位均包含于证书；零 residual 的独立半平面四宽和小于 `1e-12*total`。这些检查只验证实现，不替代上面的全域代数包含性依据。

## T2：同一原积分矩阵

57 个主积分、8 个 Core moments、864 次原实际 B/J 回溯、4 个隐藏薄层事件、8 个 weighted 场景、24 guards 全部通过。总计 **495 accepted panels、1809 evaluated panels，最大 112 panels**。原 4096 上限与所有数值参考/容差保持。

| 原场景 | 现 panels | 现 E | 原预算 |
| --- | --- | --- | --- |
| uniform shear | 112 | 5.7786364562010146e-6 | 5.859375e-6 |
| AMR shear | 71 | 2.3982779741824910e-6 | 2.44140625e-6 |
| rotating cut circle | 20 | 1.7028084049550841e-6 | 2e-6 |
| nonlinear cut swirl | 11 | 1.2253586908408132e-7 | 2e-7 |

uniform shear 的精确 `L+=0.001953125`、`T+=0.005859375` 保持；原 r1/r3 在 4096 panels 后 E=`6.6763159339198168e-6` 仍失败。旧统计的 4508 panels 仅覆盖当时成功的 56 个主积分，另有该失败的 4096 panels；现 495 panels 覆盖全部 57 个成功积分，不能把两个成功数不同的统计混作完全相同工作集。最终直接运行原积分 contract 用时 0.5271845 s。

预算/depth 拒绝检查的原常速平面现在可以被精确半平面直接积分。该 guard 夹具改用原 bilinear shear，使其具有非零剩余带，继续证明小 panels/depth 必须拒绝；未改 57 个主积分的任何场景、参考、容差或配置。

## 初始失败与修复

首轮编译因显式 rational 默认构造与聚合空初始化不兼容而失败，给矩成员加显式零初始化后通过。首轮积分矩阵已全部通过，但新独立半平面测试遇到 exact expansion overflow：outward 计算两个相同一次系数之差产生了人工 subnormal band，使本应精确经过角点的线出现极端小的裁剪边。对两个相同 binary64 系数的实数差保留精确零，消除该人工剩余；没有把真正小系数、真实余带或非零差截断。第二轮通过，随后增加可变权重与不同斜率两个检查，最终 70 场景通过。初始编译/数值失败日志均保留。

## T0、范围与归档

修改 5 production：Slab h/cpp、原 Snapshot cpp、新共享 ExactFaceTimeRational header、TemporalPhaseFlux cpp；2 test-related 为已有 affine_phase 与 bounded_spacetime_flux contracts。没有 CMake/配置/benchmark 参数迁移。原 actual velocity/RK2、原 Q4/strict sign、M/库存权重定义、reflux 与完整 driver 未改。

39 runtime 目标重建成功，70.0788984 s，无 compiler warning。Core 与两个测试的 strict/no-LTO 设置核验通过。**25 项 CTest 全通过**，163.6463864 s；22 项未改数值输出逐字相同，原 affine phase 输出完整前缀不变，原 bounded integral 从 FAIL 变为 PASS。

architecture 保持 12 targets、383 TU、5 exact legacy edges 和原依赖规则。提取共享 header 后为 908 files、3118 includes（原 907/3115），变化仅为该 Core header 与两个直接消费者的必要 include。81 frozen、3 份 PDF、r3 archive 保持 SHA 身份。旧 crossing probe 的 exit 1 与完整输出仍保持，不能用本项 PASS 重分类旧失败。

第一次最终审计通过。归档包含 5 production、2 tests、5 文档、39 executables 与原 DLL；所有 payload 单独核验 SHA256。证据：[目录](../../runs/baselines/I12a15r4_20261002/)、[验收](../../runs/baselines/I12a15r4_20261002/verification.json)、[算子与 shear](../../runs/baselines/I12a15r4_20261002/operator_verification.json)、[回归比较](../../runs/baselines/I12a15r4_20261002/regression_comparison.json)、[范围审计](../../runs/baselines/I12a15r4_20261002/source_audit.json)。

## 尚未完成

本项解决一般 owned polynomial velocity 的完整单分块支路下的相积分预算，尚未接入通用运动的原半步/整步库存与几何驱动。后续要处理所有 native faces 的真实边界/分支覆盖，保持独立 M 的严格守恒、有界性及原精度预算，完成原 rotation/vortex、时间阶、动态 AMR 和曲率验收。I12 / Gate I 未完成，不进入 Hydro。
