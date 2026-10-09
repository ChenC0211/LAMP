# I12p1 / 纯 Hermite–GALS 路线切换结果

**PASS（本任务）。当前纯路线的完整 I12 / Gate I 尚未通过。** 新接口已实现并完成原 12 个长程算例、1,036 个接受步，不再依赖独立库存或守恒相通量认证。

## 已实现

按用户授权的[设计变更](I12_Pure_Hermite_Gals_Design_Change.md)，`FirstJet` 状态要求 `phaseVolumes` 为空；原 `ConservativeVolume` 接口保留兼容默认值。新 `publishHermiteGalsFirstJetGeometry` 复用原 publication 内核及几何有效性检查，以同 revision 的 I06 mixed closure 和实际回溯中心值直接生成唯一 Q4 snapshot，不调用 bubble 体积求根或通量装配。

位置、cut、法向、曲率和体积测量共享该 snapshot。库存残差标为 NaN/NA；每个接受步记录 V0、Vn、相对初值漂移，测量结果不写回状态。原 B_h/J_Bh、顶点/边中点输运、AMR master/hanging 算法未修改。[API](../api/I12p1_Pure_Hermite_Gals_API.md)列出模式、所有权与调用语义。

六份现行规范已注明这次用户授权变更。旧守恒专项停止作为当前推进任务，冻结材料和显式历史比较保持；尚未增加 production App 字典选项或接入 Hydro。

## 验证结果

- 8 组 pure publication：XY/RZ × uniform/AMR × Q2/非多项式；88 项拒绝检查通过。覆盖无 M、中心原值、同 scalar、旧 bubble 拒绝、错误模式/identity/布局/NaN 输入和失败事务。
- 39 targets Release 构建通过（113.66 s），无 warning；28 CTest 全通过（186.98 s）。26 项旧完整输出逐字相同，包括 architecture 与原守恒 smoke；原 legacy crossing 失败输出保持。
- 另执行直接受影响的旧 state、phase-volume、publication 契约，3 项全部通过。
- 修改 5 production + 2 test 文件。81 冻结文件、3 PDF、r7 archive 哈希保持。架构仍为 12 targets / 384 translation units / 910 files / 3121 内部 includes / 5 legacy edges。

原完整运行（226.65 s）：4 个 static 各 4 步、4 个 XY/RZ uniform/AMR translation 各 40 步、2 个 rotation 各 30 步、2 个 vortex 各 400 步，全部完成。无改动原网格 n=8、速度、dt 或目标时间。静止/平移 Q2 的位置、梯度、曲率及逐 cell 几何体积仍通过原精度阈值。

| 动态算例 | 位置误差 | 相对梯度误差 | 绝对曲率误差 | 最大相对体积漂移 |
| --- | ---: | ---: | ---: | ---: |
| rotation uniform | 6.17624e-6 | 4.84536e-5 | 3.68875e-4 | 1.22700e-6 |
| rotation AMR | 6.19707e-6 | 4.83822e-5 | 3.68935e-4 | 1.22283e-6 |
| vortex uniform | 5.51434e-6 | 2.62554e-4 | 9.08527e-3 | 2.54169e-5 |
| vortex AMR | 2.29768e-6 | 1.62926e-4 | 1.68356e-2 | 9.39548e-6 |

体积漂移不是舍入守恒误差；新路线明确不守恒。AMR vortex 在该网格的最大曲率误差大于 uniform，不能用较小体积漂移推断曲率更准。下一任务须检查空间/时间/曲率收敛、跨面与 regrid，随后才能汇总 Gate I。长程跑通本身不等于精度收敛或低黏耦合通过。

审计脚本最初把新增标准库 `<limits>` 误计入 architecture 的内部 include 数，故审计预期曾失败；检查确认 architecture 只计仓库内部依赖，实际旧输出完全相同，修正审计预期后通过。没有修改 architecture 检查器、容限或数值测试。各阶段运行耗时存在主机波动，不据此声称速度提升。

证据：[验证](../../runs/baselines/I12p1_20261002/verification.json)、[原驱动逐例汇总](../../runs/baselines/I12p1_20261002/driver_verification.json)、[源码审计](../../runs/baselines/I12p1_20261002/source_audit.json)、[旧输出比较](../../runs/baselines/I12p1_20261002/regression_comparison.json)、[完整原始日志](../../runs/baselines/I12p1_20261002/first_jet_full.log)。
