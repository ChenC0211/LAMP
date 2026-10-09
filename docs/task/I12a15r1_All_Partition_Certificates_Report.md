# I12a15r1 / 所有候选分割证书报告

状态：**FAIL**。同一积分的候选证书交集已实现，原剪切反例仍未达到预设精度/工作预算。依据 [TODO](I12a15r1_All_Partition_Certificates_TODO.md) 与 [设计](I12_All_Partition_Certificates_Design.md) 归档；[a15 原失败](I12a15_Bounded_SpaceTime_Face_Integral_Report.md)保持历史身份，I12 / Gate I 仍 FAIL。

## 改动与验证

每次原有切向/时间两个完整二分候选均已计算后，将两份父积分范围求交；仍按原规则选择继续细化的树分支，但保留两者的共同证书。没有相加重复积分，也没有增加候选查询或 leaf/work/depth 预算。pending 仍由实际保留树完整追踪，因此未认证的非零通量不能因共同误差界较小而成功。

这是同一真值的严格范围交集，Core 原 `intersected` 内核、代表值、liquid<=total、所有真实 source 与 RK2 规则均保持。只改 `HermiteGalsTemporalPhaseFlux.cpp` 一个 production 文件，以及既有 `interface_bounded_spacetime_flux_contract.cpp` 一个测试文件，无新增源文件或 CMake 修改。

测试对原 8 组非单位坐标、三次速度、XY/RZ 正反面矩，新增切向与时间两个不互相嵌套的区域分割，核验它们的共同证书仍包含独立物理原函数参考，代表值合法、误差界不宽于任一候选。所有主积分场景另核验 `evaluatedPanels<=1+4*(acceptedPanels-1)`，防止把额外试算隐藏在相同 leaf 上限下。

## 反例结论

原 57 个主积分场景仍为 **56 成功、1 失败**。864 次原实际 B/J、8 组 Core moments、24 guards、四组隐藏薄层和 8 组独立库存尾量检查继续完成。

| 指标 | a15 | r1 |
| --- | --- | --- |
| uniform shear E | 6.676659722917087e-6 | 6.6763159339198168e-6 |
| 固定绝对误差预算 | 5.859375e-6 | 同左 |
| uniform shear leaf 上限/最终 leaf | 4096/4096 | 同左 |
| 深度预算 | 40 | 同左 |
| AMR shear E | 2.4408823264994990e-6 | 2.4408032004604447e-6 |
| AMR shear leaf | 3326 | 3326 |
| 成功主场景累计 leaf / evaluated | 4508 / 17864 | 同左 |

共同证书确有微小收紧：uniform 误差界降低约 0.00515%，但仍超预算约 13.94%。不能把微小改善称为该反例通过。主要矩阵的通过/失败和工作量汇总相同，完整浮点误差输出并非逐字相同；比较记录见证据中的 `counterexample_comparison.json`。新 contract 保持普通失败测试，未加 WILL_FAIL、未降低验收目标。

该证据说明，仅保留另一个已算候选的证书不足以解决现有分块精度问题。本项按预先写明的范围停止扩展算法，没有在同一任务偷偷追加前瞻、容差调整或特定 shear 快捷路径。

## T0 与归档

37 个 runtime 目标完整重建成功；随后测试文件增加上述解析参考/评价量断言，仅重建测试目标。所有构建日志无 compiler warning。相关 CTest 为 23 项：原 22 项通过，新增积分测试仍因同一反例失败；21 项旧数值输出和架构输出逐字相同。新 contract 两次直接运行及 CTest 的完整最终输出一致，构建、耗时与退出状态保存在证据中。

源码审计为 1 production + 1 test；81 个冻结文件、3 份方法 PDF、a15 archive 未变，既有旧 crossing probe 的 exit 1 与输出未变。归档包含本次交付、37 个 executable 与原 DLL，并逐条验证 SHA256。archive 完整性 PASS 与本项数值验收 FAIL 分别记录。

证据：[I12a15r1_20261002](../../runs/baselines/I12a15r1_20261002/)、[验收](../../runs/baselines/I12a15r1_20261002/verification.json)、[反例比较](../../runs/baselines/I12a15r1_20261002/counterexample_comparison.json)、[旧回归比较](../../runs/baselines/I12a15r1_20261002/regression_comparison.json)、[范围审计](../../runs/baselines/I12a15r1_20261002/source_audit.json)。API 签名与 [a15 所有权/失败规则](../api/I12a15_Bounded_SpaceTime_Face_Integral_API.md)不变。

## 下一项

按新 TODO 处理单个分块内的相交区域包围：保留空间—时间相关性，给旧 Q4 拉回相提供可验证的局部线性/余项范围，再对已确定液/气部分计量。先设计并证明误差与舍入覆盖，不用采样同号或粗细求积差替代证书，也不继续微调方向启发或增加预算。原剪切验收仍须在同一参数下通过；完整 driver、边界历史、库存可行性、rotation/vortex、时间阶、动态 AMR 和低黏/无黏耦合均未完成。