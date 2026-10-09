# H08a1 / 同一实际速度重构的 cell 分支 first jet 报告

2026-10-05，**PASS（局部运动学接口）**。[任务书](H08a1_Actual_Advector_Cell_Jet_TODO.md)、[实施前设计](H08a1_Actual_Advector_Cell_Jet_Design.md)、[API](../api/H10j0_Cch_Spatial_Advector.md)。H08a时空材料积分尚未完成，下一唯一任务[H08a2](H08a2_Piecewise_Characteristic_Source_TODO.md)，H11未启动。

## 实际改动与保持项

唯一 `CchSpatialAdvector` 新增只读 `cellVelocityJet(pressure, point)`，返回同一原速度值、实际分支的物理Jacobian及分支范围。解析微分原面trace和积分lift；RZ正确从X=r²/2变换，并给r=0正则单侧极限。没有另拟合梯度、插值权重、全局平滑速度或时间owner。

只改2 production（原h/cpp）及1个既有短契约。1376份源码/输入哈希审计中仅这3份改变；原cpp删除新增方法后与启动副本逐字相同，原测试删除新函数/调用后逐字相同。原cellVelocity、nativeVelocity、atom polynomial、材料F/R、dual/KEP及FirstJet/Core均保持。材料仍只有用户指定四硬约束，不要求R=0、quarter配分或pressure→dual精确交换。

## 验收结果

XY/RZ × uniform/2:1 × n4/8/16 × full/半域/孔洞/不连通 × affine/非线性/任意native均值，**144组合**。共69,804个分支采样、696轴点、1,008次内部knot两侧检查、23,898次共有native面检查。返回值沿用原值函数，分支范围由独立native-side枚举核对；全部四个物理导数由实际值的分片二次精确差商独立验证。

| 检查 | 最大误差 | 事前门槛 |
| --- | --- | --- |
| 同分支Jacobian / 独立差商 | 4.84531e−14 | 2e−11 |
| 可识别affine解析导数 | 5.13478e−15 | 1e−12 |
| 物理Jacobian散度 / 原D_I | 3.10475e−16 | 2e−11 |
| 共有面及内部knot法向连续性 | 4.71845e−16 | 1e−12 |

20项新拒绝和2组真实regrid重建通过，包含错pressure、非有限/越cell点、FP模式、geometry/full clock及mesh身份；副本/原对象释放后新查询仍相同。已有来源、材料高密度比、tiny/counterflow、三层几何及KE短契约保持。

标准Release通过。首轮8项CTest中7项通过；新增momentum测试最后的regrid夹具直接细化到level2造成非2:1，被原网格检查拒绝。补齐level1过渡层后，仅重编并重跑受影响momentum：**1/1通过（24.56s）**。生产代码和阈值未因此修改，原失败日志保留；最终8项所需检查均通过，不将首轮标为8/8。文档收口后架构/目录另跑**3/3通过（1.30s）**。

原5个数值契约共**4,432行确定性输出**与H10o完全一致，仅排除原ADVECTOR_LAYOUT的两个计时字段；新增输出单列。空间重构、材料/KE数值路径未改变，没有重跑H10o完整28窗口，也没有把历史v2结论扩展为新的时空验证。

## 真实跳跃与后续边界

上述字段中，cell共有面的最大切向差0.2527498，cell内部AMR knot的最大切向差0.0979139；它们是实际分片B的跳跃，不是导数错误。只共享法向不足以符合原Core tracer要求的连续速度迹。

新局部jet不包含跨跳跃事件的时间导数/传递矩阵，不能直接包装成全局连续callback。H08a2必须处理同一实际分片B的轨迹和导数，并先通过真实横穿反例。之后才可接步内FirstJet几何、时空相积分及primal/dual同一积分；本轮没有生成时空F，也没有提供R的时空收敛证据。Pure FirstJet保持独立。

检查目录：`runs/checks/H08a1/20261005_095158_314`；冻结交付：`runs/baselines/H08a1_20261005`。`verification.json`记录源码范围、旧输出逐行核对、各项误差和实际二进制哈希。momentum契约exe SHA256：`A69E83E7CE5DE52555DA59D9F672A4DAB50C8B1102D91BC44535A63DE2D4E26B`。实施前任务/设计副本与首轮失败、最终重跑日志同时保存。
