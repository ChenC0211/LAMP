# H10c / 最终native对流空间矩与共享skew系数

2026-10-05。当前扩展见[H10o报告](../task/H10o_Bounded_Final_Native_Convection_Report.md)：显式闭合入口通过有界最终矩、原反例及v2完整研究；下一H08a，材料仍只有四硬约束。下述H10c描述保留为无端口重载的契约。

H10a原API、来源/clock/axis/mesh拒绝及`-N+uR`保持。最终载荷改为`N_old+J*u`。`points/weights/centerMomentumRates`仍是原始共同trace诊断，不能散射该诊断替代最终载荷。

```cpp
const auto& pairs = convection.couplings();
const auto& moments = convection.momentReport();
```

每个`CchConvectionCoupling`保存有序native索引`first<second`及`coefficient`，表示`J(first,second)=coefficient`、反向为其负值。两索引同分量且属于实际pressure邻接生成的统一one-ring patch；该系数是最终算子的派生修正，不是第二材料通量。copy共享不可变prepared data，不改变现有primary状态。

所有active行约束J*1=0；内部行额外约束最终常数/一次及中心二次空间矩。内部判定来自真实外部momentum segment owner；外部行只保持常数/skew，外部误差仍计入完整benchmark。`interiorRows`是参与矩约束的active内部行数，不包括零率零原K行。

one-ring共享pairs共同解率尺度下的最小改动；小型局部SVD核验rank并白化1/6行约束。联合迭代只用稀疏A/A^T，补偿累积pair解，复算白化前矩残差。临时乘子/迭代向量不成为持久速度或高阶状态。存储及一次乘法为网格线性规模，迭代次数须如实报告；没有dense全局operator或AMR分支。

`CchConvectionMomentReport`给出active内部行数、约束数、迭代数，以及常数/一次/二次矩的最大相对残差。残差分母是原率gross尺度加对应J项的绝对量，生产拒绝阈值为256*double epsilon；独立短契约直接检查最终K+J矩的1e−12门槛。rank不足、非有限及残差不合格明确拒绝，不返回旧算子fallback。

J逐pair严格反对称，所以动能功代数抵消；J*1同时保证常速度和内部总动量。实际material-KE/skew、边界功和uR仍由原生产报告及独立矩阵测试检查。内部P2再现不证明整个有限域或周期端口二阶，最终裁决必须包含全部native边界自由度。

研究入口与原H10b一致：`tools/build.cmd -Preset local-benchmarks -Target LAMP_cch_convection_benchmark`，`tools/run_case.cmd benchmarks/cch_convection`。开发者`audit`保留原n=8/16 affine反例矩阵；`audit32`只检查最细四组瞬时载荷，均在`runs/checks`使用，不能替代完整时间窗。


## H10o 显式闭合端口重载

`build(mesh, material, vector<CchConvectionPeriodicPair>)` 接收实际外段索引与将second映射到first的平移。检查同分量、相反方向、几何/真测度/原F一致、非零外流量全部匹配；不允许RZ径向周期。只读 `periodicPairs()` / `representatives()` 发布声明关系；`applyPaired(mesh, values, sources, clock)` 构造该关系的中心trace，普通apply也验证同一关系和alias值。

这一路径用最终K的常数/点PDE一次矩及二次余项slab、l1球联合求解，目标来自同一actual B与真实M。周期等值类按总M形成共同加速度后分配回全部原half-dual，仍有load=-N+uR，原材料KE/动量身份保持。pairs是K相对于原trace配对的修正，不是新材料流量。

报告的 `boundedRows` / `constraints` / `iterations` / `projectionCalls` 分别记录闭合类行、约束及实际求解成本；`constantError` / `affineError` 为原始矩除原流量gross的误差。`coefficientRatio` 为sum|K|/(2s)，`quadraticBudgetRatio` 为最大二次矩/(s*h/L)，`quadraticBoundError` 为超界量/s；`maximumSupportRatio` 是实际允许对的距离/h。旧 `quadraticError` 只用于无端口H10c入口，不能将其默认0解释成闭合入口二次矩精确。

闭合入口当前验证光滑纯材料和rho1/1e3/1e6的常数/允许affine，不能把任意sharp-cut的精确一次矩当材料硬约束。不可识别/不相容或求解、最终界失败时拒绝，不回退到无端口入口。一般trace重载仍按上面的历史契约，不承诺全native有界性；一般生产边界归H14/H15。准备是联合稀疏求解，n32 AMR仍需数千迭代，没有步间缓存/更新效率保证。
