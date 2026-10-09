# H11 / 局部受控耗散

状态：2026-10-06，**PASS（规定范围）**，见[报告](../task/H11_Controlled_Dissipation_Report.md)与[原任务书](../task/H11_Controlled_Dissipation_TODO.md)。唯一 owner 仍为 `CchKepConvection`；19 项短检查、36 个研究窗口/140 gate 与原 H10b 回归通过。完整流体 App、通用物理边界和全耦合二阶未验收。原 `build`/`apply` 中心路径保留。

```cpp
CchConvectionDissipationConfig config; // strength=1, activation=.2, saturation=.5
auto op = CchKepConvection::buildControlled(mesh, dualMaterial, periodicPairs, config);
auto result = op.applyPaired(mesh, actualVelocity, actualSources, actualMidpointClock);
const auto& load = result.momentumLoads(); // -N_KEP-L+uR
const auto& L = result.dissipationRates();
const auto& samples = result.dissipationSamples();
```

无周期端口的 `buildControlled(mesh, material, config)` 配合原显式 `CchConvectionBoundary` 使用。配置必须有限，`0<=strength<=1`、`0<=activation<saturation`；非法配置直接拒绝。strength=0 不准备传感器，保持中心结果。每次 apply 消费当前给定速度重新检测；prepared 只缓存同一几何/材料的采样投影，不缓存上一次触发状态。

每个真实内部段和计一次的显式周期对使用两端原 pressure one-ring 支撑并集，按真实 dual 质量构造 Favre 加权二次采样像空间。传感器是该采样投影的加权残差范数除以同一支撑的 Favre 标准差；常数变化为零时指标为零。周期坐标来自原声明的平移关系；alias 的质量先求和，同一速度样本不重复加权。采样投影不需要欠定多项式在未采样点具有唯一值。

指标不超过 activation 时 chi=0，达到 saturation 时 chi=strength，中间为平方斜坡。传感器不读取密度跳跃标记、case 名或 AMR 类型；质量只进入加权，F 只读消费同一 H08/H09 材料 owner。原中点近似和 H08a 时空积分都使用自己的同一合法平均材料率。

新增 native 外流为 `L_i += .5*chi*abs(F)*(u_i-u_j)`、`L_j -= ...`。内部动量成对抵消，`D=sum(.5*chi*abs(F)*(u_i-u_j)^2)>=0`。原最终 KEP 矩与 R/uR 不变；周期 alias 按原真实质量分配同一加速度，`dissipationRates()` 返回分配后的真实 native L。相同 alias 间跳跃为零，该连接不需要传感器。

受控模式的 `momentumRates()` 返回总外流 `N_KEP+L`，`momentumLoads()` 已包含 `-L`；消费者不得再次扣除 `dissipationRates()`。后者用于诊断和独立功核验。

结果 samples 保存原 segment/周期 partner、两 native 索引、指标、chi、halfJump、第一端正向新增率和耗散功；它们是被动观察量。每分量 report 增加总 D、最大指标、连接数与活动连接数，继续检查真实预算：`kineticRate=materialKineticPower-boundaryKineticPower-D`。总动能还含实际非零 R 的功，不因 D≥0 就宣称任意材料条件下总动能单调。

未配对的物理外边界仍采用原显式中心 trace，没有暗设储库、ghost、backflow 或额外边界罚项。RZ 固定轴的速度约束和反力由调用方施加，径向速度0使约束功0；原始轴载荷不能隐去。完整 H14 时间接受/H15 边界及无黏两相 Euler/NS 验收不由本接口替代。

prepared/result 保留同一 immutable 材料来源、完整 clock/layout 与生命周期；错来源、外来材料 owner、mesh/regrid、非法轴、非有限或不可表示计算均按原规则拒绝。无材料库存、几何修复、隐式 donor 或最小黏度。重网格必须显式 transfer/rebuild。
