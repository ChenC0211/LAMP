# H10 / 共同位置中心重构与动能相容装配

H10a基础装配设计，2026-10-03；下文保留原推导，现行最终载荷另含[H10c联合矩修正](H10c_Native_Convection_Moments_Design.md)。前置[H09](H09b_Dual_Material_Consistency_Report.md)已通过。冻结CCH PDF §7.2要求共同物理位置重构及skew-compatible结构；PDF保持不变，材料继续采用[完全非守恒路线](H07_Pure_Gals_Hydro_Coupling_Design.md)。算子恒等式PASS不等于完整无黏benchmark PASS。

## 为什么不能只换一个平均值

记内部共同段关联列为b（negative owner +1、positive owner −1），唯一材料率为F。规则XY网格的中心trace为 `t0=(e_i+e_j)/2`，raw对流 `b F t0^T u` 满足

`u^T N = 0.5*(u_i^2-u_j^2)*F`。

2:1处两个native velocity采样点的算术中点通常不在实际共同段测度中心。用单份共同位置重构 `t^T u` 直接代替t0后，raw功变成 `F*(u_i-u_j)*(t^T u)`，一般不满足上述恒等式。分别在两侧重构再平均既不解决共同位置问题，也不能自动保持native动能。两份正确性要求必须独立检查。

本设计显式补齐AMR装配细节：共同位置中心通量仍是 `F*(t^T u)`，其raw incidence与负转置配对装配。最终native动量载荷包含下面明确写出的相容项，不将raw segment通量冒称最终载荷的唯一散射。相对于仅散射重构通量的实现，这是有意的装配差异；其数学性质和精度须分别验收，不能只靠能量恒等式证明PDE准确。

## 单份共同位置affine重构

每条内部段只构造一份t，两个owner共用。物理位置是实际边段面积测度中心：XY为中点，RZ水平段按 `2*pi*r dr` 的径向一阶矩取中心。使用两owner实际相邻cell的一环native同分量速度支持；无新全局未知量、ghost primary或另一份速度状态。

以t0为参考，解几何约束下的最小欧氏改动：`A t=(1,0,0)^T`，A行是1及以共同位置为原点、局部宽度归一化的x/y坐标。使用固定阈值1e−12的SVD求可识别约束；目标不在支持affine hull内明确拒绝。可识别的一维支持不是低阶fallback，仍核验全部affine约束。特别是单列隔离RZ slab的轴向native点共线，而真实径向测度中心偏离该列：缺少径向信息时明确拒绝，不跨不连通域借点或猜导数；短契约将该不足作为负例，正例的每个slab有两列。常数与两个坐标的再现误差按gross系数标度检查256eps。规则XY几何在t0已经满足约束时精确保留t0，不引入不必要的滤波。

此处只定义affine值functional，独立于H06的受约束P2 strain/hoopfunctional；不改写H06、不重复组装strain或压力梯度。实际外边界的半dual体即使位于规则XY网格，其部分内部共同段也不满足native采样中点对齐；同样必须重构。只有采样中点与共同物理位置对齐的段才严格恢复t0，不能把整个uniform mesh错误等同于对齐段。该几何事实由首次短契约暴露，原失败证据保留；共同位置/能量数值阈值不变。

## 相容native载荷

内部段装配

\[
N_g(u)=\frac{F_g}{2}\left[b_g(t_g^Tu)-t_g(b_g^Tu)+\operatorname{diag}(b_g)u\right].
\]

等价于 `0.5*(Nraw-Nraw^T)+0.5*diag(delta_u F_internal)`。每项来自同一F、共同t和原关联b；无另一套rho*u、全局耗散或黏度。共同位置对齐且t=t0的XY段与通常中心通量散射完全相同。每条内部贡献之和为0，且

\[
N_g(U\mathbf1)=U F_g b_g,\qquad
u^TN_g(u)=\tfrac12 F_g(u_i^2-u_j^2).
\]

由此无外力常速度与H09真实材料缺陷保持一致，内部半离散动能没有donor耗散。对称性来自同一trace与其代数转置，不能自行改一侧系数。

## 边界与非守恒材料动能

每个外部共同段使用显式同层transported-component trace `u_b`，装配 `s F u_b` 到唯一owner；端口绑定H09 owner、来源及完整clock。零面积RZ轴段率0且径向trace要求0。无缺省边界值或暗中取owner值；具体物理wall/traction/backflow策略归H15。研究driver可显式给解析或周期trace，但这不是生产周期mesh或完整边界实现。

记全部dual外流为q、质量时间变化为dotm，H09 `R=dotm+q`。无外力动量来源为 `-N+uR`。分方向动能预算是

\[
u^T(-N+uR)-\tfrac12\sum_i u_i^2\dot m_i
=\tfrac12\sum_i u_i^2R_i-P_b,
\quad P_b=\sum_{boundary} sF(u_i u_b-\tfrac12u_i^2).
\]

内部动量总和精确抵消；边界动量、动量功、材料动能项和边界动能功分别报告。零材料缺陷/封闭或相容周期功时得到半离散KEP；一般非守恒材料下不能宣称动能严格恒定。该式不是全离散端点KE恒等式，变步长predictor–corrector归H14。

## 验证和后续

H10a必须验证实际2:1中t0位置误差非零、共同affine重构准确、raw功失配可检出、配对后功恒等式和常速度同时通过；随机/高密度比、XY/RZ、selected边界、来源误用和真实regrid均独立检查。不能以恒等式代替空间精度。

H10b在正式研究目录运行光滑无黏、无界面力的完整时间窗benchmark，检查独立解析误差/空间与时间收敛、周期边界功、动能漂移及显式donor对照。若相容项的AMR一致性达不到门槛，H10保持FAIL并修正本装配；不放宽阈值或用黏度掩盖。一般生产速度Hodge/边界和完整NS时间推进仍归H14/H15。

## H10b审计结果与当前入口（2026-10-03）

[H10b报告](H10b_Smooth_Inviscid_Benchmark_Report.md) **FAIL**：28完整窗口/15,116步，四组末档空间阶1.17–1.34，尽管动能漂移约1e−12。解析外边界下的内部affine载荷反例排除了周期端口/时间步解释：XY真实2:1 n=16应为−1的导数实际−1.15061。共同trace的affine约束没有同时保证负转置项的行矩；H10a功恒等式PASS不能豁免该问题。

此为首次H10b失败记录，历史H10a/H10b报告保持；当前实现和入口以下述H10c复验为准。

## H10c联合矩与当前入口（2026-10-04）

[H10c](H10c_Native_Convection_Moments_Report.md)以最终`N=N_old+J*u`为对象，在统一one-ring共享反对称系数上约束J*1和内部一次/中心二次矩。原内部affine反例已降至7.42e−14，旧material-KE/skew及19/19回归通过；没有AMR特例、库存/bubble或trace精度补丁。

原完整H10b的28窗口/15,116步仍有10空间+4donor失败，production保持H10c；H10i actual source FAIL。[H10i1a Core只读相迹](H10i1a_Native_Quartic_Trace_Port_Report.md)已PASS：正式同源zero/两侧normal Taylor port完成，真实1,408段/40零迹/80侧区间核验通过。真实质量/clip/native出版消费者仍拒绝，材料与最终convection未迁移。下一唯一任务[H10i1b闭合裁剪矩](H10i1b_Closed_Clip_Moments_TODO.md)，随后逐项完成H10i1c–f，再复核H10i并实施H10j空间owner。材料仅四硬约束：唯一共享F、内部面反号/真实AMR子面和、来源一致、KEP/skew；quarter/精确交换仅诊断，真实R/uR与原PDE门槛保持。Pure FirstJet独立，无库存/bubble；H10未关闭，H11/H08a待最终空间通过。材料输入改善不能替代原homogeneous最终PDE空间收口。

用户指定[H08a同源时空材料积分](H08a_SpaceTime_Material_Flux_Upgrade_TODO.md)在空间离散通过之后启动：Pure FirstJet几何独立，无库存/bubble，以唯一同源积分F_rho减小R并继续R/uR与KEP。本轮只登记，未改H08源码。
