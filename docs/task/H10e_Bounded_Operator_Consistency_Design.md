# H10e / 最终native算子的有界矩设计（实现前）

2026-10-04。[任务书](H10e_Bounded_Operator_Consistency_TODO.md)只允许最终对流算子/现有契约/研究driver范围；不改H09真实材料质量、唯一F/R/uR或Interface。H10d的精确全中心P2与P2惩罚失败，不能作为fallback。**本方案数值FAIL，全部候选生产实现已退役**；下文保留实施前设计及试验补充的历史顺序，不表示现行API。见[报告](H10e_Bounded_Operator_Consistency_Report.md)。

## Taylor条件与独立数据

写闭合端口的最终quotient convection为`N=(K+diag(q)/2)u`，K共享pair反对称，q为原实际唯一F的有符号外流。对允许的周期等值采样，压力邻接通过声明的匹配段/平移闭合，原native DOF和half-dual质量均保留。重复采样class只形成同一速度导数，逐行`load=-N+uR`及原材料KE/动量预算按H10d已记录的代数关系保持；不同速度明确拒绝，不静默平均。

最终`a=(-N+u q)/M`的Taylor展开，常数要求`K*1=q/2`；一次矩来自真实有向F和segment测度中心。定义i行真实局部宽度h、最大域跨度L、原F总绝对流量s。三层材料质量/几何时钟原样消费，无另一通量或库存。

以同一全局率尺度f*归一K与s，中心二次矩余项为`r2=sum_j(K_ij/f*)*delta_a*delta_b/h²`，实际导数二次矩为`f* h² r2/M`。要求`|r2|<=(s/f*)*h/L`，可使二次导数项≤`(s*h/M)*h²/L`。同时`sum_j |K_ij|<=2*s_i`使one-ring内三次余项≤`2*(s_i*h_i/M_i)*C_patch³*h_i²`，C_patch来自真实支撑距离/h且须在加密审计中有界。不能仅约束低阶矩而允许系数增长。

此处r2指中心二次多项式`delta_a delta_b / h²`与K的内积，不是对参考解拟合。目标中心二次矩0保留为区间中心，但不要求机器精度零；宽度由Taylor阶条件在数值试验前固定。实际一次PDE再现、二次矩/h、三次矩/h²、行尺度和原完整门槛仍必须独立检查。

## 单一闭合约束最小改动

共享one-ring pair的最终K系数为唯一未知量，最小化到原raw skew/显式中心端口K的欧氏距离。约束为每个active quotient行：

- 3个硬等式：常数与两方向一次矩；
- 3个中心二次矩区间：每项绝对值≤`s_i*h_i/L`；
- 行l1界：系数绝对值总和≤`2*s_i`。

以上常数在试验前固定，所有测度/密度/AMR相同规则。没有惩罚指数、人工黏度、解析解输入或AMR特例。旧opaque接口仍只提供旧内部硬矩能力，不能冒称闭合周期；新API以明确端口数据选择闭合能力，求解失败直接拒绝，没有调用旧解的路径。

用Dykstra最小距离投影统一求解共享pair变量：one-ring内3行affine投影通过小SVD；三个slab分别一维投影；行l1球通过局部排序阈值投影。各块共享同一个K向量且保留投影修正，不能独立拟合后直接拼接破坏skew。所有变量/修正仅prepared暂存，不成为主状态；完成后只保存原raw trace诊断与最终skew correction。

每轮独立复算原硬矩、二次区间、l1界和系数变化，原256eps身份门槛不改。预设迭代预算内不满足则FAIL；不增预算、调界或fallback。s=0行须零有向作用；局部不可识别、支撑不足或全局交集不存在均明确拒绝。可识别性未证明时不宣称二阶或生产效率通过。

数值求解初始化补充：初次逐块Dykstra在uniform n=16的预算内未达到硬矩舍入精度；未改约束/容差/预算。先用同一共享one-ring A与局部3行Gram预条件计算全体C/P1的正交最小距离投影x0，再作原交集投影。因完整可行集合包含于该affine子空间，对任意可行x有`|x-Kraw|²=|x-x0|²+|x0-Kraw|²`，初始化不改变最小改动目标。全局投影的乘子仅prepared暂存，不成为physical unknown；实际硬矩/余项仍按原尺度独立复算。

联合refinement补充：逐块投影用总预算的前半段识别起作用的区间/系数界，后半段用相同one-ring的小SVD块联合解硬矩及这些active等式，检验真实全部不等式与非负KKT乘子。错号active约束释放，遗漏的真实违反约束加入；预算耗尽直接FAIL，不添加迭代额度。中心二次界在舍入误差内按原256eps gross尺度（界值、矩项绝对值及s之和）检验，不能以小区间宽度当成全部舍入尺度；另报告实际界比率及超界gross误差。硬C/P1和material-KE的原门槛保持。

该联合active试验在XY n=8/16 uniform和n=8 2:1通过，但n=16 2:1预算内失败；**该求解实现已归档退役，不能推论有界约束本身不相容**。当前改为对A^T作秩揭示SparseQR，利用Q/R直接得到同一全体硬C/P1正交投影，并在每个Dykstra周期使用此全体投影与原局部区间/l1投影，保留原共享变量、最小距离目标与预算。显式复算剩余约束以检测真实相容性，而非在失败active集合上强行解等式。SparseQR只处理prepared几何约束，无另一物理状态或稠密全局高阶未知量；存储/时间成本必须实测。

独立材料相容性观察器：对于闭合端口及沿非周期方向的admissible affine速度w，若规定advector与该方向正交，PDE要求加速度0。现有scalar R/uR与material-KE必有`sum_i m_i*w_i*a_i = sum_i q_i*w_i²/2 - P_boundary`，q来自同一实际F。因此右侧非零是精确affine PDE再现的必要条件缺口，不能靠stencil求解消失；其h阶决定是否仍允许二阶。观察器只读取真实三层FirstJet/材料、实际F与native样本，不参与系数/载荷。须保留rho=1对照与真实rho contrast、2:1/XY/RZ的加密数据，不能据一例直接改变H09或宣称二阶不可能。

## 来源、端口与验证顺序

匹配周期段/平移只是一般几何关联输入：同分量、反向orientation、端点/测度/唯一F一致、无重复使用；所有非零外部F须有关系，RZ禁止径向平移。原native/source/clock与可复制prepared owner校验保持。normal重复采样等值关系仅来自实际几何，研究driver不能投影载荷或注入解析状态。

新短契约独立组装最终dense矩阵、等值class真实质量作用及明确边界功，检查所有wrapped C/P1、余项界、one-ring和material-KE；非法关系/来源/clock/速度均拒绝。保留原124场景、1,516内部矩与真实R/uR/regrid。原内部affine和全部native逐行n=8/16/32先过，后原完整H10b28窗口/15,116步/84gate。原donor失败门槛不得改变。不先做H11或H08a。
