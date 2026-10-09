# H10h / 最终共享pair与真实PDE矩的联合设计

2026-10-04，**固定区间候选 FAIL，诊断闭合**；下文保留实施前设计，原副本归 H10h baseline。[任务书](H10h_Final_Operator_Joint_Feasibility_TODO.md)限制 2 production+1 test+1 research；实际只修改研究观察与文档。H10g 158 payload 已独立核验，production/原门槛保持，无新 QP/API；[报告](H10h_Final_Operator_Joint_Feasibility_Report.md)给出四组不可行证书及下一设计授权。H11/H08a 等待。

## 对象是最终K，不是预投影target

显式闭合端口定义真实native alias class、周期增量与geometry one-ring。所有native自由度及正质量仍保留；alias类的M、q和矩目标求和，actual apply按原M回散同一类加速度，不平均不一致的输入值。共享变量`x_e=K_ij/f*`，反向系数为−x_e，f*为原材料共享率最大绝对值。raw trace只是诊断，最终修正为K−Kraw。周期关系逐项验证真实段、相反方向、测度、唯一材料F及clock来源；非零开放F必须显式闭合，不隐造ghost。

硬常数`sum_j K_ij=q_i/2`同时保证常速度与原材料KE/skew。一次矩目标直接取实际`b_i=M_i*V_i`，V来自原H07总量率/面测度和在真实closed one-ring上可识别的affine functional；不消费解析速度或benchmark参考。周期位移显式wrapped；横向无跳跃坐标的两条必要条件只用于相应方向。

## 预先固定的二阶预算

令h_i为真实dual宽度、L为实际域跨度，s_i为同一原F的行gross，`S_i=s_i+M_i*(|V_x|+|V_y|)/h_i`。S只是一种算子尺度，不生成另一F。令`chi=|rho_l-rho_g|/max(rho_l,rho_g)`，是物理材料对比度，范围[0,1)，不是case/AMR选择器。固定约束：

- 常数矩为硬等式。
- 每个一次矩`|sum K*delta_a/h_i-b_i,a/h_i| <= chi*S_i*(h_i/L)^2`。
- 每个二次矩`|sum K*delta_a*delta_b/h_i^2| <= S_i*h_i/L`。
- 行l1范数`sum|K_ij| <= 2*S_i`。

chi=0时一次区间退化为精确affine等式，同一个约束系统保留均匀中心可行证书；没有失败后选择另一算法。固定轴径向一次矩额外为硬0，原axis/native行不删除。若f*=0且实际V=0，则唯一零算子；若材料率与实际V的目标矛盾，拒绝而不静默返回零。

从Taylor展开，常数项为0，一次矩区间给出的实际导数系数余项界为`chi*(S_i*h_i/M_i)*(h_i/L)^2`；二次矩给出同样二阶尺度，三次余项由真实l1与one-ring半径控制。**这只是候选预算的条件性一致性依据，不是已有H07相通量误差定理**。须实际证明S_i*h_i/M_i及支撑半径加密有界，并检验实际PDE误差阶；高密度比下允许预算不等于稳定/精度通过。固定系数1、2和原空间门槛在试验前规定，不提高H10f的8，不重标H10g失败，不为通过而调整。

## 先给出可验证全局可行性必要条件

横向一次余项r_i=`sum K*delta_i-bPDE_i`必须满足`sum r=c0`、`sum z*r=c1`。逐行上述区间等价于`|r_i|<=B_i=chi*S_i*h_i^3/L^2`，固定轴B=0。故任何lambda都必须满足

`|lambda0*c0+lambda1*c1| <= sum B_i*|lambda0+lambda1*z_i|`。

二维zonotope的所有法向来自lambda=(1,0)及lambda=(−z_k,1)，逐不同实际z_k检查即可。导出最大比率、实际lambda、缺口和支持函数的原始和；若>1（原gross舍入容差），这是该**固定区间**不可行的dual证书，不是任意二阶算子的普适不可能证明，也不是迭代失败。通过此必要条件仍须真实one-ring联合求解，不先计算一组修正b。

## 唯一最终最小改动求解

锚点为一般中心共享Kcenter=sign*F/2，不对raw trace作补丁。求欧氏最小改动`min |x-xcenter|^2/2`，同时满足上述硬等式和全部区间/l1。硬矩用A^T SparseQR揭示秩，验证被舍弃的依赖右端；P为其null-space正交作用。

不等式dual为`max_{mu>=0} mu*(C*x0-d)-|P*C^T*mu|^2/2`。工作集始终保留非负乘子，加入最违反的实际约束；沿新投影法向、减去现有active分量后，同时满足旧active等式。取新约束首次到等式与旧乘子首次到0中的较小步，到0即释放，不直接跳到含负乘子的active等式解。l1使用有效sign切平面，最终仍核验完整绝对值范数。最多`max(10000,4*约束数)`步，不能恢复退役Dykstra/罚函数或增迭代fallback。

检查实际K的硬矩、全部区间、支撑/l1和原gross容差；还检查mu>=0、active complementarity和最小改动stationarity。退化方向须导出可复算dual不可行证书。零误差均匀中心直接返回同一可行最优解，是凸问题的解，不是网格类型特例。

## 验证与失败处理

先真实全局支持函数预检，再实际最终QP。原opaque兼容调用数值保持；新显式端口独立契约含真实非零R/uR、拒绝、alias散射/KE和矩。原内部affine反例n8/16/32与新最终wrapped矩/actual局部affine先通过，再完整native行/48材料，随后才原H10b28窗口/15,116步/84gate。原4个donor控制不修改、不豁免。

任意预检/局部交集/求解KKT/真实apply/空间gate失败均FAIL并归档候选；不得以全局条件通过代替最终系数与精度。若本预算不可行，不在本轮调整常数；从明确证据拆下一数学任务。尚未实现生产QP。

## 首次预检后的结构证书与输入尺度审计

固定系数1的区间在 XY 2:1、rho=1e3/1e6、n8/n32 四组被全局必要条件拒绝，故不进入该候选生产 QP，不在本轮提高系数。实际 owner 与修正整数重载错误后的独立脚本均确认四组；初次临时预检漏判两组的结果已退役，见报告。最终证书最大差异1.22e−15。

两条条件的证书可以直接写在最终矩阵上：对任意实际shared pair(i,j)，硬常数行乘子`alpha_i=2*lambda0*z_i+lambda1*z_i^2`、一次矩行乘子`beta_i=h_i*(lambda0+lambda1*z_i)`，给出`alpha_i-alpha_j+beta_i*(z_j-z_i)/h_i+beta_j*(z_j-z_i)/h_j=0`。因此整个共享skew图的列被消去，剩余正是`-(lambda0*c0+lambda1*c1)`与一次余项区间的支持函数。这个证书覆盖任意pair图，不依赖小segment或one-ring的求解器；只适用横向无跳跃坐标，周期alias同z后仍成立。对真实内部段/闭合端口先做逐pair数值核验，代数展开证明覆盖全部可能pair。

在非alias、非固定轴的native行，原R/uR给出的瞬时导数矩阵为`A=M^-1*(diag(q)/2-K)`。K反对称无对角，所以A_ii=q_i/(2*M_i)，任何共享修正都不能改变这个对角；硬常数又要求K*1=q/2，故行l1至少为`|q_i|/M_i`。因此`h_i*|q_i|/M_i`是最终导数行尺度的输入强制下界，不是某个QP所得系数。记录该下界、真实相体积/dual体积、R与q、actual V；纯gas/liquid盒中仍有较大q时，下一任务需核验真实相面/quarter转移，不能直接归罪Interface或通过改moment系数掩盖。

这两个审计均不更新矩目标、F、R、状态或生产系数。前者只否定本轮固定区间；后者说明输入对算子尺度的限制，不把它自动推论为普适二阶不可能或所有时间步都不稳定。

## 后续设计授权

用户本轮解除 H08/H09 附加过约束，仅以唯一共享 F_rho、内部面反号/AMR 子面共享、来源一致、KEP/skew 为材料耦合硬约束。quarter 分配和 pressure→dual 散度/缺陷精确交换降为诊断。保持实际 R/uR、精度验收及独立 FirstJet，无库存/bubble。授权登记到下一单任务[H10i](H10i_Spatial_Material_Flux_Transfer_Design_TODO.md)，不修改本固定区间或历史 FAIL；H08a 时空积分仍待空间通过，H11 等待。
