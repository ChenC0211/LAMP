# H10f / 实现前相容性与可行证书设计

2026-10-04，**设计候选FAIL**。[任务书](H10f_Compatible_Bounded_Moments_TODO.md)限定最终native convection与现有研究/契约，不改变Interface、材料质量/F/R、旧App或输入。均匀中心系数的可行性证书通过；下述几何修正及后续率加权PDE修正均未接受，已退役，只保留被动原始PDE审计。QP未实现。见[报告](H10f_Compatible_Bounded_Moments_Report.md)，下文按实施前设计/证据补充的历史顺序保留，不表示现行生产能力。

## 先分清代数交集与实际PDE目标

H10e固定问题为最终闭合quotient skew K、硬C/P1、中心P2区间和行l1界。其迭代失败不是不可行证书。本轮先从同一实际F构造`Kcenter_ij=sign*F/2`，每个内部共享段和声明的中心周期段只装配一次；原half-dual/native DOF保留，几何等值class仅用于系数审计。Kcenter是一般中心KEP系数，不使用解析问题数据、网格类型判断或trace精度拟合。

独立逐行检查H10e原固定条件：`K*1=q/2`，`sum K*delta/h=b/h`，`|sum K*delta_a*delta_b/h²|<=s*h/L`，`sum|K|<=2*s`，并验证每个非零pair属于统一真实闭合one-ring。uniform n=16若满足即提供可复核的交集可行证书；AMR不满足亦如实记录，不能当作不存在其他系数的证明。证书模式只输出actual coefficients/targets及误差，不参与原full或状态。

空间Taylor导数为`a=(-N+u*q)/M`。几何一次目标b由实际F和segment测度中心产生，**`b/M`未必等于实际规定advector**，尤其在分相材料与非均匀采样下。真实FirstJet/H09的48组观察器将追加逐native一次导数审计：XY advector=(0.4,0)、RZ=(0,0.4)，横向affine probe的精确导数为0；记录b/M、按真实质量的全体RMS、h尺度与相容修正后的同一指标。若存在比二阶更大的局部PDE缺陷，不能靠只消除全局Psi接受该设计。

## 相容修正的候选数学依据

对于唯一周期方向的闭合端口，z为横向非周期实际native坐标。K反对称及常数矩要求`sum b+sum z*q=0`、`Psi=sum(z*b+q*z²/2)=0`。实际共享F使第一个条件成立，第二个可能非零。

设`w_i=s_i*h_i³`，权重为同一实际率与几何，零率行不产生新输运；`zbar=sum(w*z)/sum(w)`，`D=sum(w*(z-zbar)²)`。最小化`sum(delta_b²/w)`并保持`sum delta_b=0`、`sum z*delta_b=-Psi`，得到

`delta_b_i=-w_i*(z_i-zbar)*Psi/D`。

修正不会改变唯一F/q/R、常速度或skew。其实际一次导数余项满足

`|delta_b_i|/M_i <= (s_i*h_i/M_i)*(h_i²/L²)*C_aff`，

其中`C_aff=L³*|Psi|/D`。只有C_aff在加密中有界、实际b/M及最终apply的PDE误差二阶时才可接受；不以投影到range本身作为证明。固定候选识别界C_aff<=8，失败不调界。D=0且Psi非零则明确不可识别。相容修正先作为被动数据审计，不直接写入生产。

## 若可识别，最终有界QP的单一求解

最终最小改动可锚定上述一般中心KEP系数，避免对失败的trace矩做补丁；同一F与原raw trace诊断保持。硬矩A*K=t，二次slab及行l1界仍使用H10e固定预算。非均匀材料的一次目标须先通过上述实际PDE审计，不暗中切到旧解。

硬矩秩揭示QR得到唯一欧氏投影x0及null-space正交作用P，不用正常方程。对于不等式C*x<=d，dual为`max_mu>=0 mu*(C*x0-d)-|P*C^T*mu|²/2`，`x=x0-P*C^T*mu`。

从非负乘子工作集开始，每次加入实际最违反的约束；新约束的null-space方向同时保持旧active等式。步长取到新约束等式的步长与旧乘子首次到0步长的较小值；到0者释放。不存在可动方向且乘子方向无可释放项时，须给出可复算的dual不可行证书，不能只报“不收敛”。active Gram的Cholesky可增量维护，避免每次全局正常方程或一组错误active等式。行l1用实际sign支撑的有效线性cut生成，最终仍检查完整l1范数。

固定最多`max(10000,4*约束数)`次工作集步，与H10e预算同阶；最终原gross矩/所有不等式、非负乘子和最小改动KKT必须独立复算，必要的舍入refinement不增加预算。存储/准备时间实测，不预称高效。不通过则退役，不能恢复Dykstra/惩罚或增黏度fallback。

## 验证顺序与现行能力

本轮第一阶段只新增被动证书/材料一次PDE审计。上述修正和QP尚未生产实现；先得到证据，再决定是否可在本任务范围实现。后续仍须旧短契约、真实闭合端口/来源拒绝、内部affine、全部行尺度通过，再原28窗口/84gate完整研究；不进入H11/H08a。

## 第一阶段证据后的PDE目标补充（生产实现前）

中心通量证书已在XY/RZ uniform n=8/16/32满足H10e原约束，特别是XY n=16各硬矩/二次超界均为0、l1比率0.25；证明该处交集可行。真实材料审计则发现几何b/M在横向零PDE导数下产生明显局部误差：例如XY 2:1 rho=1000 n=8 RMS约0.25，消除Psi仍约0.248；RZ uniform也有分相误差。故不采用上述仅修正几何b的候选目标。

最终PDE一次目标应直接取`bPDE_i=M_i*V_i`，V从同一H07实际`totalFaceRates/faceMeasures`恢复native normal advector，再在真实one-ring用局部affine functional取得另一分量。固定轴速度0，原clock/geometry/材料来源保持；没有另一材料通量、解析速度输入或独立高阶状态。Pure uniform bPDE与几何b一致；分相材料不再以错误的flux-centroid目标代表点采样PDE。

对单周期方向的横向坐标，两个必要条件改为`c0=-sum(bPDE+z*q)`、`c1=-sum(z*bPDE+q*z²/2)`。使用`S_i=s_i+(|bPDE_x|+|bPDE_y|)/h_i`，`w_i=S_i*h_i³`，以同一weighted最小改动得到

`delta_b_i=w_i*(alpha+beta*(z_i-zbar))`，

`alpha=c0/sum(w)`，`beta=(c1-zbar*c0)/D`。

一次导数余项界为`(S_i*h_i/M_i)*(h_i²/L²)*C_aff`，`C_aff=L²*(|alpha|+L*|beta|)<=8`在生产数值试验前固定，并须实际加密证明有界及二阶。不得只让global条件为0就接受。观察器先追加真实PDE目标、两个必要条件和逐native修正数据。

S只定义算子系数尺度，不是材料通量。增加真实M*V/h项是点采样PDE目标所需尺度：相在dual内部而未穿过其边界时，真实质量/速度可非零而实际材料F很小，单以s约束K不能代表M*V；原R/uR、F与KE关系仍保持。新最终QP预算固定`sum|K|<=2*S`、`P2<=S*h/L`；这是明确的PDE目标设计变更，原H10e s预算的可行证书不改写，不能把新预算用于重标旧试验。原完整PDE数值门槛不变。
