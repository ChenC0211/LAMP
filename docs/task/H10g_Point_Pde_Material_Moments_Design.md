# H10g / 点采样材料矩的二阶误差范数设计

2026-10-04，**候选FAIL并退役**。[报告](H10g_Point_Pde_Material_Moments_Report.md)记录两个材料空间阶失败；下文保留实施前推导及被动下界审计依据，不表示生产实现通过。[任务书](H10g_Point_Pde_Material_Moments_TODO.md)的修改预算、原始门槛和失败退役规则保持。H10f冻结148个payload及当前交付均已独立校验；旧App/PDF、Interface与H08/H09不变。最终QP未实施，不进入H11/H08a。

## 从真实PDE误差而非通量率权重出发

同一H09质量和唯一材料通量给出`a=(-N+u*q)/M`。若最终`N=q*u/2+K*u`、`K^T=-K`且`K*1=q/2`，则常数保持；一次局部增量的导数为`-sum_j K_ij*delta_ij/M_i`。点采样PDE要求该一次目标为`bPDE=M*V`，V只由原H07实际总量率/面测度和真实局部functional恢复。几何flux-centroid矩不自动等于这个目标。

单周期闭合端口的横向坐标z无跳跃，其一次目标b必须满足两个必要条件：

`sum(b+z*q)=0`，`sum(z*b+q*z*z/2)=0`。

周期等值native class内，q、质量和系数须先按实际拓扑求和；逐native散射仍按原真实质量。对于本轮横向probe，各alias的z一致，因此直接native求和与quotient求和相同。轴行仍参与代数和审计；接受轴速度为零，其横向一次目标和修正强制为零，不删除该行。

设`delta_a=-delta_b/M`。以真正二阶PDE余项的质量范数定义目标：

`min 1/2 sum_i M_i*(delta_a_i/h_i^2)^2 = 1/2 sum_i delta_b_i^2/(M_i*h_i^4)`。

这唯一给出权重`W_i=M_i*h_i^4`；固定轴行W=0、delta_b=0。令`c0=-sum(bPDE+z*q)`、`c1=-sum(z*bPDE+q*z*z/2)`，`zbar=sum(W*z)/sum(W)`、`D=sum(W*(z-zbar)^2)`，解为

`alpha=c0/sum(W)`，`beta=(c1-zbar*c0)/D`，

`delta_b_i=W_i*(alpha+beta*(z_i-zbar))`。

这是指定范数下的等式最小改动，有明确KKT条件`delta_b/W=alpha+beta*(z-zbar)`，不是调节失败惩罚。sum(W)>0且D>0才可识别；用centered坐标避免Gram大数抵消，记录`D/(sum(W)*L^2)`。零Gram且非零右端明确FAIL，不用正则化或fallback。

## 必要条件通过不代表二阶或局部可行

修正实际导数为`delta_a=-h^4*(alpha+beta*(z-zbar))`。因此需观察`delta_a/h^2`，而不是把H10f率权重的C_aff=8挪到另一公式上。H10f四项失败裁决原样保留。本轮没有提高旧门槛，也不以新范数重标旧试验。

在固定两级尺度h_i~h且真实正质量覆盖非退化横向区间时，`D~sum(M)*h^4*L^2`。若两个相容缺陷为O(h^2)且质量归一化缺陷有界，公式给出O(h^2)导数余项；若缺陷慢于二阶或Gram退化，则该设计本身不能满足目标。必须从真实48组数据检验这些前提，不能用该条件性推导宣称普适二阶。

第一阶段固定审计XY/RZ、uniform/2:1、n=4/8/16/32、rho=1/1000/1e6：两条原gross等式相对误差<=1e-12、独立rank-2解及KKT；rho=1的零PDE控制<=1e-12；真实质量RMS、全部非固定native Linf、`max |delta_a|/h_i^2`及加密阶。数值成功要求原空间标准：n16→32 RMS阶>=1.7、最细RMS<=0.05且Linf下降。n4→8→16也完整报告，不把抵消或输出exit 0当作空间PASS。未通过则归档/退役这一候选，不进入最终QP和完整H10b。

## 若上述目标通过，再实现最终one-ring算子

通用中心共享系数`Kcenter=sign*F/2`为最小改动锚点，不再对segment trace作补丁。真实几何和显式端口生成唯一closed one-ring及周期alias。硬常数与相容一次矩联合装配到最终共享skew pair；周期方向用wrapped delta和正确winding项，不能把横向rank-2公式直接推广到周期坐标。全局相容条件只必要，局部rank可行性必须独立验证。

Taylor尺度预先定义`S=s+(|M*V_x|+|M*V_y|)/h`，它只用于算子预算，不改变F/q/R。行l1<=2*S，归一化二次矩绝对值<=S*h/L；另外独立记录三次矩/h^2与行l1/(M/h)的实际加密尺度。P2的O(S*h/L)配合真实M下的有界输运尺度才控制二阶误差；单独通过矩slab不足。均匀中心证书原s预算仍单独检查，不能覆盖旧证书。

约束最小改动使用硬矩SparseQR与null-space dual工作集：保持非负乘子，加入实际最违反slab/l1 sign cut，沿维持现有active等式的方向走到新约束或先释放到零乘子。最大步数`max(10000,4*constraints)`。硬矩、全部不等式、工作集乘子和stationarity/KKT必须独立复算；无可动方向须导出可验证不可行证书。禁止恢复H10e失败Dykstra/罚函数或增加迭代当作设计修复。

旧opaque API和原124场景先保持。新显式端口API及actual apply须真实来源/clock拒绝、KE/动量恒等式、内部affine n8/16/32、全部行/材料/系数尺度通过后，才原28窗口/15,116步/84gate完整重跑。原4个donor gate、KE、时间和空间门槛都不豁免。尚未实现这一阶段。

## 第一阶段预审后的独立下界与材料来源审计

质量权重候选预审在XY 2:1、rho=1000/1e6的n16→32实际RMS阶约0.80/1.11，未通过原1.7；不进入生产QP。以下新增的是被动证据，不是另一修正路线。

任意满足两条相容条件的实际一次PDE误差e，都有`sum(M*e)=-c0`、`sum(M*z*e)=-c1`。故其最小真实质量平方范数为`c0^2/sum(M)+(c1-zbar*c0)^2/sum(M*(z-zbar)^2)`；除全部native真实质量后开方得到RMS下界。以M*h^4代替M同样得到`e/h^2`范数的下界。只输出这两个标量、Gram与原缺口，并用独立2×N SVD最小范数及Cauchy–Schwarz检查核验。它们不构造新target/pair/通量/状态。固定轴不删除，等式仍含轴，其实际误差限定为0。

**有限网格下界的比值不能证明任何满足约束的算子都不可能二阶**：实际误差可以高于下界，其比值不受下界比值直接限制。非零c0/c1足以否定这组材料的“全部行精确affine+原skew”，但是否有稳定有界二阶算子仍需最终局部约束与渐近证据。

H09已冻结的逐dual恒等式为`q_u=sum_p overlapFraction_up*qPrimal_p`。定义实际映射的`xi_p=sum_{u in component0} overlapFraction_up*z_u`、`eta_p=sum overlapFraction_up*z_u^2`。两条缺口分别等于`-sum bPDE-sum xi*qPrimal`、`-sum z*bPDE-sum eta*qPrimal/2`。再将qPrimal展开为真实primal面F的带符号scatter，独立复算缺口；输出所有pressure行和与cell-center坐标的差异。该审计能区分映射空间矩与真实材料F对缺口的贡献，不改任何H08/H09定义，不把它自动归罪于Interface或非守恒材料路线。
