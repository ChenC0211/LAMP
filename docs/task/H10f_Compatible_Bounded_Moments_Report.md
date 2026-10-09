# H10f / 可行证书与真实材料PDE一次矩

2026-10-04，**H10f设计候选FAIL / H10未关闭**。均匀网格已有明确有界可行系数，H10e该处失败不能再解释为交集不可行。真实材料的flux-centroid几何一次矩则不对应所需点采样PDE；本轮几何修正和后续率加权PDE修正均未接受、已退役。生产与原短契约始终保持H10c字节，只保留被动证书/一次矩观察器。下一单任务[H10g](H10g_Point_Pde_Material_Moments_TODO.md)，不进入H11/H08a。

## 原固定有界交集的正面证据

按[实施前设计](H10f_Compatible_Bounded_Moments_Design.md)，从实际H09同一F构造一般中心KEP共享系数`K_ij=sign*F/2`，明确闭合研究周期段和native重复采样，真实half-dual质量保留。没有解析场进入系数，没有按AMR类型选择数值算法。

`certificate/certificate32`保存12组n=8/16/32、XY/RZ uniform/中央2:1的15,120个quotient行及15,424个非零pair。独立从保存的pair重组矩阵，检查原H10e常数/一次目标、中心二次区间`s*h/L`与行l1界`2*s`，与保存行数据逐值相同；source另断言每个非零pair属于真实闭合one-ring。结果见[certificate_summary.json](../../runs/baselines/H10f_20261004/certificate_summary.json)。

六组uniform均可行：XY n=16硬矩和二次超界为0，l1比率0.25；RZ最大硬矩gross误差约4.04e−17。六组AMR的中心seed不满足原一次/P2预算，不能作为其可行或不可行证书。这里证明uniform交集存在，不证明最小改动已求出、所有AMR可行或最终PDE已经通过；观察模式exit 0表示诊断完整，AMR非可行seed仍明确输出0。

## 几何矩不等于点采样PDE

真实导数`a=(-N+u*q)/M`在affine探针下为`-b/M`。以flux测度centroid装配的b并不自动等于`M*V`，即使全局Psi为0也可能有大局部偏差。研究输入仍是H10e的48组真实moving FirstJet/H09材料：n=4/8/16/32、XY/RZ uniform/2:1、rho=1/1000/1e6、dt=0.1/n。XY取w=(1+y,0)、V=(0.4,0)，RZ取w=(r,0)、V=(0,0.4)，横向导数参考为0，固定轴按原约束处理。

新增`material_affine.csv`逐native记录真实M/h、几何b、b/M、q/率gross，以及从原`totalFaceRates/faceMeasures`和可识别one-ring functional恢复的实际advector。48组共46,680行，速度恢复最大差2.78e−16；从保存行重新核算所有RMS/Linf及两个PDE相容条件，与summary相同，原H10e全部15个材料字段逐字保持。独立核验脚本只读证据，不装配/修改接受状态。

rho=1000的**几何目标对应导数**如下（不是本轮对新convection实际apply的测试）：

| 目标的质量RMS | n=4 | n=8 | n=16 | n=32 |
| --- | ---: | ---: | ---: | ---: |
| XY 2:1 | 0.299356 | 0.250173 | 0.0535031 | 0.0497217 |
| RZ uniform | 0.416000 | 0.237309 | 0.185245 | 0.102287 |

XY最细两档阶仅0.106；RZ uniform Linf由11.74增至22.51。rho=1对照近舍入精度，说明只验证均匀材料会遗漏此问题。真实质量与实际横向速度生成的`bPDE=M*V`在这些探针上为0；但它还必须满足同一q/skew的全局条件，不能直接宣称精确联合再现。

两个必要条件至少包括

`sum(bPDE+z*q)=0`，`sum(z*bPDE+q*z²/2)=0`。

当前观察器同时保存它们的缺口c0/c1，区别于几何Psi。结果见[material_summary.json](../../runs/baselines/H10f_20261004/material_summary.json)及逐行CSV。目标、实际apply、完整时间积分是三个不同证据范围；本轮只直接核验前者与原内部反例。

## 未接受的修正及退役

首先用实际`s*h³`权重消除几何Psi，XY rho=1000 n=8 RMS只从0.250173降到0.248180；RZ局部偏差同样未解决。因此这条几何目标修正不进入生产。

随后在生产实现前明确改为真实`M*V`目标，同时修正上述两个必要条件，使用率/点质量速度尺度`S=s+M*|V|_1/h`及`S*h³`权重。这个S只用于算子预算，没有修改唯一F/R/uR。候选在事前固定`C_aff<=8`下有4组失败：XY 2:1 n=8和n=32的rho=1000/1e6。n=8指标10.2004/13.0728，n=32为8.28450/10.7211；n=8 rho=1e6的导数Linf约2027.49。完整[候选裁决](../../runs/baselines/H10f_20261004/rejected_rate_weighted_candidate/verdict.json)及source/runtime关联保留于同目录。

初次n=32诊断在rho=1e6的投影身份检查中因使用小权重/variance作误差分母而拒绝；改用实际等式项的gross尺度，原1e−12门槛保持，完成12组。候选数值/预算并未因此通过；与n=4/8/16的36组构成完整48组判定。没有提高固定8，没有实现最终QP、任何生产端口API或新增配对短契约，也没有启动一般Hodge。两条修正的数值代码均从现行driver删除，失败source/binary仅在归档。

这些结果排除本轮的目标/加权选择，**不排除正确点采样PDE与skew的二阶方案**。下一步需要以真实质量及二阶PDE余项范数推导相容矩，不能再以一个全局功残差小来替代点导数与实际h阶。

## 现行验证、冻结和范围

标准Release构建通过；原短契约及架构/自测/目录共4/4通过。124场景、1,356 trace、6,144 load、24实际跨界、19拒绝、4 regrid、1,516内部矩/3,504pair的数值SUMMARY与H10e逐字相同。原1e−12载荷/动量/KE及1e−10常速度门槛保持。

原n=8/16八组合affine再次运行，内部最大1.58492e−13，dof911误差仍7.41629e−14。n=32内部6.59128e−13来自上一轮H10e的完整审计与同一生产身份，本轮没有重新运行audit32。均匀可行证书与真实材料一次矩均实际运行到n=32；不能将它们改称新operator实际apply通过。

**本轮未重跑原完整H10b**，因为相容PDE目标候选失败且没有可验收的新最终算子。最新完整研究仍是H10c的28窗口/15,116步、84gate中10空间+4donor失败；所有原门槛保持。无H11、H08a或完整CCH流体App。

冻结目录[runs/baselines/H10f_20261004](../../runs/baselines/H10f_20261004)保存原133个H10e payload重核、当前交付、全部诊断和独立verification。1,042受保护源码只有研究driver改变，346输入资产只有driver/README改变；2 production、原test、Interface、材料定义、CMake/工具、旧App/PDF保持。可逆scope核验保证Study/rhs/full/RK4/donor/gate/rowAudit/affineAudit保持，原材料字段的身份另行核验。

H10g先完成PDE材料矩的误差范数与相容性设计，再决定共享one-ring有界QP。原内部反例与真实材料/全部行前置通过之后才完整H10b；Pure FirstJet独立、无库存/bubble，空间离散通过之后才H08a。
