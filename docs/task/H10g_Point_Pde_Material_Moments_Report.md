# H10g / 真实质量误差范数与材料矩来源报告

2026-10-04，**FAIL / H10未关闭**。M*h^4相容目标通过代数和最小范数核验，但两组原空间阶未通过，已退役；没有新生产operator/API，生产与原短契约保持H10c。保留48组被动必要下界和材料来源观察器，下一单任务[H10h](H10h_Final_Operator_Joint_Feasibility_TODO.md)。[任务书](H10g_Point_Pde_Material_Moments_TODO.md)、[实施前设计](H10g_Point_Pde_Material_Moments_Design.md)。本轮是有效设计证据进展，不是外部BLOCKED。

## 候选推导与明确裁决

从真实二阶导数余项范数`sum M*(delta_a/h^2)^2`推导的唯一权重是M*h^4；它同时补偿两个必要条件缺口c0/c1，固定轴修正为0。没有提高H10f的C_aff=8，原四项失败保持。候选只在已保存的同源native数据上作rank-2解析评估，未进入convection系数或actual apply。

独立非centered 2×2 Gram解与研究观察器的centered Gram/SVD相互核对。48组、46,680个native行、22,440个pressure来源行完整；独立最大参考误差1.93e−13、相容等式误差5.96e−16、KKT误差9.43e−17，均<=1e−12。rho=1及其余横向零PDE组合在舍入尺度通过。两条必须同时成立的空间条件仍FAIL：

| XY真实2:1材料 | n16 RMS | n32 RMS | n16→32阶（须>=1.7） | n32 Linf | n32 max误差/h_i^2 |
| --- | --- | --- | --- | --- | --- |
| rho=1000 | 7.9400132e−4 | 4.5575032e−4 | **0.80090 FAIL** | 0.01316243 | 13.47833 |
| rho=1e6 | 5.5582004e−4 | 2.57013896e−4 | **1.11277 FAIL** | 0.02629664 | 26.92776 |

最细RMS<0.05和Linf下降不能豁免空间阶；12个family验收中两项失败，不是原完整84gate的重跑。全部n4/8/16/32数据都保存，不能用n8→16的较高阶或48组观察器exit 0代替空间PASS。公式、逐行目标witness和退出1裁决只保存在冻结的`rejected_mass_norm_target`。

## 不是继续寻找投影权重

研究生产观察器只输出任意相容误差必须满足的真实质量RMS及二阶缩放RMS下界。由两条约束的Gram逆/Cauchy–Schwarz得到，独立2×N SVD核验；最大参考误差6.08e−15。全部Gram可识别，最小归一化横向variance=5.72e−6，实际值完整保存，不假定高密度比下自然良态。

XY 2:1 rho=1000的最小真实RMS下界，n16为4.8002953e−4、n32为2.3605202e−4。非零缺口排除该组材料“全部native行精确affine+原skew”；**这些有限网格下界不能证明所有二阶算子都不可能**。具体候选失败与最终one-ring有界交集是否可行是不同结论；下一任务须直接在最终共享pair上检验PDE误差与矩预算，不能再预先投影一个目标然后宣称算子一致。

## 缺口已独立追溯到原材料输入

冻结H09逐dual恒等式`q_u=sum overlapFraction_up*qPrimal_p`给出实际映射xi/eta。独立D_I方向scatter与直接真实primal面owner scatter都复现两条native条件，最大来源误差4.32e−16。D_I矩阵已含面积、primal faceRates已是积分率，所以观察器只使用D_I的方向；初次审计错误地重复乘面积，已由来源断言拒绝并修正，错误版本归档，不计通过。

XY 2:1映射的一次坐标xi与pressure cell-center完全相同（最大差0）。rho=1000的c0在n8/16/32为−0.24086049、0.06874287、0.02219143，可由**原primal材料F的空间矩**完全解释；该缺口在convection系数求解前就已存在。二次坐标eta与center²的最大差为h²/32：n8/16/32分别4.8828125e−4、1.220703125e−4、3.0517578125e−5。它说明H09的交换守恒恒等式与空间材料矩是不同约束，不说明非守恒FirstJet路线应被废弃。

观察器不改Interface、材料质量/F/R、状态或接受条件。没有把来源差异变成AMR分支、解析载荷修正或另一套材料F，也没有用“全局矩可投影”代替局部最终operator可行性。

## 回归、未运行项与来源

标准Release构建通过；4/4受影响CTest通过：architecture、architecture_self_test、repository_layout、cch_kep_convection_contract。旧124场景/1,516独立矩及SUMMARY与H10f完全一致；内部n8/16 affine审计3,760行最大1.58492e−13，原XY AMR n16 dof911误差7.41629e−14。n32旧审计本轮未重跑，生产与原审计数值路径未变，历史证据保持。

前置PDE目标失败，未实施最终QP、新端口契约或重跑原完整H10b；仍以[H10c原28窗口/15,116步/84gate](H10c_Native_Convection_Moments_Report.md)的14失败为完整研究裁决。未进入H11/H08a。源预算2 production+1 test+1 research实际只改1 research；两个生产文件和原短契约逐字节不变。原Study/rhs/full/RK4/donor/gates/rowAudit/affineAudit和原材料观察字段由逆向删除新增被动代码独立核对，不调门槛。

启动时间2026-10-04 05:36:29+08，检查目录`runs/checks/H10g/20261004_053627_707`；冻结`runs/baselines/H10g_20261004`。H10f148个payload、当前18份交付、1,042份源码及346份输入启动时全部验证；本轮受保护源码只driver改变，输入只driver/README改变。原App SHA256=7F0C61DDE621C2B5967BC7ACE32DFC7E5D09516D5029718470D95587F46CBA5B；冻结CCH PDF SHA256=7E9803EE0E339AC58DBD77B0DA71BB1F0E385798CB2F7D2C6EC34EBB132F7207。最终研究exe SHA256=F52CA4BD96B01501CA5EFEF6480586BCD856E5519B17E45707DEB463BB19FDDE，绑定source与48组证据见run/manifest/verification。
