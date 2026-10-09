# H10l / 实施前 same-F dual 与缺陷分解设计

2026-10-04。前一轮progress；启动独立复核H10k的229 payload/109 current canonical及j/j0/i2/f/e/d/c归档通过。按[任务书](H10l_Same_F_Dual_Material_Migration_TODO.md)固定5 production/2既有tests，现research driver/文档，无新源码/target/tools。检查目录`runs/checks/H10l/20261004_144355_065`，before全源/输入/原文件与App/PDF已保存。

唯一材料F仍由H10k Flux的canonical atoms计算，H09的faceRates只返回原Transport.dualFaceRates引用，不另存/按面积分配/重新乘密度。原momentumGeometry/三层真实dual质量和完整clock保留；prepare先核验原三层质量及原layout/D_I/密度，再用expected midpoint.geometry精确核验原spatial材料instance，validateFor也使用它。actual advector stamp随原material传播，不贴给旧数组。

删除CchDualMaterial中Local/quarter分块、tensor猜测、SPD延拓及segmentRate。原真实side incidence计算q_dual和dt*F。由原D_I方向及同一primal F计算q_primal/gross，仅作诊断，不另建D_I。原overlapFraction映射T仅作用于massChange、integratedDefect及q_primal。

cell新增outwardRate、mappedPrimalOutwardRate及spatialExchangeRate，E_sp=q_dual-T*q_primal。G=delta_m_dual-T*delta_m_primal；r_dual=delta_m_dual+dt*q_dual。硬身份改为r_dual=T*r_primal+G+dt*E_sp，按原gross≤256eps核验。不要求E_sp/G/R=0。report删除无消费者的maxQuarterExchange，旧maxDualExchange改名maxRelativeSpatialExchange并保留实际非零诊断，另有maxAbsSpatialExchange；绝不能赋零继续冒称旧交换通过。consistencyImpulse/materialKineticChange的算法和定义不改。

原CchPhaseRate::evaluate仅余dual/KEP两tests，迁移后删除声明/实现。Data保留同一spatial对象，不再optional legacy结果，isSpatial为恒真兼容port；prepareSpatial实际数值积分逐字保持（仅Data构造变更）。MaterialTransport旧prepare只编排唯一Flux，不再legacy density-weighted face表达；private assemble因固定header预算保留签名但强制nonnullpreparedFlux，Data仅持有该Flux，无legacy rates数组/optional分支。

dual/KEP两tests helper只在原mass.geometry上准备actual B和phase spatial rate。原basis/linearity、geometry/时间/动量/能量/寿命/regrid拒绝保留。旧exact exchange与primal面积比例断言换为独立atom physical积分/实际view共享、真实incidence与含dtE_sp分解；basis参考允许读实际kinematic polynomial作为输入，不能读生产F当reference。独立可识别MMS另有声明系数参考、zero正负/零/变号和自身尺度小相检查。未识别tangent仍原56粗selected声明，不新增豁免。

研究增加只读dual模式，输出actual shared segments、三层dual几何质量、true R、mapped量/G/E_sp、uR/材料KE、sources/完整clock与gas witness；参考不参与F/operator。原full/RK4/donor/gates不改。几何/advector一致族才讨论R收敛；不同规定运动族仅作强制几何身份检查。实际8/16内部affine先重跑，若迁移新F使最终native失败，保留witness并拆final operator任务，禁止在本轮改KEP production或gate。

关闭仅H10l，完整H10b/H11/H08a仍待final native联合空间/KE验收。Pure FirstJet独立，无库存/bubble/mass修复；四材料硬约束保持，其余诊断真实报告。
