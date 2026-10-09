# H10k / 实施前唯一材料 owner 设计

2026-10-04。前一轮为 progress；H10j 的289 payload/101 canonical及 j0/i2/f/e/d/c归档已独立复核。按[任务书](H10k_Material_Spatial_Owner_Migration_TODO.md)只改现 CchMaterialTransport header/cpp、现 phase/material 两份契约及既有研究 driver/文档；其他源、输入、App/PDF保持。

原 prepare 对 spatial midpointRate 保留同一 actual advector/momentumGeometry/layout/D_I，逐 atom 唯一计算 F=rho_l*Q_l+rho_g*Q_g。两相乘积分别有限性/representability校验；primal和dual视图只按原 ordered atom indices用补偿求和。dt积分仍为原 midpoint近似。原 pressure三层几何质量、实际质量差、R和全域 balance定义保持，不把R=0/quarter/材料commuting加回约束。

新增 isSpatial/momentumGeometry/atomRates/dualFaceRates只读port，spatial专属getter在legacy结果上拒绝。CchMaterialSources增加默认无效的CchAdvectorIdentity advector；spatial来源由原 midpointRate取得，legacy暂只有0标记，不能把0当实际advector。旧调用者三geometry aggregate仍能构造legacy stamp，但spatial validateFor必须精确比较actual stamp；另加带expected momentumGeometry的validateFor以拒绝等坐标fresh实例。所有prepared复制共享原数据，无setter/source投票。

prepare仍接受任务书明确暂留的legacy midpoint输入，只执行原面表达和原回归；输入类型在构造前已确定，新spatial失败不得转入legacy。这不是两种production选项，最后旧dual/KEP两个测试在H10l迁移后删除旧入口与表达。本轮不把现dual quarter视图冒称新F，也不修改其算法或KEP。

phase契约移除旧evaluate调用，以同一actual factory核验geometry/clock/axis/lifetime及独立physical polynomial切割积分；普通geometric测度和收敛门槛保留。旧A_phase*u_face断言改为实际积分，不要求新varying trace满足旧常值算法。material原constant/axis-regular族仍有独立闭式相面积参考，basis输入保留但独立从actual polynomial和解析相交重算所有真实面/相流：零面平均不意味零材料率，因此不再断言basis仅影响两个cell。保留其真实邻接求和/全域边界身份和静态质量差0。

追加source-advector错配、expected几何实例、复制/生命周期、非有限/细相、重合迹正负/零/局部变号和拒绝事务。实际XY/RZ uniform/2:1、full/selected、rho1/1e3/1e6、正式三层移动/crossing及pure-gas21段研究。原56个不可识别tangent仍明示不称二维affine精确，但finite/source/sharing不可删行。

研究的两处材料准备helper消费原momentumGeometry上的actual advector；native advector观察直接读该owner，不再从material q重新拟合另一V。新增material观察只写实际atom/view/三层压力质量/R/来源与耗时，独立后处理重算物理积分及材料balance，不改变full/RK4/donor/gate。旧CSV若变化按实际范围解释；内部affine必须重跑，完整H10b留到最终native收口。PASS只关闭H10k，下一H10l，无H11/H08a。
