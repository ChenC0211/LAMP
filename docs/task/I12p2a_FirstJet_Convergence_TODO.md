# I12p2a / 纯 Hermite–GALS 收敛与步内跨面

目标：验证已归档 I12p1 纯路线的空间/时间/曲率精度及非对齐步内 crossing，不改 production 算式。

允许修改：tests/interface_hermite_gals_benchmark.cpp、tests/CMakeLists.txt（0 production / 2 tests），本任务文档和证据。

禁止修改：tracer、state/publication、I05/I06、几何、体积/守恒算法、原算例参考、冻结资料及其他模块；不恢复 M、体积修正或重新初始化。

数学约束：所有非初始内部值来自前一步唯一 scalar 的实际 first-jet 回溯；解析场仅初始化、真实外边界与误差参考。纯体积保持诊断，不用它调节几何。空间/时间误差按独立解析 jet/curvature 计算。

实现要求：复用 I12p1 原纯步与 publication。静态非多项式 curvature：XY/RZ uniform/AMR，n=8/16/32，原梯度/曲率阶门槛 2/1.5。非多项式平移空间研究复用原 translationSpaceStudy 的 n=8/16/32、dt=.01、4步（终止时间 .04）；保留四组 XY/RZ uniform/AMR，检查末档梯度/曲率阶 >=2/1.5。旋转时间研究使用原 n=8、t=.3、dt=.01/.005/.0025，uniform/AMR，末档位置阶 >=1.8。需明确空间误差/舍入底限，不伪造阶数。

跨面：原 48 组 XY x/y 与 RZ z、uniform/AMR、速度/相位正负、1/4 和 3/4 步内穿越 .5 面；直接纯步生成 midpoint/end 候选，逐 cell 几何体积 <=5e-12、first jet/normal <=1e-11、curvature <=1e-10，与原解析标准相同。没有守恒通量或库存 ledger 声称。

验收：T0 只有2 tests修改，旧源码/runtime及 freeze81/PDF3/p1 archive保持；新测试编译无warning，原28 CTest保持。T1/T2上述三档收敛/48crossing全执行，任何失败据实记录。动态 refinement/coarsening 属后续 I12p2b，不能以静态 AMR替代。p2a PASS 不等于完整 Gate I。

失败规则：内部重置解析场、修改容限迎合结果、删失败组或未完成三档均 FAIL；内部精度问题不标外部 BLOCKED。

交付：测试实现、原始误差与阶数、回归/审计/归档、报告及下一阻塞。

结果：**PASS**，见[报告](I12p2a_FirstJet_Convergence_Report.md)。12 静态 + 12 平移空间 + 6 旋转时间 + 48 crossing 通过；原 28 CTest/输出保持。动态 regrid 留给 I12p2b，Gate I 未通过。
