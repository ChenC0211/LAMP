# I12a16r7 / 舍入修复遵守原面通量证书

目标：解决r6中三个原算例第二步full在原舍入修复后越出面积分认证区间的问题；仍使用同一个库存装配kernel及原64epsilon预算。

允许修改：HermiteGalsTemporalPhaseFlux.h/.cpp、Detail/HermiteGalsInventoryAssembly.h（3production）；原integrated_inventory_contract及hermite_gals_benchmark（<=2tests）。禁止修改Core circulation、actual tracer、速度/相积分、独立M/几何闭合、reflux算法、原case/网格/dt/步数/容限、冻结和其他模块。

数学约束：始终使用唯一共享面通量，M_new=M_old+入-出，严格0<=M_new<=C；不clip M。新证书只缩小原舍入修复的允许变化集合，不能增加64epsilon*min(C)额度或128次修复上限。total和boundary固定；最终liquid同时满足原认证上下界与[0,total]。几何按原精度预算闭合，不恢复逐cut cell机器精度体积匹配。

实现要求：先写 [设计](I12_Certified_Roundoff_Repair_Design.md)。原private finish新增可选完整面bounds；增长分支用证书upper，减小分支只消费liquid-lower的可减余量。无证书分支保持原表达式顺序及行为。两个公开assembly入口共用原private装配/发布路径；不在benchmark复制修复算法，不增加公开可伪造积分证书。没有合法余量时明确拒绝，不按小误差删掉来源pending。

验收：T0 3production/<=2tests，39targets/26CTest、81freeze/3PDF/r6archive保持，无新增warning；旧无证书数值结果不变。T1正负方向、空/满库存、非零lower、极小通量、固定boundary、可行替代修复与无可用余量拒绝；严格证书/守恒/原舍入预算同时核验，旧integrated contract完整前缀保持。T2原四个smoke及原full30/400步实际执行，要求原三个full证书越界被解决；如后续来源认证或其他真实gate失败，原样记录，不将局部修复冒充I12通过。

失败规则：证书越界、边界/total被改、扩大舍入预算、库存clip或旧输出意外改变均FAIL。内部数值缺口不标外部BLOCKED。一个task验收归档后才开始下一项，I12完整gate未过不接Hydro。

交付：设计、代码、算子/原驱动与回归证据、报告/审计/归档、未解决问题。

结果：**PASS（本项）**，见[报告](I12a16r7_Certified_Roundoff_Repair_Report.md)。四例第二步通过；第三步来源认证仍失败，I12 未通过。按用户后续授权，归档后转纯 Hermite–GALS 非守恒路线，不再扩展此守恒通量分支。
