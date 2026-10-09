# I12a16r6 / 认证面通量区间与库存可行性

目标：解决原rotation/vortex uniform/AMR第二步半步库存越界；首先判定它是否能在各共享通量的原认证误差区间内消除，避免把实际离散不相容误当舍入问题。

允许修改：HermiteGalsTemporalPhaseFlux.h/.cpp、Core/Discrete/BoundedCirculation.h/.cpp、根CMakeLists.txt（5production）；原hermite_gals_benchmark和bounded_spacetime_flux_contract（2tests）。诊断完成后、production实现前明确此范围：区间内共享选择需要中立的有界循环流及独立严格浮点编译单元，不能在benchmark复制算子或在非严格TU内使用exact expansion。先只加基于原bounded result的诊断，再依证据决定production改动；若需要区间外的输运方法变更，另写设计/任务，不在本项扩大范围。禁止改actual tracer、速度evaluate、旧库存装配/64epsilon修复规则、I10/I11、benchmark案例/步数/网格/容限、冻结和其他模块。

数学约束：M独立守恒、严格0<=M<=C；Hermite几何只满足既定精度预算，不恢复逐cut cell机器精度匹配。每个native face只能有一个共享通量，边界history/source及coarse-fine reflux保持。任何协调必须在原certified gross bounds及liquid<=total内，原误差界仍成立；不得clip M、增加旧舍入修复预算或随意改动boundary通量。必要设计必须先于production实现。

实现要求：记录原失败cell的M、C、wl/wg、几何体积、所有相邻gross通量及其四个区间、原代表值净更新、完整区间净更新的向外上下界和原64epsilon预算。若区间不可能满足库存界，给出明确证据并保持拒绝；不得选择区间外通量冒充积分准确性。若确实可行，再实现保持sameInterval及唯一面所有权的有界选择，并验证原装配严格通过。

验收：T0范围/原39targets/26CTest、81freeze/3PDF/r5archive保持；原未涉及数值结果不变。T1诊断向外界与独立小矩阵参考一致，极端库存/方向/边界/AMR覆盖；若增加协调，验证守恒、有界、原certified bounds、确定性、source/provenance与不可行拒绝。T2原四个smoke和原full30/400步实际执行，不隐藏失败；记录新accepted时间及真实失败位置。只有原第二步越界被解决才将本项修复记PASS；若需要另一个数学机制则本项FAIL归档原因证据，继续下一项。I12原完整gate不变，未通过不进入Hydro。

失败规则：超区间修改、库存clip、虚假complete、容限放宽、旧行为意外改变均FAIL。内部实现/方法缺口不标外部BLOCKED；保留原accepted状态并继续证据驱动的必要工作。

交付：TODO、诊断/必要设计与代码、测试及原full证据、报告/审计/归档、未解决问题。

执行结果：**FAIL（完整修复未完成）**。按 [设计](I12_Certified_Inventory_Circulation_Design.md) 实现原认证区间内的共享选择，8网络/12真实积分批次/14guards通过；独立BigInteger复核41+82个诊断区间。原四个第二步half通过，只有AMR vortex第二步full通过；另外三例在原舍入修复越出极小面证书时保持拒绝。39targets及26CTest通过、旧数值结果保持，原full及I12仍FAIL。下一项给原舍入kernel添加证书约束，见 [报告](I12a16r6_Certified_Flux_Inventory_Admissibility_Report.md)。
