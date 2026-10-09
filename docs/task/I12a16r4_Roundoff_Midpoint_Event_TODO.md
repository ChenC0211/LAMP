# I12a16r4 / 原回溯的舍入分辨率中点边界事件

目标：修复原合法CFL下近壁内部点因中点不可采样而耗尽128次缩步的运行缺口；在原已有时间/坐标事件分辨率内解析物理边界命中，仍使用同一actual tracer。

允许修改：include/Core/Geometry/LevelSetBoundaryTrace.h（1production）；原hermite_gals_trace_contract与bounded_spacetime_flux_contract（2tests）。禁止修改速度evaluate、区间/相积分、库存装配/reflux、I10/I11、冻结、其他模块及已有容限。

数学约束：完整内部RK2、Jacobian和正常target-event表达式顺序不变。中点仍用原P−u(P)*.5*d，不在域外采样。当中点不在域内，且d/2<=原timeTolerance、两点最大坐标差<=原coordinateTolerance时，复用原舍入bracket的直线交点规则，返回边界事件；事件时间允许舍入到当前elapsed，但绝不接受零时间内部substep。继续原物理域、CFL、axis、activeSegment、唯一侧及边界入流速度检查。命中不返回内部Jacobian。

实现要求：先写独立设计说明，明确这是原舍入极限事件规则的扩展，不宣称连续ODE精确first-hit或从中点得到原full-step target。共享原bracket交点helper，原正常路径不变；新事件标志便于验证。原128/96循环限及容限不增加。切向、非入流、多壁、axis负半径或coverage违规不得被snap通过。

验收：T0 1+2范围、39targets/26CTest、81冻结/3PDF/r3archive保持。原已有编译警告按冻结源/原记录区分，不允许新增警告。T1 XY四边/RZ南北、非零壁坐标、nextafter/极小/亚正规距离、时间下溢、两callback路径/速度jet一致，无域外query；同一原解析常速事件参考，物理time、无边界J及有限callback计数；真实切向/多壁/CFL/axis/coverage拒绝，旧trace完整输出保持。T2原面5的两条y=1e-100实际轨迹须转为有来源的south事件；其余7hit/5interior保持，原同四个driver smoke/full仍实际运行。边界混合通量证书未接入时仍保持a16/I12 FAIL，不能将事件修复冒充完整驱动通过。

失败规则：容限外触发、域外sampling、虚构边界切线、旧物理结果变化或违规被放行记FAIL。若暴露不可解析数据仍明确拒绝；内部实现缺口不标goal BLOCKED。没有新增数学模型、边界相或内部解析重初始化。

交付：设计、代码、旧/新轨迹证据及完整驱动结果、report/审计/归档。

执行结果：**PASS（本项）**。66近壁事件+1子步后事件、12零时间命中、8类拒绝通过；原两个polynomial复现点转为south事件。原I02完整输出保持；原积分只改变最后的近壁诊断。39targets及26CTest/审计完成，25PASS+原driver真实FAIL，24旧完整输出保持。11条历史C4267已按源hash和原日志核验，新增warning0。见 [报告](I12a16r4_Roundoff_Midpoint_Event_Report.md)。混合相通量覆盖仍待实现，a16/I12保持FAIL，goal继续。
