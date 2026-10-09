# 近壁中点事件的舍入分辨率处理

依据 [开发总纲](../00_LAMP_Development_Charter.md) 的设计说明要求，先说明本次有限精度边界规则，再实现 [r4 TODO](I12a16r4_Roundoff_Midpoint_Event_TODO.md)。原 I02 的完整内部离散RK2及其Jacobian保持；本项扩展的是实际tracer的边界事件舍入极限处理。

## 原实现不能解析的情况

[r3报告](I12a16r3_Sign_Preserving_Velocity_Bounds_Report.md)已用同一polynomial source实测：合法half-cell CFL下，面5中y=1e-100的两个内部点因midpoint仍在域外而耗尽128次缩步。仅增加循环限也无法覆盖所有浮点极限：很小的正距离和较大的有限法向速度可能需要小于最小正double的事件时间。不能在域外求值、clip/reflect RK stage、伪造内部Jacobian，或让积分在未定义trace上宣告成功。

## 沿用原事件分辨率

原tracer已以128epsilon*dt的时间容限及128epsilon*domainScale的坐标容限隔离事件，并在最终舍入大小的inside/outside bracket间求直线交点。新分支沿用这两个容限，不提高它们：若原midpoint H=P−u(P)*.5*d 在物理域外，d/2不大于原时间容限，且max|H−P|不大于原坐标容限，则P与H作为舍入分辨率内的事件bracket，按原同一几何交点公式解析命中。该H没有被采样，也不作为full-step target或合法内部stage。

hitTime=elapsed+(d/2)*fraction（实现保留共有bracket表达式顺序）。时间可能舍入到elapsed；由于立即返回边界事件，无内部零时间step被接受。投影到物理壁上的点仍须通过原activeSegment、唯一侧、速度CFL/轴和真实边界入流方向校验；多壁或切向等不支持情况仍拒绝。返回命中时没有内部Jacobian。新增标志仅报告这个分支。

## 改变与可比较性

新分支只在原中点不可采样且两种既有分辨率条件同时满足时启动。所有完整no-hit路径和其B/J运算顺序、正常target事件的bracket/refinement/交点表达式保持。原正常样本结果要求逐字回归相同；此前缩步耗尽的近壁数据新增可解析事件结果。

这是数值边界事件在舍入极限的定义，**不是连续ODE first-hit精确证书**，也不把H假装成实际RK2 target；其直接界是bracket空间跨度与时间跨度，不应误标为一般连续流场的全局ODE误差。原tracer继续要求resolved、continuous的velocity traces，不能证明任意振荡callback的首次碰壁。后续相积分还须对更新后的实际branch及旧内部/显式边界来源建立完整覆盖，不能仅凭本项事件运行成功推导通量。

## 算子验收

对四边/轴向边、极小/亚正规距离、非零墙坐标、两callback路径，比较独立常速事件位置和时间；所有query记录必须在闭域内，新分支必须无内部J。时间可舍入为零的事件需有明确夹具。多壁/切向/CFL/负半径/coverage拒绝保留。原rotation/vortex near-wall两个复现点转为south事件，原12个普通点保持7hit+5interior；旧I02数值输出以及r3前缀逐字比较。完整driver仍按原gate实际运行。
