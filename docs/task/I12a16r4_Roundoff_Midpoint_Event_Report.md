# I12a16r4 / 近壁中点事件报告

**PASS（本项事件可解析性）。原a16驱动、I12 / Gate I仍FAIL。** 见 [TODO](I12a16r4_Roundoff_Midpoint_Event_TODO.md) 与实施前写明的 [设计](I12_Roundoff_Midpoint_Event_Design.md)。独立库存和几何预算继续分离；本项没有恢复逐cut cell机器精度几何体积匹配。

## 修复范围

仅1production+2tests。原不可采样midpoint已在既有时间/坐标舍入分辨率内时，用原inside/outside bracket交点规则结束边界事件，再检查activeSegment、唯一侧及边界真实入流速度。新增`roundoffMidpointEvent`标志用于识别该分支；命中不返回内部Jacobian。原128缩步、128子步和96事件细分上限、两个128epsilon事件容限都没有增大。

共享原交点helper保持原target事件表达式顺序。新midpoint没有被采样、clip或作为full-step target使用。hitTime可以舍入到当前elapsed，但不会接受零时间内部substep。完整内部RK2及J、速度源、相历史、库存装配/reflux和I10/I11均未修改。这是数值事件的舍入分辨率规则，直接界为bracket跨度；不是任意callback连续ODE first-hit精确证书。

## 原轨迹验证

66个新事件场景覆盖XY四边、RZ南北、零/非零壁坐标、nextafter/1e-16/1e-100/1e-300/最小正double距离、两种速度和时间下溢。与独立常速事件参考比较，位置/时间在原分辨率内；primal与jet两条路径的事件、采样顺序完全一致。没有域外query，新事件最多3次callback；12个事件时间舍入到零，均只返回boundary hit且无内部J。

另验证已接受一个子步后才靠壁的事件，acceptedSubsteps=1，时间不倒退。8项多壁/切向/非入流/coverage/axis负半径/CFL/非有限速度/非有限导数拒绝通过。普通target hit、起点已在壁上的hit和dt=0 identity不走新分支。原I02完整输出（含B/J差分和二阶时间收敛）逐字不变。

原rotation/vortex面5的两个y=1e-100复现点现在都成功返回south事件，物理域内查询、时间和无J均核验。原12普通点仍是7个south hit+5个interior。原bounded integral全部输出仅最后诊断行从“2个缩步拒绝”改为“2个near-wall events”，其余内容逐字不变。

## 完整驱动仍未通过

同四个driver smoke/full均实际执行。它们仍在nativeface5、半步的内部／边界混合来源认证上拒绝：rotation E=3.3260603868273389e-6、pending3780；vortex E=9.8832032391906215e-7、pending3772，均4096panels，预算仍1.5625e-8。

本项未修改相积分证书，不能以actual tracer近壁事件已修复推导完整混合通量。driver的完整输出与r3逐字一致；accepted时间0、旧identity/M保持，没有执行到半步库存装配/闭合。smoke5.201837s、full5.2763216s，均exit1；full仍请求原30/400步。

## T0与归档

首次3targets构建83.158777s；全部39targets重建259.9778984s，均exit0。共享头触发旧模块重建，出现既有11条`ImmersedGeometry_impl.h` C4267；与a13历史日志逐条相同，该源文件SHA仍638A623F0C0AF5ED30920403C14D9F3C9276055B75CDD9BCCC9B43E6A26F70EA，**新增warning为0**，未抑制或顺手修改旧代码。

trace contract首次.0886968s、最终.13512s；integral contract首次.7593105s、最终.7179811s，均exit0。26CTest165.4345949s：25PASS+原driver FAIL，无WILL_FAIL。24项旧完整输出（包括driver的真实FAIL）、原trace完整前缀、上述integral除诊断外的所有输出、旧crossing probe exit1/输出均保持。

架构12targets、383TU、908files、3118includes、5exact legacy edges不变；81frozen、3PDF和r3 archive身份保持。首次审计PASS。归档交付3源文件、TODO/设计/本文/入口、39runtime及原DLL，逐项核对payload SHA。

证据：[验收](../../runs/baselines/I12a16r4_20261002/verification.json)、[原轨迹](../../runs/baselines/I12a16r4_20261002/trace_verified.log)、[原polynomial复现](../../runs/baselines/I12a16r4_20261002/integral_verified.log)、[驱动](../../runs/baselines/I12a16r4_20261002/driver_verification.json)、[回归](../../runs/baselines/I12a16r4_20261002/regression_comparison.json)、[警告核验](../../runs/baselines/I12a16r4_20261002/warning_audit.json)、[范围](../../runs/baselines/I12a16r4_20261002/source_audit.json)。

## 下一项

建立更新后原actual tracer的内部／边界混合覆盖证书，再接入原face5及其余native faces。需同时约束合法速度/可达区域、实际事件/子步分支及有限终止；内部来自同一旧snapshot的wl/wg，边界来自同一显式history。不能将“液滴离边界较远”或有限点trace成功当作完整来源证明，也不能按小E丢弃pending。之后继续同四个first steps、原full及时间/空间/AMR/曲率gate，未进入Hydro。
