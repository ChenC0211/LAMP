# I12a16r3 / 原速度符号包围报告

**PASS（本项）。原a16驱动、I12 / Gate I仍FAIL。** 见 [TODO](I12a16r3_Sign_Preserving_Velocity_Bounds_TODO.md) 与 [说明](../api/I12a16r3_Sign_Preserving_Velocity_Bounds.md)。本项保持独立M和几何精度预算分离，未恢复逐cut cell机器精度体积匹配。

## 实现

2production+2tests。新增Core `encloseSignPreserving` 对同一slab、同一归一化/Horner/导数递推/时间混合表达式包围。利用输入相对origin的位置认证归一化坐标符号，退化坐标等于origin严格为零；复用Number同号加法、乘法及正数除法的符号单调性；theta和1-theta在[0,1]，包含原endpoint/mean/general mix全部分支。结果与原enclose相交。没有近零阈值，原evaluate、原enclose及actual tracer保持。

现有tight-time查询显式消费新界后，继续原混合舍入误差包络、CFL/axis/coverage、initial-hit、midpoint/target和r2壁面单调性认证。新API本身不宣称characteristic或phase有效。未修改历史来源、库存装配/reflux、I10/I11、积分误差预算和64eps修复限。

## 验证与真实驱动

120个速度区间、20,600次值/四导数检查通过：包含原8类流场、64个全tensor基函数/四象限/非单位且非二进制基底组合；12个对称中心精确零、70个符号/相邻可表示非零/亚正规检查、时间反转及3项输入/FP拒绝检查。新界包含实际求值与原流场的独立解析值，且是原界子集。原velocity测试的整个输出前缀逐字相同。

原导致face4失败的结构扩为16个XY四个壁面端点/rotation-vortex/uniform-AMR完整面，1,296条原actual tracer验证全部no-hit及B界，保留内部几何气相gasWeight=.25；16个面均1panel通过，未借用边界history。原57主积分、r1/r2测试的整个输出前缀逐字相同，未改变reference、budget或主积分认证统计。

原四个driver first steps均已完成face0–4的积分，随后在 **内部nativeface5（negative4、positive5）** 拒绝。该面x=.625、y∈[0,.125]存在真实南壁回溯。

| 原场景 | 半步末端 | E | panels / pending | 最后拒绝 |
| --- | --- | --- | --- | --- |
| rotation uniform/AMR | .005 | 3.3260603868273389e-6 | 4096 / 3780 | target domain unresolved |
| vortex uniform/AMR | .0005 | 9.8832032391906215e-7 | 4096 / 3772 | cannot exclude initial boundary hit |

面预算仍1.5625e-8、depth40。smoke5.2239243s、full5.1502925s，均exit1、completed0/planned4；full仍请求原rotation30步、vortex400步。实际mid/end预测与首面history证书通过；accepted时间0及旧identity/M保持，未完成半步装配/闭合。

为下一步提供直接证据，测试在相同owned endpoint coefficients、相同half arrival和面5上运行原tracer：12个点中7个确实命中south，5个为内部轨迹。另两个y=1e-100内部点均被原tracer的128次midpoint shortening上限拒绝，错误为`midpoint remained outside after bounded substep shortening`。这是**尚未修复的真实运行限制**，作为DIAGNOSTIC明确记录；本项PASS不表示这个问题已解决。不能再把面5全部视作no-hit，也不能将其pending区域直接填成气相。

## T0 / 归档

首次3targets构建30.6278258s成功；全部39targets重建69.0752852s、exit0、无warning。速度contract首次.137389s、重建后.0706067s；积分contract首次.773882s，补充实际回溯诊断后.7913955s，均exit0。

26CTest163.8015648s：25PASS+原driver FAIL，无WILL_FAIL。23旧完整输出、velocity/integral两项的完整旧前缀、旧crossing probe exit1及输出保持。架构12targets、383TU、908files、3118includes、5exact legacy edges不变。81frozen、3PDF、r2 archive SHA保持，首次审计PASS。归档交付4源文件、TODO/说明/本文/入口、39runtime及原DLL，逐项核对payload SHA。

证据：[验收](../../runs/baselines/I12a16r3_20261002/verification.json)、[真实驱动](../../runs/baselines/I12a16r3_20261002/driver_verification.json)、[原tracer诊断](../../runs/baselines/I12a16r3_20261002/integral_verified.log)、[回归](../../runs/baselines/I12a16r3_20261002/regression_comparison.json)、[范围](../../runs/baselines/I12a16r3_20261002/source_audit.json)。

## 下一项

按原程序路线处理近壁事件与混合来源。先形成明确设计：原tracer在合法CFL下的近壁内部点也可能缩步耗尽，不能仅增加积分panels、越过实际trace失败或用解析内部重初始化规避。解决其可解析性后，认证内部到旧snapshot及边界到同一history的完整覆盖；内部仍使用原wl/wg，边界使用显式物理相数据。再执行同四个first steps及原全程/收敛gate。I12未过，未进入Hydro。
