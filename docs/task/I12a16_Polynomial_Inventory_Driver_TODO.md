# I12a16 / 原半步整步的 polynomial 库存驱动接线

目标：将已通过 a15r4 的原区间相积分接入原 benchmark 的 actual first-jet prediction、库存装配/reflux、几何预算闭合和唯一发布链。对原 XY rotation/vortex、uniform/AMR 的全部 native faces 实际运行；不得把只在内部面通过的 primitive 当作完整驱动。

数学约束：半步/整步使用同一个旧 accepted snapshot、owned polynomial velocity 与原 RK2 characteristic source；full step 的线性时间重构由原 Case 的两端速度给出，midpoint view 必须是其真实子阶段而非新的 velocity owner。M 只由共享 integrated faces 装配/reflux更新，first jets/center来自原 actual tracer，I10/I11沿原 mesh-scaled geometry/Hessian预算发布。M 不重置到解析值、不clip，不提高64eps修复限制。

允许修改：只改原 tests/interface_hermite_gals_benchmark.cpp 与 tests/CMakeLists.txt（0 production+2 test-related），在原 initial/predict/finish/runCaseUsing 内显式选择新路线，默认旧路线数值与输出不变。保留旧 mid-point flux A/B；不以旧 producer 作失败 fallback。

实现要求：给原 motion=2/3 生成相同 tensor-cubic endpoint coefficients；使用 a12 polynomialArrivalView 连接实际 Midpoint 身份。原整个mesh每个 native face都调用已通过的 polynomial integral，汇集原 transfer→inventory→reflux；不可跳过边界face或从理论总质量代填face。误差预算预先固定为每面关联最小cell容量*1e-6（积分absolute budget），最大4096 panels/depth40；与原 geometry τ=max(1e-10,(hmin/L)^3) 和 h² Hessian预算分离。

验收：T0 本task仅2测试文件，benchmark重建且其余38 runtime hash保持r4，81 frozen/3PDF/r4 archive不变；原25CTest完整数值输出逐字一致。T1 verify endpoint coefficients/actualsubstage owner，默认旧路线不变；T2新增 ordinary CTest的4个完整原 first steps：rotation n8 dt.01，vortex n8 dt.001，均uniform/AMR，cx=.65 cy=.53 radius=.13，speed2/40 period.4。每个成功步骤须完成两阶段M/firstjet/闭合/发布且M严格有界、ledger/reflux<=1e-10；4个全过才本项PASS。可额外同命令full运行原30/400步，不能缩成只跑安全内部面。

失败规则：任一真实步骤在boundary/source/积分/库存/闭合被拒绝则保留新CTest真实FAIL，不用WILL_FAIL重分类、不放宽budget或强行补气相。记录native face、phase、stage、旧accepted时间/未修改M证据，下一项按实际缺口写独立设计与TODO。真正外部依赖才BLOCKED；本任务自身FAIL不是线程goal BLOCKED。

交付：新driver选择/commands、完整失败或成功日志、旧回归与runtime身份核验、report/API入口及归档。通过本first-step接线仍不等于I12全时间/空间/AMR/曲率gate通过。

执行结果：**FAIL（本项）**。四场景的实际end/midpoint预测同源通过，smoke/full均在首个物理入流面0被原no-hit分支认证拒绝，4096 panels全部pending，accepted时间仍0且M未变。原25项CTest完整输出不变，新ordinary CTest真实FAIL；0production+2tests。见 [报告](I12a16_Polynomial_Inventory_Driver_Report.md)，下一项补显式polynomial边界相历史及完整可达覆盖。I12 / Gate I仍FAIL，线程goal继续。
