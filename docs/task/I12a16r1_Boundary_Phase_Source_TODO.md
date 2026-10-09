# I12a16r1 / 显式边界相历史与立即入流分支

目标：解决一个来源事实——owned polynomial 速度下，原回溯器在到达点立即命中物理入流边界时，积分可从显式、同源的边界相历史取得内容。a16 的完整 first-step gate 保持，不把本项局部分支通过当成完整驱动通过。

允许修改：Core PolynomialVelocitySlab.h/.cpp 与 Interface HermiteGalsTemporalPhaseFlux.h/.cpp（4 production）；现有 bounded_spacetime_flux_contract 与 hermite_gals_benchmark（2 tests）。禁止修改 tracer、I10/I11、库存装配/reflux、冻结代码、Hydro/Electricity 及原误差预算。

数学约束：每边的整段物理时间历史显式声明 Gas/Liquid，不从 phi、M 或采样推断。历史保留 unbound interval 的旧状态、速度、characteristic、时间范围；interval/transfer/append/assembly 都保留该所有权。闭面上每个正持续时间的原 tracer 必须立即且唯一命中同一边；须验证原 domain、axis、CFL、无 opaque coverage，其他边不得同时命中。t=旧时间的 dt=0 原分支不改变，它是时间积分的零测度端点；必须在说明中明确。不得据此放松其余 no-hit/事件分支。

实现要求：独立于旧 constant-translation inflow capability 的常相四边历史；不授予 polynomial 恒速 capability。严格区间包围真实速度，不能拿有限采样证明分支。没有唯一立即命中证据时仍使用原 no-hit 证书；未认证分支仍 pending/拒绝，不用气相 history 越过内部旧 wl/wg。新 driver 显式提供其既有解析边界条件，全程气相声明需由场景的液滴支持范围推导，不能靠离散检查。

验收：T0 范围4+2、重建现有39 targets、26 CTest、冻结81文件/3PDF/前项archive不变；除扩展的contract和原4场景driver诊断外，旧输出保持。T1 四边/正负/XY/RZ适用边、time-subinterval、复制和独立所有权、append/assembly拒绝、缺历史/切向/多壁/CFL/axis/opaque拒绝；和原 actual tracer 比较，但采样仅作验证。T2 原四个驱动first steps全部实际运行，不跳native faces；要求此前首面0由证书完成，并记录随后完整步骤真实结果。原a16 ordinary CTest依旧要求4步全过，不设WILL_FAIL。若后续分支仍拒绝，则a16/I12保持FAIL，记录真实下一缺口。

失败规则：本项source或证书/回归失败记FAIL；下一边界事件能力未实现不假装已通过。无外部阻断则继续goal，不以内部数学缺口声明BLOCKED。

交付：代码、边界来源与零测度端点说明、原场景执行记录、回归/冻结/hash审计、报告和archive。

执行结果：**PASS（本项）**。28个边界积分、1,134条原trace、48项拒绝/来源检查通过，原四场景首面均1panel认证。完整驱动下一内部面1仍pending1/depth40拒绝，a16/I12保持FAIL。26CTest为25PASS+原driver真实FAIL；24旧完整输出及原积分全部前缀不变。见 [报告](I12a16r1_Boundary_Phase_Source_Report.md)。
