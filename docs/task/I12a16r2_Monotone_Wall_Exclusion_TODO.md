# I12a16r2 / 原RK2坐标单调性排除新壁面接触

目标：只补一个数学事实——原 target=P−W*d、d>=0 的浮点坐标单调性可以排除方向相反的新壁面接触，覆盖含壁面起点与旧时间的闭参数块；不改变原离散轨迹。

允许修改：src/Core/Geometry/PolynomialVelocitySlab.cpp（1 production），原 tests/interface_bounded_spacetime_flux_contract.cpp（1 test）。禁止修改原tracer、velocity求值、历史来源、积分算法、M装配/reflux、I10/I11、冻结与其他模块。

数学约束：已证明midpoint/target在域内、初始无hit、CFL/axis/coverage均保持。若 W.upper<=0，则RN(P−RN(W*d))>=P，内部点不可能新触及lower wall；W.lower>=0同理排除upper wall。负零/乘法下溢不破坏弱单调性。未证明的方向仍用原严格departure/arrival空间界判据。不能丢掉正测度pending或按小误差豁免事件。

实现要求：只增强noNewWall逻辑，不改任何bounds计算/返回数值；原已通过query仍相同。补四边XY、RZ南北边、常/变空间速度，含完整旧时间端点和壁面起点的闭块；对照原实际trace和全部B/J界。真实初始/后续壁面事件仍拒绝。

验收：T0 1+1范围、39targets/26CTest、81冻结/3PDF/r1archive保持；除扩展contract及原4场景driver外旧完整输出保持。T1单调闭块证书与原轨迹、零时间identity、反向/真实壁面事件拒绝，原57积分和r1全前缀保持。T2同四个driver first steps全部实际运行，要求原内部face1缺口被解开并记录随后真实结果；原a16仍4完整步全过才PASS，不设WILL_FAIL或降预算。

失败规则：本项证明或回归不满足记FAIL；driver后续缺口保持真实FAIL，未完成原I12 gate不得进入Hydro，无外部阻断不标记goal BLOCKED。

交付：补充证明、tests、原driver结果、报告/审计与归档。

执行结果：**PASS（本项）**。18闭块、4,050原trace、8真实事件拒绝通过；同四个driver均完成原face1，推进至face4的initial-hit认证缺口。39targets构建及审计通过，26CTest为25PASS+原driver FAIL，24旧完整输出及整个原contract前缀保持。见 [报告](I12a16r2_Monotone_Wall_Exclusion_Report.md)。a16/I12仍FAIL，goal继续。
