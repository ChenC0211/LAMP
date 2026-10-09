# I12a16r3 / 原速度表达式的符号与精确零认证

目标：为原tensor-cubic速度及导数表达式提供保留可证明符号/精确零的区间包围，解决一般outward界在对称中心和壁面方向判定中人为跨零的问题。保留真实回溯，不用有限采样判零。

允许修改：Core PolynomialVelocitySlab.h/.cpp（2production）；原polynomial_velocity_contract与bounded_spacetime_flux_contract（2tests）。禁止修改tracer、旧速度evaluate、相历史、库存装配/reflux、I10/I11、冻结、Hydro/Electricity。

数学约束：坐标normalize中x>=origin可认证RN((x-origin)/scale)>=0，x<=origin同理，退化x=origin严格为零；正scale不变号。复用原Number的同号加法和乘法符号约束，用原Horner操作顺序传播全部value/first derivative。时间theta和1-theta都在[0,1]，endpoint/mean/general mix原分支全部覆盖。仅用代数及单调性收紧区间，不以阈值把近零改成零。最终与原enclose界相交，公开原enclose/evaluate保持；现有tight-time消费显式改用新界。

实现要求：一个新可检验Core查询，沿用原输入/FP/非有限拒绝；全tensor-cubic、非单位basis、时间反转、单边/跨零/退化box都覆盖。认证速度、导数和原实际RK2轨迹；核心验收必须含导致原driver面4失败的rotation/vortex及其邻近反例。后续真的回溯到物理边界仍须拒绝未认证分支。

验收：T0 2+2范围、39targets/26CTest、81冻结/3PDF/r2archive保持。T1原velocity完整输出保持；新界包含原实际query及独立解析值，且是旧界子集；精确零/符号/反转/nextafter邻点/亚正规/rounding guards，原57积分references和budget不变。更紧bounds可以改变认证统计/E，但不得改变旧发布状态或用扩大预算通过。T2同四个driver first steps及原full命令实际运行，证明面4缺口是否解开，记录后续真实失败或完整通过；a16仍要求4完整steps，不设WILL_FAIL。

失败规则：新界逃逸或来源/FP/regression不满足记FAIL；真实事件尚缺能力保持原driver/I12 FAIL，不以较易的局部PASS代替完整gate，不将内部实现缺口标记goal BLOCKED。

交付：证明/API说明、tests、原driver与回归结果、report及archive。

执行结果：**PASS（本项）**。120速度区间/20,600样本、16完整面/1,296实际trace通过；23旧完整输出及velocity/integral整个旧前缀保持。原四场景越过face4，在face5的真实边界/内部混合分支拒绝。另直接诊断确认两个近壁点触发原tracer的128次缩步上限。39targets/26CTest审计通过；driver仍真实FAIL，a16/I12仍FAIL。见 [报告](I12a16r3_Sign_Preserving_Velocity_Bounds_Report.md)。
