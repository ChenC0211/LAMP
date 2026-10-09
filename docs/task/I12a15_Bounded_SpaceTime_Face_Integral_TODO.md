# I12a15 / 有误差上界的同源面—时间积分

目标：只解决 owned polynomial velocity 与同一 old Hermite/库存含量下，实际 native face × 时间区间的四个 gross 总/液体积分及绝对误差上界，产出可由原装配器消费的单份共享 transfer。结果须同时覆盖速度变号与采样间隐藏相事件。

允许修改：`PolynomialVelocitySlab.h/.cpp`（已有 strict/no-LTO Core TU）与 `HermiteGalsTemporalPhaseFlux.h/.cpp`，共 4 production 文件；新增 `tests/interface_bounded_spacetime_flux_contract.cpp` 与 `tests/CMakeLists.txt`，共 2 test-related 文件；本文、设计/API/报告、入口与证据。

禁止修改：原 RK2 代数、速度数据/求值、CFL/轴/边界/缩步规则、Q4 符号算术、库存权重、M 装配/舍入修复、bubble、几何发布、reflux、原 benchmark 参数、冻结文件及 gate。不把未知边界历史或 branch 当相=0，不放宽工作预算后声称旧配置通过。

数学约束：分别积分 `max(u_n,0)`、`max(-u_n,0)` 与它们乘同一 `rho^n(B_h)` 的四个非负量；XY 为面积、RZ 为真实 2*pi*r 体积。所有积分使用同一时间区间与 owned polynomial 数据。每个 tile 的上下界经 outward 算术累计，四项宽度之和不大于预先指定的绝对体积误差预算才成功；净总/净液量误差也由该预算控制。最终代表值须位于证实的区间内并满足 liquid+/- <= total+/-，不是对 M 作 clipping。

实现要求：先写设计。同一 polynomial 数据用解析 monomial moments 包围有符号总积分；速度符号明确的 tile 得 gross 总积分，符号未确定时以速度范围与净量关系包围。沿 a13/a14 原表达与 closed 叶覆盖获取 rho 范围；不能用粗细求积差替代证明。自适应二分 face tangent 或 time，父/子区间可作有证明的交集收紧。新增有类型的“RK2 branch 尚未认证”异常，保持原错误文本与规则；允许细分这些 box，其他输入/浮点错误直接失败。所有非零通量叶必须最终认证回溯分支，不能仅因其体积小而忽略未定义边界历史。零测度/严格零通量只在代数确证时可独立于相返回零；RZ 轴正则条件仍核验。

架构要求：私有构造只读积分结果保留真实 interval、唯一 old geometry/velocity/characteristic 与显式 inventory owner。原装配/reflux 继续消费 transfer，不在本项修改库存算法或 driver。a14 公开入口保留校验；其 closed-box 内核可在同一 TU 抽取复用，避免每个自适应 tile 重建 interval/全网格 identity。Core 算术与 Interface 数据职责保持。

验收：T0 构建、无新增 warning、范围/冻结/runtime 审计与旧相关回归输出不变。T1 纯相、正反平移、两端速度反号、平面跨越、时间节点间薄层、shear 与旋转/非线性速度，XY/RZ、uniform/AMR native face；独立解析或有明确误差包围的参考落在积分范围中，净量/代表误差满足报告；纯液/纯气与总量关系、半区间组成、重复调用、同源库存尾量、axis 零面积、细化预算收敛与有界失败均测试。外来/opaque/moved source、过期 mesh、缺少 inflow/branch 证明、非法配置与浮点环境明确拒绝，失败不改 accepted geometry/M。

失败规则：参考与积分界冲突、gross 子集关系失败、未认证 branch 返回成功、旧行为改变、work/depth/可表示性耗尽而未达预算，均 FAIL 并保存数据。不得将收敛慢或尚未支持边界包装成完整 I12 PASS。仅真正需外部输入且不能继续时报告 BLOCKED。

交付：源码、完整 T0/T1、精度/工作量记录、失败记录、API/报告与 archive。下一项再处理实际完整驱动的边界覆盖、库存可行性及原旋转/涡旋矩阵；本项不提前通过 Gate I。

实施补充（剪切反例定位后）：保留 a13 公开入口的数值输出不变；同一 Core RK2 区间表达增加面积分专用 tight-time 查询入口。它把原速度区间与 endpoint 值凸包加显式浮点混合误差带相交，消除 `(1-theta)` 与 `theta` 的重复依赖。真实速度数据、实际求值和 RK2 运算顺序不变。新增界必须覆盖 endpoint、普通混合、half-time 分支及渐进下溢；原入口回归逐字比较。自适应方向可用旧标量值变化作工作启发，任何采样均不参与相判定、积分值/界或成功条件。
执行状态：**FAIL**，不是 BLOCKED。原 57 个主积分场景中 56 成功，uniform shear 在既定 4096 panels / 40 depth、E<=5.859375e-6 下未通过，最终 E=6.676659722917087e-6。新测试与 CTest 保持非零；原 22 项回归通过。详见 [报告](I12a15_Bounded_SpaceTime_Face_Integral_Report.md)，未修改验收阈值或重分类为 PASS。