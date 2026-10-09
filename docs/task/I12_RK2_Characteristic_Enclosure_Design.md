# I12 / 实际 RK2 特征的区间覆盖

[I12a12](I12a12_Polynomial_Velocity_Interval_Source_Report.md) 已建立 owned polynomial velocity slab。本项沿[已授权的一般区间源设计](I12_General_Velocity_Interval_Source_Design.md)传播该数据，保持原 tracer；不替换为解析旋转或另一个保体积映射。

## 同一个离散表达

旧时刻 t0，实际到达时刻 t，d=t-t0。无缩短、无 hit 的原 RK2 为

`H = X - u_arrival(X) * .5 * d`，

`W = (u_arrival(H) + u_old(H)) * .5`，

`B = X - W*d`，

`J = I - d * D_W(H) * (I - .5*d*D_u_arrival(X))`。

端点函数来自同一 owned slab，但 W 必须按 Core 的实际 endpoint 平均生成；不能以直接在中间时间求值替代其浮点运算。原 double tracer 与新区间传播共用 midpoint、target 和 J 的代数表达；只改变求值数域，不改变 double 的括号、乘加顺序或物理方法。

若正的 subnormal 时长满足 `.5*d == 0`，原 tracer 会先进入 `elapsed==0` 的 arrival 端点分支。区间下界允许该情况时，W 与 D_W 额外包住 arrival 中点查询；不能一概假定每个正时长都会执行平均。零时长单点直接返回恒等映射。

包围同时覆盖空间 box 与到达时间段，包含闭区间 t=t0 的恒等映射。时间与空间依赖可导致范围膨胀；这是需要进一步细分的信号，不允许挪动实际回溯点或放宽物理域。

## 充分的分支条件

成功结果必须覆盖全部输入点与时刻：输入在 closed domain；起点不触发 outward backtrace；H 与 B 有有限域内包围；没有新的 wall contact；actual sampled old/arrival/midpoint velocities 满足原半 cell 限制；RZ 轴值/切向导数与不跨负 r 条件保持。old/arrival 查询所在的每个范围都单独核验。

CFL 用向外舍入的速度 L1 范数乘完整回溯 d 作保守充分条件，同时检查原 worst-component hypot 条件；不改变实际 tracer 的 `hypot(u)*d <= h/2`。域及尺度校验复用原 tracer 的构造检查。opaque active coverage 不能凭采样认证，明确拒绝。纯代数的已知非负时长、乘积符号和减法单调性可以收紧区间；这是范围证明，不是坐标裁剪。

如果包围跨越边界，即使某些实际样本可以推进，也不能据此声称整个 box 同属 interior 分支。本项拒绝这类区域；后续事件算子需要细分或明确的边界历史。RZ 轴线可在精确零径向速度下覆盖；包含轴邻域且无法证明非负 r 的一般范围仍须拒绝，不能把数值区间裁到 r=0。

## 来源与边界

Interface 结果保留原 `HermiteGalsPhaseFluxIntervalSource`、同一 old geometry、velocity 与 characteristic owner，并采用现有 mesh identity 校验。Core 仅负责中立范围运算；结果不建几何、不改 M、不提交候选。输入 geometry 要求覆盖 native domain，不能让部分 active cells 冒充完整旧标量支持。

现有 `NativeMeshIdentity` 比较精确、有序的网格几何/拓扑数据，不比较 C++ 对象地址。独立分配但这些数据完全相同的 mesh 可以使用；细化后的原对象或拓扑不同的另一对象都必须拒绝。

返回的是一个已核验的**单步 interior 离散分支**的 B/J 范围，不是连续 ODE 的轨迹管，也不是 boundary first-hit 证书。多项式值/导数的区间包含性、共用 RK2 代数及明确分支条件构成证明前提；密集测试只是检查实现。相符号、跨叶拓扑、事件完备、时空求积、共享库存有界性仍分别待实现。
