# Hermite–GALS V2 characteristic API（I02）

对应 [I02 TODO](../task/I02_Characteristic_Map_Jacobian_TODO.md)、[报告](../task/I02_Characteristic_Map_Jacobian_Report.md)。声明位于 [HermiteGalsTrace.h](../../include/Modules/Interface/HermiteGalsTrace.h)，命名空间为 `Quantum::Interface`。

I02 建立返回实际回溯点 `B_h(x)` 与同公式 Jacobian `J_Bh(x)` 的 V2 接口。**它只推进特征线，不推进 skeleton 或相体积。** `conservativeHermiteGals` 尚未注册为可运行方法；I03–I12 的验收继续按路线图执行。

## 一次调用的语义

```cpp
auto result = traceHermiteGalsWithJacobian(domain, arrivalPoint,
                                          previousVelocity, arrivalVelocity);
auto departure = result.departurePoint();
auto jacobian = result.interiorJacobian(); // boundary hit 时明确抛错
```

`HermiteGalsVelocityView` 持有非零 geometry identity、I01 的 `HermiteGalsTimeLayer` 和一个 paired query。query 一次返回 `VelocityJet { value, derivative }`，其中 derivative 为同一空间速度插值的

\[
D u=\begin{pmatrix}\partial_x u_x&\partial_y u_x\\
\partial_x u_y&\partial_y u_y\end{pmatrix}.
\]

不需要保存另一套全局速度梯度场，也不需要速度 Hessian。provider 必须保证 value/derivative 来自同一插值分支、可重复且在调用期间不变。view 拥有 callable；若 callable 借用了其它对象，provider 仍负责其寿命与并行只读安全。

旧速度和 arrival 速度可以属于不同几何 revision。结果保留两个 identity、两个时间层和 arrival point；tracer 不修改接受态、候选态或几何，也不控制主时钟。后续 I03/App 必须把这些 view 与自己的几何候选正确绑定，I02 不凭 callback 推断其物理来源。

## 时间区间

实际起点为 `previousVelocity.timeLayer().time()`，终点为 `arrivalVelocity.timeLayer().time()`。要求旧速度时间精确等于 arrival 时间层的 `begin`，再由两者之差确定回溯时长；没有独立的裸 `dt` 参数与它们竞争。

- `End` 的 arrival velocity 是层的 `end` 时刻速度。
- `Midpoint` 的 arrival velocity 必须是 `begin + (end-begin)/2` 时刻速度；这一调用只回溯半个区间。不得把完整步末速度标为 midpoint。
- 零区间只允许 `Initial`/restart 语义；不采样速度，返回精确点值及单位 Jacobian。
- 若 advancing/midpoint 时间在浮点数中已无法与 begin 区分，显式拒绝，不伪装成零步。

正常完成后 departure time 直接保留旧时间，避免通过 `newTime-dt` 再做一次舍入恢复。物理边界命中时返回实际 `arrivalTime-elapsed`。

## 离散公式与复用边界

固定使用显式 midpoint RK2。对时长 `delta` 的一个成功子步：

\[
x_m=x-\frac{\delta}{2}u(x,t),\qquad
B=x-\delta u(x_m,t-\delta/2),
\]

\[
J=I-\delta D u(x_m,t-\delta/2)
       \left[I-\frac{\delta}{2}D u(x,t)\right].
\]

速度与空间导数在旧/arrival 端点之间用相同的线性时间插值权重。子步的时间使用真实累计 elapsed，不能每次套用原完整区间的平均速度。多个子步的 Jacobian 按轨迹顺序左乘组合：`J_total <- J_local * J_total`。

实现直接调用通用 Core 的 [LevelSetBoundaryTrace::traceWithJacobian](../../include/Core/Geometry/LevelSetBoundaryTrace.h)。Core 的 value-only 与 tangent 版本使用同一模板控制流，因此没有新增第二套 RK、边界或时间插值 stencil。本轮没有修改 Core，也没有调用旧 `ActiveQuarticHermiteGals` predictor、parametricGals 或 CompatibleMoment。

Jacobian 的行是 departure 分量，列是 arrival 坐标；后续 first-jet pullback 应使用 `J^T*grad(phi_old)`。它是当前插值/子步选择分支内的离散导数。跨非光滑速度接缝或子步分支切换的导数不自动存在；provider/调用方必须明确这种位置的求值约定。

## 边界及有效域

保留 Core 已有的有限矩形域、active-segment、axis parity 和半单元位移检查：

- 每个实际采样速度均满足 `dt*|u| <= min(cellSize.x,cellSize.y)/2`，这里 dt 是本次完整 trace interval。
- domain 及 cellSize 必须有限、有效；尺度应由实际网格/活动域提供。I02 不要求速度来自 uniform grid，也不创建第二张网格。
- 外部物理域不采样速度；跨越 inactive coverage 明确报错；RZ 不裁剪或反射负半径。
- r=0 要求径向速度与其切向导数精确为零。轴线上即使点值不动，横向 Jacobian 仍可非平凡。

这是受上述分辨率条件约束的离散 characteristic，并非任意振荡速度连续 ODE 的 certified first-hit oracle；没有从半单元检查推导完整输运稳定性。

物理 wall/inflow 命中时，`boundaryHit()`、`boundarySide()`、`elapsed()` 和命中时间有效，`jacobian()` 为空。命中时间与 prescribed scalar history 的导数尚未定义，所以 `interiorJacobian()` 拒绝提供错误的内部 pullback。后续边界处理必须显式提供相应物理历史；不能将此状态改成单位矩阵或自动切换其它 integrator。

## 验证范围

[interface_hermite_gals_trace_contract.cpp](../tests/interface_hermite_gals_trace_contract.cpp) 验证：

- affine shear 的解析位置/Jacobian、明确的 midpoint 物理时长、零时间和静止速度；
- 3 个 dt、12 个位置、6 个扰动尺度上的四分量中央差分收敛；
- 触发两个成功子步的非交换速度梯度，独立 value-only Core map 与 tangent 路径点值及采样序列完全相同；
- time-dependent rotation 在固定总时间下的均匀和交替 dt 二阶收敛，分别检查位置和全部 Jacobian 分量；
- 四面边界事件、非平凡轴线 tangent、非法数据/时间/identity/domain/axis/CFL 的拒绝。

有限差分只在测试中使用。rotation 的角速度线性依赖时间，端点时间插值对该解析场精确，因此该项单独验证 RK2 的全局时间阶。任意 Hydro 速度数据的时间/空间误差、AMR 接缝和耦合误差不在此结论内。原始误差与阶数见 [I02 证据目录](../../runs/baselines/I02_20261001/)。
