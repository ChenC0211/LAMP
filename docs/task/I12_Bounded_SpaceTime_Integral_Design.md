# I12 / 同源 gross face×time 积分与误差界

沿[时间区间通量授权](I12_Phase_Flux_Time_Integration_Design_Change.md)，复用 [a13](I12_RK2_Characteristic_Enclosure_Design.md) 与 [a14](I12_Characteristic_Phase_Coverage_Design.md)。积分对象是同一 owned coefficients 定义的重构速度及同一离散 RK2 表达拉回的旧含量；实际浮点速度/B_h 求值也由已有范围覆盖，不用连续 ODE 替换 B_h。

## 积分与量纲

全局 +x/+y 面法向速度为 u。分别计算 `T+ = integral max(u,0) dA dt`、`T- = integral max(-u,0) dA dt`、`L+ = integral rho max(u,0) dA dt` 与 `L- = integral rho max(-u,0) dA dt`。先平均速度再取正负部分会丢失反向通量，因此禁止。XY 的 dA 为 ds，RZ 为 2*pi*r ds，所有四项都包括时间积分。

重构速度对空间是 tensor cubic、对两个 endpoint 线性时间插值。沿 native face 的有符号积分可用 monomial moments 解析求值并向外舍入：规范化切向两端为 a,b 时，平均 q^j 为 `sum(a^(j-k)*b^k,k=0..j)/(j+1)`。该式避免小区间的原函数相减；RZ 水平面多乘 r 后使用到四次矩。时间只需同一 slab 的平均 theta。pi 也用包含真实常数的浮点区间。

当速度范围同号，解析净积分就是相应 gross 积分；变号 tile 使用正/负速度范围 × 测度给上界，并由 T+−T-=net 收紧。对几何/库存含量范围 [rho_lo,rho_hi]，L+/- 属于 [rho_lo,rho_hi] × T+/-，进一步与 0<=L<=T 相交。这些界不是高低阶 quadrature 差。

## 自适应覆盖

根 tile 是完整实际 native face × 请求时间区间。每个 tile 使用原 a13 RK2 范围及 a14 closed 叶/相/权重内核，所有来源都来自同一保留的旧 state。只对 typed branch-unresolved 情形尝试细分；非法输入、非有限算术或浮点环境错误不能被伪装成普通未确定。

候选切分为切向或时间的可表示内部二分。比较两组子范围的不确定量决定保留哪个方向，平局时使用参数宽度使两个维度都能继续收缩。所有最终叶构成原矩形分割，父结点用 strict outward 加法合并；旧父范围与子范围的交集仍是有效包围。选择与队列只影响效率，不提供额外相证明。

回溯 branch 尚未认证的 tile 可以暂时用 [0,1] 含量作为工作量估计，**不得作为成功结果导出**。无论误差界看起来多小，所有非零通量叶最终都必须覆盖原完整单步分支；缺少 inflow/history 或缩步覆盖时达到预算也不能成功。严格零通量的乘积独立于相；RZ 轴面仍检查精确零径向速度及切向导数。域/尺度/opaque coverage 的原校验不因零通量被绕开。

## 返回值与误差

Core 在已有 strict/no-LTO TU 中完成积分、区间合并/交集、误差宽度及代表值选择。Interface 只负责 native face、来源、旧相/库存内核与自适应树。每个 gross 范围必须非负，并保持液/总的上下端点有序。以半端点和在区间内选代表值；如发生 subnormal 舍入，选本区间端点保持包含性。单调的共同规则保证代表液量不超过总量。这是已证明误差区间内的近似积分选择，不修改或裁剪独立 M。

四个最终区间宽度的 outward 和 E 是代表值四项总绝对误差的上界，也控制净总/净液量。只有 E<=预设 absoluteVolumeTolerance 且未认证分支数为零才返回。该预算独立于守恒、M 界限、几何精度与曲率容限；不能当作 64eps 库存舍入修复的新增额度。有限 panel/depth 上限及不可表示切分导致明确失败，不自动提高预算。

结果私有构造、只读，保留 bounds/report、同一 interval 的单份 `HermiteGalsIntegratedFaceTransfer` 和可选 inventory owner。总/液都取自这次同源时空积分，实际 coarse side 以后仍须汇总这些共享子面 transfer。公共 a13/a14 查询契约不变；同一 TU 的积分内部在入口验证 mesh/来源后复用它们的范围内核，不反复构造虚假的子时间 geometry。

本项不承诺任意边界/缩步区间都可认证，也不预先证明原 M 装配必能接受一般积分结果。完整旋转/涡旋、时间收敛、动态 AMR 和低黏/无黏耦合仍按原 gate 验证。

## 实施补充与已知限制

积分专用 `encloseInteriorRk2TightTime` 复用原 a13 表达与全部 branch 守卫，把每次速度/导数查询与 endpoint 范围凸包加 outward 浮点带相交。实际 theta 由同一正 endpoint 差分构造，依浮点单调性位于 [0,1]；普通混合的 `1-theta`、两个乘积和一次加法，以及 `(a+b)*.5` 分支，其额外误差由 `8*epsilon*maxabs + 8*denorm_min` 包围。非有限中间算术仍由原区间拒绝。原 `encloseInteriorRk2` 入口数值不变，真实速度与 RK2 查询未替换。

切分方向先比较两组子区间误差；若差不超过父误差的 1/8，使用 departure box 的 3×3 旧标量值变化选择方向，仍相同则按参数宽度。这些样本只决定先做哪项工作，不能缩窄相范围、给出积分值或认证事件缺失。1/8 是调度参数，绝不进入成功条件或误差预算。所有证书依旧只来自 strict 区间与 closed Q4 相覆盖。尝试过的一层前瞻没有解决剩余反例，最终未保留该额外复杂度；过程日志归档。

I12a15 的预设均匀网格剪切精度未在 4096 panels 内达到，因此本项 FAIL。最大分块、深度、容差以及参考均未为通过测试而修改；完整新 contract 继续返回非零。后续须在新 TODO 下改进该反例的认证积分，不能把本报告的 56 个成功场景等同于原 I12 全矩阵。