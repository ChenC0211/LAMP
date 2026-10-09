# H08a2r16 / 实施前参考设计

2026-10-05。上轮 r15 为 progress，其 82 个冻结 payload、14 份交付及 1383 源码/输入已核验。本项遵循[任务书](H08a2r16_Resolved_Nonlinear_Motion_Reference_TODO.md)，production 0/test≤1；先做隔离参考审计，若原生产 gate 失败则保留反例，不向永久契约写入一个虚假 PASS。

直接以同一 regular preparation 与原 native values 得到 actual Q4/面修正系数；冻结 r11 evaluate 中也只调用同一 apply。导出两时端系数相等性、真实 support rectangles、完整系数、public sample 的值/J及原各步数轨迹，绑定输入、运行程序与实际库。参考侧以高精度代数展开每片完整四次 tensor；用分片 Bernstein 系数界及独立物理采样检查值/J，禁止用解析 sin 替代实际系数。

先核验实际 u_y−1 及其一阶导数的全片界。仅当能通过 Gronwall 与 Hessian 界把相应 point/J 扰动预算压到 1e−12 内，才把 y 当作精确积分时钟；否则回到完整二维系统。保留实际 u_x 的全部 x/y 依赖。对这个受控静态辅助场，积分 x'=u_x、P'=∂x u_x P，并用自治流身份 `J_xy=u_x(end)−P*u_x(start)` 求敏感度，其余 J 行由 y 时钟给出。这不是假定纯剪切；若 x 依赖非零，P 必须实际积分。

在每个真实 y support 区间内使用高精度 Taylor ODE，与区间 Taylor 导数余项界、先验轨迹包围及局部 Lipschitz 放大共同控制误差；跨片在实际几何 knot 处截断步长。确认整个轨迹及误差管留在选定 x slab，不能遗漏 x 事件。记录阶数/步长细化、十进制算术舍入、系数转换、真实场到辅助场的扰动，以及最终 point/J 总预算。原 native 数据的二进制浮点值按精确值读取，不能先截成低精度十进制。

用原 40/80/160、2560 步查询与参考比较，同时保留实际离散 map 的 FD。这两项分别检验连续流准确性和对自身 map 的导数一致性，不能互相替代。参考门槛仍 1e−12，生产 point/J 加预算仍 1e−11；若出现真实生产误差，不增加步数或容差，另拆生产修复。完整原 characteristic 分组仍需运行，不能由新 witness 的首个失败遮住其他覆盖。材料四硬约束、Pure FirstJet 和当前已接受生产路径均不变。

实施补充：原 cell 原点的幂基在区间求值时出现严重依赖膨胀，内部域第一步的高阶余项包围失去用途。改为在每步区间中点平移同一多项式，再做区间 Taylor；没有删系数、缩小实际输入或放宽精度。旧失败实现与日志保留为 `reference_initial*`。

二进制系数按精确有理数展开，65 位 Decimal 计算 Taylor 系数，52 位向下/向上舍入区间包围系数及高阶余项。用全片 Bernstein 界包围速度、一阶及二阶导数；每步先包围 x 和 P 的整个运动，再用第 p+1 个 Taylor 系数的区间绝对界乘 |dt|^(p+1)。另一次退化区间求值包围截断多项式求值舍入。平移/系数转换使用 1e−48 余量；独立 100 位计算核验转换误差远低于该余量。

误差传播利用三角变分结构。记旧误差 Ex/Ep、局部截断与舍入合计 bx/bp、全片界 A=|fx|、B=|fxx|、Pg=(|P|+Ep)exp(A|dt|)，取 `Ep_new=(Ep+bp+B*Pg*|dt|*(Ex+bx))*exp(A|dt|)`，`Ex_new=(Ex+bx)*exp(A|dt|)`。最终 Jxy 用端点 fx 界、Ex 和 Ep 包围；point 从无穷范数转为原 Euclidean 门槛时保守乘 2。

对于真实 uy 到 1 的扰动，令 eps=|uy−1|、epsJ=|grad uy|、L 为速度梯度行和界、H 为 Hessian 行和界，取 `pointPert=eps*T*exp(L*T)`、`jacPert=T*exp(2*L*T)*(epsJ+H*pointPert)`。额外核验每个被跨水平连接的完整值/梯度多项式差与法向速度下界，加入保守跳跃余量；本两域的精确有理连接差均为零。整个真实/参考轨迹误差管保留在导出的共同 x slab 内，未遗漏 x 事件。预算均针对实际系数定义的连续场，public double 求值的采样差另列，不将生产 trace 返回 J 当参考。

执行 12 阶/0.0025、12 阶/0.00125、16 阶/0.0025 三次参考，细化差还须落在不含 uy 扰动的 Taylor 预算内。两个可解析 affine ODE 检验 point 和全 J；实际把余项清零的参考负控被拒绝。将实际场错当理想剪切的 Jxx=1 负控也被拒绝。可信参考没有回写任何生产结果。最终裁决见[报告](H08a2r16_Resolved_Nonlinear_Motion_Reference_Report.md)。
