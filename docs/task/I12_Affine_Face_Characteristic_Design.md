# I12a15r2 / 原浮点 RK2 的相关仿射包围设计

见 [TODO](I12a15r2_Affine_Face_Characteristic_TODO.md)、[上一反例](I12a15r1_All_Partition_Certificates_Report.md)。本项仅新增 Core 面—时间 characteristic enclosure，不改变 a15 当前积分行为或 I12 gate。

## Pointwise 语义

原 native 端点 a,b 及时间 t0,t1 解释为精确实数，切向坐标 s=(s0+s1)/2+(s1-s0)*q_s/2，t=(t0+t1)/2+(t1-t0)*q_t/2，q∈[-1,1]²。固定法向坐标不参与独立参数。实际 arrival point/time 是这个连续参数域中的 double 输入，按原程序规则回溯。

内部 Form 为 C+S*q_s+T*q_t+R(q)，C/S/T/R 均 outward interval。与仅返回矩形不同，同一 q 在所有分量及原运算中一致；R 只保存不可保持为一次项的非线性/舍入不确定性。实际值的 pointwise 包含性由运算归纳证明，不对舍入函数求导。

## 运算归纳

加法将同类系数相加；额外加入原 RN 加法舍入带。乘法保留 CA*CB、CA*SB+CB*SA、CA*TB+CB*TA；余项包含 CA*RB+CB*RA，以及 (SA*q_s+TA*q_t+RA)*(SB*q_s+TB*q_t+RB)。后一个乘积用两个变动部分各自范围作 outward 乘法。每项运算覆盖的是原真值，区间系数计算的舍入不确定性本身由 interval 表示。

为覆盖原求值顺序中的浮点操作，另在 R 中加入 eps*maxabs(原操作值范围)+denorm_min 的对称带；这是 round-to-nearest 的保守误差界，含渐进下溢。除以正的 double 常量逐系数向外除，再补同样 RN 带。精确零、加零、乘零/一和精确反号按代数保留，避免无依据地破坏轴上的精确零。任何非有限区间拒绝；不以人为大数代替溢出。

原坐标归一化、Horner component 和 RK2 shared helper 都直接消费 Form。时间 theta 使用原程序实际 double 分母 end-old；实际 theta=0/1 时普通混合等于 endpoint。theta=.5 特殊 `(a+b)*.5` 与理想半和的误差用 4*eps*maxabs(a,b)+4*denorm_min 包围，加入一般混合的 R，不把整个时间依赖抹成常数 hull。

原 dt/2 舍入为 0 时，tracer 的 arrival-only 分支与一般 mean 分支不同。该情况只有 dt<=denorm_min；使用已覆盖的 arrival/mean 速度范围和 arrival 坐标范围，给最终 B 增加 8*eps*maxabs(X)+8*denorm_min*(1+maxabs(U_new)+maxabs(U_mean)) 的对称差异带。它覆盖两个乘积/减法的差异及 RN 误差；原 tight-time witness 仍核验原规则。这是小分支舍入带，不将完整时间区域替换为 arrival-only map。

## 输出与消费

最终在 C/S/T 的区间各选一个有限代表 c/s/t，把 C-c 及 (S-s)*[-1,1]、(T-t)*[-1,1] 向外累计到 R。公开结果包含实际 face/time 元数据、原 branch witness、两个 departure forms；物理点查询由端点差向外反算参数，已知点在原域内时与[-1,1]相交，再评价 affine interval并与 witness 的 box相交。

下一项可以把旧 Q4 的权威多项式施加于同一 forms，得到 phi(B) 的一次参数项与严格剩余带；随后对一次项+带所认证的纯液/纯气子区域计量。多 native donor leaf、库存 wl/wg、RZ 真实测度和 gross 正反向量仍需独立处理。本项未实现这些消费者，也不宣称原 shear 已通过。