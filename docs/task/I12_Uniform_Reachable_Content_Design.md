# I12 / 常量可达来源的相通量求值规则

## 动机与语义

r4已修复原tracer的近壁缩步事件；原face5仍有混合内部和boundary来源。对任意非恒定含液量，仍需完整来源映射证明。对所有可达旧库存及相关边界历史都严格等于同一常量c的区域，输运后的含液量就是c，故L±=cT±；求出分支、first-hit时间或J并不提供更多信息。

这是显式新增的**常量输运量代数求值分支**，不是“可达box证明了完整tracer”或失败后的gas fallback。它在求map之前执行；同一输入所有成功actual traces都与之相容。当原map因grazing、多壁或有限迭代不能返回时，只延拓唯一的常量输运量，不延拓/伪造B或J。I-R4的first-jet传输仍完整使用原同一B/J；代码不改变tracer或其拒绝条件。新分支不支持非恒定相内容，不能解决一般boundary相输运。

## 同源可达包围

使用原owned tensor-cubic slab在整个闭物理域及old→arrival上的严格sign-preserving包围，沿用r3的endpoint-hull加舍入界时间收紧，再增加原tracer二次时间混合的8epsilon*maxabs+8denorm_min舍入界，精确零分量保持零。所有速度与一阶导数有有限界，原全步half-cell CFL的充分L1界及实际hypot检查均必须通过。不接受opaque activeSegment。RZ本分支仅支持全域严格零径向速度及径向速度切向导数；其他RZ仍走既有认证。

记D为最大实际double duration，U_i为混合后速度绝对界，A_i为域坐标绝对界，e=double epsilon（比unit roundoff更保守），eta=最小正double。原tracer至多N=128个子步。每个非末步时间累加，以及末步remaining减法的误差不超过2(eD+eta)，故所有实际duration之和被D+2(N+2)(eD+eta)包围。这也包含最后一个event trial的duration；event bisection的trial从同一个起点重算，不累加运动距离。

每个target或midpoint的坐标增量舍入误差不超过4e(A_i+U_iD)+4eta(1+D)，包含原speed*.5*duration的半乘下溢。至多N+1个这样的增量加上最终bracket插值（fraction在[0,1]）的16e(A_i+U_iD)+16eta(1+D)，可用以下向外计算的半径统一包围：

`R_i = U_i [D + 2(N+2)(eD+eta)] + [4(N+1)+16] [e(A_i+U_iD)+eta(1+D)]`。

另外保守约束参与中间加减的绝对界不超过max_double/16，避免证明中的溢出。严格零速度分量的所有坐标操作保持原double坐标，因此R_i=0。源box为arrival box±R与闭物理域的交；这里只限定来源集合，不clip轨迹。若非零分量的box接触某壁，保守纳入该壁，不依赖采样或壁上速度符号来缩小集合。这包含普通event、r4event及成功no-hit；也包含有限速连续路径，因而支持上述常量延拓，不提供完整map的终止或Jacobian证书。

时间范围使用同一old→arrival历史。新能力限old=0或0<old<=arrival<=2old，使actual duration差值及arrival-duration恢复old满足精确减法关系；物理hit time不会越过历史起点。其他时间范围不获得新能力，既有算法保留。

## 相内容与执行次序

全闭source box使用既有native leaf覆盖与authoritative Q4严格符号，带入同一独立库存的wl/wg。任何leaf交界或零点的不确定性都保留在范围中。所有可能物理侧的同源显式history加入同一含液量范围；history必须覆盖old起点，不能仅覆盖当前时间panel。只有`minimum==maximum`时可提取常数，不能增加几何预算或舍入修复预算。

原zero-normal-flux与原unique immediate boundary规则保持优先。新分支只在附有显式polynomial boundary history时尝试；若能力不具备或内容不恒定，执行原no-hit/affine/phase算法。单独记录uniformReachableQueries；不计为characteristic证书。最终transfer和inventory assembly继续绑定原interval/source/history。

## 比较与限制

原不带history的全部算例与结果应保持。原带history的immediate-hit和有非零gas库存测试也应保持。新增常量例对解析总通量乘c独立检验，并运行原four drivers全部面、原full30/400步。所有失败原样报告。相通量积分误差仍受原预算控制；未证明的非恒定混合来源与I12曲率/时空收敛仍须逐项完成。
