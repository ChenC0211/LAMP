# I12 / 认证相通量的舍入可行域

r6的精确循环流已在原认证区间内找到共享通量，四个第二步half库存装配通过。然而原舍入修复只知道[0,total]和64epsilon*min(C)，不知道有些面liquid upper只有约1e-26；它增加这些通量后可以满足库存界，却越过原积分证书。r6的最终核验正确拒绝了这些结果。

本项不给积分、几何或守恒增加容限。只在原库存装配核函数中增加可选完整的liquid±上下界。已有rate/algebraic API不带这些bounds，保持原运算路径。仅r6的真实BoundedFaceIntegral批量入口把同一证书、经[0,total]相交后的上下界交给核函数；单位仍是已积分体积，故该可选能力仅接受unitVolume=1。所有bounds须有限、有序、包含当前选中流量且不超过total。

原增加分支可用余量由total-liquid缩为upper-liquid，仍取旧舍入预算及邻居库存空间的最小值；实际nextafter或浮点相加后再严格检查upper和旧预算。若不能增加，不得退回忽略证书的路径。

对带证书的减小分支，各方向的可减库存为liquid-lower。按原incoming/outgoing约束计算保留比例f∈[0,1]，仅将可减部分乘f：L_new=lower+(L_old-lower)*f。f=1直接保持原值；浮点结果与数学区间[lower,L_old]相交，必要时向lower取一个相邻double以保证有表示变化，随后核对原64epsilon预算。每次仍重新汇总同一共享面通量；不更改cell M来满足约束。如果没有合法可减余量或没有可表示的进展，明确拒绝。原128次上限保持。

无bounds时保留原incoming/outgoing求和、比例和乘法表达式，避免改变已有legacy/普通API结果。新分支没有另一个库存owner：既有private finish仍唯一执行CFL/divergence、严格库存界、舍入修复和global ledger。两个公开入口通过同一private批量装配/发布方法进入该kernel，最终仍由r6再复核原证书、total与boundary不变。

测试首先构造旧修复会超过给定upper的空/满、正/负、非零lower和极小通量反例，要求新修复利用其他可用变化且不触碰lower/upper；完全固定的不可行流量、超64epsilon及错误bounds仍拒绝。再运行r6相同原full，确认三个证书越界消失。之后的第三步来源认证属于独立问题，不放宽本项范围或误称完整I12通过。
