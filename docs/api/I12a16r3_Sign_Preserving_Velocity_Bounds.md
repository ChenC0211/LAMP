# 原速度表达式的符号认证

见 [TODO](../task/I12a16r3_Sign_Preserving_Velocity_Bounds_TODO.md)。新增encloseSignPreserving对同一owned slab的值与四个一阶导数包围，不创建另一份速度源。原evaluate和enclose保持。

在正scale的原坐标归一化中，输入整个区间在origin一侧时，减法及正数除法的舍入保留非负/非正性；退化坐标等于origin时结果严格为零。随后仍按原component的两层Horner和解析导数递推计算，以现有Number区间乘法/加法的符号约束去掉不可能的跨零舍入邻点。没有“小于epsilon视为零”。

原theta的时刻减法/正数除法单调，查询时间在闭slab中，故theta∈[0,1]，RN(1-theta)也在[0,1]。混合包含两个endpoint、theta=.5均值及一般加权分支。两endpoint同号的分量不会因外向舍入变成异号；反向时间速度仍必须保留真实跨零范围。最终结果与原一般enclose相交。

原tight-time证书消费新界及其endpoint界，随后继续原8epsilon/8eta混合误差包络、CFL/axis、initial-hit/no-hit、midpoint/target及r2单调性检查。新query本身不认证characteristic或phase；实际初始/后续物理边界事件不会被声明成内部轨迹。
