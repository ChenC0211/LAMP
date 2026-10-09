# 原RK2的新壁面接触排除

范围见 [TODO](../task/I12a16r2_Monotone_Wall_Exclusion_TODO.md)。沿用 [r1边界来源](I12a16r1_Boundary_Phase_Source.md) 与原完整闭块no-hit证书，只加强最后的新壁面接触判断。

对任一坐标，原target的逐运算顺序是 P−(W*d)。已有证书覆盖d>=0及有限运算。W<=0时，IEEE逐步舍入的W*d<=0，随后P−(W*d)>=P；P严格大于lower wall时target也严格大于该壁，P就在该壁时则不构成原touchesNewWall的“从内部新触及”。因此W.upper<=0能排除lower-wall的新接触。同理W.lower>=0排除upper wall。零速度、负零、正duration乘积下溢为零仍遵守弱单调性。

这项证明不声称实际离开点严格远离壁面，也不丢弃t=旧时间或空间端点。原midpoint与target在域内、初始无hit、CFL/轴/覆盖验证先行；方向无法认证时仍采用原严格空间界条件。原已成功查询的数值界与实际tracer完全保持。
