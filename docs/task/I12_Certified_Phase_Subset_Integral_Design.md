# I12a15r4 / 子区域积分证书设计

见 [TODO](I12a15r4_Certified_Phase_Subset_Integral_TODO.md) 与 [r3](I12_Authoritative_Q4_Affine_Pullback_Design.md)。改变的数学事实是：不把一个相未决 tile 的全部面—时间测度都计入液体不确定性，而对严格认证的子区域计量。原守恒库存独立于 geometry 的路线不变。

对同一闭 donor 覆盖中的 forms `c_i+s_i q_s+t_i q_t+R_i`，选首个 form 的 c/s/t 为公共一次部分。每个系数差按 outward 算术移入剩余区间，取所有剩余区间的 hull。于是任一真实 donor 的 phi(B) 都位于共同 form 中；这不是平均相场或重新拟合几何。仅在完整覆盖与同一原 RK2 branch 已证明后使用。

若公共 form 为 `l(q)+[rlo,rhi]`，集合 `P_l={l+rhi<0}` 是明确液相，`P_g={l+rlo>0}` 是明确气相。两者在同一参数方形内不相交；其余保持未决。转换 `q=2u-1` 后，精确 Expansion 保存半平面常数/两个系数，沿原四条方形边用精确有理交点裁剪。复用原 Snapshot 的 rational 加减乘除、面积和第一 u 矩。非恒定线的等号边界零测度；恒等零必须返回空认证区域，不能把整个方形当作液体。

XY 的物理测度为 length*duration*area；RZ 横向面为 `2*pi*length*duration*(r0*area+(r1-r0)*moment)`，竖向面为 `2*pi*r*length*duration*area`。面积/矩/长度/时间先精确有理组合，最后对正分母商向外包围；pi 仍用原 strict 区间。没有用近似交点或 rounded polygon 面积宣称证书。

记已存在的 gross 总量为 T±，点态 gross 速度最小界为 v±，两个子区域的真实测度为 A_l/A_g。令 wl_min/wl_max/wg_min/wg_max 为所有已覆盖 donor 原库存权重的界（几何模式为 1/1/0/0）。则：

`L± >= wg_min*T± + (wl_min-wg_min)*v±*A_l`

`L± <= wl_max*T± - (wl_max-wg_max)*v±*A_g`。

所有运算向外，子区域贡献只使用非负下界。用事实 `subset gross <= total gross` 降低其表示下界到不超过 T 的下界，可保持原组件界的单调排序；这是降低一个证书下界，不改物理积分或库存。结果与原 a15 bounds 相交。速度反号使点态最小界为零时该项不提供改善，原完整分块与 gross 正负界仍继续细化；不能用均值流速抵消 counterflow。

这条不等式允许先用全体候选 donor 的共同符号和权重范围，不要求此项反解每条 donor 边界在参数域中的曲线。它可能保守，但不会漏覆盖；任一未解决分支仍 pending。误差验收仍用原四个积分区间宽度的 outward 和，不另设局部替代预算。

实现细节：共同 band 中相同 binary64 一次系数的实数差严格为零，应直接保留零，不能因通用 outward 加法产生人工 subnormal 裁剪事件。真实非零差仍向外包围。点态速度最小界用 tile 两个时间端点的空间 range hull，加原 8eps/8tiny 混合包络；同一实数线性时间插值在端点 hull 内，避免对常速时间重复依赖导致不必要的零下界。
