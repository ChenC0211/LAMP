# I12a11r12 / Corner-cut bubble 条件性算子验证

## 目标与范围

只检验一个数学事实：固定 Q4 C1 skeleton 时，cardinal center bubble 对趋近角点的小相体积的灵敏度是否按 cut 尺寸五次方退化，以及实际原 I10 求根是否吻合这一解析预言。0 production + 1 test：`tests/interface_hermite_gals_benchmark.cpp`；另含本 TODO/报告、README 与 `runs/baselines/I12a11r12_20261001/`。不改任何生产几何、M 更新、B/J、容差、实际 benchmark 或冻结文件。

## 数学与实现

4×4 均匀网格中 cell 5 的单位 chart（左下 [.25,.25]，物理 h=.25），解析场 phi=±((x+y−.5)/h−a)，a=2^-2,…,2^-8。该单元较小相 / capacity 为 a²/2；其它单元取同一平面的独立解析裁剪体积。将角点试验置于较长的真实平面链中，以满足原 complete discovery 至少三个 retained graphs 的条件；首轮外边界单个 arc 的不合格 fixture 拒绝日志保留。两种符号分别测试小液相/小气相。以下 x、y 为所选单元的单位 chart 坐标，center basis 为 B=256 x²(1−x)² y²(1−y)²；沿 x+y=a，归一化体积导数为

`dV/dc = -256 a^5 ∫_0^1 u²(1-u)²(1-a*u)²(1-a+a*u)² du`。

用测试内独立多项式乘积及逐幂原函数计算积分，与原 discovery sensitivity 比较。a→0 时导数 / a^5→−256/30；不从生产积分节点生成参考。实际 target 在解析液量上加固定 128 eps capacity，使用原 I10 工作/激活/位移预算和真实边界 contacts 求根。该扰动是算子制造量，不进入任何真实 benchmark 库存。

物理导数为上述值乘 capacity。核验初始体积、实际闭合残差、系数响应 Δc≈ΔV/(dV/dc)、固定 skeleton、center Hessian 改变 Δphi_xx=Δphi_yy=−16Δc/h²、Δphi_xy=0。报告系数/内场响应随小相缩小的增长，不把几何条件性测试通过当数值收敛通过，也不据此放宽体积预算。

## 验收与失败规则

T0：Release 无 warning；11 项直接回归输出与 r11 一致；其它 27 可执行、production、81 冻结/参考文件、3 PDF、r11 archive 不变。

T1：14 组完整；独立 derivative 相对误差≤1e-9，最后两级 sensitivity 斜率在 [4.9,5.1]；实际 root 达到原工作预算，线性响应的体积预测缺陷≤原工作预算量级，center jet 独立检查≤1e-10。任何拒绝或参考失败照实归档 FAIL，不增加求根工作量或修正幅值上限。

交付原始参数/日志/误差、源码与 runtime 哈希、报告及 archive。后续优先检查现有路线的减小步长/局部加密能否控制已定位的系数响应；改变状态分解或放宽预算另需具体设计变更。
