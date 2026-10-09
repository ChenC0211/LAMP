# H08a2r12 / 完整物理chart浮点验收设计

2026-10-05；[任务书](H08a2r12_Regular_Source_Arithmetic_Consumer_TODO.md)。单一事实：数值消费者的误差预算应针对同一物理量、真实运算尺度，不把旧低次表示的字面系数条件当作新方法硬约束。保持材料四硬约束及原物理门槛。

先按production 0、tests 2（phase、characteristic）实施。基线已完整核验，上一轮属于progress。生产准备与运动方法不更改；若独立物理精度检查证实缺陷，再按TODO设计生产修复。临时研究只在runs/checks。

对端点chart `p(s)=Σ c_k s^k, 0≤s≤L`，无量纲系数 `d_k=c_k L^k` 都是速度量纲。使用 `N_L(c)=Σ|c_k|L^k`；多项式值误差受此范数控制，一阶切向导数用 `Σk|δc_k|L^k` 控制L倍导数。平面积分≤L N，RZ积分≤2π L r_max N。双端chart分别检查，不从两个舍入全局端点相减构造tiny width。若长度单位改为α、速度单位β，系数改为β c_k/α^k、L改为αL，范数和work同乘|β|，相对预算不变。

避免只看最终被抵消的系数：对固定geometry/boundary，逐个有效native自由度查询同一来源的单位响应T_i，用 `W(v)=Σ |v_i| N_L(T_i)` 作可复算的运算尺度。它是实际浮点线性响应观察，不独立证明重构正确；独立解析/物理值/J/GL3积分和原56布局/208basis负责准确性。原phase的208basis可复用构造该work，保留旧raw系数结果为诊断。

另在XY/RZ、uniform/AMR、full/half/hole/disconnected和无墙/全墙/单墙矩阵比较现有source与regular暂存场。共同区间由原atom与实际supportKnots细分，真实值/J来自对应source。原1e-12门槛作用于量纲一致且无绝对floor的误差。使用TwoSum分离native加法舍入，再以单位响应重组，独立100位Decimal复算已导出的系数/输入/work；对零、1e±200、随机和强抵消均检查。错符号、删高项、错width及漏support必须被独立物理参考检出。

实施中补充（新增检查定位）：旧来源RZ disconnected全墙、极小幅值的一条零functional，所有单位响应终值均为0，但密集输入抵消得到约1.30e-217。故逐atom的最终响应仍不足以覆盖运算抵消；新增完整场检查采用向量范数 `||T_i||=max_(region,chart) N_L(T_i)`，`W(v)=Σ|v_i| ||T_i||`。这是整个输出向量的normwise线性检验，不是每条零functional的相对误差，也不是有证明的逐操作舍入上界；该区别须在报告中保留。它由真实列响应可复算，随输入线性缩放，没有绝对floor。原208basis/local物理门槛与独立tiny相积分保持。

实际候选常数中间轨迹(.3125→.1875)的事件数由旧1变为2，新增一个真实support边界，不是多报事件。事件计数仍严格相等，但参考改为独立遍历实际piece闭包与直线轨迹的交点（排除起点、包含精确终点），不硬编码旧分片数量。坐标、oneSided、splitEvents、身份与拒绝仍保持原意义。

完整phase候选补验又定位同类native零面矩尺度问题：RZ uniform n4，native输入u20=5.0929581789406511、其它0；面u32期望0，积分−2.90656e-17，旧分母仅残余chart work 2.86114e-16，误判相对误差0.101587。原208basis先全部构造并缓存，仍逐一检查，不增加/删减输入或重复几何准备；以这些响应的完整场列范数建立fieldWork。零/非零面均值与积分Gauss预算增加同量纲 `area*fieldWork`（cell按真实边界面积相加），门槛仍1e-12。原DI独立incidence、zero-native-gross的平均div严格为0、解析MMS、tiny相对精度及全部负控保持。此为normwise物理通量验收，不宣称每条零流量都能有有限的局部相对误差。

常数事件保持精确末点、oneSided、piece、事件计数、外边界/来源/full-clock/工作预算拒绝；仅将非零时长的J字面等式改为原解析矩阵已有1e-12误差要求，并检查全部四项及同一map FD。零时长身份映射仍严格相等。误差判定及负控应在不变生产source与已冻结r11候选上分别运行。后者通过runs/checks中的独立CMake研究目标链接捕获的四份候选源码与当前未改库，绝不安装/替换当前库；旧跳跃专属fixture仍归实际source数学能力迁移，不在本任务静默删除。

结束前：复算原145系数反例、验证单位变换与负控，原24时间族/288FD/2304分区矩阵及全部受影响契约。候选研究与当前生产证据分开，候选其它失败如实归后续部署先决，不扩大本任务通过范围。
