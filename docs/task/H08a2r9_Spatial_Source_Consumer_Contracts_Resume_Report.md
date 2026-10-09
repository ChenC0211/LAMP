# H08a2r9 / 空间来源消费者契约迁移完成

2026-10-05，**PASS（消费者契约）**。H08a2r10 关闭精度先决后，本次完成原 phase/momentum 两份消费者迁移，production 0/test 2。原[前置 FAIL](H08a2r9_Spatial_Source_Consumer_Contracts_Report.md)保留历史；本报告是恢复后的当前裁决。[任务书](H08a2r9_Spatial_Source_Consumer_Contracts_TODO.md)、[恢复设计](H08a2r9_Spatial_Source_Consumer_Contracts_Resume_Design.md)。

实际 B/J/atom、材料 F/R/uR/KEP 及 characteristic 未变。下一唯一任务为 [H08a2r11 同步部署](H08a2r11_Common_Regular_Source_Deployment_TODO.md)，本项不代表 H08a 时空材料积分或 H08a3 桥接完成。

## 迁移后的验收

| 原检查 | 当前检查与保留要求 |
| --- | --- |
| phase 单点中心差分强制 B-div=DI/V | 实际分片的物理体积分与原 native incidence/DI 比较；208基、极大/极小缩放、叠加、抵消和零场保持 |
| 独立再拟合旧 tangent/trace/lift | 删除测试中的旧拟合/重构副本；保留独立解析输入均值、原解析值/J/native中心精度；atom primitive 对实际查询的独立GL3线积分 |
| dual flux=Σ overlap.volume×cell平均div | 边界通量对同一实际场在真实 overlap 的物理散度积分；旧配分差单列诊断 |
| 三个volume坐标节点微分旧二次场 | 六个片内物理Chebyshev节点，以完整次数Lagrange导数微分值；返回J不进入参考计算 |
| 点态常散度、切向跳跃必须非零、旧存储常数 | 点态/配分/跳跃仅诊断；实际atom/view库存和原空间覆盖保持，旧rank/side-piece统计不作新表示门槛 |

XY 积分 `∂x u_x+∂y u_y`；RZ 积分 `2π[r(∂r u_r+∂z u_z)+u_r]`，轴上不除 r。每片密度最高五次，GL3×GL3 足阶。分界缓存来自实际 `pieces()`，正则暂存场来自同一准备对象的 `supportKnots`；真实 overlap 单独裁剪分片。误差尺度累加求积点上各物理项绝对贡献，避免净通量抵消掩盖尺度。

导数参考在 RZ 微分 `(r*u_r,u_z)`，可覆盖五次加权速度。非轴处恢复物理导数，轴处使用 `∂r u_r=∂rr(r*u_r)/2`，并保留独立轴奇偶检查。六节点的实际舍入物理位置用于插值，避免把理想节点当作实际查询位置。

原解析资格集合直接按冻结矩阵核对：64布局中56个适用、8个不适用，不能被未来某个rank标签缩小。解析值/J仍 **1e-12**，值采样导数仍 **2e-11**；没有靠修改容差放行。

## 原覆盖保持

phase 原 **208 个非轴 native 基、236 次更新、13,528 面检查、5,528 cell检查、4,664 零gross cell**均保持。新体积分最大尺度化差 **1.47e-16**。真实native均值、原DI/独立incidence、完整atom双端primitive、线性、1e±200缩放及抵消输入通过。原phase几何/相率闭式参考及来源/clock/wall拒绝未改。

momentum 原 **64布局/128输入、114,300 atoms、72,080 cell点、247,560视图**保持，解析值最大差 **2.23e-16**。原 **144个jet场景、69,804点、696轴点、1,008内部knot、23,898实际native面**保持；新导数参考最大差 **2.34e-13**，原解析J最大差 **5.14e-15**。原几何身份、完整时钟、生命周期/regrid、zero/tiny/极小width及错误来源拒绝保持。

旧 `cellPieces` 报告统计的是侧片段，不是矩形 `pieces().size()`。首轮迁移误把两者作等值比较，短契约报错；读明实现后移除该错误统计假设，保留真实 atom/view 库存及实际piece闭包/积分检查。失败日志与当时两测试源码保存于 `pilot_source`，未改生产或精度门槛。

## 两种来源的同一数学检查

XY/RZ × uniform/2:1 × n4/8/16 × full/half/hole/disconnected，共 **48布局**；无wall/全外部wall/单wall，原解析/非线性/随机native三种输入，共 **432组输入**。每组分别检查当前接受场与正则暂存场；同一个模板执行同样的求积、导数、原解析资格和阈值，来源标号只用于诊断记录。

| 同一数学检查（两来源合计） | 数量 | 最大尺度化差 |
| --- | ---: | ---: |
| 由实际值独立微分的J | 406,728点 | 4.70e-13 |
| cell物理体积分/原DI | 80,136 | 5.37e-15 |
| dual真实边界/overlap体积分 | 183,204 | 3.86e-15 |
| native真实面均值 | 183,204 | 4.20e-15 |
| 原允许解析族值/J/native中心 | 全部适用输入 | 5.22e-14 |

正则场的点态div与cell平均最大绝对差 **28.3631**，dual边界通量与旧平均div配分最大差 **0.195326**；上述真实物理积分仍通过。当前旧场相应差为 5.33e-15/1.81e-16。由此确认新判据确实解除旧过约束，而不是只让仍满足旧公式的场通过。

测试端三项故意扰动均被同一独立物理参考检测：normal反号、高次项截断、遗漏支持分片的尺度化差分别约 **0.0370、0.00540、0.0127**。原错来源拒绝代码保持。扰动不进入生产，也不生成备用算法。

永久短契约保留 n4 的16布局/144输入；完整矩阵由捕获源码/runtime的临时入口运行，约9.65秒，入口现已退役。原r8/r10的完整次数、C1、native面矩、缺秩低阶和guard契约继续执行。

## 回归与交付边界

最终标准 Release 及 **9/9 CTest PASS**，约197.68秒。生产、输入和未涉及源码哈希保持；1,382项库存只有两份测试变化。material/dual/KEP/characteristic和benchmark可执行哈希保持；只重建两测试runtime。

对原 **10,413行**输出逐项核对：原覆盖、未迁移数值及拒绝保持；明确单列已迁移参考误差列、诊断与计时。新增匹配矩阵另独立校验，不能把排除后的对照称为全部输出逐字相同。上一轮373冻结payload保持。完整H10b未重跑，因本项未切换实际来源或算子。

检查：`runs/checks/H08a2r9/20261005_164956_929`；冻结：`runs/baselines/H08a2r9_resume_20261005`。`audit_full`绑定完整矩阵源码/runtime；`LastTest_final.log`、`verification.json`记录原覆盖对照、扰动敏感性、生产未改和退役证明。旧FAIL冻结不被覆盖。

材料仍仅唯一共享F、内部反号/真实AMR子面和、来源一致、KEP/skew四硬约束。H08a整体、FirstJet桥接及完整Hydro App尚未完成；H11未启动。
