# H08a2r11 / 同一正则来源部署审计

2026-10-05，**FAIL（部署验收先决），候选已退役**。实际B/J/native/piece/atom和完整运动证书曾同步切换进行检查；未获准作为当前实现保留。下一唯一任务为[H08a2r12](H08a2r12_Regular_Source_Arithmetic_Consumer_TODO.md)，核清完整Q4表示的误差预算与事件断言；随后[H08a2r13](H08a2r13_Benchmark_Source_Observer_TODO.md)清理仍依赖旧patch的研究观察器，再恢复r11。不是将任一失败门槛自动改判PASS。

本轮按[任务书](H08a2r11_Common_Regular_Source_Deployment_TODO.md)和[设计](H08a2r11_Common_Regular_Source_Deployment_Design.md)工作。候选4 production、2 tests；最终production 0、tests 0，启动时1,382项源码/输入全部恢复。Pure FirstJet、材料F/R/uR/KEP、最终native算子及原H10b输入均未修改。没有保留候选分支或fallback。

## 原验收发现的缺口

原phase契约的系数线性叠加断言失败。独立临时入口重放原XY/RZ × uniform/AMR n4的两组native输入及其浮点和，共784 atom、3,920系数：**121 atom中的145系数超出原1e-12门槛**，其中三次26项、四次119项。最大原尺度化误差**8.30056e-12**，没有降低门槛或删掉输入。

最坏见证为XY AMR、atom94、owner26，法向x、区间y=[0.625,0.6875]。前两场的四次系数均为0，叠加场为−1.4551915228366852e-11。原检查把不同物理量纲的各次系数直接相加作为work，未乘区间长度的幂；高次系数的残差不能直接等同于速度误差。

补充诊断按实际区间长度L计算 `sum |δc_k| L^k / sum (|a_k|+|b_k|) L^k`。以保存的浮点系数作100位Decimal精确输入复算，全部atom最大**6.14038e-15**。它只衡量已发表多项式的叠加残差，不证明准备权重的后向误差、任意tiny/RZ积分或整个材料模块通过；本轮不以此代换原gate。下一任务须给出具有量纲一致性、单位变换不变性及独立负控的数值验收定义，再决定是否确需改进算术。

原characteristic精确常数fixture也失败，但**事件位置和语义标志未失败**：末点x=.25、oneSided=true、events=1、splitEvents=0；Jxx=0.99999999999999933、Jxy=−2.53270e-16，旧断言要求Jxx字面等于1。实际常数查询包含约1e-14的导数舍入残差。应区分拓扑/来源的严格条件与浮点解析值检查，不能把此见证描述成轨迹失稳或伪造单位Jacobian消除它。后续尚未执行到的事件拒绝不计作已通过。

此外完整benchmark编译发现遗漏的消费者：`cases/benchmarks/cch_convection/driver.cpp` 中spatial观察器读取patch.rank/support，closure观察器重算旧tangent slope。其仍调用已删除的 `patches()`，发生C2039。设计初稿“仅momentum测试读取旧patch”不完整，现明确更正。仅升级五系数CSV尚未迁移这些观察器，不能保留无消费者数值用途的旧拟合作为兼容壳。驱动本轮未修改，前置任务r13须真实迁移其输出与对应参考。

## 候选中已通过的范围

- 原momentum契约完整通过：64布局/128更新、原56适用解析门槛、144 jet场景、regular准备/低阶/完整次数/C1/native均值/Gauss及双来源消费者检查。其通过不替代失败的phase契约。
- 完整物理Q4迹证书候选采用五节点Bernstein换基，六个value/scaled-J分量，按变换绝对项累计work，原128 epsilon门槛；54项常数/三次/四次负控（含旧三节点不可见多项式、1/1e-20/1e-200尺度）通过。
- 两个原非线性跳跃位置fixture改查新C1场：独立RK4 2048/4096步参考差≤1.16e-13，实际轨迹收敛与FD通过；无非透明重启。两种场的数学能力不同，不保留“新场必须有旧跳跃”的要求。
- 原288个实际map Jacobian差分点通过，最大归一化差**4.31e-10**；来源/物理边界8组通过。
- 原2,304次分区比较通过，最大差**1.96e-15**；零/1e-20/普通幅值及两时间端点连续性通过。
- 原24时间族、原7个失败见证全部重验；最小位置阶**2.99969**、J阶**2.99883**，最细误差最大**8.60e-10**，时间门槛0失败。这是实际候选owner运行结果，不是r3暂存运动原型。

常数事件失败后，以临时独立分组入口继续收集分区与时间结果；总退出码仍为2，不能把这些部分通过包装成characteristic整体通过。完整H10o/native/affine与H10b-v2本轮未重跑：T0/原短契约未闭合，benchmark还未能编译候选。历史H10o 28窗口/84项PASS仅属于此前接受实现。

## 恢复、证据与边界

候选和临时入口已捕获到 `runs/checks/H08a2r11/20261005_172448_525`：`candidate`保存系数输入、Q4 DOF、面修正、候选源码/runtime与失败退出；`motion`保存实际24族/分区重验；`coefficient_witnesses.json`保存全部121反例，`coefficient_analysis.json`保存独立Decimal复算。恢复时刷新mtime强制重新编译，原源码哈希全部一致。最终标准构建及**9/9回归PASS，196.89秒**；排除明确计时行/列后，六契约的**10,465行**原输出保持。前次冻结72项payload校验通过。最终核验见同目录 `verification.json`，冻结为 `runs/baselines/H08a2r11_20261005`；恢复检查不计为候选部署通过。

候选的准备统计改为真实regular SVD/权重/pieces；由于部署失败，本轮没有完成新的准备/更新时间成本验收。完整Hydro App、FirstJet消费桥接、真实时空相积分及R收敛仍未完成。材料仍只要求唯一共享F、内部反号/真实AMR子面和、来源一致、KEP/skew四硬约束；不要求R为零，不增加库存、bubble或界面体积修复。H11未启动。
