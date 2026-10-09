# H08a2r12 / 正则来源浮点消费者契约

状态：**PASS（消费者契约）**，见[报告](H08a2r12_Regular_Source_Arithmetic_Consumer_Report.md)与[误差预算设计](H08a2r12_Regular_Source_Arithmetic_Consumer_Design.md)。来自[r11部署审计](H08a2r11_Common_Regular_Source_Deployment_Report.md)，属于恢复实际部署的先决，不是新数值路线。production 0、tests 2；下一唯一任务[r13](H08a2r13_Benchmark_Source_Observer_TODO.md)。

目标：判明完整Q4物理chart的线性/缩放误差与常数轨迹事件见证中哪些是实际精度缺陷、哪些是旧表示专属断言；建立可复核的单位一致预算，保留原物理精度及全部拒绝语义。不得只把失败阈值改大。

允许修改：production原则上0；确有算术精度缺陷时先设计，仅既有CchRegularAdvectorGeometry.h/.cpp（≤2）。test≤2：原phase与characteristic契约。不得在此任务切换实际source、修改材料F/R/uR、最终native算子、FirstJet或H10b输入。临时审计仅runs/checks，无长期工具或生产开关。

数学约束：材料仅唯一共享F、内部反号/真实AMR子面和、来源一致、KEP/skew四硬约束；实际R/uR保留，不要求R=0或恢复库存。实现要求见误差预算设计，区分严格拓扑/来源约束与量纲一致的数值误差。

失败规则：独立物理精度/负控或原来源拒绝退化即FAIL；真正缺少外部依赖才可BLOCKED，不以耗时、内部实现难度或表示换算问题声称外部阻塞。

T0：复现r11原784atom/3920系数/145失败及常数末点J见证，绑定原数据、源码/runtime；保持原观察值作为诊断。生产变更前定义推导和量纲。

T1：同时检查当前来源和同一owner的regular暂存场。线性/缩放覆盖完整双chart、实际support、零/1e±200/抵消、轴/墙/AMR/selected；在无量纲chart统一各次项，并独立检查物理值/J/积分。使用独立高精度重组或可复算后向误差工作量，区分native加法舍入与系数运算；禁止用近零的最终系数掩盖实际工作量，也禁止任意绝对floor。坐标/单位变换不得改变通过结论；遗漏高次、错误符号、错误width/support的反例必须被检出。原独立解析/物理1e-12等门槛和原208basis/56布局适用性不缩减。

T2：常数fixture将严格事件位置、piece选择、oneSided、事件计数、身份/full clock拒绝与数值Jacobian误差分开。原终点(.25,.43137)、初始/同时/擦边/预算拒绝不能删除或弱化；数值部分须以同一实际map FD和既有精度预算为依据，不能将J强置单位或使用第二场。原24时间族、288FD及2304分区测试仍必需，零场与极小幅值不因floor失去敏感性。

交付：误差预算推导、负控、原覆盖对照、真实算术缺陷或无缺陷的证据，标准构建及受影响契约。任何实际物理精度退化按FAIL处理；表示尺度诊断不是物理通过的替代。通过后执行[r13](H08a2r13_Benchmark_Source_Observer_TODO.md)，再恢复r11完整部署及H10o/H10b-v2。H11未启动。
