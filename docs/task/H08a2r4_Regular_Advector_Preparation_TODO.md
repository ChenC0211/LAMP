# H08a2r4 / 共用正则advector的生产几何准备

状态：**PASS（生产准备层）**，见[报告](H08a2r4_Regular_Advector_Preparation_Report.md)。144结构组、原12+6 source布局及9/9回归通过，旧4727+3364确定性行保持。下一[H08a2r5](H08a2r5_Material_Reference_Polynomial_Degree_TODO.md)准备材料/dual的完整次数参考消费者。[实施前设计](H08a2r4_Regular_Advector_Preparation_Design_Change.md)固定边界模式、所有权及文件预算。[H08a2r3可行性PASS](H08a2r3_Common_Regular_Advector_Report.md)是输入证据，不表示生产B已经更换。H08a3消费目标保留，H11未启动。

目标：把已通过的native均值P2重构、共享Q4约束及实际face均值校正的**几何线性映射**纳入现有`CchAdvectorReconstruction`，使后续速度evaluate只应用prepared weights。只解决准备层的来源、秩、C1、面矩及成本，不在本项同时更换完整field/atom/characteristic消费者。

允许修改：实施前固定≤5 production/build、≤2 test/build清单。候选范围为现有`CchSpatialAdvector.h/.cpp`、职责独立的`include/Modules/Hydrodynamics/Detail/CchRegularAdvectorGeometry.h`及对应`src/Modules/Hydrodynamics/CchRegularAdvectorGeometry.cpp`、根CMake；测试优先扩展原momentum契约，必要时仅调整tests CMake。新helper必须由现有reconstruction实际持有，是其准备实现，不另建可发布速度/材料/clock的owner。若实际部署需超预算，先拆分，禁止把临时Python留作生产依赖。

禁止修改：FirstJet primary/geometry、材料F/R/uR与最终KEP算子、压力/黏性/电场/App/字典。当前实际B及atom数值在准备层过渡任务中保持；不能半切换速度查询而材料F仍旧，也不能把准备对象通过称为H08a通过。

数学约束：四项材料硬约束不变；新几何映射来自实际native faces/selected domain，逐subface修正，XY/RZ真测度，零输入及指定wall/axis条件。共享Q4约束同源，保留缺秩证书及统一最小范数；不强制信息不足布局P2精确，不强制点态常散度或R=0。

实现要求：先核对H08a2r3的全部输入、原型raw数据和消费者表，写实施前设计，特别明确boundary law目前在evaluate绑定而新stencil受law影响的准备策略。应在几何阶段准备必要的有限约束模式/权重，由同一个显式law选择；或采用经设计的不可变law准备绑定。不得在每次evaluate重新SVD或依靠隐式fallback。各布局只能使用实际活动侧和原mesh身份；相同坐标重建不自动等同。保留一圈起始、统一扩展、rank/deficit/范数/存储报告，稳定保存实际支持分界。新准备数据不持有速度、phi、密度、时间层或mesh裸引用。

验收：T0标准Release、架构/目录及原受影响契约；原实际B、atom、F、R及原确定性数值保持，新增被动准备报告单列。T1复用原60布局与144种law/field检查：与冻结原型及Core约束数据对照；C1/native面均值/Gauss/可辨识P2的scaled缺陷≤1e−11，全部缺秩行及独立证书保留；重新验证错误mesh/来源/law、释放原对象、regrid重建、FP环境和零/轴/选定wall。T2在原12收敛布局上应用准备权重至同一native数据，source值/导数与冻结原型逐点比较，并保持原2.5/1.8阶或舍入精确；额外sin三档不能删除。记录C++准备/重复应用的完整n8/16/32成本，应用过程中SVD次数为0，权重/支持规模必须有界增长。

失败规则：来源、C1、面矩、边界、精度、未声明rank缺口、应用重复求解或最细规模不可完成，任一FAIL；失败实现归档退役，原已通过来源保持。只有真实外部依赖阻断才BLOCKED。不得借本准备任务恢复旧材料过约束。

交付：唯一owner中的准备实现、API/报告、原型逐点证据、成本与拒绝范围、source/runtime冻结；列明生产B尚未切换。随后按预算迁移旧三系数消费者并部署同一正则B/atom/J，重验H10o最终native/原affine/H10b-v2完整28窗口84项，再关闭H08a3与H08a时空材料积分。一次仍只推进一项。
