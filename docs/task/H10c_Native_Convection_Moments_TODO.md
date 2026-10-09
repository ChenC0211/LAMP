# H10c / 最终native对流算子的联合空间矩与动能约束

状态：**FAIL**，2026-10-04；[报告](H10c_Native_Convection_Moments_Report.md)。T0/T1/T2通过，原内部affine反例修复；完整T3空间门槛仍失败（28窗口/15,116步，14原失败门槛保留）。下一单任务[H10d整体一致性收口](H10d_Final_Operator_Consistency_Closure_TODO.md)，不进入H11或H08a。任务书及[实现前设计](H10c_Native_Convection_Moments_Design.md)先于源码，遵守[路线图§0](../roadmap/05_LAMP_Module_Task_Roadmap.md)。

## 目标

只解决一个数学事实：最终native convection operator同时满足内部constant/affine空间一致性及已有material-KE/skew关系，消除负转置配对造成的内部PDE矩缺陷。对象为最终载荷，不继续调整单段trace来代替载荷验收。

## 允许修改

- `include/Modules/Hydrodynamics/CchKepConvection.h`、`src/Modules/Hydrodynamics/CchKepConvection.cpp`：只读最终算子系数/矩诊断及统一one-ring联合最小改动。
- `tests/contracts/hydrodynamics/cch_kep_convection_contract.cpp`：原门槛保持，独立reference计入新系数并追加内部矩检查。
- `cases/benchmarks/cch_convection/driver.cpp`及研究README：追加最终算子诊断，既有full物理参数/门槛不变。
- 本TODO/设计/报告/API、现行H10/H07设计、路线图及索引；登记用户要求的H10c后H08同源时空相通量升级，暂不实现。
- 证据 `runs/baselines/H10c_20261003`；临时目录由start.json登记于 `runs/checks/H10c/<timestamp>`；正式完整研究仍由通用run_case产生。

预算2 production +1 research source、1 test source；不改CMake/工具/字典。

## 禁止修改

H01–H09、Interface/Core/Solvers、旧Hydro/Compatibility/App、旧算例输入、冻结PDF及历史TODO/报告。不得添加AMR分支、第二材料通量、守恒库存、bubble、解析状态修正、人工黏度或upwind fallback。H11、H08时空升级和生产H14/H15边界不在本轮。

## 数学约束

- 保留唯一H09实际F/质量/R/来源/clock，保留 `uR` 和 `0.5*u^2*R`。
- 最终N=N_old+J*u；J逐共享native pair反对称，J*1=0（机器精度），故保持内部动量抵消、常速度及KE/skew恒等式。
- one-ring由实际pressure cell邻接和同分量native faces统一生成，重叠patch共用反对称系数，不进行逐patch修正后互相覆盖；无AMR特例、跨不连通域支持或新primary未知量。
- 内部row的线性矩目标由同一F与真实segment位置构造；局部中心二次矩为0作为native点处二阶输运目标，不能以几何体积平均位置偏移冒充native点导数。外部row保留constant/skew约束、显式现有端口，边界误差仍全部进入H10b。
- 局部小型Gram/SVD求最小改动与可识别性；共享系数的联合约束须迭代收敛并显式检查，不留未收敛系数、不静默降低阶数。

## 实现要求

- [x] 核验H10b 43payload、源码/旧输入/PDF/App哈希，保存开始任务书/设计。
- [x] 不可变prepared系数、same-owner绑定和只读矩诊断；finite/尺度/收敛拒绝，无IO/clock/接受事务。
- [x] 最终算子的独立dense reference、内部affine/二次矩、旧随机动量/能量/常速度/材料跨界与regrid契约。
- [x] 先运行原内部affine反例，过后才启动既有H10b完整研究；保持全部失败行和84条原门槛。
- [x] Release、受影响契约/架构/目录/前置回归、来源与范围审计及归档。

## 验收 T0/T1/T2/T3

**T0：** 标准Release无新增warning/error；生产源唯一CCH owner，旧App/PDF/输入不变、默认CTest无长研究；只读API及误用拒绝保持。

**T1：** XY/RZ uniform/2:1/selected域，独立dense载荷/随机动量/KE原1e−12门槛保持，常速度原1e−10保持；共享pair严格反对称、同分量/one-ring支持。constant及内部线性/二次矩按gross系数误差≤1e−12；联合求解未达到该门槛明确拒绝。

**T2：** 原XY 2:1内部n=16、dof911反例应为−1，实际误差≤1e−10；n=8/16 XY/RZ uniform/2:1全部内部affine native记录≤1e−10。正式非零R/uR、密度比1/1e3/1e6及真实regrid原契约继续通过。

**T3：** T1/T2通过后重跑H10b完整28窗口/15,116步；四组末档空间/截断阶≥1.7、最细L2≤0.05、Linf下降及原动能门槛均通过，独立时间阶≥3.5。不删外边界自由度、不改变T/dt/n或原门槛。donor原控制裁决原样保留，本任务不把旧有限网格控制失败改成PASS；H10完整关闭仍需对应控制证据。

## 失败规则

矩/恒等式/常速度/来源拒绝不通过，或共享one-ring约束不可识别、不收敛即FAIL；不得保存不合格生产系数或放宽门槛。内部检查通过而完整空间研究失败仍判H10c FAIL，冻结剩余缺陷并先拆下一单任务；不进入H11、不提前做H08。真正无可推进外部依赖才BLOCKED，数学失败不是goal BLOCKED。

## 交付

最终native算子约束/系数/诊断、独立契约、原反例及完整研究结果、来源/哈希证据和明确未解决问题。用户新增顺序：H10c空间离散完成后，H08升级为Pure FirstJet独立几何上的同源时空积分相通量，减小R，继续R/uR与KEP；不恢复库存/bubble。
