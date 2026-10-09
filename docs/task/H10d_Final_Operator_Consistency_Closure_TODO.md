# H10d / 最终native算子整体空间一致性收口

状态：**FAIL，2026-10-04**。独立全部native逐行审计完成；成对周期闭合的精确中心P2候选在2:1网格不相容，硬C/P1加P2惩罚候选仍有加密系数增长且最细准备失败。失败实现已归档并从生产源码退役，现行H10c内部affine/材料KE保持。见[报告](H10d_Final_Operator_Consistency_Report.md)，下一单任务[H10e有界最终算子一致性](H10e_Bounded_Operator_Consistency_TODO.md)。H10未关闭；不进入H11/H08a。

目标：最终native convection operator连同显式边界端口满足整体二阶空间一致性，并保持同一材料F/R/uR、现有material-KE/skew恒等式。内部低阶矩精确并不足以证明高阶余项有正确h尺度，整体空间精度与边界作用必须独立定位，不继续打单段trace补丁。

允许修改：实现前先冻结最终H10c源码/证据，审计每个native row的完整smooth truncation及边界/内部来源，形成联合矩可识别性设计。实现预算仍为CchKepConvection.h/.cpp、现有短契约、现有研究driver及必要文档；需要额外文件或方法改变时先拆任务。研究数据仍由通用启动器进入时间戳runs，不增加专属工具。

禁止修改：Interface、H08/H09真实几何质量/材料来源、旧App、历史报告/冻结PDF；禁止AMR特例、解析状态修正、人工黏度、upwind fallback、删边界自由度或恢复库存/bubble。不提前做H11、H08a或完整H14/H15。

数学约束：保持H10c已验证内部常数/一次/二次矩及旧随机动量、非零R/uR、材料KE关系。外部trace数据、native采样与实际边界关联必须在同一最终operator/lift中明确；若现有端口数据不足以支持约束，先明确API/拓扑需求和可识别性，不暗中改变端口功公式或质量定义。

实现要求：统一真实几何one-ring约束最小改动，无独立AMR分支。区分“内部moment精确”“高阶系数余项具有正确h尺度”“全部native瞬时PDE再现”“完整轨迹收敛”；剩余误差应由独立分解验证，不能只凭能量恒等式或终值某一范数推断。H10c较大准备开销如实保留，性能改造另行划定目标。

验收：T0标准Release/来源/原App和输入不变；T1现有124场景、1,516内部矩、one-ring支持、1e−12载荷/KE与1e−10常速度门槛保持，补充最终外部row/端口矩的独立反例；T2原affine n=8/16及最细四组保持，无解析数据进入accepted状态；T3原H10b的28窗口/15,116步/n=8/16/32/T=1/84门槛原样重跑，整体截断/终值阶≥1.7、最细L2≤0.05、Linf下降、能量/时间门槛通过。原donor有限网格控制失败继续单列保留，不能事后调整其原门槛。

失败规则：联合矩不相容、来源/恒等式失败或完整空间gate失败判FAIL并保留证据；无旧算子fallback。仅真正外部依赖阻断才BLOCKED。H10空间离散通过后再按用户顺序启动[H08a](H08a_SpaceTime_Material_Flux_Upgrade_TODO.md)。

交付：实现前设计、最终operator/端口契约、内部与外部独立误差分解、原完整研究及明确剩余限制。不把内部PASS改称H10整体PASS。

## 本轮执行边界（实现前）

- 只改`CchKepConvection.h/.cpp`、现有`cch_kep_convection_contract.cpp`、现有`cases/benchmarks/cch_convection/driver.cpp`/README及对应文档；预算2 production +1 research +1 test，不改CMake/工具/字典。
- 先只向研究driver加入被动逐行矩阵审计：独立原K+J+现有显式周期trace矩阵与真实apply载荷核对≤1e−12；内部/外部/固定轴分类的误差平方贡献必须加回完整误差，不能删除行。报告局部wrapped一次/二次/三次矩、行系数尺度及实际误差。
- n=8/16/32、四原组合的瞬时定位放`runs/checks/H10d/<timestamp>`；开始源码/输入/旧App/PDF哈希和全部证据放`runs/baselines/H10d_20261004`。瞬时诊断不是完整窗口验收。
- 数学修正必须先写实现前设计并定义明确的来源/拓扑/端口与KE关系，再改生产源；不得让研究参考参与生产系数或轨迹。
- 新最终端口/空间矩/KE契约通过，原内部affine及材料/来源契约保持后，才重跑原28窗口/15,116步/84门槛。原空间门槛及donor失败不改变，不进入H11/H08a。

## 实际关闭

开始100个H10c payload核验保持。原算子的四组合n=8/16/32全部行审计完成，独立矩阵与真实apply一致，分类贡献加回整体误差。联合可识别性与惩罚候选均保留失败证据；T1新闭合/有界余项前置失败，所以本轮未重跑完整H10b，也未声称新端口API通过。最终2个生产文件及原测试恢复为H10c精确字节，仅研究driver增加被动rows/rows32观察器与README解释。原内部affine n=8/16复核通过，旧契约与架构/目录共4项通过；最终源码/输入/实际运行归档见报告。
