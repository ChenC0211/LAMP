# BF-01 / Regularized color 的唯一只读派生

状态：**PASS（仅BF-01）**。实际实现和完整证据见 [报告](BF-01_Report.md)：新契约123.714秒通过原180秒门槛，13/13相关回归通过，660行真实cell数据与来源hash归档。用户已解除旧档案缺失前置；旧H13 FAIL、历史缺失/timeout和全部物理gate保持。本轮只做BF-01；下一次允许仅执行BF-02，BF-02/H14/H15本轮未启动。

任务 ID / 标题：BF-01 / 从同一FirstJet publication派生CSF color物理体积平均。

目标：按BF-CSF V1 §3.2–3.6实现只读 `color_csf[K]`，正确处理有量纲distance proxy、输入标度与XY/RZ/AMR真实体积求积；不推进color、不改变sharp材料。

允许修改：精确production清单为新 `include/Modules/Hydrodynamics/CchRegularizedColor.h`、新 `src/Modules/Hydrodynamics/CchRegularizedColor.cpp` 和根 `CMakeLists.txt` 的唯一CCH注册，共3个。测试清单为新 `tests/contracts/hydrodynamics/cch_regularized_color_contract.cpp` 与 `tests/CMakeLists.txt` 注册，共2个。API/TODO/报告和方法状态索引；D:输出归 `runs/checks/BF-01/20261009_143534`。

禁止修改：Interface/Core数值实现、CchPhaseGeometry/MaterialTransport/DualMaterial、压力/鞍点/对流/黏性、旧H12/H13；不得改原CTest180秒门槛或物理阈值。

数学约束：`c=H_eps(-phi/|grad phi|)`，epsilon为明确的物理长度；力带内梯度须以本cell输入局部标度可辨识。正权物理平均，0<=c<=1，远离band纯相精确0/1，RZ乘2*pi*r一次。source/layout/mesh/full clock来自同一phase publication；不复制/重置phi，不clip、不用cell-center替代体积分，不独立输运、不修材料体积。

实现：新只读owner保留原phase及config、单元color/积分报告。使用原Q4 jet和原仿射包络只读查询；端点稳定的H/delta kernel；自适应正权Gauss、明确误差估计和工作上限、未解析薄带/退化/身份早拒绝。所有cell使用一致epsilon，报告band体积估计/包络、c与sharp alpha差异和积分成本。

验收：T0 Release无新warning、architecture/self_test/repository_layout与文件预算；T1独立kernel积分/导数/端点有界性、平面位置/标量缩放扫掠和真实体积解析平均、circle/sphere/radial near-axis、uniform/2:1 AMR/selected域、重复构造稳定与错误来源/clock/epsilon/退化拒绝；两求积预算/独立高阶积分收敛核验，旧相关压力/质量/相几何检查保持。原完整momentum旧候选timeout继续单列，不冒称已关闭。

失败规则：缺失/未解析数学条件明确FAIL或unsupported，不clipping、改phi距离、扩大epsilon、放宽误差或变更消费者换通过；同BF-ID持续修正实际错误。完整证据前不计PASS。

交付：真实color数据/epsilon、求积与输入标度报告、单一owner/API、现有回归/来源hash、未解决项及PASS/FAIL/BLOCKED。源码、临时、构建和数值输出全部位于D:。
