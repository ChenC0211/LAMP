# BF-03 / 真实 Hermite 曲率投影与一般 CSF 力

状态：**FAIL / 实施中**。本轮只执行BF-03；BF-01/BF-02封存PASS，BF-04/05/06和H14/H15不启动。

任务 ID / 标题：BF-03 / 同一FirstJet surface quadrature的局部正权曲率与实际CSF RHS。

目标：原publication实际κ_q、physicalWeight唯一来源；按固定局部正权核投影到真实native faces，构造 `S=sigma*kappa_f*C(color)`，核验常量再现、独立物理积分和真实误差随加密下降。

允许修改：既有 `include/Modules/Hydrodynamics/CchBalancedCsfLoad.h` / `src/Modules/Hydrodynamics/CchBalancedCsfLoad.cpp`；根 `CMakeLists.txt` 只给该源文件启用严格浮点/禁LTO边界，满足原Core exact expansion consumer约束，共3 production（仍在5上限内）。新增 `tests/contracts/hydrodynamics/cch_csf_curvature_contract.cpp` 和 `tests/CMakeLists.txt` 注册，2测试。API/TODO/报告、当前状态索引；D:输出 `runs/checks/BF-03/20261009_163129`。

禁止修改：Core/Interface/原actual source插值/Jacobian、BF-01、旧H12/H13、C/D/质量/F_rho/R/uR/KEP/黏性/鞍点、已有BF-02契约及原门槛。不加入生产App/时间推进、解析κ替换、曲率cap、远场补支撑、可调滤波宽度或拟合频率。

数学约束：原mesh/layout/C owner/color owner/geometry/full clock保持；仅调用唯一C得到g。固定物理cone核 `W_f(x)=max(0,1-dist(x,fCenter)/R_f)`，`R_f=epsilon+max(相邻真实leaf物理对角线)`，再乘selected域可见性指示：q在当前压力分量内，fCenter到q的闭线段必须完全在selected叶闭包内。不得跨hole/断连域借曲率，不另加硬ring裁断。

实现要求：原publication q物理权重只计一次，RZ不再乘r/2*pi；正权平均使用稳定normalized anchor-delta算术，常κ精确再现，归一化下溢不能静默丢掉非零均值。提供数学正权mean与精确segment-visible kernel，独立验证常量/有界性/正宽hole，均不能作为factory的几何/曲率覆盖输入。用物理空间桶和本地leaf覆盖检查降低无关扫描；ray参数复用原Core exact expansion有理端点，禁止浮点商相等跨hole。候选与visibility leaf查询均有独立上限/账本。暴露κ/supportCount/totalWeight/minmax/radius、boundary/axis/AMR/near-boundary与不可用状态。

力带与壁规则：所有有局部正支撑的非零g照算S，弱力不阈值归零，roundoff signal仅审计。缺支撑的有效free face或真实wall-band明确unsupported；边界邻cell已由BF-01证明transitionUpper=0且无支撑时，标纯C1边界反力项，variable工厂S=0且κ不适用，不能伪造远κ。无界面且全无band的纯相也物理零；BF-02显式constant工厂完整wall stress保持。sigma=0在来源/config合法后准确零并标投影未计算，不冒称κ=0已投影。

半径限制：phi/|grad phi| band不等同到Γ的epsilon邻域，固定核不能承诺所有FirstJet输入均覆盖band；实际缺支撑拒绝而不扩大R。固定物理epsilon的收敛目标为同hat连续积分；Sharp曲率收敛用epsilon=O(h)。边界/selected/轴截断如降阶须报告，不隐藏。

验收：T0 Release无新warning、architecture/self-test/layout、精确2/2范围，旧BF-01/BF-02及H12/H13相关回归保持。T1真实circle/sphere/ellipse/正弦XY/RZ、uniform/2:1 AMR，同pub两epsilon，三档n=8/16/32以上；原q独立权重平均核验到1e-12尺度，稳定正权kernel常量再现及bounds；真实κ/force独立连续曲线积分与误差下降/阶数完整报告（§4.2不承诺首版projection硬二阶，BF-05原>=1.8频率gate保持）。near-axis物理权重、纯相/σ0/constant channel、source/fullclock/regrid/跨孔/缺支撑/work-limit拒绝。保存所有source q与全部face κ/support/Cc/S及独立reference。

失败规则：任何来源变更、第二gradient/几何/库存、隐藏弱力过滤、错权重、非真实κ替换、跨hole/远field补支撑或误差不随加密下降均FAIL；同BF-ID修正实际实现，不改原gate。超出本任务的依赖明确BLOCKED，不能把静态projection称完整动力学。

交付：3/2源码范围与BF-02基线diff、唯一owner/API、实际构建/测试、曲率及力原始数据/收敛/边界异常与hash归档、PASS/FAIL/BLOCKED及下一次是否允许仅执行BF-04的明确判断。
