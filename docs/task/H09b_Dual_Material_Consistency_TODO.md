# H09b / 唯一 dual 材料流量与非守恒动量一致性

状态：**PASS**，2026-10-03，见[报告](H09b_Dual_Material_Consistency_Report.md)。前置H01–H08、H09a PASS，本项与H09a共同关闭H09，下一项H10；依据[路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md)、冻结CCH §5.2的流量交换条件及[完全非守恒材料设计](H07_Pure_Gals_Hydro_Coupling_Design.md)。不代表Gate H或完整两相时步。

## 目标

只解决一个数学事实：由H08唯一中点材料流量在线性几何转移下得到共享dual流量，保留真实三层dual质量及非零连续性缺陷，提供同层uR一致性项，使无外力常速度不受材料缺陷驱动。

## 允许修改

- 新增 `include/Modules/Hydrodynamics/CchDualMaterial.h`、`src/Modules/Hydrodynamics/CchDualMaterial.cpp`，root `CMakeLists.txt`唯一source注册。
- `tests/contracts/hydrodynamics/cch_dual_material_contract.cpp`、`tests/CMakeLists.txt`。
- 本TODO、API、报告、现行路线图/架构/文档索引/耦合设计；证据 `runs/baselines/H09b_20261003`，临时工作 `runs/checks/H09b/20261003_195159_436`。

预算3 production/build +2 tests；超出先拆任务。

## 禁止修改

H01–H09a实现/旧测试、Interface/Core/Solvers、旧Hydro/Compatibility、电场/App、cases、冻结PDF、历史TODO/报告。无独立材料库存、几何/体积修复、rho*u重新插值、第二primal散度、全局附加速度未知量或隐藏黏度。H10一般KEP对流、H14完整corrector、H15兼容Hodge/一般速度transfer/生产边界留在对应任务。

## 数学约束与固定局部转移

- 三个H09a质量对象绑定H08开始/中点/结束的实际publication identity、完整clock、相同layout/D_I、相同密度及同一dual geometry identity；不由本项重建Q4或拟合mass。
- 唯一primal材料rate为H08 `F_rho`。actual primal-face上的dual边段按真实面积比例分配原rate；RZ轴面积0要求原rate和分配均精确0。
- 每个实际cell仅在局部求值中用实际中心分成2×2四块，k=x+2y。这是固定的派生转移，不创建mesh、全局unknown或材料更新时间状态。边界原rate按实际子面/物理面积分配给四块，目标外流为 `w_k*delta_c F`，`w_k=V_k/V_c`。
- 四个内部率h=[X_lower,X_upper,Y_left,Y_right]，内部关联
  `B h=[X_lower+Y_left,-X_lower+Y_right,X_upper-Y_left,-X_upper-Y_right]`。
  设 `d=w*delta_c F-externalOutward`。固定参考g：x中线按径向体积分数theta插值两侧率（XY theta=1/2，RZ theta=(r_mid²-r_lo²)/(r_hi²-r_lo²)）；y中线为两侧率的1/2。
- 固定唯一转移 `h=g+W B^T lambda`，`B W B^T lambda=d-Bg`，固定lambda_3=0；W为内部段物理面积/相邻四块中心的normal距离。严格正几何权重，局部3×3 SPD Cholesky，失败明确拒绝。它选择距离g最近的流量延拓，不修改任何primal F、几何质量或缺陷；不是非守恒材料修复/接受后redistribution。
- H09a的cell内部共同段取对应内部半段率或其真面积部分；每个共享dual段仅计算一次。必须逐dual验证 `delta_u R_F F=R_m delta_c F`，`R_m(f,c)=V_fc/V_c`；交换误差门槛256*double epsilon，scaled by gross transferred rates。
- `delta_m_dual=rho_l*(V_l,end-V_l,begin)+rho_g*(V_g,end-V_g,begin)`，相容积先作差；`r_dual=delta_m_dual+sum(s*dt*F_dual)`，`R_dual=r_dual/dt`。非零缺陷不拒绝、不填零。
- 显式保留 `geometryTransferChange=delta_m_dual-R_m delta_m_cell`，检验 `r_dual=R_m*r_cell+geometryTransferChange`。不能用Rm平均cell mass或cell缺陷覆盖真实dual量。
- 同层动量一致性积分为 `u_mid*r_dual`；材料动能诊断为 `sum(0.5*u_mid²*r_dual)`，分别核算，不宣称严格动能守恒。常速度U的无外力更新为 `m_end*U=m_begin*U-U*outwardIntegral+U*r_dual`。
- 每方向signed/L1缺陷、最大R/相对缺陷、边界材料积分及代数总收支显式报告；相对缺陷分母为 `abs(delta_m)+sum(abs(I_dual))`，全零返回0。保留实际来源、完整时间区间与dt，不拥有主时钟/commit。

## 实现要求

- [x] TODO先于源码；核验H09a全部交付并冻结源码/输入/PDF/程序。
- [x] 不可变prepared对象持有H08与三层H09a质量；局部线性转移、共享积分、缺陷/转移差异、同层uR和材料动能诊断；finite/身份/clock/mesh错误拒绝。
- [x] 独立native邻接与dual矩形边积分核验，单面basis/随机率覆盖实际2:1和边界；XY/RZ、全域/半域/孔洞/不连通域。
- [x] 实际FirstJet平面/圆球平移、密度比1/1e3/1e6、非零缺陷、关闭uR对照和常速度保持；真实Mu进入H04单体参考块。
- [x] 实际Interface exact-refine及projection-coarsen后重建来源并续算；旧prepared拒绝，新来源/物理clock/误差核验。一般流体velocity transfer仍归H15。
- [x] 标准Release构建、新契约/前置回归、哈希/来源/依赖/数值审计与归档。

## 验收

**T0：** cpp唯一CCH owner、无旧Hydro/Compatibility/App依赖，公开头自足，无kernel IO/字典/主时钟；local Release无新增warning/error，help/架构/目录/App通过。

**T1：** local四块和全局dual交换、共享段/边界总收支/单面basis线性叠加≤1e−12，生产交换检查更严格按256eps；XY affine normal flux与RZ `u_r=a*r,u_z=b+c*z`内部真测度率解析再现≤1e−12。真实几何dual mass/变化/缺陷/几何转移差异独立解析核验≤1e−10；实际非零缺陷与geometryTransferChange保留。同层常速度恢复≤1e−10、关闭uR在预设crossing fixture产生>1e−5误差；真实材料Mu/H04速度误差≤1e−10、制造压力误差≤1e−8。来源/密度/layout/D_I/dual geometry/clock/overflow/mesh错误明确拒绝，不污染原对象。

**T2：** 固定mesh纯相光滑时间率 `u_x=sin(t)*x,u_y=cos(t)*y`（RZ径向为sin(t)*r），T=.3，N=4/8/16，全部dual边率时间积分与独立sin/cos原函数误差递减，末档阶≥1.8，XY/RZ uniform/2:1分别验证。actual FirstJet圆/球和平面跨面、真实exact-refine/projection-coarsen前后来源及常速度恢复通过。该项与H09a几何空间阶共同支持H09算子验收；不以光滑率时间阶冒称任意sharp crossing或完整NS二阶，完整稳定性仍归H10/H14/IH。

## 失败规则

同源、交换、预定收敛、误用拒绝或常速度控制未通过即FAIL，在本范围修复，不放宽门槛/密度/时窗。需要改Interface/H09a或超预算先拆缺失contract，H09保持未关闭。仅真正无可推进的外部依赖才BLOCKED。

## 交付

共享dual材料流量/积分、真实非零dual缺陷/转移差异、同层uR/材料动能诊断、独立短契约与API/TODO/报告及原始数值/哈希证据。只有本项与H09a均通过才关闭H09、进入H10；总目标保持active。
