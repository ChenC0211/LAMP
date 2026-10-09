# H08 / 唯一中点材料流量与非守恒连续性缺陷

状态：**PASS**，2026-10-03；验收见[报告](H08_Material_Flux_Continuity_Defect_Report.md)。前置 H01–H07 PASS；依据[路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md)、冻结 CCH §4.2/5.1 的中点时间求积和[完全非守恒材料设计](H07_Pure_Gals_Hydro_Coupling_Design.md)。本任务不代表 Gate H。

## 目标

只解决一个数学事实：在同一固定 native layout/mesh 上，由开始/中点/结束的唯一 FirstJet 几何与明确同层面速度，生成唯一共享中点材料流量及其时间积分，逐单元显式保留非零材料连续性缺陷，不把几何质量覆写为守恒更新。

## 允许修改

- 新增 `include/Modules/Hydrodynamics/CchMaterialTransport.h`、`src/Modules/Hydrodynamics/CchMaterialTransport.cpp`，正式 CCH 材料 owner。
- root `CMakeLists.txt` 显式注册唯一 source owner。
- `tests/contracts/hydrodynamics/cch_material_transport_contract.cpp`、`tests/CMakeLists.txt`。
- 本 TODO、API、报告、现行路线图/接口索引与耦合设计；证据 `runs/baselines/H08_20261003`，临时工作 `runs/checks/H08/<时间戳>`。

预算：3 production/build +2 tests；实现前超过预算须拆任务。

## 禁止修改

H01–H07实现/旧测试、Interface状态/输运/几何/求积、Core/Solvers、旧Hydro/Compatibility、电场/App、cases、冻结PDF及历史TODO/报告。禁止独立材料库存、几何/权重修复、flux redistribution、中心/顶点判相、第二散度或插值 `rho_f*u_f`。H09 dual质量/动量、H10对流、H14完整corrector、H15动态regrid/边界留待对应任务。

## 数学约束

- `m_c^k=rho_l*V_l,c^k+rho_g*V_g,c^k`，k为开始/中点/结束；只读派生量。两相密度为有限正数，无density floor；本任务不支持真空。
- `F_rho,f^(mid)=rho_l*q_l,f^(mid)+rho_g*q_g,f^(mid)`，液气率均消费H07a同一对象。直接用独立气相率，不能用两次大数相减抹去小相。`I_rho,f=dt*F_rho,f^(mid)` 是具名中点近似，不冒称精确时空相通量。
- 同一真实native face只保存一次材料流量；cell外流使用保留D_I非零项的符号，不重装配关联矩阵。2:1使用实际fine subface。
- `delta_m=rho_l*(V_l,end-V_l,begin)+rho_g*(V_g,end-V_g,begin)`，`r_int=delta_m+sum(s_cf*I_rho,f)`，`R=r_int/dt`。不填零、不修改几何质量，非零缺陷本身不触发守恒性拒绝。
- cell相对缺陷为 `abs(r_int)/(abs(delta_m)+sum(abs(I_rho,f)))`；分母为零时要求分子也零并返回0。另保存有量纲R、signed/L1缺陷、几何质量及边界材料积分；总cell缺陷与总质量变化+边界积分满足同一代数收支。
- 开始层必须是上一步Initial/End，时间等于本区间begin、step=mid.step-1；中点和结束必须属于完全相同 `{step,begin,end}`，分别为Midpoint/End。三层共用同一layout及实际mesh；geometry身份各自保留。无主时钟、accepted commit或材料更新。
- 中点求积对光滑时间面率应二阶；当界面与整面重合造成phase rate时间跳跃时，不预先宣称局部二阶或零质量残差，须保留真实缺陷并测出该局限。完整sharp-interface时间收敛由H14/IH验收，不能靠本项短MMS替代。

## 实现要求

- [x] TODO先于源码；核验H07b交付并冻结源码/输入/PDF/程序身份。
- [x] 不可变prepared对象持有三层几何、中点相率、密度、完整区间、共享材料流量和逐cell诊断；复制共享只读数据。
- [x] 实际mesh/layout/完整clock与材料校验；finite arithmetic、工作尺寸及总收支检查；拒绝不污染已有对象。
- [x] XY/RZ、uniform/2:1、选中子域、纯相/平面/圆球、小液气相、等密度/1e3/1e6、正负/零速度、实际FirstJet候选与跨面非零缺陷检查。
- [x] 标准Release构建与新短契约、前置回归；来源/数值/链接/交付哈希检查和归档。

## 验收

**T0：** local Release无新增warning/error；cpp仅属正式CCH，依赖不进入旧Hydro/Compatibility/App；公开头自足，无kernel IO/字典/主时钟。help、architecture/self_test、repository_layout与前置App通过。

**T1：** 独立解析cell/face参考核验几何材料质量、中点面材料率、积分及逐cell缺陷，归一化误差≤1e-10；共享内面抵消/选中边界与2:1真实子面求和≤1e-12；总收支恒等式≤1e-12。高密度比不改变通量来源，等密度退化正确，小相材料贡献仍非零。人为非守恒候选及实际FirstJet跨面步骤的非零缺陷必须如实保留。身份/clock/layout/mesh/材料/溢出明确拒绝；重建与寿命保持。

**T2：** 非多项式圆/球geometry的 n=4/8/16 uniform/2:1 材料mass/face-rate参考误差递减且末档阶≥2。固定mesh光滑移动圆/球的T=.1、N=4/8/16中点材料积分与独立解析时间积分误差递减，末档阶≥1.8；均匀与2:1分别验证。平面整面crossing另测解析非零缺陷，不用该不光滑情形冒称二阶。检查只调用同一生产owner；参考不进入接受条件；不是完整PDE/高密度比稳定性证明。

## 失败规则

同源、代数收支、误用拒绝或预定光滑收敛未通过即FAIL，在允许范围修复，不放宽门槛或改密度/时窗。若需要改变Interface或超预算，先拆缺失contract；H08保持未关闭。只有真正无可推进的外部依赖才BLOCKED。

## 交付

正式材料流量/缺陷prepared对象、短契约、API/TODO/报告和原始数值/回归/来源哈希证据；明确H09/H10/H14/H15、CCH App、完整两相及低黏EHD仍未完成。总目标保持active。
