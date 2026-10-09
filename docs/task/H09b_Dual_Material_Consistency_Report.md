# H09b / 共享 dual 材料流量与非守恒动量一致性验收

日期：2026-10-03。**H09b PASS / H09 PASS；下一项H10。** 按[任务书](H09b_Dual_Material_Consistency_TODO.md)执行，与[H09a](H09a_Momentum_Dual_Geometry_Report.md)共同关闭H09；总目标、Gate H/IH/E/C 尚未完成。

## 实现与数学边界

[CchDualMaterial](../api/H09b_Cch_Dual_Material.md)在正式CCH库保留H08唯一材料流量及H09a三个真实几何质量对象，绑定实际publication identities、完整clock、同一native layout/D_I、dual geometry identity及密度。材料沿用用户授权的[完全非守恒路线](H07_Pure_Gals_Hydro_Coupling_Design.md)，不恢复Interface M、不增加材料库存或质量/几何修复。

primal面上的dual段按实际面积分配原流量；cell内部采用TODO预定的局部四块线性延拓，以真实四块体积作目标外流权重、内部面积/normal距离作metric，解固定gauge的3×3 SPD系统选定唯一内部率。它仅计算派生流量，不创建mesh或全局未知量，不改变原primal F；每条共享dual段只存一次。固定几何满足冻结方法§5.2的交换条件 `delta_u R_F F=R_m delta_c F`。RZ径向权重按真实r²测度；零面积轴段流量精确0。

真正的dual相积分质量一般不等于cell质量的几何加权平均，因此独立保留

\[
\Delta m_f=\rho_l\Delta V_{l,f}+\rho_g\Delta V_{g,f},\quad
r_f=\Delta m_f+\sum_g s_{fg}I_g,\quad R_f=r_f/\Delta t,
\]

\[
G_f=\Delta m_f-R_m\Delta m_c,\qquad r_f=R_m r_c+G_f.
\]

非零R和G均不填零、不拒绝或修复；相容积先作差，避免相减两份大质量。生产逐局部四块/dual检查交换关系、缺陷分解及总收支恒等式，门槛256*double epsilon，标度为gross rates/integrals；检查代数自洽，不要求物理材料守恒。两个速度方向分别报告三层质量、边界积分、signed/L1缺陷和最大R/相对缺陷。

同一中点速度的积分一致性项是 `u_mid*r_f`；常速度U的无外力材料更新满足

\[
m_{end}U=m_{begin}U-U\sum_g s_{fg}I_g+U r_f.
\]

材料动能项单列 `sum(0.5*u_mid²*r_f)`，不是完整动能预算或一般速度下的端点KE变化。返回load不是accepted state，内核无主时钟、commit、字典或IO；完整候选速度来源仍由后续App/corrector保证。

## 独立数值验收

| 检查 | 实测 / 预定门槛 |
| --- | --- |
| prepared场景 / segment检查 / dual检查 | 273 / 49,056 / 14,996 |
| 单面basis / 线性叠加 | 208 / 704；XY/RZ uniform/2:1 |
| affine / selected非均匀率 | 16 / 12组；全域、半域、孔洞、不连通域 |
| 独立native→dual交换最大误差 | 2.7756e−16 / 1e−12 |
| XY affine及RZ线性normal物理率解析误差 | 8.3267e−17 / 1e−12 |
| 解析dual质量变化/缺陷/一致性项最大归一化误差 | 1.3637e−15 / 1e−10；含独立segment积分核验 |
| 缺陷分解最大归一化误差 | 1.9062e−16 / 1e−12 |
| 无外力常速度恢复最大误差 | 3.5822e−12 / 1e−10 |
| 真实Mu进入H04：速度 / 制造压力误差 | 7.7521e−13 / 9.2974e−11；门槛1e−10 / 1e−8 |
| 实际FirstJet平移 / 真实Mu鞍点检查 | 24 / 12；平面、XY圆与RZ球，密度比1/1e3/1e6 |
| 明确拒绝 / 实际refine与coarsen | 20 / 4，XY/RZ各一次细化和粗化后续算 |

独立参考从native face实际相邻cell枚举半矩形和物理体积分数，直接以native adjacency散射primal F；dual外流从共享segment的negative/positive owner独立散射，不读取生产转移权重。解析材料参考使用独立圆/球裁剪原函数及平面积分；actual moving几何由正式FirstJet producer推进，非静态相体积重贴时钟。粗细交界球半径0.34与H09a实际cut覆盖相同。

实际平移中最大的 `abs(G)` 为27,488.90823，关闭uR后的最大绝对速度漂移为63,999.936；两项控制只统计actual FirstJet案例，不把合成fixture算进去。它们是该组单位/尺度下的数值，说明几何转移差异和动量一致性项不能省略。RZ 2:1平面、密度比1e6的两方向L1积分材料缺陷为196,349.3445，显式保留；本项没有用强制零残差达成PASS。

另有4组明确标注的合成常速度fixture：三个静止cut几何层配非零F，验证两个速度方向处理已知不相容材料数据，不将其冒称物理GALS benchmark。H04检查用真实end Mu作对角动量块，并加非零制造pressure load验证可消费性；实际无外力保持来自前述独立材料更新，不来自这项有制造载荷的参考solve。

光滑时间积分T=.3，N=4/8/16，纯相固定mesh、`u_x=sin(t)*x,u_y=cos(t)*y`（RZ径向为sin(t)*r），全部共享段累计积分对独立sin/cos原函数：

| 序列 | N=16相对L1积分误差 | 末档阶 |
| --- | ---: | ---: |
| XY uniform | 1.4648587705e−5 | 2.00004438 |
| XY 2:1 | 1.4648587705e−5 | 2.00004438 |
| RZ uniform | 1.4648587705e−5 | 2.00004438 |
| RZ 2:1 | 1.4648587705e−5 | 2.00004438 |

12行均误差递减，超过预定阶≥1.8。这是光滑面率的中点积分时间阶；不证明sharp crossing或完整NS的时间二阶。空间材料积分精度由H09a独立三档验证；本项不报告水平集曲率。

实际Interface ExactRefinement与FirstJetProjection粗化后，在新layout/来源上重建H07a/H08/H09a/H09b，并按原物理clock推进下一段；旧fine prepared在coarse mesh明确拒绝。测试参考是常速度，未实现一般Hydro velocity transfer。此次平面fixture两次相体积转移误差为0，不据此宣称一般粗化无误差。

拒绝覆盖密度、dual geometry、三层身份、layer位置、layout、保留D_I、完整区间/step、非法轴速度、溢出/下溢和mesh坐标/拓扑变化；原流量保持不变，复制共享只读data，原mesh析构后等效mesh仍可消费。下溢fixture按Windows long double=double实际表示范围检查。

## 工程验证与证据

预算3 production/build +2 tests：新增header/cpp及root CMake注册，新增短契约及tests CMake。前置1,034文件中1,032保持，仅两个CMake更新、新增3文件；312 case资产、31 fixture资产和原LAMP.exe/冻结PDF不变。新source唯一归 `lamp_hydrodynamics_cch`，只消费正式Core/Interface/线代；无旧Hydro/Compatibility/App依赖。H01–H09a实现和旧tests未修改。

- local Release定向构建通过，无新增编译warning/error。
- 最终新契约 **1/1 PASS，39.12秒**（测试本身39.11秒），包含加强后的selected非均匀率覆盖。
- 前置/受影响回归 **17/17 PASS，402.13秒**：H01–H09a、正式phase-region、FirstJet App、旧pressure-gradient load、help/architecture/self_test/repository_layout。H01–H09a选定数值行与H09a冻结证据逐字一致，App仍12smokes/49reject。
- 文档收尾后repository_layout另行通过，本地链接、来源/依赖/交付哈希核验归档；H09a历史证据保持。
- 未运行全部81项CTest、Debug或完整流体论文benchmark；本项没有接入CCH App。

初次编译暴露测试端既有publication API用法和auto不同类型声明错误，已修正。两次独立测试因净量抵消时误用净量标度失败，已按TODO预定gross rates/integrals归一化；生产256eps与独立1e−12门槛未放宽，失败原始日志保留。

原始证据：[verification](../../runs/baselines/H09b_20261003/verification.json)、[operator detail](../../runs/baselines/H09b_20261003/operator_detail.log)、[regression detail](../../runs/baselines/H09b_20261003/regression_detail.log)、[source audit](../../runs/baselines/H09b_20261003/source_audit.json)、[交付清单](../../runs/baselines/H09b_20261003/delivery_manifest.json)。TODO-before与失败日志保留；审计脚本关闭任务时归档至baseline/engineering_scripts。H09a历史TODO/报告、源码与全部26份证据保持。

复现新契约：`tools/build.cmd -Target quantum_cch_dual_material_contract`；构建后使用 `ctest --preset local -R "^cch_dual_material_contract$" --output-on-failure`（本机CMake已在PATH或使用其完整路径）。

## 未解决范围

H10一般KEP/symmetric convection与完整材料动能预算、H11受控耗散、H12–H13毛细载荷、H14完整变步长corrector、H15兼容AMR Hodge/一般速度transfer/生产边界、CCH App和两相无黏/低黏/EHD联合验收仍待完成。高密度比常速度通过只证明本项材料—动量一致性，不证明一般高密度比流动稳定。H09a区域查询逐次验证mesh的规模成本仍须在完整算例中测量。
