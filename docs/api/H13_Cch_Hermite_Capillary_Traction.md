# H13 同源 Hermite 毛细载荷候选

最新045122：**H13 FAIL，当前API不变。** [原M条件证据](../task/H13_Phase_Cut_Kinetic_Metric_Certificate.md)和[待授权联合方法](../task/H13_Phase_Cut_Kinetic_Hodge_Method_Change.md)说明，常矩阵/伴随配对通过不能证明准确高密度响应；phase-cut物理动能Hodge及同源分相source尚未获准或接入。H14未启动。

状态：2026-10-08，**H13 FAIL，授权联合改造的开发候选**。当前接口要求同一 prepared source；局部完整弱牵引已实现，完整物理验收待完成。证据见[本次任务书](../task/H13_Weak_Traction_Acceptance_TODO.md)及[报告](../task/H13_Hermite_Capillary_Traction_Report.md)。H12 显式常应力工厂保留。

2026-10-09 / 035647：当前API/算子未修改。[精确平面直接反例](../task/H13_Exact_Plane_Cut_Response_Counterexample.md)已确认一般cut normal约一阶及高密度回复方向错误；下列完整矩阵/局部formula性质不证明准确物理响应，该候选不能采用。

## 接口与所有权

```cpp
const auto pressure = CchPressureLoad::build(phase.metrics().divergence());
const auto normal = CchInterfaceNormalLoad::capillary(
    mesh, pressure, phase, sigma, reconstruction, &boundary);
normal.validateFor(mesh, pressure, phase.geometryIdentity(), phase.timeLayer(),
                   reconstruction, &boundary);
auto rhs = normal.balanceLoads();
for (double& value : rhs) value = -value; // A*u + C*p = rhs
```

`CchInterfaceNormalLoad` 仍为唯一法向应力沉积 owner，属于 `lamp_hydrodynamics_cch`。工厂保留原 pressure/phase owner，要求相同 layout、同一 D_I、正确 mesh/source/full clock。复制共享不可变数据；不建立独立材料库存、几何、梯度、时间推进或接受状态。

只读入口包括 `isCapillary()`、`surfaceTension()`、`curvatureAverages()`、`supportMeasures()`、`stressMatrix()`、`reconstruction()`、`boundary()` 和 `constantReactions()`。毛细对象保留相同 reconstruction/boundary owner，新的验证重载拒绝不同准备或边界身份。`constantJump()` 仅用于显式常应力对象。`sigma` 必须有限且非负，零值使载荷精确消失。

## 当前离散与符号

FirstJet 的 phi<0 为液相，曲率为液→气外法向的 `+div(n)`；冻结 CCH §2 使用相反法向。因此采用 `J=-sigma*kappa_FirstJet`。

令 `a=C*alpha`，`b_qu` 为同一 prepared source 在实际界面求积点的法向虚速度权重，`w_q` 为真实正物理测度，`J_q=-sigma*kappa_q`。当前候选为

`barJ_u = sum_q(w_q*b_qu^2*J_q) / sum_q(w_q*b_qu^2)`，

`F_u = barJ_u*a_u - sum_q w_q*b_qu*(J_q-barJ_u)`。

均值以一个应力锚点的差值稳定计算；无局部虚速度支撑的行使用常量锚点，getter 的支撑曲率/测度为零。已声明不可穿透边界的虚法向自由度被排除，其常量压力反力单列。该常量通道保留原完整压力反力，不能把无支撑行全部删去。

`C*alpha` 只调用已有 pressure operator。support 为真实 source 的局部稀疏足迹，来自 `cellStencil/cellJetWeights` 的同一 Q4/master/subface algebra；不作全域 unit probes 或运行时拟合。ds / 2*pi*r*ds、曲率与法向均来自同一 publication。

完整稀疏 `J_f` 将上述既定弱力延拓到 native 行列，满足 `F=J_f*a`。每个点的局部线性 functional 满足 `r_q*a=1`，仅用于因子分解，不是另一载荷或质量逆。足迹中的 incidence 全零时明确拒绝。常应力差值精确为零，因此完整矩阵为 `J_f(J0)=J0*I`，包括 inactive rows；并非只验证在 `C*alpha` 上的作用。H12 工厂仍直接返回 `C(J0*alpha)`。

## 已知限制

当前局部矩阵、P1/P2/jump 虚功与常量平衡契约不自动证明实际 FirstJet 几何运动的能量相容性。原静滴72几何/1728解、完整无黏周期、空间/时间/AMR/低黏及全部受影响 source gate 均须按新身份验证；不得通过曲率滤波、解析替换或最小黏度掩盖失败。

以下日期条目记录历史候选与隔离研究。其中“未改变本 API/生产载荷”等表述仅对应各历史检查点；当前实现由上述授权联合候选说明。

冻结§9.1没有把支撑限定为上述两个cell，§9.2也未明定所有J_f必须对角；当前实现的反例不能升级为对所有允许映射的不可能性证明。扩大支撑的固定正二次B样条隔离参考恢复了约二阶频率，但实际FirstJet表面能仍不匹配，因此没有改变本API或生产载荷。详见[H13报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-06--084708-接续载荷支撑范围核验)。

下游尚无 H14 完整耦合时间推进；研究驱动的毛细子系统 split 不是生产 App。H13定位中已复现直接破坏H08a的运动Jacobian错误，限定修正及验证范围见[反例说明](../task/H13_H08a_Transparent_Picard_Jacobian_Counterexample.md)。毛细候选接口及R/uR、KEP、H12门槛保持；新版本完整材料28窗口及原空间/时间/R历史回归已通过，这一处限定修复已收口。

2026-10-06后续只读审计给出当前对角C*alpha通道的运动范围缺口。隔离的功配对参考虽改善能量，但实际周期频率和四档线性化阶数未达标，未采用；本API公式和生产数值实现没有因此改变。不可把该研究参考当成第二个可选生产工厂。


2026-10-06 / 093237：由冻结系数局部动能分析得到的被动载荷响应仅在runs/checks隔离验证；单一自然模态频率改善，完整初值响应及实际能量仍需验收。当前只定义平面均匀等密度参考，未形成一般曲界面/AMR方法，未增加生产API，也没有改变本工厂的公式、质量或common source。见[本轮报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-06--093237-接续最终载荷响应的隔离验证)。

2026-10-06 / 104140：只读平移检查拒绝固定反馈的通用性；Green核与最小native动能分析仍属隔离载荷设计诊断，未采用P1反馈或全局投影。本生产工厂、原质量、source与H08a均未变，不能将矩阵参考的改善当成毛细生产能力。见[本轮报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-06--104140-接续位置依赖与多余动能诊断)。

2026-10-06 / 111706：旧n64实际短窗已完成并通过原始审计。局部native环流验证了常应力平衡和正响应的代数可能性；精确体积消元仍只作平面等密度参考，n128完整初值响应改善但一般曲界面/AMR/黏性定义未完成。n256尺度复用被原查询门槛拒绝，未采用参考结果，也未重开r11。本API、原native质量和材料耦合均不变。见[最新报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-06--111706-接续局部环流与精确体积消元)。

2026-10-06 / 120753：斜界面实际endpoint更新证实，完整几何运动包含FirstJet primary、mixed closure和中心重构的导数，不能未经核验以公开速度trace替代。两方向三档真实更新及表面能二次项与该导数相符；这是载荷设计依据，尚未接入新的生产工厂。原对齐gamma不能直接推广，r11/H08a均保持收口，本API仍为未通过完整H13的原候选。见[最新报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-06--120753-接续真实firstjet几何运动的载荷配对对象)。

2026-10-06 / 124031：闭合曲界面XY/RZ、uniform/AMR的真实几何方向导数及固定分支载荷转置已核验。分片C1旧scalar使精确方向导数可能依赖实际回溯分支；本组该差随加密减小，不能把固定分支B宣称为唯一速度无关的生产算子。常应力平衡配对后保留平均曲率乘几何体积变化率的功余项；不强制其归零，也不新增来源或精确功硬约束。小网格参考未接入本API，原H13动态失败保持。见[最新报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-06--124031-接续曲界面方向导数与保留体积功余项)。

2026-10-06 / 130703：从原native正交象限的旧scalar Hessian对称极限构造了速度无关的几何测试映射B_c，XY/RZ及AMR观察再现真实方向导数的奇部分；未更改实际FirstJet，偶部分余项保留。B_c转置在原高密度静滴例外位置的144次实际求解中虽有空间衰减，却增大纯气相响应，因此仍未接入本工厂。必须联合验证原材料惯性/压力、几何运动和最终载荷，不能单凭功配对采用。见[最新报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-06--130703-接续对称几何测试映射与高密度静滴对照)。

2026-10-06 / 132735：隔离参考只用原native质量归一化缩放B_c转置的变量曲率部分，保留原常应力；240次真实静滴求解显示已测高密度气相伪流降低，24组平界面谱核验完整归一化响应。不过等密度自然频率末档仍1.17阶、曲界面真实几何功差保留，未接入本工厂。所有原门槛和H13 FAIL状态保持，H08a/r11不重开。见[报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-06--132735-接续变量载荷密度缩放与完整平界面响应)。

2026-10-06 / 135228：原质量比例D直接缩放变量几何转置的参考，在四分之一网格偏移的局部平界面法向测试中产生不随加密消失的18.54%偏差，已排除其作为通用载荷定义。真实椭圆/AMR弱功数据同时保留；14个观察均结束，未接入本API。继续同源应力积分与原压力/材料惯性/几何运动的联合响应；此反例不触发r11/H08a重开。见[报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-06--135228-接续密度缩放载荷的偏移反例)。

2026-10-06 / 211858：两侧独立外推载荷参考在18组最终回复谱中出现负/复增长，三档原native鞍点实际求解与公开查询复现代表性负分支，已拒绝。该参考未接入本API；原工厂、质量和common来源不变。继续联合压力/惯性/几何/牵引响应，保留可收敛功余项，不新增精确功硬约束。见[报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-06--211858-接续两侧载荷外推的native增长反例)。

2026-10-06 / 215208：原native质量逆的Gram配对和原Green首项导出的系数在隔离平面参考中消除增长，原完整初值与能量有所改善；等密度末档仍1.73986阶，未到原1.8门槛，未接入本API。全部原分支、初值、能量和初次诊断退出保持；一般曲面/AMR/效率及实际动态验收未完成。生产来源及r11/H08a不变。见[报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-06--215208-接续正的联合native响应与green首阶宽度)。

2026-10-07 / 224337：隔离响应的原Dirichlet边界表达修正保持完整界面回复逐位不变，等密度频率仍1.73986阶失败。完整速度历史及纯气相native环流显示原对角惯性与连续source液相支撑的失配；下一步仅作同源分相场动能的只读比较，未改变本工厂、生产M或source，未新增壁面/功硬约束。此证据不直接阻断H08a，r11/H08a不重开。见[报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-07--224337-接续牵引弱功壁面回复与气相惯性诊断)。

2026-10-07 / 091136：公开source积分量化了原气相环流的液相体动能；六个观察及18组原完整初值的被动场能量历史完成，原native能量定义不变。有效牵引环流功显示联合Gram参考在远界面壁面附近仍有载荷，下一步只检查辅助回复的局部表达；未恢复非对角生产质量，未改本API/source/F/R/uR/KEP，也未新增零功或场能量硬约束。H13 FAIL，r11/H08a不重开。见[报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-07--091136-接续同源分相体动能与有效牵引支撑)。

2026-10-07 / 093259：辅助Gram局部化保持完整界面回复逐位不变，消去远壁面环流功，却未解决原频率/气相速度失败；固定该辅助族的代数质量探针显示其有效载荷系数仍依赖全局P。停止该族的生产推广，回到冻结允许的局部应力映射；本工厂、原生产M/source、H12常应力和H13门槛均保持。初始非判别对称控制保留且原结论作废，不作为H08a反例。见[报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-07--093259-接续辅助gram的几何支撑局部化)。

2026-10-07 / 094907：重新视觉核对冻结J_f(J0)=J0*I及同层对角M，简单加性候选等密度1.557阶失败。三个实际publication/本工厂观察确认当前单元常权平均丢失native可感的应力分量；三个正权的局部应力表示只作测量，未接入本工厂。下一步先保留应力变化再构造压力相容的完整常量矩阵；原native M/C、source、F/R/uR/KEP和H13门槛保持，r11/H08a不重开。见[报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-07--094907-接续从载荷系数回到应力数据)。

2026-10-07 / 130155：局部应力偏差参考完成平面常矩阵条件，但等密度0.953阶及复增长失败；实际原native鞍点/public查询重放了代表性复耦合。该参考未接入本工厂，下一步仅隔离检验由同源局部应力构造的已知压力跳跃lifting。原M/C/source/F/R/uR/KEP及门槛保持，r11/H08a不重开。见[报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-07--130155-接续局部应力偏差映射的最终算子拒绝)。

2026-10-07 / 133530：已知压力跳跃lifting保持平面常量矩阵及原C，但等密度负刚度、0.716阶频率和增大的原能量缺陷失败；原native求解及source应力加权散度分解确认机制。固定光滑弱功收敛不能覆盖这一失败，未接入本API。下一步仅诊断局部二相压力传输的最终native/FirstJet耦合；本生产工厂、M/C/source/F/R/uR/KEP及门槛保持，r11/H08a不重开。见[报告](../task/H13_Hermite_Capillary_Traction_Report.md#2026-10-07--133530-接续已知压力跳跃-lifting-的负刚度机制)。
