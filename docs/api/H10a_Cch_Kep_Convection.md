# H10a / 共同位置中心对流与非守恒材料动能

状态：2026-10-03，**H10a PASS**，验收见[TODO](../task/H10a_Compatible_Center_Convection_TODO.md)及[报告](../task/H10a_Compatible_Center_Convection_Report.md)。数学装配详见实现前的[H10设计](../task/H10_Compatible_Center_Convection_Design.md)。H10a只验收算子，下一项H10b完整无黏benchmark；H10未关闭。

**当前状态（2026-10-06）：** 默认中心算子的规定空间/KE 范围由 [H10o](../task/H10o_Bounded_Final_Native_Convection_Report.md) 与原 H10b 验证；H08a 时空积分也已在规定范围通过。[H11](H11_Cch_Controlled_Dissipation.md) 已通过规定的局部非负耗散、材料动能与光滑收敛验证。下一单任务 H12，本轮未启动。 原 `build` 的默认中心路径保留；材料仍仅唯一共享 F、内部反号/AMR 子面和、来源一致、KEP/skew 四硬约束。当前总状态以[路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md)为准；以下 H10a 部分保留其历史 API 范围，完整流体 App、一般物理边界及全耦合二阶未验收。

公开头：[CchKepConvection.h](../../include/Modules/Hydrodynamics/CchKepConvection.h)，唯一owner为正式 `lamp_hydrodynamics_cch`。

```cpp
auto convection = CchKepConvection::build(mesh, actualDualMaterial);
auto boundary = CchConvectionBoundary::bind(mesh, actualDualMaterial,
    actualSegmentTraces, actualSources, actualMidpointClock);
auto result = convection.apply(mesh, actualMidpointVelocity, boundary,
    actualSources, actualMidpointClock);
const auto& N = result.momentumRates();
const auto& materialMomentumLoad = result.momentumLoads(); // -N + uR
const auto& budget = result.report();
```

## 来源与所有权

prepared保留H09同一不可变材料owner，包含三层真实几何质量、唯一F、实际三层identity及完整中点clock。重构只由实际几何生成，无质量更新时间、独立速度primary、字典/IO/主时钟或commit。复制共享只读data；result持有原材料和实际bound boundary。

`apply`核验mesh、实际来源、完整区间与速度layout/轴值。boundary还核验同一prepared材料owner（共用data的copy有效，重新prepare即便相同geometry/clock也拒绝），不能拿另一套F的边界冒充当前候选。裸CchNativeValues没有自己的物理clock；实际velocity与H08生成率的同层物理关系仍须由后续App/corrector保证，不能靠此API验证任意两组外部数组的物理来源。

## 重构和相容载荷

`points()`、`weights()`与H09共同segment顺序相同。内部段每份weight仅引用当前同分量native velocity；外部段weight为空。位置为实际面积测度中心，RZ水平段用径向一阶矩。支持来自两个owner的实际交集cell及一环邻接，不跨不连通区域；SVD最小改变满足常数/x/y的affine约束。可识别rank不足仍检查全部约束，不能识别目标的位置则明确拒绝。单列隔离RZ slab缺少测度中心所需径向信息是已测试拒绝，两个采样列的slab正例通过。

在native采样中点与共同位置相同的XY段上，weights精确保留两个0.5；uniform mesh外边界的部分半dual段也需要重构，不能将所有uniform段称为对齐段。RZ近轴无需新的primary反射ghost；使用实际支持、真实位置与H05轴速度约束，零面积轴材料率精确0。

每条内部段按同一trace与其代数转置配对：

\[
N_g=\tfrac12 F_g[b_g(t_g^Tu)-t_g(b_g^Tu)+\operatorname{diag}(b_g)u].
\]

`momentumRates()`是正向外流N，不是已除质量的速度导数。该装配同时保持内部动量抵消、常速度一致性和native动能关系；无sensor/limiter、人工黏度或upwind。`centerMomentumRates()`只是原始 `F*t^T*u`（边界为F*u_b）的诊断，不能将其直接散射替代最终相容载荷；共同位置重构后的raw散射在实际2:1处会破坏KE恒等式，这有独立负对照。

## 显式边界和预算

`segmentTraces`长度必须等于所有segment数量；外部段显式给transported-component trace，内部未使用位置必须0，所有值finite。径向轴外部段的径向trace必须0；零面积段不除面积。没有缺省owner trace、ghost状态或隐式物理边界模式。具体wall/traction/open/backflow和生产周期支持仍归H15，研究driver可在该数据端口显式提供解析/周期trace。

`momentumLoads()`已经包含H09同层uR：`-N+u_mid*integratedDefect/dt`，消费时不能再加第二次uR。它不包含压力、黏性、毛细/电场或其他外力，也不推进accepted velocity。

每分量的report保存：N总和/边界动量率、`u^T N`、边界动能功 `sum(sF*(u_i*u_b-0.5*u_i²))`、材料项 `sum(0.5*u_i²*R_i)`及半离散动能率 `sum(u_i*load_i-0.5*u_i²*delta_m_i/dt)`。动量收支和 `kineticRate=materialKineticPower-boundaryKineticPower` 独立核验，生产256eps且按组装前gross contributions/powers标度归一化；非零R不修复。这不是全离散端点KE变化或一般非守恒材料的严格动能保持。

[H08a2r15](../task/H08a2r15_Closed_Material_Energy_Reference_Report.md)已固化消费者参考：以真实 dual 四侧的唯一共享 F 构造 q、以几何相体积差构造 Δm，独立核验 `R=Δm/dt+q`。闭域静态材料中 `u^T N=0.5*sum(q_i*u_i²)`；native `D_I=0` 不推出所有 dual q/R 为零，也不推出这项对流功为零。墙面零通量和 KEP/skew 仍为硬约束，旧零功数值只作诊断。当前旧来源和隔离正则候选的完整 KEP 契约均通过；正则来源尚未正式部署。

非有限计算、非零乘积下溢、不可识别支持、foreign owner/layout/来源/clock、非法轴和mesh坐标/拓扑变更均明确拒绝；已有材料/weights/result保持不变。计算先乘F再乘第二份速度，并先分别乘1/2再求和/差，避免可避免的u²或平均值中间溢出；最终量或gross诊断仍不可表示时拒绝。

regrid后先做实际Interface transfer，再在新mesh/来源上重建H07a/H08/H09a/H09b及本算子；旧prepared不能重贴标签。H10a测试常速度续算，不实现一般Hydro速度transfer。

复现：`tools/build.cmd -Target quantum_cch_kep_convection_contract`；`ctest --preset local -R "^cch_kep_convection_contract$" --output-on-failure`（CMake在PATH或使用本机完整路径）。完整空间/时间收敛与无一阶耗散结论只能由H10b研究证据支持。
