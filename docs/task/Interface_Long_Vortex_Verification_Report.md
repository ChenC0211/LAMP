# Interface / 可逆差速旋涡长时间验证报告

日期：2026-10-03。**两个完整窗口均 PASS。** 此报告关闭此前用户提出的长算例请求，并回答当前曲率表现；数据作为 H 阶段接线前的 Interface 基线，不属于流体或 EHD 验收，也不新增 H01 之外的模块算法任务。

## 配置与参考

当前生产方法为纯 Hermite–GALS / FirstJet，体积只作派生诊断。两组输入位于 [长涡旋算例](../../cases/hermite_gals_vortex_long_time/README.md)，生产代码没有专用算例分支。单位 XY 域中，初态圆心 (0.65,0.53)、半径 0.13，规定速度为

\[
(u,v)=40\cos(\pi t)((x-0.5)^2+(y-0.5)^2)(-(y-0.5),x-0.5).
\]

粒子相对旋涡中心的半径不变，转角为 (40/π)sin(πt)ρ²；映射面积 Jacobian 为 1，准确液相面积始终为 π·0.13²。t=0.5/1.5 变形最大，t=1/2 回归初始圆；窗口覆盖一个完整速度周期和两次反向变形—回归过程。

两组均从中央 level-1 加密的 16×16 基础网格开始，dt=0.001、2,000 步、每 25 步输出，共 81 帧。dynamic_amr 有 15 次成功网格切换，叶单元数在 256/448 间变化；fixed_amr 保持 448 个叶单元。两组使用同一 LAMP.exe，SHA-256 为 `d67fda82f4bb33cce27ba212ca9c94fac93d207f9af6087b96629f04bdd5819b`。

[独立后处理器](../../tools/postprocess/analyze_interface_reference.py)使用逆映射的二阶链式导数计算 level-set 法向和曲率，并以 forward 参数曲线及有限差分交叉检查。参考自检最大曲率差 1.67e-13、法向差 3.78e-15、面积 Jacobian 差 4.44e-16。参考只读取保存结果，不反馈生产状态。

reference/verification.json 的五输入哈希与阈值在数值运行前确定；两组均通过来源、完整时钟、81 帧、重网格次数、拒绝状态、面积、位置、法向和曲率检查。曲率阈值原为相对 L2≤0.1，未在看到结果后收紧或放宽。

## 曲率与长时间误差

| 指标 | dynamic_amr | fixed_amr |
| --- | ---: | ---: |
| 接受步 / 成功 regrid | 2,000 / 15 | 2,000 / 0 |
| 曲率相对 L2 峰值 | **0.174186%**，t=1.700 | **0.102961%**，t=1.750 |
| t=1 回归圆的曲率相对 L2 | 0.0580188% | 0.0253185% |
| t=2 末态曲率相对 L2 | **0.0483109%** | **0.0365724%** |
| 所有保存帧的曲率最大绝对误差 | 0.123442 | 0.0728049 |
| 相对面积误差峰值 | 6.55891e-5 | 3.11725e-5 |
| 法向 L2 峰值 | 1.82036e-4 | 8.94924e-5 |
| 法向位置缺陷最大值 / R | 1.90339e-4 | 1.06534e-4 |
| 采样折线 Hausdorff 最大值 / R | 3.27317e-3 | 2.61967e-3 |
| 曲率不可用样本数（所有接受态） | 0 | 0 |
| 启动器完整运行耗时 | 639.55 s | 745.68 s |
| 接受步耗时 median / p95 | 0.3254 / 0.4161 s | 0.3734 / 0.4078 s |

曲率误差定义为

\[
E_\kappa=\left(\frac{\sum_i w_i(\kappa_h(x_i)-\kappa_{exact}(x_i,t))^2}{\sum_i w_i\kappa_{exact}(x_i,t)^2}\right)^{1/2},
\]

其中 x_i 是保存 VTP 的数值界面采样点，w_i 为相邻折线段长度各半的权重，参考是该位置的准确 level-set 曲率。这是弧长加权采样误差，包含数值形状位置的影响；不宣称完整隐式几何的连续范数，也不是最大逐点相对误差。绝对曲率误差与相对 L2 峰值可能发生在不同时间。误差按每 25 步保存的帧检查；步间仅有同源 kappa_max/rms 诊断，不能把 81 帧峰值称为全部 2,001 个状态的精度峰值。

初始圆准确曲率为 1/R=7.69230769；动态运行的同源 kappa_max 峰值约 39.97，覆盖约五倍曲率增大。初始 Q2 圆的曲率接近机器精度是多项式再现结果，不代表一般变形也达到机器精度。两组形状与曲率保持良好，但 dynamic 的粗化时段和 transfer 都会影响误差；对照差值不能全部归因于 AMR transfer。

面积误差覆盖全部 2,001 个接受态。位置缺陷是 |phi_exact|/|grad(phi_exact)|，并非准确欧氏距离；Hausdorff 是独立参考采样到保存折线的双向近似，包含折线输出离散误差。固定 AMR 此次耗时更长与其全程使用较多单元相符，单次运行不足以断言一般速度收益。

## 可追溯结果与使用

| 数据 | dynamic_amr | fixed_amr |
| --- | --- | --- |
| 程序状态 / 日志 | [run.json](../../runs/hermite_gals_vortex_long_time__dynamic_amr/20261003_113027_852_validation_11d7700e/run.json) / [run.log](../../runs/hermite_gals_vortex_long_time__dynamic_amr/20261003_113027_852_validation_11d7700e/run.log) | [run.json](../../runs/hermite_gals_vortex_long_time__fixed_amr/20261003_114341_057_control_43d9ffb7/run.json) / [run.log](../../runs/hermite_gals_vortex_long_time__fixed_amr/20261003_114341_057_control_43d9ffb7/run.log) |
| 独立验收与来源哈希 | [verification.json](../../runs/hermite_gals_vortex_long_time__dynamic_amr/20261003_113027_852_validation_11d7700e/output/analysis/reference_20261003_034203_415786/verification.json) | [verification.json](../../runs/hermite_gals_vortex_long_time__fixed_amr/20261003_114341_057_control_43d9ffb7/output/analysis/reference_20261003_035930_638769/verification.json) |
| 误差原始数据 | [accuracy_history.csv](../../runs/hermite_gals_vortex_long_time__dynamic_amr/20261003_113027_852_validation_11d7700e/output/analysis/reference_20261003_034203_415786/accuracy_history.csv) | [accuracy_history.csv](../../runs/hermite_gals_vortex_long_time__fixed_amr/20261003_114341_057_control_43d9ffb7/output/analysis/reference_20261003_035930_638769/accuracy_history.csv) |
| 形状比较图 | [PNG](../../runs/hermite_gals_vortex_long_time__dynamic_amr/20261003_113027_852_validation_11d7700e/output/analysis/reference_20261003_034203_415786/interface_evolution.png) / [SVG](../../runs/hermite_gals_vortex_long_time__dynamic_amr/20261003_113027_852_validation_11d7700e/output/analysis/reference_20261003_034203_415786/interface_evolution.svg) | [PNG](../../runs/hermite_gals_vortex_long_time__fixed_amr/20261003_114341_057_control_43d9ffb7/output/analysis/reference_20261003_035930_638769/interface_evolution.png) / [SVG](../../runs/hermite_gals_vortex_long_time__fixed_amr/20261003_114341_057_control_43d9ffb7/output/analysis/reference_20261003_035930_638769/interface_evolution.svg) |
| 误差历程图 | [PNG](../../runs/hermite_gals_vortex_long_time__dynamic_amr/20261003_113027_852_validation_11d7700e/output/analysis/reference_20261003_034203_415786/error_history.png) / [SVG](../../runs/hermite_gals_vortex_long_time__dynamic_amr/20261003_113027_852_validation_11d7700e/output/analysis/reference_20261003_034203_415786/error_history.svg) | [PNG](../../runs/hermite_gals_vortex_long_time__fixed_amr/20261003_114341_057_control_43d9ffb7/output/analysis/reference_20261003_035930_638769/error_history.png) / [SVG](../../runs/hermite_gals_vortex_long_time__fixed_amr/20261003_114341_057_control_43d9ffb7/output/analysis/reference_20261003_035930_638769/error_history.svg) |

两组 PNG 形状图与误差历程图已目视检查，保存界面与参考曲线对应，标签/单位/坐标可读。ParaView 打开各运行的 output/fields/solution.pvd。

```powershell
.\tools\run_case.cmd hermite_gals_vortex_long_time/dynamic_amr
.\tools\run_case.cmd hermite_gals_vortex_long_time/fixed_amr
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/postprocess/analyze_interface_reference.ps1 -Run runs/<算例>/<时间戳>
```

历史独立空间研究的光滑曲率阶为 **2.44–3.41**，见 [Gate I 最终审计](I12p3r2_FirstJet_Gate_Final_Report.md)与[收敛复核](I12p3r1_FirstJet_Gate_Reaudit_Report.md)。本次没有重跑空间/时间阶扫描，当前 long-window 不能单独给出阶数。

结论范围为规定光滑速度下的光滑曲线非线性变形、回归与重复 AMR。没有拓扑断裂、真实细丝颈缩、表面张力反馈、Hydro/EHD 或低黏流体动力学，因此仍需 H 阶段及其对应耦合 gate。
