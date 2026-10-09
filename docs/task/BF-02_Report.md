# BF-02 / C 共用的恒曲率 CSF 载荷报告

状态：**PASS（仅 BF-02）**。新不可变 `CchBalancedCsfLoad` 已实现恒曲率实际动量 RHS，并通过原门槛下的实际两相静态鞍点、符号、完整自由面/壁面和来源检查。本轮未实现变量曲率或时间推进；BF-03/H14/H15尚未启动，完整BF-CSF及原H13物理gate没有被本结果替代。

任务格式/边界见 [BF-02 TODO](BF-02_TODO.md)，使用说明见 [API](../api/BF-02_Cch_Balanced_Csf_Load.md)。本轮与已封存BF-01对应的当前真实源码/输入比较；所有构建、临时文件和输出在D:。

## 范围与实际方法

3个production文件：新 `include/Modules/Hydrodynamics/CchBalancedCsfLoad.h`、新 `src/Modules/Hydrodynamics/CchBalancedCsfLoad.cpp`、根 `CMakeLists.txt` 的唯一CCH注册。2个测试文件：新 `tests/contracts/hydrodynamics/cch_balanced_csf_load_contract.cpp`、`tests/CMakeLists.txt` 注册。后者将同一个可执行的core和domains两组独立注册，均保持原180秒timeout、OMP2环境和物理阈值。没有新增第三个测试文件。

唯一load owner保留原 `CchPressureLoad` 和 `CchRegularizedColor`，通过已有 typed `CchNativeValues` 调用一次 `C.apply(color)`，然后形成

\[
g=Cc_\epsilon,\qquad S_\sigma=+\sigma\kappa_0g.
\]

`momentumRhsLoad()` 已经是实际RHS，调用者直接相加。不沿用H12/H13的balanceLoads取负入口；没有第二gradient、面积或密度倍乘、虚拟粗面均分、解析pressure写入solver、材料color替代、独立界面状态或几何复制。

sigma有限非负，kappa0有限且允许正负。零sigma/零kappa0时所有RHS行精确正零，仍保留真实Cc诊断；非零输入的系数溢出/下溢为零，以及非有限面载荷明确拒绝，无floor。构造检查layout与原D owner；运行检查原pressure matrix、geometry/full clock，实际color重载进一步检查同一不可变color owner及epsilon。副本共享合法，重新构造的等值对象不能误当成已绑定owner。

闭域按每个真实压力分量使用物理体积gauge：`p*=sigma*kappa0*(c-mean_component,V(c))`。自由行平衡；固定行保存 `R_wall=Au+Cp-S`，静态参考为 `-sigma*kappa0*C(mean_component(c))`。不要求壁C1=0；纯液只有自由行精确零，固定边界的C1载荷/反力保留。

## 实际构建和验收

VS2022/MSBuild17.11.9、MSVC19.41.34120.0、CMake3.29.5-msvc4，x64 Release。构建沿已有 `build/h13_acceptance` 实际CCH库；最终 [Release日志](../../runs/checks/BF-02/20261009_160120/build_csf_analytic_curvature.log) 无新warning/error。

| 检查 | 实际结果 | 证据/原门槛 |
| --- | --- | --- |
| core注册契约 | PASS，86.6855秒 | 原180秒；48个rho/mu/几何场景、146次真实saddle、26来源/参数拒绝、50符号负控 |
| domains注册契约 | PASS，9.58667秒 | 原180秒；18个场景、54次真实saddle；孔洞、断连、annular及纯相 |
| 常曲率全row gradient/RHS | 最大相对误差1.95547e-16 | 原1e-12；独立native incidence和物理边面积，不读取C/D系数作参考 |
| 实际静态速度 | 最大相对量4.87619e-15 | 原H12 1e-11 |
| 实际pressure/gauge | 最大归一误差2.66454e-15 | 原H12 1e-10；各真连通分量物理体积mean |
| 自由行Cp-S | 最大相对量2.36633e-15 | 原H12 1e-11；全部native自由面逐行保存 |
| 壁反力 | 最大相对误差6.70182e-16 | 1e-11；Au+Cp-S与独立component mean参考 |
| 原局部continuity gate | 最大7.97833e-12 | 原saddle 1e-11；没有修改solver或接受条件 |
| 旧相关回归 | 10/10 PASS，377.38秒 | architecture/self-test/layout、pressure/saddle/真实测度、phase geometry/viscosity、原H12及BF-01 |

两新契约覆盖XY/RZ、uniform和两种实际2:1 AMR图样、真实circle/sphere；正曲率分别为独立解析1/R、2/R。密度比1、1e3、1e6，黏度0及1e-3，气黏度为其0.01倍；使用原sharp两相质量/原分相黏性，mu=0算子精确消失。还有六个确定随机平面位置及原标量7倍缩放；负kappa、零kappa、零sigma及极端零参数准确性均核验。

真实正S解的液侧压力较高；XY/RZ uniform凸滴均额外实际解反号-S，独立压力检查检出了反转。selected域包括孔洞和两个断开压力分量；RZ断连采用sphere，让两分量color均值确实不同，不用一个全域mean替代组件gauge。实际FirstJet regrid后新phase/color/C/load重建成功；旧mesh、旧identity、完整clock字段、另一个D/C/color owner和同pub不同epsilon均明确拒绝。

首次两组测试在有限mu轴面的 `reaction==0` 精确断言失败。该断言混淆了零面积C/CSF结构零与原黏性固定行的Au舍入量：旧H12只要求axis load精确零。正确核验改为S和Cp精确零、mu=0时Au精确零，并保留每个完整reaction及原wall<=1e-11。最终528条轴行的有限mu最大|Au|为4.36588e-37，全部壁门槛通过。该修正未修改任何生产算子、黏度、phi或门槛；首次失败日志原样保留。

## 真实数据和来源

最终 [CTest原始日志](../../runs/checks/BF-02/20261009_160120/csf_delivery_ctest.log) / [完整JUnit](../../runs/checks/BF-02/20261009_160120/csf_delivery_ctest.xml) 保存所有stdout；[回归日志](../../runs/checks/BF-02/20261009_160120/regression_ctest.log) / [JUnit](../../runs/checks/BF-02/20261009_160120/regression_ctest.xml) 独立保存。两次通过版本及首次失败没有覆盖。

[cell CSV](../../runs/checks/BF-02/20261009_160120/csf_cells.csv) 共9,528行，含color/sharp alpha、真实volume、component、实际pressure与独立expected/gauge；[face CSV](../../runs/checks/BF-02/20261009_160120/csf_faces.csv) 共22,002行，含真实side/area/free/axis、Cc/S/Cp/Cp-S/Au及实际/参考reaction。独立后处理对保存数据逐行复核，最大绝对pressure误差6.66134e-15，自由Cp-S 1.23323e-15，壁reaction误差1.33227e-15；零参数及528条轴面Cc/S/Cp结构零、无黏轴Au零全部成立。详见 [数据汇总](../../runs/checks/BF-02/20261009_160120/csf_summary.json)。

与BF-01的1,400份数值源码/输入逐项SHA256比较，仅两份CMake注册改变，新增类对及本契约；BF-01/Core/Interface/C-D/材料/KEP/黏性/鞍点/旧H12-H13数值源码与输入逐字保持。新树身份 `4c1356180e3d57c72494503e96703bf30aa307d8c25a62c6eb67d5edf7251e87`；见 [来源guard](../../runs/checks/BF-02/20261009_160120/source_diff.json) / [注册diff](../../runs/checks/BF-02/20261009_160120/registration.diff)。BF-01已封存archive的SHA256重新核验保持不变。

## 范围限制与下一任务

本任务只验证恒曲率代数和上述实际静态范围。BF-04的完整72几何原矩阵、BF-05无黏完整周期/能量/全谱、BF-06动态AMR/边界与H14许可尚未完成。旧H13 FAIL、旧模态完整momentum契约180秒TIMEOUT及历史缺失标记保留；本轮没有把它们改判PASS。

按已授权路线，**下一次允许仅执行BF-03**：以同一FirstJet实际surface quadrature作局部正权真实曲率投影，扩展当前load owner；必须先立同格式TODO。本轮未写BF-03代码/TODO，也未启动H14/H15；总体目标继续进行。
