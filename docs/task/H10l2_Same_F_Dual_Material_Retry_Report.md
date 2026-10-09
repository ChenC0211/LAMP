# H10l2 / 同一 F 的 dual 材料迁移验收

**PASS，限定本项迁移范围。** 2026-10-04完成实现、计算与独立数值验收，2026-10-05补齐报告并复核归档。[任务书](H10l2_Same_F_Dual_Material_Retry_TODO.md)、[实施前设计](H10l2_Same_F_Dual_Material_Retry_Design.md)的预算与门槛保持；上一轮为progress，本项goal仍active。下一唯一[H10m](H10m_Final_Native_Same_F_Convection_TODO.md)尚未实施；H10整体、原完整H10b、H11与H08a没有由本项自动通过。

## 当前实现

H09直接消费H08唯一 `CchMaterialFlux` 的dual视图。`CchDualMaterial::faceRates()`返回 `material().dualFaceRates()` 的同一vector引用；只保存dt乘积，不再重新积分、乘密度、按面积分摊primal F或求quarter。内部面反号共享、真实AMR子面和与原actual B/phase来源贯通。

真实三层dual质量仍由同一Pure FirstJet几何取得。以实际incidence计算

\[
q_d=\sum_s sF_s,\quad r_d=\Delta m_d+\Delta t q_d,\quad R_d=r_d/\Delta t,
\]
\[
E_{sp}=q_d-Tq_p,\qquad G=\Delta m_d-T\Delta m_p,
\qquad r_d=Tr_p+G+\Delta t E_{sp}.
\]

T仅为原真实overlapFraction，E_sp为率、G/r为质量。三者允许非零，不赋伪零、修复质量或反向约束界面。原consistencyImpulse、materialKineticChange数值逐字保持；KEP生产代码未改，继续使用同一F/R。

删除quarter/SPD局部材料延拓及其精确交换断言，退役最后 `CchPhaseRate::evaluate`、optional legacy相率存储和Transport legacy表达。现有 `isSpatial()`仅保留恒true兼容查询，不存在第二数值分支。相积分数值、唯一Flux owner、原D_I/C、actual advector与显式边界law均未修改。现行来源、字段、单位与调用见[H09 API](../api/H09b_Cch_Dual_Material.md)。

## 迁移前证据与范围

启动复核H10l1：323 payload、119 current canonical、实际source/runtime及九份前置档案通过。逐文件证明本项5 production、2 tests与H10l原候选启动前源码相同，才提取并审阅其迁移改动；重新编译、运行和独立核验当前H10l1上的结果，没有把旧失败候选或旧研究的通过数据当成当前验收。

本项5 production为 `CchPhaseGeometry.h/.cpp`、`CchMaterialTransport.cpp`、`CchDualMaterial.h/.cpp`；2 tests为现dual/KEP契约。5份生产文件净减少173行。原研究driver只增加两段H10l2 observer；删除这两段后逐字恢复本项启动版本，包含H10l1观察器及原full/RK4/donor/gates。

起始1,044 source成员不变，1,036份哈希不变；仅上述7份与原research driver变化。0新源码/target/tools。Core/FirstJet、phase spatial数值、Flux owner、advector/boundary、D_I/C、KEP、App/案例配置/冻结PDF保持。源码与接口更新没有扩散到额外production文件。

## 原失败回归已实际完成

H10l的3组basis拒绝由H10l1稳定face-local积分修复；这次迁移后 **208个非轴单面basis全部实际完成**，并通过704个线性组合segment检查。H10l的4组closed流函数在现KEP测试中显式列出真实boundary DOFs并绑定typed不可穿透law，逐一核验same identity、全部指定normal polynomial两端系数、phase atom和material atom严格0，以及同一dual F引用和实际boundary segment零率。

未按mean==0、算例名或AMR位置压通量。四closed场景最大无黏对流功率尺度误差 **2.776e−17**，全部后续constant/KEP检查实际运行；dual/KEP两项均输出 `MIGRATION_RESULT,0,PASS`。旧H10l FAIL的归档与历史结论保持。

| 当前短契约 | 实际覆盖 / 结果 |
| --- | --- |
| dual总体 | 537场景、95,520 segment检查、29,252 dual检查；208 basis、24实际跨面、12鞍点消费、22来源拒绝、12时间行、4 regrid |
| 独立spatial dual | 264准备、42,240 atoms；最大独立物理F误差1.928e−16 |
| 零净反向相流 / 小相自身尺度 | 32 / 136个见证，均保留非零相/材料贡献 |
| KEP总体 | 124场景、1,356 trace检查、6,144载荷检查、24实际跨面、20拒绝、4 regrid |
| 最终native原矩检查 | 1,516矩、3,504 pairs，最大尺度误差5.076e−13 |
| KEP独立载荷 / 动量 / 能量身份 | 5.451e−16 / 8.400e−16 / 5.457e−17 |

高密度跨面常速度/鞍点测试沿用原各自门槛：dual常速度最大9.702e−13、鞍点速度1.610e−15、压力4.668e−11；KEP常速度最大8.513e−11，不把这些不同归一化结果合并成统一1e−12声明。关闭uR与真实几何转移的对照仍暴露非零漂移。原H10l1 phase/material契约原文不改，208 closure、40 boundary guards、4实际counterflow及352 causal scale拒绝均在必要回归中保留。

## 独立实际研究

`dual`运行从2026-10-04 20:33:53至20:51:16，原进程正常exit0、OMP1，没有超时重启或截短。XY/RZ、uniform/真实2:1、n4/8/16/32与原四selected布局共64布局；原56个coarse-selected tangent缺口保留，无新增豁免。

三族为静止平面/MMS、移动圆球/MMS、移动圆球与同向0.4平移的coupled制造族；每族三密度1/1e3/1e6，实际192个phase factory、576个材料准备。独立脚本从声明多项式、解析相交与circle-rectangle primitive、真实native face/pressure坐标重算相率/F、三层质量、half-cell overlap、incidence和各项缺陷，不以production F/R作为物理参考。

| 独立CSV核验 | 实测 |
| --- | ---: |
| dual cell / shared segment行 | 348,930 / 1,071,360 |
| pressure cell / primal face / overlap行 | 162,180 / 348,930 / 655,200 |
| 独立相率 / 材料F最大尺度误差 | 3.775e−15 / 3.766e−15 |
| 独立三层体积 / 质量误差 | 4.613e−13 / 4.613e−13 |
| 真实R/E_sp/G分解与incidence/global balance | 5.456e−15 |
| uR/材料动能身份 | 0 |
| 三层×三密度gas dof103见证 | 63段，独立误差0 |
| coupled gas dof103的三份interval q/R | 均0 |

脚本按实际gross/work核验，零work要求精确零；没有用经验floor或抵消后的E_sp自身放大舍入误差。最大实际 |E_sp|=298,456.526、|R|=519,600.343，包含高密度比强制运动族的有量纲诊断；不宣称一般R小余项、coupled收敛或完整PDE通过。单个pure-gas见证为零也不表示整场R为零。

原内部n8/16 affine四组合（八组）重跑，3,760行中3,024内部行最大误差 **1.084e−13**，满足原1e−10门槛；外边界原始误差保留。没有运行n32 affine或原完整H10b。原28窗口/84 gates及失败历史保持，下一任务直接面向最终native算子与整体空间精度。

## 工程交付与后续

标准Release六targets构建成功；必要CTest **8/8**：architecture、architecture_self_test、repository_layout、phase、material、momentum、dual、KEP，总177.50秒。最终文档/架构检查另 **3/3**。这次没有数值编译或回归失败；工具读取中两次无效Select-Object参数及一次API文档同路径delete/add补丁拒绝均未修改数值路径，随后修正调用方式，不改变门槛或已运行结果。

`compiled_runtime.json`绑定33份source与8份实际exe/lib/DLL；research exe SHA256前缀 `FCA559FFA6C5F644`，dual/affine退出记录均匹配该构建，CTest对应同一批二进制。2026-10-05续接时确认无存活研究进程，重新核验当前source/runtime、scope、已保存实际basis/closed/affine结果，没有重启已结束研究。

检查目录 `runs/checks/H10l2/20261004_203031_209`；独立原始CSV/脚本、实现前设计、运行日志、编译source/runtime、scope及flat SHA256 manifest见[H10l2冻结基线](../../runs/baselines/H10l2_20261004/README.md)。目录日期沿用任务启动日；最终复核时间在verification.json。十份前置archive仍按各自manifest校验，不修改旧payload。

下一唯一H10m已写任务书：在同一F与真实R上联合最终native constant/affine空间矩、二阶余项及material-KE/skew；内部反例先于原完整H10b，维持所有原gate。材料仅四硬约束，Pure FirstJet独立，无库存/bubble；H11等待，H08a时空材料通量升级随后另项实施。本项材料时间积分仍为中点近似。
