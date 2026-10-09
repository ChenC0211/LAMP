# H05 / native CCH 真测度与轴正则报告

**状态：PASS，2026-10-03。** 按[任务书](H05_Rz_True_Measure_Axis_Regularity_TODO.md)完成冻结 CCH §3.2/§6.1/§12.1 的真实物理散度与轴边界规则，下一项 **H06：能量型黏性 K_mu**。本项没有实现strain/hoop黏性、RZ压力/速度PDE、AMR pressure Hodge、流体时间推进或Gate H。

## 实现与数学事实

[CchPhysicalDivergence](../api/H05_Cch_Physical_Divergence.md)保留同一权威D_I/layout，消费Core真V_c/A_f；apply调用唯一D_I得到积分净流量，再逐cell除V_c。没有另一套散度矩阵/gradient、经验1/r项或重新造网格。压力仍由积分C=-D_I^T进入[H04](H04_Single_Phase_Saddle_Report.md)，不使用物理归一化散度的负转置。

RZ轴面A_f=0并不足以保证正则，因此显式拒绝非零u_r，包含极小非零，不自动改零。具名接口为caller的正半径样本产生径向速度奇、轴向速度/pressure偶的负半径边界值，对应u_r=0与∂r u_z=0。反射只在选中实际RZ轴面的域启用，派生值不成为额外primary ghost未知量。反射实施边界条件，并不证明任意内部场光滑；轴上u_r/r的strain极限仍须由H06独立验收。

prepared测度/轴编号均只读，复制共享，同一cell-volume视图与原D_I来源保留；无借用mesh生命周期、IO/字典/主时钟。错来源/stale实际mesh、非法轴速度/反射数据及归一化非有限均拒绝，失败不污染原对象或先前返回结果。

## 独立解析与正则验证

72个测度场景=XY/RZ轴上/RZ离轴×n=4/8×depth=0/1/2×全域/左半域/两块不连通/带孔域。真实边端点与cell bounds独立计算full旋转体体积/面积，不用domain.face()/CartesianMeasure或数值矩阵构造参考。

每组比较三个解析场的准确face averages与体积平均divergence：规则多项式（径向奇/轴向偶）、正则stream无散流、径向线性流，共216组。RZ多项式u_r=.7r+.2r³-.15rz²、u_z=.45+.3r²z-.25z³，准确div=1.4+1.1r²-1.05z²；XY按其本身Cartesian微分积分。另核对真实体积加权总散度与独立选中域边界通量闭合。

| 指标 | 最大归一化缺陷 | 预先gate |
| --- | ---: | ---: |
| 72组V_c/A_f端点测度 | 0（本组网格实测） | ≤1e-13 |
| 216组解析体积平均散度 | **3.32760e-15** | ≤1e-13 |
| 全局真实体积/边界通量闭合 | 3.04691e-16 | ≤1e-13 |
| 轴极限/独立奇偶参考 | **3.01346e-15** | ≤1e-13 |

9组径向线性流以n=4/8/16、域半径1/1e-3/1e-6检查div=2a在全部cell（含轴邻域）有限，不在轴上除r。216组反射比较正/负半径解析场，轴向反射中心导数精确为零。

16项拒绝涵盖错layout、拉伸build/apply、真实顶点改变、AMR后build/apply、非零轴面及denorm极小非零、XY/离轴反射、非法分量/非有限反射、物理归一化溢出。两项关键反例分别确认：非零轴速度的原D_Iu仍精确零，但physical apply拒绝；有限极大速度的原积分散度有限，但除V_c溢出后必须拒绝。复制/输入寿命、独立等价mesh、零输入和失败后恢复通过。

## 非多项式RZ散度MMS

单位RZ域的u_r=r sin(πz)、u_z=cos(πr)sin(πz)，准确div=2sin(πz)+πcos(πr)cos(πz)。production输入native face-center采样；准确参考按full 2πr cell volume独立积分，全部cell（含轴邻域）进行物理体积加权相对L2比较，不读D_I/C来反造参考。

| n | physical divergence相对L2 |
| --- | ---: |
| 4 | 3.76625e-2 |
| 8 | 9.32336e-3 |
| 16 | **2.32500e-3** |
| 8→16观测阶 | **2.00362** |

误差递减，独立重算末档阶符合≥1.8的运行前门槛，n=16误差约0.2325%。这是连续性测度/散度MMS，不是pressure/velocity、黏性、AMR Hodge或时间精度研究；不能从准确face-average多项式的舍入级再现推广一般零误差。

## 构建、回归和修复记录

最终local Release无新增compiler warning。新契约 **1/1 PASS（0.92 s）**，另10项定向回归 **10/10 PASS（69.73 s）**：lamp_help、architecture/self_test、repository_layout、FirstJet App、H01/H02b/H03/H04、旧hydro_gradient_load。没有重跑或声称完整75项CTest、Debug、论文完整benchmark或流体/EHD。

H01原23场景/29拒绝及H02b/H03/H04完整数值汇总逐字保持。前置8份CCH源码、所有旧production、architecture checker、LAMP.exe、冻结方法与243份case输入不变。开工1,016来源中1,014未改，仅root/tests CMake变更，新增header/cpp/短契约；总预算 **3 production/build + 2 tests**。新.cpp只有CCH owner，依赖仍Core/Eigen，无旧Hydro/Compatibility。

初次构建已编译新实现，但错误使用不存在的architecture target；改为现行quantum_architecture_check。首次数值检查暴露新夹具二级加密与一级共享轴向边界、形成non-2:1接口；修正为嵌套加密，保留轴面、两级AMR与全部selection组合。另在首次运行前校正独立XY多项式微分的z²系数。以上只修正新测试夹具/构建命令，production实现、参数、样本数量和验收gate保持；原失败日志保留。

## 交付与复现

[公开头](../../include/Modules/Hydrodynamics/CchPhysicalDivergence.h)、[实现](../../src/Modules/Hydrodynamics/CchPhysicalDivergence.cpp)、[短契约](../../tests/contracts/hydrodynamics/cch_physical_divergence_contract.cpp)，以及root/tests CMake、本TODO/API/报告和现行索引/路线图。

```powershell
.\tools\build.cmd -Target quantum_cch_physical_divergence_contract
ctest --preset local -R '^cch_physical_divergence_contract$' --output-on-failure
```

CTest在可用开发shell执行；通用完整短检查仍为tools/build.cmd -Test。临时定位脚本归档[来源/数值审计](../../runs/baselines/H05_20261003/engineering_scripts/source_audit.ps1)与[交付审计](../../runs/baselines/H05_20261003/engineering_scripts/delivery_audit.ps1)，不进入生产或通用tools。

证据：[独立误差/阶](../../runs/baselines/H05_20261003/verification.json)、[新契约完整记录](../../runs/baselines/H05_20261003/operator_detail.log)、[十项回归](../../runs/baselines/H05_20261003/regression_detail.log)、[最终构建](../../runs/baselines/H05_20261003/build.log)、[契约重编译](../../runs/baselines/H05_20261003/operator_rebuild.log)、[原构建失败](../../runs/baselines/H05_20261003/initial_build.log)、[原夹具失败](../../runs/baselines/H05_20261003/initial_operator_detail.log)、[来源审计](../../runs/baselines/H05_20261003/source_audit.json)、[依赖图](../../runs/baselines/H05_20261003/architecture_graph.txt)、[文档审计](../../runs/baselines/H05_20261003/document_audit.json)、[交付清单](../../runs/baselines/H05_20261003/delivery_manifest.json)。最终文档更新后单独复跑repository_layout并保留其记录，不覆盖十项回归日志。
