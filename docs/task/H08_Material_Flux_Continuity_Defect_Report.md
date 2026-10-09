# H08 / 中点材料流量与非守恒连续性缺陷验收

日期：2026-10-03。**H08 PASS，下一项 H09。** 按[任务书](H08_Material_Flux_Continuity_Defect_TODO.md)完成唯一材料流量、中点时间积分及实际非零连续性缺陷；总目标和 Gate H 保持未完成。

## 实现与范围

[公开API](../api/H08_Cch_Material_Transport.md)在正式CCH提供不可变 `CchMaterialTransport`。开始/中点/结束三层几何共用同一layout/D_I/实际mesh，geometry identity各自保留；两相密度与完整区间绑定在同一对象。材料质量直接消费H07a几何，面材料率直接消费同一H07a液气率，不另设库存、不插值另一份 `rho_f*u_f`。

`I_rho=dt*(rho_l*q_l+rho_g*q_g)`是具名中点近似；逐cell保留几何质量变化、共享面外流积分、`r_int=delta_m+outflow`及`R=r_int/dt`。非零缺陷不填零、不修复质量/界面/通量。总cell缺陷与质量变化+选中边界积分满足同一代数收支；该恒等式不代表材料质量守恒。

## 独立数值验收

cell参考为独立解析平面/圆/球积分，face参考为独立线段/圆交截；方向由Core原生面邻接直接构造，不借生产D_I装配参考。时域参考用圆弧/球截面的独立解析原函数，参考不进入生产接受条件。

| 检查 | 实测 |
| --- | --- |
| 场景 / cell / face | 309 / 6,780 / 16,444；XY/RZ uniform/2:1、选中半域、纯相/平面/圆球 |
| 密度/速度 | 密度比1、1e3、1e6；零/正/反向面速度，RZ轴值严格零 |
| 三层材料质量最大归一化误差 | 8.1491e-16 |
| 面材料率最大归一化误差 | 2.9104e-16 |
| 质量变化/实际缺陷最大归一化误差 | 8.1491e-16 |
| 内部共享面 / 粗面聚合 | 11,676 / 1,000；真实fine subface积分聚合最大缺陷1.0701e-15 |
| 独立单面注入 | 8组；XY/RZ uniform/2:1及选中半域，各相邻cell收到相反积分，无关cell缺陷精确零 |
| 选中边界/总缺陷收支恒等式 | 最大缺陷3.5084e-16 |
| 极小液/气相 | XY/RZ各2^-40，4组均保留可观察的非零小相材料贡献 |
| 拒绝 | 28；密度、三层clock/identity/layout/D_I/实际mesh、溢出/下溢及退化中点 |

零速度/静态几何的积分缺陷精确零。失败不污染既有对象；实际regrid后旧来源拒绝，新mesh重新准备及等价mesh验证通过。副本共享只读材料数据和几何寿命。真实材料密度不设置最小值；不能表示的非零乘积/积分明确拒绝，不悄悄下溢为零。

实际 FirstJet producer 平面平移从 .49 到 .53，速度 .4、dt=.1、rho_l=1e6/rho_g=1，跨过原生 .5 面。整个面phase rate的时间跳跃使单次中点近似产生局部缺陷，结果如实保留：

| 序列 | L1积分质量缺陷 | 最大绝对R |
| --- | ---: | ---: |
| XY uniform / 2:1 | 19999.98 / 19999.98 | 24999.975 / 24999.975 |
| RZ uniform / 2:1 | 62831.79024 / 62831.79024 | 137444.54115 / 137444.54115 |

L1参考为 `.02*(rho_l-rho_g)` 乘横截面积1或pi；signed总缺陷接近零，不能以总量抵消隐藏上述局部缺陷。本项未把该不光滑情形冒称局部二阶或零残差；H09动量/动能须消费非零缺陷，H14/IH还需验收完整sharp-interface时间精度。

非多项式圆/球geometry的三档n=4/8/16材料测量，原门槛≥2保持：

| 序列 | n=16质量误差 | 质量末档阶 | 面材料率末档阶 |
| --- | ---: | ---: | ---: |
| XY uniform | 2.32243e-9 | 5.32524 | 4.74000 |
| XY 2:1 | 3.71448e-10 | 4.81560 | 5.26670 |
| RZ uniform | 2.11173e-10 | 5.41430 | 5.17448 |
| RZ 2:1 | 2.11173e-10 | 5.41384 | 5.18983 |

光滑移动圆/球穿过固定物理侧面的T=.1、N=4/8/16中点时间积分，实际生产owner准备112个区间。XY/RZ uniform/2:1均下降，末档阶分别1.99722/1.99722/2.00000/2.00000；N=16相对积分误差6.26343e-5、6.26343e-5、4.00281e-6、4.00281e-6。时间收敛使用解析制造publication隔离材料求积误差，另外的crossing使用实际FirstJet producer；两者不等同于完整两相PDE收敛。

## 工程验证与证据

3 production/build +2tests：新增header/cpp/rootCMake与新短契约/testsCMake。前置1028个文件仅两个CMake改变，1026保持，新增3个；312个case资产和31个fixture资产保持。cpp仅属正式CCH；依赖闭包不进入旧Hydro/Compatibility/App。冻结PDF和LAMP.exe身份均保持。

- local Release定向构建通过，无新增warning/error。首轮显式/隐式指针initializer-list类型不一致及混合auto声明已修正，原失败日志保留。
- 新契约 **1/1 PASS，39.91秒**；首轮小相测试的RZ径向速度为零，没有实际小相流量，已改为满足轴规则的径向线性速度，保留原失败；门槛、物理小相大小和材料参数未放宽。
- 前置回归 **14/14 PASS，106.79秒**：help、architecture/self_test、repository_layout、FirstJet App、H01–H07与旧pressure-gradient load。前置数值汇总/收敛行保持；App契约仍为12 smokes/49拒绝。
- 最终文档后通过repository_layout、本地链接与交付哈希核验。未跑全79项CTest、Debug或完整物理benchmark；没有改现有流体App行为。

证据：[verification](../../runs/baselines/H08_20261003/verification.json)、[operator detail](../../runs/baselines/H08_20261003/operator_detail.log)、[regression detail](../../runs/baselines/H08_20261003/regression_detail.log)、[source audit](../../runs/baselines/H08_20261003/source_audit.json)、[交付清单](../../runs/baselines/H08_20261003/delivery_manifest.json)。TODO先于源码，开始时核验H07b的23项交付/31项证据；原TODO/源码/输入身份、失败定位和最终审计脚本归档。

## 未解决范围

H09 dual几何质量/材料流量映射与常速度动量一致性、H10/H14对流/动能材料项与完整corrector、H15边界/Hodge/regrid、CCH App和低黏/EHD仍待后续任务。H08证明材料来源和缺陷核算以及光滑中点求积精度，不证明高密度比完整输运稳定或任意跨面时间二阶。
