# H06 / native CCH 能量型黏性验收

日期：2026-10-03。**H06 PASS；下一项 H07。** 按[任务书](H06_Energy_Viscosity_TODO.md)完成实际应变、能量黏性与无黏极限；总目标和 Gate H 保持未完成。

## 实现与范围

[生产接口](../api/H06_Cch_Energy_Viscosity.md)在独立 CCH 库提供 CchStrainOperator/CchEnergyViscosity。E_h 从实际 native 面速度、局部两层支撑构造，不接受调用者自造 E；非轴 P2 重构精确保留当前单元 normal face 值，规则 normal 导数回到中心差分。2:1 patch 使用同一缩放多项式再现思路。

RZ 显式包含 hoop；轴邻域将径向奇/轴向偶多项式限制写入系数空间，r=0 用 ∂r u_r 极限。仅验证 axis face=0 不足以防止近轴除小半径引起的抵消误差，此约束已在半径1/1e-3/1e-6的解析场通过。复用H05反射和同一D_I/layout；无新速度 primary、第二散度/压力矩阵或旧 Hydro 依赖。

H_mu 的实际行权重为2μ_q dV_q，K_mu=E_h^T H_mu E_h。单相GL2×GL2使用真V_c及RZ物理权重；调用者求积端口不推断phase。H07仍需唯一Hermite phase geometry/material identity与切分覆盖；完整no-slip/free-slip/open边界归H15。没有人工最小黏度或黏性稳定化。

## 独立数值验收

| 检查 | 实测 |
| --- | --- |
| 实际场景 / 解析与刚体场 | 21 / 55；XY/RZ轴上/离轴、uniform/2:1 depth0/1/2、半域/不连通/孔洞 |
| 多项式应变最大归一化误差 | 1.0192e-13 |
| 单相求积体积/r/z一阶矩缺陷 | 5.5511e-17 |
| 规则 normal 中心差分缺陷 | 1.4211e-14 |
| K_mu 对称相对缺陷 | 8.3693e-17 |
| 与独立逐点2μD:D求积/矩阵功的最大缺陷 | 5.6547e-13 |
| 刚体耗散最大值 | 6.7998e-28 |
| 6个谱 | PSD至舍入误差；XY仅3刚体零模，约束RZ轴值后仅1轴向平移零模 |
| 105个黏度缩放 | 0/1e-3/1e-6/1e-12/1e-300；归一化缩放缺陷4.0081e-15 |
| 6个实际K_mu进入H04单体求解 | 速度/压力/连续性/均值最大缺陷1.8652e-13 |
| 拒绝 / 轴极限 | 25 / 3；身份/AMR/尺寸/材料/点/轴/rank/构造与应用溢出 |

全部权重为零时矩阵nonzeros=0、力/耗散精确零；有限但可使E u溢出的极端速度也通过零黏度检查。正黏度下同样的overflow输入明确拒绝，不把失败当作零耗散。数据/来源失败不污染已有对象；独立等价mesh与原mesh析构后持有对象均通过。

单相μ=1的非多项式MMS，使用真实面中心值；精确strain独立解析计算，积分耗散参考使用独立GL8体积分。网格n=4/8/16，原末档门槛≥1.8保持：

| 序列 | n=16 strain相对L2 | strain末档阶 | 积分耗散末档阶 |
| --- | ---: | ---: | ---: |
| XY uniform | 0.0207348 | 2.05945 | 1.88160 |
| XY 2:1 AMR | 0.0125628 | 2.15660 | 1.84334 |
| RZ uniform | 0.0172031 | 2.21747 | 1.87242 |
| RZ 2:1 AMR | 0.0132470 | 2.35570 | 1.81103 |

误差为整个体域的strain求积范数，不是流体解的PDE误差。6个H04求解使用离散制造RHS验证真实K_mu可消费、同源C/D_I和压力均值保持；没有用离散RHS冒充连续黏性NS收敛。

## 工程验证与证据

3 production/build +2tests：新header/cpp/rootCMake与新短契约/testsCMake。前置1019个审计文件仅2个CMake修改，1017保持，新增3个；243个case输入、H01–H05全部实现/旧测试、冻结PDF、LAMP.exe与旧生产实现保持。source有唯一lamp_hydrodynamics_cch owner，依赖closure仍为Core/Eigen。

- local Release定向构建通过，无新增warning。首轮编译使用了当前Eigen不支持的接口/初始化方式及Vector2D除法，已修正；原失败日志保留。
- 新短契约 **1/1 PASS，1.44秒**。首轮无约束宽模板的strain/energy收敛未过、近轴抵消误差及不产生overflow的测试输入均保留原FAIL。通过改局部重构约束/正则空间和真实放大导数fixture解决，未降低gate、缩短窗口或放宽阈值。
- 定向回归 **11/11 PASS，70.37秒**：H01–H05、旧pressure-gradient load、App接线、architecture/self_test/repository_layout/lamp_help。H01–H05数值汇总与前置证据逐字相同。
- 最终文档后再次通过repository_layout和本地链接/来源/交付哈希检查。未跑全76项CTest、Debug或完整物理benchmark；没有改现有流体App数值行为。

原始证据：[verification](../../runs/baselines/H06_20261003/verification.json)、[operator detail](../../runs/baselines/H06_20261003/operator_detail.log)、[regression detail](../../runs/baselines/H06_20261003/regression_detail.log)、[source audit](../../runs/baselines/H06_20261003/source_audit.json)、[交付清单](../../runs/baselines/H06_20261003/delivery_manifest.json)。临时审计脚本归档到baseline/engineering_scripts；原TODO身份和H05四十项payload开始时均核验。

## 未解决范围

H07同源相分割/材料与非守恒Interface的质量—动量耦合设计、H08/H09共享质量与dual-mass通量、完整边界/Hodge、对流/时间推进/表面张力、CCH App及低黏/EHD仍待对应任务。H06证明黏性算子耗散且零黏度消失，不证明完整Euler/两相/电喷雾算法稳定。
