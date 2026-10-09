# H10e / 有界最终算子与材料空间矩复核

2026-10-04，**H10e FAIL / H10未关闭**。最终算子的硬constant/affine、中心二次余项区间、行系数界与skew联合候选未达到验收精度，已归档退役；生产和原短契约恢复H10c冻结字节。新增被动材料观察器给出真实非均匀材料的精确affine/skew相容性缺口及加密数据。下一单任务[H10f](H10f_Compatible_Bounded_Moments_TODO.md)；不进入H11/H08a，不改变Pure FirstJet或材料非守恒路线。

## 最终算子的联合有界试验

[任务书](H10e_Bounded_Operator_Consistency_TODO.md)及[实施前设计](H10e_Bounded_Operator_Consistency_Design.md)从最终闭合端口的`N=(K+diag(q)/2)u`出发。K为统一真实one-ring上的共享反对称pair，q是同一H09唯一材料F的外流。所有原native自由度/质量保留，重复周期采样使用明确几何等值关系；不依赖问题名、解析参考、AMR分支或隐藏载荷修正。

试验前固定硬C/P1、中心二次矩区间`|sum(K*delta_a*delta_b/h²)| <= s*h/L`及行l1界`sum|K| <= 2*s`；s为该行实际F的总绝对率。这使实际二次/三次Taylor项具备所需h尺度的必要控制，但**不预先证明这些约束有交集，也不证明原完整PDE通过**。

局部Dykstra、完整硬矩初始化、联合active refinement以及秩揭示SparseQR硬矩投影均只改变同一数学问题的求解实现，没有根据参考误差调界、加黏度或回退。联合active试验通过XY uniform n=8/16及2:1 n=8，但2:1 n=16在固定预算内失败。最新SparseQR/Dykstra在uniform n=8通过，独立矩阵与actual apply差≤5.02e−16；在uniform n=16达到14,784周期仍失败，常数/一次gross误差分别为9,445,469.95与36,656,296.46倍256eps。界比率接近1及行l1比率0.283464不能抵消硬矩失败。

这证明**这些求解实现未通过**，不证明均匀材料的有界凸交集不可行。没有扩大迭代预算来包装PASS，没有形成全部n=8/16/32有界候选数据，更没有新增配对全场景短契约或完整时间窗。候选source/实际binary/运行关联保存在[失败档案](../../runs/baselines/H10e_20261004/rejected_sparse_candidate/README.md)和同目录active档案；它们不属于现行生产API。

## 真实材料的独立必要条件

保留原标量质量和`R/uR`。令最终实际速度导数

`a=(load-u*dm/dt)/M=(-N+u*q)/M`。

原material-KE/skew恒等式要求

`sum(M*w*a) = sum(q*w²)/2 - P_boundary = Phi`。

取规定advector沿周期方向、w沿横向非周期坐标的admissible affine样本，PDE导数为0。Phi非零时不可能逐点机器精度满足该PDE与同一功恒等式；Cauchy不等式还给出加速度质量RMS下界`|Phi|/sqrt(sum(M)*sum(M*w²))`。这里dm与R原样消费后代数消去，不是将它们设0。

独立地，硬几何一次目标`b_i=sum(sign*F*(segment_centroid_z-z_i))`和`K*1=q/2`必有

`K*z = z*q/2+b`，因此`Psi=sum(z*b+q*z²/2)=0`。

Psi非零即与共享skew不相容，**扩大支撑也不能消除这个全局必要条件**；但改为有受控离散余项的二阶相容设计仍可能。Psi针对当前几何centroid目标，Phi针对实际PDE/功；二者不能混为同一残差。

[现有研究driver](../../cases/benchmarks/cch_convection/driver.cpp)新增被动`compatibility`/`compatibility32`：n=4/8/16/32、XY/RZ、uniform/真实中央2:1、rho=1/1000/1e6，共48组。初始圆形Pure FirstJet经真实characteristic推进产生begin/mid/end，dt=0.1/n、最细CFL=0.08；三层相几何、真实dual质量、唯一F、非零R/G均由原H07–H09产生。XY取w=(1+y,0)、advector=(0.4,0)；RZ取w=(r,0)、advector=(0,0.4)，满足轴约束及声明研究周期采样。

观察器不构造或修改convection系数。native求和与独立逐segment表达式交叉核验，48组最大gross相对差2.70e−16；rho=1功对照最大1.52e−16。真实周期测度/方向/F匹配有独立断言，观察器exit 0仅表示诊断身份有效，不是PDE PASS。最初固定dt=0.1在加密时触发既有半格CFL，被拒绝的诊断日志保留；最终证据均来自合法dt。

XY 2:1、rho=1000得到：

| n | Phi | Psi | 加速度RMS下界 |
| --- | ---: | ---: | ---: |
| 4 | 0.9858170 | −0.03379529 | 0.00402768 |
| 8 | 0.4329269 | −0.01031310 | 0.00177017 |
| 16 | −0.1002380 | −0.002780179 | 0.000409830 |
| 32 | −0.03045037 | −0.0007188339 | 0.000124499 |

Psi后两段阶1.891/1.951，RMS下界阶2.111/1.719。rho=1与本组uniform及RZ对照近舍入精度，rho=1e6的XY 2:1缺口同样明确非零。细网格缺口减小支持研究二阶余项，**不证明二阶不可能，也不证明已构造二阶算子**；几何与dt同时加密，该48组只核验必要条件，不替代固定原H10b研究。完整数值见[material_summary.json](../../runs/baselines/H10e_20261004/material_summary.json)及最终两个CSV。

## 恢复、验收与限制

H10e全部候选从生产退役，CchKepConvection.h/.cpp和原短契约为H10c精确字节。现行没有新增配对API、bounded report字段或SparseQR依赖。研究driver仅增加上述被动观察器与两个模式；原Study、rhs、full/RK4/donor/gate及rowAudit保持H10d版本，可逆scope核验记录在冻结目录。

标准Release构建通过。原短契约124场景、1,356 trace、6,144 load、24实际跨界、19拒绝、4 regrid及1,516空间矩/3,504pair通过，数值SUMMARY与H10c/H10d一致；原1e−12载荷/动量/KE与1e−10常速度门槛保持。架构、自测、目录3/3通过，加短契约共4/4受影响CTest。前次19/19回归是历史证据，未冒称本轮全部重跑。

原内部affine n=8/16/32、四组合共15,344行复验完成：内部最大6.59128e−13≤1e−10，原XY AMR n=16 dof911仍为−1.0000000000000742、误差7.41629e−14。外部误差始终保留，不能将内部反例通过改称全部空间PASS。见[affine_summary.json](../../runs/baselines/H10e_20261004/affine_summary.json)；全部实际观察日志保存于[runs/checks/H10e/20261004_024250_022](../../runs/checks/H10e/20261004_024250_022)。

**本轮未重跑原完整H10b**：新的有界算子前置失败。最新完整研究仍是H10c原28窗口/15,116步、84gate中10空间+4donor失败；原门槛未修改。无H11、H08a、一般Hodge改造、CCH流体App或高密度比低黏稳定性新结论。

## 冻结与下一步

开始独立复核H10d的121个payload。结束冻结源文件/输入/原App/PDF与实际交付关联；受保护生产、工具、CMake、Interface和材料定义保持，仅现有研究driver/README为长期修改。冻结身份及独立verification放[runs/baselines/H10e_20261004](../../runs/baselines/H10e_20261004)，工程PASS与方法FAIL分开。

H10f先对均匀材料有界交集给出可行解或不可行证书，再依据真实材料相容缺口设计有二阶依据的最终空间矩；不能继续强加不相容目标、仅加迭代或把一次数值失败解释为一般方法不可能。保持用户指定的最终native operator设计方向与Pure FirstJet独立、无库存/bubble，空间离散通过后才启动H08a。
