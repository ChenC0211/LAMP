# H10k / 唯一空间材料通量与来源迁移报告

**2026-10-04当前收口：** [H10l2](H10l2_Same_F_Dual_Material_Retry_Report.md)已PASS；旧quarter/legacy evaluate退役，H09直接共享H08原dual F，真实R/E_sp/G与uR/材料KE通过。现行接口见[H09 API](../api/H09b_Cch_Dual_Material.md)。下一唯一[H10m](H10m_Final_Native_Same_F_Convection_TODO.md)联合最终native空间矩与material-KE/skew，内部affine先于原完整H10b。Pure FirstJet独立；R/E_sp/G不要求零，H11/H08a继续等待。

2026-10-04，**PASS / engineering PASS**。仅关闭 H10k，goal 保持 active，本轮为 progress。按[任务书](H10k_Material_Spatial_Owner_Migration_TODO.md)、[实施前设计](H10k_Material_Spatial_Owner_Migration_Design.md)和[实施前瞬时 owner 补充](H10k_Instantaneous_Material_Flux_Addendum.md)执行；[实际 API](../api/H10k_Cch_Spatial_Material_Flux.md)。最终 native/H10b 未关闭，H11/H08a 未开始。

## 实际改动

固定2 production：现 `CchMaterialTransport.h/.cpp`；2既有 tests：phase/material 契约；现研究 driver/README。0新源码、target、tools。Core/FirstJet、phase/advector生产、D_I/C、dual/KEP生产及其两份测试、App/字典/PDF均保持。

同一文件内 `CchMaterialFlux` 唯一保存 `rho_l*Q_l+rho_g*Q_g` 的 canonical atoms，primal/dual port 只从这份 F 求和。它接受真实 Initial/Midpoint/End layer，无 dt/R/库存；`CchMaterialTransport` 保留原 Flux storage，计算原三层质量差、中点积分和实际 R。来源新增实际 advector stamp；Flux/Transport均核验原 momentumGeometry instance、layout/D_I、完整 clock 和 mesh。copies 与区间包装的三个率数组地址相同。

瞬时 owner 补充解决了原区间对象不能接受真实 Initial/End witness 的边界；没有重贴 snapshot 时钟、观察器另乘密度或第二材料通量算法。原设计保留，补充先于 typed Flux 实现。

phase/material 两份契约不再调用旧 evaluate；研究三处 helper 使用原 momentumGeometry 上实际 spatial factory，native advector 观察直接读同一实际 B，删除第二次拟合。两个新增只读 observer blocks 输出材料研究。逆向恢复三处 helper 并删除两个 observer 后，原 driver 全文相同；原 full/RK4/donor/gate 未改。

## 实际验收

标准 Release；最终必要 CTest **8/8**，包括architecture/自检/目录与phase/material/momentum/dual/KEP；最终文档检查3/3。既有 dual/KEP 回归仍属旧算法范围，不能称它们的新 F 已迁移。

phase 原几何/时钟/寿命/轴/收敛检查保留，face-constant 假设换为独立 physical polynomial 积分。material原309场景、6,780 cell/16,444 face检查保留；8个basis单面输入按实际 reconstructed polynomial 重算所有真实面及邻接，零面平均不再被误当作零材料率。新增264材料准备、42,240 atom核验，含12个零净 counterflow witness及136个自身尺度小相 witness，最大独立误差1.928e−16；复制、三真实时间层、不同actual stamp/instance/fullclock拒绝均通过。原phase 88场景/14,080 atom和8种斜切/多分量形状保持。

| 独立最终 CSV 核验 | 实测 |
| --- | ---: |
| XY/RZ uniform/2:1 full/selected布局 | 64 |
| actual phase factories / material preparations | 128 / 384 |
| canonical atom / primal+dual view / pressure cell行 | 685,800 / 946,860 / 108,120 |
| 内部primal / dual视图记录 | 204,180 / 633,240 |
| 独立相积分 / 材料率误差 | 2.304e−15 / 2.244e−15 |
| 同源phase closure / 几何重组F视图差 | 1.946e−16 / 0 |
| 独立三层体积 / 质量误差 | 2.068e−13 / 2.068e−13 |
| incidence/R/全域balance身份误差 | 6.669e−15 |
| 三真实层×三密度gas witness | 63行，误差0 |

每个布局保留 stationary plane 和正式 `advanceFirstJet` 的移动 circle/sphere，rho=1/1e3/1e6。独立脚本只读实际 CSV，按声明物理 MMS 闭式积分、解析相交和 circle-rectangle chord primitive 重算 phase/F/三层质量；按坐标重组 canonical atom、检查无gap/overlap和原数量；用真实 native owners重算selected内部反号和boundary/R/balance。全部56个预列 coarse selected tangent仍有限/同源/共享；只对其不适用的二维affine精确参考作原声明，没有新豁免。

原XY AMR n16 dof103的Initial/Midpoint/End各7段，三密度均消费真实层上的新Flux；液率严格0、气率及F=.4*A或0，signed outward sums为0，原三geometry/三advector identities及合法完整clock保持。它只证明新F，尚不能改判旧H09 quarter R。

384次材料准备累计1.385s（本机单次，不含先前phase/advector成本，不作性能门槛）。源Factory复用128份实际phase，没有第二轮速度/几何准备。

## 边界与失败记录

几何按规定平移或静止，Hydro advector为不同的MMS场，属于强制几何下的空间材料率/代数验证，**不是耦合PDE或R小余项验证**。实际最大|R|=442,287.79、最大relative R=0.99998759，作为有量纲/尺度诊断保留，不用R=0作接受条件。时间仍中点近似，一般时空精度归H08a；Core trace estimate不改称严格bounds。

首轮phase迁移在±1e−300下触发既有Core“exact product tail ... unsupported arithmetic scale”，未修改生产或容差。既有fixture把可用因果侧的明确 arithmetic rejection记为支持边界；最终352条此类拒绝，其余正常±流速/零/变号、小相正常验收。原始exit1和随后exit0均归档。后处理的首次C#编译/布局几何计数错误已修正独立脚本；没有据此更改数值输出或门槛。

原source链20份CSV字节保持；affine的4份中3份字节保持，`affine.csv` 3,688数值字段变化，最大scaled差3.290e−13。这是允许的atomic spatial material准备变化，字段/行/原门槛不变，未宣称全部字节相同。实际8/16内部affine四组合最大误差1.235e−13，原counterexample检查通过。n32与原完整H10b本轮未重跑；原28窗口/84gate与14失败历史保持。

启动前H10j 289 payload/101 canonical与j0/i2/f/e/d/c归档通过。1,044源成员保持，仅上述5源改变，1,039不变；所有cases配置保持。compiled sources、CCH library/exe/DLL/五测试二进制、实际退出/CSV、实施前设计和独立核验冻结于[H10k baseline](../../runs/baselines/H10k_20261004/README.md)。临时检查归`runs/checks/H10k/20261004_141203_693`；旧归档不变。

## 唯一下一任务

[H10l](H10l_Same_F_Dual_Material_Migration_TODO.md)，固定5 production/2既有tests：dual消费同一真实F，保留实际R/uR与材料KE，定义非零E_sp，删除quarter算法、最后旧evaluate和material legacy表达。当前旧入口只由dual/KEP两份契约消费，新spatial失败不回落旧路线。

材料硬约束仅唯一共享F、内部反号/实际AMR子面和、来源一致、KEP/skew；quarter、精确material commuting、R=0只诊断。Pure FirstJet完全独立，无库存/bubble。H10l后仍须联合收口最终native constant/affine与material KE，先内部affine再原H10b全部28窗口/84gate，H11/H08a等待。
