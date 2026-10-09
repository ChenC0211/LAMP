# H10l1 / 实际速度场闭合与显式不可穿透来源验收

2026-10-04，**PASS，限定本项范围**。本轮仅完成 [H10l1任务书](H10l1_Spatial_Advector_Closure_And_Boundary_TODO.md)，goal 保持 active，分类 progress。[实施前设计](H10l1_Spatial_Advector_Closure_And_Boundary_Design.md)原文保持。下一唯一 [H10l2](H10l2_Same_F_Dual_Material_Retry_TODO.md) 尚未实施；H10l 的失败归档不改判，当前 H09 仍旧 quarter，完整 H10/H11/H08a 尚未通过。

## 修复及来源

仅修改 `CchSpatialAdvector.h/.cpp` 两份 production。原 local one-ring 权重/支持、原 D_I、真实 AMR faces/atoms 和 cell lift 保持；同一 face 使用局部体积坐标的零均值切向基函数，前缀积分为 `w*(b+k*(w-L)/2)`。完整 face 的 slope 贡献精确为零，避免先分段积出正负大数再相消。实际 field divergence 仍由实际积分计算，原 D_I 是独立参考；生产门槛仍为 `1e-12*native gross`，去掉经验 1e−300 floor，gross=0 时要求精确零。

新增只读 `CchAdvectorBoundary::impermeable`，绑定原 momentumGeometry instance 与调用者明确列出的 selected-domain boundary DOFs。指定 native normal mean 必须严格零，该 face 的整个 normal polynomial 为零；其他面包括 mean=0 的内部/一般边界继续使用原重构。law 没有第二份 velocity/clock/phase 或接受状态，actual B 保留其不可变共享副本、原完整 clock 和独立 stamp。全部 native/cell velocity、primal/dual atoms 与 phase/material F 消费同一实际场，没有在 H09 后置压通量。接口及拒绝规则见 [API](../api/H10l1_Cch_Advector_Boundary_Closure.md)。

## 原反例定位

先只增加实际 cell、field divergence、原 D_I 与 native gross 的拒绝上下文，不修改数值；保存该诊断版 source/exe/lib/DLL 后，运行 n4 XY/RZ、uniform/2:1 的全部 **208** 非轴单面 basis。**88** 次实际 evaluate 拒绝；独立脚本按原片段算术重算 **4,912 cells**，识别 **226** 个失败 cell，全部位于 **4,552** 个零 native gross cell 中。脚本预测的失败 scene 集合与实际拒绝完全相同。

原最大误差/实际 trace L1 为 **6.890e−16**；高精度 decimal 算术的零 gross divergence 最大 **4.3e−27**，稳定全 face mean 与 D_I 的最大尺度误差 **2.180e−16**。这证明该组失败源自零均值切向积分相消的舍入，而非依据 gate 的 floor 猜测放宽容差。典型 XY AMR u3/cell5：旧 actual divergence=−5.551e−17，原 D_I=0，native gross=0。

修改后同一 **208/208** basis 全部接受。既有 phase 契约还运行 **236** 次 basis/组合/零场/近相消/tiny/large 更新，核验13,528个物理face、5,528个cell和4,664个零gross cell。独立 face GL3、原 incidence、actual cellVelocity 在真实piece内的导数与多项式线性核验，最大误差分别 **1.352e−16、2.220e−16、1.601e−15、6.109e−16**。tiny/large 输入尺度1e−200/1e200；全零精确零。新增边界/source/instance/copy/lifetime/clock/regrid/非零denormal均值拒绝检查 **40** 项通过。

## 显式边界与反向相流

原 closed streamfunction 在零 native boundary mean 下仍有变号 trace；此次不把零均值自动当作不可穿透条件。研究覆盖 XY/RZ、uniform/2:1、n4/8/16/32、full/left-selected 共 **32布局**，绑定/未绑定共 **64 actual revisions**。独立解析 polynomial 物理积分、流函数输入均值、原 incidence、actual cellVelocity 散度及原真实 atom/view 分割核验结果：

| 边界研究 | 实测 |
| --- | ---: |
| atom / face / cell / segment 行 | 140,280 / 46,980 / 22,440 / 143,760 |
| 指定 boundary atoms / dual subsegments 精确零 | 3,480 / 5,220 |
| 未绑定对照的非零 boundary atoms | 3,000 |
| 独立 F / physical mean 最大尺度误差 | 4.348e−16 / 4.324e−16 |
| 原 D_I / actual cellVelocity divergence | 2.479e−16 / 1.571e−15 |
| 内部共享/实际子面和身份 | 4.348e−16 |

既有 material 契约追加 pure-liquid/平面两种实际 FirstJet、三密度1/1e3/1e6、full/selected 的 **48** 次准备：1,392 boundary atoms、2,088 boundary segments 精确零，1,200未绑定非零对照，4,488内部段；独立 phase/F、真实 incidence/R 和 global closed boundary 身份最大 **4.719e−15**。

额外四种布局只指定一个零均值边界 face，其他面保留 **15** 个零均值非零切向迹，真实平面 causal phase 保留 **4** 个液气反向流见证。XY AMR 的正负相流分属真实共同侧的两个子面，求和 `Q_l=0.00390625`、`Q_g=−0.0039062500000000009`、`F=3906.24609375`（rho=1e6/1）。逐 atom 独立积分和真实子面求和同时核验，没有按 AMR 位置修改生产算法。

原 closed 研究场本身没有内部 mean=0 变号 DOFs（该统计为0）；不能把它冒称 counterflow 覆盖，该要求由上述 selective-law 契约独立满足。原小相自身尺度仍保留：phase 88准备/14,080 atoms/68 tiny witnesses及8形状/1,280 atoms；material **264准备/42,240 atoms/136小相见证**，最大独立误差 **1.928e−16**。原352个 causal scale 拒绝仍明确执行。

## 既有研究和工程回归

64布局/128 actual advector 更新仍保留72 rank-deficit、56 coarse-selected tangent 缺口，无新增豁免。独立核验77,540 face、228,600 atom、36,040 cell：物理均值3.084e−16、affine2.220e−16、atom系数6.281e−14、闭式积分3.156e−16、原divergence1.188e−14，均在原1e−12门槛内。

原材料研究 **64布局/128 phase factories/384材料准备** 独立核验685,800 atoms、946,860 views、108,120 cells。最大 phase/F 误差 **2.304e−15 / 2.244e−15**，解析三层体积/质量 **2.068e−13**，incidence/balance **6.669e−15**，三层×三密度 pure-gas63行误差0。最大实际 |R|=442,287.7885、相对R=0.9999875868为规定几何与MMS速度不同的诊断，**不表示 R 小余项或耦合 PDE 收敛**。

原内部 n8/16 affine 四组合、共八组运行通过，最大内部误差 **1.235e−13**（原gate1e−10）；外边界误差原样保存，不计入内部通过结论。没有运行 n32 affine 或原完整 H10b；原28窗口/84 gates及失败历史保持。删除 driver 两段 H10l1 observer 后，与启动前原文相同；原full/RK4/donor/gates未改。

标准 Release 六 targets 构建成功，最终必要 CTest **8/8**：architecture、architecture_self_test、repository_layout、phase、material、momentum、dual、KEP。当前旧 dual/KEP 测试只证明原 H10k 范围，不能替代下一项 same-F 迁移验收。最终文档检查另 **3/3**。source成员1,044不变，其中1,039哈希不变；仅2 production、2 tests、1研究driver变化。无新源码/targets/tools；Core/FirstJet、phase/material/dual/KEP数值、App/字典/PDF保持。

## 失败记录与归档

新增测试初版选择domain corner的零迹plane，被原Core正确拒绝；改为已有内部coincident plane fixture。XY AMR 初版只找单segment的零净反向流而漏掉真实共同侧，保留失败后补充真实子面和，逐atom断言保持。RZ单相内变号多项式积分应为0，独立GL3留下4.337e−19舍入；初版用signed净值归一化而误失败，改为实施前设计规定的实际系数work界（无floor，1e−12及原tiny-phase own-scale门槛不变）。这些修改只在两份允许测试内，生产phase/material数值未变。

首次完整必要回归7/8时的日志与source/runtime、candidate1/3/4失败测试及candidate5成功退出均保留；最终重新标准构建和8项回归。原closure诊断版另存source/runtime；最终五个研究运行绑定同一个研究exe（SHA256前缀`EC3B8C1B7B864F14`），最终31份compiled source/8 binaries与当前canonical核验一致。临时PowerShell执行策略及UTF-8无BOM解析问题仅改变调用/脚本编码，不改数值。

实际CSV、失败/最终日志、独立核验脚本、source/runtime和flat SHA256 manifest见 [H10l1冻结目录](../../runs/baselines/H10l1_20261004/README.md)；检查目录 `runs/checks/H10l1/20261004_193638_664`。启动前H10l的328 payload/115 canonical及全部前置归档复核通过；关闭时旧归档仍按原manifest校验。

下一项仅H10l2：5 production/2 dual/KEP tests迁移同一材料F并删除quarter，closed fixture显式消费已验收law。其后联合最终native constant/affine空间矩与material-KE/skew，先内部反例再原完整H10b；H11/H08a等待。四材料硬约束继续限定为唯一共享F、内部反号/实际AMR子面和、来源一致、KEP/skew；Pure FirstJet独立，无库存/bubble。
