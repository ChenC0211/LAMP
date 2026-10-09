# H10l2 / 实施前 same-F dual 迁移设计

2026-10-04；上一轮 progress。启动复核 H10l1 的323 payload/119 current canonical、实际编译源码/runtime及九份前置归档通过；当前1,044 source成员。按[任务书](H10l2_Same_F_Dual_Material_Retry_TODO.md)仅做5 production/2既有dual/KEP tests，不进入下一数学任务。

## 迁移依据

逐文件比较确认本项5 production及2 tests与H10l候选启动前原文哈希相同；H10l1只修改actual advector及其phase/material tests和两个被动研究观察块。本次可从H10l最终失败候选提取已审阅的same-F消费与legacy退役改动，但必须对当前H10l1重新编译、复核范围并重新跑必要回归，不能沿用旧研究的通过数据。

原失败有两类：basis的零native gross假散度由H10l1的真实face-local前缀积分修复；closed流函数均值0并未给出pointwise不可穿透条件，由已验收typed law显式绑定。后一项只在现KEP closed fixture传入明确selected-domain boundary DOFs，不在生产按mean0、算例名或AMR位置推断墙面。指定actual polynomials、phase atoms、material dual segments均逐一核验严格0；未指定面的反向相流由原phase/material tests及迁移dual小相/counterflow tests保留。

## 唯一F与真实缺陷

对同一momentumGeometry每个segment直接消费 `material.dualFaceRates()`，`CchDualMaterial::faceRates()`返回原vector引用，不复制或重新求F。仅保存 `dt*F` 积分；真实三层质量来自原FirstJet/H09a。移除原quarter/SPD/local extension及其硬交换断言。

真实incidence给出 `q_d=sum(s*F)`、`r_d=delta_m_d+dt*q_d`。原D_I符号给出pressure q，原真实overlapFraction定义T，仅用于诊断。显式保存 `E_sp=q_d-T*q_p`（率）、`G=delta_m_d-T*delta_m_p`（质量），核验 `r_d=T*r_p+G+dt*E_sp`。E_sp/G/R不要求零；归一化使用真实gross dual+mapped work，不能用抵消后诊断自身作舍入尺度。原consistencyImpulse和materialKineticChange数值保持逐字相同；KEP生产文件不改。

完整三层source/clock、实际advector stamp、momentumGeometry instance、原layout/D_I与密度一致性保留。先检查原具体来源以保留拒绝原因，再检查material与mass同一momentumGeometry。copies共享原owner，重网格失效/不同geometry/时钟/nonfinite仍拒绝。静态law无新增clock；本项无接受事务。

无调用者后删除 `CchPhaseRate::evaluate` 及optional legacy存储；现prepareSpatial数值逐字保持（仅Data构造从optional改直接owner），原MaterialFlux数值owner不改，Transport只删除legacy分支。现公开spatial接口保持，避免新增第6production文件。

## 验证与交付边界

先标准Release并运行8项必要CTest。原208单面basis/linearity、四closed场景（显式law）及原其他dual/KEP场景均必须实际完成；任何失败保留日志，不忽略或放宽gate。原136 tiny-phase own-scale见证、counterflow、真实三层跨面与uR/materialKE/source拒绝保留。两份允许tests内必要新断言可以加强，不能修改已验收H10l1 tests。

既有research driver只追加两个H10l2 observer块；删除这两块必须逐字恢复本项启动时driver（含H10l1观察器）。复用原64布局、56 coarse-selected tangent缺口、三密度与三族实际GALS：静止平面/MMS、移动圆球/MMS为强制几何，移动圆球/同向0.4平移为coupled制造族。独立CSV重算相交物理积分、三层质量、真实half-cell overlap、q/R/E_sp/G及uR身份；gas dof103三层21段×三密度保留。coupled族也仅验证实际身份，不凭有限点宣称R收敛。

重跑原n8/16内部affine，门槛1e−10、原外边界数据保留；此次不跑完整H10b或n32，不改full/RK4/donor/gates。若最终operator受新F影响失败，冻结反例并另拆最终operator任务，不绕回segment trace补丁。

最终2tests/5production与1044成员预算、原空间phase/Flux及uR数值、App/config/PDF/CMake保持核验；源码/真实exe/lib/DLL、失败/最终退出、独立审计与flat manifest绑定。只有本项PASS后才在下一任务联合final native constant/affine spatial moments与material-KE/skew，内部通过后原H10b28窗口/84gate，H11/H08a等待。
