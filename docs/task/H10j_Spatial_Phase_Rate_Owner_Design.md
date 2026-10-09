# H10j / 实施前同源空间相率设计

2026-10-04。前一轮为 progress，H10j0 193 payload/97 canonical 与前置归档已独立复核。按[任务书](H10j_Spatial_Phase_Rate_Owner_TODO.md)只改现 CchPhaseGeometry header/cpp 和一个现 phase contract；research只观察 actual factory。旧 evaluate 与 H08/H09/KEP、Core/FirstJet/CCH advector、字典、App、PDF保持。

prepareSpatial 比较原 phase/advector/layout/mesh/full clock 和唯一 D_I，保留原 immutable advector与publication。H10j0 canonical atom 是唯一率 owner，primal/dual views仅作真实索引求和。报告原rank/deficit，不重新拟合V。空间owner getters在legacy evaluate结果上明确拒绝，不伪造advector来源。

每个 atom 用原 owning geometry view alias 调用 NativeQuarticTrace。ordinary ranges分别积分 signed normal polynomial，液/气单独积累。正长度zero范围按 actual velocity 的正部分取minus侧，负部分取plus侧；每侧的相ranges独立消费，不把两侧root ordinal或opaque chart坐标合并。Core已有 integratePositiveUnitPolynomialRange 提供正部分和误差估计；密度多项式只由同一actual normal coefficients与正XY/RZ面积组成。offset/width仍为原输入chart，绝不重造全局root点。

selected availability属于Hydro实际域：即使原scalar覆盖域外，zero流入侧不在layout中也拒绝。缺少侧时对该zero interval检查对应正部分，无法证明无流入即拒绝；不读空间平均、clip迹或邻居材料。轴面积零且radial velocity零得到严格零率。净面速度/率为零不代表normal polynomial恒零。

源查询的每个尚有宽度的非point-only事件用当前fresh inventory作同源refinement；每轮重算两相率，最多4轮，只是factory有限工作量、不改变Core隔离预算。停止要求两相估计变化分别<=1e-13各自工作尺度；不能认证/无可用精度则拒绝。positive primitive的请求误差为1e-13原range物理工作尺度，报告真实error；terminal ranges和相邻两轮变化是误差估计，不能称严格bounds。完整atom phase closure与原actual polynomial total按1e-12 gross核验，不修正两相以强制closure。

验收先独立closed-form平面/斜切/圆球/小相多项式积分，再64真实full/selected布局、16重合源、positive/negative/zero/counterflow、完整clock/来源/失败事务与生命周期。外部zero流入为明确拒绝范围。原H10h pure-gas dof103的三层21段必须液率严格0；只观察新相率，不生成材料F/R。必要phase/material/dual/KEP等回归和旧研究数值保持。PASS只关闭H10j，下一H10k固定4 production/2 tests迁移入口与来源stamp，再H10l/final native/H10b；H11/H08a等待。
