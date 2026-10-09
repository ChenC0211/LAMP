# H10j / 同源空间相通量 prepared owner

状态：**PASS / engineering PASS**，2026-10-04，见[报告](H10j_Spatial_Phase_Rate_Owner_Report.md)和[实际API](../api/H10j_Cch_Spatial_Phase_Rate.md)。2 production/1既有test，64布局/128factory及8/8必要回归通过；只关闭空间phase owner。下一唯一[H10k](H10k_Material_Spatial_Owner_Migration_TODO.md)，材料F与最终native未迁移，H10/H11/H08a顺序保持。

目标：在 Hydro 侧用同一 FirstJet 和实际 compatible advector，准备 primal/dual 的同源空间相通量视图，使后续唯一 F_rho owner 可直接消费真实相区间，不再以材料散度/面积分配代替材料积分。只实现这一个事实。

允许修改：`include/Modules/Hydrodynamics/CchPhaseGeometry.h`、`src/Modules/Hydrodynamics/CchPhaseGeometry.cpp`（2 production）；现 `tests/contracts/hydrodynamics/cch_phase_geometry_contract.cpp`（1 test）；既有 `cases/benchmarks/cch_convection/driver.cpp`/README 和对应设计/API/报告。使用原 CMake owner 与标准 target，不新增工具或变体文件。临时目录 `runs/checks/H10j/<时间戳>`，关闭证据归 baselines。

禁止修改：Interface primary/出版/Q4、库存/bubble/质量修复、D_I/C、黏性/电场、材料密度定义、旧 H08/H09 F/R/uR、旧 App/字典/PDF、原 H10b gate。无独立高阶全局状态、AMR 特例、解析输入修正、界面重构或隐式 fallback。本阶段不把新 phase-rate 拿去替代未迁移的 material owner。

数学约束：材料耦合硬约束仅唯一共享 F_rho、内部面反号共享/真实 AMR 子面求和、来源一致、KEP/skew；本阶段消费H10j0已验收的实际advector构造同源phase-rate，不要求材料commuting/quarter/R=0。H10j0的运动学总量`Q=A*u`及`div B=(D_I*u)/V`保持，不重新强制pressure材料散度分配给dual。

实现要求：

- 现CchPhaseRate增加`prepareSpatial(mesh, phase, actualSpatialAdvector)`，直接比较advector持有的geometry/full clock、原layout/D_I/mesh，保留同一momentumGeometry instance和actual prepared owner。新增只读`spatialAdvector/momentumGeometry/dualFaceRates/spatialReport/advectorIdentity`等必要port；不再接裸速度另造B，不裸露mutable cache/率setter/snapshot shared_ptr。
- 原 evaluate 暂为当前消费者的历史面常值入口；新 factory 不调用旧材料延拓，也不在失败后转入 evaluate。下一 H10k 完成入口迁移并删除旧空间算法；本任务不声称 production 路线已替换。禁止将新 factory 的 API PASS 冒称 H08/H09 或 H10 PASS。
- one-ring physical-mean SVD、transfinite/zero-normal lift及native同点B均直接消费H10j0；本任务不重做重构或更改权重。RZ轴正则、未识别tangent及实际divergence scope保持，不能在phase积分里补另一算法。
- 复用H10j0的canonical atom/有向primal-dual视图（包含kinematic piece断点）；每个正向atom只作一份同源NativeQuarticTrace查询。使用publication owning alias，不复制scalar或穿透Detail。ordinary/zero side与实际速度正负分片使用原endpoint-relative范围/shifted polynomial；root的opaque chart坐标不能作为global atom参数。
- 正物理相权重，sign 来自 actual velocity。候选 trace 的三点 Gauss精度范围按 H10i；细相不可丢弃/clip/吸收，轴零面积率严格0。非有限/分解/几何query失败明确拒绝。
- 原CchAdvectorIdentity来自H10j0实际evaluate，本factory只保留它；下游H10k加入CchMaterialSources。所有点查询消费同一advector，不制造接受状态或H14事务证明。
- 不从u_f=0或total Q=0推断各相率为零；actual polynomial可以变号。两相分别积累signed正物理求积，不能total减大相丢小相。内部zero取actual causal inflow；外部/selected zero流入缺正式相数据拒绝，无空间平均/内部迹fallback。
- phase准备和更新按原actual atom规模记录工作量，复用H10j0几何storage/SVD；一次registry排序成本单列，禁止持续膨胀global DOF。

验收：T0 标准 Release、architecture/self-test/repository layout、原 App/输入/PDF哈希保持。T1 新 spatial factory 来源拒绝、同面反号/AMR atomic求和/physical closure、原 D_I total恢复和 volume-coordinate divergence均≤1e−12 gross；所有原时间/identity/轴/正负与零速度/小相/regrid拒绝场景保留。

16 full布局（XY/RZ uniform/2:1 n4/8/16/32）及原三类 selected：分别记录识别/不识别 tangent functional。H10i实测38,770面中72 rank deficient、56 tangent未识别（仅n4 selected2/3），不得把这56行从共享/有限/来源检查删掉，也不得把它们标成二维affine精确。full16和n≥8 selected必须再现可识别的affine，≤1e−12；极粗信息不足者按同一个最小范数构造验代数自洽并明确精度限制。

独立解析测试必须区分面平均与点样本：XY全affine、RZ `u_r=a*r,u_z=b+c*z`及可识别轴向even-r² physical-average family；带真实斜切/曲面/小相时，phase polynomial积分与独立closed-form或shifted moment核验≤1e−12相尺度。不得用新kernel自身做参考；RZ point/mean差异按O(h²)报告，原H10b点参考不改。

真实Pure FirstJet三层移动/crossing/regrid来源、H10h pure gas dof103全部三层21边段、正物理XY/RZ相测度、同物理段与父/子段moment加回、prepared复制和source寿命验证；测试新 spatial实际输出，无密度库存或材料时间update。原 phase/material/dual/KEP回归全部记录。本任务前置失败时不重跑完整H10b；不进入H11/H08a。

失败规则：任何空间构造/可识别矩/真实source/共享/细相/正规性或必要回归失败均FAIL；冻结诊断并修正或拆接口任务，无旧evaluate fallback。极粗selected不可识别是已预先列出的精度限制，不能在试验后借此豁免full/n≥8精度。只有真实外部依赖才BLOCKED。

交付：实现/API、实际new factory短契约/源段研究结果、所有不支持或未验证项、明确候选/旧consumer状态、平面manifest及protected核验；下一H10k单独迁移phase/material入口。H10k后H10l迁移唯一dual F/R/E_sp，再独立收口最终native PDE和原H10b，不能仅以材料输入改善关闭H10。

补充正式零迹契约：[H10i1设计](H10i1_Coincident_Phase_Trace_Clip_Contract_Design.md)。新factory只用同一typed NativeQuarticTrace的两侧迹与实际advector，按a正负分片取因果流入相；clip内部迹、无速度对称诊断测度均不得代替材料率。外部/selected零迹流入无正式相数据须拒绝，axis率严格0。完成后验原16重合场景、正负零速度及变号、normal反向/子面求和、wrong source与生命周期，并用独立平移右时间导数核验；不得恢复half分配、R=0或库存。H10i1f已迁移H07 trace kind与legacy拒绝，H10j预算只实施spatial owner，不隐含Core/Interface修改。

H10i1a正式Core NativeQuarticTrace已验收；H10i1b–f实际clip/volume/publication/H07/dual mass链通过H10i2独立复核。本工厂已消费H10j0 actual advector并通过实际验收。Core point-only value incidence仅保留endpoint来源，不是一条正长度零曲线；ordinary与side事件不能从rounded global参数合并。

实施后调用库存更正退役时序：旧evaluate仍由四份phase/material/dual/KEP契约和既有driver消费，一次删除超过原2-test预算。下一H10k固定2 production/2既有tests，先迁移唯一材料owner/source及phase/material fixtures；H10l再迁移dual/KEP并删除最后旧入口/quarter数值算法。新factory无旧fallback，不能把阶段保留称为全面迁移。详见报告与下一任务书。
