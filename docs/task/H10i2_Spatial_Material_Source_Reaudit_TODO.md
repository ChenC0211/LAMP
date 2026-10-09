# H10i2 / 同源空间材料设计与来源链独立复核

状态：**设计/来源 PASS，engineering PASS**，2026-10-04；[报告](H10i2_Spatial_Material_Source_Reaudit_Report.md)。16实际来源、864face/864dual/1472overlap、39原CSV、7/7必要回归及独立公式核验通过；0 production/0 tests。前置f210/e268/d189/c175归档与f87 canonical已核验，[实施前复核方法](H10i2_Spatial_Material_Source_Reaudit_Design.md)保持。下一唯一任务[H10j0 actual advector准备](H10j0_Spatial_Advector_Reconstruction_TODO.md)，然后H10j；[现行API/固定预算](../api/H10i2_Spatial_Material_Owner_Contract.md)。原H10i设计FAIL与原H10b 14失败保留，新材料F/最终native未实施。

目标：从当前实际 FirstJet publication 经 typed trace/cache、H07 spatial kind、RegionSource/真实 dual mass，独立复核全部来源和能力；按现行接口重审 H10i 的空间 owner 数学设计，明确 H10j 可以实施的范围与固定预算。

允许修改：本 TODO/独立设计复核/报告，H10i/H10j 现行设计与 API 状态说明、现行路线图/架构入口；如必要，仅既有 cch_convection research driver/README 的被动证据。0 production / 0 tests / 0 新 tools/targets。开始前核验 H10i1f flat manifest、canonical、前置 archive，保存源码/输入 catalog 和原数值入口；临时证据归 runs/checks/H10i2，结束归 runs/baselines。

禁止修改：Core geometry、Interface primary/Q4/transport/publication/clock、H07–H10 production、材料 F/R/uR/KEP、压力/散度、冻结 PDF/App/字典、原研究容限。无库存/bubble、质量修复、几何扰动、half 分配、AMR 特例、解析输入或数值 fallback。本任务不构造新 phase rate，也不迁移 material consumer。

数学约束：材料硬约束仅唯一共享 F、内部面反号共享/真实 AMR 子面和、来源一致、KEP/skew。保留真实质量与真实 R/uR；quarter、精确 pressure→dual commuting、R=0 只诊断。固定 mesh 的实际定义为 `R_d=Delta_m_d/dt+q_d=T*R_p+G/dt+E_sp`，其中 `E_sp=q_d-T*q_p`；不能重新强制 E_sp/G 为零。Pure FirstJet 完全独立。其他模块原压力伴随/真实物理测度等约束保持。

实现要求：逐项核查原 Snapshot/revision/full time/actual domain/AMR subface，区分实际 publication 与独立 seed、几何 spatial average 与 causal material trace。内部 zero face 的新率必须共用实际 advector 的因果流入侧；外部/selected zero 流入缺少正式相数据应拒绝，不能用已知内侧或对称平均填补。ordinary、point-only incidence、side support 和完整 region 证书不得混用。审查 canonical atom overlay、physical-mean one-ring advector、RZ 体积坐标/轴正则、positive phase quadrature、生命周期及工作预算，尤其确认仅 source 完整不足以证明运动学重构或 PDE 正确。

T0：Release architecture/self-test/layout 与必要原短契约，保护源码/输入/App/PDF 哈希、原内部 affine 与原完整 H10b 的固定门槛；不得为独立复核扩大测试或生产文件预算。

T1：独立重新检查 i1f 保存的原 16 实际 publication/H07/dual mass、两方向 clip/region/XY–RZ 物理参考及侧 availability；原 64 普通布局、157,810 查询、56 极粗 selected 切向未识别行保留。核查 missing exterior、ordinary→side、多 component、source/time/domain、FP、lifetime、拒绝与失败事务的实际记录。不能将未认证 RZ 域角点当作支持范围，不能把非零 causal rate 拒绝误报为几何来源失败。必要追加观察只能使用原实际来源和独立参考。

T2：重审 H10h pure-gas dof103 的三层 21 实际边段与 H10j 接口；列出唯一 owner 的输入/来源、可识别 affine 与信息不足的区别、共享/AMR 求和及 R/空间误差的独立 gate。重新确认 H10j <=5 production/<=2 tests 的固定文件清单；若现 2 production/1 test 预算不足，先拆单一接口任务，不能隐含 Core/Interface/material 改造。身份 <=1e−12 gross、常速度 <=1e−10；最终内部 affine 精确/二阶先通过，随后原完整 H10b 的空间阶 >=1.7、最细 L2<=0.05、Linf下降及全部原 KE/时间/donor gate 保持。

失败规则：来源/API/数学设计不闭合即 FAIL，冻结具体缺口后登记唯一下一接口任务；只在真实外部依赖时 BLOCKED。不能以解除材料 commuting 约束替代 PDE 精度验收。

交付：独立H10i2裁决、旧到现接口核对、支持/不支持/未实施项、来源和source scope证据、后续固定预算及flat manifest。本次复核按“原两文件预算不足先拆单一接口”要求，明确前置H10j0 actual advector几何复用，随后原H10j。H10k/H10l与最终native/PDE原H10b分别验收；H11/H08a仍等待，不在本任务执行。
