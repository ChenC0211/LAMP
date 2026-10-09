# I12p3f / 同源曲率接受态时序诊断

目标：关闭 I12p3r1 中输出规范 §9 的 RMS/unavailable count 缺口。只读当前 accepted publication，补齐全部接受态曲率统计，不改变数值方法。

允许修改：InterfaceFirstJetOutput.cpp（1 production）、interface_first_jet_output_contract.ps1（1 test），本项 API/TODO/报告、当前状态入口及 runs 证据脚本。

禁止修改：Core/Modules/App producer/typed config/reader、几何/曲率/求积/AMR 算法、输入/容限、历史报告/归档/PDF；不实施性能优化或 Hydro。

数学约束：kappa 使用 publication 原 physical curvature（XY meridional；RZ meridional+azimuthal）。RMS=sqrt(sum physicalWeight*kappa²/sum physicalWeight)，权重来自同一个 publication.quadrature；不使用 VTP 显示点替代求解几何。按原 report.maxAbsCurvature 缩放平方避免无必要的大数平方，无物理状态/acceptance/输出调度反馈。

实现要求：conservation.csv 追加 kappa_rms、curvature_quadrature_samples、curvature_unavailable_samples、curvature_measure。空界面 RMS 为 NA、measure=0、counts=0；任何不可用采样显式计数，RMS 为 NA，不能用零或部分平均冒充完整统计。当前 publication 已拒绝非法曲率/权重，实际接受态预期 unavailable=0。保持既有字段次序/内容、stdout、VTK/PVD/frames。计时开销归入原 IO stage。

验收：T0 39 targets/no new warning/architecture；T1 18 频率组每个初始/接受态有上述字段且 schedule-independent，样本数正、measure正、unavailable=0、0<=RMS<=kappa_max；Q2 XY/RZ 静止/平移 RMS 对照解析 1/R、2/R。原输出/失败/原子重试、p3e 15 精确范围/6 合法边界保持。T2 实际 12 App smoke/stdout/既有series/48帧VTK全部与 p3e 比较；数学 kernel/38其它 runtime 保持，复用 p3d 长窗口/收敛/规模证据，不声称本项重跑这些实验。

失败规则：另拟合曲率、隐藏无效采样、改变原数值/容限/合法输入、schedule改变诊断或解即 FAIL。内部问题不记外部 BLOCKED。

交付：1 production/1 test、同源具名统计与局限、实际验证/身份归档。随后只读最终 Gate I 复核；不自动宣告低黏耦合或 broad goal 完成。

执行结果：**PASS**，见[报告](I12p3f_FirstJet_Curvature_Diagnostics_Report.md)。90 个频率实验初始/接受态统计、30实际组/150行旧字段、96VTK字节、原App/范围/失败/重试保持；1 production/1 test，无数值 kernel 修改。下一项 I12p3r2 只读最终复核。
