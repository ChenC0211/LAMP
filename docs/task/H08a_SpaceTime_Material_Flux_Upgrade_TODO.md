# H08a / 同源时空相积分材料通量

状态：2026-10-06，**PASS：本任务规定的同源时空积分、材料一致性与小余项收敛范围**，见[报告](H08a_SpaceTime_Material_Flux_Report.md)。r11 基础设施保持封口，必要 FirstJet 运动消费已并入本任务完成；未另立 source/Jacobian/低秩精修任务。Pure FirstJet 独立，真实 R/E_sp/G/uR 与 KEP/skew 保留，历史 FAIL 不改判。下一单任务为 H11 controlled dissipation，本轮未启动；完整流体 App、全耦合二阶及效率尚未验收。

目标：保持Pure FirstJet/Hermite–GALS几何完全独立，在Hydro侧从同一实际步内几何和特征运动积分共享face的液/气相输运，构造唯一材料通量F_rho，使材料连续性缺陷R成为可收敛的小离散余项。

允许修改：沿用 r11 已部署的同源实际场、characteristic 和完整 clock，在 Hydro owner 内实现步内几何消费、共享相通量积分及 primal/dual 唯一材料积分；具体文件预算在实现前设计中固定。必要的 FirstJet 运动消费接入纳入积分目标，保持 Interface 不反向依赖 Hydro。旧 H08a1/source 准备顺序为历史，不能继续以精修基础设施替代本任务。完整积分与 R 收敛验收归本任务。

禁止修改：Interface primary、几何发布/曲率、库存、bubble、逐cell体积拟合、材料修复或另一套rho*u通量。不以积分通量更新一个守恒质量状态；真实材料质量仍直接来自唯一GALS几何。

数学约束：液气相积分共享同一时空几何，材料积分为rho_l*Q_l+rho_g*Q_g；内部共享face和2:1真实subfaces保持唯一来源与代数抵消。逐cell和dual的R必须从真实几何质量变化及该唯一积分重新核算，不能设为0。继续现有uR、|u|²R/2及KEP相容对流，不反向约束界面。

实现要求：先定义同源步内几何/characteristic来源、积分精度预算、时间层和拒绝条件。仅三个独立snapshot做插值或midpoint面积乘dt不能自动称为同源时空相积分；不得重新启用历史守恒Interface driver作为生产owner。

实施进度（2026-10-06）：已增加中立已解析运动的 FirstJet 消费与 `CchMaterialTransport::prepareSpaceTime`，mid/end/时间求积节点均从同一原始 FirstJet 和完整 source slab 生成；primal/dual 共同消费已积分 atom 材料通量，原 RK2 和中点入口保留原语义。24 项受影响短检查均已通过；包括原 12 组移动平面及真实粗细边界非零相通量切换 6 组（24 个共享子面检查）、8 个端点、4 个真实 transfer/rebuild/resume、8448 个同源稠密场比较、3268 个实际内部骨架轨迹及独立材料/KE 恒等式。证据为 `runs/checks/H08a/20261005_232800/current_short_verification.json`。原材料断言和 CMake 门槛已用任务前哈希及原文本复核，未放宽。

最新部署已完成全部 28 个曲面时空窗口，独立复核通过；四布局的空间通量/几何/primal/dual 缺陷率误差下降，同网格时间加密的真实逐步缺陷率历史与终态几何自误差也下降。最新 H10b 完成 28 窗口/15116 步/84 PASS，四旧 donor 仍 FAIL，8 表与 r11 逐字相同。总复核见 `runs/checks/H08a/20261005_232800/current_complete_verification.json`；最新运行是 `runs/benchmarks__cch_spacetime_material/20261006_012544_866_H08a_corner_547e7f08`。原 Green 算术首轮仅作历史比较。实施范围、原失败和预设门槛见[设计](H08a_SpaceTime_Material_Flux_Implementation_Design.md)。时间自收敛以同网格 8 步为参考，固定 h 的绝对 R 可趋于非零空间余项，不能据此宣称完整耦合二阶或 R=0。

积分对象直接取真实 canonical atom 上的液气相时空积分 `Q_phase = integral_dt integral_atom chi_phase(phi_h(x,t)) B(x,t)·n dS dt`。Hydro 只在此唯一结果上做一次密度加权 `F_rho = rho_l Q_l + rho_g Q_g`，primal/dual 视图复用并按既有 incidence 求和；对流所需平均率来自同一积分除以 dt。随后以真实几何质量差加该积分核算 R，继续既有 uR 和材料动能项。这是本积分任务的接入顺序，不增加守恒库存、交换恒等式或新的 source 精修先决任务。

验收：启动前在实现前设计中固定T0/T1/T2数值门槛及完整研究窗口。至少包括真实moving/crossing相面、XY/RZ uniform/2:1、共享face求和、光滑与穿越情况下R的空间/时间收敛、密度比1/1e3/1e6、既有常速度/uR/KE预算、真实regrid来源重建及全部前置回归。不得以严格R=0为目标，不能在结果出来后改阈值。

R 收敛证据须区分步积分缺陷与缺陷率，并使用固定物理观察窗口和明确的质量/体积尺度；不能仅因 dt 变小、单步积分缺陷自然减小就宣称收敛。分别说明时间积分误差、空间几何/重构余项与舍入平台；这些是误差诊断，不变成 R=0 或 quarter commuting 等附加硬约束。

失败规则：不同源、第二材料通量、不可收敛或来源拒绝失败判FAIL；无旧midpoint或守恒库存fallback。真正外部依赖阻断才BLOCKED。

交付：唯一积分API/实现、R及材料功的独立完整证据、明确误差和剩余缺陷。H08历史报告保持原裁决，不把旧中点研究改标为时空积分。

## H08/H09 解除过约束（2026-10-04 用户授权）

材料耦合硬约束仅四项：唯一共享材料通量、内部面反号共享/AMR 子面守恒、来源一致、KEP/skew。旧 quarter 分配、`q_dual=T*q_primal` 和无空间交换差的缺陷关系不再强制。固定 mesh 的实际缺陷显式包含 `E_sp=q_dual-T*q_primal`，与几何质量转移差分别诊断；不要求严格 R=0，精度和 R 收敛仍需验收。空间 owner 沿用 H10i 后续通过的设计，不能恢复旧过约束。

**基础设施收口（2026-10-05 用户授权）：** r11 是 common regular source/characteristic 的最后一道关口。除非出现可复现且直接阻断 H08a 的数学反例，不再拆 source/Jacobian/低秩精修任务，不新增硬约束、特殊 case，也不放宽原门槛。r11 通过后立即回到 H08a 的真正同源时空相通量积分；必要的 FirstJet 运动消费接入属于该积分任务，不再单列桥接先决任务。目标仅为 R 成为可收敛的小离散余项，保持 Pure FirstJet 独立、无 M_K/bubble、无 R=0 或 quarter commuting 强制。材料硬约束仍仅为唯一共享材料通量、内部面反号共享/AMR 子面守恒、来源一致和 KEP/skew。
