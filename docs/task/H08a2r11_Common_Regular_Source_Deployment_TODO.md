# H08a2r11 / 同一正则来源同步部署

状态：**第二次恢复 PASS**，见[本次部署报告](H08a2r11_Common_Regular_Source_Deployment_Resume2_Report.md)。原两次 FAIL 保持历史裁决。本次四 production / 两 test 已部署并通过全部原验收，不另拆基础设施精修任务；下一唯一任务立即回到 H08a 真正同源时空相通量积分，H11 未启动。

目标：在同一Hydro owner内一次切换实际 B/J、native查询、piece/atom及运动连续性证书，使材料与运动消费同一已验证正则场，关闭“只有暂存准备可用、实际来源仍旧”的部署缺口。

允许修改：production ≤5，核心为 `CchSpatialAdvector.h/.cpp`、`CchCharacteristic.h/.cpp`；若需第五个既有内部文件，在实施前设计中指明职责。test ≤2，优先原momentum/characteristic契约。无新生产变体、长期脚本、第二owner或运行时方法开关；超预算先按真实依赖拆解。

禁止修改：Pure FirstJet状态/库存/几何、material/dual的F/R/uR/KEP关系、最终native对流算子、原H10b输入及门槛。不得只平滑characteristic用速度而保留旧材料atom，也不得将旧路径留下作为隐式fallback。

数学约束：唯一共享F、内部反号/真实AMR子面和、来源identity/full clock一致、KEP/skew仍为材料四硬约束。真实native物理均值、cell积分散度、C1、原可辨识低阶与完整P2精度保持；cellDivergence是平均量，点态和dual平均配分差只诊断。无R=0、quarter库存或界面bubble。

实现要求：先列明当前trace/lift/piece/normal和运动证书消费者，形成部署前设计。实际field持有同一次regular apply的系数及原velocity/geometry/boundary/full-clock身份；所有值/J/完整normal多项式从这些系数取得。按真实support分片构建唯一primal/dual overlay，独立width与双端chart完整保留。native值须说明同一场的物理查询语义。退役无消费者的旧拟合/lift实现，prepare完成全部几何工作，更新不重复SVD。

运动连续性必须改为实际物理次数的完整值/J迹证书，不能继续用三点volume-coordinate二次拟合；明确换基/舍入工作量、轴正则与实际同源查询。轨迹与Jacobian采用[r17](H08a2r17_Nonlinear_Crossing_Jacobian_Accuracy_Report.md)已验证的同源三阶分片Picard内部积分/三次dense事件公式，保持原clock和事件语义；同步更新CchCharacteristic.h及现行API中的公式注释，不能保留旧“仅Ralston3端点”的描述。不伪造Jacobian或任选有跳跃的分支。旧“必须真实跳跃”的fixture应依据当前生产场数学能力改写为完整C1/负向证书检查，并保留能检测不连续/漏高项的独立反例，不能仅删除断言。

验收：T0标准构建、架构/目录和全部受影响契约；预算外源码/原input保持，旧路线确实退出生产。T1 r9的原64/144矩阵、208basis、原56适用精度门槛及r8/r10完整C1/face/Gauss/次数/缺秩检查通过，身份/完整clock/axis/wall/tiny/反号/AMR共享拒绝不减。T2当前实际运动重验原24时间族、分区比较、Jacobian差分、零/极小幅值和完整次数连续性负控，以及r16/r17两个原2560步point/J可信参考和point-only负控；不能用暂存原型代替实际source。然后重跑原内部affine反例、最终native约束和H10b-dissipation-v2原28窗口/84门槛，旧4项donor失败仍单列。完整研究依旧走既有benchmark入口，临时审计只放runs/checks。

失败规则：混用来源、遗漏高次/支持分片、界面反向修复、任何原精度/稳定性门槛失败或来源拒绝退化即FAIL，不能调容差或缩减矩阵；失败候选退役并保留证据。真正外部依赖才可BLOCKED。

交付：实际来源/运动证书改动、退役清单、完整物理与运动/最终算子研究证据、准备/更新成本及未完成范围。通过后立即返回H08a真正同源时空相通量积分与R收敛，必要的FirstJet消费接入并入该任务。H11未启动。

**基础设施收口（2026-10-05 用户授权）：** r11 是 common regular source/characteristic 的最后一道关口。除非出现可复现且直接阻断 H08a 的数学反例，不再拆 source/Jacobian/低秩精修任务，不新增硬约束、特殊 case，也不放宽原门槛。r11 通过后立即回到 H08a 的真正同源时空相通量积分；必要的 FirstJet 运动消费接入属于该积分任务，不再单列桥接先决任务。目标仅为 R 成为可收敛的小离散余项，保持 Pure FirstJet 独立、无 M_K/bubble、无 R=0 或 quarter commuting 强制。材料硬约束仍仅为唯一共享材料通量、内部面反号共享/AMR 子面守恒、来源一致和 KEP/skew。
