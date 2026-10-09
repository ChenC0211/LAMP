# I12p3r1 / 当前纯路线 Gate I 复核

目标：只读核实 I12p3 首次审计的三个交付缺口是否已由 p3a–p3d 关闭，并按现行六份规范和用户授权的非守恒设计变更裁决完整 Gate I。

允许修改：本 TODO/报告、当前规范入口中的状态及必要的历史/当前措辞、runs/baselines/I12p3r1_20261002 下审计脚本/证据。production/test/case 修改预算均为 0。

禁止修改：所有源码、测试、输入、容限、旧报告/归档、三份 PDF；不实施 Hydro，不恢复相库存/通量认证，不重跑已经有当前二进制身份和完整输出证明的确定性计算。

数学约束：FirstJet 无 M/体积拟合，中心来自同一实际旧 scalar 回溯，B/J 同源，I06 mixed 和 master/hanging 保持，唯一 authoritative Q4 geometry；体积漂移必须可见。AMR exact/projection 语义、原位置/梯度/曲率/时间门槛不变。后续质量—动量耦合须单独设计。

实现要求：建立规范要求→当前 API/源码→实际原始证据→裁决的矩阵，区分 Interface 规定速度验收与尚未开始的 Hydro/EHD。直接核对 producer/App/reader/output/AMR 接线，核对旧数值输出、完整窗口、收敛、体积漂移、结构/资源和输出证据。逐条校验 p1/p2a/p2b/p3/p3a/p3b/p3c/p3d archive 及内部清单哈希，当前 production/test/case 与最后相关归档相同；69 硬冻结、12 接线参考和 PDF 分开审计。

验收：T0 当前源码/39 runtime 和冻结身份匹配，37 CTest 结果属于当前实现；T1 原代数/拒绝/ownership 契约和三个缺口逐项复核；T2/T3 原 24 App 完整组、三档空间/曲率/时间及动态 AMR、48 crossing、输出/规模研究均完整且无放宽门槛。基于已验证原始输出审计，不以绿色汇总文件代替证据。纯 Interface 的 tip/jet/low-viscosity 耦合验收留给 H/IH/E/C，不能宣称已经完成。

失败规则：有当前必需项未实现或没有有效证据即 FAIL 并列出下一单任务；只有实际缺失外部接口才记 BLOCKED。不得因为投入多或子任务全部 PASS 就自动通过。原 I12p3 FAIL 及历史守恒 FAIL 不改写。

交付：零源码修改审计、完整 gate 矩阵、准确性/漂移/性能及局限说明、当前状态更新、可重现脚本和逐项哈希归档。若 Gate I PASS，当前独立 pure Hermite–GALS 才成为候选，下一路线图任务 H01；宽泛电喷雾目标仍 active。

执行结果：**FAIL（两项实际诊断缺口）**，审计完整性 PASS，见[报告](I12p3r1_FirstJet_Gate_Reaudit_Report.md)。当前数学/数值、独立 App/正式输出主接线与效率证据通过；字典范围错误上下文及曲率 RMS/unavailable count 仍须补齐。零源码/测试修改，无外部阻塞。
