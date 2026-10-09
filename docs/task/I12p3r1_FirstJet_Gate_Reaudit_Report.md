# I12p3r1 / 当前纯路线 Gate I 复核报告

**Gate I：FAIL，两项实际诊断缺口；数值与此前三个主要交付缺口已通过。** 没有外部阻塞。按[TODO](I12p3r1_FirstJet_Gate_Reaudit_TODO.md)执行只读复核，production/test/case 修改均为 0。当前纯 Hermite–GALS 路线有效，原守恒路线不恢复；下一项先补字典范围错误上下文，再补同源曲率时序诊断。

本次直接检查六份现行规范、用户授权设计变更、当前 producer/App/reader/output/state/publication/mesh-transfer 接线及原始 CTest/实际 App 输出。p1、p2a、p2b、初次 p3、p3a–p3d 共 8 个 archive 外部 SHA 与每项内部 payload SHA 全部核验；1,625 个当前源码/测试/工具/输入/构建文件与最后归档重建出的身份相同，39 runtime 与实际测量二进制相同。69 硬冻结实现、11 未改接线参考、3 PDF 保持；允许扩展的 reader 参考不误算为硬冻结，旧 reader 主体证明仍在 p3d 已校验的证据中。历史 FAIL 报告及 archive 未改写。

| 现行要求 | 实际证据 | 裁决 |
| --- | --- | --- |
| I-R1 唯一高阶 scalar/geometry | 当前 state/publication、原 owner/revision 拒绝契约 | PASS |
| I-R2 first jet、直接中心、无 M/体积拟合 | producer 实际旧 scalar departure 查询；FirstJet state 拒绝库存；所有 full 输出 NA/drift | PASS |
| I-R4 同一实际 B/J | 原 trace/vertex/edge 完整 CTest 与当前模块调用 | PASS |
| I-R5 master/hanging、exact/project AMR | 16 exact/8 projection/45 guards、12 粗化三档组、36 完整 regrid | PASS |
| I12 原静止/平移/旋转/涡旋与收敛/crossing | 24 实际组/2,072 步、原空间/曲率/时间阶、48 crossing | PASS（规定范围） |
| 正式模块 producer 与唯一 App 事务 | p3a 16 场景/80 guards；App 在物理步及计划转移均准备成功后接受 | PASS |
| 显式支持字典、默认与配置追踪 | 当前 strict reader，实际 preflight、effective config、input hash/trace | PASS |
| 字典规范 §9 完整范围错误上下文 | 实际 dt/resolution/regrid 范围错误缺 file:line/限定字段/raw value | **FAIL** |
| 同源只读 VTK/accepted/rejected/series/meta | p3c ParaView 实读；p3d 完整 VTK 字节相同；18 频率组、真实 IO 失败/原子重试 | PASS |
| 输出规范 §9 曲率 max/RMS/unavailable count | 当前 header/source 仅 kappa_max，无 RMS 和不可用采样计数 | **FAIL** |
| 实际路径阶段计时、资源与 native 状态 | p3d 18×3 固定/4 动态规模、24 full 结构及真实计时包含关系 | PASS（已测范围） |

原 36 CTest 和新增 App 契约共 37 项通过，属于当前实现；本审计不重复运行这些确定性计算。直接从 CTest detail 抽取数值，再与保存的完整输出逐字比较。三个空间研究的原 gate 均为末档 gradient 阶 ≥2、curvature 阶 ≥1.5，旋转时间位置阶 ≥1.8，未放宽。

| 三档研究（n=8/16/32） | XY uniform | XY AMR | RZ uniform | RZ AMR |
| --- | --- | --- | --- | --- |
| 静止非多项式 gradient / curvature 阶 | 3.958 / 3.304 | 4.095 / 3.413 | 4.015 / 2.970 | 4.082 / 2.970 |
| 非多项式平移 gradient / curvature 阶 | 3.520 / 2.817 | 4.283 / 3.191 | 3.503 / 2.442 | 3.503 / 2.460 |
| 实际粗化后继续平移 gradient / curvature 阶 | 3.504 / 3.033 | 3.721 / 3.034 | 3.594 / 2.888 | 3.503 / 2.442 |

旋转 t=.3，dt=.01/.005/.0025 的末档位置阶为 uniform 2.00417、AMR 2.00468。48 crossing 的原 volume ≤5e−12、jet ≤1e−11、curvature ≤1e−10 逐项通过。Q2 静止/平移精度与静止漂移原要求保留；三档研究与传输测试中的非多项式数据避免只验证可精确再现的多项式。

四个完整 vortex 都运行 400 步，dt=.001，t=.4。下面保留最大界面/梯度/曲率误差和整个运行的体积漂移峰值，不能把最终漂移替代峰值。

| full 配置 | 位置 | 梯度 | 曲率 | peak ΔV/V0（绝对值） |
| --- | ---: | ---: | ---: | ---: |
| fixed uniform | 5.5143e−6 | 2.6255e−4 | 9.0853e−3 | 2.5417e−5 |
| fixed AMR | 2.2977e−6 | 1.6293e−4 | 1.6836e−2 | 9.3955e−6 |
| dynamic，初始 uniform | 3.6854e−6 | 4.6136e−4 | 3.2678e−2 | 2.0661e−5 |
| dynamic，初始 AMR | 6.1109e−6 | 2.4299e−4 | 8.3749e−3 | 2.1376e−5 |

dynamic 每组交替三次，最终网格与初始不同，不能凭表声称 AMR 总是更准。完整误差、全局解析体积误差与逐步漂移见机器证据。FirstJet 没有独立库存或相通量 ledger；取消相体积守恒不等于所有体积误差为零。

两项缺口均有当前源码和已归档真实运行证明：

1. `resolution=3`、`dt=0`、非递增/最终步 regrid 等输入经 typed config 的 validate 拒绝，进程确实失败且不生成解，但错误只有 `[Interface FirstJet App] ...`。相同 reader 中类型、未知字段和重复字段通过 Node.fail 已带完整上下文。修复应留在 reader，把范围检查映射回已有 Node；typed config 保留独立防御校验。任务 I12p3e 不改变合法输入、默认或数值路径。
2. `conservation.csv` 当前记录 cut/graph 数和 `kappa_max`；输出规范 §9 还要求曲率 RMS 与 unavailable sample count。当前 writer 和 header 都缺这两项，不能以存在 VTP 点值代替全部接受步诊断。任务 I12p3f 应只读当前 publication 的同源 quadrature，定义 RMS 权重和空界面语义，明确记录不可用数量，不能另拟合曲率或改变 acceptance。

p3d 在已测三档未出现稳定二次或更高增长，也没有返回同网格后的 primary/face/center 累积；publication 占 n=32 App step 约 66.4%–91.5%。这是有限网格/窗口、Windows 采样内存的观察，非渐近复杂度或异机吞吐保证。方法只验证六类规定速度独立 Interface，尚未接 Hydro，也没有任意入口、真实 jet/tip/neck breakup 或低黏稳定证据。后续 H07/H-R3 的质量—动量方案仍须单独设计。

证据：[I12p3r1_20261002](../../runs/baselines/I12p3r1_20261002/)。gate_matrix/verification 为当前裁决，archive_identity/source_review 保存身份与检查范围，convergence/full_metrics 来自原始数据独立抽取；input_diagnostic_gap/output_diagnostic_gap 固定缺口证据。package_evidence.ps1 只归档本审计与文档，引用已校验的前序运行包；最终完整性标识见 archive_verification.json。

下一顺序：I12p3e → I12p3f → 只读再复核 Gate I。两项补齐并验证后才转 H01。原 I12p3 的主要缺口已关闭，其首次 FAIL 仍作为历史审计保留；当前 Gate I FAIL 是上述两项必要诊断，而不是原数值方法未通过。
