# I12p3r2 / 当前纯 Hermite–GALS Gate I 最终报告

**Gate I：PASS，当前显式 `hermiteGals` / FirstJet 独立规定速度路线。** I12p3r1 的两个必需诊断缺口已由 p3e/p3f 关闭。按[TODO](I12p3r2_FirstJet_Gate_Final_TODO.md)完成零源码/测试/输入修改复核，下一路线图任务为 H01。宽泛电喷雾目标仍 active；此 gate 不宣告 Hydro、低黏/无黏耦合或 jet/neck/breakup 完成。

| 当前必要项 | 证据与裁决 |
| --- | --- |
| 唯一 Q4 scalar/geometry 与 first jet | **PASS**：当前 state/publication/producer，原 owner/revision/拒绝契约；无第二份几何 |
| 无 M/无体积拟合/真实中心预测 | **PASS**：FirstJet phaseVolumes 为空；中心从同一旧 scalar 的实际 departure 查询，几何体积只作 drift/error 诊断 |
| 实际 B/J、I06/master/hanging | **PASS**：原 trace/vertex/edge/AMR closure 完整数值契约和当前模块接线 |
| exact/project 动态 AMR | **PASS**：16 exact/8 projection/45 guards、12 三档粗化组、36 完整 regrid，原 source/time 保持 |
| 静止/平移/旋转/涡旋/收敛/crossing | **PASS（原范围）**：24 完整组/2,072 步、三档空间/曲率/时间、48 crossing；误差门槛保持 |
| production 模块与唯一 App 事务 | **PASS**：模块 producer 16 场景/80 guards；App 物理步及计划转移均成功后接受 |
| 支持字典/默认/追踪 | **PASS**：实际 strict reader、effective/configuration trace/input hashes |
| 完整输入错误上下文 | **PASS**：p3e；当前15精确范围拒绝、6合法边界，字段/原始值/行号/method/range 完整 |
| 同源只读输出与失败语义 | **PASS**：当前18频率组/48帧/96VTK、数值拒绝、真实IO失败、已有输出保护和原子替换重试；前序ParaView实读 |
| 曲率 max/RMS/unavailable count | **PASS**：p3f；当前90初始/接受态统计，XY/RZ解析 RMS、物理求积权重和频率独立，原数值/VTK保持 |
| 阶段成本/结构/资源 | **PASS（已测范围）**：p3d 18配置×3重复和4动态组；native状态界、完整400步vortex结构及原计时包含关系 |
| legacy/数学/证据身份 | **PASS**：11前序archive、69硬冻结/11未改接线参考/3PDF、当前1,625文件/39runtime身份与诊断修复证明 |

当前数学 kernel 与 p3d 相同。p3d 后的完整源码差异仅为 p3e 的 typed error 包装/reader 定位、p3f 的只读曲率统计和两个既有 App 测试。p3e 去掉异常包装后，Config 判据/顺序、App 其余主体、状态布局、reader 合法输入和旧分支逐字相同；p3f 去掉统计块及追加列后 writer 逐字相同。当前 App/reader 与 p3e 验证版本、writer 与 p3f 验证版本精确匹配。

验收证据的复用范围明确如下：

- 35 个数值测试：数学源码及38其它运行文件与 p3d 精确相同，复用原完整数值输出。
- 当前 App/architecture：p3f 实际执行2个CTest，二进制 SHA 与当前一致；App包含原12 smoke/48步/36regrid、24非法输入、15新范围/6合法边界、CFL及legacy，再执行18输出频率/真实失败/重试契约。
- 24长窗口/2,072步/36regrid、空间/曲率/时间和48crossing、58正式规模运行：由上述未改数学实现和合法 App/诊断同一性证明继承 p3d/p2/r1。没有声称最终App二进制重新执行这些长窗口或规模测量。

原 spatial gate 为末档 gradient阶≥2、curvature阶≥1.5，旋转时间位置阶≥1.8；实际末档范围分别为 gradient 3.503–4.283、curvature 2.442–3.413，旋转约2.004。三档完整数据与vortex误差/体积漂移表保留在[I12p3r1报告](I12p3r1_FirstJet_Gate_Reaudit_Report.md)，没有将体积漂移隐去或当作守恒残差。原四个完整400步vortex全部完成，最大相对V0漂移约9.40e−6至2.54e−5；dynamic最终网格不同，不能据此保证AMR总是更准确。

效率结论是有限三档/窗口的结构与增长检查。p3d 主成本为唯一 publication；p3f 只在IO中遍历既有只读求积节点，没有新增持久几何或全局物理自由度。当前IO包含这项新统计开销；前序绝对wall/采样内存不包装成最终App的重新测量或跨硬件吞吐保证。

共11份前序 archive 的外部SHA和**11,438项内部payload SHA**逐项通过；其中初次p3与p3r1的FAIL历史保持。当前1,625个源码/测试/工具/输入/构建文件从p3f归档重建并验证，39runtime均匹配当前文件，69硬冻结实现、11未改接线参考及3PDF保持。允许扩展的reader参考单列，不把其误算为硬冻结。当前 QUANTUM.exe SHA-256：`77D2369F4E8D8D5888570EA1A56C727F4D3D5DEFAC7C082A5ABB73D192C4CB1F`。

此次PASS仅适用于当前规定速度独立Interface候选。部分 native chain 仍需显式 boundary contacts，正长度coincident zero等未支持几何保持明确拒绝；没有任意入口几何、真实tip/jet/neck破碎或Hydro/EHD验收证据。H01先建立native face/cell布局；H07/H-R3接入前须另设计非守恒Interface下的密度更新、唯一质量—动量通量和误差/稳定性验收。后续压力相容、黏性PSD、mu=0、电荷/Gauss和统一受力要求保持。

证据：[I12p3r2_20261002](../../runs/baselines/I12p3r2_20261002/)。gate_matrix为当前完整裁决；verification记录本次实际检查和继承范围，archive_identity为全部前序内部完整性，source_before/runtime_hashes固定当前实现。package_evidence.ps1只归档本审计与当前状态文档，不重复包装旧数值运行；最终标识见archive_verification.json。
