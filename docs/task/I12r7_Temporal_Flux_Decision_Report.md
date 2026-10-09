# I12r7 / 时间通量设计变更的实施前审查报告

**状态：PASS（设计审查范围）；I12：FAIL；Gate I：未通过。** 没有修改生产源码、既有测试或冻结方法。下一步是用户决定是否批准 [时间通量设计变更](I12_Phase_Flux_Time_Integration_Design_Change.md)，之后才可实施相应子任务。

按 [TODO](I12r7_Temporal_Flux_Decision_TODO.md) 执行。上一项 I12r6 属于 progress；本次以新的事件遗漏反例和现有 API 审查，修正实施条件和任务拆分，没有把诊断通过当作数值修复。

## 1. 当前状态复核

当前 1579 个源码/配置/算例文件与 I12r6 最终身份一致，19 个 runtime 文件一致。69 个冻结源码、12 个 integration reference、3 份方法 PDF 哈希一致；I12r6 归档 SHA256 仍为 `A3C2B79D439E0C013595BB1D467E93D0DD948C364A75E596A4600AB923D0FA17`。

当前程序实际执行两个独立入口：

| 命令 | 实际结果 | 含义 |
| --- | --- | --- |
| `quantum_interface_hermite_gals_benchmark.exe --crossing-probe` | 退出 1，两次拒绝 | q=1/4 库存虽有界但必须制造未标记气泡；q=3/4 超容 2%，仍不可接受 |
| `quantum_interface_hermite_gals_benchmark.exe --gross-flux-contract` | 退出 0 | I12r1 空间子集积分修复仍有效，不能继续将该旧缺陷列为未解决的时间问题 |

q=1/4 的 midpoint 跨面液量为 0.005，解析值 0.00375，上游库存 0.06125；q=3/4 分别为 0 和 0.00125，上游库存 0.06375；单元容量为 0.0625。I08/I10 的拒绝保留，未加 `WILL_FAIL` 或替换返回码。

本轮没有重建、没有重跑全量 CTest/完整 benchmark：源码/runtime 没有变化，已有完整结果按原字节哈希复核。I12r6 完整 benchmark 的 SHA256 为 `FC4483C3E9E94594CF6304995698E784FC7EB5CCFDB0EA557ED05D924289BB25`，其中仍有 10 项失败。这是已验证的历史完整运行，不冒充本轮新执行。

## 2. 新增独立解析证据

新增脚本使用 Python 标准库 `fractions.Fraction`，不调用生产几何或求积代码生成参考量。共检查 62 个 q、9 档 dt、正反方向、XY 与 4 个 RZ 轴向环带（含轴），即 5580 个**解析参考组合**。它们不是 5580 次生产求解。

若面测度为 A，流向宽度为 h，速度大小为 1，则

`M_old=A*(h-q*dt)`，`V_in=A*dt`，`V_exact=A*dt*(1-q)`。

精确更新恒为 `A*h`，而单点中点液量为 `A*dt` 或 0，因此缺陷比例是 `(dt/h)*min(q,1-q)`。RZ 使用 z 向平移及真实环带面积，共同 pi 精确消去。扫描 CFL 从 0.08 至 0.0003125，最小缺陷仍超过既有舍入预算约 `3.44e8` 倍。q=1/2 是当前不支持的中点整面零迹，明确排除于此矩阵并单独记录；没有用排除它来声称完整覆盖。

这一族初始界面位置随 dt 调整，证明固定 CFL 不能普遍保障全部跨面相位；不据此否认任意固定初值的时间细化或自适应算法。

另一个精确反例改变了提案的事件检测要求：在固定面上令 `phi(tau)=(tau-1/16)*(tau-1/8)`，它对应匀速平移薄层的面迹。液相仅在 `(1/16,1/8)` 存在。中点规则与两半区间中点规则都得 0，二者之差也是 0，而真实时间积分为 1/16。因此“粗细求积一致”不能单独认证没有漏掉相区间；必须另有同源事件覆盖证据。这不表示该薄层已能通过当前完整几何发布。

## 3. 可批准的变更及边界

更新后的提案请求将守恒相通量定义为同源时间区间积分除以 dt，保留 first-jet/M 分解、唯一共享实际面液量、reflux、严格库存、bubble 拓扑/幅值门槛和旧 A/B。

本次审查确认不能仅在 benchmark 插入多个中点样本：I08 的 velocity/result 绑定单个 midpoint identity，I09 的 source 类型固定。实施应按顺序拆为区间源身份、事件覆盖/单面积分、库存装配/I09 消费，最后接半步/全步 benchmark；每项另写 TODO，各自遵守文件预算。时间包围能力、正长度零迹和未解析事件是明确接口缺口，不以零通量或旧算法后备填补。

已修正提案中关于 vortex 空间 gross 缺陷的过时描述，并加入 I12r4 精确几何中点仍产生负库存的证据。时间积分不会自动解决所有 bubble 证书不足、空间重构、动态 regrid 或拓扑变化问题。若实施需进一步改变 characteristic/状态/时间几何方法，必须另列数学变更，不能借本提案默许。

## 4. 复现与交付

证据目录为 [I12r7_20261001](../../runs/baselines/I12r7_20261001/)，入口：

- [验证脚本](../../runs/baselines/I12r7_20261001/verify_temporal_decision.py)、[身份核验](../../runs/baselines/I12r7_20261001/identity_audit.json)、[综合状态](../../runs/baselines/I12r7_20261001/verification.json)。
- [当前跨面原始日志](../../runs/baselines/I12r7_20261001/crossing-probe.log)、[当前空间通量日志](../../runs/baselines/I12r7_20261001/gross-flux-contract.log)、[运行结果](../../runs/baselines/I12r7_20261001/probe_results.json)。
- [有理数矩阵](../../runs/baselines/I12r7_20261001/exact_crossing_matrix.csv)、[解析摘要及遗漏事件](../../runs/baselines/I12r7_20261001/exact_reference_summary.json)。

本机使用已安装 ParaView 的 Python 解释器运行标准库脚本，不需要安装包或启动 GUI：

```powershell
& 'C:\Users\impul\AppData\Local\Programs\ParaView-6.1.1\bin\pvpython.exe' --disable-registry runs/baselines/I12r7_20261001/verify_temporal_decision.py
```

首次读取旧 PowerShell 日志时误按 UTF-8 解码，因其 UTF-16LE BOM 失败；已记录 [初次错误](../../runs/baselines/I12r7_20261001/initial_audit_failure.json)，改为按 BOM 只读解码后完整脚本退出 0。旧日志原字节及哈希未变。这是证据读取修复，不是数值修复。

## 5. 待决定及剩余目标

请求批准上述时间通量设计变更。确认要求直接来自 [开发总纲 §2.1](../00_LAMP_Development_Charter.md) 的“不得由 Codex 在实现过程中自行改变方法”；当前仅有推进目标的授权，没有更换冻结时间离散的明确确认。本审查已把提案、证据、实施边界及验收准备完成，不再以旁支任务回避该决定。

本轮审查完成，尚未启动 I12a；用户决定未到之前不执行依赖此批准的生产修改。总目标继续 active；这不是连续三轮无法推进的 blocked 判定。I12 的全部失败、Gate I 后续要求，以及 Hydro、Electricity、低黏/无黏电喷雾验收均未完成。
