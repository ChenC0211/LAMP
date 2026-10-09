# I12p3a / 可复用 FirstJet 时间步结果

**PASS（本任务）；Gate I 尚未通过。** [任务](I12p3a_FirstJet_Step_Producer_TODO.md) / [API](../api/I12p3a_FirstJet_Step_API.md)。新增 2 个 production 文件，修改 benchmark/测试注册 2 个文件；既有 production 文件全部保持。

新 `prepareHermiteGalsFirstJetStep` 接收 accepted FirstJet publication、实际 trace domain、前后速度 view 和真实外边界 first jet，组合原 vertex/edge transport、实际中心回溯、I06 与纯 publication，返回新 publication 和来源/DOF 计数。没有新增离散公式、独立 M、体积校正或 accepted commit。

原测试中的 `finishFirstJet` 已删除；所有纯固定网格、收敛、48 crossing、粗化和动态网格时间步改调新 API。测试只提供规定速度、真实外边界数据、时钟和配置。历史 `predict/finish` 留作 ConservativeVolume 比较，明确禁止被纯路径调用，没有 fallback。

新契约覆盖 XY/RZ × uniform/AMR × Q2/非多项式 × midpoint/end 共 16 组，检查时间/来源、新 identity、单一 scalar、实际内外边界 DOF 计数和原 Q2 精度。80 项拒绝覆盖旧守恒源、错误来源/时间/domain/网格、无效 publication 配置、缺失/非有限/抛错 boundary provider、moved-from、原 CFL 超限，以及真实内部 edge/center 边界事件。边界 provider 不能填补内部回溯。

测试拒绝条件复核保留了真实过程：首版宽泛的 `boundary` 字符匹配也能匹配 `LevelSetBoundaryTrace` 的 CFL 异常，因此该首版不作为最终“内部命中”证据。将检查细化到实际事件消息后，原大 dt 夹具被 CFL 提前拒绝；后续 center 剪切夹具也因速度系数算错被 CFL 拒绝。最终分别设置独立的超 CFL、合法半 cell CFL 上限的内部 edge 命中、零水平边速度的 center 命中，16 组/80 拒绝全部通过。两次修正期间的失败日志保留，未放宽 tracer、CFL 或数值门槛。

39 个目标编译通过（41.04 s），无新 warning；36 项 CTest 全通过（369.17 s）。原 34 项数值测试完整输出逐字保持；架构检查仅新增 1 TU、2 文件、6 个仓内 include，5 条 legacy edge 保持。历史 crossing 的原拒绝和完整输出也保持。

固定网格 full（223.87 s）与 p1 档案完整输出逐字相同，动态网格 full（232.40 s）与 p2b 档案完整输出逐字相同。共 24 个完整算例、2,072 个接受步及动态组的 36 次网格切换，原物理参数、网格、dt、终时和所有误差/体积漂移均保持。此次 wall time 不用于声称性能优化；规模研究仍属 I12p3d。

范围审计确认既有 production 全部不变，仅新增 producer 两个文件；81 freeze、3 PDF、p2b/p3 archives 保持。证据：[验收汇总](../../runs/baselines/I12p3a_20261002/verification.json)、[算子与完整驱动核验](../../runs/baselines/I12p3a_20261002/step_verification.json)、[范围审计](../../runs/baselines/I12p3a_20261002/source_audit.json)、[旧输出比较](../../runs/baselines/I12p3a_20261002/regression_comparison.json)。

此项仅完成模块时间步归属；正式独立 App/字典入口为下一项 I12p3b，正式输出和效率研究随后完成。不把当前模块接口声明为已接通 Hydro 或已完成完整 Gate I。
