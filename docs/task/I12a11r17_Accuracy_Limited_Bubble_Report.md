# I12a11r17 / 几何精度预算下的局部 bubble 闭合报告

**PASS（局部算子）；完整驱动未切换，I12 / Gate I 仍 FAIL。** 按 [TODO](I12a11r17_Accuracy_Limited_Bubble_TODO.md) 修改 2 个 production 文件和 2 个 test-related 文件。[API](../api/I12a11r17_Accuracy_Limited_Bubble_API.md) 与[已授权设计变更](I12_Geometry_Inventory_Accuracy_Design_Change.md)说明调用及边界。

## 数值变化

原 I10 内核增加显式 `GeometryAccuracy` 契约：cut cell 以用户请求的 `relativeVolumeTolerance`（原默认 `1e-10`）验收几何体积差，不再暗中缩到库存 `16/64*epsilon`。误差已达标时保留预测几何；需要校正时，将 Core 的拓扑/位移证书与逐 cell 的 Hessian 扰动预算相交，复用原 safeguarded root。区间端点向内处理舍入，单向扩展仍受原始精度预算限制，最终实际存储值再次检查。

新增预算由调用者按几何 cell 顺序提供，表示物理 scalar Hessian 任一分量的最大新增变化。全单元界为 `abs(delta_c)*max(32/hx^2,32/hy^2,16/(hx*hy))`。这控制新增导数扰动，不声称已经控制总曲率误差或确定完整推进的最优预算。

独立 M、first jets、mixed 和候选 identity 不变；pure 仍使用原 roundoff 预算，拓扑拒绝及原迭代/求积工作量保持。旧 `InventoryRoundoff` 契约保留数值行为及过渡期默认，用于现有驱动和显式 A/B；没有隐式 fallback。

## 算子结果

14 组 corner cut：正反相 × 7 档 cut 尺度，原 `128*epsilon*capacity` 的库存差为 `2.8421709430404007e-14` 相对容量。新契约全部保留该差并如实报告，center 完全不变，Hessian 改动界为零。最弱一例还通过现有 I11 默认容限发布，保留原 M 和 identity；这是单次发布检查，不是接受时间步。

16 组实际校正：XY/RZ × uniform/AMR × 正反相 × center 正负变化。用同源目标制造约 `0.02` 的 Hessian 分量扰动，再以明确的 `0.04` 预算求解；该预算只属于此算子夹具，不作为生产校准。独立显式 bubble 多项式在每组 9×9 点检查 value/gradient/Hessian 与边界 C1。

| 量 | 最大实测值 | 门槛 |
| --- | ---: | ---: |
| active cut 最终几何体积差 / 容量 | `1.6574223840070993e-14` | 请求的 `1e-10` |
| 全单元 Hessian 分量变化界 | `0.02000000000479929` | 输入预算 `0.04` |
| 与独立 bubble 公式的 jet 差 | `7.2442052356791464e-15` | `5e-11` |
| 边界 value/gradient 变化 | `0` | `1e-12` |

本次 active 根残差偶然远小于请求值，不意味着新契约仍要求机器精度。9 类拒绝通过：缺少预算、错误 cell 顺序、负/NaN 预算、旧契约携带新预算、非法契约值、零/过小 active 预算及 pure 超原预算。重复结果一致，拒绝不改变原 M/identity。

## 构建与回归

31 个 executable 全部构建成功，命令耗时 `45.3573196 s`。12 项 CTest 全部通过，命令耗时 `249.602423 s`；新增测试 `16.74 s`，首次直接执行 `17.0828498 s`。11 项原输出（含 architecture）与 r15 逐字一致；legacy crossing 原 exit=1 及输出保持。

最终审阅补充了区间端点的乘积检查，确保除法舍入后的端点仍满足物理 Hessian 预算。随后重新构建全部目标，耗时 `35.5917959 s`，并直接复验受影响的新契约；标准输出与此前 CTest 逐字一致。原契约数值分支未变，保留两轮执行的实际日志，不把它们合并为一次原始 CTest。

源码审计确认 2 production + 2 test-related，构建无 warning；81 个冻结/参考文件、3 PDF 以及 r15/r16 archive 不变。新文件未改变原构建模块依赖，architecture 输出保持。

## 未完成边界

本项实现了用户要求的局部几何匹配放宽及导数保护，但尚未把新契约用于完整时间推进。旧几何通量排空单元时可能把非机器精度的 M−V 差带成负库存；必须先完成同源共享通量的库存协调。不能用放大 roundoff 修复常数掩盖，也不能将新算子 PASS 当作非多项式曲率收敛通过。后续仍需原 h/dt、rotation/vortex、AMR crossing 与曲率矩阵，未通过前不进入 Hydro。

证据：[目录](../../runs/baselines/I12a11r17_20261001/)、[验证摘要](../../runs/baselines/I12a11r17_20261001/verification.json)、[30 组数据](../../runs/baselines/I12a11r17_20261001/accuracy_comparison.json)、[原输出比较](../../runs/baselines/I12a11r17_20261001/regression_comparison.json)、[源码审计](../../runs/baselines/I12a11r17_20261001/source_audit.json)。同目录保存 runtime 哈希、执行脚本、交付清单、`evidence.zip` 和 archive 验证结果。
