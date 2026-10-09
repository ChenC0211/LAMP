# I12p2b / 纯 FirstJet 动态 AMR 结果

**PASS（本任务）；完整 Gate I 仍未通过。** [任务](I12p2b_FirstJet_AMR_Transfer_TODO.md) / [设计](I12_FirstJet_AMR_Transfer_Design.md) / [API](../api/I12p2b_FirstJet_AMR_Transfer_API.md)。本任务修改 4 个 production 文件、2 个测试文件，实现零时间 first-jet 网格传递，未引入 M、体积拟合或相库存 reflux。

`ExactRefinement` 直接发布 Core exact Q4 restriction 的 snapshot，primary/closure 仅为其派生记录；`FirstJetProjection` 从旧 scalar 采样并用原 I06/纯 publication 重构。两个操作显式选择，失败不切换。共同的测量、拓扑、cut、法向和曲率检查提取为 private publication 方法，外部不能传任意 snapshot 绕过 producer。

16 组 exact 传递覆盖 XY/RZ、uniform/AMR、同网格/细化和 Q2/非多项式。最大 scaled value/gradient/Hessian 差为 7.561e-13，trace 差 2.665e-15，同点物理曲率差 4.113e-13，父叶体积加和差/容量 1.791e-15，均满足任务原阈值。8 组细化后投影返回源网格，其中全局 Q2 保持原精度；非多项式的误差如实记录，不声称普遍精确。

45 项拒绝检查包括：exact 模式禁止粗化、非法模式、错误源网格、非法配置、moved-from、非嵌套叶、变域、旧守恒模式、局部支持，以及粗化丢失已解析小液滴。源 publication 保持可用。分量及无方向端点类型检查只能识别相应拓扑变化，不能代替形状/曲率误差界。

动态 full 沿用原 12 组 n=8、物理参数与结束时间，完成 1,036 个实际接受步、36 次网格切换。每组在 floor(N/4)、floor(N/2)、floor(3N/4) 完成步后交替 uniform/AMR；每次切换保留完整时间层，并继续原 characteristic first-jet 输运和 I06，没有解析重置。smoke 为每组 4 步，也实际发生 36 次切换。

| 算例/初始网格 | 结束位置误差 | 结束梯度误差 | 结束曲率误差 | 峰值相对液体体积漂移 |
| --- | ---: | ---: | ---: | ---: |
| rotation / uniform | 6.19710e-6 | 4.84508e-5 | 3.72592e-4 | 1.22514e-6 |
| rotation / AMR | 6.17614e-6 | 4.84508e-5 | 3.69529e-4 | 1.22367e-6 |
| vortex / uniform | 3.68539e-6 | 4.61361e-4 | 3.26784e-2 | 2.06606e-5 |
| vortex / AMR | 6.11089e-6 | 2.42993e-4 | 8.37487e-3 | 2.13759e-5 |

标签表示初始网格，三次切换后最终网格相反。表中并非固定网格精度比较，不能据此声称 AMR 总能改善曲率。vortex 没有逐 cell 解析体积参考，原始 `cellVolume=nan` 在结构化证据中记为 null；总体几何体积和初值漂移仍实际测量。库存残差和 flux ledger 始终为 NA。

非多项式粗化收敛覆盖 XY/RZ × uniform 细网格/局部 AMR 来源。每档先细网格真实平移 2 步，投影到 n=8/16/32，再粗网格推进 2 步，dt=.01、终止 t=.04。

| 度量/来源 | 末档梯度阶 | 末档曲率阶 |
| --- | ---: | ---: |
| XY / uniform | 3.50361 | 3.03292 |
| XY / AMR | 3.72056 | 3.03409 |
| RZ / uniform | 3.59352 | 2.88817 |
| RZ / AMR | 3.50307 | 2.44213 |

全部三档误差递减，满足原 >=2 / >=1.5 门槛。这是实际粗化后继续推进的验证，不借用固定网格结果冒充动态传递。

首次接口测试的局部支持夹具缺少原 Core 要求的外边界接触点；补充接触点后，x=.5 平面又触发已有 coincident-zero-trace 拒绝。最终改为具有实际外边界接触点的 x=.53 平面，进入并通过预期的“禁止部分支持传递”检查。两次失败输出保留，未修改 Core 或放宽检查。另一次错误可执行路径未启动测试，单独记为 NOT_RUN。回归脚本误选全部 119 个注册测试，人工终止后保留 INTERRUPTED 记录，并改为预定 35 项范围；不把被中断运行计为通过或数值失败。

39 个目标编译通过（45.34 s），无新 warning；35 项 CTest 全通过（365.27 s），另直接执行 state/volume/publication/exact-refinement 四个受影响契约通过。原 27 项数值测试及 p2a 四项研究完整输出逐字保持；架构检查仅增加 1 TU、2 文件、6 个仓内 include，仍为 5 条显式 legacy edge。历史 crossing 拒绝及其完整输出保持。最终动态 full 为 238.31 s，1,036 步、36 regrid 的完整输出与首次 full 一致；smoke、粗化研究也与最终 CTest 输出一致。

独立审计重新计算了全部粗化阶数、每步体积漂移及每组峰值，核验每次网格切换的步号/时间/显式方法序列；没有把 wall time 用作数值分支。81 个冻结文件、3 PDF 与 p2a archive 哈希保持。

证据：[验收汇总](../../runs/baselines/I12p2b_20261002/verification.json)、[完整传递/粗化/动态数据](../../runs/baselines/I12p2b_20261002/transfer_verification.json)、[范围审计](../../runs/baselines/I12p2b_20261002/source_audit.json)、[旧输出比较](../../runs/baselines/I12p2b_20261002/regression_comparison.json)。

当前完整 Gate I 仍未通过，后续 I12p3 必须汇总纯路线准确性、体积漂移、效率及 production 入口/诊断缺口；不能据本项宣告低黏流体耦合或项目目标完成。
