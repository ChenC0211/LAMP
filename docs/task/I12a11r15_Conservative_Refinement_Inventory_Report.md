# I12a11r15 / 精确细化的守恒库存转移报告

**状态：PASS（受限库存转移算子）；I12 / Gate I 仍为 FAIL。** 本项不是 accepted regrid，也没有完成动态加密推进。按 [TODO](I12a11r15_Conservative_Refinement_Inventory_TODO.md) 新增 2 个 production 文件及 2 个 test-related 文件；接口见 [API](../api/I12a11r15_Conservative_Refinement_Inventory_API.md)。

## 实现与适用范围

`refineHermiteGalsInventory` 拥有原不可变状态、旧新网格身份、r14 精确细化 Snapshot、目标库存与原始几何测量。未变叶复制原库存；细化叶在同一父叶内，以实际子几何液量为参考，把原父库存差额按可用增减量分配。目标始终为原始 `M_P`，没有几何重置、跨父叶借量或事后 clipping。

输入及子库存与同源几何之差受本任务原 `64*epsilon` 容量预算限制。逐父/全局守恒和实际存储后的严格容量界限分别检查；无可行分配时明确拒绝。结果仅为准备对象，不能冒充 I11 发布或接受步。pure/cut 测量沿用原相位证书及局部 Green 路线。

## 验收结果

60 个场景覆盖 XY/RZ、正反相、未变叶、均匀细化、uniform→AMR、AMR→更深层和跨两层细化；另含正负库存扰动、纯相、重复转移、copy/move 与外部来源/网格销毁。10 项拒绝覆盖错误身份、域变化、coarsening、部分覆盖、非法库存、超原几何预算、父叶无可行分配及 moved-from。

| 量 | 实测最大值 | 本任务门槛 |
| --- | ---: | ---: |
| 父子库存残差 / 父容量 | `1.884773144570973e-16` | `4*epsilon` |
| 原库存与原几何偏差 / 容量 | `9.423865722854864e-16` | `64*epsilon` |
| 子库存与子几何偏差 / 容量 | `3.788394020587656e-15` | `64*epsilon` |
| 原父库存与子几何和的差 / 父容量 | `2.085030291181639e-15` | 诊断 |

累计检查 720 个细化父叶、384 个未变父叶及 4,032 个目标叶。输出库存与再次从同一结果 Snapshot 积分得到的测量一致；扰动库存没有被静默替换为几何体积。

30 个 executable 构建成功，完整构建命令耗时 `31.996725 s`。11 项相邻 CTest 全部通过，命令耗时 `431.5216561 s`；9 项原输出逐字不变。architecture 变化仅为新增文件及 include：12 targets / 381 translation units / 903 files / 3102 includes。legacy crossing 仍以原输出退出 1；未以本项代替完整 I12 重跑。

最后修正了父级不可行分配异常中的 child 上下文：父级判断显示 `child=-1`，逐子判断保留实际 child。最终重新构建全部目标并直接复验本契约，60 场景 / 10 拒绝通过，耗时 `26.798178 s`；完整标准输出与此前 CTest 逐字一致。测试不要求某个固定父叶先失败，而是验证有效父叶及正确的父级上下文。

全部构建无 warning；81 个冻结/参考文件、3 份 PDF 以及 r14 archive 哈希保持不变。保存了初次 CMake glob 仅再生成、测试 lambda 编译错误、partial fixture 错误调用全域测量、错误硬编码先失败父叶等中间记录；均未通过放宽数值预算修复。

## 用户新决定与后续顺序

用户在本项收尾时明确授权：如果逐 cell 守恒和曲率精度难以兼顾，放弃每个 cut cell 的 Hermite 几何必须机器精度匹配守恒体积。该决定改变后续优先级。本项按原受限契约归档，其 `64*epsilon` 几何适用域不再被推广为后续生产方案的硬要求。

下一步优先形成已授权的几何—库存误差分离设计变更，保留独立 `M_K` 的共享通量守恒、严格界限与唯一 Hermite 几何；单独控制几何体积差和 bubble 对导数/曲率的扰动。仅调大闭合容差会重现已有原始负库存反例，必须与通量有界性协调。原计划中的 accepted regrid / 局部细化重试不再是恢复机器精度闭合的先决路线。I12 原精度与收敛门槛不变，未通过前不进入 Hydro。

证据：[目录](../../runs/baselines/I12a11r15_20261001/)、[验证摘要](../../runs/baselines/I12a11r15_20261001/verification.json)、[60 场景](../../runs/baselines/I12a11r15_20261001/inventory_comparison.json)、[回归比较](../../runs/baselines/I12a11r15_20261001/regression_comparison.json)、[源码审计](../../runs/baselines/I12a11r15_20261001/source_audit.json)。同目录保存 `evidence.zip`、交付清单、失败日志及 archive 核验结果。
