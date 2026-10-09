# I06b / AMR mixed closure 交付报告

**结论：PASS；I06a + I06b 完成父任务 I06。** 执行日期：2026-10-01。下一项 I07：守恒 cell phase volume。Gate I、低黏/无黏与完整电喷雾目标仍未完成。

依据：[实施前 TODO](I06b_AMR_Mixed_Closure_TODO.md)、[API 与公式](../api/I06b_Hermite_Gals_AMR_Closure_API.md)。开工核验 I06a 的 10 个交付文件及 archive 哈希，均一致；保存 1,559 文件的开工清单。上一目标轮 I06a 实际完成接口和验证，归为 progress。

## 交付范围

| 文件 | 变化 |
| --- | --- |
| [HermiteGalsAmrClosure.h](../../include/Modules/Interface/HermiteGalsAmrClosure.h) | 同一 primary revision 的 mixed + exact skeleton 只读结果与诊断 |
| [HermiteGalsAmrClosure.cpp](../../src/Modules/Interface/HermiteGalsAmrClosure.cpp) | master 覆盖、连通五点、约束依赖消元、RZ、直接公式代回和失败隔离 |
| [Core skeleton header](../../include/Core/Geometry/ActiveQuarticHermiteSkeleton.h) / [cpp](../../src/Core/Geometry/ActiveQuarticHermiteSkeleton.cpp) | nodeSensitivity 公有入口，保留严格浮点实现边界 |
| [Core constraint engine](../../include/Core/Geometry/Detail/LevelSetHermiteQuarticLeafData.h) | 复用原 cell pullback 的同一 DAG 累加段 |
| [AMR contract](../../tests/interface_hermite_gals_amr_contract.cpp) / [注册](../../tests/CMakeLists.txt) | Q4、MMS、C1、独立参考、真实输运、support/identity/拒绝路径 |

预算：5 production + 2 test；4 个既有文件改变、3 个新文件。I01–I05、mesh、Apps/reader、Hydro/Electricity、case 和 S08 专有文件不变。Core 原 reverse accumulation 段仅去空白后逐字一致；I06a 已有约束输出保持不变。

## 数学验收

所有数值容限和空间阶门槛保持 TODO 预设值。证据目录：[runs/baselines/I06b_20261001](../../runs/baselines/I06b_20261001/)，最终原始日志 `amr_final.log`，各 CSV 由该日志提取。

| 检查 | 实际覆盖 | 结果 |
| --- | --- | --- |
| Q4 | shifted/anisotropic XY/RZ；uniform、单格、相邻区域、三级 2:1；单格上 26 个 XY / 16 个 RZ 场，其余混合 Q4 | 最大 master/all-node mixed 误差 3.41061e-13 |
| uniform 对照 | XY/RZ 与 I05 公有入口逐项比较 | 位级相同 |
| constant gradients | 三级 dyadic AMR，49 master vertices | mixed 精确为零 |
| 非解析 master first jet | 49 master、24 hanging、48 fine subedges，每条 coarse/fine 段 9 点 | C1 最大绝对差 3.55271e-15 |
| 独立五点公式 | Vandermonde 小系统，未调用 production 求导函数 | 最大绝对差 1.23457e-13 |
| node sensitivity | 292 次查询；全部 master 参数叠加扰动，mixed 扰动量级 1e7 | 最大缩放误差 2.58742e-15；绝对误差 4.65661e-10 |
| 支持完整的 AMR subset | 48 active leaves、57 masters | Q4 误差 1.47105e-15，C1 为零 |
| disconnected support | 两个独立条带共 24 line components；另一输入一侧只有三点 | 有效条带闭合通过；不足五点明确拒绝，未跨 gap |
| identity/isolation | 未绑定 candidate、mesh mutation、primary replacement、已发布候选、缺失/多余/fine 数据、axis、极值溢出 | 拒绝符合预期；源和候选不变 |

局部 first jet 来自同一 candidate；源 first jet 与库存逐项保持原值。hanging 仅进入 Core 代数约束和采样，不进入 primary/独立 transport target。

在本次合法 2:1 网格上，完整 hanging jet 有 40 项非零 mixed 依赖，但闭合选取的切向 samples 不依赖 mixed，因此 compatibility matrix 为单位阵（49 行、49 nonzeros）。实现保留方程 `(I-LC)q=La` 的所有非零依赖，未用任意 mixed seed 作为返回结果；但本报告**不将单位阵测试冒充非对角闭合系统验证**。代数残差和最终原公式代回残差分别检查，实测均为零。

## Master-only 输运接入

正 Q4 的真实旧 snapshot 经 I03/I04 非线性速度回溯，两个输运结果放入同时间层新 candidate，再调用 AMR closure。M=0 是该正 scalar 的真实负相库存；此处未更新 phase volume 或发布新 geometry。

| 测试 | master vertex / edge | hanging / fine subedge | 独立闭合公式最大差 | C1 差 |
| --- | --- | --- | --- | --- |
| XY 三级 AMR | 49 / 100 | 24 / 48 | 3.88578e-16 | 0 |
| RZ 三级 AMR | 55 / 106 | 18 / 36 | 1.11022e-16 | 0 |

结果绑定新 primary identity 和 End 时间 t=0.05；旧 source 的 snapshot 地址/identity 和 first jet 不变，新 candidate 未绑定 geometry。replacePrimary 产生新 identity 后，旧闭合结果没有静默重新绑定。

## 网格收敛

base N=8,16,32,64，内部固定区域 level-1 refinement。XY 为 `exp(0.7x+0.4y)`，RZ 为 `cos(0.7r)*exp(0.4z)`。最细网格有 7,233 master、128 hanging vertices。

| MMS | N=64 master error | N=64 hanging error | N=32 master 阶 | N=64 master 阶 |
| --- | --- | --- | --- | --- |
| XY | 1.33237e-10 | 2.98450e-12 | 3.97482 | 3.29657 |
| RZ | 3.62534e-11 | 8.86902e-13 | 3.96864 | 3.52505 |

最后两级均满足预设 `[2.7,4.3]`。AMR tangent trace 再求导可能引入三阶项；本轮没有将观察到的部分四阶区间推广为所有 AMR 点四阶。所有 axis mixed 精确零。均匀网格四阶由保持不变且再次通过的 I05 contract 单独证明；完整曲率和 Interface benchmark 仍需 I12。

## 回归与证据完整性

- Release 主程序及全部目标构建通过；首次全目标 33.938 s，无新增 compiler warning。
- CTest 10/10 PASS：help、architecture、旧 Interface state、I01–I05、I06a、I06b，总计 10.46 s；AMR contract 1.31 s。
- 旧 parametric control 与 Core GALS control 均退出 0；GALS 完整输出逐字不变。
- 原 native geometry/AMR contract 完整输出、I06a skeleton contract 完整输出，与 I06a 归档逐字一致；I03/I04 的 23 行数值记录亦相同。
- 1,555 个既有源码/输入/工具文件哈希不变；S08 69/69、方法 PDF 3/3 一致。

失败记录保留：最初测试错误地要求“闭合矩阵必须有非对角项”，实际合法网格给出单位阵；用于扩大覆盖的 N=3 嵌套 refine-box 组合被拒绝为非 2:1，继续细化的试验又生成不合法 Cartesian side。这些 fixture 没有进入最终验收，不作为有效网格覆盖。最终保留原合法三级网格，记录实际依赖，不降低 Q4/C1/空间阶/残差门槛。极值输入实际先触发 Core expansion overflow，测试改为核对这个上游拒绝，未改生产防御行为。

`verification.json`、逐文件 delivery/runtime hashes、原始日志、source audit 和 evidence archive 给出可追溯记录。Git 元数据不可用。

## I06 完成审计与后续

I06a 已证明 master trace 的精确 hanging/fine restriction；I06b 已实现同 revision 的非均匀 mixed closure，并验证 Q4、MMS、coarse/fine C1 与 I03/I04 master-only 接线。路线图 I06 的两项唯一验收点均有代码与执行证据，故父任务 I06 PASS。

下一项按路线图为 I07。I07–I10 的守恒库存/唯一 flux/reflux/bubble、I11 唯一 geometry publication、I12 完整 Interface benchmark 以及后续低黏/无黏 Hydro/EHD 验收均未完成。当前测试不证明 dynamic regrid/coarsening 事务、contact/sector、physical inflow history 或最终低黏稳定性；S01 已知 `hydro_moment_xy` 首步失败仍保留。
