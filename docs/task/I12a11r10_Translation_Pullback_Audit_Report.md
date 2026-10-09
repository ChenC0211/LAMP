# I12a11r10 / 平移 predictor 拉回与重构审计报告

**PASS（只读诊断）；I12 仍 FAIL。** 对应 [TODO](I12a11r10_Translation_Pullback_Audit_TODO.md)。第三步 predictor 大误差主要来自新 Q4 重构；仅替换 mixed derivative 的私有反事实仍保留大部分误差。

## 实现与检查

0 production + 1 test，新增 `--translation-pullback-audit` 和 `--translation-pullback-control`。在实际 predictor 界面求积节点，查询同一旧 accepted scalar 在真实 B 的二阶 jet，核验实际 owned constant velocity/map 的 B=x-uΔt、J=I。实际 vertex/edge first jets 及 center 与该拉回逐位一致。candidate/old identity 和 M 均不变；诊断异常不遮蔽原推进。

私有反事实 patch 只把 master vertex dxy 换为旧 scalar 在实际 B 的 dxy，其余 transported 数据不动；boundary hit 保持 actual mixed（n=32 有 33 个，控制 n=8 有 9 个）。这份 patch 只在实际 predictor 的同一组节点求值，不发现其界面、不发布、不推进。旧场共享边界的二阶导数按原 native scalar 的确定性选边查询，因此该反事实只用于排除“仅 mixed closure 导致全部误差”的解释，不是可直接采用的输运方法。

每个阶段记录四类误差的独立最大值，以及四个最差点的 actual / old pullback / private mixed / analytic 全部二阶 jets。所有 32 个最差点的 departure 均只有一个实际覆盖单元。逐点带符号曲率分解的闭合缺陷为零；不能将下表不同位置的最大值直接相减作为分解。

## 同一节点集上的最大曲率误差

| 步 / 阶段 | actual − analytic | old pullback − analytic | actual − old | private mixed − old |
|---|---:|---:|---:|---:|
| 1 / midpoint | 8.0572e-6 | 9.6089e-6 | 9.4460e-6 | 9.4479e-6 |
| 1 / end | 8.7629e-6 | 9.6008e-6 | 7.9277e-6 | 7.9315e-6 |
| 2 / midpoint | 1.2739e-4 | 5.0301e-5 | 1.2825e-4 | 1.4742e-4 |
| 2 / end | 2.0692e-4 | 4.7996e-5 | 2.0673e-4 | 2.0015e-4 |
| 3 / midpoint | 5.2019e-3 | 7.4011e-4 | 5.1173e-3 | 4.4924e-3 |
| 3 / end | 1.4937e-2 | 6.4749e-4 | 1.4919e-2 | 1.2882e-2 |
| 4 / midpoint | 8.3431e-3 | 4.7579e-3 | 6.8133e-3 | 6.6199e-3 |
| 4 / end | 1.3911e-2 | 4.5589e-3 | 1.3814e-2 | 1.1834e-2 |

8 个 actual predictor 曲率最大值与 r9 原始日志逐位一致，最终四步 BENCHMARK 行也逐字一致。原空间收敛失败仍存在。

## 最差点的带符号分解

第三步 end 最差点位于 cell 650，departure 也在 650：actual − analytic = 1.4936763e-2，其中 old pullback − analytic = 1.7684529e-5，actual − old = 1.4919079e-2。后者又分为 actual − private mixed = 2.2477696e-3 和 private mixed − old = 1.2671309e-2。这证明该点主要在新 Q4 重构时放大，而单独变更 mixed derivative 无法消除主要误差。

第二步 end 全局最大系数修正的身份已找到：cell 650，bounds 为 [.3125,.34375] × [.625,.65625]，容量 9.765625e-4，c 从 −7.7970379442e-3 变到 −7.7901356890e-3，Δc=6.9022552e-6。初始相对体积残差仅 6.2061467e-14，最终为零，但最终 dV/dc 只有 −8.7807392e-12。该单元随后成为第三步两个 predictor 的最差单元。

第三步 end 另一个最大修正位于 cell 363，Δc=−9.9896362e-6，初始相对体积残差 −1.1235457e-13，dV/dc=−1.0983509e-11；它随后成为第四步 end predictor 的最差单元。同点误差 −1.3910988e-2 中，旧场继承项仅 −9.6716105e-5，重构项 −1.3814272e-2。

这条单元关联支持继续检查弱灵敏度 bubble 修正的来源，但尚不证明初始体积残差来自库存、积分舍入还是重构几何差额，也未证明某种限制系数或步长调整能同时保持 M 与原精度 gate。本项没有做这些生产变更。

## 验证及交付

Release 构建 3.81 秒，无 warning。Q2 控制 2.12 秒、退出 0，两阶段的最大拉回/重构曲率误差为 1.5632e-13，低于 1e-10；B/J、transport 和逐点分解缺陷均为零，无 bubble 激活。

原非多项式四步 131.11 秒、退出 0；8 个 predictor、8 个 root set、389 个 corrected-cell report、32 个最差点完整，无诊断异常。11/11 直接回归通过，46.69 秒，各项输出与 r9 逐字一致；旧 midpoint crossing probe 保持预期退出 1。回归与长诊断有并行时段，耗时仅记录，不作严格性能比较。

其它 27 个可执行、全部 production、81 项冻结/参考文件、3 PDF 和 r9 archive 不变。审计脚本直接使用原始十进制日志作 Double 解析对照，避免 PowerShell 5 的 JSON Decimal 自动类型影响逐位数值比较。原始日志、字段表、阶段/单点分解、全部 root 数据、源码/runtime 哈希、交付清单与 archive 在 `runs/baselines/I12a11r10_20261001/`。

下一任务针对 cell 650/363 的实际闭合前体积残差作守恒分解，区分 old M 与旧几何差额、实际 sweep/修复/组装差额、精确平移区域与 predictor 的几何差额，再决定最小修复。I12/Gate I 仍 FAIL，不进入 Hydro。
