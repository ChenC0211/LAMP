# I12r2 / bubble 拒绝的可达性审计

## 任务 ID / 标题

I12r2：量化 I12 原始运动 benchmark 的 bubble 拒绝，区分已用安全区间不足、可证明的安全外部根与尚无证据的情形。前置 I12r1 PASS 已核验；上一 goal 轮为 progress，父任务 I12 仍 FAIL。

**结束状态：PASS（审计完成）；父任务 I12 仍 FAIL。** 全部 8 个失败已复现，4 个存在通过同伦证书的外部根、2 个 pure mismatch、2 个区间外探测仍未解决；生产代码未改。详见 [报告](I12r2_Bubble_Rejection_Audit_Report.md)。

## 目标

只解决一个事实：原始失败候选的 M 在既有 bubble 区间中的位置，以及是否存在同 skeleton、满足原体积/位移门槛且有既有 Core 同伦证书的外部根。诊断结果用于决定后续 I12 修复范围，不把测试探测当作 production 求解。

## 允许修改

- `tests/interface_hermite_gals_benchmark.cpp`：保留原步进及拒绝，增加显式 CLI 的只读失败探针，复用公共 Core/Interface 评估和证书 API。
- 必要时 `tests/interface_hermite_gals_bubble_contract.cpp`：独立复现证书比较；最多 2 tests。
- 本 TODO/报告、docs README、相关 API 的验证边界；`runs/baselines/I12r2_20261001/`。

## 禁止修改

所有 production、Core 证书、I10 配置/容限、时间取样、旧 frozen paths、方法 PDF、App/reader。不实施待确认的 phase-flux 时间积分设计变更。禁止 M clipping、解析参考回灌推进、接受区间外试算、隐式 fallback 或将父 benchmark FAIL 改成 PASS。

## 数学约束

1. 使用真实失败候选的 first jets、mixed、predicted center 和守恒 M；不重置为解析解。
2. 只变动一个失败 cell 的 interior center，边界及相邻 skeleton 固定。沿用同一 zero-center cardinal gauge。
3. 记录 `M(c)` 非增、原证书、端点体积、误差/capacity、moment estimate/capacity；不得把积分估计当严格真误差界。
4. 正常步进保持旧 accepted identity/M；探针得到的 geometry 不返回 production，不推进时钟。
5. `nativeQuarticChangeCertificate` 成功只证明其声明的 regularity/displacement；失败为 inconclusive，不证明不可达，也不自动替代 retainedDerivativeFraction 判据。
6. 圆/球平移与旋转使用独立解析 cell-volume reference；vortex 无独立逐 cell 真值时显式 NaN。

## 实现要求

显式 probe 模式在原 `finish` 失败后检查候选，保留原 exception。报告 cell bounds、时间层、pure/cut、target/predicted/reference volume、预测系数、原安全区间与端点体积。对 cut cell 仅作有上限的单方向扩区和二分诊断，同时检查 topology discovery、同伦证书及未改变的位移/体积门槛；任何证书不足/拓扑发现失败都明确报告。记录默认运行无 probe 时的原结果，防止诊断改变数值路径。

## 验收

- T0：Release benchmark 与 architecture 构建/检查；production 哈希不变，无新增 warning，文件预算符合。
- T1：原失败目标位于测得 bracket 外或 pure mismatch 与原 exception 一致；独立圆/球参考、端点单调性、试算误差、证书字段可复核；安全根只在实际证书与体积门槛通过时记录。
- T2：覆盖 I12r1 全部 8 个原始 moving failures（XY/RZ、uniform/AMR、rotation/vortex），保持 accepted time、失败阶段与 accepted M；先前 plane crossing 控制不变。
- 回归：默认 benchmark 的数值行与 I12r1 相同（不以运行时间作数值比较）；既有可执行文件哈希不变的前序算子不重复运行冒充新验证。

## 失败规则

没有完整证书不能宣称 guard 过严或目标不可能；诊断本身失败须记录为 unresolved。若发现可修复 production 缺陷，形成下一独立任务，不在本项越界改方法/容限。任何时间离散设计变更仍需明确决策。

## 交付

显式探针、逐失败分类数据、原行为不变证据、源码/runtime 哈希、完整成功/失败日志、报告与后续建议。本任务结束状态仅表示审计是否完成，不表示 I12/Gate I 完成。
