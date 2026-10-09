# I12a10r1 / AMR 平面对称性丢失定位

## 目标：一个事实

定位 I12a10 的 12 个 AMR 正向平移失败中，exact outer symmetry 在“旧 Q4 → 原 B/J first-jet → I06 mixed → 新 Q4”哪一层丢失。只做可复现诊断，不先假定是 roundoff、BC 或 mixed 算法错误。

## 范围与禁止项

0 production；最多 2 测试相关文件：benchmark 与 CMake。另含本 TODO/报告、`runs/baselines/I12a10r1_20261001`。禁止改变现有推进结果、BC 声明、任何 gate、源 scalar、accepted M/first jets；禁止用诊断修补正式候选。

## 实现

- 显式 `--plane-symmetry-audit` 重放 12 个失败配置，并含负向/均匀对照。检查旧 geometry 与预测 candidate 身份/M 不变。
- 在真实 transverse 外边界检查 vertex/edge first jets 和 mixed 值；对原回溯 departure、同源 scalar jet、J_B 逐项记录。Q4 的对应 5 个边界 normal-derivative DOF 全为零当且仅当该边界导数恒零，使用公开 nativeQuarticLocalDofs 核验，不用容差冒充 exact symmetry。
- 允许独立 private 反事实快照只用于诊断“边界 normal first jets 置零”“边界 mixed 置零”各自影响；不发布，不返回给 advance，不改变 interior/原候选，不拿反事实成功代替正式推进验收。只有证据确定需要时才实施该对照。
- 确认原拒绝来自哪些实际边界/单元，避免把原异常上下文 cell=27 当作已定位故障单元。

## 验收与交付

新诊断可执行且解释全部 12 个拒绝，给出与真实标量/候选同源的数值证据；原 crossing publication 仍 36/48，旧直接回归输出保持。所有 production、未改可执行、81 冻结/参考、3 PDF 与 I12a10 archive 哈希保持。诊断 PASS 不将 I12a10/I12 标 PASS。若证据不足即 FAIL，不能直接修改规则。报告只提出由证据支持的最小下一项修复，保留日志/源码/runtime 归档。
