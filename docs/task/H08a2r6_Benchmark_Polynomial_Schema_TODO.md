# H08a2r6 / benchmark 法向多项式导出 schema

状态：**PASS（导出格式）**，见[报告](H08a2r6_Benchmark_Polynomial_Schema_Report.md)及[实施前设计](H08a2r6_Benchmark_Polynomial_Schema_Design.md)。原34表/716228行保持（两耗时列单列），18条三/五项回读及3/3架构目录通过。下一[H08a2r7](H08a2r7_Complete_Normal_Polynomial_TODO.md)升级完整生产表示。输入为 [H08a2r5](H08a2r5_Material_Reference_Polynomial_Degree_Report.md)，材料/dual 参考已接受完整四次；生产 B/atom/J 仍旧，H08a/H11 尚未完成或启动。

目标：关闭完整多项式表示升级前的 benchmark 列数障碍。当前三处 header 固定 `c0,c1,c2`（一处另含 `e0,e1,e2`），数据却遍历实际 coefficient array；直接升级数组会造成 header/row 错位。只解决真实导出 schema 随系数库存一致这一事实。

允许修改：production 0；benchmark/test 1：`cases/benchmarks/cch_convection/driver.cpp`，对应 README 属文档。临时审计归本任务 `runs/checks`，结束时移除临时入口。禁止扩张成新的调试程序或长期脚本。

禁止修改：数值生产代码、接受状态、Interface、材料 F/R/uR、KEP、原 benchmark 模式/gate/物理数据。不能填造高次系数或用旧 header 截断实际输出。历史冻结及其分析脚本只读。

数学约束：四材料硬约束保持；导出仅观察原 owner。三类数据为 `trace-contract` 的 `spatial_zero_rate.csv`、`advector` 的 atom 表、`boundary` 的 atom 表。统一 header helper 按实际 coefficient 容器长度生成完整 `c/e` 列；当前三项格式和值保持，无新速度或多项式 owner。

实现要求：先登记这三处 header 与 row 及剩余 phase/momentum 固定次数消费者。完整系数列由实际类型给出，不由 case 名或临时参数选择。用同一个 writer 的三/五项合成记录核对列名、顺序和末尾物理字段；合成数据只在可退役审计中使用。生产表示升级到五项后应直接输出 c0…c4/e0…e4。

验收：T0 标准 benchmark Release、架构/目录及相关短检查；当前生产/source 未改的哈希。T1 原三个导出模式在独立时间戳目录完整运行，所有实际行字段数与 header 一致、无重复列，原物理/来源数据保持；旧计时列按既有规则单列。T2 对三/五项共同 writer 的独立 parser 检查，所有系数及尾部字段按名称回读完全相同。无物理变化，不把本项代替 H10o/H10b-v2 复验。

失败规则：任何字段错位、截断、输出物理量变化或来源丢失均 FAIL；修正或退役候选，不放宽检查。只有真实外部依赖才 BLOCKED。

交付：一个通用 schema writer、原三模式和合成五项的可复算证据、source/runtime 冻结、剩余表示迁移清单。完成后应直接推进完整 `CchNormalPolynomial`/phase rate 表示及 phase/momentum 两份参考（预算 production ≤5/test ≤2），再同步部署新 B/atom/J。正式来源切换后，原 native/affine/H10b-v2 28 窗口 84 项仍为前置。
