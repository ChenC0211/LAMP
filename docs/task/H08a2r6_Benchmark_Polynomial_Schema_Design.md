# H08a2r6 / 导出次数一致的实施前设计

2026-10-05。生产 0，benchmark 1：只修改 `cases/benchmarks/cch_convection/driver.cpp`，文档单列。上一轮 H08a2r5 的 52 个冻结 payload、12 个当前交付文件与 1381 个源码/输入已复核。

实际障碍为三处：`trace-contract/spatial_zero_rate.csv` 的 c 列、`advector/advector_atoms.csv` 的 c/e 列、`boundary/boundary_atoms.csv` 的 c 列。数据已有逐系数遍历，header 却写死三项。新增一个按 coefficient 容器类型的 `tuple_size` 生成 header 的 helper 和对应逐值 writer，三处共同调用。类型直接取自 `CchNormalPolynomial::beginCoefficients/endCoefficients`；不创建虚假的高次速度或任何物理 owner。当前数组长度为 3 时，列名、顺序和数值字节保持；未来数组为 5 时同一源码生成全部 c0…c4/e0…e4。

先对当前源码标准构建 benchmark，复制实际 exe 与运行依赖，在新检查目录分别完整运行 trace-contract/advector/boundary，作为当次 before。三个模式原配置不缩短、不增添参数。完成 helper 后通过临时审计入口调用同一 header/data writer，输出三/五项、单端/双端及不同尾部字段的记录；独立 parser 按原真实 header 扩展期望列，并逐项核对非零 c3/c4/e3/e4 和尾部字段。临时入口随后归档移除。

最终源码再标准构建并完整重跑原三模式。所有 CSV 校验 header 唯一、行字段数一致、before/after 文件集合一致；除既有 `advector_summary.csv` 的 prepare/evaluate 两个计时列，所有输出逐字保持。完整来源字段也比较，不能只比较速度数值。记录标准架构/目录检查；生产源码和六个现行短契约 exe 哈希保持，不以格式任务重复完整 H10b。

本项不迁移 production 多项式、不改变 phase/材料/KEP 或旧 benchmark gate。剩余固定次数包括 `CchNormalPolynomial` 数组/value/integral、`CchPhaseGeometry::ratePolynomial` 和 phase/momentum 两份参考；material/dual 已由 H08a2r5 关闭。下一任务直接升级这些表示消费者，再统一部署实际正则 B/atom/J 与 characteristic 证书。H08a/H11 状态不因 CSV 通过而改变。
