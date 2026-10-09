# I12a10 / 任意步内平面跨面的完整发布

## 目标：一个数学事实

原 q=1/4、3/4 的跨面反例必须在新的同源 half/full 通量 → 独立 M → I09 → 原 I10/I11 完整链路通过；算子代数通过不能代替实际 first-jet 与发布。I12a9 常速度圆/球 40 步范围已由 I12a9r2 验收归档。

## 允许与禁止修改

0 production；最多 2 测试相关文件：interface_hermite_gals_benchmark.cpp、tests/CMakeLists.txt。文档及 `runs/baselines/I12a10_20261001`。禁止改 B/J、相积分、bubble、任何容限、冻结代码/PDF 或 Hydro。保留原 --crossing-probe 为显式 midpoint 失败对照。

## 实现与数学约束

- 为现有制造解 Case 增加明确的轴对齐平面与 XY 两个平移方向，默认圆/球和旧输出保持。初始 I07 初始化 M 一次；后续 interior first jets 和 center 仅消费原 accepted scalar 经原 characteristic。
- 平面 transverse 物理边界使用已有 symmetry-plane 声明，RZ 轴沿用原真实 axis 规则；同一配置传递给 I07、region、I10 和 I11。不提供任何 interior 解析重置。
- 实际入流 BC 从明确的平面外边界历史给出：整面相态与整个 base time window 常值且有严格同号证明，允许 Liquid/Gas。仍共享 I12a8 的 immutable BC owner。
- 复用完整 advanceTranslation：midpoint 和 end 都必须原样闭合发布，M 与共享通量保持逐项相同。失败检查旧 accepted identity/M 不变。

## 验收

- T1/T2：XY x/y、RZ z × uniform/AMR × 正反速度 × 正反相 × q=1/4、3/4，共 48 组。n=4、dt=.02、速度幅值 1，实际穿越原 native face .5，不缩小 dt 对齐事件。核对真实半/全时间、解析逐单元 M、first jets、正常/曲率、全区间 crossing face 四个 gross volume 与边界 ledger。
- q=1/2 的已知 midpoint 整面零迹不以此次非对齐反例代替支持；零迹和其它拓扑项目保留独立 gate。
- T0：Release benchmark，新增 crossing publication、constant translation smoke、circle reference、直接区间/库存/reflux/边界合同、architecture/CLI；受影响旧测试输出逐字比较，production 和其它 runtime hash、81 冻结/参考 + 3 PDF、I12a9r2 archive 不变。

## 失败规则与交付

任一 I10/I11、库存或几何 gate 失败均保留真实原因，不更改 gate 或用解析 M 覆盖。缺失 Core 能力另列有证据的单项修复，不混入本项。交付最小测试接线、48 组原始数据、哈希审计/归档与下一项；I12 整体必须继续完成非恒定速度与其余验收。

## 执行结果

FAIL：36/48 组完整 half/full 发布通过；12 个 AMR 正向平移在中点的 exact outer symmetry 检查拒绝。9 项旧回归输出逐字不变，production 未改。详见 [报告](I12a10_Planar_Crossing_Publication_Report.md)，继续定位边界 first-jet/mixed/Q4 的零法向导数丢失。
