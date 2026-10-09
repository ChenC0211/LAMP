# I12a7 / 常速度 curved Q4 穿面扫掠积分

## 一个数学事实

对无切向分量的常速度平移，单面在真实时间区间的穿面相体积等于同一旧 scalar 在矩形扫掠区域的相体积。XY 两方向与 RZ 轴向成立。利用 I12a6 的完整曲线分解覆盖时间事件，不按中点相指示猜测无事件。

## 范围

最多 4 个 production 文件：Core LevelSetHermiteGeometry.h/.cpp、Interface HermiteGalsTemporalPhaseFlux.h/.cpp；最多 2 个测试相关文件：新增 curved sweep contract 与 tests/CMakeLists.txt。另含本 TODO/API/报告及 runs/baselines/I12a7_20261001。I12a6 已通过并归档。

## 合同

- 新入口显式接收同一旧 accepted state 的 prepared region source；几何指针与 primary identity 均匹配原 interval。仍使用 owned constant velocity 和原 RK2 characteristic identity。拒绝 opaque coverage callback、非零切向速度、未支持的 RZ 径向运动；零法向速度的穿面量严格为零。
- Core strict TU 从原 double 速度和时间的精确 expansion 构成扫掠位置，并输出包含精确区域的外矩形、被其包含的内矩形、精确几何 total volume 的浮点值以及坐标包围带体积。不能先将两个全局 departure 舍入后相减并忽略误差。此入口暂要求完整 active 原域，外矩形越界/无可分辨内矩形明确拒绝；不外推 BC。
- 内外矩形使用同源完整 Q4 测量。坐标包围带是几何界定；两次浮点液量本身不是严格误差界。用两次测量的平均表达相液量，并从较小相稳定形成另一相；纯相保持与独立 total volume 一致。记录坐标带和 GL8 估计，默认相对坐标/积分估计门槛 1e-10；它们不能代替或放宽 I12a3 的 64 eps 库存门槛。无法解析/严格 gross 越界直接拒绝。
- 每面四个 integrated volumes 只生成一次；rate 仅派生，不改库存/I09、midpoint/end 发布或旧 affine 算子。半段与全段保留各自 interval，不制造 quadrature-node snapshot。
- 本项不处理变化速度、带切向平移的平行四边形、域外边界历史或完整 driver。必须如实标注支持范围。

## 验收

- T0：主程序/直接回归、architecture、81 冻结/参考+3 PDF、4+2 文件预算、前序 archive。
- T1：独立解析圆/球交集，XY/RZ、均匀/AMR、正反法向、二进制和十进制时间、区间切分，严格 total/liquid 与旧 snapshot 不变。
- T1：完整时间窗包含小组件而中点/四分点均为气相的漏事件反例，完整扫掠必须得到非零相体积；不以 coarse/fine 时间采样一致代替事件覆盖。
- T1：同源/异源、错误 mesh、opaque/nonconstant/tangential/RZ-radial、域外扫掠和不可分辨区间拒绝。测量误差与坐标误差分开记录。
- 原 I12a1–a6、Core 几何和 I07–I11 受影响回归及旧 crossing probe；本项通过不代表 I12 benchmark PASS。

## 失败规则与后续

如库存/拓扑实际接线暴露其它不一致，保留当前容限并单独拆项；不由本项暗改 B/J 或 primary reconstruction。随后增加真实 BC 历史身份和消费，再推进完整半/全步与变化速度矩阵。

## 执行结果

**PASS。** 4 production + 2 test 文件；66 组检查、12 类拒绝、6 个完整漏事件反例与 24 项 CTest 通过，22 项旧输出一致。单叶片未解析闭合小环仍明确拒绝，未将其伪报为支持。见 [报告](I12a7_Curved_Translation_Sweep_Report.md)。
