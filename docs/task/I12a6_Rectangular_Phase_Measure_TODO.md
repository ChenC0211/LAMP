# I12a6 / 同源 native Q4 矩形子区域液量

## 一个数学事实

从同一旧 native Q4 snapshot 的完整界面分解测量任意受支持矩形与负相的交集，不能通过重新插值出另一份 scalar 或以采样遗漏小组件。此项只提供空间测量，不宣称已解决非二进制扫掠端点或一般时空事件。

## 范围

最多 2 个 production 文件：新增 Modules/Interface/HermiteGalsPhaseRegion.h/.cpp；最多 2 个测试相关文件：新增 region contract 和 tests/CMakeLists.txt。另含本 TODO、API/报告、runs/baselines/I12a6_20261001。I12a5 已 PASS 并归档。

准备对象持有原 state/view 的共享 immutable owner，按已有 complete NativeComponents 发现全部 graph；query 必须验证精确 mesh identity。仅支持 native enriched Q4 与真实 XY/RZ 度量，接受原来的显式边界 contacts/symmetry planes。原 M 不用于测量、更不覆盖。

## 实现合同

1. 矩形端点是明确的有限物理 double，宽高严格正、完全位于 active scalar 支撑内。检查覆盖不能用浮点总面积接近代替拓扑无缝覆盖；空洞/域外明确拒绝。无外推。
2. 将矩形与原实际叶片相交，纯相使用同源 certificate。切割叶片保留原 complete graphs；按矩形切线将整个叶片划分为目标与补集的凸区域，不生成新 scalar。给所有外部区域边界供应实际 retained cuts，不遗漏曲线 owner。
3. 使用 I12a5 局部原点及同原点 fixed-boundary moments；双 Green 公式、负/正相容量与误差估计诊断。确定性地按矩形宽高选择主矩公式，较小相的直接测量用于稳定互补表达；不得裁剪负液量、修改原库存或用另一几何 fallback。
4. Core 请求根细化时更新指定 retained cut 并重做所有相关 fixed/free moments，使用有界工作量。未解析事件、正长度零迹或工作上限明确拒绝；不把粗细求积接近当拓扑完整证明。
5. 默认相对 ledger/浮点误差估计门槛 1e-10，与现有 I07 分开记录；不是严格误差包围，更不是允许 I12a3 放宽舍入预算。若严格容量或控制失败，留存证据并另拆修复。

## 验收

- T0：Release 主程序、相关 targets、architecture、81 冻结/参考和 3 PDF、前序 archive、2+2 文件预算。
- T1：XY/RZ、均匀/AMR、矩形跨叶片/完全包含/薄条、纯相、圆/环形及反相；独立解析矩形交集，分割可加性和同源 identity/生命期、原 M 不变。
- T1：非法/越界/空洞、错误与过期 mesh、moved-from handles、已有未支持零迹明确拒绝；失败隔离。
- 直接回归 I07、I10/I11、I12a1–a5、architecture/CLI；原 crossing probe 保留真实失败。未接 driver，不宣称 I12 PASS。

## 后续

本项通过后继续常速度 curved Q4 扫掠积分与边界历史身份，再推进变化速度与完整半/全步消费。不能以此局部算子替代 I12 benchmark。

## 执行结果

**PASS。** 2 production + 2 test 文件；212 组测量、11 类拒绝及 11 项直接回归通过，9 项旧输出相同。partial active domain 的纯场在 prepare 阶段保留原 Core 显式接触要求，未伪造接触点以绕过。见 [报告](I12a6_Rectangular_Phase_Measure_Report.md)。
