# I11 / 唯一权威 Hermite 几何发布交付报告

**状态：PASS。日期：2026-10-01。下一项：I12。** 总体 goal 与 Gate I 仍未完成。

## 完成内容

新增 `publishHermiteGalsGeometry` 与不可伪造的只读 Publication。入口消费同 revision 的 I10 bubble 结果，构造并验收唯一 native Q4 snapshot，向体积、完整零集、normal/curvature、半面 cuts 和 I08 midpoint phase flux 提供相同 scalar owner。保留 primary first jet/M、time、measure、geometry identity；不产生第二份 production geometry、不裁剪 M、不提交 coupled accepted step。

发布在私有同 revision candidate 副本上完成。所有构造和检查完成后才以 noexcept 步骤消耗输入候选；失败保留原候选供重试，旧 accepted owner 不变。RZ 使用同源子午曲率加 `n_r/r`，真实 axis endpoint 采用偶对称极限，无拟合或近轴阈值。

修改 **4 production + 2 tests**：新增 Publication h/cpp；State h/cpp 仅增加私有事务入口及注释；新增 contract 和 CMake 注册。移除三段新增私有定义后，State.cpp 全文与 I01 交付文件在忽略空白后完全一致。Core、I02–I10 算式、App 和所有冻结实现不改。

任务格式见 [TODO](I11_Authoritative_Geometry_Publication_TODO.md)，接口语义见 [API](../api/I11_Hermite_Gals_Publication_API.md)。

## 验收结果

原始证据位于 `runs/baselines/I11_20261001/`，主要数据为 `publication.csv`、`publication_final.log`、`numerical_summary.json`、`ctest_details.log`。

| 项目 | 结果 |
| --- | --- |
| Release 构建 | 主程序、前序和本项目标全部成功；完整回归目标构建 36.90 s，无 compiler warning |
| CTest | 15/15 PASS，总计 51.63 s；I11 8.32 s |
| 解析矩阵 | XY/RZ × uniform/2:1 AMR × 正反 phase 共 8 组圆/球；另有 XY/RZ AMR 双分量，共 10 组 |
| 同源采样 | 1424 个 GL8 节点，另逐 graph 两端及 fraction=0.37 抽样；12 个真实 RZ axis endpoint |
| 最大解析 normal 误差 | `2.7577367830154345e-15`，门槛 `2e-10` |
| 最大解析 meridional/physical curvature 误差 | `2.4336088699783431e-13`，门槛 `2e-8` |
| 最大解析总量误差 / 总容量 | `1.1102230246251565e-16`，门槛 `1e-10` |
| 非零 bubble | 18 cut cells；重新测量的最大逐 cell 库存相对 residual `4.6045389723303742e-11`，门槛 `1e-10` |
| 纯相 | XY/RZ 全液/全气 4 组，空 topology、无伪曲率 |
| I08 链路 | XY/RZ × uniform/AMR 4 组 published midpoint；Q/M 结果保留同一 snapshot，零速度 M 精确不变 |
| 事务/ownership | 输入深拷贝、source 指针恒等、结果 copy/move、输入/mesh 析构后 surface 可读、stale revision/mesh、重复发布、严格晚期体积 gate 失败后重试均通过 |
| 拒绝记录 | 共 76 次，含多组重复的 consumed/moved-from/sample index/fraction；及 config、mesh、revision、晚期体积失败 |

单圆半径 0.27，XY 中心 `(0.47,0.51)`、RZ 中心 `(0,0.51)`；双圆半径 0.13，XY 中心 `(0.31,0.49)/(0.73,0.49)`、RZ 中心 `(0,0.31)/(0,0.73)`。双分量以两圆 level-set 乘积构造全 Q4。解析几何测试由 I07 测量 fixture 初始化独立 M，再经 I06/I10/I11 发布；解析法向、曲率和总体积独立于测量代码。

非零 bubble 测试以 I10 安全区间 0.3 幅度交替改变 cut-cell 系数，测量目标 M，再从原 center prediction 闭合。先用 `1e-18` 严格体积门槛强制晚期失败，确认候选未消耗、未绑定、primary/M 不变；同一候选使用预先规定的默认 `1e-10` 门槛成功。并非因测试失败而放宽任务验收门槛。

## 回归与审计

- I07、I08、I09、I10、I06a 与旧 native state 的完整输出和 I10 归档逐字相同；I03/I04 的 23 行数值指标相同。
- 旧 GALS 与 parametric control 均退出 0，完整文件哈希不变；GALS SHA256 为 `C6D8D43979F6A5998E9C54B292EA7EE38EA627A6352E66F9CB0BDBA779AAA5D1`。
- 前置 1574 个文件中 1571 个未变，3 个既有文件按允许范围修改，另新增 3 个源码/测试文件。69 个 S08 专用源与 3 PDF 哈希不变，无越界代码改动。
- 首次构建成功；首次运行在 `validateMesh` 检出 mapped cache 的 domain 地址不匹配。根因是 cache 在临时 domain 上绑定后随 Data 移动。已改为先在最终 Data 中存放 domain，再原位构造 partitions；首轮失败日志保留，专项及最终回归通过。
- 新增采样生命周期测试验证网格析构后的非空 surface，另验证已发布状态对 changed mesh 的拒绝。当前增量构建没有 warning；既往冻结 Hydro C4244 不视为已修复。

## 未完成范围

此项完成唯一权威几何发布，不证明整条 Interface 时间推进或真实颈缩破碎。I12 尚须完成 static、translation、rotation、vortex、AMR crossing 和 curvature convergence。

Core 的半面内部精确相切/重根、正长度零迹与未解析 closed subcell cycle 限制仍在；I10 的充分 signed-partial/位移 guard 可能拒绝可行修正。任意步内跨面通量的时间精度和有界性仍未证实。派生缓存复用和局部重建性能优化未开展。S01 `hydro_moment_xy` 首步已知失败、后续 Hydro/Electricity/EHD 及低黏/无黏电喷雾稳定性均未解决。

## 交付归档

`verification.json`、`source_audit.json`、`state_arithmetic_comparison.json`、`regression_comparison.json` 汇总检查；`delivery_manifest.json` 记录交付文件和 archive SHA256；`evidence_inventory.csv` 与 `archive_verification.json` 验证全部 zip 条目、live/staged 源码和运行程序。Git commit 为 `unavailable`。本项通过后才进入 I12，不提前接 Hydro 新主线。
