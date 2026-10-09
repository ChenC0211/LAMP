# I10 / 单元内部 bubble 体积闭合交付报告

**状态：PASS。日期：2026-10-01。下一项：I11。** 总体 goal 与 Gate I 仍未完成。

## 完成内容

新增 V2 `closeHermiteGalsBubble`：固定候选 primary first jet 和 I06 mixed closure，仅求派生 bubble c，使每个单元几何体积在具名容限内闭合到其独立 M。采用明确的 zero-center cardinal extension 基准；该基准下 c 才等于 Q4 center scalar。候选与预测输入不变，结果没有生产 geometry/accepted state 所有权。

Core 已有 `nativeQuarticCenterInterval` 为整个修正区间提供 Bernstein 零集管道 signed-partial 与位移证书；V2 使用该安全区间、实测端点体积括区及 safeguarded Newton/二分。所有迭代都复核同源体积、导数、拓扑、求积估计，超界/不可达/预算耗尽明确拒绝。没有 M clipping、全局 KKT、legacy fallback 或超预算的另一个数值任务。

改动为 **4 production + 2 tests**：新 `HermiteGalsBubble.h/.cpp`；共用 `InterfaceCellMeasure.h/.cpp` 新增完整 Q4 bubble 体积/导数入口和逐 cell moment estimate；新 contract 与 CMake 注册。Core、mesh、I01–I09 数值公式及三条冻结实现不改。

旧 `discoverVolumes` 的拓扑/积分算式只提取为共用模板；false 实例保留原路径，true 实例使用既有 sensitivity 积分。移除新增分支、逐 cell 估计汇总和 wrapper 后，整个 cpp 与 I07 交付文本在忽略空白后完全一致，见 `measure_arithmetic_comparison.json`。

任务格式见 [TODO](I10_Interior_Bubble_Closure_TODO.md)，详细语义见 [API](../api/I10_Hermite_Gals_Bubble_API.md)。

## 数值验收

原始数据位于 `runs/baselines/I10_20261001/`，包括 `closure.csv`、`pure.csv`、`bubble_final.log`、`ctest_details.log`。

| 验收项 | 结果 |
| --- | --- |
| T0 Release | 全部目标构建成功，完整构建 76.65 s，无 compiler warning |
| CTest | 14/14 PASS，44.17 s；I10 contract 23.82 s |
| T2 闭合矩阵 | XY/RZ × uniform/AMR × 安全区间 0.3/0.95 幅度，共 8 组、138 个 cut-cell 实例 |
| 初始解析体积 | 圆 `pi*R^2` / 球 `4*pi*R^3/3`，R=0.27，误差低于预设 `1e-10` 总容量尺度 |
| 最终逐 cell 相对 residual | 最大 `9.654373478435894e-11`，门槛 `1e-10` |
| 独立解析截面积分 residual | 最大 `9.6544111738987855e-11`，门槛 `1e-10` |
| 导数中心差分归一误差 | 最大 `1.0841947893358564e-8`，门槛 `2e-5` |
| 单调扫描 | 每种 XY/RZ/网格组合在 5 个系数位置扫描，M 非增 |
| 所有试算的相对 moment estimate | 最大 `2.5202479699273835e-12`，门槛 `1e-11` |
| 位移上界 / 最小 cell 尺寸 | 最大 `0.11867111471400293`，不超过 0.125 |
| safeguarded root | 388 个 Newton、128 个二分、2 个端点步；最多 7 次更新批次 |
| 纯相与 ownership | 4 组全液/全气 XY/RZ；系数精确不变，约 `1e-12` 的库存 tail 保留且报告；closure 在源析构后可读 |
| 拒绝路径 | 21 次（含各纯相组重复）：错误尺寸/非有限值/revision/mesh/config、容量越界、不可达 M、迭代耗尽、纯相造分量、moved-from、无 signed-partial 证书及未解析 closed cycle |

导数误差定义为 `abs(dM_dc_FD-dM_dc_integral)*(c_upper-c_lower)/capacity`；FD 步长为初始证书区间宽度的 `0.001`。单调扫描位置为预测点及正负方向的 0.25/0.5 区间幅度。两种目标幅度都在不同 cell 上交替选择正/负方向，不仅测试单侧增减。

独立 reference 直接积分解析 Q2 圆/球加已知 `delta_c*psi`，不采样 Snapshot、不使用 Core roots/moments。选取证书方向做截面二分求根，再 adaptive Simpson；以解析边界交点分段防止遗漏小角切区域。最终还通过独立构造的 V2 fixture 调用 I07 体积 API，检查几何—库存 mismatch。该 fixture 是测试用发布，未实现 I11 production publisher。

结构上全部 master nodes/edges/mixed 输入固定，bubble 基函数及其梯度在边界恒零；数值上在每 cell 四边抽样对比 value/完整 gradient，尺度误差小于 `2e-13`。AMR 使用同一 hanging constraint engine，RZ 保留轴对称边界 trace。M 数组、candidate identity、未绑定 geometry 状态均保持。

## 回归、审计与失败记录

- I07 volume、I08 flux、I09 reflux、I06a skeleton、旧 native state 的完整输出与 I09 归档逐字相同；I03/I04 共 23 行数值一致。
- 旧 GALS 与 parametric state control 完整输出相同，退出 0；GALS SHA256 仍为 `C6D8D43979F6A5998E9C54B292EA7EE38EA627A6352E66F9CB0BDBA779AAA5D1`。
- 额外运行旧 `interface_conservative_contract --gals-study`，原 bubble/shared closure 控制通过；这是现版本回归通过，不声称存在未保存的完整数值 A/B。
- 基线 1571 文件中 1568 不变、3 个既有文件变化，另新增 3 个源码/测试文件。69 个 S08 专用源、3 份 PDF 哈希不变；无越界代码修改。
- 本轮增量构建无 warning；既往冻结 Hydro 的 C4244 仍属历史已知项，不声称已修复。S01 `hydro_moment_xy` 的已知首步失败也未重分类。
- 初次编译因测试误用 `ScalarJet2.gradient` 而失败，改为真实 `dx/dy` 字段后通过。
- 初次完整测试试图先用 I07 测量单 cell 内的小闭合圈以制造 target，提前触发 Core 的 `closed subcell cycle` 拒绝。随后改用解析 target，使测试直接验证 V2 的明确拒绝；另外直接验证已有 Core center guard 无 signed-partial 证书。生产实现未因该失败缩小容限或加入后备算法。早期日志完整保留。

## 未完成范围

本任务只闭合已给定 skeleton/预测与 M，不证明整条 GALS 时间推进或真实拓扑变化。guard 是充分安全条件，不能关闭所有数学上存在的体积目标；过大修正、未解析小闭合圈或证书不足需要调用方缩步/细化。当前仍按批次构建全域临时几何，局部缓存复用/性能工程未开展。

Core 的半面内部精确相切/重根、正长度零迹及未解析 closed subcell cycle 限制仍在。I09 暴露的任意步内整面穿越时间求积问题仍需 I12 收敛/有界性验证。I11 唯一权威几何发布、I12 完整 Interface benchmark、后续 Hydro/Electricity/EHD 及低黏/无黏电喷雾稳定性均未完成。

## 交付归档

`verification.json`、`source_audit.json`、`regression_comparison.json` 汇总验收；`delivery_manifest.json` 记录交付与 archive SHA256；`evidence_inventory.csv` 和 `archive_verification.json` 校验所有归档条目以及 live/staged 源码、运行程序一致性。保留构建/测试/早期失败日志，Git commit 为 `unavailable`。
