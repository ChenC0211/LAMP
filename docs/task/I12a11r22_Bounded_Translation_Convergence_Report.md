# I12a11r22 / 有界平移非多项式收敛报告

**状态：FAIL。固定的几何/Hessian 预算未使原空间矩阵全部完成，I12 / Gate I 仍 FAIL。** 依据 [TODO](I12a11r22_Bounded_Translation_Convergence_TODO.md)，未改变 production 算式、r21 预算、原网格/dt/终点或精度 gate。

全部计划运行与回归已于 2026-10-02 完成；不存在以未完成运行代替结论的情况。

## 执行范围

仅修改 `tests/interface_hermite_gals_benchmark.cpp`，让原空间研究及 r13 时间控制选择 bounded 契约；原入口保持默认旧契约。新增命令 `--bounded-translation-space-study`、`--bounded-translation-dt-half`、`--bounded-translation-dt-quarter`。

几何容限固定 1e-10，每叶 Hessian 扰动预算仍为 `min(hx,hy)^2`。初态为原 `exp(q)-1`，first jets/center 与 M 分别沿原输运及新共享库存通量推进，半步/全步由原 I10/I11 发布。没有解析重置、clipping、临时缩步或放宽预算。

## 空间矩阵

原 12 组全部执行，6 组到达 t=.04，6 组拒绝，四条三档序列均不完整。空间入口 exit 1，用时 1122.2157035 s；“单组完成”与“空间收敛通过”分别记录。

| 场景 | n | 接受终点 | 结果 |
| --- | --- | --- | --- |
| XY uniform | 8 | .04 | 完成 |
| XY uniform | 16 | 0 | 首步 end closure 拒绝，cell 89 |
| XY uniform | 32 | .04 | 完成 |
| XY AMR | 8 | 0 | 首步 end closure 拒绝，cell 58 |
| XY AMR | 16 / 32 | .04 | 两组完成 |
| RZ uniform | 8 | .01 | 第 2 步 midpoint closure 拒绝，cell 49 |
| RZ uniform | 16 | 0 | 首步 midpoint closure 拒绝，cell 82 |
| RZ uniform | 32 | .04 | 完成 |
| RZ AMR | 8 | .01 | 第 2 步 midpoint closure 拒绝，cell 33 |
| RZ AMR | 16 | 0 | 首步 midpoint closure 拒绝，cell 74 |
| RZ AMR | 32 | .04 | 完成 |

6 次拒绝均发生在 bubble 的 `volume bracket` 阶段，标签为 `target outside Hessian accuracy/geometry volume budgets`。这是原几何证书与新增导数预算交集中的可达性拒绝；现有日志没有逐一分离 Core 证书边界和 Hessian 边界的贡献，不能仅凭标签断言每次都由同一边界主导。失败保留了原 accepted identity/M 和真实已接受时刻。未出现库存装配越界拒绝。

相同 n=32、dt=.01、终点 .04 的曲率最大误差比较：

| 场景 | r8 旧严格闭合 | 本项 bounded 闭合 |
| --- | --- | --- |
| XY uniform | 1.0386418829e-2 | 7.0989368016e-5 |
| XY AMR | 1.1499912514e-3 | 6.8005510796e-6 |
| RZ uniform | 8.1841629536e-5 | 3.4115159258e-5 |
| RZ AMR | 8.1841629536e-5 | 3.4115159258e-5 |

XY uniform 曲率误差减少约 146.31 倍，梯度误差由 2.6939273783e-5 降为 1.1316375575e-7。RZ AMR 的梯度误差则由 1.1170564870e-7 略增至 1.1190483182e-7；不能声称所有指标都改善。粗网格拒绝使原逐档递减和末档梯度/曲率阶数要求尚未满足，不跳过失败档计算通过阶数。旧 r8 路径虽然完成 12 组，但 XY 空间收敛同样 FAIL。

## 时间控制与验收记录

同一 XY uniform n=32、终点 .04，dt=.005 的 8 步完成，耗时 226.6230944 s。曲率误差 7.2923855922e-5，相对本项 dt=.01 的结果增加约 2.7%；梯度误差 1.3666982026e-7 也增加。旧同参数严格闭合路径的曲率误差约 .4572，较旧粗时间步增加约 44 倍；当前结果显著改善，但仍不能宣布时间收敛。

dt=.0025 的 16 步也完成到 .04，耗时 476.4704196 s。旧路径在 t=.025 后拒绝，本项通过了那个位置；两者终点不同，不计算错误的同终点改善倍数。最小 dt 的曲率误差较本项粗 dt 增加 8.24%，梯度误差增加 33.87%。

| dt | 步数 | 接受终点 | gradient | curvature |
| --- | --- | --- | --- | --- |
| .01 | 4 | .04 | 1.1316375575e-7 | 7.0989368016e-5 |
| .005 | 8 | .04 | 1.3666982026e-7 | 7.2923855922e-5 |
| .0025 | 16 | .04 | 1.5149305516e-7 | 7.6841402082e-5 |

恒定平移的离散回溯与区间积分本身精确，此时间控制检验重构/闭合跨步行为，不替代原一般运动时间阶 >=1.8 的验收。

116 条成功闭合记录的最大几何残差 9.99255565445334e-11、pure 残差 8.91017121196214e-11；实际 Hessian 扰动仍在所选预算内，最大 .015625。108 条已生成区间面记录的库存 ledger 最大绝对值 1.38777878078145e-17，reflux 为 0；原 roundoff 面修正最大相对值 2.22044604925031e-16，单阶段最多修复 2 面。守恒约束与几何误差分别验证。

首次构建通过，4.2571775 s，无 warning。4 项 CTest 全通过，64.2385587 s，全部输出与 r21 逐字相同；旧 crossing probe 保持相同输出和 exit 1。回归与 dt_quarter 独立进程重叠执行，未重建运行中的程序，耗时不用于性能比较。production 无修改，其余 32 个运行目标哈希保持。81 个冻结文件、3 份 PDF 和 r21 archive hash 不变。

证据目录：[I12a11r22_20261002](../../runs/baselines/I12a11r22_20261002/)。完整 space、dt_half、dt_quarter 日志记录实际过程；`convergence_comparison.json` 对照 r8/r13 并区分相同终点与提前拒绝，禁止把提前停止的误差当作更准确结果。`numerical_verification.json` / `verification.json` 为任务 FAIL，T0/来源完整性检查单独为 PASS；`delivery_manifest.json` 和 `archive_verification.json` 核验 `evidence.zip`，归档完整性 PASS 不等于数值任务 PASS。

## 后续依据

取消机器精度匹配并保持库存有界，已经改善部分细网格曲率；固定 1e-10 体积差门槛与 O(h²) scalar Hessian 预算仍不能使完整矩阵闭合。下一项应检查这些预算与网格截断误差、实际可行区间的关系，先明确数学依据和验收，再实施。不能把本次失败改写为通过，也不在本项增大预算。一般 rotation/vortex 区间积分、完整时间阶和动态 AMR 仍未完成。
