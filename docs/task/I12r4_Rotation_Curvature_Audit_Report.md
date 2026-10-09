# I12r4 / 旋转曲率与库存误差来源审计

**状态：PASS（诊断审计）；父任务 I12：FAIL；Gate I 未通过。** 2026-10-01。

按 [TODO](I12r4_Rotation_Curvature_Audit_TODO.md) 完成来源定位：AMR rotation 首步大曲率来自满足当前库存目标的 interior bubble 修正；在已测节点，曲率公式和 bubble 的独立微分均核对通过。库存目标的主要误差来自瞬时 midpoint 时间取样。不过，将目标换为精确体积仍有一个候选无法通过现有证书，不能承诺时间积分修复就足以完成 I12。

本次没有修改 production 方法、容限或接受条件。只新增显式 `--rotation-curvature-audit` 测试入口；原完整 benchmark 输出逐字不变，仍以退出码 1 报告 10 个失败。

## 验证设置与隔离

原 XY rigid rotation：初始圆心 `(0.65,0.53)`、半径 `0.13`、角速度 `2`，基础网格 `n=8`，分别 uniform / 2:1 AMR，`dt=0.01`。另有 AMR `dt=0.005` 首步对照。实际路径仍用 I02–I11 原 API，记录 initial、quarter predictor、midpoint/end 的 predictor 与 actual closure。

每个半步/全步复制候选，仅将副本 M 换为独立 circle–rectangle 解析相交量；first jet、mixed closure、预测 c 相同。该反事实结果不返回给真实步进、不发布为 accepted。诊断前后核对实际 candidate identity、M 和未发布状态；后续拒绝核对 accepted identity/time/M 保持。

对真实及反事实闭合，在各自 retained graph 的 GL8 点独立求导 `delta_c * 16 xi²(1-xi)² * 16 eta²(1-eta)²`，核对 actual−predictor 的 value、gradient、Hessian；最大归一化偏差 `5.333491e-14`，门槛 `1e-10`。独立 jet 曲率公式与原 frame 曲率符合 `1e-11*(1+abs(kappa))` 门槛。这是采样核验，不是整条曲线连续最大值的区间证明。

6 次成功半步/全步的面通量，用实际 face cuts 和独立线性速度原函数重算，最大 `abs(dt*(Q_check-Q_I08))/capacity` 为 `2.220446e-18`。解析圆的中点通量直接用圆与直线相交区间计算，不使用 Hermite face cuts；仅用于误差分解。

## 首步曲率来源

下表为每个阶段自身曲线 GL8 节点上的最大曲率绝对误差。精确 M 对照保持同一预测 skeleton；其失败不是缺失值被当作 0。

| 网格、dt | 阶段/物理时间 | predictor | actual closure | 仅替换精确 M |
| --- | --- | ---: | ---: | ---: |
| uniform，0.01 | midpoint，0.005 | 1.508494e-6 | 9.960438e-2 | 3.115042e-3 |
| uniform，0.01 | end，0.01 | 1.206795e-5 | 1.106278 | 2.970493e-2 |
| AMR，0.01 | midpoint，0.005 | 1.508434e-6 | 4.403673e-1 | 6.309356e-2 |
| AMR，0.01 | end，0.01 | 1.206834e-5 | **32.942469** | **cell 72 证书区间不足，拒绝** |
| AMR，0.005 | midpoint，0.0025 | 1.885720e-7 | 1.784444e-1 | 1.618788e-2 |
| AMR，0.005 | end，0.005 | 1.508434e-6 | 4.403673e-1 | 6.309356e-2 |

AMR 大误差点位于 cell 72、`(0.5624734063,0.4363038937)`：实际曲率 `40.63478537`，解析曲率 `7.69231608`，实际梯度范数 `0.26030755`，没有接近零。首步 actual−predictor 的 bubble Hessian 分量在全部已测节点的最大变化为 `20.23780247`；独立导数核验确认其来自中心 bubble。

cell 72 为 level 1，容量 `0.00390625`，液量很小：

| 量 | 数值 |
| --- | ---: |
| predictor 几何 M | 8.2086099570e-7 |
| 实际守恒目标 M | 8.1720738213e-7 |
| 精确圆相交 M | 8.2107424740e-7 |
| 精确几何 + 单点 midpoint 更新 M | 8.1715613740e-7 |
| predictor c | 0.013116375390625 |
| actual 闭合 c | 0.379115923517020 |
| 初始 dM/dc | -1.0125222840e-8 |
| 最终 dM/dc | -9.8426651612e-9 |

体积对中心系数很不敏感，小库存差触发约 `0.366` 的系数修正。I12r3 的正则性/位移证书允许它，体积残差也通过，但二阶导数及曲率并未因此获得精度保证。证据不支持修改曲率公式、输出平滑，或继续单纯扩大系数搜索来解决这个精度问题。

## 时间取样与几何通量误差分解

逐单元核对 `(M_actual-M_exact)/capacity = (M_exact_geometry_midpoint-M_exact)/capacity + (M_actual-M_exact_geometry_midpoint)/capacity`。第一项只保留该步的 midpoint 时间取样误差；第二项含预测/闭合后的面几何、空间积分差及实际初始 M 的舍入差。

| 网格、dt、阶段 | 第一项最大绝对值 | 第二项最大绝对值 |
| --- | ---: | ---: |
| uniform，0.01，midpoint | 4.716604e-7 | 4.793186e-10 |
| uniform，0.01，end | 3.722094e-6 | 7.609270e-9 |
| AMR，0.01，midpoint | 1.714256e-6 | 1.863412e-9 |
| AMR，0.01，end | 1.349660e-5 | 3.060095e-8 |
| AMR，0.005，midpoint | 2.160341e-7 | 1.149350e-10 |
| AMR，0.005，end | 1.714256e-6 | 1.863412e-9 |

两列为各自的单元最大范数，不保证来自同一 cell。曲率峰值 cell 72 的 end 阶段两项分别为 `-1.0030361604e-6`、`+1.3118650966e-8`，时间取样误差也占主导。

6 个精确 M 反事实中 5 个闭合成功、曲率误差减小；AMR `dt=0.01` end 在 cell 72 仍以 `directed search cannot extend bracket` 拒绝。本任务没有该支路的合格修正后曲率，也没有证明精确目标不存在拓扑/位移可接受根。证书不足、可达性和连续方程误差应分开判断。

减半 dt 后单步误差下降，但终点分别为 0.005 和 0.01，不能据此报告固定时间收敛阶，或以缩短算例替代原完整 benchmark。

## 下一半步库存失败

从真实 accepted `t=0.01` 计算到 `t=0.015`，quarter geometry 时间为 `0.0125`。独立实际面内积分给出 cell 72：

| 量 | 数值 |
| --- | ---: |
| accepted M | 8.1720738213e-7 |
| face 82 流出液量 | 3.8890181087e-7 |
| face 210 流出液量 | 4.2833508439e-7 |
| face 76、190 液量 | 0 |
| 原始新 M | **-2.9513122344e-11** |
| 精确圆新 M | 5.0990819215e-10 |
| 精确旧 M + 精确几何中点通量的新 M | **-1.9777403857e-10** |

实际负库存约为容量的 `-7.555359e-9`，是 I08 原 `64*eps*capacity` 舍入预算的约 `531661` 倍，不能用舍入修复或 clipping 接受。原 `advance` 确实在 `half_flux`、cell 72 拒绝，accepted state 保持。

相对容量的总误差 `-1.3809185651e-7` 分解为：继承库存误差 `-9.8991750941e-7`、解析中点时间误差 `-1.8116665106e-7`、本步面几何/通量差 `+1.0329923040e-6`。这些项有抵消，不能当作正的“责任百分比”。去掉继承库存和几何误差后，解析中点更新本身仍为负；这是原旋转算例中的时间取样反例，补充了既有平面跨面证据。

## 构建、回归与交付

- Release benchmark 和 architecture 构建通过，无新增 compiler warning。初次测试编译引用了不存在的 time-layer 成员，已改为实际 `begin/end` 接口，初始失败日志保留。
- 显式审计最终退出 0，20.618 秒；3 配置、23 个曲率阶段、6 个面积分对照、6 个精确 M 支路及下一步原始库存重现。此 PASS 只表示审计完成。
- architecture CTest 通过。只改 1 个测试文件、0 production；1578 个受审计文件中 1577 个不变，无新增源码。69 个 S08 frozen、12 个 integration reference 和 3 份方法 PDF 哈希不变。
- 移除新增诊断函数和 CLI 入口后，原 benchmark 源码与 I12r3 完全一致。当前完整运行退出 1，255.143 秒，仍 10 个失败；原始输出 SHA256 与 I12r3 同为 `FC4483C3E9E94594CF6304995698E784FC7EB5CCFDB0EA557ED05D924289BB25`。
- 18 个 runtime 中只有 benchmark 重编译变化，其余 17 个哈希不变。本次没有重跑全部 15 项旧 CTest；旧通过记录属于 I12r3。上轮归档与交付文件在开始时核验，归档保持。

复现入口为 `build/windows/tests/Release/quantum_interface_hermite_gals_benchmark.exe --rotation-curvature-audit`，不带参数运行原完整 benchmark。环境沿用 I12r3：MSVC `19.40.33812.0`、CMake `3.28.3-msvc11`、OMP 线程及上限 2、OMP_DYNAMIC=false、BLAS/MKL 各 1；见 [verification.json](../../runs/baselines/I12r4_20261001/verification.json)。

证据：[审计日志](../../runs/baselines/I12r4_20261001/curvature_audit_final.log)、[数值汇总](../../runs/baselines/I12r4_20261001/numerical_summary.json)、[逐 cell 分解](../../runs/baselines/I12r4_20261001/cell_decomposition.csv)、[下一步 face 明细](../../runs/baselines/I12r4_20261001/next_face_contributions.csv)、[默认 benchmark](../../runs/baselines/I12r4_20261001/benchmark_final.log)、[回归比较](../../runs/baselines/I12r4_20261001/regression_comparison.json)、[源码审计](../../runs/baselines/I12r4_20261001/source_audit.json)。归档哈希/清单见 [delivery_manifest.json](../../runs/baselines/I12r4_20261001/delivery_manifest.json) 和 [archive_verification.json](../../runs/baselines/I12r4_20261001/archive_verification.json)。归档完整性 PASS 不等于方法 PASS。

## 下一步边界

继续停留在 I12。证据支持优先处理 [时间通量设计提案](I12_Phase_Flux_Time_Integration_Design_Change.md) 所指缺口；提案仍待用户明确决策，本轮未实施。新通量还必须验证与固定 skeleton 的可达性、小 cut cell 曲率精度及原失败序列，不能直接将本任务解析 M 替换当作 production 算法。

在方法决策之外，I12 仍有独立的几何相切/零迹、未解析闭合子单元、动态 regrid 与配置/输出验收缺项，可按单任务约束推进。本任务不宣称修复曲率、不进入 H01，项目低黏/无黏电喷雾总目标保持 active。
