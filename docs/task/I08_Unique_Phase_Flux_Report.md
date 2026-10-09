# I08 / 唯一相面通量与有界约束报告

**状态：PASS；2026-10-01。** 路线图 I08 的内部面抵消和有界通量约束通过；下一项 I09。完整运动几何、低黏/无黏与电喷雾目标未完成。

## 实现

按 [TODO](I08_Unique_Phase_Flux_TODO.md) 修改 **4 production + 2 tests**：新 `HermiteGalsPhaseFlux.h/.cpp`，Core `UnitPolynomialRoots.h` / `LevelSetHermiteGeometry.cpp`，新 flux contract / tests CMake。无 Hydro、Apps、reader 或冻结专有路径修改。

本次 midpoint scalar 的所有 actual native half-face cuts 来自同一 geometry。输入 normal velocity 为拥有独立存储的局部 polynomial，显式绑定 mesh、geometry identity、时间层和 boundary kind。XY/RZ metric 统一加入；Core 新 range port 用 exact expansion 平移/缩放后复用原 positive integrator，保留小区间的 offset/width。

每个实际 face 只返回一条 `Q+、Q-、Q_l+、Q_l-` 及其净值。`dt*Q_l` 生成一次，两个内部 owner 一加一减。gross CFL 包含正反流，即使净值为零。old M 独立参与更新，midpoint M 不覆盖它；结果保留 old/midpoint identities 和原 scalar 的只读所有权。

更新后严格要求 `0<=M<=capacity`。有限幅值越界拒绝；只有 `64*eps(double)*capacity` 内的舍入 defect 可缩减内部共享液相流出/流入。boundary 和 total flux 不改；重组装 M 后再次检查严格 bounds。没有 M clipping、upwind/FCT fallback 或隐藏时间推进。

## 数值证据

8 组解析 midpoint 输运覆盖 XY/RZ、uniform/2:1 AMR、正反速度。每组额外验证全液/全气通量恒等式及零速度 M 严格不变。XY 为斜平面沿 x 平移，RZ 为横向平面沿 z 平移；逐 face 与独立解析 primitive、逐 cell 与解析新体积比较。

| 指标 | 实测 |
| --- | --- |
| 解析输运最大 cell M 误差 | 1.3877787808e-17 |
| 解析输运最大 face rate 误差 | 2.7755575616e-17 |
| 最大全局 liquid ledger 绝对 closure | 5.5511151231e-17 |
| XY 半面内 polynomial counterflow 最大误差 | 6.2612675461e-18 |
| RZ 半面内 polynomial counterflow 最大误差 | 4.3097036356e-18 |
| 微修复 | 1 次 sweep，4 个内部面 |
| 微修复最大相对 face correction | 1.3877787808e-17 |
| 微修复后全局 closure | 0 |

每条内部 face 的正反 owner contribution 逐条检查，`(-liquidVolume)+liquidVolume` 严格为零。全局求和另按 `1e-12` gate 检查；没有用全局小残差替代逐面一致性。

零净通量 counterflow fixture 有 20 个双向面，gross outflow CFL=0.125，M 不变。增大 dt 后，即使 signed net 仍为零也明确拒绝 gross CFL>1。多根 phase trace 在 XY/RZ 两种 measure 下均保留两个相界根；正反 liquid directions 同时非零的 faces 各 7 个。

Core 直接测试 `offset=0.5, width=2^-60`，常数积分为 `2.6020852139652106e-18`，位于 root 附近的线性正部积分为 `3.76158192263132e-37`；两种 fromEnd 方向均不丢失，返回误差估计为 0。另验证 degree-10、正部 crossing、旧 signed-range 和旧 positive entry；小区间有独立相对误差检查。

## Face-cut 收敛

沿流向均匀、因此库存不变的测试中，曲面位置仅取决于 x/r：XY 使用 `exp(0.7*x)-exp(0.7/3)`，RZ 使用 `exp(r*r)-exp(1/9)`；速度沿 y/z。它们的精确截断位置是 1/3。输入 mixed 来自 I06；这里刻意隔离 face cut 表示误差，库存作为独立 operator fixture 提供，不冒充 I11 的体积闭合结果。

| measure | N=8 face flux error | N=16 error / 阶 | N=32 error / 阶 |
| --- | --- | --- | --- |
| XY | 4.9653209322e-10 | 1.5800916131e-11 / 4.97381 | 4.8933079810e-13 / 5.01305 |
| RZ | 2.8079973013e-7 | 9.8443710961e-9 / 4.83410 | 2.9011409941e-10 / 5.08461 |

误差下降，最后阶满足原定 `>=3`。没有用这个局部约五阶结果声称完整时间推进或曲率五阶；I12 的全步 benchmark 尚未开始。

## 失败隔离和支持范围

最终 18 个拒绝检查覆盖：非法 mapped range、gross CFL、宏观库存越界、boundary-only 微修复无内部容量、缺/重复 face、NaN、错误 boundary kind、零净值但非恒零的不可渗透速度、非 midpoint、foreign geometry、mesh/time 错配、非法容限、divergence、old M 超容量，以及 XY/RZ 面内精确相切根。

外部系数修改不影响 velocity；结果在局部输入 State 析构后仍保留同一 scalar。失败不改变 old/midpoint identities 或旧 M。成功修复只改变内部 liquid flux，boundary 不改；无内部容量时明确拒绝。

**现有 Core 的面内精确重根/相切根限制仍然存在。** 补充测试 `phi=y-(x-r)^2` / RZ `phi=z-(r²-r0²)^2` 时，面内 root 被现有 isolator 以 near-critical/multiple root 拒绝；actual half-face 端点 root 在 XY/RZ 下通过，且整液相通量严格等于 total flux。未修改 root policy、删除零点或虚构 phase fraction。positive-length coincident zero trace 同样仍拒绝。

这些不是完整运动几何的完成证据；面内相切限制必须在 I11/I12 几何与运动验收中继续处理。当前 I08 在既有可解析 face partition 合同内实现唯一通量，并对不支持的分段显式失败。

## 构建、回归与失败记录

- Release 主程序和全部目标构建通过，完整构建 241.234 s；Core 仍为 `/fp:strict /GL-`。
- 初次完整 CTest 11/12 通过，新增探索性面内相切检查失败；修正该检查对既有 Core 支持范围的错误预期后，I08 单独 CTest 1/1 通过（0.88 s）。生产代码在这次测试修正前后没有变化。合计 12 个目标均有通过证据，不把初次失败日志覆盖成绿色。
- 原 native geometry、I06a skeleton、I07 volume 完整输出与 I07 归档逐字相同；I03/I04 23 行数值相同。旧 GALS 完整输出哈希相同，parametric state control 退出 0。
- Core signed-range 的 endpoint arithmetic / integral arithmetic 提取前后逐段相同，旧 positive integrator 未改变；没有把旧积分替换成近似重拟合。
- 1565 个既有源码/输入/工具文件中 1562 个不变；69/69 S08 专有文件与 3/3 PDF 一致。
- 构建出现的 2 个 C4244 warning location 均在冻结的 `HydroCompatibleMomentSplineTransport.cpp:115`；10 行完整诊断与 S08 build.log 相同。无新增 warning，未越权改冻结 Hydro 来清除历史警告。

早期日志保留：首次编译缺少 Snapshot 完整头，补 include 后通过；首次数值测试对 native face 0 的内外属性作了错误假设，改为按实际 owner 选择内部面；额外相切测试明确区分端点支持与面内拒绝。预先写定的误差、收敛阶、严格 bounds 与守恒门槛均未放宽。

## 交付和下一项

[API](../api/I08_Hermite_Gals_Phase_Flux_API.md) 说明输入单位、时间层、gross bounds、修复预算与错误语义。`runs/baselines/I08_20261001/` 保存所有原始日志、数值 CSV、Core 提取/strict flags/冻结/旧输出审计、源码和 runtime 哈希及逐条验证的 ZIP；manifest 给出归档 SHA-256，Git 为 unavailable。

下一项 **I09：2:1 reflux，coarse flux=sum(fine subflux)**。本轮只对实际 native subfaces 做唯一组装，未实现该聚合。I10 bubble、I11 authoritative geometry、I12 完整 Interface、Hydro/IH/E/C gate 继续未完成；S01 `hydro_moment_xy` 已知首步失败保持原分类。
