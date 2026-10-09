# I12r1：空间 gross phase-flux bounds 修复

**状态：PASS（本项空间积分一致性修复）；父任务 I12 仍为 FAIL。** 总体目标保持 active，Gate I 未通过，不进入 H01。执行范围见 [TODO](I12r1_Gross_Phase_Flux_Bounds_TODO.md)。上一轮属于 progress；本轮接续完成代码、测试和证据交付。

## 结果与数学范围

修复了 I12 均匀网格 vortex 首个 half-step 中，独立液相子集积分略大于同向 full integral 导致的拒绝。仍使用同一 midpoint scalar、同一速度 polynomial、同一 XY/RZ metric 和原有时间取样公式。总通量保持由速度独立决定，最终液相通量在 actual face 上唯一生成并共享给两侧 cell。

修复后，原 vortex 均匀网格与 AMR 的完整第一步均通过。完整 I12 重跑中，均匀网格 vortex 从 `t=0` 首步拒绝推进到 `t=0.031`（31 个 accepted steps），之后在 midpoint bubble volume bracket 拒绝；AMR 仍在 `t=0.034`（34 步）拒绝。两者要求的终点均为 `t=0.4`，因此不能据此声明 vortex 或 I12 完成。

## 原始缺陷与误差尺度

复现使用原 I12 vortex：XY、`n=8`、圆心 `(0.65,0.53)`、半径 `0.13`、速度系数 `40*cos(pi*t/0.4)`、`dt=0.001`，quarter-time scalar 为 `t=0.00025`。失败面为 native face 33，其 half 1 为 `x=0.75, y∈[0.4375,0.5]`。局部参数 `s∈[0,1]` 的速度系数依次为：

```text
0.16601530497890299, -0.18554651732936217,
0.029296818525688762, -0.0097656061752295874
```

metric factor 为 `0.0625`。gas range 为 begin-anchored `offset=0,length=0.15288541563591540`；liquid range 为 end-anchored `offset=0,length=0.84711458436408460`。

| 该 half-face 的反向量 | 原始积分 | primitive 返回的 absolute error |
| --- | ---: | ---: |
| full | `2.9815559743351415e-19` | `2.9815559743351415e-19` |
| liquid | `3.0496459334343837e-19` | `3.0496459334343837e-19` |
| gas | `0` | `0` |

子集超量为 `6.8089959099242166e-21`，小于 full/liquid 合计误差估计 `6.0312019077695252e-19`。以该半面正反向总 gross rate 归一化为 `1.3522279277000121e-18`，低于未改变的 `64*eps(double)=1.4210854715202004e-14`。

根因是 positive-part primitive 对不同区间独立给出带误差估计的结果，并不保证这些估计满足域的包含单调性。不能以接近零的反向估计作为唯一误差尺度，也不能把这些估计当成精确积分。Core primitive 未修改；前后原始 probe 的系数、区间和独立积分相同。完整数据见 [修复前 probe](../../runs/baselines/I12r1_20261001/gross_probe_before.log) 和 [最终复现](../../runs/baselines/I12r1_20261001/gross_contract_final.log)。

## 实现与诊断

仅当某一 half-face/direction 的液相估计 `L>Q` 时，另积分同 polynomial 的气相补集 `G`，取 `L_consistent=Q*L/(L+G)`；`G=0` 时直接复用固定总量 `Q`。总量不随相几何改变。

初始超量与完整调整都必须同时受到既有积分误差估计及 `64*eps` half-face gross 预算限制；`liquid 原误差 + 调整` 必须不超过既有 half-face integration tolerance。每个 face 的单方向调整乘 dt 还受 `64*eps*min(adjacent capacities)` 限制。任何不足均拒绝，未增大预算、裁剪 M、改时间公式或引入旧算法后备。

`integratedLiquidPositiveRate/NegativeRate` 保存独立积分原值；原 `rawLiquid...` 字段表示空间一致性处理后、库存微修复前的值；最终 `liquid...` 字段仍是共享库存更新依据。新增报告记录总调整数、零气相恒等式恢复数和最大相对调整，调整量计入积分误差估计。空间一致性处理与原有 M-level flux repair 的诊断分开。详细合同见 [phase flux API](../api/I08_Hermite_Gals_Phase_Flux_API.md)。

最初只处理零气相方向的方案，在相界面紧贴速度根且气相非零的 AMR 配置失败，保留于 `flux_matrix_initial.log`。最终方案覆盖该情形。后续首次全量 CTest 又因新增测试漏统计互补相侧的调整而失败；补齐两侧 ledger、总通量不变性和预算断言后 I08 通过。失败日志均保留，未删除或重分类。

## 验收证据

证据根目录：[I12r1_20261001](../../runs/baselines/I12r1_20261001/)。Release 使用 MSVC 19.40.33812、CMake 3.28.3-msvc11；OpenMP 2 threads，动态线程关闭，BLAS/MKL 1 thread。无 Git 信息，记录为 unavailable。

| 检查 | 结果与边界 |
| --- | --- |
| 主程序及相关 contracts 构建 | PASS，无新增 compiler warning |
| 前序 15 项回归 | 最终全部 PASS：`ctest_details.log` 的其他 14 项通过；测试统计修正后 `ctest_flux_final_details.log` 的 I08 通过。首次全量耗时 53.60 s，I08 最终复跑 3.82 s |
| 新空间积分矩阵 | 48 配置 = XY/RZ × uniform/AMR × 两流向 × 六切分位置，每项含液相/互补相/纯气相总量对照 |
| 两类调整覆盖 | 共 88 次 half-face/direction 调整，其中 70 次 gas 值及其误差为零，18 次 gas 非零；每侧均检查共享 ledger、原始值/调整与 capacity 预算 |
| 最大调整/half-face gross | `1.7625900290109935e-16`，低于原有 `64*eps` 预算 |
| 最大互补和误差/face gross | `1.7225311647152888e-16`；full 正反向和净通量对不同 phase geometry 完全相同 |
| 既有 I08 验证 | pure phase、counterflow、微修复、多根、tiny ranges、符号、收敛、失败和生命周期测试仍通过 |
| 原 vortex gross 复现 | `--gross-flux-contract` exit 0，原始超量保留，实际 half flux PASS、ledger closure 0 |
| 原 vortex 首步 | `--vortex-smoke` exit 0：uniform/AMR 均达到 `t=0.001`；AMR 确实含 coarse/fine cut cells 2/13 |
| 时间离散不变控制 | `--crossing-probe` exit 1；7 行 CROSSING 原始数据/拒绝与 I12 完全相同，保留 q=1/4 和 3/4 的时间取样失败 |
| 完整 I12 | 重跑 209.81 s、exit 1：4 行静止、12 行静态曲率与 I12 逐字相同；8 个运动算例与 2 个跨面反例仍失败 |
| legacy controls | GALS extraction 与 parametric state 控制 exit 0，输出 SHA256 均与 I11 相同 |

完整 I12 的均匀 vortex 在最后 accepted time `t=0.031` 的位置/梯度/曲率误差分别为 `6.9176588684894183e-05`、`1.1307707837784004e-02`、`2.4727450666732764e-01`；全局体积误差为 0，独立逐 cell vortex 体积参考仍未建立，字段为 NaN。全局守恒并不证明几何精度或完整时间推进通过。

完整运动失败的前后阶段、accepted time 和原因逐项保存在 [moving_comparison.csv](../../runs/baselines/I12r1_20261001/moving_comparison.csv)。除均匀 vortex 外，其余 accepted error 行相同；后续对象 identity 因增加 accepted steps 而变化。旧 14 项完整测试输出及 I08 既有输出前缀/44 行数值指标均相同，见 [regression_comparison.json](../../runs/baselines/I12r1_20261001/regression_comparison.json)。

可重复关键命令（Release）：

```text
ctest --test-dir build/windows -C Release -R "^(architecture|quantum_help|interface_state_contract|interface_hermite_gals_(state|trace|vertex|edge|mixed|skeleton|amr|volume|flux|reflux|bubble|publication)_contract)$" --output-on-failure
build/windows/tests/Release/quantum_interface_hermite_gals_benchmark.exe --gross-flux-contract
build/windows/tests/Release/quantum_interface_hermite_gals_benchmark.exe --vortex-smoke
build/windows/tests/Release/quantum_interface_hermite_gals_benchmark.exe --crossing-probe
build/windows/tests/Release/quantum_interface_hermite_gals_benchmark.exe
```

最后两条当前仍返回 1，不应通过 `WILL_FAIL` 或忽略退出码将父 gate 变绿。本轮没有复跑未改变设置的 refinement sweep；先前 I12 的 refinement FAIL 仍有效，不能声明时间/网格收敛。

## 交付与冻结检查

- Production 2 文件：[header](../../include/Modules/Interface/HermiteGalsPhaseFlux.h)、[implementation](../../src/Modules/Interface/HermiteGalsPhaseFlux.cpp)。
- Tests 2 文件：[flux contract](../../tests/interface_hermite_gals_flux_contract.cpp)、[benchmark](../../tests/interface_hermite_gals_benchmark.cpp)。没有新增 production/test 文件。
- 文档：本报告、任务 TODO、phase flux API、docs README；时间积分提案仍保持未批准状态。
- 1578 项源码/配置/case 基线中 1574 项未变，改动严格为上述 4 文件。S08 三条路径 69 个 frozen 文件、12 个 integration/reference 文件及 3 份方法 PDF 哈希均匹配；Core、Hydro、Electricity、App、reader 与 CMake 未改。

归档包含修改文件、运行程序、全部成功/失败日志、源码审计和数值对比；最终文件清单、archive SHA256 和解压流逐项验证结果见 [delivery_manifest.json](../../runs/baselines/I12r1_20261001/delivery_manifest.json) 与 [archive_verification.json](../../runs/baselines/I12r1_20261001/archive_verification.json)。归档完整性 PASS 不等于父任务 I12 PASS。

## 未完成事项

下一工作仍属于 I12。任意步内跨面的时间通量误差未解决；冻结方法的时间积分变更见 [独立提案](I12_Phase_Flux_Time_Integration_Design_Change.md)，需要明确决策后实施（[总纲 §2.1](../00_LAMP_Development_Charter.md)：“不得由 Codex 在实现过程中自行改变方法。”）。本项不改变该决策状态。

bubble guard 的可达性、Core 内部精确相切/多重根、正长度零迹、未解析 closed subcell cycle、动态 regrid、完整移动界面收敛与强变形 A/B 仍未通过。V2 reader/accepted VTK 输出、Hydro/Electricity/EHD 及低黏/无黏目标亦未完成；ParaView 安装本身不构成这些输出或数值验收。
