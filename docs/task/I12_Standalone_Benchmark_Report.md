# I12 / standalone benchmark 交付与失败报告

**状态：FAIL。日期：2026-10-01。** I12 和 Gate I 未通过；不进入 H01。总体 goal 保持 active。已提交通量时间离散设计变更提案，未获批准前保持生产公式不变。

## 本轮完成的工作

按照 [TODO](I12_Standalone_Benchmark_TODO.md) 建立可执行的完整一步验证驱动：从同一旧 accepted snapshot 回溯 vertex/edge first jet 与 center prediction，形成私有 quarter-time geometry，使用 I08 更新半步 M，再经 I06/I10/I11 形成 midpoint；同源 midpoint 通量与 I09 reflux 给出全步 M，最后与全步 skeleton 闭合发布。只有全部成功后才替换 accepted handle。

预测、FV、closure、publication 均调用现有公开模块 API，没有复制这些算式，也没有逐步用解析 geometry/M 覆盖数值解。解析 first jet 仅用于初始化，以及真实外边界上有 boundary-hit 的已知边界数据。内部 master points 全部走 I03/I04；hanging 仍由 master constraints 导出。失败前后检查原 accepted identity 与 M 精确不变。

新增 `tests/interface_hermite_gals_benchmark.cpp` 并在 `tests/CMakeLists.txt` 注册。**0 production + 2 tests**；生产公式、验收容限和冻结方法均未修改。新 benchmark 失败时返回非零，未设置 `WILL_FAIL`、禁用测试或略过失败用例。

## 数值结果与证据解释

- XY/RZ × uniform/AMR 的四组静止界面完成 4 步，`dt=0.005`，终时 0.02；保持库存及 first jet，无系统 drift。
- 静态非多项式 `phi=exp(r^2-R^2)-1` 的 8/16/32 网格 curvature study 实际执行；gradient/curvature 分别报告末档阶数。它验证空间几何，不能替代移动界面的时间收敛。
- 平移、旋转、vortex、AMR crossing 均未全部完成目标时间。完整日志保留失败阶段、原始异常和最后 accepted time；没有把拒绝候选的时间计入进度。
- 平移另外执行 XY/RZ 的 h/dt 三档细化矩阵，目标终时统一为 0.08；无法在共同终时比较的结果不计算伪收敛阶。
- 原来 RZ R=0.19 不触及 `[0.25,0.75]^2` 的细化区。复核后将最终 RZ 组统一设为 R=0.27，并增加粗/细 cut-cell 同时非空的检查；早期数据保留但不作为 RZ AMR crossing 证明。
- 初次 vortex 的 `dt=0.005` 超出 I02 半个最小单元的 trace CFL。最终采用 `dt=0.001`、相同速度/目标终时，不改变数值容限；仍出现 phase gross flux 界限拒绝或 bubble 目标不可达。前者的误差大小/根因尚未定位，不能假定它只是舍入问题。

最终数据以 `runs/baselines/I12_20261001/benchmark_verified.log`、`refinement_final.log`、`numerical_summary.json` 为准。`benchmark_initial.log`、`benchmark_final.log` 与旧 `refinement.log` 保留排查过程，后缀 final 不自动代表最后版本，以 verification/manifest 记录为准。

误差定义：position 为 surface quadrature 点上的 `abs(phi_exact)/abs(grad_exact)`；gradient 为同点 relative gradient norm；curvature 为同源发布曲率与解析 Hessian 曲率的最大绝对差；global volume 为库存总量与解析圆/球体积的绝对差。另独立解析积分每个矩形内的圆/球部分，报告 `max(abs(M-M_exact)/capacity)`；该 reference 的全域求和也须通过解析总量检查。vortex 没有独立逐 cell reference，字段明确为 NaN，不声称此项完成。

解析 differential-rotation vortex 的速度为 `omega(t)*r^2*(-y,x)`（坐标相对 `(0.5,0.5)`），`omega=40*cos(pi*t/0.4)`；其连续精确逆映射为按 `-r^2*integral(omega)` 旋转。测试自带二阶解析自动微分，仅用于独立 gradient/Hessian reference，不调用数值 tracer 或 Snapshot 的导数作为真解。

| 最终证据 | 数值 / 状态 |
| --- | --- |
| 静止 4 步最大 first-jet/M drift | `2.2204460492503131e-16` |
| 静止最大逐 cell 解析 M 误差 / capacity | `9.9920072216264089e-15` |
| 静态 curvature study | 4 组 × 3 网格，共 12 行；末档 gradient 阶数最低 3.958，curvature 阶数 2.970–3.413 |
| AMR 初始 cut 覆盖（coarse/fine） | XY 静止/平移 2/21，RZ 静止/平移 8/4，XY rotation/vortex 2/13 |
| 移动 benchmark | 8 组均未完成目标时间；最终 AMR vortex 接受至 0.034，目标 0.4 |
| 平移 h/dt study | 12 次运行，10 个不同配置，均未达到共同终时 0.08 |
| 前序回归 | 15/15 PASS，完整输出逐字等于 I11；17 个前序 runtime 文件哈希不变 |
| 新 benchmark | 退出 1，8 个移动失败 + 2 个解析跨面失败；CTest 返回 8 |
| 运行时间 | 16 项完整 CTest 194.43 s；修正 AMR fixture/逐 cell reference 后最终专项 141.34 s，最终细化矩阵 65.45 s |

完整 CTest 后仅修改 benchmark fixture/reference，再执行该目标和细化矩阵；未改动或重编前序 15 项的数值代码。前序 runtime 哈希相同，故不为汇集日志重复运行已通过的测试。

## 决定下一步的解析反例

在 4×4 网格中以 `u=(1,0)` 平移平面界面，`dt=0.02`，CFL 仅 0.08。穿越 `x=0.5` 的时刻分别取 dt 的 1/4 与 3/4。

| 穿越时刻 | 原始 midpoint 过面液量 | 解析过面液量 | 上游原始库存 / 容量 | 结果 |
| --- | --- | --- | --- | --- |
| dt/4 | 0.005 | 0.00375 | 0.06125 / 0.0625 | I08 有界但 I10 拒绝：纯液单元库存缺口会要求内部造气泡 |
| 3dt/4 | 0 | 0.00125 | 0.06375 / 0.0625 | I08 拒绝：2% 容量超界，远超舍入预算 |

该反例使用解析中点/末端 skeleton，排除了第一步预测器误差。bubble 在 face 上值和梯度为零，不能把整面中点的全液/全气指示改成时间平均的部分液相。增大 bubble guard、clipping M 或只选步端穿越都不能作为修复。

现有 strict bound/topology guard 正确拒绝了不相容输入；问题是单点时间取样不能普遍产生同时有界、且与末端 skeleton 拓扑相容的库存。详见 [设计变更说明](I12_Phase_Flux_Time_Integration_Design_Change.md)，其中提出事件分段的时间平均相通量及保持的离散约束、风险和后续 I12 子任务。

## 验收结论

构建与前序 I01–I11 回归通过，新 standalone benchmark 明确失败。前序 operator PASS 不扩张为完整算法 PASS。源审计确认前置 1577 文件中只有 tests/CMakeLists.txt 修改，另新增一个 benchmark 文件；69 个 S08 专用源、3 PDF 不变，无 production 改动或新增 compiler warning。

完整 Gate I 还缺少：动态 refinement 的 exact polynomial transfer、coarsening 的严格 M 汇总与新 authority 事务、moving-interface 时间/空间收敛、完整强变形长时新旧 A/B，以及 V2 reader/配置追踪/accepted VTK diagnostics 接线。已有 I06 仅证明 hanging/master restriction，不证明动态 regrid。搜索与当前 API/前序报告均无足够证据证明这些要求已实现。

Core 面内精确相切/重根、正长度零迹、未解析 closed subcell cycle；I10 充分 guard 的可达性限制；S01 `hydro_moment_xy` 首步失败与历史冻结 Hydro C4244 均仍保留。后续低黏/无黏 Hydro/EHD、电喷雾射流和破碎目标尚未达到。

## 交付与下一步

交付包含源码、可执行、初次失败、最终结果、refinement matrix、逐项 Gate I 审计和全部哈希。`archive_verification.json` 的 PASS 只代表归档字节/哈希核验通过；任务及 `verification.json` 的数值状态仍为 FAIL，不能混淆。

下一步仍在 I12 内：用户确认设计变更后，按受控文件预算推进同源时空 phase flux 子任务；若不批准，先比较其它满足冻结关系的方案。本轮没有实施提案中的时间离散变更。
