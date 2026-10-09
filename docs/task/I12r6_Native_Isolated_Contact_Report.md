# I12r6 / native 孤立接触的完整界面连接报告

**状态：PASS（本任务范围）；父任务 I12：FAIL；Gate I：未通过。**

同一 native Q4 标量的正则孤立相切现在可进入完整界面 discovery、I07 体积、I10 闭合及 I11 发布。此前圆形界面的 `isolated noncorner contact is unsupported` 已解除；真正奇异零点和正长度零迹仍拒绝。完整移动基准仍有原来的 10 项失败，总项目目标没有完成。

本任务先按 [TODO](I12r6_Native_Isolated_Contact_TODO.md) 固定单一数学事实和变更预算。上一轮 I12r5 完成交付并提供新的相区间能力，归类为 progress。本轮只继续 I12，不进入 Hydro。

## 1. 接触的数学分类

边上的偶重根是二维零集的切向接触，并不意味着二维梯度为零。一个不同根在同一单元内可能有 0、1 或 2 条曲线关联；不能把“根数”等同于“曲线端点数”。

对原 scalar 的边限制保留所有不同根及包围；原标量在根两侧的 exact phase 决定是否变号。内部变号根保持一条关联。对于同号接触和边端点，限制原标量的法向导数到**完整根包围**，只有其确实非零才可认证隐式曲线的进入方向。法向朝单元内的导数符号与邻接边相位相同，表示向外接触、关联数为 0；相反时，内部根有 2 条关联，角点的单侧根有 1 条关联。角点可由另一侧提供证书，但两侧已得到的证书必须一致。

内部 chart 的根也复用 I12r5 的 exact square-free 认证。该因子只提供原零集身份，不能提供原 scalar 的法向、曲率或普通形状导数。用于积分的 transverse section 仍只接受原标量的简单根。边根细化保留同一零集和 opaque root identity，未引入第二份几何。

同一点的两条向内关联仍使用同一个 node。完整界面在该 node 连续；负相片则需区分两个局部扇区。当同一接触点同时连接负相实际边区间和界面弧时，每条负相片边界在“界面弧—实际边区间”之间转接，不能把两条液相片通过一个零测度点并成一片。跨单元合并仍仅依据共享的正长度负相区间。

已有内部 chart 分区的单元直接消费局部认证和实际边身份匹配，不再额外要求整条粗边可由一个浮点参数隔离。该安排保留实际 subedge/chart ratio 的分辨能力；图 regularity、非零跨度、完整 node incidence 及负相闭合检查均保留。

## 2. 最小修改

| 类别 | 文件与职责 |
| --- | --- |
| production 1 | `include/Core/Geometry/Detail/LevelSetMappedInterface.h`：边根相位、仅 native Q4 边上的重根认证/细化、原法向导数包围证书 |
| production 2 | `include/Core/Geometry/Detail/LevelSetNativeInterfaceChain.h`：0/1/2 条关联、完整曲线配对、相切处负相片连接 |
| production 3 | `include/Core/Geometry/LevelSetInterfacePanels.h`：公开能力和限制注释 |
| test 1 | `tests/interface_hermite_gals_flux_contract.cpp`：原完整圆形拒绝改为闭合曲线及解析曲率检查 |
| test 2 | `tests/interface_hermite_gals_publication_contract.cpp`：完整发布、体积、曲率、材料面板及失败隔离矩阵 |

最终为 **3 production + 2 tests**，未新增生产文件。没有修改 first-jet + M 状态、通量时间公式、M/bubble/publication 容限、物理网格或其它模块的数值实现。

## 3. T0/T1 验收

环境为 Windows / Release / MSVC 19.40.33812.0 / CMake 3.28.3-msvc11。`OMP_NUM_THREADS=OMP_THREAD_LIMIT=2`、`OMP_DYNAMIC=false`，OpenBLAS/MKL 线程均为 1。共 18 个可执行目标（含主程序）构建成功，无新增 compiler warning。

最终 16 项 CTest 全部通过，耗时 106.05 秒。其中 14 项完整程序输出与 I12r5 相同，I08 原有 73 条非相切数值记录相同，I11 原验证直到 transaction 检查的完整输出相同。新增测试没有替换原覆盖。

| 新增验证 | 结果 |
| --- | --- |
| 圆形 / RZ 球面 | 20 组精确相切配置，覆盖 XY/RZ、uniform/真实 AMR、正反相、水平/竖直边及网格角点 |
| 近相切两侧 | 圆心沿横向分别偏移 `±2^-12`，共 8 组，保留新出现或消失的细小相区间 |
| 非二进制二/四重接触 | `phi=(x-cx)^2+(3y-1)^m-r²`，`m=2,4`，XY/RZ × uniform/AMR × 正反相，共 16 组；接触纵坐标为 `1/3` |
| 完整消费 | I07 初始化 → I10 bubble → I11 publication，验证唯一 source/time/identity、分量/闭合/方向、负相区域、M 未被改写 |
| 解析测度 | 28 组圆形/球面逐单元及总相体积；椭圆/四次卵形的独立解析总量，包含 RZ 真体积 |
| 法向与曲率 | 圆形/球面用独立公式核验端点、内部样本和所有发布 GL8 节点；非二进制矩阵用原多项式梯度/Hessian 核验端点、内部及 416 个接近端点的样本 |
| 界面长度 / 面积 | 同一 retained chain 的既有受控 material-panel 积分，圆弧长度和球表面积均通过 `2e-9` 绝对误差检查 |
| 重根细化的失败隔离 | 2 个 `2^-60` 端点邻近采样先细化同一重根，最终仅因 binary64 chart 无法区分内部位置而明确拒绝；原源仍可用 |
| 真正不支持的对象 | `phi=(x-1/2)^2+(y-1/2)^2` 的孤立奇异零点拒绝；正长度零迹仍拒绝 |

44 组矩阵的实际最大误差：圆形/球面逐单元体积 `7.6327832942979512e-17`；总相体积 `4.4408920985006262e-16`；法向 `1.9641850382783467e-15`；物理曲率 `1.2256862191861728e-13`。这些是所测位置/配置的最大值，不是任意连续零集的全域误差界。

圆形逐单元参考值由解析 disk/rectangle 交面积分给出，不使用生产根或图。非二进制椭圆使用仿射圆面积/球体积；四次卵形的 XY 面积为 `r^(3/2) B(1/4,3/2)/3`，RZ 体积为 `8*pi*r^(5/2)/15`。曲率始终来自原多项式，特别检查四重接触处不能借用无重因子的导数。

三个 legacy 对照均退出 0。extraction 和 state 控制的完整日志相同；conservative GALS 的程序输出及历史 OMP Warning #96 相同，只排除了 PowerShell 调用行号/源码片段的差异。

## 4. 中间失败与复核记录

- [原拒绝复现](../../runs/baselines/I12r6_20261001/reproduction.log)：先把原圆形改成完整拓扑验收，在修改生产代码前重现 `isolated noncorner contact`。
- [初始 publication 检查](../../runs/baselines/I12r6_20261001/publication_initial.log)：错误地对固定单条 GL8 规则要求自适应弧长精度，出现约 `1e-6` 的弧长差。发布接口没有这项误差保证；改为消费已有受控材料面板验证同一 chain，未改求积公式或放宽容限。固定 publication GL8 的这一限制在 API 中明确记录。
- [一次链接失败](../../runs/baselines/I12r6_20261001/build_local_charts.log)：上一版 benchmark 运行期间，其 Windows 可执行文件仍被占用，重链接报 `LNK1104`。确认原进程终止后重建成功；失败日志保留。随后对最终程序重新运行完整回归、legacy 和 benchmark。
- `ctest_final*` / `benchmark_final*` 是前一实现迭代的原始记录；交付源码对应的是 **`build_verified*`、`ctest_verified*`、`benchmark_verified*` 和 `legacy_*_verified*`**。

## 5. 完整 I12 结果与范围审计

最终完整 benchmark 实际运行 **241.43 秒**，退出码 **1**，仍为 **10 项失败**。完整输出与 I12r5 逐字节相同，SHA256 为 `FC4483C3E9E94594CF6304995698E784FC7EB5CCFDB0EA557ED05D924289BB25`。这些运行时间仅作证据记录；并行验证和系统负载不同，未进行性能改进声明。

1579 个源码/配置/算例文件中，1574 个不变，仅上述 5 个文件修改。S08 的 69 个冻结文件、12 个 integration reference 及 3 份方法 PDF 哈希不变。前序 I12r5 归档哈希保持 `AFBFD6581443E455006D77DB6DF4143453A0C548642CB8898A6AB1B3F44F8238`。

证据入口：

- [最终 CTest](../../runs/baselines/I12r6_20261001/ctest_verified_details.log)、[完整 I12 失败](../../runs/baselines/I12r6_20261001/benchmark_verified.log)。
- [接触矩阵](../../runs/baselines/I12r6_20261001/isolated_contact_matrix.csv)、[非二进制矩阵](../../runs/baselines/I12r6_20261001/nonbinary_contact_matrix.csv)、[数值摘要](../../runs/baselines/I12r6_20261001/numerical_summary.json)。
- [源码审计](../../runs/baselines/I12r6_20261001/source_audit.json)、[回归比较](../../runs/baselines/I12r6_20261001/regression_comparison.json)、[验收状态](../../runs/baselines/I12r6_20261001/verification.json)、[归档完整性](../../runs/baselines/I12r6_20261001/archive_verification.json)。

复现使用 `cmake --build --preset local` 的对应目标，以及 `ctest --preset local -R '^(architecture|quantum_help|interface_state_contract|interface_jet_transport_contract|interface_hermite_gals_(state|trace|vertex|edge|mixed|skeleton|amr|volume|flux|reflux|bubble|publication)_contract)$' -V`。独立运行 `quantum_interface_hermite_gals_benchmark.exe` 得到完整 I12 结果。证据审计脚本为 [audit.ps1](../../runs/baselines/I12r6_20261001/audit.ps1)；当前 PowerShell 策略下用 `& ([ScriptBlock]::Create((Get-Content -Raw -Encoding UTF8 runs/baselines/I12r6_20261001/audit.ps1)))` 执行，不修改系统策略。

## 6. 未解决问题

该认证仍有限且是充分条件，不能保证任意 Q4 零集都能解析。正长度零迹、真正奇异零点、未解析 closed subcell cycle、动态 regrid 和实际拓扑变化仍有缺口。面重根处的普通负相变分继续拒绝；固定 publication GL8 不具有自适应弧长精度保证。

下一任务仍属于 I12，须处理剩余时间通量、库存/闭合与几何失败，保持原完整 benchmark 验收。[时间通量设计变更提案](I12_Phase_Flux_Time_Integration_Design_Change.md) 仍未批准、未实施，不自行改变冻结方法。V2 配置/输出、Hydro、Electricity 及低黏/无黏电喷雾联合验收尚未完成。总目标保持 active，不进入 H01。
