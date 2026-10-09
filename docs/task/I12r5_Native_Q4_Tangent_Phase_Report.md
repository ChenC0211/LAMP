# I12r5 / native Q4 孤立相切的面相区间报告

**状态：PASS（本任务范围）；父任务 I12：FAIL；Gate I：未通过。**

本任务解除 native Q4 线段相测度对代数精确孤立重根的拒绝。共享面液相通量可消费这些相区间；偶重根不翻相，奇重根按原标量翻相。完整界面拓扑、时间通量和低黏/无黏电喷雾目标尚未完成。

任务范围先由 [TODO](I12r5_Native_Q4_Tangent_Phase_TODO.md) 固定，遵守 [05 路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md) 的单任务顺序。上一轮 I12r4 产生曲率/库存来源证据，归类为 progress；本轮继续 Interface 的 I12 修复，不进入 H01。

## 1. 数学事实与实现

线段上的切向导数为零，并不自动意味着二维界面梯度奇异。对于同一 retained native Q4 标量的精确限制多项式 `P(t)`，孤立重根两侧的相区间仍可定义。原先 strict isolator 把内部重根作为无法认证的 critical root 拒绝，导致实际 I08 面通量也被拒绝。

新增有限 exact-expansion 多项式认证：以 division-free pseudo remainder 求 `G=gcd(P,P')`，仅精确零余式可证明公因子。伪除法所得 `Q` 再与 `G` 重组，逐系数交叉乘积验证 `QG` 与原 `P` 精确成比例。所有归一化都是统一的正二进制幂缩放，不删除小系数、不做近似 GCD、不合并近根。

原 strict root isolation 成功时完全沿用原路径；仅 native Q4 需要额外认证时，以 `Q` 的简单根隔离得到 `P` 的全部不同零点。相区间符号仍来自原 `P` 的 exact sign，不能用 `Q` 的符号，也不能无条件翻相。细化使用同一原限制的无重因子零集，收紧 retained root enclosure；源 query 和源 snapshot 保持不可变。

含此类认证 chart 的负相变分接口明确拒绝，避免在原重根导数为零时返回伪导数。这个拒绝是整个 query 的保守合同，未增加相切处的广义形状导数。算术范围或 bounded isolation 无法证明时仍可拒绝；已有同标量 dyadic subdivision 保留，未引入另一份近似几何或物理算法。

## 2. 修改范围

| 类别 | 文件与职责 |
| --- | --- |
| production 1 | `include/Core/Geometry/Detail/PolynomialSquareFree.h`：精确公因子和重组核验，新文件 |
| production 2 | `include/Core/Geometry/Detail/LevelSetMappedSegment.h`：native Q4 接线、原标量相位及根细化 |
| production 3 | `include/Core/Geometry/Detail/LevelSetMappedSegmentSensitivity.h`：非横截 chart 变分拒绝 |
| production 4 | `include/Core/Geometry/LevelSetMappedSegment.h`：公开能力/限制注释 |
| test 1 | `tests/interface_hermite_gals_flux_contract.cpp`：真实 I08 通量、AMR/方向/近根/拒绝边界 |
| test 2 | `tests/interface_jet_transport_contract.cpp`：更新旧重根拒绝断言；真实 native query 的高次因子/缩放矩阵 |

共 4 个 production、2 个 tests。相通量空间/时间公式、M bounds、bubble/publication 门槛以及 Hydro/Electricity/App/reader 均未修改。临时直接算术测试所在 `polynomial_root_enclosure_contract.cpp` 已逐字节恢复，原 strict isolator 拒绝合同不变；最终验证通过真实 native query 消费新路径。

## 3. T0/T1 证据

环境：Windows，Release，MSVC 19.40.33812.0，CMake 3.28.3-msvc11；`OMP_NUM_THREADS=OMP_THREAD_LIMIT=2`、`OMP_DYNAMIC=false`，OpenBLAS/MKL 线程均为 1。

| 验证 | 结果与边界 |
| --- | --- |
| Release 构建 | 主程序、architecture、I01–I11、benchmark、legacy 与 jet 相关目标成功，无新增 compiler warning |
| XY/RZ 面通量 | uniform/AMR × 端点/面内相切共 8 组通过；实际通量保持纯液恒等式、内部抵消与库存约束 |
| 有理重根 | `(3x-1)^m+y`，`m=2,3,4`，两相符号、正反向、uniform/AMR 共 24 组通过；验证根库存、相测度、细化不修改源查询及变分拒绝 |
| 无理重根 | `(x²-x+1/8)²+y` 两方向通过，保留两个根，偶重根不翻相，细化后拒绝普通相变分 |
| 最高八次线段限制 | 3 个 native Q4 tensor 制造多项式 × 正负号 × `2^-400,1,2^400`，共 18 组通过；验证解析根位置、全部区间长度、相符号及细化 |
| 近重根控制 | 非二进制根的 `±2^-20` 扰动、dyadic 根的 `±2^-40` 扰动，保留 0/2 根及解析负相宽度；简单根变分仍可用 |
| 失败边界 | 正长度零迹继续拒绝；圆形完整拓扑仍在 isolated noncorner contact 处明确拒绝 |
| 最终回归集合 | 16 项各自最终结果通过：首批 15 项，加仅受后续 test 改动影响的 jet 重新编译/运行 |

高次制造多项式使用 `A(z)=z²-z+1/8`：`A(x)² A(y)²` 的对角线为八次、具有两个四重无理根；`(3x-1)³(5y-4)^4` 的对角线为七次、具有三重和四重根；`x²(x-1/2)²(y-1)²` 的对角线为六次、包含两个输入端点和内部重根。测试输入是可精确表示的 tensor Q4 数据，参考根/区间直接由这些因子给出。

原 I08 的 73 条非相切数值记录逐行相同，包括 48 组 directional case；方向子集摘要仍为 `70,88,872,1.7625900290109935e-16,1.7225311647152888e-16`。除 architecture（新增一个头文件）及 I08 扩展输出外，原 13 项完整程序输出相同，包含 I10 已有的 108 个定向扩区单元。

三个 legacy 对照 exit 0：extraction/state 全日志相同；conservative GALS 的程序数据及历史 OMP Warning #96 相同，只有 PowerShell 的调用行号/源码片段不同。未把这个既有运行时线程警告计为新增编译警告。

## 4. 保留的失败与 T2

- `flux_final.log` 是中间失败：新增 full topology 断言误以为错误仍发生于 near-critical root，实际已推进到后续 `isolated noncorner contact is unsupported`。仅修正测试预期，未放开拓扑接受条件；之后 I08 通过。
- `ctest_final.log` 名称沿用启动时计划，但该首批结果确有 1 项失败：旧 jet 断言要求 native Q4 重根必须拒绝。更新合同并增加上述解析矩阵后，`ctest_jet_verified_details.log` 为最终 jet 证据。未删改首批失败日志，也未将首批总退出码 8 改写为 0。
- `polynomial_initial.log` 的 30 组临时直接算术检查只是中间证据；该测试代码已恢复，不作为最终源码覆盖声明。

完整 I12 benchmark 实际运行 269.04 秒，退出码 1，仍为 **10 项失败**。其完整输出与 I12r4 逐字节相同，SHA256：`FC4483C3E9E94594CF6304995698E784FC7EB5CCFDB0EA557ED05D924289BB25`。因此本轮几何面相区间改进不能宣称解决旋转/vortex/跨面库存、曲率准确性或完整 I12。

## 5. 范围审计与复核入口

前序 1578 个源码/配置/算例文件中，1573 个不变、5 个修改，新增 1 个头文件；当前总计 1579。S08 的 69 个冻结文件、12 个 integration reference 以及 3 份方法 PDF 哈希保持。I12r4 归档保持 `D6668425CC16E25C8CB45AE3DB12F023BE407EB26F3458A4CFBDEA32A04B7BD6`。

原始证据位于 [I12r5 目录](../../runs/baselines/I12r5_20261001/)：

- [源码审计](../../runs/baselines/I12r5_20261001/source_audit.json)、[回归对照](../../runs/baselines/I12r5_20261001/regression_comparison.json)、[数值覆盖](../../runs/baselines/I12r5_20261001/numerical_summary.json)、[验收状态](../../runs/baselines/I12r5_20261001/verification.json)。
- [首批 CTest 原始结果](../../runs/baselines/I12r5_20261001/ctest_final_details.log)、[最终 jet 结果](../../runs/baselines/I12r5_20261001/ctest_jet_verified_details.log)、[完整 benchmark 失败](../../runs/baselines/I12r5_20261001/benchmark_final.log)。
- [归档逐项哈希复核](../../runs/baselines/I12r5_20261001/archive_verification.json)；归档完整性不代替数值验收。

复现：使用 `cmake --build --preset local` 构建上述目标；`ctest --preset local -R '^(architecture|quantum_help|interface_state_contract|interface_jet_transport_contract|interface_hermite_gals_(state|trace|vertex|edge|mixed|skeleton|amr|volume|flux|reflux|bubble|publication)_contract)$' -V` 运行 16 项集合。独立运行 `quantum_interface_hermite_gals_benchmark.exe` 得到完整 I12 结果。当前环境 PowerShell 禁止直接执行脚本时，可用 `& ([ScriptBlock]::Create((Get-Content -Raw -Encoding UTF8 runs/baselines/I12r5_20261001/audit.ps1)))` 重做证据审计，无需修改执行策略。

## 6. 剩余工作

下一项应继续 I12，单独处理 native full graph 的 isolated contact 连接与拓扑认证，不能从面测度通过推断 I07/I11 的完整发布已支持相切。正长度零迹、未解析 closed subcell cycle、动态 regrid、完整移动 benchmark 与 V2 配置/输出仍有缺口。

[时间通量设计变更提案](I12_Phase_Flux_Time_Integration_Design_Change.md) 仍未批准、未实施；依总纲 §2.1 不自行改变冻结方法。后续需要同时解决时间库存误差及剩余几何/闭合问题，保持 I12 原验收。总目标继续 active，低黏/无黏电喷雾及 Gate I/H/IH/E/C 均没有完成声明。
