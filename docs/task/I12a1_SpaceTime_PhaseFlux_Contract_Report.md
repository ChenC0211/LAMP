# I12a1 / Space–Time Phase Flux 区间契约与来源身份报告

TASK: I12a1  
STATUS: **PASS（仅区间/来源/积分体积契约）**

本任务按用户 [任务书](I12a1_SpaceTime_PhaseFlux_Contract_TASK.md)、后续“中心思想优先、冲突细节遵循现有程序路线”的补充和先行 [TODO](I12a1_SpaceTime_PhaseFlux_Contract_TODO.md) 完成。旧 I12 仍 FAIL，未声明 Gate I 通过。

## MODIFIED FILES

| 类别 | 文件 | 内容 |
| --- | --- | --- |
| production 1 | `include/Modules/Interface/HermiteGalsTemporalPhaseFlux.h` | 新的区间/来源/积分体积只读契约 |
| production 2 | `src/Modules/Interface/HermiteGalsTemporalPhaseFlux.cpp` | 来源绑定、结构校验、派生 rate 和受检区间组合 |
| test 1 | `tests/interface_temporal_flux_interval_contract.cpp` | 所有新契约的正向、错误来源、非法输入、生命周期及失败隔离检查 |
| test 2 | `tests/CMakeLists.txt` | 注册新 contract target/test |
| 文档 | `docs/I12a1_Hermite_Gals_Temporal_Phase_Flux_API.md` | API、单位、来源身份、使用前提和边界 |
| 文档 | 本报告、对应 TODO、`docs/README.md`、既有时间通量设计提案 | 任务记录、文档入口及本轮授权范围 |
| 证据 | `runs/baselines/I12a1_20261001/` | 原始编译/测试/跨面日志、审计脚本及身份清单 |

范围为 **2 个 production + 2 个测试相关文件（含 CMake）**。没有修改用户任务书。新增 V2 文件沿用 `HermiteGals*` 前缀；已存在的 `InterfaceTemporalPhaseFlux.*` 承载旧数值积分，因此保持原样，避免把新契约加入旧算法类中。

## CONTRACT IMPLEMENTED

- **authoritative integrated quantities:** 四个有限、非负的 total/liquid positive/negative integrated volume 是唯一保存的面传输数值。liquid 逐方向不得超过 total；零测度面只能传递零体积；`-0` 接受，负微量不裁剪、不引入新容差。实际 native face 及相邻 cell 来自现有 `NativeLeafDomain`，只采用 canonical +x/+y 方向。
- **interval/source identities:** 区间保留旧 `HermiteGalsState` 原 handle 和原 snapshot，以及真实 `[begin,end]`、现有 measure 和精确 `NativeMeshIdentity`。速度来源持有实际两端 `HermiteGalsVelocityView`；characteristic 来源持有现行 RK2 的实际 Domain/coverage 及 mesh identity。来源复制共享不可变 owner；独立构造不视为同一来源，几何/时间相同不足以认定相同速度。区间身份由旧 geometry、两个来源 owner 和端点共同确定，half/full 不混淆。
- **derived average-rate API:** signed volume 按正负分量作差，`averageTotalFlux()` / `averageLiquidFlux()` 现场除以所属区间 dt，不存储 rate，不接受外部另传 dt。派生值溢出明确失败，原积分量保持。
- **interval-composition rule:** `left.appended(mesh,right)` 仅允许同源、同面、同方向且精确首尾相接的两个区间；四个体积分量相加，返回 `[a,b]`。不同旧 geometry、velocity/characteristic owner、无效 mesh、不同 face、倒序、重叠或间隙都拒绝。不执行求积、reflux 或 M 更新，任一失败不改变输入。

来源身份沿用现有共享不可变 owner 路线（如 `InterfaceConservativePhase::Revision` 和 `InterfaceNormalFluxView::scaled`），没有重新发明数值 ID、全局计数器、地址哈希或任意标签。既有 velocity view 只含几何身份的问题通过持有**真实 view 对象的共享 owner**解决，不修改 `HermiteGalsTrace` 或 Core 公式。

mesh revision 的含义沿用现有精确 geometry/connectivity 内容校验：同内容 clone 可用，实际坐标或拓扑改变后拒绝；完全恢复原内容后又相容，不额外要求单调递增 revision。provider 对任意 callback capture 的有效性/不可变性责任与现有 tracer 相同，不声称能验证外部引用的不可变性。App 负责选择 accepted Initial/End handle，本契约不替 App 接受状态。

## TESTS RUN

环境：Windows / Release / MSVC 19.40.33812.0，CMake 3.28.3-msvc11。OpenMP 线程数/上限 2，动态线程关闭，OpenBLAS/MKL 线程数 1。

1. **command:** `cmake --build --preset local --target <20 个已有/新增目标> --parallel 1 -- /nr:false`。
   **result:** 退出 0，25.30 秒，含主程序、新 interval contract、旧 time-flux contract 和既有 Interface targets；无新增 compiler warning。完整目标数组见 [build_verified_result.json](../../runs/baselines/I12a1_20261001/build_verified_result.json)。不是声称所有仓库测试目标都已构建。
2. **command:** `ctest --preset local -R '^(architecture|quantum_help|interface_state_contract|interface_jet_transport_contract|interface_time_flux_contract|interface_temporal_flux_interval_contract|interface_hermite_gals_(state|trace|vertex|edge|mixed|skeleton|amr|volume|flux|reflux|bubble|publication)_contract)$' -V`。
   **result:** 18/18 PASS，108.58 秒。新增测试涵盖 XY/RZ × 1×1/2×2 网格、共享内部面的实际 owner、169 次预期拒绝、非等长区间组合、旧几何与 callback 生命周期、restart 子区间，以及对保留源调用原 tracer 的可用性。构造和组合本身没有调用 velocity callback。
3. **command:** `quantum_interface_hermite_gals_benchmark.exe --crossing-probe`。
   **result:** 原始退出码 1；q=1/4 和 q=3/4 两次失败保留，日志换行规范化后与 I12r7 完全相同。没有将这次运行标为数值 PASS，也没有 `WILL_FAIL`。
4. **command:** `& ([ScriptBlock]::Create((Get-Content -Raw -Encoding UTF8 runs/baselines/I12a1_20261001/audit.ps1)))`。
   **result:** 范围、冻结文件、测试输出和旧证据身份审计通过。无需修改系统 execution policy。

具体 T1 覆盖：

| 任务书要求 | 证据 |
| --- | --- |
| T1.1 正常 interval 与不复制几何 | 原 snapshot 地址及 geometry ID 保持，两个来源的只读对象共享，dt 正确 |
| T1.2 非法 interval | 等长、倒序、NaN/Inf、有限端点差溢出、越出 velocity 时间范围拒绝；M/geometry 不变 |
| T1.3 volume/rate | `5,2,3,1` 与 dt=0.5 给出 total/liquid signed volume `3,2`，rate `6,4`；修改调用方输入无效 |
| T1.4 gross | 四分量逐个测试负数、最小负非正规数、Inf、NaN；正反液量越过 total 拒绝；-0 接受 |
| T1.5 组合 | 四分量分别相加，合并区间等于同源 full；等长及非等长分段都通过 |
| T1.6 错误组合 | 不同旧 geometry、同 geometry/time 的不同 velocity、不同 characteristic、不同网格/陈旧坐标、不同面/方向、重叠/倒序/间隙拒绝 |
| T1.7 half/full | 同旧 authority 的 `[0,.5]` 与 `[0,1]` 不同，绝不以一个虚假 midpoint identity 代替 |
| 生命周期/强失败隔离 | 输入 state/mesh/domain/view 释放后区间仍能读旧 scalar；query capture 生命周期保留；moved-from、加法/派生 rate 溢出拒绝 |

## REGRESSION

- **old midpoint numerical path changed? NO.** 15 项既有测试的完整程序输出与 I12r6 相同，其中包括 I08 phase-flux；第 16 项 architecture 正常反映新增文件：translation units 378→379、source files 894→896、includes 3068→3074，仍为 12 个 build-graph targets、5 条 exact legacy edges，依赖约束通过。
- **dictionary changed? NO.** 字典、reader、case 和配置源码不变。
- **inventory/bubble/reflux changed? NO.** 没有改这些源码或消费类型，原 crossing 失败仍在。
- 1579 个原源码/配置/算例文件中仅修改测试 CMake；新增 3 个文件为 2 production + 1 test。69 个冻结文件、12 个 integration reference 和 3 份方法 PDF 均不变。
- I12r7 旧归档 SHA256 保持 `074BAC7A46EE42885E23C5F53D095485651DBEA00858B92B695C46B9858F125E`。本任务没有重跑完整 I12 benchmark，不借数据契约测试替代其验收。

首次审计错误地要求 architecture 文件计数输出也逐字不变；检查确认全部差异正是新增两个 production 文件带来的计数变化。已改为单独核验这组明确差异，其余 15 项仍逐字比较。记录见 [audit_initial_failure.json](../../runs/baselines/I12a1_20261001/audit_initial_failure.json)，没有为此改数值代码或容限。

首次打包时，Windows PowerShell 的 JSON 数组被额外嵌套，导致校验器把所有 runtime 行作为一个对象比较并拒绝；改为直接反序列化后逐项检查，未发现真实 runtime 改动。记录见 [package_initial_failure.json](../../runs/baselines/I12a1_20261001/package_initial_failure.json)。

证据入口：[验收状态](../../runs/baselines/I12a1_20261001/verification.json)、[源码审计](../../runs/baselines/I12a1_20261001/source_audit.json)、[回归比较](../../runs/baselines/I12a1_20261001/regression_comparison.json)、[CTest 原始输出](../../runs/baselines/I12a1_20261001/ctest_verified_details.log)、[API](../api/I12a1_Hermite_Gals_Temporal_Phase_Flux_API.md)。

## KNOWN LIMITATIONS

- event detection: **NOT IMPLEMENTED**。
- space-time quadrature: **NOT IMPLEMENTED**。
- crossing benchmark: **NOT EXPECTED TO PASS YET**，本次确认仍拒绝。
- 人工填入满足结构约束的 integrated volume 并不构成数值积分证书；尚无时空事件完整性、CFL/库存有界性或几何闭合保证。
- 来源 owner 为进程内身份；不证明 arbitrary callback 等价/不可变，不提供 restart 序列化，不声称跨任意不同加法顺序逐 bit 一致。

## STOP

**I12a1 已结束。Do not start I12a2.**
