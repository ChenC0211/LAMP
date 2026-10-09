# I02 / 回溯映射与同公式 Jacobian 交付报告

**结论：PASS。** 执行日期：2026-10-01。本轮只完成 I02；下一任务为 I03。Gate I 与低黏/无黏电喷雾完整目标仍未完成。

任务定义：[I02 TODO](I02_Characteristic_Map_Jacobian_TODO.md)。接口、公式、时间层和边界语义见 [API 说明](../api/I02_Hermite_Gals_Characteristic_API.md)。开工前复核 I01 的 8 个交付文件及归档 SHA-256，确认前序未变；三份冻结方法 PDF 仍匹配 S01。

## 交付与数学事实

| 文件 | 变化 |
| --- | --- |
| [HermiteGalsTrace.h](../../include/Modules/Interface/HermiteGalsTrace.h) | timed paired velocity view、只读带来源的 B/J 结果及 V2 tracer 声明 |
| [HermiteGalsTrace.cpp](../../src/Modules/Interface/HermiteGalsTrace.cpp) | 时间/身份校验、共用 Core midpoint 调用、边界 tangent 缺失语义及错误上下文 |
| [interface_hermite_gals_trace_contract.cpp](../../tests/interface_hermite_gals_trace_contract.cpp) | 独立解析、数值微分、子步、时间收敛和失败路径检查 |
| [tests/CMakeLists.txt](../../tests/CMakeLists.txt) | 注册新 tracer contract |

本次是 V2 对通用 Core 算子的明确接入与独立验收，**没有重新实现第二套 tracer**。Core 原已具备共用控制流的 value-only/tangent midpoint。新接口从明确端点时间推导实际 trace interval，成对消费同一插值分支的速度值与一阶导数，返回具有来源 geometry identities/time layers 的结果。

已验证 `J_Bh` 与实际 `B_h` 的离散导数一致；多个成功子步按实际顺序组合。生产路径没有空间差分 Jacobian、速度 Hessian、反射/裁剪、RK4 自动回退或旧 Interface state 共享。遇到物理边界时保留事件信息，Jacobian 为空；调用内部 pullback accessor 会明确失败。

## 原始数值结果

证据目录：[runs/baselines/I02_20261001](../../runs/baselines/I02_20261001/)。所有阈值均在实现前 TODO 中指定，首次正式构建与测试即通过，未调整容限。

nonlinear map 的中央差分检查覆盖 12 个位置和全部四个矩阵分量。每个 dt 使用 `epsilon=0.02,...,0.000625` 共六级扰动，以下为最细一级；完整 18 行见 `finite_difference.csv`。

| 请求 dt | 最大分量误差 | 最后一级差分收敛阶 |
| --- | --- | --- |
| 0.037 | 4.57821e-10 | 2.000083 |
| 0.113 | 1.38928e-9 | 2.000030 |
| 0.210 | 2.55796e-9 | 2.000009 |

两成功子步案例使用不同端点的非交换梯度；value-only Core 与 V2 tangent 路径的采样序列和映射值完全相同，四分量数值微分误差为 **1.64074e-12**，见 `substep_difference.csv`。affine shear 的位置和矩阵符合 `2e-14*(1+|reference|)` 的既定阈值。轴线静止点保留了非单位 Jacobian，没有将“位移为零”误当作单位映射导数。

对 `omega(t)=0.6+0.4*t` 的解析旋转，在 `[0,1]` 组合回溯并分别检查位置和全 Jacobian。`N=8,16,32,64`，交替序列的相邻 dt 比为 1:2。以下为 N=64 的结果；完整 8 行见 `time_convergence.csv`。

| 时间步序列 | 位置误差 | J 最大分量误差 | 位置阶 | J 阶 |
| --- | --- | --- | --- | --- |
| uniform | 5.70918e-6 | 2.95097e-5 | 2.001928 | 2.000892 |
| alternating | 7.64209e-6 | 3.94895e-5 | 2.007506 | 2.005890 |

两种序列的最后两级阶均在预定 `[1.85,2.15]` 内。这证明所测 tracer 对已知时间插值场的二阶时间精度，**不证明完整 Interface/Hydro 耦合时间阶或曲率收敛**。

## 构建、回归与范围

| 检查 | 结果 |
| --- | --- |
| Release：QUANTUM、新 tracer contract、I01 contract、原 Interface contract、architecture | PASS，退出 0，10.312 s；无新增 compiler warning |
| CTest：quantum_help / architecture / I01 state / I02 tracer | 4/4 PASS，0.09 / 0.74 / 0.08 / 0.08 s |
| 原 `quantum_interface_state_contract --parametric-state-control` | PASS，退出 0，0.373 s |
| 原有代码树/输入/工具哈希 | 1,541 文件未变；唯一既有改动为 tests/CMakeLists.txt，另增 3 个源码/测试文件 |
| S08 冻结专有实现 | 69/69 哈希一致 |

预算：2 个 production 文件、2 个测试文件。Core、I01 类型、旧数值实现、Apps、reader、默认配置和 case 输入均未修改。SuperLU 既有 `/LTCG` 自动重链接信息保留。Git 元数据不可用，没有伪造 commit。

测试原始日志保存在 `ctest_details.log` 与 `legacy_interface.log`；构建命令、编译器、线程及结果见 `verification.json` / `build_result.json`。`runtime/` 和 `runtime_hashes.json` 保存实际执行文件身份。`evidence.zip`、`evidence_inventory.csv`、`delivery_manifest.json` 保存交付副本和逐文件哈希。

## 尚未完成

I03/I04 skeleton transport、I05 closure、AMR、守恒通量、bubble 求根、I11 同源几何发布和 I12 benchmark 均未在本轮开始。`conservativeHermiteGals` 仍未进入 reader 或生产 App。

边界命中时间/边界历史的 tangent 仍需专门物理边界处理；本接口明确暴露该缺口，不伪造内部 Jacobian。任意 callback 的 value 与 derivative 是否真为同一插值仍是 provider 前置条件；I02 没有构造 Hydro 速度插值或全局梯度场。Core 的半单元采样位移、active-domain、非光滑分支及 RZ parity 约束在 API 说明中明确列出。

S01 已知的 `hydro_moment_xy` 首步失败仍保留；本次短回归不能将它记为通过。本轮结束后可按路线图进入 I03。
