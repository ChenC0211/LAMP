# I12p3e / FirstJet 字典范围错误上下文报告

**PASS（本项）。** I12p3r1 中字典规范 §9 的范围错误定位缺口已关闭。Gate I 还需 I12p3f 曲率 RMS/unavailable count 及最终复核。范围见 [TODO](I12p3e_FirstJet_Input_Diagnostics_TODO.md)，异常字段及 reader 定位语义见 [API](../api/I12p3e_FirstJet_Input_Diagnostics_API.md)。

修改 3 production：FirstJetBenchmark.h/.cpp、BenchmarkInput.cpp；1 test：既有 App CMake 契约。Config 原数值条件与检查顺序保持，只将既有异常携带的字段类型化；reader 捕获后调用现有 Node.fail，补上 file:line、限定字段、raw value、期望范围及 method/consumer。没有复制范围条件，没有用错误文本猜字段。typed config 脱离字典仍能独立拒绝非法值。

例如实际 resolution=3 的新日志：

```text
system/interface.dict:4: Interface.Benchmark.resolution = '3'; expected for method=hermiteGals / Interface.FirstJet: [Interface FirstJet App] resolution must be in [4,256]
```

15 个实际新增拒绝涵盖 resolution/steps 的下界、上界、负值，dt=0/负值/有限 dt 导致 total-time 溢出，writeEverySteps 下界/上界/负值，以及 regrid 乱序/最终步/超出最终步。各字段置于不同输入行，逐项核对实际行号、限定路径、完整 raw value、范围/method、非零退出及无 output/cache。6 个合法边界 preflight 涵盖 resolution=4/256、steps=1/1,000,000、正的极小 dt、频率上界；原频率 1 及其它正常配置由原有效组覆盖。

| 验证 | 结果 |
| --- | --- |
| 39 targets Release build | PASS，42.579 s，无新增 warning |
| 实际 App CTest | PASS，62.230 s |
| 新范围/边界 preflight | 15 精确拒绝 + 6 合法输入 |
| 原 App 接线契约 | 12 smoke/48 步/36 regrid，24 原 invalid，1 CFL，3 legacy preflight 保持 |
| 原输出契约 | 18 频率组、48 帧/96 VTK/6,264 同源样本；真实 IO 失败、数值拒绝、旧输出保护和文件锁重试通过 |
| 合法 App 与 p3d 对照 | 12 组完整 stdout 仅去掉独立运行路径 banner 后逐字相同；conservation/solver/events、accuracy、配置/追踪和 frames 字节保持 |
| architecture CTest | PASS，.856 s；12 targets/388 TU/918 files/3149 includes/5 legacy edges 保持 |

从已校验 p3d/p3c archive 提取原 App/header/reader，与任务开始 SHA 匹配；剥除 typed error 包装后，Config 判据及顺序、App 其余完整主体、primary config 布局、reader 除异常映射外的完整内容逐字相同。Core、Modules、output writer 和数值测试源码未改；39 runtime 中仅 QUANTUM.exe 改变，38 个其它可执行 SHA 保持。

本项实际执行 **2 个 CTest**（App 与 architecture）。原 35 个数值测试的完整输出及对应未改源码/二进制身份复用 p3d 证据，没有声称本项重跑 37 项、24 长窗口或规模研究。有效 App 本次实跑了上述短窗及输出频率组。它们与源码/运行文件身份共同证明错误定位改动没有进入时间步或几何算法。

69 硬冻结实现、11 未改接线参考、3 PDF、p3c/p3d/p3r1 archive SHA 保持。唯一 reader 参考扩展属于此前明确允许的纯入口，旧 reader 主体未动。最终 QUANTUM.exe SHA-256：`A7C916AD335B0F210504F13A1A71B2DD44DEFE8913C11B727981588B7ABB453F`。

证据：[I12p3e_20261002](../../runs/baselines/I12p3e_20261002/)。verification、source_audit、numeric_identity 给出判据；保留实际输入/错误日志/退出码、所有有效输出、build/CTest 及可重现脚本。独立新目录使 stdout 的 Case banner 不同，比较只剔除此一行，未剔除数值或 geometry identity。package_evidence.ps1 归档并逐项校验 payload/runtime/doc links，最终标识见 archive_verification.json。

下一项 I12p3f 只补同源曲率时序诊断。当前纯路线没有独立 M/通量 ledger；本项没有声称 Hydro、低黏稳定或电喷雾耦合已通过。
