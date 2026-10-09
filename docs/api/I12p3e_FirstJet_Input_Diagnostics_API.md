# I12p3e / FirstJet 配置错误定位

合法输入、范围、默认和 consumer 保持 [App API](I12p3b_FirstJet_App_API.md) / [输出 API](I12p3c_FirstJet_Output_API.md) 原定义。

`InterfaceFirstJetBenchmarkConfig::validate()` 对既有范围错误抛出 `InterfaceFirstJetConfigError`，仍派生自 `std::invalid_argument`，原 `what()` 含义保持。其只读 `field` 是 `Config::Field::{Resolution,Steps,Dt,WriteEverySteps,RegridSteps}`；没有字典 Node、路径或行号。直接 typed App 调用仍执行同一防御检查。

reader 捕获该类型后，通过 field 找到实际输入 Node，调用已有 `Node.fail()` 输出：

```text
system/interface.dict:4: Interface.Benchmark.resolution = '3'; expected for method=hermiteGals / Interface.FirstJet: [Interface FirstJet App] resolution must be in [4,256]
```

dt 字段承担原 `finite(dt) && dt>0 && finite(dt*steps)` 检查；总时长溢出会定位 dt 并说明 total time 要求。Regrid 的严格递增与 `< steps` 检查定位 `Interface.Benchmark.Regrid.atSteps`。类型、未知字段、重复字段等继续使用原 strict Node 诊断。禁用空间场时显式未使用的 schedule 仍按原规则拒绝。

异常字段映射由类型完成，不读取错误文本猜字段，也不在 reader 复制范围判据。没有新增 schema、数值开关或 fallback。`--check-case` 的拒绝保持非零退出、不构造网格、不创建 output/cache。执行范围见 [TODO](../task/I12p3e_FirstJet_Input_Diagnostics_TODO.md)。
