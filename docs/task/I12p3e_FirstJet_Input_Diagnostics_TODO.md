# I12p3e / FirstJet 字典范围错误上下文

目标：关闭 I12p3r1 中字典规范第 9 节的实际缺口：范围/关联检查失败也必须有 file:line、限定字段、原始值、有效范围和 method/consumer。

允许修改：InterfaceFirstJetBenchmark.h/.cpp、InterfaceBenchmarkInput.cpp（3 production）；interface_first_jet_app_contract.cmake（1 test）。本任务 API/TODO/报告、当前规范状态及 runs 证据脚本。

禁止修改：Core、Interface 数学内核、output writer、合法输入/默认/schema、容限、步长、旧 reader 主体/历史报告/冻结实现/PDF；不混入曲率诊断。

数学约束：原 config.validate 数值判据和合法输入返回值不变，仍能脱离字典执行防御校验；只让其已存在范围错误带类型化 Field 标识，reader 使用现有 Node.fail 添加输入上下文，不复制另一套数值条件。kernel 不接收 Node/path。失败发生在 check-case/生成网格与输出之前。

实现要求：类型化异常覆盖 resolution、steps、dt/total-time、writeEverySteps、Regrid.atSteps；保留原错误含义，reader 将 field 映射回实际 Node，并说明 hermiteGals。没有通过消息字符串猜字段的映射。旧 reader 分支逐字保持。

验收：T0 39 targets 构建、无新增 warning、架构依赖通过；T1 实际 App 输入契约扩展所有上下界、负值/零、总时长溢出、乱序/最终步 regrid、frequency，逐项断言路径/行号/字段/raw/range/method 与非零退出、无 output/cache。原 12 有效 smoke/48 步/36 regrid、24 原 invalid、CFL 拒绝、legacy preflight、输出频率/真实 IO/原子重试保持。T2 数值 kernel/test executable 身份保持，合法 App 结果与 p3d 逐字对照；不因错误消息任务重跑已有同一数值实现的长窗口或规模实验。

失败规则：合法输入路径/数值、config 独立防御或旧入口改变即 FAIL；少一个范围上下文也 FAIL。仅外部缺接口才 BLOCKED。

交付：3 production/1 test 以内修改、API/报告、实际输入/合法输出比较和范围/身份归档。本项通过后再做 I12p3f；Gate I 仍需最终复核。

执行结果：**PASS**，见[报告](I12p3e_FirstJet_Input_Diagnostics_Report.md)。15 精确范围拒绝/6 合法边界、原 App/输出及架构通过；合法数值、typed 判据、旧 reader 和全部数值 kernel 保持。实际执行 2 CTest，原 35 数值测试按未改源码/运行文件身份复用。
