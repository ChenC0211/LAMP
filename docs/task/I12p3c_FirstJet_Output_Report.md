# I12p3c / 纯路线正式输出报告

**PASS。** 纯 Hermite–GALS 独立 App 已接通只读 VTK XML、accepted-frame、series/meta、配置追踪和失败诊断。完整 Gate I 仍缺 I12p3d 效率/规模证据；不进入 Hydro。使用方式及字段含义见 [API](../api/I12p3c_FirstJet_Output_API.md)，范围见 [TODO](I12p3c_FirstJet_Output_TODO.md)。

生产修改为新 `InterfaceFirstJetOutput.h/.cpp` 及既有 FirstJetBenchmark.h/.cpp、InterfaceBenchmarkInput.cpp，共 5 文件。测试修改为原 App CMake 契约及新 PowerShell 输出契约，共 2 文件。Core/Modules、数学算法、容限、旧 benchmark 和冻结输入均未改。

`writeFields=true` 开启真实叶网格 VTU 和同源曲线 VTP，可用 `writeEverySteps` 控制频率；`false` 仍保存全部接受态的标量时序及元数据。显示点从当前 publication.sample 获取法向/组合曲率，显示折线不反哺计算。所有库存/通量 ledger 明确 NA/null。初始和实际接受帧进入 PVD，拒绝步可定位阶段/原因/接受时钟。IO 失败保持非零退出，已有 output 不允许覆盖。

| 验证 | 结果 |
| --- | --- |
| 39 targets 最终 Release 构建 | PASS，23.006 s；无新增 warning；只有 QUANTUM.exe 改变 |
| 原 36 CTest | PASS，368.502 s；35 个完整数值输出逐字相同 |
| 最终 App CTest | PASS，约 61.96 s；原输入/运行契约与扩展输出契约均通过 |
| 原 App 基础契约 | 12 短算例/48 步/36 regrid、24 输入拒绝、1 原 CFL 拒绝、3 旧入口预检 |
| 输出频率对照 | 六问题 × 关闭/每步/稀疏，共 18 组，数值及几何 identity 完全相同 |
| 输出内容检查 | 48 帧、96 VTK 文件、6,264 个界面采样，独立解析 XML/CSV/JSON；体积、时钟、identity 与 stdout 保持 |
| 失败与重试 | 数值拒绝、已有输出保护、接受后真实写入失败、实际文件锁触发原子替换重试、4 个非法频率均通过 |
| 24 实际完整算例 | PASS，2,072 接受步/36 regrid，461.072 s（含预检）；116 帧/232 VTK、2,096 初始及接受时序行 |
| 完整数值对照 | 每个接受步、最终误差、动态转移、几何 identity 与 p3b 逐字相同；原 p1/p2b 模块数值 token 也保持 |
| ParaView 6.1.1 实读 | XY/RZ、非多项式初值、固定/动态完整涡旋及轴对称平移 PVD/VTU/VTP 均成功读取 |

完整矩阵保持 p3b 的 resolution=8、原六问题、fixed/dynamic × uniform/AMR、全部 dt/终时；VTK 频率为 max(1,floor(steps/3))，未缩短窗口。每步仍写守恒诊断、solver、performance，初始/最终帧保证输出。最终 24 组没有触发 IO 重试；受控文件锁测试单独实际触发 retry=1，释放锁后正常退出 0。接受后写入失败测试实际退出 3，accepted_step=3、正常帧截至 step=2；没有用缺失退出码冒充失败。

输入 SHA-256 用系统 Get-FileHash 独立校验。字段长度、native quad 连接、轴对称容量 π、同源体积和正向法线、Q2 静止/平移曲率、元数据时钟、NA 量及输出频率均有机器断言。两个额外非多项式 AMR XY/RZ 案例各接受 2 步并粗化一次，生成字段并由 ParaView 实读。复杂完整涡旋文件也由真实 ParaView reader 读取，不只检查扩展名或 XML 合法性。

架构统计为 **12 targets、388 TU、918 文件、3148 内部 includes、5 条原 legacy edges**。S08 的 69 冻结实现、11 个未改接线参考、3 PDF 和 p3b archive 哈希保持；新增纯 reader 功能仍位于允许扩展的第 12 个接线参考中，旧 reader 主体与共用验证逻辑单独核验。历史 crossing probe 仍按原样退出 1。38 个其它运行文件 SHA 不变。

开发过程中保留的失败/修正：

- 首次构建成功但有两处 int64 step 到 int 的 C4244；输出时钟和拒绝接口改用 int64，后续最终构建零 warning。
- 第一次完整矩阵在 fixed AMR vortex 的 step=313 被 Windows 拒绝替换 manifest 临时文件，数值已接受，正常帧截至 266。前 313 步与 p3b 完全相同，failure.json 和非零退出均保留。原子替换现对 Windows permission-denied 做最多 6 次同一完整 payload 的有界重试，总退避 630 ms；不删除旧文件、不隐藏持续错误。随后从新目录完整重跑 24 组通过。
- 首次新增文件锁 fixture 实际重试并成功运行，但 Windows PowerShell 未保留子进程 handle，ExitCode 为 null，成功断言失败。现显式保留 handle，并要求失败 fixture 的退出码非 null；最终受控重试及真实 IO 失败分别确认 0/3。失败日志、旧目录和修正后目录全部保存。
- 首次 ParaView pilot 尝试保存用户设置时提示写入不可用；改用 `--disable-registry` 后读取不访问该设置文件，后续验证通过。

证据根目录：[I12p3c_20261002](../../runs/baselines/I12p3c_20261002/)。`verification.json` 为数值/范围结论，`output_verification.json` 和 `output_full_comparison.json` 为输出审计，`first_full_failure_analysis.json` 定位上述部分窗口失败；`paraview_*.log` 保存实际 reader 结果。五次 App 契约目录、原部分完整运行、新完整运行、程序/DLL、输入及所有日志/场文件由 package_evidence.ps1 归档，逐条检查 SHA；最终标识见 archive_verification.json。

当前性能文件是实际 App 阶段的诊断，不能当作效率/复杂度已通过的证明。下一项 I12p3d 记录规模增长、阶段成本与实际状态/内存，再复核 I12p3。没有声称非守恒路线已经解决质量—动量耦合或低黏/无黏电喷雾。
