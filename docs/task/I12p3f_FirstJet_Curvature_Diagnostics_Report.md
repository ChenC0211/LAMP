# I12p3f / 同源曲率接受态诊断报告

**PASS（本项）。** 输出规范 §9 要求的曲率 RMS 和 unavailable count 已补齐，全部初始/接受态保存，与空间帧频率无关。范围见 [TODO](I12p3f_FirstJet_Curvature_Diagnostics_TODO.md)，权重/空界面/NA 语义见 [API](../api/I12p3f_FirstJet_Curvature_Diagnostics_API.md)。下一项只读复核完整 Gate I。

修改 1 production（InterfaceFirstJetOutput.cpp）和 1 test（既有输出 PowerShell 契约）。conservation.csv 保留原字段，再追加 kappa_rms、curvature_quadrature_samples、curvature_unavailable_samples、curvature_measure。数据直接读取当前 publication.quadrature，XY 用 ds、RZ 用 2πr ds 和原组合曲率。使用原 report 最大曲率缩放平方；无新几何、重拟合或可写状态，计算开销归入原 IO stage。

当前 publication 已检查 sample/weight 有限且有效，全部实际接受态 unavailable=0，显式计数而非占位。空界面/任何不可用样本的 RMS 为 NA；空界面 measure/counts=0，定义不把缺失数据伪装成零曲率。本项六类 App 问题都有界面，未声称新增空界面 App 基准。

| 验证 | 结果 |
| --- | --- |
| 39 targets Release build | PASS，25.926 s，无新增 warning |
| 实际 App CTest + architecture CTest | PASS，App 62.305 s、architecture .791 s；原接线/精确范围/输出与架构均保持 |
| 新曲率统计 | 18 频率组、90 初始/接受态；字段完整、数量/measure 正、unavailable=0、RMS 在 [0,max] |
| 解析 RMS | XY/RZ Q2 静止与平移对照 1/R、2/R，误差 <1e−10 |
| 频率独立 | 关闭/每步/稀疏三种输出下，统计及解/geometry identity 逐字相同 |
| 与 p3e 比较 | 12 原 smoke +18 频率组，共30组/150初始及接受态；全部旧 conservation 列、solver/events/accuracy/配置/frames 保持 |
| VTK/PVD 内容 | 48帧、96VTK及所有 PVD 字节完全相同；6,264 原同源显示样本验证保持 |
| 失败语义 | 原 CFL/数值拒绝、接受后真实 IO 失败、已有输出保护及受控文件锁原子重试通过 |
| p3e 精确输入诊断 | 15 精确范围拒绝/6 合法边界保留；原24输入拒绝、3 legacy preflight 保持 |

将新增统计块和四列 header/row 去除后，writer 完整内容与任务开始前 p3d 归档逐字相同。Core/Modules/App producer/typed config/reader 未改；38 其它 runtime SHA 保持，仅 QUANTUM.exe 改变，最终 SHA-256 为 `77D2369F4E8D8D5888570EA1A56C727F4D3D5DEFAC7C082A5ABB73D192C4CB1F`。实际执行 2 个 CTest，35 个原数值测试按未改源码/二进制身份复用；p3d 的完整窗口、收敛及规模证据继续有效，本项没有重跑或重新计为这些实验。

架构仍为12 targets/388 TU/918 files/3149 includes/5 legacy edges；69硬冻结、11未改接线参考、3PDF及 p3e/p3d archive保持。空界面 NA 分支及不可用计数属 IO 定义，没有修改 publication 接受条件或允许非法曲率进入物理 state。

证据：[I12p3f_20261002](../../runs/baselines/I12p3f_20261002/)。verification/source_audit/output_comparison 保存原始比较和范围判据；build/CTest、实际输入、错误/退出码、全部输出、脚本与runtime一并归档，最终 SHA 见 archive_verification.json。后续 I12p3r2 核对 p3e/p3f 对此前两项缺口的关闭及当前身份，才裁决完整 Gate I。Hydro/低黏/电喷雾耦合尚未验证。
