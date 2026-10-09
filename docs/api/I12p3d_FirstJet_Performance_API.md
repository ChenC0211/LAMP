# I12p3d / FirstJet 被动性能诊断

独立 App 输入保持 [p3c](I12p3c_FirstJet_Output_API.md) 不变。新增计时只读 steady_clock，在原 producer 阶段边界记录，不参与状态、步长、几何验收或输出调度；模块仍不写文件。

[HermiteGalsFirstJetStepReport](../../include/Modules/Interface/HermiteGalsFirstJetStep.h) 保留原 transported/prescribed/center 数量，新增 `timings` 及源 publication 的 cells/nativeFaces/cutCells/graphs 数量。timings 单位为墙钟秒，均为诊断：

| report 字段 / performance stage | 范围 |
| --- | --- |
| `input` / `interface_input` | 原输入验证、源网格/时层检查、primary 副本及规模记录 |
| `boundary` / `interface_boundary` | 真实外边界 first jet 查询、master vertex/edge 分类 |
| `vertices` / `interface_vertices` | 原 vertex transport 与 primary 写入 |
| `edges` / `interface_edges` | 原 edge transport 与 primary 写入 |
| `centers` / `interface_centers` | 旧 scalar 的实际 center departure 查询、报告计数及 candidate 构造 |
| `mixed` / `interface_mixed_closure` | 原 mixed AMR closure 和 center prediction 容器 |
| `publication` / `interface_publication` | 原唯一 publication，包括原几何/测量/验收 |
| `total` / `interface_producer_total` | 上述完整 preparation 及计时/诊断开销，不含返回 wrapper 分配 |

App 保留原 `interface_step`，包括 velocity/characteristic 构造与 producer 返回处理；`interface_mesh_transfer`、`io`、`step_total` 也保持。七个 producer 子阶段可相加，`interface_producer_total` 与 `interface_step` 是包含关系，不能再与这些子阶段相加后冒充更大的总耗时。初始态 step=0 没有 producer 调用，只保存原四个初始化/App 阶段；随后每个接受步有 12 行 performance。

重网格步需要区分规模：`interface_step` 和 producer 的八行使用**该次输运源 publication** 的 native cells/faces、cut/graph、primary vertex/edge 与 center 数量；mesh transfer、IO、step_total 使用**最后接受目标网格**的数量。源数量取自原 source 及原 transport 报告，不能把转移后的细网格规模用于解释粗网格上完成的物理输运成本。

没有新增字典开关、时间层或可写物理状态。primary 模式、identity、返回 publication 及原拒绝语义不变。源代码审计从原 p3a 归档去除新增计时/只读数量行后，producer 主体逐字相同；原数值是否保持仍由实际回归和完整窗口核验。

规模测量脚本在任务证据目录，使用同一个最终可执行及固定 OMP2/BLAS1/MKL1，串行独立进程、每配置三次。其 Windows 外部进程内存采样不是新的 production state；保存原始时序和请求间隔，只能报告采样峰值。三个网格档的增长指数是已测范围的观察，不是渐近复杂度证明，也不是其它机器上的绝对吞吐承诺。具体测量范围及验收见 [TODO](../task/I12p3d_FirstJet_Performance_Scale_TODO.md)。
