# 应用层导航

独立模块与耦合工作流遵守同一 [App 组织合同](../../docs/api/App_Workflow_Contract.md)。头文件位于 `include/Apps`，本目录保存对应实现；数值算法归 Core / Solvers / Modules。

| 位置 | 职责 |
| --- | --- |
| `lamp_main.cpp` / `RunCase.cpp` | CLI、case 目录、统一进程资源与工厂分派 |
| 根目录 `AppRegistry`、`IApp`、`RunConfig`、`RunControls` | 所有 App 的入口合同与类型化调度值 |
| `IO` | 通用字典设施、被动来源记录、纯流体输入/输出组件 |
| `Electric` | 独立电场输入、运行、同源输出 |
| `Hydrodynamics` | 纯流体七文件；流体事务、唯一接受时钟、调度及自身 Output |
| `Electrohydrodynamics` | 独立八文件；电配置/Q/Maxwell、联合接受事务及自身 Output，组合纯流体运行服务 |
| `Interface` | IO 读配置，Setup 准备类型化 provider/网格，App 拥有实际输运循环 |
| `Pic` | IO 读文件和写报告，Setup 准备单位/规范化，RunAdapter 拥有 session 与 accepted bind；显式拒绝预检 |

仅保留五个同名 App 目录和共享 IO；不设 Runtime、Coupling、Verification 生产目录。短验证归 `tests/contracts`，完整研究归 `cases/benchmarks`。新增 App 同步注册工厂、目录/依赖检查与输入输出合同，不创建空职责文件。

旧工作流与输入保持退役，见[清单](../../docs/task/App_Retired_Input_Inventory_20261014.md)。当前数值调度见[主线](../../docs/roadmap/LAMP_Two_Phase_Solver_Mainline.md)；结构整理不改变历史数值验收。
