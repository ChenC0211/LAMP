# 应用层导航

[代码组织规范](../../docs/QUANTUM_Code_Organization_Standard.md)是职责和依赖的唯一规则入口；[模块化主计划](../../docs/QUANTUM_Modular_Consolidation_Plan.md)维护进度，目录中不另抄历史待办。`include/Apps` 与本目录对应。

| 入口 | 职责 |
| --- | --- |
| `quantum_main.cpp`、`RunCase.cpp` | CLI、具体 App 注册、输入来源和共享服务绑定；不拥有物理时间循环 |
| `Runtime` | 公共运行协议、注册/执行和运行控制 reader；不包含具体 App 或物理模块头文件 |
| `SingleModule` | 独立模块注册及 PIC 会话循环；Plasma 模块接收物理输入，提供整步和检查点完成服务 |
| `Hydrodynamics` | Hydro 主循环、Hydro/EHD 共用且无主循环的 `HydroCase`；绑定公共域、界面及载荷 |
| `Electrohydrodynamics` | EHD 主循环、电学时层接线及跨模块诊断 |
| `Electrospray` | 喷嘴主循环；Setup 绑定相/边界，Coupling 交换模块数据，IO 输出，Compatibility 单向读取旧输入 |
| `Electricity`、`Interface` | 正式独立运行入口；使用与耦合运行相同的模块服务，明确各自支持范围 |
| `Verification/Pic`、`Verification/Mesh` | 验证专有输入和基线处理，调用 Plasma/Solvers/Core 原服务；不把验证开关塞回物理模块 |

方程装配、界面输运/体积修正、电荷迭代、粒子推进和物理历史由所属模块负责。App 保留单位/几何身份绑定、唯一主时钟、模块间接受条件及输出调度。新增能力同时更新所属 reader、实际消费者和[字典参考](../../docs/QUANTUM_Dictionary_Reference.md)，并运行受影响的现有检查。
