# H02a / 原生 CCH 编译入口

状态：2026-10-03，H02a PASS；[报告](../task/H02a_Cch_Build_Boundary_Report.md)。本项为 H01 已有布局建立真实编译所有权；其后 [H02b 积分散度](H02b_Cch_Integral_Divergence.md)、[H03 压力载荷](H03_Cch_Pressure_Load.md)与 [H04 封闭域参考求解](H04_Cch_Single_Phase_Saddle.md)已 PASS。[H05](H05_Cch_Physical_Divergence.md)已通过，[H06](H06_Cch_Energy_Viscosity.md)已 PASS，下一项 H09；可运行 CCH 流体 App 尚未实现。

## 消费与所有权

公开 API 仍为 [CchDofLayout / CchNativeValues](H01_Cch_Native_Dof_Layout.md)。CchDofLayout::build 已从头内原样移到 `src/Modules/Hydrodynamics/CchDofLayout.cpp`，其唯一 CMake owner 为 static library `lamp_hydrodynamics_cch`；其余短容器操作保持头内实现。

新消费者链接真实库：

```cmake
target_link_libraries(my_hydro_service PRIVATE lamp_hydrodynamics_cch)
```

该库 PUBLIC 导出 quantum_core 及其 include/Eigen 使用要求，PRIVATE 使用统一编译设置。H07a 新增正式只读 quantum_interface 依赖，用于消费唯一 publication；闭包包含其 quantum_linalg/SuperLU。没有旧 Hydro、Interface Compatibility、App 或 Runtime 依赖；Core 内部的 quantum_hermite_geometry object 属于已有通用几何。

root CMake 的 `_lamp_cch_sources` 显式列出新路线来源，并从旧 Hydro glob 移除同一列表。新增 CCH .cpp 时只在该列表添加，不编译进旧 quantum_hydrodynamics。H01 原契约现经新库链接，仍核对原生 DOF/AMR/identity 与寿命。

## 边界检查

architecture 将 `Modules/Hydrodynamics/Cch*` 和 `Modules/Hydrodynamics/Detail/Cch*` 的头/实现归为 HydrodynamicsCch，新 target 是该编译层的唯一 owner。旧 Hydro 源码继续归原 owner，不移动、不新增兼容边。

- CCH 不能链接旧 quantum_hydrodynamics、lamp_interface_compatibility 或 App/Runtime；也不能包含旧 Hydro header 后间接获得其兼容依赖。
- CCH 可按模块层级使用 Core/Solvers 和正式 Interface；H02a 历史交付仅依赖 Core，当前 H07a 明确新增正式 Interface 并验证真实构建闭包。该依赖不代表完整流体耦合或时步已实现。
- 错 source owner、重复编译、反向模块依赖、App include directory 继续被拒绝。

H02a 迁移保持 H01 的构造顺序、native 编号、Nu+Np、RZ 轴位置、layout/mesh identity、字段复制与非法输入语义。该迁移本身未引入配置开关、边界消元、D_I、C、Hodge 或时间层；D_I 由后续 H02b 单独实现。

## 验证与后续

local Release 构建和 7/7 定向 CTest 通过；布局仍为 23 场景/29 拒绝，架构自测 133 项及端到端 31 个拒绝均通过。实际依赖闭包见 [source_audit.json](../../runs/baselines/H02a_20261003/source_audit.json)。

后续 [H02b](../task/H02b_Integral_Divergence_Report.md)已使用该库实现唯一 D_I=B A_f，并通过 random flux、内部抵消、selected 边界与 2:1 局部闭合；H02a 本身不代替该数学验收。[H03](H03_Cch_Pressure_Load.md)已从同一矩阵生成 C=-D_I^T 并通过压力功伴随，[H04](H04_Cch_Single_Phase_Saddle.md)已完成单相封闭域参考求解与均匀 XY 连续 MMS，当前下一单任务 H09（[H06](H06_Cch_Energy_Viscosity.md)已 PASS）。
