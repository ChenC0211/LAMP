# H02a / CCH 编译所有权与依赖边界报告

**结论：PASS。** 日期：2026-10-03。[TODO](H02a_Cch_Build_Boundary_TODO.md)按默认文件预算完成；H02 的积分散度仍未实现，下一单任务为 H02b。本项不代表 Gate H、流体 MMS 或低黏稳定性通过。

## 改动与数学保持

将 H01 已有 CchDofLayout::build 从 [公开头](../../include/Modules/Hydrodynamics/CchDofLayout.h)原样迁至对应 [CchDofLayout.cpp](../../src/Modules/Hydrodynamics/CchDofLayout.cpp)。建立真实 static library `lamp_hydrodynamics_cch`，其 PUBLIC 依赖为 quantum_core，PRIVATE 依赖为统一编译设置；已有布局契约改为通过该库消费。

root CMake 显式统一 CCH 来源列表，并从旧 Hydro glob 排除同一列表。architecture 单独分类 Cch / Detail/Cch 来源，防止旧 Hydro 编译新 CCH、CCH 编译旧 Hydro、重复 owner、兼容层或 App 反向依赖。旧 compatibility allowlist 保持字节不变，新消费者不能经旧 Hydro header 获得其兼容依赖。详见 [API](../api/H02a_Cch_Build_Boundary.md)。

[迁移审计](../../runs/baselines/H02a_20261003/source_audit.json)逐行比较迁移前后的 build body，除缩进外完全相同；移除 build body 并替换其声明/owner 注释后，公开头的其余 contract 完全相同。没有修改 native 编号、拓扑、RZ 测度、身份、输入拒绝、值隔离、寿命或默认 App 路线。

拆分依据：同时建立新 library/architecture 分类并增加散度契约会修改 tests/architecture、tests/CMake 和新契约三个测试文件，超过默认两文件预算。因此本次只完成 H02a 的编译架构事实；H02b 后续单独实现/验收散度，不提前将 H02 判 PASS。

## 验收证据

标准 local preset Release 构建成功，无新增 warning。定向 CTest **7/7 PASS，0 失败，总耗时 67.30 s**；[构建](../../runs/baselines/H02a_20261003/build.log)、[最终 architecture 编译](../../runs/baselines/H02a_20261003/architecture_build.log)、[CTest](../../runs/baselines/H02a_20261003/ctest.log)及[完整断言](../../runs/baselines/H02a_20261003/ctest_detail.log)已归档。

| 层级 | 检查 | 结果 |
| --- | --- | --- |
| T0 | lamp_help / repository_layout | PASS；标准程序入口与源码布局保持 |
| T0/T1 | architecture | PASS；14 个实际 production/build targets，371 个翻译单元各一个 owner；旧兼容边未扩展 |
| T1 | architecture_self_test | PASS；133 项 policy/scanner 检查，31 个端到端非法变更拒绝 |
| T1 | cch_dof_layout_contract | PASS；经独立库链接，仍为原 23 场景 / 29 拒绝；测试源码保持原哈希 |
| T2 | hydro_gradient_load_contract | PASS；原 XY/RZ native 载荷、积分与拒绝保持 |
| T2 | interface_first_jet_app_contract | PASS；12 smokes / 49 rejections，五输入来源、正式输出及接受/拒绝检查保持 |

新增端到端拒绝包括 CCH→旧 Hydro、CCH→Compatibility、CCH→App，旧 Hydro 编译 Cch 来源、CCH 编译 Core/旧 Hydro 来源、CCH 同来源重复 owner，以及 Cch 源码包含旧 Hydro header。通过项是实际拒绝 fixture，不仅是静态 allow/deny 判定。

[实际生成的构建图](../../runs/baselines/H02a_20261003/architecture_graph.txt)与独立依赖闭包审计确认，CCH 仅到达 `lamp_hydrodynamics_cch / quantum_compile_options / quantum_core / quantum_hermite_geometry / Eigen3::Eigen`。CchDofLayout.cpp 只有 lamp_hydrodynamics_cch 一个 owner；旧 Hydro 没有它。Core 内部原有几何 object 不构成 Interface Compatibility 依赖。

## 来源与预算

[开工清单](../../runs/baselines/H02a_20261003/source_before.json)包含 893 个旧来源/方法/程序文件；最终 889 个保持哈希，四个合法改动是 CchDofLayout.h、root CMake、tests/CMake 和 architecture checker；新增来源只有 CchDofLayout.cpp。旧数值算法与冻结 CCH PDF 保持哈希。

[算例输入清单](../../runs/baselines/H02a_20261003/inputs_before.json)中的 **243 份原 .dict/.json 输入全部保持**。LAMP.exe 保持 SHA-256 `d67fda82f4bb33cce27ba212ca9c94fac93d207f9af6087b96629f04bdd5819b`。H01 的历史报告/归档不改写；当前 source/API 的编译入口以本项为准。

预算为 **3 个 production/build 文件 + 2 个测试文件**，符合 TODO。测试文件仅 tests/CMakeLists.txt、tests/architecture/check.cpp；没有复制布局契约或新增占位 solver/App。临时审计在 runs/checks/H02a，结果冻结到本基线，不成为运行依赖。

复现定向验证：

```powershell
& ./tools/build.ps1 -Target @('LAMP','lamp_hydrodynamics_cch','quantum_cch_dof_layout_contract','quantum_architecture_check','quantum_hydro_gradient_load_contract')
ctest --preset local -R '^(lamp_help|architecture|architecture_self_test|repository_layout|interface_first_jet_app_contract|cch_dof_layout_contract|hydro_gradient_load_contract)$' --output-on-failure
```

CMake/CTest 不在 PATH 时使用本机 VS 随附程序；通用入口保持 ./tools/build.cmd。没有运行全部 CTest、Debug 配置、流体求解或完整物理 benchmark；本项没有改变其数学算法，不能扩大验收结论。

下一单任务为 H02b：在已有真实库中实现唯一 D_I=B A_f，并独立验证 local/global closure、random flux、selected 边界及 coarse–fine subfaces。H03 压力负转置和流体推进仍按路线图随后执行。
