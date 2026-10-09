# H01 / 原生面速度与单元压力自由度布局报告

**结论：PASS。** 完成日期：2026-10-03。按 [TODO](H01_Native_Face_Cell_Dof_Layout_TODO.md)完成 CCH 原生编号与数据所有权；下一单任务为 H02。尚未装配散度、压力算子或流体时间推进，不能将本项通过解释为 Gate H、低黏或 EHD 通过。

## 实现与范围

新增公开头 [CchDofLayout.h](../../include/Modules/Hydrodynamics/CchDofLayout.h)，包含两个紧凑值类型，详见 [API](../api/H01_Cch_Native_Dof_Layout.md)：

- CchDofLayout：每个参与 selected 域的真实 native face 一个正坐标方向速度位置，每个 selected leaf cell 一个压力位置。2:1 接口使用实际 fine subfaces，两侧共享编号，没有额外 coarse-face unknown。
- CchNativeValues：拥有 Nu 个面速度与 Np 个单元压力的副本，仅提供只读访问；错尺寸、非有限值、不同布局及过期网格明确拒绝。

布局保留现有 NativeLeafDomain 的 compact metadata 和同一不可变 Core mesh identity。pressure 编号直接复用 Core，velocity 紧凑双向映射由 native faces 一次生成；不复制 mesh/face geometry，不依赖 C1/Q8 或 Interface Compatibility。每次 build 产生独立 layout identity；复制共享 identity；输入 vector、原 domain 和原 mesh 的销毁不影响保存数据。内容完全等价的独立 mesh 可通过 Core 身份验证，同尺寸但坐标/拓扑改变的 mesh 必须拒绝。

本项是容器与编号 contract，头内只实现容器操作，几何验证继续属于 Core；不为此增加算法翻译单元或占位 CCH library。新测试只链接 quantum_core。RZ 轴面即使真实 Af=0，仍保留 native 数据位置；轴条件、边界消元和 pressure gauge 留给后续任务，不把此位置宣称为独立可求解速度。

## 验收

标准 local preset 的 Release 配置和构建成功，无新增 warning。定向 CTest **7/7 PASS，0 失败，67.76 s**；[构建日志](../../runs/baselines/H01_20261003/build.log)、[CTest 日志](../../runs/baselines/H01_20261003/ctest.log)与[断言输出](../../runs/baselines/H01_20261003/ctest_detail.log)已归档。

| 层级 | 检查 | 结果 |
| --- | --- | --- |
| T0 | lamp_help、architecture、architecture_self_test、repository_layout | 4/4 PASS；公开头作为测试首个 include 独立编译 |
| T1 | cch_dof_layout_contract | PASS；23 网格场景，29 非法输入拒绝 |
| T2 | interface_first_jet_app_contract | PASS；12 smokes / 49 rejections，包含五输入来源、输出与拒绝检查 |
| T2 | hydro_gradient_load_contract | PASS；现有 XY/RZ native 几何/压力载荷与来源拒绝保持 |

T1 覆盖 XY/RZ 的 2、4、8、16、32 档均匀网格、4×8 非方形域，以及 4、8、16 档初始 2:1 AMR 的全域/左半 selected 域。FaceSegmentGraph 的唯一真实边集合独立提供计数与 cell incidence，核对内部面共享、coarse side 的两片 fine subfaces、轴类型及压力双射。均匀全域严格满足 Np=n²、Nu=2n(n+1)；所有 selected 场景满足 Nu≤8Np，primary 标量数严格为 Nu+Np。

拒绝检查包括 native/DOF 索引越界、空 domain、短/长 vector、NaN/±Inf、同尺寸但不同 selection/RZ/layout、拉伸网格、修改坐标及 AMR 改变。合法零值、有限正负值、输入隔离、layout 复制和 mesh 销毁后的寿命均通过。

复现受影响检查：

```powershell
& ./tools/build.ps1 -Target @('LAMP','quantum_cch_dof_layout_contract','quantum_architecture_check','quantum_hydro_gradient_load_contract')
ctest --preset local -R '^(lamp_help|architecture|architecture_self_test|repository_layout|interface_first_jet_app_contract|cch_dof_layout_contract|hydro_gradient_load_contract)$' --output-on-failure
```

CMake/CTest 不在 PATH 时使用本机 VS 附带的同名程序；通用完整构建/短检查入口仍为 `./tools/build.cmd -Test`。

## 预算与身份

H01 修改预算为 **1 个 production 文件 + 2 个测试文件**：上述新头、新 [cch_dof_layout_contract.cpp](../../tests/contracts/hydrodynamics/cch_dof_layout_contract.cpp)、[tests/CMakeLists.txt](../../tests/CMakeLists.txt)。没有改动旧 Core/Solvers/Modules/App 代码或 root CMake。

[开工清单](../../runs/baselines/H01_20261003/source_before.json)记录 890 个旧来源/程序文件；[最终审计](../../runs/baselines/H01_20261003/source_audit.json)确认其中 889 个保持，唯一变更是测试 CMake 注册。冻结 CCH 方法 PDF 和 LAMP.exe 保持原哈希；新增 production 只有该头。程序 SHA-256 为 `d67fda82f4bb33cce27ba212ca9c94fac93d207f9af6087b96629f04bdd5819b`。

对照 G04 预检运行保存的来源，[原输入审计](../../runs/baselines/H01_20261003/original_inputs_audit.json)确认 24 个原模板的 120 份字典哈希保持。新长涡旋算例是在此前用户请求下构造/启动的独立 Interface 基线，其结果在[单独报告](Interface_Long_Vortex_Verification_Report.md)记录，不计入 H01 流体验收。

最终文档/索引更新后再次通过 repository_layout；[文档检查](../../runs/baselines/H01_20261003/document_audit.json)确认 15 份相关文档的本地链接均可解析，两组新算例可由标准 -List 入口列出。

本次没有运行全套 CTest、流体 MMS、低黏流体 benchmark 或耦合 gate；H01 没有改动这些数值行为。下一项 H02 应从本布局与同一 Core face geometry 装配唯一 integral divergence，并单独写 TODO 后执行。
