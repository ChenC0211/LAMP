# H02a / CCH 编译所有权与依赖边界

任务依据：[现行路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md) H02，前置 [H01 PASS](H01_Native_Face_Cell_Dof_Layout_Report.md)。H02 拆为 a（真实编译边界）和 b（唯一积分散度）：新增库需要更新 architecture；若与散度实现/测试同做，将涉及 architecture、tests/CMake 和新契约三个测试文件，超过默认两文件预算。本项只解决编译架构事实，H02 的数学验收仍由 H02b 完成。

## 目标

让已实现的 CCH 原生布局拥有独立 production 编译入口及唯一翻译单元 owner，依赖 Core，避免通过旧 quantum_hydrodynamics 传递链接 Interface Compatibility。保持 H01 编号、值隔离、身份、拒绝和寿命语义。

## 允许修改

- `include/Modules/Hydrodynamics/CchDofLayout.h`、`src/Modules/Hydrodynamics/CchDofLayout.cpp`：将已有 build 实现原样移入对应翻译单元，其余值 contract 保持。
- root `CMakeLists.txt`：新增真实 static library `lamp_hydrodynamics_cch`，显式列出 CCH 源码并从旧 Hydro glob 移除，链接 Core 和编译设置。
- `tests/CMakeLists.txt`：已有 cch_dof_layout_contract 改为链接上述库。
- `tests/architecture/check.cpp`：分类新 target、检查 Cch 翻译单元唯一 owner、阻止 CCH 链接旧 Hydro/Compatibility；更新对应实际 build graph、自测与拒绝 fixture。
- 本 TODO、报告、当前 API/索引/路线图；身份/日志在 `runs/baselines/H02a_20261003`，临时检查在 `runs/checks/H02a/<时间戳>`。

预算：3 个 production/build 文件、2 个测试文件。不新增占位算法或 executable；新库包含原 H01 的真实实现。

## 禁止修改

- Core/Solvers、旧 Hydro 数值实现、Interface/Electricity/Plasma/App、输入 reader、case、方法 PDF 和现有数值容限。
- 不装配 H02b 的 D_I，不实现 H03 负转置压力载荷、H04 saddle solve 或流体推进。
- 不增加兼容 include allowance，不链接旧 quantum_hydrodynamics 或 lamp_interface_compatibility，不改变默认 App 路线。

## 数学约束

- H01 的原生 face/cell 双射、Nu+Np、实际 fine-subface 共享、RZ 轴位置、selected 域和 exact mesh/layout 身份保持。
- build 算法原样迁移；不改变索引顺序、输入拒绝顺序或数值。
- CCH 源码只有一个 CMake owner；旧 Hydro 拥有其原来源，不能把 CCH 新实现编译两次或归入 Core。
- 新 target 的实际依赖闭包不得包含旧 Hydro、Interface Compatibility、Apps 或 Runtime。

## 实现要求

- [x] 保存开工前 production、CMake、受影响检查、冻结 PDF、程序与原输入身份；复核 H01 前置。
- [x] 将 CchDofLayout::build 原样移入 .cpp，建立 Core-only 独立 static library。
- [x] root CMake 显式统一 CCH 源码列表，旧 Hydro 排除该列表；既有布局契约通过新库消费。
- [x] architecture 分类与唯一 owner 检查覆盖 CCH；拒绝 Compatibility、旧 Hydro、App 反向链接和错 owner，不扩大旧兼容边。
- [x] 标准构建和 T0/T1/T2 通过，核对迁移前后源算法、程序、原字典和已有测试输出。
- [x] 报告预算、来源、实际依赖及未实现项；H02 入口更新到 H02b。

## 验收

- **T0：** local preset Release 构建；lamp_help、architecture、architecture_self_test、repository_layout PASS；无新增 warning。实际 CMake graph 包含新库，每个 Cch .cpp 恰有一个 owner。
- **T1：** 原 cch_dof_layout_contract 经新库链接，仍为 23 场景/29 拒绝；自测拒绝新库→旧 Hydro/Compatibility/App、旧 Hydro 编译 Cch 源、CCH 编译旧 Hydro 源及重复 owner；原 compatibility allowances 逐条有效且未增加。
- **T2：** 原 hydro_gradient_load_contract 与 interface_first_jet_app_contract 保持；源码审计确认旧数值算法、原 case、冻结方法和运行程序未变。构建图独立验证新库实际依赖闭包只含编译设置与 Core/Eigen。

## 失败规则

本项编译/所有权/身份拒绝错误判 FAIL 并限域修复。若隔离必须修改旧数值路线或超出预算，记 BLOCKED，说明缺失接口并再拆任务，不引入 fallback。不得以 H01 契约通过宣称 H02 的散度已实现。

## 交付

独立真实库、迁移实现、编译边界拒绝检查、更新后的 API/索引、TODO/报告、原/新来源及检查日志。完成后只允许下一单任务 H02b；整体 H 阶段目标保持。

**最终状态：PASS。** 7/7 定向 CTest、133 项架构自测、31 个端到端拒绝；布局原 23 场景/29 拒绝保持，迁移 body/其余 contract、旧算法/程序/PDF 和 243 份输入保持。详见 [报告](H02a_Cch_Build_Boundary_Report.md)。H02b 的数学散度尚未实现。
