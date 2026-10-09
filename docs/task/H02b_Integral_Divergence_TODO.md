# H02b / 唯一原生积分散度 D_I

**状态：PASS，2026-10-03；H02 数学验收通过。** 结果见[报告](H02b_Integral_Divergence_Report.md)，下一单任务为 H03。

任务依据：[现行路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md) H02、H-R1/H-R2、冻结 CCH 方法 §6.1/V1。前置 [H01](H01_Native_Face_Cell_Dof_Layout_Report.md)和 [H02a](H02a_Cch_Build_Boundary_Report.md)已 PASS。本项只解决一个数学事实：D_I=B A_f 是唯一连续性算子，实际 fine subfaces 与真实面测度决定其系数。

## 目标

在独立 CCH 库中准备不可变、绑定 H01 layout/mesh identity 的稀疏积分散度矩阵。通过同一矩阵计算各 pressure cell 的净外流量，验证 random flux 与 2:1 AMR 局部闭合，为 H03 直接生成 -D_I^T 提供唯一来源。

## 允许修改

- `include/Modules/Hydrodynamics/CchIntegralDivergence.h`、`src/Modules/Hydrodynamics/CchIntegralDivergence.cpp`：公开不可变 prepared operator、同一矩阵的 apply 与身份/有限性拒绝。
- root `CMakeLists.txt`：将新 .cpp 添加到已有 CCH 唯一源码列表；旧 Hydro 自动排除同一列表，不增加新 library。
- `tests/contracts/hydrodynamics/cch_integral_divergence_contract.cpp`、`tests/CMakeLists.txt`：独立 geometry/flux 参考、解析积分、AMR 和拒绝契约，链接实际 CCH 库。
- 本 TODO、API、报告、当前文档索引/路线图。证据位于 `runs/baselines/H02b_20261003`，临时检查位于 `runs/checks/H02b/<时间戳>`。

预算：3 个 production/build 文件 + 2 个测试文件。

## 禁止修改

- H01 layout/value 实现、Core/Solvers、旧 Hydro、Interface/Electricity/Plasma/App、字典、case、冻结 PDF、现有容限与默认运行路线。
- 不实现 H03 压力载荷、H04 saddle solve、H05 轴边界消元/正则、H06 黏性或质量/动量/界面力/时间推进。
- 不保存第二套 production face geometry、B、normalized-divergence stencil 或 coarse-face unknown；不增加参考修正、边界猜测、clipping 或 fallback。

## 数学约束

\[
(D_I u)_c=\sum_{f\subset\partial c}s_{cf}A_f u_f,\quad D_I=B A_f.
\]

- u_f 按 H01 正坐标方向编号；Core negative cell 的贡献为 +A_f，positive cell 为 -A_f。selected 域外的 cell 不进入矩阵；所选域边界仅有一个贡献。
- 同一内部面两项系数严格相反；实际 coarse side 使用两个 fine subfaces 的贡献之和，不装配替代粗面通量。
- XY 面测度为真实边长；RZ 使用 Core 的 full 2πr 测度。轴面 A_f=0，保留原 velocity 列，其全部系数精确为零。
- D_I 为 Np×Nu 的积分算子，不除以 cell volume，也不乘 dt、density、phase fraction 或 pressure。后续归一化只可由它派生。
- 只保存 prepared 稀疏系数及已有 layout handle；非零项≤2Nu、≤8Np。构造与 apply 随 native mesh/稀疏规模线性增长，身份只在批次入口验证。
- matrix 只读、复制共享；apply 必须验证实际 mesh 和值的 layout identity。错误来源和非有限输出明确拒绝，原字段/矩阵保持。

## 实现要求

- [x] 保存当前源码、相关测试、PDF、程序及所有 case 输入的开工身份；复核冻结定义和 H02a 前置。
- [x] 每个参与 native face 只查询一次 Core 几何，构造唯一 compressed sparse D_I；保留零轴列，不存另一个 stencil。
- [x] apply 仅使用该矩阵及同布局 CchNativeValues，输出 cell 积分净流量；拒绝 stale mesh、错 layout 和有限输入导致的溢出。
- [x] 以独立 FaceSegmentGraph 的逐 cell 有向边积分核对 random flux、内部列抵消、selected 边界与真实 coarse–fine 聚合。
- [x] 验证 XY/RZ（轴上及离轴）、非方域、非二进制精确坐标、五档均匀、三档 AMR、部分/不连通域、解析非零散度和多项式流函数局部零散度；网格重建后重新绑定并拒绝旧来源。
- [x] 标准 CMake 构建与 T0/T1/T2 通过；报告范围和未运行项，裁决 H02b / H02，下一项 H03。

## 验收

- **T0：** local Release 构建；公开头独立作为首 include；architecture、自测、repository_layout、lamp_help PASS；无新增 warning；新来源仅有 CCH 一个 owner。
- **T1：** 矩阵尺寸、只读/寿命、每列 incidence/符号/测度、压缩和 nnz 规模正确；32 组随机速度/场景与独立有向图积分相符。local defect=|DIu−reference|/Σ|A_fu_f|、global defect=|ΣDIu−boundary_flux|/Σ_all oriented |A_fu_f| 均≤1e-13；零尺度必须精确为零。内部列系数和精确为零，粗面由实际子面求和。解析多项式净外流和流函数局部闭合使用同样 1e-13 相对通量尺度，不放宽容限。
- **T2：** 重建均匀/细化/粗化网格后，重新准备的 CCH 算子通过同一解析闭合，旧绑定明确拒绝；这是算子重建验证，不声称物理 velocity transfer。H01 布局契约、现有 hydro_gradient_load_contract、interface_first_jet_app_contract 保持。旧代码、程序、方法和输入哈希保持；本项不宣称流体求解或低黏 gate。

## 失败规则

若散度符号、局部闭合、零轴、来源/寿命或有限性拒绝不通过，判 FAIL 并在范围内修复，不降精度。若必须更改 Core 几何、H01 layout 或超过预算才可继续，记 BLOCKED 并说明缺失 contract，再拆任务；不引入第二套通量/散度。

## 交付

唯一积分矩阵及 apply、独立契约、API、TODO/报告、标准构建与检查日志、开工/最终哈希与实际依赖图。H03 尚未通过，H02 通过也不代表 Gate H 或可运行流体模块。
