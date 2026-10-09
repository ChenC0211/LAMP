# I06a / AMR master–hanging 精确 skeleton 约束接口

## 任务 ID / 标题

I06a：开放 Core 既有 quartic master/hanging 约束，独立于 cell interior 与 geometry publication。I06 的唯一最终验收仍是 coarse/fine C1 trace 一致、hanging 不独立推进。

## 目标

只解决一个数学/架构事实：在不构造 cell polynomial、bubble 或第二份 geometry 的情况下，由同一套现有精确 master 约束生成 hanging node 与 fine-subedge 数据，并查询同一个 value/normal trace。

完整 I06 至少需要 Core 公有 header、独立实现 cpp、现有约束引擎、严格浮点 CMake 配置、Interface header/cpp 共 6 个 production 文件，超过路线图的 5 文件预算。因此按架构规范 §8 拆分：先 I06a 数据/算子接口，再 I06b 非均匀 mixed closure 与 I03/I04 接线。不得以本子任务 PASS 宣称 I06 或 Gate I 完成。

## 允许修改

- 新 `include/Core/Geometry/ActiveQuarticHermiteSkeleton.h`、`src/Core/Geometry/ActiveQuarticHermiteSkeleton.cpp`。
- `include/Core/Geometry/Detail/LevelSetHermiteQuarticLeafData.h`：只开放 skeleton 构造/查询，复用既有约束和 exact arithmetic。
- 根 `CMakeLists.txt`：新实现加入已有 strict FP、no-LTO geometry object target。
- 新 `tests/interface_hermite_gals_skeleton_contract.cpp`、`tests/CMakeLists.txt`：注册独立公开 API contract。
- 本 TODO、对应 API/报告、`docs/README.md`、`runs/baselines/I06a_20261001/` 证据。

预算：4 production + 2 test。Core 的通用算子不属于 S08 冻结专有实现。

## 禁止修改

不改 I01–I05 算法、first-jet transport、mesh、Hydro、Electricity、Apps/reader/case、冻结方法 PDF、S08 69 个专有文件；不接入生产方法或输出帧。不引入 mixed closure 取点规则、体积/通量/bubble、regrid/coarsening 算法或新的全局未知量。

## 数学约束

1. 只在 regular master vertices 的 `(phi,phi_x,phi_y,phi_xy)` 和 maximal shared edges 的 `(phi_mid,phi_n_mid)` 注册独立数据；这里 phi_xy 是调用方已经闭合的派生量，Core 不推进它。
2. hanging jet 与 fine-subedge midpoint 为 `D_hang=C_AMR D_master`；fine value 和 normal derivative 分别精确限制同一 coarse quartic trace。输入 hanging/fine payload 不参与计算。
3. 新查询复用原有 rational quartic basis 与 exact homogeneous linear combination，不从舍入后的 hanging values 重建独立 trace。
4. 仅最后的数值查询转为 double；公开 rounded node/edge 值不能替代 exact owner 创建第二份 production geometry。
5. RZ 轴线 phi_r、phi_rz、edge normal derivative 严格零；不自动修正不合格 master 数据。
6. 构造整批成功后发布只读 owner；输入/mesh 不被修改，无 borrowed mutable state、cell centers、cut geometry 或 IO。

## 实现要求

新增中立 Core skeleton input/owner；明确 native vertex、actual leaf-side 与 master parent 元数据。允许 full mesh 和具备完整 coarse 约束支持的 active subset。构造使用现有 QuarticLeafData 的 skeleton-only 路径，不运行 cell interior/sector/geometry 构造。保留旧默认构造与运算次序。所有 exact 运算处于 strict FP TU；并发只读查询不注册新数据。

## 验收

- **T0：** Release 主程序、新 contract、I01–I05/architecture 编译；无新增 warning；公有头独立 include，结果无默认未绑定状态，数据只读。新 cpp 实际编译选项含 `/fp:strict /GL-`。
- **T1 Q4：** shifted/anisotropic XY、RZ even Q4；uniform、单格局部细化、相邻细化、多层 2:1 网格。全部 used node / actual side midpoint / interior edge query 与解析 first/mixed jet 相比 `<=3e-11*(1+|reference|)`；RZ axis 精确零。
- **T1 C1：** 任意非解析 master DOF（固定种子），至少两种 coarse/fine 拓扑；每条 coarse/fine 重合段取至少 9 点，两边 phi、dx、dy、dxy 相差 `<=3e-11*(1+|reference|)`。独立 1D quartic Hermite 公式核对约束，覆盖水平/竖直、perpendicular master edge 端点为 hanging 的情况。
- **T1 ownership：** hanging/fine 输入置 NaN 或改变不影响结果；独立 DOF 数等于 `4*N_master_vertex+2*N_master_edge`，不含 hanging 和 cell center；源 input 修改/析构、mesh 析构后 owner 查询稳定；复制的只读 owner 和并发查询一致。
- **T1 拒绝：** 缺少 master edge、重复/反向/非 side identity、非法 active mask、缺 coarse support、非有限独立数据、非法 RZ parity、域外/非有限 query 明确失败；错误包含 Core skeleton 阶段或查询对象；原输入不变。
- **T2 专项回归：** 原有 native Q4/AMR geometry contract、I01–I05 和旧 Interface state control 通过；I03/I04、旧 Core GALS 数值输出与 I05 证据逐字一致；S08 69 文件和三份方法 PDF 哈希不变。

这些测试只证明 skeleton 约束，不证明 AMR mixed 空间阶、完整输运、曲率、守恒或破碎。

## 失败规则

错误或 C1/Q4/ownership 不满足即 FAIL，在范围内修正且不放宽阈值。若原有 Core 引擎缺失必要数学能力、需要改变冻结方法或超出本子任务预算，说明具体接口缺口并返回 BLOCKED；不加入另一套几何、低阶 stencil、clipping 或 fallback。

## 交付

4 个 production 文件、2 个测试文件、公开 API 说明、验收报告、原始指标/日志、可执行与源文件哈希、归档。仅 I06a 可以记为 PASS；下一项是 I06b，I06 和 Gate I 仍未完成。

## 执行结论

**PASS（仅 I06a），2026-10-01。** 见 [报告](I06a_AMR_Skeleton_Constraint_API_Report.md) 和 [API](../api/I06a_AMR_Skeleton_Constraint_API.md)。Q4、C1、ownership、拒绝、旧 geometry/transport 回归、strict FP 及冻结哈希检查通过。多层测试 fixture 和拒绝消息检查的初始失败原文已保留；生产数学和容限未改。下一项 I06b，I06 整体尚未验收。
