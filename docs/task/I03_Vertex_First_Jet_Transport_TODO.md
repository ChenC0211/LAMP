# I03 / master vertex first-jet 输运

任务依据：[05 路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md) I03、[01 数值约束](../01_LAMP_Numerical_Invariants.md) I-R1/I-R2/I-R4，冻结方法 §4.2 的式 (20)–(21)。执行日期：2026-10-01；开工复核 I02 交付及归档。

## 目标

只解决一个数学事实：对选定的 independent master vertex，用同一旧 Hermite scalar 和 I02 离散回溯得到 `phi_new=phi_old(B_h)`、`grad_new=J_Bh^T grad_old(B_h)`。

## 允许修改

- 新建 `include/Modules/Interface/HermiteGalsVertexTransport.h`、`src/Modules/Interface/HermiteGalsVertexTransport.cpp`。
- 新建 `tests/interface_hermite_gals_vertex_contract.cpp`，修改 `tests/CMakeLists.txt` 注册该测试。
- 本 TODO、`I03_Vertex_First_Jet_Transport_Report.md`、`docs/I03_Hermite_Gals_Vertex_API.md`、`docs/README.md` 的任务索引。
- 证据目录 `runs/baselines/I03_20261001/` 及正常构建产物。

预算：2 个 production 文件、2 个测试文件。

## 禁止修改

- 不修改 I01/I02、Core 几何/回溯公式、S08 冻结实现、Hydro/Electricity、Apps、reader、case 或默认参数。
- 不开始 I04 edge 输运、I05 mixed closure、I06 AMR transfer、I07–I10 守恒闭合或 I11 新几何发布。
- 不重拟合第二份 scalar/gradient，不重新归一化梯度，不输运 Hessian/phi_xy，不裁剪守恒量，不加入 legacy fallback。
- 不把部分 vertex 结果冒充完成的 Interface state 或耦合接受步。

## 数学约束

- phi 和 gradient 必须在同一 departure point 从旧 `HermiteGalsState` 的只读 snapshot 一次取出；Jacobian 来自 I02。
- `J` 行为 departure 分量、列为 arrival 坐标，使用转置拉回；点位置不动不代表 `J=I`。
- 显式目标列表是已存在的 master vertex native indices；严格有序、无重复、不得含 hanging/未知节点。允许 active transport band 选择子集；完整闭合不在本任务内。
- 旧 velocity 的 geometry identity/time layer 必须匹配源状态；arrival 时间层明确。mesh、domain、axis 和网格尺度必须一致。
- physical boundary hit 无 boundary-history tangent 时明确失败；departure 不在旧活动几何内不得外推。
- 批量结果在所有顶点成功后才返回；失败不能修改源状态、库存、edge 或旧 snapshot。

## 实现要求

- [x] 复核前序与六份规范，记录 source baseline 和开工 TODO。
- [x] 新建只读 vertex result，携带源 identity/time、arrival velocity identity/time 及实际更新的 native vertex IDs。
- [x] 从真实 mesh 取 arrival 坐标；复用 I02 tracer 和 Core native scalar query，仅计算 first-jet pullback。
- [x] 入口校验 mesh/domain/时间/身份、非空目标集合及其 master 归属；含不支持 corner-sector 一侧 jet 的来源明确拒绝。
- [x] 失败信息包含 source/arrival identity、阶段及 vertex index；不修改输入，不预先发布部分结果。
- [x] 完成解析 Q4 + affine velocity、复合 scalar 数值微分、全顶点 XY/RZ 非线性解析速度的检验及前序回归。

## 验收

- **T0：** Release 构建主程序、I01/I02/I03 contract、architecture；header 可独立首个 include；无新增 warning；architecture/help 通过。
- **T1：** 正确构造的真实 native Q4 snapshot 与 affine shear/translation 的解析解相比，phi/两个 gradient 分量误差 `<=2e-12*(1+|reference|)`；非对称 J 检测转置方向。
- **T1：** 内部 master vertices 上，对实际复合 scalar `phi_old(B_h(x))` 做独立中央差分，gradient 最大分量误差 `<2e-8`。
- **T1：** 零步/静止速度不漂移；固定位置但非平凡 J 的边界/轴线点正确更新梯度；来源及结果只读、不改变 phi_xy、edge、M 或源几何。
- **T1：** 未知/重复/hanging target、错误时间/identity/mesh/domain/axis/尺度、boundary hit、离开活动几何、callback 中途失败均明确拒绝，无部分状态修改。
- **T2（vertex 专项）：** 使用全网格顶点与真实 Q4 scalar，XY logistic 速度和 RZ 奇对称 cubic radial 速度均有解析 flow map。单步 `dt=0.2,0.1,0.05,0.025` 的 phi/gradient 最后两级局部误差阶在 `[2.8,3.2]`，对应 midpoint RK2 局部三阶误差。不能据此宣称多步完整 Interface 收敛。
- **回归：** I01/I02 contracts 通过；原 `quantum_interface_state_contract --parametric-state-control` 通过；旧源码/输入（除测试注册）与冻结实现哈希保持不变。

## 失败规则

- 公式、来源、转置或失败隔离不符合要求判 `FAIL`，仅在本任务范围修复；不改容限配合结果。
- 如需修改 Core/前序契约或超预算，先说明最小必要原因并拆分；无法在范围内完成则 `BLOCKED`，不改用解析 callback 冒充 production source。
- 尚未实现的边界历史/单侧 sector 物理处理显式拒绝，不反射、外推或假造 tangent。

## 交付

vertex kernel、独立 contract、TODO/报告/API 说明、原始误差与阶数、构建/回归/哈希证据。只做 I03；结束状态为 `PASS / FAIL / BLOCKED`，下一项为 I04。

## 执行结论

**PASS，2026-10-01。** 见 [交付报告](I03_Vertex_First_Jet_Transport_Report.md) 与 [API 说明](../api/I03_Hermite_Gals_Vertex_API.md)。Release 构建通过，CTest 5/5 和旧 Interface 状态回归通过；下一项 I04，Gate I 尚未通过。

静止验收的来源解释：初版测试误将“无输运漂移”写成与原始 primary 导数逐位相等。Core 的已有 normalized metric cache 在一个节点产生 `2.77556e-17` 的舍入差，结果与旧 snapshot 的同点采样逐位相同。最终检查对同一 snapshot 保持逐位相等，并以 `8*epsilon*(1+|primary|)` 独立限制缓存舍入误差。没有修改生产计算、Core、解析误差容限或收敛阶范围；原 TODO、首次失败和定位日志均保留在证据目录。
