# I04 / master edge-midpoint first-jet 输运

日期：2026-10-01。依据：[05 路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md) I04、[01 数值硬约束](../01_LAMP_Numerical_Invariants.md) I-R1/I-R2/I-R4/I-R5，冻结 Hermite–GALS 方法 §4.2 的 edge 更新公式。开工复核 I03 交付与归档，前一轮归类为 progress。

## 目标

只解决一个数学事实：在真实 master leaf side 的 midpoint 上，从同一旧 Hermite snapshot 计算 `phi_new=phi_old(B_h)`、`g_new=J_Bh^T grad(phi_old)(B_h)`，保存 `dn_new=n_e dot g_new`。

## 允许修改

- 新建 `include/Modules/Interface/HermiteGalsEdgeTransport.h`、`src/Modules/Interface/HermiteGalsEdgeTransport.cpp`。
- 新建私有 `src/Modules/Interface/HermiteGalsFirstJet.h`，将 I03 的来源校验/单点 first-jet pullback 提取为共用代码；相应修改 `src/Modules/Interface/HermiteGalsVertexTransport.cpp`，保持原行为和原测试不变。
- 新建 `tests/interface_hermite_gals_edge_contract.cpp`，修改 `tests/CMakeLists.txt` 注册测试。
- 本 TODO、同目录 I04 报告、`docs/I04_Hermite_Gals_Edge_API.md`、`docs/README.md` 与 `runs/baselines/I04_20261001/` 证据。

预算：4 个 production 文件、2 个测试文件。

## 禁止修改

不改 Core、I01/I02、公有 I03 API、I03 测试或冻结实现；不开始 I05 closure、I06 AMR restriction/transfer、I07–I10 守恒闭合、I11 geometry publication；不改 Hydro/Electricity、reader、App、case 或默认参数。

## 数学约束

- edge 以原生端点 indices 标识；沿 +x/+y 排列。水平边 normal 为 +y，竖直边 normal 为 +x，不使用 outward cell normal 或液体界面 normal。
- edge 的两个量来自相同 departure point、同一旧 snapshot、同一 I02 map/Jacobian。必须先拉回完整 gradient，再取 normal 分量。
- 显式目标为非空、严格字典序、无重复的 source master edge 集合；每条必须是实际 native leaf side，不能为更长 coarse side 的受约束子边。允许合法 master 边以 hanging vertex 为端点；判据是边的包含关系，不是简单禁止 hanging endpoint。
- midpoint 从真实端点生成，要求可表示在边的严格内部；不吸附到端点或虚构位置。
- 来源 identity/time、mesh/domain/axis/最小 leaf 尺度规则与 I03 一致。固定点仍使用实际 J。
- 只返回不可变 edge stage result；失败无部分发布、不改变 source vertices/edges/M/closure/snapshot。

## 实现要求

- [x] 核验六份规范、冻结公式、I03 哈希并保存开工 TODO。
- [x] 提取 I03 已有单点公式与来源校验，让 vertex/edge 使用同一执行链，不改变 I03 结果。
- [x] 只读带 source/arrival identity/time 的 edge result，完整明确的 normal 与目标语义。
- [x] native topology 核验、来源列表核验及 midpoint 可表示性检查；在进入 velocity callback 前检查全部目标。
- [x] 错误包含 source/arrival identity、edge endpoints 和阶段；boundary-history tangent 缺失、sector、活动几何覆盖不足与非有限结果明确拒绝。
- [x] 完成算子数值验证、来源/拓扑负面契约和前序回归，交付文档与原始证据。

## 验收

- **T0：** Release 构建主程序及 I01/I02/I03/I04/旧 Interface/architecture；新头文件首个 include；无新增 warning；architecture/help 通过。
- **T1：** 真实 native Q4 scalar 配合平移、两个方向非对称 affine shear，水平/竖直边各自的 phi/dn 解析误差 `<=2e-12*(1+|reference|)`，覆盖完整内部 edge 集合。
- **T1：** 对实际复合 scalar `phi_old(B_h(x))` 沿 edge normal 中央差分，覆盖两种方向，最大误差 `<2e-8`；检测转置、法向与旧切向梯度混合。
- **T1：** 零步/静止结果与旧 snapshot midpoint 采样逐位相等；与原 DOF 的差异限制为 `8*epsilon*(1+|DOF|)`。固定 midpoint 具有非单位 J 时仍更新 dn；RZ axis 上径向 dn 保持零。
- **T1：** empty/duplicate/reversed/unknown/non-native/受约束 subedge、错误来源/时间/mesh/domain、boundary hit、离开旧活动几何及 callback 中途失败明确拒绝；完整 source primary/closure/identity/snapshot 不变。
- **T2（edge 算子专项）：** 全网格 edge midpoints 的 XY logistic 与 RZ 奇对称 cubic radial 解析速度，单步 `dt=0.2,0.1,0.05,0.025`，phi/dn 最后两级局部阶在 `[2.8,3.2]`；保留 RZ axis 与 physical boundary 边。只证明 RK2 局部精度，不替代 I12 多步界面收敛。
- **回归：** I01/I02/I03 contracts 与原 `quantum_interface_state_contract --parametric-state-control` 通过；I03 affine/FD/局部阶原始指标与已归档结果逐位一致；冻结文件及方法 PDF 哈希不变。

## 失败规则

不满足公式、正常方向、native master 身份或失败隔离即 FAIL，只在上述范围修复；不改容限配合结果。若需修改 Core/前序数学契约或超预算，说明最小原因并拆任务。无法在范围完成则 BLOCKED，不以 legacy fallback、解析 callback source 或伪造 tangent 替代。

## 交付

edge kernel、共享 first-jet 代码、独立 contract、TODO/报告/API 说明、解析误差/收敛原始数据、构建回归及 SHA-256 证据。结束状态为 `PASS / FAIL / BLOCKED`；本轮只做 I04，下一项 I05。

## 执行结论

**PASS，2026-10-01。** 见 [交付报告](I04_Edge_First_Jet_Transport_Report.md) 和 [API 说明](../api/I04_Hermite_Gals_Edge_API.md)。Release 构建通过；CTest 6/6 与旧 Interface 状态回归通过。I03 原测试未改，11 行 17 位精度数值输出与 I03 归档逐字一致。

初次编译有测试声明错误，已拆开不同类型的 `auto`；首次 CTest 的所有数值项目通过，但不可表示 midpoint 的测试网格在 Core 构造阶段就被拒绝。最终检查该上游拒绝，不伪造非法 source 来覆盖防御性分支。预设数值容限与阶数范围未改；两次初始失败和原 TODO 均已留档。下一项 I05，Gate I 尚未通过。
