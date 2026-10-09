# Interface V2 state/API（I01）

**2026-10-02 当前扩展：** 依用户授权新增 `HermiteGalsPrimaryState.transport`，见[纯路线设计](../task/I12_Pure_Hermite_Gals_Design_Change.md)。显式 `FirstJet` 状态的 `phaseVolumes` 必须为空，局部 Q4 系数由实际回溯中心值直接提供；几何体积仅在 publication 的只读测量结果中存在。历史 `ConservativeVolume` 默认值维持旧接口行为。下文 I01 原始说明中的独立 M/体积闭合条款仅适用于后者；当前接口实现进展见 [I12p1 TODO](../task/I12p1_Pure_Hermite_Gals_Publication_TODO.md)。

本说明对应 [I01 TODO](../task/I01_Interface_V2_State_API_TODO.md) 与 [交付报告](../task/I01_Interface_V2_State_API_Report.md)。接口位于 [HermiteGalsState.h](../../include/Modules/Interface/HermiteGalsState.h)，命名空间为 `Quantum::Interface`。

**当前只完成状态、所有权和结构校验。** 没有实现新的输运、闭合或几何生成算法，`conservativeHermiteGals` 尚未注册到 reader，也不能作为 production 方法运行。I02–I12 的验收仍按路线图逐项执行。

## 数学对象与存储

| 类型/字段 | 含义 | 归属 |
| --- | --- | --- |
| `HermiteGalsPrimaryState.vertices` | master vertex 的 `phi, dx, dy`，物理一阶导数 | 输运主量 |
| `HermiteGalsPrimaryState.edges` | master edge 的精确参数中点处 `phi, dn` | 输运主量 |
| `HermiteGalsPrimaryState.phaseVolumes` | native cell 的负相守恒体积 `M_K`，含纯相单元 | 独立守恒主量 |
| `HermiteGalsPrimaryState.geometryCells` | 存有高阶多项式的活动 native cell 集合 | 几何支持域 |
| `HermiteGalsClosure.vertexMixedDerivatives` | 与 master vertex 顺序相同的 `phi_xy` | I05/I06 派生闭合量 |
| `HermiteGalsClosure.cellBubbleCoefficients` | 与 geometryCells 顺序相同的 `c_K` | I10 派生守恒闭合量 |
| `HermiteGalsGeometryView` | 保有唯一 Core native enriched-quartic snapshot 的只读视图 | Interface 发布；其它模块消费 |

所有整数 ID 都是当前 mesh 的 native **index**，不是 label。vertex、cell 数组按 index 严格递增；edge 用两个端点 index 唯一标识，不能用在 AMR 后含义不清的数组偏移代替几何身份。edge 端点沿物理正 x（水平）或正 y（垂直）排序；`dn` 分别沿 **+y** 或 **+x**，不是单元外法向。edge 全梯度是 I04 的输运中间量，状态只保留方法指定的值与法向导数。

master side 可以是跨两个 fine native faces 的完整 coarse leaf side。hanging vertex 不允许作为独立 vertex jet。完整 master 集合、fine subedge 约束、边的最大性及 exact restriction 由 I06 的拓扑/闭合实现产生和验证；I01 没有偷偷生成 hanging 值或重新拟合多项式。

`geometryCells` 必须被 inventory 覆盖；`phaseVolumes` 可以覆盖更大的物理域。因此不能仅在 cut/active band 内保存全部液相库存。容器不会从零水平集重算或覆盖 `M_K`，也不会 clipping。`measure` 明确区分 `dx dy` 的平面面积与 `2*pi*r dr dz` 的轴对称体积，不自动换算。

`c_K` 指 `phi_boundary + c_K*q4(xi)*q4(eta)` 中的系数，存储语义与旧 patch 的绝对 `cellCenters` 不同。I10/I11 应按明确的 boundary extension 将它转换为最终多项式数据，不能把两种字段名直接互换。`phi_xy` 和 bubble 均不出现在 first-jet 主量中。

## 生命周期和唯一 owner

1. App 提供当前物理时间层，Interface 数值流程生成主量；`HermiteGalsCandidate(mesh, time, primary)` **复制**主量。外部数组之后的修改不会改变候选。
2. I05/I06/I10 以候选主量计算闭合，I11 从同一份闭合主量生成 snapshot。`bindGeometry(mesh, sourceIdentity, closure, Snapshot&&)` 是接收这些结果的所有权接口。
3. `std::move(candidate).finish()` 消费已绑定候选，产生只读 `HermiteGalsState`。**这不是数值接受检查或 coupled commit。** App 可以先把该完成态用于 midpoint corrector、跨模块验收，所有检查通过后才替换自己持有的 accepted handle。
4. 拒绝候选时销毁其对象即可。没有写入旧 accepted state、修改主时钟、替换网格或登记输出帧的副作用。

`HermiteGalsState` 可复制；复制共享同一份不可变 authority，不能通过任一副本修改 primary、closure 或 snapshot。geometry view 同样保有整个状态寿命；原 state handle、原输入数组和原 mesh 释放后，仍可查询同一 snapshot。网格修改后的复用通过 `validateMesh()` 检查精确 native geometry/connectivity identity。

API 不接受外部 `shared_ptr<Snapshot>`。绑定会把 Core snapshot 移入内部 `const` 对象，原 snapshot 进入 moved-from 状态，避免外部持有可写别名并替换已经发布的对象。每个主量 revision 只能绑定一次，不能为同一状态安装第二份法向/曲率/体积几何。

## 身份、时间与失效

每次构造候选以及每次成功 `replacePrimary()` 都分配新的非零 `HermiteGalsIdentity`。相同时间、相同数组内容也不会复用身份；这是进程内 revision identity，**不是** restart ID、mesh hash 或 accepted-step counter。时间层由 App 传入，容器只校验，不推进时钟。

`HermiteGalsTimeLayer` 包含 `step, begin, end, position`：

- `Initial`：`begin == end`，可以表达非零时间的 restart；
- `Midpoint`：`begin + (end-begin)/2`；
- `End`：`end`；
- advancing interval 必须有限、严格递增且 `step > 0`。

geometry view 同时带 identity 与 time layer。corrector 的两个候选即使采样时间相同，也应通过不同身份隔离几何相关 cache。后续 cache 仍须纳入 material/boundary/topology 等自身依赖，不能仅使用矩阵尺寸或本身份代替所有依赖。

`replacePrimary()` 成功后同时清除 closure 和 snapshot。修改主量必须经过该入口或新建候选，没有可写数组 getter；旧身份计算出的闭合结果会被拒绝。替换失败时保留旧候选。绑定失败在检查阶段不移动传入的 snapshot；全部校验完成后才接管所有权。moved-from 或已消费的 handle 查询明确报错。

## I01 校验的范围

已实施：

- 时间、measure、native index 范围、唯一性/顺序、edge 方向、hanging vertex 角色；
- 有限的一阶数据和闭合数据，有限且非负的 inventory；
- geometry support 必须有 inventory；closure 数量与主量布局对应；
- 精确 native mesh identity；同一个 primary revision；禁止重复绑定；
- snapshot 必须为 native enriched-quartic，cell 数量、活动集合及原生角点边界一致。

**未实施且不能从绑定成功推导的性质：** master skeleton 完整性与 AMR C1、一阶输运正确性、`phi_xy` polynomial reproduction、`0 <= M_K <= V_K` 的完整物理测度约束、唯一通量/reflux、bubble root/topology guard、snapshot 的 polynomial/volume 与主量逐项一致。这些是 I02–I11 的数值生产者前置条件，必须由相应任务实现并验证。I11 必须从同源数据构造 snapshot，不能拿任意同尺寸几何冒充该候选。I01 的 `finish()` 只证明所有权装配完成。

`alpha` 和 `Q_l` 的数值发布接口留待 I07/I08/I11，不为其创建空的 production 实现或默认零值。新状态也没有与旧 `InterfaceState` 共享可写数组、进入 parametricGals 或 KKT fallback 的入口。

## 测试和使用边界

`interface_hermite_gals_state_contract` 使用单元内解析平面及真实 Core snapshot 验证所有权，不以 mock 替代几何对象。它覆盖输入副本、接受/拒绝隔离、时间/身份、原对象释放后的读取、失效和失败原子性，并包含小数坐标的原生边界绑定。这里的解析赋值不是 V2 初始化或几何发布算法。

Gate I 仍需 I02–I12 全部通过；在此之前不接 Hydro 新主线，不开始电喷雾算例调参，也不把状态契约测试当作低黏/无黏稳定性的证据。
