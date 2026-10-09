# I12r6 / native 孤立接触的完整界面连接

## 任务 ID / 标题

I12r6：在同一 native Q4 标量上区分接触点的单元内 0/1/2 条曲线关联，并保留正确的液相连通片。上一轮 I12r5 为 progress；父任务 I12 仍 FAIL。

## 目标

只解决一个拓扑事实：正则二维零集与实际单元边或内部几何 chart 边孤立相切时，不应因奇数个不同边根或单个接触根而拒绝，也不能错误连接仅点接触的液相片。将 I12r5 的根身份接入完整 discovery、体积与 publication 消费。

## 允许修改

- `include/Core/Geometry/Detail/LevelSetNativeInterfaceChain.h`：单元/chart 的接触关联和负相片边界连接。
- `include/Core/Geometry/Detail/LevelSetMappedInterface.h`：同源边根的重根认证、相位、细化与图端点消费。
- 必要时新增一个专用 `Detail/LevelSetNativeContact.h`，以及 `Detail/LevelSetHermiteQuarticLeafData.h` 的最小证书接口；公开 `LevelSetInterfacePanels.h` 只改能力注释。production 总计不超过 5 个文件。
- `tests/interface_hermite_gals_flux_contract.cpp`（更新原完整相切拒绝）和 `tests/interface_hermite_gals_publication_contract.cpp`（完整解析拓扑/体积/法向/曲率矩阵），不超过 2 tests。
- 本 TODO/报告、I08/I11 API、docs README、`runs/baselines/I12r6_20261001/`。

## 禁止修改

冻结方法/PDF/S08 源、时间通量公式、first-jet + M 状态、bubble 可达性/幅值/位移门槛、publication 容限、Hydro/Electricity/App/reader。不得平移几何绕开相切、删根、容差合并近根、用第二份几何或 clipping M。

## 数学约束

1. 接触必须来自原标量的完整、精确认证根库存；局部正则性与进入单元的方向必须有同源证据，不靠一个近似梯度或浮点阈值猜测拓扑。
2. 真正 transverse crossing 保持原结果；偶重接触向单元外有 0 条关联、向单元内有 2 条关联；角点/AMR hanging 身份不能按全局坐标近似合并。
3. 同一几何点保留唯一根身份，所有界面关联、edge phase intervals、graph endpoints 均指向它。内部 chart 不成为物理通量/体积边。
4. 界面连接与液相连通片是不同问题：液相仅跨正长度负相面区间连通，不能穿过一个零测度点。所有 graph 仍须具备原 regularity 证书和非零坐标跨度。
5. graph 端点的重根不等于 transverse integration section 的简单根；不能借无重因子的导数替代原标量的法向或敏感度。
6. 源几何与接受态不可变；失败必须明确拒绝，保留原源。正长度零迹、真正奇异零点及尚未认证的闭子单元环不以本任务放行。

## 实现要求

先复现原圆形 `isolated noncorner contact` 拒绝，认证相切的单元关联及负相连接；覆盖实际边与 chart 边、正反相、uniform/AMR、XY/RZ。沿既有 retained roots/graphs 消费，不新增全局曲线状态，不放宽容限。先以测试发现缺口，再迭代代码。

## 验收

- T0：Release 主程序及相关目标、architecture、无新增 warning、变更预算、冻结/PDF 与前序归档哈希。
- T1：完整曲线数量/闭合/方向、根关联、相片连通、解析逐单元与总相体积、法向/曲率、根细化和 source identity；真实 I07/I11 消费。近相切两侧控制及奇异/零迹拒绝。
- 回归：I01–I11、jet/polynomial geometry、legacy 对照。完整 I12 benchmark 必须实际运行，记录所有改善或仍有失败，不更改验收。
- 本任务 PASS 不代表 Gate I 或电喷雾 production；低黏/无黏与后续耦合关卡仍按路线图进行。

## 失败规则

缺少正则性、根身份、相片连接或同源消费证明时返回 FAIL/BLOCKED 并给出具体缺口，不以减少验证范围替代目标。若预算确需超过 5 production + 2 tests，先说明并拆分，不能继续扩张。时间积分设计变更仍待用户决定，本轮不实施。

## 交付

最小代码和完整解析/回归证据，任务报告、源码/runtime/归档哈希、明确未解决的 I12 缺口；总目标保持 active。

## 执行结果

**PASS（本任务范围）；父任务 I12 仍 FAIL。** 最终 3 production + 2 tests，无新增生产文件。16 项最终回归通过，44 组完整接触/近接触发布矩阵、2 个细化后表示限制/源隔离检查及 3 条 legacy 对照通过。完整 I12 的原 10 项失败和输出保持。全过程及限制见 [报告](I12r6_Native_Isolated_Contact_Report.md)。
