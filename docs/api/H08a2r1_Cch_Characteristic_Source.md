# H08a2r1 / 同源 characteristic 来源

现行接口：[r11 第二次恢复](../task/H08a2r11_Common_Regular_Source_Deployment_Resume2_Report.md) **PASS**，共同正则场与同源运动已部署。r11 为最后基础设施关口；[H08a](../task/H08a_SpaceTime_Material_Flux_Report.md) 已完成规定范围的同源时空积分、材料一致性与小余项收敛验收。本 API 只提供同源运动，完整流体耦合仍未验收。

```cpp
auto geometry = CchCharacteristicGeometry::prepare(mesh, reconstruction);
auto source = geometry.linear(mesh, previousActualB, arrivalActualB);
source.validateFor(mesh, previousActualB.identity(), arrivalActualB.identity());
auto trace = source.traceBackward(point, arrivalTime, departureTime, {maxStep});
```

`prepare` 持有同一 reconstruction，按实际 cell/support 闭包构造 BVH 与重叠面邻接，不另造运动网格。`pieceVelocityJet` 查询给定片的原物理值/J，闭包两端可查询。getter 不保存 mesh 裸引用；每批消费前验证真实 mesh/regrid。

`linear` 要求两端 B 属于同一 prepared reconstruction 和 momentum geometry，boundary law 同一实例或均为空。previous 为 Initial/End，arrival 为下一 step 的 End，arrival.begin 精确等于 previous.time 且区间正长。保留实际 B 身份和完整 clock；坐标相同不代表实例相同，标签也不证明 native 值的物理时间或 H14 提交。

`sample(piece,x,time)` 在线性物理时钟上插值同一端点的值/J，同时返回同源时间导数。只接受实际片闭包和完整 slab 内时间；它不是任意 Core 速度 callback，也没有 FirstJet 接受事务。

`traceBackward` 使用固定物理时间槽。完整内部槽由二阶预测曲线上的一次 Picard 积分更新，保持三阶时间精度；沿真实 support 分片，各片九点 Gauss 完整积分 Q4 tensor/线性时间的 17 次被积函数。节点位置、区间宽度、连接根和时间都随实际离散 map 求导。最终三次 dense 路径采用同一端点修正并重查事件；命中真实非 C1/物理边界的预测路径继续按三次事件分支处理，不跨域取样。公式与负控见 [r17](../task/H08a2r17_Nonlinear_Crossing_Jacobian_Accuracy_Design.md)。

只有两端都通过完整物理 Q4 迹核验的连接才透明：五个 quarter 节点检查值和按片宽缩放的四项梯度，固定 Bernstein 换基计入全部绝对工作量，预算 128 epsilon，无幅值 floor。零场为零预算。它是浮点解析证书，非任意病态输入的有向区间证明，没有用户容差开关。

舍入量级的积分修正仅在点及 J 两者都落入各自工作量预算时共同保留原公式，保持精确 terminal 连接的原语义；真实修正仍执行。此分支不放宽精度门槛，不用独立 J 代替实际点映射。原非定常 affine 独立离散多项式及 r16 原两可信连续流参考分别核验离散与连续准确性。

`CchCharacteristicConfig` 要求正 `maxStep`，默认最多 4096 子步、256 实际几何事件；超预算拒绝。`events()` 统计真实横穿，`splitEvents()` 统计结束 dense 子步的不光滑/外部事件，二者不等同。`continuityReport()` 给出总连接、透明连接及最大相对迹差，只用于观察。

完整内部结果含 point、time、initial/final piece、子步和事件及可选 J。J 是固定物理起止时间的同一离散映射导数。起点/终点位于分片面时 `oneSided()` 保留几何骨架标记。H08a 的直接消费反例要求零场/轴向/对角常平移在真实 vertex/edge/center 可用：自动选择仅允许通过完整 C1 证书和实际面连接相连的局部片；断开的点接触、未证连续连接、真正物理边界及非法显式分支仍拒绝。C1 骨架与同时穿连接不需要另拟合 J 或混合单侧值。

原 64 epsilon 根并合与 128 epsilon 连接舍入带保持。并合根附近的终点/求积点保持实际坐标，选择包含该点的真实 C1 连通片，不把点夹回旧片或外推其多项式；下一槽重新选取该点的内向 C1 片。原时间/point/J/证书/工作预算门槛不变。该修正及 3268 个真实内部骨架轨迹的解析检查属于 H08a 消费接入，不另立 source 精修任务，详见[实施前设计](../task/H08a_SpaceTime_Material_Flux_Implementation_Design.md)。

外部/selected-domain 命中返回实际法向、片和时间，没有内部 J。轴保留原零法向及单侧正则，不反射或 clip。步末根沿用 64 epsilon 舍入带，透明连接仍从包含实际点的分支取样；零时长为精确恒等。越 slab、非法片、nonfinite、非标准 FP、错误来源/时钟/regrid 明确拒绝，失败不返回部分内部轨迹。

当前共同 C1 来源已关闭旧切向跳跃的精度缺陷。真实 vertex/edge/AMR/axis DOF 的运动消费与 pullback 接入属于 [H08a 时空相通量积分](../task/H08a_SpaceTime_Material_Flux_Upgrade_TODO.md)，不再单列桥接或 source/Jacobian 精修任务。只有可复现且直接阻断 H08a 的数学反例才允许例外。本 API 不发布步内几何、相通量或材料提交；不能把三个独立 snapshot 插值当成同源时空路径。
