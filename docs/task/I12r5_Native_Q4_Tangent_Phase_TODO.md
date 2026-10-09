# I12r5 / native Q4 孤立相切的面相区间

## 任务 ID / 标题

I12r5：保留 native Q4 线段的精确重根并正确构造相区间，解除 I08 在孤立面内接触处的拒绝。上一轮 I12r4 为 progress；父任务 I12 仍 FAIL。

## 目标

只解决一个几何事实：同一 native Q4 标量在线段上具有代数精确重根时，相测度仍有定义。保留全部不同根、根包围与相区间；偶重根不翻相，奇重根翻相，不将面切向导数为零等同于二维界面梯度奇异。

## 允许修改

- 新 `include/Core/Geometry/Detail/PolynomialSquareFree.h`：有限 exact-expansion 多项式最大公因子/无重因子商，供同一标量的根身份认证，不新增物理表示。
- `Detail/LevelSetMappedSegment.h`：仅 native Q4 接入、同一查询根细化与精确相位判定。
- `Detail/LevelSetMappedSegmentSensitivity.h`、`LevelSetMappedSegment.h`：非横截根的变分明确拒绝及文档，最多 4 production 文件。
- `tests/interface_jet_transport_contract.cpp`、`tests/interface_hermite_gals_flux_contract.cpp`：最多 2 tests。原 jet regression 将 native Q4 的内部重根拒绝写成预期，须同步更新。新增因子/缩放验证通过实际 native mapped query 执行，原 `polynomial_root_enclosure_contract.cpp` 保持其 strict isolator 拒绝合同。
- 本 TODO/报告、I08 API、docs README、`runs/baselines/I12r5_20261001/`。

## 禁止修改

冻结方法 PDF/S08 源、时间通量公式、M/bubble/publication 接受门槛、完整 graph discovery、Hydro/Electricity/App/reader。不得将有限近根误差当重根、合并近根、抹去系数尾项或 clipping M。时间积分设计变更尚未批准，本轮不实施。

## 数学约束

1. 原 regular simple-root 路径及旧 cubic/uniform 行为保持；新认证仅适用于 native Q4。
2. 仅 exact expansion 的零余式可证明公因子，使用 division-free 多项式运算及统一正 2 的幂缩放。商与原多项式的代数关系必须核验；near-critical 不能被误判为 exact multiple。
3. distinct-root enclosure 仍由原 bounded isolator 给出；所有根保留，重根只去代数重数，不删几何根。相区间符号来自原多项式的精确有符号样本，不能沿用无条件翻相。
4. arithmetic/隔离不能证明时继续保留明确拒绝；原 dyadic subdivision 仍只细分同一标量。不得换为近似多项式或另一份几何。
5. 根细化必须收紧同一根的包围并保持原查询不变。相切处负相测度通常无普通一阶导数，变分接口必须拒绝，不能用小分母返回巨大伪导数。
6. 共享 face flux、XY/RZ 测度、内部抵消和容量容限保持；positive-length zero trace 仍单独拒绝。完整拓扑/发布的相切能力不由本任务部分结果推断。

## 实现要求

先在精确系数上验证重复线性/二次因子、非二进制有理根、无理根、多重数、端点、正负缩放、近重根；再通过 native mapped segment 和 I08 XY/RZ/AMR 实际消费。报告新接受范围及仍拒绝边界。不得增加仅供测试、不被实际面相区间消费的 production 分支。

执行中范围整理：先前临时直接算术测试已证明 30 个精确因子/缩放配置；为保持 2-test 预算，将其恢复并把覆盖转至真实 native Q4 查询，同时更新上述既有 jet 断言。该整理不改变 production 算式或扩大预算；初始日志保留。

## 验收

- T0：Release 主程序/相关 targets、architecture、无新增 warning、变更预算/frozen/PDF hash；上轮归档与当前交付身份核验。
- T1：独立制造因子根/相区间、反向线段、root refinement、变分拒绝、正长度零迹与近根控制；I08 原 interior tangent 拒绝变为正确液相通量，保留库存/守恒验证。
- 回归：I01–I11、polynomial/jet geometry 合同、legacy 对照；原完整 I12 benchmark 仍实际运行并记录全部失败/差异。
- T2 不能用本任务通量级接触处理代替完整拓扑发布、moving benchmark 或 Gate I。

## 失败规则

算术位宽、根包围或相位不确定即拒绝且保留源状态。完整相切修复若需要超过文件预算，保留明确分层边界、后续单独接线；不得顺手改其他几何或方法。若通量层本身未达到以上验收，本任务不得 PASS。

## 交付

最小数值代码、根/相位/通量合同、原始通过和失败证据、前后源码/runtime/归档哈希、报告与后续缺口。目标保持 active，不进入 H01。

## 执行结果

**PASS（本任务范围），父任务 I12 仍 FAIL。** 4 production + 2 tests；最终 16 项回归各自通过，8 组实际相切通量、24 组有理重根相区间、2 组无理重根及 18 组最高八次/正负/幅值缩放的 native 查询通过。完整 benchmark 原 10 项失败和输出保持。初始失败日志保留，临时直接算术测试文件已逐字节恢复；最终以真实 mapped query 验证新增算术路径。见 [报告](I12r5_Native_Q4_Tangent_Phase_Report.md)。
