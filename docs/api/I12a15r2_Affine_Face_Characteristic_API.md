# I12a15r2 / 面—时间仿射 characteristic API

见 [TODO](../task/I12a15r2_Affine_Face_Characteristic_TODO.md)、[设计](../task/I12_Affine_Face_Characteristic_Design.md)、[报告](../task/I12a15r2_Affine_Face_Characteristic_Report.md)。本项新增 Core 包围原语，尚未接入 a15 相积分；原剪切预算失败和 I12 / Gate I 不因本 API 通过而改变。

`PolynomialVelocitySlab::encloseAffineFaceRk2(domain, a, b, first, last)` 要求有限、升序、正长度的轴对齐线段及 slab 内正时间区间。Core 不持有 mesh；Interface 消费者以后须从原 native face/interval 提供输入。入口先调用原 `encloseInteriorRk2TightTime`，只有同一区域的完整单步、无 hit 分支认证后才构造新 form；不支持的 CFL、边界、轴、opaque coverage 按原规则拒绝。

返回 `PolynomialAffineRk2FaceBounds` 是独立 Core 值记录，不借用 mesh/geometry/velocity 的指针，不是 Interface provenance capability。包含 faceBegin/faceEnd、timeBegin/timeEnd、原 branch witness 和两个 departure forms。公开可复制字段仅用于数学原语；调用者自行改字段不能建立新的认证来源。

每个 `PolynomialAffineCoordinate` 满足 pointwise 包含关系：

`B_k(q_s,q_t) ∈ center + tangent*q_s + time*q_t + remainder`

其中 q_s、q_t∈[-1,1]，由原 double 端点解释为精确实数后的仿射参数定义；两个分量共用同一对 q。center/tangent/time 是有限 double，remainder 是严格 outward 上下界，包含非线性项、系数代表值转换误差以及原浮点求值/分支舍入。time 系数对应归一化 q_t，不是物理 dB/dt；本 API 也不以浮点映射可微为前提。

`proof.encloseDeparture(point,time)` 检查物理查询位于保留的线段/时间中，用 outward 端点除法反算参数，再评价两个 forms，并与原 branch departure box 相交。参数计算避免 subnormal 分母的倒数溢出；不是简单用整个 B box 代替局部查询。非法、非有限、非 nearest 或 FTZ/DAZ 环境拒绝。

原 a13/a15 包围入口、实际速度和 RK2 tracer 保持数值不变；M、geometry、phase、库存权重、reflux 与 driver 均不由本 API 修改。此处没有相分类、gross flux 或误差预算验收；后续须把同一权威 Q4 施加于这些 forms，并对明确液/气子区域计量。