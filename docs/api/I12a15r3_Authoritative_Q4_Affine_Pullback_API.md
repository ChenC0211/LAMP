# I12a15r3 / 权威 Q4 的仿射拉回 API

见 [TODO](../task/I12a15r3_Authoritative_Q4_Affine_Pullback_TODO.md)、[设计](../task/I12_Authoritative_Q4_Affine_Pullback_Design.md)、[报告](../task/I12a15r3_Authoritative_Q4_Affine_Pullback_Report.md)。本接口只包围同一 scalar；尚未接入 a15 相积分，不能据此宣布原剪切或 I12 通过。

`Snapshot::nativeQuarticAffineEnclosure(index, physicalCoordinates)` 要求原 snapshot 是 active native Q4，指定原叶的 native index。调用处需同时包含 `LevelSetHermiteGeometry.h` 与 `PolynomialVelocitySlab.h`。两个物理坐标采用 r2 的 `PolynomialAffineCoordinate`，共用 `q_s,q_t∈[-1,1]`；返回一个同类型的 scalar form：

`phi_cell(B(q_s,q_t)) ∈ center + tangent*q_s + time*q_t + remainder`。

其来源是原叶的精确齐次 Hermite 数据，包括共享/悬挂约束及原 exact refinement 的限制数据。先做精确 Bernstein→power 基变换，再对系数商向外包围；不取 `nativeQuarticLocalDofs` 的 rounded 诊断缓存，不创建第二几何、不改变源数据。实际 cell 宽度也由真实端点的实数差向外包围。

**叶外结果只是该叶多项式的数学延拓。** 接口有意允许 forms 的 box 越过该叶，供未来完整 closed donor 覆盖使用；只有真实 B 属于该叶的参数点，才能将其符号解释为实际 phase。该接口不选择 donor、不提供 physical boundary 历史，也不验证不同参数点的完整覆盖；这些仍由 Interface 消费者负责。跨叶可分别使用其子域证书，或采用对所有候选叶均成立的共同符号区域，不能平均系数重新拟合界面。

`PolynomialAffineCoordinate::enclose(qsRange, qtRange)` 接受 [-1,1] 内有限有序的闭参数区间，返回 scalar 的向外范围。这个入口只做参数域求值，不接受物理坐标，也不与原 branch box 相交。物理回溯查询仍使用 r2 的 `proof.encloseDeparture(point,time)`。

Core 通用原语 `encloseTensorQuarticAffine(coefficients, cellLower, cellUpper, physicalCoordinates)` 使用 25 对有限、有序的 chart 幂系数上下界，索引 `5*y+x`，要求两个正 cell 宽度。它没有 mesh/geometry provenance，不能把任意传入系数当作原相场；Snapshot 入口才使用原齐次数据构造这组界限。

无效 index、inactive/nonquartic/moved source、非有限或反向区间、不可表示的商/中间范围、非 nearest 舍入及 FTZ/DAZ 环境拒绝。没有小系数截断、容差拟合或放宽 64eps 舍入修复。M、库存权重、共享通量与几何发布保持原职责。
