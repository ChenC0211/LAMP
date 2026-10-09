# H08a2r9 / 恢复迁移设计

2026-10-05。上一轮为 progress；H08a2r10 的373冻结payload、14当前文件与1382源码/输入已复核。原六反例关闭后恢复同一任务，production 0/test 2；原历史FAIL报告不改写。检查归 `runs/checks/H08a2r9/20261005_164956_929`。

沿用[原迁移清单](H08a2r9_Spatial_Source_Consumer_Contracts_Design.md)，具体参考如下：

1. 对实际 owner 的 pieces 建一次按cell的坐标分界缓存；正则暂存场使用同一 owner 的 supportKnots。只从几何来源取得分片，不引入另一速度/clock owner。
2. 物理GL3×GL3积分 `div u`（XY）或 `2π[r(∂r u_r+∂z u_z)+u_r]`（RZ），分别积到整cell及真实dual overlap。原场的加权散度密度与正则Q4场的密度在每片最高五次，故求积足阶。误差尺度用各物理项绝对贡献，不用相互抵消后的净量。
3. atom 的完整primitive与同一实际速度查询的独立GL3线积分比较；原解析族继续用闭式值/J/native中心核验。新旧来源应用相同数学参考与阈值；同源性检查不能替代解析精度。
4. 对 `(r*u_r,u_z)`（RZ）或 `(u_x,u_y)`（XY），用六个片内物理Chebyshev节点的Lagrange导数参考，覆盖四次速度/五次加权速度。在轴上使用 `∂r u_r=∂rr(r*u_r)/2`，不除以零；轴奇偶继续独立检查。原点、片端点、实际native面和原144场景覆盖不减。
5. 点态div与cell均值、dual边界与平均div配分之差仅诊断。piece/atom库存按真实数据自洽核验，旧固定存储常数与“跳跃必须非零”不再作正确性条件；原解析资格集合不缩小。
6. 正则未发布场通过同一求积/导数核验，并用符号反转、截断高项、遗漏支持分片和错源拒绝检查证明敏感性。所有扰动只存在测试，不能成为生产开关。保留原zero/tiny/wall、来源identity/full clock/AMR共享及原最终native/affine/H10b门槛。

完成后只为统一来源部署准备验收，尚不切换实际B/characteristic。若参考独立性、原精度或覆盖无法保持，即FAIL，不删断言掩盖；H08a和H11不因本项通过而完成。

执行结果：PASS，见[报告](H08a2r9_Spatial_Source_Consumer_Contracts_Resume_Report.md)。原208基、64布局/128输入、144 jet场景保持；48布局432输入的新旧来源应用同一参考及阈值通过，三项故意扰动被检测，9/9回归。首轮将旧cellPieces侧片段计数误比为矩形piece数量的假设已按实际库存纠正，未改数值门槛或生产。
