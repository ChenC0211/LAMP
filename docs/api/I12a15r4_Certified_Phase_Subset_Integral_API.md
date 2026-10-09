# I12a15r4 / 认证相子区域积分 API

见 [TODO](../task/I12a15r4_Certified_Phase_Subset_Integral_TODO.md)、[设计](../task/I12_Certified_Phase_Subset_Integral_Design.md)、[报告](../task/I12a15r4_Certified_Phase_Subset_Integral_Report.md)。本项将 r2/r3 接入原 a15 积分，未改变其预算、来源契约或输出类型。

`PolynomialVelocitySlab::encloseNormalFaceIntegralWithAffinePhase(a,b,begin,end,radial,forms,wlLower,wlUpper,wgLower,wgUpper)` 是 Core 数学原语。默认权重为 1/1/0/0。所有 forms 必须使用传入 face/time 的同一 `q_s,q_t∈[-1,1]`，并覆盖每个实际参数点可能使用的全部 donor scalar；Core 本身不持有 mesh/snapshot，也不验证这些来源，不能用任意 caller forms 宣称已认证相。

权重范围要求 `0<=wgLower<=wgUpper<=1`、`0<=wlLower<=wlUpper<=1` 且两个 gas 端点分别不大于相应 liquid 端点。相同 donor 中原关系 `wg<=wl` 的全局 extrema 满足这些约束。空 forms、非有限/反向 form 区间、非法 face/time/weights、FP 环境或无法表示的精确算术均明确拒绝。

返回仍为 `PolynomialFaceIntegralBounds` 的四个 gross 界。原总正/负通量保持原定义；用所有 scalar forms 的共同 band 确定严格负与严格正半平面，再以精确有理面积及第一径向矩得到物理测度。子区域贡献采用点态 gross 速度的非负下界；反向速度导致下界为零时，原分块法继续控制误差，不用净均值抵消 gross。最终由同一四宽和 `errorBound()` 控制验收。

Interface 的 `integrateHermiteGalsPolynomialFace` 公共签名、配置及结果不变。原完整 branch、closed donor 覆盖与 inventory owner 检查后，仅对相未决 tile 计算原 snapshot 的 Q4 affine forms，按原 donor 权重取 extrema，调用上述 Core 原语并与已有 tile 证书相交。原 pending、细分、最大 panels/depth、来源身份、integrated transfer/reflux 职责全部保留。

`Detail/ExactFaceTimeRational.h` 只提取原 Snapshot 的 rational 运算和 unit-square 半平面矩，供旧常速平面入口与新证书共同使用。旧入口仍用原 `faceTimeValue` 表示结果；新入口先组合精确物理矩再向外包围，因此没有把旧浮点输出误当严格区间。恒零 form 的严格相子区域为空，其余非恒定直线的等号边界为零测集。

本项不推进独立库存 M，不修改几何，也不自动接通完整旋转/涡旋 driver 或 Hydro；I12 的其余验收仍须完成。
