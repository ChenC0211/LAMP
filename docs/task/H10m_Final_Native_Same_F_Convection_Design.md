# H10m / 实施前约束与独立验收检查

2026-10-05。当前材料 owner 为 H10l2，唯一 F、真实 R、Pure FirstJet 独立及用户规定的四项材料硬约束保持。本项不改变 H10b 的任何参数或门槛。

对内段，实际兼容散射可写为 `N = (K + diag(q)/2) u`，K 为共享反对称 native pair；外段另外贡献显式 trace。载荷为 `-N + u R`。闭合周期端口后，常数再现要求 `K 1 = q/2`；几何一次矩不是一般非均匀材料的点 PDE 目标。周期重复 native 采样及各自真实 half-dual 质量必须一起处理。精确二次矩在旧闭合 one-ring 上已有不可相容反例，不能照搬成硬约束。

准备先复核当前 source/runtime、原 affine 及逐 native 矩阵，检查 H10b 各门槛是否依赖待修改的中心算子。原 donor 分支只使用固定的 F、质量、R 和研究周期配对，不使用中心算子的 weights/couplings/apply；因此，仅修改 CchKepConvection 不可能改变 donor 的终值或验收结果。

原 donor 门槛要求完整 T=1 的累计能量损失在 n16→32 时呈至少 0.7 阶。独立的均匀周期上风 Fourier 模式满足

`lambda_h = -(1-exp(-i*2*pi/n))/h`，`E(T)/E(0) = exp(2*T*Re(lambda_h))`。

采用研究的 RK4 和 `steps=32*n` 时，精确离散预测为 `|P4(dt*lambda_h)|^(2*steps)`，其中 `P4(z)=1+z+z²/2+z³/6+z⁴/24`。累计损失 `1-E(T)/E(0)` 在有限粗网格可饱和，不能直接将其任意两档比率当作渐近 O(h) 耗散率。此处是实施前可检验假设，尚非实际四组 gate 裁决。

在现有研究 driver 增加 `donor-audit`，直接复用原 Study/rhs/integrate、原四组 n8/16/32、全部 native DOF、T=1 和 32n 步。记录原 loss gate，另以均匀网格横向分量的 Fourier 多项式逐 DOF 交叉验证 donor。瞬时能量衰减率只作诊断，不替换原门槛，不进入状态或系数。所有原 Study/full/affine/rowAudit 保持。先明确这项与中心算子无关的必要验收条件；若固定 gate 已失败，保留 FAIL，不以修改中心算子、通量、质量或门槛掩盖，下一任务须明确处理门槛的数学定义。

新增内容只属于既有研究 driver/README，生产和短测试暂不修改；输出进入本任务 runs/checks。该检查不等于最终算子完成，不等于完整 H10b 重跑，也不能豁免当前空间阶/有界性问题。
