# H10j0 / 实际 compatible advector

2026-10-04，**PASS / engineering PASS**。本任务关闭可复用、与 phi/密度无关的 actual native advector；没有实现 phase/material F、最终 native convection 或 H14 时间推进。长期 goal 保持 active，本轮为 progress；下一唯一任务 [H10j](H10j_Spatial_Phase_Rate_Owner_TODO.md)。[固定任务书](H10j0_Spatial_Advector_Reconstruction_TODO.md)、[实施前设计](H10j0_Spatial_Advector_Reconstruction_Design.md)、[实际 API](../api/H10j0_Cch_Spatial_Advector.md)。

## 改动及唯一 owner

生产固定3文件：新 `CchSpatialAdvector.h/.cpp` 与现 CMake 显式 CCH 源登记；仅扩展既有 `cch_momentum_geometry_contract.cpp`。没有新 target/lib/tools。新 cpp 在 `lamp_hydrodynamics_cch` 编译一次，从 legacy Hydro glob 排除；不进入 Core/Interface。

`CchAdvectorReconstruction::prepare` 持有原 momentum geometry/layout/D_I，准备实际选中 closed one-ring 的 physical-mean 三列 SVD 最小范数 tangent weights、cell side pieces 与 canonical atoms。`evaluate` 复制 actual native values，只更新系数，生成私有 fresh advector identity，保留 caller geometry/full clock。两个 actual update 和 copied preparation 的 weights/atoms 是同一 immutable storage，没有再次分解或重建 registry。声明来源不是 H14 commit 证明。

共享 face trace 保持 A_f*u_f；实际 cell 使用体积坐标 transfinite 加显式零 normal lift。对边 piece knots 参与积分与 overlay，返回场的局部 divergence 等于原 D_I/V。canonical atom 不含 component/material key，来源恰为原 primalVelocity 或 ownerPressure；所有原 primal/dual 段映射为 ordered atom views。内部反号直接来自同一值，AMR 只求实际子段和。

两端分别保存 endpoint-relative normal polynomial，保留 independently supplied offset/width。RZ radial axis 固定零，不除 r；native point 的 transverse 分量按真实 half-dual overlap 容量权重读取邻接 cell B。只有运动学总量相容，不把材料 commuting 或 quarter 加回硬约束。

## 数值证据

标准 local-benchmarks Release，OMP=1。必要 CTest **7/7通过**：architecture/self-test/layout 与现 phase geometry/momentum geometry/dual material/KEP。现 momentum 契约保留原全部检查，并增加64布局/128更新：XY/RZ × uniform/实际2:1 × n4/8/16/32 × full/三类 selected。

| 核验 | 实测 |
| --- | ---: |
| actual canonical atoms / 两更新原子表 | 114,300 / 228,600行 |
| 测试局部点参考 / 实际 dual view 使用 | 72,080 / 247,560 |
| 研究 actual face / cell 记录 | 77,540 / 36,040行 |
| 契约 face 物理均值误差 | 1.8464e−16 |
| 可识别 cell/native affine 再现误差 | 2.2205e−16 |
| 独立非仿射 lift/atom 积分误差 | 2.9184e−16 |
| 原 D_I 与实际返回场局部散度误差 | 2.3641e−15 |
| CSV 独立 atom coefficient 误差 | 6.2811e−14 |
| CSV 独立 closed integral 误差 | 3.1559e−16 |

短契约参考从实际 native faces、选中 cell 邻接和 GL3 物理面均值取得；非仿射参考独立做 N×3 最小二乘 fit，随后按真实 knots 的独立 GL3 积分构造局部 lift，不读 owner slopes。cell 返回场用同一 piece 内 volume-coordinate 中心割线验证实际 divergence；实际 dual 边界 atom 积分加回原 overlap 的 divergence 积分。轴、normal反向、真实 subinterval/ordered views、所有原 native 行均保留。

64布局仍有原72 rank deficient /56 tangent unidentified，全部未识别行仅 n4 selected2/3；full及n>=8 selected为0。这些极粗行仍检查有限值、总量/divergence、sharing、来源与复用；只保留原明示的精度缺口，没有 fallback 或删行。

既有 research driver 新增独立只读 `advector` 模式，实际生成16 geometry revision /128 advector revision，保存 full midpoint clock。独立 PowerShell 后处理从实际坐标及原 affine 数据重算容量、速度、系数和积分；不使用 expected 列。其最大 original divergence 绝对误差1.2213e−14，仍低于固定1e−12真实尺度。

RZ `u_z=b+c*z+d*r²` 的实际 input mean 与 owner 同点差异为 `.31*dr²/8`；uniform和AMR各四档加密最大差异为 2.421875e−3、6.0546875e−4、1.513671875e−4、3.7841796875e−5，比值约4，闭式残差≤4.89e−17。保留原 point 参考，不冒称一般 RZ point/physical mean 同时 exact。

## 来源、复用及拒绝

两次 evaluate 保留原 input/geometry storage，fresh identities 不同；copied owner 与原 handles 释放后的 owner 均可用。wrong layout、fresh momentum instance、wrong geometry/full clock、真实 mesh 变更/regrid、invalid clock、nonfinite/overflow、axis 冲突、rounding/FTZ、越界点与非法 shifted range 明确拒绝，失败后原 owner 系数/identity 保持。Initial/Midpoint/End 的合法声明全部检查。width=1e−25 即使无法表示新的 global endpoint，正 shifted integral 仍保持。

零 native 输入的所有实际 polynomial 与积分严格为零；另取 XY/RZ 实际零面均值输入，仍保留非恒定、变号的 normal polynomial 与非零 gross，不由净率零推断局部速度恒零。

报告逐布局保存原 face/segment、SVD/weights、piece/atom/view数量及实际准备/更新耗时，验证局部库存界。全部64布局共392,134 weights，准备耗时合计0.1206s，两次 evaluate合计0.0602s（本机一次观察，不是计时性能门槛）。最大 n32 AMR为3,680 faces/11,072 atoms；准备约11.4–11.7ms、两次更新约4.46–4.47ms。geometry复制共享storage；一次 line sorting仍 O(N log N)，局部存储及系数更新为线性。

首轮 guard 因 expected rejection 文本错误停止；补测 secant 最初在 RZ axis cell 取到了负体积坐标，测试据实际 knots/边界选择 piece 内步长后通过。失败日志保存，未修改生产方法或数值门槛。

## 原路径与退出边界

H10i2的188 payload/current canonical 在开始前独立复核；implementation-before design 与 source/input目录已保存。删除新 driver 的三段 observer block 后原源码全文保持。原 source-chain 与 affine研究重新运行，23份保存 CSV字节相同；旧完整 ordinary source 证据保留前置archive，本轮未重跑 segments/full H10b。Core/FirstJet/Q4/transport/clock、H07–H10原生产文件、F/R/uR/KEP/D_I/C、App/PDF和字典保持。

证据冻结于 [H10j0 baseline](../../runs/baselines/H10j0_20261004/README.md)，包含 actual运行/exe/DLL、compiled source/短测试、独立验证、失败日志和 flat manifest。最终 docs architecture/self-test/layout另行核验。

只关闭H10j0。四材料硬约束仍为唯一共享 F、内部反号/真实AMR子面和、来源一致、KEP/skew；quarter/commuting/R=0仅诊断。Pure FirstJet完全独立、无库存/bubble/质量修复。它是分片 Hdiv normal 场，未证明切向C1或已接Interface characteristic/Jacobian。下一H10j同源causal phase factory，再H10k/H10l和最终native PDE；内部affine先通过，再原H10b全部28窗口/84gate。原14失败保持，H10尚未关闭，H11/H08a等待。
