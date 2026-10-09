# I12a11r9 / 平移曲率逐阶段审计

## 目标及范围

只定位一个事实：r8 的 XY uniform n=32、dt=.01、四步推进中，梯度/曲率误差首次显著放大是在 predictor，还是实际 bubble closure。0 production + 1 test：`tests/interface_hermite_gals_benchmark.cpp`；另含本 TODO/报告、README 和 `runs/baselines/I12a11r9_20261001/`。不改 M、B/J、mixed closure、bubble 方程/求根、全部误差预算、拓扑证书或解析参考。

## 实现与数学约束

新增显式 `--translation-curvature-stage-audit`，复用同一原四步 driver。只读记录四次 accepted 起点以及八个 midpoint/end 的 predictor、closed scalar。复用已有独立解析曲率与 cardinal bubble value/gradient/Hessian 检查，记录真实最大曲率误差的 cell/point。

逐闭合记录最差单元的 actual root report：是否激活、容量、c 修正、实际初始/最终体积残差、最终 dV/dc、体积/系数割线斜率、是否位于证书端点，以及该批次迭代、Newton/二分/端点次数。系数除以 h² 仅为量纲参考，不宣称曲率误差界。只读 private scalar 不发布、不用于推进；M 与 actual identity/绑定状态必须保持。

不重复 r8 已确认的整单元双测量，只记录实际报告。原工作精度诊断默认仍只测首步，不改变旧 CLI 输出。诊断异常独立计数，不能遮蔽真实推进结果。

## 验收及失败规则

T0：Release 无 warning；11 项直接回归输出逐字等同 r8；其它 27 可执行、全部 production、81 冻结/参考文件、3 PDF 与 r8 archive 不变。

T1：完成原四步，最终 BENCHMARK 行与 r8 同案例逐字一致；4 个 accepted 起点、8 次实际 closure、20 条几何阶段记录、8 条 root stats，诊断异常为零。给出首次误差增长与后续阶段传播的证据，不从全程总误差直接推断某个模块为根因。

若证据不完整或诊断扰动 actual solution，本任务 FAIL；只读诊断 PASS 仍不改变 I12 的数值 FAIL。交付原始数据、逐阶段对照、源码/runtime 哈希、未解决问题和 archive，之后再按证据写最小修复任务。
