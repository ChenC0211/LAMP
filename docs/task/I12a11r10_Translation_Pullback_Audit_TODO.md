# I12a11r10 / 平移 predictor 拉回与重构审计

## 目标及允许修改

只解决一个事实：r9 中第三步 predictor 曲率增长主要继承自旧 accepted 场，还是来自新 Q4 重构；并区分 compact mixed closure 的影响。0 production + 1 test：`tests/interface_hermite_gals_benchmark.cpp`；另含本 TODO/报告、README、`runs/baselines/I12a11r10_20261001/`。禁止改变生产模块、M、B/J、bubble 方程/系数、混合导数生产闭合、解析参考或原验收门槛。

## 数学约束与实现要求

新增显式只读入口，沿原 XY uniform n=32、dt=.01、四步平移推进。在每个真实 predictor 的 native interface quadrature 节点，使用同一 owned velocity/characteristic 与实际旧 accepted snapshot，核验 B=x-uΔt、J=I，比较旧场在 B 的二阶 jet 与新 predictor 在该节点的二阶 jet。两者与同位置的解析参考分别比较；逐点 signed curvature 差应闭合，不能把不同位置的最大误差直接相减作分解。

仅在私有测试 patch 中，将 master vertex mixed derivative 替换为旧场在实际 B 的 dxy，其余实际 transported first jets、edge 数据和 predicted center 保持；边界 hit 保留原 mixed。该反事实只在同一组 predictor 节点求值，不发布、不推进、不改生产的 I-R2 compact closure，也不声称是可用替代算法。记录被保留的边界数、master mixed 改变量及反事实剩余重构误差。

记录实际闭合的所有 corrected cell 的系数/体积残差/灵敏度及 bounds，确定全局最大系数修正的身份。候选与旧状态的 identity/M 保持，诊断异常单独记录且继续原推进。

## 验收与失败规则

T0：Release 无 warning；11 项直接回归与 r9 输出逐字一致；其它 27 可执行、全部 production、81 冻结/参考文件、3 PDF、r9 archive 不变。

T1：原非多项式四步的最终 BENCHMARK 行逐字等同 r9，8 份 predictor 分解与 8 份闭合记录完整、诊断异常为零；补充 Q2 常平移控制（n=8、一步），同样的逐点拉回/重构误差应接近舍入（1e-10）。检查实际 transport 首量与旧场拉回一致、B/J 关系，以及逐点误差分解。保存最大误差点的全部 jets/source-cell 候选身份与数据字段表。

不满足不变性、覆盖或参考检查则 FAIL；只读定位 PASS 不代表 I12 数值 PASS。遇到不可用外部依赖则记录具体阻塞；数值不收敛照实 FAIL，不增加容限或绕开条件。交付源码、原始日志/耗时、可执行与冻结哈希、分析报告及 archive；依据证据另立最小修复任务。
