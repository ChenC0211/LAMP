# I12a11r13 / 非多项式平移步长控制试验

## 目标及范围

只检验现有方法路线中的一个方向：固定 XY uniform n=32 和终点 t=.04，减小实际时间步能否减轻 r11/r12 已证实的弱灵敏度闭合及后续曲率放大。0 production + 1 test：`tests/interface_hermite_gals_benchmark.cpp`；另含本 TODO/报告、README 和 `runs/baselines/I12a11r13_20261001/`。禁止改 M/flux/B/J/I06/I10/I11、体积/拓扑/位移/工作预算、旧 benchmark 配置或原收敛门槛。

## 实现与数学约束

新增两个显式控制入口，使用同一初始非多项式场、同一 actual constant-translation driver，分别 dt=.005 / 8 步和 dt=.0025 / 16 步。原 dt=.01 / 4 步参考直接采用 r11 已归档原始 BENCHMARK 行，保留该失败事实。各配置的实际步数与物理时间明确记录，不暗中细分旧 benchmark。

每个 accepted end 只读记录实际位置/梯度/曲率/体积误差与 ledger/reflux；每次实际 closure 记录最大系数修正及所属单元、体积残差和迭代统计，不另建几何参与推进。保存每个配置独立日志、真实退出码与耗时；一次拒绝不妨碍执行另一配置。

## 验收与失败规则

T0：最终 Release 无 warning；11 项直接回归完整输出与 r12 相同；其它 27 可执行、全部 production、81 冻结/参考文件、3 PDF 和 r12 archive 不变。

T1：两个控制均达到 t=.04，所有 accepted/closure 记录完整，ledger/reflux 原门槛保持；与原 .01 参考比较最终梯度、曲率误差。只有两者随 .01→.005→.0025 均严格下降，才将本项“减小步长是有效修复方向”判为 PASS；任何推进拒绝或误差不下降则本项数值 FAIL，照实归档。固定空间的此试验不宣称独立时间收敛阶或 I12 通过。

交付源码、三步长逐字段对照、实际阶段统计、测试/冻结/runtime 哈希及 archive。若 FAIL，不能据此提高容差或开启未声明 fallback；下一项依据原方法允许的空间分辨率/局部加密方向或具体设计变更推进。
