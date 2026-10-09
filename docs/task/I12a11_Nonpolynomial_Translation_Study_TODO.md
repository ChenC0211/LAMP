# I12a11 / 非多项式界面完整平移的空间收敛

## 目标

只检验一个事实：已通过多项式 circle/sphere 的区间通量完整推进，对一般 Q4 数据是否仍保持可达库存、准确梯度/曲率和空间收敛。使用原 `exp(q)-1` 初始 profile（q 为 circle/sphere 二次式），零面与解析圆/球一致，之后完全由旧 authority/B/J 推进。

## 允许与禁止

0 production，最多 1 test：既有 benchmark；新增显式 opt-in CLI，不修改默认历史 benchmark。另含本 TODO/报告、README、`runs/baselines/I12a11_20261001/`。不改 M、B/J、任何 gate、closure/时空积分公式或解析参考；不替换 interior 为真解，不放宽失败期望，不因常速度时间误差处于舍入底限伪造二阶。

## 实现要求与数学约束

- 复用 advanceTranslation 和原 runCaseUsing。允许 runCaseUsing 返回其已计算的最终 Errors（失败为空），不重复构造几何或误差，也不改变原模式输出/行为。
- XY/RZ × uniform/AMR，n=8/16/32，固定 dt=.01、4 步至 t=.04。固定物理速度 .2；RZ 半径沿原配置 .27，其余 .19。保留实际粗细切割覆盖和 accepted-state 失败隔离。
- 打印全部 12 组的状态、最终误差、M ledger/reflux 与 h 方向阶数。发生失败仍运行其它配置；不能从已接受中间时间的数据计算最终收敛阶。
- 本项只验空间收敛，不代替独立非恒定速度时间阶或长时 vortex 验收。

## 验收

T0：Release 无 warning，architecture/直接前序回归输出不变；production、其它 27 可执行、81 冻结文件、3 PDF、前序 archive 哈希保持。

T2/T3：12 组全部完成目标时间；逐单元容量和原 I10/I11 guard 不变，ledger/reflux <=1e-10。每种配置梯度/曲率误差随 h 递减，末档阶数分别 >=2 / >=1.5。明确记录 position、phase-volume 与逐单元 M 误差，不把守恒 residual 当解误差。

## 失败规则与交付

任一推进或阶数条件不满足即 FAIL，保留真实 accepted time、stage 和异常，不进入生产修复；另开一个有证据支持的定位任务。报告区分“解析参考”“真实旧 Q4 平移”和“新 first-jet/Q4 重构”，不能预先断言失败来自时间求积。交付源码/日志/报告/runtime/hash archive；I12 继续保持 FAIL，按用户授权持续推进。
