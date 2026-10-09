# I12a11r2 / 闭合工作精度与库存预算衔接

## 目标

只修复 I12a11r1 已量化的事实：I10 以 1e-10 容量残差停止后，剩余差额可超过下一步库存允许的 64 eps，几何液量流空时造成负库存。使同一 bubble 方程求解到与消费预算相容的工作精度；不改变 M 或通量。

## 范围

最多 4 production：新增 `Detail/HermiteGalsVolumeRoundoff.h`（唯一预算定义）、原 `Detail/HermiteGalsInventoryAssembly.h`、HermiteGalsBubble.h/.cpp。最多 2 tests：bubble contract 和 benchmark。另含 TODO/报告、bubble API/README、`runs/baselines/I12a11r2_20261001/`。不改 CMake 数值目标或冻结方法。

## 数学约束与实现

- 原 64 eps 库存预算原值、原运算顺序移至共享私有 header，I08 行为保持。I10 的相对工作残差取 `min(用户 relativeVolumeTolerance, 共享预算/4)`，为后续体积/通量数值误差留下余量；这不是对全部求积误差的严格界证明。
- 同一工作精度用于初始跳过、纯相检查、括区判定、Newton/二分终止和最终核验；不只改最后 Report。Report 记录实际工作精度，算法不读取 Report。
- 原 scalar/bubble、M、first jets、B/J、正则性/位移/拓扑证书、64 次 root 工作限制及 moment gate 保持。用户的 1e-10 仍是请求上限，不提高任何门槛。
- 原允许的 1e-12 纯相尾量超过库存预算，新增拒绝验证；仍验证预算内小尾量被如实报告、M 不变，不制造几何分量或裁剪。
- 新诊断模式随修复改为验证原 4 步配置真实完成，不以“成功观测拒绝”替代推进成功。历史失败源码/日志保存在原 archive。

## 验收

T0：Release 无 warning、预算/冻结/前序 archive 核验；architecture 解释新 header/include。完整前序 30 项 CTest，逐项分析因更严格 root 精度产生的数值/迭代变化，未受影响算子输出应保持。

T1：原 XY/RZ、uniform/AMR、directed bracket 的 bubble 制造解达到预算内残差，first jets/M/约束/拒绝隔离保持；纯相预算内/外对照。严格精度导致的停滞、证书或求积失败不能掩盖。

T2：原 XY AMR n=32、dt=.01 的 4 步到达 t=.04；旧库存差额及同源扫掠诊断不再触发越界。四组 40 步多项式平移及静止重跑。非多项式 12 组矩阵重跑，库存与曲率条件分别记录；若曲率问题独立保留，本项库存精度验收和父任务 I12a11/I12 状态分开报告。

## 失败规则与交付

任一新停滞、精度不可达、未解释前序回归或原库存拒绝仍在则本项 FAIL，先保留证据再拆根因；不得增加预算或放宽 geometry gate。交付源码/API/报告/原始日志/runtime/hash archive。即使本项修复通过，I12 其余空间/时间、旋转/vortex 等未完成条件仍继续推进，不进入 Hydro。
