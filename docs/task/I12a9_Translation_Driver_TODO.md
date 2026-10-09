# I12a9 / 同源常速度半步与整步驱动

## 目标：一个事实

常速度完整推进的 half/full M 必须消费同一旧 scalar、速度、characteristic 与 BC owner 的实际区间积分；中点几何继续真实发布，不能作为相通量单点时间求积。I12a8 已验收归档，本项验证实际 first-jet → I06 → I10 → I11 链路。

## 允许修改与禁止修改

最多 2 production 文件：HermiteGalsTemporalPhaseFlux.h/.cpp；最多 2 测试相关文件：interface_hermite_gals_benchmark.cpp、tests/CMakeLists.txt。文档及 `runs/baselines/I12a9_20261001` 不计生产预算。禁止改 B/J、reconstruction、bubble 证书、库存容限、冻结路径、Hydro 和 PDF。旧 midpoint driver 保留为明确的历史对照，不做自动 fallback。

## 数学约束与最小实现

1. owned constant velocity 增加 actual candidate/time 的 arrival view；保持原不可变 provider 生命期，只允许原真实时间范围内的同一步正向时间。opaque provider 不具备此能力。不复制另一套速度插值，不改固定 RK2。
2. 整步 prediction 建立真实速度和 characteristic 来源；半步 prediction 从同一 source 派生实际中点 view。interior vertices/edges/center 仅消费旧 snapshot 经原 B/J 拉回；解析 first jet 仅用于实际入流外边界，解析误差参考不能写入 interior。
3. 一次准备旧 phase-region；显式声明完整入流面的常气相历史，并检验制造解在整个时间段的入流面保持气相。half/full interval 共享该 history 与速度/特征 owner。实际所有面经 I12a8 → I12a3 → I09；M 只来自最终共享通量。
4. 中点和步端均经原 I10/I11。失败不接受 candidate，不修改原 accepted identity/M。没有 quarter 几何或 midpoint 求积；也不通过减小 dt 对齐事件。

## 验收

- T0：Release 主程序与受影响合同，architecture，81 冻结/参考 + 3 PDF，I12a8 archive 和文件预算；旧合同及 crossing 对照的完整输出检查。
- T1：派生 view 的真实身份/时间、B/J、生命周期，非法时间/opaque/moved-from 拒绝；完整 half/full 同源身份及最终 M 消费检查。
- T2：显式 `--translation-sweep-smoke` 跑 XY/RZ、uniform/AMR 的静止四步及平移首步；`--translation-sweep` 跑原 n=8、dt=.01、40 步的四个完整平移算例，以及静止四步。原 stationary drift、ledger/reflux 和闭合 gate 原样保留，记录位置/梯度/曲率/逐单元库存误差。
- 新 mode 不能替代 rotation/vortex/收敛/其它未完成 I12 验收；I12 仅在全部条件实际通过后改状态。

## 失败规则与交付

任何原数值 gate 失败即记录真实 stage/cell/accepted time 和旧状态不变证明；本项记 FAIL，并在新 TODO 中只修复已证实的单一缺陷。若必须更改冻结数学路线且现有授权不覆盖，记 BLOCKED 并提交证据。交付代码、可复现运行日志/哈希归档、结果报告与未解决问题。

## 执行结果

FAIL：26 项回归通过、25 项旧输出不变；RZ 两组完成 40 步，XY 两组在第五步的横向根包围归属检查拒绝。已保留失败与旧状态不变证据，见 [报告](I12a9_Translation_Driver_Report.md)。继续单项修复，不降低数值 gate。
