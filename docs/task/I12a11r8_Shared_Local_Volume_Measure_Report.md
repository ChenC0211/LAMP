# I12a11r8 / 闭合与扫掠共享局部体积测量报告

**FAIL：12/12 配置完成推进，但两个 XY 组未通过空间收敛；I12 仍 FAIL。** 对应 [TODO](I12a11r8_Shared_Local_Volume_Measure_TODO.md)。两个原 XY 拒绝已消失，但梯度/曲率精度问题尚未解决。

## 已实现

4 production + 2 tests。把既有 PhaseRegion 的单叶矩形积分提取到私有 `HermiteGalsRegionMeasure.h`；与 I12a6 archive 比较，除 header 所需 inline 修饰外积分主体逐字一致。保留局部原点、retained cuts、两相直接积分、较小相补量、求积及 256 次细化工作上限。I07/I10/I11 的完整 native Q4 discovery 对 cut cell 采用相同 helper，原 sensitivity、拓扑证书、M、B/J、root 工作/验收预算与位移上限均不变。

全局 fixed/free moments 仍按共同原点计算，保留共享面账本；局部体积和其 moment estimate 单独写入消费字段。两条代数等价路径的舍入结果允许不同，不把局部体积伪写成全局 free moment。legacy evaluateVolumes 未启用此变化。接口语义见 [I07 API](../api/I07_Hermite_Gals_Phase_Volume_API.md) 和 [I10 API](../api/I10_Hermite_Gals_Bubble_API.md)。

## 已取得的验证

首轮构建 56.22 秒，无 warning。区域测试 61.44 秒 PASS：130 个 XY/RZ × uniform/AMR、circle/ring、液气互换的 cut cell，其 I07 rawLiquid 与同一 scalar 的 whole-cell PhaseRegion 液量逐位相等；原 212 个解析/加性/owner 案例和 11 个拒绝检查保持通过，M 不变。

两组原 r7 失败配置全部到达 t=.04，耗时 196.20 秒：

| 配置 | 位置误差 | 梯度误差 | 曲率误差 | ledger / reflux |
|---|---:|---:|---:|---:|
| XY uniform n=32 | 7.8330e-8 | 2.6939e-5 | 1.0386e-2 | 1.3878e-17 / 0 |
| XY AMR n=16 | 7.8023e-8 | 2.6833e-5 | 1.0242e-2 | 0 / 0 |

原 cell 531/341 的纯相拒绝未重现；上述曲率误差仍大，推进完成不代表空间收敛。由于保留全局账本并额外计算局部体积，本实现增加计算成本；bubble 回归本轮约 119.30 秒，r6 首轮约 75.01 秒。尚未做局部缓存或免重复求积优化。

## 回归与冻结审计

完整 Release 重建 89.71 秒，无 warning；30/30 CTest 通过，294.66 秒。23 项完整输出逐字不变；区域测试除新增 130 项比较汇总外原输出逐字不变。其余变化逐字段审计：I07 的 2809 个 cell 保持身份、capacity 与解析参考；bubble 的场景/迭代计数不变（Newton/bisection/endpoint 仍 503/128/0）；I11 仅体积相关列改变，几何导数不变；48 组 crossing 仅库存/通量误差列改变，first jets 与曲率不变。smoke 的 RZ AMR 平移曲率误差从 1.5898e-11 降至 6.2172e-14，其余场景曲率不变，体积列随新测量变化。

源码范围 4 production + 2 tests，81 项冻结/参考文件、3 PDF 和 r7 archive 不变。架构为 12 targets、380 translation units、901 files、3092 includes、5 legacy edges。旧 midpoint crossing probe 仍退出 1 且输出不变。

最细 XY AMR n=32 四步诊断完成，440.66 秒、退出 0；无库存越界、闭合拒绝或诊断异常，ledger/reflux 为零。最终位置/梯度/曲率误差为 2.1657e-9、1.8602e-6、1.1500e-3。cell 1408 的 predictor/closed 两种消费测量在首步 midpoint/end 均逐位相等，r6 首步 end 的 2.51 个预算差额消失；闭合后的端点体积残差为零。历史 `WORKING_MEASURE` 字段名中的 global 仍指原 I10 discovery 入口，本轮该入口已经采用局部路径，不应再将该列解读为零原点 free+fixed 账本。

## 长平移验证

四组静止（各 4 步）及四组平移（各 40 步、t=.4）全部通过，250.64 秒、退出 0。最大 ledger 残差为 1.3878e-17，reflux 为零。XY uniform/AMR、RZ uniform/AMR 四组运动的最终曲率误差分别为 7.9048e-14、1.4033e-13、5.6843e-14、6.1284e-14；完整数据见 `full_precheck.json`。这些多项式场景的通过不代替非多项式空间收敛。

## 已确认的空间收敛失败

XY uniform 三个级别均完成四步，但不满足逐级下降：

| n | 梯度误差 | 曲率误差 |
|---:|---:|---:|
| 8 | 9.2154e-6 | 7.2089e-4 |
| 16 | 1.9760e-5 | 2.9548e-3 |
| 32 | 2.6939e-5 | 1.0386e-2 |

最后一级梯度/曲率阶为 −0.4471 / −1.8135，程序报告 ORDER_FAILURE。XY AMR 的曲率误差为 2.9628e-3、1.0242e-2、1.1500e-3；最后一级梯度/曲率阶虽为 3.8505/3.1547，但 n=8→16 误差增长，仍不满足原逐级下降条件。两个独立 XY 案例及最细 AMR 的最终 BENCHMARK 行均与矩阵逐字一致。

RZ uniform 曲率误差为 1.0264e-2、9.8473e-4、8.1842e-5，最后一级梯度/曲率阶 5.0329/3.5888；RZ AMR 对应 1.0263e-2、9.8472e-4、8.1842e-5，阶数 5.9057/3.5888。两组 RZ 均满足下降及阶数要求。完整矩阵耗时 1324.21 秒，12/12 推进完成，0 次推进拒绝、2 次 order failure，退出 1。

## 初始静态对照

为区分初始表示与推进放大，追加运行原 `--curvature` 入口，使用相同 shape/grid 的 t=0 几何，不改源码或求积参数。12/12 初始配置均满足逐级下降；最后一级梯度/曲率阶分别为 XY uniform 3.9582/3.3036、XY AMR 4.0945/3.4125、RZ uniform 4.0149/2.9696、RZ AMR 4.0816/2.9696。耗时 10.21 秒，与运动矩阵并行运行，因此本次完整矩阵耗时不作为严格的串行性能对比。

XY uniform n=32 的最大曲率误差从初始 9.7002e-6 增至四步后的 1.0386e-2（约 1071 倍），梯度误差约放大 890 倍；XY AMR n=32 的曲率误差约放大 243 倍。该对照将问题定位到推进过程，尚不能单独归因于某次 bubble 或 skeleton 更新。原始数据与逐配置比值见 `static_translation_comparison.json`。

## 结论与交付

共享局部体积解决了已证实的两种同源测量差额及两个纯相拒绝，所有推进配置完成，但未解决运动曲率收敛。下一项只读追踪 XY uniform n=32 四步中的 accepted、predictor 和 bubble-closed 几何，保留同一 M、B/J、全部求根证书与门槛，以定位误差首次明显放大的阶段。

原始日志、静态/运动逐配置对照、字段级回归审计、源码/runtime 哈希、交付清单及 archive 位于 `runs/baselines/I12a11r8_20261001/`。本项数值 FAIL 与归档完整性 PASS 分开记录；不进入 Hydro。
