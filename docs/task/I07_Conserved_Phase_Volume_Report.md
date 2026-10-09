# I07 / 守恒单元相体积状态与初始化报告

**状态：PASS；2026-10-01。** 按 [TODO](I07_Conserved_Phase_Volume_TODO.md) 完成路线图 I07。下一项 I08；本轮未开始通量、reflux、bubble 或 production 几何发布。总体目标继续 active。

## 实现和所有权

新增 `HermiteGalsPhaseVolume.h/.cpp`，在 `InterfaceCellMeasure.h/.cpp` 开放完整拓扑发现后的既有 moments-only 积分。总计 **4 production + 2 test**，符合预算。共用 `InterfaceCellMeasure.cpp` 的全部旧函数保持不变，仅追加入口；旧算法没有修改。

从唯一 immutable native Q4 snapshot 获取所有零集 component、pure certificate 与真实 half-face 分段；XY 使用面积，RZ 使用完整 `2*pi*r`。结果保留原始 Green integral、cell capacity、source identity/time 与 owned readonly geometry。纯相初值严格为 0/capacity，cut 使用严格有界的 raw integral，无 clipping。

`initializeCandidate` 仅适用于 step-0/time-0 Initial 模板，生成保存 M 的新 primary revision，first jet 不变、几何尚未绑定。源 M 不参与积分且不被改写；Midpoint/End 的测量只能用作几何诊断，restart/advancement 初始化拒绝。报告显式区分原 conserved inventory、raw geometry volume 与两者差值。

初始模板是 I01 的模块 handle，不是已接受的数值解。测试以零 M 占位、解析 center samples 与 I06 mixed closure 构造模板，再由 I07 初始化 M。初始模板到最终 authoritative geometry 的 production 接线和数值一致性仍属 I11；没有把此处测试夹具称为完整 publisher。

## 解析与 AMR 验收

执行 **28 组、2,809 条逐单元记录**：XY/RZ 全液、全气、平面 cut、圆/球的内外相、两个分离组件、2:1 AMR、离轴 torus、平移且各向异性单元、极薄液膜。圆形逐 cell 参考值由独立解析 rectangle/disk primitive 积分，RZ 同时解析积分半径权重；不使用 production graph、根或求积。

| 指标 | 实测 |
| --- | --- |
| 最大解析总量绝对误差 | 4.4408920985e-16 |
| 最大逐 cell 初始 M 绝对误差 | 1.7347234760e-16 |
| 最大归一化求积估计差 | 9.5159990998e-14 |
| XY coarse cell 与 fine 子单元求和最大差 | 1.3183898417e-16 |
| RZ coarse cell 与 fine 子单元求和最大差 | 1.0234868508e-16 |
| 薄膜总面积 `2^-30` | 9.3132257462e-10，误差 0 |

直线/纯相满足预设 `3e-12`，圆形/多组件和逐 coarse cell 可加性满足 `1e-9`，未放宽门槛。薄膜另有相对误差 `<1e-10` 的检查，避免仅凭绝对容限放过丢失小相体积。

粗细网格验证的是初始几何体积可加性，不是动态 remap、FV 更新或 I09 reflux。

## 非多项式收敛

水平集为 `exp((x-cx)^2+(y-cy)^2-R^2)-1`，R=0.231；XY cx=0.47，RZ cx=0，cy=0.53。解析零集为同一圆/球，表示不再是可精确再现的 Q4。使用 I06 mixed closure、解析 first jets/center samples，分别在 N=8,16,32 测量。

| measure | N=8 总量误差 | N=16 总量误差 / 阶 | N=32 总量误差 / 阶 |
| --- | --- | --- | --- |
| XY | 1.2209825123e-8 | 1.1900497532e-9 / 3.35895 | 5.3854060100e-11 / 4.46582 |
| RZ | 1.9302788153e-8 | 2.7002126923e-10 / 6.15959 | 2.2375060138e-11 / 3.59311 |

误差单调下降，最后一级均满足预设 `>=3`。这些阶是本算例的体积误差实测值，不是任意几何、曲率或长时间输运的阶数证明。GL8 细分差仅为误差估计，没有声称严格上界。

## 隔离、拒绝与回归

11 个拒绝检查覆盖非法容限、变更 mesh、XY/RZ 不一致、inventory 缺 scalar 支持、Midpoint/End/restart 重置、退化零集、一段图的未解析 cap 以及结果初始化时 mesh 错配。结果/geometry 生命周期保留、初始化产生新 identity、旧 M 不变、独立库存修改均通过。

Release 主程序和全部目标构建通过：首次完整构建 99.049 s，最终局部修订后构建 14.512 s；无新增编译 warning。**CTest 11/11 PASS，总计 12.77 s，I07 contract 2.60 s**。旧 native geometry、I06a skeleton 完整输出与 I06b 相同，I03/I04 的 23 行指标相同；旧 GALS 完整输出哈希相同，parametric state control 退出 0。

1562 个既有源码/输入/工具文件中 1559 个不变，只修改共用 measure 的 h/cpp 和 tests CMake；新建 V2 h/cpp 与测试。69 个 S08 专有冻结文件、3 份方法 PDF 均不变。旧 state contract 中已有的 OpenMP 请求 12 线程、限为 2 的提示保留；它也出现在旧输出中，不是本轮新增编译 warning。S01 `hydro_moment_xy` 已知首步失败仍保留，未因这些短回归重新判为通过。

## 调试记录与已知边界

保留四次早期运行日志：

1. 试验性 `r_min>0` fixture 在既有径向 snapshot 构造处拒绝（要求 r=0）。最终验证完整 r>=0 域上的 even-Q4 torus，不声称环形截断域支持。
2. 一个仅占 corner cell 的小 cap 被 Core 完整拓扑以“fewer than three retained graphs”明确拒绝。保留该拒绝测试，另验证拥有八段图但极小体积分数的薄膜。未删除界面、未增加 fallback、未更改 Core 最小分辨率要求。
3. 缺 scalar 支持的拒绝 fixture 初次仍含已 inactive cell 的独占 edge midpoint，先被 Core 拒绝。修正测试输入的 active-side inventory 后，覆盖预期的 I07 拒绝位置。
4. 独立解析 primitive 最初使用 `asin(u/r)`，在圆端点舍入导致约 1e-10 的参考误差；改为数学等价且稳定的 `atan2(u,sqrt(r²-u²))` 后，参考与 production 一致到约 1e-16。生产积分及验收门槛未因此改变。

当前不支持不具 scalar 支持的外延 inventory 初始化，不提供环形截断域径向 snapshot，也不解决低于 Core 拓扑分辨率的 component。初始化接口明确要求完整且可解析的原生 Q4 几何；这些限制不等于最终 jet/neck/breakup 目标已完成。

## 证据和下一项

`runs/baselines/I07_20261001/` 保存 TODO 先行的 preflight/source inventory、构建与全部原始运行日志、解析/逐单元/收敛/可加性 CSV、范围/冻结/旧输出审计、源码及运行文件哈希和逐项验证的 ZIP。`delivery_manifest.json` 给出归档 SHA-256，Git 为 unavailable。

I07 的解析初始化/积分、M 独立保存与失败隔离已有代码和执行证据，故 PASS。下一项 **I08：unique phase-face flux**。I08–I12、Gate I、Hydro 低黏/无黏及最终电喷雾联合验收仍未完成。
