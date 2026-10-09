# I12a12 / 同源多项式速度区间源报告

**PASS（本项）；I12 / Gate I 仍 FAIL。** 按 [TODO](I12a12_Polynomial_Velocity_Interval_Source_TODO.md) 与先写定的[设计](I12_General_Velocity_Interval_Source_Design.md)完成 owned 速度数据、连续区间包围及原 RK2 view 接线。所有计划测试和来源审计于 2026-10-02 完成；一般相通量和完整旋转/涡旋推进尚未实现。

## 交付与保持

新增中立 Core `PolynomialVelocitySlab`，使用同一组 tensor cubic 空间系数生成速度值和一阶导数，时间上采用已有 endpoint 线性插值。其区间求值复用既有 outward arithmetic，加入原 strict/no-LTO geometry 对象；没有复制到 Interface/App 的舍入内核，没有增加全局速度未知量。

`HermiteGalsTemporalVelocitySource::polynomial` 从同一个不可变 slab 生成实际 old/end callbacks；`polynomialArrivalView` 为真实 candidate/time 生成同源 End/Midpoint 视图，保留 owner 生命周期。opaque 和 constant source 不自动获得新 capability，原平移 producer 仍要求 owned constant；未设置失败后切换。

原 RK2、B/J、边界事件、CFL/轴条件、M、bubble、几何发布与 reflux 算式均未改。范围是 5 production/configuration + 2 test-related 文件，详见 [API](../api/I12a12_Polynomial_Velocity_Interval_Source_API.md)。

## T1 结果

8 类来源覆盖：定常/时间变号线性旋转、定常/时间变号 cubic differential rotation、非单位坐标基底的一般 cubic 场、RZ 线性应变、RZ 轴向变号平移和零速度。共 5832 个独立解析值/导数样本；每个样本同时核验整个 box 范围、退化点范围和实际浮点求值。解析内部极值、时间反向、缩小 box 和零多项式另作直接检查。

| 检查 | 最大误差/宽度 | 固定验收上限 |
| --- | --- | --- |
| 多项式求值/导数 vs 独立解析表达 | 5.5511151231e-16 | 2e-14 |
| 实际 RK2 B/J vs 原 callback A/B | 2.3852447795e-18 | 5e-14 |
| 实际离散 B 的中心差分 vs 返回 J | 1.1810067247e-11 | 2e-9 |
| 点 box 的最大包围宽度 | 1.7763568394e-14 | 2e-12 |

32 次实际 B/J A/B 包含三个 End 时间和真正的 Midpoint 时间层。45 项拒绝/契约检查覆盖非法 scale/coefficient/box/time、overflow、舍入模式和 FTZ、opaque/constant 能力边界、moved-from、step/extent/candidate identity、原 CFL 与 RZ 轴规则。输入修改不影响已构造 source，拷贝与 interval 保持同一 owner，retained views 在局部 source 销毁后仍可使用；accepted geometry identity、snapshot 和 M 未改变。

包围结论来自既有区间运算的逐步包含性及严格浮点前提，采样测试只核验实现。范围可能保守；未把密集采样解释为一般速度/相事件证明。

## 原 RK2 的特征时间收敛

定常旋转速度 omega=.8，固定终点 .2，实际原 RK2 重复推进特征点：

| dt | 步数 | 位置误差 | 阶数 |
| --- | --- | --- | --- |
| .05 | 4 | 8.7956990147e-6 | — |
| .025 | 8 | 2.1989733286e-6 | 1.99996813 |
| .0125 | 16 | 5.4974639114e-7 | 1.99999197 |

此处只验证特征算子。fixture 为速度/时间身份提供结构性 state，未推进相界面，不以该测试代替 I12 的完整旋转、曲率或相体积时间阶。一般时空相通量尚未消费这些数据。

## T0 与归档

首次构建通过，8.8127124 s；首次直接测试通过，0.1755039 s。之后补齐真实 Midpoint、CFL 和轴拒绝覆盖，重建全部 34 个运行目标通过，57.0104332 s。最终 20 项 CTest 全通过，164.2322786 s；18 项旧数值输出与 r14/r19/r20/r21/r23 的最近归档逐字相同。旧 crossing probe 保持相同输出与 exit 1。

无 compiler warning。architecture 仍为 12 个主构建目标，源码统计由 382→383 个 TU、905→907 个文件、3110→3114 个 include，变化对应新 Core header/TU；测试源码另由改动清单审计。所有旧依赖方向检查通过。81 个冻结文件、3 份 PDF 和 r23 archive hash 不变。

证据目录：[I12a12_20261002](../../runs/baselines/I12a12_20261002/)。`operator_verification.json` 保存数值及特征阶数；`regression_comparison.json` 保存 18 项旧输出对照；`source_audit.json` 和 `strict_compile_evidence.txt` 记录改动范围及实际生成的 Strict/no-LTO 配置。`delivery_manifest.json` / `archive_verification.json` 核验 `evidence.zip` 的交付源码、34 个程序、运行库与证据。所有构建完成后才执行最终回归，没有重建运行中的程序。

## 下一项与未完成范围

下一项需要将这些速度范围沿**同一个实际 RK2 运算**传播，给出空间—时间回溯 box、实际阶段与物理域覆盖的可验证条件。还要处理跨叶相区间、速度变号、边界历史和事件完整性，才可积分共享相通量。

Core 的 endpoint 混合和直接 slab 求值在实数层面表示同一线性重构，浮点表达次序仍可能不同；下一项不能跳过实际混合/分支。cosine vortex 的时间重构误差必须相对原连续真解验证。相通量有界性、完整旋转/涡旋时间阶、动态 AMR 和 Gate I 都未由本项通过，不进入 Hydro。
