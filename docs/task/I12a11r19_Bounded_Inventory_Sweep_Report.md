# I12a11r19 / 有界库存共享面传递报告

**状态：PASS（2026-10-02）。父任务 I12 / Gate I 仍 FAIL。**

范围：同一旧 M/geometry 的库存区域表示接入单轴恒定平移区间通量，随后消费原库存装配及 reflux。依据 [TODO](I12a11r19_Bounded_Inventory_Sweep_TODO.md)；接口见 [API](../api/I12a11r19_Bounded_Inventory_Sweep_API.md)。完整时间驱动、I10/I11 与原 I12 gate 未切换。

## 实现

2 个 production 文件：`HermiteGalsTemporalPhaseFlux.h/.cpp`。增加显式库存来源 overload，与原几何来源共用扫掠/边界内核。旧来源身份及同 snapshot、时间区间、owned constant velocity、实际 map、覆盖、RZ 轴向限制、显式入流历史均保留。原总面量不变，库存相量严格位于 gross 总量界限内。

新增诊断从同一次区域查询获得原几何液量、库存修正和 L1；未建立第二份生产通量。每个实际 native face 的唯一缓存结果进入原库存装配和 reflux。报告仅输出；原代数 Transfer 仍不认证外部积分算法。原 M、原 64 epsilon 修复预算、Core 求积和 I10/I11 均未修改。

2 个 test-related 文件：新增 `interface_inventory_sweep_contract.cpp` 并注册 CTest。

## 算子验证

84 场景、4536 个实际面，包括 XY 两方向与 RZ 轴向、uniform/AMR、正反相、正反速度及静止。48 组平面排空/充满使用两种切分时刻；24 组移动曲线和 12 组静止曲线含几何纯相中的微量库存差。制造的相对容量差为 4e-11，使用默认 1e-10 几何门槛。

48 组平面中，旧几何通量均在原库存装配处明确拒绝，新库存通量全部通过同一装配。测试保留旧库存、几何 snapshot 和时间身份，未在推进后解析重置 M。移动 CFL 为 uniform .375、AMR fine .75。

| 测量 | 最大值 |
| --- | --- |
| 面更新与 departure-region 库存差 / 容量 | 1.11022302462516e-15 |
| 逐 cell 共享面余额缺陷 / 容量 | 1.02805807885689e-16 |
| 同源半区间拼接差 / 总传递量 | 1.25825276124184e-15 |
| 单面库存液量相对旧几何参考的绝对修正 | 1.37444708392387e-11 |
| 全局边界 ledger 绝对缺陷 | 2.49800180540660e-16 |
| coarse=sum(fine) 液量/总量残差 | 0 |

库存严格有界；前两项均低于原 `64*epsilon`，区间切分低于原 1e-10。原共享面舍入修复仍可能工作，单场景最多修复 3 面，没有增大预算。静止 M 保持逐值相同。

departure-region 参考使用同一库存区域算子，并以独立物理容量公式处理域外全相入流；它检验扫掠分割和面更新的一致性，不构成独立几何求积精度证明。独立常量积分及区域权重验证见 r18。7 类拒绝覆盖 moved-from、错误 mesh、缺失入流历史、不同旧 state、对角速度、不透明速度回调和非法容限。

## 初始失败与证据

初次编译因测试中重复的 `Interval` 别名失败，移除未使用别名后通过。首次执行在第一组 RZ 平面失败：测试将 r=0 轴端点误声明为 non-axis 边界接触点。修正为仅声明外径端点，轴端点由原轴对称几何契约识别；未改变生产边界规则。原失败日志保留。修复后的直接运行 17.7143676 s，84 组通过。

33 个目标完整构建通过，48.3853944 s，无 warning。17 项相邻 CTest 全通过，339.5393196 s；15 项旧数值输出与 r18 归档逐字相同。架构保持 12 targets / 382 translation units / 905 files，仅新增所需头文件依赖后为 3110 includes。旧 `--crossing-probe` 仍以相同输出 exit 1，未掩盖原失败。

81 个冻结文件、3 份方法 PDF 及 r14/r18 archive hash 未改变。证据目录：[I12a11r19_20261002](../../runs/baselines/I12a11r19_20261002/)。其中 `operator_first.log` 与 `build_first.log` 保留初始错误，`operator_verified.log` 为最终矩阵；`verification.json`、`source_audit.json`、`regression_comparison.json` 和 `sweep_comparison.json` 给出验收。可复核 `evidence.zip` 的哈希和逐文件核验见 `archive_verification.json` / `delivery_manifest.json`。

## 未解决事项

纯相闭合目前仍按旧库存 roundoff 接受 M−V；r18/r19 可运输更大的预算内差额，但不能据此绕过 I10/I11 发布。下一项应明确几何纯相与微量库存差的语义并验证发布不制造成分、不重置 M，然后接入半步/全步和原 h/dt/曲率矩阵。一般 rotation/vortex 区间积分与原 I12 门槛仍未通过；**父任务 I12 / Gate I 仍 FAIL**。
