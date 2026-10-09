# I12r1 / gross phase-flux bounds

## 任务 ID / 标题

I12r1：定位并修复 I12 vortex 中的液相 gross flux 超过同向总 flux 问题。此项属于 I12 的独立算子修复，不实施尚未获准的 I12a 时间积分设计变更。前置 I12 FAIL 归档已核验；上一轮为 progress。

**结束状态：PASS（仅本项）；I12 仍 FAIL。** 2 production + 2 tests；原 vortex 首步、48 组算子配置和前序 15 项回归通过。完整 benchmark 的均匀 vortex 推进至 t=0.031 后仍因 bubble bracket 拒绝。详细交付及证据见 [报告](I12r1_Gross_Phase_Flux_Bounds_Report.md)。

## 目标

使同一 midpoint scalar/速度的正反向液相空间积分，在既有积分误差与舍入预算下满足 `0<=Q_l^±<=Q^±`，保留总通量和共享库存更新的权威定义。先量化现有失败，再确定最小修复；不能以放宽 gate 代替修复。

## 允许修改

- `src/Modules/Interface/HermiteGalsPhaseFlux.cpp`，必要时对应 header：同一空间积分结果的一致性及诊断。
- 若明确定位到共用 positive polynomial primitive，允许 `src/Core/Geometry/LevelSetHermiteGeometry.cpp`、`include/Core/Geometry/UnitPolynomialRoots.h` 内该 primitive 的必要修改，先验证不改变其余几何路径；预算最多 4 production。
- `tests/interface_hermite_gals_benchmark.cpp`、`tests/interface_hermite_gals_flux_contract.cpp`：原始复现与算子回归，最多 2 tests，不增加第三个测试文件。
- 本 TODO/报告、相关 API/README，`runs/baselines/I12r1_20261001/`。

## 禁止修改

phase-flux 时间离散、I02/I03/I04、I10 guard/容限、三份 PDF、S08 frozen sources、Hydro/Electricity/App/reader。禁止 M clipping、增大 roundoff budget、放宽 gross bound、忽略 finite 校验、隐式 fallback 或将 I12 其它失败重分类为成功。

## 数学约束

1. 同一 scalar 的 retained half-face ranges 与同一法向速度 polynomial；RZ 使用相同真实 `2*pi*r` 度量。
2. 总通量由给定速度决定，不因 phase 几何改变；正反向分解满足已有正部积分语义。液相是同一正部测度的子集。
3. 同一个实际内部 face 的最终液量只计算/写入一次，反号写入相邻 cell；coarse 汇总实际 fine flux。
4. 仅允许既有精度预算内、具备误差依据的算术一致性处理；原始结果/改变量不得隐藏，不能把实际大误差当舍入修复。
5. M 不裁剪；失败不污染旧 state。原始 midpoint 跨面反例仍须明确失败，除非后续已授权独立时间离散变更。

## 实现要求

先保存原 vortex uniform 第一个 half-step 的 face/half/range/velocity coefficients、正反向 full/liquid/gas 积分、误差估计、缺陷量与尺度。以测得根因确定修复，保留最小可复现算子测试和真实 vortex 首步测试；无 production debug 开关或文件写入。

## 验收

- **T0：** Release 主程序及相关 contracts，无新增 warning，architecture、文件预算/冻结 hash 通过。
- **T1：** 原失败的 gross bound 通过；pure-liquid/gas、counterflow、速度根与相根相对位置、XY/RZ、2:1 AMR、小 phase interval、符号翻转均验证；总通量的 phase 独立性和共享 ledger 保持。
- **T2：** 原始 vortex 首步不再在 gross subset bound 提前失败；若在 bubble/其它阶段失败，保留原因和 accepted time，不宣称 vortex/I12 通过。
- **回归：** I08/I09/I11 及前序数值测试、legacy GALS/parametric controls；输出变化逐项说明，不要求有意修正的最后舍入位与旧结果逐字相同。
- **不变控制：** I12 平面 q=1/4 与 3/4 的时间取样反例维持原始失败；其修复仍等待设计变更决策。

## 失败规则

若缺陷超过既有数值预算，须定位真实公式/区间错误；禁止直接截断。若修复需改变冻结时间离散或超出预算，记录具体依赖并拆分，不跨越未批准变更。本项 PASS 不等于父任务 I12/Gate I PASS。

## 交付

代码/算子与原始 vortex 复现、前后误差明细、测试/失败日志、源码/runtime/归档哈希与未解决事项。总体目标保持 active。
