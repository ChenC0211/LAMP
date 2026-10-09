# I07 / 守恒单元相体积状态与初始化

## 任务 ID / 标题

I07：conserved cell phase volume state；唯一验收点为 analytic geometry initialization / volume integration。

## 目标

从同一 native quartic Hermite 零集积分得到初始液相体积 `M_K=∫_K H(-phi_h) dV`，形成可独立保存的 primary 库存；后续几何测量不覆盖库存。

## 允许修改

- 新 `include/Modules/Interface/HermiteGalsPhaseVolume.h`、`src/Modules/Interface/HermiteGalsPhaseVolume.cpp`。
- `InterfaceCellMeasure.h/.cpp`：为 V2 开放同一完整零集发现与既有 moments-only 积分，不依赖旧 InterfaceState，不复制积分公式。
- 新 `tests/interface_hermite_gals_volume_contract.cpp`、`tests/CMakeLists.txt`。
- 本 TODO、API、报告、docs README、`runs/baselines/I07_20261001/`。

预算：4 production + 2 test；不扩张后续任务。

## 禁止修改

I01–I06 离散公式、mesh、Hydro/Electricity、App/reader/case、冻结 PDF、S08 专有实现。不得进行通量更新、reflux、bubble 求解、I11 几何发布或 coupled acceptance；不增加生产方法入口。

## 数学约束

1. XY 为面积，RZ 为含完整 `2*pi*r` 的体积；使用真实 native cell/half-face，AMR coarse/fine 不重叠计数。
2. 复用同一 snapshot 的完整曲线 inventory、phase certificate、实际边界负相分段和 free-curve Green moments。禁止 PL/中心判相、第二份零集标量和体积 clipping。
3. 保留 raw geometric integral。经证明的纯相初始化为严格 0 / 实际 cell capacity；cut cell 使用严格有界的 raw integral，越界明确拒绝。
4. 测量结果保留 source identity/time 和只读源。初始化只允许 step=0、Initial、time=0 的未接受初始模板；返回新的未发布 candidate，不能原地覆盖已输运库存。初始化后不自动复用旧 revision 的 geometry/closure。
5. I01 的 completed State 是模块 handle，不等于 coupled accepted。初始模板只供 I07 测量，模板内 M 不作为测量输入；I11 仍须验证并发布最终 primary 对应的唯一几何。源 handle 和已有 M 全程不变。
6. 本任务要求每个库存 cell 都有本次 scalar 支持；不从无 scalar 的 cell 猜测纯相。AMR/restart 库存迁移不属于初始化。

## 实现要求

只读 Result 提供按 native index 排序的 M、capacity、raw integral、phase certificate 及具名诊断。积分误差估计与边界 ledger 闭合容限分开，前者不是严格误差上界。对 mesh、measure、activity、时间层、source 一致性、非有限值、缺 inventory、越界和未解析零集显式拒绝；异常带模块和 source identity。失败不污染源或 partial candidate。纯几何积分 API 不接 velocity/物性/文件路径。

## 验收

- **T0：** Release 主程序、新 contract、I01–I06、architecture；头文件自包含，无新增 warning，只读结果。
- **T1 解析：** XY/RZ 全液、全气、直线 cut、圆/球、多个分离界面；uniform 与 2:1 AMR。总量及逐 cell 与独立解析 primitive 比较，归一化绝对误差 `<=1e-9`；直线和纯相 `<=3e-12`。记录每个 cut cell、capacity、raw、M、证书及误差。
- **T1 守恒库存：** 初始化新 candidate 的 M 与积分一致，first jets 不变；原模板 M 不变、identity 不变；新 identity 且无 geometry。测量带相同 geometry identity/time；非初始/restart 初始化拒绝。
- **T1 隔离：** 错 mesh、XY/RZ 不一致、缺 scalar 支持、非法容限、未解析/退化几何等失败；无 clipping/fallback。同源实例生命周期检查。
- **T2：** 同一解析曲线 coarse/refined 测量逐 coarse cell 与 fine 子单元之和相符 `<=1e-9`，多组件和小 cut 不丢失；非多项式平滑圆形水平集 N=8,16,32 的体积误差下降并记录实测阶，最后一阶 `>=3`。其误差来自 Hermite 表示与求积，不能称严格解析几何。
- **回归：** I01–I06、旧 native geometry、旧 GALS 和 parametric control、architecture；69 个 S08 frozen files 和 3 PDF 哈希不变。共用旧积分数值回归，无任务外源码变化。

## 失败规则

数学验收失败记 FAIL 并在范围内修正；需要改冻结方法、超出预算或缺少必要接口时记 BLOCKED 并说明证据。不得放宽既定容限来配合结果，不把局部积分 PASS 宣称为 I08–I12 或完整目标完成。

## 交付

修改清单、API、执行报告、原始构建/测试/失败与数值 CSV、源码和运行文件哈希及证据归档。所有验收通过后才记 I07 PASS，下一任务 I08；总体 goal 保持 active。

## 执行结论

**PASS，2026-10-01。** 见 [报告](I07_Conserved_Phase_Volume_Report.md) 与 [API](../api/I07_Hermite_Gals_Phase_Volume_API.md)。28 组解析几何、2809 条逐 cell 数据、XY/RZ 空间收敛和 AMR 可加性、11 个拒绝、11/11 CTest、旧数值回归及冻结哈希通过。未解析的一段图小 cap 保留显式拒绝，环形截断径向域未实现；未降低任何验收门槛。下一项 I08，Gate I 未完成。
