# H08a2r5 / 材料与 dual 参考积分次数迁移报告

2026-10-05，**PASS（参考消费者）**。两份现有 material/dual contract 已接收一至五系数，计算全部高次项，并处理三/四次速度在共面界面上的变号分段。生产修改 0，测试修改 2；当前 B/J/atom/F 与几何路线未变。

[任务书](H08a2r5_Material_Reference_Polynomial_Degree_TODO.md)、[实施前设计](H08a2r5_Material_Reference_Polynomial_Degree_Design.md)。本项只关闭完整多项式迁移的两个参考消费者，不表示时空材料输运或 H08a 整体通过；H11 未启动。

## 实现

`referenceAtom` 改为显式支持 constant 至 quartic normal 的模板。通用 Horner 和 work-scale 包含全部输入系数。GL3 积分核保留旧累加顺序，足以积分四次速度乘 RZ 线性半径后的五次多项式。material 中原两处三项 work-scale 也已迁移；旧 MMS 本身的二次定义保持。

共面面的相分类使用 canonical normal 速度，normalSign 只给共享积分取反。测试参考用导数根划分单调区间并括区二分求三/四次根；低次沿用原解析公式。它不消费生产 phase/F 或认证相根结果，不是新的 production root authority，也不声称能认证任意病态近重根。

显式 `(offset,width)` 参考先平移系数，再在独立 width 内使用同一 GL3 核。反向 chart 同时改变多项式坐标和径向测度方向。没有从两个 rounded 全局端点相减重建极短 width。原完整 atom 的解析相交包装仍使用浮点物理端点，未将该测试包装宣称为生产认证切分。

## 验证

每份原 contract 末尾新增以下检查；不增加 CTest、App 或长期研究入口。

| 每份契约的检查 | 数量 | 独立参考 |
| --- | ---: | --- |
| 0–4 次单项式、混合四次；完整/部分区间、双端 chart、XY/RZ 两方向与正负 normal | 216 | 闭式幂积分 |
| 四重根处的极短非零五/六阶矩，width=2^-60/2^-80 | 24 | 已知四重根的 width^5/5、width^6/6 |
| 四个简单根、重复根、端点根；真实切分/共面因果分相、两种 phase sign | 96 | 预先指定因子根定义的子区间与闭式积分 |
| 原二次 shifted 参考与当前生产 integral | 24 | 当前生产原积分及独立 GL3 |

两份合计 **720 项**。实际新增检查分别约 **1.76 ms / 1.79 ms**，没有增加完整 benchmark。整体既有 material/dual 契约分别约 56.74 / 57.89 秒，主要为其原场景运行。

实际输出的 range/tiny/phase **672 条记录**另用 **240 位 Decimal** 独立重算：相根来自已声明因子，积分直接做高精度原始幂函数端点差；不调用 C++ 的平移、Horner、GL3 或根求解。最大 scaled 缺陷 **1.70e-16**，低于固定 1e-12 门槛。两种极短 width 均满足以实际全局坐标类型 double 计算的 `offset+width==offset`，输出仍非零；四次项及 RZ 五次被积项未丢失。检查显式使用 double，避免误要求其他平台的扩展精度 long double 也必须舍入到同一点。

标准 Release 构建、架构/架构自检/目录和原 material/dual 契约 **5/5 通过**。原 material **457 行**、dual **697 行**确定性输出逐行相同，新增结果只在原输出末尾；原数值断言及容差未修改。

启动 **1,381 项**源码/输入核验仅有预算内两份测试变化，没有新增长期源码。其余四份 phase/momentum/KEP/characteristic 可执行程序哈希与 H08a2r4 相同。没有改生产重构、相材料来源、真实 R/uR、KEP/skew、FirstJet 或字典，没有库存/bubble/体积修复，也没有强制 R=0。

## 剩余事项

未运行完整 H10o/H10b-v2，因为本项没有改数值生产路径。实际来源切换后仍必须重跑原最终 native、内部 affine 和 H10b-v2 全部 28 窗口/84 项；不能用本参考测试代替 PDE 验收。

下一唯一任务 [H08a2r6](H08a2r6_Benchmark_Polynomial_Schema_TODO.md)解除三个 benchmark 导出 header 的固定三项假设。随后迁移 `CchNormalPolynomial` 的完整次数、`CchPhaseGeometry::ratePolynomial`、phase/momentum 参考，再统一部署正则 B/atom/J 与 characteristic 证书，返回 H08a3/H08a。旧三系数表示和运动骨架拒绝范围本轮仍在。

检查：`runs/checks/H08a2r5/20261005_141614_102`；冻结：`runs/baselines/H08a2r5_20261005`。`LastTest_candidate.log` 保存完整实际输出；`decimal_rows.json` 保存独立高精度逐条结果；`verification.json` 绑定源码、实际 runtime、原输出与未改 exe；临时实施脚本留在检查归档，不是项目运行依赖。
