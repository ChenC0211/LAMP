# H08a2r7 / 完整法向多项式表示报告

2026-10-05，**PASS（完整表示）**。生产法向多项式现在保存全部五系数，积分包括四次速度和 RZ 五次密度项；生产相率转换与两份既有参考已同步升级。实际速度重构仍为原来源，尚未部署正则 B/atom/J。

[任务书](H08a2r7_Complete_Normal_Polynomial_TODO.md)、[实施前设计](H08a2r7_Complete_Normal_Polynomial_Design.md)。预算 production 4、test 2；没有修改 CMake、benchmark、材料、characteristic 或 Interface。

## 实现

`CchNormalPolynomial` 的 begin/end 数组由三项扩为五项。`value` 用完整 Horner；`integral` 先在原 endpoint chart 的 offset 处计算移位系数，再对独立 width 积分。RZ 同时积分线性半径，反向 chart 的半径斜率取负，normalSign 独立取反。低次仍走同一公式，没有次数 fallback、微小高次截断或全局 endpoint 相减。

原私有 `ratePolynomial` 移至同模块 `Detail::normalRatePolynomial`，声明在 `Detail/CchNormalRatePolynomial.h`，实现仍在 `CchPhaseGeometry.cpp`。生产 `prepareSpatial` 与 phase 契约直接调用这一个转换；它遍历完整库存，保留既有 Core expansion 和长度/测度分离项，不创建新相根或速度 owner。

phase 的 GL3 参考、绝对系数尺度、线性检查使用全部系数；momentum 的零数组与非恒定检查使用真实库存。当前旧 B 的二次分支导数参考仍适用于旧来源，本次未将其冒充正则场的运动证书。

## 验证

新增检查放在两份原契约末尾，无新测试程序或永久研究入口。

| 检查 | 数量 | 独立参考 |
| --- | ---: | --- |
| normal 整段/子段积分与 chart 点值 | 576 行 | 物理坐标原始幂函数积分、直接幂求值 |
| normal 四重零点极短积分 | 48 行 | width^5/5 与 RZ width^6/6 闭式矩 |
| 实际 rate 转换的整段/子段积分 | 288 行 | 物理坐标 GL3 |
| 转换后 Core 正部积分 | 96 行 | 已声明简单/重复/端点根划分、独立 GL3 |
| 转换后极短普通/正部积分 | 24 行 | 四重零点闭式矩 |

覆盖 0–4 次单项式及混合四次、两种长度、XY/RZ 固定与变化半径含轴、双端 chart、正负 normal。转换的逐点密度及非零高次库存也检查；显式核验 RZ 五次项存在。宽度 **2^-60 / 2^-80** 在非零 offset 处满足 double 全局端点重合，独立积分仍非零；轴起点的 width^6 首项亦通过。

上述 **1,032 条记录中的 1,632 个积分/点值**另用 **240 位 Decimal** 复算。高精度参考直接使用物理原始多项式、给定因子根及原始幂函数端点差，不调用 C++ 的移位、Horner、GL3、转换或求根。最大尺度化误差 **3.23e-16**，低于预定 1e-12；极短矩按真实非零积分尺度核验。新增 normal / rate 检查分别约 **5.80 ms / 31.73 ms**。

标准 Release 构建、架构/架构自检/目录及 phase/material/momentum/characteristic/dual/KEP **9/9 通过**，总计约 194.62 秒。原六契约 **8,813 行确定性输出保持**；只剥离原明确计时行/列，新记录在旧输出之后，原 gate 和阈值未改。

同一既有 benchmark 源码重编后完整运行 trace-contract/advector/boundary：**34 份 CSV / 716,228 行**按列名对比，原物理/来源字段逐项相同（两计时列单列）；新增 **1,215,440 个 c3/c4/e3/e4 值严格为零**，与实际旧二次场一致。其他表逐字保持。没有伪造已发布高次速度。

启动 1,381 项源码/输入仅有预算内五个旧文件修改及一个 Detail 声明新增，当前清单 1,382 项；上一轮 141 个冻结 payload 保持。

## 剩余事项

当前完成的是表示与转换，尚非同源时空材料输运。下一唯一任务 [H08a2r8](H08a2r8_Regular_Advector_Evaluation_TODO.md)将已验证的正则系数落实为共同解析求值/导数/atom 提取内核；随后以同一内核同步切换实际 B/atom/J 与完整次数运动证书，清理旧来源专属断言并按预算迁移消费契约。不能只平滑 characteristic、保留旧材料 atom。

来源切换后必须重新通过原 native、内部 affine、H10b-v2 全部 28 窗口/84 项，再返回 H08a3 桥接及 H08a 同源时空相积分。本轮未跑完整 H10b，历史 FAIL 保持，H11 未启动。

检查：`runs/checks/H08a2r7/20261005_145304_426`；冻结：`runs/baselines/H08a2r7_20261005`。`LastTest_final.log`、`numeric_verification.json`、`mode_verification.json`、`verification.json` 分别保存完整原测试、逐值高精度复算、导出比较及源码/runtime 身份；临时实施/核验脚本仅在检查归档，不是运行依赖。
