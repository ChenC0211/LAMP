# H08a2r4 / 共用正则 advector 生产准备报告

2026-10-05，**PASS（生产准备层）**。现有 `CchAdvectorReconstruction` 已实际持有 native P2 → 共享 Q4 → 实际子面均值校正的不可变几何映射。边界有限模式在准备时建立，应用不再求解。生产 B/J/atom、F/R/uR 和 KEP 仍为原来源；H08a 时空材料积分尚未完成，H11 未启动。

[任务书](H08a2r4_Regular_Advector_Preparation_TODO.md)、[实施前设计](H08a2r4_Regular_Advector_Preparation_Design_Change.md)、[API](../api/H08a2r4_Cch_Regular_Advector_Preparation.md)。四项材料硬约束保持：唯一共享 F、内部反号/实际 AMR 子面和、来源一致、KEP/skew。没有添加库存、界面 bubble、体积拟合、强制 commuting 或 R=0。

## 实现与来源

修改预算为 **5 production/build + 1 test**：现有 advector header/实现、Detail 几何 header/实现、根 CMake，以及原 momentum contract。没有新 App、字典或长期研究脚本。

helper 只有自有几何、实际 layout/geometry 身份、权重和 master DAG。`apply` 返回 Q4 tensor 数据及 face 修正幅度，没有第二份可发布 velocity owner。未来接受的速度必须仍由同一个 source 持有原 values、boundary 和完整 clock。

每个 datum/component 最多预制无 wall/零 normal 两种模式；只有可被实际外边界影响的数据才增加模式。RZ 轴奇偶条件一直有效。显式 law 选择模式，非零指定 wall 或轴输入拒绝；求值不拟合、不改变输入。应用调用图只包含验证、已有权重乘法、共享 Q4 装配和 face mean；没有 SVD/Eigen 调用。准备报告为不可变库存。

逐个核对了原 **60 个实际 native 布局**的 cell/face/vertex 编号、坐标、测度、相邻关系、full side/master 和修正宽度，完全相同。冻结 Core 的 60 个共享 Q4 probe 继续为零差；C++ 的实际系数另与经过该 probe 检查的原型逐项比较，不把原型当作生产依赖。

## T1：几何、边界与完整多项式

原 XY/RZ × uniform/2:1 × n4/8/16 × 4 selected/disconnected 布局，分别无 wall/全外壁/单个实际 wall，**144 组全部通过**。输入沿用冻结 native 数组，包含零场、所有相容 P2 单项式和非多项式 native 数据。

| 检查 | 最大缺陷 |
| --- | ---: |
| C++ Q4 数据对原型 | 5.07e-14 |
| face 修正幅度对原型 | 8.41e-15 |
| 独立 Bernstein/de Casteljau C1 系数检查 | 1.14e-13 |
| 独立真实面均值 | 3.61e-16 |
| Gauss cell 积分散度 | 8.95e-16 |
| 可辨识 P2 的完整场值/导数误差界，含全部校正 | 2.26e-12 |
| 全部准备权重对独立最小范数求解 | 1.76e-14 |
| 缺秩 null polynomial 的真实面矩残差 | 9.69e-16 |
| 零输入、指定 wall、RZ 轴条件 | 0 |

门槛均为原定 1e-11，未调整。C1 检查覆盖整个面上的多项式系数；P2 界包括所有 compact 修正，不能只以 datum 再现代替完整场精度。原 48 个无 wall 布局中 **38 个可辨识 P2**，其余 10 个仍如实标识缺信息。

所有准备模式独立复核 rank/权重，保存 **1,526 个独立缺秩模式证书**。这与前项按三种 law 重复列出的 3,841 次缺秩行出现次数不同；没有删除缺秩布局。采用同一最小范数定义和统一圈扩展，最多支持 24 个实际 native 面均值，最大 coefficient norm 约 26.76。

长期短契约覆盖 XY/RZ、原对象与 mesh 释放、身份复制、不同 layout/geometry/law、regrid 旧准备拒绝与 level-2 新准备、FP 模式、wall/axis 共 **15 个拒绝检查**，并检查零场、物理 P2 和重复应用完全相同。

## T2：source 精度与成本

原 12 个 n8/16/32 布局和额外 6 个 sine 布局全部运行；在原实际 FirstJet DOF 坐标逐点对冻结原型比较，值/导数最大 scaled 差 **1.69e-13**。

| 数据 / 布局 | source 值阶 | source 导数阶 |
| --- | ---: | ---: |
| XY y² / uniform、2:1 | 舍入精确 | 舍入精确 |
| RZ r⁴ / uniform | 4.03124 | 1.99849 |
| RZ r⁴ / 2:1 | 3.55389 | 2.08147 |
| sine / uniform | 3.86111 | 1.83542 |
| sine / 2:1 | 3.63919 | 2.07590 |

原 source 门槛 2.5/1.8 保持。未在本准备任务重跑原 96 个运动时间窗或改变 characteristic；有限时间与 J 的前项可行性证据仍不能代替正式来源部署后的复验。

下表为本机 Release 实测。准备时间包含同一 reconstruction 中旧来源与新准备；应用包括实际 mesh/layout 验证及系数装配，排除文件写入。每个收敛输入连续应用 20 次，结果逐项相同，XY 另有 sine 输入；应用 SVD 次数为 0。计时只用于规模观察，不作为跨机器速度保证。

| 布局 | n | 完整 prepare 秒 | 每次 apply 毫秒 | 存储权重数 | 最大支持 |
| --- | ---: | ---: | ---: | ---: | ---: |
| XY / uniform | 8 | 0.005415 | 0.036617 | 27456 | 16 |
| XY / uniform | 16 | 0.018130 | 0.137463 | 103552 | 16 |
| XY / uniform | 32 | 0.066886 | 1.428820 | 400128 | 16 |
| XY / 2:1 | 8 | 0.012202 | 0.100890 | 65104 | 20 |
| XY / 2:1 | 16 | 0.043692 | 0.395695 | 250432 | 22 |
| XY / 2:1 | 32 | 0.173979 | 3.583950 | 983552 | 22 |
| RZ / uniform | 8 | 0.005522 | 0.036345 | 25852 | 16 |
| RZ / uniform | 16 | 0.019178 | 0.132480 | 100476 | 16 |
| RZ / uniform | 32 | 0.067889 | 0.861960 | 394108 | 16 |
| RZ / 2:1 | 8 | 0.012746 | 0.090620 | 63692 | 20 |
| RZ / 2:1 | 16 | 0.045107 | 0.357870 | 247404 | 22 |
| RZ / 2:1 | 32 | 0.179342 | 3.810575 | 977532 | 22 |

最细组全部完成，权重/支持有界。当前准备仍同时保留旧来源，不能据此宣布完整动态耦合效率已验收。

## 回归、退役与后续

标准 Release 构建和 **9/9 CTest** 通过（架构/架构自检/目录，phase/material/momentum/dual/KEP/characteristic）。原五数值契约 **4,727 行**及原 characteristic **3,364 行**确定性输出与 H08a2r1 完全相同；只剔除原 `ADVECTOR_LAYOUT` 的两个计时列。momentum 末尾新增准备契约输出单列，原数值前缀逐行保持。

临时完整审计入口已移除，只保留短契约。启动 1,379 项源码/输入中 4 项按预算修改，新增 2 个实现文件，现有输入全集为 **1,381 项**。实际审计 source/exe、原生输入、每组 raw 系数、全模式权重/证书、完整日志与交付 source/runtime 分别绑定哈希，不能混用临时与最终可执行程序。

本项未重跑 H10o/H10b-v2：接受的 B/atom/F 未变且全部原输出保持。**实际来源切换时必须重跑**最终 native 算子、原内部 affine counterexample 及 H10b-v2 完整 28 窗口/84 项。

下一唯一任务 [H08a2r5](H08a2r5_Material_Reference_Polynomial_Degree_TODO.md)先解除 material/dual 参考消费者的固定三系数类型。之后按预算升级 phase/momentum 与 benchmark 消费者、完整 atom 次数和同一实际 B/J/characteristic 证书，再回到 H08a3/H08a。不得截断 quartic 项，也不得只给 Interface 换速度。

检查目录：`runs/checks/H08a2r4/20261005_134029_873`；冻结目录：`runs/baselines/H08a2r4_20261005`。`numeric_verification.json` 为完整独立数值证据，`verification.json` 为最终生产源码/回归，`data/rows_*.json` 包含全部模式的真实权重和证书，`data/raw_*.bin` 保存逐次输出，读取器为归档 `verify_numeric.py`。
