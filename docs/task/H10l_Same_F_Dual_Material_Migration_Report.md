# H10l / same-F dual 迁移与失败反例

2026-10-04，**FAIL；恢复后的 engineering PASS**。仅完成本项候选、必要回归及失败收口，goal 保持 active，本轮为 progress。[任务书](H10l_Same_F_Dual_Material_Migration_TODO.md)、[实施前设计](H10l_Same_F_Dual_Material_Migration_Design.md)未放宽数值门槛。候选已归档，八份源码恢复到 H10k；当前生产仍含旧 H09 quarter，不能声称 same-F dual 迁移完成。下一唯一任务为 [H10l1](H10l1_Spatial_Advector_Closure_And_Boundary_TODO.md)，本轮没有实施它。

## 候选范围与实际身份

候选仅修改五份 production：`CchPhaseGeometry.h/.cpp`、`CchMaterialTransport.cpp`、`CchDualMaterial.h/.cpp`；两份既有测试：dual/KEP；原研究 driver 增加两段只读 observer。无新源码、target、tools、输入配置或 App 路径。Core/FirstJet、actual advector、phase spatial 积分、唯一 Flux 数值、D_I/C、真实 dual 几何/质量及 KEP 生产算法保持；原 uR/材料动能数值实现逐字不变。

候选 `faceRates()` 返回原 `material.dualFaceRates()` 的同一引用，只有这份 F；删除 quarter/SPD/局部材料延拓及最后 legacy phase/material 表达。保持真实三层质量和完整来源，实际 incidence 计算

\[
q_d=\sum sF,\quad r_d=\Delta m_d+\Delta t q_d,\quad
E_{sp}=q_d-Tq_p,\quad G=\Delta m_d-T\Delta m_p,
\qquad r_d=Tr_p+G+\Delta t E_{sp}.
\]

T 只使用原真实 overlapFraction。E_sp 是率，G/r 是质量；三者可以非零，报告没有填零、旧交换硬约束或界面质量修复。候选字段和新行为仅存在于归档中，不是当前生产 API。

## 必要回归为何拒绝

首轮完整必要 CTest **6/8通过、2失败**。精确来源拒绝信息的先后次序已修正后，只重跑受影响 dual/KEP 两项，仍全部失败，最终真实失败为三组 basis 和四组 closed。测试打印具体见证、继续覆盖其余场景，最后非零退出，没有忽略失败。

| 原单面 basis 场景（n4） | 首个输入 | 原拒绝 |
| --- | --- | --- |
| XY 2:1 | u3，area=0.25，value=4 | compatible field differs from original D_I |
| RZ uniform | u2，area=π/4，value=1.2732395447351628 | 同上 |
| RZ 2:1 | u2，同上 | 同上 |

拒绝发生于原 `CchSpatialAdvector`，先于 phase/material/H09。原检查按 native `A*u` gross 核验实际空间 divergence；此处尚未独立定位为几何/系数问题还是零净流量的浮点 work 标度问题，不能先放宽容差或把实际 divergence 贴成 D_I 期望值。

| 原封闭流函数场景 | segment / native u | 实际 boundary F | 原 native mean / primal F |
| --- | --- | ---: | ---: |
| XY uniform | 32 / 0 | 0.0013671875 | 0 / 0 |
| XY 2:1 | 72 / 0 | 0.0016961348684210518 | 0 / 0 |
| RZ uniform | 41 / 4 | 0.0072304960636530684 | 0 / 0 |
| RZ 2:1 | 81 / 4 | 0.0093881634531974498 | 0 / 0 |

零 native 面均值没有给出 pointwise 不可穿透条件；actual B 重构保留了边界变号 normal trace，使子段 F 非零。不能在 H09 压掉这些 F，也不能按 mean==0 自动清除斜率，否则会破坏必须保留的内部零净 counterflow。需要下一任务提供明确、通用的 typed 不可穿透来源，并独立核验同一 B、D_I 和所有实际边界 atoms。

## 通过的有限范围

候选 dual 短契约新增 **264 次准备/42,240 atoms**，含64个 counterflow、136个自身尺度小相见证；独立最大误差 **1.928e−16**。真实 R、含 dt·E_sp 的分解、共享引用、uR/材料动能及来源拒绝通过。其他未异常终止的 dual/KEP 身份、能量、constant 与 spatial moment 场景完成；失败 basis/closed 函数中拒绝点之后的检查没有执行，包括四组 closed 功率检查，不能计作通过。

研究使用真实 Initial/Midpoint/End GALS 和同一 actual B：平面静止+MMS、移动圆球+MMS两族属于强制几何；移动圆球+同向0.4平移为 coupled manufactured 族。每族三密度1/1e3/1e6，保留原 full/selected、XY/RZ、uniform/2:1、n4/8/16/32与56个预列 coarse selected tangent 缺口，没有新增豁免。

| 独立实际 CSV 核验 | 实测 |
| --- | ---: |
| 布局 / phase factories / material preparations | 64 / 192 / 576 |
| dual cells / shared segments | 348,930 / 1,071,360 |
| pressure cells / primal faces / overlap 行 | 162,180 / 348,930 / 655,200 |
| 独立 phase / F 最大尺度误差 | 3.775e−15 / 3.766e−15 |
| 独立三层 dual volume / mass | 4.613e−13 / 4.613e−13 |
| R/E_sp/G/incidence/global balance 身份 | 5.457e−15 |
| uR/材料动能身份 | 0 |
| 三真实层×三密度 gas 见证 | 63行，误差0 |
| coupled pure-gas u103 的三份 interval q/R | 均为0 |

独立脚本按声明物理多项式、解析相交、circle-rectangle primitive、真实 native face/pressure 坐标重建半 cell overlap、三层质量和真实 incidence，不用 production F/R 当参考。共享内部面按反号求和、真实 AMR 子段和与 component boundary/balance同时核验。最大实际 |E_sp|=298,456.53、|R|=519,600.34，只是含强制运动族的有量纲诊断；未声称一般 R 小余项、coupled 收敛或完整 PDE 通过。Pure-gas 的零 q/R 也不代表整场 R 为零。

原内部 n8/16 affine 四组合重跑通过，最大误差 **1.075e−13**。原 driver 删除两段 observer 后逐字等于 H10k，full/RK4/donor 和 gate 均不变；没有运行 n32 affine 或原完整 H10b。原28窗口/84 gate及14失败历史保留，H10/H11/H08a尚未关闭。

## 归档、恢复与交付

启动前独立复核 H10k：229 payload/109 current canonical，及 j/j0/i2/f/e/d/c 历史归档。候选 scope 验证5 production/2 tests、1,044成员/1,036不变；全部八份候选 source 及 runtime 已保留，随后八份 canonical source 恢复，**1,044份源码与启动前哈希全部相同**。案例配置、CMake/App/PDF保持。

候选初版 compiled source/runtime 对应完整 CTest、`affine`、`dual` 研究。最后候选仅修正来源检查顺序及细化测试见证，candidate2 source/runtime 对应最后两项失败 CTest；没有拿初版研究 exe 冒称后版编译产物。恢复版另存 compiled source/runtime，标准 Release 重建六 targets，必要 CTest **8/8**、最终架构/目录检查 **3/3**。

一次 C++ auto 多类型声明编译失败、PowerShell target-array/脚本解析和文档脚本 alias 失败均保留，修正不改变数值。首次恢复复制保留了旧源码时间戳，增量 MSBuild 未重编译，导致 CTest 实际仍执行失败候选；保留该次日志后，更新八份恢复文件时间戳并重建，再以恢复版二进制复测，源码哈希保持。独立后处理初次把近零最大 E 差归一化为 E 自身，错误放大 signed cancellation；最终使用实际 `gross_dual+mapped_gross` work 核验 E 诊断，与逐 cell 身份相同，1e−12门槛保持。初稿脚本保留，没有修改 CSV。圆球独立质量参考在最高细度仍满足原门槛。

全部实际日志/CSV、候选与恢复版 source/binaries、实施前设计、scope/独立核验和 flat SHA256 manifest 见 [H10l冻结目录](../../runs/baselines/H10l_20261004/README.md)；检查目录 `runs/checks/H10l/20261004_144355_065`。失败候选未留在生产，也不作为 fallback。

下一唯一 [H10l1](H10l1_Spatial_Advector_Closure_And_Boundary_TODO.md)：2 production（actual advector h/cpp）/2既有 phase/material tests，先定位 basis closure，再明确 pointwise 不可穿透来源。它通过后另起任务按5 production/2 tests重试 H10l；再联合 final native constant/affine 与 material-KE/skew，内部 affine 先通过，然后原完整 H10b，H11/H08a等待。材料硬约束仍仅唯一共享 F、内部反号/实际 AMR 子面和、来源一致、KEP/skew；quarter/严格 commuting/R=0仅诊断，Pure FirstJet独立，无库存/bubble。
