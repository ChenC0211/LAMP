# H08a2r11 第二次恢复 / 同一正则来源同步部署报告

2026-10-05，**PASS**。本次接受共同 regular source/characteristic；原 r11 两次 FAIL、r16 原 Jacobian FAIL 和四项旧 donor gate FAIL 均保留历史裁决。实现前范围见[部署设计](H08a2r11_Common_Regular_Source_Deployment_Resume2_Design.md)。

## 改动及来源

四个生产文件为 `CchSpatialAdvector.h/.cpp`、`CchCharacteristic.h/.cpp`，两个测试为原 momentum/characteristic 契约。四个生产文件采用 r17 冻结候选，仅 characteristic header 同步公式注释；未再精修算法。momentum 只迁移 prepared owner 身份检查；characteristic 采用 r17 已验证九组、同一可信参考及原门槛，清理未调用的旧 RK4 参考。1383 份源码/输入中仅上述 6 份变化，其余 1377 份（含 Interface、材料/dual/KEP、研究 driver/输入、工具及 CMake）保持。

`prepare` 持有唯一 regular 权重、真实 support 分片与 primal/dual overlay；每次 `evaluate` 只 apply 一次，不重复 SVD。实际值/J、native 物理面中心值和完整五系数 normal polynomial 全部来自同一私有不可变系数及原 velocity/geometry/boundary/full-clock 身份。native 面中心值与 native 面均值语义明确区分；cellDivergence 保留 D_I/V 的积分平均。

旧 `CchAdvectorPatch`、patches/faceB/cellLine 及旧 trace/lift 实际来源退出生产，无 fallback。运动采用 r17 同源三阶分片 Picard 内部积分及三次 dense 事件；积分节点、宽度、连接时刻同步求导，完整 Q4 值/J 迹证书替代旧低次数假设。没有修改 Pure FirstJet、材料 F/R/uR 或最终 native 对流算子。

## 原验收结果

| 核验 | 结果 |
| --- | --- |
| 标准 Release、架构/目录、受影响契约 | 9/9 PASS，454.74 s |
| 实际部署 characteristic 原矩阵 | 24 时间族；最小 point/J 阶 2.998745 / 2.998115；288 FD、2304 分区比较、54 完整次数负控通过 |
| 原两条非线性跨片 2560 步可信参考 | point 和 J 原 1e−11 门槛均通过；J 误差加参考预算（含 double 转换余量）7.368e−13 / 2.595e−13 |
| 实际代码 point-only 变异 | 原 affine J、uniform J、internal AMR J 三组明确 FAIL，负控有效 |
| 材料动能与 R/uR | 原 124 场景/6144 DOF，最大身份缺陷 1.862e−16，原四 closed 非零材料功保留 |
| 同源物理观察 | 原 64 布局/128 revision/56 适用解析布局；77540 faces、297800 atoms、36040 cells；208 basis、64 边界 revision 全部通过 |
| 内部 affine 反例 | 12 组/15344 行，最大误差 8.371082e-13 ≤ 1e−10 |
| 最终 native 独立矩阵 | 12 组/15344 行，最大作用误差 1.204426e-15 ≤ 1e−12 |
| H10b-dissipation-v2 原完整研究 | 28 窗口、15116 RK4 步、928 时序行、84/84 gate PASS；四旧 donor gate 仍 FAIL |

独立物理参考复算原完整双端多项式及实际采样：native 物理均值最大相对误差 3.423e−16，原平均散度误差 1.188e−14，原 affine 点/atom 误差 2.332e−15。原 schema/覆盖保持；r14 修复过低阶权重，故不要求当前实际场与旧冻结场逐值相同。来源/完整 clock/regrid、axis/wall、tiny/反号及负向次数证书保留。外部 affine 最大误差 4.661670 继续单列诊断，未将其冒充内部二阶或添加新门槛。

原 H10b-v2 n16→32 空间结果（原门槛：两阶数 ≥1.7，最细 L2≤0.05）：

| 组合 | 截断阶 | 解误差阶 | n32 L2 |
| --- | --- | --- | --- |
| XY / uniform | 1.991655 | 1.990209 | 0.04029228 |
| XY / 2:1 AMR | 1.825562 | 2.003767 | 0.03103141 |
| RZ / uniform | 1.991655 | 1.990270 | 0.03657952 |
| RZ / 2:1 AMR | 1.769790 | 1.958653 | 0.03154876 |

24 份 donor 功快照由原 native/segment 表独立重组：同一 F、真实 q/R、uR、耗散、材料功和边界功满足原 1e−12 门槛。该固定材料光滑对流研究仍不是完整 Euler/NS 或一般 sharp-cut 两相验收。

## 成本、复现及交付

advector 原 64 布局累计 prepare 2.084129 s、128 次更新累计 evaluate 0.203375 s；n32 affine/native 审计分别 759.16/595.40 s（两模式并行，非独占计时），原 full 916.53 s。9 项短契约从 r17 旧来源的 196.88 s 增至本次 454.74 s；场复杂度增加，动态准备效率未验收，不能宣称性能通过，也不据此另开 source 精修任务。

完整研究走原通用入口：

```powershell
.\tools\run_case.cmd benchmarks/cch_convection -Threads 1 -Label H08a2r11_resume2
```

本次[标准运行](../../runs/benchmarks__cch_convection/20261005_230436_872_H08a2r11_resume2_d4e57c2f/run.json)保存原输入/程序/完整输出；[检查记录](../../runs/checks/H08a2r11_resume2/20261005_223648_021/actual_verification.json)与[冻结基线](../../runs/baselines/H08a2r11_resume2_20261005/README.md)绑定本次实际源码、测试、运行库、负控和独立复算。基线保留上一 r17 的 102 份 payload 身份；临时脚本只归检查/冻结目录。后续文档/目录检查与最终哈希复核记录在 verification.json。

## 基础设施收口与下一项

按用户 2026-10-05 指令，**r11 为 common regular source/characteristic 最后一道关口**。除可复现且直接阻断 H08a 的数学反例外，不再拆 source/Jacobian/低秩精修，不新增硬约束、特殊 case 或放宽门槛。若需要例外，必须先保存原输入、实际来源、独立参考和直接阻断积分的证据，普通误差/成本诊断不能自动产生先决任务。

下一唯一任务立即回到 [H08a 真正同源时空相通量积分](H08a_SpaceTime_Material_Flux_Upgrade_TODO.md)，必要 FirstJet 运动消费接入并入该目标。现有材料时间仍为中点近似，步内几何消费、真正积分及 R 收敛尚未完成。只使 R 成为可收敛的小余项；Pure FirstJet 独立，无 M_K/bubble、无 R=0 或 quarter commuting 强制；继续同一 F、内部反号/真实 AMR 子面和、来源一致、R/uR 与 KEP/skew。H11 未启动。
