# H10i1b / 同源闭合裁剪矩验收

2026-10-04，**PASS / engineering PASS**；本轮仅关闭 H10i1b，goal 为 progress/active。H10 最终空间与原 H10b 仍未关闭，完整流体/EHD 未验收。下一唯一任务 [H10i1c 正规 native 零边拓扑](H10i1c_Native_Zero_Edge_Topology_TODO.md)。[TODO](H10i1b_Closed_Clip_Moments_TODO.md)、[实施前设计](H10i1b_Closed_Clip_Moments_Design.md)、[API](../api/H10i1b_Native_Q4_Region_Moments.md)。

前置已核验 H10i1a 冻结的 **191 payload / 44 当前交付**。临时检查 `runs/checks/H10i1b/20261004_100116_372`，最终基线 `runs/baselines/H10i1b_20261004`；独立保存实施前设计、目录快照、实际编译源码、研究二进制/依赖、退出记录、原始 CSV 和验证脚本。最终平面 manifest 对 archived payload 及当前 canonical 交付逐项核验。

实现显式 Core geometry-only region port，并由真实 PhaseRegionSource 消费。fixed 边界取盒内侧同源相迹；原 regular graph 的 root/ordinal、owner enclosure、bounded refinement 与 controlled GL 保留，图方向矩积分后用 retained endpoints 的 shifted primitive 恢复另一方向，无修正 Jacobian。仅在 exact zero support 与原完整 zero-set normal partial 证书共同证明唯一正规 branch 后，free/fixed 重合才只计 fixed。严格 FP owner 内部 bridge 读取原 NativeQuarticTrace，不从 terminal ranges 重造 Q4/root。

不同 regionId 的内部零边两侧相迹允许不同，新账本显式为 `sum(fixed)=exterior_fixed+coincident_interface_moment`。同一 regionId 的人工 seam 排除 fixed/coincident 项，内部相界面由一条 free 曲线恢复。测试包含合并前后体积加回，防止重复计数；旧 assemblyClosure 和材料反号共享语义没有重解释。真实 canonical/transverse root refinement 与外部 typed event retry 分别计入诊断。原非负/capacity/短相及工作预算保留，无 clipping、half 边界相或质量修复。

| 固定验收 | 实际结果 |
| --- | --- |
| H10i 已出版 dual 内部位置：XY/RZ × uniform/AMR × 正负相 | **8/8** 正式 publication、H07 phase、CchMomentumMass 成功 |
| 真实 PhaseRegionSource 两侧裁剪盒 | **16/16**；液/气及两方向 Green 矩相对实际 capacity 误差 **0** |
| 真实 dual mass 独立平面参考 | **432 行**；phase 容量误差 **0**，mass/(1001·capacity) 最大误差 **9.6404e−17** |
| 新既有短契约：平面/正非仿射 Q4、XY/RZ/AMR/正负相 | **16 来源 / 208 盒 / 16 direct Core 查询**；最大容量相对误差约 **1.41358e−16**，显式 fixed/exterior/coincident 账本误差 **0** |
| direct Core 扩展 | 反向输入、共同非零 origin、same-owner 合并、缺边/重复边/坏矩形/错源/moved-from/失败后重测、临界 crossing 拒绝通过 |
| full/selected 真实质量短契约 | **16 场景 / 668 dual / 1,104 actual overlap**；独立 phase/mass 参考最大误差 **0** |
| 普通来源回归 | **64 layouts / 157,810 actual queries**；原 source CSV 保持，物理/reference 原门槛保持 |
| 最终 native 内部 affine | **3,760 行逐字节相同**；内部最大 **1.58492e−13**，原 AMR dof 911 反例 **7.41629e−14** |
| 受影响标准 Release CTest | **16/16 PASS**，含 architecture/self-test/layout、root/shifted Green/FirstJet/publication/region、phase/material/momentum/dual/KEP |

普通来源的 22 份 CSV 中 **16 份逐字节相同**，其中包括完整 spatial-source segment/patch/summary/witness 和 material source 数据；另外六份 `material_affine`、`material_compatibility`、`material_moment_bounds` 随 geometry-only 积分的浮点计算路径变化。保存逐列差异：最大绝对变化 **8.36736e−10**（高密度比 R 诊断），`max(1,abs(old),abs(new))` 归一化最大 **7.10543e−14**。这个归一化只作差异诊断，不是新增或放宽验收容限；研究程序的原独立 material expression/SVD/source/advector **1e−12** gate 全部保持并通过。原研究数值函数代码未修改。旧 KEP SUMMARY、NATIVE_MOMENTS、32 CAP、已冻结 extended NATIVE_TRACE 输出一致；不声称这六份诊断表逐字节相同，也不要求 R=0。

范围核验：**5 production / 2 既有 tests / 1 既有 research driver**，0新生产工具/target。source catalog **1,042：1,034 不变、8 变更**；inputs **346：344 不变，仅研究 driver/README 变更**，两者成员完全保持。production 为两份 MappedInterface 头、strict Core geometry TU、RegionMeasure 头、PhaseRegion cpp；Hydro F/R/uR/KEP、Interface primary/Q4/transport/clock、App/字典/PDF 与其他输入保持。最终 research driver SHA256 `B7C5058BE432DA722B7D969A20330FCCDA6D260F051413F82FCA654B21BEB509`，exe `DE093488777BC8FB029F6CB2FCE0225CA61FA6012469E010371DE7EA7F3A8A6B`，全部最终观察对应这个实际二进制。

实施过程保留原始 preflight：旧“恒零 clip 必须拒绝”测试随正式能力迁移改为真实闭式验收；新测试最初误用另一 AMR fixture，令 dual offset 落在尚未迁移的 native 零边上，随后改为原 H10i 已出版的实际 mesh fixture，原八组 native-face 拒绝仍完整保留。补充同 owner seam 不重复计数后重新构建、重跑最终 16 CTests / 四个研究模式；此前观察标为 preflight，不作为最终交付。验证脚本一次 PowerShell 系统变量名错误已修正，原错误日志保留，数值 gate 未改。

限制：本 port 支持完整 native leaf 内的轴对齐矩形分区和原正规 graph；不支持斜 polygon、一般临界界面或未认证退化。原 native-face 位置 **8/8 仍在旧 publication 入口拒绝**，这不是本任务失败或已经修复的能力；H10i1c/d/e/f 逐项迁移完整来源链。没有新 material F、空间 factory 或时空通量；原完整 H10b 未重跑，历史 14 项失败保持，H11/H08a/H10j 未启动。

材料硬约束仍只有**唯一共享材料 F、内部反号/真实 AMR 子面和、来源一致、KEP/skew**。quarter/压力到 dual 精确交换及 R=0 只作诊断；继续实际 R/uR，不把几何体积拟合或守恒库存加回 Pure FirstJet。
