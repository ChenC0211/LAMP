# H10j / 同源空间相通量 owner

2026-10-04，**PASS / engineering PASS**。只关闭 H10j；长期 goal 保持 active，本轮为 progress。新空间 phase factory 已实现，H08/H09 的材料 F 与最终 native convection 仍待迁移；没有重跑完整 H10b，不进入 H11/H08a。[任务书](H10j_Spatial_Phase_Rate_Owner_TODO.md)、[实施前设计](H10j_Spatial_Phase_Rate_Owner_Design.md)、[实际 API](../api/H10j_Cch_Spatial_Phase_Rate.md)。

## 实际改动

固定 2 production：现 `CchPhaseGeometry.h/.cpp`；只扩展现 `cch_phase_geometry_contract.cpp`。原 CMake owner/target 保持。既有 research driver 增加四段只读 observer：`spatial` 模式观察新工厂；`trace-contract` 原 16 来源链增加实际 causal rate 表。移除四段 observer 后，原 driver 全文保持。

`prepareSpatial(mesh, phase, actualAdvector)` 保留同一 actual advector、momentumGeometry、原 layout/D_I 和 owning FirstJet source。原 canonical atom 是唯一相率 owner，primal/dual 仅作真实索引求和。ordinary 液/气分别积分原 normal polynomial；zero 使用 actual velocity 正负部分选择原独立两侧流入相，不读取几何空间平均。没有第二套速度重构、scalar 副本、库存、bubble 或质量修复。

使用原 Core 的 endpoint-relative signed/positive polynomial primitive 和 same-source trace refinement。相邻估计变化按各相工作尺度控制在 1e−13，最多四轮；所有有限性/来源/预算失败拒绝，缺少外部或 selected 流入相 port 的 zero 也明确拒绝。原 A_f*u_f/D_I 总量及 H10j0 运动学 divergence 不变。closure 仅核验，不调整液/气通量。误差估计没有冒称严格 bounds。

## 数值与来源证据

标准 local-benchmarks Release，OMP=1。最终必要 CTest **8/8通过**：architecture、自检、目录及 phase/material/momentum/dual/KEP 契约；文档完成后再执行三项架构/目录检查。

现 phase 契约原检查完整保留，新增 88 场景 / 14,080 atom：XY/RZ × uniform/2:1 的普通平面、圆/球、纯相、16 零迹位置/方向的四种 actual velocity，以及宽度 2⁻⁴⁰ 的小相。六个零净面保持两个非零相率，68 个小相按自身实际相率尺度核验。再补八个 XY 斜切、重合直段接曲线、多分量场景 / 1,280 atom；独立解析积分最大误差 5.56e−17。真实 dual.sides 的相邻 owner/反号 incidence 和实际 atom addback 同时核验。

| 独立最终研究核验 | 实测 |
| --- | ---: |
| full/selected 布局 / 实际 factory | 64 / 128 |
| actual atom / primal+dual view 表 | 228,600 / 315,620 行 |
| shared interior dual view 表 | 211,080 行 |
| 独立液/气物理积分最大误差 | 2.304e−15 |
| 所有 atom/view 的总量 closure | 1.946e−16 |
| 坐标重组的实际 view addback | 1.847e−16 |
| 独立 affine total / 原 A*u | 2.873e−16 |
| 原 16 零迹来源 / actual revisions | 16 / 64 |
| 零迹 atom / 重合 atom 记录 | 10,240 / 480 行 |
| 独立 backward right-time 因果相率误差 | 1.508e−16 |
| 独立零迹 actual polynomial 系数误差 | 1.555e−15 |
| actual source refinement 最大变化 | 6.159e−15 |

128 工厂分别消费原 stationary plane 和由正式 `advanceFirstJet` 产生的真实移动 circle/sphere midpoint；完整 step/begin/end/position、原 scalar/actual advector 共享均保存。独立后处理仅读实际 CSV 坐标、声明 affine 场与解析形状，自行重算闭式相交和 XY/RZ 物理积分，不读 kernel 的 expected 值。primal/dual 的几何 atom 区间在独立脚本中重组，检查无缝隙/重叠及原数量；全部行核验有限性、容量、closure 和来源。

原 56 个不可识别 tangent 仍仅 n4 selected2/3，所有行参与 finite/sharing/source；只有它们的二维 affine 精确参考不适用，full16/n≥8 selected 无豁免或新缺口。

H10h 的 XY AMR n16 pure gas dof103，由正式三层 initial/midpoint/end 的同一 Q4 和 actual uniform advector 消费新工厂：三层共 21 段液率**严格0**，气率为独立 .4*A，三个 signed outward sums 为0。不同 geometry/advector revisions 保留；本轮没有生成材料 F/R，不能把该 phase 证据称为材料更新或 PDE PASS。

wrong geometry/full clock/layout、同 layout 的不同 D_I、selected/exterior zero 流入、source 寿命/copy、FP 环境、真实改坐标/regrid 均检查明确拒绝与原结果保持。原始 Initial/Midpoint/End 声明和 regrid 拒绝保留；未实施 Hydro 的 regrid transfer 或正式外部 inflow boundary port。

128 工厂按实际 atom 查询一次，未再次准备 advector SVD/registry；最终实际 phase 工厂耗时合计约 4.134s（本机单次观察，非性能门槛）。查询/refinement/普通与重合范围逐布局保存；不声称含所有 Core 查询的端到端成本严格线性。H10j0 的局部存储和一次排序成本声明保持。

## 回归、失败记录及后续边界

开始前独立核验 H10j0 193 payload / 97 canonical 和 i2/f/e/d/c 前置归档。源码/输入全库存、implementation-before design、实际 executable/DLL、编译源码、退出记录及独立核验冻结于 [H10j baseline](../../runs/baselines/H10j_20261004/README.md)。最终运行均使用该最终 executable，原 source/affine 23 份 CSV 字节相同；所有 cases 配置、旧 App/PDF、Core/FirstJet/advector/F/R/uR/KEP/D_I/C 保持。

补充斜切 fixture 的首个非二进制端点被正式 publication 拒绝“undeclared active boundary”；改为明确且可精确表示的斜切边界接触点后通过，没有修改生产出版条件。第一次独立 zero 后处理把常速度误赋给另一个零速度分量，修正独立公式后通过。失败记录与初轮研究保存；终值未靠放宽门槛。首次短运行的 executable 路径输错，在创建运行目录前拒绝，随后使用标准 CMake 的实际路径。

原 H10b 28窗口/84gate、14失败历史和 point reference 保持，内部原 affine 审计仍完整，未宣称新最终 operator 通过。新 phase owner 尚未成为当前 H08/H09 consumer 的唯一输入。

调用库存证明旧 evaluate 由 phase/material/dual/KEP 四份既有契约和研究 driver 使用，一次全部删除会突破原 2-test 预算。因此 [H10k](H10k_Material_Spatial_Owner_Migration_TODO.md) 先以 2 production / 2 tests 迁移唯一材料 owner、actual advector stamp 和 phase/material fixtures，显式保留未迁移旧消费者；H10l 再迁移 dual/KEP 的既有 fixtures 并删除最后旧入口和 quarter 数值转移（预核预算≤5 production/2 tests）。不得在其中隐含修改最终 convection。此处修正 H10j 草案“下一任务立即删全部旧入口”的时序，不改变用户四硬约束或生产路线。

随后独立联合设计最终 native spatial moments 和 material KE/skew，内部 affine counterexample 先通过，再原完整 H10b；不把 phase/material 输入改善当成 H10 关闭。Pure FirstJet 完全独立，无库存/bubble，quarter/严格材料 commuting/R=0 只诊断。
