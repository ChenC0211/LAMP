# H10i1 / 重合相迹与裁剪测量的同源契约设计

状态：**PASS（仅数学/API契约设计）/ engineering PASS**，2026-10-04；[报告](H10i1_Coincident_Phase_Trace_Clip_Contract_Report.md)。定义了clip内部迹、材料因果迹、exact normal Taylor认证、重合边界唯一归属及完整分阶段文件预算；本轮0 production/0 test，实际source能力仍FAIL。16场景：8个dual位置publication成功但mass拒绝、8个native位置publication拒绝；下一唯一任务[H10i1a Core trace port](H10i1a_Native_Quartic_Trace_Port_TODO.md)。H10i的132 payload/32当前交付已独立核验；上一goal turn为progress。H10/H11/H08a顺序保持。

目标：对同一保留Q4在真实段上恒为零的情况，明确**裁剪区域的几何边界迹**与**共享材料输运的面迹**各自含义，设计可认证的只读来源接口及后续实施清单。只解决这一测量契约事实；不在本任务实现新相通量、dual转移、最终convection或时间推进。

允许修改：本TODO、实施前设计/报告及当前入口状态；必要时仅扩展既有`cases/benchmarks/cch_convection/driver.cpp`的被动来源观察。**0 production / 0 test / 1 research**。可读取Core/Interface正式源代码核验契约，不能把其Detail实现复制到Hydro或脚本里。临时证据归`runs/checks/H10i1/<时间戳>`，关闭证据归baselines。先写设计再观察，每项后续实现须单独TODO，≤5 production/2 test，超出则按接口与消费者拆任务。

材料耦合仍仅四硬约束：唯一共享F_rho、内部面反号共享/真实AMR子面求和、来源一致、KEP/skew。Pure FirstJet保持独立，无库存/bubble、逐cell体积拟合或质量修复。禁止改变primary/Q4 publication、D_I/C、原H08/H09/F/R/uR、App/字典、冻结PDF和原H10b门槛；不得以微扰界面、epsilon采样、任意half分配、材料平均或特殊AMR分支处理重合面。未闭合前保留明确拒绝。

## 设计必须先回答

1. 几何裁剪：如`phi=x-a`，矩形`[a,b]`内部为气相、`[c,a]`内部为液相，正长度零边界不改变其体积。边界积分须使用该裁剪区域的单侧内部相迹/取向，并处理free interface与fixed clip边界的重合归属，不能双计。说明一般Q4零段、横穿/相切、端点与chart分片、轴和selected区域的认证或明确拒绝条件。
2. 材料面迹：在正测度零段上，严格`phi<0`和`phi>0`不能直接覆盖全物理面。需从离散输运模型推导一个共享且可加的规范，并说明正/负/零法向率、两侧相迹和时间crossing的意义。不得把裁剪体积的内部迹直接当作两邻域各自材料通量；共享面只有一个值。候选中心或单侧规则须有数学理由，不能为通过反例任意指定；空间迹也不冒充H08a时空积分。
3. 来源与认证：公开typed测量port须保留原geometry identity、native chart、方向、零区间和必要单侧符号/normal jet认证；不得从diagnostic rounded local dofs重建第二Q4。审计现有Core/Interface是否有足够精确的restriction/零迹/单侧符号能力；缺失就拆正式接口任务。owning view生命周期和同一geometry消费者的来源拒绝须明确。
4. 对普通短相/近零但非恒零trace，保持原endpoint-relative ranges与严格来源，不能扩大零阈值以吞掉小相。真正无法认证的退化Q4明确拒绝并计入能力限制；不能把所有重合平面继续列为不支持以收口。
5. 后续消费路径：同源port分别服务真实clipped phase region/momentum mass与共享phase face-rate，声明唯一owner与无循环依赖。列出Core接口、Interface派生测量及Hydro消费的实际文件和测试预算，必要时分开实施。H10j现有2 production预算不能隐含夹带Core/Interface修改。

## 固定验收

T0：标准Release、architecture/self-test/repository layout与受影响短契约通过；本轮protected生产/测试/App/PDF及原研究数值路径保持。H10i最终归档先独立核验。

T1：数学契约覆盖H10i的XY/RZ uniform/2:1四个实际反例。后续实施必须让这些publication的真实dual质量成功，按独立平面裁剪参考核验液/气体积非负、capacity与物理XY/RZ积分≤1e−12实际测度；两侧纯相质量正确且无重合边界双计。补原生面重合、非仿射Q4零迹、相切/短相/近零非零、正负零速度、axis/selected/AMR分段，以及来源与生命周期拒绝。一般退化范围与明确不支持项须事前列出。

共享材料面规则须同面反号、分片加回、同源物理测度闭合≤1e−12 gross，仍可定义唯一F_rho及实际R；KEP/skew与R/uR使用该实际来源自洽，不新增R=0或pressure→dual commuting硬等式。几何clip迹与材料迹的区别须用独立闭式反例核验。时间跳跃或大R如实记录，归后续空间精度/H08a验收，不能用相迹选择隐藏。

本设计任务PASS只表示上述规则、来源port与分阶段实现预算闭合，不表示能力已实现。实施后必须保留H10i 64布局常规测度/GL核验、原phase/material/dual/KEP回归及内部affine反例；原完整H10b仍在最终空间前置通过后复验。任何缺少认证、归属、共享或必要数学证明的设计均FAIL并列明确接口任务；内部困难不标BLOCKED。

交付：单一契约设计/报告、旧→新port与消费者表、真实反例证据、固定实施文件/门槛、平面SHA256 manifest及当前唯一下一任务。H10i历史FAIL保持，关闭缺口后的复核单独记录；H10j及最终native精度均不能提前报PASS。
