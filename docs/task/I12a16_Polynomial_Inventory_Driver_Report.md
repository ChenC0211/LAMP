# I12a16 / Polynomial 库存驱动接线报告

**状态：FAIL（本项）。I12 / Gate I 仍 FAIL。** 原 a15r4 的积分算子仍 PASS；本项完整驱动在真实物理入流面缺少相历史/边界分支证书，4 个 first-step 场景均未推进 accepted 时间。见 [TODO](I12a16_Polynomial_Inventory_Driver_TODO.md) 与 [驱动说明](../api/I12a16_Polynomial_Inventory_Driver.md)。这是已定位的内部实现缺口，不是外部依赖，不标记线程 goal BLOCKED。

## 已实现并执行的接线

只修改原 benchmark 与 tests/CMakeLists，0 production+2 test-related。显式 polynomialIntervalInventory 路线从原 Case 的 rotation/vortex 端点速度构造同一 tensor-cubic source：rotation 是一次空间场，vortex 是三次空间场、cosine 两端的线性时间重构。程序核验端点值与四个导数的原表达式；没有伪造恒定速度 capability。

end 与 midpoint prediction 都使用同一旧 accepted scalar、owned velocity 与 actual RK2 characteristic owner。实际 midpoint candidate 身份通过原 a12 polynomialArrivalView 提供。两阶段预测可独立于库存积分先完成，之后才允许库存、闭合和发布。四个场景的实际两阶段预测均通过：uniform 每阶段 81 个 regular vertices、144 条 edges、64 个 centers；AMR 为121、232、112，初始界面确实切过2个coarse和13个fine cells。

新 advance 函数已连接所有 native face → bounded polynomial transfer → 原 integrated inventory assembly → reflux → 原 finish/I10/I11 两阶段发布；实际执行在首个面拒绝，因此**不能声称本项已经验证库存装配或任何完整新运动步**。旧 producer 仅供独立旧命令运行，不是失败 fallback。边界 first-jet 的原显式条件仍由原 predictor 处理，它不能代替缺失的整个面—时间相历史。

## 原步骤的真实拒绝

smoke 的4个场景采用原 n=8、cx=.65、cy=.53、radius=.13；rotation speed=2、dt=.01，vortex speed=40、period=.4、dt=.001，各跑1步。另实际运行 full 命令，请求原 rotation30步、vortex400步；四者均在第一步相同位置拒绝。

| 场景 | 首个半步时间末端 | 失败面及邻接 | E | 预算 |
| --- | --- | --- | --- | --- |
| rotation uniform/AMR | .005 | native0，vertical，negative=-1、positive=0 | 5.4687500000000265e-4 | 1.5625e-8 |
| vortex uniform/AMR | .0005 | native0，vertical，negative=-1、positive=0 | 4.8705679135140431e-4 | 1.5625e-8 |

面预算始终为关联最小 cell capacity*1e-6；max panels=4096、depth40。所有失败均为 `polynomial_half_integrated_flux`、`panels=4096 pending=4096`，根因 `RK2 enclosure cannot exclude initial boundary hit`。四个 accepted time 仍为0，原 runCaseUsing 的拒绝检查确认 accepted identity 和每个 M 未变；没有发布新状态、填充气相通量、clipping 或提高64eps修复限制。

它不是 a15r4 剪切的积分精度问题：这里全部参数块都缺少有效的边界分支/历史，单纯细分不能把存在的物理边界碰撞变成无碰撞。原旧恒速 inflow source 绑定常速 translation，不可借用它替本 polynomial provider 建立边界来源。尽管这些解析验证算例的液滴远离域边界，也不能据此越过来源契约直接给所有边界面写零。

首次 smoke 3.7075645 s、exit1；完成 midpoint 先行预测记录后，最终 smoke 3.8028988 s、full 3.7553884 s，均exit1、completed0/planned4/failures4。日志中的 vortex `cellVolume=nan` 是原 harness 明确未提供逐cell vortex 解析参考的标记，不是本次边界拒绝导致的状态NaN；位置/梯度/曲率等旧 accepted 参考仍有限。

## T0 与归档

首次 configure/build 成功；最终 benchmark 重建4.5161314 s、exit0、无 compiler warning。其余38个 runtime executable 保持 r4 SHA，原 production 完全未改。新增 `interface_polynomial_inventory_driver_smoke` 是 ordinary CTest，未设 WILL_FAIL。

26项相关CTest用时161.198711 s：**25旧测试PASS，1新driver FAIL**。原25项完整输出（包含architecture、a15全部57个积分和r4子区域检查）逐字相同。架构仍12 targets、383 TU、908 files、3118 includes、5 exact legacy edges。旧 crossing probe 的exit1与输出保持。

81 frozen文件、3份方法PDF、r4 archive身份保持。第一次审计遇到PowerShell将ConvertFrom-Json数组保留为单个外层对象的计数问题；显式枚举后重新核验39个runtime，38个未改程序SHA完全一致，原失败记录保留。它只影响审计记录，不改变测试或构建结果。最终审计完成，但task状态仍按真实步骤结果记FAIL。

证据：[目录](../../runs/baselines/I12a16_20261002/)、[验收](../../runs/baselines/I12a16_20261002/verification.json)、[真实driver缺口](../../runs/baselines/I12a16_20261002/driver_verification.json)、[旧回归](../../runs/baselines/I12a16_20261002/regression_comparison.json)、[范围与runtime审计](../../runs/baselines/I12a16_20261002/source_audit.json)。归档交付两个测试文件、TODO/report/命令说明/入口、39程序及原DLL，并逐项核验payload SHA。

## 下一项

为原 polynomial velocity 增加显式 physical boundary phase history，以及可验证的 interior/exterior characteristic 覆盖。边界历史与旧 snapshot/velocity/time extent保持真实所有权；即使声明全边界气相，也须证明真实旧内部可达区域的相/库存范围，不能把几何气相等同于 M=0，不能用外部速度延拓或解析内部重初始化绕过来源约束。

先解决同一四个 first steps 的全部 native faces，再检验同一驱动的严格库存界、几何预算、全程rotation/vortex及时间/空间/AMR收敛。I12及后续Hydro均未通过，继续原路线。
