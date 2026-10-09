# LAMP 架构规范

**版本：** 2026-09-30  
**目的：** 保证 V2 方法落地时程序保持模块化、可测试、可维护，避免数值逻辑重新散落到 App、IO 和 case 中。

**2026-10-03 输入现状：** 独立 Interface 使用 [G04 五文件 schema](api/G04_Interface_Dictionary_Schema.md)。App 汇总 run/geometry/mesh/fields/interface 的已验证输入；初始/给定数据不按 problem 名硬编码，时钟/网格/域不放模块字典。下述 p3b 为历史验收记录，不代表旧 Benchmark 布局仍可运行。

**2026-10-02 Interface 主线：** [纯 Hermite–GALS](task/I12_Pure_Hermite_Gals_Design_Change.md) 使用现有 HermiteGals state/identity 和唯一 snapshot，以显式 FirstJet 模式区分历史 ConservativeVolume。新状态无 M，几何体积为只读派生诊断，现有 owner 表中的 M/Q_l 只适用于历史守恒分支；不能复制 geometry owner 或在失败后切换路线。I12p3b 已验证独立规定速度 App：`Interface.Benchmark.method hermiteGals`，见 [API](api/I12p3b_FirstJet_App_API.md)。App 持有唯一 accepted mesh/publication/clock；I12p3c 已接通同源只读[正式输出](api/I12p3c_FirstJet_Output_API.md)，输出只持有派生诊断，不拥有第二份几何；[I12p3d](task/I12p3d_FirstJet_Performance_Scale_Report.md) 已补齐被动阶段计时、三档重复规模及动态结构证据，[Gate I 复核](task/I12p3r1_FirstJet_Gate_Reaudit_Report.md)发现的字典诊断已由[I12p3e](task/I12p3e_FirstJet_Input_Diagnostics_Report.md)关闭，[I12p3f](task/I12p3f_FirstJet_Curvature_Diagnostics_Report.md)也已补齐同源曲率时序诊断，[I12p3r2](task/I12p3r2_FirstJet_Gate_Final_Report.md)已裁决Gate I PASS，下一项H01；当前独立路径尚未接Hydro。

**当前 Hydro 状态（2026-10-07）：** [H13 同源Hermite毛细载荷](task/H13_Hermite_Capillary_Traction_Report.md) **FAIL，仍为唯一当前任务**。原静滴七组伪流例外、无黏波能量/频率未收口；所有隔离载荷参考均未进入生产；Gram族及等密度1.557阶的简单加性候选均未采用；三个实际观察确认当前常权应力平均的信息损失；局部应力偏差及已知压力lifting均未采用；后者的等密度负刚度、0.716阶与能量升高已由原native/source分解复现；二相传输81组及九次FirstJet方向更新已完成；共享一圈应力矩候选末档频率阶-0.435、原能量峰值缺陷随加密增至0.452，独立矩阵指数/功核算复现，已拒绝；18组光滑传输及原实际算子复现确认native压力/加速度近二阶、最终source界面运动约一阶，固定二次应力profile预测后余量二阶以上；局部normal最小代价消差虽可达，但邻段扰动1.506倍，逐段叠加留下全部误差，已排除；完整P2压力矩差=tangent+jump/h已由原全矩阵核验，实际可达jump端口相对尺度约0.335且不随加密消失；仅改弱牵引虽频率二阶但能量缺陷约5.56；假设面通量矩匹配运动/牵引的线性参考频率1.987阶、原能量舍入闭合；具体共同source与P2弱牵引联合改造已获用户明确授权；真实source候选契约通过，完整弱牵引及全物理门槛待完成，候选未验收；一般曲面/AMR/实际动态及联合精度仍待解决。定位发现[直接阻断H08a的实际Jacobian反例](task/H13_H08a_Transparent_Picard_Jacobian_Counterexample.md)，已按用户限定例外在一个现有文件中先代数消去C1内部边界相消项；25/25受影响检查通过，新版本完整H08a的28窗口及原空间/时间/R历史门槛通过；H13单周期复测仍失败。r11其他来源与原门槛保持，不另拆精修任务。历史H08a/H10/H11/H12规定范围结果保留；这一处限定正确性修复已收口。Pure FirstJet独立，无M_K/bubble、R=0或quarter commuting强制；唯一共享F、真实R/uR、KEP/skew保持。H14未启动。

## 1. 依赖层级

**2026-10-03 Hydro V2 现状：** [H01](task/H01_Native_Face_Cell_Dof_Layout_Report.md)已完成原生面速度/单元压力的紧凑不可变[布局与数据 contract](api/H01_Cch_Native_Dof_Layout.md)；[H02a](task/H02a_Cch_Build_Boundary_Report.md)将其 build 实现迁入独立真实库 `lamp_hydrodynamics_cch`，直接依赖 Core，不链接旧 Hydro/Interface Compatibility。[H02b](task/H02b_Integral_Divergence_Report.md)已在该库实现唯一只读积分散度 D_I=B A_f，真实子面与 RZ 测度闭合通过，H02 PASS。[H03](task/H03_Negative_Transpose_Pressure_Load_Report.md)已保留同一权威散度并直接生成 C=-D_I^T，最大压力功伴随缺陷7.39e-16。[H04](task/H04_Single_Phase_Saddle_Report.md)已在此库完成封闭域单体参考求解、normal velocity 消元及逐连通域真体积压力均值，均匀 XY 连续 MMS 通过。A_u由物理owner提供，矩阵、factor和结果均由typed CCH对象只读持有。Cch / Detail/Cch来源有独立编译层与唯一owner，不能通过旧Hydro header/link绕回兼容层。[H05](task/H05_Rz_True_Measure_Axis_Regularity_Report.md)保留同一D_I并以真V_c归一化，验证零面积轴面的零径向速度和奇偶规则，无新primary ghost或第二散度矩阵。[H06](task/H06_Energy_Viscosity_Report.md)由局部native P2 strain和非负物理求积构造K_mu，轴奇偶/hoop正则、PSD、无黏及算子收敛通过。[H07a](task/H07a_Same_Geometry_Phase_Rate_Report.md)已绑定同源相测度/瞬时率，[H07b](task/H07b_Phase_Quadrature_Viscosity_Report.md)从同一FirstJet Q4生成相内正权体积求积并接入H06黏性，H07关闭；[H08](task/H08_Material_Flux_Continuity_Defect_Report.md)绑定三层几何质量、唯一中点材料流量与实际非零连续性缺陷。[H09a](task/H09a_Momentum_Dual_Geometry_Report.md)提供native速度的真实dual交集/共享段及同一Q4材料质量，[H09b](task/H09b_Dual_Material_Consistency_Report.md)从唯一H08流量作局部线性几何转移，保留真实dual缺陷/几何转移差异，提供同层uR积分项和单列材料动能项，两项共同关闭H09。prepared对象持有来源/clock/layout/原D_I，无新库存或全局自由度，不更新接受状态。下一项H10；尚无可运行CCH App、完整时间推进或流体联合验收。

固定依赖方向：

```text
Core
  ↑
Solvers
  ↑
Modules/{Interface,Hydrodynamics,Electricity}
  ↑
Apps/{Interface,Hydrodynamics,Electricity,Electrohydrodynamics,Electrospray}
```

允许 `Solvers` 依赖 `Core`；Modules 可依赖 `Core + Solvers`。App 负责跨模块接线。

**禁止：**

- `Core` 认识具体物理模块或 App；
- 数值 kernel 读取 case 路径、字典或输出目录；
- Electricity 为方便直接依赖 Electrospray；
- Interface 内部调用 Hydro 的时间推进；
- App 复制一套 NS、Poisson、interface transport 算法。

跨模块数据优先使用只读 view/value object；如果必须建立共享 contract，应放在中立层或 App/Coupling，不以实现类互相穿透。

## 2. 对象角色固定

新增代码遵守以下命名和职责：

| 类型 | 职责 | 不允许 |
| --- | --- | --- |
| `Config` | 已解析、已校验的不可变配置 | 保存运行时状态、打开文件 |
| `State` | accepted/candidate 物理状态 | 重新解析字典、写输出 |
| `Operator/Assembly` | 离散算子和 prepared data | 管主时钟、决定输出 |
| `Result` | 一次求解的权威结果 | 暴露可写内部场 |
| `Report` | 诊断、残差、迭代统计 | 参与数值分支 |
| `Writer/IO` | 写文件和格式化控制台 | 推进方程、改变接受条件 |

一个运行过程只有一个主时钟。candidate/accepted state 的所有权必须显式区分。

## 3. 唯一 owner

**2026-10-03 材料路线更新：** [完全非守恒材料耦合](task/H07_Pure_Gals_Hydro_Coupling_Design.md) 不设 Interface 或 Hydro 守恒材料库存。下表 `M/Q_l` 守恒状态仅归历史比较；当前几何相测度由 Interface 唯一 publication 提供，Hydro 只读消费并生成绑定 geometry/time 的瞬时材料率及后续时间积分结果。正式 `quantum_interface` 可被 CCH 只读消费；禁止进入 Compatibility 或反向依赖 Hydro。App 仍拥有唯一接受事务。

下列对象在 production 路径中只能有一个 authoritative owner：

| 对象 | Owner | 其他模块如何使用 |
| --- | --- | --- |
| Hermite geometry snapshot | Interface | 只读 geometry handle/view |
| cell phase volume `M/alpha` | Interface | Hydro/Electricity 只读消费 |
| phase-face flux `Q_l` | Interface | Hydro 直接消费，不重算 |
| native velocity `u_f` | Hydro | Interface tracer 使用明确时间层的只读速度 view |
| divergence `D_I` | Hydro topology operator | pressure load 由它生成 |
| mass flux `F_rho` | Hydro transport | momentum 与 diagnostics 共享 |
| electrostatic skeleton/operator | Electricity | charge/current/Maxwell 共享 |
| panel-integrated charge `Q` | Electricity | EHD 只接收结果/traction |
| accepted-step transaction | App/Coupling | 各模块只准备 candidate，不各自私自 commit |

任何新代码想创建第二份同类 production 对象，必须先停下来说明为何现有 owner 无法承担职责。

## 4. V2 路线与 legacy 路线隔离

### 4.1 生产路线

- Interface：纯 Hermite–GALS / `FirstJet`；`conservativeHermiteGals` 保留为历史显式比较；
- Hydro：`cch` native staggered；
- Electricity：`sharpHHO2`。

### 4.2 Legacy/experimental

`parametricGals`、大型 CompatibleMoment、`sharpHermiteP3` 等只允许：

- A/B benchmark；
- 旧 case 回归；
- 为新方法提供可比较输入/结果。

禁止：

- production 新路径失败后自动回退；
- 新旧路径共享可写 state；
- 为保持旧测试通过而让新方法继承旧算法内部约束；
- 在 legacy 类中继续加入 V2 新功能。

## 5. 文件布局原则

2026-10-02 用户授权的目录整理将头文件统一放入 `include/`，实现翻译单元放入 `src/`。内部实现头位于所属层/模块的 `detail/`，不因此成为公开物理接口；消费者使用层/模块限定 include，不用 `../src/...`。`cases` 只保存配置与静态输入，运行快照/日志/结果统一写入 `runs/<算例>/<时间戳_标识>`。固定运行、构建与调试入口见 [07 / 目录与运行指南](guides/Run_And_Debug.md)；路径迁移与原内容身份单列，不重写历史冻结证据。

不新建新的顶层目录。V2 文件继续落在已有模块下，并以清楚前缀区分职责。例如：

```text
include/Modules/Interface/
    HermiteGals*.h
    ConservativePhaseFlux*.h

include/Modules/Hydrodynamics/
    CchTopology*.h
    CchMassTransport*.h
    CchMomentum*.h
    CchViscosity*.h

include/Modules/Electricity/HHO/
    HhoGeometry.h
    HhoBasis.h
    HhoLocalOperator.h
    HhoCondensation.h
    HhoFluxRecovery.h
    HhoElectrostaticSolver.h
    SurfaceChargeSchurSolver.h
```

具体拆分随任务进行，不提前创建空壳。`.h/.cpp` 职责必须单一；不得新增 `Utils/Common/Manager` 式大杂烩。

## 6. API 规则

1. 数值函数参数中使用类型化对象，不传 `DictionaryInput::Node`、case path 或裸字符串 method name。
2. `initialize` 只初始化；`prepare/assemble/solve/advance/evaluate/commit` 按真实语义使用。
3. prepared operator 的 cache key 必须包含影响矩阵的 geometry/topology/material/boundary identity；不能只看矩阵尺寸。
4. rejected candidate 后，任何 cache 若依赖 rejected geometry/state 必须失效。
5. interface/flux/traction view 必须带 geometry identity 和 time layer；不靠注释猜时间层。
6. numerical kernel 不调用 writer；writer 失败不改变已经接受的物理状态，但必须让 App 报告 IO failure。
7. 异常不得直接穿越 OpenMP 并行区；并行段汇合后统一传播。
8. 结果对象默认只读；兼容视图不得反向参与 production 计算。

## 7. 代码卫生

每个任务完成前必须满足：

- 源文件首先 include 自己的接口；
- 头文件自行声明依赖，不依赖偶然 include 顺序；
- 无新增 compiler warning；
- 无未使用 production 分支、临时 debug 开关、case-specific magic number；
- 无复制粘贴的第二套单位/几何/边界转换；
- 失败路径有明确错误消息，包含模块、阶段、关键 identity；
- 删除已被本任务真正替代且没有 legacy 价值的 dead code；否则标明 legacy ownership，不“半删半留”。

## 8. 单任务变更预算

Codex 单任务默认预算：

- `<= 5` 个 production 文件；
- `<= 2` 个测试文件；
- 一个数学目标；
- 一组明确验收。

若实现前就确定会超过预算，先拆为“数据结构/API任务”和“数值实现任务”；不得边做边扩张范围。

任务禁止同时包含：

- 目录大搬迁 + 数值格式修改；
- 新算法 + 性能重写；
- 新算法 + 大规模字典迁移；
- 修三个不相关模块；
- 修改容限来配合新结果。

## 9. 测试位置与边界

- operator/unit test 验证代数性质，不依赖完整 App；
- module benchmark 通过公开模块 API；
- App test 只验证接线、时钟、字典、输出和跨模块 transaction；
- production case 不作为 unit test；
- 测试可使用 analytic/mock velocity/geometry，但必须明确它验证哪一层，不得把 mock 结果冒充完整耦合验收。

## 10. 与旧架构文档的关系

[QUANTUM_Code_Organization_Standard.md](history_todo/QUANTUM_Code_Organization_Standard.md)仍用于查询既有文件、兼容入口和历史职责。本文件是 V2 新开发的简明硬规则；两者冲突时，**新代码以本文件为准，旧代码按任务逐步迁移，不做一次性全仓重写。**

G02 更新：当前 Interface 只注册 FirstJet App；旧 Hydro/EHD 依赖归 lamp_interface_compatibility，其 71 条既有 include 边逐条冻结。新消费者不允许目录级豁免。头文件归 include，翻译单元归 src，测试辅助头仅归 tests/support。工程细则见 [代码与测试](06_LAMP_Code_And_Test_Rules.md)。
