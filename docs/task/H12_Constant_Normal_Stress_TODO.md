# H12 / 解析常量法向应力载荷

2026-10-06，实施前登记的唯一任务；现已按原范围 **PASS**，见[报告](H12_Constant_Normal_Stress_Report.md)。前置 H11 的 161 个冻结文件、1387 个源码/输入及 19 份现行交付已逐一复核，完全匹配。检查目录：`runs/checks/H12/20261006_051700`；上一 goal turn 为实际进展。总体 Hydro/EHD goal 保持 active，H13 本轮不启动。

## 唯一数学目标

按冻结 [CCH 方法](../methods/QUANTUM_CCH-EHD_Method.pdf) §9.2、性质9.1及 V9/V11/V13，给定常数法向总应力 J0，只读消费同一 Pure FirstJet 派生的 `alpha=V_l/V_c` 和现有唯一 `C=-D_I^T`，生成 `F_Gamma,n=J0 C alpha`，使 `p=-J0 alpha` 与载荷在舍入误差内平衡。浮点装配采用等价的 `C(J0 alpha)`，直接调用现有压力算子，不再装配梯度或几何法向沉积。

必须进一步把真实载荷接入原封闭域单体鞍点系统，求出接近线性/舍入误差的速度及正确压力跳跃；不能用解析零速度替代求解，也不能仅靠两份同公式相消作为验收。此处 J0 是显式给定总应力，不计算 Hermite 曲率或电场；将来一般法向应力仍归同一载荷 owner。

## 允许修改及预算

- 新 `CchInterfaceNormalLoad.h/.cpp`：一个只读压力相容法向载荷 owner，常数工厂、来源绑定、验证和原生载荷；生产头在 include，实现位于 src。
- 根 CMake：将新实现纳入现有 `lamp_hydrodynamics_cch`，从旧 Hydro glob 排除；不建新求解库。
- 如真实零速平衡揭示原鞍点线性求解精度不足，仅允许在 `CchSinglePhaseSaddle.h/.cpp` 对**同一增广系统、同一 factor**进行有界残差校正并记录次数；不得改矩阵、物理方程、压力 gauge 或任何原残差门槛。必须先保存具体失败。上述总生产预算最多5文件。
- 两份测试文件预算：一个新短契约 `cch_constant_normal_stress_contract.cpp` 及既有 tests/CMakeLists 注册；不修改旧数学断言。架构归属沿用既有 Cch 分类。
- 一个 `cases/benchmarks/cch_constant_normal_stress` 标准研究驱动/静态 descriptor/README及研究目标注册；独立参考和数据导出只属于 tests/cases。临时审计归本次 runs/checks。
- TODO/API/报告/路线图和现行索引；不覆盖旧报告或基线。

## 不变项

Pure FirstJet 与唯一 geometry、相测度、H08a 时空积分、唯一 F、真实 R/uR、H10/H11 对流、H06/H07 黏性数学及接受事务保持。r11 source/characteristic/Jacobian 基础设施继续封口；无库存/M_K/bubble、R=0、quarter commuting、几何修复、最小黏度、特殊 case/AMR 生产分支。不得把恒定应力专用解析答案写进线性求解器。

## 离散、符号和边界契约

载荷保存同一 `CchPressureLoad` 与 `CchPhaseGeometry`，核对其同一 D_I owner/layout，保留完整 geometry identity 和 `{step,begin,end,position}`。不复制另一份界面或材料状态；alpha 始终来自保留 phase object 的只读测量。J0 必须有限；非法 mesh/layout/D_I/source/time/regrid/nonfinite 明确拒绝。复制共享不可变数据，等价内容 mesh 生命周期按既有规则验收。

载荷是上述 `Cp+F=0` 中的**平衡项**。原鞍点契约为 `A_u u + Cp = r`，静态检查传入 `r=-F`，禁止混淆符号。参考 `A_u=M/dt+K_mu` 使用既有实际 dual 几何质量和分相能量黏性；不从解析压力构造 RHS。

全原生矩阵的外边界行一般 `C 1 != 0`。无 gauge 的 `p=-J0 alpha` 检查完整行相消；原封闭域求解则在每个连通域取真实体积均值零，参考为 `p=J0 mean_V(alpha)-J0 alpha`，自由面满足 `Cp+F=0`，固定 normal velocity 的剩余载荷为显式壁面反力、功为零。RZ轴零面积载荷及径向速度严格零。不得通过删除外边界 DOF 或强迫边界反力零来通过。selected/断开域的检查属于代数边界契约，不代表润湿、接触线或一般开放边界已经实现。

## T0：构建、来源与回归

标准 CMake Release；architecture、architecture_self_test、repository_layout通过。新短契约、原 H03/H04/H05/H06/H07a/H07b/H09a 受影响检查通过；必要的既有零 RHS、MMS、压力伴随、真实质量/黏性和拒绝检查保留。逐文件证明 H08a/r11/H10/H11 数学实现与输入未改；仅链接变化不触发再次重跑无关的102分钟H11研究。保存实际源码/配置/可执行/库哈希和退出码。

## T1：独立算子与真实静平衡

- 覆盖 XY、含轴及离轴 RZ；uniform、真实2:1/嵌套、selected含洞/断开域；平面与曲面、pure两相、切单元及粗细接口；正/负/零 J0。
- 从真实 native face/cell 拓扑、外向号及独立 XY/RZ 物理面积重组载荷，参考不得读取 C 的系数；alpha来源与保留 publication 的相体积逐项核对。共享段只引用真实子面。工作尺度相对误差≤1e-12；全行无gauge平衡及离散功关系同样≤1e-12。
- 用实际 dual质量和分相黏性组装原 `A_u`，原鞍点求解 `r=-F`。原四个线性/物理残差门槛保持；自由面 `Cp+F` 的全局工作尺度缺陷≤1e-11。速度用 `U_load=max_i |F_i|/A_ii` 归一化，max|u|/U_load≤1e-11；压力对逐连通域真体积均值参考的 max误差/|J0|≤1e-10。J0=0须返回精确零，不使用人为分母下限。固定轴/边界速度保持规定零，原始反力另行保存。
- 高密度比1/1e3/1e6，非零、低黏及零黏；检验零黏矩阵消失，不靠黏性维持静止。真实 FirstJet regrid 后显式 transfer/rebuild，新结果继续通过，旧来源拒绝；错误输入不能污染已返回结果。
- 反号或错位 alpha 的独立负控必须被载荷/压力参考发现；不能仅用零速度判对错（反号常力同样可能被错误压力平衡）。

## T2/T3：标准静滴与 AMR crossing 研究

固定界面的解析常应力研究，不冒充 H14 时间推进。四布局 XY/RZ × uniform/2:1，n8/16/32；三个位置使同一圆/轴对称球的界面依次跨过、接触及离开 refinement 边界。XY圆半径0.19、圆心x为0.40/0.44/0.48且y=0.5；RZ球的轴向圆心取同三位置，径向中心0。XY中心加密盒[0.25,0.75]^2；RZ加密盒r∈[0,0.5]、z∈[0.25,0.75]。

每个几何扫描 rho_l/rho_g=1/1e3/1e6、mu_l=1/1e-3/1e-6/0且mu_g=0.01mu_l、J0=-1.75/2.5；dt=0.125仅定义静态质量块，rho_g=1。共36个几何、432个质量/黏性块、864个静态求解。沿用 T1 门槛，各参数独立过关；不以平均结果掩盖粗细面或高密度比失败。

保存压力/速度/alpha/实际载荷/质量/反力、真实拓扑与物理测度、每次求解残差及性能；独立复算必须检查每个求解和每个参数组合。此研究验证解析常应力的静平衡与无黏连续极限；不声称 Hermite曲率误差、动态毛细波、完整Euler/NS、全耦合时间阶或生产效率已通过。标准入口仍为通用 build/run_case；临时探查只能在 runs/checks。

## 失败与交付

任何原门槛未过即 FAIL，保留实际输入/程序/原始失败。不通过放宽容差、归零速度、重贴来源、跳过接口或切换算法求通过；仅缺失外部条件才 BLOCKED。必要的原系统线性残差校正属于本静平衡接线，不拆新 source/Jacobian/低秩任务。

交付生产/API、短契约、完整静态研究、独立审计、报告及来源/程序/数据冻结。全部证据完成后才能 PASS，并转入下一单任务 H13；当前总目标未完成。

执行结果：13/13受影响检查，36短场景/108静态解/17拒绝/124负控，以及完整36几何/432块/864解与独立复算通过。原一次LU及普通残差校正的失败保留；正式改动仅在登记范围内对同一增广系统补偿计算残差、复用同一factor，原门槛未动。离轴RZ按已有含轴标量支撑+Hydro环形子域契约实现。没有重新开启r11基础设施或新增材料硬约束。
