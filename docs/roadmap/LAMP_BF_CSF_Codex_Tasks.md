# LAMP BF-CSF / Codex 执行任务书

配套唯一方法规范：`LAMP_BF_CSF_Method_V1.md`。**先阅读方法规范全文**，然后**每次只执行一个 BF-ID**。BF-01只读color及BF-02恒曲率载荷已验收，完整BF-CSF尚未验收；原 H13 FAIL 不得改写为 PASS，H14/H15 不得自动启动。

**本轮BF-02 PASS，下一次唯一入口BF-03**：见 [BF-02 TODO](../task/BF-02_TODO.md) / [报告](../task/BF-02_Report.md)。新两组契约86.7/9.6秒、200次实际saddle及10/10回归通过，完整cell/face归档。本轮未启动BF-03/H14/H15；BF-01保持封存PASS。用户已解除旧档案缺失前置；[BF-00报告](../task/BF-00_Report.md)登记当前源码/输入新A/B基线、旧缺失/FAIL与180秒timeout，原物理和数值门槛保持。完整BF-CSF尚未验收。

## 15. 可直接交给 Codex 的分阶段任务

### 15.0 严格执行规则

每次**只交给 Codex 一个 task ID**，每个任务保留：目标、允许修改、禁止修改、数学约束、实现、验收、失败规则、交付。遵守当前项目默认的 `<=5` 个 production 文件、`<=2` 个测试文件预算；超过先拆任务。状态只允许 `PASS / FAIL / BLOCKED`。**BF-00 是新方法选择的立法/文档变更，不得把原 H13 的失败改成通过。**

#### BF-00：方法切换确认与冻结基线（先做，0 生产文件）

**目标**：确定新 BF-CSF 与原 H13 的并行 A/B 关系；冻结输入、原动态失败与验收门槛。

- **允许修改**：现行数值方法索引、roadmap 中的**方法变更注记**及本设计/新任务文件。
- **禁止修改**：`src/`、`include/`、Interface 及 H13 已有结果与 FAIL 判定；不启动 H14。
- **数学约束**：承认 `C=-D_I^T`、Pure FirstJet、唯一 `F_rho`/真实 `R/uR`；明确 `c_eps` 仅是派生力场。
- **验收**：记录仓库身份/版本、现有 H12 与 H13 原报告、当前受影响 CTests、原 long-wave FAIL 和全套 gate，确认 `sigma>0` 物理 RHS 符号。
- **失败规则**：任何引用文件缺失/版本不明记 `BLOCKED`，不得以模型推测代替源码。**用户明确例外（2026-10-09）：** 已确认丢失的旧H12/H13历史payload不阻断新路线，按BF-00冻结的当前真实源码/输入建立新A/B基线，旧FAIL/缺失和全部物理/数值gate保持；不得把新运行冒称恢复旧数据。
- **交付**：`BF-00_Report.md`、新 task 入口、现有状态索引引用；0 数值变化。

#### BF-01：Regularized color 的唯一只读派生

**目标**：从当前 FirstJet publication 生成 `color_csf[K]` 与 `epsilon`/求积报告。

- **允许修改**：一个 color `.h/.cpp` 类对与至多两个对应测试/注册文件；确需 CMake 配合先在预算内核对。
- **禁止修改**：Interface 产出、`CchPhaseGeometry` 材料质量、`CchMaterialTransport`、对流、鞍点和旧 H13。
- **数学约束**：式 (3.2)–(3.6)、真物理体积、boundedness、纯相精确常值、同一 geometry/clock；不独立推进 `c`。
- **实现**：实现 `d_h` 的有限/退化审计、正权体积分、RZ/XY 和 selected 域身份绑定。
- **验收**：解析 `H_eps/delta_eps`、平界面网格位置扫掠、圆/球、RZ near-axis、uniform/AMR 只读结果；旧质量/压力等测试不变。
- **失败规则**：不允许通过 clipping、强行以 `phi` 作 signed distance 或直接 cell-center 取值放宽误差。
- **交付**：真实 `c` 数据、求积收敛、测试日志与 `PASS/FAIL/BLOCKED`。

#### BF-02：C 共用的恒曲率 CSF 载荷

**目标**：实现 `g=C c_eps` 与 `S_sigma=sigma*kappa0*g`，建立准确符号和 BF 代数门槛。

- **允许修改**：一组新 CSF load `.h/.cpp` 和 1–2 份独立合同/注册测试。
- **禁止修改**：现有 H12/H13 工厂、`C`、`D_I`、`H_rho`、`F_rho`、`R`、KEP 或黏性。
- **数学约束**：`S_sigma` 是实际 RHS；`sigma=0` 精确零；对自由面 `C(sigma*kappa0*c)-S_sigma=0`；压力 gauge / wall reaction 正确。
- **实现**：copy H12 的 owner/identity 检查结构，不复制其负号定义。支持 XY/RZ、真 AMR 子面、压力分量和 selected 域。
- **验收**：随机强界面对齐、同一 convex drop `p_l>p_g`、静态实际鞍点、误配来源拒绝。
- **失败规则**：出现错误符号、不相容边界行或重建第二套 gradient 即 FAIL。
- **交付**：新的常曲率工厂/测试、完整 free-row residual、旧 H12 回归。

#### BF-03：真实 Hermite 曲率投影与一般 CSF 力

**目标**：将 §4 的 `kappa_q` 正权投影到 native force faces，形成真正的 variable-curvature V1 候选。

- **允许修改**：仅 BF-02 的新 CSF load 类及相关测试；超预算分任务，不改旧 H13。
- **禁止修改**：原 actual source 的插值/Jacobian 规则、FirstJet、质量/KEP/黏性、原 baseline。
- **数学约束**：常 `kappa` 精确再现、局部真实支撑、同一 `C c_eps`、真实 RZ/AMR 物理权重、零支撑拒绝。
- **验收**：真实 Hermite 圆/球/椭圆/正弦曲线、两种 epsilon、uniform/AMR、曲率与载荷独立积分复核。
- **失败规则**：不可拟合频率、默认过滤局部高曲率、扩远程支撑绕过缺失，或偷换解析曲率。
- **交付**：真实曲率误差/积分记录/载荷误差及 PASS/FAIL。

#### BF-04：静滴全范围与旧 H13 A/B

**目标**：在**实际** `CchSinglePhaseSaddle` 和 `CchPhaseViscosity` 上比较 CSF、旧 H12 常跳跃和原 H13 候选。

- **允许修改**：研究 driver/算例输入/基准测试与报告，不修改压力–速度/材料生产模块。
- **禁止修改**：旧 H13 数据、H14/H15、FirstJet、`R/uR`、boundary 契约。
- **数学约束**：实际 `Cp-S_sigma` 静态 residual，真实 `rho/mu`、正交力/壁 reaction，真实 null-space gauge。
- **验收**：完整原 H13 静态矩阵（72 geometries × densities × viscosities），特别保留原七组例外及 RZ/AMR。
- **失败规则**：某些旧失败不改善时不得以解析曲率、人工黏度或更宽 `epsilon` 替换正式门槛。
- **交付**：A/B 表、全部输入/输出、原始失败与新状态，不能提前宣称生产能力。

#### BF-05：无黏毛细波与能量–谱完整 gate

**目标**：只更换冻结旧 H13 研究驱动中的毛细 RHS 为 BF-03；保留原真实 native/source/FirstJet、质量/连续性算法和完全相同的初值。

- **允许修改**：隔离研究驱动、测试参数/结果/新报告；如需修改非毛细生产代码则先暂停并另立设计。
- **禁止修改**：现有 H14、实际 source/Jacobian、材料/动量算子、旧 gate 与已保存失败数据。
- **数学约束**：两相真实 `rho,mu`、`D_Iu=0`、真实 `E_sigma`、`P_sigma`、`R/uR`、完整壁面响应、实际共时几何。
- **验收**：完整一个周期、三至四网格、至少两 `dt`、所有设定 `epsilon/h`、密度比序列；末档频率原目标、完整初值全部谱、能量残差随加密下降、低黏极限无人工 floor。
- **失败规则**：即使频率达二阶，若能量/谱仍 FAIL 则**新方案 FAIL**；禁止仅报告最优 `epsilon` 或使用对照振子替代真实运动。
- **交付**：逐时刻物理账本、全谱/全初值、旧方案 A/B，是否允许进入 BF-06 的明确决定。

#### BF-06：AMR/RZ/接触线与 H14 接入许可（仅 BF-05 全 PASS 后）

**目标**：在真 AMR regrid、轴线、固定接触点及长时间实际界面运动中确认 BF-CSF 没有引入新的不相容性；决定是否冻结为 H13 的替代方法。

- **允许修改**：本模块必要的局部 AMR/边界力带支持和对应测试；任何 H15 职责先拆单任务。
- **禁止修改**：修改 Pure FirstJet 的物理语义、重写 pressure-gradient stencil、隐瞒 `R`、绕过已有 wall reaction。
- **验收**：XY/RZ 动态 AMR、尖端/轴线、pinned contact、断连域、全部报告和身份恢复；通过后出具给 H14 的**明确输入/输出/中点时序契约**。
- **失败规则**：AMR/H15 所需的高阶 metric 或 boundary 超出本任务时，显式标记 `BLOCKED` 并给出依赖项，不把研究驱动误称为 production。
- **交付**：方法冻结版、独立全 gate 报告、下一 H14 所需接口和遗留风险。


---

以下独立任务执行时，每次必须输出：允许修改文件清单、与基线 diff、实际构建/测试、未达标数据、`PASS/FAIL/BLOCKED` 和下一任务是否获准的明确判断。
