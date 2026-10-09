# I12 / phase flux 时间积分设计变更申请

**状态：I12a1 区间契约已完成；用户已授权“开始 I12a2”，并进一步要求“I12a2 后持续推进，直到 I12 通过”。2026-10-01。** 按任务书中心思想和现有程序路线逐项实施、验收；不重复请求已获得的授权。此文不修改冻结 PDF，不提前把 I12 或 Gate I 标为通过。I12a2 先实现常速度、同源平面的单面事件积分，库存、一般界面及完整步进各自拆项。

I12r7 已完成批准前复核，见 [TODO](I12r7_Temporal_Flux_Decision_TODO.md) 与 [审查报告](I12r7_Temporal_Flux_Decision_Report.md)。本次请求的决定是：**允许以同源时间区间平均相通量替换守恒推进中的单点中点相通量，并按下文拆分实施和验收。** 这不批准放宽库存/拓扑门槛、重置解析几何、增加黏度或跳过 Gate I。

## 触发证据

I12 完整链路在平移、旋转、AMR crossing 的半步/全步 bubble 闭合中出现不可达目标。独立平面反例不依赖预测器、曲率、近轴处理或 bubble 区间宽度，直接证明单点 midpoint phase flux 与严格库存及拓扑保持不能普遍兼容。

可执行命令：

```powershell
build/windows/tests/Release/quantum_interface_hermite_gals_benchmark.exe --crossing-probe
```

原始数据：`runs/baselines/I12_20261001/crossing_probe.log`。现版本退出 1，明确暴露待修问题；不设置 `WILL_FAIL` 来使 benchmark 显示为通过。

## 解析反例

单位方形的 4×4 均匀 XY 网格，液相 `x<a(t)`，`u=(1,0)`，`dt=0.02`，`h=0.25`，gross CFL 为 0.08。考虑内部面 `x=0.5`，面长度 `L=0.25`。令 `a(0)=0.5-q*dt`，则穿越发生于 `q*dt`。

真实步内液量为

`DeltaM_exact = L*(1-q)*dt`。

当前 I08 使用同一中点几何的瞬时通量，因此

`DeltaM_midpoint = L*dt`（q<1/2），否则为 0。

| 穿越时刻 q | 原始 midpoint 穿越量 | 解析穿越量 | 上游单元原始 M_new | 上游容量 |
| --- | --- | --- | --- | --- |
| 1/4 | 0.005 | 0.00375 | 0.06125 | 0.0625 |
| 3/4 | 0 | 0.00125 | 0.06375 | 0.0625 |

- q=1/4：I08 接受有界库存，但解析 end skeleton 在上游单元处处为负；库存却缺少 0.00125。若强行 bubble 闭合，必须在全负边界内部制造气泡。I10 正确拒绝 `pure-cell mismatch would create an unmarked component`。
- q=3/4：上游超容 0.00125，即容量的 2%；I08 正确拒绝 `inventory violation exceeds roundoff budget`。不能把 2% 误差算作舍入修复。

对任意更小正 dt，都能选取同样 q 构造反例；减小 CFL 不能给出对所有穿越相位的可靠性保证。这并不证明所有固定初值的自适应时间算法都不可能工作，但排除了以固定 CFL 门槛证明该离散普遍有界的做法。精确事件对齐还会遇到当前 Core 对正长度零迹的拒绝。

该反例中的 midpoint skeleton 是精确平移的平面；bubble 在 face 上恒为零，不能改变这一整面的瞬时相指示。因此它不能通过更精确的 midpoint bubble 预测消除。

**证据更新，保留历史结论的范围：** 初次提案时尚未定位的 vortex `gross phase flux outside total bounds` 已由 I12r1 确认为空间子集积分的舍入缺陷并修复，I12r7 用当前程序重跑 `--gross-flux-contract` 退出 0。它不能再作为未修复时间缺陷的证据。I12r4 独立分解表明，被测 AMR rotation 的下一半步即使用精确旧库存与精确中点几何，仍产生约 `-1.98e-10` 的原始负库存；同样不能据此承诺时间积分能解决全部 bubble/几何问题。I12r6 完整 benchmark 仍为 10 项失败。

I12r7 进一步以精确有理数检查 62 个跨面相位、9 档 dt、正反方向、XY 与 4 个 RZ 轴向环带（含轴），共 5580 个**解析参考组合**，未把解析状态送入生产推进。若 `A` 为面测度、`h` 为流向单元宽度，原始上游库存缺陷满足

`abs(M_midpoint - A*h)/(A*h) = (dt/h)*min(q,1-q)`。

本次 CFL 从 `0.08` 至 `0.0003125`，最小缺陷仍为既有 `64*eps*capacity` 的 `8589934592/25` 倍。RZ 使用沿 z 的不可压平移，`A=pi*(r1^2-r0^2)`；共同 pi 被精确消去。`q=1/2` 对应中点整面零迹，单独保留为未支持事件，未计入矩阵。此反例族随 dt 调整初值，不是固定初值的时间收敛实验。

## 建议变更

将守恒更新中的 `Q_l,f^(n+1/2)` 明确定义为该时间区间的平均相通量：

\[
\bar Q_{l,f}=\frac1{\Delta t}\int_{t_n}^{t_{n+1}}\int_f
H[-\phi_h(x,t)]\,u_h(x,t)\cdot n_f\,dS\,dt.
\]

保留 `M_new=M_old-dt*sum(sign*Qbar)`、唯一实际 face owner、coarse=sum(fine)、独立 M、first-jet skeleton 与标量 bubble 结构。时间段内几何必须由同一旧 authority、同一 velocity/map 及明确预测状态产生；不允许各 quadrature node 由解析真解重新初始化，也不建立第二份同时间 production 界面。

采用事件分段/受检自适应时间求积，显式处理 half-face 相区间出现、消失和跨 endpoint 的事件。共享的是一次计算的 **integrated face liquid volume**；rate 由它除以 dt，仅作对应区间的表达。总通量、正反向 gross flux、RZ 度量和所有时空误差估计也须一致。保留严格容量及守恒 gate，不提高 roundoff budget。

**事件检测是独立前提，不能只靠两个求积值接近。** 例如归一化时间内，一个匀速平移薄液层在固定面上给出 `phi ∝ (tau-1/16)*(tau-1/8)`。中点 `1/2` 和两个半区间中点 `1/4,3/4` 的液相指示全为 0，粗细求积差也为 0，真实积分却为 `1/16`。I12r7 的精确脚本保留这一漏事件反例。未来实现须从同一时空源提供事件完整性的证据（例如受控时空包围/事件隔离）；单靠任意黑盒时间采样，不能宣称已认证未采到的区间为空。无法解析、出现正长度零迹或达到细化工作上限时明确拒绝，不返回伪造的零通量。

在上表反例中，新目标应分别得到 0.00375 与 0.00125，原始上游 M 恰为 0.0625；不需要裁剪或造气泡。一般光滑/AMR 几何的可达性仍须实际证明，不能以该平面算例代替完整 I12。

## 改变与保持的离散性质

改变的是 phase indicator 的时间求积及通量结果的区间语义：从一个中点 snapshot 上的瞬时空间积分，变为受控的时空积分。它超出本轮只接线/验证的 I12 benchmark，且影响冻结方法式 (32)–(33) 与后续中点耦合中 Q 的解释。

保持空间高阶 Hermite、I02 coherent B/J、I03/I04 first jet、I06 master/hanging、I09 共享 reflux、I10 topology guard、I11 单一发布、M 不裁剪、不加全局 upwind/fallback。半步预测也使用同一时间积分规则，不保留隐藏的旧 midpoint 更新分支。

## 拟定实施范围与验收

获准后先将 I12 拆为按顺序执行的子任务，仍不进入 Hydro：

1. **I12a：** 同源时间区间/phase-flux contract 与单面事件积分。原估算为新的 Interface temporal-flux h/cpp，加现有 PhaseFlux h/cpp 的必要共用入口，共 4 production + 最多 2 tests；这只覆盖接口和单面算子，不能包办完整步进。实施前进一步拆成一次一个数学事实的任务：a1 区间源与几何身份；a2 同源事件覆盖及单面四方向积分；a3 共享库存装配与 I09 消费。每个子任务先写 TODO，各自遵守文件预算；若需要 Core 新事件接口，单独列入顺序依赖，不能挤入同一任务。
2. **I12b：** 完整步进的半步、全步消费通过验收的新通量，复跑当前保留的全部失败用例及 h/dt matrix，不改变门槛。若仍有不可达 M，单独研究空间/时间一致性，不能放宽 bubble guard。
3. **I12 后续验收：** 几何相切/零迹、动态 regrid、完整长时强变形 A/B、配置与输出接线等现有缺项分别形成有范围的子任务；全部证据齐全才宣布 Gate I。

I12a 的算子验收至少包括：上述 q=1/4、3/4 及连续穿越相位扫描；XY/RZ；速度变号；AMR 四方向；单面共享/内部抵消；纯相 total/liquid 恒等；时间积分细化；失败不改变 accepted state。库存和拓扑验收使用当前容限。时空 quadrature tolerance 与守恒 gate 分开，不把求积估计当严格误差界。

### 当前 API 对实施顺序的限制

- `HermiteGalsFluxVelocity` 只接受 `Position::Midpoint`，`HermiteGalsPhaseFluxResult` 保留 old/midpoint 并以 `.geometryIdentity()` / `.geometry()` 返回单个中点 authority。新的区间结果必须明确保留 `previousIdentity`、真实 `[begin,end]`、mesh revision、measure 与时空源身份，不能把任意时间样本伪装成同一个 midpoint，也不能把某个 quadrature node 的 ID 冒充整段通量来源。
- 中点几何仍可作为同一推进候选的受力/曲率视图保留，但必须与通量的区间来源分清；多个时间层的只读预测样本不等于第二套同时间 production 几何。每个样本 first jet 必须来自同一旧 accepted state、同一实际 B/J 和 master/hanging closure。面上 bubble 为零，不能为得到面通量反复把 M 当解析几何重建。
- 当前 `advanceHermiteGalsPhaseVolumes` 将空间积分、gross 约束、库存装配和舍入修复集中在一处。必须先明确共用装配的区间单位：四个非负量 `V_total+/-`、`V_liquid+/-` 是时间积分；有符号共享量 `V_liquid=V_liquid+ - V_liquid-` 只存一次；`Qbar=V/dt` 是派生值。总/液相都采用同一时间规则，不能将瞬时 total 与时间平均 liquid 混合。既有空间原始值/修复值诊断应保留可追溯语义。
- `HermiteGalsRefluxResult::source()` 的类型目前固定为瞬时算子结果；a3 必须显式解决结果所有权和消费类型。不得在各 coarse side 重新时间积分；继续汇总 actual subface 的同一个缓存液量，且不二次修正 M。I09 的接口/实现修改计入 a3 自己的预算。
- 半区间与全区间有各自真实边界及源身份；区间切分一致性、相同源的重复调用、过期 mesh、混合 geometry/time、未覆盖事件和失败隔离列为 T1。时间变号时对瞬时正负部分分别积分，不能先平均速度再取正负部分。
- 泛型 velocity callback 和离散 tracer 当前不提供时间包围证书。a1/a2 须明确可验证的时空源表示及事件覆盖能力；若必须额外改变 characteristic 公式、first-jet 状态或加入不同的时间几何重构，则超出本提案，仅有本提案的批准不能默许这些额外方法改变。

## 旧结果的可比性与风险

I01–I11 原始归档与程序哈希保留，原瞬时 flux 算子仅可作为显式 operator A/B，不能成为新路径失败后的 fallback。新旧差异应在相同 mesh/dt/velocity/initial state 下比较，记录原始 Q/M、误差和性能。

新增时间分段将增加求积成本；先做数值正确性，不同时进行缓存/并行优化。一般时空拓扑和符号判定可能要求额外 Core 能力；在未实现前不能保证所有强变形测试通过。该提案只修复目前已证实的时间取样缺口，不替代其余 I12 或耦合验收。

## 授权依据

[开发总纲 §2.1](../00_LAMP_Development_Charter.md) 要求先形成设计变更说明，并明确：**“不得由 Codex 在实现过程中自行改变方法。”** 用户随后明确要求按其 I12a1 任务书执行，并补充接口细节以现有程序路线为主；继而明确授权 I12a2 及持续推进直至 I12 通过。上述时间区间通量路线及其必要的顺序接线据此执行，不再把 I12a1 的历史停止边界当作当前限制。仍不得放宽库存/拓扑门槛、改动冻结 PDF、用解析重置替代推进或跳过验收。I12r7 及更早归档保留当时未批准的历史状态。
