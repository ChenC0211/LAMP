# I12a1：Space–Time Phase Flux 区间契约与来源身份

**任务类型：** Interface V2 / I12 设计修正的第 1 个实施任务  
**状态：** 可实施  
**任务边界：** 本任务只建立“时间区间相通量”的数学/API 契约与来源身份；**不实现事件检测、不实现时空积分、不接入库存推进、不修改 bubble。**  
**完成后必须停止，不得自行进入 I12a2。**

---

## 1. 背景

当前 I12 已经用解析 crossing 反例证明：把 `Q_l,f^(n+1/2)` 实现为**单个 midpoint geometry 上的瞬时空间积分**，在界面于一个时间步内穿过 face 时，不能普遍同时满足：

- 严格守恒库存；
- `0 <= M_K <= capacity_K`；
- end-of-step Hermite skeleton 的拓扑语义；
- cell-local interior bubble 不制造未标记新相分量。

因此，后续 Interface V2 将 phase flux 的 fully-discrete 语义明确为**同源时间区间上的 space–time integrated phase volume**。原有 first-jet GALS、compact `phi_xy` closure、独立守恒 `M_K`、scalar interior bubble、共享 face ownership、AMR reflux 等结构保持不变。

本任务只把这个新语义做成清晰、不可误用的代码契约。

---

## 2. 本任务唯一目标

建立一个 authoritative 的时间区间 phase-flux 数据契约，使后续代码能够明确区分：

1. **时间区间内穿过 face 的积分体积**（authoritative quantity）；
2. 由积分体积除以 `dt` 得到的**区间平均通量率**（derived view）；
3. 该结果究竟来自哪一个 accepted Hermite authority、哪一个 mesh revision、哪一个真实时间区间、哪一个 velocity/characteristic source。

本任务结束时，代码中必须能够表达：

\[
\mathcal V_{\ell,f}^{[t_a,t_b]}
=
\int_{t_a}^{t_b}\int_f
H[-\phi_h(x,t)]\,u_h(x,t)\cdot n_f\,dS\,dt,
\]

以及仅作为派生量

\[
\bar Q_{\ell,f}^{[t_a,t_b]}
=
\frac{\mathcal V_{\ell,f}^{[t_a,t_b]}}{t_b-t_a}.
\]

**程序内部后续守恒推进应消费 `integrated volume`，而不是先除以 `dt` 再乘回来。** 本任务暂不接库存推进。

---

## 3. 不允许改变的数学结构

以下内容在本任务中全部冻结：

- transported state 仍为 first-jet Hermite skeleton；
- `phi_xy` 仍是 compatibility closure，不变成 transported second jet；
- `M_K` 仍是独立守恒相体积；
- interior bubble 仍只做 cell-local scalar volume closure；
- bubble 在 cell face 上仍为零，不允许用 bubble 修 face flux；
- 同一实际 face 后续仍只能有一个 authoritative phase transfer；
- AMR coarse flux 后续仍必须等于 actual fine-subface transfers 之和；
- 不 clipping `M_K`；
- 不扩大 topology guard；
- 不增加人工黏性/upwind/fallback；
- 不把解析真解或重新初始化几何作为 production 时间样本。

---

## 4. Authoritative quantity 与单位

后续 temporal phase-flux 结果的基础量必须使用**时间积分后的体积**，而不是瞬时 rate。

对一个有方向的 native face，至少要能表达四个非负积分量：

```text
integratedTotalPositiveVolume
integratedTotalNegativeVolume
integratedLiquidPositiveVolume
integratedLiquidNegativeVolume
```

语义：

```text
V_total+   = integral max(u.n, 0) dS dt
V_total-   = integral max(-u.n, 0) dS dt
V_liquid+  = integral H_liquid * max(u.n, 0) dS dt
V_liquid-  = integral H_liquid * max(-u.n, 0) dS dt
```

它们的单位是：

```text
XY:  area swept per unit depth / physical 2D volume convention used by current solver
RZ:  true axisymmetric physical volume
```

必须沿用当前 `InterfaceCellMeasure` / face physical measure 的既有约定，不得在本任务里重新定义 RZ 测度。

派生有符号量：

\[
\mathcal V_f = V_{total+}-V_{total-},
\qquad
\mathcal V_{\ell,f}=V_{liquid+}-V_{liquid-}.
\]

派生 rate：

\[
\bar Q_f=\mathcal V_f/\Delta t,
\qquad
\bar Q_{\ell,f}=\mathcal V_{\ell,f}/\Delta t.
\]

**禁止把 `Qbar` 存成第二份独立 authority。** 如果 API 暴露 rate，只能由 integrated volume 和 interval length 现场计算。

---

## 5. 时间区间来源身份契约

新的 interval source/result 必须明确保存或可查询以下语义：

```text
previous accepted Hermite geometry identity
mesh identity / mesh revision
timeBegin
timeEnd
interval length dt
measure mode (XY / RZ，复用现有类型)
velocity source identity
characteristic / tracer source identity
```

如当前仓库已有更合适的 identity 类型，必须复用，不要重新发明平行 ID 系统。

### 5.1 重要约束

- `timeEnd > timeBegin`；否则构造/调用明确失败。
- 整个 interval 必须来自**同一个旧 accepted Hermite authority**。
- 未来在 interval 内生成的预测样本可以有各自 sample identity，但**不得把任意一个 quadrature node 的 geometry identity 冒充整个 interval 的 source identity**。
- midpoint geometry 可以继续作为同一候选步的曲率/受力视图，但不能再代表整个 temporal phase flux 的来源。
- half interval `[t_n,t_{n+1/2}]` 与 full interval `[t_n,t_{n+1}]` 是两个不同 interval object；它们可以共享旧 accepted source identity，但不能共用虚假的 midpoint identity。
- stale mesh / incompatible geometry-time source 必须可检测；如果当前基础类型无法表达这一点，本任务应返回 `BLOCKED`，不得用静态全局计数器或地址哈希伪造 identity。

---

## 6. 区间可加性必须进入 API 设计

后续 space–time integration 会对子区间分段，因此本任务建立的数据结构必须天然支持：

\[
\mathcal V_f^{[a,c]}
+
\mathcal V_f^{[c,b]}
=
\mathcal V_f^{[a,b]},
\]

液相和 gross positive/negative quantities 同理。

建议给 integrated result 提供**显式累加/组合接口**，但只做纯代数合并，不做任何求积或事件判断。例如概念上允许：

```cpp
resultAB += resultAC;
resultAB += resultCB;
```

合并前必须验证：

- face identity / orientation 一致；
- mesh revision 一致；
- old accepted geometry source 一致；
- velocity/characteristic source 一致；
- 两区间首尾连续或符合当前 API 允许的组合语义。

不满足时明确失败，不能静默相加。

---

## 7. 建议文件范围

优先只修改以下文件：

```text
include/Modules/Interface/InterfaceTemporalPhaseFlux.h
src/Modules/Interface/InterfaceTemporalPhaseFlux.cpp
```

测试优先使用：

```text
tests/interface_time_flux_contract.cpp
```

如该测试当前职责过重，可新增一个很小的 contract test，例如：

```text
tests/interface_temporal_flux_interval_contract.cpp
```

并仅在：

```text
tests/CMakeLists.txt
```

注册。

### 文件预算

- production：原则上 **2 个文件**；
- tests：最多 **2 个测试相关文件**（含 CMake 注册）；
- 若必须修改 `InterfacePhaseFlux.h/.cpp` 才能消除旧结果类型与新 interval contract 的类型冲突，**先停止并报告 BLOCKED**，说明需要什么接口变化；不要在 I12a1 顺手扩大到 a3。

---

## 8. 明确禁止修改的区域

本任务禁止修改：

```text
InterfaceQuarticBubble.*
InterfaceQuarticClosure.*
InterfaceConservativeGals.* 的推进算法
AMR/reflux 实现
Hydro
Electricity
Core characteristic 数值公式
first-jet transport
字典/config
时间步控制/CFL
容差/roundoff budget
topology guard
ParaView/IO 输出
```

也禁止：

- 实现事件检测；
- 实现 adaptive temporal quadrature；
- 实现 space–time subdivision；
- 实现 crossing root solve；
- 把 midpoint operator 删除；
- 把旧 midpoint operator 设为新路径失败后的 fallback；
- 修改 I12 benchmark 让当前失败变成“期望失败而通过”。

旧 midpoint operator 目前仅保留为显式 A/B reference。

---

## 9. API 设计要求

不要机械照抄以下名字；如果仓库现有命名体系已有对应类型，优先适配现有风格。必须满足的是**语义**。

建议有两个清晰层次：

### 9.1 Interval source descriptor

概念上类似：

```cpp
struct HermiteGalsPhaseFluxIntervalSource {
    double timeBegin;
    double timeEnd;
    // existing geometry/mesh/source identities
    // existing XY/RZ measure mode
};
```

其职责只回答：

> “这一个 temporal flux interval 从什么 accepted state / mesh / tracer / velocity source 来？”

### 9.2 Integrated face-transfer result

概念上类似：

```cpp
struct HermiteGalsIntegratedFaceTransfer {
    double totalPositiveVolume;
    double totalNegativeVolume;
    double liquidPositiveVolume;
    double liquidNegativeVolume;

    double signedTotalVolume() const;
    double signedLiquidVolume() const;
    double averageTotalFlux(double dt) const;
    double averageLiquidFlux(double dt) const;
};
```

若当前 result 已有 face/geometry identity，则不要复制出一套不一致字段；通过组合/嵌套清晰表达即可。

### 9.3 禁止使用模糊命名

新 authoritative integrated quantity 不要继续只叫：

```text
flux
phaseFlux
midpointFlux
```

至少从类型名或字段名中能看出：

```text
integrated / interval / sweptVolume / transferVolume
```

避免未来再次把 rate 与 integrated transfer 混用。

---

## 10. 本任务 T1 验收

至少新增/更新以下 contract tests。

### T1.1 正常 interval

构造合法 `[t0,t1]`：

- `t1 > t0`；
- `dt()` 返回 `t1-t0`；
- identity 与输入保持一致；
- 不复制/重建 Hermite geometry authority（若现有 API 用共享只读对象，则验证 pointer/source identity）。

### T1.2 非法 interval

以下必须明确失败：

```text
t1 == t0
t1 < t0
non-finite time
```

失败不得改变 accepted interface state。

### T1.3 integrated quantity / rate 关系

人工指定：

```text
V_liquid+ = 3
V_liquid- = 1
V_total+  = 5
V_total-  = 2
dt = 0.5
```

验证：

```text
signedLiquidVolume = 2
signedTotalVolume  = 3
averageLiquidFlux  = 4
averageTotalFlux   = 6
```

rate 不能有独立存储状态。

### T1.4 非负 gross contract

positive/negative integrated components 不允许为显著负值。

本任务不要自行设计新的“大容差”。若只需处理构造输入的浮点 `-0`，使用现有 roundoff 语义；任何超过现有 roundoff contract 的负值必须拒绝。

### T1.5 区间组合

建立 `[a,c]`、`[c,b]` 两个同源结果，验证四个 integrated components 逐项相加，并得到 `[a,b]` 的 interval identity/extent。

### T1.6 禁止错误组合

至少覆盖：

- 不同 mesh revision；
- 不同 old geometry source；
- 不同 velocity/characteristic source；
- 不连续或方向错误的 interval（按最终 API 语义）。

应明确拒绝，不能静默合并。

### T1.7 half/full interval 身份

同一个 accepted source 上：

```text
[t_n, t_n+1/2]
[t_n, t_n+1]
```

应保留同一 old authority，但 interval identity/extent 必须不同。不得都伪装为 `Position::Midpoint`。

---

## 11. 本任务不要求通过的测试

以下当前失败不属于 I12a1 的验收目标，**不要为了让它们通过而改算法**：

- q=1/4、3/4 crossing 数值积分；
- continuous crossing phase scan；
- AMR crossing；
- vortex / deformation；
- bubble reachability；
- I12 full benchmark；
- Gate I。

这些从 I12a2 起逐步处理。

---

## 12. 编译与回归要求

Codex 应使用仓库现有构建方式完成：

1. 修改文件可编译；
2. temporal flux contract test 通过；
3. 与本任务无关的 Interface 基础 contract 不出现新失败；
4. 不改变旧 midpoint phase-flux 数值输出；
5. 不改变任何 production dictionary 行为。

如果平台上的完整构建被已有、与本任务无关的问题阻塞，必须：

- 给出阻塞文件和编译错误；
- 尽可能完成目标 translation unit / test 的最小编译验证；
- 状态报告为 `BLOCKED` 或带明确限制的 `PASS`，不得假称完整 build 已通过。

---

## 13. 完成判据

只有同时满足以下条件，I12a1 才能标记为 `PASS`：

- [ ] integrated face volume 被定义为 authoritative temporal quantity；
- [ ] average flux 只是 `integratedVolume / dt` 的派生视图；
- [ ] interval 明确保存真实 `[begin,end]`，不再用单个 midpoint position 代表整个来源；
- [ ] old accepted geometry / mesh / velocity / characteristic source identity 可追溯；
- [ ] half/full interval 不混淆；
- [ ] integrated result 支持受检的区间组合；
- [ ] 非法/不兼容组合明确失败；
- [ ] 没有实现事件检测或 temporal quadrature；
- [ ] 没有改变库存、bubble、AMR、Hydro、Electricity 或容差；
- [ ] 指定 contract tests 通过；
- [ ] 旧 midpoint operator 数值行为未改变。

---

## 14. BLOCKED 规则

遇到以下任何一种情况，停止并报告 `BLOCKED`，不要自己扩大任务：

1. 当前 geometry/mesh/tracer/velocity API 没有可复用的稳定 identity；
2. 必须修改 Core characteristic 算法才能表达 interval source；
3. 必须修改 `InterfacePhaseFlux` / reflux / inventory consumer 才能完成本任务；
4. 现有类型把“midpoint geometry identity”硬编码为唯一 source，无法在 2 个 production 文件预算内安全拆开；
5. 发现需要决定新的数值方法，而不是单纯建立契约。

报告需要说明：

```text
BLOCKED AT:
WHY:
MINIMAL REQUIRED INTERFACE CHANGE:
FILES THAT WOULD NEED TO CHANGE:
WHY THIS EXCEEDS I12a1:
```

---

## 15. Codex 最终报告格式

完成后只报告本任务，不提出或实施下一步。

```text
TASK: I12a1
STATUS: PASS / FAIL / BLOCKED

MODIFIED FILES:
- ...

CONTRACT IMPLEMENTED:
- authoritative integrated quantities:
- interval/source identities:
- derived average-rate API:
- interval-composition rule:

TESTS RUN:
- command:
  result:
- command:
  result:

REGRESSION:
- old midpoint numerical path changed? YES/NO
- dictionary changed? YES/NO
- inventory/bubble/reflux changed? YES/NO

KNOWN LIMITATIONS:
- event detection: NOT IMPLEMENTED
- space-time quadrature: NOT IMPLEMENTED
- crossing benchmark: NOT EXPECTED TO PASS YET

STOP:
Do not start I12a2.
```

---

## 16. 任务后的下一步（仅说明依赖，不实施）

I12a1 通过后，下一张独立任务才是 **I12a2：常速度 + 平面界面 + 单一 native face 的 crossing-event exact integral**。

I12a2 必须以本任务建立的 interval/integrated-volume contract 为输入，不得绕开它直接在旧 midpoint result 上继续打补丁。
