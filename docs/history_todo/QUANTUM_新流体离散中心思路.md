> **状态：historical/reference（历史资料，仅供追溯）。**
> 本文原始正文保留为历史记录，其中“当前”“下一步”、TODO、旧规则与通过记录均指写作当时；不再自动产生新任务，也不能证明 V2 gate 通过。
> 现行任务唯一入口：[分模块任务清单](../roadmap/05_LAMP_Module_Task_Roadmap.md)。文档权威见[开发总纲](../00_LAMP_Development_Charter.md)，完整入口见[文档索引](../README.md)。
> 旧字典/职责说明可供查询既有实现；V2 新增或修改仍遵守现行六份规范。历史正文中的链接保留原样，目标缺失时不能作为本机已验证证据。
> 本目录中的旧设计、TeX/PDF 和索引同属历史参考；总纲指定的三份冻结方法 PDF 仍以 `docs/` 下的权威文件为准，不因本目录的历史标记而降级。

# QUANTUM 新流体离散中心思路（Codex 实施版）

## 1. 总体原则

最终路线统一为：

**Compatible-Moment Hydro + Conservative Hermite/GALS + consistent mass flux + compatible electrocapillary traction + midpoint coupling**。

核心思想只有一句话：

> **cut geometry, not cut topology**：界面切割材料积分和界面力积分，但不切碎原生 Cartesian/RZ 流体拓扑。

因此，压力单元、速度面、AMR 叶单元拓扑保持固定；Hermite 界面只负责提供真实几何、相体积、相通量、曲率和界面应力积分。

## 2. 界面与两相质量

Hermite/GALS 给出 sharp interface Γ_h。由它计算每个单元的液相体积

    V_l,c = ∫_{Ω_c ∩ Ω_l} dV,
    α_c = V_l,c / V_c.

α 不是 diffuse phase field，而是 cut-volume 的真实体积分数。

界面运动时，不再通过 cell-center 符号重新分类材料，而是计算共享的液相面通量 Q_l,f，并更新

    V_l,c^{n+1} = V_l,c^n - Δt Σ_f s_cf Q_l,f.

总质量通量必须由同一份相通量得到

    F_ρ,f = ρ_l Q_l,f + ρ_g (Q_f - Q_l,f).

**动量输运必须使用同一 F_ρ，不能重新插值 ρ_f。**

## 3. 压力—速度离散

以 Compatible-Moment 空间作为最终 Hydro 核心。离散散度 D 表示原生控制体上的守恒通量平衡；压力载荷与散度必须满足

    C = -D^T

（严格按离散内积定义时可理解为负伴随关系）。

AMR 粗细接口不允许为了“提高压力梯度精度”单独修改 C 破坏该关系。AMR 精度通过局部 Hodge/moment reconstruction 提升。

最终压力—速度耦合直接解鞍点系统

    [ A   C ] [u] = [b]
    [ C^T 0 ] [p]   [0]

projection 可以作为 Schur 预条件器，但不再定义最终物理解。

## 4. 动量对流

默认基准不是 donor-cell，而是与 F_ρ 一致的中心/KEP 型动量通量，例如

    F_mom,g = F_ρ,g * u_adv,g.

光滑区域使用兼容中心通量，目标是无粘极限下不产生人工动能耗散；仅在欠分辨区域加入受控耗散。

原则：

- 质量怎么过面，动量就跟同一份质量过面；
- 不允许 moving variable-density 情况下重新构造一套独立的 ρ_f u_f 通量；
- TraceUpwind 可以作为稳定化/低阶备份，但不应成为默认全局耗散格式。

## 5. 黏性项

黏性采用能量形式组装，而不是逐项手写不兼容 stencil：

    K_μ = E_h^T H_μ E_h.

H_μ 由 Hermite cut geometry 对液相和气相分别积分得到，因此 μ_l、μ_g 的跳变通过真实 cut-volume 积分进入。

必须保持

    K_μ = K_μ^T,    K_μ >= 0.

当 μ = 0 时对应部分应自然消失，不使用最小人工黏度维持稳定。

## 6. 界面力：统一 electro-capillary traction

不要把表面张力和 Maxwell 法向力作为两套互不相关的 body-force deposition。

先在 Hermite 界面上形成统一法向应力

    J_Γ = σ κ - [[ n · T_E · n ]].

理想导体电喷雾中，它就是“毛细压力 - Maxwell 电压力”。

法向界面力必须尽量进入与压力相同的离散 range。对于局部常值 J_Γ，采用

    F_Γ,n = J_Γ C α

或完全等价的 compatible-gradient 写法，使静态平衡满足代数消去：

    C p + F_Γ,n = 0.

若 J_Γ 沿界面变化，则用 Hermite interface quadrature 积分变系数 traction；切向 Maxwell traction 单独做 sharp interface quadrature，不强行写成 pressure-like 项。

## 7. 时间耦合

最终采用单步 midpoint/Picard 耦合，而不是旧的“步初界面力 + CN/AB2”。

同一个 u^{n+1/2} 应同时决定：

- 半步/整步 Hermite 界面；
- Q_l 和 F_ρ；
- cut-material mass/Hodge；
- μ 的 cut integration；
- κ 与 Maxwell traction；
- momentum equation。

即所有两相材料和界面力尽量在 n+1/2 一致评价，避免 force/material lag。

## 8. 必须保留的结构

1. 原生 Cartesian/RZ topology 不因自由界面切割而改变。
2. Hermite/GALS 是唯一几何描述；α/V_l 是由 Hermite 几何导出的守恒库存，不再另建第二套界面。
3. D 是守恒散度；C 必须与 D 保持负伴随兼容关系。
4. phase mass flux 与 momentum mass flux 必须是同一份 F_ρ。
5. 黏性矩阵必须保持对称非负，并允许 μ=0。
6. 法向 electro-capillary traction 与 pressure operator 必须离散兼容。
7. Hydro → Interface 保留 flux-compatible velocity reconstruction。

## 9. 不允许回退的旧做法

最终新路径中不要再出现：

- `phi(cell center) -> 整格 liquid/gas -> rho_c, mu_c`；
- moving two-phase convection 中重新平均/插值得到独立 `rho_f`；
- AMR 接口通过修改压力梯度导致 `C != -D^T`；
- 默认一阶 donor-cell 作为全局动量对流；
- `F_interface^n -> u^{n+1}` 的步初界面力滞后；
- 为了低黏稳定性引入非物理最小黏度；
- capillary 与 Maxwell normal traction 使用彼此不一致的界面几何或 pressure support。

## 10. Codex 实施时的主线

不要继续发展新的并行求解路线。应把现有模块汇合为一条：

    Conservative Hermite/GALS
        -> conservative phase flux/inventory
        -> compatible mass/momentum transport
        -> Compatible-Moment Hydro
        -> compatible electro-capillary traction
        -> midpoint coupled EHD step.

旧 `FittedRZ` 路径只保留 regression/reference，不再作为新算法的承载体。
