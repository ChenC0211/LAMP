# I12r4 / 旋转曲率与库存误差来源审计

## 任务 ID / 标题

I12r4：定位 I12r3 后 AMR rotation 首步曲率误差 32.942469 与下一步库存越界的来源。前轮 I12r3 为 progress，父任务 I12 仍 FAIL。

## 目标

只审计一个事实：在真实 rotation first-step 的相同预测 skeleton 上，体积目标误差如何通过唯一 interior bubble 放大到曲率，并与下一步库存失效关联。给出可执行、独立解析对照，不凭大曲率日志推断算法错误或批准方法变更。

## 允许修改

- `tests/interface_hermite_gals_benchmark.cpp`：仅新增显式 `--rotation-curvature-audit` 诊断，默认路径保持；0 production + 1 test。
- 本 TODO/报告、`docs/README.md`、`runs/baselines/I12r4_20261001/`。

## 禁止修改

所有 production、冻结方法、I08 时间公式、I10 接受门槛及 I11 publication。禁止将解析 M/geometry 接入真实推进、放宽体积/位移容限、修改曲率或隐藏失败。时间积分提案仍未批准。本轮不实现修复、不进入 H01。

## 数学约束

1. 实际半步/全步使用原 first-jet、compact mixed closure、原唯一 face flux 与原 bubble；解析支路仅为不返回给 advance 的测试对照。
2. 在完全相同的预测 skeleton 上，仅将复制候选的 M 换为独立 circle–rectangle 解析相交量，检验体积目标的因果作用；拒绝时如实记录。
3. actual 与 predictor 的差应等于 `delta_c*q4(xi)*q4(eta)`。独立微分该 bubble，逐点核对 value/gradient/Hessian 差；独立 jet 曲率公式与原 frame 曲率一致。
4. 解析 rigid-rotation 圆与面相交区间、线性法向速度积分给出独立瞬时 midpoint 通量；分解 `M_actual-M_exact = (M_midpoint_exact_geometry-M_exact) + (M_actual-M_midpoint_exact_geometry)`，不把空间几何误差和时间取样混同。
5. 后续失败用实际 accepted state/quarter prediction 重现，记录原始未裁剪库存和每个相邻 native face 的贡献；独立面内积分先在成功步与 I08 结果交叉验证。
6. 所有诊断及反事实分支不改变实际 accepted state、candidate 或原默认 benchmark 输出。

## 实现要求

XY rotation 原 n=8、dt=0.01 的 uniform/AMR 两组；记录 initial、quarter、midpoint/end predictor 与闭合态的 position/gradient/curvature、最大误差位置/cell、c 修正/体积灵敏度/原始 M。额外 dt=0.005 的 AMR 首步对照仅诊断敏感性，不冒充完整时间收敛。保留原默认全部算例与容限。

## 验收

- T0：benchmark Release 构建/architecture，无新增 warning，源码范围和 frozen/PDF hash；上轮归档与交付文件身份核验。
- T1：独立解析体积总量、面通量积分、bubble jet 分解、曲率公式与 snapshot 对照；对照分支隔离。
- T2：原两组旋转首步及 AMR 后续库存失败有完整诊断；默认完整 I12 输出与 I12r3 逐字一致、退出码仍为 1，现存失败不能变绿。
- 若诊断支持生产缺陷或方法缺口，在报告提出下一单一任务；诊断不足则明确局限，不推断修复成立。

## 失败规则

若曲线/求积/反事实闭合无法解析，保留错误及缺失证据；不改 production 或容限绕过。若需要新数学方法，另写设计变更且不在本任务实施。

## 交付

显式审计入口、原始日志/数据、误差来源结论及限制、源码/runtime/归档 hash、报告。父任务和全项目状态按完整证据保留。

## 执行结果

**PASS（诊断审计），父任务 I12 仍 FAIL。** 0 production + 1 test；独立 bubble 导数/曲率及 6 次面通量核验通过，定位了小 cut cell 的大系数修正与中点时间取样误差。5 个精确 M 对照成功，1 个因证书区间不足拒绝；下一步原始负库存约为舍入预算的 53 万倍。默认完整 benchmark 输出逐字不变、退出码仍为 1，所有失败保留。详见 [报告](I12r4_Rotation_Curvature_Audit_Report.md)。本轮不启动下一任务。
