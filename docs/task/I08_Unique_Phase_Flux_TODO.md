# I08 / 唯一相面通量与有界约束

## 任务 ID / 标题

I08：unique phase-face flux；唯一验收点为 internal face cancellation + bounded flux constraints。前序 I07 已核验 PASS。

## 目标

同一 midpoint Hermite geometry 与同时间层 native normal velocity 在每个实际 face 上只产生一份液相通量，供两侧 cell 反号使用；有限体积更新保持相库存有界且内部严格成对抵消。

## 允许修改

- 新 `include/Modules/Interface/HermiteGalsPhaseFlux.h`、`src/Modules/Interface/HermiteGalsPhaseFlux.cpp`：不可变 face velocity 输入、唯一通量和独立 M 更新结果。
- `include/Core/Geometry/UnitPolynomialRoots.h`、`src/Core/Geometry/LevelSetHermiteGeometry.cpp`：给既有 positive-polynomial integrator 增加保留 offset/length 的 mapped-range 入口；旧 signed/positive 算法不改。
- 新 `tests/interface_hermite_gals_flux_contract.cpp`、`tests/CMakeLists.txt`。
- TODO/API/报告、docs README、`runs/baselines/I08_20261001/` 证据。

预算：4 production + 2 tests。Core 保持现有 strict-FP/no-LTO 编译环境。

## 禁止修改

I01–I07 公式、mesh、Hydro/Electricity、Apps、reader/case、冻结 PDF/S08 专有实现。不得调用旧 upwind/FCT fallback；不实现 I09 coarse-face 聚合/reflux、I10 bubble 或 I11 publisher。不增加另一份零集几何或独立 cell-face 通量。

## 数学约束

1. `Q_l,f = integral_f H(-phi_mid) u_n,mid dA`；XY 真长度，RZ 真 `2*pi*r ds`。负相分段来自同一 native Q4 snapshot；沿实际两个 half-face chart 保留 root offset/width，不以端点相减丢弃小区间。
2. 每个 actual native face 一个 owner/result entry，朝 +x/+y。`M_new = M_old - dt sum(s_Kf Q_l,f)`。每个内部条目严格以反号进入两侧；边界净出量单独记录。
3. normal velocity 是给定的局部 polynomial trace（每半面最高 degree 10，RZ metric 后最高 11），不是新增 Hydro unknown 或全局 C1 速度状态。输入深拷贝，绑定 midpoint geometry/time/mesh，不依赖任意积分回调的一致性承诺。
4. 分别保留 `Q+ / Q- / Q_l+ / Q_l-`，要求非负且 `Q_l+<=Q+、Q_l-<=Q-`。空间 counterflow 的 gross outflow 必须参与 CFL，即使 signed net=0；不以有符号分数误判有界性。
5. old State 必须对应 midpoint interval 起点与前一步；midpoint geometry 和 velocity 的 identity/time 一致。old M 是更新输入，midpoint M 不覆盖库存；输出记录两种来源身份。
6. 不可渗透 face 与 RZ axis 的 normal polynomial 必须恒零，不能只检查净积分为零。外边界流入相位由显式给定的 midpoint scalar trace 提供。
7. 要求总通量的离散 divergence 在独立容限内，gross-outflow CFL<=1，最终每 cell `0<=M<=capacity`。有限幅值越界明确拒绝；仅机器舍入尺度可经内部共享 face 的微小通量调整修复。修复不得改 prescribed boundary、total flux 或直接 clipping M；失败/预算耗尽仍拒绝，不切到 upwind/FCT。

## 实现要求

Core mapped positive primitive 用 exact expansion 进行同一 polynomial 的平移/缩放，再调用已有 bounded positive integration；不重建 rounded range endpoints。纯相与零速度路径保持精确零/完整通量。Result 只读，含唯一 faces、updated M、identity/time 和具名 Report；原 old/midpoint States 不变，失败不发布 partial result。

微修复预算固定为 `64*eps(double)*cell capacity`，只用于舍入级 defect，不作为用户调参。具名 flux integration / divergence / conservation 容限分别校验并报告；求积估计不冒充严格 PDE 精度。Report 不参与数值分支。生产内核不读字典、不写文件。

## 验收

- **T0：** Release 主程序、新 contract、I01–I07、architecture；无新增编译 warning，Core strict-FP/no-LTO 保持。
- **T1 Core：** positive mapped integrals 与解析常数/线性/高次 polynomial 比较，near-endpoint 和 sub-ULP phase width 不丢失；旧 signed-range 与旧 positive 调用回归。尺度归一化误差 `<=3e-12`，小区间另检查相对误差。
- **T1 唯一通量：** XY/RZ、uniform/actual AMR faces、静止/全液/全气/解析平面、多根 phase trace、正反速度与同面 counterflow；逐 face 对独立解析 primitive，归一化误差 `<=3e-11`。内部成对贡献严格相反，global ledger 相对 defect `<=1e-12`。
- **T1 有界：** 正反 gross phase bounds、gross CFL、zero-net counterflow、不可渗透 polynomial 恒零；合法 raw flux 不改。构造机器舍入级越界验证仅内部共享通量修复、boundary 不改、M 严格有界、全局守恒；宏观越界/CFL/非散度自由明确拒绝，无 fallback。
- **T1 ownership：** 外部输入修改、旧 M 与 midpoint M 不同、source 生命周期、错误 ID/time/mesh/缺 face/重复 face/非有限项等；拒绝不改源。
- **T2：** 公开 API 完成平面界面 XY/RZ 的单步解析 midpoint 输运，液量及逐 cell 更新 `<=3e-11`，零速度库存严格不变。单独光滑 face cut 的 N=8/16/32 空间误差下降，最后阶 `>=3`；不把此局部验证称作 I12 完整输运。
- **回归：** I01–I07、旧 native geometry/GALS/parametric control 数值不变，69 个 S08 文件与 3 PDF 哈希不变。Core extraction 前后数学实现保持；所有异常/早期失败日志保留。

## 失败规则

数值/有界/ownership 失败记 FAIL 并在范围内修正；缺必要数学接口、需改冻结方法或超预算记 BLOCKED 并说明证据。不放宽预设验收容限，不把通量 gate 的拒绝冒充稳定长时间运行。不开始 I09。

## 交付

源码/测试、API/TODO/报告、原始数值与构建/失败日志、源码/runtime 哈希、逐项校验的证据归档。只有全部验收通过才记 I08 PASS；总体 goal 与 Gate I 保持未完成。

## 执行结论

**PASS，2026-10-01。** 见 [报告](I08_Unique_Phase_Flux_Report.md) 与 [API](../api/I08_Hermite_Gals_Phase_Flux_API.md)。唯一实际面、gross bounds/CFL、内部反号抵消、微修复、8 组解析输运、XY/RZ face-cut 收敛、18 项拒绝、前序数值及冻结回归通过。12 个 CTest 目标最终均有通过证据；初次额外相切支持检查失败及重测保留。

实际半面端点相切通过；Core 对半面内部精确相切/重根仍拒绝，这项已知几何限制没有被宣称解决，须在 I11/I12 完整运动几何 gate 前继续处理。2 条冻结 Hydro 警告与 S08 一致，无新增 warning。下一项 I09，总体 goal 继续 active。
