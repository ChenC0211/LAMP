# I06b / 非均匀 mixed closure 与 master 输运接线

## 任务 ID / 标题

I06b：完成 I06a 后的 AMR master/hanging skeleton 闭合。父任务仍是路线图 I06，不开始 I07。

## 目标

同一未发布 primary revision 的 master first jet 决定 compact mixed derivatives，并由同一 Core exact master constraints 生成 hanging/fine 数据；coarse/fine C1 trace 一致，hanging 不独立推进。

## 允许修改

- 新 `include/Modules/Interface/HermiteGalsAmrClosure.h`、`src/Modules/Interface/HermiteGalsAmrClosure.cpp`。
- Core `ActiveQuarticHermiteSkeleton.h/.cpp`、`Detail/LevelSetHermiteQuarticLeafData.h`：开放 node 对 master 输入的既有约束链式导数，不另写 restriction/stencil。
- 新 `tests/interface_hermite_gals_amr_contract.cpp`、`tests/CMakeLists.txt`。
- 本任务 TODO、报告、API、docs README 与 `runs/baselines/I06b_20261001/` 证据。

预算：5 production + 2 test。复用 I05 `CompactQuarticDerivative` 和 I06a 的严格浮点 TU，不改编译配置。

## 禁止修改

I01–I05 的 primary/transport/uniform 公式、mesh、Hydro/Electricity、Apps、reader/case、冻结 PDF、S08 专有实现；不做 phase volume/flux/bubble 更新、geometry publication 或输出帧。不增加独立 mixed/hanging 动力学状态。

## 数学约束

1. 每条 connected master-edge 直线上按真实物理坐标排序端点与 master midpoint，不跨越 inactive/gap。每个 regular target 用连续五点的 quartic interpolant 求导；内部尽量居中，端部单侧；不足五点明确失败。
2. x 线使用 phi_y，y 线使用 phi_x；master midpoint 读取 primary.dn，hanging 端点由 Core trace 限制取得相应 tangential component。固定 `q_v=(Dx(phi_y)+Dy(phi_x))/2`；RZ axis q=0。
3. hanging first jet 可能线性依赖待闭合的 master q，不能把未知 mixed 任意置零后当成结果。用 Core 的同一依赖图写成 `g_h(q)=a_h+C_h q`，代入上述固定闭合，求解 `(I-LC)q=La`。这是对 frozen compact closure 与 hanging constraints 的代数消元，不改变方程，不是额外时间推进状态或 KKT 几何修正。
4. 在 q=0 处只装配 affine 常量项；该临时算子不向外发布。求解后重新绑定最终 exact skeleton，并直接核查原闭合公式残差。禁止省略非零依赖、用任意 seed 代替解或通过低阶 stencil 避免约束。
5. primary、M、source geometry、mesh 不变；输出只读 master mixed + exact skeleton，绑定 candidate identity/time。不将 rounded hanging data 重新注册为独立量。

## 实现要求

从 primary.geometryCells 构造支持集，验证完整 master vertex/edge 覆盖和真实 owner。支持完整 2:1 网格及 coarse/五点支持完整的 active subset；缺支持明确失败。沿用已有 quartic constraint transpose 传播，仅新增 node 查询入口。局部五点模板、稀疏 transient compatibility system；不持久化矩阵/新增全局状态。求解与直接代回残差分别记录并检查 `<=1e-12` 的缩放误差；diagnostic disagreement 不改变平均权重。

## 验收

- **T0：** Release 主程序、新 contract、I01–I06a、architecture；无新增 warning，头自包含，结果只读。
- **T1 数学：** uniform 与 I05 一致（`<=3e-11*(1+|reference|)`）；XY 全 Q4、RZ even Q4 在 isolated/adjacent/multilevel 2:1、shifted/anisotropic 域再现；常数梯度 mixed 严格零。Core node sensitivity 与独立扰动/线性叠加吻合；冻结现有约束结果不变。
- **T1 AMR：** 任意 master first jet 经闭合后，每条 coarse/fine 重合段至少 9 点，phi/gradient/dxy trace 差 `<=3e-11*(1+|reference|)`；覆盖 horizontal/vertical 和端点为 hanging 的 perpendicular master edges。矩阵规模、实际支持、hanging/fine 数目和残差记录。
- **T1 失败隔离：** 缺 master、fine 伪装 master、缺 coarse/五点支持、geometry 已绑定、mesh 身份变化、axis 非零、非有限计算明确失败；candidate 内容/identity 不变。
- **T2 接线：** 真实正 Q4 source snapshot → I03/I04 master-only 非线性输运 → 新 candidate → AMR closure；独立五点权重代回、C1 和来源身份检查通过。不查询或推进独立 hanging target；旧 source 和新 primary 不改。
- **T2 收敛：** XY `exp(0.7x+0.4y)`、RZ `cos(0.7r)*exp(0.4z)`，base N=8,16,32,64，加局部 level-1 refinement。记录 master/hanging mixed 最大误差、最后两级 master 阶需在 `[2.7,4.3]`；uniform 保持 I05 四阶。AMR tangential trace 求导的 O(h^4) 误差再沿 skeleton 求导可能产生 O(h^3) 项，故不能预称所有 AMR 点四阶。
- **回归：** I01–I06a、architecture、旧 native geometry/AMR 与 parametric control；I03/I04 数值、旧 GALS 输出不变；69 frozen files、3 PDFs 哈希不变。

## 失败规则

不改变 frozen 方法、不放宽已定容限/阶数。公式/残差/覆盖/C1/Q4 失败则 FAIL 并在范围内修正；缺数学或架构支持、预算无法满足则 BLOCKED，记录具体证据。不得以 C1 单项通过替代 mixed 或接线验收。

## 交付

5 production + 2 tests、API/TODO/报告、原始指标/构建/失败与回归日志、运行文件与源码哈希及归档。只有所有条目通过才记 I06b 与 I06 PASS；Gate I、I07–I12 及完整目标继续保持未完成。

## 执行结论

**PASS，2026-10-01；结合 I06a 完成父任务 I06。** 见 [报告](I06b_AMR_Mixed_Closure_Report.md) 和 [API](../api/I06b_Hermite_Gals_AMR_Closure_API.md)。Q4、MMS、C1、独立公式代回、真实 master-only 输运接线、supported/disconnected subset、拒绝隔离、10/10 CTest 与冻结哈希通过。观察到的 compatibility matrix 为单位阵，未声称覆盖非对角系统；试验性非法网格及消息检查修正有原始失败日志。数值公式、误差和阶数门槛未放宽。下一项 I07，Gate I 未完成。
