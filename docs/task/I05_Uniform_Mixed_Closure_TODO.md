# I05 / 均匀网格 phi_xy 紧凑闭合

日期：2026-10-01。依据：[05 路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md) I05、[01 数值硬约束](../01_LAMP_Numerical_Invariants.md) I-R1/I-R2，冻结 Hermite–GALS 方法 §5.2 式 (27)–(28)。前一轮为 progress；开工复核 I04 交付和归档。

## 目标

只解决一个数学事实：从同一 primary revision 的已输运 vertex gradients 和 edge normal derivatives 恢复 `phi_xy=0.5*(Dx(phi_y)+Dy(phi_x))`，作为派生 Hermite closure，不推进 Hessian 或改变 first jet。

## 允许修改

- 新建 `include/Modules/Interface/HermiteGalsMixedClosure.h`、`src/Modules/Interface/HermiteGalsMixedClosure.cpp`。
- 新建 `include/Core/Geometry/CompactQuarticDerivative.h`，把 `src/Core/Geometry/ActiveQuarticHermiteGals.cpp` 已有五点四次插值求导函数原样提取；旧调用保留数学/取点/求和顺序和错误语义。不调用旧 predictor 推进 V2。
- 新建 `tests/interface_hermite_gals_mixed_contract.cpp`，修改 `tests/CMakeLists.txt`；该目标同时编译已有、未修改的 `interface_gals_contract.cpp`，用于提取前后回归。
- 本 TODO、I05 报告、`docs/I05_Hermite_Gals_Mixed_Closure_API.md`、`docs/README.md`，证据目录 `runs/baselines/I05_20261001/`。

预算：4 个 production 文件，2 个测试文件。共用 Core 求导段不属于 S08 冻结专有实现；提取前后必须核对旧 predictor，不能扩展旧方法。

## 禁止修改

不改 I01–I04、Core mesh/snapshot/tracer、S08 冻结专有实现、Hydro/Electricity、Apps/reader/case/默认配置。不开 I06 AMR closure/transfer、I07–I10 守恒体积/bubble 或 I11 新几何发布；不引入第二份 scalar、独立 phi_xy 输运或 legacy fallback。

## 数学约束

- 内点沿一条 skeleton 线用 `{-h,-h/2,0,h/2,h}` 的四次插值导数，等价于 `[f(-h)-8f(-h/2)+8f(h/2)-f(h)]/(6h)`。
- x 方向读取 vertex.dy 与水平 edge.dn；y 方向读取 vertex.dx 与竖直 edge.dn。两个方向固定等权平均，不使用 disagreement 调节数值公式。
- 物理边界使用相同四次插值的一侧五点 `{0,h/2,h,3h/2,2h}`（上边界镜像）；保持 Q4 再现，不读取域外值。
- 采样点使用实际 native 坐标及端点 midpoint；在局部平移/缩放坐标求权重，减去 target 值以使常数导数严格为零。
- I05 只接受至少每方向两单元的完整同层 uniform leaf skeleton；缺失、非 native、AMR/hanging 或不一致输入明确拒绝。支持各向异性网格与全局同层 refinement。
- 输入为尚未发布 geometry 的 I01 candidate，输出只读 mixed vector，顺序与 primary.vertices 相同，带 primary identity/time。返回后 candidate 改写会使旧结果 identity 失效。
- RZ min-r 非负；若包含 r=0，axis vertex.dx 和竖直 edge.dn 必须严格零。even scalar 的轴线 closure 由 parity 得到严格零，不“修复”不合格 first jet。

## 实现要求

- [x] 复核六份规范/冻结公式，保存开工 TODO 和 source baseline。
- [x] 在改 Core 前运行现有 GALS contract，保留原始日志/指标；纯求导段只保留一个实现。
- [x] 实现 uniform mixed closure，校验来源/拓扑/覆盖/axis，结果无可写别名；错误包含 primary identity、阶段及 vertex。
- [x] 输出 directional disagreement 诊断，不能参与接受、混合权重或 fallback。
- [x] 验证 full Q4 reproduction、四阶收敛、常量导数严格零、两方向平均、失败隔离及 I03/I04 输运结果接入。
- [x] 完成前序/旧 predictor 回归、文件哈希核验和交付文档。

## 验收

- **T0：** Release 构建主程序、I01–I05 contracts、旧 Interface 与 architecture；新 header 首个 include；无新增 warning；architecture/help 通过。
- **T1：** 单点求导对次数 0..4 的多项式，在内部/左右端点且平移/非二进制尺度坐标上误差 `<=2e-12*(1+|reference|)`；常数严格零，非法/重复支持点拒绝。
- **T1：** uniform native skeleton 对全部 25 个 `x^a*y^b, 0<=a,b<=4` 及其非平凡线性组合，在含边界的各向异性/平移 mesh 上再现 dxy，误差 `<=3e-11*(1+|reference|)`；全局同层 refinement 同样通过。
- **T1：** 人为不同的两个 directional estimate 按固定平均得到正确标量，并保留 disagreement 诊断；phi 与 edge.phi 不参与 dxy 重构。
- **T1：** RZ 偶 Q4 及 axis parity，缺失 vertex/edge、假 edge、AMR、过小 mesh、已发布候选、mesh 变化、非有限结果拒绝；source primary/identity/geometry 状态不被修改。
- **T1 接线：** 真实 Q4 snapshot 经 I03/I04 的全 skeleton 非线性输运后，组装未发布 candidate 并闭合；与独立冻结中心/一侧公式相比 `<=3e-11*(1+|reference|)`；来源 identity/time 由新 primary revision 绑定。
- **T2 专项：** XY `exp(0.7x+0.4y)` 与 RZ `cos(0.7r)*exp(0.4z)`，Nx=Ny=8,16,32,64，内部/边界/全域 dxy 最大误差分别记录；最后两级全域收敛阶 `[3.7,4.3]`，axis error 严格零。只证明 closure 空间阶，不替代完整 Interface/曲率收敛。
- **回归：** I01–I04、旧 `--parametric-state-control` 与提取前后的现有 GALS contract 通过；旧 GALS 数值输出逐字一致。S08 69 文件、三份方法 PDF 哈希不变。

## 失败规则

公式、Q4 再现、固定平均或来源隔离失败即 FAIL，在范围内修正；不放宽容限以适配结果。需要超预算、改变冻结方法或前序数值契约时说明最小缺口并拆任务；无法在范围内完成则 BLOCKED，不用回退、额外全局未知量或二阶速度导数替代。

## 交付

I05 kernel/只读结果、共用 Core 求导、独立 contract、TODO/报告/API、误差和阶数原始证据、构建/回归及哈希归档。仅完成 I05；结束状态 `PASS / FAIL / BLOCKED`，下一项 I06。

## 执行结论

**PASS，2026-10-01。** 见 [交付报告](I05_Uniform_Mixed_Closure_Report.md) 和 [API 说明](../api/I05_Hermite_Gals_Mixed_Closure_API.md)。Release 构建、CTest 7/7、旧 Interface 与旧 GALS 提取前后回归通过；Q4 再现和四阶空间收敛满足全部预设阈值。

首次编译修复了新测试中的不可赋值 DictionaryMesh 和多余括号；运行测试前也修正了临时结果的引用生命周期。首次 CTest 的负半径用例由 I01 构造器提前拒绝，最终检查这个真实上游拒绝；未改变生产公式或验收容限。旧 GALS 的 55 行输出、I03/I04 的 23 行数值输出保持一致。下一项 I06，Gate I 尚未完成。
