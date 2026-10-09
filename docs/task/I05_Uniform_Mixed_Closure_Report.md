# I05 / 均匀网格 phi_xy 紧凑闭合交付报告

**结论：PASS。** 执行日期：2026-10-01。本轮只完成 I05，下一项为 I06。Gate I 与完整低黏/无黏电喷雾目标仍未完成。

依据和预设阈值见 [TODO](I05_Uniform_Mixed_Closure_TODO.md)，公式、单侧边界与 axis 语义见 [API 说明](../api/I05_Hermite_Gals_Mixed_Closure_API.md)。开工核验 I04 的 10 个交付文件及归档，保存 1,552 文件的 source baseline 和开工 TODO；前一轮归类为 progress。

## 交付与数学关系

| 文件 | 变化 |
| --- | --- |
| [HermiteGalsMixedClosure.h](../../include/Modules/Interface/HermiteGalsMixedClosure.h) | 未发布 candidate 输入、只读 mixed result、primary identity/time、诊断 |
| [HermiteGalsMixedClosure.cpp](../../src/Modules/Interface/HermiteGalsMixedClosure.cpp) | complete uniform skeleton 校验、两个方向固定平均、边界与 RZ parity、整批隔离 |
| [CompactQuarticDerivative.h](../../include/Core/Geometry/CompactQuarticDerivative.h) | 纯五点四次插值求导的共用实现 |
| [ActiveQuarticHermiteGals.cpp](../../src/Core/Geometry/ActiveQuarticHermiteGals.cpp) | 只把既有求导段改为共用调用，旧公式/取点/求和/错误上下文不变 |
| [mixed contract](../../tests/interface_hermite_gals_mixed_contract.cpp) / [注册](../../tests/CMakeLists.txt) | 独立解析/收敛/来源检验，并编译未改的已有 GALS contract 做前后对照 |

预算：4 个 production、2 个测试文件。I01–I04、Core mesh/snapshot/tracer、69 个冻结专有文件、Hydro/Electricity、Apps、reader、case 和默认配置均未改。

新闭合从 vertex.dy/水平 edge.dn 得到 Dx(phi_y)，从 vertex.dx/竖直 edge.dn 得到 Dy(phi_x)，固定取平均。所有输入属于同一 candidate primary revision；不读取旧 closure，不改 phi/gradient/edge/M，不推进 Hessian，也不发布 snapshot。Core 共享段不是旧算法 fallback，V2 没有调用旧 predictor。

## 多项式和接线验证

证据目录：[runs/baselines/I05_20261001](../../runs/baselines/I05_20261001/)。原始测试见 `ctest_details.log`，表格数据来自相应 CSV。

| 检查 | 覆盖 | 最大绝对误差 |
| --- | --- | --- |
| 一维 degree 0..4 | 25 个内部/左右边界组合，平移非二进制点 | 1.59872e-14 |
| 全 Q4 与线性组合 | 5×3 各向异性平移域，26 个场、全部顶点 | 2.48690e-14 |
| 全 Q4 与线性组合 | 全局细化后的 6×4 actual leaves，26 个场 | 3.21965e-15 |
| RZ even Q4 | 25 顶点，含 axis | 3.55271e-15，axis 严格零 |
| I03/I04 输运后闭合 | 81 vertex + 144 edge，与独立中心/单侧显式公式比较 | 2.44249e-15 |

常数 gradients 闭合严格零。人工设置 Dx(phi_y)=2、Dy(phi_x)=4 得到 mixed=3，relative disagreement=0.5；改变 phi/edge.phi 不改变结果，表明闭合使用已输运的一阶数据。candidate replacePrimary 后旧 mixed result 仍保留旧 identity；无自动重新绑定。

输运接线用真实正 Q4 snapshot，经 I03/I04 的非线性速度回溯，将两个阶段结果装入新的 End candidate，再由 I05 闭合。全节点核对独立冻结公式，旧 source 和新 primary 均不修改，也没有提前发布新 geometry。

Q4/收敛测试直接构造未发布 first-jet samples，不声称它们已有几何相体积闭合；接线测试中的正 scalar 具有真实零负相库存。此次只验证 closure 算子，不以构造样本替代完整 Interface benchmark。

## 四阶空间收敛

XY 使用 `exp(0.7x+0.4y)`，RZ 使用 `cos(0.7r)*exp(0.4z)`；Nx=Ny=8,16,32,64。内部、物理边界和全域最大 dxy 误差分别记录于 `spatial_convergence.csv`。最细一级为：

| 场 | 内部最大误差 | 边界/全域最大误差 | 最后一级全域阶 |
| --- | --- | --- | --- |
| XY | 1.36453e-11 | 8.24512e-11 | 3.988949 |
| RZ | 4.35818e-12 | 2.63802e-11 | 3.983710 |

最后两级全部在预定 `[3.7,4.3]` 内，所有 RZ axis error 精确为零。数值容限未调整：单点导数为 `2e-12*(1+|reference|)`，Q4/接线为 `3e-11*(1+|reference|)`。该结果证明所测 mixed closure 的四阶空间精度，不证明完整输运的时间阶、曲率收敛或低黏稳定性。

## 拒绝路径与回归

缺失 vertex/edge、假长边、AMR/hanging、每方向不足两单元、已发布候选、mesh 身份变化、非有限求导结果、RZ axis 不合格均明确拒绝，primary 数据、identity 与 geometry 状态保持原样。负半径候选实际由 I01 构造器提前拒绝，测试核对该上游行为；不声称直接覆盖 I05 的重复防御分支。

| 检查 | 结果 |
| --- | --- |
| Release 主程序、I01–I05、旧 Interface、architecture | PASS；全目标构建 9.974 s，无新增 compiler warning |
| CTest：help / architecture / I01 / I02 / I03 / I04 / I05 | 7/7 PASS；0.02 / 0.72 / 0.01 / 0.01 / 0.05 / 0.08 / 0.37 s |
| 旧 `--parametric-state-control` | PASS，退出 0，0.380 s |
| Core GALS 提取前后 | 原有 contract 通过；55 行完整输出逐字一致；计算体仅去空白后完全一致 |
| I03/I04 数值输出 | 23 行与 I04 归档逐字一致 |
| 既有源码/输入/工具 | 1,550 文件哈希不变；仅旧 Core 求导所在 cpp、tests/CMakeLists.txt 改动，另增 4 文件 |
| 冻结专有实现 / 方法 PDF | 69/69、3/3 哈希一致 |

提取前的程序与日志保存在 `legacy_gals_before.exe/log`，提取后独立模式 `--legacy-extraction-control` 的日志与对照见 `legacy_gals_after.log`、`legacy_extraction_comparison.json`。I03/I04 对照见 `I03_I04_numeric_comparison.json`。

首次编译的测试文件错误是 DictionaryMesh 不支持赋值以及多余括号；已修正。运行前检查也修正了临时结果引用的生命周期。首次 CTest 为 6/7，负半径测试在 I01 构造阶段提前拒绝；最终检查正确的拒绝位置。所有初始失败日志保留，生产公式与验收阈值未因此改变。

`verification.json` 记录命令、线程与摘要，`runtime/` 及 `runtime_hashes.json` 记录实际执行文件；交付副本、逐文件哈希和归档见 `delivery/`、`delivery_files.json`、`evidence_inventory.csv`、`evidence.zip`、`delivery_manifest.json`。Git 元数据不可用，没有伪造 commit。

## 未完成

I06 AMR master/hanging closure、I07–I10 守恒体积/通量/bubble、I11 同源几何发布和 I12 完整 benchmark 均未在本轮开始。当前入口只支持 complete uniform skeleton；局部 active-band 缺失、contact/sector、物理边界输运历史及 regrid 仍需后续实现。

V2 未接入 reader/生产 App 或新 Hydro 主线。S01 已知 `hydro_moment_xy` 首步失败继续保留，本轮短回归不能重新分类它。下一项按路线图执行 I06。
