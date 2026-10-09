# H08a2r4 / 同一 advector 的正则几何准备

Owner 为 `lamp_hydrodynamics_cch`。现有 `CchAdvectorReconstruction::prepare` 实际持有 [Detail header](../../include/Modules/Hydrodynamics/Detail/CchRegularAdvectorGeometry.h) 与 [实现](../../src/Modules/Hydrodynamics/CchRegularAdvectorGeometry.cpp)。公开 `regularGeometry()` 只读访问准备数据；复制保留身份。[r11 第二次恢复](../task/H08a2r11_Common_Regular_Source_Deployment_Resume2_Report.md) 已将实际 `evaluate`、B/J、atom 和 characteristic 同步切换为这份准备的唯一系数。

```cpp
const auto owner = CchAdvectorReconstruction::prepare(mesh, momentumGeometry);
const auto& prepared = owner.regularGeometry();
const auto data = prepared.apply(mesh, nativeValues, &explicitBoundary);
```

`apply` 生成代数数据，不发布速度 owner、clock 或 Interface callback。`cellDofs[p][component][5*y+x]` 为 Q4 tensor 数据；一维顺序为端点 0/1 值、端点 0/1 无量纲导数、中点值。法向和混合导数分别乘实际 cell 宽度、宽度乘积。`faceCorrections[u]` 为 canonical normal 方向的 `15/8*(nativeMean-baseMean)`；轴零测度 face 的幅度严格为零。H08a2r8 为这些系数附带私有原几何 shared identity；它只标识准备层，仍不具有可发布 velocity identity 或 clock。正式 `CchSpatialAdvector` 已将这份系数作为私有不可变数据，与同一 native values、原 geometry、boundary 和完整 clock 绑定；Detail 暂存数组本身不能替代正式 source。

`cells/edges/faces` 保存真实 active full sides、master/hanging 依赖、实际子面、相邻压力索引、半 cell 法向支持宽度及支持分界。`rows()` 给独立 vertex jet、master edge midpoint 值/法向导数和 cell center 的 P2 权重与审计数据。所有坐标均自有；mesh 参数只用于来源验证，不被借用。

P2 基为物理平移/尺度化 `[1,X,Y,X²,XY,Y²]`；真实 XY/RZ 面均值用充分阶数积分。incident cells 加一圈真实面邻居开始，统一扩到可辨识、四圈或分区耗尽。满秩使用原最小范数逆 A=M⁺；[H08a2r10](../task/H08a2r10_Identifiable_Lower_Space_Preparation_Report.md)在缺秩时改用 `B=A+(L−A M L)(M L)⁺`，优先保持可辨识低阶物理空间，同时保持观测最小二乘投影。L 为 XY `{1,x,y}`、RZ 径向 `{r}`/轴向 `{1,r²,z}` 与同一实际 axis/wall 条件的交集，采用正交归一系数基；物理 r² 在局部坐标的线性项不可省略。

[H08a2r14](../task/H08a2r14_Identifiable_Constant_Subspace_Design.md)将上述修正限于 ML 满列秩；ML 本身缺秩时，按物理次数先保常数，再保可观测变化量，XY 线性 x/y 联合求解。设 H 为分级低阶观测逆，则 `B=A+(I−A M)H`，仍保留原最小二乘投影。投影残差沿用原 ML 的绝对秩尺度，层级与原 rank 不一致明确拒绝。原 120532 条满秩和 1310 条完整低阶可辨识权重逐位保持；不能再以整组缺信息跳过相容常数精度。

原 rank、允许列、rings、support、coefficient norm、deficit 和 unit null polynomial 保留；新增 `priorityDimension/priorityRank/priorityDeficit` 报告优先空间的实际信息。coefficientNorm 是修正后的逆范数；优先缺秩时不声明全部再现，完整P2缺信息仍不能称P2精确。原H08a2r9的六个低阶回退已关闭，原56适用布局均通过1e-12门槛。准备层的历史任务不单独代表实际 source 切换；现行部署已由 [r11](../task/H08a2r11_Common_Regular_Source_Deployment_Resume2_Report.md) 验收，原失败记录仍保留。下一任务为 H08a 真正时空相通量积分，不再扩张准备层精修。

每个 datum/component 最多准备两种边界模式：未指定 wall，以及零 normal trace。仅实际外边界可触发第二模式，RZ 轴奇偶条件始终有效。显式 boundary 仍由原 `CchAdvectorBoundary` 绑定；`apply` 只选择已有模式并乘权重，没有 Eigen/SVD 或拟合调用。`report()` 为不可变准备库存，区别于旧 `CchAdvectorPreparationReport`。

拒绝错误 layout/momentum geometry/mesh、regrid 后沿用旧准备、非零 wall/axis 输入、非有限结果、非 nearest 或 FTZ/DAZ。符合既有 `NativeMeshIdentity` 的等价网格可验证；重新创建的 momentum geometry/layout 不自动等于原身份。原 mesh、准备局部变量、调用方数组释放不影响持有对象；regrid 必须创建新 geometry 和 reconstruction。失败不改变原对象或 FP 环境。

这是 H08a2r3 的准备层设计及 r11 的实际来源。它保持 native 面矩和 cell **积分**散度，不要求点态常散度，不引入材料守恒库存、界面 bubble、逐 dual commuting 或 R=0。完整次数 atom、相率消费者与实际 characteristic 证书已随 r11 验收；本接口不代表 H08a 时空材料积分已经通过。

H08a2r8 新增同一准备对象的解析消费（需包含 `CchSpatialAdvector.h` 以取得返回值类型）：

```cpp
const auto jet = prepared.cellJet(data, pressure, physicalPoint);
const auto normal = prepared.normalPolynomial(data, pressure, normalAxis, begin, end);
```

两者只解析原 Q4/真实面均值修正，不另建速度 owner。`cellJet` 返回物理值/J 与原 supportKnots 的闭包，内部 knot 选正侧、cell 上端选最后片；C1 保证两侧值/一阶导数一致。`normalPolynomial` 接受同一 cell 内 canonical 正切向段，严格内部 support knot 必须先切分；两端独立提取完整四次，RZ 半径方向交给既有 normal integral。不是采样拟合。

查询检查原几何 shared identity、数组库存、pressure/axis、区间及实际使用的有限系数，拒绝跨 knot、非法 FP 模式和非有限结果；只做局部计算，无 SVD/动态拟合。查询不全网格重验，调用者先以 `validateFor` 检查当前 mesh。mutable Detail 数据只用于代数暂存和检查，正式 owner 已持有私有不可变来源绑定。原 prepare/apply 数值代数保持。

[验证报告](../task/H08a2r8_Regular_Advector_Evaluation_Report.md)：144 结构组的 C1、真实面矩/Gauss 与完整 atom 通过，当时尚未切换实际 B/J/atom；现行 r11 已完成切换、运动证书和原 H10o/H10b-v2 复验。
