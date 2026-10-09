# I12a15r3 / 权威 Q4 的仿射拉回包围报告

**本项 PASS；父 I12 / Gate I 仍 FAIL。** 新算子与 T0 验收通过，原 uniform shear 积分失败保持。见 [TODO](I12a15r3_Authoritative_Q4_Affine_Pullback_TODO.md)、[设计](I12_Authoritative_Q4_Affine_Pullback_Design.md) 和 [API](../api/I12a15r3_Authoritative_Q4_Affine_Pullback_API.md)。用户授权的几何—库存分离继续保留；本项不改变 M、库存权重或几何发布策略。

## 同一个权威多项式

新 Snapshot 入口从原 native leaf 的 `local.numerator[25] / local.denominator` 开始，使用原 exact Bernstein 变换和精确有限差分转到 chart power basis。两轴 Bernstein 变换共清除 144，因此分母取 `144*原 denominator`。精确为零的高次幂系数仍为零，没有取 rounded local DOF 缓存、拟合、截断或新几何 DOF。

在共同幂二尺度平衡后，对每个分子和正分母给出 outward 区间。四个端点商的最小/最大值分别向外取一个相邻 binary64 数，包含正确舍入除法及 gradual underflow；不可表示为有限范围时明确拒绝。严格编译的原 affine engine 继续包围 physical→chart 变换和二维 degree-4 Horner。cell 宽度采用真实 double 端点的精确实数差的区间，不能换成舍入的网格 spacing。

每个返回 form 为 `c+s*q_s+t*q_t+R`，同一组 q 同时作用于 B 的两个分量。区间系数的不确定性及全部非线性积保留在原仿射算术里，输出时再将系数代表值误差移入 R。其包含关系来自基变换恒等式与逐操作区间归纳；下列样本用于检测实现，不替代全域证明。

叶外求值只包围该叶多项式的延拓。将符号用于实际相场时，仍须证明真实 B 位于对应叶，且所有参数点有完整 closed donor 覆盖。接口不生成 Interface 来源 capability，不处理外部边界历史，也不以这个延拓修改原 geometry。

## T1

新 contract 通过 **232 个区域、33,498 次实际正时长回溯、1,998 次旧时刻恒等映射、57,470 次闭叶检查、1,157 个 Core 多项式检查和 20 类拒绝检查**。范围包含性容差为 0。实际 tracer 仍逐次核验单 substep、无 boundary hit、Jacobian 可用；对每个实际 B，遍历每个包含 B 的闭叶验证该叶 scalar，同时核验原 authoritative negative predicate。

| 组别 | 覆盖 |
| --- | --- |
| XY uniform/AMR | 两种实际 face 方向；zero、axial、shear、rotation、反向 rotation、cubic swirl；plane、Q2、包含 x⁴y⁴ 的完整 Q4、薄层、内部液泡、相切 |
| RZ uniform/AMR | 两种实际 face 方向；zero、axial、strain；满足轴正则性的相同六类 scalar |
| Core 独立参考 | 双参数 affine 输入及非零 residual 两端组合；完整 biquartic 解析幂表达式；subnormal 系数保留 |
| 精确细化 | 非 dyadic 平移域，uniform 与 AMR 两种 target；原 source、exact refinement、再次 exact refinement；不重新采样解析场初始化 target |
| 隐藏事件 | 三个节点都为气相的薄层穿越；内部液泡角点/中心都为气相；精确相切的零值保留；薄层局部缩小后可认证负值 |
| 拒绝 | invalid index、inactive/nonquartic/moved source、非有限系数、反向范围、非法 bounds/parameters、非 nearest 舍入和 FTZ |

测试中的 `nativeQuarticLocalDofs` 仅用于查询前后不变性比较，不作为任何证书或解析参考的输入。非 dyadic source 的原样 scalar 值用于精确细化比较；没有将非 dyadic 采样初始化误差误称为全局解析 Q4 的精确重现。

精确 plane 经 Q4→power 转换后的 remainder 为 `3.6776137690708355e-16`，没有因 Bernstein 依赖性制造曲率。相同中心的 face/time 半宽同时减半五级：

| 级别 | shear + plane 的 R | shear + Q4 的 R |
| --- | --- | --- |
| 0 | 3.9672851562558043e-4 | 4.5567811964090292e-4 |
| 1 | 9.5367431641186603e-5 | 1.0846511577172012e-4 |
| 2 | 2.3365020752506282e-5 | 2.6451258694376784e-5 |
| 3 | 5.7816505437620687e-6 | 6.5307155706980453e-6 |
| 4 | 1.4379620557535454e-6 | 1.6224801733634691e-6 |

最小相邻比 `4.0207253874527789`。这是局部证书随参数域缩小的性质，不是完整 Interface 的时间阶或原 shear 积分已通过。

## 初始失败

首轮 production 编译暴露同一 `auto` 声明中类型不一致，以及私有 interval namespace 在不同 TU 的可见性/类型不匹配。分别拆开声明，并用正分母的 outward 四端点除法直接表示精确系数商；没有引入另一套根/相位判据。

第一轮数值运行的 XY uniform/AMR 已通过，但 RZ 夹具的平面/偏心曲线给出了非零轴径向导数，原 factory 正确拒绝。改为轴偶对称多项式，保持原轴硬约束。第二轮完成所有数值主体后，legacy cubic-grid 拒绝包含隐藏闭环的 guard 夹具；该 guard 只需验证 nonquartic source 类型，改用合法纯相 cubic 夹具。第三轮和最终直接运行通过；所有原失败日志保留，未放宽算法、容差或验收预算。

## T0 与交付

修改范围为 5 production（Slab h/cpp、Snapshot h/cpp、原 QuarticLeafData header）与 2 test-related（新 contract 与 tests/CMakeLists）。新增的唯一 production include 是 QuarticLeafData 对 Core PolynomialVelocitySlab.h 的依赖，供同一仿射值类型与原语使用；它不改变库依赖方向或原 scalar owner。architecture 通过，include 计数从 3114 增至 3115；12 targets、383 TU、907 files、5 exact legacy edges 及依赖规则保持。单列这一接口所需的显式变化，不声称 architecture 输出逐字不变。

最终直接运行 0.8181616 s、exit 0。全部 39 个 runtime 目标重建成功，398.4343939 s，所有本项构建日志无 compiler warning；Core 与新 contract 的 strict/no-LTO 设置已核验。初始构建记录因 PowerShell 嵌套数组误记 targets=2；原日志逐项确认全部 39 个目标完成后修正计数，原 JSON 与日志 SHA 保留，未重写实际构建结果。

25 项 CTest 用时 162.4706916 s，**24 PASS、1 原有 FAIL**。原 23 项数值测试的完整输出逐字相同，新增 affine phase contract 通过；唯一失败仍为 `interface_bounded_spacetime_flux_contract` 的同一 uniform shear：E=`6.6763159339198168e-6` 大于原 `5.859375e-6`，4096 panels、pending=0。旧 crossing probe 的 exit 1 与完整输出也相同。

81 个 frozen 文件、3 份方法 PDF 及 r2 原 archive 保持 SHA 身份。第一次最终审计即通过。归档包含 5 production、2 test、5 份文档、39 executables 与原 DLL，逐项核验 archive/live payload SHA256。

证据：[目录](../../runs/baselines/I12a15r3_20261002/)、[验收](../../runs/baselines/I12a15r3_20261002/verification.json)、[算子](../../runs/baselines/I12a15r3_20261002/operator_verification.json)、[旧回归](../../runs/baselines/I12a15r3_20261002/regression_comparison.json)、[范围审计](../../runs/baselines/I12a15r3_20261002/source_audit.json)、[构建目标核验](../../runs/baselines/I12a15r3_20261002/build_target_verification.json)。本项没有修改实际 velocity/RK2、原 strict sign、M/库存权重、reflux 或 benchmark 配置。

## 下一项

在完整 closed donor 覆盖下，用本项 scalar 的一次项及严格余带认证分块内明确液/气子区域，计算其 XY 面—时间面积及 RZ 径向矩，进而收紧原 gross 液体积分界。保持同一 source/库存权重、原 panel/depth/误差预算，再验同一 uniform shear。完整 driver、时间阶、动态 AMR 和后续 Hydro 仍须各自通过验收。
