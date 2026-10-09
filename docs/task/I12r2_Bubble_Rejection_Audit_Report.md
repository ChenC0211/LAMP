# I12r2：bubble 拒绝可达性审计

**状态：PASS（审计完成）；父任务 I12 仍 FAIL。** 本轮只修改 1 个测试文件，未修改生产算法、容限或冻结方法。上一 goal 轮为 progress；本轮新增了决定下一修复方向的数值证据。任务格式见 [TODO](I12r2_Bubble_Rejection_Audit_TODO.md)。

## 结论

原 I12r1 的 8 个运动失败不能全部归因于通量时间取样，也不能全部靠扩大 bubble 修正来解决。

- 旋转与 vortex 的 uniform/AMR 共 **4 个失败单元**，存在原对称系数区间外的根；其体积 residual 通过现有 `1e-10*capacity` 门槛，且既有 Core 同伦证书证明整个预测→根路径无临界零、位移上界小于 `0.125 cell`。这证明“原区间不能括住目标”不等于“目标无正则的小位移解”。
- RZ 平移 uniform/AMR 的 **2 个失败**是严格纯液几何与欠填守恒 M 不一致。预测几何及独立解析体积均为满容量，通量更新后的欠量为容量的 `1.2738421241237034e-5`；不可通过允许纯液 cell 无标记造气泡来接受。
- 剩余 **2 个 XY 平移失败**仍未证明存在符合门槛的外部根：uniform 扩区预算耗尽；AMR 的一个外部试算触发 Core 的 closed-subcell-cycle 几何发现拒绝。这些结果为 inconclusive，不证明不存在解。

上述外部根仅为诊断，不返回给步进、不发布、不接受。本轮没有使任何原始 moving benchmark 通过，也不认为 volume/regularity/displacement 通过就证明曲率准确或长期输运稳定。

## 探针及真实状态

新增显式命令：

```text
build/windows/tests/Release/quantum_interface_hermite_gals_benchmark.exe --bubble-probe
```

命令在原 I12 完整步进的 `finish` 失败后只读检查真实 candidate：原 first jets、I06 mixed、predicted center、独立守恒 M 均保留。它复用公共 `discoverBubbleVolumes`、`nativeQuarticCenterInterval` 和 `nativeQuarticChangeCertificate`，未复制生产 tracer/flux/closure 求解器。

探针记录失败单元的物理 box、时间层、容量、pure/cut、预测系数、体积/导数、moment estimate、25 个显示用 local DOF、原区间和端点体积。local DOF 是 rounded 诊断值，未用于重建权威 scalar。纯相按同源证书读取容量/零，原 raw measure 留在 Core 评估中。

对 cut cell，只改变该 cell 的 center：沿目标方向以原单侧区间宽度的 `2,4,...,256` 倍试探；一旦括住目标，最多 40 次二分诊断。每次成功试算检查原来的相分量/phase、严格容量、非正导数及 `moment estimate/capacity<=1e-11`，并记录同伦证书。几何发现失败、扩区或求根预算耗尽均保留为 unresolved。证书失败只是不能证明安全，不能当作反证。

解析 circle/sphere cell volume 仅用于输出误差分解，来自独立截面积分；不反馈 M、skeleton 或 center。vortex 的独立逐 cell 真值仍缺失，reference 与相应误差字段为 NaN。

## 已认证的区间外根

表中时间为**被拒绝 candidate 的时间**，不是 accepted time。相对 residual 为 `(M_trial-M_target)/capacity`；区间倍率为 `abs(c_root-c_pred)/abs(c_original_endpoint-c_pred)`，沿实际所需修正方向取端点。

| case / cell | candidate time | 外部根 c | 相对 residual | 位移证书 / cell | 最小 physical partial | 原区间倍率 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| rotation uniform / 37 | 0.015 | -0.018482662131393587 | -1.33307e-11 | 0.0665143 | 0.110357 | 1.40417 |
| vortex uniform / 37 | 0.0315 | -0.021964765490400620 | 9.25982e-12 | 0.0722186 | 0.126160 | 1.18750 |
| rotation AMR / 72 | 0.010 | 0.37909467603515645 | 5.60897e-11 | 0.0925155 | 0.190988 | 37.0625 |
| vortex AMR / 94 | 0.0345 | -0.016187066832573569 | -7.00671e-11 | 0.0825638 | 0.163183 | 1.53125 |

原 `centerInterval` 搜索预测系数两侧的**对称区间**，要求同一 Cartesian partial 在整个可能零集管道内保持方向与下界。它先找到可证区间，并不寻找最大可达区间。现有 `changeCertificate` 则覆盖指定两份 Q4 间的整个线性同伦，允许不同子盒用不同非零 partial 证明 regularity，并证明整体位移上界。

两种证书的合同不同。本表不证明它们可以无条件互换，也不把 `retainedDerivativeFraction` 或原报告的 partial lower bound 自动解释成已满足。后续实现须明确并验证定向搜索使用的证书，保持原物理位移预算、体积/拓扑门槛和失败隔离。AMR rotation 的 center 改变量较大，尤其需要单独检查曲率与梯度；本轮未将它声明为可直接使用的生产修正。

## 另外四个失败

| case | candidate time / cell | 已观察事实 | 分类边界 |
| --- | --- | --- | --- |
| translation XY uniform | 0.03 / 29 | M_target=0；预测体积 `3.8630190218880145e-6`；解析体积 `3.8630190114595439e-6` | 预测体积误差/capacity 仅 `6.67422e-13`，M 误差/capacity 为 `-2.47233e-4`。扩区至 256 倍仍未括住目标，未证明一般不可达 |
| AMR crossing XY | 0.005 / 101 | 预测体积与解析值差/capacity 为 `3.33067e-16`；M 误差/capacity 为 `-6.96386e-8` | 超原 bracket；区间外试算出现 unresolved closed subcell cycle。不能把该失败当成所有区间外解均不安全 |
| translation RZ uniform | 0.005 / 40 | 解析与预测均纯液；M 欠量/capacity 为 `1.27384e-5` | 原 pure-cell mismatch 拒绝保持 |
| AMR crossing RZ | 0.005 / 28 | 与上项同 physical box、同容量与欠量 | 原 pure-cell mismatch 拒绝保持 |

XY AMR 和 RZ 两项发生在首个 midpoint，预测几何已与独立解析体积一致到舍入尺度；其 M 缺陷不能归咎于求根迭代没有收敛。XY uniform 的 M=0 也明显丢失了真实的非零小液相。时间通量提案仍是独立待决事项，本审计没有实施或证明该提案能自动解决所有情况。

## 验证与前后对比

证据目录：[I12r2_20261001](../../runs/baselines/I12r2_20261001/)。Release 使用既有 MSVC 19.40.33812 / CMake 3.28.3-msvc11；OpenMP 2 threads、动态线程关闭、BLAS/MKL 1 thread。两个 benchmark 进程各自使用同配置并行运行，不据此作性能优化比较。

| 项目 | 结果 |
| --- | --- |
| benchmark / architecture target 构建 | PASS，无新增 compiler warning |
| architecture CTest | PASS，0.75 s |
| probe | 205.94 s、exit 1，`BUBBLE_AUDIT_SUMMARY,8,0,8`：8 个原失败、0 个探针错误；64 次成功记录的体积/证书试算 |
| 默认完整 benchmark | 206.74 s、exit 1、10 个数值失败，完整输出与 I12r1 逐字相同；默认路径没有 probe 输出 |
| probe 行为对照 | 8 个 case 的 accepted time、失败阶段/原因、accepted error 与原运行相同，仅对象 ID 因省略静态算例而重编号 |
| 原跨面反例 | 默认完整 benchmark 内 7 行 CROSSING 与 I12r1 相同；没有 WILL_FAIL 或忽略退出码 |
| 生产与前序回归 | production 0 变更；除 benchmark 外的 17 个原 runtime 文件 SHA256 与 I12r1 相同，既有算子结果沿用已归档验证；architecture 另行执行如上 |
| 源码/冻结 | 1578 项中 1577 未变，仅 benchmark.cpp 改变；69 个 frozen 文件、12 个 integration/reference 文件、3 份 PDF 哈希一致 |

探针保留原 numerical rejection，所以 exit 1 是本次观察到的数值失败，不是 I12 成功；探针基础设施或覆盖断言失败会进入 harness error、exit 2。最初手动调用 architecture executable 漏掉所需参数，返回 usage；随后通过正式 CTest 注册命令验证通过，失败日志保留于 `architecture.log`。

原始数据和可复核表见 [bubble_probe.log](../../runs/baselines/I12r2_20261001/bubble_probe.log)、[rejected_cells.csv](../../runs/baselines/I12r2_20261001/rejected_cells.csv)、[initial_brackets.csv](../../runs/baselines/I12r2_20261001/initial_brackets.csv)、[trials.csv](../../runs/baselines/I12r2_20261001/trials.csv)、[classifications.csv](../../runs/baselines/I12r2_20261001/classifications.csv) 和 [audit.ps1](../../runs/baselines/I12r2_20261001/audit.ps1)。

## 交付与下一步

修改文件为 [benchmark.cpp](../../tests/interface_hermite_gals_benchmark.cpp)、本报告、任务 TODO、[bubble API](../api/I10_Hermite_Gals_Bubble_API.md) 和 [docs README](../README.md)。归档保存对应测试源码、运行程序、原始成功/失败日志及哈希；最终 archive SHA256 和逐条完整性校验见 [delivery_manifest.json](../../runs/baselines/I12r2_20261001/delivery_manifest.json)、[archive_verification.json](../../runs/baselines/I12r2_20261001/archive_verification.json)。

下一独立 I12 子任务可针对已发现的 **定向 bubble 安全区间搜索不足**：保留单一 cell-local bubble、固定 skeleton/M、原位移/守恒门槛，在接受任何新增区间前给出整个区间的明确证书，并验证曲率、失败隔离和完整 benchmark。此项有实际反例支持，不必等待时间积分提案后才能研究；不能简单扩大配置中的位移数值或取消判据。

[时间积分变更](I12_Phase_Flux_Time_Integration_Design_Change.md) 仍待明确决策；依 [总纲 §2.1](../00_LAMP_Development_Charter.md) 不擅自改变冻结离散。原 I12 的运动收敛、拓扑/零迹限制、动态 regrid、accepted VTK/reader、强变形及后续 Hydro/Electricity/EHD、低黏/无黏目标均未完成。Gate I 保持未通过。
