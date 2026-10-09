# H10a / 共同位置中心对流与材料动能契约验收

日期：2026-10-03。**H10a PASS；H10未关闭，下一项H10b完整无黏benchmark。** 按[TODO](H10a_Compatible_Center_Convection_TODO.md)完成一个算子事实；总体电喷雾目标、Gate H/IH/E/C保持未完成。

## 实现与数学范围

[CchKepConvection](../api/H10a_Cch_Kep_Convection.md)消费同一H09不可变材料owner，保留三层真实质量/身份、唯一共享F、完整中点clock及native layout。每条内部共同段仅一份实际面积测度中心和同分量affine重构；RZ水平段使用径向一阶矩，one-ring支持不跨不连通区域。SVD最小改动满足常数/x/y约束，对齐XY段精确保留两个0.5；无新全局速度、primary ghost、另一套rho*u或压力散度。

实现前的[H10装配设计](H10_Compatible_Center_Convection_Design.md)明确了共同位置与native动能之间的差异：直接把重构中心通量散射到两个owner并不自动满足KEP。最终内部载荷采用同一raw矩阵与其负转置配对，逐段为

\[
N_g=\tfrac12F_g[b_g(t_g^Tu)-t_g(b_g^Tu)+\operatorname{diag}(b_g)u].
\]

该相容项是明确的装配差异，不能将公开raw中心通量冒称最终载荷的直接散射；正确性和精度分别验收。共同位置对齐时恢复通常中心散射；内部动量总和抵消、常速度与H09缺陷一致，并保持半离散native动能关系。冻结PDF字节不变，原方法§7.2的共同位置/中心/skew-compatible目标保持；没有把仅满足功恒等式视为AMR空间精度证明。

外部段采用显式同层component trace，bound boundary验证同一H09 owner/身份/完整clock。result持有材料与边界，只读返回N和已经含uR的 `-N+uR`；不推进接受状态，内核无IO/字典/主时钟。报告单列边界动量、对流功、材料动能项和边界KE功，核验

\[
\dot K=\sum_i[u_i(-N_i+u_iR_i)-\tfrac12u_i^2\dot m_i]
=\tfrac12\sum_i u_i^2R_i-\sum_{boundary}sF(u_i u_b-\tfrac12u_i^2).
\]

真实非零R/G保留，材料继续采用[完全非守恒路线](H07_Pure_Gals_Hydro_Coupling_Design.md)；无库存、几何/材料修复、人工黏度、upwind或sensor。此预算是半离散材料动能关系，不是一般全离散端点KE守恒。

## 独立数值证据

| 检查 | 实测 / 预定门槛 |
| --- | --- |
| 场景 / 共同trace / native载荷检查 | 124 / 1,356 / 6,144 |
| 实际coarse–fine trace / 精确保留的对齐XY段 | 284 / 188 |
| 常数/x/y共同trace最大归一化误差 | 3.9026e−16 / 1e−12 |
| 独立dense矩阵载荷最大误差 | 4.3694e−16 / 1e−12 |
| 内部动量抵消/边界总收支误差 | 9.9597e−17 / 1e−12 |
| 独立非守恒材料KE预算误差 | 1.0793e−16 / 1e−12 |
| 常速度材料更新恢复误差 | 1.0780e−11 / 1e−10 |
| 同一冻结F下的载荷线性误差 | 7.0506e−15 / 1e−12 |
| actual FirstJet平移 / 明确拒绝 / 实际regrid | 24 / 19 / 4 |

16组代数几何覆盖XY/RZ、uniform/2:1以及全域/半域/孔洞/不连通域。独立参考用公开的单份trace系数建立dense raw incidence矩阵，再做独立矩阵转置和对角项乘法，未复用生产逐段载荷scatter。物理centroid由独立面积一阶矩核验，affine采样值来自实际native face坐标；总动量与材料动能由独立共享边关联和真实质量变化计算。归一化按gross贡献/功，生产检查256eps。

真实粗细交界段的native算术中点位置最大偏差为 **0.38034188*h**；仅散射重构通量的raw功失配为 **0.04475344**（归一化），分别超过预定负对照1e−3*h/1e−5。该控制只取两owner实际cell层级不同的段，未把外边界半控制体偏移冒称AMR证据。配对装配同时通过共同位置、动量、功和常速度检查。

24组actual几何由正式FirstJet producer推进平面、XY圆/RZ球，密度比1/1e3/1e6；最高密度比的真实非零dual积分缺陷最大绝对值27,488.90823，几何转移差异最大绝对值27,488.90823，未抹去。常速度测试直接用本项实际 `momentumLoads()` 作无外力材料更新；另用非均匀transported-component探针检查冻结F下的预算。这些探针是算子响应检查，不冒称同一速度驱动的完整非线性流体轨迹。

四组封闭无黏瞬时检查以独立streamfunction端点差构造实际native体积率，RZ使用真实2*pi测度；边界F精确0，离散散度在舍入误差内闭合。任意非均匀速度探针的归一化对流功分别为XY uniform **4.77e−18**、XY 2:1 **3.17e−17**、RZ uniform **2.08e−17**、RZ 2:1 **2.86e−17**。它们是瞬时算子检查，不是周期完整时间窗或流体PDE收敛。

实际Interface ExactRefinement和FirstJetProjection粗化后，重新生成H07a/H08/H09a/H09b/本算子，在原物理clock的step2、[.1,.2]续算，旧fine prepared在coarse mesh拒绝。这里只以常速度参考重建；一般Hydro velocity transfer仍归H15。

拒绝覆盖错误维数/nonfinite/未使用内部trace、foreign材料owner、三层identity、完整区间/step、native layout、轴速度/trace、overflow/underflow、不可识别支持及坐标/拓扑变更。原F不变，复制共享data，原mesh析构后等效mesh可用。数值范围正例用速度1e160、纯液密度1e−290验证最终有限的边界功，避免先构造超范围u²；这不是物理极端密度benchmark。

## 已发现并处理的实际限制

- uniform XY外边界的半dual段也可能不与native采样中点对齐。最初测试把整个uniform mesh都当成普通中心散射，几何前提错误；已修正为实际对齐段逐段核验，全部共同位置/动量/能量门槛保持，失败原日志保留。
- 单列隔离RZ slab的轴向native点共线，而实际面积测度中心偏离该列，无法识别径向affine值。生产明确拒绝，未跨不连通区域借点或猜梯度；负例保留，正例每个slab有两列。其他正例及actual平移采用n=4；该支持限制必须由后续生产输入预检暴露。
- 高密度比常速度净载荷接近零时，独立测试初次误用净量标度；已按TODO的gross贡献核验。生产诊断也保留组装前各项gross功，防止零边界时内部抵消隐藏实际数值尺度；256eps/1e−12门槛不变。
- actual regrid fixture初次漏传既有publication开放边界声明；已按正式Interface API显式传入，没有修改Interface或放宽拒绝条件。

## 工程验收与归档

预算3 production/build +2 tests：新header/cpp/root CMake与新短契约/tests CMake。前置1,037文件中1,035不变，只更新两份CMake、新增3文件；312 case与31 fixture资产、H01–H09/旧tests/Interface/Core/Solvers/旧Hydro/Compatibility/App/PDF及原LAMP.exe保持。新source唯一归CCH，依赖闭包只有正式Core/Interface/Eigen/线代。

- local Release定向构建通过，无新增warning/error；最初测试helper与std::bind的ADL冲突/auto类型/类型转换编译错误已修正并保留日志。
- 最终新契约 **1/1 PASS，17.31秒**（test本身17.30秒）；包括封闭域检查与数值范围正例。
- 前置/受影响回归 **18/18 PASS，243.03秒**：H01–H09、正式phase-region、FirstJet App、旧pressure-gradient load、help/architecture/self-test/layout。H01–H09选定数值行与H09b冻结证据逐字一致，App仍12smokes/49reject。
- 文档收尾后repository_layout另行通过，来源/依赖、本地链接和交付哈希另行核验归档；H09b全部24份证据和历史TODO/报告保持。
- 未运行全部82项CTest、Debug或H10b完整时间窗/空间与时间收敛/论文benchmark；未接入CCH App或完整NS求解。

一次CTest在构建尚持有exe时未启动，资源锁日志保留；待已确认的原构建session完成后才重新测试，没有重新启动构建或以Not Run计PASS。其他失败日志、TODO-before、design-before和最终原始数据均保留。

证据：[verification](../../runs/baselines/H10a_20261003/verification.json)、[operator detail](../../runs/baselines/H10a_20261003/operator_detail.log)、[regression detail](../../runs/baselines/H10a_20261003/regression_detail.log)、[source audit](../../runs/baselines/H10a_20261003/source_audit.json)、[交付清单](../../runs/baselines/H10a_20261003/delivery_manifest.json)。临时审计脚本归档于baseline/engineering_scripts。

## 下一步与未完成范围

下一单任务 **H10b**：在 `cases/benchmarks` 调用本生产算子，执行光滑无黏/无界面力完整时间窗、独立解析空间/时间收敛、边界KE功/动能漂移和donor对照；特别验证相容项在AMR/边界的空间一致性。只有对应证据通过才能关闭H10。H11受控耗散、H12–H13毛细、H14corrector、H15兼容Hodge/生产边界/速度transfer以及IH/EHD仍未完成；本项不报告水平集曲率。
