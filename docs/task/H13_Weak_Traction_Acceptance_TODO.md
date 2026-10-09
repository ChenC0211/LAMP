# H13 / 授权联合候选的完整弱牵引与原物理门槛接续验证

**当前授权081136（2026-10-09）：** 用户明确回复“授权”，批准[phase-aware source/物理动能Hodge联合方案](H13_Phase_Cut_Kinetic_Hodge_Method_Change.md)。下面原四文件/M冻结和授权等待段保留历史；现行范围以本次接续格式为准。仍仅H13，FAIL；H14/H15未启动。

**当前修订授权：** 用户已明确授权独立切向模态的初值、演化及动能/弱力/黏性/材料能量联合耦合。q作为同一source的组合速度状态，覆盖下文旧“无独立切向状态”限制；24生产/10已有契约预算和其余原约束保持。1030单格证明不能推出全邻域阻断，实际支持四组合rank2证据已形成。

## 2026-10-09 / 独立切向模态授权后执行格式

本接续的对流执行记录：仍同一H13，目标为唯一材料flux owner的具名物理通量视图与完整[U,q]体对流双线性形式。允许仅使用原清单中的 `CchMaterialTransport.h/.cpp`、`CchKepConvection.h/.cpp` 及原momentum-geometry契约；不增加production文件。物理flux必须由原midpoint PhaseRate的同一source与真实density取得，不另设F。原始体对流矩阵不强行反对称化；独立curl carrier的散度、界面normal和外边界flux为零时，真实积分应自行给出零对流功，并在完整H的midpoint演化中核验。非零边界/材料R及移动基的完整配对仍须随后接入，不能用此子步判定H13 PASS。输出归103518目录。

任务 ID / 标题：同一H13 / 独立切向模态与完整动能块演化。
目标：以显式q初值和演化承载内部切向运动，联合验证同源物理场、H完整块、法向弱功、黏性及移动材料能量；最终仍关闭全部原H13精度/稳定/效率gate。
允许修改：原24文件内本步先使用 `CchSpatialAdvector.h/.cpp`、`CchMomentumGeometry.h/.cpp`、`CchSinglePhaseSaddle.h/.cpp` 和现有momentum-geometry/saddle短契约；其他已批准文件按必要消费者逐步接入。不另开task或新增production文件。
禁止修改：native U/p定义和D/C、FirstJet/characteristic公式、物理阈值及无关模块保持；q不能被默认置零、丢弃或静态伪装成U派生值。禁止仅接H却继续旧能量/材料配对。
数学约束：`x=(U,q)`、`R_phys=R_base(U)+Tq`；T是同Q4分相curl、native面均值和interface normal皆0；`C_ext=(C,0)^T`。`H`含所有交叉块且SPD；黏性同源非负，exact mu=0零；固定几何midpoint独立动力模态保持/功耗散恒等式，移动H的metric变化与唯一F/R/uR及模式转移余项显式保留。
实现要求：先完成实际邻域rank纠正，再构造带真实phase绑定、显式q初始化及解析XY/RZ first-jet的内部模态子空间和H/压力saddle全块。内部bubble只是第一子空间，共享trace模态与phase-P2/P3最终再现不能由其替代。
验收：同源/clock/模态数量拒绝、面均值/normal/散度/轴/边界jet、真实体积分H块与交叉能量、固定几何无黏保持/强迫功/黏性耗散；随后原T0/T1全部与T2/T3，未完成项不计PASS。
失败规则：H13全范围未通过保持FAIL，数学或实现困难继续同任务；不再申请已授权模态方法，不启动H14/H15。
交付：组合状态/API、初值与方程、真实分相H/弱功/黏性/材料诊断、完整回归和未完成项，输出归103518 D:目录。

修订授权的必要接线文件：独立q的物理phase atom rate必须进入唯一 `CchPhaseRate` owner。当前仍至多24生产文件，将原清单中尚未修改且本方案直接复用的 `CchEnergyViscosity.h/.cpp` 两个名额改为 `CchPhaseGeometry.h/.cpp`；同源bulk黏性接线放在已允许的PhaseViscosity owner，不复制EnergyViscosity算法。此为用户明确要求“材料能量一并耦合”的直接依赖；预算不增加，原EnergyViscosity文件不修改。其余原清单及10既有契约保持。

## 2026-10-09 / 103518 同一H13接续执行格式

任务 ID / 标题：同一H13 / 既有邻域观测下的分相source局部可识别性复核及已授权构造。
目标：先核验1030单格反例是否覆盖实际prepared邻域；若只证明单格欠定，则修正过宽结论，按081136原授权继续相内多项式/分相curl空间构造，不以未经证明的阻断要求新增状态。
允许修改：现有24生产/10契约清单保持；先只修改既有 `cch_momentum_geometry_contract.cpp` 的直接source支持/独立分相矩秩断言及H13文档，定位证据归 `runs/checks/H13/20261009_103518`。后续正常构造使用原批准清单，不另开task。
禁止修改：原T1/T2/T3、native U/p数量、accepted FirstJet/characteristic/Jacobian及唯一F/R/uR等原约束保持；未经证据和明确方法决定不新增独立切向状态。
数学约束：区分仅四面观测矩阵与实际邻域完整观测矩阵；采用真实XY/RZ面测度及AMR子面，不从单格零均值推断全局U为零；必须给出rank、支持外非零观测及独立常场系数复现。
实现要求：保存上轮原有理数脚本/证据；新增实际prepared support独立观测检查，做四面局部负控与完整邻域正控。仅证据支持后修正roadmap/报告中的方法阻断结论，继续原授权source构造。
验收：XY/RZ × uniform/2:1 AMR四组合的实际source支持复核与独立分相常场观测矩阵；完整原momentum-geometry短契约及后续直接相关检查。支撑可识别不等于source精度、Hodge配对或H13全物理PASS。
失败规则：无法识别的配置提供确切source支持及null证书；可识别配置继续已授权实现。保持H13 FAIL直至全部原gate通过，H14/H15不启动。
交付：实际支持/秩/常场系数/负控证据、必要纠正文档与原批准方案的下一步实现；所有输出在D:。上轮改变了源码并产生新证据，分类为progress；本轮复核其推论，尚无经全邻域验证的外部阻塞。

## 2026-10-09 / 081136授权后接续执行格式

任务 ID / 标题：同一H13 / 分相物理source—局部动能Hodge—弱力及材料能量联合改造。
目标：关闭一般cut的准确normal、真实分相动能与应力回复之间的失配，同时实现原无黏/低黏精度、稳定性及效率验收。
允许修改：以下24个精确production文件为本次预算上限，未需要的文件不改。依次为 `include/Modules/Hydrodynamics/Detail/CchRegularAdvectorGeometry.h`、`src/Modules/Hydrodynamics/CchRegularAdvectorGeometry.cpp`、`include/Modules/Hydrodynamics/CchSpatialAdvector.h`、`src/Modules/Hydrodynamics/CchSpatialAdvector.cpp`、`include/Modules/Hydrodynamics/CchMomentumGeometry.h`、`src/Modules/Hydrodynamics/CchMomentumGeometry.cpp`、`include/Modules/Hydrodynamics/CchInterfaceNormalLoad.h`、`src/Modules/Hydrodynamics/CchInterfaceNormalLoad.cpp`、`include/Modules/Hydrodynamics/CchPhaseViscosity.h`、`src/Modules/Hydrodynamics/CchPhaseViscosity.cpp`、`include/Modules/Hydrodynamics/CchEnergyViscosity.h`、`src/Modules/Hydrodynamics/CchEnergyViscosity.cpp`、`include/Modules/Hydrodynamics/CchDualMaterial.h`、`src/Modules/Hydrodynamics/CchDualMaterial.cpp`、`include/Modules/Hydrodynamics/CchMaterialTransport.h`、`src/Modules/Hydrodynamics/CchMaterialTransport.cpp`、`include/Modules/Hydrodynamics/CchKepConvection.h`、`src/Modules/Hydrodynamics/CchKepConvection.cpp`、`include/Modules/Hydrodynamics/CchSinglePhaseSaddle.h`、`src/Modules/Hydrodynamics/CchSinglePhaseSaddle.cpp`、`include/Modules/Hydrodynamics/Detail/HydroMomentumAssembly.h`、`include/Modules/Hydrodynamics/HydrodynamicsMomentum_impl.h`、`include/Modules/Hydrodynamics/HydrodynamicsCoupled_impl.h`、`include/Modules/Hydrodynamics/HydroCoupledBalance.h`。新局部数据类型先置于其现有owner，不增加编译边界文件。
测试范围：`tests/contracts/hydrodynamics/`下现有 `cch_momentum_geometry_contract.cpp`、`cch_capillary_traction_contract.cpp`、`cch_characteristic_contract.cpp`、`cch_single_phase_saddle_contract.cpp`、`cch_phase_viscosity_contract.cpp`、`cch_energy_viscosity_contract.cpp`、`cch_dual_material_contract.cpp`、`cch_material_transport_contract.cpp`、`cch_kep_convection_contract.cpp`、`cch_phase_geometry_contract.cpp`，上限10份。仅直接必要断言；完整研究归现有 `cases/benchmarks/cch_capillary`。文档及 `runs/checks/H13/20261009_081136`证据不新增task。
禁止修改：Interface accepted字段、characteristic/Jacobian积分公式、电场/无关模块、字典语义、原物理阈值；不引入第二primary速度/相库存、case专用公式、密度或曲率滤波、最小黏度、R=0及fallback。
数学约束：原native面均值与D/C负转置；source唯一owner，分相物理normal与C1几何normal一致；局部H对称正定、真实分相积分与常场总质量；完整J0矩阵及同源弱力；移动H的时间能量项与唯一F/真实R/uR配对；mu=0精确零。旧K_diag/失败单列，不将新能量定义冒称旧gate通过。
实现要求：先冻结当前文件/原build身份，维护两行KEP与两组可信reference；随后构造/验证source局部基、Hodge与材料变化配对，最终接入同一实际推进。prepared辅助局部消元，不加全球高阶DOF。不能只接H并沿用未经核验的旧配对。依赖不足按实际具体原因处理，不重复申请已批准范围。
验收：T0标准Release/architecture/repository_layout，无新warning；T1真实face/分相多项式、AMR/RZ轴及跨相、H正定/真实体积分/物理功、C/完整常应力、全部来源拒绝；T2原36平面响应>=1.8阶、72几何/1728静滴、H08a28窗口与H10b84原范围；T3原完整周期/三档空间时间/AMR/RZ/high-density/strict mu=0/低黏和效率。退出0与部分绿检查不代替物理PASS。
失败规则：数学失败保留H13 FAIL并继续同一任务；真实外部阻塞方可BLOCKED，不将困难、耗时或未完成当阻塞。全部完整H13范围通过前不启动H14/H15。
交付：真实修改、分相/惯性/能量诊断、原短契约及完整研究、来源归档、报告/API与未解决问题。所有工作/临时/输出位于D:。

## 2026-10-09 / 1023 H13接续执行记录：phase-cut Hodge求积核

任务 ID / 标题：同一H13 / 分相物理source—局部动能Hodge—弱力及材料能量联合改造；本记录仅覆盖其中局部切分相矩求积子步，不另开task。
目标：让局部phase-cut Hodge的分相矩积分在平面/曲面、AMR/RZ切面上具有可审查的界面事件处理、误差账本及有界工作量。
允许修改：仅本授权清单中的 `src/Modules/Hydrodynamics/CchMomentumGeometry.cpp` 与 `tests/contracts/hydrodynamics/cch_momentum_geometry_contract.cpp`；本轮未新增production/test文件，日志写入 `runs/checks/H13/20261009_081136`。
实现：将切分限制在单一原生Q4图表并按规则支撑结点切盒；按可认证偏导选切片轴；用Q4盒符号检测界面事件。根包络的chart方向由两种候选物理映射的原生Q4符号证书判定，避免把chart局部参数误按输入边方向映射。界面不穿过的子盒由同一Q4符号认证；GL8对最高双次数8的162个相矩（RZ乘径向Jacobian后x次数至9）精确积分，舍入预算显式加入；cut子盒仍比较GL8/GL16并递归，单盒65,536次slice查询上限和原 `1e-11` 几何预算不变。对单点Q4盒符号为0的求积点不强行赋相；相归属取自同源slice的保留phase range，若点盒符号已解析且与range矛盾仍拒绝。
验证：Release增量构建 `lamp_hydrodynamics_cch` 和 `quantum_cch_momentum_geometry_contract` 成功，构建日志未出现warning/error；新增Hodge-only八组合（平面/圆形 × XY/RZ × uniform/AMR）全部PASS，最大报告分相体积差 `4.44e-15`，均通过对称、正定、真实总质量常场动能及稀疏矩阵作用检查。完整 `cch_momentum_geometry_contract` 退出0，耗时54.95秒；日志见 [完整契约输出](../../runs/checks/H13/20261009_081136/h13_cch_momentum_geometry_contract_full.log)，构建日志为 `h13_cch_roottable_build.log` 与 `h13_contract_cch_roottable_build.log`。根映射单例也以 `--kinetic-hodge-circle-amr-only` 独立退出0。
当前状态：**H13仍FAIL / 联合候选未完成**。本子步中的H矩仍以现有 `regular.cellJetWeights` 在两相共用同一C1基；尚未实现phase-enriched物理source及界面切向跳跃，也未把H接入生产动量/弱力/真实材料通量与黏性能量配对。完整T0 architecture/repository_layout、原物理T2/T3及全部联合验收均未完成；不得将本子步PASS误报为H13 PASS，H14/H15仍不启动。

## 2026-10-09 / 1030 同一H13的source可识别性审计

任务 ID / 标题：同一H13 / 分相物理source—局部动能Hodge—弱力及材料能量联合改造；此项是已授权方法的必要可识别性检查，不另开task。
目标：核验仅由现有native normal face means和同一几何能否唯一恢复含切向跳跃的分相速度场。
允许修改：本轮仅在授权的D:证据目录新增精确有理数审计；无production/永久测试修改，无新task。
方法与结果：单元方格、界面 `x=3/8`，两相常速度 `(0,5/8)` 与 `(0,-3/8)`。四个native面法向均值全为0，与零速度场完全相同；两相散度皆0，界面法向迹皆0，但切向跳跃为1。精确观测矩阵存在null向量 `(0,0,5/8,-3/8)`，对应解析及独立有理数脚本见 [可识别性审计](../../runs/checks/H13/20261009_081136/phase_source_identifiability_audit.md) 和 [复算脚本](../../runs/checks/H13/20261009_081136/phase_source_identifiability_audit.py)。这是P0（故属于P2/P3）且可置于周期y方向、无穿透x边界的静止无黏涡片。
结论：任何只以 `U` 和几何为输入的确定性分相source（包括固定 `curl(phi_+ W)` 补空间）都无法区分这两个不同的物理场。若要保持T1切向跳跃再现，缺失信息必须由额外切向模态/状态或明确的含历史局部演化方程提供；仅把该系数静态消去成 `U` 的函数不能保留自由无黏涡片。当前批准方案的“只用native normal速度、无额外切向状态”与T1一般切向跳跃再现要求相冲突。
状态与下一步：**H13仍FAIL**。不降低T1/T2/T3、不将切向跳跃测试弱化为可观测子空间。生产source/Hodge耦合需先审查并明确如何携带该模态；在原批准方案内不能诚实构造唯一系数，故当前不臆造来源或接入真实动量方程。H14/H15不启动。


任务 ID / 标题：H13 / 同源 native 面矩运动—完整 Hermite 毛细弱牵引候选验证（2026-10-08 23:35，Asia/Shanghai）。

目标：只解决一个数学事实——当前共同 source 与完整弱牵引是否在原压力/材料质量/FirstJet 运动下同时满足常应力矩阵、真实虚功及原无黏能量门槛。仍属路线图 H13，不启动 H14。

允许修改：原授权的四 production 文件 CchRegularAdvectorGeometry.h/.cpp、CchInterfaceNormalLoad.h/.cpp；仅扩展已有 source/traction 两份契约；现有 cch_capillary benchmark 必要验证绑定、H13 TODO/报告/API/状态索引；证据归 runs/checks/H13/20261008_233555。先验证当前已存在候选，仅在可复现错误确定后作最小修正。

禁止修改：Interface accepted state/边界回调、characteristic/Jacobian 积分、M/C/D_I、F/R/uR、KEP、黏性、其他模块、既有物理门槛。不得另开 source 精修、加最小黏度、曲率滤波、case 分支或 fallback。旧报告和旧失败结论不改写。

数学约束：同一 prepared source / publication / layout / clock；C=-D_I^T；完整 J_f(J0)=J0 I，F=J_f C alpha；同源局部 basis 与 actual apply/cellJet 一致；常量反力显式保留，真实 ds/2*pi*r ds 与 AMR 子面通量；mu=0 严格为零；Pure FirstJet 无独立库存，实际材料缺陷不隐去。

实现要求：保存当前候选源码/输入/可执行与哈希；旧 runs 不存在如实记录，不凭文档声称复核历史归档。独立 D: 构建，复用已安装 D: 依赖；先 T0/T1，再实际无黏周期和静滴原扫描。若短契约或原周期出现实质失败，保留失败和接受前缀，定位最小反例后裁决，禁止用退出0替代物理验收。

验收：T0 CMake Release 编译无新增 warning、architecture/architecture_self_test/repository_layout；T1 原 source/traction、完整常量矩阵/P1/P2/jump、来源拒绝、压力/材料/特征线受影响短契约；T2 原实际 n8 无黏周期（256步），原能量定义及频率不变，独立保存即时功、质量缺陷和全部状态；通过后再完整空间/时间/AMR/低黏及静滴72几何/1728解。末档频率阶>=1.8，原静滴/能量 gate 保留。

失败规则：数学或 gate 不通过记 FAIL；外部缺失接口才记 BLOCKED。任一必要前置失败，不将后续未运行记通过；H13 全范围未通过前仍 FAIL，不采用候选，不启动 H14。

交付：本 TODO、实际改动（如有）、构建/测试和物理数据、独立核算、来源哈希、H13 报告与剩余问题。所有工作/临时/计算输出存 D:。

授权前045122检查点状态：**FAIL**。36组条件端口与36组原M投影完成；准确normal与伴随弱力下，真实native对角动能的一阶缺失给出响应下界。详见[动能证据](H13_Phase_Cut_Kinetic_Metric_Certificate.md)及[联合方法方案](H13_Phase_Cut_Kinetic_Hodge_Method_Change.md)。该检查点production/permanent test零修改，原一般cut反例及六组静滴回升保持FAIL；永久characteristic/KEP两项仍未关闭，原完整静滴及n8周期数据保留。完整三档/AMR/低黏/材料验收未闭合，H14不启动。当前授权后进展见上方1023记录。


## 2026-10-09 / 003903接续执行格式

任务 ID / 标题：同一H13 / 原高密度常速度反例的系数—求和因果分解。
目标：只解决一个数学事实——失败来自新source规则、prepared KEP系数残差，还是相同算子的浮点求和，确定保持原方程/门槛的最小修复范围。
允许修改：runs/checks/H13/20261009_003903中隔离重放/高精度原系数核算/明确标识的旧source语义对照；已有H13 TODO/报告/具体设计说明。production/test先零修改。
禁止修改：当前production source/KEP/F/R/uR/M/C、characteristic/Jacobian/FirstJet、接受状态、原输入及阈值；不启动H14、不增加fallback或原子级常速度特判。
数学约束：逐原共享segment通量及原skew配对；真实begin/end mass、非零材料缺陷、原常速度更新限值1e-10保留。高精度只复算保存的实际二进制系数，不参与接受状态。旧source对照仅移除新normal-port消费，保留原P2行/父子约束，不能冒称旧字节档案恢复。
实现要求：核验上轮封存SHA与当前四production身份；导出actual F、weights、couplings及所有DOF；独立Decimal分离prepared行残差与apply求和误差。仅证据支持时形成具体算术修复提案并隔离验证，不由本轮自动解除四production/两契约及冻结KEP范围。
验收：T0来源/构建与架构；T1原204行/三密度失败逐值重放，高精度原系数核算、旧source单变量对照和负控；如形成候选，原KEP/材料/能量门槛全部保持，不能仅凭一行通过采用。
失败规则：任何预算/模型扩展先登记具体方案；原FAIL保存，不调阈值或选择依赖版本换通过。仅外部接口缺失记BLOCKED。
交付：实际原系数、因果分项、可复算程序/结果及报告；H13完整未通过前仍FAIL，目标active。

003903接续补充：80位Decimal原二进制系数核算证明prepared残差占最差常速度误差约98.50%，普通运行时补偿求和不足。旧source语义控制通过旧门槛但残差仍大；原KEP两行内部残差精修隔离候选通过原204行及完整124场景/6144 DOF/24 crossing/57耗散场景，生产采用仍待第五文件范围例外。同一实际Q4的独立100位30阶Taylor/Gronwall/Cauchy参考、旧source正控、34阶及两负控完成；仅同步两个reference tuple的暂存原完整characteristic契约通过，永久数据同步另须第三测试文件预算补充。原完整静滴72/1728扫描以当前实际source/载荷独立诊断运行，数据完成仍需独立收敛审计，不预先迁旧PASS。H14不启动。

003903原物理验证范围：在等待两个明确预算例外期间，当前source/载荷和tracer均未修改。原72几何/1728静滴扫描已自然退出0，下一步读取原始数据独立核算/收敛，不由退出0宣布PASS；同一当前程序重新执行n8原256步完整无黏周期独立诊断，原研究流程不调用被冻结KEP，source/载荷/边界/初值/原参数保持。该运行依据已确认终态的前次86步停止记录新建目录，不因观察超时重启活进程。两个待批准文件例外不在此运行中采用；不把T2数据反向覆盖原T1失败或称H14通过。

003903保存数据source审计执行格式：仍为同一H13；目标为判断六组静滴回升是否来自实际弱牵引/source支持或保存数据/装配误差。允许只在本证据目录新增source-only只读重放，逐一核对保存的mesh cell/face序号、坐标及连接后，以保存quadrature导出实际局部basis；不构造第二份phi/curve/相体积，不解新接受状态。独立Python重算72几何的变应力载荷，并用保存的27组XY uniform非多项式无黏场检验稀疏basis与actual apply/cellJet一致，定位最差面、质量及虚功。禁止production/permanent test修改、力过滤、M/C改动和新门槛；T0只读构建，T1载荷/原场重放与扰动负控，T2沿用已完成原扫描；任何原非收敛保留FAIL。n8完整256步已自然退出0，独立原能量/功/频率复算完成，尚无完整三档频率收敛结论。

## 2026-10-09 / 035647接续执行格式

任务 ID / 标题：同一H13 / 六组静滴回升的原载荷分项响应。
目标：一个数学事实——以原M/C/鞍点线性映射判断回升由压力相容平均项、曲率偏差弱项还是两项叠加产生。
允许修改：runs/checks/H13/20261009_035647只读诊断、已有H13任务/报告/状态索引；保存003903原场及basis只读输入。先零production/permanent test修改，不采用待批两补充。
禁止修改：source、曲率/界面quadrature、M/C、KEP/F/R/uR、characteristic/FirstJet、阈值；不得将分项或解析去常量研究力当production候选。
数学约束：F=-kbar*a+sum W*b*(k-kbar)；只为因果诊断减去已知圆曲率常量压力通道，三分项之和恢复实际F，原strict mu=0时A=8*M、C=-D_I^T、真实边界、压力均值与实际气相质量保持。
实现要求：同一9几何/27组保存原无黏场；用原public saddle对平均偏差、弱项、实际完整载荷分别求响应，重放并记录叠加/原场差、峰值位置/质量、分项随n的响应。实际load二进制舍入差保留，不改写六个原失败。若分项不能复原原场先定位数据或数值问题。
验收：T0实际来源哈希和只读构建；T1完整27组原场复现/线性叠加及故意污染负控；T2使用原72/1728结论，不重复扫长运行或迁为PASS。
失败规则：仅方法证据可推动下一步；超出授权的实质方法变更先登记具体离散性质、可比性及原验收，自动续跑不是授权。完整H13仍FAIL，H14未启动。
交付：分项力、原求解器结果与独立因果统计、原门槛、归档及报告；全部存D:。

035647结果：9几何/27原无黏场、7种原力分项共189次原public saddle求解自然完成，完整原场逐位复现，两个分解叠加误差最大1.61e-20。原μ=0 A=8M、全部public saddle门槛不变。n16的共同应力偏差响应与局部体积余项强相消，n32相消减弱；不能把六组回升归为浮点、独立几何或删除弱项就能修复。只保留平均项在rho=1e3仍末档-1.758阶；共同参考诊断即使单调，也不代表可替换局部support方法。0 production/permanent test修改，原两项待批预算补充保持，H13仍FAIL。

035647顺序接续执行格式：仍为H13，前一载荷分项诊断已终态；只检验一个事实——精确平面/余弦应力的原M/C/native response与实际共同source的界面normal精度是否受cut位置和密度比影响。只允许本目录制造解与原public算子诊断；原source/load/压力/真实相质量/边界均保持，不构造接受的未来形状。平面phi=y-H严格FirstJet，J=-cos(2*pi*x)为制造应力，解析液/气势流满足上下零normal、两侧零normal、界面normal连续及p_liquid-p_gas=cos。H=.5+theta/n，theta=0/.25/.5；rho=1/1e3/1e6；n8/16/32/64。不将制造应力冒称实际曲率或原物理波验收；导出真实原生面、相体积、求积、实际/真值source query，分别量化native误差与source插值误差。T0只读来源/构建；T1解析PDE/边界/跳跃独立核验与完整输入矩阵/扰动负控；T2仍使用原完整H13 gate，任何新失败保留。0 production/permanent test修改，H14不启动。

035647平面结果：36组自然完成并独立核验，精确alpha/M逐位吻合；真native face经source查询的off-grid normal末档约0.995阶，内部及模态同为约一阶。rho=1e6时theta=.5的n64实际normal模态/解析值=-10.4379，pressure不做功，表面梯度与加速度正配对。n8/.5/1e6已单独约1秒复现，与原case9逐字一致。此反例需要共同source一般内部normal与原M/C/弱力相容的方法推导；不能靠待批的两项算术/参考补充、曲率过滤或最小黏度关掉。具体实质方法扩展须先形成可复核方案并按总纲§2.1登记，不自动扩大原normal端点范围。

## 2026-10-09 / 045122接续执行格式

任务 ID / 标题：同一H13 / 两相连续normal矩与原质量的条件端口推导。
目标：一个数学事实——在上轮精确plane反例中，保留真实M/C和沿面通量迹，用连续normal、两侧P2及原逆质量代价能否同时恢复normal插值与应力回复的正确密度/位置精度。
允许修改：runs/checks/H13/20261009_045122隔离参考/原public saddle诊断、已有H13 TODO/报告/方法证据；修正文档的重复最小复现段。0 production/permanent test修改，两个待批补充保持。
禁止修改：当前共同source/内部规则、实际界面/曲率/quad、M/C、F/R/uR/KEP、characteristic/Jacobian/FirstJet、原阈值；不得把条件Γ trace冒称实际全source或production修复。
数学约束：沿面迹为原native FV P3端点+Q4面均值；法向P=[1,d_-,d_-²,d_+,d_+²]；P^T*w=e0，最小化w^T*M^-1*w，由SVD计算而非拟合响应。应力与normal使用同一条件B及真实quad；常应力完整I/gauge由明确公共常量通道处理，是否满足一般局部support/全C1尚未预先证明。
实现要求：只读035647完整36组原字段；独立FV多项式重组与原aligned source query核对；生成normal矩权重及法向测试B，原public saddle验证同参数矩阵。来源/负控/全部失败保留；phase/M绑定超出原geometry-only规则，参考通过也不自动实施。
验收：T0来源及只读编译；T1矩再现、原沿面迹一致、正确回复与位置/密度收敛、原nativeM/pressure/saddle及恒等式，扰动矩/交换测试B负控；T2/T3和真实C1/face mean/characteristic/AMR/RZ/曲面/原72静滴及完整波尚未实施，保持H13 FAIL。
失败规则：若条件参考失败先保存/定位，不切算法、调参数或隐式fallback；任何原method/API范围扩展先形成具体设计并由总纲§2.1决定。当前执行不增加第二production字段。
交付：推导、条件权重、原求解与独立audit、方法可行性/限制/来源及封存；所有工作、temp、输出存D:。

045122结果：36条件响应与36原M正交投影共72次public saddle自然退出0；真native的off-grid normal约三阶，响应仍约一阶，参考FAIL且phase跨越权重不连续，未采用。80位独立解析范数复算证明半格n64原对角动能漏约4.85%；同源准确normal/伴随力及原M的Rayleigh下界为一阶，明确仅限该假设族。已写独立证据和同一H13的phase-aware source/局部物理动能Hodge联合变更提案，涉及原M/能量冻结，需要方法授权；新提案至多24生产文件/10已有短契约，原全部物理门槛及历史FAIL保留。此前两项小维护待批可在批准后的联合预算内吸收，不额外开任务。所有模拟终态，无后台接受步；H13仍FAIL。

授权依赖复核：联合方法决定连续三轮仍未收到直接用户答复；具体设计、证据、原验收/可比性及封存已完成，现有授权不允许实施该M/动能改造。目标推进满足blocked审计条件，等待既有授权问题的明确答复；这不改写H13数学FAIL或整体目标，不启动H14/H15，不把自动续跑当作方法授权。审计记录在 `runs/checks/H13/20261009_045122/method_authorization_blocker_audit.json`，不改动已封存数值证据。
