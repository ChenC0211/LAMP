# I12p3b / 纯路线正式独立 App 与字典入口

目标：让主程序从严格字典配置实际运行当前 FirstJet 模块时间步和动态网格，独立于测试可执行；不重写数值内核。

允许修改：新增 Apps/Interface/Verification/InterfaceFirstJetBenchmark.h/.cpp；既有 InterfaceGeometryBenchmark.h、InterfaceBenchmarkInput.cpp、InterfaceGeometryBenchmark.cpp（5 production）。新增 tests/interface_first_jet_app_contract.cmake、修改 tests/CMakeLists.txt（2 tests）。任务/API/报告及 runs 证据。

禁止修改：Core/Modules、旧 benchmark 算式/参数、冻结 case/资料、Hydro/Electricity、全局 CMake 编译选项、数值容限；不在 App 复制 tracer/J、pullback、mixed closure、publication 或 mesh-transfer 算法。

数学约束：显式 method hermiteGals 表示 FirstJet 非守恒；无 M/体积闭合/相库存 reflux。只初始化一次解析 first jet；后续调用 prepareHermiteGalsFirstJetStep 和 transferHermiteGalsFirstJetMesh。外边界解析数据只经已验收的 boundary query，内部不重置。App 是唯一 accepted clock/mesh owner；失败不发布新状态，网格和 publication 成对替换。所有诊断消费同一 publication。

实现要求：沿用 system/run.dict 的 Run.app=Interface 和 interface.dict 的 Interface.Benchmark，新增独立 typed config 和 ConfigurationTrace。纯分支严格拒绝 legacy/no-op 选项，保留旧 reader/dispatcher 语义；--check-case 不生成网格或输出。命名六个原基准物理问题：staticCircleXY、staticSphereRZ、translationCircleXY、translationSphereRZ、solidRotationXY、polynomialVortexXY；原单位域/形状/规定速度显式报告。配置 resolution、dt、steps、initialMesh uniform/amr、nonpolynomial、可选 Regrid.atSteps（严格递增、位于完成前，交替 exact refinement/显式 projection）。不假装这是任意几何/速度或 Hydro 耦合入口。

初值 master/hanging 归属调用现有 Core skeleton，不复制 fine-side 筛选。初始 primary/center 由规定解析场采样，再调用 I06/纯 publication，不建立第二个几何 owner。时间推进用原 constant/owned polynomial velocity provider，沿用原浮点时钟加法顺序。stdout 输出生效配置、接受步体积/漂移/曲率和重网格信息；正式文件输出属 I12p3c，本项必须显式 writeFields=false，true 暂时报“不支持”，不可忽略。

验收：T0 39 targets无新warning，原36CTest/34数值及producer完整输出保持（仅架构预期新增项变化），69个S08冻结实现/3PDF/p3a archive保持；12个接线参考中允许本项reader新增独立分支，原分支需单独核验（S08原报告不冻结reader其它分支）。T1 主可执行的 --check-case/实际运行/拒绝输入：六类问题、XY/RZ、两种网格、非多项式、真实regrid；未知/重复key、非法类型/数值/时钟/网格/步序、未接通输出、旧无效方法选项必须失败；运行时CFL拒绝不能产生虚假accepted步。保留代表旧App入口回归。T2 与原12组fixed、12组dynamic同resolution=8、dt、步数、速度及终时，通过主可执行运行全部24组/2,072步；核对逐步体积漂移、最终位置/gradient/curvature与原档案，动态36次切换，不降低原阈值。检查配置值确有consumer，不复跑无关模块。

失败规则：解析内场重置、隐式legacy回退、未消费配置、错误accepted时间/网格配对、过早声称VTK/production/Hydro已接通，或原数值门槛/基准范围被缩小均FAIL。内部实现问题不标外部BLOCKED。

交付：正式App入口、严格typed字典和配置追踪、可复现输入/日志/原结果对照、API说明/审计/归档；之后进入I12p3c，完整Gate I仍待输出与效率证据。

结果：PASS，见[报告](I12p3b_FirstJet_App_Entry_Report.md)。37项CTest、实际主程序24完整组/2,072步/36regrid通过；逐步及最终数值token与原模块档案相同。初次失败及输入全部归档；完整Gate I仍待输出/效率证据。
