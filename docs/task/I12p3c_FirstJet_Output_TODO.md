# I12p3c / 纯路线正式输出与诊断

目标：为已验收的独立 FirstJet App 接通规范目录的只读 VTK XML、accepted-frame、标量时序和配置/运行元数据；所有几何来自当前唯一 publication，输出频率不改变数值路径。

允许修改：新增 Apps/Interface/IO/InterfaceFirstJetOutput.h/.cpp；修改 InterfaceFirstJetBenchmark.h/.cpp、InterfaceBenchmarkInput.cpp（共5 production）；修改 tests/interface_first_jet_app_contract.cmake，新增 tests/interface_first_jet_output_contract.ps1（共2 test相关）。任务/API/报告与 runs 证据。只在既有 App CTest 内调用新输出检查，不增第三个测试配置文件。

禁止修改：Core/Modules 数值算法、mixed/trace/publication/transfer、旧方法 reader 主体、冻结实现/输入/PDF、数值阈值、全局构建选项；不借输出重建第二份曲线或维护独立相库存。

数学/所有权：VTU 使用真实叶网格、同源几何体积分数/中心标量；VTP 点从 publication.sample 查询零面位置、法向、组合曲率，采样连线仅为显示，不进入求解。XY/RZ 标明测度。只有初始接受态及实际 accepted step 能进 PVD/frames。失败记录 attempted/accepted 时钟、原因与阶段；数值失败不得提交新状态，IO失败不得吞掉或写成成功。

实现要求：沿用 method hermiteGals 和既有 writeFields；false 只关闭空间帧，series/meta 仍写。true 开启字段，新增可选 writeEverySteps（默认1；false时显式提供此键应拒绝）。初始/最后接受帧始终输出。--check-case 不建网格/目录。正式运行要求全新 output，禁止覆盖已有/符号链接目录。输出规范 fields/volume/*.vtu、fields/interface/*.vtp、fields/solution.pvd、series/{conservation,solver,performance,events}.csv、meta/{frames,effective_config,configuration_trace,manifest,build_info,failure}；尚未求解的库存/质量/通量/能量/求解器量标 NA 或明确不适用，不填零冒充有效数据。无需新增空 energy 文件。

元数据：保存完整生效配置和已有 ConfigurationTrace 原始/默认/单位/consumer；case路径、独立run_id、UTC开始时间、实际编译器/Release标识、无Git时unavailable、输入哈希、初始/最终geometry identity与网格规模、接受步/终时/退出状态。哈希须由独立工具校验。App 层实际 step/mesh-transfer/IO 耗时写入performance；更细核内阶段和规模研究归后续p3d，计时不参与接受。

T0：39 targets，无新增warning；原36 CTest完整数值输出保持，新App契约因输出能力更新按明确预期改变；69冻结实现、11未改接线参考、3PDF/p3b归档保持，reader旧分支仍单独核验。

T1：实际主程序六类问题、uniform/AMR/重网格/非多项式、XY/RZ。解析所有VTK/PVD/CSV/JSON，核对网格拓扑/数据长度/有限数值/正向法线/曲率参考/时间/identity/volume与stdout同源。字段关闭、每步、稀疏三种输出必须保持每步数值及最终误差；输入追踪与哈希、fresh-output拒绝、原数值CFL拒绝、接受后IO失败及失败manifest通过。用本机现有ParaView读取代表性PVD/VTU/VTP（若无可用工具，明确记录缺口，不冒充可读性验证）。

T2：与p3b同输入的24实际完整组/2,072步/36regrid；正式series开启，VTK按稀疏频率输出，核对所有旧数值token不变和最终帧/终时。小规模输出频率对照只做必要测试，不替代完整原数值窗口。

失败规则：输出派生量反哺数值、重拟合法向曲率、rejected出现在正常帧、输出频率改变结果、IO失败被成功退出掩盖、虚构未有物理量或配置consumer、缩短原完整窗口均FAIL。内部实现问题不标外部BLOCKED。

交付：只读App/IO、可追溯配置/时序/VTK、输入/运行证据、独立输出审计和归档；之后执行I12p3d，完整Gate I待复核。

结果：PASS，见[报告](I12p3c_FirstJet_Output_Report.md)。37CTest、18频率对照、真实数值/IO失败与有界重试、ParaView实读及24完整组/2,072步通过；116正式帧，数值与几何identity保持。完整Gate I仍待I12p3d。
