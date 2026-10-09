# I12p3d / 纯路线阶段成本与规模报告

**PASS（本项）。** 实际 FirstJet App 的被动阶段计时、三档网格重复测量和动态结构检查完成。完整 Gate I 随后独立复核，本项不宣告 Hydro 或低黏耦合通过。范围见 [TODO](I12p3d_FirstJet_Performance_Scale_TODO.md)，计时包含关系及重网格源/目标规模见 [API](../api/I12p3d_FirstJet_Performance_API.md)。

修改 5 个 production 文件和 1 个既有测试，无新增 production/test 文件。producer 只增加 steady_clock 阶段边界和只读数量记录，App/IO 在接受后消费；未改 Core、字典、数学表达式、容限或步长。去掉新增诊断行后，producer 主体与 p3a 归档逐字相同。七个子阶段属于 producer total，后者属于 App interface_step，禁止重复相加。

| 验证 | 结果 |
| --- | --- |
| 39 targets 最终 Release 构建 | PASS，43.673 s，无新增 warning |
| 原 36 CTest | PASS，373.357 s；35 个完整数值输出逐字保持 |
| 最终 App 契约 | PASS，62.689 s；阶段非负、包含关系、源/目标数量、频率独立及真实 IO 失败/原子重试通过 |
| 原 24 实际完整组 | PASS，2,072 步/36 regrid，470.495 s；逐步与最终数值、geometry identity 与 p3c 相同 |
| 完整输出对照 | 116 帧/232 VTK、2,096 初始及接受行；全部 VTK、PVD、frames、conservation/solver/events、accuracy、有效配置与追踪字节保持 |
| 固定规模 | 18 配置 × 3 独立重复，216 步；同配置数值及 identity 完全相同 |
| 动态规模 | 4 组/80 步/16 regrid；返回同网格时 primary/face/center 数量不累积 |
| 完整窗口结构 | 24 组全部满足 native 状态界，含四个完整 400 步 vortex |

固定规模采用 n=8/16/32、uniform/AMR、XY 平移/RZ 平移/XY 多项式涡旋。每次 4 步，平移 dt=.005，涡旋 dt=.0005；每类问题跨网格保持相同物理终时。动态组 n=8/16、XY/RZ 平移各 20 步，dt=.005，在 4/8/12/16 步交替网格。规模实验关闭空间场输出，保留 series/meta，使用最终同一可执行、串行进程、OMP2/BLAS1/MKL1。58 次正式运行 wall time 合计 119.958 s；独立 pilot 另存，没有替代正式重复。

以下是 n=32 的每次四个物理步平均耗时，再取三次独立运行的中位数；初始化另列在原始汇总。private 为每次运行采样峰值的三次中位数。

| 问题 | 网格 | 叶数 N | primary 标量数 | 最大 graph 数 | App step 秒 | publication 秒 | private MiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| XY 平移 | uniform | 1024 | 8515 | 48 | .469391 | .428516 | 67.80 |
| XY 平移 | AMR | 1792 | 14499 | 93 | .846104 | .774303 | 96.47 |
| RZ 平移 | uniform | 1024 | 8515 | 35 | .363210 | .301892 | 67.90 |
| RZ 平移 | AMR | 1792 | 14499 | 44 | .507890 | .400883 | 97.08 |
| XY 涡旋 | uniform | 1024 | 8515 | 34 | .468748 | .315082 | 68.07 |
| XY 涡旋 | AMR | 1792 | 14499 | 58 | .785451 | .521344 | 96.62 |

每叶一个中心，独立 primary 标量数 3V+2E+C ≤ 21N，native faces ≤ 4N。没有新增独立 M 或全局二阶物理状态。完整 fixed vortex 的最大 graph 数为 14/19，dynamic 两组均为 19；几何表示没有随 400 步持续累积。该结论针对实际已测窗口，不能推及任意拓扑复杂度。

相邻档增长指数按 log(成本比)/log(叶数比) 计算，保留三次 min/median/max。六种问题/网格组合的 App step 中位增长指数为 .036–.594，producer total 为 .035–.593，采样 private 为 .194–.544，working set 为 .150–.484。预设“两个相邻区间均稳定不低于二次”的标记为空。小网格固定开销明显，三档实验不证明渐近复杂度，也不保证其它硬件吞吐。n=32 的 publication 占 App step 约 66.4%–91.5%，是已测主要成本；本项没有实施性能优化。

内存由同一 Windows 子进程计数器采样，保存 WorkingSet64、PrivateMemorySize64、PeakWorkingSet64、CPU/wall 和非空退出码。请求间隔 20 ms，3,315 个实际间隔为 min 25.957 / median 32.022 / max 51.039 ms；报告采样峰值，不冒充严格全程最大值。环境为 Windows NT 10.0.26200、PowerShell 5.1、32 logical processors、Release/local/jobs1，完整环境见 measurement_environment.json。

架构为 12 targets、388 TU、918 文件、3149 内部 includes、5 原 legacy edges。69 冻结实现、11 未改接线参考、3 PDF 和 p3c archive 保持；第 12 个 reader 参考是 p3b 已允许的扩展，本项未改 reader，旧主体继续逐字验证。Interface 库重建导致 39 运行文件中 33 个 SHA 改变、6 个保持；不能把本项误写成仅 App 二进制改变。最终 QUANTUM.exe SHA-256：`102AC6E612D62B7A3E81E881080BB28B4516B1813DCCAC18AF87D82571A01FA8`。

证据：[I12p3d_20261002](../../runs/baselines/I12p3d_20261002/)。verification/source_audit/producer_numeric_identity/regression_comparison/output_verification/full_structure/scale_verification 给出机器判据；scale_details、scale_summary、scale_growth 和 dynamic_scale_summary 保留全部重复和原始量。所有输入、输出、内存时序、初次及最终契约、脚本和运行文件归档，最终 SHA 见 archive_verification.json。先行计时版本通过后补齐了 regrid 源/目标规模区分，最终构建、契约及全部运行均使用修订后版本。

图表：[规模与资源](../../runs/baselines/I12p3d_20261002/scale_overview.png) / [PDF](../../runs/baselines/I12p3d_20261002/scale_overview.pdf)，[七阶段成本](../../runs/baselines/I12p3d_20261002/stage_costs.png) / [PDF](../../runs/baselines/I12p3d_20261002/stage_costs.pdf)。两张图已实际打开检查，图中误差范围和采样含义与原始数据一致。下一项对当前完整纯路线执行只读 Gate I 复核。
