# I12p3d / 纯路线阶段成本与规模证据

目标：补齐 Gate I 对实际独立 App 路径的阶段耗时、状态/几何数量和内存增长证据；本项只测量，不改变算法或调参。

允许修改：HermiteGalsFirstJetStep.h/.cpp（只增只读计时诊断）、InterfaceFirstJetBenchmark.cpp、InterfaceFirstJetOutput.h/.cpp（共5 production）；tests/interface_first_jet_output_contract.ps1（1 test）。任务/API/报告、runs 下可复现测量与审计脚本。不改现有求解顺序或数值表达式。

禁止修改：Core/几何/trace/transport/mixed/publication/transfer 算法、字典与物理输入、网格策略、容限、冻结实现/PDF、旧数值门槛；不混做性能优化，不缩短原完整窗口。

数学约束：计时使用 steady_clock，只在原阶段边界读钟，绝不参与接受、步长或输出频率选择。FirstJet 无 M、唯一 accepted geometry、实际特征/J 及 AMR exact/projection 语义保持；原逐步与最终数值及 identity 必须逐字保持。

实现：producer report 增加输入验证/边界数据/vertex transport/edge transport/center prediction/mixed closure/unique publication/total 的实际秒数；模块不写文件，App/IO 在接受后消费报告。performance.csv 保留原 App inclusive 阶段，并追加清楚命名的 producer 子阶段，不能把 inclusive total 与子阶段重复相加当作总成本。初始态没有 producer 步，只保留既有初始化/IO阶段。原 transported/prescribed/count报告不变。

T0/T1：39 targets，无新增 warning；原36 CTest完整数值输出保持，App契约核查每阶段有限非负、子阶段和不超过producer总计/所属App时间，且原输出频率独立、IO失败/原子重试保持。69冻结实现、11未改接线参考、旧reader分支、3PDF/p3c archive保持。

T2：再次运行原24实际完整组/2,072步/36regrid，与p3c及原模块数值/identity/输出诊断对照，不因测量改变结果。

规模研究：同一个最终主可执行，固定线程（OMP2、BLAS/MKL1），串行进程；n=8/16/32、uniform/AMR、translationCircleXY/translationSphereRZ/polynomialVortexXY，共18配置，每配置3次独立新目录、4实际接受步。平移dt=.005、vortex dt=.0005，保持相同物理终时与配置之间的参数；关闭空间帧但保留series/meta。报告初始化及每阶段每步的median/min/max、CPU/wall与采样内存，并保存采样频率/退出码/进程handle证据。内存采用同一子进程的Windows计数器，不把采样峰值谎称严格全程最大值。

动态规模补充：n=8/16、XY与RZ平移，各20步、dt=.005，在(4 8 12 16)交替网格；同一进程记录接受步结构规模，检查返回同一网格时primary/face/center数量不累积。既有完整400步vortex记录另用来检查长步数的表示数量保持。

判据：每叶center=1，独立primary标量数=3V+2E+C，按四边形拓扑不超过21*N；native faces不超过4*N，无第二份独立M/高阶全局物理状态。列出cut/graph随网格和时间的实际增长，不能用绘图点代替求解几何计数。测量CPU/内存对叶数的增长指数和每叶成本，不声称有限三档证明渐近复杂度；若两相邻档均呈稳定不低于二次的增长，或返回同网格时状态/几何数量持续累积，记FAIL并另立定位任务。明显噪声/不可用计数器必须明确，不能挑最快的一次代替三次记录。

失败规则：计时影响数值、减小原误差/物理窗口要求、阶段重复计费或伪造内存/规模数据、资源结构失控未解释即宣称效率通过，均FAIL。内部问题不标外部BLOCKED。

交付：被动阶段诊断、完整回归、三档多次测量与动态结构证据、局限和瓶颈说明、API/报告/归档；之后回到I12p3复核完整Gate I。

执行结果：**PASS**，见[报告](I12p3d_FirstJet_Performance_Scale_Report.md)。37 CTest、24 完整组/2,072 步和 58 次正式规模运行通过，完整数值与 VTK 内容保持。Gate I 待后续独立复核。
