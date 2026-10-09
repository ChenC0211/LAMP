# H10h / 最终native算子的PDE矩联合可行性

状态：**FAIL，诊断闭合**，2026-10-04；固定最终矩区间在 XY 2:1、rho=1e3/1e6、n8/n32 四组不可行，见[报告](H10h_Final_Operator_Joint_Feasibility_Report.md)。只改研究观察与文档，无生产 QP/新 API；工程 4/4 通过。H10 未关闭，用户本轮解除 H08/H09 附加过约束，下一单任务[H10i](H10i_Spatial_Material_Flux_Transfer_Design_TODO.md)。H11/H08a 等待；下文保留启动要求，历史副本和设计/失败证据归 H10h baseline。

目标：直接以最终共享one-ring skew pair为设计变量，联合constant/affine空间矩、真实材料点采样PDE误差、material-KE/skew及高阶/行范数预算。均匀材料中心K已存在原有界可行证书，内部affine须先精确/二阶通过；非均匀材料两条非零缺口不允许全部行精确affine，若采用二阶余项必须在最终算子里联合求出并证明实际尺度，不能先任选一个修正target再强制匹配它。

启动：独立校验H10g冻结链和protected哈希，先写本任务数学设计、固定预算与失败判据，再源码。允许CchKepConvection.h/.cpp、原短契约、现有研究driver/README和对应文档，2 production+1 test+1 research；临时runs/checks/H10h/时间戳、关闭runs/baselines。禁止改Interface、H08/H09的M/F/R定义、原App/输入/PDF、CMake/通用工具；无库存/bubble、AMR特例、人工黏度/donor fallback、解析状态或载荷修正，不删边界/轴行。若证明必须一般Hodge/材料矩重设计，拆出明确数学任务，不偷偷改变KE或提前H08a。

实施要求：

- 先区分均匀内部affine硬矩、显式闭合端口和真实分相材料PDE余项。显式端口需实际拓扑/alias/clock归属；周期方向使用wrapped空间增量及winding，不把横向rank-2条件错误推广。
- 最终K^T=-K、K*1=q/2保持唯一F/R/uR；原opaque API和124场景先保持。真实actual advector只能来自原H07率，不能由参考数据提供。
- 原primal F已产生的全局一次/二次缺口作为约束数据进入最终可行性判断。以最终pair的PDE一次余项为联合未知，均匀原硬矩保持；若分相需二阶区间/残差预算，从Taylor展开、真实M和同源F误差预先推导。**不得仅提高H10f的8、放宽原空间验收、故意增加粗网格误差制造收敛阶，或把下界比值当作普适不可能证明。**
- 用已知必要条件预检任意固定逐行余项区间：对任意lambda，缺口投影必须不超过该区间对应的支持函数。不可行时导出可独立复算的dual证书；可行只代表全局必要条件，仍检查真实one-ring局部rank、l1、P2/P3和实际apply。
- 一般中心K作为最小改动锚点；正确约束工作集/QR、可行系数、非负乘子和KKT均需实际证据。禁止退役Dykstra/罚函数fallback、单segment trace补丁或仅增迭代预算。若适合只约束必要二阶预算，不强制已证明不相容的完整P2等式。
- 最终每个native行系数/M、实际二次/三次余项及全部PDE误差贡献须独立复算。保持研究driver单一最终operator调用，不内藏第二套算法。

验收顺序：T0标准Release/架构/目录/protected哈希；T1原1e−12 KE/载荷、1e−10常速度、真实非零R/uR/跨界/regrid，加新显式端口/材料矩/拒绝契约；T2首先原内部affine反例n8/16/32四组合，之后全native行与48组真实材料、系数和PDE余项尺度，原空间阶>=1.7、最细L2<=0.05、Linf下降；T3前置全部通过才原H10b28窗口/15,116步/84gate，原KE、时间和四个donor控制全部保留。不得把观察器exit 0当作PDE通过。

裁决：任何目标/全局或局部可行性/实际系数/余项/恒等式/完整研究失败均FAIL并退役数值候选；保留可复算证据，不标外部BLOCKED，不进入H11。交付最终联合可行性设计、实际唯一operator/API/独立契约（若通过）、全部原反例/材料/行/完整研究或明确未运行项，以及平面manifest/哈希验证。本任务诊断已闭合；新空间材料设计尚未启动。
