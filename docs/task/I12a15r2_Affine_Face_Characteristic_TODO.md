# I12a15r2 / 面—时间仿射 RK2 包围

目标：只解决原实际完整单步 RK2 departure 在 face×time 上的相关性包围，返回中心+切向/时间线性项+严格残差，供下一项计算分块内的相区域。原 a15r1 的矩形 [0,1] 含量积分未在 4096 panels 下通过 shear；本项不将该 FAIL 改名为已解决。

允许修改：`PolynomialVelocitySlab.h/.cpp`（2 production）；新 `tests/interface_affine_characteristic_contract.cpp` 与 `tests/CMakeLists.txt`（2 test-related）；TODO/设计/API/报告/入口/证据。继续使用既有 strict/no-LTO Core TU，不改 root CMake。

禁止修改：原 velocity 数据及实际求值、原 RK2 表达/CFL/axis/boundary/substep，a13/a15 既有范围入口，Q4 相与库存、M 装配、几何发布/reflux、现有测试参数/容差与全部 gate。

数学约束：单位参数 q_s,q_t∈[-1,1] 用真实 double 端点作为实数定义原 face/time 的精确仿射坐标。每个 departure 分量必须 pointwise 属于 c + a_s*q_s + a_t*q_t + R，其中 c/a_s/a_t 为返回的有限 double，R 为 strict outward 余项。必须覆盖原有限 round-to-nearest/gradual-underflow 运算、endpoint/half-time 混合分支，以及极小 dt 的 mean/arrival 分支差异；不能把浮点 map 当光滑函数后漏掉舍入。原 a15 tight-time branch witness 在入口完整验证同一区域，无新分支豁免。

实现要求：先写设计并独立推导误差。使用带两个一次参数项的 affine arithmetic；常量/一次项用 outward 区间，所有非线性乘积进 residual，逐原运算加入显式 RN 误差。复用同一 polynomial component/Horner 与 rk2Midpoint/rk2Target 表达。最后把系数选择为区间内代表值，系数余宽移入 R。来源/面/时间元数据保留，物理 point 查询用 outward 参数反算，结果与既有全 box witness 相交。零/恒定量保留代数精确简化，任何非有限状态和非法输入显式失败。

验收：T0 无新增 warning、<=2 production+2 test、全 runtime 重建、81 frozen/3 PDF/r1 archive不变；原23相关测试保持22 PASS+1相同shear FAIL，全部原输出不变。本项新独立 contract T1：XY/RZ uniform/AMR 及 native face/轴、常速/反号/旋转/shear/cubic，多时空点与原 actual tracer B 对比；同源真实 polynomial/RK2 解析参考和独立高精度/有界误差参考；区间细化残差按预期收缩，单个 affine form 不能以全 departure box 伪装；非单位basis、非dyadic endpoints、subnormal速度/dt与half-time原分支；moved/无关mesh由既有Interface守卫，本Core验证time/face/domain/FP/nonfinite/branch拒绝。不得以采样代替 pointwise 包含性的代数证明。

失败规则：实际或实数代数参考逃出范围、未认证branch返回、残差漏舍入、旧输出改变，均 FAIL并保留。仅真实外部阻塞才BLOCKED。工作中必要修正限定本数学事实和文件预算，不扩展到相积分。

交付：源代码、数学证明、完整 T0/T1、API/report与 archive。之后另立任务将 Q4 相沿 affine departure 拉回，并认证已知液/气子区域；原 shear 和 I12全矩阵继续保持未通过直到实际验收。
执行状态：**PASS（本项）**。新 affine characteristic contract 的 63 区域、9,129 actual traces、514 identity、18 guards 通过，原 23 项结果与完整输出不变。相积分未消费本原语，原 uniform shear、I12/Gate I 仍 FAIL。详见 [报告](I12a15r2_Affine_Face_Characteristic_Report.md)。