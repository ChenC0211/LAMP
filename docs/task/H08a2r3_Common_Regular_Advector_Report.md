# H08a2r3 / 共用正则advector可行性

2026-10-05，**PASS（构造可行性，尚未部署）**。[任务书](H08a2r3_Common_Regular_Advector_TODO.md)、[实施前设计及补充门槛](H08a2r3_Common_Regular_Advector_Design_Change.md)。完成原144组结构检查、32族96窗口及新增sin非多项式收敛；生产速度仍为H08a2r1的原分片B，FirstJet消费桥接、H08a时空材料积分和H11均未完成。

## 验证的是哪一个field

临时C++导出模式从实际CCH布局输出60份native cell/face/dual段、实际主FirstJet DOF与Core Q4约束probe。NumPy原型仅用这些几何和native法向均值：局部缩放P2最小范数重构共享jet → 同一Q4 tensor/master约束 → 每个实际face的紧支撑均值校正。所有位置/J查询及面积分都用这一个候选field；没有用解析参考修正接受值，也没有给J单独恢复另一个速度导数。

均值校正为δu_f=(15/8)(u_f−mean_f(u_0·n)) b(xi)g(s)n，g=16s²(1−s)²，b=(1−xi²)²；系数由native运动学均值残差唯一决定。g/b在支持端的值和一阶导数为零，b在面上的法向导数为零；两侧支持宽度不同仍保持C1。XY与真实RZ测度都有mean(g)=8/15。此项不携带材料库存，也不是界面内部bubble或体积拟合。

新的field允许点态散度随位置变化，只保持cell积分散度等于原native面通量和。全局C1与任意逐cell不同的点态常散度不能同时强制；本项没有改变D_I及C=−D_I^T。四项材料硬约束及真实R/uR保持，未引入pressure→dual精确commuting、quarter分配或R=0。

该原型刻意更改原B的空间重构，不能与H08a2r1视作同一来源。它不发布C++ advector identity、FirstJet candidate、phase/material F或新的接受状态；完整clock/typed source接线仍须在生产部署时验收。临时Python不是另一个生产求解器。

## T1：C1、native面矩及边界

48个XY/RZ × uniform/central-box真实2:1 × n4/8/16 × 全域/半域/内孔/分离域布局，各检验无指定wall、全外部wall、单个实际外部wall，共**144组**。保留全部零场、P2基、随机native输入、缺秩和边界记录。

| 独立证据 | 结果 |
| --- | ---: |
| 60份Q4局部25数据与Core原约束probe | 逐值相同，最大差0 |
| 同向native面校正的交叉支持检查 | 2,220,298对，无非对角正支持 |
| 最终field的内部面检查（含三种law） | 23,898面 |
| 独立de Casteljau/Bernstein C1缺陷 | 1.14e−13，门槛1e−11 |
| native面均值：实际query / 独立Bernstein积分 | 6.78e−16 / 6.01e−16 |
| 1,335,600个足阶体积求积点的Gauss缺陷 | 8.52e−16 |
| 零输入、指定wall法向值、轴u_r及轴奇偶所需一阶条件 | 精确零 |
| 可辨识P2的完整field值/导数界，含全部面校正 | ≤8.86e−12 |

C1独立检查直接从保存的Hermite数据转为Bernstein controls并限制到实际subface；未使用“几个样点相同”代替系数检查。体积分在每个实际cell及校正支持分界上分片，RZ积分2π[r(∂_r u_r+∂_z u_z)+u_r]，不在轴上除以r。随机输入的点态散度确实不是cell常数，积分仍满足原D_I。

48个无wall布局中**38个支持可辨识完整P2**。其余10个极窄/分离布局保留统一最小范数定义，不能标成P2精确。跨三种law共有**3,841条缺秩stencil证书**；独立用闭式真测度面矩核验非零null多项式，最大scaled残差1.17e−14。它们仍通过全部C1/面矩/边界检查，且未替代原收敛矩阵；T2全部支持可辨识。这里不恢复“所有不足信息的布局都必须P2精确”的过约束。

本批60布局均满足Core现行活动粗owner条件；因此与Core的逐值比较全部执行。设计中对活动侧独立定义约束的更一般selected情形，**不能由这60例自动声称已验证**；生产部署须保留明确支持/拒绝范围。

## T2：同一个运动及其梯度

沿用H08a2r2的全部**32族、96窗口、56,992条真实主vertex/edge/center记录**。XY u=(y²,0)，RZ u=(0,r⁴)，uniform/横向2:1条带、n8/16/32、正负方向、静态/0.7→1.3线性时变；完整slab为0→0.25，目标时刻0.125/0.25。native均值用闭式积分。

位置/J来自同一固定步RK4阶段及其链式微分，64/128/256子步全部运行。最细两档时间差与空间误差之比最大2.97e−6，低于10%门槛；XY舍入误差单列，不对舍入噪声计算阶。原CCH事件tracer尚未切换，这里的RK4仅用于证明新field可形成可微的同源运动。

| 布局 | source值阶 | source导数阶 | 有限时间phi阶 | 有限时间grad阶 |
| --- | ---: | ---: | ---: | ---: |
| XY uniform / 2:1 | 舍入精度再现 | 舍入精度再现 | 舍入精度再现 | 舍入精度再现 |
| RZ uniform | 4.0312 | 1.9985 | 4.0780–4.0844 | 2.0066–2.0115 |
| RZ 2:1 | 3.5539 | 2.0815 | 3.5380–3.5422 | 2.1010–2.1097 |

原门槛source值/导数为2.5/1.8，运动phi/grad为1.8/1.8。对**288个实际骨架点**另做普通双侧差分，最大scaled差1.31e−7，小于5e−7；参考解与差分只用于检查。原始point/J、native数据、时间细化和每条DOF均保存，独立分析器重算RMS与全部观察阶。

## 原反例及非多项式补充

XY uniform n8顶点(0.5,0.5)、T=0.25：新y²来源的departure x=0.4375，J_xy=−0.25；四个incident cell给相同value/J。原sin面均值来源改用新重构后，departure x=0.494869610646，J_xy=−0.007760939391；左右实际映射的差分随扰动缩小趋于同一个J，1e−7扰动时最坏差3.88e−7。没有平均原B的两个单侧J。

**粗网格sin的物理梯度误差仍明显。** 唯一导数本身不证明高精度，因此在主矩阵之后、运行前固定了六个额外空间窗口：XY uniform/条带AMR、n8/16/32全部原DOF、正向静态T=0.25，共**3,562条DOF**。native均值来自cos差闭式积分，沿用原精度门槛；128/256时间差也检查。

| sin布局 | source值阶 | source导数阶 | phi阶 | grad阶 |
| --- | ---: | ---: | ---: | ---: |
| uniform | 3.8611 | 1.8354 | 3.9150 | 1.8365 |
| 2:1 | 3.6392 | 2.0759 | 3.6901 | 2.0586 |

补充检查通过，但uniform末档尚不能报成精确二阶，更不表示任意欠分辨率速度都准确。后续完整H10b-v2仍不可省略。

## 成本、退役与生产部署边界

T1支持最多24个native均值，系数范数最大26.76；原T2最细RZ/2:1准备有968,892个存储权重、单datum最多22个native均值，Python准备约1.67秒（n16约0.43秒）。存储/准备随实际网格规模增长，最细组完整完成；这些是原型成本，**不是C++动态耦合效率已通过**。生产须在几何准备时完成SVD/约束选择，evaluate不重复求解。

首轮pilot把两份已舍入单项式转换相减并用系数绝对和估计P2误差，界过于宽松，四组报警；该版本和日志保留在`pilot_certificate`。正式检查直接对同一Hermite数据差作Bernstein界，独立检查另用de Casteljau限制；原1e−11门槛未改变，最终field的全部均值校正也计入独立P2误差界。

临时C++导出入口已退役，**1379份启动源码/输入哈希全部恢复，生产和长期测试改动均为0**。恢复后标准Release及架构/目录/momentum契约**4/4通过**，原momentum的**2,024行确定性输出**与H08a2r1相同（只排除原两项计时列）；其他五个数值exe哈希保持。文档定稿后另跑三项架构/目录检查。未重跑H10o/H10b-v2，因为本轮未更换生产B；部署新B时必须重跑。

已核对的消费者如下；不能只改局部velocity query留下旧atom。

| 现行位置 | 部署时必须处理的假设 |
| --- | --- |
| `CchSpatialAdvector.h/.cpp` | 三个法向系数、二次value/integral、旧trace/lift/piece；改为同一正则field与实际支持分片 |
| `CchPhaseGeometry.cpp::ratePolynomial` | 固定k<3；新quartic normal与RZ测度最高五次必须完整传递 |
| `UnitPolynomialRoots.h` | Core已支持≤11次，不需要增加另一套相根/积分引擎 |
| `CchCharacteristic.cpp` | 旧volume-coordinate二次C1证书、初始/同时事件拒绝；须适配新field的真实次数及等价分支 |
| 原phase/material/dual/momentum契约与benchmark导出 | 多处固定三系数参考/列格式；应按依赖和文件预算迁移 |
| `cellDivergence`消费者 | 旧点态常散度语义变为积分平均；点态散度不可伪称相同 |

下一唯一任务[H08a2r4](H08a2r4_Regular_Advector_Preparation_TODO.md)将已验证的重构移入同一Hydro owner的几何准备层。随后按消费者预算升级atom次数、切换唯一生产B及轨迹证书，再重验最终native算子/原affine/H10b-v2并回到H08a3。旧任务报告不改判；H08a整体仍未通过。

检查：`runs/checks/H08a2r3/20261005_124715_641`；冻结：`runs/baselines/H08a2r3_20261005`。`audit_verification.json`含独立证据，`t1/raw_*.npz`及`t2/raw_*.npz`保存实际候选数据/轨迹，`sine.json`记录补充观察阶，`verification.json`记录恢复。导出exe SHA256：`96ABB7FF537850EF232AC012F45A03E83F3382E4A85F4739BD4A466B69C8290D`。
