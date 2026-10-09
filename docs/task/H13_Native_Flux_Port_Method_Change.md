# H13 / native面通量矩一致的运动—弱牵引联合改造提案

状态：2026-10-08，**用户已明确回复“同意”，授权在同一H13内实施本提案；候选尚未采用/验收**。H13仍FAIL，H14未启动。

授权登记：2026-10-08 09:36（Asia/Shanghai检查点）。用户回复对应本提案的共同source面通量矩端口与完整P2弱牵引联合改造；不授权放宽门槛、改M/C/F/R/uR/KEP或自动fallback。此前“待授权”状态是历史检查点；以下实现进入候选验证，仍不得宣称H13 PASS。

本提案只属于同一H13任务。冻结的r11/common source准备和H13载荷都没有改动；以下面端口是隔离的假设算子，不是实际发布速度场或FirstJet推进。依据与程序在 `runs/checks/H13/20261008_084706`。

## 1. 为什么提出这项变更

原准确native压力/加速度近二阶，但当前共同source把切向压力梯度与横向加速度变化带入界面normal trace：normal值误差约一阶，切向导数相对误差在n16/32/64约1.267/1.204/1.179，不趋零。原C1和face平均正确不等于界面运动梯度正确。这解释了为何正确压力矩、局部消差和频率通过仍不能保证原表面能。

本轮用独立、完全由面平均矩推导的端口作两条对照。保留原运动而只换弱牵引，原完整初值频率末档2.084阶，却仍有5.56倍初始能量的峰值缺陷。将假设运动端口与同一弱牵引匹配，在三密度比上频率末档1.987阶、原native动能加表面能的峰值缺陷最大5.96e-12。未改变原M/C、初值、频率定义、动能或表面能，全部几何零分支也保留。这不是所有原方法都不可能的证明，只证明这套已复算方案若保留当前运动端口，不能通过原能量验收。

三个独立C++运行以原M/C/鞍点和同源制造应力重新装配牵引并求初始加速度，最大相对差1.04e-13。条件端口normal与切向导数在n64误差0.000819/0.000830；它们明确标为条件公式。独立矩阵指数重放18组全周期参考，原能量差最大5.75e-11。没有实际新毛细周期或general-curve/AMR证明。

## 2. 具体端口，不从频率反算系数

均匀面长h，U_i为原native normal face平均。共享端点的值和切向导数由有限体积P3矩再现确定：

`v_i=7/12*(U_(i-1)+U_i)-1/12*(U_(i-2)+U_(i+1))`，

`d_i=[15/12*(U_i-U_(i-1))-1/12*(U_(i+1)-U_(i-2))]/h`。

每个面的Q4 normal trace由 `v_i,v_(i+1),h*d_i,h*d_(i+1)` 和其准确native平均U_i唯一确定。保留共享C1端点与原face通量，不增加全球velocity DOF。原Q4边基和零端点face平均bubble已经能表达此端口；原cell内部和横向/mixed jet仍需完整核验。

平面参考用同一trace的弱牵引 `F=-integral J*B_face^T ds`，加原常应力 `J0*C*alpha` 反力通道，J来自原同源完整P2曲率，不作曲率投影/滤波。原压力矩、切向及端点jump不能丢掉。参考constant balance作用缺陷约2.8e-17；生产必须用锚点差实现完整 `J_f(J0)=J0*I`，仅C*alpha上的常量作用还不够。

参考的两端偶对称闭合对应这次闭域波研究的已声明对称条件；不得因此给所有impermeable边界暗加free-slip。一般闭域端点必须使用当前typed边界条件和一致单侧矩闭合，不能用case名称选公式。

## 3. 请求授权的实现范围

请求允许在**同一H13内**联合修改：

1. 共同CCH source的master顶点normal-component值及切向导数准备，从混合横向P2观察改为共面的原native flux矩。均匀内部退化为上式；非均匀/AMR/RZ用真实face区间与物理测度的同一有限体积矩构造。现有共享Q4、master/hanging约束和原face平均修正保持；不增加global DOF或每步拟合。
2. 该共同source对应的局部、稀疏弱牵引/完整P2应力矩装配。原常应力完整矩阵、原C负转置、材料M和原native动能保留；保留一般曲面的几何运动/体积功余项，不强制精确功或R=0。

预计生产范围不超过四文件：`CchRegularAdvectorGeometry.h/.cpp`（prepared顶点权重及同一局部native basis/jet只读查询）、`CchInterfaceNormalLoad.h/.cpp`（同源完整应力弱装配）。测试仅扩展现有regular/traction两份短契约；完整研究归benchmarks。必要跨API预算须先登记，不能边做边扩scope。

这不是开启独立r11/source精修任务。它改变此前关闭的共同source准备规则，并与H13能量/精度同时验证。不得修改characteristic积分/Jacobian公式、容差、材料F/R/uR、KEP或accepted FirstJet状态；新速度与Jacobian仍来自同一个common source，不能留两条production主线或失败后fallback。

## 4. 新的算子与完整验收

- T0：标准构建无新warning，architecture/repository_layout和来源身份。
- T1：共面P2/P3矩、共享C1/value/derivative、全部真实face平均、RZ真测度/轴奇偶、2:1/master限制；同一query/normalPolynomial/actual characteristic map与FD-Jacobian；完整常应力矩阵、原压力负转置、局部弱力与实际FirstJet方向功，曲率jump保留。
- T2/T3：原r11/H08a全部受影响短契约与28原窗口/空间/时间/R门槛重跑，保持原断言；原H10/H11对流/耗散规定范围回归。它们的历史PASS不预先代表新source通过。
- H13原全范围：原静滴1728解/七组伪流例外、三个空间/时间层级、完整实际波周期、strict mu=0及低黏、high-density/AMR、原频率>=1.8与原能量/效率；不以本线性平面参考宣称PASS，不进入H14。

一般曲面内部normal trace、相分割物理速度梯度、RZ/AMR/开边界和实际非线性功尚未证明；这些是明确剩余工作，不用uniform参考覆盖。任何原gate失败保留数据并继续H13，不能放宽断言或自动回退。

## 5. 旧结果的可比性和权限

旧程序/source/输入/所有失败候选和原物理门槛都封存，新结果用新身份单列。对旧source的上述对照仍可逐文件复算，材料M、C、原初值/频率/能量定义完全可比。源规则一旦变更，不能声称新轨迹与旧值逐字相同，也不能改写旧PASS或FAIL。

总纲§2.1要求先形成设计变更并写明离散性质、验收、可比性，且规定“不得由 Codex 在实现过程中自行改变方法”。用户现已明确授权本联合提案，允许在本提案范围内解除H13共同source冻结；不扩大到独立source精修任务。本提案不是部署或H13通过证明。

## 6. 实施前2:1/RZ兼容性补充

092915只读复核完成48个非均匀/径向真测度P3端点场景、6个轴相容场景及5个parent/subface控制。值/缩放切向导数的多项式再现误差约4e-15；RZ以真实径向权重积分，不对零面积轴面作除法。

单一parent Q4无法同时满足四个端点jet与两个独立子面平均，五个控制的不可相容残差0.087–0.150必须保留。现有CCH的子面 `16*s²*(1-s)²` bubble在端点值与导数均为零，其XY/RZ归一化平均均为8/15；按原15/8修正可闭合真子面平均，且保持parent真通量，误差4.45e-16。

因此AMR仍保留现有piecewise子面bubble与supportKnots，不能把完整source父面trace宣称为单一Q4，不能把子面独立均值当成新的global DOF。此项只证明代数表示兼容，还需真实source/master/normalPolynomial/characteristic/FirstJet测试。


## 7. 093638真实source候选状态

授权后两个production文件已实施面矩端点准备，旧横向/mixed/edge/center与parent/subface逻辑保留。完整四面使用P3物理面平均；不足四面用可观察低阶矩加原P2高阶补空间（cubic seed为零）。这是一条确定的准备规则，不是运行时fallback。轴/typed wall保持原物理约束，一般端点使用单侧真实矩，不暗加free-slip。真实cellJet新增7136端点及920零trace检查通过，原source契约保持。旧弱载荷尚未替换，候选尚未验收；第1节隔离参考为历史证据，不代表实际联合周期通过。


## 8. 163744实际单侧端点反例

局部stencil已与实际apply/cellJet一致，九次原鞍点/固定弱力通过独立解析矩参考。一般单侧闭合的实际normal末档1.738阶、切向导数0.528阶。完整18算子及独立full-state指数确认：actual联合九组均有两个增长分支，原能量随加密严重放大。全部image功差等于两个接触端点[B*g]边界项，预测缺陷不超过1.12e-14；初始g端点为零的单次功测试不足以证明动态稳定。

第2节条件偶对称参考保持历史事实，但不能证明当前一般单侧闭合通过。候选未验收、生产弱牵引尚未采用。后续须在同一授权联合方法中明确当前geometry/source/weak traction端点相容条件，不能自动给所有impermeable边界添加free-slip或反射，也不能删除两个增长分支。完整常应力矩阵及原物理门槛保留。
