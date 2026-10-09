# H08a2r3 / 候选设计：一个兼容native面均值的正则advector

2026-10-05，**实施前候选设计**；后续[可行性报告](H08a2r3_Common_Regular_Advector_Report.md)为PASS，尚未生产部署。以下保留运行前设计及按发生顺序记录的补充门槛。[任务书](H08a2r3_Common_Regular_Advector_TODO.md)。[H08a2r2](H08a2r2_Skeleton_Motion_Selection_Report.md)证明原分片B的单侧选择不能达到预设二阶骨架梯度；改动应前移到运动与材料共用的空间来源。本项先验证构造可行性，不同时部署所有消费者。

## 一份来源与明确放弃的旧性质

拟由实际native normal均值构造一个C1分片多项式物理速度u_h。只有这一个field及其解析导数产生characteristic的point/J和primal/dual atom上的法向积分。保留同一native状态、prepared geometry、boundary identity、advector identity及完整clock；不创建全局求解的C1速度未知量，也不让Interface拥有另一份速度。

保留native面均值是该候选的运动学相容选择，确保原D_I及C=−D_I^T不变，不增加材料代数硬约束。当前`cellDivergence()`既是D_I/V又等于各分片点态散度的性质必须显式更改：一般C1速度的散度连续，而任意native输入的D_I/V逐cell不同，不能同时要求它等于不连续的分片常数。新候选仅通过Gauss保持**cell积分散度**与原native通量和；点态散度单独查询/诊断，不能假称恒定。不可压缩native约束与这种点态离散余项须分别检查。

材料仍仅四硬约束：唯一共享F、内部反号/实际AMR子面和、来源一致、KEP/skew。真实R、E_sp、G及uR保留；不强制逐dual commuting或R=0。几何状态仍为独立Pure FirstJet，没有库存、材料质量更新或界面体积拟合。

## 候选构造

先从真实面均值functionals作局部缩放P2重构，给共享master的速度值/导数，再用同一AMR master/hanging约束形成C1分片多项式base u_0。可以复用Core中中性的Hermite基函数/几何约束，不能借用Interface已接受phi状态作为速度owner。XY在矩可辨识的支持上至少精确再现所有总次数≤2速度多项式；RZ以物理正则基为准，先固定包含u_r=a*r、u_z=b+c*z+d*r²，以及本次r⁴剪切的梯度收敛研究。RZ轴满足u_r=0及正确奇偶导数。rank与系数范数必须可审计；不能把缺秩点静默降阶后计作P2通过。

这个base如何在真实2:1与selected/disconnected边界上保持共享C1、可辨识矩精度和有界准备成本，是下一项首先实测的事实，**尚未有实现保证**。若one-ring缺秩，只按统一实际邻接规则扩大支持并记录；不得为具体算例或AMR位置写特例。对极小disconnected分区，即使全分区也可能不能辨识P2，不能把“所有几何上P2精确”作为新硬约束：候选采用统一的缩放最小范数定义，并完整记录不可辨识函数及误差，仍须验证C1、原面矩和边界；不另加运行时降阶分支。这类行不能计入P2精度通过，且不能取代有足够支持的原收敛矩阵。

再对每个真实native面f校正其均值，避免向最终材料F施加修复。设s∈[0,1]为该实际面切向参数，xi为法向距离除以两侧各自不超过相邻leaf宽度一半的支持宽度，取

\[
g(s)=16s^2(1-s)^2,\quad b(\xi)=(1-\xi^2)^2\quad(|\xi|\le1),
\]

支持外为零。g在切向两端、b在支持端均为零且一阶导数为零；b在面上为1且法向导数为零，因此两侧支持宽度不同仍C1。以canonical normal e_n定义

\[
\delta u_f=\frac{15}{8}(u_f-\overline{u_0\cdot e_n}_f)b(\xi)g(s)e_n.
\]

在XY及RZ真实测度下都有∫_f g dA=(8/15)A_f：RZ切向r为线性函数，g关于中点对称；径向面r恒定。因此该项恰好校正本native面的均值。其余平行native面在支持外，同平面的其他subface不在切向支持内，垂直面上的法向分量为零；需用实际全face交叉矩阵核实这一cardinal性质。AMR每个实际subface独立使用同一构造，粗侧仅取真实subface和，不引入quarter分配。

这是由速度面矩残差唯一决定的运动学基函数系数；不引入守恒库存、界面内部bubble或材料量修复。零输入必须生成精确零field。

若base在光滑场下value误差O(h³)，均值校正幅度O(h³)、导数O(h²)，可保留目标梯度阶。base若仅O(h²)，校正导数可能O(h)，会重现H08a2r2的问题；不能只证明校正后面均值准确。

## 边界和旧API的影响

selected impermeable law要求整条面的法向field恒为零，不能只令其均值为零。base须先满足此点态边界；对应残差及校正系数应精确为零。角点的多法向约束需相容；外部历史缺失仍拒绝，不造内部J。轴面A_f=0不能作均值除法，径向轴点由正则base决定；其余校正须检查轴奇偶性，不能仅检验有限值。

现有`CchNormalPolynomial`固定3系数；候选校正含四次切向函数，不能装入原三系数后截断。atom overlay必须含实际base及校正的piece边界，两端相对坐标、极短正区间积分、RZ真测度继续保留。`CchPhaseGeometry.cpp::ratePolynomial`及数个测试参考也固定旧次数，下一项须逐消费者登记；不能只改`cellVelocityJet`让atom仍旧。

现行`CchCharacteristic`的C1核验使用旧volume-coordinate次数，且自动零场/同时角点仍有H08a3拒绝；即使新速度C1，也必须相应修订证书与等价分支/事件规则。不能把局部速度样本连续当成这些门槛已经关闭。

因此先在一个可退役审计中验证base+面矩校正，再根据实际依赖划分≤5 production/≤2 test的部署任务。正式切换时，位置/J、atom、phase/material、primal/dual必须同时指向同一个新来源；不能在半迁移状态宣布H08a通过。原H10o最终native算子、内部affine counterexample、H10b-v2完整28窗口/84项须重验，不能沿用旧B的PASS。

## 下一项的固定检查范围

T0核对H08a2r2冻结与1379源码/输入，记录实际测试source/exe绑定；生产改动0。候选P2系数和校正系数仅在审计中生成，不进入接受状态。

T1在XY/RZ、uniform/真实2:1、n4/8/16，以及既有selected/disconnected/axis/impermeable fixture上验证：每个实际连接两侧value/J，全部native面矩交叉影响，Gauss积分与原D_I，矩可辨识支持的恒定/仿射/规定多项式再现，零场及边界。设计目标scaled缺陷≤1e−11，缩放须按实际运算工作量事先固定。缺秩行完整保留，给出独立不可辨识证书及统一最小范数结果，不计P2精度通过；所有布局仍接受C1、面矩、边界检查，失败不能调大阈值或删布局。若原规则收敛矩阵本身缺秩并达不到目标，仍FAIL。

T2对H08a2r2完整32族/96窗口原DOF点与两种方向重做空间精度，先检查局部source value/J，再用同一候选field微分实际离散map作有限时间phi/grad验证。XY/RZ均保留n8/16/32；RZ r⁴不要求P2精确，但收敛不能豁免。候选value目标≥2.5、derivative≥1.8；同一map的phi/grad目标≥1.8。时间细化需确认积分误差低于最细空间误差的10%，并独立差分同一map；不能以另行恢复的J代替。原sin跳跃输入改为新来源后，其新左右field/map应一致，并明确已更换空间重构。

记录准备时间、稀疏系数/atom数量和每次evaluate工作量。n8→16→32应有有限支持；支持数量上界与权重范数由实际几何给出并报告，不接受随h恶化而无法完成最细档。成本失败也不能只交最粗矩阵。上述阶数是候选精度验收，不是材料第五硬约束。

## 实施细化（2026-10-05，数值运行前）

H08a2r2的75 payload/11文档/1379源码输入及上一冻结46 payload已复核；上一goal turn为实测进展。临时扩展原momentum-geometry测试一个导出模式，输出同一生产native cell/face/dual段及实际FirstJet DOF、完整布局标识。NumPy原型只读取这些导出，独立构造候选；它不是生产advector，也不接入接受状态。无SciPy依赖。临时源码、exe、输入与Python源码全部冻结绑定。

Core Q4约束当前要求活动挂点的粗owner也活动；selected-domain可能不满足。原型按**实际活动full sides**定义master/hanging依赖，只约束存在活动粗侧的挂点；未活动侧不生成速度未知量。沿用同一五项一维Q4基及其tensor数据顺序，全部full-domain局部25数据与Core同一约束引擎独立probe比较；selected域另作完整系数C1核验。该支持差异必须在未来Core/运动学owner部署时明确处理，不能偷偷给inactive几何填接受数据。

每个独立datum以实际incident cells及一圈面邻接开始，逐圈扩展至可辨识或四圈/分区用尽；原面均值矩使用物理平移/尺度化P2基{1,X,Y,X²,XY,Y²}。同一缩放SVD最小范数用于所有布局。零normal边界对master侧上的该分量施加整条线为零的多项式条件；轴u_r=0、轴u_z的dx/dxy=0。反射只用于核验延拓的C1奇偶性，不用于实际轨迹越轴。全部系数/支持及缺秩证书记录。

T1的两种AMR为原central box和uniform，48个XY/RZ×n4/8/16×4种selection布局；各检验无wall、全部外部wall与单个实际外部wall三种声明。T2用原横向条带12布局及原32族96窗口。精确P2的空间误差若全部≤1e−11，标记roundoff reproduction，不对舍入噪声计算阶；RZ r⁴保留真实收敛。sin反例另外保留，并检查新field下的双侧map/J。

所有C1/矩/再现缺陷用abs(diff)/(1+实际绝对项和)缩放，目标1e−11；RMS继续沿用H08a2r2的abs(error)/(1+abs(reference))。C1以分片迹的多项式系数/Bernstein界核验，不能以几个样点相等代替。面矩与Gauss使用逐实际分片的足阶积分并与闭式校正矩对照。T2采用同一field及解析导数的固定步RK4离散map；64/128/256子步细化确认时间误差，普通双侧差分目标5e−7。若任何前提失败，保留全矩阵和真实失败；不会把原型改称生产同源时空积分。

## 非多项式补充复核（主矩阵完成后、补充运行前）

原sin反例在新field上已具有唯一普通导数，但n8该导数与连续sin的物理解仍有明显差距。为区分欠分辨率与不一致重构，新增XY uniform/条带AMR、n8/16/32全部原DOF、正向静态T=0.25的sin空间检查；原96窗口全部保留。native均值用cos差的闭式积分，source和实际map均对连续sin解检查。仍使用已有source value/derivative≥2.5/1.8、phi/grad≥1.8门槛，128/256时间差小于最细空间误差的10%；不改已有阈值，也不能用此补充代替原RZ矩阵。
