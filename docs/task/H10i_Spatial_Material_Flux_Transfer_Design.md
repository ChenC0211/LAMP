# H10i / 同一空间相通量 owner 的设计

**2026-10-04现行接口复核：** 以下保留H10i原未实施设计及FAIL历史。i1a–f已补齐实际几何消费，[H10i2](H10i2_Spatial_Material_Source_Reaudit_Report.md)独立复核后将几何复用/实际advector明确拆为前置[H10j0](H10j0_Spatial_Advector_Reconstruction_TODO.md)，再原H10j相通量。新工厂签名、kinematic piece overlay及固定预算以[现行契约](../api/H10i2_Spatial_Material_Owner_Contract.md)为准；不改下述冻结历史数字与裁决。

2026-10-04，**FAIL：同源重合相迹/裁剪契约不完整；engineering PASS**。此文件保留未实施的空间设计，实际验收见[报告](H10i_Spatial_Material_Flux_Transfer_Report.md)。[任务书](H10i_Spatial_Material_Flux_Transfer_Design_TODO.md)限制本轮 0 production / 0 test / 1 research；本轮只测量当前 source 和固定后续实现。本设计依据用户指定四项材料硬约束；不改 FirstJet、质量、D_I/C、时间积分、冻结 PDF 或原研究门槛。H10h 的 115 payload/28 当前交付已独立核验；原 FAIL 保留。

**实施前缺口：** 64 个移动圆形/球形、selected 布局支持下述普通相区间与物理测度；不能由此推论全部合法几何可用。4 个实际 Pure FirstJet 平面 publication 成功，但界面落在 cell 内部的 dual 面上时，现有 Core 查询拒绝正长度 `phi=0` 段，真实 `CchMomentumMass::build` 也同因失败。704 段中20段触发该拒绝。下述严格 `phi<0` / `phi>0` 积分在这种面上没有覆盖全测度的分相定义，不能用任意 half trace、微扰或平均材料补齐。下一[H10i1](H10i1_Coincident_Phase_Trace_Clip_Contract_TODO.md)先区分真实裁剪体积的单侧边界迹与共享材料面通量的零相迹，冻结公开typed契约及实现预算；H10j/k/l均等待。此缺口是内部数学/API问题，不是外部BLOCKED，也不否定下述运动学lift。

## 1. 空间 owner，而非标量材料率的几何延拓

H10h 的纯气体 dual 有真实质量变化 0，却有旧 q=R=12499.9875；其原因是材料散度按 pressure 体积强制分配，不能靠 skew 矩修正消除。新设计先准备**与密度无关**的 compatible advector，再用同一 FirstJet 相区间定义空间相积分。唯一材料 owner 只作 `rho_l*Q_l+rho_g*Q_g`，primal/dual 是它在不同真实面的积分视图。

**材料硬约束仅四项：** 单一 F_rho owner；同面反号共享/实际 AMR 子面求和；几何/advector/layout/完整 clock 来源一致；KEP/skew 与真实 R/uR。下面的局部运动学、矩再现和求积阶是该候选的构造与精度验证，不额外要求材料 commuting、逐 quarter 散度或 R=0。

## 2. 保留 D_I 的总量运动学，补充面内信息

native 面标量 u_f 对应既有总量 `Q_f=A_f*u_f`；不变更接受速度数组、压力伴随或 D_I。为每个真实 native 子面构造一份线性切向 normal trace，保持其物理面积平均为 u_f。切向斜率从同分量实际 native 面的一般 one-ring patch 的三列最小二乘获得，物理 functional 使用真实面积平均，不用解析解；几何矩阵独立于速度/密度，prepared 后复用。

XY 拟合空间为 `{1,x,y}`。RZ 采用体积坐标 `s=r^2/2`：轴向拟合 `{1,s,z}`；径向拟合 `u_r=r*(a+b*s+c*z)`，正半径面提供 u_r/r 数据，零面积轴面只施加 u_r=0，不做除零。以实际 physical-face mean 为插值 functional，最后在该面把切向一次项减去真实 weighted mean；normal 总积分仍为 A_f*u_f。

RZ 的面积平均一般不等于面中心点值。本设计不同时强制“点值 affine 精确”和“真面积平均精确”这两个可能冲突的条件：XY affine exact；RZ 原轴正则 affine family `u_r=a*r, u_z=b+c*z` exact，偶径向轴向项用真实面积平均 functional 验证。原 benchmark 的点采样参考和 H10 精度门槛不修改；一般 RZ 点采样与面平均差异须在后续 MMS 中按实际 O(h²) 余项报告，不能重标参考来制造通过。

one-ring 是实际选中压力邻域及其真实面，不跨非选中区域造 ghost。按归一化坐标和真实测度对同一三列矩阵统一求 Moore–Penrose 最小范数解；rank deficient 时仍是同一优化问题，不选择另一算法或扩大解析数据。逐面报告 rank、切向 functional 的投影缺口及权重范数：可识别的 affine trace 才声明 exact；极粗 selected 的单列或细连接区域可能缺少切向信息，不能强行宣称其二维 affine 再现。这个缺口是精度诊断，不增加第五项材料硬约束；对应极粗 selected 场景仍须共享/KEP/有限值自洽，完整精度须在足够几何信息的加密序列或后续真实边界 port 上验证。非有限数据、轴冲突和数值分解失败明确拒绝，无隐式 fallback。

## 3. cell 内的显式 divergence lift

用直角体积坐标 `(X,Y)`：XY 为 `(x,y)`；RZ 为 `(s,z)`。对应通量向量 B：XY 为 `(u_x,u_y)`，RZ 为 `(r*u_r,u_z)`；RZ 所有积分统一带 2*pi。每个实际 pressure cell 的两侧共有 normal trace 已确定，不同子面只作为 piecewise trace 的实际区间。

先作 normal 方向的线性 transfinite 延拓 B0，切向按实际子面断点分片。令 `d_c=(D_I*u)_c/V_c`，在体积坐标上 `f=div(B0)-d_c`；边界总积分给出 `integral f=0`。以下是显式局部零 normal lift，a/b 为 X 界、c/d 为 Y 界：

\[
\bar f(Y)=\frac1{b-a}\int_a^b f(X,Y)dX,
\]
\[
H_X(X,Y)=\int_a^X f(t,Y)dt-(X-a)\bar f(Y),\qquad
H_Y(X,Y)=\int_c^Y\bar f(t)dt.
\]

`div H=f`，Hx 在 X 两侧为零，Hy 在 Y 两侧为零；`B=B0-H` 保留所有 primal normal trace 并有 `div B=d_c`。分片多项式按实际区间积分，不解全局方程、不创建新 primary 或高阶 unknown。在 r=0，Hx=O(s)，B_r=O(s)，故 u_r=B_r/r=O(r)；轴值取其正规极限 0，不能直接除 r=0。normal trace 在相邻实际 cell 间共享，所以 B 的 normal 分量相容，切向不要求跨 cell 连续。

这个约束作用于**纯运动学总量**。材料 `div(rho B)` 随真实相区间变化，不被设置为体积权重乘 pressure 材料散度；旧 quarter SPD 材料延拓完全不属于新构造。

XY 全局 affine 输入的 trace 与 B0 exact、f=0；RZ 上述轴正则 affine family同样 exact。光滑无散度输入的 native 总量满足 D_I*u=0 时，lift 维持局部零物理散度。一般 smooth 输入的误差还依赖面 functional、patch 条件数和相求积；形式再现不是完整 PDE 二阶证明，须实际测试。

## 4. canonical 段、相积分与所有权

所有 primal 真子面和两个 velocity component 的 dual segment 先按同一物理直线的实际端点做共同 overlay。正向 key 为 `(normalAxis,plane,tangentBegin,tangentEnd)`；输入段由一组有向原子段索引表达。相同原子段只查询/保存一份液/气积分，粗面与各 component 的视图只求和；不能按几何面积拆一份全 cut 面材料标量。

primal 原子段使用共享 normal trace；cell-interior 原子段使用该实际 cell 的 B。来源必须明确为 primalVelocity 或 ownerPressure 二者之一，与现有 CchMomentumGeometry 一致。不同 component 在同一物理位置使用同一求值；内部 segment 两个外流号相反。selected 边界仍是真实显式边界，不能暗自把相反端口配对成生产周期边界。

相区域查询使用当前 publication 的 `state().geometry()` owning view，aliasing shared_ptr 仅为了调用既有 Core `LevelSetMappedSegment`，同时保留 view/publication 的 identity/clock；不复制 snapshot、不访问 Interface Detail。既有 `HermiteGalsPhaseQuadrature.cpp` 已采用相同 owning alias 模式。mapped ranges 的 `fromEnd/offset/length` 保留短相，不通过两个全局 root 点相减重建区间。

\[
Q_{l,g}=\int_{g\cap\{\phi<0\}} v\cdot n\,dA,
\qquad Q_{g,g}=\int_{g\cap\{\phi>0\}} v\cdot n\,dA.
\]

对每个正长度相 range 在其 endpoint-relative chart 积分局部 polynomial trace。XY 候选 trace 次数不超过 2；RZ axial trace 在 r 中不超过 4，乘 r 后不超过 5，三点 Gauss 在每个实际相 range 足够，或使用相同 shifted polynomial moments；权重为正真实面积，速度/通量可带符号。区间极小而无可表示内部全局点时采用 chart 多项式积分，不能 snap、丢弃短相或改 scalar。全部 root/degenerate query 拒绝规则沿用 Core，不静默近似。

由同一源上的 Ql/Qg 定义 F_rho；不单独拟合 density，也不加入材料散度目标。纯气体 box 的每条面若在同一源上为纯气体，其 F=rho_g*Q，局部散度来自 compatible kinematics；H10h 的常速度 witness 因此应 q=0。该性质在实现前冻结为相尺度验收，不能使用邻近 rho_l 的 gross 归一化掩盖泄漏。

## 5. API 与后续拆分

本轮不新增production类。H10i1来源契约关闭并复核本设计后，H10j才在现有`CchPhaseGeometry.h/.cpp`给同一`CchPhaseRate`增加显式`prepareSpatial(mesh,phase,momentumGeometry,actualVelocity,identity,fullClock)`工厂。prepared数据持有原geometry/clock、实际advector、与密度无关的reconstruction、canonical atom registry及原native/segment顺序的primal/dual phase-rate视图。公开`momentumGeometry()/dualFaceRates()/spatialReport()`只读typed查询，不接受caller自行重贴来源的裸候选率数组；header forward declare CchMomentumGeometry避免include cycle。该实现预算需按H10i1落定的正式来源port复核，不能隐含修改Core/Interface。

prepared owner 还持有不可伪造的 `CchAdvectorIdentity`（每次实际输入准备生成新 identity，不是新的时间推进状态）。H10k 把它加入 `CchMaterialSources`，H08/H09/H10 验证链核验该 identity，拒绝同一几何/clock下另一份 advector revision。CchNativeValues 本身只有 layout，没有物理 clock；调用者提供实际 provenance 的责任保留，新 identity 不能冒称已验证 H14 的 coupled commit。

最终 native PDE 所需 actual advector point 查询也来自该 B：normal 分量取同一 primal trace，切向分量取实际相邻 cell 的 B 在该点的几何交集加权平均（不按材料密度加权）。轴值用正规极限；同点所有视图共用这一约定。后续移除研究中的第二份跨面恢复，不从解析参考或 moment target 修正产生 V。
迁移拆为三个任务，避免同时改六个 production 文件或三个以上短测试：

| 任务 | 单一事实与允许文件 | 消费与退役 |
| --- | --- | --- |
| H10j | 新 prepared 空间相通量；2 production（CchPhaseGeometry.h/.cpp）+1 现短契约 cch_phase_geometry_contract.cpp +既有研究 driver | 独立验收新 factory；旧 evaluate 仅维持现消费者的历史面常值测量，不是失败后备。新 factory 从不调用旧材料 quarter 延拓。 |
| H10k | 将原 phase/material 消费入口迁移到新空间 owner；4 production（CchPhaseGeometry.h/.cpp、CchMaterialTransport.h/.cpp）+2 现短契约（phase_geometry/material_transport）+研究 driver | 加入 advector identity 来源检查，删除旧 A_l*u_f 空间算法，evaluate 委托同一 prepareSpatial；无 geometry 输入的入口构造自己的同源 momentum geometry。此阶段不重设计 dual 材料转移，旧 H09 仍冻结，H10 未关闭。 |
| H10l | 唯一 F_rho 的 dual 消费及真实 E_sp/R；4 production（CchMaterialTransport.h/.cpp、CchDualMaterial.h/.cpp）+2 现短契约（dual_material/kep_convection）+研究 driver | H08 同时拥有 primal/dual 材料率与 dt 积分，H09 只读取该材料 owner 的 dual 视图；删除 localFlux/quarter SPD/面积材料分配。 |

H10l 的真实 dual 消费者必须在 prepareSpatial 时传入与三层 CchMomentumMass **同一 instance** 的 momentumGeometry；迁移两份短契约 helper 和研究调用者，不把不同 geometry 重新贴 identity。仅用 primal 材料时可用自建 geometry 的 evaluate 入口；后续构建 dual mass 时应消费 `rate.momentumGeometry()`。所有消费者经此明确路径绑定，不能以“坐标看起来一样”替代 identity/完整 clock 检查。

唯一材料 F_rho 由 H08 对两个视图统一密度乘积：H09 不另乘密度、不另积分，不拥有另一组独立材料率。真实质量三层仍为既有 CchMomentumMass；中点 dt 近似、consistencyImpulse/materialEnergy 保持。固定 mesh 上 `q_d=incidence(F_d)`、`E_sp=q_d-T*q_p`、`G=Delta_m_d-T*Delta_m_p`，`R_d=T*R_p+G/dt+E_sp`。旧 quarter 精确交换字段退役或明确改名为诊断，不能让消费者继续假设 E_sp=0。

每项先写 TODO 和固定 gate 再实施。CMake/工具/字典不改；旧数学报告与冻结 PDF 保持，后续 production 路线唯一，不保留无消费者的旧 numerical 代码。空间精度前置仍必须先过；H11/H08a 等待。
## 6. 实现后固定验收与本轮限制

H10j：T0标准 Release/architecture/layout；T1来源与同面共享 1e−12 gross、D_I 总量恢复和 affine 再现 1e−12、正物理相测度闭合 1e−12、三点 Gauss/shifted moment 独立核验1e−12。XY/RZ uniform/2:1 n4/8/16/32，以及原 selected/轴、正负/零速度、相 crossing、小相、重网格来源拒绝，全部实际 native/segment 保留。MMS包含 physical-average 输入与 RZ point/mean O(h²) 差异，不宣称全部 trace 点值 exact。

H10k：上述 phase gate 保持，新增同几何/clock不同 advector identity 的拒绝、实际来源/质量/中点材料缺陷1e−12 gross、原光滑时间二阶及真实 crossing 非零缺陷。原 phase 参考按新空间 functional 更新，不能删除来源、微小相、正负速率、实际GALS三层或拒绝场景；完整phase/material/dual/KEP回归必须记录。

H10l：原真实质量/缺陷定义、uR/材料 KE/边界功 1e−12 gross，常速度 1e−10；H10h纯 gas witness 的 `abs(q)/max(rho_g*sum(abs(Q)),rho_g*|V|*h)` ≤1e−12（真实RZ面积替换 h），三层 phase/mass/q/R 来源逐项核验。rho=1/1e3/1e6、16布局/48材料组。混相瞬时空间缺陷 R_sp=dm_geo/dt|mid+q 与 phase face 矩按独立解析平移/平面及光滑 solenoidal 参考收敛（dm_geo/dt 只作瞬时参考，不用端点差代替）：全部误差达1e−12相尺度时 exact PASS，否则 n16→32 的质量 RMS误差阶 ≥1.7、最细相对量 ≤0.05、Linf下降；此 gate 针对有明确参考的空间算子，不把中点时间 crossing 误差混入或要求所有非守恒R严格为0。原 moving circle 的 R 与 E_sp/G 则完整观察，并在H08a验证时间积分。

随后最终 native convection 使用新实际 F/R；内部 affine counterexample首先 ≤1e−10 或达到二阶，再全部 native 材料/PDE矩及真实系数尺度，前置通过才原 H10b 28窗口/84gate，空间阶≥1.7、最细L2≤0.05、Linf下降、原KE/时间及4 donor控制门槛不变。不同材料下所有 cut 行 exact affine 与 skew 的不相容条件不重新列为硬等式；应报告实际二阶余项。

本轮只验证 source 能支撑该设计和公式闭合；没有构造新的材料 F/最终 K 或接受状态。若数据 port、functional 处理、RZ正规性、canonical共享或解析 lift不闭合，本设计 FAIL 并拆明确接口任务，不通过增预算常数或改原研究目标收口。

三个空间材料任务不等于 H10 的整体空间闭合：rho=1、常速度时新旧材料率可以一致，H10c 原 homogeneous benchmark 的误差问题仍需最终 native operator 的独立收口。不能以相积分/材料恒等式 PASS 替代原 28窗口/84gate，或只选择真实混相研究作为新成功目标。旧 donor 失败保留，若后来审计其判据必须另列设计变更，不能在材料升级中隐去。
