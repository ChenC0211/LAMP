# H10i1 / 同源零迹的几何与材料契约

2026-10-04，实施前设计；实际能力审计后补齐以下文件预算，最终设计裁决见[报告](H10i1_Coincident_Phase_Trace_Clip_Contract_Report.md)。本轮0 production/0 test，只扩展既有研究driver的被动能力观察；不构造新通量或更新接受状态。[TODO](H10i1_Coincident_Phase_Trace_Clip_Contract_TODO.md)。原材料四硬约束、Pure FirstJet独立、无库存/bubble、原H10b门槛保持。初稿与修改前TODO独立保存在baseline，不将后补清单冒称试验前已经完整。

## 1. 两种迹及材料规则的推导

设真实轴对齐面g的规范normal为正坐标方向n；物理normal速度a=V·n由唯一advector给出。定义同一Q4的两侧相迹`chi_minus/chi_plus = lim(epsilon→0+) 1[phi(x ∓ epsilon*n)<0]`，在正长度零迹上使用实际两侧native叶的多项式。它们不使用采样epsilon；离散接口只返回已认证的相符号。普通非零面迹仍按phi<0判相，两侧相迹在几乎处处一致。

**几何clip迹**取裁剪区域内部侧：上/右边界取minus，下/左取plus。`phi=x-a`时，[a,b]内部气相，[c,a]内部液相；边界phi=0不改变体积。Green积分只取真实液体区域的约化边界：在clip边界上的零曲线由fixed clip项唯一计入，不再加一份free curve；内部free curve按原liquid→gas normal积分。两侧同相的偶阶零迹不是相边界。有限孤立点不贡献线测度，但保留来源/事件记录。

**材料面迹**取因果输运迹：

\[
Q_l(g)=\int_g\bigl(a^+\chi_-+a^-\chi_+\bigr)dA,
\qquad Q_g(g)=\int_g a\,dA-Q_l(g),
\]

其中a+=max(a,0)，a−=min(a,0)。这不是为反例选half分配：固定瞬时实际V，解`phi_t+V·grad(phi)=0`的右时间极限；重合平面上切向phi为0，backward characteristic从流入侧到达该面。一般正规零面`phi(x,t+delta)=−delta*a*partial_n(phi)+O(delta²)`，给出相同选择。高阶法向零迹用完整normal Taylor的首个非零项，backward normal位移的奇偶号给出两侧相。当a=0时液气率都0，不需发明面密度；时间端点的此右极限约定不是已经实现的时空积分。此推导要求同一actual advector提供一致的共享normal trace；若Interface实际特征速度不同，R如实保留该差异，不能通过面迹反推或改变几何。

该材料选择只在两侧迹不同的正长度零面改变普通phi<0积分；不是momentum donor耗散，也不改变KEP/skew装配或Interface输运。速度在零段内变号时按**实际advector多项式**的正/负区间积分，与相迹事件共同overlay；不能用整段平均速度选相。没有可用流入侧（真实外边界/selected支持之外）时，非零流入要求正式外部相数据port，否则拒绝；内部面不得改用另一侧、中心平均或ghost。axis物理面积0且径向速度0，所有材料率严格0。

同一规范面只计算一份Ql/Qg，两邻域取反号；normal反转同时交换minus/plus和a符号，率严格反号。actual子面区间可加，故AMR粗面为子面和。F_rho仅由一个owner作rho_l*Ql+rho_g*Qg；来源保留同一geometry/advector/layout/full clock。KEP对任意这样的共享F保持`N=diag(q)/2+K, K^T=−K`，仍用实际`R=Delta_m/dt+q`与uR；不附加R=0、quarter或精确commuting。

## 2. 一般Q4的认证与退化边界

Core现有`QuarticLeafData::restrictLocal(index,x,y,dx,dy)`支持同一homogeneous Q4的0–4阶精确偏导，物理尺度的denominator严格正；外部rounded local dofs不得替代它。对轴对齐段保留原rational chart cuts，value restriction全部系数精确0才是恒零迹；小值/区间包0都不证明恒零。正长度native incidence的value trace需逐一验证相同零迹，不能只选一个叶。

每个实际法向侧σ按native闭盒与段的精确关系选择叶：段在叶内部时同一叶服务两侧；段在叶边时分别选择真实incident叶。非选中或inactive侧标记unavailable，不延拓。计算`p_k(t)=partial_n^k phi|g`；沿该侧的第一条**不是恒零多项式**的p_k，k≤4。由原endpoint-relative有界隔离器保留p_k的全部根（包括exact multiple），在各正长度开区间认证其符号。该侧相符号为`sign(p_k)*σ^k`。孤立p_k根点的高阶退化不决定积分值，也不能消灭其事件记录；精确全零normal Taylor表示该叶在二维开区域恒零，明确拒绝。

导数的root隔离、refine与dyadic restriction必须调用**同一阶偏导**的restriction functor；不能复用原value-only fallback把恒零phi再次当作p_k。positive clearing factors不改变符号。两侧事件以retained chart作精确overlay，不能用rounded global端点排序/合并；预算耗尽、算术尺度超界、缺支持或未认证符号均拒绝。Range保持fromEnd/offset/length；逆序只调整chart取向，不以1−root或snap重造细相。

非仿射正规示例`phi=(x-a)*(1+y+y²)`具有整段零迹，但normal一阶系数严格正，两个侧迹分别液/气。偶阶`(x-a)²*(1+y²)`两侧气，三阶`(x-a)³*(1+y²)`两侧异相：Core测量port可表达，既有publication若因临界梯度拒绝仍如实拒绝，不能为此降低正规曲线要求。`(x-a)*(y-b)`需保留b事件并分相，交叉点的完整界面publication能力另由原正规拓扑规则决定。近零非零restriction走普通根/相算法；零阈值不扩大。

## 3. typed来源port

拟在现有`LevelSetMappedSegment.h`定义独立只读`NativeQuarticTrace`结果，由同一LevelSetMappedSegment owner的显式factory准备，持有owning immutable Snapshot。禁止给原bool negative Range塞零相或把旧constructor改成隐式half规则。普通constructor、roots/variation/integrate的历史契约继续明确拒绝恒零。

新typed结果记录normalAxis、原始方向、ordinary phase ranges、coincident ranges及每个侧的availability/nativeCell/chart/order/sign事件，提供usesSource、root endpoint身份与同阶refined查询。用于clip时只传region取向选内部迹；用于material时只给same-source/advector owner读两侧信息，不把速度/密度或Hydro时钟塞Core。孤立事件不是一条free interface，不冒充旧Root有限零库存；正长度零range不能转为一串采样点。

几何query要区分scalar active支持与Hydro selected支持：Core只认证源的active leaf，不知道selected网格；Hydro的typed owner须以实际选中域将外侧标记为外部port，不能因Snapshot在外侧有数据就暗自消费。裁剪一个真实源叶的内部侧总有实际支持，inactive正面积交集仍拒绝。Core正常数阶数最多4，局部限制/分片只增长实际查询事件；按现有root/refinement预算处理，不新增全局DOF。

Core另提供geometry-only `integrateNativeQuarticRegionMoments`：输入同一source的实际graphs、上述typed cuts、真实regions、XY/RZ和共同Green origin；按region输出两相axial/radial closed Green moments、error estimate和真实arc/refinement计数，无修正Jacobian。复用原regular graph root/moment engine，ordinary相仍走同一restriction/controlled GL路径；不另造Fubini几何或第二Q4。

coincident cut若对应实际free graph：以相同native/chart、exact zero restriction及原graph唯一正规分支证明贴合，在该region的外clip边界上只计fixed项；不得按浮点坐标相近跳过曲线。内部零曲线保持free项；相邻region各取内部迹，区域和恢复全体积。事件支撑和图分支不明确时请求原typed root refinement，不能凭GL点owner投票。legacy敏感度接口对拓扑/边界零迹仍拒绝，不能把缺导数标成0。

闭合账本不能继续假定几何内部迹为一份反号共享材料迹。重合内部面两侧液相迹可不同，逐region fixed项和应满足`sum(fixed)=exterior_fixed+coincident_interface_moment`；assembly误差应减去这项真实内部相界面矩。原`BoundaryMoments::assemblyClosure=sum(fixed)-exterior_fixed`只描述旧共享迹账本，不能静默重解释或强制新相几何满足0。新`DiscoveredTraceVolumes`显式输出该几何clip账本和coincident项，FirstJet派生volume报告增加该项/trace说明；capacity总量仍按真实盒闭合。材料内部面反号抵消则始终由唯一F保证，与此几何账本无冲突。

## 4. 原生面重合与publication范围

H10i反例的publication已成功，仅dual内部clip失败；但要求补原生面重合不能假定其也能发布。现有MappedFacePartitions只持旧LevelSetMappedSegment，native component discovery也有coincident edge拒绝。本轮将实际核验两个位置×正负phi；若原生面publication拒绝，需要独立的**派生拓扑/测量**迁移任务，不能隐藏在H10j的两文件预算。

正规native zero edge的拟规则：精确value零迹和非零normal一阶认证证明真实相界面支撑，按实际native edge片的共同端点/incident来源去重；canonical edge arc只有一个几何owner，两个侧phase分别保留。native-side-attached graph需独立typed支撑（side、tangent retained extent、native owner/incident），不得伪造无限根的rootOrdinal。regular interior graph继续原规则。phase相同的零边不造界面，inactive侧不补延拓。frame/normal只在canonical真实叶的保留chart求值，保持原严格正规性；不改primary、Q4、transport或accepted clock。

native整叶内部phase证书须在完整zero set上保持原严格normal partial认证：每条normal截线上所有零的穿越方向相同，零边已经给一个零，因此不能有第二个内部零；连续性与内部侧相迹认证整叶内部为该相。此处“phase纯”指二维内部，可有零边；不能仅因零边normal为正就忽略其他内部环/分支。否则继续原完整regular graph inventory或明确拒绝。canonical edge端点用实际subedge/vertex来源与ordinary exact端点incidence绑定，正规component保持原2-valence/declared boundary规则；边附着graph使用真实闭叶chart及验证过的支撑，不取消普通graph的interior sample要求。旧公开手工Graph不能仅设置一个flag绕过source零迹/正规认证。

只读面测度在零相面没有严格phi<0全覆盖，须区别于material rate：若需要无速度的诊断面测度，采用有明确定义的**对称空间迹**(chi_minus+chi_plus)/2，即两侧相迹在对称normal mollification下的极限，并显式记录trace kind；它不能被旧face-constant rate当作新材料Ql。在新spatial factory迁移前，旧evaluate遇此kind应拒绝；material仅消费第1节因果迹。缺侧的边界几何诊断使用明确inside-domain迹并记录boundary kind；流入材料仍需外部数据。

## 5. 固定实施顺序与预算

本轮16个实际场景中：dual内部位置8个publication成功但mass拒绝，共1,408段/40零段；其两侧16个clip盒均拒绝。native-face位置8个publication在geometric inventory阶段同因拒绝。因此不能只实施clip后就跳到H10j。新增源能力按下表分开，均复用原库/CMake/短契约，无新工具或变体算法。带Core exact数据的实现仍归既有strict geometry TU；legacy敏感度/Conservative入口的拒绝与原证据保持。

| 单一任务 | production文件（完整预算） | 短契约 | 消费与边界 |
| --- | --- | --- | --- |
| **H10i1a：下一唯一任务** | `include/Core/Geometry/LevelSetMappedSegment.h`、`include/Core/Geometry/Detail/LevelSetMappedSegment.h`、`include/Core/Geometry/Detail/LevelSetMappedCellLayout.h`、`src/Core/Geometry/LevelSetHermiteGeometry.cpp`（4） | 既有`interface_root_enclosure_contract.cpp`（1） | 准备同源NativeQuarticTrace/normal Taylor认证；现有LeafData的0–4阶restrictLocal已足够，不复制其实现；未迁移消费者 |
| H10i1b：闭合clip moments | `include/Core/Geometry/LevelSetMappedInterface.h`、对应`Detail/LevelSetMappedInterface.h`、`src/Core/Geometry/LevelSetHermiteGeometry.cpp`、`include/Modules/Interface/Detail/HermiteGalsRegionMeasure.h`、`src/Modules/Interface/HermiteGalsPhaseRegion.cpp`（5） | 既有`interface_phase_region_contract.cpp`、`cch_momentum_geometry_contract.cpp`（2） | 新geometry-only region moments和真实clip消费；先验收已能发布的dual内部8例，native publication缺口仍保持；正规free arc/clip归属不计算Jacobian |
| H10i1c：正规native零边拓扑 | `include/Core/Geometry/LevelSetMappedInterface.h`、`include/Core/Geometry/LevelSetInterfacePanels.h`、`include/Core/Geometry/Detail/LevelSetNativeInterfaceChain.h`、`Detail/LevelSetMappedInterface.h`、`src/Core/Geometry/LevelSetHermiteGeometry.cpp`（5） | 既有`interface_root_enclosure_contract.cpp`、`interface_hermite_gals_publication_contract.cpp`（2） | 显式`discoverNativeTraceComponents`；typed native-side graph支撑、唯一正规edge arc与两侧phase/端点incidence；Core frame/sampling/geometry moments支持此来源，旧default不改 |
| H10i1d：只读半面trace cache | `include/Modules/Interface/InterfacePhaseGeometry.h`与对应`.cpp`（2） | 既有`interface_hermite_gals_flux_contract.cpp`（1） | `MappedFacePartitions::prepareNativeTrace`和`traceHalfFace`显式只读cache；普通`halfFace`仅ordinary legacy view，零迹访问拒绝；同源refinement后重新计算派生量 |
| H10i1e：FirstJet派生出版测量 | `include/Modules/Interface/InterfaceCellMeasure.h`、对应`.cpp`、`include/Modules/Interface/HermiteGalsPhaseVolume.h`、对应`.cpp`、`src/Modules/Interface/HermiteGalsPublication.cpp`（5） | 既有`first_jet_contract.cpp`、`interface_phase_region_contract.cpp`（2） | 独立typed `discoverNativeTraceVolumes`/几何clip账本；FirstJet选择同源新派生测量，拓扑/面cache/frame来源统一；primary/Q4/transport/clock/正规接受要求不变，Conservative/敏感度不fallback |
| H10i1f：同源几何消费者 | `src/Modules/Interface/HermiteGalsPhaseRegion.cpp`、`include/Modules/Hydrodynamics/CchPhaseGeometry.h`与对应`.cpp`（3） | 既有`cch_phase_geometry_contract.cpp`、`cch_momentum_geometry_contract.cpp`（2） | Region准备new complete topology；H07读typed普通/重合对称诊断面迹并显式trace kind，旧evaluate遇重合kind拒绝；16场景的publication/真实质量闭合，仍没有新材料F |

每项先写自己的TODO/设计/gate，不依赖本轮TODO作为跨任务无限预算。H10i1b以后单独冻结完整ordinary回归；H10i1f完成后重验64普通来源布局、全部16实际重合场景和真实identity/clock/selected。随后才复核H10i设计并实施H10j的空间factory：材料取因果迹，不把对称诊断测度乘u当作Ql。其后H10k/H10l移植唯一材料owner与实际R/E_sp，最后最终native PDE/原H10b收口；来源PASS不关闭H10。

Core trace新port的恒等式/测度/源拒绝为1e−12实际尺度；非零值与normalTaylor符号必须exact认证，无epsilon零判据。真实clip phase volume按独立平面/非仿射closed form≤1e−12实际capacity，不裁负值或扩大账本容差。原4反例及后续native-face位置、正负相、axis/selected/AMR实际共享/子面核验均必须包含，原64普通source布局保持。

causal flux按独立平移weak transport右极限核验正负/零法向率与速度变号，区间和/反向严格自洽≤1e−12 gross；常速度/uR/KEP原gate保持。新phase工厂的全部空间矩、PDE精度及原H10b门槛仍归H10j以后。一般临界publication、不支持外部流入/未认证退化和时间跳跃如实记录，不能把设计PASS称为这些能力已可用。

## 6. 独立闭式反例与身份检查

令phi=x-a-v*t、h>0，两盒B−=[a-h,a]×[c,d]、B+=[a,a+h]×[c,d]，液体在左，单位面面积A=d−c。v>0时右时间质量导数为`dm_-/dt=0`、`dm_+/dt=(rho_l-rho_g)*v*A`，共享面F=rho_l*v*A；v<0时为`dm_-/dt=(rho_l-rho_g)*v*A`、`dm_+/dt=0`，共享面F=rho_g*v*A。远侧纯相流率分别rho_l*v*A和rho_g*v*A，故独立continuum右导数与第1节因果迹相容。若用对称half材料迹，将各盒引入相反的`(rho_l-rho_g)*v*A/2`瞬时缺口，不能据对称几何诊断选材料。负phi号互换液气；RZ轴向平移把A替换为真实pi*(r_hi²-r_lo²)，推导相同。这里只是规范的闭式参考，不要求一般GALS实际R为0。

几何参考另由phi=x-a：Vl(B−)=capacity，Vl(B+)=0；XY Green potentials为(x-x0)e_x与(y-y0)e_y，RZ meridional potentials为pi*(r²-r0²)e_r与2*pi*r*(z-z0)e_z，分别的div都是正确volume measure。clip内部迹配合不重复free arc得到此结果；边界点值的half不能作为clip内迹。重新合并两盒时，原内clip上的相界面转为一个内部free arc，其Green矩与原fixed矩相同，保证体积可加。平面乘正的非仿射factor不改变相区域，提供与新kernel独立的Q4参考。

normal反向：令n'=-n、a'=-a、chi'_-=chi_+、chi'_+=chi_-，直接代入得Ql'=-Ql、Qg'=-Qg。分片加回由同一有符号积分的可加性给出。对任意该F与实际q/R，skew双线性`w^T K w=0`仍成立，材料动能项由q/R原恒等式给出；没有借用旧pressure→dual交换。以上证明选择规范与代数相容性，不是完整sharp-interface PDE、时间阶或无黏稳定性证明。
