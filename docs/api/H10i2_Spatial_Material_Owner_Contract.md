# H10i2 / 空间材料 owner 的现行实施契约

状态：H10i2为实施前设计历史；H10j0/H10j/H10k/H10l1及[H10l2](../task/H10l2_Same_F_Dual_Material_Retry_Report.md)已完成实际来源与同一F迁移。quarter/legacy evaluate已退役，当前[H09 API](H09b_Cch_Dual_Material.md)为准；当前任务见[接口索引](README.md)，H10m原目标FAIL，旧FAIL不改判。

## 1. 现有来源已经具备的能力

| 现行入口 | 已证明的能力 | 不能据此声称 |
| --- | --- | --- |
| publication.state().geometry() | 原 owning Snapshot/revision/full clock；唯一 FirstJet | 已有 Hydro coupled commit |
| NativeQuarticTrace | ordinary ranges、zero charts、实际侧/normal阶数、endpoint-relative宽度与独立event/refinement | 内部zero面已有共享因果材料率 |
| traceHalfFace / faceTraces | 同源native子面缓存、selected owner与可用侧、空间 kind/物理迹 | 对称空间平均乘速度可作为材料率 |
| RegionSource / CchMomentumMass | 实际rectangle/dual/overlap，两方向Green、非负容量与原source | 材料质量可由压力库存/通量更新替代 |
| CchMomentumGeometry | 实际segment的唯一primalVelocity或ownerPressure来源、共享/边界位置与sameIdentity | 等坐标的fresh instance可自动重贴来源 |
| CchNativeValues | 相同layout的原native标量，不含独立梯度/库存 | 数组本身已证明物理clock与H14事务 |

## 2. 先准备实际 compatible advector：H10j0

H10j0已新增 `CchSpatialAdvector.h/.cpp`，其中两个阶段属于同一个运动学 owner：

1. `CchAdvectorReconstruction::prepare(mesh, CchMomentumGeometry)`：持有原geometry/layout/D_I，准备与phi/密度/速度无关的实际one-ring physical-mean functional、同一scaled SVD最小范数权重、cell piece布局、canonical atom与有向primal/dual视图。copies共享不可变data；fresh准备有自己的identity，不使用全局裸cache。
2. `reconstruction.evaluate(mesh, actualValues, velocityGeometry, fullClock)`：返回immutable `CchSpatialAdvector`，只用原native速度生成面trace/cell B/同点查询与fresh `CchAdvectorIdentity`；保留同一重构storage、原输入值和声明来源，不重做SVD/atom拓扑。没有setter、时间推进、接受状态或密度。caller必须传实际provenance，不能重贴旧数据；H14实际事务核验仍未实现。

只读advector给出 `reconstruction()/momentumGeometry()/values()/geometryIdentity()/timeLayer()/identity()`、`validateFor(mesh, expectedMomentumGeometry, velocityGeometry, fullClock)`、actual cell/atom normal polynomial及native position的两分量velocity。合法时间字段在evaluate校验，已有owner的完整clock/source/instance在validateFor精确比较；caller声明并不替代H14事务。phase consumer不能自行重新拟合第二份V。normal polynomial以原atom的endpoint-relative局部参数表达，保留原fromEnd/offset/length；不用rounded global root点。private identity只能由实际evaluate生成。

primal normal trace的真物理平均等于原u_f，故总量仍 A_f*u_f；XY affine可识别时精确。RZ以s=r²/2、B=(r*u_r,u_z)处理，axis u_r=0，不除零；原轴正则 `u_r=a*r, u_z=b+c*z` 与可识别even-r² axial族精确。一般RZ face物理平均与原中心点采样有O(h²)差异，保留原benchmark点参考，不能重标 primary 数据或声称所有 transformed polynomial 精确。

one-ring仅由选中真实压力邻域与子面组成。full16/n>=8 selected可识别affine的门槛<=1e−12；原n4 selected2/3的56切向不可识别行保留同一最小范数问题与显式缺口，不fallback、删行或假标exact。

cell B采用原 transfinite extension 和显式零normal lift。piecewise布局包含左右面在Y的实际子面断点和上下在X/s的断点。axis对齐段上的normal多项式在XY为分片一次、RZ axial为分片一次s（至多二次r），其有符号速度积分可用现行三点Gauss/shifted moments；不需要新的Core scalar/root API。RZ径向正半径面只在固定r除r，轴零面积直接为零。这是Hdiv normal相容场，不保证跨cell切向C1；尚未接入FirstJet characteristic/Jacobian，后续耦合须单独验证实际离散map与其导数。

## 3. 规范atom与共享来源

按实际normalAxis/plane/有向tangent区间overlay primal和两个component dual的真实segment；还须切开上述kinematic piece断点，不能只叠加输入segment端点。key不含速度component或材料；primal atom用原共享trace，interior atom用原ownerPressure的B。同位置的来源必须一致；若有相矛盾的primal/owner来源，拒绝而非owner投票。每条actual segment映射到ordered signed atom indices，真实AMR粗视图仅求实际子面和。

overlay、one-ring及几何图只准备一次；只查询所选domain，边界没有隐式周期配对或ghost。局部数据/速度更新随实际face/atom线性；一次全域line排序可为O(N log N)，不得冒称全部准备严格线性。native同点的normal取该shared trace，transverse取实际邻接cell B按真实half-dual overlap的几何权重平均，不按密度平均；这仅定义运动学点view，不恢复材料quarter约束。

## 4. 同源causal phase factory：H10j

已实现 `CchPhaseRate::prepareSpatial(mesh, phase, actualSpatialAdvector)`，工厂直接比较advector的声明geometry/full clock与phase，核验原layout/D_I/mesh和momentumGeometry instance；不再另收裸速度数组或创建另一套advector。只读 `spatialAdvector()/momentumGeometry()/dualFaceRates()/spatialReport()/advectorIdentity()` 与原primal率视图绑定同一owner。

每个canonical atom从phase的原owning view调用NativeQuarticTrace，普通相区间与原两侧zero区间不混用。ordinary与速度/kinematic piece交集、zero的两侧事件与实际normal polynomial正负分片均保留短范围；两侧event没有全局ordinal对应，不能按近似root坐标相近合并。Root.parameter/lower/upper属于其原opaque chart，不是全局atom参数；terminal ranges是独立offset/width的浮点近似，不是认证bounds。使用原same-source refinement/有限工作预算和shifted积分，不扩大隔离容限或重造source；保留求积/误差估计与无法认证的拒绝。zero上

`Q_l=integral(a^+*chi_minus+a^-*chi_plus)dA`。

**u_f=0或total Q=0不意味着各相率为零**：该面trace可在内部变号。只有实际normal polynomial恒零才无需流入相。两相分别积累signed polynomial积分，再以gross核验总量，不能用total减去大相而吸收小相或零净量中的相反流动。正物理求积权重保留，速度可以有符号。内部用actual causal side；外部/selected zero流入缺少正式相数据时拒绝，已知内侧或额外scalar支撑不能被静默当作外部boundary law。

正反normal同时交换side与速度号，得到Q反号；atom单份值提供真实共享/AMR可加性。zero的空间平均、clip内部迹、Interface Green ledger均不能替代它。H10l2已删除旧evaluate，当前仅有spatial factory，无失败fallback。

## 5. 固定后续预算和退出条件

| 任务 | 唯一事实 | 固定文件/测试 |
| --- | --- | --- |
| **H10j0：PASS** | 实际compatible advector与可复用几何准备，不查phi/不生成phase F | 新CchSpatialAdvector.h/.cpp + CMakeLists.txt现lamp_hydrodynamics_cch单一源登记，共3文件；仅既有cch_momentum_geometry_contract.cpp，1 test；0新target/tools |
| **H10j：PASS** | 上述advector上的同源causal phase owner | 既有CchPhaseGeometry.h/.cpp，2 production；既有cch_phase_geometry_contract.cpp，1 test；既有研究driver |
| **H10k：PASS** | 唯一Flux/source stamp，Transport复用原F；phase/material fixtures已迁移 | 2 production/2 tests，384实际准备及8/8必要回归 |
| **H10l2：PASS** | dual消费同一F、真实R/E_sp；quarter/legacy已退役 | 5 production/2既有tests，KEP production未改 |

CMakeLists登记是H10j0唯一编译边界改动：新cpp只属于现CCH源表，同时从legacy Hydro glob排除，不新增库/target/另一实现。当前H10i2不修改这些生产文件。

H10j0只证明运动学可识别再现/同源/共享/divergence/复用；H10j证明实际causal phase积分及来源。唯一材料F由H08乘密度，H09不能另积分/另乘密度。H10l的 `R_d=Delta_m_d/dt+q_d=T*R_p+G/dt+E_sp` 只核验真实定义，不令E_sp/G/R为零。

材料仅唯一共享F、内部反号/真实AMR子面和、来源一致、KEP/skew四类硬约束。不同材料所有cut行exact affine与skew不重新列硬等式；最终native常数/可识别affine空间矩及其二阶余项仍与material KE一起验收，不能只验segment trace。原内部affine先通过，再原完整H10b空间阶>=1.7、最细L2<=0.05、Linf下降和全部KE/时间/donor gate。H11/H08a等待。

2026-10-04 H10j调用库存核对：原四份契约均有旧evaluate调用，故早期H10k立即删全部旧入口的时序改为H10k迁移材料owner和两份fixtures、H10l迁移最后消费者后删除。新phase工厂不回落历史入口，实际材料F尚未实现。

2026-10-04 H10k实施后：唯一材料F已实现，瞬时Flux与区间Transport共享同一storage，实际接口以H10k API为准；上文早期拟定路径和未实施说明只作历史。剩余旧consumer仅dual/KEP两份fixtures，H10l迁移后删除。

2026-10-04 H10l2收口：上述H10j/H10k调用库存为历史快照；当前旧evaluate及quarter算法库存为0，真实same-F dual、R/E_sp/G、8项必要回归通过。完整PDE精度不由来源迁移自动证明。
