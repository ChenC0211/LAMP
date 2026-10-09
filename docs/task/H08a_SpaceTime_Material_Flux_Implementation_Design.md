# H08a / 同源时空相积分实施前设计

2026-10-05 实施前登记；当时 H08a 尚未验收。上一轮 r11 PASS 的 196 payload、19 份交付、1383 源码/输入已逐一复核。本设计不另立 source/Jacobian/低秩任务。2026-10-06 最终状态为规定范围 **PASS**，见[报告](H08a_SpaceTime_Material_Flux_Report.md)；下述原门槛与失败记录保留。

## 直接消费反例与修正边界

任务开始时 r11 冻结库上的 `runs/checks/H08a/20261005_232800/witness.log` 复现 XY uniform n8、实际 FirstJet 主顶点 (0.5,0.5)、t=0.25→0、maxStep=0.0125。零场、u=(1,1)、u=(0,1) 三个来源均 480/480 连接通过完整 C1 证书，却分别因多分片起点/同时穿面拒绝。解析结果分别为到达点自身 (0.5,0.5)、(0.25,0.25)、(0.5,0.25)，J 均为 I，轨迹完全在域内。因此这是直接阻断实际 FirstJet 消费的反例，不是提高 source 精度的新目标。

仅对已有完整 C1 证书连通的局部骨架星形邻域，允许等价分片选择和同时连接穿越；跨真实物理边界、未证连续连接和不连通点接触不因此允许。C1 保证值/J 唯一，零宽过渡的连接积分一阶边界项相消，不需要任意单侧混合或第二套 J。保持原三阶公式、数值容差、时间矩阵、可信 point/J 和全部证书负控；原合法常数运动的“应拒绝”能力检查改为原输入解析通过，其数字精度门槛不放宽。若实现无法保持原精度，仍为 FAIL。

## 同一积分的数据流

### 直接消费反例补充：端点的局部 Green 算术

2026-10-06，新增的原定区间起止 crossing 检查在 RZ 对齐平面复现旧测量器的严格容量拒绝；证据为本检查目录 `ctest_endpoints.log`、隔离的 `endpoint_probe/run.log` 与 `run_core.log`。实际 cell 6 的 liquid=0.24543692606170261494，capacity=0.24543692606170258719，差一个浮点单位；gas=0，原 Green ledger 相对差 1.13e-16。另一个方向直接相减产生 −1.39e-17 的小相。这个问题在 publication 测量阶段直接阻断 H08a，不是 common source/Jacobian 的精度反例，也不能用 source 精修、体积裁剪或提高黏度处理。

隔离试验保持同一 scalar/graphs/roots/regions/求积规则及全部原断言，只调整既有 Green 势的自由原点。先以同一 scalar 的 box-center FirstJet 作一次法向坐标估计，再选取离估计最近的既有 box 角点；梯度退化时取 lower 角点。角点选择避免在固定边的整容量和中引入额外舍入的内部偏移。这个坐标只进入 Green 势，不写回 phi、梯度、根、几何或曲率，没有第二几何。原点平移的数学恒等性不依赖此估计成为真根。

此前直接采用连续内部估计点的候选虽通过 8 个端点，却使原 phase-region 严格容量检查失败，故未被接受；失败留在 `ctest_gauge.log`。角点候选的隔离复核 `endpoint_probe/run_corner_endpoints.log`、`run_corner_regions.log` 同时通过全部 8 个端点（通量误差最大 1.77e-15）及原 213 场景/123 拒绝。正式部署仍须完整受影响回归，不能以隔离试验代替验收。

为关闭这个直接消费缺口，文件范围增加一个现有内部算术头 `include/Modules/Interface/Detail/HermiteGalsRegionMeasure.h`，总计十个 production；不另拆任务。部署前后必须复验原 phase-region/shifted-Green、FirstJet publication/transfer、H07/H09 和本次端点原失败。保留 0≤V≤capacity 的严格检查、原 ledger/误差预算及 tiny-phase 门槛；若任何原门槛失败，仍不能判 H08a PASS。完整曲面研究的旧运行保留为改动前证据，新算术影响的结果须重新验证，不能直接改标。

同一原始 Pure FirstJet 和完整 endpoint characteristic source 定义整个步内几何路径。Interface 新增中立、携带源 identity/源 clock/目标 clock 的已解析运动消费入口，规定速度 RK2 入口仍使用原公式；Hydro 传入实际 CCH 运动及同一映射 J。共用 FirstJet 值/梯度与 center 更新、mixed closure 和唯一 publication，不新增几何 owner、库存或 bubble，也不让 Interface 认识 Hydro。

任意采样 t 的几何从原始 FirstJet 直接生成，不从另一个采样层递推或在三个独立 snapshot 之间插值。中间积分节点使用明确的局部候选 [t0,t] 时钟；完整父 slab 和实际 endpoint source 一直由 Hydro 区间 owner 保留。t0、完整 midpoint、end 的质量仍测量唯一几何。

边界消费显式区分原 RK2 的入流命中历史与给定整条外边界 FirstJet。H08a 若传入 boundary provider，其含义是每个采样时刻的完整外边界数据（含切向/出流边与角点）；无 provider 时仍由实际运动决定是否需要历史。触边平面的对称性由该明确标注的制造边界数据满足，不能把未施加边界条件的微小非零梯度强行清零，也不改变几何发布的精确对称性检查。内部轨迹若出域而无完整历史仍拒绝。

同一 actual B 的步内场由两端已冻结系数直接线性组合，atom 双端完整多项式也组合原系数；不重新拟合或把重算 native 数组冒充相同实际来源。每个时刻调用既有同源空间相积分，时间上用正权 Gauss 3/5 对及细分误差估计，按真实 canonical atoms 积分液气两相。该估计是数值误差指标，不宣称严格区间包络。显式初始分段、最大深度/采样量与失败拒绝确保工作有界；不足精度时不回退中点。

唯一积分 `F_atom=rho_l Q_l+rho_g Q_g` 只密度加权一次；primal/dual 均复用这份 atom 积分和既有有序 incidence。对流平均率取 F/dt，midpoint 即时率保留独立诊断语义。R 为实际几何质量差加唯一通量积分，dual 继续真实 E_sp/G、uR、材料动能及 KEP/skew，不强制为零或恢复 quarter commuting。

## 文件范围

本任务是一个贯穿现有 owner 的积分事实，按用户要求不再机械拆成基础设施先决任务。实施前将默认五文件预算明确扩为九 production：`CchCharacteristic.cpp`（上述直接消费反例）、`CchSpatialAdvector.h/.cpp`（同一系数时间消费）、`HermiteGalsFirstJetStep.h/.cpp`（中立运动消费）、`CchMaterialTransport.h/.cpp`（同源路径、唯一积分及真实 R）、`CchDualMaterial.h/.cpp`（共享积分消费）。不动 Interface primary/几何发布算法、regular 准备权重、最终 KEP 算法、压力或黏性。

先扩展原 characteristic/material 两份短测试。新增材料检查按 XY/RZ×uniform/AMR 四布局在 `tests/CMakeLists.txt` 独立注册，保留每组原 180 秒超时及全部三密度/数值断言；旧材料契约独立运行，不把同一 CTest 总耗时的超时误作数值通过。完整收敛研究放在 `cases/benchmarks/cch_spacetime_material/driver.cpp`，注册同目录 benchmark.json/README 与现有 benchmarks/CMakeLists；该驱动只调用生产积分入口，输出独立解析通量、质量与 R 的空间/时间研究，不增加生产问题分支。临时定位只归本次 runs/checks。不得将仅局部消费或单一平面通过改称 H08a 完成。

验收覆盖复核补充：原平面位置 0.48137 的穿面主要在细网格内部，不能单独证明真实粗细边界的非零相通量切换。保留全部原场景，另在同一测试入口取 0.23137，使相面在原时间窗内穿过 0.25 的实际粗细边界；XY/RZ、三密度仍沿用原 1e-7 通量预算与 1e-12 材料/KE 身份。检查实际两侧 level 不同、面相从全气变全液以及共享材料积分，不以切向零通量充当该覆盖。这是原 AMR 子面要求的证据补齐，未增加生产问题分支、结构硬约束或 source 精修任务。

## 预设验收

T0：标准 Release、架构/目录、原九项 Hydro 契约及受影响 FirstJet step/vertex/edge 契约；源码/输入变化符合本设计。

T1：原三个骨架反例 point/J 按 1e-12 解析误差检查，XY/RZ uniform/2:1 的真实主顶点/edge/center 均纳入消费；保留原 24 时间族、288 FD、2304 分区和两条 1e-11 可信参考。wrong source/full clock/regrid、预算耗尽、失败不改变原 FirstJet 和不可支持边界明确拒绝。同源采样场对两端公开值/J 的线性组合按原物理工作尺度 1e-12 检查。共享 atom/primal/dual、质量差/R/uR/KEP 身份沿用既有 1e-12 门槛。

T2：真实移动/穿面平面与曲面，XY/RZ、uniform/2:1，rho 比 1/1e3/1e6。平面解析时空通量按配置积分预算 1e-7 检查，事件落在区间内部/端点均保留；对比原中点误差但不以中点为真值。空间 n=4/8/16、独立时间加密及固定物理观察窗口都报告缺陷率和步积分缺陷；曲面相对连续参考、几何投影余项与时间参考自收敛分开。超舍入平台的误差必须随加密下降，细化后误差相对可解析物理输运尺度不超过 1%；不把严格 R=0 或 quarter commuting 变为门槛。原 H10b-v2 28窗口/84项最后复验，旧四 donor 失败保持。

最终状态（2026-10-06）：真正积分入口、24 项受影响短检查、当前角点算术的完整 28 个空间/时间窗口及原 H10b 28 窗口/84 项均通过；独立复核见 `runs/checks/H08a/20261005_232800/current_complete_verification.json`。原始 Green 算术上的 28 窗口仅作比较。新增重网格测试曾误构造非 2:1 目标，被原规则拒绝；改用原 FirstJet 契约中的平衡嵌套布局后四项通过，生产算法和门槛未变。历史失败、当前程序和来源绑定均保留。PASS 限于本任务积分、材料一致性与小余项收敛范围；同网格 8 步时间参考不证明完整流体二阶，完整 App、一般 sharp-cut PDE 和效率未验收。H11 本轮未启动。
