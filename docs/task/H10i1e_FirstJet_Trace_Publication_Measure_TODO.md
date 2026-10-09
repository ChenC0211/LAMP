# H10i1e / FirstJet 同源 trace 派生出版测量

状态：**PASS / engineering PASS**；[报告](H10i1e_FirstJet_Trace_Publication_Measure_Report.md)。H10i1d 的 189 个冻结 payload 与 69 个当前交付在实施前已核验。上一轮为 progress；本轮仅关闭 H10i1e。实际 16 组发布/352 个单元、原 native 8 组/40 side graphs/80 zero halves 通过，17/17 回归、原 30 CSV/affine/82 flux 数值保持；原 native H07 8 组仍拒绝，下一唯一任务 [H10i1f](H10i1f_Native_Trace_Phase_Consumers_TODO.md)。[总契约](H10i1_Coincident_Phase_Trace_Clip_Contract_Design.md)、[Core拓扑](../api/H10i1c_Native_Zero_Edge_Topology.md)、[半面cache](../api/H10i1d_Readonly_HalfFace_Trace_Cache.md)、[实施前设计](H10i1e_FirstJet_Trace_Publication_Measure_Design.md)、[具体身份/边界核对](H10i1e_Trace_Publication_Scope_Note.md)。

目标：让正规 native zero-edge 的 Pure FirstJet 出版消费同一 immutable source 的 typed topology/half cache及派生几何测量，闭合完整几何/边界账本；不引入 primary 库存或改变 phi/grad/center/Q4 来匹配体积。

允许修改：`include/Modules/Interface/InterfaceCellMeasure.h`及对应`.cpp`、`include/Modules/Interface/HermiteGalsPhaseVolume.h`及对应`.cpp`、`src/Modules/Interface/HermiteGalsPublication.cpp`（5 production），既有 `tests/contracts/interface/interface_hermite_gals_publication_contract.cpp`、`tests/contracts/interface/interface_phase_region_contract.cpp`（2 tests）；既有 research/docs 可更新。0新工具/target；超过预算先拆任务。实施前核对纠正：原 native-side 8 组拒绝事务实际位于 publication 测试，故用它替代初稿误列的 first_jet 测试；后者四种回归保持必跑。开始前核验 H10i1d 平面 manifest、保存实施前设计/目录/原 gate。

禁止修改：Core Q4/root/regularity/trace引擎、Interface primary/transport/clock/trace中心预测、Hydro phase/material/convection/压力/散度，App/字典/PDF及原容限。禁止 inventory/bubble 闭合、逐cell几何体积拟合、质量修复、source perturbation、half相分配、AMR特例、失败后fallback或新的 scalar/root 来源。历史 Conservative/敏感度入口不自动迁移。

数学约束：typed topology、half traces、raw frame/sample、派生 volume 均使用出版事务中唯一原 Snapshot/revision/time。regular zero-edge incident leaf 的内部相使用原完整证书，不凭沿边法向采样排除内部零分支。canonical geometry owner 与实际液侧可不同；两侧 availability 和无内部分支的 interior phase 语义保持。真实外部/轴/声明边界的固定与free贡献不能重复或遗漏；negativeRegion不存在时不得制造外部液体。

实现要求：先核对现有 measurement/volume/publication 事务，为 FirstJet 显式接入 `discoverNativeTraceVolumes`（最终签名在设计确定）、Core complete trace topology 与 typed cache。纯 interior-phase leaf 即使有 side-attached boundary graph 也按同源相证书测量；一般内部 regular arcs保留原精确根/controlled积分，不把attached graph送到尚不支持它的sensitivity或region consumer。两方向Green账本、capacity/volume诊断和端点/采样完整性显式记录。publication失败保留原candidate/accepted状态/identity/clock；成功只提交原无M FirstJet数据与其派生来源。

T0：标准Release、architecture/self-test/layout及原FirstJet/publication/phase-region/half-cache/flux/H07/momentum/material/dual/KEP必要回归；H10i1b clip、c topology、d half cache、64普通source/原内部affine保持。默认legacy路径与冻结历史PASS/FAIL不改判。

T1：原8个native-face场景的实际FirstJet publication（XY/RZ×uniform/真实2:1×正负相）通过；Core来源、raw frame、half-cache/source/time/revision一致，正非仿射Q4及ordinary→side/多component/axis/selected/真实边界、wrong source、生命周期/失败事务分别核验。独立平面和几何闭式体积/两方向Green贡献 ≤1e−12实际capacity，严格非负/capacity/regularity/FP与原预算保持；actual shared/subedge几何加回正确。一般临界/全零/未认证退化仍拒绝。

H10i1e 不提前把 H07 phase consumer、RegionSource或真实 native-face dual mass全链算作通过：这些消费尚归H10i1f。本任务若必要链需要额外Core/H07修改，保存具体反例并先拆分，而非绕过source或门槛。

失败规则：任一同源/正则性/边界账本/非负容量/事务/必要回归失败，冻结反例后修正或拆任务；不扩大容限或隐藏未支持场景。

交付：API/实施前设计/报告、实际出版和独立参考、标准短契约、实际编译源码/二进制/退出/CSV、完整source/input范围核验与平面SHA256 manifest；单独登记下一H10i1f。来源链闭合后才H10i复核/H10j；最终native内部affine先通过，再原H10b，H10未关闭，H11/H08a继续等待。

材料仅四硬约束：唯一共享F、内部反号/真实AMR子面和、来源一致、KEP/skew；quarter/压力到dual精确交换与R=0仅诊断，继续真实R/uR。Pure FirstJet几何完全独立。
