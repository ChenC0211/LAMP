# H10i1e / 实施前设计

日期：2026-10-04。实施前 H10i1d 189 payload、69 canonical delivery 与 c 的 175 archived payload 核验通过。原文件/目录目录表、允许变更文件、TODO 与保护 App/PDF hash 保存到本任务 `runs/checks`；本设计在数值修改之前保存副本，后续不可追改。

唯一事实：Pure FirstJet 的实际发布可以消费既有 complete native trace topology、实际 half cache 与原 Q4 的派生体积，primary/闭包值/clock 不改变。普通敏感度、ConservativeVolume、H07/RegionSource 仍采用原入口；本任务不宣布完整 Hydro native 来源消费链通过。

公开测量入口：`discoverNativeTraceVolumes(mesh, domain, owningSnapshot, contacts, symmetry)` 负责原 Core trace discovery 和 typed half preparation；另有 `measureNativeTraceVolumes(mesh, domain, partitions, topology)` 复用发布对象的同一只读 cache/topology。后者验证原 Snapshot 地址、布局、chain 与 side support，不接受用户 phase mask。结果增加两个 Green 方向的 fixed/free/volume/capacity 及闭合诊断，不暴露敏感度。

零边账本：fixed phase 仅包括原普通负相区间；正长度零迹不被赋予 half phase。完整拓扑的 side graph 自原 `nativeMoments` 积分，并且仅给其真实存在的液侧 leaf 添加 free moment，canonical graph owner 可以在气侧。非零普通内部分支仍给其原 leaf，沿用原 controlled/root/local integration。零边 incident leaf 的完整 signed-normal 证书给出纯内部相，体积用真实 CV capacity 或 0；fixed + free 的两方向几何账本必须与证书一致。对外部/轴/selected 支撑亦执行同样的 availability 规则；不存在液侧则不得制造外部液体。shared fixed face 始终一次积分、两侧反号，side free 只计一次；不能把 side graph 送到旧敏感度或尚未支持它的 region 消费。

普通路径保真：没有零支撑的 typed cache 仍返回原 `ordinary()`，轴向 fixed 复用原 boundaryMoments；原 ordinary geometry-only native evaluator 的 root/GL8/local integration 保持。只有显式 trace 测量新增 radial 诊断与 side ledger。原 1e-10 内部 ledger/estimate gate 与配置容限不放宽；独立平面验收仍 1e-12 actual capacity。局部 strict 非负/capacity 不裁剪，临界/全零/不完整认证保持拒绝。

发布事务：最终 Data 按 state→domain→topology→partitions→volumes 初始化，domain 到达 final allocation 后才绑定 cache。FirstJet volume 显式从 Data 的相同 topology/cache 得到；普通独立 volume 测量按 primary transport 显式分派。所有 sample/quadrature/half/source/measure 检查在最后 consumePublication 前完成。任意失败保留原 candidate、accepted source、identity、primary、clock。成功返回原 FirstJet 数据与只读派生值；无库存/bubble、体积拟合、phi/grad/center 修正。

预算核对：5 production；2 tests 使用 publication 契约和 phase-region 契约。原 8 native-side 拒绝测试实际在 publication 文件，初稿误列 first_jet 测试已在实施前纠正；FirstJet 四种回归仍必跑。既有 research 仅补充来源/发布的被动证据，0 新工具/target。原普通 64 layout、内部 affine、clip、topology、half-cache 必跑；旧 legacy 拒绝仍可单独验证，不把旧失败历史改写为新成功。

完成只关闭 H10i1e，登记独立 H10i1f。之后来源链复核/H10j、最终 native 内部 affine 与原完整 H10b；H11/H08a 继续等待。材料仅唯一共享 F、内部反号/真实 AMR 子面和、来源一致、KEP/skew 四类硬约束；quarter、精确交换与 R=0 仅诊断，真实 R/uR 保持。
