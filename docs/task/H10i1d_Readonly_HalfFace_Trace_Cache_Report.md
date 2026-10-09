# H10i1d / 同源只读半面 trace cache 验收

2026-10-04，**PASS / engineering PASS**；本轮只关闭H10i1d，goal保持active/progress。H10、原完整H10b和完整流体/EHD未关闭。下一唯一任务 [H10i1e FirstJet派生出版测量](H10i1e_FirstJet_Trace_Publication_Measure_TODO.md)，本轮未启动它。[TODO](H10i1d_Readonly_HalfFace_Trace_Cache_TODO.md)、[实施前设计](H10i1d_Readonly_HalfFace_Trace_Cache_Design.md)、[API](../api/H10i1d_Readonly_HalfFace_Trace_Cache.md)。

前置独立核验H10i1c的 **175 archived payload / 61当前交付**，保存实施前设计及1,042 source/346 input目录快照；先用标准Release重建并运行原flux短契约保存对照。检查目录 `runs/checks/H10i1d/20261004_111021_237`，冻结目录 `runs/baselines/H10i1d_20261004`；包含编译源码、实际研究二进制/依赖、各模式真实退出记录、原始CSV、标准回归日志和独立验证脚本。

`prepareNativeTrace`为现有move-only cache新增显式native Q4来源。每个真实row/half唯一owns原Core NativeQuarticTrace；原顶点/face.center、normalAxis和selected活动校验保持。普通`halfFace`直接返回同一typed查询的legacy view，zero support明确拒绝。默认constructor继续legacy恒零拒绝，未被自动换到新来源。first-active caller FP检查、原OpenMP异常收集/顺序及未完成cache不发布的事务保持；未复制Q4/root引擎。

typed count明确分为value incidences和全部value/normal events。`refineTraceEvent`用原immutable same-root API，成功替换唯一query，依赖测度与typed/ordinary引用失效并须重算；失败保留query/source/count。typed ordinary `refineRoot`检查当前legacy root与对应value event后委托同一细化，不保留两份不同enclosure。

首轮实际value根暴露wrapper错把opaque Root chart ID当作稳定身份。原Core细化会追加同根child chart，公开parameter/lower/upper也属于新chart，不能跨chart比较。已检查原引擎并保存[身份核对](H10i1d_Refinement_Chart_Identity_Note.md)，wrapper保留source/native cell/event index/order/side/coincidentChart与incidence核验，移除错误chart相等要求；同根认证仍由原exact parent bracket完成，没有放宽门槛。

| 验收项 | 实际结果 |
| --- | --- |
| 原16场景的真实half-cache来源 | **16/16、1,728 actual halves**；**80 zero halves / 576 value events / 0 normal events**（正规平面）；独立XY/RZ两侧/普通相物理面积最大误差 **0** |
| 既有flux短契约扩展 | **16 layouts / 1,728 halves / 160 zero halves / 64 actual hanging halves**；两种normal、XY/RZ、正负相与正非仿射Q4，physical measure/反向/完整实际子边加回通过，原1e−12 actual-scale gate保持 |
| refinement | **3次成功**：独立非二进制normal根、ordinary value根与同一value的重复细化；child chart/唯一legacy view一致，重新计算两侧相测度与独立1/3参考通过 |
| source/事务 | wrong family/domain/activity/source、moved-from/lifetime/禁copy、invalid/solid slot、改网格identity、strict FP、invalid event及有效exact-point事件的晚期失败保留query/count均通过 |
| 几何边界/来源种类 | **64 physical exterior halves / 8 selected zero halves / 36 solid halves**；RZ轴normalOrder=2和2/3/4阶Taylor保留；ordinary→zero在相邻真实half/subedge拼接，未制造非实际合并face；全零仍拒绝 |
| OpenMP cache | **496 actual halves / 248 faces**，跨原`>128`并行分支；`OMP_NUM_THREADS=2`真实运行通过，与单线程短契约数值相同，标准target启用`/openmp:llvm` |
| 原flux契约 | 实施前标准构建的 **82行原数值输出逐条相同**；新增cache观察在原流程之后，不改历史库存/flux参考 |
| 原来源链 | c的 **7份trace/clip/mass/topology/support CSV逐字节相同**：dual内部8组/16 clip/432 mass、Core16组/native8组/40支持及旧native publication8组拒绝保持 |
| 原普通来源/材料诊断 | **64 layouts / 157,810 queries**，**22/22 CSV逐字节相同**，原独立物理/reference门槛保持 |
| 最终native内部affine | **3,760行逐字节相同**；内部最大 **1.5849193387254136e−13**，原AMR dof911反例 **7.41629e−14** |
| 标准Release受影响CTest | **17/17 PASS**（原16加既有flux）；原KEP SUMMARY/NATIVE_MOMENTS、CAP、NATIVE_TRACE、COINCIDENT_CLIP/MASS以及cside/FirstJet失败事务保持 |

RZ水平half的物理面积分别为各自 `pi*(r_hi²-r_lo²)`，不等同于完整face.area/2；只在真实竖直恒r面才按等长各半得到等面积。fixture最初的边界/point-root判据也按原API核对：pointOnly不是所有exact endpoint root的标志，精确enclosure的晚期细化拒绝另行验证。preflight日志全部保存；修正错误fixture/元数据断言，没有改变生产几何或数值容限。

范围 **2 production / 1既有test / 1既有research driver**，0新工具/target。只修改InterfacePhaseGeometry public header/cpp和既有flux短契约，research只追加被动half缓存观察；原数值函数和c旧Core/source观察去除新增观察后源码逐字符相同（仅规范化换行）。source catalog **1,042：1,038不变、4变更**，input **346：344不变，仅driver/README变更**，成员保持。Core原Q4/root、Interface primary/transport/clock/volume/publication、Hydro F/R/uR/KEP、App/字典/PDF保持。

实际final driver SHA256 `1F80056BF30F6D2AF78D9DFE9D7ABF127F6C8451DC9C92671432C367A8974EA9`，research exe `CB604FE96BB209728028C515FDB9740AB640966A2AD2AF187367172A14DD4088`；四个最终研究模式对应该二进制，双线程短契约对应最终flux目标hash。`half_trace_cache.csv`/summary清楚标明独立FirstJet seed identity，不是accepted publication。平面manifest核验最终archived payload/current canonical delivery与全部175前置archived payload。

限制：typed cache支持原Core normal Taylor/事件，不意味着临界界面可以出版。新factory尚未接入实际FirstJet派生volume/publication/H07消费者，native出版8组仍拒绝；H10i1e/f逐项迁移来源链后才H10i复核/H10j。没有新材料率、convection factory或时空积分，原完整H10b未重跑、14项历史失败保持。内部affine先通过再原H10b；H11/H08a等待最终空间通过。

材料硬约束仍只有 **唯一共享F、内部反号/真实AMR子面和、来源一致、KEP/skew**。quarter/压力到dual精确交换及R=0仅诊断，继续真实R/uR；Pure FirstJet完全独立，无primary库存/bubble/几何体积拟合/质量修复。
