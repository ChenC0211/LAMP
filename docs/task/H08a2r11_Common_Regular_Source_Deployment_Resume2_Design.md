# H08a2r11 第二次恢复 / 实施前部署设计

2026-10-05。只恢复 [r11 同步部署任务](H08a2r11_Common_Regular_Source_Deployment_TODO.md)。上一完成轮 r17 为 progress：102 个冻结 payload、11 份交付和 1383 源文件/输入已核验；原 r11 两次失败记录保持。r14 常数/低阶准备、r15 真实材料动能消费者、r16 可信连续流参考和 r17 跨片 point/J 修复均已闭合。

预算为四个生产文件 `CchSpatialAdvector.h/.cpp`、`CchCharacteristic.h/.cpp`，两个既有测试 `cch_momentum_geometry_contract.cpp`、`cch_characteristic_contract.cpp`。生产四文件采用 r17 冻结候选，只更新 header 的积分公式注释；不另加算法。momentum 只迁移原共享准备身份检查；characteristic 采用 r17 全九组检查，清理归档中不再调用的未解析 RK4 参考，保留完整原输入、阈值及可信 point/J 常量与预算。禁止顺手调整材料、最终 native 算子或研究输入。

owner/data-flow：prepare 持有唯一 regular 权重与真实 support knots/primal-dual overlay；每次 evaluate 只 apply 一次。实际值/J、native 物理中心值、cell/piece 和完整 normal 多项式均查询同一系数；cellDivergence 保留原 D_I 平均量。旧 CchAdvectorPatch、tangent 拟合、faceB/cellLine/lift 实现退出生产；只剩 momentum 一处旧存储身份消费者需迁移。既有研究观察器已在 r13 迁移，无需修改 driver。

运动采用五点完整 Q4 值/J 迹证书；内部时间槽为 r17 三阶分片 Picard，同一求导包括积分节点/宽度/连接时刻。真实事件仍按三次 dense 处理；舍入带内点/J共同保留原公式。公开注释同步说明，禁止只更换运动场而保持材料 atom 来源不变。原透明连接、constant/tiny/axis/clock/source 身份拒绝全部重验。

先标准构建并通过原九项 CTest；再运行实际部署程序的 advector/closure/boundary 原 64 布局、208 basis、32 边界布局及独立完整次数物理参考。原 regular 144 结构组、低阶原矩阵及算术拒绝仍由完整 momentum/phase 等契约覆盖。point-only 负控将基于当前部署代码/库另行构建隔离变异，确认原可信 J gate 拒绝；不写生产开关。

短契约通过后，依次核验 audit/audit32 原内部 affine、rows/rows32 的最终 native 独立矩阵；全部通过才由原通用 run_case 入口运行 full 原 28 窗口/84 gate。runtime、原 inputs、source identity、最终窗口数和完整 gate 表必须逐项核验；旧四 donor 失败单列。若任何原门槛失败，按 r11 规则保存候选并退役，不改输入或容差；若全部通过，接受统一来源并立即返回 H08a 真正同源时空相通量积分，将必要的 FirstJet 运动消费接入并入该任务，H11 未启动。

**基础设施收口（2026-10-05 用户授权）：** r11 是 common regular source/characteristic 的最后一道关口。除非出现可复现且直接阻断 H08a 的数学反例，不再拆 source/Jacobian/低秩精修任务，不新增硬约束、特殊 case，也不放宽原门槛。r11 通过后立即回到 H08a 的真正同源时空相通量积分；必要的 FirstJet 运动消费接入属于该积分任务，不再单列桥接先决任务。目标仅为 R 成为可收敛的小离散余项，保持 Pure FirstJet 独立、无 M_K/bubble、无 R=0 或 quarter commuting 强制。材料硬约束仍仅为唯一共享材料通量、内部面反号共享/AMR 子面守恒、来源一致和 KEP/skew。
