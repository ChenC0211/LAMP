# I12a11r20 / 几何纯相与预算内库存差的发布契约

目标：只解决一个事实：在已验证的 r18/r19 有界库存输运前提下，GeometryAccuracy 允许经过原几何认证的纯相叶保留同一几何预算内的 M−V 差，发布不改变 M、不生成新几何成分。

允许修改：`include/Modules/Interface/HermiteGalsBubble.h`、`src/Modules/Interface/HermiteGalsBubble.cpp`（最多 2 production）；扩展 `tests/interface_accuracy_bubble_contract.cpp`（1 test-related）；本文、Report/API、README、库存表示设计补充及本项证据。

禁止修改：Core 相位/成分认证、Hessian 预算、共享面通量、库存权重公式、库存装配、roundoff 常数、reflux、I11 发布实现、完整时间驱动、冻结文件及原 I12 gate。

数学约束：几何纯相由同一 scalar 认证，与严格库存端点 M=0/C 区分；`abs(V-M)/C<=relativeVolumeTolerance`，库存仍严格 `[0,C]`。纯相不迭代 center，不为微量 M 建立未解析成分；超过几何预算则拒绝，不改 M 或跳过拓扑认证。InventoryRoundoff 保持原 pure roundoff 门槛。GeometryAccuracy 的 cut/pure 共用同一显式几何容限，Report 单独记录实际纯相差。

实现要求：复用原纯相分支及最终验收，只调整 GeometryAccuracy 的纯相误差契约，先在设计说明中解释离散语义和授权关联。没有新方法开关或失败后 fallback。原 I11 再独立测量几何误差；结果必须能直接供 r18/r19 消费。

验收：T0 无 warning；旧 InventoryRoundoff 输出不变。T1 保留 r17 的 14 弱 cut、16 校正、9 拒绝；其中 GeometryAccuracy 的 pure 超预算拒绝夹具应改为超过显式几何容限，保留其拒绝语义。新增 XY/RZ×uniform/AMR×正反相 8 场景，纯相相对差 4e-11；GeometryAccuracy 接受、旧契约拒绝。验证原 first jets/M/identity 保留，pure/cut center 不动、scalar jets 不变、component 数目不变、解析曲率、I11 唯一 owner、r18 整格恢复、下一次 r19 面→库存→reflux 链路严格有界并守恒。追加几何超预算、M 超界、发布容限收紧拒绝和失败隔离；未解析成分继续由原发现器拒绝。

失败规则：纯相身份、几何误差、库存或现有回归断言失败即 FAIL，保留原日志；缺失必要 API 才 BLOCKED。不得把本项当作 cut→pure 完整时间驱动验收。

交付：源码/契约说明、测试矩阵与旧契约反例、回归/耗时、冻结 hash 审计、报告和可复核 archive。完成后才接入半步/全步驱动。
