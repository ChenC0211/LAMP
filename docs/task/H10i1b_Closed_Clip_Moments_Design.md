# H10i1b / 实施前闭合裁剪矩设计

2026-10-04；H10i1a 的 191 个冻结 payload 和 44 个当前交付已核验。实施证据位于 `runs/checks/H10i1b/20261004_100116_372`。五 production、两个既有 test 的上限和 TODO 的原 gate 保持。

Core 新 port 只接受同一 native Q4、真实 native cell 内的轴对齐矩形分区、完整正规 graphs、实际 NativeQuarticTrace cuts、XY/RZ 与共同 Green origin。当前真实 consumer 就是这类矩形；斜四边形明确拒绝，不以它冒称一般 polygon 能力。严格 FP TU 通过既有 Segment owner 的内部 bridge 访问 retained value/normal 查询，公开 terminal ranges 不参与 root 身份或重合判定。

fixed 边界按区域内部法向侧读取普通 phase 与 zero chart 的同源 normal phase。所有区间在 retained rational 参数中裁剪后才形成物理长度和 RZ 一阶矩。外边界和同 regionId 的内部 seam 分别记录；不同 regionId 的 zero seam 两侧迹可以不同，其 fixed 和等于真实 coincident interface moment，不能强制几何内部迹反号共享。

free 积分复用原 CellIntegrator 的 retained cut 分割、root ordinal、有界 refinement、正规 graph 和 controlled GL。geometry-only 模式只积分图方向的 Green 矩，另一方向用同一 retained 端点的 shifted primitive 分部积分恢复，不生成 25 列 Jacobian。区域 owner 和 enclosure 检查仍保留。恒零 cut 的实际 chart 端点参加分割；只有同 cell、exact zero support 和整个原 graph subbox 的原 zero-set normal-partial 证书共同证明唯一正规分支时，才可把该段视为 fixed/free 重合。区域外边界上重合段仅计 fixed；同 owner seam 不删除内部 free 段。

新结果分别输出 liquid/gas 两方向 closed moment、fixed/free 项、coincident 项、capacity、实际 arcs/refinement/error。真实 PhaseRegionSource 消费此 port，仍沿原严格非负/capacity/短相容限选较小相及补相；不改变 topology publication。预算耗尽、缺侧、foreign source、未认证退化均拒绝且不修改输入。

固定验收：真实 dual 内部位置八组出版来源的 16 盒与八个 momentum mass；独立平面及正非仿射 factor 闭式参考，XY/RZ 正负相、axis、AMR、反向、partition、生命周期/失败事务；误差 ≤ 1e−12 实际容量。原 native-face 八组仍记录 publication 拒绝。受影响短契约、64 普通来源、内部 affine 先验证；本轮不运行尚未满足前提的完整 H10b，不进入 H11/H08a/H10j，不生成新材料 F。
