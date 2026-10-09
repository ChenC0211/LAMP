# I12a11r20 / 几何纯相库存差发布报告

**状态：PASS（2026-10-02）。I12 / Gate I 仍 FAIL。**

依据 [TODO](I12a11r20_Pure_Geometry_Inventory_TODO.md) 与[设计补充](I12_Bounded_Inventory_Transport_Design.md)，将已授权的几何—库存分离扩展至经原 scalar 认证的纯相叶。API 见[纯相契约说明](../api/I12a11r20_Pure_Geometry_Inventory_API.md)。本项不切换完整时间驱动。

## 改动

2 个 production 文件为 `HermiteGalsBubble.h/.cpp`：GeometryAccuracy 的纯相误差门槛采用同一 `relativeVolumeTolerance`，默认 1e-10；预测纯相 scalar/center 保留，误差超标拒绝。InventoryRoundoff 继续使用原 `min(relativeVolumeTolerance,64*epsilon)`。没有新开关、M 重置或隐藏成分生成。原拓扑发现、cut 求根/Hessian 预算、库存权重、面通量、舍入修复和 I11 实现均不变。

1 个测试文件 `interface_accuracy_bubble_contract.cpp`：保留 r17 的 14 弱 cut、16 有界修正和 9 类拒绝。其中纯相拒绝夹具的库存差由 128 epsilon 改为显式几何门槛的两倍，保持“超过所选契约即拒绝”的目的；原弱 cut/修正数值输出未改变。

## 新验证

XY/RZ × uniform/AMR × 正反相共 8 场景、126 个认证纯相叶。每个纯相叶保留相对容量 4e-11 的库存差，曲线 cut 保留原测得 M，所有 center 的额外 Hessian 预算为零。

- GeometryAccuracy 全部接受，原 InventoryRoundoff 对同样 8 场景全部拒绝。
- first jets、mixed、center 和原 M 保留；发布前后在每叶 3×3 点比较 scalar 值、一二阶导数，全部逐值相等；component 数仍为 1。
- I11 保持 identity、唯一 snapshot、cut/volume owner 和独立 M。曲率与解析圆/球值比较，最大绝对误差 2.75335310107039e-14。
- 发布后的同一 state 直接进入 r18→r19→原库存装配→reflux。整格库存恢复差在本矩阵为 0，下一次库存更新严格有界，边界 ledger / 总容量最大 2.22044604925031e-16，粗细液量/总量残差为 0，未触发舍入面修复。
- 新增共 11 次拒绝：上述 8 个旧契约反例，加上收紧几何容限、超容量 M、收紧 I11 容限。拒绝后候选的 M 和 identity 不变，I11 失败不消费候选。旧 bubble 回归继续覆盖隐藏成分拒绝。

纯相最大相对差 4.00000949553956e-11；I11 独立 raw 积分最大相对差 4.00009174018587e-11，均低于同一默认几何门槛。后者不是库存守恒残差。几何中不存在的微量库存作为有界相含量误差记录，未被解释为新物理微滴。

首次构建和直接测试通过；直接测试 21.5003621 s。33 个目标完整构建通过，42.0890457 s，无 warning；12 项 CTest 全通过，235.9855832 s。其余 11 项完整输出（含 architecture）与 r19 逐字相同；扩展测试中的原 r17 输出前缀也逐字相同。旧 `--crossing-probe` 保持相同输出及 exit 1。

81 冻结文件、3 方法 PDF 及 r14/r19 archive hash 不变。证据目录：[I12a11r20_20261002](../../runs/baselines/I12a11r20_20261002/)。`verification.json`、`pure_comparison.json`、`regression_comparison.json`、`source_audit.json` 记录验收，`delivery_manifest.json` 与 `archive_verification.json` 给出 `evidence.zip` 的逐文件完整性核验。

## 限制

本矩阵检验“认证纯相误差保留→I11 发布→下一次共享通量”，未实施下一次 first-jet 更新或几何发布，不能视作完整 cut→pure 时间步测试。后续将 r17/r20 与 r18/r19 共同接入原半步/全步驱动，使用预先固定的导数预算并重跑原 I12 h/dt/曲率矩阵。一般 rotation/vortex 时空积分与完整收敛仍未通过，**I12 / Gate I 仍 FAIL**。
