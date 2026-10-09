# I12a11r19 / 库存加权共享时空面传递 API

`integrateHermiteGalsSweptFace` 和 `integrateHermiteGalsTranslationFace` 增加接受 `const HermiteGalsInventoryRegionSource&` 的显式 overload。调用者在步骤开始时用同一旧 state 准备 r18 库存区域 owner，再为每个实际 native face 生产一次区间传递。原几何 overload 保留用于比较；没有失败后切换。

两种来源共用原扫掠矩形、坐标内外包围、相积分检查和全相入流内核。来源必须同时具有旧 state identity 与同一个 snapshot 地址；支持 owned 单轴恒定平移，XY 或 RZ 轴向。非零法向扫掠要求零切向速度、完整旧几何覆盖和可验证的特征域。入流边界仍要求绑定于原区间的显式全相历史。没有 ghost 几何外推。

库存分支消费同一次区域查询的加权 liquid/gas，沿用 smaller-phase/complement 与严格 gross 面量界限。生成的 `HermiteGalsIntegratedFaceTransfer` 直接交给原 `assembleHermiteGalsIntegratedInventory`，再由原 reflux 消费同一缓存传递；原 M、舍入修复预算及库存装配均未改变。

新增诊断如下：

| 字段 | 含义 |
| --- | --- |
| `TranslationFaceResult.inventoryWeighted` | 显式 producer 类型，包含 prescribed inflow 和零速度结果 |
| `SweptFaceReport.inventoryWeighted` | 当前扫掠是否消费库存区域来源 |
| `geometricLiquidVolume` | 同次查询的原几何液量参考，非负，无第二次几何积分 |
| `signedInventoryCorrection` | 库存液量相对几何参考的差，带全局 +x/+y 面方向 |
| `l1InventoryCorrection` | 同次区域查询的 L1 修正估计 |

原 `innerLiquid/outerLiquid/innerGas/outerGas` 对库存来源记录库存相量；`momentErrorEstimate` 仍是原几何求积估计。区间切分复用相同旧 M/geometry、速度、map 和历史。诊断报告不控制数值分支。

Transfer 本身继续是纯代数契约，不是任意外部求积的认证对象；其构造器不能证明调用者使用了哪种积分器。生产调用者必须明确选用库存 overload，并将该唯一结果用于库存及 reflux。此项没有切换完整时间驱动，也未调整 I10/I11 的纯相发布语义。

参见 [TODO](../task/I12a11r19_Bounded_Inventory_Sweep_TODO.md)、[报告](../task/I12a11r19_Bounded_Inventory_Sweep_Report.md) 和[库存表示设计](../task/I12_Bounded_Inventory_Transport_Design.md)。
