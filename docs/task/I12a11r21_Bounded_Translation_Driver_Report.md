# I12a11r21 / 有界库存平移驱动报告

**状态：PASS（2026-10-02）。I12 / Gate I 仍 FAIL。**

依据 [TODO](I12a11r21_Bounded_Translation_Driver_TODO.md)，把 r17/r20 几何精度闭合与 r18/r19 库存通量共同接入原 I12 平移驱动。固定配置、命令和过程字段见[驱动说明](../api/I12a11r21_Bounded_Translation_Driver.md)。本项不宣告 I12 通过。

## 接线

0 production、2 test-related 文件：原 benchmark 和 CTest 注册。复用原 `initial/predict/finish/advanceTranslation`；显式 bounded 契约在旧 accepted state 上准备一个库存区域 owner，半步和全步共享它及原速度、实际 B/J 来源与入流历史。每个阶段的 M 经原共享面装配/reflux 更新，first jets/center 经原 characteristic 输运，随后原 I10/I11 成功才更新 accepted handle。

几何容限固定为原默认 1e-10。单位方域无量纲测试的 Hessian 分量校正预算事先取每叶 `min(hx,hy)^2`，系数 1；没有从真解或失败结果调节。旧契约保留用于同函数 A/B，无失败 fallback。新增 producer 类型检查和逐阶段诊断，Report 不控制算法。

## 完整 Q2 平移

四组静止各 4 步，四组平移各 40 步至 t=.4，共 176 个完整步骤通过；每一步完成 midpoint 和 end 的真实 first-jet/M 更新及发布。AMR 初始界面均同时切过 coarse/fine 叶。

| 40 步场景 | position | gradient | curvature | cell M 解析差 / 容量 |
| --- | --- | --- | --- | --- |
| XY uniform | 2.00863e-16 | 1.68628e-15 | 7.90479e-14 | 2.44249e-15 |
| XY AMR | 1.09561e-16 | 2.41809e-15 | 1.40332e-13 | 4.10783e-15 |
| RZ uniform | 1.79897e-16 | 1.92935e-15 | 5.68434e-14 | 8.48148e-16 |
| RZ AMR | 1.79897e-16 | 1.92935e-15 | 6.12843e-14 | 2.92140e-15 |

四个静止场景的最大 first-jet/M drift 为 2.22044604925031e-16。全矩阵最大库存 ledger 1.39417183748739e-17，reflux 为 0，严格库存容量检查保持。包括初始化共 360 条闭合记录、352 条区间面记录，未发生拒绝或 M 重置。

Q2 全程几何差最大 4.88498130835069e-15，所有预测已满足请求几何容限，bubble 修正数及新增 Hessian 扰动均为零。单面库存修正最大 3.86555142996954e-18；原 roundoff 面修正最大相对值 2.25514051876985e-17，单阶段最多修复 1 面。这里通过的是完整接线和 Q2 再现，不能由零校正样例证明非多项式预算选择或曲率收敛。

首次构建通过，6.1666426 s；完整 8 场景执行 271.0078807 s。结果与所有过程记录保存在 [I12a11r21_20261002](../../runs/baselines/I12a11r21_20261002/) 的 `translation_full.log`。

## 跨面与回归

48 组平面任意步内跨面全部通过原门槛：最大相对库存误差 2.22044604925031e-16、面量误差 2.20871852879411e-15、first-jet/normal 误差 6.10622663543836e-16、绝对曲率 5.74128228499174e-14。新增 smoke 的 8 场景及 20 个完整步骤也通过。

11 项 CTest 全通过，120.9750231 s；9 项旧完整输出（含 architecture、旧 smoke 和旧平面跨面）与归档逐字相同。旧解析反例 `--crossing-probe` 保持相同输出及 exit 1。没有 production 修改，因此只重建 benchmark；其余 32 个运行目标哈希保持 r20 版本。无新增 warning。

81 个冻结文件、3 份 PDF 和 r20 archive hash 不变。`driver_comparison.json` 汇集 8+8+48 场景及 552 条闭合/488 条区间面诊断，`verification.json`、`regression_comparison.json` 和 `source_audit.json` 记录验收；`delivery_manifest.json` / `archive_verification.json` 给出 `evidence.zip` 完整性核验。

## 尚未解决

当前能力为单轴恒定平移，XY 或 RZ 轴向。下一项使用同一预定预算执行原非多项式 h/dt/曲率矩阵，并同历史严格闭合结果比较；失败须保留门槛和日志，不能扩大几何或 Hessian 预算以通过。一般 rotation/vortex 真正区间积分、完整动态 AMR 与原 Gate I 清单仍未完成，**I12 / Gate I 仍 FAIL**。
