# I12a11r21 / 有界库存平移驱动配置

现有 standalone benchmark 的 `Case.boundedInventory` 显式选择 r17/r20 几何精度闭合与 r18/r19 库存面通量。默认保留旧契约，供原测试和归档比较。两种契约复用原 `initial/predict/finish/advanceTranslation`，没有失败后切换或第二份推进公式。

新命令为 `--bounded-translation-smoke`、`--bounded-translation-full`、`--bounded-crossing-publication`。前两者分别执行同一 XY/RZ×uniform/AMR Q2 测试的 1 步/40 步平移，并都执行原 4 步静止；第三者执行原 48 组平面跨面解析验证。新增 CTest 仅注册 smoke 与平面跨面。

验证配置在执行前固定：几何容限沿用默认 1e-10；每个叶的物理 Hessian 分量新增变化预算为 `min(hx,hy)^2`，系数 1，单位方域的无量纲基准。预算随加密为 O(h²)，不取自真解误差或运行结果。此处只规定 benchmark 配置，不宣告生产应用默认预算或曲率误差上界。

每次完整步骤仅从旧 accepted state 准备一次库存区域 owner。半步/全步的 actual face 传递均消费此 owner，并共享旧 scalar、owned velocity、同一 B/J 来源及明确入流历史。midpoint/end 的 first jets 和 center 仍按原 characteristic 路线推进；M 由对应旧时刻至阶段时刻的唯一共享通量更新。原 I10/I11 发布检查与 accepted 状态替换规则保留。

新增过程输出字段：

- `BOUNDED_CLOSURE`：case、XY/RZ、AMR、step、真实时刻、阶段、几何容限、最大请求 Hessian 预算、初始/最终几何差、pure 差、实际最大 Hessian 扰动、校正叶数、有符号/L1 几何差。
- `BOUNDED_FLUX`：case、XY/RZ、AMR、step、区间两端、旧几何差最大值/L1、最大单面库存修正、所查询面区域的 L1 修正之和、原 roundoff 面修正最大值与面数、库存 ledger 和 reflux 残差。

面区域 L1 求和仅为该区间面查询的诊断，不是互不重叠的全域误差积分。Report 不控制数值分支；实际 producer 类型在驱动中检查，不允许混用旧几何传递。

本接线能力限单轴恒定平移；一般 rotation/vortex 尚无该区间积分器，原非多项式 h/dt/曲率矩阵仍需独立通过。见 [TODO](../task/I12a11r21_Bounded_Translation_Driver_TODO.md) 与[报告](../task/I12a11r21_Bounded_Translation_Driver_Report.md)。
