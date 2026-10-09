# I12a11r3 / 首步纯相闭合差额审计

## 目标与范围

只定位 I12a11r2 首步末 cell 1513 的纯相残差拒绝，不预设它只是 16 eps 与 64 eps 的差别。0 production + 1 test（既有 benchmark），另含本 TODO/报告、docs/README 与 `runs/baselines/I12a11r3_20261001/`。所有生产容限及失败行为保持 r2。

## 同源分解

新增显式 `--translation-pure-gap-audit`，重放原 XY AMR n=32、非多项式、dt=.01 的 4 步请求。在原 half/full 库存返回、I10 消费之前，只读构造同一预测 scalar，按现有纯相证书查找超过实际 16 eps 工作残差的单元。保留真正 I10 拒绝及 accepted 状态隔离。

对命中的实际单元记录 old M、旧 published volume、旧区域 volume、原始与修复后共享 face exchange、目标 M、同源平移后矩形量、预测 certified/raw volume 与容量。将 `target-predicted` 分解为旧 M/区域差额、扫掠差额、原有限舍入修复、装配舍入、平移区域/预测几何差额。逐 face 输出引用同一实际 transfer，验证分解恒等式并保留浮点缺陷。moment estimates 如实列出，不冒充严格界。

不以解析库存替代 M，不造气泡，不调容限，不接受诊断候选、不修改 actual flux，不增加第二套生产 authority。若纯相判据或下游测量本身失败，报告诊断失败。

## 验收与交付

T0：Release 无 warning；production、其它可执行、81 冻结/参考文件、3 PDF、前序 archive 不变。重跑既有 11 项直接回归，与 r2 对应输出逐项比较，architecture 只允许新增共享预算 header 的 include。

T1：原 first end_closure/cell 1513 拒绝和 accepted error 精确重现；命中单元满足大于实际工作残差。分解量和 native face 证据足以区分精度分配、旧库存误差、扫掠和空间重构的贡献，不将此诊断 PASS 写成原推进 PASS。

完成后归档源码/runtime/日志/哈希与报告，再根据实测根因拆最小修复，继续推进 I12。原 I12a11r2 与父 I12 均保留 FAIL。
