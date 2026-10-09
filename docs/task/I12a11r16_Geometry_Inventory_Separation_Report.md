# I12a11r16 / 几何精度与守恒库存分离报告

**PASS（授权设计与文档）；I12 / Gate I 仍 FAIL。** [TODO](I12a11r16_Geometry_Inventory_Separation_TODO.md) 与[设计变更](I12_Geometry_Inventory_Accuracy_Design_Change.md) 已落实用户明确授权：取消逐 cut cell 的机器精度几何体积硬匹配，保留独立库存的局部共享通量守恒和严格界限。01 对应条款、05 的当前 I12 顺序与 README 入口已更新。

设计明确了 `V(phi)-M` 与守恒残差的区别、既有弱 cut 条件性证据、导数扰动预算、纯相/拓扑边界以及单纯调大容差仍会产生负库存的风险。后续顺序改为精度预算闭合算子、共享通量有界性协调、半步/全步接线和原 I12 验收。r14/r15 保留为受限算子能力，不再把加密恢复机器精度闭合作为主线前提。

本任务没有改动数值代码或测试源码；不重复运行数值测试，不声称新推进已通过。文档链接、源码哈希、81 个冻结/参考文件、3 份 PDF 和 r15 archive 均由[审计脚本](../../runs/baselines/I12a11r16_20261001/audit.ps1)检查。原 PDF 和历史归档不变，新的例外由用户授权及独立变更说明表达。

证据：[目录](../../runs/baselines/I12a11r16_20261001/)。`verification.json`、`manifest.json`、`archive_verification.json` 和 `evidence.zip` 记录源码不变、链接及文档交付哈希。下一项先实现局部闭合契约；库存协调未通过前不切换完整驱动，也不进入 Hydro。
