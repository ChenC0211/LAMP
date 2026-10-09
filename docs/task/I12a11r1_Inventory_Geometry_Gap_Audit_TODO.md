# I12a11r1 / 库存、旧几何与平移扫掠差额诊断

## 目标

只定位 I12a11 的 XY AMR n=32、t=.03 后半步库存越界：原始 M 更新中，旧库存与旧几何的差额、几何扫掠与平移区域测量的差额各占多少。不是生产修复，也不预先认定 I10 的 1e-10 残差是原因。

## 范围

0 production，1 test：既有 benchmark，仅显式 `--translation-inventory-audit` 启用只读探针。另含本 TODO/报告、README 和 `runs/baselines/I12a11r1_20261001/`。不新增 CTest 长程默认任务。不修改任何 primary/M、B/J、积分、闭合、修复预算或验收 gate，不更换运行参数。

## 数学约束与实现

- 重放同一 exp(q)-1、XY AMR n=32、dt=.01、4 步；仍由原 advanceTranslation 完成前三步并重现半步拒绝。
- 在真实 face transfer 已生成、原 assembler 调用之前，独立只读累加各 actual face 的有符号液量，列出所有超过原 64 eps 容量预算的原始库存。
- 对这些实际单元记录旧 M、I11 published raw/phase volume、同一 PhaseRegionSource 的原单元量及名义平移后矩形量；分解 `Mnew - Vshift = (Mold - Vold) + (Vold + exchange - Vshift)`，同时记录相关实际面四个 gross volumes。这里的区域测量是浮点估计，不伪称严格体积界；保持源 identity/snapshot 一致。
- `Vold + exchange` 只作为纸面反事实数值输出，不送入 assembler、candidate 或 advance，不作为修复后的库存。
- 原 assembler 仍真实运行并抛出原拒绝；诊断异常单独记录，不替代原异常。验证 accepted identity/M 不变，新增模式的诊断 PASS 不把原推进失败改成成功。

## 验收与失败规则

T0：Release 无 warning，直接旧回归输出相同；所有 production、27 其它可执行、81 冻结、3 PDF、I12a11 archive 哈希不变。

T1：准确重现同一 accepted time/失败阶段/实际 cell，完成上述来源一致的差额分解并记录误差尺度；若仍不能区分来源，诊断按 FAIL 归档并继续更窄的测量。不得提高舍入预算、删除异常或用解析 M 替换 actual M。

交付完整日志、源码/runtime/hash archive 与有证据支撑的下一修复建议。I12a11/I12 仍 FAIL；曲率非单调问题另行验收。
