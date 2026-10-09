# I12a11r5 / 求根工作精度与同源体积测量审计

## 目标及范围

定位 r4 的 cell 1454 可达区间拒绝，并量化同一严格闭合对首步曲率的影响。只检验一个事实：实际工作残差与同源体积测量差额的尺度是否一致。0 production + 1 test（既有 benchmark），另含本 TODO/报告、README 与 `runs/baselines/I12a11r5_20261001/`。不更改任何 gate、M、B/J、scalar 求根或解析参考。

## 实现

显式 `--translation-working-precision-audit` 重放原 XY AMR n=32、dt=.01、四步请求。成功的 midpoint/end 闭合只读比较实际 predictor 与 closed scalar 的位置/梯度/曲率；复用已有独立 bubble value/gradient/Hessian 核验，记录最大曲率单元的真实求根报告和目标。

对成功最大曲率单元的 predictor/closed scalar，以及失败 cell 1454 的 predictor、原证书上下端点，分别比较原 `discoverBubbleVolumes` 与已有局部 Green `PhaseRegionSource.measure`。局部测量所需 state 是显式测试 fixture 的私有绑定，DOF 必须与被测 scalar 完全一致；它不推进、不发布 I11、不会返回给 actual driver，M 也不作为几何输入。记录两种量、moment estimates、容量、目标残差及 16/64 eps 预算；不把不同求积的接近或差异当作认证误差界。

保留实际 I10 拒绝与接受时间；失败诊断自身异常须单独记录并使审计失败。所有 actual candidate/accepted identity、M、几何绑定状态保持不变。原模式不开启诊断，输出不应变化。

## 验收及失败规则

T0：Release 无 warning，11 项直接回归逐字等同 r4；所有 production、其它 27 可执行、81 冻结/参考文件、3 PDF 和 r4 archive 不变。

T1：原 t=.01 accepted error 与 cell 1454 midpoint 拒绝重现（只允许显式私有测试 fixture 带来的 identity 序号差异）。曲率实际 Hessian 变化与同一 cardinal bubble 独立公式一致；original/certified endpoints 的完整数据可判断后续是否需改变测量数值精度。不得在没有数据时预先选定全局消减、求积或空间重构为根因。

本项只诊断，诊断 PASS 不改变 I12 数值 FAIL。完成源码/runtime/原始数据/哈希归档后，再按证据写最小修复 TODO 并继续。
