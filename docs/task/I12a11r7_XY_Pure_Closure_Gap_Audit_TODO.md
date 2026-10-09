# I12a11r7 / XY 新纯相闭合拒绝差额审计

## 目标及范围

只判断 r6 两个新拒绝的主要差额来自旧库存与几何测量、通量扫掠，还是预测几何：XY uniform n=32 的 cell 531，以及 XY AMR n=16 的 cell 341，均在接受 t=.03 后的 midpoint closure。0 production + 1 test（既有 benchmark），另含本 TODO/报告、README 和 `runs/baselines/I12a11r7_20261001/`。不修改 M、预算、门槛、B/J、闭合或推进路线。

## 实现

新增显式 `--translation-xy-pure-failure-audit`，重放以上两个原四步请求。在 step=3 的真实 midpoint candidate 上复用只读纯相差额审计，分别限定 cell 531/341。对真实同源量分解：旧 M 减旧区域体积、区域扫掠余差、raw/final flux 修复、库存组装舍入、平移区域与预测纯相的差额。保留原面通量、测量误差估计和分解残差，不将估计当认证界，不把私有反事实用于推进。

诊断异常单独计数并继续原闭合；最终必须仍由原 I10 拒绝。检查 actual accepted/candidate identity、M 和几何绑定状态未被改变。原 CLI 模式不开启此诊断。

## 验收及失败规则

T0：Release 无 warning；11 项直接回归逐字等同 r6；production、其它 27 可执行、81 冻结/参考文件、3 PDF 及 r6 archive 保持不变。

T1：两个原拒绝和最后 accepted error 重现，仅允许独立重放导致的 primary 序号差异。每个目标单元有完整差额分解及面通量，失配超过原 64 eps 验收预算，分解余差不超过该预算；诊断异常为零。未获得完整数据则本诊断 FAIL，不推断根因。

交付：修改文件、原始测试日志、差额表、未解决问题、源码/runtime 哈希和 archive。诊断 PASS 不等于 I12 通过；归档后按数据选最小修复 TODO，继续推进原门槛。
