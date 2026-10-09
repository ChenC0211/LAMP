# H13 定位中发现的 H08a 直接阻断反例

2026-10-06；检查目录 `runs/checks/H13/20261006_074226`。这是用户“r11最后关口”规则的限定例外核验，不开启新的source精度、Jacobian阶数或低秩任务。H13仍为唯一任务，H14不启动。

## 实际输入与失败

原生产源码、库和测试保持070108归档身份。均匀n8/n16、闭域native无散速度、原impermeable边界与同一regular reconstruction；速度来自H13对角毛细载荷的负谱分支。输入逐面值、原程序和全部原始查询保存在检查目录。这个速度只是可复现的输入，不进入生产问题选择。

原 `traceBackward` 在n8、速度倍率1e-6、到达点(0.875,0.5)、dt=0.001返回

`departure=(0.875,0.49999999992515193)`，

`J=[[9.6771001923, 2.7974e-11], [12.1519353954, 1.0000000000]]`。

实际速度约(-6.72e-25,7.48e-8)，场梯度约1e-7；三个差分距离独立复算同一个原点映射的导数均接近单位阵。偏差不是新的精度目标：它违反既有“返回同一点映射的Jacobian”以及原FD 5e-7门槛。全部分片连接仍通过原完整C1证书；该轨迹只有一次透明连接、零次物理/非C1分步。

H08a `CchMaterialTransport::prepareSpaceTime` 的 `sampleGeometry` 原样把该映射和J传入 `prepareHermiteGalsFirstJetStep`，后者以 `J^T grad(phi)` 更新独立FirstJet。实际endpoint调用中，n8倍率1e-4/1e-5/1e-6的标量时间差分相对误差为3.53e5/2.31e6/8.21e7；倍率1e-7直接在唯一publication处因伪闭合subcell cycle拒绝。n16也复现大误差。完整原材料运行两次被显式中断，标作未完成，**不以耗时或中断本身充当反例**；直接阻断证据是原运动的错误导数及H08a必经FirstJet发布。

## 数学原因与限定修正

内部Picard积分沿同一二次预测曲线X(theta,a)，在真实分片区间[a_j(a),b_j(a)]用GL9精确积分同一Q4/线性时钟多项式。旧实现分别求节点和区间宽度导数，产生与穿面法向速度倒数成比例的大数，最后才数值抵消。近切向穿片时法向量仅为舍入量，导致严重相消。

对通过原C1证书的内部连接，Leibniz边界项在相邻片间为

`(w_left(X_j)-w_right(X_j))*dtheta_j/da = 0`。

因此可先代数消去这些内部边界项，再在每个原区间积分固定theta的导数：

`dI/da = sum_j integral [Dw * dX/da + partial_t(w)*dt/da] dtheta`。

这仍是同一Picard点映射的导数；多项式次数仍在原GL9的精确范围内。只适用于既有 `predicted.side<0` 分支的已认证连续内部片；外边界、非C1事件和原dense事件导数不改。不引入小速度阈值、裁剪、Jacobian替换、不同来源、额外求积或fallback。

## 隔离结果与生产部署核验

`prototype/CchCharacteristic_stable_reference.cpp` 只做上述代数重排，链接原其他生产库。n8/n16、三个速度倍率的原点/J独立FD重新核验；八个真实FirstJet发布完成。n8原样完整H08a调用的四个微幅区间均正常完成，实际标量变化恢复相应微小尺度，等密度材料L1 R约2.7e-15至5.1e-15；R未被置零。

后续只允许修改现有 `CchCharacteristic.cpp` 的这一段及扩展既有characteristic回归（保留实际输入fixture），沿用原FD/精度/身份/物理事件/预算门槛；原H08a材料与受影响契约必须重跑。原完整H08a研究尚未因局部参考通过而自动继承新版本PASS。完成后立即回到H13原静滴、频率和能量问题，不能据此更改regular权重、材料通量、R/uR、KEP、质量或Pure FirstJet。

该限定生产修改已实施。实际fixture的原库负控误差12.1519并明确失败；新库25/25原CCH及架构/目录检查通过（337.56秒）。独立66条轨迹的点值修正前后精确相同，FD误差最大2.62e-9；八次FirstJet发布及生产库四个完整H08a材料区间通过。详见[独立复算](../../runs/checks/H13/20261006_074226/counterexample_audit.json)。标准28窗口H08a研究现已正常退出0，XY/RZ与uniform/AMR的四组空间及四组固定窗口时间门槛、实际R历史自收敛均通过；逐单元缺陷恒等式保持。完整独立审计见[084708检查点](../../runs/checks/H13/20261006_084708/h08a_complete_audit.json)。这一处等价求导修复收口，继续按用户限定冻结common regular source/characteristic；本结果不授权其他source/Jacobian/低秩精修。原毛细n8一周期仍有45.2792峰值能量缺陷，H13保持FAIL。

历史r11/H08a规定范围的已接受结果保留；这个新输入拓展了已知失效范围，报告不得宣称原程序无此问题。H13负谱分支和未通过的隔离动能/正权方案仍保留，不能用本修复替代H13验收。
