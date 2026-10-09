# H08a1 / 实施前设计：实际分片速度的局部微分

2026-10-05，代码前固定。[任务及数值门槛](H08a1_Actual_Advector_Cell_Jet_TODO.md)。本项不改变H10j0重构方程或用户授权的四项材料硬约束。

## 为什么先拆这一项

当前 H08 对唯一空间 F 使用中点率乘 dt。真正的时空相积分还需要一份可追溯的步内几何及其实际运动源。现有 Hydro 只有 `CchSpatialAdvector::cellVelocity`；原 FirstJet tracer 需要同一速度分支的 value/Jacobian，且明确要求连续的已解析 interior/boundary trace。原 Hdiv 场只共享法向，切向可跳跃，不能把它随意接成光滑 callback。

完整接线至少涉及 advector、phase、material 各一对 h/cpp 和 dual cpp（7 production），超出默认预算，并且跨跳跃运动是独立数学契约。因此先给同一owner补齐局部jet与精确分支范围，下一项处理实际运动/事件，再实现时空相积分和共同材料积分接线。局部jet是后续消费者的实际能力；不创建另一条速度重构，也不宣布时空来源已经闭合。

## 原场的解析导数

令 XY 中 X=x、B=(u_x,u_y)；RZ 中 X=r²/2、B=(r*u_r,u_z)。原 cell 场在实际 native-side knots 划分的矩形分支上为分片二次函数。

记 A(Y)=(B_R(Y)−B_L(Y))/ΔX，C(X)=(B_T(X)−B_B(X))/ΔY。原 lift 是对应对边差的积分减其均值乘坐标，所以：

```
∂X Bx = A(Y) − C(X) + meanC
∂Y By = C(X) − A(Y) + meanA
∂Y Bx = 原 cellLine(normal=0) 返回的 tangent slope
∂X By = 原 cellLine(normal=1) 返回的 tangent slope
```

这些量直接来自原 faceB、cellLine、lift均值和相同 sideFace 选择，没有独立 gradient fit。原 D_I/V 不重新计算，作为独立验收权威。

RZ 对 r>0：`Jrr=∂X Bx−Bx/r²`，`Jrz=∂Y Bx/r`，`Jzr=r*∂X By`，`Jzz=∂Y By`。在实际轴 cell 的 r=0，由原零轴法向及 Bx(X=0)=0 得 `Jrr=∂X Bx/2`，`Jrz=Jzr=0`，`Jzz=∂Y By`；真实散度极限是 `2*Jrr+Jzz`。这是原场的正则极限，不是任意零除保护。

## 查询范围与来源

`cellVelocityJet(pressure, physicalPoint)` 返回 `CchAdvectorCellJet`：Core 的 `VelocityJet`（value和physical Jacobian）及所选矩形分支的 `lower/upper` 物理坐标。pressure由调用者指定，不做全局 cell 猜测。内部 knot 选正方向，在cell最上/右端选最后分支，与既有 sideFace 一致。包围盒给分支闭包；上边界内部knot的另一侧jet需另行选择其分支，不能据闭包断言场在两侧一致。

分支内、原cell外边界上的返回值必须逐位等于原 `cellVelocity`。保留原值计算次序；Jacobian只增加只读计算。元数据和采样不持mesh裸引用，复制保留原owner identity、native values、geometry、full clock；消费者仍先调用原 `validateFor`。不把局部sample升级为新的时间来源，也不伪造publication identity。

## 不可跨越的范围

同一cell内部 AMR knot 也可能有切向跳跃，不能只检查原cell外边界。切向跳跃的流映射可能需要事件时间导数/传递矩阵；局部 J 的RK求积不自动包含这些项。下一项须在实际分片B上定义运动及其导数、覆盖/regrid和拒绝条件，并先验证横穿反例，再接FirstJet；不得通过平滑B、换一套拟合速度或给任意callback贴身份标签规避同源约束。

本项不修改Core tracer契约，不改变原材料 F/R 或 KEP/skew，不引入强制连续性作为第五个材料硬约束。H08a仍须在XY/RZ、2:1、密度比及真实moving/crossing/regrid上证明时空积分与R收敛。
