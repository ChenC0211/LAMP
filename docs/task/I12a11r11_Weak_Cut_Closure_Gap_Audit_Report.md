# I12a11r11 / 弱灵敏度 cut-cell 闭合残差来源报告

**PASS（只读定位）；I12 数值仍 FAIL。** 对应 [TODO](I12a11r11_Weak_Cut_Closure_Gap_Audit_TODO.md)。三个被选中的大系数修正，其闭合前体积差几乎全部来自精确平移区域与 predictor 重构区域之间的差额。

## 实现及不变性

0 production + 1 test。`--translation-weak-cut-gap-audit` 复用原四步 driver，仅在 end 的 step 1 / cell 648、step 2 / cell 650、step 3 / cell 363 记录守恒分解。数据来自实际 old PhaseRegion、raw/final integrated faces、assembled M 和实际 predictor DOFs。private predictor fixture 仅适配原只读测量 API，DOFs 与实际 predictor 逐位一致，未发布或反馈推进。

3 份预测体积、目标 M 与随后的真实 root 输入逐位相同；实际 c、预测 c 和最终 dV/dc 与 r10 原始记录逐位一致。四步最终 BENCHMARK 行也逐字一致，位置/梯度/曲率误差仍为 7.8330e-8 / 2.6939e-5 / 1.0386e-2。old/candidate identity 与 M 不变。

## 残差分解

符号约定为 target M − predictor volume；五个加数依次为 old M − old geometry、old geometry + raw exchange − shifted geometry、final − raw exchange、assembled M 的舍入项、shifted geometry − predictor geometry。

| step / cell | 实际体积差 | 相对 capacity | old stock | sweep | repair / assembly | reconstruction |
|---|---:|---:|---:|---:|---:|---:|
| 1 / 648 | +1.0649414e-16 | +1.0905000e-13 | 0 | +2.5411e-20 | 0 / 0 | +1.0646873e-16 |
| 2 / 650 | −6.0606901e-17 | −6.2061467e-14 | 0 | 0 | 0 / 0 | −6.0606901e-17 |
| 3 / 363 | +1.0972126e-16 | +1.1235457e-13 | 0 | 0 | 0 / 0 | +1.0972126e-16 |

逐项重组缺陷全部为零。cell 648 的 reconstruction 占 99.9761%；另两个为 100%。因此，本次残差并非旧 M 与旧几何不一致，也不是 flux 修复或库存组装引入。三个残差绝对值分别为原 64 eps capacity 预算的 7.6737 / 4.3672 / 7.9063 倍，不能按原预算直接忽略。

## 灵敏度及实际闭合

| cell | 预测较小相 / capacity | 初始 dV/dc | 实际 Δc | ΔV / 初始 dV/dc |
|---|---:|---:|---:|---:|
| 648 | 1.6439305e-4（液） | −1.9103734e-9 | −5.5745199e-8 | −5.5745199e-8 |
| 650 | 1.8984978e-5（气） | −8.7807390e-12 | +6.9022552e-6 | +6.9022552e-6 |
| 363 | 2.0166697e-5（气） | −1.0983509e-11 | −9.9896362e-6 | −9.9896362e-6 |

实际修正与初始线性预测一致，且三个最终体积均等于目标 M。本轮没有显示 Newton 走错根或异常过冲；大系数修正与体积对 center 系数的弱灵敏度一致。r10 已记录这些修正所在单元在后续重构中成为最大曲率误差位置，但本项仍不把两者的关联提升为完整稳定性证明。

cell 650 的 shifted/predictor 气量为 1.8540017155e-8 / 1.8540017095e-8；cell 363 为 1.9694040014e-8 / 1.9694040124e-8。预测两种 Green 气量之差分别为 1.5551e-22 / 9.5953e-23，远小于需要修正的体积差；对应 moment estimates 为 0 / 3.3087e-24。cell 648 的双 Green 液量差为 −7.9409e-23。这些是独立方向一致性与求积估计，不能当作严格积分误差上界；没有据此调整任何容差。

## 验证与下一项

Release 构建 3.93 秒，无 warning。完整四步诊断 135.48 秒、退出 0，3 gaps / 3 roots / 9 region measures / 12 face records / 3 组各 25 个 DOFs 齐全，异常为零。11/11 直接回归通过（约 47.2 秒），完整输出与 r10 相同；旧 midpoint crossing probe 仍预期退出 1。回归与四步诊断并行，耗时不作严格性能比较。

其它 27 可执行、全部 production、81 冻结/参考文件、3 PDF 和 r10 archive 不变。原始日志、逐项/比例分解、DOFs、字段表、源码/runtime 哈希、交付清单及 archive 位于 `runs/baselines/I12a11r11_20261001/`。

下一项以独立可解析的 corner-cut 算子试验检验 Q4 cardinal bubble 灵敏度随小相尺寸的退化，区分数值实现缺陷与表示/闭合本身的条件性。任何改变 first-jet、单个 cell bubble 或允许误差预算的方案都需要具体方法论证；本项没有实施这些变化。I12/Gate I 仍 FAIL，不进入 Hydro。
