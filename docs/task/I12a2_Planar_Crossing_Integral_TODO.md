# I12a2 / 常速度平面跨面事件积分

## 任务 ID / 标题

I12a2：从旧 accepted native Q4 与常速度来源，计算单一 native face 的真实时间区间积分。用户已明确要求 I12a2 后持续推进直至 I12 通过；每个子任务仍单独定义、验收和归档。

## 目标

只解决一个事实：平面在常速度下穿越单一面的相体积必须按整个时间区间积分，不能用单点 midpoint 相指示替代。

## 允许修改

最多 5 个 production 文件：现有 HermiteGalsTemporalPhaseFlux.h/.cpp、LevelSetHermiteGeometry.h/.cpp、Detail/LevelSetHermiteQuarticLeafData.h。最多 2 个测试相关文件：新 planar temporal contract 与 tests/CMakeLists.txt。另含本 TODO、报告、API 说明及 runs/baselines/I12a2_20261001 证据。

## 禁止修改

不改库存、bubble、driver、reflux、first-jet、Hydro/Electricity、字典、容差或 CFL；不改旧 midpoint 数值行为。本子任务不声称 I12 已通过。

## 数学约束

- 积分消费 I12a1 的旧 accepted authority、同一实际 velocity/characteristic source 和真实区间，输出四个非负 integrated gross volumes。
- 常速度由持有值的来源工厂生成实际 endpoint callbacks；任意 callback 即使少量采样相等也不获常速度证明。
- 从原 native Q4 的 exact homogeneous DOFs 证明扫掠覆盖中的多项式为同一 affine scalar；不以 rounded diagnostic DOFs 拟合，不容差降阶，不重建解析状态。
- 常速度下现有固定 RK2 即平移。将 face × time 上的 affine 相指示分割，解析积分；XY 使用长度，RZ 使用 2*pi*r 的真实测度。保留正反方向，内部两侧共享同一结果。
- RZ 只支持本任务物理有效的常轴向平移；扫掠必须全部处于旧 active geometry 与 domain 内。未知 coverage callback、非平面、径向常平移或需外部边界历史明确拒绝。
- 事件边界的瞬时零迹为时间零测集；覆盖正时空测度的恒零 scalar 明确拒绝。失败不改变 accepted state。

## 实现要求

几何证明和精确谓词留在现有 strict Core TU；接口模块只绑定来源和方向。使用解析时空半平面裁剪（等价于完整端点 crossing-event 分割），不使用未证明事件完备的采样或自适应求积。数值最终输出仅有表示舍入，不声称普通 double 返回值为符号精确数。

## 验收

- T0：Release 构建、architecture、文件预算与冻结文件身份检查。
- T1：原 q=1/4、3/4 的 0.00375/0.00125 反例；多相位、方向、无穿越/端点穿越、斜平面、时间组合；XY/RZ 与独立闭式参考一致，half/full 同源且真实区间明确。
- T1：非平面/未认证速度/越界/覆盖与 mesh 失配明确拒绝；接受态 identity、snapshot、M 不变。
- 回归：I12a1 的 18 项检查及原 crossing 探针；本任务未接 driver，旧探针仍失败并保留证据。

## 失败规则

未通过精确平面或覆盖证明的输入返回明确错误，绝不隐式退回 midpoint。若必须超出文件预算或改变其它数学事实，先拆出后续任务。

## 交付

代码、算子测试、API/报告、原始构建测试日志与 SHA256 归档。本子任务验收后按用户授权继续推进 I12 的下一项必要工作。

## 执行结果

**PASS。** 5 production + 2 test 文件；741 组有效检查、8 类拒绝和 19 项 CTest 通过。17 项前序输出一致，冻结文件未变。详见 [报告](I12a2_Planar_Crossing_Integral_Report.md)。I12 整体仍 FAIL，按最新授权继续推进下一独立子任务。
