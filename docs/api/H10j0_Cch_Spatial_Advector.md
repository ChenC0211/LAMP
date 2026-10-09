# H10j0 / 同源空间速度来源

现行接口：[r11 第二次恢复](../task/H08a2r11_Common_Regular_Source_Deployment_Resume2_Report.md) **PASS**，共同正则场与同源运动已部署。r11 为最后基础设施关口，下一立即回到 H08a 真正同源时空相通量积分；本 API 不代表该积分或 R 收敛已通过。

Owner：`lamp_hydrodynamics_cch`；实现 [CchSpatialAdvector.cpp](../../src/Modules/Hydrodynamics/CchSpatialAdvector.cpp)，公开 [header](../../include/Modules/Hydrodynamics/CchSpatialAdvector.h)。几何准备和实际速度 revision 共用一个 owner；无字典、IO、phi、密度、材料库存或接受事务。

```cpp
const auto reconstruction = CchAdvectorReconstruction::prepare(mesh, momentumGeometry);
const auto advector = reconstruction.evaluate(mesh, actualValues, velocityGeometry, fullClock);
advector.validateFor(mesh, momentumGeometry, velocityGeometry, fullClock);
```

`prepare` 持有原 momentum geometry/layout、唯一 D_I 和不可变 [regular 准备](H08a2r4_Cch_Regular_Advector_Preparation.md)。真实 cell/support knots 构成 `pieces()`；原 primal/dual 物理段共同切分成唯一 canonical `atoms()`，`primalAtoms()/dualAtoms()` 只保存有序索引。内部符号由原 incidence 给出，AMR 使用真实子面和。旧 `patches()` 和 trace/lift 路线已退役。

每次 `evaluate` 复制原 native values，并只调用一次同一准备的 `apply`。实际场持有所得 Q4/面均值修正系数，所有值/J/normal 查询都从它取得；更新不重复 SVD。私有构造的 velocity identity 每次 revision 新建；geometry identity、完整 clock 和显式 boundary identity 保持，消费者精确比较。时钟标签是调用方声明，不证明 native 数组的物理时刻或 H14 耦合提交。

`cellVelocity`、`cellVelocityJet` 查询指定真实 cell 内的物理速度及解析 Jacobian；`CchAdvectorCellJet {jet, lower, upper}` 给出实际 support 闭包。内部 knot 默认正侧，cell 上端取最后片。`pieceVelocityJet` 明确选择当前 reconstruction 的片索引，允许闭包两端；同一 C1 场在连接处共享值与一阶导数，不平均邻接值。`nativeVelocity(u)` 返回该实际场在 native 面中心的完整向量，**不是**把 native 面均值当点值。`cellDivergence()` 保留原 D_I/V 平均量；点态及 dual 分配差仅作诊断，不强求逐点常散度或 R=0。

`atomPolynomials()` 同时涵盖 primal 及 cell 内 dual normal 段。每个 atom 已切开所有严格内部 support knot。`CchNormalPolynomial` 在 begin/end 两端分别保存沿正距离的完整五系数，可含非零三/四次项；RZ 物理测度由线性半径乘入积分。它是 metric 数据，不是相几何/root authority。

`value(fromEnd, distance)` 和 `integral(fromEnd, distance, width, normalSign)` 使用物理切向距离。Core 的 0–1 offset/width 须分别乘原 atom length；不能先重构舍入后的全局端点再作差。fromEnd 系数与半径方向已反向，测度仍为正；normalSign 只取 ±1，邻侧率复用同一积分取反。独立 width 保留无法表示全局 endpoint 的正短范围。

支持 XY/RZ、真实 2:1、selected/disconnected、axis 和合法 Initial/Midpoint/End clock。P2 及可辨识低阶由同一准备的层级重构和原完整矩阵验收；缺秩仍公开报告，不能称所有稀疏布局 P2 精确。native 物理均值、真实 cell 积分散度和共享面由同一场核验。显式 impermeable law 只约束调用方选定的真实外部面，轴奇偶条件始终有效。

副本保持原 prepared 身份，不保存 mesh 裸引用；调用方数组/局部变量释放后仍持有完整数据。每批消费前检查当前 mesh，regrid 必须新建准备。wrong layout/momentum/source/boundary/mesh、未绑定几何、非法时钟、非零 axis/wall 输入、nonfinite、非法 FP 模式、越 cell/piece 或非法 shifted range 明确拒绝。失败不修改 owner 或 FP 模式。

同源 [characteristic](H08a2r1_Cch_Characteristic_Source.md) 消费该场，材料相率也消费同一 atom。FirstJet 运动桥接尚未完成；现有材料时间仍为中点近似，本空间 API 不冒充同源时空材料积分。保持唯一共享 F、内部反号/真实 AMR 子面和、来源与完整 clock 一致、KEP/skew 四硬约束，继续显式处理真实 R/uR。

浮点消费者使用完整物理 chart 的 `sum |c_k| L^k` 和实际 native 列响应定标；零期望积分不以残余噪声自身作尺度，零 work 要求零残差。原 tiny/counterflow、独立物理积分、来源拒绝和完整次数负控保持；详见 [r12](../task/H08a2r12_Regular_Source_Arithmetic_Consumer_Report.md)、[r13](../task/H08a2r13_Benchmark_Source_Observer_Report.md) 和当前部署报告。
