# I12a16r1 / 显式边界相历史报告

**PASS（本项来源与立即命中证书）。a16 完整驱动、I12 / Gate I 仍 FAIL。** 见 [TODO](I12a16r1_Boundary_Phase_Source_TODO.md) 和 [说明](../api/I12a16r1_Boundary_Phase_Source.md)。取消逐 cut cell 几何体积的机器精度硬匹配后，独立库存、几何预算、原守恒装配与64eps修复限保持。

## 实现与证据边界

4 production+2 tests。新增常相四边 history 保留未绑定 interval 的旧状态、owned polynomial velocity、原 characteristic 和完整声明时间；不同 history 即使四边值相同也不是同一来源。interval、transfer append、inventory assembly 全链保留并核验该来源。未给 polynomial 伪造恒速能力，也未从旧 phi 或 M 推断边界气相。

Core 的初始命中证书覆盖整个闭空间块和每个正持续时间，验证实际速度严格区间、原构造/CFL/轴条件，且唯一 outward-backtrace side 必须相同。原 tight-time 界提取成同一私有函数供复用，运算顺序保持。无证明时继续要求原完整 no-hit 证书；没有 exterior velocity、替代回溯或待定分块误差豁免。dt=0 原 identity 不变，时间积分只明确排除该零测度时间切片；未放松空间角点或正测度事件。

benchmark 从其既有 rotation/vortex 解析物理条件显式声明四边气相：径向支持上界 hypot(cx-.5,cy-.5)+radius 小于.5。该证明只提供边界数据。测试中特意令内部几何气相 M=.25C，确认即使声明边界气相，内部通量仍含原 gasWeight=.25；没有以几何纯相清零库存。

## T1 / T2

扩展原积分测试：28个边界积分场景，包含四边、正负流向、XY/RZ适用边、uniform/AMR、气/液历史、空间和时间变化速度；独立多项式原函数与RZ真实径向矩作参考。1,134次原实际 tracer 比较包含正持续时间立即命中，以及零持续时间 identity。48项拒绝/来源检查覆盖不同历史的append、assembly丢失history、旧状态/速度/map/time/mesh、移动对象、缺历史、轴/CFL/opaque/切向/多壁/反转等。全部通过。整项contract首次.6601404 s，重建后.6841406 s。

原四个driver first steps均实际执行，首个入流 native face0 由4096 pending改为 **1个认证panel**。rotation的四宽误差界1.6263032587282568e-18，vortex为2.0599841277224588e-18，预算仍1.5625e-8。

随后四场景均在 **内部native face1（negative0、positive1）** 拒绝：41panels、pending1、depth40，`RK2 enclosure new wall contact unresolved`。rotation E=5.7006035634153762e-16，vortex E=4.4593211213733091e-16，虽低于预算仍不得发布。实际mid/end预测保持同源通过，accepted时间0、旧identity/M保持；没有完成半步装配/闭合。smoke 3.0294733 s、full 2.9928453 s，均exit1；full仍请求原30/400步，不减少为已认证面。

## T0 / 归档

首次构建28.0426029 s成功；全部39targets重建71.1131083 s、exit0、无warning。26CTest用时162.5154204 s：25PASS、原完整driver FAIL。24旧完整输出逐字不变；扩展contract中原57主积分、864条tight-time轨迹及原guards的完整输出前缀也逐字不变。旧crossing probe仍exit1、输出相同。

架构12targets、383TU、908files、3118工程includes、5exact legacy edges保持。81冻结文件、3PDF、a16 archive SHA保持。首次审计误把新增标准库`<array>`当成工程依赖计数+1；核对scanner实际输出后修正审计为工程依赖数不变，保留首次审计失败记录。编译和数值测试无需重跑；最终审计PASS。

证据：[验收](../../runs/baselines/I12a16r1_20261002/verification.json)、[驱动](../../runs/baselines/I12a16r1_20261002/driver_verification.json)、[回归](../../runs/baselines/I12a16r1_20261002/regression_comparison.json)、[范围](../../runs/baselines/I12a16r1_20261002/source_audit.json)。归档含6源文件、本文/TODO/说明/入口、39runtime及原DLL，逐项核验payload SHA。

## 下一项

补原离散RK2的壁面接触排除证明：现判据在闭块含壁面起点/零持续时间极限时，仅要求 departure 整块严格远离该壁，未利用速度符号导致的坐标单调性。用原浮点表达式可证明的单调性排除不可能的新接触，再运行同一四个first steps；不能按E很小丢掉pending。随后仍须处理真实后来命中、完整rotation/vortex及时间/空间/AMR/曲率gate，未进入Hydro。
