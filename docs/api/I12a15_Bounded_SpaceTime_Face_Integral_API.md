# I12a15 / 有界面—时间积分 API

状态：**FAIL（预设矩阵尚有 1 项超工作预算）**。见 [TODO](../task/I12a15_Bounded_SpaceTime_Face_Integral_TODO.md)、[设计](../task/I12_Bounded_SpaceTime_Integral_Design.md)、[报告](../task/I12a15_Bounded_SpaceTime_Face_Integral_Report.md)。实现已可供显式验证调用，尚未接入完整 driver，不能视为 I12 通过。

## 入口与所有权

`integrateHermiteGalsPolynomialFace(mesh, interval, nativeFace, config)` 与显式 `inventory` 重载返回私有构造的 `HermiteGalsBoundedFaceIntegral`。要求 owned polynomial velocity、完整旧 native Q4、同一 mesh，以及可认证的原完整 RK2 单步无 hit 分支；库存重载另核验 old identity 与 snapshot 指针。结果保留 interval/source、old state、可选库存 source，复制共享只读 owner；moved-from 和过期 mesh 拒绝。

`transfer()` 是原 `HermiteGalsIntegratedFaceTransfer`，四个量依次为 total+/total-/liquid+/liquid-；全部按同一真实 face 和时间区间积分，方向为全局 +x/+y，RZ 含 2*pi*r。`bounds()` 提供四个量的上下界；`report().absoluteErrorBound` 是四项宽度的 outward 和；`config()` 保留原 controls；`inventorySource()` 在显式库存入口非空。

默认 `absoluteVolumeTolerance=1e-8`、`maximumPanels=4096`、`maximumDepth=40`。成功要求绝对体积误差界不超过该容差，并且所有非零通量叶的 branch 都已认证。容差不提供额外库存修复额度。最大分块是当前树的叶数；`evaluatedPanels` 包括为选择方向而试算后未保留的候选；本实现不超过 `1+4*(acceptedPanels-1)` 次完整分块计算。`phaseQueries`、`branchRefusals` 同样包括候选。`maximumDepthReached` 包括候选最大深度，`unresolvedPhasePanels` 为最终叶中的未确定相数，不代表分支未认证。

预算耗尽、不可表示的继续二分、未知边界历史、非法数据/环境均抛异常；没有返回部分积分的成功句柄，也没有修改 M、geometry 或 accepted state。严格零通量的积分可为零，但仍检查原 domain/opaque coverage，轴面另查精确径向正则条件。分支未知不会因所占误差较小而被吞掉。

## Core

`PolynomialVelocitySlab::encloseNormalFaceIntegral(a,b,t0,t1,radial,rhoLo,rhoHi)` 用同一多项式 monomial moments 包围有符号积分，结合速度符号范围得到 gross 上下界。它只消费独立给定的含量范围，不持有相或状态。端点必须是升序、正长度、轴对齐的线段，时间为 slab 内正长度区间。

`PolynomialFaceIntegralBounds` 提供 outward `added`、`intersected`、`withContentRange`、`errorBound` 与位于范围内且 liquid<=total 的 `representative`。这些原语在既有 strict/no-LTO TU 执行。`PolynomialRk2EnclosureUnresolved` 标识可细分的分支尚未认证，沿用原错误文字；非有限算术、非法参数、浮点环境不作为可细分分支吞掉。

`encloseInteriorRk2TightTime` 是原 RK2 区间表达的附加收紧入口，覆盖真实求值的时间混合舍入；原 `encloseInteriorRk2` 与实际求值规则保持。详见设计中的误差带。样本用于检查与工作方向选择，成功证书不依赖采样密度。

## 使用边界

成功得到 transfer 不代表原库存装配一定可接受它；库存严格界限、散度/守恒、边界历史、reflux 与完整几何推进仍需各自验收。当前没有一般物理 inflow/history 或多缩步证书；本任务不解决完整旋转/涡旋输运、时间收敛、动态 AMR 和低黏/无黏 Hydro。