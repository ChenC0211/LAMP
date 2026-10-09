# I12 / 纯 Hermite–GALS 非守恒路线设计变更

**现行授权，2026-10-02。** 用户要求：“如果当前路线还是很难实现，就转为‘纯 Hermite–GALS，不守恒’的路线”。r7 已消除第二步守恒装配证书越界，但四个原 rotation/vortex 算例第三步仍被边界来源认证阻塞。依据这一证据执行路线切换，停止把守恒库存/通量认证作为当前 Interface 主线。

本变更优先于六份规范中与 Interface 独立库存、体积闭合、相通量守恒相关的旧条款；冻结 PDF、既有代码与结果保留为显式历史比较，不能隐式回退。项目的强变形/曲率/低黏目标不变。

## 新离散状态与发布

沿用现有 HermiteGals state/identity/time/唯一 snapshot，增加显式 `FirstJet` 模式，与历史 `ConservativeVolume` 模式区分。FirstJet 状态的 `phaseVolumes` 必须为空；禁止把每步几何体积写回该字段伪装成守恒 M。

顶点直接输运 phi 和梯度；边中点输运 phi 与法向导数；B_h 与 J_Bh 仍来自同一实际 RK2 map。phi_xy 仍由 I05/I06 紧凑 compatibility closure 派生，hanging data 仍从 master 生成。局部 Q4 内部值采用同一旧 snapshot 在实际回溯中心处的 phi 值，保留现有零中心 cardinal 表示；不再求体积约束根或调整中心值来守恒。它是本步重构输入，不引入独立守恒体积或独立二阶导数状态。

新增明确的纯几何 publication 入口：输入未发布 FirstJet candidate、同 revision 的 I06 mixed result 及中心预测；复用原发布核函数生成并验证唯一 Q4 snapshot、拓扑、cut、法向、曲率、XY/RZ 测量。失败保留 candidate，成功消费，identity/time 不变。保留原表面残差、法向、有限曲率、几何积分/拓扑检查。取消几何与 M 比较，因为该模式没有 M。

几何相体积 V(phi_h) 是派生诊断，报告初始/当前值、相对初值漂移和解析误差。原 inventory residual/ledger 在该模式为“不适用”，不能输出零并声称守恒。液体体积的时间守恒、phase-flux cancellation 和 reflux 不属于纯界面接受条件；几何测量自身的闭合/有效性仍须验证。旧独立 M 接口及约束只适用于显式历史守恒模式。

## 验证与比较

I12p1 验证状态、发布、拒绝事务及原实际输运接线，并实跑原 XY/RZ static/translation、XY rotation 30 步和 vortex 400 步及 uniform/AMR。没有 half/full 守恒积分；给定速度的完整 RK2 first-jet 步由旧接受状态直接生成候选。原静止/多项式精度、几何合法性要求保持，体积漂移作为诊断。

I12p2 验证空间/时间收敛、curvature convergence 和 AMR crossing/regrid，保持原解析参考及可比较网格/时间参数。体积漂移不能隐去，也不能把取消守恒要求写成零体积误差。实际失败继续修复，不以换路线自动标记 I12 PASS。

I12p3 汇总准确性、体积漂移、AMR 敏感性及效率 gate；仅通过当前纯路线的完整独立验证后才进入 Hydro。旧守恒 I12 的 FAIL 保留历史意义，不能混同两条路线的验收。

## 后续耦合约束

Hydro/Electricity 的压力相容、黏性半正定、零黏极限、统一几何与电荷/Gauss 要求保持。原 H-R3/H07 中“由守恒 Q_l 自动保证密度/动量一致”的前提已不成立，必须在开始耦合前单独设计密度更新和质量—动量误差诊断/稳定性验收。此次仅解除 Interface 相体积守恒，不擅自改变流体/电场离散方程或声称已解决低黏耦合。

执行入口：[I12p1 TODO](I12p1_Pure_Hermite_Gals_Publication_TODO.md)。
