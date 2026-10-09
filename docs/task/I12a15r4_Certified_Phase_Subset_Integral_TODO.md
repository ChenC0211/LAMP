# I12a15r4 / 分块内认证相子区域的 gross 积分界

目标：在原 a15 完整分块积分里使用 r2/r3 的相关回溯与权威 Q4，认证每个 face-time tile 内明确液/气的子区域，收紧同一 gross 液体积分。必须用原 uniform shear 的 4096 panels、depth 40、5.859375e-6 误差预算验证，不提高预算或换算例。

数学约束：原闭 donor 覆盖与 single/no-hit RK2 分支先认证；所有候选叶的 scalar forms 转成一个包含全部原叶的共同 affine band。上包络<0 的半平面是明确液相，下包络>0 的半平面是明确气相；只能作用于已经完整覆盖的参数域。用精确有理裁剪取得 unit-square 面积与第一径向矩，最终向外转换，边界直线为零测集；恒零多项式不得认证任何正测度相。库存权重取原同源候选 donor 的 wl/wg 范围。既有总 gross 界保持，液体界按已认证区域改进，再与原界相交。

允许修改：Slab h/cpp（Core 区域积分与复用 strict arithmetic）、原 Snapshot cpp（提取已有 exact rational/半平面矩以复用）、新 Detail/ExactFaceTimeRational.h（原算法原序搬移及共享裁剪）、TemporalPhaseFlux cpp（已认证来源内接入），最多5 production；已有 affine_phase_contract 与 bounded_spacetime_flux_contract 两个测试。先写设计。无 root CMake/配置迁移，不另起第二套相场或根求解器。

禁止修改：原 velocity/RK2、Q4 原始数据/constraints/sign、M/库存权重定义、reflux/driver、benchmark参数、误差与 panel/depth 预算、64eps 修复限制或任何 gate。样本只能检测或选工作，不能决定事件/相积分。不能以新 primitive PASS 替代原 shear。

实现要求：复用 Snapshot 原精确有理运算与裁剪表达式，原入口保持运算顺序与数值结果。Core 新入口仅数学值接口，不生成 Interface provenance。Interface 使用同一 old snapshot、r2 forms 与原完整 closed coverage，不在 tile 内制造新 interval/mesh。对未证明分支保持 pending/refinement；数学范围溢出应明确失败，不静默降级到采样规则。

验收：T0 <=5 production+2 tests、全39 runtime 重建、无新增warning，81 frozen/3 PDF/r3 archive身份不变。原未改测试数值输出逐字相同；改动的两测试与必要 Core include 计数单列。T1 精确 half-plane XY/RZ 两面向、正反 gross、wl/wg、非零 residual、跨叶共同 band、恒零/tangent、非法输入/FP环境；T2 原57主积分全通过，尤其原 shear 与旋转/非线性/薄层/AMR/weighted/反向速度及拒绝规则，真实四宽误差 E 不超过原预算。

失败规则：来源混用、漏 donor、以叶外延拓代替 actual phase、恒零当纯相、库存重置/clip、提高预算、解析参考逃出界或原 shear 不过均FAIL保留证据。只有真实外部依赖才BLOCKED。不接入 Hydro。

交付：实现、共享代数说明、T0/T1/T2 与失败记录、API/report/入口和 SHA archive。通过后仍按原 I12 的完整驱动/时间阶/AMR/曲率要求推进，不宣告 I12 已完成。

执行结果：**PASS（本项）**。70 独立子区域场景及原 57 主积分全部通过，uniform shear 为112 panels、E=5.7786364562e-6<=5.859375e-6。39 runtime 重建、25 CTest全通过，22项未改数值输出及原 affine phase 输出前缀不变，81 frozen/3PDF/r3 archive不变。见 [报告](I12a15r4_Certified_Phase_Subset_Integral_Report.md)。父 I12 / Gate I 仍 FAIL。
