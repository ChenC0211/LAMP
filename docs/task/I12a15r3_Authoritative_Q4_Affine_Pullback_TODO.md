# I12a15r3 / 权威 native Q4 的仿射拉回包围

目标：只解决同一 old native Q4 多项式作用于 r2 相关 departure forms 后，phi(B) 的一次参数项和严格剩余带。使用 snapshot 的精确齐次数据；公开 rounded local DOFs 不得进入新证书。为下一项分块内液/气子区域积分提供数学输入。

允许修改：`PolynomialVelocitySlab.h/.cpp`（复用既有 strict affine engine）、`LevelSetHermiteGeometry.h/.cpp`（snapshot入口）、`Detail/LevelSetHermiteQuarticLeafData.h`（权威齐次数据适配），共5 production；新 `tests/interface_affine_phase_contract.cpp` 与 `tests/CMakeLists.txt`，共2 test-related；TODO/design/API/report/入口和证据。

禁止修改：实际速度/RK2、r2/a13/a15已有入口行为、原Q4 basis/constraints/root/sign/发布、原M/库存权重/reflux、benchmark参数、误差/panel/depth预算和任何gate。不添加第二几何或拟合多项式，不靠复制rounded缓存声称权威相。未接入相积分前，原剪切FAIL保持。

数学约束：原 native leaf 的25项齐次 DOF 经已有 exact Bernstein变换，再用精确有限差分转成同一张量Q4的幂系数；common positive denominator=144*原denominator。只做精确幂二尺度平衡、outward商包围，零系数仍精确零。与两个共享q_s/q_t的physical forms复合后，返回 c+s*q_s+t*q_t+R 包含该权威多项式的实际值。该证书属于指定叶的多项式；在叶外只表示同一多项式的数学延拓，不赋予物理donor/边界历史。调用者后续必须根据真实B所属叶使用它，并做完整closed覆盖。

实现要求：先写设计。Core既有affine算术增加固定degree4二维Horner与真实physical->cell chart的区间除法；不另写一套affine引擎。snapshot只取既有active native quartic权威数据；原state/source owner不变，无新增cache或重构。Core值原语不生成Interface capability。精确变换的系数若不可表示为有限区间则明确拒绝，不能删除小系数或改容差。

验收：T0 <=5 production+2 test；全runtime重建、无新增warning；81 frozen/3 PDF/r2 archive不变；原24项测试保持23 PASS+1原shear FAIL，完整原输出逐字相同。T1 新contract覆盖XY/RZ uniform/AMR、plane/Q2/Q4/薄层/内部液泡/相切、原polynomial实际RK2；每个B落入的closed native叶都按原snapshot核验，phi范围包含权威值和独立可精确表示的解析参考，保留hidden event而不按samples认证。空间/时间缩小时非线性残差二次缩小，平面不从roundedcache重建；跨叶/hanging及exact refinement来源正确。invalid/inactive/nonquartic/moved/source外输入与FP环境拒绝；调用不改geometry/M。后续相积分前仍需完整覆盖/来源契约。

失败规则：参考逃出证书、利用诊断cache取代精确scalar、未授权physical延拓、旧行为改变、预算/舍入规则放宽，均FAIL并保留。仅真正外部依赖才BLOCKED。新primitive通过不能重分类原剪切或I12状态。

交付：源代码、代数证明、T0/T1/失败记录、API/report与archive。之后单独任务用phi的一次项及严格余带确定分块液/气子区域，按原gross XY/RZ积分预算复验同一shear。

验收记账澄清：新 Q4 适配需要在原 QuarticLeafData header 显式包含 Core PolynomialVelocitySlab.h，因此 architecture 的 include 计数预期增加 1。原 23 项数值输出（含旧 shear FAIL）仍逐字比较；architecture 则要求原 targets/TU/files/5 legacy edges 及依赖规则保持，仅该 Core 内部 include 可变。报告单列这一变化，不宣称 architecture 文本逐字相同，不放宽任何数值或依赖 gate。

执行结果：**PASS（本项）**。232 区域、33,498 次实际回溯及 57,470 次闭叶检查通过，最小 residual 缩小比 4.0207；39 runtime 目标重建、25 CTest 为 24 PASS+1 原 shear FAIL。23 项旧数值输出逐字相同，81 frozen/3 PDF/r2 archive 不变。见 [报告](I12a15r3_Authoritative_Q4_Affine_Pullback_Report.md)；I12 / Gate I 仍 FAIL，继续下一项分块内相积分证书。
