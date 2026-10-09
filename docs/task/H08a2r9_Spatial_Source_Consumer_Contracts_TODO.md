# H08a2r9 / 空间来源消费者契约迁移

状态：**PASS（恢复后的消费者契约）**，见[当前报告](H08a2r9_Spatial_Source_Consumer_Contracts_Resume_Report.md)与[恢复设计](H08a2r9_Spatial_Source_Consumer_Contracts_Resume_Design.md)。原[精度前置FAIL](H08a2r9_Spatial_Source_Consumer_Contracts_Report.md)由H08a2r10关闭并保持历史；本次完成下述production 0/test 2迁移。下一唯一任务[H08a2r11](H08a2r11_Common_Regular_Source_Deployment_TODO.md)同步部署实际来源/运动证书；当前B/J/atom仍旧。

目标：使现有 phase/momentum 空间契约验收同源物理量与准确性，而非强制已经授权放弃的旧重构公式、点态 cell 常散度或 pressure→dual 交换；为随后原子式切换实际来源准备正确验收。

允许修改：production 0；test 2：`cch_phase_geometry_contract.cpp`、`cch_momentum_geometry_contract.cpp`。实施前逐条登记旧来源专属断言、对应仍保留的物理要求及独立替代检查。无新长期脚本/入口，不增加生产方法开关。

禁止修改：生产 B、regular 准备/查询、characteristic、FirstJet、F/R/uR/KEP、benchmark、原最终 native/affine/H10b-v2 门槛及历史证据。不得简单移除或放宽失败断言；不得用当前待测值构造一个恒等于自身的参考。

数学约束：材料仍仅唯一共享 F、内部反号/AMR 子面和、来源一致、KEP/skew 四硬约束。空间准确性仍须独立检查；geometry/advector/full clock、轴/wall、真实 native 物理面矩及 Gauss 恒等式保留。点态散度可变化；dual 通量不必等于 pressure cell 平均散度按 overlap 体积配分，但应等于同一实际场在该真实区域的物理散度积分。

实现要求：将 phase 的点态 B-div=原 DI/V 改为足阶分片体积分核验；momentum 中旧 `AdvectorReference` trace/lift 重构专属对照与 dual 常散度配分，替换为当前实际查询/J、独立物理求积及解析族的闭式参考。积分分片来自真实 owner 的 cell/piece/atom 支撑，不从 rounded 极短端点差重建宽度。原 affine/允许 RZ 解析族、零输入、tiny、符号、边界、内部共享及错误来源拒绝必须保持或增强。旧缺信息统计保留历史，不将旧一次 tangent rank/存储常数当成未来 P2/Q4 场的硬门槛；实际信息不足必须继续显式报告并限制精确性声明。

验收：T0 标准构建、架构/目录及受影响契约，生产和未涉及源码哈希保持。T1 替代物理积分参考具有独立闭式/足阶依据，原全部布局/输入/边界/tiny 覆盖不缩小；当前接受场通过，故意错误的 normal 符号、遗漏高项/支撑或错误来源可被对应检查识别。T2 H08a2r8 的已知非恒定散度正则场在同一数学验收下亦通过，且旧点态/dual 配分差仍被如实诊断；不能靠算法枚举选择宽松判据。明确区分同源性、精度和 PDE gate，不以此任务替代实际来源切换后的完整验证。

失败规则：任何独立性/覆盖/精度丢失、仅删除旧断言、把参考回灌生产即 FAIL。无法在当前事实内完成则保留具体剩余项，不私自越过 H08a；仅真实外部依赖可 BLOCKED。

交付：两份原消费者契约、断言迁移表、独立数值证据与生产未改证明。随后在一次来源部署中同步迁移实际 owner、piece/atom、完整次数 C1 运动证书和相应 characteristic/benchmark 消费者，按预算执行，再重验 H10o/H10b-v2 并返回 H08a3/H08a；H11 未启动。
