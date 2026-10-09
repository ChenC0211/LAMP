# BF-00 / 方法切换确认与冻结基线

**当前裁决：PASS（仅立法/冻结任务）。** 用户明确要求忽略旧档案丢失的前置、立即切入新路线，故该档案阻断已解除。全部原FAIL/缺失/timeout和原物理gate保持；下面BLOCKED为此前检查点，不作为继续实现的障碍。新任务为 [BF-01](BF-01_TODO.md)。

状态：**BLOCKED / 原payload丢失**。用户已明确答复原long-wave数据“丢失”。本轮只做BF-00，已冻结1,397数值文件与现存归档、复核物理RHS符号和25项基线检查；24项通过，完整momentum契约在原180秒门槛两次TIMEOUT，失败保留。本轮0数值改动，完整结果/待审查恢复方案见 [报告](BF-00_Report.md)。依据 [BF任务书](../roadmap/LAMP_BF_CSF_Codex_Tasks.md)、[BF-CSF V1方法](../roadmap/LAMP_BF_CSF_Method_V1.md)，原H13保持FAIL，BF-01/H14/H15不启动。

任务 ID / 标题：BF-00 / 方法切换确认与冻结基线。

目标：确定BF-CSF与原H12/H13的独立A/B关系，冻结当前实际源码、输入、原失败证据及全套门槛，核实物理毛细RHS符号。

允许修改：本TODO、BF-00报告、现行方法/状态索引（docs/README.md、00开发总纲、01数值硬约束、05路线图中的方法变更注记）和新BF任务入口。只读冻结/核验程序、归档、构建和检查输出置于 `runs/checks/BF-00/20261009_135037`；0 production、0注册测试。

禁止修改：src、include、Interface、tests、cases、CMake及旧H13报告/结果；不继续旧切向模态候选的开发，不改判H13，不进入BF-01/H14/H15。

数学约束：Pure FirstJet为唯一界面；native U/p、`C=-D_I^T`、sharp材料质量、唯一F_rho和真实R/uR、KEP/skew及mu=0保留。`c_eps`仅为从同一publication派生的力场，不能成为材料库存。BF-CSF实际RHS为 `+sigma*kappa*C(c_eps)`，旧balanceLoads必须取负号才进入RHS；凸滴液内压力高。

实现：读取完整新方法和现有源码；按文件SHA256标识无Git的当前checkout。分别冻结当前代码（含未采用的旧H13模态候选）、当前native毛细driver/input及原动态/静态payload，不将不同方法身份混为同一基线。保存受影响CTest注册、Release工具链、实际构建与结果。查找缺失历史payload，若无法恢复明确登记缺失而不推测数值。

验收：T0当前源码/输入身份与前后零数值diff、来源归档可独立复验；T1用实际C/载荷/saddle源码及现有符号断言确认物理RHS，记录现有H12/H13原报告、受影响CTest结果；T2冻结原long-wave FAIL、72/1728静滴、完整周期/空间时间/epsilon/密度/低黏/谱/能量/AMR/RZ gate及全部未达标项。报告不以短契约通过判BF-CSF物理PASS。

失败规则：任务书明确规定“任何引用文件缺失/版本不明记BLOCKED，不得以模型推测代替源码”。先搜索实际文件与现存归档；缺失原失败payload如无法恢复则BF-00 BLOCKED，不放宽此条件。goal不因首次遇到缺失立即标blocked；仍可完成其余冻结/检查工作。

交付：BF-00_Report.md、版本/文件/hash清单、冻结输入和现存原失败数据、回归日志、符号证据、缺失项/明确后续许可。BF-00全部门槛通过才允许下一轮BF-01。
