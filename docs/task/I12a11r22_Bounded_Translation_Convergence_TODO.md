# I12a11r22 / 有界平移的非多项式 h/dt 验证

目标：只检验 r21 同一固定配置在原非多项式平移矩阵上的准确性，不修改数值算式或以参数调节获得通过。

允许修改：`tests/interface_hermite_gals_benchmark.cpp`（0 production + 1 test-related），只向原空间/时间控制入口传递 bounded 契约；本文、报告、README、必要配置说明和本项证据。

禁止修改：production 算子、r21 Hessian 预算 `min(hx,hy)^2`、默认 1e-10 几何容限、旧 gate、初态、mesh/dt/终点、Core 工作量或 roundoff 常数、冻结资料、rotation/vortex 路线及其它物理模块。

数学约束：原初始非多项式 `exp(q)-1`，同一 accepted first-jet/M 输运与半步/全步发布，XY/RZ、uniform/AMR；不重置解析几何或 M。三档网格 8/16/32、固定 dt=.01、终点 .04；曲率/梯度必须逐档递减且末档阶数 >=1.5/2。原界面/M/ledger/reflux 门槛不变。

实现要求：复用原 `translationSpaceStudy` 和 r13 时间控制编排，增加明确 bounded 命令，不复制数值驱动。记录原所有误差及 r21 实际预算/几何差/面修正/闭合诊断。空间 12 运行全部执行，即使部分失败也保留其它组证据。补充 XY uniform n=32、终点 .04 的 dt=.01/.005/.0025 控制；首点来自同矩阵，不重复。恒定平移的特征及区间积分本身精确，此时间控制用于检查跨步稳定性/重构累计，不代替原一般运动时间阶 >=1.8 的 Gate I 要求。

验收：T0 原 geometry 契约 smoke/平面跨面和已通过 bounded smoke 输出相同、architecture 通过、无 warning。T1/T2 12 组空间运行完成，满足原梯度/曲率空间 gate；两个新增小 dt 运行完成，并报告相对 dt=.01 的同终点误差，误差增大不得隐藏。完整原 I12 时间阶/rotation/vortex/dynamic AMR 未执行部分继续 FAIL。

失败规则：空间 gate 或实际推进失败即 FAIL；时间细化不能减少已接受步或导致运行拒绝。完整误差矩阵和失败阶段归档，再按具体原因另写任务，不在本项放宽容限或试参。缺明确 API 才 BLOCKED。

交付：最小入口接线、原/新误差比较、所有计划运行及失败日志、回归/冻结/hash 审计、报告与 archive。
